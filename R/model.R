# Model construction: compile a Stan program, bind data, and expose the
# log density, its gradient and the parameter transforms.

#' Check a Stan program
#'
#' Parses and type checks a Stan program without binding data.
#'
#' @param code Stan program as a character string (or vector of lines).
#' @param file Path to a `.stan` file, used when `code` is `NULL`.
#' @param include_paths Directories searched, in order, for files named in
#'   `#include` directives. The directory of `file` is searched last.
#' @return The typed intermediate representation, invisibly. Errors carry the
#'   classes `bayescomp_syntax_error`, `bayescomp_semantic_error` or
#'   `bayescomp_unsupported_error`.
#' @examples
#' bayescomp_check("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' @export
bayescomp_check <- function(code = NULL, file = NULL, include_paths = NULL) {
  src <- read_stan_source(code, file, include_paths)
  invisible(stan_check(stan_parse(src)))
}

read_stan_source <- function(code, file, include_paths = NULL) {
  if (is.null(code)) {
    if (is.null(file)) stop("supply either 'code' or 'file'", call. = FALSE)
    code <- readLines(file, warn = FALSE)
    include_paths <- c(include_paths, dirname(file))
  }
  if (!is.character(code)) stop("'code' must be a character string", call. = FALSE)
  expand_includes(unlist(strsplit(paste(code, collapse = "\n"), "\n", fixed = TRUE)), include_paths, character())
}

# replace each #include line by the named file, searched in include_paths
expand_includes <- function(lines, paths, seen) {
  re <- "^[[:space:]]*#include[[:space:]]*[<\"]?([^>\"[:space:]]+)[>\"]?[[:space:]]*$"
  hit <- grep(re, lines)
  if (!length(hit)) return(paste(lines, collapse = "\n"))
  out <- vector("list", length(lines))
  out[] <- as.list(lines)
  for (i in hit) {
    name <- sub(re, "\\1", lines[i])
    cand <- file.path(paths, name)
    f <- cand[file.exists(cand)][1L]
    if (is.na(f)) stop_syntax(sprintf("#include file '%s' not found in the include paths (%s)", name,
                                      if (length(paths)) paste(paths, collapse = ", ") else "none given; use 'include_paths'"),
                              list(line = i, col = 1L))
    f <- normalizePath(f)
    if (f %in% seen) stop_syntax(sprintf("#include file '%s' includes itself", name), list(line = i, col = 1L))
    out[[i]] <- expand_includes(readLines(f, warn = FALSE), c(paths, dirname(f)), c(seen, f))
  }
  paste(unlist(out), collapse = "\n")
}

#' Compile a Stan program and bind data
#'
#' @param code Stan program as a character string.
#' @param data Named list of data, one element per variable in the `data`
#'   block. Arrays use R arrays with the array dimensions first, so
#'   `array[N] vector[K] y` is an `N` by `K` matrix.
#' @param file Path to a `.stan` file, used when `code` is `NULL`.
#' @param include_paths Directories searched for `#include` files, as in
#'   [bayescomp_check()].
#' @param seed Seed for random number generation in `transformed data`.
#' @param backend Execution backend, see [bayescomp_backends()]. Native
#'   backends compute the log density and its gradient.
#' @param functions Named list of R functions implementing functions that
#'   the program declares without a body (Stan's external C++ functions,
#'   such as `vector CODOLS(matrix X, vector y);` in rstanarm). Each is called
#'   with the Stan arguments; arguments that depend on parameters arrive as
#'   automatic differentiation values, so functions applied to parameters
#'   must be built from operations that support them.
#' @param allow_fallback When the backend lacks a capability the model needs
#'   (constraining transforms, generated quantities or random numbers), use
#'   the R implementation for it. In this version only the log density and
#'   its gradient are ever taken from a native backend, so a native backend
#'   always needs `allow_fallback = TRUE` for sampling. The default `FALSE` raises an error
#'   instead, so implementations are never mixed silently. The model's
#'   `provenance` element records which backend provides each capability.
#'
#'   As in rstan, a data variable `int CHAIN_ID` that is declared but not
#'   supplied is set to the chain number: 1 for the model itself, and
#'   1, 2, ... in each chain of the samplers.
#' @return An object of class `bayescomp_model`.
#' @examples
#' m <- bayescomp_model("
#'   data { int N; vector[N] y; }
#'   parameters { real mu; real<lower=0> sigma; }
#'   model { y ~ normal(mu, sigma); }
#' ", data = list(N = 3, y = c(1, 2, 3)))
#' m
#' bayescomp_log_prob(m, c(0, 0))
#' @export
bayescomp_model <- function(code = NULL, data = list(), file = NULL, seed = NULL, backend = "r",
                            allow_fallback = FALSE, functions = list(), include_paths = NULL) {
  src <- read_stan_source(code, file, include_paths)
  if (!is.list(functions) || (length(functions) && (is.null(names(functions)) || !all(vapply(functions, is.function, TRUE)))))
    stop("'functions' must be a named list of R functions", call. = FALSE)
  ir <- stan_check(stan_parse(src), external = names(functions))
  if (is.null(data)) data <- list()
  if (!is.list(data)) stop("'data' must be a named list", call. = FALSE)
  be <- backend_get(backend)
  m <- compile_r(ir)
  m$external_impl <- functions
  m$chain_id_auto <- "CHAIN_ID" %in% vapply(top_decls(ir, "data"), `[[`, "", "name") && is.null(data$CHAIN_ID)
  if (m$chain_id_auto) data$CHAIN_ID <- 1L
  m$seed <- seed
  inst <- with_seed(seed, instantiate(m, data))
  caps <- c("log_density", "gradient", "constrain", "generated_quantities", "rng")
  # this version takes only the log density and its gradient from a native
  # backend; the other capabilities always run in R, whatever is declared
  used <- intersect(be$capabilities, native_capabilities)
  prov <- stats::setNames(as.list(ifelse(caps %in% used, be$name, "r")), caps)
  out <- structure(list(code = src, ir = ir, compiled = m, env = inst$env, decls = inst$decls,
                        params = inst$params, dim = inst$dim, data = data, backend = be$name,
                        provenance = prov, native = NULL),
                   class = "bayescomp_model")
  if (be$name == "r") return(out)
  core <- setdiff(c("log_density", "gradient"), be$capabilities)
  if (length(core))
    stop(sprintf("backend '%s' does not provide %s", be$name, paste(core, collapse = " and ")), call. = FALSE)
  other <- setdiff(caps, used)
  if (length(other) && !isTRUE(allow_fallback))
    stop(sprintf("backend '%s' does not provide %s; set allow_fallback = TRUE to use the R implementation for %s",
                 be$name, paste(other, collapse = ", "), if (length(other) > 1L) "them" else "it"), call. = FALSE)
  native <- be$compile(bayescomp_ir(out), out)
  if (!is.list(native) || !is.function(native$log_density_gradient))
    stop(sprintf("backend '%s' did not return a log_density_gradient function", be$name), call. = FALSE)
  out$native <- native
  out
}

# log density and gradient on the unconstrained scale from the model's
# backend: list(value, grad). Errors propagate with their classes.
model_lpg <- function(model, jacobian, propto, print = FALSE) {
  nat <- model$native
  if (!is.null(nat)) return(function(u) nat$log_density_gradient(u, jacobian, propto))
  f <- model_fun(model, "log_prob")
  function(u) {
    ctx <- new_ctx(jacobian = jacobian, print = print, propto = propto)
    ad_value_grad(function(x) f(x, ctx), u)
  }
}

#' Check a native backend against the R implementation
#'
#' Evaluates the log density and gradient of a model compiled with a native
#' backend and of the same model in R at random unconstrained points.
#'
#' @param model A `bayescomp_model` with a native backend.
#' @param n Number of points.
#' @param radius Points are uniform on `(-radius, radius)`.
#' @param seed Random seed.
#' @return A data frame with the absolute and relative differences of the
#'   log density and the largest ones of the gradient at each point, for
#'   each combination of `jacobian` and `propto`.
#' @export
bayescomp_validate_backend <- function(model, n = 20L, radius = 2, seed = 1L) {
  check_model(model)
  if (is.null(model$native)) stop("the model uses the R backend; nothing to compare", call. = FALSE)
  ref <- model
  ref$native <- NULL
  U <- with_seed(seed, matrix(stats::runif(n * model$dim, -radius, radius), n, model$dim))
  combos <- expand.grid(jacobian = c(TRUE, FALSE), propto = c(FALSE, TRUE))
  rows <- lapply(seq_len(nrow(combos)), function(k) {
    jac <- combos$jacobian[k]; pro <- combos$propto[k]
    fa <- model_lpg(model, jac, pro); fb <- model_lpg(ref, jac, pro)
    do.call(rbind, lapply(seq_len(n), function(i) {
      a <- fa(U[i, ]); b <- fb(U[i, ])
      data.frame(point = i, jacobian = jac, propto = pro, lp_abs = abs(a$value - b$value),
                 lp_rel = abs(a$value - b$value) / max(1, abs(b$value)),
                 grad_abs = max(abs(a$grad - b$grad), 0),
                 grad_rel = max(abs(a$grad - b$grad) / pmax(1, abs(b$grad)), 0))
    }))
  })
  do.call(rbind, rows)
}

with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
}

# ---- compilation ---------------------------------------------------------------

top_decls <- function(ir, block) {
  b <- ir$blocks[[block]]
  if (is.null(b)) return(list())
  Filter(function(st) st$kind == "decl" && isTRUE(st$top), b$body)
}
block_body <- function(ir, block) ir$blocks[[block]]$body %||% list()

check_calls <- function(decls, block) {
  lapply(decls, function(d) {
    cl(rt_check, vsym(d$name), cl(as.name("[["), quote(.decls), d$name), cg_cons_args(d$type), block)
  })
}

compile_r <- function(ir) {
  ctxl <- function(b) list(block = b)
  funs <- lapply(Filter(function(fd) !is.null(fd$body), ir$blocks$functions$body %||% list()), cg_fundef)
  # transformed data: statements, then validation of the declared variables
  td <- c(cg_block(block_body(ir, "transformed data"), ctxl("transformed data")), list(quote(.stmt$loc <- NULL)))
  td_checks <- check_calls(top_decls(ir, "transformed data"), "transformed data")
  pars <- top_decls(ir, "parameters")
  par_code <- lapply(pars, function(d) {
    cl(as.name("<-"), vsym(d$name), cl(rt_param, quote(.u), cl(as.name("[["), quote(.decls), d$name),
                                        cg_cons_args(d$type), CTX))
  })
  tp <- c(cg_block(block_body(ir, "transformed parameters"), ctxl("transformed parameters")),
          list(quote(.stmt$loc <- NULL)),
          check_calls(top_decls(ir, "transformed parameters"), "transformed parameters"))
  mb <- cg_block(block_body(ir, "model"), ctxl("model"))
  gq <- c(cg_block(block_body(ir, "generated quantities"), ctxl("generated quantities")),
          list(quote(.stmt$loc <- NULL)),
          check_calls(top_decls(ir, "generated quantities"), "generated quantities"))
  out_list <- function(decls) cl_list(as.name("list"), stats::setNames(lapply(decls, function(d) vsym(d$name)), vapply(decls, `[[`, "", "name")))
  tps <- top_decls(ir, "transformed parameters")
  gqs <- top_decls(ir, "generated quantities")
  # parameter transforms run before any statement: no location
  reset <- list(quote(.stmt$loc <- NULL))
  log_prob_body <- blk(c(reset, par_code, tp, mb, list(cl(rt_target_get, CTX))))
  write_body <- blk(c(reset, par_code,
    list(cl(as.name("<-"), quote(.out), cl(as.name("list"), out_list(pars), NULL, NULL))),
    list(cl(as.name("if"), quote(.tp || .gq), blk(c(tp, list(cl(as.name("<-"), quote(.out[[2L]]), out_list(tps)))))),
         cl(as.name("if"), quote(.gq), blk(c(gq, list(cl(as.name("<-"), quote(.out[[3L]]), out_list(gqs)))))),
         quote(.out))))
  unc_code <- lapply(seq_along(pars), function(i) {
    d <- pars[[i]]
    blk(list(
      cl(as.name("<-"), vsym(d$name), cl(rt_read_value, quote(.vals), cl(as.name("[["), quote(.decls), d$name))),
      cl(as.name("<-"), cl(as.name("[["), quote(.parts), i),
         cl(unconstrain_variable, cl(as.name("val"), vsym(d$name)), cl(as.name("[["), quote(.decls), d$name),
            cl(as.name("lapply"), cg_cons_args(d$type), as.name("val"))))))
  })
  unc_body <- blk(c(list(quote(.parts <- vector("list", 0L))), unc_code, list(quote(unlist(.parts)))))
  list(ir = ir, funs = funs, tdata = td, tdata_checks = td_checks, log_prob = log_prob_body, write = write_body,
       unconstrain = unc_body, pars = pars, tps = tps, gqs = gqs,
       data_decls = top_decls(ir, "data"), td_decls = top_decls(ir, "transformed data"))
}

decl_info <- function(d, env) {
  t <- d$type
  sh <- cg_decl_shape(t)
  adims <- as.integer(eval(sh$adims, env))
  eshape <- as.integer(eval(sh$eshape, env))
  if (any(c(adims, eshape) < 0L)) stop_data(sprintf("variable '%s' has a negative size (%s)", d$name, paste(c(adims, eshape), collapse = ", ")))
  info <- list(name = d$name, kind = decl_kind(t), base = decl_base_type(t$base), stan_type = t$base,
               adims = adims, eshape = eshape)
  if (info$kind == "cholesky_factor_cov" && length(eshape) == 2L && eshape[1L] < eshape[2L])
    stop_data(sprintf("cholesky_factor_cov '%s' needs at least as many rows as columns", d$name))
  info$free_size <- free_size_variable(info)
  info
}

instantiate <- function(m, data) {
  ns <- environment(rt_decl)
  fenv <- new.env(parent = ns)
  for (f in m$funs) assign(f$name, make_function(f$formals, f$body, fenv), envir = fenv)
  for (x in m$ir$externals) local({
    impl <- m$external_impl[[x$name]]
    # the generated code passes the evaluation context last
    assign(x$mangled, function(...) { a <- list(...); do.call(impl, a[-length(a)]) }, envir = fenv)
  })
  env <- new.env(parent = fenv)
  decls <- list()
  env$.decls <- decls
  env$.ctx <- new_ctx(jacobian = FALSE)
  # data
  for (d in m$data_decls) {
    info <- decl_info(d, env)
    decls[[d$name]] <- info
    v <- read_data_value(data, info)
    assign(paste0("V.", d$name), v, envir = env)
    env$.decls <- decls
    args <- lapply(d$type$cons, function(e) eval(cg_expr(e), env))
    msg <- check_variable(v, info, args)
    if (!is.null(msg)) stop_data(sprintf("variable '%s' %s", d$name, msg))
  }
  extra <- setdiff(names(data), vapply(m$data_decls, `[[`, "", "name"))
  # transformed data
  # sizes of transformed data variables may depend on earlier transformed
  # data, so their metadata is built after the statements run
  if (length(m$tdata)) eval(blk(m$tdata), env)
  for (d in m$td_decls) decls[[d$name]] <- decl_info(d, env)
  env$.decls <- decls
  if (length(m$tdata_checks)) eval(blk(m$tdata_checks), env)
  # parameters, transformed parameters, generated quantities
  off <- 0L
  for (d in m$pars) {
    info <- decl_info(d, env)
    info$offset <- off
    off <- off + info$free_size
    decls[[d$name]] <- info
  }
  for (d in c(m$tps, m$gqs)) decls[[d$name]] <- decl_info(d, env)
  env$.decls <- decls
  env$.ignored_data <- extra
  list(env = env, decls = decls, dim = off,
       params = list(pars = vapply(m$pars, `[[`, "", "name"), tps = vapply(m$tps, `[[`, "", "name"),
                     gqs = vapply(m$gqs, `[[`, "", "name")))
}

read_data_value <- function(data, info) {
  nm <- info$name
  D <- c(info$adims, info$eshape)
  if (!(nm %in% names(data))) {
    if (length(D) && prod(D) == 0L) return(rt_decl(info$base, info$adims, info$eshape))
    stop_data(sprintf("variable '%s' is declared in the data block but was not supplied", nm))
  }
  x <- data[[nm]]
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.list(x)) {
    # arrays given as lists of elements
    x <- tryCatch(simplify_data_list(x), error = function(e) stop_data(sprintf("variable '%s': %s", nm, conditionMessage(e))))
  }
  if (is.logical(x)) x <- as.integer(x)
  if (!is.numeric(x)) stop_data(sprintf("variable '%s' must be numeric", nm))
  if (anyNA(x) && !all(is.nan(x[is.na(x)]))) stop_data(sprintf("variable '%s' contains NA values", nm))
  xd <- if (is.null(dim(x))) length(x) else dim(x)
  ok <- if (length(D) == 0L) length(x) == 1L
        else if (length(D) == 1L) length(x) == D && (is.null(dim(x)) || sum(dim(x) > 1L) <= 1L)
        else identical(as.integer(xd), as.integer(D)) || (prod(D) == 0L && length(x) == 0L)
  if (!ok) stop_data(sprintf("variable '%s' has size %s but the program declares %s", nm,
                             paste(xd, collapse = " x "), if (length(D)) paste(D, collapse = " x ") else "a scalar"))
  if (info$base == "int") {
    if (any(is.nan(x)) || any(abs(x) > .Machine$integer.max) || any(x != round(x)))
      stop_data(sprintf("variable '%s' is declared int but contains non integer values", nm))
    x <- as.integer(x)
  } else {
    x <- as.double(x)
  }
  dim(x) <- if (length(D) > 1L) D else NULL
  x
}

simplify_data_list <- function(x) {
  if (!length(x)) return(numeric(0))
  parts <- lapply(x, function(e) if (is.list(e)) simplify_data_list(e) else e)
  d1 <- if (is.null(dim(parts[[1L]]))) length(parts[[1L]]) else dim(parts[[1L]])
  n <- length(parts)
  flat <- unlist(lapply(parts, as.vector))
  if (length(parts[[1L]]) == 1L) return(flat)
  arr <- array(flat, dim = c(d1, n))
  aperm(arr, c(length(d1) + 1L, seq_along(d1)))
}

# read a constrained value supplied by the user for unconstraining
rt_read_value <- function(vals, info) {
  nm <- info$name
  if (!(nm %in% names(vals))) stop_data(sprintf("value for parameter '%s' is missing", nm))
  v <- read_data_value(vals, utils::modifyList(info, list(base = "real")))
  v
}

# ---- evaluation ----------------------------------------------------------------

check_upars <- function(model, upars) {
  if (!is.numeric(upars) || length(upars) != model$dim)
    stop(sprintf("expected %d unconstrained parameter values, got %d", model$dim, length(upars)), call. = FALSE)
  as.numeric(upars)
}

model_fun <- function(model, which) {
  key <- paste0(".fn_", which)
  f <- model$env[[key]]
  if (is.null(f)) {
    body <- model$compiled[[which]]
    fm <- switch(which,
      log_prob = alist(.u = , .ctx = ),
      write = alist(.u = , .ctx = , .tp = , .gq = ),
      unconstrain = alist(.vals = ))
    f <- make_function(fm, body, model$env)
    assign(key, f, envir = model$env)
  }
  f
}

# value and gradient of the log density at unconstrained values
lp_grad <- function(model, upars, jacobian = TRUE, propto = TRUE, print = TRUE) {
  model_lpg(model, jacobian, propto, print)(upars)
}

#' Log density of a model
#'
#' @param model A `bayescomp_model`.
#' @param upars Numeric vector of unconstrained parameter values.
#' @param jacobian Include the log absolute Jacobian determinant of the
#'   constraining transforms.
#' @param propto Drop constant terms of sampling statements (`~`) and
#'   `_lupdf`/`_lupmf` calls, as the samplers do.
#' @return The log density (a number). `bayescomp_grad_log_prob()` returns the
#'   gradient with the log density in attribute `"log_prob"`.
#' @export
bayescomp_log_prob <- function(model, upars, jacobian = TRUE, propto = FALSE) {
  lp_grad(model, check_upars(model, upars), jacobian, propto)$value
}

#' @rdname bayescomp_log_prob
#' @export
bayescomp_grad_log_prob <- function(model, upars, jacobian = TRUE, propto = FALSE) {
  r <- lp_grad(model, check_upars(model, upars), jacobian, propto)
  structure(r$grad, log_prob = r$value)
}

#' Map between constrained and unconstrained parameters
#'
#' @param model A `bayescomp_model`.
#' @param pars Named list of constrained parameter values.
#' @param upars Numeric vector of unconstrained values.
#' @param include_tparams,include_gqs Also return transformed parameters and
#'   generated quantities.
#' @param seed Seed for random number generators in generated quantities.
#' @return `bayescomp_unconstrain()` returns a numeric vector;
#'   `bayescomp_constrain()` returns a named list of arrays.
#' @export
bayescomp_unconstrain <- function(model, pars) {
  if (!is.list(pars)) stop("'pars' must be a named list", call. = FALSE)
  f <- model_fun(model, "unconstrain")
  old <- .ad$tape
  .ad$tape <- NULL
  on.exit(.ad$tape <- old)
  out <- f(pars)
  if (is.null(out)) numeric(0) else as.numeric(out)
}

#' @rdname bayescomp_unconstrain
#' @export
bayescomp_constrain <- function(model, upars, include_tparams = TRUE, include_gqs = TRUE, seed = NULL) {
  upars <- check_upars(model, upars)
  r <- with_seed(seed, write_array(model, upars, include_tparams, include_gqs))
  c(r[[1L]], r[[2L]], r[[3L]])
}

write_array <- function(model, upars, tp = TRUE, gq = TRUE, print = TRUE) {
  f <- model_fun(model, "write")
  old <- .ad$tape
  .ad$tape <- NULL
  on.exit(.ad$tape <- old)
  ctx <- new_ctx(jacobian = FALSE, print = print)
  f(upars, ctx, tp, gq)
}

# flat names in Stan's column major order: a, b[1], m[1,2]
flat_names <- function(name, dims) {
  if (!length(dims)) return(name)
  if (prod(dims) == 0L) return(character(0))
  g <- do.call(expand.grid, lapply(dims, seq_len))
  paste0(name, "[", do.call(paste, c(unname(as.list(g)), sep = ",")), "]")
}

#' Parameter names of a model
#'
#' @param model A `bayescomp_model`.
#' @param include_tparams,include_gqs Include transformed parameters and
#'   generated quantities.
#' @return Character vector of flattened names in Stan's output order.
#' @export
bayescomp_param_names <- function(model, include_tparams = TRUE, include_gqs = TRUE) {
  nm <- c(model$params$pars, if (include_tparams) model$params$tps, if (include_gqs) model$params$gqs)
  unlist(lapply(nm, function(n) {
    d <- model$decls[[n]]
    flat_names(n, c(d$adims, d$eshape))
  }))
}

#' Number of unconstrained parameters
#' @param model A `bayescomp_model`.
#' @return The number of unconstrained parameters, an integer.
#' @export
bayescomp_dim <- function(model) model$dim

# flatten a write_array result to a named numeric vector
flatten_draw <- function(lst) {
  unlist(lapply(lst, function(v) as.numeric(v)), use.names = FALSE)
}

#' @export
print.bayescomp_model <- function(x, ...) {
  cat("<bayescomp_model>\n")
  cat("  unconstrained dimension:", x$dim, "\n")
  show <- function(lbl, nms) if (length(nms)) cat(sprintf("  %-24s %s\n", lbl, paste(nms, collapse = ", ")))
  show("parameters:", x$params$pars)
  show("transformed parameters:", x$params$tps)
  show("generated quantities:", x$params$gqs)
  cat("  backend:", x$backend, "\n")
  if (length(x$env$.ignored_data)) cat("  unused data:", paste(x$env$.ignored_data, collapse = ", "), "\n")
  invisible(x)
}
