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
#' @param pedantic Print the warnings of [bayescomp_pedantic()] (stanc's
#'   `--warn-pedantic`).
#' @param warn_uninitialized Print warnings about variables used before they
#'   are assigned (stanc's `--warn-uninitialized`).
#' @return The typed intermediate representation, invisibly. Errors carry the
#'   classes `bayescomp_syntax_error`, `bayescomp_semantic_error` or
#'   `bayescomp_unsupported_error`.
#' @examples
#' bayescomp_check("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' @export
bayescomp_check <- function(code = NULL, file = NULL, include_paths = NULL, pedantic = FALSE, warn_uninitialized = FALSE) {
  src <- read_stan_source(code, file, include_paths)
  ir <- stan_check(stan_parse(src))
  if (isTRUE(pedantic)) print(ped_result(c(ped_checks(ir), ped_uninit(ir))))
  else if (isTRUE(warn_uninitialized)) print(ped_result(ped_uninit(ir)))
  invisible(ir)
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

# TRUE for lines that start inside a /* */ comment
in_block_comment <- function(lines) {
  inside <- FALSE
  out <- logical(length(lines))
  for (i in seq_along(lines)) {
    out[i] <- inside
    m <- regmatches(lines[i], gregexpr('/\\*|\\*/|//|"', lines[i]))[[1L]]
    str <- FALSE
    for (t in m) {
      if (inside) { if (t == "*/") inside <- FALSE }
      else if (str) { if (t == '"') str <- FALSE }
      else if (t == '"') str <- TRUE
      else if (t == "//") break
      else if (t == "/*") inside <- TRUE
    }
  }
  out
}

# replace each #include line by the named file, searched in include_paths.
# With includes, attribute "linemap" gives the file (NA for the main
# program) and line of every line of the result, for error messages.
expand_includes <- function(lines, paths, seen, file = NA_character_) {
  re <- "^[[:space:]]*#include[[:space:]]*[<\"]?([^>\"[:space:]]+)[>\"]?[[:space:]]*$"
  hit <- grep(re, lines)
  # stanc ignores directives inside block comments
  hit <- hit[!in_block_comment(lines)[hit]]
  if (!length(hit) && is.na(file)) return(paste(lines, collapse = "\n"))
  out <- vector("list", length(lines))
  out[] <- as.list(lines)
  map <- lapply(seq_along(lines), function(i) list(file = rep(file, 1L), line = i))
  inc <- character()
  for (i in hit) {
    name <- sub(re, "\\1", lines[i])
    cand <- file.path(paths, name)
    f <- cand[file.exists(cand)][1L]
    if (is.na(f)) stop_syntax(sprintf("#include file '%s' not found in the include paths (%s)", name,
                                      if (length(paths)) paste(paths, collapse = ", ") else "none given; use 'include_paths'"),
                              list(line = i, col = 1L))
    f <- normalizePath(f)
    if (f %in% seen) stop_syntax(sprintf("#include file '%s' includes itself", name), list(line = i, col = 1L))
    sub_src <- expand_includes(readLines(f, warn = FALSE), c(paths, dirname(f)), c(seen, f), name)
    inc <- c(inc, f, attr(sub_src, "included"))
    out[[i]] <- as.vector(sub_src)
    map[[i]] <- attr(sub_src, "linemap")
  }
  txt <- paste(unlist(out), collapse = "\n")
  structure(txt, included = unique(inc),
            linemap = list(file = unlist(lapply(map, `[[`, "file")), line = unlist(lapply(map, `[[`, "line"))))
}

#' Compile a Stan program and bind data
#'
#' @param code Stan program as a character string.
#' @param data Named list of data, one element per variable in the `data`
#'   block, or the name of a CmdStan JSON or R dump file (see
#'   [bayescomp_read_data()]). Arrays use R arrays with the array dimensions first, so
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
  if (is.null(functions)) functions <- list()
  if (!is.list(functions) || (length(functions) && (is.null(names(functions)) || !all(vapply(functions, is.function, TRUE)))))
    stop("'functions' must be a named list of R functions", call. = FALSE)
  ast <- stan_parse(src)
  ir <- stan_check(ast, external = names(functions))
  for (w in ped_int_division(ir))
    warning(sprintf("line %s, column %s: %s", w$loc$line %||% "?", w$loc$col %||% "?", w$msg), call. = FALSE)
  # loops that only touch their own element run as statements on slices
  if (vectorize_enabled()) ir <- vectorize_program(ast, ir, external = names(functions))
  if (is.null(data)) data <- list()
  data <- data_arg(data)
  if (!is.list(data)) stop("'data' must be a named list or the name of a JSON or R dump file", call. = FALSE)
  be <- backend_get(backend)
  m <- compile_r(ir)
  m$external_impl <- functions
  m$chain_id_auto <- "CHAIN_ID" %in% vapply(top_decls(ir, "data"), `[[`, "", "name") && is.null(data$CHAIN_ID)
  if (m$chain_id_auto) data$CHAIN_ID <- 1L
  # transformed data rebuilt per chain for CHAIN_ID must draw the same random
  # numbers in every chain, so such a model always has a seed
  if (m$chain_id_auto && is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
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
# value only log density: the same code without a tape
model_lp <- function(model, jacobian, propto, print = FALSE) {
  nat <- model$native
  if (!is.null(nat)) {
    if (!is.null(nat[["log_density"]])) return(function(u) nat[["log_density"]](u, jacobian, propto))
    return(function(u) nat$log_density_gradient(u, jacobian, propto)$value)
  }
  f <- model_fun(model, "log_prob")
  # with propto, which terms are kept depends on which arguments are
  # parameters, so the tape still records (only the backward sweep is
  # skipped); without it the code runs with no tape at all
  if (propto) return(function(u) {
    ctx <- new_ctx(jacobian = jacobian, print = print, propto = TRUE)
    old <- tape_start()
    on.exit(tape_restore(old))
    as.numeric(val(f(ad_leaf(as.numeric(u)), ctx)))
  })
  function(u) {
    ctx <- new_ctx(jacobian = jacobian, print = print, propto = FALSE)
    old <- .ad$tape
    .ad$tape <- NULL
    on.exit(.ad$tape <- old)
    as.numeric(val(f(as.numeric(u), ctx)))
  }
}

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
#' @examples
#' # a backend that evaluates the standard normal density in closed form
#' bayescomp_register_backend("example_normal", c("log_density", "gradient"), ir_versions = "1.0.0",
#'   compile = function(ir, model) list(log_density_gradient = function(u, jacobian, propto)
#'     list(value = -0.5 * sum(u^2) - if (propto) 0 else 0.5 * length(u) * log(2 * pi), grad = -u)),
#'   description = "closed form normal density")
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }",
#'                      backend = "example_normal", allow_fallback = TRUE)
#' bayescomp_validate_backend(m, n = 5)
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

# a second seed derived from a seed, valid for every integer seed
next_seed <- function(seed) if (is.null(seed)) NULL else as.integer(as.numeric(seed) %% (.Machine$integer.max - 1) + 1)

with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  old_kind <- RNGkind()
  on.exit({
    RNGkind(old_kind[1L], old_kind[2L], old_kind[3L])
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  })
  # the same stream whatever generator the caller (or a chain) has set
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
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
  h <- hoist_ir(ir)
  hoist <- lapply(h$hoisted, function(x) list(name = x$name, code = cg_expr(x$expr)))
  ir_orig <- ir
  ir <- h$ir
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
  out_list <- function(decls) cl_list(as.name("list"), stats::setNames(lapply(decls, function(d)
    if (is_cx_base(d$type$base)) cl(cx_out, vsym(d$name))
    else if (identical(d$type$base, "tuple")) cl(tuple_out, vsym(d$name), cl(as.name("[["), quote(.decls), d$name))
    else vsym(d$name)), vapply(decls, `[[`, "", "name")))
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
  list(ir = ir_orig, hoist = hoist, funs = funs, tdata = td, tdata_checks = td_checks, log_prob = log_prob_body, write = write_body,
       unconstrain = unc_body, pars = pars, tps = tps, gqs = gqs,
       data_decls = top_decls(ir, "data"), td_decls = top_decls(ir, "transformed data"))
}

decl_info <- function(d, env) {
  t <- d$type
  if (identical(t$base, "tuple")) return(tuple_info(t, env, d$name))
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

instantiate <- function(m, data, generate = FALSE) {
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
  gen <- list()
  for (d in m$data_decls) {
    info <- decl_info(d, env)
    decls[[d$name]] <- info
    if (generate) {
      args <- eval(cg_cons_args(d$type), env)
      g <- generate_value(info, args, d$name)
      # a container of one element stays a JSON array
      if (!identical(info$kind, "tuple") && length(c(info$adims, info$eshape)) && is.null(dim(g)) && length(g) == 1L) g <- I(g)
      gen[[d$name]] <- g
      v <- read_data_value(stats::setNames(list(g), d$name), info)
    } else v <- read_data_value(data, info)
    assign(paste0("V.", d$name), v, envir = env)
    env$.decls <- decls
    args <- eval(cg_cons_args(d$type), env)
    msg <- check_variable(v, info, args)
    if (!is.null(msg)) stop_data(sprintf("variable '%s' %s", d$name, msg))
  }
  if (generate) return(gen)
  extra <- setdiff(names(data), vapply(m$data_decls, `[[`, "", "name"))
  # transformed data
  # sizes of transformed data variables may depend on earlier transformed
  # data, so their metadata is built after the statements run
  if (length(m$tdata)) eval(blk(m$tdata), env)
  for (d in m$td_decls) decls[[d$name]] <- decl_info(d, env)
  env$.decls <- decls
  if (length(m$tdata_checks)) eval(blk(m$tdata_checks), env)
  # data only subexpressions: computed now; one that fails (say an index
  # out of range in a branch that never runs) is left to fail where used
  for (hx in m$hoist) local({
    v <- tryCatch(eval(hx$code, env), error = function(e) e)
    # an active binding re-signals the error on every use (a promise warns
    # "restarting interrupted promise evaluation" from the second use on)
    if (inherits(v, "error")) {
      code <- hx$code
      makeActiveBinding(hx$name, function() eval(code, env), env)
    } else assign(hx$name, v, envir = env)
  })
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
  if (identical(info$kind, "tuple")) return(tuple_read_data(x, info))
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is_cx_base(info$stan_type %||% "")) return(read_cx_data(x, info))
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
  if (!is.list(x[[1L]])) return(flat)
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
    g <- make_function(fm, body, model$env)
    # data dependent factorizations made during the call belong to this
    # instance (see normal_glm_prep())
    store <- model$env$.glm_store
    if (is.null(store)) store <- model$env$.glm_store <- new.env(parent = emptyenv())
    f <- function(...) {
      prev <- .glm_cache$cur
      .glm_cache$cur <- store
      on.exit(.glm_cache$cur <- prev)
      g(...)
    }
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
#' @param upars Unconstrained parameter values: a numeric vector, or a
#'   matrix with one row per point (as CmdStan's `log_prob` method with
#'   `unconstrained_params`). Constrained values are accepted as a named list
#'   or the name of a JSON or R dump file (CmdStan's `constrained_params`).
#' @param jacobian Include the log absolute Jacobian determinant of the
#'   constraining transforms.
#' @param propto Drop constant terms of sampling statements (`~`) and
#'   `_lupdf`/`_lupmf` calls, as the samplers do.
#' @return The log density (a number, or a vector for a matrix of points).
#'   `bayescomp_grad_log_prob()` returns the gradient (a matrix with one row
#'   per point for a matrix) with the log density in attribute `"log_prob"`.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' bayescomp_log_prob(m, 0.5)
#' bayescomp_grad_log_prob(m, 0.5)
#' @export
bayescomp_log_prob <- function(model, upars, jacobian = TRUE, propto = FALSE) {
  f <- model_lp(model, jacobian, propto, TRUE)
  upars <- as_upars(model, upars)
  if (is.matrix(upars)) return(apply(upars, 1L, function(u) f(check_upars(model, u))))
  f(check_upars(model, upars))
}

#' @rdname bayescomp_log_prob
#' @export
bayescomp_grad_log_prob <- function(model, upars, jacobian = TRUE, propto = FALSE) {
  upars <- as_upars(model, upars)
  if (is.matrix(upars)) {
    r <- lapply(seq_len(nrow(upars)), function(i) lp_grad(model, check_upars(model, upars[i, ]), jacobian, propto))
    g <- matrix(unlist(lapply(r, `[[`, "grad")), nrow(upars), model$dim, byrow = TRUE)
    return(structure(g, log_prob = vapply(r, `[[`, 0, "value")))
  }
  r <- lp_grad(model, check_upars(model, upars), jacobian, propto)
  structure(r$grad, log_prob = r$value)
}

# constrained values (a list or a file) to unconstrained ones
as_upars <- function(model, upars) {
  if (is.character(upars)) upars <- bayescomp_read_data(upars)
  if (is.list(upars)) upars <- bayescomp_unconstrain(model, upars)
  upars
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
#' @examples
#' m <- bayescomp_model("parameters { real<lower=0> sigma; } model { sigma ~ exponential(1); }")
#' u <- bayescomp_unconstrain(m, list(sigma = 2))
#' u
#' bayescomp_constrain(m, u)
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
#' @examples
#' m <- bayescomp_model("parameters { vector[2] b; } model { b ~ normal(0, 1); }
#'   generated quantities { real s = sum(b); }")
#' bayescomp_param_names(m)
#' bayescomp_param_names(m, include_gqs = FALSE)
#' @export
bayescomp_param_names <- function(model, include_tparams = TRUE, include_gqs = TRUE) {
  nm <- c(model$params$pars, if (include_tparams) model$params$tps, if (include_gqs) model$params$gqs)
  unlist(lapply(nm, function(n) {
    d <- model$decls[[n]]
    if (identical(d$kind, "tuple")) return(tuple_flat_names(n, d))
    if (is_cx_base(d$stan_type)) return(cx_flat_names(n, c(d$adims, d$eshape)))
    flat_names(n, c(d$adims, d$eshape))
  }))
}

#' Number of unconstrained parameters
#' @param model A `bayescomp_model`.
#' @return The number of unconstrained parameters, an integer.
#' @examples
#' m <- bayescomp_model("parameters { simplex[3] p; real mu; }
#'   model { mu ~ normal(0, 1); }")
#' bayescomp_dim(m)
#' @export
bayescomp_dim <- function(model) model$dim

# flatten a write_array result to a named numeric vector
flatten_draw <- function(lst) {
  # complex values: real and imaginary parts of each element, as in CmdStan
  unlist(lapply(lst, flatten_value), use.names = FALSE)
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

#' Call user defined Stan functions from R
#'
#' Returns the functions of a program's `functions` block as R functions, as
#' stanc's `--standalone-functions` and CmdStanR's `expose_functions()` do.
#' Arguments are R values in the shapes the package uses for data (numbers,
#' vectors, matrices, and lists or matrices for arrays). Overloaded functions
#' pick the overload by the number and types of the arguments. `_rng`
#' functions use R's random number generator, so [set.seed()] makes them
#' reproducible. Functions ending in `_lp` need a log density accumulator and
#' are not exposed.
#'
#' @param x A `bayescomp_model`, Stan code, or a `.stan` file name.
#' @param ... Passed to [bayescomp_model()] when `x` is code or a file.
#' @return A named list of R functions.
#' @examples
#' f <- bayescomp_expose_functions("functions { real twice(real x) { return 2 * x; } }")
#' f$twice(21)
#' @export
bayescomp_expose_functions <- function(x, ...) {
  m <- if (inherits(x, "bayescomp_model")) x
       else if (is.character(x) && length(x) == 1L && grepl("\\.stan$", x) && file.exists(x)) bayescomp_model(file = x, ...)
       else bayescomp_model(code = x, ...)
  fds <- Filter(function(fd) !is.null(fd$body) && !isTRUE(fd$sig$lp), m$ir$blocks$functions$body %||% list())
  fenv <- parent.env(m$env)
  ctx <- m$env$.ctx
  by_name <- split(fds, vapply(fds, function(fd) fd$name, ""))
  lapply(by_name, function(ov) {
    impls <- lapply(ov, function(fd) get(fd$sig$mangled, envir = fenv))
    types <- lapply(ov, function(fd) lapply(fd$args, `[[`, "type"))
    rets <- lapply(ov, function(fd) fd$sig$ret)
    nm <- ov[[1L]]$name
    function(...) {
      a <- list(...)
      # the viable overload with the fewest promotions, as in Stan
      cost <- vapply(seq_along(ov), function(k) {
        if (length(types[[k]]) != length(a) || !all(mapply(r_value_fits, types[[k]], a))) return(Inf)
        sum(mapply(r_value_cost, types[[k]], a))
      }, 0)
      if (!any(is.finite(cost))) stop(sprintf("no signature of '%s' takes these arguments", nm), call. = FALSE)
      k <- which.min(cost)
      args <- mapply(r_to_stan_value, types[[k]], a, SIMPLIFY = FALSE)
      out <- do.call(impls[[k]], c(unname(args), list(ctx)))
      if (is_cx_ty(rets[[k]])) cx_value(out) else out
    }
  })
}

# Stan type of an R value's elements: int, real or complex
r_elem_base <- function(v) {
  if (is.complex(v)) "complex"
  else if (is.integer(v) || is.logical(v)) "int"
  else "real"
}
# promotions needed to pass an R value as a Stan argument
r_value_cost <- function(type, v) {
  if (is.list(v)) return(0)
  to <- if (is_cx_ty(type)) "complex" else if (type$base == "int") "int" else "real"
  from <- r_elem_base(v)
  if (from == "real" && to == "int") return(0.5)
  promo_cost(ty(to), ty(from))
}
# an R complex argument in the package's storage of complex values
r_to_stan_value <- function(type, v) {
  if (!is_cx_ty(type) || is.list(v)) return(v)
  z <- as.complex(v)
  cx_from_z(z, if (type$nd == 0L && type$base == "complex") integer(0) else z_dims(z))
}

# does an R value fit a Stan argument type (for overload selection)
r_value_fits <- function(type, v) {
  nd <- type$nd %||% 0L
  if (nd > 0L) return(is.list(v) || is.numeric(v) || is.complex(v))
  switch(type$base,
    int = (is.numeric(v) || is.logical(v)) && length(v) == 1L && is.null(dim(v)) && v == round(v),
    real = (is.numeric(v) || is.logical(v)) && length(v) == 1L && is.null(dim(v)),
    complex = (is.numeric(v) || is.complex(v)) && length(v) == 1L,
    vector = , row_vector = , complex_vector = , complex_row_vector = (is.numeric(v) || is.complex(v)) && is.null(dim(v)),
    matrix = , complex_matrix = is.matrix(v),
    TRUE)
}

#' Model information
#'
#' The information stanc's `--info` option prints: the type (`int`, `real`,
#' `complex` or a list for tuples) and number of dimensions of every data
#' variable, parameter, transformed parameter and generated quantity, the
#' library functions and distributions the program uses, and the included
#' files.
#'
#' @param x A `bayescomp_model`, Stan code, or a `.stan` file name.
#' @param include_paths Directories searched for `#include` files.
#' @return A named list with `inputs`, `parameters`,
#'   `transformed parameters`, `generated quantities`, `functions`,
#'   `distributions` and `included_files`.
#' @examples
#' str(bayescomp_info("data { int N; vector[N] y; } parameters { real mu; }
#'                     model { y ~ normal(mu, 1); }"))
#' @export
bayescomp_info <- function(x, include_paths = NULL) {
  if (inherits(x, "bayescomp_model")) { ir <- x$ir; src <- x$code }
  else {
    file <- if (is.character(x) && length(x) == 1L && grepl("\\.stan$", x) && file.exists(x)) x
    src <- read_stan_source(if (is.null(file)) x, file, include_paths)
    ir <- stan_check(stan_parse(src))
  }
  vars <- function(block) {
    ds <- top_decls(ir, block)
    stats::setNames(lapply(ds, function(d) info_type(d$type)), vapply(ds, `[[`, "", "name"))
  }
  funs <- character(); dists <- character()
  walk <- function(n) {
    if (!is.list(n)) return()
    if (identical(n$kind, "tilde") && is.character(n$fname) && !isTRUE(n$sig$user)) dists <<- c(dists, n$fname)
    if (identical(n$kind, "call") && is.list(n$sig) && !isTRUE(n$sig$user)) {
      f <- n$resolved %||% n$name
      if (grepl("_(lpdf|lupdf|lpmf|lupmf|cdf|lcdf|lccdf)$", f)) dists <<- c(dists, f) else funs <<- c(funs, f)
    }
    # stanc lists print, reject and fatal_error, and operators called by name
    if (is.character(n$kind) && n$kind %in% c("print", "reject", "fatal_error")) funs <<- c(funs, n$kind)
    if (is.character(n$fun_name)) funs <<- c(funs, n$fun_name)
    for (e in n) if (is.list(e) && !is.function(e)) walk(e)
  }
  walk(ir$blocks)
  list(inputs = vars("data"), parameters = vars("parameters"),
       `transformed parameters` = vars("transformed parameters"),
       `generated quantities` = vars("generated quantities"),
       functions = sort(unique(funs)), distributions = sort(unique(dists)),
       included_files = as.character(attr(src, "included") %||% character()))
}

info_type <- function(t) {
  if (identical(t$base, "tuple"))
    return(list(type = lapply(t$elems, info_type), dimensions = length(t$adims)))
  inner <- if (grepl("vector|simplex|ordered|positive_ordered|sum_to_zero_vector|unit_vector", t$base)) 1L
           else if (grepl("matrix|cov|corr|cholesky|stochastic", t$base)) 2L else 0L
  base <- if (t$base == "int") "int" else if (startsWith(t$base, "complex")) "complex" else "real"
  list(type = base, dimensions = length(t$adims) + inner)
}

#' Stan version implemented
#'
#' The version of the Stan language and library that bayescompR follows
#' (stanc's `--version`).
#'
#' @return A character string such as `"2.40.0"`.
#' @examples
#' bayescomp_stan_version()
#' @export
bayescomp_stan_version <- function() stan_signatures_version

#' Random data and initial values
#'
#' `bayescomp_generate_data()` draws a random value for every variable of the
#' data block, in declaration order so that sizes can depend on earlier
#' integers, as stanc's `--debug-generate-data` does: integers are drawn
#' within their bounds (1 to 6 when unbounded), other values are drawn
#' uniformly on (-2, 2) on the unconstrained scale and transformed to their
#' declared constraints. `bayescomp_generate_inits()` draws initial values of
#' the parameters the same way (`--debug-generate-inits`).
#'
#' @param x A `bayescomp_model`, Stan code, or a `.stan` file name; for
#'   `bayescomp_generate_inits()`, a `bayescomp_model`.
#' @param seed Random seed.
#' @param file If given, the values are also written to this JSON file.
#' @return A named list, invisibly when `file` is given.
#' @examples
#' code <- "data { int<lower=1> N; vector[N] y; } parameters { real mu; }
#'          model { y ~ normal(mu, 1); }"
#' d <- bayescomp_generate_data(code, seed = 1)
#' m <- bayescomp_model(code, data = d)
#' bayescomp_generate_inits(m, seed = 1)
#' @export
bayescomp_generate_data <- function(x, seed = NULL, file = NULL) {
  ir <- if (inherits(x, "bayescomp_model")) x$ir else {
    f <- if (is.character(x) && length(x) == 1L && grepl("\\.stan$", x) && file.exists(x)) x
    stan_check(stan_parse(read_stan_source(if (is.null(f)) x, f)))
  }
  m <- compile_r(ir)
  out <- with_seed(seed, instantiate(m, list(), generate = TRUE))
  if (is.null(file)) return(out)
  bayescomp_write_data(out, file)
  invisible(out)
}

#' @rdname bayescomp_generate_data
#' @export
bayescomp_generate_inits <- function(x, seed = NULL, file = NULL) {
  check_model(x)
  out <- with_seed(seed, bayescomp_constrain(x, stats::runif(x$dim, -2, 2), include_tparams = FALSE, include_gqs = FALSE))
  if (is.null(file)) return(out)
  bayescomp_write_data(out, file)
  invisible(out)
}

generate_value <- function(info, args, name) {
  if (identical(info$kind, "tuple")) {
    # a tuple as bayescomp_read_data() gives it: elements named "1", "2", ...
    one <- function() {
      r <- lapply(seq_along(info$elems), function(k) {
        e <- info$elems[[k]]
        a <- args[[k]] %||% list()
        if (!identical(e$kind, "tuple")) a <- lapply(a, function(x) if (length(val(x)) > 1L) val(x) else x)
        generate_value(e, a, name)
      })
      stats::setNames(r, as.character(seq_along(r)))
    }
    if (!length(info$adims)) return(one())
    return(lapply(seq_len(prod(info$adims)), function(i) one()))
  }
  n_elem <- prod(c(info$adims, info$eshape, 1L))
  if (identical(info$base, "int")) {
    lo <- if (!is.null(args$lower)) val(args$lower) else 1
    hi <- if (!is.null(args$upper)) val(args$upper) else lo + 5
    if (!is.null(args$upper) && is.null(args$lower)) lo <- hi - 5
    lo <- rep_len(lo, n_elem); hi <- rep_len(hi, n_elem)
    v <- as.integer(lo + floor(stats::runif(n_elem) * (hi - lo + 1)))
    return(if (length(info$adims) > 1L) array(v, info$adims) else if (!length(info$adims)) v[1L] else v)
  }
  r <- constrain_variable(stats::runif(info$free_size, -2, 2), info, args)$x
  r <- val(r)
  if (!startsWith(info$base %||% "", "complex")) return(r)
  # complex values are stored with a trailing [re, im] dimension
  d <- dim(r)
  if (is.null(d)) return(complex(real = r[1L], imaginary = r[2L]))
  k <- length(d)
  h <- prod(d[-k])
  z <- complex(real = r[seq_len(h)], imaginary = r[h + seq_len(h)])
  if (k > 2L) dim(z) <- d[-k]
  z
}
