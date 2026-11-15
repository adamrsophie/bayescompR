# The typed intermediate representation shared by all backends.
#
# bayescomp_ir() returns the semantic analyzer's output as plain data: every
# node is a named list with a `kind`, expression nodes carry their Stan
# `type` (base and array dimensions) and `data_only` flag, and calls carry
# the resolved signature (name, argument patterns, user or library). R
# function objects are removed, so the IR can be serialized to JSON and read
# by backends written in other languages. The schema is documented in the
# vignette "Backends and the intermediate representation".

#' Typed intermediate representation of a Stan program
#'
#' @param x A `bayescomp_model`, or Stan code as a string.
#' @param json Return a JSON string instead of a list.
#' @return A list with `ir_version`, `stan_version` and `blocks` (one
#'   entry per program block present, including `functions`), or its JSON
#'   serialization.
#' @examples
#' ir <- bayescomp_ir("parameters { real y; } model { y ~ std_normal(); }")
#' ir$ir_version
#' cat(substr(bayescomp_ir("parameters { real y; } model { y ~ std_normal(); }", json = TRUE), 1, 80))
#' @export
bayescomp_ir <- function(x, json = FALSE) {
  ir <- if (inherits(x, "bayescomp_model")) x$ir else stan_check(stan_parse(read_stan_source(x, NULL)))
  blocks <- strip_ir(ir$blocks)
  # blocks is always a JSON object, also for an empty program
  if (!length(blocks)) blocks <- structure(list(), names = character(0))
  out <- list(ir_version = ir_version, stan_version = "2.40.0", blocks = blocks)
  if (json) to_json(out) else out
}

strip_ir <- function(x) {
  if (is.function(x)) return(NULL)
  if (is.list(x)) {
    if (!is.null(x$impl) || !is.null(x$args) && !is.null(x$ret) && !is.null(x$name) && !is.null(x$loc) && isTRUE(x$user)) {
      # a resolved signature
      return(list(name = x$name, user = isTRUE(x$user), mangled = x$mangled,
                  args = if (isTRUE(x$user)) lapply(x$args, strip_ir) else as.list(unlist(x$args)),
                  ret = if (is.function(x$ret)) "computed" else strip_ir(x$ret)))
    }
    nm <- names(x)
    out <- lapply(x, strip_ir)
    # vector valued fields stay arrays in JSON even with one element
    if (is.character(out$kinds)) out$kinds <- as.list(out$kinds)
    keep <- !vapply(out, is.null, TRUE) | vapply(x, is.null, TRUE)
    out <- out[keep]
    if (!is.null(nm)) names(out) <- nm[keep]
    return(out)
  }
  x
}

drop_loc <- function(x) {
  if (!is.list(x)) return(x)
  x$loc <- NULL
  lapply(x, drop_loc)
}

json_str <- function(s) {
  s <- gsub("\\\\", "\\\\\\\\", s)
  s <- gsub("\"", "\\\\\"", s)
  s <- gsub("\n", "\\\\n", s)
  s <- gsub("\t", "\\\\t", s)
  s <- gsub("\r", "\\\\r", s)
  # any other control character as a \u escape, so the JSON stays valid
  ctl <- grepl("[\\x01-\\x1f]", s, perl = TRUE)
  if (any(ctl)) s[ctl] <- vapply(s[ctl], function(x) {
    cs <- strsplit(x, "")[[1L]]
    i <- grepl("[\\x01-\\x1f]", cs, perl = TRUE)
    cs[i] <- sprintf("\\u%04x", vapply(cs[i], utf8ToInt, 0L))
    paste(cs, collapse = "")
  }, "", USE.NAMES = FALSE)
  paste0("\"", s, "\"")
}

# numbers as JSON, vectorized; %.17g round trips doubles exactly
num_json <- function(x) {
  v <- if (is.integer(x)) as.character(x) else sprintf("%.17g", x)
  bad <- !is.finite(x)
  if (any(bad)) v[bad] <- ifelse(is.nan(x[bad]), "\"NaN\"", ifelse(x[bad] > 0, "\"Infinity\"", "\"-Infinity\""))
  v
}

json_array <- function(x) paste0("[", paste(num_json(x), collapse = ","), "]")

to_json <- function(x) {
  if (is.null(x)) return("null")
  if (is.list(x)) {
    if (!is.null(names(x)) && all(nzchar(names(x)))) {
      parts <- vapply(names(x), function(n) paste0(json_str(n), ":", to_json(x[[n]])), "")
      return(paste0("{", paste(parts, collapse = ","), "}"))
    }
    return(paste0("[", paste(vapply(x, to_json, ""), collapse = ","), "]"))
  }
  if (is.character(x)) v <- json_str(x)
  else if (is.logical(x)) v <- ifelse(is.na(x), "null", ifelse(x, "true", "false"))
  else if (is.numeric(x)) v <- num_json(x)
  else v <- json_str(as.character(x))
  if (length(x) == 1L) v else paste0("[", paste(v, collapse = ","), "]")
}

#' Generated quantities for existing draws
#'
#' Reruns the transformed parameters and generated quantities blocks for
#' every draw of a fit, for example after changing the generated quantities
#' or the data of a model with the same parameters.
#'
#' @param model A `bayescomp_model` with the same parameters as the fit.
#' @param fit A `bayescomp_fit`, the names of Stan CSV files (as CmdStan's
#'   `fitted_params`), or draws of the constrained parameters: a posterior
#'   draws object, a 3-D array (iterations, chains, variables) or a matrix or
#'   data frame with one row per draw. Columns are named as in
#'   [bayescomp_param_names()] or in CmdStan's style (`theta.1`); other
#'   columns are ignored.
#' @param seed Random seed.
#' @return A `bayescomp_fit` whose draws hold the new quantities.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }
#'   generated quantities { real y_rep = normal_rng(mu, 1); }")
#' fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100,
#'                         seed = 1, refresh = 0)
#' gq <- bayescomp_generate_quantities(m, fit, seed = 2)
#' mean(bayescomp_extract(gq, "y_rep"))
#' @export
bayescomp_generate_quantities <- function(model, fit, seed = NULL) {
  check_model(model)
  if (is.character(fit)) fit <- bayescomp_read_csv(fit, model = model)
  if (!inherits(fit, "bayescomp_fit")) fit <- fit_from_draws(model, fit)
  if (model$dim != ncol(fit$unconstrained[[1L]]))
    stop("the model and the fit have different numbers of unconstrained parameters", call. = FALSE)
  # the unconstrained draws only mean the same thing under the same ordered
  # parameter declarations (names, shapes, transforms and their bounds)
  if (!is.null(fit$model)) {
    sig <- function(m) lapply(m$params$pars, function(n) {
      d <- m$decls[[n]]
      list(n, d$kind, d$adims, d$eshape,
           deparse(drop_loc(strip_ir(Filter(function(st) identical(st$name, n), top_decls(m$ir, "parameters"))[[1L]]$type))))
    })
    if (!identical(sig(model), sig(fit$model)))
      stop("the model's parameter declarations differ from those of the fit (names, order, shapes or constraints)", call. = FALSE)
    # bounds may depend on data, so the saved parameter values (not the
    # unconstrained coordinates) are carried over to the new model
    if (!identical(model, fit$model)) {
      fit$unconstrained <- lapply(fit$unconstrained, function(q) {
        if (!nrow(q)) return(q)
        out <- vapply(seq_len(nrow(q)), function(i) {
          v <- write_array(fit$model, q[i, ], FALSE, FALSE, print = FALSE)[[1L]]
          tryCatch(bayescomp_unconstrain(model, v), error = function(e)
            stop(sprintf("draw %d of the fit is outside the constraints of the model: %s", i, conditionMessage(e)), call. = FALSE))
        }, numeric(model$dim))
        matrix(out, nrow(q), model$dim, byrow = TRUE)
      })
    }
  }
  chains <- length(fit$unconstrained)
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, chains, seed, 1L, function(id) {
    qs <- fit$unconstrained[[id]]
    dg <- matrix(0, nrow(qs), 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    if (!is.null(fit$draws)) dg[, "lp__"] <- fit$draws[seq_len(nrow(qs)), id, "lp__"]
    list(q = qs, gq = chain_write_draws(model, qs, TRUE), diag = dg, step_size = NA_real_, inv_metric = NULL,
         n_warmup_saved = fit$n_warmup_saved %||% 0L, time = c(warmup = 0, sampling = 0))
  })
  out <- make_fit(model, res, fit$config, "generate_quantities", seed, proc.time()[["elapsed"]] - t0)
  out
}

# a fit from draws of the constrained parameters
fit_from_draws <- function(model, x) {
  if (inherits(x, "draws")) {
    if (!requireNamespace("posterior", quietly = TRUE)) stop("the posterior package is required for draws objects", call. = FALSE)
    x <- unclass(posterior::as_draws_array(x))
  }
  if (is.data.frame(x)) x <- as.matrix(x[, vapply(x, is.numeric, TRUE), drop = FALSE])
  if (is.matrix(x)) x <- array(x, c(nrow(x), 1L, ncol(x)), dimnames = list(NULL, NULL, colnames(x)))
  if (!is.array(x) || length(dim(x)) != 3L || is.null(dimnames(x)[[3L]]))
    stop("'fit' must be a bayescomp_fit, Stan CSV files or named draws of the parameters", call. = FALSE)
  vn <- dimnames(x)[[3L]]
  vn <- ifelse(grepl("[.]", vn) & !grepl("\\[", vn), from_cmdstan_names(vn), vn)
  pn <- bayescomp_param_names(model, FALSE, FALSE)
  miss <- setdiff(pn, vn)
  if (length(miss)) stop(sprintf("the draws have no values for %s", paste(utils::head(miss, 5L), collapse = ", ")), call. = FALSE)
  idx <- match(pn, vn)
  n <- dim(x)[1L]
  un <- lapply(seq_len(dim(x)[2L]), function(c) {
    P <- matrix(x[, c, idx, drop = FALSE], n)
    matrix(vapply(seq_len(n), function(i) bayescomp_unconstrain(model, unflatten_params(model, P[i, ])), numeric(model$dim)),
           n, model$dim, byrow = TRUE)
  })
  structure(list(draws = NULL, unconstrained = un, model = model, n_warmup_saved = 0L, config = list(), method = "draws"),
            class = "bayescomp_fit")
}
