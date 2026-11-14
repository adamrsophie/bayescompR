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
#' @return A list with `ir_version`, `stan_version`, `blocks` and
#'   `functions`, or its JSON serialization.
#' @examples
#' ir <- bayescomp_ir("parameters { real y; } model { y ~ std_normal(); }")
#' ir$ir_version
#' cat(substr(bayescomp_ir("parameters { real y; } model { y ~ std_normal(); }", json = TRUE), 1, 80))
#' @export
bayescomp_ir <- function(x, json = FALSE) {
  ir <- if (inherits(x, "bayescomp_model")) x$ir else stan_check(stan_parse(read_stan_source(x, NULL)))
  out <- list(ir_version = ir_version, stan_version = "2.40.0", blocks = strip_ir(ir$blocks))
  if (json) to_json(out) else out
}

strip_ir <- function(x) {
  if (is.function(x)) return(NULL)
  if (is.list(x)) {
    if (!is.null(x$impl) || !is.null(x$args) && !is.null(x$ret) && !is.null(x$name) && !is.null(x$loc) && isTRUE(x$user)) {
      # a resolved signature
      return(list(name = x$name, user = isTRUE(x$user), mangled = x$mangled,
                  args = if (isTRUE(x$user)) lapply(x$args, strip_ir) else unlist(x$args),
                  ret = if (is.function(x$ret)) "computed" else strip_ir(x$ret)))
    }
    nm <- names(x)
    out <- lapply(x, strip_ir)
    keep <- !vapply(out, is.null, TRUE) | vapply(x, is.null, TRUE)
    out <- out[keep]
    if (!is.null(nm)) names(out) <- nm[keep]
    return(out)
  }
  x
}

json_str <- function(s) {
  s <- gsub("\\\\", "\\\\\\\\", s)
  s <- gsub("\"", "\\\\\"", s)
  s <- gsub("\n", "\\\\n", s)
  s <- gsub("\t", "\\\\t", s)
  paste0("\"", s, "\"")
}

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
  else if (is.numeric(x)) v <- ifelse(is.finite(x), format(x, digits = 17, scientific = NA, trim = TRUE),
                                      ifelse(is.nan(x), "\"NaN\"", ifelse(x > 0, "\"Infinity\"", "\"-Infinity\"")))
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
#' @param fit A `bayescomp_fit`.
#' @param seed Random seed.
#' @return A `bayescomp_fit` whose draws hold the new quantities.
#' @export
bayescomp_generate_quantities <- function(model, fit, seed = NULL) {
  check_model(model)
  if (!inherits(fit, "bayescomp_fit")) stop("'fit' must be a bayescomp_fit", call. = FALSE)
  if (model$dim != ncol(fit$unconstrained[[1L]]))
    stop("the model and the fit have different numbers of unconstrained parameters", call. = FALSE)
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
