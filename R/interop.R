# Methods for the generics of the Stan ecosystem (posterior, bayesplot, loo,
# rstantools),
# registered when those packages are loaded, so tools written for rstan or
# cmdstanr fits accept a bayescomp_fit.

#' Methods for posterior, bayesplot, loo and rstantools
#'
#' A `bayescomp_fit` works with the generics that rstan and cmdstanr fits
#' support:
#' * `posterior::as_draws()`, `as_draws_array()`, `as_draws_df()` and
#'   `as_draws_matrix()` return the draws (warmup excluded) of all variables,
#'   or of `variable`;
#' * `bayesplot::nuts_params()` returns the sampler diagnostics in the long
#'   format of `bayesplot::mcmc_nuts_divergence()` and related plots, and
#'   `bayesplot::log_posterior()` the `lp__` draws;
#' * `loo::loo()` computes PSIS-LOO from the pointwise log likelihood stored
#'   in the generated quantity `variable` (default `"log_lik"`), with
#'   relative efficiencies from the chains; without a `log_lik` generated
#'   quantity it uses [bayescomp_log_lik()];
#' * `rstantools::log_lik()` and `rstantools::posterior_predict()` return
#'   draws by elements matrices (chains stacked) of the generated
#'   quantities `variable` (defaults `"log_lik"` and `"y_rep"`), or of
#'   [bayescomp_log_lik()] and [bayescomp_posterior_predict()] when the fit
#'   has no such variable.
#' @param x,object A `bayescomp_fit`.
#' @param variable Variable names; for `loo()` the name of the pointwise log
#'   likelihood.
#' @param r_eff For `loo()`, `TRUE` to compute relative efficiencies from
#'   the chains, or a vector of them.
#' @param ... Passed on (to [loo::loo.array()] for `loo()`).
#' @return Objects of the classes the generics document.
#' @examples
#' m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; }
#'   model { y ~ normal(mu, 1); }
#'   generated quantities {
#'     vector[N] log_lik;
#'     for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu, 1);
#'   }",
#'   data = list(N = 10, y = seq(-1, 1, length.out = 10)))
#' fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100,
#'                         seed = 1, refresh = 0)
#' if (requireNamespace("posterior", quietly = TRUE))
#'   print(posterior::summarise_draws(posterior::as_draws_df(fit, "mu")))
#' if (requireNamespace("loo", quietly = TRUE)) print(loo::loo(fit))
#' @name bayescomp_interop
NULL

fit_draws <- function(x, variable = NULL) bayescomp_draws(x, variable, format = "array")

#' @rdname bayescomp_interop
#' @exportS3Method posterior::as_draws
as_draws.bayescomp_fit <- function(x, variable = NULL, ...) posterior::as_draws_array(fit_draws(x, variable))

#' @rdname bayescomp_interop
#' @exportS3Method posterior::as_draws_array
as_draws_array.bayescomp_fit <- function(x, variable = NULL, ...) posterior::as_draws_array(fit_draws(x, variable))

#' @rdname bayescomp_interop
#' @exportS3Method posterior::as_draws_df
as_draws_df.bayescomp_fit <- function(x, variable = NULL, ...) posterior::as_draws_df(posterior::as_draws_array(fit_draws(x, variable)))

#' @rdname bayescomp_interop
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.bayescomp_fit <- function(x, variable = NULL, ...) posterior::as_draws_matrix(posterior::as_draws_array(fit_draws(x, variable)))

post_warmup <- function(a, x) {
  w <- x$n_warmup_saved %||% 0L
  if (w > 0L) a[-seq_len(w), , , drop = FALSE] else a
}

long_format <- function(a) {
  it <- dim(a)[1L]; ch <- dim(a)[2L]; nm <- dimnames(a)[[3L]]
  data.frame(Chain = rep(rep(seq_len(ch), each = it), length(nm)),
             Iteration = rep(seq_len(it), ch * length(nm)),
             Parameter = factor(rep(nm, each = it * ch), levels = nm),
             Value = as.vector(a))
}

#' @rdname bayescomp_interop
#' @exportS3Method bayesplot::nuts_params
nuts_params.bayescomp_fit <- function(object, ...) {
  a <- post_warmup(object$sampler_diagnostics, object)
  keep <- intersect(c("accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__"), dimnames(a)[[3L]])
  long_format(a[, , keep, drop = FALSE])
}

#' @rdname bayescomp_interop
#' @exportS3Method bayesplot::log_posterior
log_posterior.bayescomp_fit <- function(object, ...) {
  d <- long_format(fit_draws(object, "lp__"))
  d[c("Chain", "Iteration", "Value")]
}

#' @rdname bayescomp_interop
#' @exportS3Method loo::loo
loo.bayescomp_fit <- function(x, variable = "log_lik", r_eff = TRUE, ...) {
  # no such generated quantity: take the terms from the model block
  a <- if (identical(variable, "log_lik") && !has_variable(x, variable)) bayescomp_log_lik(x, "array")
       else unclass(fit_draws(x, variable))
  if (!dim(a)[3L]) stop(sprintf("the fit has no variable '%s'", variable), call. = FALSE)
  if (isTRUE(r_eff)) {
    # relative efficiency does not change when each observation's likelihood
    # is rescaled, so shift by its maximum to avoid underflow
    mx <- apply(a, 3L, function(z) { z <- z[is.finite(z)]; if (length(z)) max(z) else 0 })
    r_eff <- loo::relative_eff(exp(sweep(a, 3L, mx)))
  }
  else if (isFALSE(r_eff)) r_eff <- NULL
  loo::loo(a, r_eff = r_eff, ...)
}

# draws x elements matrix of one variable, chains stacked
draws_matrix_of <- function(fit, variable) {
  a <- unclass(fit_draws(fit, variable))
  if (!dim(a)[3L]) stop(sprintf("the fit has no variable '%s'", variable), call. = FALSE)
  m <- matrix(a, ncol = dim(a)[3L])
  colnames(m) <- dimnames(a)[[3L]]
  m
}

#' @rdname bayescomp_interop
#' @exportS3Method rstantools::log_lik
log_lik.bayescomp_fit <- function(object, variable = "log_lik", ...) {
  if (identical(variable, "log_lik") && !has_variable(object, variable)) return(bayescomp_log_lik(object))
  draws_matrix_of(object, variable)
}

#' @rdname bayescomp_interop
#' @exportS3Method rstantools::posterior_predict
posterior_predict.bayescomp_fit <- function(object, variable = "y_rep", ...) {
  if (identical(variable, "y_rep") && !has_variable(object, variable)) return(bayescomp_posterior_predict(object))
  draws_matrix_of(object, variable)
}
has_variable <- function(fit, variable) {
  nm <- dimnames(fit$draws)[[3L]]
  any(nm == variable | sub("\\[.*$", "", nm) == variable)
}
