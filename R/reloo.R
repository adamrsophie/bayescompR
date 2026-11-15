# Exact leave-one-out refits for the observations where PSIS-LOO is
# unreliable (Vehtari, Gelman and Gabry 2017, section 2.1).

#' Refit the model for observations with large Pareto k
#'
#' PSIS-LOO is unreliable for observations whose Pareto k exceeds
#' `k_threshold`. For each of them, `bayescomp_reloo()` fits the model to the
#' data without the observation and computes its exact log predictive
#' density, \eqn{\log \frac{1}{S}\sum_s p(y_i \mid \theta^{(s)}_{-i})}, then
#' updates the LOO estimates.
#'
#' The held out log likelihood is evaluated with the full data model at the
#' draws of the refit, from the generated quantity `variable` when the fit
#' has it and otherwise from [bayescomp_log_lik()], so the parameters must
#' not depend on the number of observations.
#'
#' @param fit A `bayescomp_fit` from [bayescomp_sample()].
#' @param data_minus Function of the observation index `i` (a column of the
#'   pointwise log likelihood) returning the data list without that
#'   observation.
#' @param loo_object The result of `loo::loo(fit)`; computed when `NULL`.
#' @param k_threshold Pareto k above which an observation is refitted.
#' @param variable Name of the pointwise log likelihood generated quantity.
#' @param ... Passed to [bayescomp_sample()] for the refits; by default they
#'   use the chains, warmup and sampling iterations of `fit`.
#' @return The updated `loo` object, with Pareto k set to 0 for the
#'   refitted observations and their indices in element `refits`.
#' @references Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical
#'   Bayesian model evaluation using leave-one-out cross-validation and WAIC.
#'   Statistics and Computing 27, 1413-1432.
#' @examples
#' \donttest{
#' if (requireNamespace("loo", quietly = TRUE)) {
#'   y <- c(-0.3, 0.1, 0.4, -0.2, 6)
#'   m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; }
#'     model { mu ~ normal(0, 10); y ~ normal(mu, 1); }", data = list(N = 5, y = y))
#'   fit <- bayescomp_sample(m, chains = 2, iter_warmup = 300, iter_sampling = 300,
#'                           seed = 1, refresh = 0)
#'   bayescomp_reloo(fit, function(i) list(N = 4, y = y[-i]), k_threshold = 0.5, refresh = 0)
#' }
#' }
#' @export
bayescomp_reloo <- function(fit, data_minus, loo_object = NULL, k_threshold = 0.7, variable = "log_lik", ...) {
  if (!requireNamespace("loo", quietly = TRUE)) stop("the loo package is required", call. = FALSE)
  if (!inherits(fit, "bayescomp_fit") || !identical(fit$method, "sample"))
    stop("'fit' must be a fit from bayescomp_sample()", call. = FALSE)
  if (!is.function(data_minus)) stop("'data_minus' must be a function of the observation index", call. = FALSE)
  lo <- loo_object %||% loo.bayescomp_fit(fit, variable = variable)
  idx <- which(lo$diagnostics$pareto_k > k_threshold)
  lo$refits <- idx
  if (!length(idx)) return(lo)
  args <- utils::modifyList(list(chains = length(fit$unconstrained), iter_warmup = fit$config$iter_warmup,
                                 iter_sampling = fit$config$iter_sampling, seed = fit$seed), list(...))
  gq <- has_variable(fit, variable)
  full_ll <- if (gq) draws_matrix_of(fit, variable) else bayescomp_log_lik(fit)
  for (i in idx) {
    m <- bayescomp_model(fit$model$code, data = data_minus(i), seed = fit$model$compiled$seed)
    if (m$dim != fit$model$dim) stop("the model without an observation has different parameters", call. = FALSE)
    refit <- do.call(bayescomp_sample, c(list(m), args))
    ll <- if (gq) draws_matrix_of(bayescomp_generate_quantities(fit$model, refit), variable)[, i]
          else { refit$model <- fit$model; bayescomp_log_lik(refit)[, i] }
    elpd <- log_mean_exp(ll)
    lpd <- log_mean_exp(full_ll[, i])
    lo$pointwise[i, c("elpd_loo", "mcse_elpd_loo", "p_loo", "looic")] <- c(elpd, 0, lpd - elpd, -2 * elpd)
    if ("influence_pareto_k" %in% colnames(lo$pointwise)) lo$pointwise[i, "influence_pareto_k"] <- 0
    lo$diagnostics$pareto_k[i] <- 0
  }
  n <- nrow(lo$pointwise)
  for (nm in c("elpd_loo", "p_loo", "looic")) {
    v <- lo$pointwise[, nm]
    lo$estimates[nm, ] <- c(sum(v), sqrt(n * stats::var(v)))
    lo[[nm]] <- sum(v)
    lo[[paste0("se_", nm)]] <- sqrt(n * stats::var(v))
  }
  lo
}

log_mean_exp <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
