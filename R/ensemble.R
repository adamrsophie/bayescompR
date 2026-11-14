# Ensemble samplers: the affine invariant stretch move (Goodman and Weare
# 2010) and differential evolution Markov chain Monte Carlo (ter Braak
# 2006). Both move a population of walkers using the positions of the
# others, so they adapt to linear correlations without gradients or tuning.
# Each walker is reported as a chain. Walkers are not independent, so R-hat
# and ESS computed across walkers are optimistic; run several ensembles
# (`chains`) for a fair convergence check.

#' Ensemble Markov chain Monte Carlo
#'
#' @param model A `bayescomp_model`.
#' @param method `"stretch"` (affine invariant ensemble sampler) or `"de"`
#'   (differential evolution MCMC with occasional mode jumping steps).
#' @param walkers Number of walkers (default `2 * dim + 2`).
#' @param iter_warmup,iter_sampling Iterations per walker.
#' @param thin Keep every `thin`th iteration.
#' @param a Stretch scale for `"stretch"`.
#' @param gamma Jump scale for `"de"`, default `2.38 / sqrt(2 * dim)`; every
#'   tenth iteration uses 1 to jump between modes.
#' @param noise Standard deviation of the small Gaussian added to `"de"`
#'   proposals.
#' @param init,seed As in [bayescomp_sample()].
#' @return A `bayescomp_fit` with one chain per walker and element
#'   `acceptance` (mean acceptance per walker).
#' @references Goodman, J. and Weare, J. (2010). Ensemble samplers with
#'   affine invariance. *Communications in Applied Mathematics and
#'   Computational Science* 5, 65 to 80. ter Braak, C. J. F. (2006). A Markov
#'   chain Monte Carlo version of the genetic algorithm differential
#'   evolution. *Statistics and Computing* 16, 239 to 249.
#' @export
bayescomp_ensemble <- function(model, method = c("stretch", "de"), walkers = NULL, iter_warmup = 1000L,
                               iter_sampling = 1000L, thin = 1L, a = 2, gamma = NULL, noise = 1e-4,
                               init = 2, seed = NULL) {
  check_model(model)
  method <- match.arg(method)
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  W <- if (is.null(walkers)) 2L * n + 2L else as_count(walkers, "walkers", 3L)
  if (method == "stretch" && W %% 2L == 1L) W <- W + 1L
  nw <- as_count(iter_warmup, "iter_warmup")
  ns <- as_count(iter_sampling, "iter_sampling", 1L)
  thin <- as_count(thin, "thin", 1L)
  if (is.null(gamma)) gamma <- 2.38 / sqrt(2 * n)
  t0 <- proc.time()[["elapsed"]]
  out <- run_chains(model, 1L, seed, 1L, function(id) {
    lpf <- lp_value_fun(model)
    X <- t(vapply(seq_len(W), function(w) find_init(model, init, w), numeric(n)))
    if (n == 1L) X <- matrix(X, ncol = 1L)
    lp <- vapply(seq_len(W), function(w) lpf(X[w, ]), 0)
    keep <- seq_len(ns)[(seq_len(ns) - 1L) %% thin == 0L]
    Q <- array(NA_real_, c(length(keep), W, n))
    LP <- matrix(NA_real_, length(keep), W)
    acc <- numeric(W)
    row <- 0L
    for (it in seq_len(nw + ns)) {
      if (method == "stretch") {
        half <- list(seq(1L, W, by = 2L), seq(2L, W, by = 2L))
        for (h in 1:2) {
          S <- half[[h]]; O <- half[[3L - h]]
          for (j in S) {
            k <- O[sample.int(length(O), 1L)]
            z <- ((a - 1) * stats::runif(1) + 1)^2 / a
            y <- X[k, ] + z * (X[j, ] - X[k, ])
            lpy <- lpf(y)
            la <- (n - 1) * log(z) + lpy - lp[j]
            ok <- is.finite(la) && log(stats::runif(1)) < la
            if (ok) { X[j, ] <- y; lp[j] <- lpy }
            if (it > nw) acc[j] <- acc[j] + ok
          }
        }
      } else {
        g <- if (it %% 10L == 0L) 1 else gamma
        for (j in seq_len(W)) {
          r <- sample(setdiff(seq_len(W), j), 2L)
          y <- X[j, ] + g * (X[r[1L], ] - X[r[2L], ]) + stats::rnorm(n, 0, noise)
          lpy <- lpf(y)
          la <- lpy - lp[j]
          ok <- is.finite(la) && log(stats::runif(1)) < la
          if (ok) { X[j, ] <- y; lp[j] <- lpy }
          if (it > nw) acc[j] <- acc[j] + ok
        }
      }
      if (it > nw && ((it - nw - 1L) %% thin == 0L)) {
        row <- row + 1L
        Q[row, , ] <- X
        LP[row, ] <- lp
      }
    }
    list(Q = Q, LP = LP, acc = acc / ns)
  })[[1L]]
  res <- lapply(seq_len(W), function(w) {
    q <- matrix(out$Q[, w, ], ncol = n)
    dg <- matrix(0, nrow(q), 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- out$LP[, w]
    dg[, "energy__"] <- NA_real_
    list(q = q, gq = chain_write_draws(model, q, FALSE), diag = dg, step_size = NA_real_, inv_metric = NULL,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  fit <- make_fit(model, res, list(algorithm = method, walkers = W), "ensemble", seed, proc.time()[["elapsed"]] - t0)
  fit$acceptance <- out$acc
  fit
}
