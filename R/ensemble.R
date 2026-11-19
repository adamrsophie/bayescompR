# Ensemble samplers: the affine invariant stretch move (Goodman and Weare
# 2010), the side move (Chen et al. 2025), ensemble slice sampling
# (Karamanis and Beutler 2021) and differential evolution Markov chain Monte
# Carlo (ter Braak 2006). All move a population of walkers using the
# positions of the others, so they adapt to linear correlations without
# gradients.
# Each walker is reported as a chain. Walkers of one ensemble are not
# independent, so summaries compute R-hat from the walker means of the
# independent ensembles (`chains`) and not across walkers.

#' Ensemble Markov chain Monte Carlo
#'
#' @param model A `bayescomp_model`.
#' @param method `"stretch"` (affine invariant ensemble sampler), `"side"`
#'   (side move: a Gaussian step along the difference of two walkers of the
#'   other half, which scales better with dimension than the stretch move),
#'   `"slice"` (ensemble slice sampling along such differences; no rejections
#'   and a scale tuned during warmup) or `"de"` (differential evolution MCMC
#'   with occasional mode jumping steps).
#' @param walkers Number of walkers (default `2 * dim + 2`).
#' @param iter_warmup,iter_sampling Iterations per walker.
#' @param thin Keep every `thin`th iteration.
#' @param a Stretch scale for `"stretch"`. The step scale of `"side"` starts at
#'   `1.6869 / sqrt(dim)` and the slice scale of `"slice"` at 1; both are tuned
#'   during warmup.
#' @param gamma Jump scale for `"de"`, default `2.38 / sqrt(2 * dim)`; every
#'   tenth iteration uses 1 to jump between modes.
#' @param noise Standard deviation of the small Gaussian added to `"de"`
#'   proposals.
#' @param chains Number of independent ensembles.
#' @param init,seed,cores As in [bayescomp_sample()].
#' @return A `bayescomp_fit` with one chain per walker (the walkers of
#'   ensemble 1 first) and element `acceptance` (mean acceptance per walker;
#'   1 for `"slice"`, which reports `evals`, the mean log density evaluations
#'   per walker update, instead).
#'   Its summary reports R-hat computed from the per iteration walker means of
#'   the ensembles (split R-hat when `chains = 1`) and effective sample sizes
#'   computed over all walkers, as in emcee; run several ensembles for a
#'   reliable convergence check.
#' @references Goodman, J. and Weare, J. (2010). Ensemble samplers with
#'   affine invariance. *Communications in Applied Mathematics and
#'   Computational Science* 5, 65 to 80. ter Braak, C. J. F. (2006). A Markov
#'   chain Monte Carlo version of the genetic algorithm differential
#'   evolution. *Statistics and Computing* 16, 239 to 249. Karamanis, M. and
#'   Beutler, F. (2021). Ensemble slice sampling. *Statistics and Computing*
#'   31, 61. Chen, Y. et al. (2025). New affine invariant ensemble samplers
#'   and their dimensional scaling. arXiv:2505.02987.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_ensemble(m, iter_warmup = 100, iter_sampling = 100, seed = 1)
#' mean(bayescomp_extract(fit, "mu"))
#' @export
bayescomp_ensemble <- function(model, method = c("stretch", "side", "slice", "de"), walkers = NULL, iter_warmup = 1000L,
                               iter_sampling = 1000L, thin = 1L, a = 2, gamma = NULL, noise = 1e-4,
                               chains = 1L, init = 2, seed = NULL, cores = getOption("mc.cores", 1L)) {
  check_model(model)
  method <- match.arg(method)
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  W <- if (is.null(walkers)) 2L * n + 2L else as_count(walkers, "walkers", 3L)
  if (method != "de" && W %% 2L == 1L) W <- W + 1L
  # the stretch move keeps wrong means and variances for long in high
  # dimensions (Huijser, Goodman and Brewer 2022)
  if (method == "stretch" && n > 50L)
    warning(sprintf("the stretch move mixes poorly in %d dimensions; consider method = \"slice\" or bayescomp_sample()", n), call. = FALSE)
  # the walkers span the affine hull of their starting points: with W - 1 < n
  # the ensemble never leaves a lower dimensional subspace
  if (W <= n + 1L)
    stop(sprintf("'walkers' must exceed the number of parameters plus one (%d); use at least %d", n + 1L, 2L * n), call. = FALSE)
  if (W < 2L * n) warning(sprintf("%d walkers for %d parameters; at least %d are recommended", W, n, 2L * n), call. = FALSE)
  nw <- as_count(iter_warmup, "iter_warmup")
  ns <- as_count(iter_sampling, "iter_sampling", 1L)
  thin <- as_count(thin, "thin", 1L)
  chains <- as_count(chains, "chains", 1L)
  if (is.null(gamma)) gamma <- 2.38 / sqrt(2 * n)
  t0 <- proc.time()[["elapsed"]]
  rc <- run_chains(model, chains, seed, cores, function(id) {
    lpf <- lp_value_fun(model, propto = TRUE)
    X <- t(vapply(seq_len(W), function(w) find_init(model, init, (id - 1L) * W + w, gradient = FALSE), numeric(n)))
    if (n == 1L) X <- matrix(X, ncol = 1L)
    lp <- vapply(seq_len(W), function(w) lpf(X[w, ]), 0)
    keep <- seq_len(ns)[(seq_len(ns) - 1L) %% thin == 0L]
    Q <- array(NA_real_, c(length(keep), W, n))
    LP <- matrix(NA_real_, length(keep), W)
    acc <- numeric(W)
    row <- 0L
    half <- list(seq(1L, W, by = 2L), seq(2L, W, by = 2L))
    log_sigma <- log(1.6869 / sqrt(n))
    mu <- 1; tuning <- TRUE; settled <- 0L; evals <- 0
    pick2 <- function(O) O[sample.int(length(O), 2L)]
    for (it in seq_len(nw + ns)) {
      if (method == "side") {
        na <- 0L
        for (h in 1:2) {
          S <- half[[h]]; O <- half[[3L - h]]
          for (j in S) {
            r <- pick2(O)
            y <- X[j, ] + exp(log_sigma) * (X[r[1L], ] - X[r[2L], ]) * stats::rnorm(1)
            lpy <- lpf(y)
            ok <- is.finite(lpy) && log(stats::runif(1)) < lpy - lp[j]
            if (ok) { X[j, ] <- y; lp[j] <- lpy }
            na <- na + ok
            if (it > nw) acc[j] <- acc[j] + ok
          }
        }
        # Robbins-Monro on the log scale toward acceptance 0.5 during warmup
        if (it <= nw) log_sigma <- log_sigma + (na / W - 0.5) / it^0.6
      } else if (method == "slice") {
        nexp <- 0L; ncon <- 0L
        for (h in 1:2) {
          S <- half[[h]]; O <- half[[3L - h]]
          for (j in S) {
            r <- pick2(O)
            d <- 2 * mu * (X[r[1L], ] - X[r[2L], ])
            x0 <- X[j, ]
            f <- function(t) lpf(x0 + t * d)
            y0 <- lp[j] - stats::rexp(1)
            L <- -stats::runif(1); R <- L + 1
            J <- floor(1e4 * stats::runif(1)); K <- 1e4 - 1 - J
            while (J > 0 && f(L) > y0) { L <- L - 1; J <- J - 1; nexp <- nexp + 1L; evals <- evals + 1 }
            while (K > 0 && f(R) > y0) { R <- R + 1; K <- K - 1; nexp <- nexp + 1L; evals <- evals + 1 }
            evals <- evals + 2
            repeat {
              t <- L + stats::runif(1) * (R - L)
              lt <- f(t); evals <- evals + 1
              if (is.finite(lt) && lt > y0) break
              ncon <- ncon + 1L
              if (t < 0) L <- t else R <- t
              if (R - L < 1e-12) { t <- 0; lt <- lp[j]; break }
            }
            X[j, ] <- x0 + t * d; lp[j] <- lt
            if (it > nw) acc[j] <- acc[j] + 1
          }
        }
        # the scale doubles the share of expansions until it is about a half
        # (Karamanis and Beutler 2021, algorithm 4)
        if (tuning && it <= nw) {
          r <- max(1L, nexp) / (max(1L, nexp) + ncon)
          mu <- mu * 2 * r
          if (abs(r - 0.5) < 0.05) settled <- settled + 1L
          if (settled > 5L) tuning <- FALSE
        }
      } else if (method == "stretch") {
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
    list(Q = Q, LP = LP, acc = acc / ns, evals = evals / ((nw + ns) * W))
  })
  # generated quantities draw from a stream derived from the run's seed
  res <- with_seed(attr(rc, "seed") %% 2147483646L + 1L, lapply(seq_len(chains * W), function(cw) {
    out <- rc[[(cw - 1L) %/% W + 1L]]
    w <- (cw - 1L) %% W + 1L
    q <- matrix(out$Q[, w, ], ncol = n)
    dg <- matrix(0, nrow(q), 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- out$LP[, w]
    dg[, "energy__"] <- NA_real_
    list(q = q, gq = chain_write_draws(model, q, FALSE), diag = dg, step_size = NA_real_, inv_metric = NULL,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  }))
  fit <- make_fit(model, res, list(algorithm = method, walkers = W, ensembles = chains), "ensemble", seed,
                  proc.time()[["elapsed"]] - t0)
  fit$acceptance <- unlist(lapply(rc, `[[`, "acc"))
  if (method == "slice") fit$evals <- mean(vapply(rc, `[[`, 0, "evals"))
  fit
}
