# Adaptive tempered sequential Monte Carlo (Del Moral, Doucet and Jasra 2006).
#
# Stan programs define one joint density, not a separate prior and
# likelihood, so the bridge runs from a normalized Gaussian reference g on
# the unconstrained space to the posterior p:
#
#   pi_beta(u) proportional to g(u)^(1 - beta) p(u)^beta,  0 = beta_0 < ... < beta_K = 1.
#
# Each step picks the next beta by bisection so that the effective sample
# size of the incremental weights is `ess_threshold` times the number of
# particles, resamples systematically, and moves every particle with
# `num_mcmc_steps` random walk Metropolis steps targeting pi_beta, whose
# covariance is the scaled particle covariance. Because g is normalized and
# p is evaluated with all constants, the product of the mean incremental
# weights estimates the normalizing constant of the model's density on the
# unconstrained scale, which is the marginal likelihood when the program
# specifies a proper generative model with normalized densities.

#' Tempered sequential Monte Carlo
#'
#' @param model A `bayescomp_model`.
#' @param num_particles Number of particles.
#' @param reference Gaussian reference distribution: `"pathfinder"` (the
#'   best Pathfinder approximation with its covariance inflated by
#'   `reference_scale^2`) or `"standard"` (independent normals with sd
#'   `reference_scale`, centered at zero).
#' @param reference_scale Scale factor for the reference.
#' @param ess_threshold Target relative effective sample size used to pick
#'   each temperature.
#' @param num_mcmc_steps Metropolis steps per particle and temperature.
#' @param max_stages Maximum number of temperatures.
#' @param seed Random seed.
#' @param cores Workers for the particle moves (forked, not on Windows).
#' @return A `bayescomp_fit` with method `"smc"` and elements
#'   `log_evidence` (log normalizing constant estimate), `betas`
#'   (temperatures) and `acceptance` (mean Metropolis acceptance per stage).
#' @references Del Moral, P., Doucet, A. and Jasra, A. (2006). Sequential
#'   Monte Carlo samplers. *Journal of the Royal Statistical Society B* 68,
#'   411 to 436.
#' @examples
#' \donttest{
#' m <- bayescomp_model("data { real y; } parameters { real mu; }
#'   model { mu ~ normal(0, 1); y ~ normal(mu, 1); }", data = list(y = 1))
#' fit <- bayescomp_smc(m, num_particles = 200, seed = 1)
#' fit$log_evidence  # exact: dnorm(1, 0, sqrt(2), log = TRUE)
#' }
#' @export
bayescomp_smc <- function(model, num_particles = 1000L, reference = c("pathfinder", "standard"),
                          reference_scale = 1.5, ess_threshold = 0.5, num_mcmc_steps = 10L,
                          max_stages = 200L, seed = NULL, cores = 1L) {
  check_model(model)
  reference <- match.arg(reference)
  N <- as_count(num_particles, "num_particles", 10L)
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  if (!(ess_threshold > 0 && ess_threshold < 1)) stop("'ess_threshold' must be in (0, 1)", call. = FALSE)
  t0 <- proc.time()[["elapsed"]]
  info <- NULL
  res <- run_chains(model, 1L, seed, 1L, function(id) {
    lp <- lp_value_fun(model)
    mu <- rep(0, n)
    L <- diag(reference_scale, n)
    if (reference == "pathfinder") {
      pf <- tryCatch({
        lpg <- safe_lp(model, jacobian = TRUE)
        pathfinder_best_normal(lpg, find_init(model, 2, 1L), 5L, 1000L, 25L, 1e-12, 1e-8)
      }, error = function(e) NULL)
      if (!is.null(pf)) { mu <- pf$mu; L <- reference_scale * pf$L }
    }
    log_det_L <- sum(log(abs(diag(L))))
    log_g <- function(U) {
      Z <- t(forwardsolve(L, t(U) - mu))
      -0.5 * rowSums(Z^2) - 0.5 * n * log(2 * pi) - log_det_L
    }
    U <- sweep(matrix(stats::rnorm(N * n), N, n) %*% t(L), 2L, mu, "+")
    lg <- log_g(U)
    lpv <- eval_rows(lp, U, cores)
    beta <- 0
    log_z <- 0
    betas <- 0
    accs <- numeric(0)
    lw <- numeric(N)
    scale <- 2.38 / sqrt(n)
    for (stage in seq_len(max_stages)) {
      dl <- lpv - lg
      dl[!is.finite(dl)] <- -Inf
      ess_at <- function(b) {
        w <- (b - beta) * dl
        if (all(w == -Inf)) return(0)
        w <- exp(w - max(w))
        sum(w)^2 / sum(w^2) / N
      }
      nb <- if (ess_at(1) >= ess_threshold) 1 else {
        lo <- beta; hi <- 1
        for (k in 1:50) { mid <- (lo + hi) / 2; if (ess_at(mid) >= ess_threshold) lo <- mid else hi <- mid }
        max(lo, beta + 1e-8)
      }
      inc <- (nb - beta) * dl
      mx <- max(inc)
      if (!is.finite(mx)) stop("SMC: no particle has a finite log density", call. = FALSE)
      log_z <- log_z + mx + log(mean(exp(inc - mx)))
      w <- exp(inc - mx)
      beta <- nb
      betas <- c(betas, beta)
      idx <- systematic_resample(w / sum(w))
      U <- U[idx, , drop = FALSE]; lg <- lg[idx]; lpv <- lpv[idx]
      # random walk Metropolis moves targeting pi_beta
      C <- stats::cov(U)
      Lc <- tryCatch(t(chol(C + diag(1e-10 * mean(diag(C)) + 1e-12, n))), error = function(e) diag(sqrt(pmax(diag(C), 1e-12)), n))
      acc_n <- 0
      for (k in seq_len(num_mcmc_steps)) {
        P <- U + scale * matrix(stats::rnorm(N * n), N, n) %*% t(Lc)
        lgp <- log_g(P)
        lpp <- eval_rows(lp, P, cores)
        la <- beta * (lpp - lpv) + (1 - beta) * (lgp - lg)
        la[!is.finite(la)] <- -Inf
        ok <- log(stats::runif(N)) < la
        U[ok, ] <- P[ok, ]; lg[ok] <- lgp[ok]; lpv[ok] <- lpp[ok]
        acc_n <- acc_n + sum(ok)
      }
      acc <- acc_n / (N * num_mcmc_steps)
      accs <- c(accs, acc)
      scale <- scale * exp(acc - 0.234)
      if (beta >= 1) break
    }
    if (beta < 1) warning(sprintf("SMC stopped at beta = %.3g after %d stages; increase 'max_stages'", beta, max_stages), call. = FALSE)
    info <<- list(log_evidence = log_z, betas = betas, acceptance = accs)
    gq <- chain_write_draws(model, U, FALSE)
    dg <- matrix(0, N, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- lpv
    list(q = U, gq = gq, diag = dg, step_size = NA_real_, inv_metric = NULL,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  fit <- make_fit_extra(model, res, list(num_particles = N, reference = reference), "smc", seed, proc.time()[["elapsed"]] - t0)
  fit$log_evidence <- info$log_evidence
  fit$betas <- info$betas
  fit$acceptance <- info$acceptance
  fit
}

eval_rows <- function(f, U, cores = 1L) {
  if (cores > 1L && .Platform$OS.type != "windows") {
    ch <- split(seq_len(nrow(U)), cut(seq_len(nrow(U)), cores, labels = FALSE))
    unlist(parallel::mclapply(ch, function(ix) vapply(ix, function(i) f(U[i, ]), 0), mc.cores = cores), use.names = FALSE)
  } else {
    vapply(seq_len(nrow(U)), function(i) f(U[i, ]), 0)
  }
}

systematic_resample <- function(w) {
  N <- length(w)
  u <- (stats::runif(1) + 0:(N - 1L)) / N
  pmin(findInterval(u, cumsum(w)) + 1L, N)
}
