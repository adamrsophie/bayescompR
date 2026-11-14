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
#' @param num_mcmc_steps MCMC steps per particle and temperature.
#' @param kernel Rejuvenation kernel: `"metropolis"` (random walk, no
#'   gradients), `"mala"` or `"hmc"` (`hmc_steps` leapfrog steps), all
#'   preconditioned by the particle covariance and tuned toward their
#'   optimal acceptance rates between stages.
#' @param hmc_steps Leapfrog steps per HMC move.
#' @param resampling `"systematic"`, `"stratified"`, `"residual"` or
#'   `"multinomial"`.
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
                          kernel = c("metropolis", "mala", "hmc"), hmc_steps = 5L,
                          resampling = c("systematic", "stratified", "residual", "multinomial"),
                          max_stages = 200L, seed = NULL, cores = 1L) {
  check_model(model)
  reference <- match.arg(reference)
  kernel <- match.arg(kernel)
  resampling <- match.arg(resampling)
  resample <- switch(resampling, systematic = systematic_resample, stratified = stratified_resample,
                     residual = residual_resample, multinomial = multinomial_resample)
  N <- as_count(num_particles, "num_particles", 10L)
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  if (!(ess_threshold > 0 && ess_threshold < 1)) stop("'ess_threshold' must be in (0, 1)", call. = FALSE)
  t0 <- proc.time()[["elapsed"]]
  info <- NULL
  res <- run_chains(model, 1L, seed, 1L, function(id) {
    lp <- lp_value_fun(model)
    lpg_full <- model_lpg(model, jacobian = TRUE, propto = FALSE)
    ref <- gaussian_reference(model, reference, reference_scale)
    mu <- ref$mu; L <- ref$L; log_g <- ref$log_g
    U <- sweep(matrix(stats::rnorm(N * n), N, n) %*% t(L), 2L, mu, "+")
    lg <- log_g(U)
    lpv <- eval_rows(lp, U, cores)
    beta <- 0
    log_z <- 0
    betas <- 0
    accs <- numeric(0)
    lw <- numeric(N)
    scale <- switch(kernel, metropolis = 2.38 / sqrt(n), mala = 1 / n^(1 / 6), hmc = 1 / n^(1 / 4))
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
      idx <- resample(w / sum(w))
      U <- U[idx, , drop = FALSE]; lg <- lg[idx]; lpv <- lpv[idx]
      # random walk Metropolis moves targeting pi_beta
      C <- stats::cov(U)
      Lc <- tryCatch(t(chol(C + diag(1e-10 * mean(diag(C)) + 1e-12, n))), error = function(e) diag(sqrt(pmax(diag(C), 1e-12)), n))
      acc_n <- 0
      if (kernel != "metropolis") {
        mv <- smc_gradient_moves(U, lg, lpv, beta, mu, L, Lc, scale, kernel, hmc_steps, num_mcmc_steps, lpg_full, cores)
        U <- mv$U; lg <- mv$lg; lpv <- mv$lpv
        acc <- mv$acc
        accs <- c(accs, acc)
        target_acc <- if (kernel == "mala") 0.574 else 0.65
        scale <- scale * exp(acc - target_acc)
        if (beta >= 1) break
        next
      }
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
    complete <- beta >= 1
    if (!complete) warning(sprintf("SMC stopped at beta = %.3g after %d stages; the draws and log_evidence describe the tempered bridge, not the posterior; increase 'max_stages'", beta, max_stages), call. = FALSE)
    info <<- list(log_evidence = if (complete) log_z else NA_real_, partial_log_normalizer = log_z,
                  betas = betas, acceptance = accs, complete = complete)
    gq <- chain_write_draws(model, U, FALSE)
    dg <- matrix(0, N, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- lpv
    list(q = U, gq = gq, diag = dg, step_size = NA_real_, inv_metric = NULL,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  fit <- make_fit_extra(model, res, list(num_particles = N, reference = reference), "smc", seed, proc.time()[["elapsed"]] - t0)
  fit$log_evidence <- info$log_evidence
  fit$betas <- info$betas
  fit$complete <- info$complete
  fit$partial_log_normalizer <- info$partial_log_normalizer
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

stratified_resample <- function(w) {
  N <- length(w)
  u <- (stats::runif(N) + 0:(N - 1L)) / N
  pmin(findInterval(u, cumsum(w)) + 1L, N)
}
multinomial_resample <- function(w) sample.int(length(w), length(w), replace = TRUE, prob = w)
residual_resample <- function(w) {
  N <- length(w)
  k <- floor(N * w)
  idx <- rep(seq_len(N), k)
  R <- N - sum(k)
  if (R > 0L) {
    r <- N * w - k
    idx <- c(idx, sample.int(N, R, replace = TRUE, prob = r))
  }
  idx
}

# MALA or HMC moves for every particle targeting
# pi_beta = g^(1 - beta) p^beta, preconditioned by the particle covariance
# C = Lc Lc^T (the inverse metric).
smc_gradient_moves <- function(U, lg, lpv, beta, mu, L, Lc, eps, kernel, hmc_steps, steps, lpg, cores) {
  N <- nrow(U)
  C <- tcrossprod(Lc)
  Sinv_apply <- function(x) backsolve(t(L), forwardsolve(L, x))
  log_g1 <- function(u) {
    z <- forwardsolve(L, u - mu)
    -0.5 * sum(z^2) - 0.5 * length(u) * log(2 * pi) - sum(log(abs(diag(L))))
  }
  tgt <- function(u) {
    r <- tryCatch(lpg(u), bayescomp_fatal_error = function(e) stop(e), bayescomp_error = function(e) NULL)
    if (is.null(r) || !is.finite(r$value) || any(!is.finite(r$grad))) return(list(V = Inf, lp = -Inf, lg = log_g1(u), g = rep(0, length(u))))
    lgu <- log_g1(u)
    grad <- beta * r$grad - (1 - beta) * Sinv_apply(u - mu)
    list(V = -(beta * r$value + (1 - beta) * lgu), lp = r$value, lg = lgu, g = -grad)
  }
  one <- function(i) {
    u <- U[i, ]
    st <- tgt(u)
    acc_sum <- 0
    for (k in seq_len(steps)) {
      if (kernel == "mala") {
        m0 <- u - 0.5 * eps^2 * as.vector(C %*% st$g)
        un <- m0 + eps * as.vector(Lc %*% stats::rnorm(length(u)))
        sn <- tgt(un)
        a <- 0
        if (is.finite(sn$V)) {
          m1 <- un - 0.5 * eps^2 * as.vector(C %*% sn$g)
          lf <- -0.5 * sum(forwardsolve(Lc, un - m0)^2) / eps^2
          lb <- -0.5 * sum(forwardsolve(Lc, u - m1)^2) / eps^2
          la <- st$V - sn$V + lb - lf
          a <- if (is.nan(la)) 0 else min(1, exp(la))
        }
      } else {
        p <- backsolve(t(Lc), stats::rnorm(length(u)))
        H0 <- st$V + 0.5 * sum(p * as.vector(C %*% p))
        sn <- st; un <- u
        for (l in seq_len(hmc_steps)) {
          p <- p - 0.5 * eps * sn$g
          un <- un + eps * as.vector(C %*% p)
          sn <- tgt(un)
          if (!is.finite(sn$V)) break
          p <- p - 0.5 * eps * sn$g
        }
        H1 <- sn$V + 0.5 * sum(p * as.vector(C %*% p))
        a <- if (is.finite(H1)) min(1, exp(H0 - H1)) else 0
      }
      if (stats::runif(1) < a) { u <- un; st <- sn }
      acc_sum <- acc_sum + a
    }
    c(u, st$lg, st$lp, acc_sum / steps)
  }
  rows <- if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(N), one, mc.cores = cores)
  } else lapply(seq_len(N), one)
  M <- do.call(rbind, rows)
  n <- ncol(U)
  list(U = M[, seq_len(n), drop = FALSE], lg = M[, n + 1L], lpv = M[, n + 2L], acc = mean(M[, n + 3L]))
}

# A normalized Gaussian reference on the unconstrained space: the best
# Pathfinder normal with its covariance scaled by scale^2, or N(0, scale^2 I).
gaussian_reference <- function(model, reference, scale) {
  n <- model$dim
  mu <- rep(0, n)
  L <- diag(scale, n)
  if (reference == "pathfinder") {
    pf <- tryCatch({
      lpg <- safe_lp(model, jacobian = TRUE)
      pathfinder_best_normal(lpg, find_init(model, 2, 1L), 5L, 1000L, 25L, 1e-12, 1e-8)
    }, error = function(e) NULL)
    if (!is.null(pf)) { mu <- pf$mu; L <- scale * pf$L }
  }
  log_det_L <- sum(log(abs(diag(L))))
  log_g <- function(U) {
    U <- matrix(U, ncol = n)
    Z <- t(forwardsolve(L, t(U) - mu))
    -0.5 * rowSums(Z^2) - 0.5 * n * log(2 * pi) - log_det_L
  }
  list(mu = mu, L = L, log_g = log_g)
}
