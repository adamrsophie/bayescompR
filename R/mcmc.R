# Further MCMC kernels sharing the HMC chain driver (warmup windows, dual
# averaging and metric estimation in run_hmc_chain):
#
# * random walk Metropolis, block or componentwise. The proposal covariance
#   is step_size^2 times the inverse metric, so the windowed metric
#   estimation gives the adaptive Metropolis covariance of Haario, Saksman
#   and Tamminen (2001), restricted to warmup so the sampling phase is a
#   fixed Markov kernel. Dual averaging tunes the scale to `adapt_delta`
#   (default 0.234 for blocks, 0.44 componentwise).
# * the Metropolis adjusted Langevin algorithm (Roberts and Tweedie 1996)
#   with the same preconditioning (default target acceptance 0.574).
# * slice sampling with stepping out and shrinkage (Neal 2003), one update
#   per direction of the Cholesky factor of the inverse metric. Slice
#   sampling has no acceptance statistic; the bracket width is step_size
#   times the metric scale and only the metric is adapted.
#
# All kernels report the log density in lp__, the acceptance statistic in
# accept_stat__ and the number of density evaluations in n_leapfrog__.

# log density without the gradient: the tape is disabled. Without a tape
# every value is a constant, so propto would drop every sampling statement;
# the full density is evaluated and, in the chain driver, shifted by a
# constant so that lp__ matches the propto value used by HMC. Samplers pass
# propto = TRUE: the two conventions differ by more than a constant when a
# program reads target() or drops a data only statement in a parameter
# dependent branch, so samplers use Stan's (constants dropped; values that
# depend on parameters keep their class, but nothing is recorded). Evidence and model probability code needs
# the normalized density and keeps the default.
lp_value_fun <- function(model, print = FALSE, jacobian = TRUE, propto = FALSE) {
  if (!is.null(model$native)) {
    nat <- model$native
    fv <- nat[["log_density"]] %||% function(u, jacobian, propto) nat$log_density_gradient(u, jacobian, propto)$value
    return(function(q) {
      # the same conditions as the R path below: only domain errors and
      # rejections give -Inf, everything else aborts
      r <- tryCatch(fv(q, jacobian, propto), bayescomp_fatal_error = function(e) stop(e),
                    bayescomp_domain_error = function(e) -Inf,
                    bayescomp_reject_error = function(e) -Inf)
      if (length(r) != 1L || is.nan(r)) -Inf else r
    })
  }
  f <- model_fun(model, "log_prob")
  g <- if (propto) propto_value_fun(model, f, jacobian, print) else function(q) {
    old <- .ad$tape
    .ad$tape <- NULL
    on.exit(.ad$tape <- old)
    as.numeric(val(f(q, new_ctx(jacobian = jacobian, print = print, propto = FALSE))))
  }
  function(q) {
    r <- tryCatch(g(q), bayescomp_fatal_error = function(e) stop(e),
       bayescomp_domain_error = function(e) -Inf,
       bayescomp_reject_error = function(e) -Inf)
    if (length(r) != 1L || is.nan(r)) -Inf else r
  }
}

# lower triangular square root of the inverse metric applied to z
metric_chol_mult <- function(m, z) {
  if (m$type == "diag_e") sqrt(m$inv) * z else as.vector(crossprod(m$U, z))
}
metric_chol_solve <- function(m, x) {
  if (m$type == "diag_e") x / sqrt(m$inv) else as.vector(backsolve(m$U, x, transpose = TRUE))
}
metric_direction <- function(m, i) {
  if (m$type == "diag_e") { d <- numeric(length(m$inv)); d[i] <- sqrt(m$inv[i]); d } else m$U[i, ]
}

mh_out <- function(q, lp, acc, evals) {
  list(q = q, lp = lp, accept_stat = acc, treedepth = 0L, n_leapfrog = evals, divergent = FALSE, energy = NA_real_)
}

metropolis_transition <- function(st, q) {
  m <- st$metric
  lp0 <- st$lp_cache(q)
  n <- length(q)
  if (st$block) {
    qn <- q + st$eps * metric_chol_mult(m, stats::rnorm(n))
    lp1 <- st$lp(qn)
    acc <- if (is.finite(lp1)) min(1, exp(lp1 - lp0)) else 0
    if (stats::runif(1) < acc) { q <- qn; lp0 <- lp1 }
    st$last <- list(q = q, lp = lp0)
    return(mh_out(q, lp0, acc, 1L))
  }
  sd <- if (m$type == "diag_e") sqrt(m$inv) else sqrt(diag(m$inv))
  accs <- numeric(n)
  for (i in seq_len(n)) {
    qn <- q
    qn[i] <- q[i] + st$eps * sd[i] * stats::rnorm(1)
    lp1 <- st$lp(qn)
    accs[i] <- if (is.finite(lp1)) min(1, exp(lp1 - lp0)) else 0
    if (stats::runif(1) < accs[i]) { q <- qn; lp0 <- lp1 }
  }
  st$last <- list(q = q, lp = lp0)
  mh_out(q, lp0, mean(accs), n)
}

mala_transition <- function(st, q) {
  m <- st$metric
  e <- st$eps
  # the retained state of the previous transition is where this one starts
  cached <- identical(st$last_q, q)
  r0 <- if (cached) st$last_r else st$pot(q)
  drift <- function(r, q) q - 0.5 * e^2 * metric_dtau_cov(m, r$g)
  mu0 <- drift(r0, q)
  qn <- mu0 + e * metric_chol_mult(m, stats::rnorm(length(q)))
  r1 <- st$pot(qn)
  acc <- 0
  if (is.finite(r1$V)) {
    mu1 <- drift(r1, qn)
    lq_fwd <- -0.5 * sum(metric_chol_solve(m, qn - mu0)^2) / e^2
    lq_bwd <- -0.5 * sum(metric_chol_solve(m, q - mu1)^2) / e^2
    la <- (r0$V - r1$V) + lq_bwd - lq_fwd
    acc <- if (is.nan(la)) 0 else min(1, exp(la))
  }
  if (stats::runif(1) < acc) { q <- qn; r0 <- r1 }
  st$last_q <- q
  st$last_r <- r0
  mh_out(q, -r0$V, acc, if (cached) 1L else 2L)
}
# covariance (inverse metric) times a vector
metric_dtau_cov <- function(m, g) if (m$type == "diag_e") m$inv * g else as.vector(m$inv %*% g)

slice_transition <- function(st, q) {
  m <- st$metric
  lp0 <- st$lp_cache(q)
  n <- length(q)
  evals <- 0L
  for (i in seq_len(n)) {
    d <- metric_direction(m, i)
    logy <- lp0 - stats::rexp(1)
    w <- st$eps
    L <- -w * stats::runif(1)
    R <- L + w
    # stepping out, at most max_steps expansions in total (Neal 2003, Fig. 3)
    J <- floor(st$max_steps * stats::runif(1))
    K <- st$max_steps - 1L - J
    while (J > 0 && st$lp(q + L * d) > logy) { L <- L - w; J <- J - 1; evals <- evals + 1L }
    while (K > 0 && st$lp(q + R * d) > logy) { R <- R + w; K <- K - 1; evals <- evals + 1L }
    repeat {
      t <- stats::runif(1, L, R)
      qn <- q + t * d
      lp1 <- st$lp(qn)
      evals <- evals + 1L
      if (lp1 > logy) { q <- qn; lp0 <- lp1; break }
      if (t < 0) L <- t else R <- t
      if (R - L < 1e-12) break
    }
  }
  st$last <- list(q = q, lp = lp0)
  mh_out(q, lp0, 1, evals)
}

# reuse the log density of the current state between transitions
lp_cache_fun <- function(st) {
  function(q) {
    if (!is.null(st$last) && identical(st$last$q, q)) st$last$lp else st$lp(q)
  }
}

default_adapt_delta <- function(algorithm, block) {
  switch(algorithm, nuts = 0.8, hmc = 0.8, mala = 0.574, metropolis = if (block) 0.234 else 0.44, NA_real_)
}

# Generalized elliptical slice sampling (Murray, Adams and MacKay 2010;
# Nishihara, Murray and Adams 2014): the target is written as
# N(q; c, Sigma) * L(q) with L = p / N, where the Gaussian has the adapted
# inverse metric as covariance and the warmup mean as center. Elliptical
# slice sampling on L is exact for any target; it is efficient when the
# Gaussian matches the posterior. It needs no gradient and no step size.
elliptical_transition <- function(st, q) {
  m <- st$metric
  c0 <- st$center
  logN <- function(x) -0.5 * sum(metric_chol_solve(m, x - c0)^2)
  lp0 <- st$lp_cache(q)
  ll0 <- lp0 - logN(q)
  logy <- ll0 + log(stats::runif(1))
  nu <- metric_chol_mult(m, stats::rnorm(length(q)))
  th <- stats::runif(1, 0, 2 * pi)
  lo <- th - 2 * pi
  hi <- th
  evals <- 0L
  repeat {
    qn <- c0 + (q - c0) * cos(th) + nu * sin(th)
    lp1 <- st$lp(qn)
    evals <- evals + 1L
    if (is.finite(lp1) && lp1 - logN(qn) > logy) { q <- qn; lp0 <- lp1; break }
    if (th < 0) lo <- th else hi <- th
    if (hi - lo < 1e-12) break
    th <- stats::runif(1, lo, hi)
  }
  st$last <- list(q = q, lp = lp0)
  mh_out(q, lp0, 1, evals)
}
