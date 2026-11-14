lgssm <- function() {
  set.seed(1)
  Tn <- 40
  x <- numeric(Tn); x[1] <- rnorm(1)
  for (t in 2:Tn) x[t] <- 0.8 * x[t - 1] + rnorm(1, 0, 0.5)
  y <- x + rnorm(Tn, 0, 0.7)
  kf <- function(phi, q, r = 0.7) {
    m <- 0; P <- 1; ll <- 0
    for (t in seq_len(Tn)) {
      if (t > 1) { m <- phi * m; P <- phi^2 * P + q^2 }
      S <- P + r^2
      ll <- ll + dnorm(y[t], m, sqrt(S), log = TRUE)
      K <- P / S; m <- m + K * (y[t] - m); P <- (1 - K) * P
    }
    ll
  }
  ssm <- bayescomp_ssm(y,
    init = function(N, th) rnorm(N),
    transition = function(x, t, th) th[["phi"]] * x + rnorm(length(x), 0, th[["q"]]),
    log_obs = function(yt, x, t, th) dnorm(yt, x, 0.7, log = TRUE),
    log_init = function(x, th) dnorm(x, log = TRUE),
    log_transition = function(xn, xo, t, th) dnorm(xn, th[["phi"]] * xo, th[["q"]], log = TRUE))
  list(ssm = ssm, kf = kf)
}

test_that("the particle filter likelihood matches the Kalman filter", {
  s <- lgssm()
  set.seed(2)
  ll <- replicate(20, bayescomp_particle_filter(s$ssm, c(phi = 0.8, q = 0.5), 300)$log_lik)
  expect_equal(mean(ll), s$kf(0.8, 0.5), tolerance = 0.5)
  pf <- bayescomp_particle_filter(s$ssm, c(phi = 0.8, q = 0.5), 200, resampling = "stratified", ess_threshold = 0.5)
  expect_equal(dim(pf$filtered_mean), c(40L, 1L))
  expect_true(all(pf$ess > 0))
})

test_that("particle MCMC recovers the parameter posterior", {
  s <- lgssm()
  m <- bayescomp_model("parameters { real<lower=-1, upper=1> phi; real<lower=0> q; }
    model { phi ~ normal(0, 1); q ~ normal(0, 1); }")
  g <- expand.grid(phi = seq(-0.99, 0.99, length.out = 80), q = seq(0.01, 1.5, length.out = 80))
  lp <- mapply(function(p, q) s$kf(p, q) + dnorm(p, log = TRUE) + dnorm(q, log = TRUE), g$phi, g$q)
  w <- exp(lp - max(lp)); w <- w / sum(w)
  for (meth in c("pmmh", "pgibbs")) {
    fit <- bayescomp_pmcmc(m, s$ssm, meth, num_particles = 80, chains = 2, iter_warmup = 200,
                           iter_sampling = 500, seed = 1)
    expect_equal(mean(bayescomp_extract(fit, "phi")), sum(w * g$phi), tolerance = 0.12)
    expect_equal(mean(bayescomp_extract(fit, "q")), sum(w * g$q), tolerance = 0.12)
    expect_equal(dim(fit$states), c(500L, 40L, 1L))
  }
  bad <- bayescomp_ssm(1:3, function(N, th) rnorm(N), function(x, t, th) x, function(y, x, t, th) 0)
  expect_error(bayescomp_pmcmc(m, bad, "pgibbs"), "log_transition")
})
