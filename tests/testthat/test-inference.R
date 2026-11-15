test_that("NUTS recovers a known posterior", {
  skip_on_cran()
  m <- bayescomp_model("parameters { vector[2] x; real<lower=0> s; } model { x ~ normal(1, 2); s ~ gamma(3, 2); }")
  n <- if (long_tests()) 1000L else 300L
  fit <- bayescomp_sample(m, chains = 2, iter_warmup = n, iter_sampling = n, seed = 11, refresh = 0)
  s <- summary(fit)
  mx <- s[s$variable == "x[1]", ]
  expect_lt(abs(mx$mean - 1), 5 * 2 / sqrt(mx$ess_bulk))
  ms <- s[s$variable == "s", ]
  expect_lt(abs(ms$mean - 1.5), 5 * sqrt(3) / 2 / sqrt(ms$ess_bulk))
  expect_true(all(s$rhat < 1.05))
  expect_equal(bayescomp_diagnose(fit, quiet = TRUE)$num_divergent, 0L)
  # results do not depend on the global seed
  set.seed(99)
  fit2 <- bayescomp_sample(m, chains = 2, iter_warmup = 50, iter_sampling = 20, seed = 11, refresh = 0)
  set.seed(1)
  fit3 <- bayescomp_sample(m, chains = 2, iter_warmup = 50, iter_sampling = 20, seed = 11, refresh = 0)
  expect_identical(fit2$draws, fit3$draws)
})

test_that("dense metric, static HMC and fixed parameter models run", {
  skip_on_cran()
  m <- bayescomp_model("parameters { vector[2] x; } model { x ~ multi_normal([0, 0]', [[1, 0.8], [0.8, 1]]); }")
  fit <- bayescomp_sample(m, chains = 1, iter_warmup = 200, iter_sampling = 200, metric = "dense_e", seed = 2, refresh = 0)
  expect_equal(dim(fit$inv_metric[[1L]]), c(2L, 2L))
  fit <- bayescomp_sample(m, chains = 1, iter_warmup = 100, iter_sampling = 100, algorithm = "hmc", seed = 2, refresh = 0)
  expect_equal(dim(fit$draws)[1L], 100L)
  g <- bayescomp_model("generated quantities { real y = normal_rng(0, 1); }")
  fit <- bayescomp_sample(g, chains = 1, iter_sampling = 50, seed = 1, refresh = 0)
  expect_equal(dim(fit$draws)[1L], 50L)
})

test_that("optimization, Laplace, ADVI and Pathfinder agree with the analytic answer", {
  skip_on_cran()
  y <- seq(-1, 3, length.out = 20)
  m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; real<lower=0> sigma; }
    model { y ~ normal(mu, sigma); }", data = list(N = 20, y = y))
  o <- bayescomp_optimize(m, seed = 1)
  expect_equal(unname(o$par["mu"]), mean(y), tolerance = 1e-5)
  expect_equal(unname(o$par["sigma"]), sqrt(mean((y - mean(y))^2)), tolerance = 1e-4)
  on <- bayescomp_optimize(m, algorithm = "newton", seed = 1)
  expect_equal(unname(on$par), unname(o$par), tolerance = 1e-4)
  la <- bayescomp_laplace(m, seed = 1, draws = 2000)
  expect_equal(mean(bayescomp_extract(la, "mu")), mean(y), tolerance = 0.05)
  vi <- bayescomp_variational(m, seed = 1)
  expect_equal(mean(bayescomp_extract(vi, "mu")), mean(y), tolerance = 0.1)
  pf <- bayescomp_pathfinder(m, seed = 1)
  expect_equal(mean(bayescomp_extract(pf, "mu")), mean(y), tolerance = 0.1)
})

test_that("diagnostics agree with the posterior package", {
  skip_on_cran()
  skip_if_not_installed("posterior")
  set.seed(5)
  x <- matrix(stats::rnorm(800), 200, 4) + rep(c(0, 0.1, 0, 0.3), each = 200)
  x[, 1] <- cumsum(x[, 1]) / 10
  expect_equal(bayescomp_rhat(x), posterior::rhat(x), tolerance = 1e-8)
  expect_equal(bayescomp_ess_bulk(x), posterior::ess_bulk(x), tolerance = 1e-8)
  expect_equal(bayescomp_ess_tail(x), posterior::ess_tail(x), tolerance = 1e-8)
  expect_equal(bayescomp_ess_basic(x), posterior::ess_basic(x), tolerance = 1e-8)
  expect_equal(bayescomp_mcse_mean(x), posterior::mcse_mean(x), tolerance = 1e-8)
  expect_equal(bayescomp_mcse_sd(x), posterior::mcse_sd(x), tolerance = 1e-6)
  expect_equal(bayescomp_mcse_quantile(x, c(0.05, 0.5, 0.95)), unname(posterior::mcse_quantile(x, c(0.05, 0.5, 0.95))), tolerance = 1e-8)
})
