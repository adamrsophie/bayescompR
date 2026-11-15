test_that("simulation based calibration returns ranks in range", {
  skip_on_cran()
  code <- "data { int N; vector[N] y; } parameters { real mu; } model { mu ~ normal(0, 1); y ~ normal(mu, 1); }"
  sim <- function(i) { mu <- rnorm(1); list(parameters = list(mu = mu), data = list(N = 3, y = rnorm(3, mu))) }
  x <- bayescomp_sbc(code, sim, n_sims = 8, max_rank = 19, seed = 1, chains = 1, iter_warmup = 100, iter_sampling = 100, refresh = 0)
  expect_s3_class(x, "bayescomp_sbc")
  expect_equal(nrow(x$ranks), 8L)
  expect_true(all(x$ranks$rank >= 0 & x$ranks$rank <= 19))
  expect_equal(x$summary$variable, "mu")
  expect_output(print(x), "8 simulations")
  expect_error(bayescomp_sbc(code, function(i) list(parameters = list(b = 2), data = list(N = 1, y = 0)), n_sims = 1, max_rank = 9,
                             chains = 1, iter_warmup = 20, iter_sampling = 20, refresh = 0), "do not match")
  expect_error(bayescomp_sbc(code, sim, n_sims = 1, max_rank = 99, chains = 1, iter_warmup = 20, iter_sampling = 20, refresh = 0), "fewer than")
})
