test_that("multivariate ESS is near the draw count for independent draws and lower under autocorrelation", {
  set.seed(1)
  iid <- array(rnorm(8000), c(1000, 4, 2))
  expect_equal(bayescomp_ess_multi(iid)$mess, 4000, tolerance = 0.3)
  ar <- array(0, c(1000, 4, 2))
  for (ch in 1:4) for (j in 1:2) ar[, ch, j] <- stats::arima.sim(list(ar = 0.9), 1000)
  r <- bayescomp_ess_multi(ar)
  expect_lt(r$mess, 1000)
  expect_false(r$enough)
  # Vats, Flegal and Jones (2019): about 6146 for one variable
  expect_equal(bayescomp_ess_multi(iid[, , 1, drop = FALSE])$min_ess, 6146, tolerance = 0.01)
})

test_that("nested R-hat flags superchains that disagree", {
  set.seed(1)
  sc <- rep(1:8, each = 8)
  expect_lt(bayescomp_rhat_nested(matrix(rnorm(640), 10, 64), sc), 1.05)
  y <- matrix(rnorm(640), 10, 64) + rep(rep(c(0, 3), each = 32), each = 10)
  expect_gt(bayescomp_rhat_nested(y, sc), 1.5)
  expect_error(bayescomp_rhat_nested(matrix(0, 1, 3), 1:3 > 1), "equal size")
})

test_that("local R-hat flags chains that differ only in shape", {
  set.seed(1)
  x <- cbind(rnorm(2000), runif(2000, -sqrt(3), sqrt(3)))
  # Moins et al. (2023): about 1.018 in the limit for this pair
  expect_gt(bayescomp_rhat_local(x), 1.01)
  expect_lt(bayescomp_rhat_local(matrix(rnorm(4000), 2000)), 1.01)
})

test_that("kernel Stein discrepancy grows when the draws are biased", {
  skip_on_cran()
  m <- bayescomp_model("parameters { vector[2] mu; } model { mu ~ normal(0, 1); }")
  f <- bayescomp_sample(m, chains = 1, iter_sampling = 500, seed = 1, refresh = 0)
  a <- bayescomp_ksd(f)
  f$unconstrained[[1L]] <- f$unconstrained[[1L]] + 0.5
  expect_gt(bayescomp_ksd(f), 3 * a)
})

test_that("R* is near 1 for mixed chains and larger when one chain differs", {
  skip_if_not_installed("ranger")
  set.seed(1)
  x <- array(rnorm(8000), c(1000, 4, 2))
  expect_lt(bayescomp_rstar(x), 1.15)
  x[, 4, 1] <- 2 * x[, 4, 1]
  expect_gt(bayescomp_rstar(x), 1.15)
  expect_length(bayescomp_rstar(x, uncertainty = TRUE, nsim = 50), 50)
})
