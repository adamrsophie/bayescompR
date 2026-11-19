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
