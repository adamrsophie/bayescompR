test_that("pointwise log likelihood matches a generated quantity", {
  skip_on_cran()
  set.seed(1); N <- 8; x <- rnorm(N); y <- 1 + x + rnorm(N); k <- rpois(N, 3)
  code <- "data { int N; vector[N] x; vector[N] y; array[N] int k; } parameters { real a; real<lower=0> s; real<lower=0> lam; }
    model { a ~ normal(0, 5); s ~ exponential(1); for (n in 1:N) y[n] ~ normal(a + x[n], s); target += poisson_lpmf(k | lam); }
    generated quantities { vector[N] log_lik; for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | a + x[n], s) + poisson_lpmf(k[n] | lam); }"
  m <- bayescomp_model(code, data = list(N = N, x = x, y = y, k = k))
  fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 50, seed = 1, refresh = 0)
  ll <- bayescomp_log_lik(fit)
  expect_equal(dim(ll), c(100L, 2L * N))
  expect_equal(colnames(ll)[c(1, N + 1)], c("y[1]", "k[1]"))
  expect_equal(unname(ll[, 1:N] + ll[, N + 1:N]), unname(draws_matrix_of(fit, "log_lik")))
  expect_equal(dim(bayescomp_log_lik(fit, "array")), c(50L, 2L, 2L * N))
  yr <- bayescomp_posterior_predict(fit, seed = 2)
  expect_equal(dim(yr), dim(ll))
  expect_true(all(yr[, N + 1:N] == round(yr[, N + 1:N])))
  expect_equal(yr, bayescomp_posterior_predict(fit, seed = 2))
})

test_that("models without data sampling statements give an error", {
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
  fit <- bayescomp_sample(m, chains = 1, iter_warmup = 20, iter_sampling = 5, seed = 1, refresh = 0)
  expect_error(bayescomp_log_lik(fit), "no sampling statement on data")
})

test_that("loo, log_lik and posterior_predict fall back to the model block", {
  skip_on_cran()
  skip_if_not_installed("loo")
  m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; } model { y ~ normal(mu, 1); }",
                       data = list(N = 4, y = c(-1, 0, 1, 2)))
  fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100, seed = 1, refresh = 0)
  lo <- suppressWarnings(loo::loo(fit))
  expect_equal(nrow(lo$pointwise), 4L)
  expect_equal(log_lik.bayescomp_fit(fit), bayescomp_log_lik(fit))
  expect_equal(dim(posterior_predict.bayescomp_fit(fit)), c(200L, 4L))
  expect_error(loo.bayescomp_fit(fit, variable = "ll"), "unknown variable")
})
