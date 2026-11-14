norm_model <- function(s, y = c(0.5, 1.2, 2, 0.8, 1.6)) {
  bayescomp_model(sprintf("data { int N; vector[N] y; } parameters { real mu; }
    model { mu ~ normal(0, %g); y ~ normal(mu, 1); }", s), data = list(N = length(y), y = y))
}
exact_lz <- function(s, y = c(0.5, 1.2, 2, 0.8, 1.6)) {
  S <- diag(length(y)) + s^2
  as.numeric(-0.5 * determinant(2 * pi * S)$modulus - 0.5 * sum(y * solve(S, y)))
}

test_that("nested sampling estimates the evidence", {
  fit <- bayescomp_nested(norm_model(2), num_live = 200, seed = 1)
  expect_equal(fit$log_evidence, exact_lz(2), tolerance = 4 * fit$log_evidence_se + 0.05)
  expect_equal(mean(bayescomp_extract(fit, "mu")), 6.1 / 5.25, tolerance = 0.1)
})

test_that("reversible jump recovers posterior model probabilities", {
  r <- bayescomp_rjmcmc(list(narrow = norm_model(1), wide = norm_model(10)), iter_sampling = 3000, seed = 1)
  z <- c(exact_lz(1), exact_lz(10))
  expect_equal(r$probabilities$posterior, exp(z - max(z)) / sum(exp(z - max(z))), tolerance = 0.05)
  expect_s3_class(r$fits$narrow, "bayescomp_fit")
  expect_output(print(r), "Reversible jump")
})

test_that("EP and SG-MCMC fit a linear regression", {
  set.seed(1)
  N <- 120
  X <- cbind(1, rnorm(N))
  y <- as.vector(X %*% c(1, 2) + rnorm(N))
  m <- bayescomp_model("data { int N; matrix[N, 2] X; vector[N] y; } parameters { vector[2] b; }
    model { b ~ normal(0, 10); y ~ normal(X * b, 1); }", data = list(N = N, X = X, y = y))
  P <- crossprod(X) + diag(1 / 100, 2)
  mu <- as.vector(solve(P, crossprod(X, y)))
  ep <- bayescomp_ep(m, sites = 3, seed = 1)
  expect_equal(ep$mean, mu, tolerance = 1e-3)
  expect_equal(sqrt(diag(ep$cov)), sqrt(diag(solve(P))), tolerance = 1e-2)
  for (meth in c("sgld", "sghmc")) {
    f <- bayescomp_sgmcmc(m, meth, batch_size = 20, chains = 1, iter_warmup = 200, iter_sampling = 1500, seed = 1)
    expect_equal(colMeans(bayescomp_extract(f, "b")), mu, tolerance = 0.1)
  }
  expect_error(bayescomp_ep(m, n_name = "M"), "observation count")
})
