# Reference values from CmdStan 2.40 log_prob for the same program.

lap_data <- list(N = 4L, y = c(1L, 3L, 0L, 2L), b = c(1L, 0L, 0L, 1L), x = c(0.1, 0.7, 1.4, 2.2), idx = 1:4)

test_that("embedded Laplace approximations match CmdStan log_prob", {
  code <- "
functions {
  matrix kf(array[] real x, real a, real l) { return add_diag(gp_exp_quad_cov(x, a, l), 1e-6); }
  real lik(vector th, array[] int y, real s) { return poisson_log_lpmf(y | th * s); }
}
data { int N; array[N] int y; array[N] int b; array[N] real x; array[N] int idx; }
transformed data {
  print(laplace_latent_solve(lik, (y, 1.1), 1, kf, (x, 0.9, 0.8)));
}
parameters { real<lower=0> a; real<lower=0> l; real<lower=0> s; vector[N] m; }
model {
  a ~ normal(1, 0.5); l ~ normal(1, 0.5); s ~ normal(1, 0.2); m ~ normal(0, 1);
  target += laplace_marginal(lik, (y, s), 1, kf, (x, a, l));
  target += laplace_marginal_tol(lik, (y, s), 2, kf, (x, a, l), generate_laplace_options(N));
  y ~ laplace_marginal_poisson_log(idx, m, 1, kf, (x, a, l));
  target += laplace_marginal_bernoulli_logit_lpmf(b | idx, m, 1, kf, (x, a, l));
}
"
  out <- utils::capture.output(m <- bayescomp_model(code, data = lap_data))
  expect_identical(out, "([0.328628,0.40177,0.13211,0.225154],[[0.525263,0,0,0],[0.245683,0.409601,0,0],[-0.0146408,0.319172,0.417807,0],[-0.0128182,-0.0069292,0.289263,0.489739]])")
  u <- bayescomp_unconstrain(m, list(a = 0.9, l = 0.8, s = 1.1, m = c(0.2, -0.1, 0.3, 0.05)))
  g <- bayescomp_grad_log_prob(m, u, jacobian = TRUE, propto = TRUE)
  expect_equal(attr(g, "log_prob"), -26.211422, tolerance = 1e-7)
  expect_equal(as.numeric(g), c(-3.1074909, 0.12774017, -4.6359313, -0.43539394, 1.3707912, -2.2605004, 1.0999047),
               tolerance = 1e-5)
})

test_that("the negative binomial Laplace form has consistent gradients and latent draws", {
  code <- "
functions { matrix kf(array[] real x, real a, real l) { return add_diag(gp_exp_quad_cov(x, a, l), 1e-6); } }
data { int N; array[N] int y; array[N] real x; array[N] int idx; int E; }
parameters { real<lower=0> a; vector<lower=0>[E] eta; vector[N] m; }
model {
  a ~ normal(1, 0.5); eta ~ lognormal(0, 1); m ~ normal(0, 1);
  target += laplace_marginal_neg_binomial_2_log_lpmf(y | idx, eta, m, 1, kf, (x, a, 0.8));
}
generated quantities { vector[N] th = laplace_latent_neg_binomial_2_log_rng(y, idx, eta, m, 1, kf, (x, a, 0.8)); }
"
  m <- bayescomp_model(code, data = list(N = 4L, y = c(1L, 3L, 0L, 2L), x = c(0.1, 0.7, 1.4, 2.2), idx = 1:4, E = 4L))
  u <- c(0.3, 0.2, 0.4, 0.6, 0.8, 0.2, -0.1, 0.3, 0.05)
  g <- bayescomp_grad_log_prob(m, u)
  fd <- vapply(seq_along(u), function(k) {
    h <- 1e-5; up <- u; dn <- u; up[k] <- up[k] + h; dn[k] <- dn[k] - h
    (bayescomp_log_prob(m, up) - bayescomp_log_prob(m, dn)) / (2 * h)
  }, 0)
  expect_equal(as.numeric(g), fd, tolerance = 1e-5)
  draws <- vapply(1:400, function(i) bayescomp_constrain(m, u, seed = i)$th, numeric(4))
  expect_equal(dim(draws), c(4L, 400L))
  expect_true(all(abs(rowMeans(draws)) < 2))
})

test_that("Laplace calls are checked as in stanc", {
  f <- "functions { matrix kf(array[] real x, real a) { return add_diag(gp_exp_quad_cov(x, a, 1.0), 1e-6); }
                    real lik(vector th, array[] int y) { return poisson_log_lpmf(y | th); } }
        data { array[2] int y; array[2] real x; } parameters { real<lower=0> a; } "
  expect_error(bayescomp_check(paste0(f, "transformed parameters { tuple(vector[2], matrix[2, 2]) r = laplace_latent_solve(lik, (y,), 1, kf, (x, a)); }")),
               "data only", class = "bayescomp_semantic_error")
  expect_error(bayescomp_check(paste0(f, "model { vector[2] t = laplace_latent_rng(lik, (y,), 1, kf, (x, a)); }")),
               "random number generator", class = "bayescomp_semantic_error")
  expect_error(bayescomp_check(paste0(f, "model { target += laplace_marginal(lik, (y,), 1, kf, (x, x)); }")),
               "covariance function", class = "bayescomp_semantic_error")
  expect_silent(bayescomp_check(paste0(f, "model { target += laplace_marginal_tol(lik, (y,), 1, kf, (x, a), [0, 0]', 1e-8, 100, 1, 100, 1); }")))
})

test_that("short MCMC runs flag a poor approximation and pass an exact one", {
  skip_on_cran()
  m <- bayescomp_model("parameters { real<lower=0> s; real mu; } model { s ~ gamma(1.5, 1); mu ~ normal(0, s); }")
  r <- bayescomp_check_approximation(bayescomp_laplace(m, seed = 1), chains = 4, iter = 150, seed = 1)
  expect_true(any(r$flag))
  m2 <- bayescomp_model("parameters { vector[2] x; } model { x ~ normal(1, 2); }")
  r2 <- bayescomp_check_approximation(bayescomp_laplace(m2, seed = 1), chains = 4, iter = 150, seed = 1)
  expect_false(any(r2$flag))
})
