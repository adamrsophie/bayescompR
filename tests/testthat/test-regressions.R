# regressions for numerical and semantic edge cases

test_that("log_diff_exp, boundary truncation and tail CDF gradients are finite", {
  m <- bayescomp_model("generated quantities { real a = log_diff_exp(0, negative_infinity()); }")
  expect_equal(bayescomp_constrain(m, numeric(0))$a, 0)
  tr <- bayescomp_model("data { real y; } parameters { real<lower=0> lam; }
                         model { y ~ exponential(lam) T[0, 1]; lam ~ gamma(2, 1); }", data = list(y = 0.4))
  fd <- (bayescomp_log_prob(tr, 0.2 + 1e-6) - bayescomp_log_prob(tr, 0.2 - 1e-6)) / 2e-6
  expect_equal(as.vector(bayescomp_grad_log_prob(tr, 0.2)), fd, tolerance = 1e-5)
  tl <- bayescomp_model("parameters { real<lower=0> l; real<lower=0, upper=1> p; }
                         model { target += poisson_lcdf(0 | l) + binomial_lcdf(0 | 1000, p); }")
  g <- as.vector(bayescomp_grad_log_prob(tl, c(log(1000), stats::qlogis(0.9))))
  # chain rule through the transforms: d/du = -lambda + 1 and -1000 theta (1 - theta) / (1 - theta) + Jacobian
  expect_equal(g, c(-999, -900.8), tolerance = 1e-6)
})

test_that("repeated assignment indices give the adjoint to the last write only", {
  m <- bayescomp_model("parameters { vector[2] x; }
    model { vector[2] z = rep_vector(0, 2); z[{1, 1}] = x; target += z[1] - 0.5 * dot_self(x); }")
  expect_equal(as.vector(bayescomp_grad_log_prob(m, c(0, 0))), c(0, 1))
})

test_that("shapes: transformed data sizes, singleton arrays and size mismatches", {
  m <- bayescomp_model("transformed data { int K = 2; vector[K] v = rep_vector(1, K); }
                        parameters { real a; } model { a ~ normal(sum(v), 1); }")
  expect_equal(bayescomp_log_prob(m, 2), stats::dnorm(0, log = TRUE))
  a <- bayescomp_model("transformed data { array[2] vector[1] a = {rep_vector(1, 1), rep_vector(2, 1)}; }
                        generated quantities { real s = a[2][1]; }")
  expect_equal(bayescomp_constrain(a, numeric(0))$s, 2)
  expect_error(bayescomp_model("transformed data { vector[1] x = rep_vector(1, 1); vector[2] z = x + rep_vector(1, 2); }"),
               class = "bayescomp_shape_error")
})

test_that("CDFs check their parameters and impossible bounds are rejected", {
  m <- bayescomp_model("generated quantities { real a = normal_lcdf(0 | 0, -1); }")
  expect_error(bayescomp_constrain(m, numeric(0)), "Scale parameter", class = "bayescomp_domain_error")
  ok <- bayescomp_model("generated quantities { real b = poisson_lccdf(-1 | 3); }")
  expect_equal(bayescomp_constrain(ok, numeric(0))$b, 0)
  ib <- bayescomp_model("parameters { real<lower=positive_infinity()> a; } model { }")
  expect_error(bayescomp_log_prob(ib, 0), class = "bayescomp_domain_error")
})

test_that("diagnostics match the posterior package", {
  expect_equal(bayescomp_ess_basic(rep(c(-1, 1), 10)), 10)
  set.seed(3)
  y <- matrix(stats::rnorm(2000), 500)
  expect_equal(bayescomp_rhat(y * 1e-20), bayescomp_rhat(y))
  skip_if_not_installed("posterior")
  x <- y + matrix(rep(cumsum(stats::rnorm(500)) * 0.05, 4), 500)
  z <- matrix(stats::rbinom(2000, 1, 0.3), 500)
  a <- matrix(rep(c(-1, 1), 1000), 500)
  for (f in c("rhat", "ess_bulk", "ess_tail", "ess_basic", "mcse_mean", "mcse_sd"))
    for (m in list(x, z, a))
      expect_equal(get(paste0("bayescomp_", f))(m), get(f, asNamespace("posterior"))(m), label = f)
})

test_that("unequal minibatches give an unbiased stochastic gradient", {
  # N = 3, batch size 2: batches of sizes 2 and 1 used equally often
  m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; }
                        model { mu ~ normal(0, 1); y ~ normal(mu, 1); }", data = list(N = 3, y = c(0, 0, 3)))
  f <- bayescomp_sgmcmc(m, "sgld", batch_size = 2, step_size = 0.01, chains = 1, iter_warmup = 500,
                        iter_sampling = 8000, seed = 1)
  expect_equal(mean(bayescomp_extract(f, "mu")), 0.75, tolerance = 0.15)
})

test_that("particle filter increments ignore zero weight particles", {
  expect_equal(pf_increment(c(-Inf, 0), c(-Inf, 0)), 0)
})

test_that("warmup only fits and impossible diffusion horizons", {
  m <- bayescomp_model("parameters { real x; } model { x ~ normal(0, 1); }")
  f <- bayescomp_sample(m, chains = 1, iter_warmup = 20, iter_sampling = 20, seed = 1, refresh = 0)
  f$n_warmup_saved <- dim(f$draws)[1L]
  expect_length(kept(f), 0L)
  expect_error(bayescomp_diffusion(m, horizon = -1), "horizon")
  d <- bayescomp_diffusion(m, num_samples = 20, steps = 100, horizon = 4, inner_samples = 8, seed = 1)
  expect_true(all(is.finite(bayescomp_extract(d, "x"))))
})

test_that("generated quantities refuse a fit with different parameters", {
  m <- bayescomp_model("parameters { real a; real b; } model { a ~ normal(0, 1); b ~ normal(0, 1); }")
  f <- bayescomp_sample(m, chains = 1, iter_warmup = 20, iter_sampling = 20, seed = 1, refresh = 0)
  swapped <- bayescomp_model("parameters { real b; real a; } model { } generated quantities { real s = a; }")
  expect_error(bayescomp_generate_quantities(swapped, f), "parameter declarations")
})

test_that("Gibbs slice kernel leaves a bounded uniform invariant", {
  m <- bayescomp_model("parameters { real<lower=0, upper=1> x; } model { }")
  f <- bayescomp_gibbs(m, blocks = list(list(pars = "x", kernel = "slice")), chains = 2,
                       iter_warmup = 100, iter_sampling = 1500, seed = 1)
  x <- bayescomp_extract(f, "x")
  expect_equal(mean(x), 0.5, tolerance = 0.06)
  expect_equal(stats::var(x), 1 / 12, tolerance = 0.15)
})

test_that("gradient free samplers initialize without gradients", {
  # svd_U of a parameter dependent matrix has a value but no gradient
  m <- bayescomp_model("parameters { real x; } model { matrix[2, 2] A = [[x, 1], [0, 2]];
                        target += -0.5 * square(x) + 0 * svd_U(A)[1, 1]; }")
  expect_error(bayescomp_sample(m, chains = 1, iter_warmup = 5, iter_sampling = 5, seed = 1, refresh = 0),
               class = "bayescomp_unsupported_error")
  f <- bayescomp_sample(m, algorithm = "slice", chains = 1, iter_warmup = 20, iter_sampling = 20, seed = 1, refresh = 0)
  expect_equal(dim(f$draws)[1:2], c(20L, 1L))
})

test_that("state space parameter names line up with the values", {
  m <- bayescomp_model("parameters { matrix[2, 2] M; real s; } model { to_vector(M) ~ normal(0, 1); s ~ normal(0, 1); }")
  u <- c(1, 2, 3, 4, 5)
  v <- ssm_theta(m, u)
  expect_equal(unname(v[c("M[2,1]", "M[1,2]", "s")]), c(2, 3, 5))
})

test_that("the compatibility matrix is well formed and its test files exist", {
  cm <- bayescomp_compatibility()
  expect_true(all(cm$status %in% c("implemented", "partial", "not_implemented", "not_applicable")))
  files <- unique(unlist(strsplit(cm$tests[nzchar(cm$tests)], ";")))
  expect_true(all(file.exists(test_path(files))))
})
