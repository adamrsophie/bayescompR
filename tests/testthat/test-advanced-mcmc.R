target2 <- "parameters { vector[2] x; } model { x ~ multi_normal([1, -1]', [[1, 0.7], [0.7, 1]]); }"

moments_ok <- function(fit, tol = 0.25) {
  x <- bayescomp_extract(fit, "x")
  expect_equal(colMeans(x), c(1, -1), tolerance = tol)
  expect_equal(apply(x, 2, sd), c(1, 1), tolerance = tol)
  expect_equal(stats::cor(x)[1, 2], 0.7, tolerance = tol)
}

test_that("elliptical slice sampling targets the posterior", {
  m <- bayescomp_model(target2)
  fit <- bayescomp_sample(m, chains = 2, iter_warmup = 200, iter_sampling = 600,
                          algorithm = "elliptical_slice", metric = "dense_e", seed = 1, refresh = 0)
  moments_ok(fit, 0.15)
})

test_that("ensemble samplers target the posterior", {
  m <- bayescomp_model(target2)
  for (meth in c("de", "stretch")) {
    fit <- bayescomp_ensemble(m, method = meth, iter_warmup = 400, iter_sampling = 400, seed = 2)
    moments_ok(fit, 0.3)
    expect_true(all(fit$acceptance > 0.05))
  }
})

test_that("parallel tempering mixes between modes and estimates the evidence", {
  mix <- bayescomp_model("parameters { real x; }
    model { target += log_mix(0.3, normal_lpdf(x | -5, 0.5), normal_lpdf(x | 5, 0.5)); }")
  fit <- bayescomp_parallel_tempering(mix, num_temps = 8, chains = 1, iter_warmup = 300, iter_sampling = 600,
                                      reference = "standard", reference_scale = 8, seed = 1)
  expect_equal(mean(bayescomp_extract(fit, "x") < 0), 0.3, tolerance = 0.15)
  expect_equal(fit$log_evidence, 0, tolerance = 0.3)
  b <- fit$betas[[1L]]
  expect_equal(range(b), c(0, 1))
  expect_true(all(diff(b) > 0))
  two <- bayescomp_parallel_tempering(bayescomp_model(target2), num_temps = 2, chains = 1, iter_warmup = 50,
                                      iter_sampling = 100, seed = 1)
  expect_equal(fit_dims <- dim(two$draws)[1:2], c(100L, 1L))
  s <- pt_update_schedule(seq(0, 1, length.out = 5), c(0.8, 0.1, 0.1, 0.1))
  expect_lt(s[2], 0.25)
  # zero density regions: a standard normal truncated at 1 by reject()
  tr <- bayescomp_model("parameters { real x; } model { if (x > 1) reject(\"out\"); x ~ normal(0, 1); }")
  ft <- bayescomp_parallel_tempering(tr, num_temps = 8, chains = 1, seed = 2, reference = "standard")
  expect_lt(abs(ft$log_evidence - log(stats::pnorm(1))), 0.1)
  expect_lte(max(ft$draws[, , "x"]), 1)
})

test_that("Gibbs sampling detects conjugate conditionals", {
  y <- c(1, 2, 0.5, 1.5, 2.5)
  m <- bayescomp_model("data { int N; vector[N] y; }
    parameters { real mu; real<lower=0> tau; real<lower=0, upper=1> p; vector[2] b; real<lower=0> s; }
    model {
      mu ~ normal(0, 10); tau ~ gamma(2, 1); y ~ normal(mu, inv_sqrt(tau));
      p ~ beta(2, 3); 7 ~ binomial(10, p);
      b ~ normal(0, 1); [1, 2]' ~ normal(b, 1);
      s ~ student_t(3, 0, 1);
    }", data = list(N = 5, y = y))
  fit <- bayescomp_gibbs(m, chains = 2, iter_warmup = 100, iter_sampling = 300, seed = 1)
  expect_equal(fit$blocks$family[1:4], c("gaussian", "gamma1", "beta", "gaussian"))
  expect_equal(fit$blocks$kernel[5], "slice")
  expect_true(all(fit$blocks$acceptance[1:4] > 0.999))
  expect_equal(mean(bayescomp_extract(fit, "p")), 9 / 15, tolerance = 0.05)
  expect_equal(colMeans(bayescomp_extract(fit, "b")), c(0.5, 1), tolerance = 0.15)
  fit2 <- bayescomp_gibbs(m, blocks = list(list(pars = c("mu", "tau"), kernel = "metropolis"),
                                           list(pars = c("p", "b", "s"), kernel = "mala")),
                          chains = 1, iter_warmup = 300, iter_sampling = 300, seed = 2)
  expect_equal(mean(bayescomp_extract(fit2, "p")), 9 / 15, tolerance = 0.1)
  expect_error(bayescomp_gibbs(m, blocks = list("mu")), "every parameter")
  one <- bayescomp_model("parameters { real mu; } model { mu ~ normal(1, 2); }")
  f1 <- bayescomp_gibbs(one, chains = 2, iter_warmup = 50, iter_sampling = 200, seed = 1)
  expect_equal(nrow(f1$blocks), 1L)
  expect_equal(mean(bayescomp_extract(f1, "mu")), 1, tolerance = 0.4)
})

test_that("Riemannian samplers and MCLMC target the posterior", {
  m <- bayescomp_model(target2)
  fit <- bayescomp_riemannian(m, "smmala", chains = 2, iter_warmup = 150, iter_sampling = 400, seed = 1)
  moments_ok(fit, 0.25)
  fit <- bayescomp_mclmc(m, chains = 2, iter_warmup = 500, iter_sampling = 1500, seed = 1)
  moments_ok(fit, 0.2)
  # the preconditioned step size must not depend on the scale of the target
  for (sc in c(100, 0.01)) {
    ms <- bayescomp_model("data { real s; } parameters { vector[5] x; } model { x ~ normal(0, s); }", data = list(s = sc))
    fs <- bayescomp_mclmc(ms, chains = 1, iter_warmup = 500, iter_sampling = 1000, seed = 2)
    expect_equal(mean(apply(fs$unconstrained[[1L]], 2, stats::sd)) / sc, 1, tolerance = 0.2)
  }
  skip_if_not(long_tests())
  fit <- bayescomp_riemannian(m, "rmhmc", chains = 1, iter_warmup = 100, iter_sampling = 300, seed = 1)
  moments_ok(fit, 0.35)
  expect_gt(mean(fit$sampler_diagnostics[, , "accept_stat__"]), 0.6)
})

test_that("piecewise deterministic samplers target the posterior", {
  m <- bayescomp_model(target2)
  for (meth in c("zigzag", "bps")) {
    fit <- bayescomp_pdmp(m, meth, chains = 2, iter_warmup = 150, iter_sampling = 500, seed = 1)
    moments_ok(fit, 0.25)
    expect_true(all(fit$events > 0))
  }
})

test_that("normalizing flow VI fits a correlated Gaussian", {
  m <- bayescomp_model("parameters { vector[2] x; } model { x ~ multi_normal([1, -1]', [[1, 0.8], [0.8, 1]]); }")
  fit <- bayescomp_flow_vi(m, layers = 2, iter = 500, learning_rate = 0.05, seed = 2)
  x <- bayescomp_extract(fit, "x")
  expect_equal(colMeans(x), c(1, -1), tolerance = 0.15)
  expect_equal(cor(x)[1, 2], 0.8, tolerance = 0.1)
  expect_lt(fit$pareto_k, 0.7)
  expect_lt(abs(mean(utils::tail(fit$elbo, 50))), 0.2)
})

test_that("diagnostic plots draw for every type", {
  m <- bayescomp_model("parameters { real x; } model { x ~ normal(0, 1); }")
  fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 200, seed = 1, refresh = 0)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (ty in c("trace", "density", "rank", "ess_local", "ess_quantile")) expect_identical(plot(fit, type = ty), fit)
})
