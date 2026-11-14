test_that("native backends are used, validated and never mixed silently", {
  code <- "parameters { real mu; } model { mu ~ normal(1, 2); } generated quantities { real z = normal_rng(mu, 1); }"
  calls <- new.env()
  calls$n <- 0L
  bayescomp_register_backend("fake", c("log_density", "gradient"), ir_versions = "1.0.0",
    compile = function(ir, model) {
      expect_equal(ir$ir_version, "1.0.0")
      list(log_density_gradient = function(u, jacobian, propto) {
        calls$n <- calls$n + 1L
        list(value = if (propto) -(u - 1)^2 / 8 else stats::dnorm(u, 1, 2, log = TRUE), grad = -(u - 1) / 4)
      })
    }, description = "test backend")
  on.exit(rm("fake", envir = .backends))
  expect_true("fake" %in% bayescomp_backends()$name)
  expect_error(bayescomp_model(code, backend = "fake"), "allow_fallback")
  m <- bayescomp_model(code, backend = "fake", allow_fallback = TRUE)
  expect_equal(m$provenance$gradient, "fake")
  expect_equal(m$provenance$generated_quantities, "r")
  expect_equal(bayescomp_log_prob(m, 0.5), stats::dnorm(0.5, 1, 2, log = TRUE))
  expect_gt(calls$n, 0L)
  v <- bayescomp_validate_backend(m, n = 5)
  expect_lt(max(v$lp_abs, v$grad_abs), 1e-10)
  fit <- bayescomp_sample(m, chains = 1, iter_warmup = 100, iter_sampling = 100, seed = 1, refresh = 0)
  expect_equal(dim(fit$draws)[3L], 3L)
  expect_error(bayescomp_model(code, backend = "missing"), "not available")
})
