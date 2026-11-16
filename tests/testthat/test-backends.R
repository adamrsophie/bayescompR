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
  # capabilities beyond the density are never taken from a native backend
  bayescomp_register_backend("fake5", c("log_density", "gradient", "constrain", "generated_quantities", "rng"),
    ir_versions = "1.0.0", compile = function(ir, model) list(log_density_gradient = function(u, jacobian, propto)
      list(value = -(u - 1)^2 / 8, grad = -(u - 1) / 4)))
  on.exit(rm("fake5", envir = .backends), add = TRUE)
  expect_error(bayescomp_model(code, backend = "fake5"), "allow_fallback")
  m5 <- bayescomp_model(code, backend = "fake5", allow_fallback = TRUE)
  expect_equal(m5$provenance$constrain, "r")
  expect_equal(m5$provenance$rng, "r")
  # native inputs: IR and data in Stan order
  md <- bayescomp_model("data { array[2] vector[3] x; } parameters { real mu; } model { mu ~ normal(0, 1); }",
                        data = list(x = matrix(1:6, 2, 3)))
  inp <- bayescomp_native_inputs(md)
  expect_equal(inp$dim, 1L)
  expect_match(inp$data, "\"values\":\\[1,3,5,2,4,6\\]")
  expect_error(bayescomp_native_error("E_REJECT: no"), class = "bayescomp_reject_error")
  expect_error(bayescomp_native_error("E_DOMAIN: no"), class = "bayescomp_domain_error")
  expect_error(bayescomp_native_error("E_SHAPE: no"), class = "bayescomp_shape_error")
  expect_error(bayescomp_native_error("E_INTERNAL: boom"), "boom")
})

test_that("value only samplers use a backend that provides only the gradient function", {
  bayescomp_register_backend("fake_ldg", c("log_density", "gradient"), ir_versions = "1.0.0",
    compile = function(ir, model) list(log_density_gradient = function(u, jacobian, propto)
      list(value = -0.5 * sum((u - 1)^2), grad = -(u - 1))))
  on.exit(rm("fake_ldg", envir = .backends))
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(1, 1); }", backend = "fake_ldg", allow_fallback = TRUE)
  expect_equal(bayescomp_log_prob(m, 0.5), -0.125)
  fit <- bayescomp_sample(m, algorithm = "metropolis", chains = 1, iter_warmup = 50, iter_sampling = 50, seed = 1, refresh = 0)
  lp <- as.data.frame(fit$draws)[[1L]]
  expect_true(all(is.finite(lp)))
})
