interrupt_at <- function(text) function(msg) {
  if (grepl(text, msg, fixed = TRUE)) signalCondition(structure(class = c("interrupt", "condition"), list(message = "", call = NULL)))
}

test_that("an interrupt keeps the finished chains or the partial chain", {
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(1, 2); }")
  local_mocked_bindings(chain_message = interrupt_at("Chain 2: Iteration: 150"))
  expect_warning(fit <- bayescomp_sample(m, chains = 3, iter_warmup = 100, iter_sampling = 100, seed = 1, refresh = 50, cores = 1),
                 "the 1 chains that finished")
  expect_true(fit$interrupted)
  expect_equal(dim(fit$draws)[1:2], c(100L, 1L))
  local_mocked_bindings(chain_message = interrupt_at("Chain 1: Iteration: 150"))
  expect_warning(fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100, seed = 1, refresh = 50, cores = 1),
                 "the 50 draws of chain 1")
  expect_equal(dim(fit$draws)[1:2], c(50L, 1L))
  local_mocked_bindings(chain_message = interrupt_at("Chain 1: Iteration:  50"))
  expect_error(suppressWarnings(bayescomp_sample(m, chains = 1, iter_warmup = 100, iter_sampling = 100, seed = 1, refresh = 50)),
               "before any draw")
})

test_that("resume continues each chain without warmup", {
  m <- bayescomp_model("parameters { vector[2] mu; } model { mu ~ normal(1, 2); }")
  fit <- bayescomp_sample(m, chains = 2, iter_warmup = 200, iter_sampling = 50, seed = 1, refresh = 0, metric = "dense_e")
  f2 <- bayescomp_sample(m, resume = fit, iter_sampling = 300, seed = 2, refresh = 0)
  expect_equal(dim(f2$draws)[1:2], c(300L, 2L))
  expect_equal(f2$step_size, fit$step_size)
  expect_equal(f2$inv_metric, fit$inv_metric)
  expect_equal(f2$config$metric, "dense_e")
  expect_lt(abs(mean(f2$draws[, , "mu[1]"]) - 1), 0.6)
  expect_error(bayescomp_sample(m, resume = list()), "fit from bayescomp_sample")
})

test_that("divergence locations are stored with the chain and iteration", {
  skip_on_cran()
  m <- bayescomp_model("data { int J; array[J] real y; array[J] real<lower=0> sigma; } parameters { real mu; real<lower=0> tau; array[J] real theta; }
    model { mu ~ normal(0, 5); tau ~ cauchy(0, 5); theta ~ normal(mu, tau); y ~ normal(theta, sigma); }",
    data = list(J = 8, y = c(28, 8, -3, 7, -1, 1, 18, 12), sigma = c(15, 10, 16, 11, 9, 11, 10, 18)))
  fit <- bayescomp_sample(m, chains = 2, iter_warmup = 300, iter_sampling = 300, seed = 3, refresh = 0)
  nd <- sum(fit$sampler_diagnostics[, , "divergent__"])
  skip_if(nd == 0)
  d <- fit$divergences
  expect_equal(nrow(d), nd)
  expect_equal(names(d), c("chain", "iteration", "mu", "tau", sprintf("theta[%d]", 1:8)))
  for (k in seq_len(nrow(d))) expect_equal(unname(fit$sampler_diagnostics[d$iteration[k], d$chain[k], "divergent__"]), 1)
  expect_true(all(d$tau > 0))
})
