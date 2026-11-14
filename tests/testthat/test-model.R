code8 <- "
data { int<lower=0> J; array[J] real y; array[J] real<lower=0> sigma; }
parameters { real mu; real<lower=0> tau; vector[J] eta; }
transformed parameters { vector[J] theta = mu + tau * eta; }
model { mu ~ normal(0, 5); tau ~ cauchy(0, 5); eta ~ std_normal(); y ~ normal(theta, sigma); }
generated quantities { array[J] real y_rep = normal_rng(theta, sigma); }"
d8 <- list(J = 8, y = c(28, 8, -3, 7, -1, 1, 18, 12), sigma = c(15, 10, 16, 11, 9, 11, 10, 18))

test_that("log density, gradient and transforms agree with direct computation", {
  m <- bayescomp_model(code8, data = d8)
  expect_equal(m$dim, 10L)
  set.seed(3)
  u <- stats::rnorm(10)
  mu <- u[1]; tau <- exp(u[2]); eta <- u[3:10]
  ref <- stats::dnorm(mu, 0, 5, log = TRUE) + stats::dcauchy(tau, 0, 5, log = TRUE) + sum(stats::dnorm(eta, log = TRUE)) +
    sum(stats::dnorm(d8$y, mu + tau * eta, d8$sigma, log = TRUE))
  expect_equal(bayescomp_log_prob(m, u, jacobian = FALSE), ref)
  expect_equal(bayescomp_log_prob(m, u), ref + u[2])
  g <- bayescomp_grad_log_prob(m, u)
  expect_equal(as.vector(g), num_grad(function(x) bayescomp_log_prob(m, x), u), tolerance = 1e-6)
  p <- bayescomp_constrain(m, u, seed = 1)
  expect_equal(p$theta, mu + tau * eta)
  expect_length(p$y_rep, 8)
  expect_equal(bayescomp_unconstrain(m, p[c("mu", "tau", "eta")]), u)
  bad <- p[c("mu", "tau", "eta")]; bad$tau <- -1
  expect_error(bayescomp_unconstrain(m, bad), "tau", class = "bayescomp_data_error")
  expect_equal(bayescomp_param_names(m)[1:4], c("mu", "tau", "eta[1]", "eta[2]"))
})

test_that("propto drops constants in sampling statements only", {
  m <- bayescomp_model("data { real y; } parameters { real mu; } model { y ~ normal(mu, 1); target += normal_lpdf(mu | 0, 1); }", data = list(y = 0.5))
  full <- bayescomp_log_prob(m, 0.2)
  prop <- bayescomp_log_prob(m, 0.2, propto = TRUE)
  expect_equal(full - prop, -0.5 * log(2 * pi))
})

test_that("data are validated", {
  expect_error(bayescomp_model(code8, data = d8[-1]), class = "bayescomp_data_error")
  bad <- d8; bad$sigma[2] <- -1
  expect_error(bayescomp_model(code8, data = bad), "sigma", class = "bayescomp_data_error")
  bad <- d8; bad$J <- 7
  expect_error(bayescomp_model(code8, data = bad), class = "bayescomp_data_error")
  expect_error(bayescomp_model("data { int n; } model { }", data = list(n = 1.5)), class = "bayescomp_data_error")
})

test_that("arrays, matrices and indexing follow Stan semantics", {
  m <- bayescomp_model("
    data { array[2, 3] real a; matrix[2, 3] M; array[2] vector[3] v; }
    transformed data {
      array[3] real r = a[2];
      row_vector[3] mr = M[2];
      vector[2] mc = M[:, 3];
      real s = v[2, 3];
      array[2] int idx = {2, 1};
      matrix[2, 3] M2 = M[idx];
      matrix[2, 2] B = M[1:2, 2:3];
      vector[0] e = M[1, 3:2]';
    }
    generated quantities { array[3] real rr = r; row_vector[3] mm = mr; vector[2] cc = mc; real ss = s; matrix[2,3] m2 = M2; matrix[2,2] b = B; }",
    data = list(a = matrix(1:6, 2), M = matrix(1:6, 2), v = matrix(1:6, 2)))
  p <- bayescomp_constrain(m, numeric(0))
  expect_equal(p$rr, c(2, 4, 6))
  expect_equal(p$mm, c(2, 4, 6))
  expect_equal(p$cc, c(5, 6))
  expect_equal(p$ss, 6)
  expect_equal(p$m2, matrix(1:6, 2)[2:1, ])
  expect_equal(p$b, matrix(1:6, 2)[, 2:3])
})

test_that("runtime errors carry classes", {
  m <- bayescomp_model("data { vector[3] x; } generated quantities { real y = x[4]; }", data = list(x = 1:3))
  expect_error(bayescomp_constrain(m, numeric(0)), class = "bayescomp_domain_error")
  m <- bayescomp_model("parameters { real x; } model { reject(\"no \", x); }")
  expect_error(bayescomp_log_prob(m, 1), "no 1", class = "bayescomp_reject_error")
})

test_that("user functions, recursion and closures over data work", {
  m <- bayescomp_model("
    functions {
      int fib(int n) { if (n < 2) return n; return fib(n - 1) + fib(n - 2); }
      real my_lpdf(real y, real mu) { return normal_lpdf(y | mu, 1); }
    }
    parameters { real x; }
    model { x ~ my(1); }
    generated quantities { int f10 = fib(10); }")
  expect_equal(bayescomp_constrain(m, 0)$f10, 55L)
  expect_equal(bayescomp_log_prob(m, 0.3), stats::dnorm(0.3, 1, 1, log = TRUE))
})

test_that("the IR serializes to JSON", {
  j <- bayescomp_ir("parameters { real y; } model { y ~ std_normal(); }", json = TRUE)
  expect_match(j, "\"ir_version\":\"1.0.0\"")
  expect_match(j, "std_normal")
})
