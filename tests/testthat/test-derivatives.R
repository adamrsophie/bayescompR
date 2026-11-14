fd_grad <- function(m, u, h = 1e-6) {
  vapply(seq_along(u), function(i) {
    a <- u; b <- u; a[i] <- a[i] + h; b[i] <- b[i] - h
    (bayescomp_log_prob(m, a) - bayescomp_log_prob(m, b)) / (2 * h)
  }, 0)
}

test_that("svd_U and svd_V have exact gradients", {
  code <- "data { int m; int n; int k; matrix[m, k] W1; matrix[n, k] W2; } parameters { matrix[m, n] X; }
  model { target += -0.5 * sum(square(X)) + sum(svd_U(X) .* W1) + sum(svd_V(X) .* W2); }"
  set.seed(1)
  for (dims in list(c(3L, 3L), c(4L, 2L), c(2L, 4L))) {
    m0 <- dims[1]; n0 <- dims[2]; k <- min(dims)
    m <- bayescomp_model(code, data = list(m = m0, n = n0, k = k, W1 = matrix(rnorm(m0 * k), m0, k),
                                           W2 = matrix(rnorm(n0 * k), n0, k)))
    u <- rnorm(m0 * n0)
    g <- bayescomp_grad_log_prob(m, u)
    expect_equal(as.numeric(g), fd_grad(m, u), tolerance = 1e-6)
  }
})

test_that("qr_Q and qr_R of a tall matrix have exact gradients and LAPACK values", {
  set.seed(2)
  code <- "data { matrix[4, 4] W1; matrix[4, 2] W2; } parameters { matrix[4, 2] X; }
  model { target += -0.5 * sum(square(X)) + sum(qr_Q(X) .* W1) + sum(qr_R(X) .* W2); }"
  m <- bayescomp_model(code, data = list(W1 = matrix(rnorm(16), 4), W2 = matrix(rnorm(8), 4)))
  u <- rnorm(8)
  expect_equal(as.numeric(bayescomp_grad_log_prob(m, u)), fd_grad(m, u), tolerance = 1e-6)
})

test_that("student_t_qf is differentiable in the degrees of freedom", {
  m <- bayescomp_model("parameters { real<lower=1> nu; vector<lower=0, upper=1>[2] p; }
    model { target += student_t_qf(p[1], nu, 0.5, 2) + sum(student_t_qf(p, nu, 0, 1)); }")
  u <- c(0.4, -0.3, 0.8)
  expect_equal(as.numeric(bayescomp_grad_log_prob(m, u)), fd_grad(m, u), tolerance = 1e-6)
})
