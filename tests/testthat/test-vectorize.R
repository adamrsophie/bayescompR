same_as_loop <- function(code, data, u, loops = NULL) {
  old <- options(bayescompR.vectorize = FALSE)
  on.exit(options(old))
  m0 <- bayescomp_model(code, data = data)
  options(bayescompR.vectorize = TRUE)
  m1 <- bayescomp_model(code, data = data)
  if (!is.null(loops)) expect_equal(attr(m1$ir, "vectorized_loops") %||% 0L, loops)
  for (propto in c(FALSE, TRUE))
    expect_equal(bayescomp_grad_log_prob(m1, u, propto = propto), bayescomp_grad_log_prob(m0, u, propto = propto), tolerance = 1e-10)
  invisible(m1)
}

test_that("independent loops are vectorized with the same density and gradient", {
  set.seed(1)
  N <- 40; X <- matrix(rnorm(N * 3), N); g <- sample(5, N, TRUE)
  d <- list(N = N, X = X, g = g, y = rnorm(N), x = rnorm(N), k = rpois(N, 2), b01 = rbinom(N, 1, 0.5))
  dd <- "data { int N; matrix[N, 3] X; array[N] int g; vector[N] y; vector[N] x; array[N] int k; array[N] int b01; }"
  same_as_loop(paste(dd, "parameters { real a; vector[3] b; vector[5] u; real<lower=0> s; }
    model { vector[N] mu; for (n in 1:N) mu[n] = a + X[n] * b + u[g[n]]; for (n in 1:N) y[n] ~ normal(mu[n], s); }"),
    d, seq(-0.4, 0.4, length.out = 10), loops = 2L)
  same_as_loop(paste(dd, "parameters { real a; real b; } model { for (n in 1:N) b01[n] ~ bernoulli_logit(a + b * x[n]); }"), d, c(0.1, -0.2), 1L)
  same_as_loop(paste(dd, "parameters { real a; real<lower=0> s; } model { for (n in 1:N) target += normal_lpdf(y[n] | a + X[n, 2] * exp(a), s); }"),
               d, c(0.1, 0.2), 1L)
  same_as_loop(paste(dd, "parameters { real a; } model { for (n in 1:N) k[n] ~ poisson_log(a + log1p(square(x[n])) - x[n] / 2); }"), d, 0.3, 1L)
  same_as_loop(paste(dd, "parameters { real a; } model { vector[N] mu = rep_vector(a, N); for (n in 1:N) mu[n] += 2 * x[n]; y ~ normal(mu, 1); }"), d, 0.3, 1L)
  same_as_loop(paste(dd, "parameters { real a; } transformed parameters { vector[N] m; for (n in 1:N) m[n] = inv_logit(a + y[n]); } model { y ~ normal(m, 1); }"), d, 0.3, 1L)
  same_as_loop(paste(dd, "parameters { real a; real<lower=0> s; } model { for (j in 1:2) for (n in 1:N) y[n] ~ normal(a * j - x[n] / s, s); }"), d, c(0.1, 0.2), 1L)
  same_as_loop(paste(dd, "parameters { real a; } model { for (n in 1:N) { real m = a * x[n]; y[n] ~ normal(m, 1); } }"), d, 0.2, 0L)
})

test_that("loops that are not independent stay loops", {
  d <- list(N = 4L, y = c(1, 2, 3, 4), k = c(1L, 2L, 3L, 4L))
  dd <- "data { int N; vector[N] y; array[N] int k; }"
  same_as_loop(paste(dd, "parameters { real a; } model { vector[N] mu; mu[1] = a; for (n in 2:N) mu[n] = mu[n - 1] * 0.5 + a; y ~ normal(mu, 1); }"), d, 0.2, 0L)
  same_as_loop(paste(dd, "parameters { real<lower=0> s; } model { for (n in 1:N) y[n] ~ normal(0, s) T[-5, ]; }"), d, 0.2, 0L)
  same_as_loop(paste(dd, "parameters { real a; } model { for (n in 1:N) k[n] ~ poisson(exp(a + k[n] / 2)); }"), d, 0.2, 0L) |> suppressWarnings()
  same_as_loop(paste(dd, "parameters { real a; } model { for (n in 1:N) y[n] ~ normal(a * n, 1); }"), d, 0.2, 0L)
})

test_that("empty and out of range loops behave as loops", {
  same_as_loop("data { int N; vector[N] y; } parameters { real a; } model { for (n in 1:N) y[n] ~ normal(a, 1); }", list(N = 0L, y = numeric(0)), 0.1, 1L)
  same_as_loop("data { int N; vector[N] y; } parameters { real a; } model { for (n in 3:2) y[n] ~ normal(a, 1); }", list(N = 3L, y = 1:3), 0.1, 1L)
  m <- bayescomp_model("data { int N; vector[N] y; } parameters { real a; } model { for (n in 1:(N + 1)) y[n] ~ normal(a, 1); }", data = list(N = 3L, y = 1:3))
  expect_error(bayescomp_log_prob(m, 0.1), "index out of range")
})

test_that("element writes in a loop that cannot be vectorized use linear memory", {
  code <- "data { int N; vector[N] y; } parameters { real<lower=-1, upper=1> r; real a; }
    model { vector[N] mu; mu[1] = a; for (n in 2:N) mu[n] = r * mu[n - 1] + a; y ~ normal(mu, 1); }"
  grad_mb <- function(N) {
    m <- bayescomp_model(code, data = list(N = N, y = rep(0.1, N)))
    gc(reset = TRUE)
    bayescomp_grad_log_prob(m, c(0.2, 0.1))
    sum(gc()[, "max used"] * c(56, 8)) / 2^20
  }
  skip_on_cran()
  small <- grad_mb(500L)
  big <- grad_mb(8000L)
  # quadratic memory needed over 1000 MB more
  expect_lt(big - small, 300)
})

test_that("a data covariance is factored once and gives the same density", {
  set.seed(3)
  K <- 4; S <- crossprod(matrix(rnorm(K * K), K)) + diag(K); Y <- matrix(rnorm(10 * K), 10)
  m <- bayescomp_model("data { int K; matrix[10, K] Y; cov_matrix[K] S; } parameters { vector[K] mu; } model { for (n in 1:10) Y[n]' ~ multi_normal(mu, S); }",
                       data = list(K = K, Y = Y, S = S))
  m2 <- bayescomp_model("data { int K; array[10] vector[K] Y; cov_matrix[K] S; } parameters { vector[K] mu; } model { Y ~ multi_normal(mu, S); }",
                        data = list(K = K, Y = Y, S = S))
  u <- rep(0.2, K)
  expect_equal(bayescomp_grad_log_prob(m, u, propto = FALSE), bayescomp_grad_log_prob(m2, u, propto = FALSE))
  R <- chol(S)
  ref <- sum(apply(Y, 1, function(y) -0.5 * K * log(2 * pi) - sum(log(diag(R))) - 0.5 * sum(backsolve(R, y - u, transpose = TRUE)^2)))
  expect_equal(bayescomp_log_prob(m2, u, propto = FALSE), ref)
})

test_that("fused correlation and bound transforms have exact gradients", {
  m <- bayescomp_model("data { vector[3] lo; } parameters { corr_matrix[4] R; cholesky_factor_corr[3] Lc; real<lower=-1, upper=2> a;
      vector<lower=0, upper=5>[4] b; vector<lower=lo, upper=lo + 3>[3] c; real<lower=a, upper=a + 1> d; }
    model { R ~ lkj_corr(2); Lc ~ lkj_corr_cholesky(1.5); target += R[1, 2] * 0.3 + sum(Lc); a ~ normal(0, 1); b ~ normal(1, 2); d ~ normal(0, 1); }",
    data = list(lo = c(-1, 0, 1)))
  set.seed(2); u <- rnorm(m$dim) * 2
  for (jac in c(TRUE, FALSE)) {
    f <- function(v) bayescomp_log_prob(m, v, jacobian = jac)
    num <- vapply(seq_along(u), function(i) { e <- replace(numeric(length(u)), i, 1e-6); (f(u + e) - f(u - e)) / 2e-6 }, 0)
    expect_equal(as.numeric(bayescomp_grad_log_prob(m, u, jacobian = jac)), num, tolerance = 1e-6)
  }
  R <- bayescomp_constrain(m, u)$R
  expect_equal(diag(R), rep(1, 4))
  expect_true(all(eigen(R)$values > 0))
})
