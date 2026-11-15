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
