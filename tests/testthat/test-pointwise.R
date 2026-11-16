test_that("pointwise log likelihood matches a generated quantity", {
  skip_on_cran()
  set.seed(1); N <- 8; x <- rnorm(N); y <- 1 + x + rnorm(N); k <- rpois(N, 3)
  code <- "data { int N; vector[N] x; vector[N] y; array[N] int k; } parameters { real a; real<lower=0> s; real<lower=0> lam; }
    model { a ~ normal(0, 5); s ~ exponential(1); for (n in 1:N) y[n] ~ normal(a + x[n], s); target += poisson_lpmf(k | lam); }
    generated quantities { vector[N] log_lik; for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | a + x[n], s) + poisson_lpmf(k[n] | lam); }"
  m <- bayescomp_model(code, data = list(N = N, x = x, y = y, k = k))
  fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 50, seed = 1, refresh = 0)
  ll <- bayescomp_log_lik(fit)
  expect_equal(dim(ll), c(100L, 2L * N))
  expect_equal(colnames(ll)[c(1, N + 1)], c("y[1]", "k[1]"))
  expect_equal(unname(ll[, 1:N] + ll[, N + 1:N]), unname(draws_matrix_of(fit, "log_lik")))
  expect_equal(dim(bayescomp_log_lik(fit, "array")), c(50L, 2L, 2L * N))
  yr <- bayescomp_posterior_predict(fit, seed = 2)
  expect_equal(dim(yr), dim(ll))
  expect_true(all(yr[, N + 1:N] == round(yr[, N + 1:N])))
  expect_equal(yr, bayescomp_posterior_predict(fit, seed = 2))
})

test_that("models without data sampling statements give an error", {
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
  fit <- bayescomp_sample(m, chains = 1, iter_warmup = 20, iter_sampling = 5, seed = 1, refresh = 0)
  expect_error(bayescomp_log_lik(fit), "no sampling statement on data")
})

test_that("loo, log_lik and posterior_predict fall back to the model block", {
  skip_on_cran()
  skip_if_not_installed("loo")
  m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; } model { y ~ normal(mu, 1); }",
                       data = list(N = 4, y = c(-1, 0, 1, 2)))
  fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100, seed = 1, refresh = 0)
  lo <- suppressWarnings(loo::loo(fit))
  expect_equal(nrow(lo$pointwise), 4L)
  expect_equal(log_lik.bayescomp_fit(fit), bayescomp_log_lik(fit))
  expect_equal(dim(posterior_predict.bayescomp_fit(fit)), c(200L, 4L))
  expect_error(loo.bayescomp_fit(fit, variable = "ll"), "unknown variable")
})

test_that("fused GLMs give one pointwise term per row as the unfused model does", {
  set.seed(3); N <- 6; X <- matrix(rnorm(2 * N), N); z <- rnorm(N); k <- rpois(N, 2)
  hd <- "data { int N; matrix[N, 2] X; vector[N] z; array[N] int k; } parameters { real a; vector[2] b; real<lower=0> s; }
    model { a ~ normal(0, 1); b ~ normal(0, 1); s ~ lognormal(0, 0.5); "
  draws <- function(body) {
    m <- bayescomp_model(paste(hd, body, "}"), data = list(N = N, X = X, z = z, k = k))
    fit <- bayescomp_sample(m, chains = 1, iter_warmup = 30, iter_sampling = 4, seed = 1, refresh = 0)
    list(ll = bayescomp_log_lik(fit), pp = bayescomp_posterior_predict(fit, seed = 2))
  }
  for (p in list(c("z ~ normal_id_glm(X, a, b, s);", "z ~ normal(a + X * b, s);"),
                 c("k ~ poisson_log_glm(X, a, b);", "k ~ poisson_log(a + X * b);"),
                 c("k ~ neg_binomial_2_log_glm(X, a, b, s);", "k ~ neg_binomial_2_log(a + X * b, s);"))) {
    g <- draws(p[1]); u <- draws(p[2])
    expect_equal(dim(g$ll), c(4L, N))
    expect_equal(g$ll, u$ll, tolerance = 1e-8)
    expect_equal(g$pp, u$pp, tolerance = 1e-8)
  }
})

test_that("pointwise terms share non elementwise arguments and multivariate predictions keep their shape", {
  set.seed(4); N <- 5; y <- c(1L, 3L, 2L, 3L, 1L)
  m <- bayescomp_model("data { int N; array[N] int y; } parameters { real a; } model { a ~ normal(0, 1); y ~ categorical_logit([a, 0, -a]'); }",
                       data = list(N = N, y = y))
  fit <- bayescomp_sample(m, chains = 1, iter_warmup = 30, iter_sampling = 4, seed = 1, refresh = 0)
  a <- as.numeric(draws_matrix_of(fit, "a"))
  ref <- t(vapply(a, function(ai) { e <- c(ai, 0, -ai); (e - log(sum(exp(e))))[y] }, numeric(N)))
  expect_equal(unname(bayescomp_log_lik(fit)), ref)
  expect_equal(dim(bayescomp_posterior_predict(fit, seed = 1)), c(4L, N))
  Y <- matrix(rnorm(2 * N), N)
  m2 <- bayescomp_model("data { int N; array[N] vector[2] Y; } parameters { vector[2] mu; } model { mu ~ normal(0, 1); Y ~ multi_normal(mu, diag_matrix(rep_vector(1, 2))); }",
                        data = list(N = N, Y = Y))
  fit2 <- bayescomp_sample(m2, chains = 1, iter_warmup = 30, iter_sampling = 4, seed = 1, refresh = 0)
  expect_equal(dim(bayescomp_log_lik(fit2)), c(4L, 1L))
  pp <- bayescomp_posterior_predict(fit2, seed = 1)
  expect_equal(dim(pp), c(4L, 2L * N))
  expect_equal(colnames(pp)[c(1, N + 1)], c("Y[1,1]", "Y[1,2]"))
})

test_that("one cutpoint vector per observation gives one term per observation", {
  set.seed(5); N <- 4; y <- c(1L, 3L, 2L, 3L); C <- t(replicate(N, sort(rnorm(2))))
  m <- bayescomp_model("data { int N; array[N] int y; array[N] vector[2] c; } parameters { real a; }
    model { a ~ normal(0, 1); y ~ ordered_logistic(rep_vector(a, N), c); }", data = list(N = N, y = y, c = C))
  fit <- bayescomp_sample(m, chains = 1, iter_warmup = 30, iter_sampling = 3, seed = 1, refresh = 0)
  a <- as.numeric(draws_matrix_of(fit, "a"))
  ref <- t(vapply(a, function(ai) vapply(1:N, function(i) log(-diff(c(1, plogis(ai - C[i, ]), 0))[y[i]]), 0), numeric(N)))
  expect_equal(unname(bayescomp_log_lik(fit)), ref)
  expect_true(all(bayescomp_posterior_predict(fit, seed = 1) %in% 1:3))
})
