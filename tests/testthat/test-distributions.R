# Each row: Stan density, R reference, argument generator. Values are
# compared with the normalized density and gradients with finite
# differences of the AD implementation.

dist_cases <- list(
  normal = list(function(y, a, b) sum(stats::dnorm(y, a, b, log = TRUE)), c(0.3, -1.1), 0.5, 1.4),
  student_t = list(function(y, n, a, b) sum(stats::dt((y - a) / b, n, log = TRUE) - log(b)), c(0.3, 2), 4.5, 0.5, 1.4),
  cauchy = list(function(y, a, b) sum(stats::dcauchy(y, a, b, log = TRUE)), c(0.3, 2), 0.5, 1.4),
  double_exponential = list(function(y, a, b) sum(-log(2 * b) - abs(y - a) / b), c(0.3, 2), 0.5, 1.4),
  logistic = list(function(y, a, b) sum(stats::dlogis(y, a, b, log = TRUE)), c(0.3, 2), 0.5, 1.4),
  gumbel = list(function(y, a, b) sum(-log(b) - (y - a) / b - exp(-(y - a) / b)), c(0.3, 2), 0.5, 1.4),
  lognormal = list(function(y, a, b) sum(stats::dlnorm(y, a, b, log = TRUE)), c(0.3, 2), 0.5, 1.4),
  chi_square = list(function(y, n) sum(stats::dchisq(y, n, log = TRUE)), c(0.3, 2), 3.5),
  exponential = list(function(y, b) sum(stats::dexp(y, b, log = TRUE)), c(0.3, 2), 1.4),
  gamma = list(function(y, a, b) sum(stats::dgamma(y, a, b, log = TRUE)), c(0.3, 2), 2.5, 1.4),
  inv_gamma = list(function(y, a, b) sum(a * log(b) - lgamma(a) - (a + 1) * log(y) - b / y), c(0.3, 2), 2.5, 1.4),
  weibull = list(function(y, a, b) sum(stats::dweibull(y, a, b, log = TRUE)), c(0.3, 2), 2.5, 1.4),
  beta = list(function(y, a, b) sum(stats::dbeta(y, a, b, log = TRUE)), c(0.3, 0.8), 2.5, 1.4),
  uniform = list(function(y, a, b) sum(stats::dunif(y, a, b, log = TRUE)), c(0.3, 0.8), -1, 2),
  skew_normal = list(function(y, a, b, c) sum(log(2) + stats::dnorm(y, a, b, log = TRUE) + stats::pnorm(c * (y - a) / b, log.p = TRUE)), c(0.3, 2), 0.5, 1.4, 0.7),
  pareto = list(function(y, m, a) sum(log(a) + a * log(m) - (a + 1) * log(y)), c(1.3, 2), 0.5, 1.4),
  rayleigh = list(function(y, s) sum(log(y) - 2 * log(s) - y^2 / (2 * s^2)), c(0.3, 2), 1.4))

test_that("continuous densities match reference implementations", {
  for (nm in names(dist_cases)) {
    cs <- dist_cases[[nm]]
    ref <- cs[[1L]]
    args <- cs[-1L]
    f <- lib_fun(paste0(nm, "_lpdf"))
    expect_equal(val(do.call(f, args)), do.call(ref, args), tolerance = 1e-10, label = nm)
    # gradient with respect to all scalar arguments and the outcome
    lens <- lengths(args)
    x0 <- unlist(args)
    g <- function(x) {
      parts <- split(seq_along(x), rep(seq_along(lens), lens))
      a <- lapply(parts, function(ix) ad_gather(x, ix))
      do.call(f, unname(a))
    }
    expect_grad_ok(g, x0, tol = 1e-5, label = paste(nm, "gradient"))
  }
})

test_that("discrete mass functions match reference implementations", {
  k <- c(0L, 3L, 5L)
  expect_equal(lib_fun("poisson_lpmf")(k, 2.5), sum(stats::dpois(k, 2.5, log = TRUE)))
  expect_equal(lib_fun("poisson_log_lpmf")(k, log(2.5)), sum(stats::dpois(k, 2.5, log = TRUE)))
  expect_equal(lib_fun("binomial_lpmf")(k, 6L, 0.3), sum(stats::dbinom(k, 6, 0.3, log = TRUE)))
  expect_equal(lib_fun("binomial_logit_lpmf")(k, 6L, 0.4), sum(stats::dbinom(k, 6, plogis(0.4), log = TRUE)))
  expect_equal(lib_fun("bernoulli_lpmf")(c(0L, 1L, 1L), 0.3), sum(stats::dbinom(c(0, 1, 1), 1, 0.3, log = TRUE)))
  expect_equal(lib_fun("bernoulli_logit_lpmf")(c(0L, 1L, 1L), c(-1, 0.2, 2)), sum(stats::dbinom(c(0, 1, 1), 1, plogis(c(-1, 0.2, 2)), log = TRUE)))
  expect_equal(lib_fun("neg_binomial_2_lpmf")(k, 2.5, 3), sum(stats::dnbinom(k, size = 3, mu = 2.5, log = TRUE)))
  expect_equal(lib_fun("neg_binomial_2_log_lpmf")(k, log(2.5), 3), sum(stats::dnbinom(k, size = 3, mu = 2.5, log = TRUE)))
  expect_equal(lib_fun("neg_binomial_lpmf")(k, 2, 0.5), sum(stats::dnbinom(k, size = 2, prob = 0.5 / 1.5, log = TRUE)))
  expect_equal(lib_fun("beta_binomial_lpmf")(k, 6L, 2, 3),
               sum(lchoose(6, k) + lbeta(k + 2, 6 - k + 3) - lbeta(2, 3)))
  expect_equal(lib_fun("hypergeometric_lpmf")(2L, 5L, 4L, 6L), stats::dhyper(2, 4, 6, 5, log = TRUE))
  expect_equal(lib_fun("categorical_lpmf")(c(1L, 3L), c(0.2, 0.3, 0.5)), log(0.2) + log(0.5))
  expect_equal(lib_fun("multinomial_lpmf")(c(1L, 2L, 0L), c(0.2, 0.3, 0.5)), stats::dmultinom(c(1, 2, 0), prob = c(0.2, 0.3, 0.5), log = TRUE))
  th <- c(0.2, 0.7, 0.4)
  pb <- function(y) { p <- 1; for (t in th) p <- stats::convolve(p, rev(c(1 - t, t)), type = "open"); p[y + 1] }
  expect_equal(exp(lib_fun("poisson_binomial_lpmf")(1L, th)), pb(1), tolerance = 1e-12)
  ys <- function(y, a) log(a) + lbeta(y, a + 1)
  expect_equal(lib_fun("yule_simon_lpmf")(c(1L, 4L), 2.5), sum(ys(c(1, 4), 2.5)))
  expect_equal(sum(exp(vapply(1:2000, function(y) lib_fun("yule_simon_lpmf")(y, 2.5), 0))), 1, tolerance = 1e-4)
})

test_that("cumulative distribution functions match", {
  expect_equal(lib_fun("normal_lcdf")(0.3, 0.5, 1.4), stats::pnorm(0.3, 0.5, 1.4, log.p = TRUE))
  expect_equal(lib_fun("normal_lccdf")(0.3, 0.5, 1.4), stats::pnorm(0.3, 0.5, 1.4, lower.tail = FALSE, log.p = TRUE))
  expect_equal(lib_fun("gamma_lcdf")(0.3, 2, 1.5), stats::pgamma(0.3, 2, 1.5, log.p = TRUE))
  expect_equal(lib_fun("beta_lccdf")(0.3, 2, 1.5), stats::pbeta(0.3, 2, 1.5, lower.tail = FALSE, log.p = TRUE))
  expect_equal(lib_fun("poisson_lcdf")(3L, 2.5), stats::ppois(3, 2.5, log.p = TRUE))
  expect_equal(lib_fun("student_t_lcdf")(0.3, 3, 0.5, 1.4), stats::pt((0.3 - 0.5) / 1.4, 3, log.p = TRUE))
  expect_equal(lib_fun("exponential_lccdf")(0.3, 2), -0.6)
  expect_equal(lib_fun("weibull_lcdf")(0.3, 2, 1.4), stats::pweibull(0.3, 2, 1.4, log.p = TRUE))
  expect_grad_ok(function(x) lib_fun("gamma_lcdf")(ad_select(x, 1L), ad_select(x, 2L), ad_select(x, 3L)), c(0.7, 2, 1.5), tol = 1e-4)
  expect_grad_ok(function(x) lib_fun("normal_lccdf")(ad_select(x, 1L), ad_select(x, 2L), ad_select(x, 3L)), c(0.7, 0.2, 1.5))
})

test_that("multivariate densities match", {
  S <- matrix(c(2, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1.5), 3)
  y <- c(0.5, -0.2, 1)
  mu <- c(0, 0.1, 0.3)
  ref <- -1.5 * log(2 * pi) - 0.5 * log(det(S)) - 0.5 * sum((y - mu) * solve(S, y - mu))
  expect_equal(lib_fun("multi_normal_lpdf")(y, mu, S), ref)
  expect_equal(lib_fun("multi_normal_cholesky_lpdf")(y, mu, t(chol(S))), ref)
  expect_equal(lib_fun("multi_normal_prec_lpdf")(y, mu, solve(S)), ref)
  a <- c(2, 3, 1.5); th <- c(0.2, 0.5, 0.3)
  expect_equal(lib_fun("dirichlet_lpdf")(th, a), lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * log(th)))
  # LKJ normalizing constant: the density integrates to one for K = 2
  f <- function(r) exp(lib_fun("lkj_corr_lpdf")(matrix(c(1, r, r, 1), 2), 2.5))
  expect_equal(stats::integrate(Vectorize(f), -1, 1)$value, 1, tolerance = 1e-6)
  W <- matrix(c(2, 0.5, 0.5, 1), 2)
  ref_w <- function(W, nu, S) {
    p <- nrow(W)
    (nu - p - 1) / 2 * log(det(W)) - sum(diag(solve(S, W))) / 2 - nu * p / 2 * log(2) - nu / 2 * log(det(S)) -
      (p * (p - 1) / 4 * log(pi) + sum(lgamma((nu + 1 - seq_len(p)) / 2)))
  }
  expect_equal(lib_fun("wishart_lpdf")(W, 4, diag(2)), ref_w(W, 4, diag(2)))
})

test_that("rng functions draw from the right distributions", {
  set.seed(1)
  x <- vapply(1:4000, function(i) lib_fun("normal_rng")(1, 2), 0)
  expect_equal(mean(x), 1, tolerance = 0.1)
  expect_equal(stats::sd(x), 2, tolerance = 0.05)
  x <- vapply(1:4000, function(i) lib_fun("von_mises_rng")(0.5, 2), 0)
  expect_equal(atan2(mean(sin(x)), mean(cos(x))), 0.5, tolerance = 0.05)
  x <- vapply(1:4000, function(i) lib_fun("skew_double_exponential_rng")(0, 1, 0.3), 0)
  expect_equal(mean(x <= 0), 0.3, tolerance = 0.03)
})

test_that("edge cases of special functions and densities", {
  m <- bayescomp_model("generated quantities {
    real h1 = hypergeometric_2F1(1, 1, 2, -3);
    real h2 = hypergeometric_2F1(1, 1, 2, 0.4);
    real ie = inv_erfc(0.3);
  }")
  v <- bayescomp_constrain(m, numeric(0))
  # Stan 2.40 requires a positive gamma variate
  expect_error(lib_fun("gamma_lpdf")(0, 1, 2), "must be positive")
  expect_equal(v$h1, log(4) / 3, tolerance = 1e-10)
  expect_equal(v$h2, -log(0.6) / 0.4, tolerance = 1e-10)
  expect_equal(2 * stats::pnorm(-v$ie * sqrt(2)), 0.3)
  bad <- bayescomp_model("generated quantities { real h = hypergeometric_2F1(1, 1, 2, 1.5); }")
  expect_error(bayescomp_constrain(bad, numeric(0)), class = "bayescomp_domain_error")
  expect_error(bayescomp_model("parameters { real x; } model { target += lower_bound_jacobian(x, 0); }"),
               class = "bayescomp_semantic_error")
})
