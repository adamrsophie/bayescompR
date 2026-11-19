target_code <- "parameters { vector[2] x; real<lower=0> s; }
model { x ~ multi_normal([1, -1]', [[1, 0.7], [0.7, 1]]); s ~ gamma(3, 2); }"

check_moments <- function(fit) {
  s <- summary(fit)
  get <- function(v) s[s$variable == v, ]
  for (v in c("x[1]", "x[2]")) {
    r <- get(v)
    expect_lt(abs(r$mean - c("x[1]" = 1, "x[2]" = -1)[[v]]), 5 / sqrt(r$ess_bulk), label = v)
  }
  r <- get("s")
  expect_lt(abs(r$mean - 1.5), 5 * sqrt(0.75) / sqrt(r$ess_bulk))
  expect_true(all(s$rhat < 1.1))
}

test_that("MALA, Metropolis and slice sampling target the posterior", {
  skip_on_cran()
  m <- bayescomp_model(target_code)
  n <- if (long_tests()) 1000L else 300L
  for (a in c("mala", "barker", "metropolis", "slice")) {
    k <- if (a == "metropolis") 4L * n else n  # Metropolis steps are cheap and mix slowly
    fit <- bayescomp_sample(m, chains = 2, iter_warmup = k, iter_sampling = k, algorithm = a, seed = 7, refresh = 0)
    check_moments(fit)
  }
  fit <- bayescomp_sample(m, chains = 2, iter_warmup = 2L * n, iter_sampling = 2L * n, algorithm = "metropolis",
                          block = FALSE, seed = 7, refresh = 0)
  check_moments(fit)
})

test_that("value only density matches the gradient density up to a constant", {
  m <- bayescomp_model(target_code)
  f <- lp_value_fun(m)
  u1 <- c(0.1, 0.2, 0.3); u2 <- c(-1, 2, -0.5)
  d_full <- f(u1) - f(u2)
  d_prop <- bayescomp_log_prob(m, u1, propto = TRUE) - bayescomp_log_prob(m, u2, propto = TRUE)
  expect_equal(d_full, d_prop)
  expect_equal(f(u1), bayescomp_log_prob(m, u1))
})

test_that("Pathfinder fits initialize chains", {
  skip_on_cran()
  m <- bayescomp_model(target_code)
  pf <- suppressWarnings(bayescomp_pathfinder(m, num_paths = 2, num_draws = 50, seed = 1))
  fit <- bayescomp_sample(m, chains = 2, iter_warmup = 50, iter_sampling = 20, init = pf, seed = 1, refresh = 0)
  expect_equal(dim(fit$draws)[2L], 2L)
})

test_that("tempered SMC estimates the marginal likelihood", {
  skip_on_cran()
  y <- c(0.5, 1.2, 2, 0.8, 1.6)
  m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; }
    model { mu ~ normal(0, 2); y ~ normal(mu, 1); }", data = list(N = 5, y = y))
  S <- diag(5) + 4
  exact <- as.numeric(-0.5 * determinant(2 * pi * S)$modulus - 0.5 * sum(y * solve(S, y)))
  for (r in c("pathfinder", "standard")) {
    fit <- bayescomp_smc(m, num_particles = 400, reference = r, seed = 2)
    expect_equal(fit$log_evidence, exact, tolerance = 0.05, label = r)
    expect_equal(mean(bayescomp_extract(fit, "mu")), sum(y) / 5.25, tolerance = 0.1)
  }
  w <- systematic_resample(c(0.5, 0.25, 0.25))
  expect_length(w, 3L)
  expect_true(all(w %in% 1:3))
})

test_that("tempered SMC finds both modes of a mixture", {
  skip_on_cran()
  m <- bayescomp_model("parameters { real x; }
    model { target += log_mix(0.5, normal_lpdf(x | -4, 0.5), normal_lpdf(x | 4, 0.5)); }")
  fit <- bayescomp_smc(m, num_particles = 400, reference = "standard", reference_scale = 6, seed = 1)
  expect_equal(mean(bayescomp_extract(fit, "x") > 0), 0.5, tolerance = 0.15)
})

test_that("simulation based calibration of NUTS is uniform", {
  skip_if_not(long_tests())
  # Talts et al. (2018): ranks of prior draws among posterior draws
  m_code <- "data { int N; vector[N] y; } parameters { real mu; real<lower=0> sigma; }
    model { mu ~ normal(0, 1); sigma ~ lognormal(0, 0.5); y ~ normal(mu, sigma); }"
  set.seed(2024)
  reps <- 100L
  ranks <- matrix(NA_integer_, reps, 2L)
  for (r in seq_len(reps)) {
    mu <- stats::rnorm(1); sigma <- stats::rlnorm(1, 0, 0.5)
    y <- stats::rnorm(10, mu, sigma)
    m <- bayescomp_model(m_code, data = list(N = 10, y = y))
    fit <- bayescomp_sample(m, chains = 1, iter_warmup = 200, iter_sampling = 495, thin = 5, seed = r, refresh = 0)
    ranks[r, ] <- c(sum(bayescomp_extract(fit, "mu") < mu), sum(bayescomp_extract(fit, "sigma") < sigma))
  }
  for (j in 1:2) {
    counts <- tabulate(cut(ranks[, j], breaks = seq(-0.5, 99.5, length.out = 11), labels = FALSE), 10)
    expect_gt(stats::chisq.test(counts)$p.value, 0.001)
  }
})

test_that("unseeded runs advance the global stream and record their seed", {
  m <- bayescomp_model("parameters { real x; } model { x ~ normal(0, 1); }")
  run <- function(...) bayescomp_sample(m, chains = 1, iter_warmup = 30, iter_sampling = 10, refresh = 0, ...)
  set.seed(11); f1 <- run(); f2 <- run()
  set.seed(11); f3 <- run()
  expect_false(identical(f1$draws, f2$draws))
  expect_identical(f1$draws, f3$draws)
  expect_true(is.numeric(f1$seed))
  set.seed(5); a <- stats::runif(1)
  set.seed(5); run(seed = 3); b <- stats::runif(1)
  expect_identical(a, b)
})

test_that("a short NUTS run works everywhere, including on CRAN", {
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(1, 0.5); }")
  f <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 200, seed = 3, refresh = 0)
  expect_equal(dim(f$draws)[1:2], c(200L, 2L))
  expect_equal(mean(bayescomp_extract(f, "mu")), 1, tolerance = 0.3)
})
