test_that("INLA matches the exact eight schools hyperparameter posterior", {
  y <- c(28, 8, -3, 7, -1, 1, 18, 12); s <- c(15, 10, 16, 11, 9, 11, 10, 18)
  m <- bayescomp_model("data { int J; vector[J] y; vector[J] s; }
    parameters { real mu; real<lower=0> tau; vector[J] theta; }
    model { mu ~ normal(0, 5); tau ~ normal(0, 5); theta ~ normal(mu, tau); y ~ normal(theta, s); }",
    data = list(J = 8, y = y, s = s))
  fit <- bayescomp_inla(m, latent = c("theta", "mu"), seed = 1)
  tg <- seq(1e-4, 30, length.out = 3000)
  stat <- vapply(tg, function(tau) {
    V <- s^2 + tau^2; P <- sum(1 / V) + 1 / 25; m1 <- sum(y / V) / P
    c(dnorm(tau, 0, 5, log = TRUE) - 0.5 * sum(log(V)) - 0.5 * sum(y^2 / V) + 0.5 * P * m1^2 - 0.5 * log(P), m1)
  }, numeric(2))
  w <- exp(stat[1, ] - max(stat[1, ])); w <- w / sum(w)
  expect_equal(mean(bayescomp_extract(fit, "tau")), sum(w * tg), tolerance = 0.3)
  expect_equal(mean(bayescomp_extract(fit, "mu")), sum(w * stat[2, ]), tolerance = 0.3)
  expect_equal(sum(fit$grid$weight), 1)
  expect_error(bayescomp_inla(m, latent = "nope"), "not a parameter")
})

test_that("neural posterior estimation recovers a normal posterior", {
  set.seed(3)
  w <- mdn_init(2, 2, 3, 5); X <- matrix(rnorm(20), 10); Tm <- matrix(rnorm(20), 10)
  r <- mdn_loglik(w, X, Tm)
  for (nm in c("W1", "Wa", "Wm", "Ws", "bs")) {
    w2 <- w; w2[[nm]][1] <- w2[[nm]][1] + 1e-6
    expect_equal((mdn_loglik(w2, X, Tm, FALSE)$value - r$value) / 1e-6, r$grad[[nm]][1], tolerance = 1e-3)
  }
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 2); }")
  sim <- function(th) mean(rnorm(10, th[["mu"]], 1))
  fit <- bayescomp_npe(m, sim, observed = 1, num_simulations = 1500, epochs = 150, seed = 1)
  post_var <- 1 / (10 + 1 / 4)
  expect_equal(mean(bayescomp_extract(fit, "mu")), 10 * post_var, tolerance = 0.1)
  expect_equal(sd(bayescomp_extract(fit, "mu")), sqrt(post_var), tolerance = 0.1)
  expect_equal(dim(fit$posterior(0, 5)), c(5L, 1L))
})

test_that("reverse diffusion Monte Carlo approximates a mixture", {
  mix <- bayescomp_model("parameters { real x; }
    model { target += log_mix(0.3, normal_lpdf(x | -3, 1), normal_lpdf(x | 3, 1)); }")
  fit <- bayescomp_diffusion(mix, num_samples = 150, steps = 30, inner_samples = 24, reference = "standard",
                             reference_scale = 3.5, seed = 1)
  x <- bayescomp_extract(fit, "x")
  expect_equal(mean(x < 0), 0.3, tolerance = 0.12)
  expect_gt(fit$ess, 1)
})
