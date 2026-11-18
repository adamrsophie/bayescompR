test_that("BUGS densities integrate to one", {
  D <- asNamespace("bayescompR")$.bugs_dist
  chk <- function(nm, ..., lo = -Inf, hi = Inf) {
    d <- D[[nm]]
    expect_equal(stats::integrate(function(x) exp(d$lpdf(x, ...)), lo, hi)$value, 1, tolerance = 1e-5, label = nm)
  }
  chk("dnorm", 1, 4); chk("dlnorm", 0, 2, lo = 0); chk("dt", 1, 2, 5); chk("dlogis", 0, 3)
  chk("ddexp", 1, 2); chk("dweib", 2, 3, lo = 0); chk("dpar", 3, 2, lo = 2)
  chk("dgen.gamma", 2, 1.5, 0.7, lo = 0); chk("dggamma", 2, 1.5, 0.7, lo = 0); chk("dgamma", 2, 3, lo = 0)
})

test_that("a normal model with conjugate priors has the exact posterior", {
  y <- c(1.2, 0.8, 1.9, 1.1, 0.4, 1.6)
  code <- "model { for (i in 1:N) { y[i] ~ dnorm(mu, 4) }  mu ~ dnorm(0, 1) }"
  f <- jags(list(y = y, N = 6), NULL, "mu", code, n.chains = 2, n.iter = 4000, jags.seed = 3, DIC = FALSE)
  prec <- 1 + 4 * 6
  m <- f$BUGSoutput$sims.list$mu
  expect_equal(mean(m), 4 * sum(y) / prec, tolerance = 0.02)
  expect_equal(stats::sd(m), 1 / sqrt(prec), tolerance = 0.05)
})

test_that("discrete nodes are enumerated exactly", {
  # z ~ dbern(0.3); y ~ dnorm(2 z, 1) with y = 1.5: P(z = 1 | y) in closed form
  code <- "model { z ~ dbern(0.3)  m <- 2 * z  y ~ dnorm(m, 1) }"
  f <- jags(list(y = 1.5), NULL, "z", code, n.chains = 2, n.iter = 6000, jags.seed = 4, DIC = FALSE)
  p1 <- 0.3 * stats::dnorm(1.5, 2); p0 <- 0.7 * stats::dnorm(1.5, 0)
  expect_equal(mean(f$BUGSoutput$sims.list$z), p1 / (p1 + p0), tolerance = 0.03)
})

test_that("missing data are sampled and bounds are respected", {
  code <- "model { for (i in 1:3) { y[i] ~ dnorm(0, 1)T(0, ) } }"
  f <- jags(list(y = c(0.5, NA, NA)), NULL, "y", code, n.chains = 1, n.iter = 4000, jags.seed = 5, DIC = FALSE)
  s <- f$BUGSoutput$sims.list$y
  expect_true(all(s[, 2:3] > 0))
  expect_true(all(s[, 1] == 0.5))
  expect_equal(mean(s[, 2]), sqrt(2 / pi), tolerance = 0.06)
})

test_that("models written as R functions and OpenBUGS output work", {
  m <- function() {
    for (i in 1:N) { y[i] ~ dpois(lambda) }
    lambda ~ dgamma(1, 1)
  }
  y <- c(2, 4, 3, 5)
  b <- bugs(list(y = y, N = 4), NULL, "lambda", m, n.chains = 2, n.iter = 4000, seed = 6)
  expect_s3_class(b, "bugs")
  # gamma(1 + sum y, 1 + N) posterior
  expect_equal(b$mean$lambda, 15 / 5, tolerance = 0.03)
  expect_true(is.finite(b$DIC))
  expect_equal(dim(b$sims.array)[2L], 2L)
})

test_that("unsupported features stop with a clear error", {
  expect_error(jags(list(), NULL, "x", "model { x ~ dfoo(1) }", n.iter = 10), "dfoo", class = "bayescomp_unsupported_error")
  expect_error(jags(list(), NULL, "x", "model { x ~ dnorm(y, 1)  y ~ dnorm(x, 1) }", n.iter = 10), "cycle")
  expect_error(jags(list(x = 1), NULL, "x", "model { x <- 2 }", n.iter = 10), "logical node")
})
