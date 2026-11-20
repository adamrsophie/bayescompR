# Posterior and interface checks with exact answers

run_jags <- function(code, data, par, n.iter = 6000, seed = 1, ...)
  jags(data, NULL, par, code, n.chains = 2, n.iter = n.iter, jags.seed = seed, DIC = FALSE, ...)$BUGSoutput$sims.list

test_that("loop-mode statements get their own parameters", {
  s <- run_jags("model { x ~ dnorm(mean(v[]), 1) T(0, 4)  z ~ dcat(p[]) }", list(v = c(1, 3), p = c(0.2, 0.8)), c("x", "z"))
  expect_equal(mean(s$x), 2, tolerance = 0.04)
  expect_equal(mean(s$z == 2), 0.8, tolerance = 0.03)
})

test_that("a child of one block member does not change the others", {
  s <- run_jags("model { for (i in 1:2) { z[i] ~ dbern(0.5) }  y ~ dnorm(2 * z[1], 1) }", list(y = 2), "z")
  p1 <- stats::dnorm(2, 2) / (stats::dnorm(2, 2) + stats::dnorm(2, 0))
  expect_equal(colMeans(s$z), c(p1, 0.5), tolerance = 0.03)
})

test_that("the slice sampler leaves a uniform target invariant", {
  set.seed(2)
  z0 <- stats::runif(40000, 0, 200)
  f <- function(z) ifelse(z >= 0 & z <= 200, 0, -Inf)
  z1 <- asNamespace("bayescompR")$bugs_slice_core(f, z0, rep(1, length(z0)), m = 50L)
  expect_equal(mean(z1 < 10), 0.05, tolerance = 0.12)
})

test_that("NIMBLE parameters are converted before samplers are chosen", {
  m <- function() { tau ~ dgamma(2, 1); y ~ dnorm(0, sd = tau) }
  s <- jags(list(y = 0), NULL, "tau", m, n.chains = 2, n.iter = 8000, jags.seed = 3, DIC = FALSE)$BUGSoutput$sims.list
  expect_equal(mean(s$tau), 1, tolerance = 0.05)
  expect_error(jags(list(y = 0), NULL, "mu", function() { mu ~ dnorm(0, 1); y ~ dnorm(mu, sdd = 2) }, n.iter = 10),
               "no parameter 'sdd'")
  M <- bugs_compile(bugs_from_r(function() { y ~ dt(mu = 3, sigma = 2, df = 5) }), list(y = 3))
  lp <- M$stmts[[1L]]$fn$lp(1L)
  expect_equal(lp, stats::dt(0, 5, log = TRUE) - log(2), tolerance = 1e-8)
})

test_that("cut() stops feedback", {
  s <- run_jags("model { x ~ dnorm(0, 1)  z <- cut(x)  y ~ dnorm(z, 1) }", list(y = 10), "x")
  expect_equal(mean(s$x), 0, tolerance = 0.08)
  expect_error(jags(list(y = 1), NULL, "x", "model { x ~ dnorm(0, 1)  y ~ dnorm(2 * cut(x), 1) }", n.iter = 10), "whole right side")
})

test_that("truncation normalizers are accurate in the tails and for integers", {
  M <- bugs_compile(bugs_parse("model { x ~ dnorm(0, 1) T(10, )  k ~ dpois(3) T(0.5, ) }"), list(x = 11, k = 2))
  expect_equal(M$stmts[[1L]]$fn$lp(1L), stats::dnorm(11, log = TRUE) - stats::pnorm(10, lower.tail = FALSE, log.p = TRUE), tolerance = 1e-8)
  expect_equal(M$stmts[[2L]]$fn$lp(1L), stats::dpois(2, 3, log = TRUE) - log(1 - stats::dpois(0, 3)), tolerance = 1e-8)
  s <- run_jags("model { x ~ dpois(1) T(3, 4) }", list(), "x", n.iter = 2000)
  expect_true(all(s$x %in% 3:4))
})

test_that("censoring follows OpenBUGS for observed values", {
  M <- bugs_compile(bugs_parse("model { y ~ dnorm(0, 1) C(1, ) }"), list(y = 0), "openbugs")
  expect_equal(M$stmts[[1L]]$fn$lp(1L), stats::dnorm(0, log = TRUE))
})

test_that("graphs: cycles, latent indices and observed values", {
  expect_error(bugs_compile(bugs_parse("model { x ~ dnorm(x, 1) }"), list()), "cycle")
  s <- run_jags("model { k ~ dbern(0.5)  y ~ dnorm(mu[k + 1], 1) }", list(mu = c(0, 1), y = 1), "k")
  p1 <- stats::dnorm(1, 1) / (stats::dnorm(1, 1) + stats::dnorm(1, 0))
  expect_equal(mean(s$k), p1, tolerance = 0.03)
  s <- run_jags("model { y ~ dbin(p, n)  n <- x + y  p ~ dbeta(1, 1) }", list(x = 4, y = 3), "p")
  expect_equal(mean(s$p), 4 / 9, tolerance = 0.03)
  s <- run_jags("model { k ~ dbin(0.001, 101)  y ~ dnorm(mu[k + 1], 1) }", list(mu = rep(0, 102), y = 0), "k", n.iter = 1000)
  expect_true(all(s$k >= 0))
})

test_that("data transformations make observations", {
  code <- "model { a <- x + 1  z <- a  z ~ dnorm(mu, 1)  mu ~ dnorm(0, 1) }"
  b <- bugs(list(x = 9), NULL, "mu", 6000, code, n.chains = 2, bugs.seed = 4, DIC = FALSE)
  expect_equal(b$mean$mu, 5, tolerance = 0.03)
  M <- bugs_compile(bugs_parse("data { x[2] <- x[1] + 1 } model { y ~ dnorm(x[2], 1) }"), list(x = c(3, NA), y = 1))
  expect_equal(get("x", envir = M$S)[2L], 4)
})

test_that("partly observed multivariate normals are completed", {
  P <- solve(matrix(c(1, 0.5, 0.5, 1), 2))
  s <- run_jags("model { y[1:2] ~ dmnorm(mu[], P[, ]) }", list(y = c(1, NA), mu = c(0, 0), P = P), "y", n.iter = 8000)
  expect_equal(mean(s$y[, 2]), 0.5, tolerance = 0.06)
  expect_equal(stats::var(s$y[, 2]), 0.75, tolerance = 0.06)
  expect_true(all(s$y[, 1] == 1))
})

test_that("vectorized ifelse keeps per-row values", {
  s <- run_jags("model { z ~ dbern(0.5)  for (i in 1:2) { w[i] <- ifelse(z, i, 0)  y[i] ~ dnorm(w[i], 1) } }", list(y = c(0, 2)), "z")
  l1 <- stats::dnorm(0, 1) * stats::dnorm(2, 2); l0 <- stats::dnorm(0, 0) * stats::dnorm(2, 0)
  expect_equal(mean(s$z), l1 / (l1 + l0), tolerance = 0.03)
})

test_that("functions follow JAGS and OpenBUGS", {
  M <- bugs_compile(bugs_parse("model { a <- round(2.5)  b[1:2] <- rank(v[])  c <- rank(d[], 3)  m <- max(v[], 7) }"),
                    list(v = c(5, 5), d = c(0.5, -1.2, 2)))
  expect_equal(get("a", envir = M$S)[1L], 3)
  expect_equal(as.vector(get("b", envir = M$S)), c(1, 2))
  expect_equal(get("c", envir = M$S)[1L], 3)
  expect_equal(get("m", envir = M$S)[1L], 7)
  expect_error(bugs_compile(bugs_parse("model { z <- 2 %in% a[] }"), list(a = 1:3)), "operator")
})

test_that("Dirichlet elements with zero parameters stay at zero", {
  s <- run_jags("model { p[1:3] ~ ddirch(a[])  x ~ dcat(p[]) }", list(a = c(1, 0, 1), x = 1), "p")
  expect_true(all(s$p[, 2] == 0))
  expect_equal(mean(s$p[, 1]), 2 / 3, tolerance = 0.03)
})

test_that("latent multinomial counts keep their total", {
  s <- run_jags("model { x[1:3] ~ dmulti(p[], 10) }", list(p = c(0.2, 0.3, 0.5)), "x", n.iter = 4000)
  expect_true(all(rowSums(s$x) == 10))
  expect_equal(colMeans(s$x), c(2, 3, 5), tolerance = 0.06)
})

test_that("conjugate normal updates give the regression posterior", {
  set.seed(5)
  x <- c(-1, 0, 1, 2, 3); y <- c(0.1, 1.2, 1.9, 3.2, 3.8)
  code <- "model { for (i in 1:5) { mu[i] <- a + b * x[i]  y[i] ~ dnorm(mu[i], 4) }  a ~ dnorm(0, 0.01)  b ~ dnorm(0, 0.01) }"
  f <- jags(list(x = x, y = y), NULL, c("a", "b"), code, n.chains = 2, n.iter = 6000, jags.seed = 6, DIC = FALSE)
  X <- cbind(1, x); P <- diag(0.01, 2) + 4 * crossprod(X)
  m <- solve(P, 4 * crossprod(X, y))
  expect_equal(c(f$BUGSoutput$mean$a, f$BUGSoutput$mean$b), as.vector(m), tolerance = 0.03)
  expect_equal(c(f$BUGSoutput$sd$a, f$BUGSoutput$sd$b), unname(sqrt(diag(solve(P)))), tolerance = 0.06)
})

test_that("dloglik and dinterval", {
  s <- run_jags("model { mu ~ dnorm(0, 0.01)  z ~ dloglik(-0.5 * (mu - 2)^2) }", list(z = 0), "mu")
  expect_equal(mean(s$mu), 2 / (1 + 0.01), tolerance = 0.05)
  s <- run_jags("model { t ~ dnorm(0, 1)  y ~ dinterval(t, c[]) }", list(y = 1, c = c(0, 1)), "t", n.iter = 4000)
  expect_true(all(s$t > 0 & s$t <= 1))
})

test_that("R2OpenBUGS and R2jags arguments and results", {
  code <- "model { x ~ dnorm(0, 1)  y ~ dnorm(x, 1) }"
  b <- bugs(list(y = 1), NULL, "x", 100, code, n.chains = 2, n.burnin = 20, n.thin = 4, DIC = FALSE)
  expect_equal(b$n.keep, 80L)
  f <- jags(list(y = 1), NULL, c("x", "deviance"), code, n.chains = 2, n.iter = 200, DIC = FALSE)
  expect_true("deviance" %in% names(f$BUGSoutput$sims.list))
  f <- jags(list(y = c(1, 2)), NULL, "x[2]", "model { for (i in 1:2) { x[i] ~ dnorm(0, 1)  y[i] ~ dnorm(x[i], 1) } }",
            n.chains = 1, n.iter = 200, DIC = TRUE)
  expect_equal(colnames(f$BUGSoutput$sims.matrix), c("deviance", "x[2]"))
  dev <- f$BUGSoutput$sims.array[, , "deviance"]
  expect_equal(f$BUGSoutput$pD, stats::var(as.vector(dev)) / 2)
  expect_error(bugs_parse_values("y[]\n4\ntypo\nEND\n"), "not a number")
  expect_error(bugs_parse_values("list(x = rep(0, 1e12))"), "too long")
})

test_that("the rjags functions keep the chain state", {
  m <- jags.model(textConnection("model { y ~ dnorm(mu, 1)  mu ~ dnorm(0, 0.01) }"), data = list(y = 1.5), n.chains = 2,
                  n.adapt = 200, seed = 7)
  update(m, 200)
  s <- coda.samples(m, "mu", n.iter = 3000)
  expect_s3_class(s, "mcmc.list")
  expect_equal(attr(s[[1L]], "mcpar"), c(401, 3400, 1))
  expect_equal(mean(unlist(s)), 1.5 / 1.01, tolerance = 0.05)
  j <- jags.samples(m, c("mu", "deviance"), n.iter = 100)
  expect_equal(dim(j$mu), c(1L, 100L, 2L))
})

test_that("the NIMBLE functions run nimbleCode models", {
  code <- nimbleCode({
    for (i in 1:N) y[i] ~ dnorm(mu, sd = 1)
    mu ~ dnorm(0, sd = 10)
  })
  y <- c(1.1, 0.4, 1.8, 0.9, 1.3)
  out <- nimbleMCMC(code, constants = list(N = 5), data = list(y = y), niter = 4000, nburnin = 1000, nchains = 2,
                    setSeed = 1, summary = TRUE)
  expect_named(out$samples, c("chain1", "chain2"))
  expect_equal(colnames(out$samples$chain1), "mu")
  expect_equal(out$summary$all.chains["mu", "Mean"], sum(y) / (5 + 0.01), tolerance = 0.03)
  m <- nimbleModel(code, constants = list(N = 5), data = list(y = y))
  s <- runMCMC(compileNimble(buildMCMC(configureMCMC(m))), niter = 500, setSeed = 2)
  expect_equal(dim(s), c(500L, 1L))
})

test_that("mixture parameters chosen by latent indices get conjugate updates", {
  code <- "model {
    for (i in 1:4) { z[i] ~ dcat(w[])  y[i] ~ dbin(theta[z[i]], 10) }
    for (k in 1:2) { theta[k] ~ dbeta(1, 1) }
  }"
  y <- c(3, 5, 4, 6)
  M <- bugs_compile(bugs_parse(code), list(y = y, w = c(1, 0)))
  st <- bugs_chain_init(M, NULL, "typical")
  expect_true(any(vapply(st$blocks, function(b) b$sampler == "conj_dbeta" && isTRUE(b$dyn), TRUE)))
  s <- run_jags(code, list(y = y, w = c(1, 0)), "theta")
  expect_equal(mean(s$theta[, 1]), (1 + sum(y)) / (2 + 40), tolerance = 0.02)
  expect_equal(mean(s$theta[, 2]), 0.5, tolerance = 0.04)
  # cells of a stochastic node given by logical nodes
  s <- run_jags("model { Y[1] <- 1  Y[2] <- 0  Y[1:2] ~ dmulti(p[], 1)  p[1:2] ~ ddirch(a[]) }", list(a = c(1, 1)), "p")
  expect_equal(mean(s$p[, 1]), 2 / 3, tolerance = 0.03)
})

test_that("a mixture reaches Dirichlet nodes that fill part of a row", {
  # error[2, 1:2] is the node; error[2, ] is read whole through error[z[i], ]
  code <- "model {
    for (i in 1:4) { z[i] ~ dcat(w[])  y[i, 1:3] ~ dmulti(error[z[i], ], 5) }
    error[2, 1:2] ~ ddirch(a[])
  }"
  d <- list(y = rbind(c(1, 4, 0), c(2, 3, 0), c(0, 5, 0), c(1, 4, 0)), w = c(0, 1), a = c(1, 1),
            error = rbind(c(1, 0, 0), c(NA, NA, 0)))
  s <- run_jags(code, d, "error")
  expect_equal(mean(s$error[, 2, 2]), 17 / 22, tolerance = 0.02)
})

test_that("metropolized Gibbs keeps finite discrete conditionals and moves more often", {
  z <- run_jags("model { z ~ dcat(p[]) }", list(p = c(0.2, 0.5, 0.3)), "z", n.iter = 21000, n.thin = 1)$z
  expect_equal(as.numeric(prop.table(table(factor(z, 1:3)))), c(0.2, 0.5, 0.3), tolerance = 0.03)
  # plain Gibbs moves with probability 1 - sum(p^2) = 0.62
  expect_gt(mean(z[-1] != z[-length(z)]), 0.75)
})

test_that("multivariate normal nodes get elliptical slice updates with the exact posterior", {
  set.seed(2)
  t <- seq(-1, 1, length.out = 20); y <- 1 + 2 * t + rnorm(20, 0, 0.5)
  s <- run_jags("model { b[1:2] ~ dmnorm(m[], O[,])  for (i in 1:N) { y[i] ~ dnorm(b[1] + b[2] * t[i], 4) } }",
                list(y = y, t = t, N = 20, m = c(0, 0), O = diag(0.1, 2)), "b")$b
  X <- cbind(1, t); V <- solve(diag(0.1, 2) + 4 * crossprod(X))
  expect_equal(unname(colMeans(s)), as.vector(V %*% (4 * crossprod(X, y))), tolerance = 0.03)
  expect_equal(unname(apply(s, 2, sd)), unname(sqrt(diag(V))), tolerance = 0.1)
})

test_that("correlated regression coefficients get a joint normal update with the exact posterior", {
  skip_on_cran()
  set.seed(1); n <- 60; x1 <- rnorm(n); x2 <- x1 + 0.1 * rnorm(n); y <- 1 + x1 - x2 + rnorm(n)
  code <- "model { for (i in 1:n) { y[i] ~ dnorm(mu[i], 1); mu[i] <- b0 + b1 * x1[i] + b2 * x2[i] }
    b0 ~ dnorm(0, 0.001); b1 ~ dnorm(0, 0.001); b2 ~ dnorm(0, 0.001) }"
  m <- jags.model(textConnection(code), data = list(n = n, x1 = x1, x2 = x2, y = y), n.chains = 1, quiet = TRUE, seed = 3)
  expect_true("joint_normal" %in% vapply(m$state$chains[[1L]]$blocks, `[[`, "", "sampler"))
  s <- do.call(rbind, coda.samples(m, c("b0", "b1", "b2"), n.iter = 3000))
  X <- cbind(1, x1, x2); V <- solve(crossprod(X) + diag(0.001, 3))
  expect_equal(unname(colMeans(s)), as.vector(V %*% crossprod(X, y)), tolerance = 0.05)
  expect_equal(unname(apply(s, 2, sd)), unname(sqrt(diag(V))), tolerance = 0.1)
})

test_that("hidden Markov chains are updated by forward filtering backward sampling", {
  skip_on_cran()
  set.seed(2); T <- 40; P <- matrix(c(0.95, 0.05, 0.1, 0.9), 2, byrow = TRUE); mu <- c(0, 1.5)
  z <- integer(T); z[1] <- 1; for (t in 2:T) z[t] <- sample(2, 1, prob = P[z[t - 1], ]); y <- rnorm(T, mu[z], 1)
  code <- "model { z[1] ~ dcat(p0[]); for (t in 2:T) { z[t] ~ dcat(P[z[t-1], ]) }
    for (t in 1:T) { y[t] ~ dnorm(mu[z[t]], 1) } }"
  m <- jags.model(textConnection(code), data = list(T = T, P = P, mu = mu, y = y, p0 = c(0.5, 0.5)), n.chains = 1, quiet = TRUE, seed = 3)
  expect_identical(unique(vapply(m$state$chains[[1L]]$blocks, `[[`, "", "sampler")), "ffbs")
  # exact smoothing probabilities by the forward backward algorithm
  la <- matrix(0, T, 2); la[1, ] <- log(0.5) + dnorm(y[1], mu, 1, log = TRUE)
  for (t in 2:T) for (j in 1:2) la[t, j] <- dnorm(y[t], mu[j], 1, log = TRUE) + log(sum(exp(la[t - 1, ]) * P[, j]))
  lb <- matrix(0, T, 2)
  for (t in (T - 1):1) for (i in 1:2) lb[t, i] <- log(sum(P[i, ] * dnorm(y[t + 1], mu, 1) * exp(lb[t + 1, ])))
  pz <- exp(la + lb); pz <- pz[, 2] / rowSums(pz)
  s <- do.call(rbind, coda.samples(m, "z", n.iter = 3000))
  expect_lt(max(abs(colMeans(s) - 1 - pz)), 0.05)
})
