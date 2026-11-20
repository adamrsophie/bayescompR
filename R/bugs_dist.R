# Distributions of the BUGS language (JAGS and OpenBUGS), in their BUGS
# parameterizations: normal and related distributions take a precision,
# gammas a rate, dweib(v, lambda) has density v lambda x^(v - 1)
# exp(-lambda x^v), dnegbin(p, r) counts failures before the r-th success.
#
# Each entry has
#   lpdf(x, ...)   log density, vectorized over the elements of a scalar
#                  distribution (returns one value per element);
#   rng(n, ...)    n draws;
#   support        "real", "pos" (x > 0), "unit" (0 < x < 1), "int" (any
#                  integer >= 0), "finite" (an integer in lo..hi),
#                  "simplex", "spd" (positive definite matrix), "mvreal";
#   bounds(...)    for supports with parameter dependent limits, a list of
#                  the lower and upper limits;
#   dim            "scalar" or "vector" / "matrix" for multivariate ones;
#   typical(...)   the "typical value" JAGS starts from when no initial
#                  value is given.
# Multivariate distributions take whole vectors and matrices.

.bugs_dist <- new.env(parent = emptyenv())

bdist <- function(name, lpdf, rng, support, typical, bounds = NULL, dim = "scalar", npar, aliases = character(), dialect = "both", discrete = FALSE) {
  d <- list(name = name, lpdf = lpdf, rng = rng, support = support, typical = typical, bounds = bounds, dim = dim,
            npar = npar, dialect = dialect, discrete = discrete)
  for (nm in c(name, aliases)) assign(nm, d, envir = .bugs_dist)
}

sd_of <- function(tau) 1 / sqrt(tau)
ints_ok <- function(x) all(is.finite(x)) && all(x == round(x))
# outside the support or with invalid parameters the log density is -Inf;
# a NaN from R's density functions means invalid parameters
fin <- function(v) { v[is.nan(v)] <- -Inf; v }

# ---- discrete -----------------------------------------------------------------------
bdist("dbern", function(x, p) fin(stats::dbinom(x, 1, p, log = TRUE)), function(n, p) stats::rbinom(n, 1, p),
      "finite", function(p) as.numeric(p >= 0.5), bounds = function(p) list(0, 1), npar = 1L, discrete = TRUE)
bdist("dbin", function(x, p, n) fin(stats::dbinom(x, n, p, log = TRUE)), function(m, p, n) stats::rbinom(m, n, p),
      "finite", function(p, n) round(n * p), bounds = function(p, n) list(0, n), npar = 2L, aliases = "dbinom", discrete = TRUE)
bdist("dpois", function(x, lambda) fin(stats::dpois(x, lambda, log = TRUE)), function(n, lambda) stats::rpois(n, lambda),
      "int", function(lambda) round(lambda), npar = 1L, discrete = TRUE)
bdist("dnegbin", function(x, p, r) fin(stats::dnbinom(x, r, p, log = TRUE)), function(n, p, r) stats::rnbinom(n, r, p),
      "int", function(p, r) round(r * (1 - p) / p), npar = 2L, aliases = "dnbinom", discrete = TRUE)
bdist("dgeom", function(x, p) fin(stats::dgeom(x - 1, p, log = TRUE)), function(n, p) stats::rgeom(n, p) + 1,
      "int", function(p) round(1 / p), npar = 1L, dialect = "openbugs", discrete = TRUE)
bdist("dgeom0", function(x, p) fin(stats::dgeom(x, p, log = TRUE)), function(n, p) stats::rgeom(n, p),
      "int", function(p) round((1 - p) / p), npar = 1L, dialect = "openbugs", discrete = TRUE)
# dcat(pi[]): probabilities that need not sum to one
bdist("dcat", function(x, p) { p <- p / sum(p); if (length(x) != 1L || x < 1 || x > length(p) || x != round(x)) -Inf else log(p[x]) },
      function(n, p) sample.int(length(p), n, replace = TRUE, prob = p), "finite", function(p) which.max(p),
      bounds = function(p) list(1, length(p)), npar = 1L, dim = "scalar_vecpar", discrete = TRUE)

# ---- continuous ---------------------------------------------------------------------
bdist("dnorm", function(x, mu, tau) fin(stats::dnorm(x, mu, sd_of(tau), log = TRUE)), function(n, mu, tau) stats::rnorm(n, mu, sd_of(tau)),
      "real", function(mu, tau) mu, npar = 2L)
bdist("dlnorm", function(x, mu, tau) fin(stats::dlnorm(x, mu, sd_of(tau), log = TRUE)), function(n, mu, tau) stats::rlnorm(n, mu, sd_of(tau)),
      "pos", function(mu, tau) exp(mu), npar = 2L)
bdist("dt", function(x, mu, tau, k) fin(stats::dt((x - mu) * sqrt(tau), k, log = TRUE) + 0.5 * log(tau)),
      function(n, mu, tau, k) mu + stats::rt(n, k) / sqrt(tau), "real", function(mu, tau, k) mu, npar = 3L)
bdist("dlogis", function(x, mu, tau) fin(stats::dlogis(x, mu, 1 / tau, log = TRUE)), function(n, mu, tau) stats::rlogis(n, mu, 1 / tau),
      "real", function(mu, tau) mu, npar = 2L)
bdist("ddexp", function(x, mu, tau) fin(log(tau / 2) - tau * abs(x - mu) + ifelse(tau > 0, 0, NaN)),
      function(n, mu, tau) mu + ifelse(stats::runif(n) < 0.5, -1, 1) * stats::rexp(n, tau), "real", function(mu, tau) mu, npar = 2L)
bdist("dunif", function(x, a, b) fin(stats::dunif(x, a, b, log = TRUE)), function(n, a, b) stats::runif(n, a, b),
      "interval", function(a, b) (a + b) / 2, bounds = function(a, b) list(a, b), npar = 2L)
bdist("dbeta", function(x, a, b) fin(stats::dbeta(x, a, b, log = TRUE)), function(n, a, b) stats::rbeta(n, a, b),
      "unit", function(a, b) a / (a + b), npar = 2L)
bdist("dgamma", function(x, r, mu) fin(stats::dgamma(x, r, mu, log = TRUE)), function(n, r, mu) stats::rgamma(n, r, mu),
      "pos", function(r, mu) r / mu, npar = 2L)
bdist("dexp", function(x, lambda) fin(stats::dexp(x, lambda, log = TRUE)), function(n, lambda) stats::rexp(n, lambda),
      "pos", function(lambda) 1 / lambda, npar = 1L)
bdist("dchisqr", function(x, k) fin(stats::dchisq(x, k, log = TRUE)), function(n, k) stats::rchisq(n, k),
      "pos", function(k) k, npar = 1L, aliases = "dchisq")
bdist("dweib", function(x, v, lambda) fin(ifelse(x > 0 & v > 0 & lambda > 0, log(v) + log(lambda) + (v - 1) * log(x) - lambda * x^v, -Inf)),
      function(n, v, lambda) stats::rweibull(n, v, lambda^(-1 / v)), "pos", function(v, lambda) lambda^(-1 / v) * gamma(1 + 1 / v), npar = 2L, aliases = "dweibull")
bdist("dpar", function(x, alpha, c) fin(ifelse(x >= c & alpha > 0 & c > 0, log(alpha) + alpha * log(c) - (alpha + 1) * log(x), -Inf)),
      function(n, alpha, c) c * stats::runif(n)^(-1 / alpha), "interval", function(alpha, c) c * 2^(1 / alpha),
      bounds = function(alpha, c) list(c, Inf), npar = 2L)
bdist("dflat", function(x) rep(0, length(x)), function(n) stats::rnorm(n), "real", function() 0, npar = 0L, dialect = "openbugs")
# JAGS: dgen.gamma(r, lambda, b), density b lambda^(b r) x^(b r - 1) exp(-(lambda x)^b) / Gamma(r)
bdist("dgen.gamma", function(x, r, lambda, b) fin(ifelse(x > 0 & r > 0 & lambda > 0 & b > 0, log(b) + b * r * log(lambda) + (b * r - 1) * log(x) - (lambda * x)^b - lgamma(r), -Inf)),
      function(n, r, lambda, b) stats::rgamma(n, r)^(1 / b) / lambda, "pos", function(r, lambda, b) r^(1 / b) / lambda, npar = 3L, dialect = "jags")
# OpenBUGS: dggamma(r, mu, beta), density beta mu^(beta r) x^(beta r - 1) exp(-(mu x)^beta) / Gamma(r)
bdist("dggamma", function(x, r, mu, beta) fin(ifelse(x > 0 & r > 0 & mu > 0 & beta > 0, log(beta) + beta * r * log(mu) + (beta * r - 1) * log(x) - (mu * x)^beta - lgamma(r), -Inf)),
      function(n, r, mu, beta) stats::rgamma(n, r)^(1 / beta) / mu, "pos", function(r, mu, beta) r^(1 / beta) / mu, npar = 3L, dialect = "openbugs")
bdist("df", function(x, n, m) fin(stats::df(x, n, m, log = TRUE)), function(k, n, m) stats::rf(k, n, m), "pos", function(n, m) 1, npar = 2L, dialect = "jags")
bdist("dnchisqr", function(x, k, delta) fin(stats::dchisq(x, k, delta, log = TRUE)), function(n, k, delta) stats::rchisq(n, k, delta),
      "pos", function(k, delta) k + delta, npar = 2L, dialect = "jags")
bdist("dnt", function(x, mu, tau, k) fin(stats::dt(x * sqrt(tau), k, ncp = mu * sqrt(tau), log = TRUE) + 0.5 * log(tau)),
      function(n, mu, tau, k) stats::rt(n, k, ncp = mu * sqrt(tau)) / sqrt(tau), "real", function(mu, tau, k) mu, npar = 3L, dialect = "jags")

# ---- multivariate -------------------------------------------------------------------
lmvn_prec <- function(x, mu, Omega) {
  R <- tryCatch(chol(Omega), error = function(e) NULL)
  if (is.null(R)) return(-Inf)
  z <- R %*% (x - mu)
  sum(log(diag(R))) - 0.5 * sum(z * z) - 0.5 * length(x) * log(2 * pi)
}
rmvn_prec <- function(n, mu, Omega) {
  R <- chol(Omega)
  as.vector(mu + backsolve(R, stats::rnorm(length(mu))))
}
bdist("dmnorm", function(x, mu, Omega) lmvn_prec(x, mu, Omega), function(n, mu, Omega) rmvn_prec(n, mu, Omega),
      "mvreal", function(mu, Omega) mu, dim = "vector", npar = 2L)
bdist("dmnorm.vcov", function(x, mu, Sigma) { O <- tryCatch(solve(Sigma), error = function(e) NULL); if (is.null(O)) -Inf else lmvn_prec(x, mu, O) },
      function(n, mu, Sigma) as.vector(mu + t(chol(Sigma)) %*% stats::rnorm(length(mu))), "mvreal", function(mu, Sigma) mu, dim = "vector", npar = 2L, dialect = "jags")
bdist("dmt", function(x, mu, Omega, k) {
  p <- length(x); R <- tryCatch(chol(Omega), error = function(e) NULL)
  if (is.null(R)) return(-Inf)
  z <- R %*% (x - mu)
  lgamma((k + p) / 2) - lgamma(k / 2) - p / 2 * log(k * pi) + sum(log(diag(R))) - (k + p) / 2 * log1p(sum(z * z) / k)
}, function(n, mu, Omega, k) rmvn_prec(1, rep(0, length(mu)), Omega) / sqrt(stats::rchisq(1, k) / k) + mu,
      "mvreal", function(mu, Omega, k) mu, dim = "vector", npar = 3L)
# elements with a zero parameter are fixed at zero (JAGS)
bdist("ddirch", function(x, alpha) {
  if (any(alpha < 0) || all(alpha == 0) || any(x < 0) || abs(sum(x) - 1) > 1e-6 || any(x[alpha == 0] != 0)) return(-Inf)
  a <- alpha > 0
  lgamma(sum(alpha[a])) - sum(lgamma(alpha[a])) + sum(ifelse(alpha[a] == 1, 0, (alpha[a] - 1) * log(x[a])))
}, function(n, alpha) { g <- stats::rgamma(length(alpha), alpha); g[alpha == 0] <- 0; g / sum(g) }, "simplex", function(alpha) alpha / sum(alpha),
      dim = "vector", npar = 1L, aliases = "ddirich")
bdist("dmulti", function(x, p, N) {
  p <- p / sum(p)
  if (sum(x) != N || any(x < 0) || !ints_ok(x)) return(-Inf)
  fin(stats::dmultinom(x, N, p, log = TRUE))
}, function(n, p, N) as.vector(stats::rmultinom(1, N, p)), "multi", function(p, N) {
  # largest remainders, so that the counts sum to N
  e <- N * p / sum(p); x <- floor(e)
  k <- N - sum(x)
  if (k > 0) x[order(e - x, decreasing = TRUE)[seq_len(k)]] <- x[order(e - x, decreasing = TRUE)[seq_len(k)]] + 1
  x
}, dim = "vector", npar = 2L, discrete = TRUE)
# dwish(R, k): the density of a precision matrix with mean k inverse(R)
bdist("dwish", function(x, R, k) {
  p <- nrow(x)
  cx <- tryCatch(chol(x), error = function(e) NULL)
  if (is.null(cx)) return(-Inf)
  cr <- chol(R)
  ldx <- 2 * sum(log(diag(cx))); ldr <- 2 * sum(log(diag(cr)))
  k / 2 * ldr + (k - p - 1) / 2 * ldx - 0.5 * sum(R * x) - k * p / 2 * log(2) - (p * (p - 1) / 4 * log(pi) + sum(lgamma((k + 1 - seq_len(p)) / 2)))
}, function(n, R, k) stats::rWishart(1, k, solve(R))[, , 1], "spd", function(R, k) k * solve(R), dim = "matrix", npar = 2L)

# OpenBUGS: dloglik(phi) adds phi to the log likelihood (the value of the
# node is ignored); JAGS: dinterval(t, c[]) is the number of cut points
# below t, as used for censored data
bdist("dloglik", function(x, phi) phi, function(n, phi) rep(0, n), "real", function(phi) 0, npar = 1L, dialect = "openbugs")
bdist("dinterval", function(x, t, cut) if (length(x) == 1L && x == sum(t > cut)) 0 else -Inf,
      function(n, t, cut) sum(t > cut), "finite", function(t, cut) sum(t > cut), bounds = function(t, cut) list(0, length(cut)),
      npar = 2L, dim = "scalar_vecpar", discrete = TRUE, dialect = "jags")
# OpenBUGS evaluates the Poisson and binomial likelihoods for non integer
# observations too
.bugs_dist[["dpois@openbugs"]] <- utils::modifyList(.bugs_dist$dpois, list(lpdf = function(x, lambda)
  fin(ifelse(x >= 0 & lambda >= 0, ifelse(x == 0, 0, x * log(lambda)) - lambda - lgamma(x + 1), -Inf))))
.bugs_dist[["dbin@openbugs"]] <- utils::modifyList(.bugs_dist$dbin, list(lpdf = function(x, p, n)
  fin(ifelse(x >= 0 & x <= n & p >= 0 & p <= 1, lgamma(n + 1) - lgamma(x + 1) - lgamma(n - x + 1) +
             ifelse(x == 0, 0, x * log(p)) + ifelse(n == x, 0, (n - x) * log1p(-p)), -Inf))))

bugs_dist_get <- function(name, dialect = "jags", loc = NULL) {
  if (dialect == "openbugs" && !is.null(.bugs_dist[[paste0(name, "@openbugs")]])) return(.bugs_dist[[paste0(name, "@openbugs")]])
  d <- .bugs_dist[[name]]
  if (is.null(d)) stop_unsupported(sprintf("the distribution '%s' is not supported", name), loc)
  d
}
