suppressMessages(library(brms))
out <- commandArgs(TRUE)[1]
set.seed(1)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:6, each = 10)), t = rep(1:10, 6))
d$y <- 1 + 0.5 * d$x + rnorm(6)[as.integer(d$g)] + rnorm(n)
d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x))
d$bin <- rbinom(n, 1, plogis(0.2 + d$x))
d$trials <- 10L; d$succ <- rbinom(n, 10, plogis(-0.3 + 0.5 * d$x))
d$ord <- cut(d$y, c(-Inf, 0, 1, 2, Inf), labels = FALSE)
d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
d$zi <- ifelse(runif(n) < 0.3, 0L, d$cnt)
d$pos <- exp(d$y / 3)
d$prop <- plogis(d$y / 2)
d$cens <- ifelse(d$pos > 2, 1L, 0L); d$posc <- pmin(d$pos, 2)
d$y2 <- d$y + rnorm(n)
mods <- list(
  gaussian = bf(y ~ x),
  gaussian_ri = bf(y ~ x + (1 | g)),
  gaussian_rs = bf(y ~ x + (1 + x | g)),
  student = bf(y ~ x),
  poisson = bf(cnt ~ x),
  negbin = bf(cnt ~ x + (1 | g)),
  bernoulli = bf(bin ~ x + z),
  binomial = bf(succ | trials(trials) ~ x),
  cumulative = bf(ord ~ x),
  categorical = bf(cat ~ x),
  zip = bf(zi ~ x),
  hurdle_pois = bf(zi ~ x),
  lognormal = bf(pos ~ x),
  gamma = bf(pos ~ x),
  weibull_cens = bf(posc | cens(cens) ~ x),
  beta = bf(prop ~ x),
  spline = bf(y ~ s(x, k = 5)),
  gp = bf(y ~ gp(x, k = 5, c = 5/4)),
  ar1 = bf(y ~ x + ar(time = t, gr = g, p = 1)),
  sigma_model = bf(y ~ x, sigma ~ z),
  nonlinear = bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
  mvgauss = bf(mvbind(y, y2) ~ x) + set_rescor(TRUE),
  mo = bf(y ~ mo(ord)),
  me = bf(y ~ me(x, 0.1)),
  mi = bf(y ~ x)
)
fams <- list(student = student(), poisson = poisson(), negbin = negbinomial(), bernoulli = bernoulli(),
  binomial = binomial(), cumulative = cumulative(), categorical = categorical(), zip = zero_inflated_poisson(),
  hurdle_pois = hurdle_poisson(), lognormal = lognormal(), gamma = Gamma(link = "log"), weibull_cens = weibull(), beta = Beta())
pri <- list(nonlinear = prior(normal(1, 2), nlpar = "b1") + prior(normal(0, 2), nlpar = "b2"))
for (nm in names(mods)) {
  fam <- if (is.null(fams[[nm]])) gaussian() else fams[[nm]]
  r <- try({
    code <- make_stancode(mods[[nm]], data = d, family = fam, prior = pri[[nm]])
    sd <- make_standata(mods[[nm]], data = d, family = fam, prior = pri[[nm]])
    saveRDS(list(code = as.character(code), data = unclass(sd)), file.path(out, paste0(nm, ".rds")))
  }, silent = TRUE)
  cat(nm, if (inherits(r, "try-error")) paste("GEN-ERR", substr(r, 1, 100)) else "ok", "\n")
}
