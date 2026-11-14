suppressMessages(library(rstanarm))
OUTDIR <- commandArgs(TRUE)[1]
CAPNAME <- NULL
trace(rstan::sampling, signature = "stanmodel", tracer = quote({
  saveRDS(list(code = object@model_code, data = data), file.path(OUTDIR, paste0(CAPNAME, ".rds"))); stop("captured")
}), print = FALSE, where = asNamespace("rstan"))
cap <- function(name, expr) {
  CAPNAME <<- name
  r <- try(expr, silent = TRUE)
  cat(name, file.exists(file.path(OUTDIR, paste0(name, ".rds"))), if (!file.exists(file.path(OUTDIR, paste0(name, ".rds")))) substr(r, 1, 120), "\n")
}
set.seed(1)
cap("lm", stan_lm(mpg ~ wt + hp, data = mtcars, prior = R2(0.5, what = "mean")))
cap("polr", stan_polr(factor(gear) ~ wt + hp, data = mtcars, prior = R2(0.3, what = "mean")))
cap("betareg", stan_betareg(qs ~ wt, data = transform(mtcars, qs = qsec / 30)))
cap("glm_binomial2", stan_glm(cbind(carb, 10L - carb) ~ wt, data = mtcars, family = binomial()))
cap("glm_invgauss", stan_glm(mpg ~ wt, data = mtcars, family = inverse.gaussian()))
cap("aov", stan_aov(mpg ~ factor(cyl), data = mtcars, prior = R2(0.5, what = "mean")))
