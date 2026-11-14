suppressMessages(library(multinma))
out <- commandArgs(TRUE)[1]
CAPNAME <- NULL
trace(rstan::sampling, signature = "stanmodel", tracer = quote({
  saveRDS(list(code = object@model_code, data = data), file.path(OUTDIR, paste0(CAPNAME, ".rds"))); stop("captured")
}), print = FALSE, where = asNamespace("rstan"))
OUTDIR <- out
cap <- function(name, expr) {
  CAPNAME <<- name
  r <- try(expr, silent = TRUE)
  cat(name, file.exists(file.path(out, paste0(name, ".rds"))), if (inherits(r, "try-error")) substr(r, 1, 400), "\n")
}

# binary outcome, fixed and random effects (smoking cessation)
net <- set_agd_arm(smoking, study = studyn, trt = trtc, r = r, n = n, trt_ref = "No intervention")
cap("smk_fe", nma(net, trt_effects = "fixed", prior_intercept = normal(scale = 100), prior_trt = normal(scale = 100)))
cap("smk_re", nma(net, trt_effects = "random", prior_intercept = normal(scale = 100), prior_trt = normal(scale = 100), prior_het = normal(scale = 5)))
# continuous contrast data (Parkinson's)
pk <- set_agd_contrast(parkinsons, study = studyn, trt = trtn, y = diff, se = se_diff, sample_size = n)
cap("park_fe", nma(pk, trt_effects = "fixed", prior_trt = normal(scale = 100)))
# rate outcome (Poisson, dietary fat)
df <- set_agd_arm(dietary_fat, study = studyc, trt = trtc, r = r, E = E, trt_ref = "Control", sample_size = n)
cap("fat_re", nma(df, trt_effects = "random", prior_intercept = normal(scale = 10), prior_trt = normal(scale = 10), prior_het = half_normal(scale = 5)))
sv <- set_ipd(ndmm_ipd[ndmm_ipd$studyf %in% c("Attal2012", "McCarthy2012"), ], study = studyf, trt = trtf, Surv = survival::Surv(eventtime / 12, status))
cap("surv_mspline", nma(sv, likelihood = "mspline", n_knots = 3, prior_intercept = normal(scale = 10), prior_trt = normal(scale = 10)))
cap("surv_weibull", nma(sv, likelihood = "weibull", prior_intercept = normal(scale = 10), prior_trt = normal(scale = 10), prior_aux = half_normal(scale = 10)))
# ordered categorical (hta_psoriasis)
hp <- set_agd_arm(hta_psoriasis[hta_psoriasis$studyc %in% unique(hta_psoriasis$studyc)[1:6], ], study = paste(studyc, year), trt = trtc, r = multi(r0 = sample_size - rowSums(cbind(PASI50, PASI75, PASI90), na.rm = TRUE), PASI50, PASI75, PASI90, inclusive = FALSE, type = "ordered"))
cap("ordered_fe", nma(hp, trt_effects = "fixed", link = "probit", prior_intercept = normal(scale = 10), prior_trt = normal(scale = 10), prior_aux = flat()))
# ML-NMR, Bernoulli with numerical integration (plaque psoriasis, small)
pso_ipd <- plaque_psoriasis_ipd[plaque_psoriasis_ipd$studyc %in% c("UNCOVER-1"), ]
pso_agd <- plaque_psoriasis_agd[plaque_psoriasis_agd$studyc == "FIXTURE", ]
pso_ipd <- transform(pso_ipd, bsa = bsa / 100, weight = weight / 10, durnpso = durnpso / 10)
pso_ipd <- pso_ipd[complete.cases(pso_ipd[, c("bsa", "weight", "durnpso", "pasi75")]), ]
pso_agd <- transform(pso_agd, bsa_mean = bsa_mean / 100, bsa_sd = bsa_sd / 100, weight_mean = weight_mean / 10, weight_sd = weight_sd / 10, durnpso_mean = durnpso_mean / 10, durnpso_sd = durnpso_sd / 10)
pnet <- combine_network(set_ipd(pso_ipd, study = studyc, trt = trtc, r = pasi75), set_agd_arm(pso_agd, study = studyc, trt = trtc, r = pasi75_r, n = pasi75_n), trt_ref = "PBO")
pnet <- add_integration(pnet, durnpso = distr(qgamma, mean = durnpso_mean, sd = durnpso_sd), bsa = distr(qlogitnorm, mean = bsa_mean, sd = bsa_sd), weight = distr(qgamma, mean = weight_mean, sd = weight_sd), n_int = 16)
cap("mlnmr", nma(pnet, trt_effects = "fixed", link = "probit", likelihood = "bernoulli2", regression = ~(durnpso + bsa + weight) * .trt, prior_intercept = normal(scale = 10), prior_trt = normal(scale = 10), prior_reg = normal(scale = 10), QR = TRUE))
