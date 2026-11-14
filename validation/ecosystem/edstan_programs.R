# edstan's shipped IRT programs with small data: Rscript edstan_programs.R <outdir>
suppressMessages(library(edstan))
out <- commandArgs(TRUE)[1]
ex <- system.file("extdata", package = "edstan")
sp <- spelling[1:200, ]
ag <- aggression[aggression$person <= 60, ]
bin <- irt_data(response_matrix = sp[, 2:5], covariates = sp[, "male", drop = FALSE], formula = ~ male)
poly <- irt_data(y = ag$poly, ii = ag$item, jj = ag$person, covariates = ag[, c("male", "anger")], formula = ~ male + anger)
for (f in list.files(ex, "\\.stan$")) {
  d <- if (grepl("^(rasch|2pl)", f)) bin else poly
  if (!grepl("latent_reg", f)) d <- d[setdiff(names(d), c("K", "W"))]
  saveRDS(list(code = paste(readLines(file.path(ex, f)), collapse = "\n"), data = d), file.path(out, sub("\\.stan$", ".rds", f)))
}
