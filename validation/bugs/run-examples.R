# Runs the OpenBUGS Examples (volumes 1 to 3) with bugs() and reports which
# compile and sample. The example files are not part of the package: point
# BAYESCOMPR_BUGS_EXAMPLES at the doc/Examples folder of the OpenBUGS 3.2.3
# sources.
#   Rscript validation/bugs/run-examples.R [n.iter]
suppressMessages(pkgload::load_all(quiet = TRUE))
dir <- Sys.getenv("BAYESCOMPR_BUGS_EXAMPLES")
if (!nzchar(dir)) stop("set BAYESCOMPR_BUGS_EXAMPLES")
args <- commandArgs(TRUE)
n_iter <- if (length(args)) as.integer(args[1L]) else 1000L
models <- sort(list.files(dir, "model\\.txt$"))
res <- lapply(models, function(mf) {
  nm <- sub("model\\.txt$", "", mf)
  df <- file.path(dir, paste0(nm, "data.txt"))
  inf <- list.files(dir, paste0("^", nm, "inits[0-9]*\\.txt$"), full.names = TRUE)
  t0 <- Sys.time()
  out <- tryCatch({
    b <- bugs(if (file.exists(df)) df else NULL, if (length(inf)) inf else NULL, NULL, n_iter, file.path(dir, mf),
              n.chains = 2, bugs.seed = 1, DIC = TRUE)
    rh <- max(b$summary[, "Rhat"], na.rm = TRUE)
    sprintf("ok   max Rhat %.2f, DIC %s", rh, if (is.null(b$DIC)) "NA" else format(round(b$DIC, 1)))
  }, error = function(e) paste("FAIL", conditionMessage(e)))
  line <- sprintf("%-14s %6.1fs  %s", nm, as.numeric(difftime(Sys.time(), t0, units = "secs")), substr(gsub("\n", " ", out), 1, 150))
  writeLines(line); flush(stdout())
  line
})
cat(sprintf("\n%d of %d models ran\n", sum(grepl("  ok ", unlist(res))), length(models)))
