# Compares bugs() with the results OpenBUGS 3.2.3 publishes for its Examples
# (validation/bugs/reference.csv, transcribed from the result tables of the
# example pages). A posterior mean passes when it is within 4 combined Monte
# Carlo standard errors of the reference mean, or within a tenth of the
# reference posterior sd. The example files are not part of the package:
# point BAYESCOMPR_BUGS_EXAMPLES at the doc/Examples folder of the OpenBUGS
# sources.
# Each model's rows are saved in BAYESCOMPR_BUGS_CACHE (default
# validation/bugs/cache) as soon as it finishes, and a rerun skips the
# models found there.
#   Rscript validation/bugs/validate-examples.R [n.keep] [cores] [models...]
suppressMessages(pkgload::load_all(quiet = TRUE))
dir <- Sys.getenv("BAYESCOMPR_BUGS_EXAMPLES")
if (!nzchar(dir)) stop("set BAYESCOMPR_BUGS_EXAMPLES")
args <- commandArgs(TRUE)
n_keep <- if (length(args) >= 1L) as.integer(args[1L]) else 5000L
cores <- if (length(args) >= 2L) as.integer(args[2L]) else 1L
ref <- utils::read.csv("validation/bugs/reference.csv", stringsAsFactors = FALSE)
models <- if (length(args) >= 3L) args[-(1:2)] else unique(ref$model)
cache <- Sys.getenv("BAYESCOMPR_BUGS_CACHE", "validation/bugs/cache")
dir.create(cache, showWarnings = FALSE, recursive = TRUE)
cached <- function(nm) file.path(cache, paste0(nm, ".rds"))

check <- function(nm) {
  r <- ref[ref$model == nm, ]
  df <- file.path(dir, paste0(nm, "data.txt"))
  inf <- list.files(dir, paste0("^", nm, "inits[0-9]*[.]txt$"), full.names = TRUE)
  t0 <- Sys.time()
  fit <- tryCatch(bugs(if (file.exists(df)) df, if (length(inf)) inf, unique(r$node), n.iter = 1000L + n_keep,
                       model.file = file.path(dir, paste0(nm, "model.txt")), n.chains = 2, n.burnin = 1000L,
                       DIC = FALSE, bugs.seed = 1),
                  error = function(e) e)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(fit, "error")) return(data.frame(model = nm, node = NA, ref = NA, mean = NA, z = NA, pass = FALSE,
                                                note = conditionMessage(fit), secs = secs))
  out <- lapply(seq_len(nrow(r)), function(k) {
    x <- fit$sims.array[, , r$node[k]]
    m <- mean(x)
    ess <- bayescomp_ess_bulk(x)
    mcse <- stats::sd(as.vector(x)) / sqrt(max(ess, 1))
    d <- abs(m - r$mean[k])
    z <- d / sqrt(mcse^2 + r$mc_error[k]^2)
    data.frame(model = nm, node = r$node[k], ref = r$mean[k], mean = signif(m, 4), z = round(z, 1),
               pass = is.finite(z) && (z < 4 || d < 0.1 * r$sd[k]), note = "", secs = secs)
  })
  do.call(rbind, out)
}

run <- function(nm) {
  if (file.exists(cached(nm))) return(invisible())
  out <- tryCatch(check(nm), error = function(e) data.frame(model = nm, node = NA, ref = NA, mean = NA, z = NA,
                                                            pass = FALSE, note = conditionMessage(e), secs = NA))
  saveRDS(out, cached(nm))
}
todo <- models[!file.exists(vapply(models, cached, ""))]
if (cores > 1L) parallel::mclapply(todo, run, mc.cores = cores, mc.preschedule = FALSE) else lapply(todo, run)
res <- do.call(rbind, lapply(models, function(nm) if (file.exists(cached(nm))) readRDS(cached(nm))
                             else data.frame(model = nm, node = NA, ref = NA, mean = NA, z = NA, pass = FALSE,
                                             note = "no result", secs = NA)))
utils::write.csv(res, "validation/bugs/validation-results.csv", row.names = FALSE)
for (nm in models) {
  x <- res[res$model == nm, ]
  bad <- x[!x$pass, ]
  cat(sprintf("%-13s %6.0fs  %s\n", nm, x$secs[1L],
              if (!nrow(bad)) sprintf("pass (%d nodes)", nrow(x))
              else paste("FAIL", paste(sprintf("%s: %s vs %s (z %s)%s", bad$node, bad$mean, bad$ref, bad$z,
                                               ifelse(nzchar(bad$note), paste0(" ", bad$note), "")), collapse = "; "))))
}
cat(sprintf("\n%d of %d models pass\n", sum(tapply(res$pass, res$model, all)), length(models)))
