# Stan CSV files: the output format of CmdStan

#' Write and read Stan CSV files
#'
#' `bayescomp_write_csv()` writes a fit in CmdStan's CSV format: comment
#' lines with the configuration, a header of column names in CmdStan's style
#' (`theta.1`, `z.real`, `t:1`), one row per draw, and for NUTS and HMC the
#' adaptation results and timing. MCMC fits give one file per chain; other
#' methods give one file. The files can be read by CmdStan's `stansummary`,
#' `cmdstanr::read_cmdstan_csv()` and `bayescomp_read_csv()`.
#'
#' `bayescomp_read_csv()` reads CSV files written by CmdStan (any method) or
#' by `bayescomp_write_csv()` into a `bayescomp_fit`. With `model`, the
#' unconstrained draws are recovered, so the fit can be used with
#' [bayescomp_generate_quantities()] and as initial values.
#'
#' @param fit A `bayescomp_fit`.
#' @param dir Directory for the files.
#' @param basename File name stem; chains add `_1`, `_2`, ... .
#' @param save_metric For NUTS and HMC fits, also write the adapted step
#'   size and inverse metric of each chain as JSON (`<file>_metric.json`, as
#'   CmdStan's `save_metric`); such a file can be given to
#'   `bayescomp_sample(inv_metric = )`. Fits sampled with
#'   `save_diagnostics = TRUE` also get `<file>_diagnostic.csv` (CmdStan's
#'   `diagnostic_file`). Programs with `profile` blocks also get
#'   `<file>_profile.csv` (CmdStan's `profile_file`) from element `profile`
#'   of the fit.
#' @param save_config Also write the run configuration as
#'   `<basename>_config.json` (CmdStan's `save_cmdstan_config`).
#' @param sig_figs Significant digits of the values (CmdStan's default is
#'   8; 17 or more round trips exactly).
#' @param files Paths of CSV files, one per chain.
#' @param model Optional `bayescomp_model` the draws belong to.
#' @return `bayescomp_write_csv()` returns the file paths, invisibly;
#'   `bayescomp_read_csv()` returns a `bayescomp_fit`.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_sample(m, chains = 1, iter_warmup = 100, iter_sampling = 100,
#'                         seed = 1, refresh = 0)
#' f <- bayescomp_write_csv(fit, tempdir())
#' fit2 <- bayescomp_read_csv(f, model = m)
#' all.equal(bayescomp_draws(fit), bayescomp_draws(fit2), tolerance = 1e-7)
#' @export
bayescomp_write_csv <- function(fit, dir = ".", basename = "output", sig_figs = 8L, save_metric = FALSE,
                                save_config = FALSE) {
  if (!inherits(fit, "bayescomp_fit")) stop("'fit' must be a bayescomp_fit", call. = FALSE)
  sig_figs <- as_count(sig_figs, "sig_figs", 1L)
  if (sig_figs > 18L) stop("'sig_figs' must be at most 18", call. = FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  chains <- dim(fit$draws)[2L]
  files <- file.path(dir, if (chains > 1L) sprintf("%s_%d.csv", basename, seq_len(chains)) else paste0(basename, ".csv"))
  for (c in seq_len(chains)) writeLines(csv_lines(fit, c, sig_figs), files[c])
  stem <- sub("\\.csv$", "", files)
  if (isTRUE(save_metric) && identical(fit$method, "sample") && length(fit$inv_metric[[1L]])) {
    for (c in seq_len(chains))
      bayescomp_write_data(list(stepsize = fit$step_size[c], metric_type = fit$config$metric %||% "diag_e", inv_metric = fit$inv_metric[[c]]),
                           paste0(stem[c], "_metric.json"))
  }
  if (!is.null(fit$profile)) {
    for (c in seq_len(chains)) {
      p <- fit$profile[fit$profile$chain == c, -1L, drop = FALSE]
      if (nrow(p)) utils::write.csv(p, paste0(stem[c], "_profile.csv"), row.names = FALSE, quote = FALSE)
    }
  }
  if (!is.null(fit$single_paths)) {
    for (i in seq_along(fit$single_paths)) {
      sp <- fit$single_paths[[i]]
      writeLines(csv_lines(sp, 1L, sig_figs), paste0(stem[1L], "_path_", i, ".csv"))
      bayescomp_write_data(list(path = i, elbo = I(sp$elbo_trace)), paste0(stem[1L], "_path_", i, ".json"))
    }
  }
  if (!is.null(fit$diagnostics) && !is.null(fit$model)) {
    un <- unconstrained_names(fit$model)
    for (c in seq_len(chains)) {
      sd <- fit$sampler_diagnostics[, c, , drop = TRUE]
      if (is.null(dim(sd))) sd <- matrix(sd, 1L, dimnames = list(NULL, dimnames(fit$sampler_diagnostics)[[3L]]))
      X <- cbind(lp__ = fit$draws[, c, "lp__"], sd, fit$unconstrained[[c]], fit$diagnostics[[c]]$p, fit$diagnostics[[c]]$g)
      body <- apply(X, 1L, function(r) paste(csv_num(r, sig_figs), collapse = ","))
      writeLines(c(csv_header(fit, c), paste(c(colnames(X)[seq_len(1L + ncol(sd))], un, paste0("p_", un), paste0("g_", un)), collapse = ","), body),
                 paste0(stem[c], "_diagnostic.csv"))
    }
  }
  if (isTRUE(save_config)) write_config_json(fit, file.path(dir, paste0(basename, "_config.json")), files, sig_figs)
  invisible(files)
}

# the run configuration in the layout of CmdStan's save_cmdstan_config
write_config_json <- function(fit, path, files, sig_figs) {
  v <- as.integer(strsplit(stan_signatures_version, ".", fixed = TRUE)[[1L]])
  cfg <- Filter(function(x) is.atomic(x) && length(x) <= 1L && !is.null(x), fit$config %||% list())
  method <- fit$method %||% "sample"
  m <- list(value = method)
  m[[method]] <- c(cfg, list(num_chains = dim(fit$draws)[2L]))
  out <- list(stan_major_version = as.character(v[1L]), stan_minor_version = as.character(v[2L]),
              stan_patch_version = as.character(v[3L]), model_name = fit$model$name %||% "bayescompR_model",
              start_datetime = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
              method = m, id = 1L, random = list(seed = fit$seed %||% -1L),
              output = list(file = files[1L], sig_figs = as.integer(sig_figs), save_cmdstan_config = TRUE),
              num_threads = 1L, stanc_version = paste0("bayescompR ", utils::packageVersion("bayescompR")))
  jsonlite::write_json(out, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
  invisible(path)
}

csv_num <- function(x, sig_figs) {
  out <- formatC(x, digits = sig_figs, format = "g")
  out[is.nan(x)] <- "nan"
  out[is.infinite(x)] <- ifelse(x[is.infinite(x)] > 0, "inf", "-inf")
  gsub(" ", "", out)
}

# CmdStan's column names from the package's: A[1,2] is A.1.2, z[real] is
# z.real, t:2[1] is t:2.1
to_cmdstan_names <- function(x) gsub("[[,]", ".", gsub("]", "", x, fixed = TRUE))

# CmdStan's names of the unconstrained parameters (theta.1, theta.2 for a
# simplex[3] theta)
unconstrained_names <- function(model) {
  unlist(lapply(model$params$pars, function(n) {
    d <- model$decls[[n]]
    k <- d$free_size
    if (k == 1L && !length(c(d$adims, d$eshape))) n else paste0(n, ".", seq_len(k))
  }))
}

# the inverse (as tuple_flat_names): indices after the first dot go in
# brackets, so at.1:2.real is at[1:2,real]
from_cmdstan_names <- function(x) {
  vapply(x, function(s) {
    p <- regexpr(".", s, fixed = TRUE)
    if (p < 0L) return(s)
    paste0(substr(s, 1L, p - 1L), "[", gsub(".", ",", substr(s, p + 1L, nchar(s)), fixed = TRUE), "]")
  }, "", USE.NAMES = FALSE)
}

# constrained parameter values from one flat draw (CmdStan's column order)
flat_count <- function(d) {
  if (identical(d$kind, "tuple")) {
    one <- sum(vapply(d$elems, flat_count, 0))
    return(one * prod(d$adims))
  }
  prod(c(d$adims, d$eshape))
}
unflatten_value <- function(v, d) {
  if (identical(d$kind, "tuple")) {
    one_tuple <- function(w) {
      k <- vapply(d$elems, flat_count, 0)
      ends <- cumsum(k)
      lapply(seq_along(d$elems), function(i) unflatten_value(w[(ends[i] - k[i]) + seq_len(k[i])], d$elems[[i]]))
    }
    if (!length(d$adims)) return(one_tuple(v))
    if (length(d$adims) > 1L) stop("reading multidimensional arrays of tuple parameters is not supported", call. = FALSE)
    per <- length(v) / d$adims
    return(lapply(seq_len(d$adims), function(i) one_tuple(v[(i - 1L) * per + seq_len(per)])))
  }
  dims <- c(d$adims, d$eshape)
  if (length(dims) > 1L) array(v, dims) else v
}
unflatten_params <- function(model, r) {
  out <- list()
  off <- 0L
  for (n in model$params$pars) {
    d <- model$decls[[n]]
    k <- flat_count(d)
    out[[n]] <- unflatten_value(r[off + seq_len(k)], d)
    off <- off + k
  }
  out
}

optim_return_code <- function(msg) {
  pats <- c("line search" = -1L, "absolute change in objective" = 20L, "relative change in objective" = 21L,
            "relative gradient" = 31L, "gradient norm" = 30L, "change in parameters" = 10L, "maximum number" = 40L,
            "no parameters" = 0L)
  for (p in names(pats)) if (grepl(p, msg, fixed = TRUE)) return(pats[[p]])
  30L
}

csv_header <- function(fit, chain) {
  cfg <- fit$config %||% list()
  m <- fit$method
  kv <- function(k, v, depth) sprintf("#%s%s = %s", strrep("  ", depth), k, paste(format(v), collapse = ","))
  lines <- c("# stan_version_major = 2", "# stan_version_minor = 40", "# stan_version_patch = 0",
             "# model = bayescompR_model", paste0("# start_datetime = ", format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"), " UTC"),
             paste0("# method = ", m), paste0("#   ", m))
  if (m == "sample") {
    alg <- cfg$algorithm %||% "nuts"
    lines <- c(lines, kv("num_samples", cfg$iter_sampling %||% dim(fit$draws)[1L], 2L),
               kv("num_warmup", cfg$iter_warmup %||% 0L, 2L), kv("save_warmup", tolower(isTRUE(cfg$save_warmup)), 2L),
               kv("thin", cfg$thin %||% 1L, 2L), "#     adapt", kv("engaged", tolower(isTRUE(cfg$adapt_engaged)), 3L),
               kv("gamma", cfg$adapt_gamma %||% 0.05, 3L), kv("delta", cfg$adapt_delta %||% 0.8, 3L),
               kv("kappa", cfg$adapt_kappa %||% 0.75, 3L), kv("t0", cfg$adapt_t0 %||% 10, 3L),
               kv("init_buffer", cfg$adapt_init_buffer %||% 75L, 3L), kv("term_buffer", cfg$adapt_term_buffer %||% 50L, 3L),
               kv("window", cfg$adapt_window %||% 25L, 3L),
               kv("algorithm", if (alg %in% c("nuts", "hmc")) "hmc" else alg, 2L))
    if (alg %in% c("nuts", "hmc"))
      lines <- c(lines, "#       hmc", kv("engine", if (alg == "nuts") "nuts" else "static", 4L),
                 if (alg == "nuts") c("#           nuts", kv("max_depth", cfg$max_treedepth %||% 10L, 6L))
                 else c("#           static", kv("int_time", cfg$int_time %||% 6.28319, 6L)),
                 kv("metric", cfg$metric %||% "diag_e", 4L), kv("stepsize", cfg$step_size %||% 1, 4L),
                 kv("stepsize_jitter", cfg$step_size_jitter %||% 0, 4L))
    lines <- c(lines, kv("num_chains", 1L, 2L))
  } else {
    for (k in names(cfg)) if (is.atomic(cfg[[k]]) && length(cfg[[k]]) == 1L && !is.na(cfg[[k]]))
      lines <- c(lines, kv(k, if (is.logical(cfg[[k]])) tolower(cfg[[k]]) else cfg[[k]], 2L))
  }
  c(lines, kv("id", chain, 0L), kv("random", "", 0L), kv("seed", fit$seed %||% 0, 1L),
    "# output", kv("sig_figs", "", 1L), "# stanc_version = bayescompR")
}

csv_lines <- function(fit, chain, sig_figs) {
  m <- fit$method
  d <- fit$draws[, chain, , drop = TRUE]
  if (is.null(dim(d))) d <- matrix(d, nrow = dim(fit$draws)[1L], dimnames = list(NULL, dimnames(fit$draws)[[3L]]))
  vars <- colnames(d)
  pnames <- setdiff(vars, c("lp__", "log_p__", "log_g__", "log_q__"))
  P <- d[, pnames, drop = FALSE]
  n <- nrow(d)
  cols <- switch(m,
    sample = {
      sd <- fit$sampler_diagnostics[, chain, , drop = TRUE]
      if (is.null(dim(sd))) sd <- matrix(sd, nrow = n, dimnames = list(NULL, dimnames(fit$sampler_diagnostics)[[3L]]))
      cbind(lp__ = d[, "lp__"], sd)
    },
    optimize = cbind(lp__ = d[, "lp__"], converged__ = c(rep(0L, n - 1L), optim_return_code(fit$message %||% ""))),
    variational = cbind(lp__ = 0, log_p__ = d[, "log_p__"], log_g__ = d[, "log_g__"]),
    pathfinder = cbind(lp_approx__ = d[, "log_q__"], lp__ = d[, "log_p__"], path__ = fit$path %||% rep(1, n)),
    laplace = cbind(log_p__ = d[, "log_p__"], log_q__ = d[, "log_g__"]),
    cbind(lp__ = d[, intersect("lp__", vars)]))
  if (m == "variational" && !is.null(fit$vi$mean) && !is.null(fit$model)) {
    # CmdStan's first row is the mean of the approximation
    mean_row <- chain_write_draws(fit$model, matrix(fit$vi$mean, 1L), FALSE)
    cols <- rbind(c(lp__ = 0, log_p__ = 0, log_g__ = 0), cols)
    P <- rbind(mean_row[1L, colnames(P)], P)
  }
  X <- cbind(cols, P)
  body <- apply(X, 1L, function(r) paste(csv_num(r, sig_figs), collapse = ","))
  hdr <- csv_header(fit, chain)
  hdr[hdr == "#   sig_figs = "] <- paste0("#   sig_figs = ", sig_figs)
  out <- c(hdr, paste(c(colnames(cols), to_cmdstan_names(colnames(P))), collapse = ","))
  nw <- fit$n_warmup_saved %||% 0L
  if (m == "sample" && isTRUE(fit$config$adapt_engaged) && fit$config$algorithm %in% c("nuts", "hmc")) {
    im <- fit$inv_metric[[chain]]
    adapt <- c("# Adaptation terminated", paste0("# Step size = ", csv_num(fit$step_size[chain], 6L)),
               if (is.matrix(im)) c("# Elements of inverse mass matrix:", apply(im, 1L, function(r) paste0("# ", paste(csv_num(r, 6L), collapse = ", "))))
               else c("# Diagonal elements of inverse mass matrix:", paste0("# ", paste(csv_num(im, 6L), collapse = ", "))))
    out <- c(out, body[seq_len(min(nw, n))], adapt, body[setdiff(seq_len(n), seq_len(nw))])
  } else {
    out <- c(out, body)
  }
  if (m == "sample" && !is.null(fit$time$chains)) {
    tm <- fit$time$chains[chain, ]
    out <- c(out, "# ", sprintf("#  Elapsed Time: %s seconds (Warm-up)", format(tm[[1L]])),
             sprintf("#                %s seconds (Sampling)", format(tm[[2L]])),
             sprintf("#                %s seconds (Total)", format(sum(tm))), "# ")
  }
  if (m == "optimize") out <- c(out, paste0("# Optimization terminated: ", fit$message %||% ""))
  out
}

#' @rdname bayescomp_write_csv
#' @export
bayescomp_read_csv <- function(files, model = NULL) {
  if (!length(files) || !all(file.exists(files))) stop("every file in 'files' must exist", call. = FALSE)
  if (!is.null(model)) check_model(model)
  parts <- lapply(files, read_csv_one)
  method <- parts[[1L]]$method
  if (!all(vapply(parts, function(p) identical(p$method, method), TRUE))) stop("the files come from different methods", call. = FALSE)
  cn <- colnames(parts[[1L]]$values)
  if (!all(vapply(parts, function(p) identical(colnames(p$values), cn), TRUE))) stop("the files have different columns", call. = FALSE)
  diag_cols <- grepl("__$", cn)
  ours <- if (!is.null(model)) {
    nm <- bayescomp_param_names(model)
    map <- stats::setNames(nm, to_cmdstan_names(nm))
    ifelse(cn[!diag_cols] %in% names(map), map[cn[!diag_cols]], from_cmdstan_names(cn[!diag_cols]))
  } else from_cmdstan_names(cn[!diag_cols])
  vi_skip <- method == "variational"
  rows <- function(p) if (vi_skip) p$values[-1L, , drop = FALSE] else p$values
  niter <- min(vapply(parts, function(p) nrow(rows(p)), 0L))
  chains <- length(parts)
  lead <- switch(method,
    sample = , optimize = "lp__",
    variational = c("lp__", "log_p__", "log_g__"),
    pathfinder = c("lp__", "log_p__", "log_q__"),
    laplace = c("lp__", "log_p__", "log_g__"),
    intersect("lp__", cn))
  draws <- array(NA_real_, c(niter, chains, length(lead) + length(ours)),
                 dimnames = list(iteration = NULL, chain = NULL, variable = c(lead, ours)))
  sdn <- if (method == "sample") setdiff(cn[diag_cols], "lp__") else character(0)
  sdiag <- if (length(sdn)) array(NA_real_, c(niter, chains, length(sdn)), dimnames = list(iteration = NULL, chain = NULL, variable = sdn)) else NULL
  for (c in seq_len(chains)) {
    v <- rows(parts[[c]])[seq_len(niter), , drop = FALSE]
    L <- switch(method,
      pathfinder = cbind(v[, "lp__"], v[, "lp__"], v[, "lp_approx__"]),
      laplace = cbind(v[, "log_p__"], v[, "log_p__"], v[, "log_q__"]),
      v[, lead, drop = FALSE])
    draws[, c, seq_along(lead)] <- L
    draws[, c, -seq_along(lead)] <- v[, !diag_cols, drop = FALSE]
    if (!is.null(sdiag)) sdiag[, c, ] <- v[, sdn, drop = FALSE]
  }
  hdr <- parts[[1L]]$config
  cfg <- list(algorithm = hdr$engine %||% hdr$algorithm, iter_sampling = as.integer(hdr$num_samples %||% niter),
              iter_warmup = as.integer(hdr$num_warmup %||% 0L), save_warmup = identical(hdr$save_warmup, "true") || identical(hdr$save_warmup, "1"),
              thin = as.integer(hdr$thin %||% 1L), max_treedepth = as.integer(hdr$max_depth %||% 10L),
              adapt_engaged = !identical(hdr$engaged, "false") && !identical(hdr$engaged, "0"), metric = hdr$metric)
  if (identical(cfg$algorithm, "static")) cfg$algorithm <- "hmc"
  nw <- if (method == "sample" && cfg$save_warmup) parts[[1L]]$n_warmup else 0L
  fit <- structure(list(
    draws = draws, sampler_diagnostics = sdiag, method = method, model = model, n_warmup_saved = nw,
    step_size = vapply(parts, function(p) p$step_size %||% NA_real_, 0),
    inv_metric = lapply(parts, `[[`, "inv_metric"), unconstrained = NULL,
    time = list(chains = do.call(rbind, lapply(parts, `[[`, "time")), total = NA_real_),
    config = cfg, seed = suppressWarnings(as.numeric(hdr$seed))), class = "bayescomp_fit")
  if (!is.null(model)) {
    np <- length(bayescomp_param_names(model, FALSE, FALSE))
    fit$unconstrained <- lapply(seq_len(chains), function(c) {
      P <- matrix(draws[, c, -seq_along(lead), drop = FALSE], niter)[, seq_len(np), drop = FALSE]
      matrix(vapply(seq_len(niter), function(i) bayescomp_unconstrain(model, unflatten_params(model, P[i, ])), numeric(model$dim)),
             niter, model$dim, byrow = TRUE)
    })
  }
  if (method == "optimize") {
    fit$par <- draws[niter, 1L, -1L]
    fit$value <- draws[niter, 1L, 1L]
  }
  fit
}

read_csv_one <- function(file) {
  ln <- readLines(file, warn = FALSE)
  com <- startsWith(ln, "#")
  hdr_i <- which(!com & nzchar(trimws(ln)))[1L]
  if (is.na(hdr_i)) stop(sprintf("'%s' has no header line", file), call. = FALSE)
  cfg_lines <- ln[seq_len(hdr_i - 1L)]
  kv <- regmatches(cfg_lines, regexec("^#\\s*([A-Za-z_0-9]+) = ([^(]*?)\\s*(\\(Default\\))?$", cfg_lines))
  kv <- kv[lengths(kv) > 0L]
  config <- list()
  for (x in kv) if (is.null(config[[x[2L]]])) config[[x[2L]]] <- x[3L]
  method <- sub(" .*", "", config$method %||% "sample")
  data_i <- which(!com & nzchar(trimws(ln)))[-1L]
  cn <- strsplit(gsub(" ", "", ln[hdr_i]), ",", fixed = TRUE)[[1L]]
  vals <- if (length(data_i)) {
    tx <- gsub(" ", "", ln[data_i])
    v <- matrix(suppressWarnings(as.numeric(unlist(strsplit(tx, ",", fixed = TRUE)))), ncol = length(cn), byrow = TRUE)
    v
  } else matrix(numeric(0), 0L, length(cn))
  colnames(vals) <- cn
  # warmup rows come before "# Adaptation terminated"
  ad <- which(ln == "# Adaptation terminated")
  n_warmup <- if (length(ad)) sum(data_i < ad[1L]) else 0L
  step <- NULL; im <- NULL
  if (length(ad)) {
    s <- grep("^# Step size = ", ln, value = TRUE)
    if (length(s)) step <- as.numeric(sub("^# Step size = ", "", s[1L]))
    j <- grep("inverse mass matrix:$", ln)
    if (length(j)) {
      k <- j[1L] + 1L
      rows <- character(0)
      while (k <= length(ln) && grepl("^# [-0-9.eE+naif, ]+$", ln[k]) && !grepl("Elapsed", ln[k])) { rows <- c(rows, ln[k]); k <- k + 1L }
      mat <- do.call(rbind, lapply(rows, function(r) as.numeric(strsplit(sub("^# ", "", r), ",\\s*")[[1L]])))
      im <- if (grepl("Diagonal", ln[j[1L]])) as.numeric(mat) else mat
    }
  }
  tm <- function(lab) {
    x <- grep(sprintf("seconds \\(%s\\)", lab), ln, value = TRUE)
    if (length(x)) as.numeric(sub(".*?([0-9.eE+-]+) seconds.*", "\\1", x[1L], perl = TRUE)) else NA_real_
  }
  list(method = method, config = config, values = vals, n_warmup = n_warmup, step_size = step, inv_metric = im,
       time = c(warmup = tm("Warm-up"), sampling = tm("Sampling")))
}
