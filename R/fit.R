# Fit objects returned by the inference functions.

make_fit <- function(model, res, cfg, method, seed, elapsed) {
  nm <- bayescomp_param_names(model)
  chains <- length(res)
  niter <- nrow(res[[1L]]$gq)
  dvars <- colnames(res[[1L]]$diag)
  draws <- array(NA_real_, c(niter, chains, 1L + length(nm)),
                 dimnames = list(iteration = NULL, chain = NULL, variable = c("lp__", nm)))
  sdiag <- array(NA_real_, c(niter, chains, length(dvars) - 1L),
                 dimnames = list(iteration = NULL, chain = NULL, variable = dvars[-1L]))
  for (c in seq_len(chains)) {
    draws[, c, 1L] <- res[[c]]$diag[, "lp__"]
    if (length(nm)) draws[, c, -1L] <- res[[c]]$gq
    sdiag[, c, ] <- res[[c]]$diag[, -1L, drop = FALSE]
  }
  nw <- res[[1L]]$n_warmup_saved
  structure(list(
    draws = draws, sampler_diagnostics = sdiag, method = method, model = model,
    n_warmup_saved = nw, step_size = vapply(res, `[[`, 0, "step_size"),
    inv_metric = lapply(res, `[[`, "inv_metric"),
    unconstrained = lapply(res, `[[`, "q"),
    time = list(chains = do.call(rbind, lapply(res, `[[`, "time")), total = elapsed),
    config = cfg, seed = if (is.null(seed)) attr(res, "seed") else seed,
    diagnostics = if (!is.null(res[[1L]]$diagnostics)) lapply(res, `[[`, "diagnostics"),
    profile = profile_of(res)),
    class = c("bayescomp_fit"))
}

# methods whose draws are Markov chains, for which R-hat and ESS apply
mcmc_methods <- c("sample", "parallel_tempering", "ensemble", "gibbs", "rmhmc", "mclmc", "pdmp",
                  "sgmcmc", "pmcmc", "rjmcmc")

kept <- function(fit, inc_warmup = FALSE) {
  nw <- fit$n_warmup_saved %||% 0L
  n <- dim(fit$draws)[1L]
  if (inc_warmup || nw == 0L) seq_len(n) else if (nw >= n) integer(0) else (nw + 1L):n
}

#' Draws from a fit
#'
#' @param x,fit A `bayescomp_fit`.
#' @param variables Variable names: base names such as `"theta"` select all
#'   elements. `NULL` selects everything.
#' @param inc_warmup Include saved warmup draws.
#' @param format `"array"` (iterations by chains by variables), `"matrix"`
#'   (draws by variables) or `"df"` (data frame with `.chain`, `.iteration`,
#'   `.draw`).
#' @param ... Unused.
#' @return Draws in the requested format.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100,
#'                         seed = 1, refresh = 0)
#' dim(bayescomp_draws(fit))
#' head(bayescomp_draws(fit, format = "df"))
#' head(as.matrix(fit))
#' @export
bayescomp_draws <- function(fit, variables = NULL, inc_warmup = FALSE, format = c("array", "matrix", "df")) {
  format <- match.arg(format)
  a <- fit$draws[kept(fit, inc_warmup), , select_vars(dimnames(fit$draws)[[3L]], variables), drop = FALSE]
  switch(format,
    array = a,
    matrix = draws_to_matrix(a),
    df = {
      m <- draws_to_matrix(a)
      d <- as.data.frame(m, check.names = FALSE)
      d$.chain <- rep(seq_len(dim(a)[2L]), each = dim(a)[1L])
      d$.iteration <- rep(seq_len(dim(a)[1L]), dim(a)[2L])
      d$.draw <- seq_len(nrow(d))
      d
    })
}

draws_to_matrix <- function(a) {
  m <- matrix(a, dim(a)[1L] * dim(a)[2L], dim(a)[3L])
  colnames(m) <- dimnames(a)[[3L]]
  m
}

select_vars <- function(all, variables) {
  if (is.null(variables)) return(all)
  base <- sub("\\[.*$", "", all)
  sel <- all %in% variables | base %in% variables
  miss <- setdiff(variables, c(all, base))
  if (length(miss)) stop(sprintf("unknown variable%s: %s", if (length(miss) > 1L) "s" else "", paste(miss, collapse = ", ")), call. = FALSE)
  all[sel]
}

#' Draws of one variable with its declared shape
#'
#' @param fit A `bayescomp_fit`.
#' @param variable Base name of a parameter, transformed parameter or
#'   generated quantity.
#' @return An array whose first dimension indexes draws (chains stacked) and
#'   whose remaining dimensions are the variable's declared dimensions.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100,
#'                         seed = 1, refresh = 0)
#' mu <- bayescomp_extract(fit, "mu")
#' c(mean(mu), sd(mu))
#' @export
bayescomp_extract <- function(fit, variable) {
  m <- bayescomp_draws(fit, variable, format = "matrix")
  d <- fit$model$decls[[variable]]
  dims <- if (is.null(d)) NULL else c(d$adims, d$eshape)
  if (!length(dims)) return(as.vector(m))
  array(m, c(nrow(m), dims))
}

#' Posterior summary
#'
#' @param object A `bayescomp_fit`.
#' @param variables Variables to summarize (base names allowed).
#' @param probs Quantiles to report.
#' @param ... Unused.
#' @return A data frame with mean, median, sd, mad, quantiles, R-hat and bulk
#'   and tail effective sample sizes (the last three for MCMC fits only).
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100,
#'                         seed = 1, refresh = 0)
#' summary(fit)
#' @export
summary.bayescomp_fit <- function(object, variables = NULL, probs = c(0.05, 0.95), ...) {
  a <- bayescomp_draws(object, variables)
  vars <- dimnames(a)[[3L]]
  mcmc <- object$method %in% mcmc_methods
  rows <- lapply(vars, function(v) {
    x <- matrix(a[, , v], dim(a)[1L], dim(a)[2L])
    xv <- as.vector(x)
    fin <- xv[is.finite(xv)]
    q <- if (length(fin)) stats::quantile(fin, probs, names = FALSE) else rep(NA_real_, length(probs))
    r <- c(mean = mean(xv), median = stats::median(xv), sd = stats::sd(xv), mad = stats::mad(xv), q)
    if (mcmc) {
      # walkers of an ensemble are dependent: R-hat compares the walker means
      # of the independent ensembles instead
      xr <- if (identical(object$method, "ensemble") && !is.null(object$config$ensembles)) {
        W <- object$config$walkers
        vapply(seq_len(object$config$ensembles), function(e) rowMeans(x[, (e - 1L) * W + seq_len(W), drop = FALSE]), numeric(nrow(x)))
      } else x
      r <- c(r, rhat = bayescomp_rhat(matrix(xr, nrow(x))), ess_bulk = bayescomp_ess_bulk(x), ess_tail = bayescomp_ess_tail(x))
    }
    r
  })
  out <- as.data.frame(do.call(rbind, rows))
  names(out)[4L + seq_along(probs)] <- paste0("q", formatC(100 * probs, format = "g"))
  out <- cbind(variable = vars, out, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  class(out) <- c("bayescomp_summary", "data.frame")
  out
}

#' @export
print.bayescomp_summary <- function(x, digits = 3, ...) {
  y <- x
  class(y) <- "data.frame"
  num <- vapply(y, is.numeric, TRUE)
  for (j in which(num)) {
    nmj <- names(y)[j]
    y[[j]] <- if (nmj %in% c("ess_bulk", "ess_tail")) formatC(round(y[[j]]), format = "d", big.mark = "")
              else if (nmj == "rhat") formatC(y[[j]], digits = 3, format = "f")
              else formatC(signif(y[[j]], digits), format = "g", flag = " ")
  }
  print(y, row.names = FALSE, right = TRUE)
  invisible(x)
}

#' @export
print.bayescomp_fit <- function(x, variables = NULL, max_rows = 10L, ...) {
  a <- dim(x$draws)
  cat(sprintf("<bayescomp_fit> method = %s, %d chain%s, %d draws per chain\n", x$method, a[2L],
              if (a[2L] > 1L) "s" else "", length(kept(x))))
  s <- summary(x, variables)
  if (is.null(variables) && nrow(s) > max_rows) {
    print(s[seq_len(max_rows), ])
    cat(sprintf(" # ... with %d more variables\n", nrow(s) - max_rows))
  } else {
    print(s)
  }
  if (x$method == "sample") {
    d <- bayescomp_diagnose(x, quiet = TRUE)
    if (d$num_divergent > 0L) cat(sprintf("Warning: %d divergent transitions after warmup.\n", d$num_divergent))
    if (d$num_max_treedepth > 0L) cat(sprintf("Warning: %d transitions hit the maximum tree depth.\n", d$num_max_treedepth))
    if (any(d$ebfmi < 0.3, na.rm = TRUE)) cat("Warning: low E-BFMI in chains", paste(which(d$ebfmi < 0.3), collapse = ", "), "\n")
  }
  invisible(x)
}

#' Sampler diagnostics
#'
#' Checks what CmdStan's `diagnose` utility checks: transitions that
#' saturated the maximum tree depth, divergent transitions, the energy
#' Bayesian fraction of missing information (E-BFMI, over all chains as
#' CmdStan computes it, and per chain), effective sample sizes below 0.001
#' per draw, and rank normalized split R-hat above 1.01, and prints
#' CmdStan's report.
#'
#' @param fit A `bayescomp_fit` from [bayescomp_sample()], or the names of
#'   Stan CSV files (for example from CmdStan).
#' @param quiet Do not print the report.
#' @return A list with `num_divergent`, `num_max_treedepth`, `ebfmi` (per
#'   chain), `ebfmi_all`, `low_ess` and `high_rhat` (parameter names), and
#'   `problems`
#'   (`TRUE` when any check failed), invisibly.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100,
#'                         seed = 1, refresh = 0)
#' bayescomp_diagnose(fit)
#' @export
bayescomp_diagnose <- function(fit, quiet = FALSE) {
  if (is.character(fit)) fit <- bayescomp_read_csv(fit)
  sd <- fit$sampler_diagnostics[kept(fit), , , drop = FALSE]
  nt <- dim(sd)[1L] * dim(sd)[2L]
  div <- sum(sd[, , "divergent__"], na.rm = TRUE)
  mtd <- fit$config$max_treedepth %||% 10L
  hit <- if (identical(fit$config$algorithm, "nuts")) sum(sd[, , "treedepth__"] >= mtd, na.rm = TRUE) else 0L
  ebfmi <- apply(matrix(sd[, , "energy__"], dim(sd)[1L]), 2L, function(e) {
    if (length(e) < 2L || anyNA(e)) return(NA_real_)
    sum(diff(e)^2) / length(e) / stats::var(e)
  })
  e_all <- as.vector(sd[, , "energy__"])
  ebfmi_all <- if (length(e_all) > 1L && !anyNA(e_all)) mean(diff(e_all)^2) / stats::var(e_all) else NA_real_
  a <- bayescomp_draws(fit)
  vars <- setdiff(dimnames(a)[[3L]], "lp__")
  vars <- vars[!grepl("__$", vars)]
  chains <- dim(a)[2L]
  low <- character(0); high <- character(0)
  for (v in vars) {
    # CmdStan's definitions; undefined (NaN) values fail no check
    d <- cmdstan_mcmc_summary(matrix(a[, , v], dim(a)[1L]))
    n_eff <- if (isTRUE(d[["ess_bulk"]] < d[["ess_tail"]])) d[["ess_bulk"]] else d[["ess_tail"]]
    if (isTRUE(n_eff / length(a[, , v]) < 0.001)) low <- c(low, v)
    if (isTRUE(d[["rhat"]] > 1.01)) high <- c(high, v)
  }
  bad_e <- isTRUE(ebfmi_all < 0.3)
  r <- list(num_divergent = as.integer(div), num_max_treedepth = as.integer(hit), ebfmi = ebfmi, ebfmi_all = ebfmi_all,
            low_ess = low, high_rhat = high,
            problems = div > 0 || hit > 0 || bad_e || length(low) > 0 || length(high) > 0)
  if (!quiet) {
    pct <- function(k) sprintf("%d of %d (%s%%)", as.integer(k), as.integer(nt), format(100 * k / nt, digits = 6))
    cat("Checking sampler transitions treedepth.\n")
    if (hit > 0) cat(sprintf("%s transitions hit the maximum treedepth limit of %d, or 2^%d leapfrog steps.\nTrajectories that are prematurely terminated due to this limit will result in slow exploration.\nFor optimal performance, increase this limit.\n\n", pct(hit), mtd, mtd))
    else cat("Treedepth satisfactory for all transitions.\n\n")
    cat("Checking sampler transitions for divergences.\n")
    if (div > 0) cat(sprintf("%s transitions ended with a divergence.\nThese divergent transitions indicate that HMC is not fully able to explore the posterior distribution.\nTry increasing adapt delta closer to 1.\nIf this doesn't remove all divergences, try to reparameterize the model.\n\n", pct(div)))
    else cat("No divergent transitions found.\n\n")
    cat("Checking E-BFMI - sampler transitions HMC potential energy.\n")
    if (bad_e) {
      cat(sprintf("The E-BFMI, %s, is below the nominal threshold of 0.3 which suggests that HMC may have trouble exploring the target distribution.\nIf possible, try to reparameterize the model.\n\n", format(ebfmi_all, digits = 6)))
    } else cat("E-BFMI satisfactory.\n\n")
    if (length(low)) cat(sprintf("The following parameters had fewer than 0.001 effective draws per transition:\n  %s\nSuch low values indicate that the effective sample size estimators may be biased high and actual performance may be substantially lower than quoted.\n\n", paste(low, collapse = ", ")))
    else cat("Rank-normalized split effective sample size satisfactory for all parameters.\n\n")
    if (length(high)) cat(sprintf("The following parameters had rank-normalized split R-hat greater than 1.01:\n  %s\nSuch high values indicate incomplete mixing and biased estimation.\nYou should consider regularizing your model with additional prior information or a more effective parameterization.\n\n", paste(high, collapse = ", ")))
    else cat("Rank-normalized split R-hat values satisfactory for all parameters.\n\n")
    cat(if (r$problems) "Processing complete.\n" else "Processing complete, no problems detected.\n")
  }
  invisible(r)
}

# CmdStan 2.40's definitions (stan/analyze/mcmc): MCSE of the mean from the
# ESS of unsplit chains, tail ESS from x <= q05 and x >= q95 with one used
# when the other is undefined, and NaN when a chain is constant or all chains
# start at the same value (is_finite_and_varies)
cmdstan_varies <- function(m) {
  if (nrow(m) < 2L || !all(is.finite(m))) return(FALSE)
  flat <- function(v) diff(range(v)) <= 1e-12 * max(abs(v))
  if (any(apply(m, 2L, flat))) return(FALSE)
  !(ncol(m) > 1L && flat(m[1L, ]))
}
cmdstan_ess <- function(m) {
  if (nrow(m) < 4L || na_diag(m)) NaN else ess_raw(m)
}
cmdstan_mcmc_summary <- function(m) {
  out <- c(mcse = NaN, ess_bulk = NaN, ess_tail = NaN, rhat = NaN)
  if (nrow(m) >= 4L && cmdstan_varies(m)) out[["mcse"]] <- stats::sd(as.vector(m)) / sqrt(cmdstan_ess(m))
  s <- split_chains(m)
  if (!cmdstan_varies(s)) return(out)
  out[["rhat"]] <- max(rhat_raw(z_scale(s)), rhat_raw(z_scale(abs(s - stats::quantile(s, 0.5, names = FALSE)))))
  if (nrow(s) < 4L) return(out)
  out[["ess_bulk"]] <- ess_raw(z_scale(s))
  e <- c(cmdstan_ess((s <= stats::quantile(s, 0.05, names = FALSE)) + 0),
         cmdstan_ess((s >= stats::quantile(s, 0.95, names = FALSE)) + 0))
  out[["ess_tail"]] <- if (all(is.nan(e))) NaN else min(e, na.rm = TRUE)
  out
}

#' Posterior summary in the style of CmdStan's stansummary
#'
#' Prints and returns the table of CmdStan's `stansummary` utility: mean,
#' Monte Carlo standard error, standard deviation, median absolute
#' deviation, percentiles, bulk and tail effective sample sizes, bulk ESS per
#' second of sampling and R-hat, for `lp__`, the sampler diagnostics and
#' every variable.
#'
#' @param x A `bayescomp_fit` or the names of Stan CSV files.
#' @param percentiles Percentiles to report (stansummary's `-p`).
#' @param sig_figs Significant figures (`-s`).
#' @param include_param Variables to include (`-i`); base names select all
#'   elements. `NULL` includes all.
#' @param autocorr A chain number: also print the autocorrelations of that
#'   chain's draws (`-a`).
#' @param csv_file File to write the table to (`-c`).
#' @return The table as a data frame, invisibly.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100,
#'                         seed = 1, refresh = 0)
#' bayescomp_stansummary(fit)
#' @export
bayescomp_stansummary <- function(x, percentiles = c(5, 50, 95), sig_figs = 2L, include_param = NULL,
                                  autocorr = NULL, csv_file = NULL) {
  fit <- if (is.character(x)) bayescomp_read_csv(x) else x
  if (!inherits(fit, "bayescomp_fit")) stop("'x' must be a bayescomp_fit or Stan CSV file names", call. = FALSE)
  if (any(percentiles < 0 | percentiles > 100)) stop("'percentiles' must be in [0, 100]", call. = FALSE)
  a <- bayescomp_draws(fit, include_param)
  sd <- if (!is.null(fit$sampler_diagnostics)) fit$sampler_diagnostics[kept(fit), , , drop = FALSE]
  all_vars <- c(if ("lp__" %in% dimnames(a)[[3L]]) "lp__", if (!is.null(sd) && is.null(include_param)) dimnames(sd)[[3L]],
                setdiff(dimnames(a)[[3L]], "lp__"))
  get <- function(v) if (v %in% dimnames(a)[[3L]]) matrix(a[, , v], dim(a)[1L]) else matrix(sd[, , v], dim(sd)[1L])
  mcmc <- fit$method %in% mcmc_methods
  samp_time <- if (!is.null(fit$time$chains)) sum(fit$time$chains[, "sampling"], na.rm = TRUE) else NA_real_
  rows <- lapply(all_vars, function(v) {
    m <- get(v)
    xv <- as.vector(m)
    q <- stats::quantile(xv, percentiles / 100, names = FALSE, na.rm = TRUE)
    d <- if (mcmc) cmdstan_mcmc_summary(m) else c(mcse = NaN, ess_bulk = NaN, ess_tail = NaN, rhat = NaN)
    c(Mean = mean(xv), MCSE = d[["mcse"]], StdDev = stats::sd(xv), MAD = stats::mad(xv), q,
      ESS_bulk = d[["ess_bulk"]], ESS_tail = d[["ess_tail"]], `ESS_bulk/s` = d[["ess_bulk"]] / samp_time,
      R_hat = d[["rhat"]])
  })
  tab <- as.data.frame(do.call(rbind, rows), check.names = FALSE)
  names(tab)[4L + seq_along(percentiles)] <- paste0(formatC(percentiles, format = "g"), "%")
  rownames(tab) <- all_vars
  chains <- dim(a)[2L]
  cat(sprintf("Inference for Stan model: %s\n", fit$model$name %||% "bayescompR_model"))
  cat(sprintf("%d chains: each with iter=%d; warmup=%d; thin=%d; %d iterations saved.\n\n", chains,
              as.integer(fit$config$iter_sampling %||% dim(a)[1L]), as.integer(fit$config$iter_warmup %||% 0L),
              as.integer(fit$config$thin %||% 1L), dim(a)[1L]))
  if (!is.null(fit$time$chains)) {
    tw <- fit$time$chains[, "warmup"]; ts <- fit$time$chains[, "sampling"]
    cat(sprintf("Warmup took (%s) seconds, %s seconds total\n", paste(signif(tw, 2), collapse = ", "), signif(sum(tw), 2)))
    cat(sprintf("Sampling took (%s) seconds, %s seconds total\n\n", paste(signif(ts, 2), collapse = ", "), signif(sum(ts), 2)))
  }
  shown <- tab
  for (j in seq_along(shown)) shown[[j]] <- ifelse(is.nan(tab[[j]]), "nan", formatC(signif(tab[[j]], sig_figs), digits = sig_figs, format = "g"))
  for (j in c("ESS_bulk", "ESS_tail")) shown[[j]] <- ifelse(is.finite(tab[[j]]), formatC(round(tab[[j]]), format = "d"), "nan")
  print(shown, right = TRUE, quote = FALSE)
  cat("\nSamples were drawn using", fit$config$algorithm %||% fit$method, "\n")
  cat("MCSE and ESS were computed with rank normalization; R_hat is the rank normalized split R-hat.\n")
  if (!is.null(autocorr)) {
    ch <- as.integer(autocorr)
    if (ch < 1L || ch > chains) stop("'autocorr' must be a chain number", call. = FALSE)
    vs <- setdiff(dimnames(a)[[3L]], "lp__")
    ac <- vapply(vs, function(v) stats::acf(a[, ch, v], lag.max = min(20L, dim(a)[1L] - 1L), plot = FALSE, na.action = stats::na.pass)$acf[, 1, 1],
                 numeric(min(20L, dim(a)[1L] - 1L) + 1L))
    cat(sprintf("\nDisplaying the autocorrelations for chain %d:\n\n", ch))
    print(round(cbind(Lag = seq_len(nrow(ac)) - 1L, ac), sig_figs + 1L), row.names = FALSE)
  }
  if (!is.null(csv_file)) utils::write.csv(cbind(name = all_vars, tab), csv_file, row.names = FALSE, quote = TRUE)
  invisible(tab)
}

#' @rdname bayescomp_draws
#' @export
as.matrix.bayescomp_fit <- function(x, ...) bayescomp_draws(x, format = "matrix")

#' @rdname bayescomp_draws
#' @param row.names,optional Unused.
#' @export
as.data.frame.bayescomp_fit <- function(x, row.names = NULL, optional = FALSE, ...) bayescomp_draws(x, format = "df")

#' @rdname bayescomp_draws
#' @export
as.array.bayescomp_fit <- function(x, ...) bayescomp_draws(x, format = "array")

#' Convert a fit to a posterior draws object
#'
#' Requires the posterior package.
#' @param fit A `bayescomp_fit`.
#' @param inc_warmup Include saved warmup draws.
#' @return A `draws_array`.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100,
#'                         seed = 1, refresh = 0)
#' d <- bayescomp_as_draws(fit)
#' str(d)
#' @export
bayescomp_as_draws <- function(fit, inc_warmup = FALSE) {
  if (!requireNamespace("posterior", quietly = TRUE)) stop("the posterior package is required", call. = FALSE)
  posterior::as_draws_array(bayescomp_draws(fit, inc_warmup = inc_warmup))
}

#' Diagnostic plots
#'
#' `"trace"` and `"density"` show the draws per chain. `"rank"` shows rank
#' histograms: the pooled draws are ranked and each chain's ranks are
#' binned, so well mixed chains give flat histograms (Vehtari et al. 2021).
#' `"ess_local"` shows the effective sample size of the indicators of small
#' quantile intervals and `"ess_quantile"` the effective sample size of
#' quantile estimates, both against the probability; dashed lines mark the
#' bulk ESS, tail ESS and the ESS of the folded draws |x - median|, which
#' measures the efficiency for scale.
#'
#' @param x A `bayescomp_fit`.
#' @param variables Variables to plot (default: up to four parameters).
#' @param type `"trace"`, `"density"`, `"rank"`, `"ess_local"` or
#'   `"ess_quantile"`.
#' @param bins Number of bins for rank histograms and of intervals for the
#'   local ESS.
#' @param ... Passed to the plotting functions.
#' @return `x`, invisibly.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100,
#'                         seed = 1, refresh = 0)
#' plot(fit)
#' plot(fit, type = "rank")
#' @export
plot.bayescomp_fit <- function(x, variables = NULL, type = c("trace", "density", "rank", "ess_local", "ess_quantile"),
                               bins = 20L, ...) {
  type <- match.arg(type)
  if (is.null(variables)) variables <- utils::head(setdiff(dimnames(x$draws)[[3L]], c("lp__", "log_p__", "log_g__")), 4L)
  a <- bayescomp_draws(x, variables)
  vars <- dimnames(a)[[3L]]
  nc <- dim(a)[2L]
  if (type == "rank") {
    op <- graphics::par(mfrow = c(length(vars), nc), mar = c(2, 2, 1.5, 0.5))
  } else {
    op <- graphics::par(mfrow = c(length(vars), 1L), mar = c(3, 4, 1.5, 1))
  }
  on.exit(graphics::par(op))
  cols <- grDevices::hcl.colors(max(2L, nc), "Dark 3")
  for (v in vars) {
    xm <- matrix(a[, , v], dim(a)[1L])
    if (type == "trace") {
      graphics::matplot(xm, type = "l", lty = 1, col = cols, ylab = v, xlab = "", main = "", ...)
    } else if (type == "density") {
      ds <- apply(xm, 2L, function(col) stats::density(col))
      yl <- range(unlist(lapply(ds, `[[`, "y")))
      xl <- range(unlist(lapply(ds, `[[`, "x")))
      plot(ds[[1L]], col = cols[1L], main = "", xlab = v, xlim = xl, ylim = yl, ...)
      for (k in seq_along(ds)[-1L]) graphics::lines(ds[[k]], col = cols[k])
    } else if (type == "rank") {
      r <- matrix(rank(xm, ties.method = "average"), nrow(xm))
      br <- seq(0, length(r), length.out = bins + 1L)
      cnt <- apply(r, 2L, function(rc) graphics::hist(rc, breaks = br, plot = FALSE)$counts)
      yl <- c(0, max(cnt) * 1.1)
      for (k in seq_len(nc)) {
        graphics::barplot(cnt[, k], space = 0, col = cols[k], border = NA, ylim = yl,
                          main = sprintf("%s, chain %d", v, k), cex.main = 0.8, ...)
        graphics::abline(h = nrow(xm) / bins, lty = 2)
      }
    } else {
      ref <- c(bulk = bayescomp_ess_bulk(xm), tail = bayescomp_ess_tail(xm),
               folded = bayescomp_ess_basic(abs(xm - stats::median(xm))))
      if (type == "ess_local") {
        pr <- seq(0, 1, length.out = bins + 1L)
        qs <- stats::quantile(xm, pr, names = FALSE)
        e <- vapply(seq_len(bins), function(k) {
          ind <- (xm > qs[k] & xm <= qs[k + 1L]) + 0
          if (k == 1L) ind <- (xm <= qs[2L]) + 0
          bayescomp_ess_basic(ind)
        }, 0)
        px <- pr[-1L] - 0.5 / bins
      } else {
        px <- seq(0.025, 0.975, by = 0.025)
        e <- bayescomp_ess_quantile(xm, px)
      }
      yl <- c(0, max(c(e, ref), na.rm = TRUE) * 1.1)
      plot(px, e, type = "b", pch = 19, cex = 0.6, ylim = yl, xlim = c(0, 1), xlab = "probability",
           ylab = sprintf("ESS (%s)", v), main = "", ...)
      graphics::abline(h = ref, lty = 2, col = c("grey30", "firebrick", "steelblue"))
      graphics::abline(h = 400, lty = 3, col = "grey60")
      graphics::legend("bottomright", c(sprintf("%s %.0f", names(ref), ref), "400"), lty = c(2, 2, 2, 3),
                       col = c("grey30", "firebrick", "steelblue", "grey60"), bty = "n", cex = 0.7)
    }
  }
  invisible(x)
}

# profile tables of the chains, with a chain column
profile_of <- function(res) {
  p <- lapply(seq_along(res), function(i) if (!is.null(res[[i]]$profile)) cbind(chain = i, res[[i]]$profile))
  p <- Filter(Negate(is.null), p)
  if (length(p)) do.call(rbind, p) else NULL
}
