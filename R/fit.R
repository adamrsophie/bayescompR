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
    config = cfg, seed = if (is.null(seed)) attr(res, "seed") else seed),
    class = c("bayescomp_fit"))
}

# methods whose draws are Markov chains, for which R-hat and ESS apply
mcmc_methods <- c("sample", "parallel_tempering", "ensemble", "gibbs", "rmhmc", "mclmc", "pdmp",
                  "sgmcmc", "pmcmc", "rjmcmc")

kept <- function(fit, inc_warmup = FALSE) {
  nw <- fit$n_warmup_saved %||% 0L
  n <- dim(fit$draws)[1L]
  if (inc_warmup || nw == 0L) seq_len(n) else (nw + 1L):n
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
      r <- c(r, rhat = bayescomp_rhat(x), ess_bulk = bayescomp_ess_bulk(x), ess_tail = bayescomp_ess_tail(x))
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
#' Counts divergent transitions and transitions that saturated the maximum
#' tree depth, and computes the energy Bayesian fraction of missing
#' information (E-BFMI) per chain.
#'
#' @param fit A `bayescomp_fit` from [bayescomp_sample()].
#' @param quiet Do not print messages.
#' @return A list with `num_divergent`, `num_max_treedepth` and `ebfmi`.
#' @export
bayescomp_diagnose <- function(fit, quiet = FALSE) {
  sd <- fit$sampler_diagnostics[kept(fit), , , drop = FALSE]
  div <- sum(sd[, , "divergent__"], na.rm = TRUE)
  mtd <- fit$config$max_treedepth %||% 10L
  hit <- if (identical(fit$config$algorithm, "nuts")) sum(sd[, , "treedepth__"] >= mtd, na.rm = TRUE) else 0L
  ebfmi <- apply(matrix(sd[, , "energy__"], dim(sd)[1L]), 2L, function(e) {
    if (length(e) < 2L || anyNA(e)) return(NA_real_)
    sum(diff(e)^2) / length(e) / stats::var(e)
  })
  r <- list(num_divergent = as.integer(div), num_max_treedepth = as.integer(hit), ebfmi = ebfmi)
  if (!quiet) {
    cat(sprintf("Divergent transitions: %d\n", r$num_divergent))
    cat(sprintf("Transitions at maximum tree depth: %d\n", r$num_max_treedepth))
    cat(sprintf("E-BFMI per chain: %s\n", paste(formatC(ebfmi, digits = 3, format = "f"), collapse = ", ")))
  }
  invisible(r)
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
