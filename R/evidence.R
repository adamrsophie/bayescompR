# Evidence and model choice: nested sampling (Skilling 2006) and reversible
# jump MCMC across a list of models (Green 1995) with independence
# proposals from Gaussian approximations of each model's posterior.

#' Nested sampling
#'
#' Estimates the marginal likelihood by shrinking a set of live points
#' through nested likelihood contours (Skilling 2006). New points are drawn
#' from the prior restricted to the current contour by hit-and-run slice
#' sampling along random directions whitened by the covariance of the live
#' points (Handley, Hobson and Lasenby 2015; Yallup et al. 2026), which needs
#' no step size and rejects nothing. The prior must be
#' proper; it is the model evaluated on an empty data set, as in
#' [bayescomp_sgmcmc()].
#'
#' @param model A `bayescomp_model` with an observation count in its data.
#' @param num_live Number of live points.
#' @param mcmc_steps Slice sampling steps per replacement; the default is
#'   three times the number of parameters, at least 5.
#' @param dlogz Stop when the evidence left in the live points is below
#'   this on the log scale.
#' @param max_iter Maximum number of replacements.
#' @param n_name,subset As in [bayescomp_sgmcmc()].
#' @param output_samples Posterior draws resampled from the weighted dead
#'   points.
#' @param seed Random seed.
#' @return A `bayescomp_fit` with method `"nested"` and elements
#'   `log_evidence`, `log_evidence_se` and `information`.
#' @references Skilling, J. (2006). Nested sampling for general Bayesian
#'   computation. *Bayesian Analysis* 1, 833 to 859. Handley, W. J., Hobson,
#'   M. P. and Lasenby, A. N. (2015). PolyChord: next-generation nested
#'   sampling. *Monthly Notices of the Royal Astronomical Society* 453, 4384
#'   to 4398. Yallup, D. et al. (2026). Nested slice sampling.
#'   arXiv:2601.23252.
#' @examples
#' \donttest{
#' m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; }
#'   model { mu ~ normal(0, 10); y ~ normal(mu, 1); }",
#'   data = list(N = 20, y = seq(-1, 3, length.out = 20)))
#' fit <- bayescomp_nested(m, num_live = 100, seed = 1)
#' fit$log_evidence
#' }
#' @export
bayescomp_nested <- function(model, num_live = 400L, mcmc_steps = NULL, dlogz = 0.01, max_iter = 50000L,
                             n_name = "N", subset = NULL, output_samples = 1000L, seed = NULL) {
  check_model(model)
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  sp <- split_setup(model, n_name, subset)
  nl <- as_count(num_live, "num_live", 10L)
  mcmc_steps <- if (is.null(mcmc_steps)) max(5L, 3L * n) else as_count(mcmc_steps, "mcmc_steps", 1L)
  S <- as_count(output_samples, "output_samples", 1L)
  t0 <- proc.time()[["elapsed"]]
  pseed <- next_seed(seed)
  pr_fit <- bayescomp_sample(sp$prior, chains = 1L, iter_warmup = 500L, iter_sampling = 4L * nl, seed = pseed, refresh = 0L)
  info <- NULL
  res <- run_chains(model, 1L, seed, 1L, function(id) {
    lprior <- lp_value_fun(sp$prior)
    lfull <- lp_value_fun(model)
    logL <- function(u) { a <- lfull(u); b <- lprior(u); if (is.finite(a) && is.finite(b)) a - b else -Inf }
    P <- pr_fit$unconstrained[[1L]][seq(4L, 4L * nl, by = 4L), , drop = FALSE]
    LL <- apply(P, 1L, logL)
    LP <- apply(P, 1L, lprior)
    dead_u <- list(); dead_l <- numeric(0); dead_L <- numeric(0); converged <- FALSE
    log_z <- -Inf
    log_x <- 0
    for (i in seq_len(max_iter)) {
      w <- which.min(LL)
      Lstar <- LL[w]
      log_x_new <- -i / nl
      log_wt <- Lstar + log_x + log1p(-exp(log_x_new - log_x))
      log_z <- log_sum_exp_v(c(log_z, log_wt))
      dead_u[[i]] <- P[w, ]; dead_l[i] <- log_wt; dead_L[i] <- Lstar
      log_x <- log_x_new
      # replace the worst point with a constrained prior draw
      C <- stats::cov(P)
      Lc <- tryCatch(t(chol(C + diag(1e-12, n))), error = function(e) diag(sqrt(pmax(diag(C), 1e-12)), n))
      k <- sample.int(nl, 1L)
      while (k == w && nl > 1L) k <- sample.int(nl, 1L)
      x <- P[k, ]; lx <- LL[k]; px <- LP[k]
      for (s in seq_len(mcmc_steps)) {
        # slice sampling of the prior along a whitened random direction,
        # restricted to L > Lstar, with a width of three live point scales
        z <- stats::rnorm(n)
        d <- as.vector(Lc %*% (z / sqrt(sum(z^2))))
        lev <- px - stats::rexp(1)
        inside <- function(t) {
          y <- x + t * d; py <- lprior(y)
          if (!is.finite(py) || py <= lev) return(NULL)
          ly <- logL(y)
          if (ly > Lstar) list(y = y, py = py, ly = ly)
        }
        lo <- -3 * stats::runif(1); hi <- lo + 3; j <- 0L
        while (j < 50L && !is.null(inside(lo))) { lo <- lo - 3; j <- j + 1L }
        j <- 0L
        while (j < 50L && !is.null(inside(hi))) { hi <- hi + 3; j <- j + 1L }
        repeat {
          t <- stats::runif(1, lo, hi)
          v <- inside(t)
          if (!is.null(v)) { x <- v$y; px <- v$py; lx <- v$ly; break }
          if (t < 0) lo <- t else hi <- t
          if (hi - lo < 1e-12) break
        }
      }
      P[w, ] <- x; LL[w] <- lx; LP[w] <- px
      if (max(LL) + log_x < log_z + log(expm1(dlogz))) { converged <- TRUE; break }
    }
    # the remaining live points share the final prior volume
    lw_live <- LL + log_x - log(nl)
    log_z_final <- log_sum_exp_v(c(log_z, lw_live))
    U <- rbind(do.call(rbind, dead_u), P)
    lw <- c(dead_l, lw_live) - log_z_final
    # information H = sum p_i log L_i - log Z over dead and live points
    pw <- exp(lw); Li <- c(dead_L[seq_along(dead_l)], LL)
    H <- sum(pw[pw > 0] * Li[pw > 0]) - log_z_final
    idx <- sample.int(nrow(U), S, replace = TRUE, prob = exp(lw))
    qs <- U[idx, , drop = FALSE]
    if (!converged) warning("nested sampling stopped at max_iter before the evidence tolerance was reached", call. = FALSE)
    info <<- list(log_z = log_z_final, se = sqrt(max(H, 0) / nl), H = H, iter = i, converged = converged)
    dg <- matrix(0, S, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- apply(qs, 1L, lfull)
    list(q = qs, gq = chain_write_draws(model, qs, FALSE), diag = dg, step_size = NA_real_, inv_metric = NULL,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  fit <- make_fit_extra(model, res, list(num_live = nl), "nested", seed, proc.time()[["elapsed"]] - t0)
  fit$log_evidence <- info$log_z
  fit$log_evidence_se <- info$se
  fit$information <- info$H
  fit$iterations <- info$iter
  fit$converged <- info$converged
  fit
}

log_sum_exp_v <- function(x) { m <- max(x); if (!is.finite(m)) m else m + log(sum(exp(x - m))) }

#' Reversible jump MCMC across models
#'
#' Samples the joint posterior over a finite set of models and their
#' parameters. Within a model the parameters move by adaptive random walk
#' Metropolis; between models the sampler proposes the other model's
#' parameters from an independent Gaussian approximation of its posterior
#' (Pathfinder, widened by `proposal_scale`), which keeps the reversible
#' jump acceptance ratio simple and exact for any dimensions.
#'
#' @param models A list of `bayescomp_model` objects (named for the
#'   output).
#' @param prior_prob Prior model probabilities (default uniform).
#' @param iter_warmup,iter_sampling Iterations.
#' @param proposal_scale Scale applied to the Gaussian approximations.
#' @param seed Random seed.
#' @return An object of class `bayescomp_rjmcmc` with `probabilities`
#'   (posterior model probabilities and their Monte Carlo standard errors),
#'   `fits` (a `bayescomp_fit` per model with the draws visited there),
#'   `trace` (model index per iteration) and `jump_rate`.
#' @references Green, P. J. (1995). Reversible jump Markov chain Monte Carlo
#'   computation and Bayesian model determination. *Biometrika* 82, 711 to
#'   732.
#' @examples
#' \donttest{
#' y <- c(0.5, 1.2, 2, 0.8, 1.6)
#' mk <- function(s) bayescomp_model(sprintf("data { int N; vector[N] y; } parameters { real mu; }
#'   model { mu ~ normal(0, %g); y ~ normal(mu, 1); }", s), data = list(N = length(y), y = y))
#' r <- bayescomp_rjmcmc(list(narrow = mk(1), wide = mk(10)), iter_sampling = 1000, seed = 1)
#' r$probabilities
#' }
#' @export
bayescomp_rjmcmc <- function(models, prior_prob = NULL, iter_warmup = 500L, iter_sampling = 5000L,
                             proposal_scale = 1.2, seed = NULL) {
  if (inherits(models, "bayescomp_model") || !length(models)) stop("'models' must be a list of models", call. = FALSE)
  for (m in models) check_model(m)
  J <- length(models)
  nms <- names(models) %||% paste0("model", seq_len(J))
  if (any(!nzchar(nms))) nms[!nzchar(nms)] <- paste0("model", which(!nzchar(nms)))
  pp <- if (is.null(prior_prob)) rep(1 / J, J) else prior_prob / sum(prior_prob)
  nw <- as_count(iter_warmup, "iter_warmup"); ns <- as_count(iter_sampling, "iter_sampling", 1L)
  t0 <- proc.time()[["elapsed"]]
  rc <- run_chains(models[[1L]], 1L, seed, 1L, function(id) {
    lps <- lapply(models, lp_value_fun)
    refs <- lapply(models, function(m) gaussian_reference(m, "pathfinder", proposal_scale))
    us <- lapply(refs, function(r) r$mu)
    cur_lp <- vapply(seq_len(J), function(j) lps[[j]](us[[j]]), 0)
    sc <- vapply(models, function(m) 2.38 / sqrt(m$dim), 0) / proposal_scale
    j <- 1L
    trace <- integer(ns)
    draws <- lapply(seq_len(J), function(k) list())
    jumps <- 0; tries <- 0
    for (it in seq_len(nw + ns)) {
      # within model
      for (s in 1:3) {
        y <- us[[j]] + sc[j] * as.vector(refs[[j]]$L %*% stats::rnorm(models[[j]]$dim))
        ly <- lps[[j]](y)
        a <- if (is.finite(ly)) min(1, exp(ly - cur_lp[j])) else 0
        if (stats::runif(1) < a) { us[[j]] <- y; cur_lp[j] <- ly }
        if (it <= nw) sc[j] <- sc[j] * exp((a - 0.3) / sqrt(it))
      }
      # between models
      if (J > 1L) {
        others <- setdiff(seq_len(J), j)
        k <- others[sample.int(length(others), 1L)]
        v <- refs[[k]]$mu + as.vector(refs[[k]]$L %*% stats::rnorm(models[[k]]$dim))
        lv <- lps[[k]](v)
        la <- lv + log(pp[k]) + refs[[j]]$log_g(us[[j]]) - cur_lp[j] - log(pp[j]) - refs[[k]]$log_g(v)
        a <- if (is.finite(la)) min(1, exp(la)) else 0
        if (it > nw) { tries <- tries + 1; jumps <- jumps + a }
        if (stats::runif(1) < a) { us[[k]] <- v; cur_lp[k] <- lv; j <- k }
      }
      if (it > nw) {
        i <- it - nw
        trace[i] <- j
        draws[[j]][[length(draws[[j]]) + 1L]] <- c(cur_lp[j], us[[j]])
      }
    }
    list(trace = trace, draws = draws, jump = jumps / max(tries, 1))
  })
  out <- rc[[1L]]
  # generated quantities draw from a stream derived from the run's seed
  fits <- with_seed(attr(rc, "seed") %% 2147483646L + 1L, lapply(seq_len(J), function(k) {
    d <- out$draws[[k]]
    if (!length(d)) return(NULL)
    M <- do.call(rbind, d)
    q <- M[, -1L, drop = FALSE]
    dg <- matrix(0, nrow(q), 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- M[, 1L]
    r <- list(q = q, gq = cbind(log_p__ = M[, 1L], chain_write_draws(models[[k]], q, FALSE)), diag = dg,
              step_size = NA_real_, inv_metric = NULL, n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
    make_fit_extra(models[[k]], list(r), list(), "rjmcmc", seed, 0)
  }))
  names(fits) <- nms
  p <- tabulate(out$trace, J) / ns
  # batch means standard error of the model indicators
  nb <- max(2L, floor(sqrt(ns)))
  bm <- vapply(seq_len(J), function(k) {
    b <- tapply(out$trace == k, ceiling(seq_len(ns) / (ns / nb)), mean)
    stats::sd(b) / sqrt(length(b))
  }, 0)
  structure(list(probabilities = data.frame(model = nms, prior = pp, posterior = p, mcse = bm),
                 fits = fits, trace = out$trace, jump_rate = out$jump,
                 time = proc.time()[["elapsed"]] - t0),
            class = "bayescomp_rjmcmc")
}

#' @export
print.bayescomp_rjmcmc <- function(x, digits = 3, ...) {
  cat(sprintf("Reversible jump MCMC over %d models (%d draws, jump acceptance %.2f)\n",
              nrow(x$probabilities), length(x$trace), x$jump_rate))
  print(format(x$probabilities, digits = digits), row.names = FALSE)
  invisible(x)
}
