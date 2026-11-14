# Normalizing flow variational inference (Rezende and Mohamed 2015): the
# approximation is a diagonal Gaussian pushed through a stack of planar
# flows f(z) = z + u tanh(w'z + b), with u reparameterized so that every
# layer is invertible (u'w >= -1). The ELBO gradient comes from reverse mode
# AD through the flow combined with the model gradient at the transformed
# draws. The fitted density is exact, so importance weights and Pareto k
# diagnose the fit as for Pathfinder.

#' Normalizing flow variational inference
#'
#' @param model A `bayescomp_model`.
#' @param layers Number of planar flow layers.
#' @param iter Adam iterations.
#' @param grad_samples Monte Carlo draws per gradient.
#' @param learning_rate Adam step size.
#' @param output_samples Number of approximate posterior draws returned.
#' @param psis_resample Resample the output with Pareto smoothed importance
#'   weights.
#' @param init,seed As in [bayescomp_sample()].
#' @param refresh Print the ELBO every `refresh` iterations (0 for none).
#' @return A `bayescomp_fit` with method `"flow"`, `elbo` (trace), and
#'   `pareto_k`.
#' @references Rezende, D. and Mohamed, S. (2015). Variational inference
#'   with normalizing flows. *Proceedings of the 32nd International
#'   Conference on Machine Learning*, 1530 to 1538.
#' @export
bayescomp_flow_vi <- function(model, layers = 6L, iter = 2000L, grad_samples = 8L, learning_rate = 0.02,
                              output_samples = 1000L, psis_resample = TRUE, init = 2, seed = NULL, refresh = 0L) {
  check_model(model)
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  K <- as_count(layers, "layers")
  iter <- as_count(iter, "iter", 1L)
  M <- as_count(grad_samples, "grad_samples", 1L)
  S <- as_count(output_samples, "output_samples", 1L)
  t0 <- proc.time()[["elapsed"]]
  info <- NULL
  res <- run_chains(model, 1L, seed, 1L, function(id) {
    lpg <- safe_lp(model, jacobian = TRUE, propto = FALSE)
    # parameter layout: mu (n), omega (n), then per layer u_hat (n), w (n), b (1)
    np <- 2L * n + K * (2L * n + 1L)
    th <- numeric(np)
    th[seq_len(n)] <- find_init(model, init, id)
    lay <- function(k) {
      o <- 2L * n + (k - 1L) * (2L * n + 1L)
      list(u = o + seq_len(n), w = o + n + seq_len(n), b = o + 2L * n + 1L)
    }
    for (k in seq_len(K)) th[lay(k)$w] <- stats::rnorm(n, 0, 0.1)
    # flow forward pass; works on numbers and on AD values alike
    flow <- function(th, eps) {
      mu <- ad_select(th, seq_len(n)); om <- ad_select(th, n + seq_len(n))
      z <- mu + exp(om) * eps
      ld <- sum(om)
      for (k in seq_len(K)) {
        ix <- lay(k)
        uh <- ad_select(th, ix$u); w <- ad_select(th, ix$w); b <- ad_select(th, ix$b)
        wu <- sum(w * uh)
        u <- uh + (ad_log1p_exp(wu) - 1 - wu) * w / (sum(w * w) + 1e-10)
        h <- tanh(sum(w * z) + b)
        ld <- ld + log(abs(1 + (1 - h * h) * sum(w * u)))
        z <- z + u * h
      }
      list(z = z, ld = ld)
    }
    # ELBO surrogate: its gradient is the reparameterization gradient
    elbo_grad <- function(th) {
      E <- matrix(stats::rnorm(M * n), M, n)
      gs <- vector("list", M); vals <- numeric(M)
      for (m in seq_len(M)) {
        f <- flow(th, E[m, ])
        r <- lpg(as.numeric(f$z))
        if (!is.finite(r$value)) return(NULL)
        vals[m] <- r$value + f$ld + 0.5 * n * (1 + log(2 * pi))
        gs[[m]] <- r$grad
      }
      out <- ad_value_grad(function(t) {
        acc <- 0
        for (m in seq_len(M)) { f <- flow(t, E[m, ]); acc <- acc + sum(gs[[m]] * f$z) + f$ld }
        acc / M
      }, th)
      list(elbo = mean(vals), grad = out$grad)
    }
    m1 <- numeric(np); m2 <- numeric(np)
    trace <- numeric(0)
    best <- th; best_e <- -Inf; run_e <- NA_real_
    for (it in seq_len(iter)) {
      r <- elbo_grad(th)
      if (is.null(r) || any(!is.finite(r$grad))) {
        th <- best
        m1[] <- 0; m2[] <- 0
        next
      }
      run_e <- if (is.na(run_e)) r$elbo else 0.95 * run_e + 0.05 * r$elbo
      trace <- c(trace, r$elbo)
      if (it > 50L && run_e > best_e) { best_e <- run_e; best <- th }
      m1 <- 0.9 * m1 + 0.1 * r$grad
      m2 <- 0.999 * m2 + 0.001 * r$grad^2
      lr <- learning_rate * min(1, it / 100)
      th <- th + lr * (m1 / (1 - 0.9^it)) / (sqrt(m2 / (1 - 0.999^it)) + 1e-8)
      if (refresh > 0L && it %% refresh == 0L) message(sprintf("iter %5d ELBO %.3f", it, run_e))
    }
    # the iterate with the best running ELBO, not the last one
    if (is.finite(best_e)) th <- best
    E <- matrix(stats::rnorm(S * n), S, n)
    qs <- matrix(NA_real_, S, n); log_g <- numeric(S)
    for (s in seq_len(S)) {
      f <- flow(th, E[s, ])
      qs[s, ] <- f$z
      log_g[s] <- -0.5 * sum(E[s, ]^2) - 0.5 * n * log(2 * pi) - f$ld
    }
    log_p <- apply(qs, 1L, function(q) lpg(q)$value)
    lw <- log_p - log_g
    lw[!is.finite(lw)] <- -Inf
    ps <- psis_smooth(lw)
    if (psis_resample) {
      w <- exp(ps$lw - max(ps$lw))
      idx <- sample.int(S, S, replace = TRUE, prob = w)
      qs <- qs[idx, , drop = FALSE]; log_p <- log_p[idx]; log_g <- log_g[idx]
    }
    info <<- list(elbo = trace, pareto_k = ps$k, theta = th)
    dg <- matrix(0, S, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- log_p
    list(q = qs, gq = cbind(log_p__ = log_p, log_g__ = log_g, chain_write_draws(model, qs, FALSE)), diag = dg,
         step_size = learning_rate, inv_metric = NULL, n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  fit <- make_fit_extra(model, res, list(layers = K), "flow", seed, proc.time()[["elapsed"]] - t0)
  fit$elbo <- info$elbo
  fit$pareto_k <- info$pareto_k
  if (is.finite(info$pareto_k) && info$pareto_k > 0.7)
    warning(sprintf("Pareto k = %.2f > 0.7: the flow approximation is poor", info$pareto_k), call. = FALSE)
  fit
}
