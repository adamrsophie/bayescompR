# Variational inference: automatic differentiation variational inference
# (Kucukelbir et al. 2017, as implemented in Stan) and Pathfinder (Zhang,
# Carpenter, Gelman and Vehtari 2022) with Pareto smoothed importance
# resampling.

#' Automatic differentiation variational inference
#'
#' Fits a normal approximation on the unconstrained space by stochastic
#' gradient ascent on the evidence lower bound (ELBO), with Stan's adaptive
#' step size sequence and convergence checks on the relative ELBO change.
#'
#' @param model A `bayescomp_model`.
#' @param algorithm `"meanfield"` (diagonal) or `"fullrank"`.
#' @param iter Maximum iterations.
#' @param grad_samples Monte Carlo draws per gradient.
#' @param elbo_samples Monte Carlo draws per ELBO estimate.
#' @param eta Step size scale; `NULL` selects it adaptively.
#' @param adapt_iter Iterations per candidate `eta` during adaptation.
#' @param tol_rel_obj Convergence tolerance on the relative ELBO change.
#' @param eval_elbo Evaluate the ELBO every `eval_elbo` iterations.
#' @param output_samples Number of approximate posterior draws.
#' @param init,seed As in [bayescomp_sample()].
#' @param refresh Report progress every `refresh` ELBO evaluations (0 quiet).
#' @return A `bayescomp_fit` with method `"variational"`, extra columns
#'   `log_p__` and `log_g__`, and element `vi`: the ELBO trace (`elbo`), the
#'   step size `eta`, `iterations`, `converged` (whether the relative ELBO
#'   change fell below `tol_rel_obj`) and a `status` message.
#' @export
bayescomp_variational <- function(model, algorithm = c("meanfield", "fullrank"), iter = 10000L,
                                  grad_samples = 1L, elbo_samples = 100L, eta = NULL, adapt_iter = 50L,
                                  tol_rel_obj = 0.01, eval_elbo = 100L, output_samples = 1000L,
                                  init = 2, seed = NULL, refresh = 0L) {
  check_model(model)
  algorithm <- match.arg(algorithm)
  t0 <- proc.time()[["elapsed"]]
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  trace <- NULL
  res <- run_chains(model, 1L, seed, 1L, function(id) {
    lpg <- safe_lp(model, jacobian = TRUE)
    mu0 <- find_init(model, init, 1L)
    full <- algorithm == "fullrank"
    vp0 <- list(mu = mu0, omega = rep(0, n), L = diag(n))
    draw_z <- function(k) matrix(stats::rnorm(k * n), k, n)
    transform <- function(vp, z) {
      if (full) sweep(z %*% t(vp$L), 2L, vp$mu, "+") else sweep(sweep(z, 2L, exp(vp$omega), "*"), 2L, vp$mu, "+")
    }
    entropy <- function(vp) {
      0.5 * n * (1 + log(2 * pi)) + if (full) sum(log(abs(diag(vp$L)))) else sum(vp$omega)
    }
    elbo <- function(vp) {
      th <- transform(vp, draw_z(elbo_samples))
      lp <- apply(th, 1L, function(q) lpg(q)$value)
      # as Stan's calc_ELBO: draws whose log density is not finite are
      # dropped, and the ELBO fails only when every draw is dropped
      if (!any(is.finite(lp))) return(-Inf)
      mean(lp[is.finite(lp)]) + entropy(vp)
    }
    grad <- function(vp) {
      z <- draw_z(grad_samples)
      th <- transform(vp, z)
      gs <- matrix(0, grad_samples, n)
      for (i in seq_len(grad_samples)) {
        r <- lpg(th[i, ])
        if (!is.finite(r$value) || any(!is.finite(r$grad)))
          stop("the gradient of the log density is not finite during ADVI; the model may be ill conditioned", call. = FALSE)
        gs[i, ] <- r$grad
      }
      gmu <- colMeans(gs)
      if (full) {
        gL <- crossprod(gs, z) / grad_samples
        gL[upper.tri(gL)] <- 0
        diag(gL) <- diag(gL) + 1 / diag(vp$L)
        list(mu = gmu, L = gL)
      } else {
        list(mu = gmu, omega = colMeans(gs * z) * exp(vp$omega) + 1)
      }
    }
    sgd <- function(vp, eta, iters, check = FALSE) {
      s <- NULL
      cb <- numeric(0)
      cb_size <- max(as.integer(0.1 * iter / eval_elbo), 2L)
      elbo_prev <- NA_real_
      tr <- numeric(0)
      for (k in seq_len(iters)) {
        g <- grad(vp)
        gv <- unlist(g)
        s <- if (is.null(s)) gv^2 else 0.1 * gv^2 + 0.9 * s
        step <- eta * k^(-0.5 + 1e-16) / (1 + sqrt(s))
        flatp <- if (full) c(vp$mu, vp$L) else c(vp$mu, vp$omega)
        flatp <- flatp + step * gv
        vp$mu <- flatp[seq_len(n)]
        if (full) { vp$L <- matrix(flatp[-seq_len(n)], n, n); vp$L[upper.tri(vp$L)] <- 0 }
        else vp$omega <- flatp[-seq_len(n)]
        if (any(!is.finite(flatp))) return(list(vp = vp, ok = FALSE, trace = tr))
        if (check && k %% eval_elbo == 0L) {
          e <- elbo(vp)
          tr <- c(tr, e)
          if (refresh > 0L && (length(tr) %% refresh == 0L)) message(sprintf("iter %5d ELBO %.3f", k, e))
          if (!is.na(elbo_prev)) {
            rel <- abs((e - elbo_prev) / e)
            cb <- utils::tail(c(cb, rel), cb_size)
            if (mean(cb) < tol_rel_obj || stats::median(cb) < tol_rel_obj)
              return(list(vp = vp, ok = TRUE, converged = TRUE, trace = tr, iter = k))
          }
          elbo_prev <- e
        }
      }
      list(vp = vp, ok = TRUE, converged = FALSE, trace = tr, iter = iters)
    }
    if (is.null(eta)) {
      best <- -Inf
      for (cand in c(100, 10, 1, 0.1, 0.01)) {
        r <- tryCatch(sgd(vp0, cand, adapt_iter), error = function(e) list(ok = FALSE))
        e <- if (isTRUE(r$ok)) elbo(r$vp) else -Inf
        if (is.finite(e) && e > best) { best <- e; eta <- cand }
      }
      if (is.null(eta)) stop("ADVI: all step size candidates diverged; try supplying 'eta'", call. = FALSE)
    }
    r <- sgd(vp0, eta, iter, check = TRUE)
    if (!r$ok) stop("ADVI diverged; try a smaller 'eta'", call. = FALSE)
    if (!r$converged)
      warning(sprintf("ADVI reached the maximum of %d iterations without meeting tol_rel_obj; the approximation may be poor", iter), call. = FALSE)
    vp <- r$vp
    z <- draw_z(output_samples)
    qs <- transform(vp, z)
    log_g <- -0.5 * rowSums(z^2) - 0.5 * n * log(2 * pi) - (entropy(vp) - 0.5 * n * (1 + log(2 * pi)))
    log_p <- apply(qs, 1L, function(q) lpg(q)$value)
    gq <- chain_write_draws(model, qs, FALSE)
    trace <<- list(elbo = r$trace, eta = eta, iterations = r$iter, converged = r$converged,
                   status = if (r$converged) "converged: relative ELBO change below tol_rel_obj" else "stopped: maximum iterations reached",
                   mean = vp$mu,
                   sd = if (full) sqrt(rowSums(vp$L^2)) else exp(vp$omega), L = if (full) vp$L else NULL)
    dg <- matrix(0, output_samples, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    list(q = qs, gq = cbind(log_p__ = log_p, log_g__ = log_g, gq), diag = dg, step_size = eta,
         inv_metric = NULL, n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  fit <- make_fit_extra(model, res, list(algorithm = algorithm), "variational", seed, proc.time()[["elapsed"]] - t0)
  fit$vi <- trace
  fit
}

# ---- Pathfinder -----------------------------------------------------------------------

# L-BFGS inverse Hessian approximation as a dense matrix, with initial
# diagonal `alpha` and curvature pairs (columns of S and Y).
lbfgs_inv_hessian <- function(alpha, S, Y) {
  n <- length(alpha)
  H <- diag(alpha, n)
  for (j in seq_len(ncol(S))) {
    s <- S[, j]; y <- Y[, j]
    rho <- 1 / sum(y * s)
    Hy <- as.vector(H %*% y)
    H <- H - rho * (outer(s, Hy) + outer(Hy, s)) + (rho^2 * sum(y * Hy) + rho) * outer(s, s)
  }
  0.5 * (H + t(H))
}

pathfinder_path <- function(lpg, u0, history, max_iter, num_elbo_draws, num_draws, tol_rel_obj, tol_grad) {
  n <- length(u0)
  x <- u0
  r <- lpg(x)
  if (!is.finite(r$value)) stop("Pathfinder: the log density is not finite at the initial point", call. = FALSE)
  S <- matrix(0, n, 0L); Y <- matrix(0, n, 0L)
  alpha <- rep(1, n)
  best <- list(elbo = -Inf)
  evaluate <- function(x, g, alpha, S, Y) {
    H <- tryCatch(lbfgs_inv_hessian(alpha, S, Y), error = function(e) NULL)
    if (is.null(H)) return(NULL)
    L <- tryCatch(t(chol(H)), error = function(e) NULL)
    if (is.null(L)) return(NULL)
    mu <- x + as.vector(H %*% g)
    z <- matrix(stats::rnorm(num_elbo_draws * n), num_elbo_draws, n)
    th <- sweep(z %*% t(L), 2L, mu, "+")
    lq <- -0.5 * rowSums(z^2) - 0.5 * n * log(2 * pi) - sum(log(diag(L)))
    lp <- apply(th, 1L, function(q) lpg(q)$value)
    e <- mean(ifelse(is.finite(lp), lp, -Inf) - lq)
    list(elbo = if (is.nan(e)) -Inf else e, mu = mu, L = L)
  }
  for (it in seq_len(max_iter)) {
    # L-BFGS direction on -log p
    d <- -lbfgs_two_loop(-r$grad, alpha, S, Y)
    step <- 1
    f0 <- r$value
    gd <- sum(r$grad * d)
    if (!(gd > 0)) { d <- r$grad; gd <- sum(r$grad^2) }
    repeat {
      xn <- x + step * d
      rn <- lpg(xn)
      if (is.finite(rn$value) && rn$value >= f0 + 1e-4 * step * gd) break
      step <- step / 2
      if (step < 1e-12) break
    }
    if (step < 1e-12) break
    s <- xn - x
    y <- r$grad - rn$grad
    if (sum(s * y) > 2.2e-16 * sum(y * y)) {
      # diagonal update of the initial inverse Hessian (Pathfinder, Algorithm 3)
      a <- sum(y * alpha * y); b <- sum(y * s); c <- sum(s * s / alpha)
      an <- 1 / (a / (b * alpha) + y^2 / b - a * s^2 / (b * c * alpha^2))
      if (all(is.finite(an)) && all(an > 0)) alpha <- an
      S <- cbind(S, s); Y <- cbind(Y, y)
      if (ncol(S) > history) { S <- S[, -1L, drop = FALSE]; Y <- Y[, -1L, drop = FALSE] }
    }
    x <- xn
    rel <- abs(rn$value - f0) / max(abs(rn$value), abs(f0), 1)
    r <- rn
    if (ncol(S) > 0L) {
      ev <- evaluate(x, r$grad, alpha, S, Y)
      if (!is.null(ev) && ev$elbo > best$elbo) best <- ev
    }
    if (rel < tol_rel_obj || sqrt(sum(r$grad^2)) < tol_grad) break
  }
  if (!is.finite(best$elbo)) stop("Pathfinder: no usable normal approximation was found along the optimization path", call. = FALSE)
  z <- matrix(stats::rnorm(num_draws * n), num_draws, n)
  th <- sweep(z %*% t(best$L), 2L, best$mu, "+")
  lq <- -0.5 * rowSums(z^2) - 0.5 * n * log(2 * pi) - sum(log(diag(best$L)))
  lp <- apply(th, 1L, function(q) lpg(q)$value)
  list(q = th, log_p = lp, log_q = lq, elbo = best$elbo, mu = best$mu, L = best$L)
}

lbfgs_two_loop <- function(g, alpha, S, Y) {
  m <- ncol(S)
  q <- g
  a <- numeric(m)
  rho <- if (m) 1 / colSums(S * Y) else numeric(0)
  for (j in rev(seq_len(m))) { a[j] <- rho[j] * sum(S[, j] * q); q <- q - a[j] * Y[, j] }
  r <- alpha * q
  for (j in seq_len(m)) { b <- rho[j] * sum(Y[, j] * r); r <- r + S[, j] * (a[j] - b) }
  r
}

# Pareto smoothed importance sampling (Vehtari, Simpson, Gelman, Yao and
# Gabry 2024): returns smoothed log weights and the Pareto k estimate.
psis_smooth <- function(lw) {
  S <- length(lw)
  lw <- lw - max(lw)
  M <- min(ceiling(0.2 * S), ceiling(3 * sqrt(S)))
  ord <- order(lw)
  if (M < 5L || S < 10L) return(list(lw = lw, k = NA_real_))
  tail_idx <- ord[(S - M + 1L):S]
  cutoff <- lw[ord[S - M]]
  x <- exp(lw[tail_idx]) - exp(cutoff)
  fit <- gpd_fit(sort(x))
  if (is.finite(fit$k)) {
    p <- (seq_len(M) - 0.5) / M
    qs <- exp(cutoff) + gpd_quantile(p, fit$k, fit$sigma)
    smoothed <- log(pmin(qs, 1))
    lw[tail_idx[order(lw[tail_idx])]] <- smoothed
  }
  list(lw = lw, k = fit$k)
}

# Zhang and Stephens (2009) generalized Pareto fit with the weakly
# informative prior used in PSIS
gpd_fit <- function(x) {
  n <- length(x)
  m <- 30L + floor(sqrt(n))
  prior <- 3
  xstar <- x[floor(n / 4 + 0.5)]
  theta <- 1 / x[n] + (1 - sqrt(m / (seq_len(m) - 0.5))) / prior / xstar
  lx <- vapply(theta, function(t) { k <- -mean(log1p(-t * x)); log(t / k) + k - 1 }, 0) * n
  w <- 1 / vapply(lx, function(l) sum(exp(lx - l)), 0)
  theta_hat <- sum(theta * w)
  k <- mean(log1p(-theta_hat * x))
  sigma <- -k / theta_hat
  k <- (n * k + 0.5 * 10) / (n + 10)
  list(k = k, sigma = sigma)
}
gpd_quantile <- function(p, k, sigma) sigma * expm1(-k * log1p(-p)) / k

#' Pathfinder variational inference
#'
#' Runs `num_paths` L-BFGS optimization paths, forms a normal approximation
#' from the L-BFGS inverse Hessian estimate at every iterate, keeps the one
#' with the largest ELBO on each path, and combines the paths by Pareto
#' smoothed importance resampling.
#'
#' @param model A `bayescomp_model`.
#' @param num_paths Number of paths.
#' @param num_draws Number of draws returned.
#' @param history_size L-BFGS history size.
#' @param max_lbfgs_iters Maximum L-BFGS iterations per path.
#' @param num_elbo_draws Draws per ELBO estimate.
#' @param psis_resample Resample `num_draws` draws with Pareto smoothed
#'   importance weights. With `FALSE`, all `num_paths * num_draws` draws of
#'   the paths are returned, as in Stan.
#' @param tol_rel_obj,tol_grad L-BFGS convergence tolerances.
#' @param init,seed As in [bayescomp_sample()].
#' @param cores Paths to run in parallel; defaults to `getOption("mc.cores", 1L)`.
#' @return A `bayescomp_fit` with method `"pathfinder"`, extra columns
#'   `log_p__` and `log_q__`, and element `pareto_k`.
#' @export
bayescomp_pathfinder <- function(model, num_paths = 4L, num_draws = 1000L, history_size = 5L,
                                 max_lbfgs_iters = 1000L, num_elbo_draws = 25L, psis_resample = TRUE,
                                 tol_rel_obj = 1e-12, tol_grad = 1e-8, init = 2, seed = NULL, cores = getOption("mc.cores", 1L)) {
  check_model(model)
  if (model$dim == 0L) stop("the model has no parameters", call. = FALSE)
  t0 <- proc.time()[["elapsed"]]
  paths <- run_chains(model, num_paths, seed, cores, function(id) {
    lpg <- safe_lp(model, jacobian = TRUE)
    pathfinder_path(lpg, find_init(model, init, id), history_size, max_lbfgs_iters,
                    num_elbo_draws, num_draws, tol_rel_obj, tol_grad)
  })
  q <- do.call(rbind, lapply(paths, `[[`, "q"))
  lp <- unlist(lapply(paths, `[[`, "log_p"))
  lq <- unlist(lapply(paths, `[[`, "log_q"))
  k <- NA_real_
  res <- run_chains(model, 1L, if (is.null(seed)) NULL else seed + 1L, 1L, function(id) {
    sel <- seq_len(nrow(q))
    if (psis_resample) {
      lw <- ifelse(is.finite(lp), lp - lq, -Inf)
      ps <- psis_smooth(lw)
      k <<- ps$k
      w <- exp(ps$lw - max(ps$lw))
      sel <- sample.int(nrow(q), num_draws, replace = TRUE, prob = w)
    }
    # without resampling every path's draws are returned, as in Stan
    qs <- q[sel, , drop = FALSE]
    gq <- chain_write_draws(model, qs, FALSE)
    dg <- matrix(0, length(sel), 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    list(q = qs, gq = cbind(log_p__ = lp[sel], log_q__ = lq[sel], gq), diag = dg, step_size = NA_real_,
         inv_metric = NULL, n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  fit <- make_fit_extra(model, res, list(num_paths = num_paths), "pathfinder", seed, proc.time()[["elapsed"]] - t0)
  fit$pareto_k <- k
  fit$elbo <- vapply(paths, `[[`, 0, "elbo")
  if (is.finite(k) && k > 0.7) warning(sprintf("Pareto k = %.2f: the importance weights are unreliable", k), call. = FALSE)
  fit
}

# the best normal approximation along one Pathfinder path
pathfinder_best_normal <- function(lpg, u0, history, max_iter, num_elbo_draws, tol_rel_obj, tol_grad) {
  r <- pathfinder_path(lpg, u0, history, max_iter, num_elbo_draws, 1L, tol_rel_obj, tol_grad)
  list(mu = r$mu, L = r$L)
}
