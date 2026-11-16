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
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_variational(m, output_samples = 200, seed = 1)
#' mean(bayescomp_extract(fit, "mu"))
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
    # the ELBO and the output draws need values only, not gradients
    lv <- attr(lpg, "value")
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
      lp <- apply(th, 1L, lv)
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
    # as CmdStan: log_g__ unnormalized in the standard normal draws, log_p__
    # with all constants
    log_g <- -0.5 * rowSums(z^2)
    log_p <- apply(qs, 1L, lv)
    log_p <- log_p + propto_shift(model, qs, log_p)
    gq <- chain_write_draws(model, qs, FALSE)
    trace <<- list(elbo = r$trace, eta = eta, iterations = r$iter, converged = r$converged,
                   status = if (r$converged) "converged: relative ELBO change below tol_rel_obj" else "stopped: maximum iterations reached",
                   mean = vp$mu,
                   # CmdStan's first output row, drawn in this run's RNG stream
                   mean_draw = chain_write_draws(model, matrix(vp$mu, 1L), FALSE),
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

# L-BFGS inverse Hessian approximation in compact form (Byrd, Nocedal and
# Schnabel 1994), H = diag(alpha) + B G t(B) with B = [diag(alpha) Y, S],
# and the normal approximation it defines, drawn and evaluated in O(n m)
# per draw as in Pathfinder (Zhang et al. 2022, Algorithm 4); no n x n
# matrix is formed.
lbfgs_normal <- function(x, g, alpha, S, Y) {
  m <- ncol(S)
  SY <- crossprod(S, Y)
  Rm <- SY
  Rm[lower.tri(Rm)] <- 0
  Rinv <- backsolve(Rm, diag(m))
  B <- cbind(alpha * Y, S)
  G <- rbind(cbind(matrix(0, m, m), -Rinv),
             cbind(-t(Rinv), t(Rinv) %*% (diag(diag(SY), m) + crossprod(Y, alpha * Y)) %*% Rinv))
  mu <- x + alpha * g + as.vector(B %*% (G %*% crossprod(B, g)))
  q <- qr(B / sqrt(alpha))
  Q <- qr.Q(q)
  R <- qr.R(q)[, order(q$pivot), drop = FALSE]
  M <- diag(nrow(R)) + R %*% G %*% t(R)
  Lt <- t(chol(0.5 * (M + t(M))))
  if (!all(is.finite(mu)) || !all(is.finite(Lt))) stop("non finite approximation")
  sa <- sqrt(alpha)
  list(mu = mu, Q = Q, Lt = Lt, sa = sa,
       logdet = sum(log(sa)) + sum(log(diag(Lt))),
       # draws theta = mu + sqrt(alpha) (u + Q (Lt - I) t(Q) u), u ~ N(0, I)
       draw = function(U) sweep(sweep(U + (U %*% Q) %*% t(Lt - diag(ncol(Q))) %*% t(Q), 2L, sa, "*"), 2L, mu, "+"))
}

# the covariance of an lbfgs_normal as a dense matrix
lbfgs_normal_cov <- function(a) {
  A <- diag(length(a$sa)) + a$Q %*% (a$Lt - diag(ncol(a$Q))) %*% t(a$Q)
  A <- a$sa * A
  tcrossprod(A)
}

# the same inverse Hessian as a dense matrix (reference for tests)
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

pathfinder_path <- function(lpg, u0, ctl) {
  history <- ctl$history_size; num_elbo_draws <- ctl$num_elbo_draws; num_draws <- ctl$num_draws
  eps <- .Machine$double.eps
  n <- length(u0)
  lv <- attr(lpg, "value") %||% function(q) lpg(q)$value
  x <- u0
  r <- lpg(x)
  if (!is.finite(r$value)) stop("Pathfinder: the log density is not finite at the initial point", call. = FALSE)
  S <- matrix(0, n, 0L); Y <- matrix(0, n, 0L)
  alpha <- rep(1, n)
  best <- list(elbo = -Inf)
  trace <- numeric(0)
  evaluate <- function(x, g, alpha, S, Y) {
    a <- tryCatch(lbfgs_normal(x, g, alpha, S, Y), error = function(e) NULL)
    if (is.null(a)) return(NULL)
    z <- matrix(stats::rnorm(num_elbo_draws * n), num_elbo_draws, n)
    th <- a$draw(z)
    lq <- -0.5 * rowSums(z^2) - 0.5 * n * log(2 * pi) - a$logdet
    lp <- apply(th, 1L, lv)
    e <- mean(ifelse(is.finite(lp), lp, -Inf) - lq)
    list(elbo = if (is.nan(e)) -Inf else e, mu = a$mu, approx = a)
  }
  for (it in seq_len(ctl$max_iter)) {
    # L-BFGS direction on -log p
    d <- -lbfgs_two_loop(-r$grad, alpha, S, Y)
    step <- if (it == 1L) ctl$init_alpha else 1
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
    df <- abs(rn$value - f0)
    rel <- df / max(abs(rn$value), abs(f0), 1)
    rel_grad <- abs(sum(rn$grad * lbfgs_two_loop(rn$grad, alpha, S, Y))) / max(abs(rn$value), 1)
    r <- rn
    ev <- NULL
    if (ncol(S) > 0L) {
      ev <- evaluate(x, r$grad, alpha, S, Y)
      if (!is.null(ev)) {
        trace <- c(trace, ev$elbo)
        if (ev$elbo > best$elbo) best <- ev
      }
    }
    if ((ctl$refresh %||% 0L) > 0L && it %% ctl$refresh == 0L)
      message(sprintf("Path [%d] :Iter %6d  log prob %14.6g  ||dx|| %10.4g  ||grad|| %10.4g  ELBO %12.6g",
                      ctl$path_id %||% 1L, it, r$value, sqrt(sum(s^2)), sqrt(sum(r$grad^2)), ev$elbo %||% NA_real_))
    if (df < ctl$tol_obj || rel < ctl$tol_rel_obj * eps || sqrt(sum(r$grad^2)) < ctl$tol_grad ||
        rel_grad < ctl$tol_rel_grad * eps || sqrt(sum(s^2)) < ctl$tol_param) break
  }
  if (!is.finite(best$elbo)) stop("Pathfinder: no usable normal approximation was found along the optimization path", call. = FALSE)
  z <- matrix(stats::rnorm(num_draws * n), num_draws, n)
  th <- best$approx$draw(z)
  lq <- -0.5 * rowSums(z^2) - 0.5 * n * log(2 * pi) - best$approx$logdet
  lp <- if (isTRUE(ctl$calculate_lp)) apply(th, 1L, lv) else rep(NaN, num_draws)
  list(q = th, log_p = lp, log_q = lq, elbo = best$elbo, mu = best$mu, approx = best$approx, trace = trace)
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
#' smoothed importance resampling, as CmdStan's `pathfinder` method does.
#'
#' @param model A `bayescomp_model`.
#' @param num_paths Number of paths.
#' @param num_draws Number of draws from each path's approximation.
#' @param num_psis_draws Number of draws returned after importance
#'   resampling.
#' @param history_size L-BFGS history size.
#' @param max_lbfgs_iters Maximum L-BFGS iterations per path.
#' @param num_elbo_draws Draws per ELBO estimate.
#' @param psis_resample Resample `num_psis_draws` draws with Pareto
#'   smoothed importance weights. With `FALSE`, all `num_paths * num_draws`
#'   draws of the paths are returned, as in Stan; so they are with one path.
#' @param calculate_lp Evaluate the log density of the returned draws. With
#'   `FALSE`, `log_p__` is `NaN` and no resampling is done, as in Stan.
#' @param save_single_paths Keep each path's draws (element `single_paths`,
#'   a list of fits, which [bayescomp_write_csv()] writes as
#'   `<file>_path_<i>.csv`) and ELBO trajectory.
#' @param init_alpha First line search step size.
#' @param tol_obj,tol_rel_obj,tol_grad,tol_rel_grad,tol_param L-BFGS
#'   convergence tolerances, as in [bayescomp_optimize()] (the relative
#'   tolerances are multiples of the machine epsilon).
#' @param init,seed As in [bayescomp_sample()].
#' @param cores Paths to run in parallel; defaults to `getOption("mc.cores", 1L)`.
#' @param refresh Report each path's L-BFGS iteration, log density, step
#'   length, gradient norm and ELBO every `refresh` iterations, and its best
#'   ELBO at the end (0, the default, is quiet).
#' @return A `bayescomp_fit` with method `"pathfinder"`, extra columns
#'   `log_p__` and `log_q__`, and elements `pareto_k`, `elbo` (the best ELBO
#'   of each path) and `path` (the path of each draw).
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_pathfinder(m, num_paths = 2, num_draws = 200, num_psis_draws = 200, seed = 1)
#' mean(bayescomp_extract(fit, "mu"))
#' @export
bayescomp_pathfinder <- function(model, num_paths = 4L, num_draws = 1000L, num_psis_draws = 1000L, history_size = 5L,
                                 max_lbfgs_iters = 1000L, num_elbo_draws = 25L, psis_resample = TRUE, calculate_lp = TRUE,
                                 save_single_paths = FALSE, init_alpha = 1e-3, tol_obj = 1e-12, tol_rel_obj = 1e4,
                                 tol_grad = 1e-8, tol_rel_grad = 1e7, tol_param = 1e-8,
                                 init = 2, seed = NULL, cores = getOption("mc.cores", 1L), refresh = 0L) {
  check_model(model)
  if (model$dim == 0L) stop("the model has no parameters", call. = FALSE)
  refresh <- as.integer(refresh)
  num_paths <- as_count(num_paths, "num_paths", 1L)
  ctl <- pathfinder_control(history_size = as_count(history_size, "history_size", 1L),
                            max_iter = as_count(max_lbfgs_iters, "max_lbfgs_iters", 1L),
                            num_elbo_draws = as_count(num_elbo_draws, "num_elbo_draws", 1L),
                            num_draws = as_count(num_draws, "num_draws", 1L), calculate_lp = isTRUE(calculate_lp),
                            init_alpha = init_alpha, tol_obj = tol_obj, tol_rel_obj = tol_rel_obj, tol_grad = tol_grad,
                            tol_rel_grad = tol_rel_grad, tol_param = tol_param, refresh = refresh)
  num_psis_draws <- as_count(num_psis_draws, "num_psis_draws", 1L)
  t0 <- proc.time()[["elapsed"]]
  paths <- run_chains(model, num_paths, seed, cores, function(id) {
    lpg <- safe_lp(model, jacobian = TRUE)
    p <- pathfinder_path(lpg, find_init(model, init, id), utils::modifyList(ctl, list(path_id = id)))
    if (refresh > 0L) message(sprintf("Path [%d] :Best ELBO %.6g", id, p$elbo))
    p
  })
  q <- do.call(rbind, lapply(paths, `[[`, "q"))
  lp <- unlist(lapply(paths, `[[`, "log_p"))
  # CmdStan's lp__ has all constants
  shift <- propto_shift(model, q, lp)
  lp <- lp + shift
  for (i in seq_along(paths)) paths[[i]]$log_p <- paths[[i]]$log_p + shift
  lq <- unlist(lapply(paths, `[[`, "log_q"))
  pid <- rep(seq_len(num_paths), vapply(paths, function(p) nrow(p$q), 0L))
  resample <- isTRUE(psis_resample) && ctl$calculate_lp && num_paths > 1L
  k <- NA_real_
  pick <- NULL
  res <- run_chains(model, 1L, next_seed(seed), 1L, function(id) {
    sel <- seq_len(nrow(q))
    if (resample) {
      lw <- ifelse(is.finite(lp), lp - lq, -Inf)
      ps <- psis_smooth(lw)
      k <<- ps$k
      w <- exp(ps$lw - max(ps$lw))
      sel <- sample.int(nrow(q), num_psis_draws, replace = TRUE, prob = w)
    }
    pick <<- sel
    pathfinder_result(model, q[sel, , drop = FALSE], lp[sel], lq[sel])
  })
  fit <- make_fit_extra(model, res, list(num_paths = num_paths, num_draws = ctl$num_draws, num_psis_draws = num_psis_draws,
                                         psis_resample = resample, calculate_lp = ctl$calculate_lp, history_size = ctl$history_size,
                                         init_alpha = init_alpha, tol_obj = tol_obj, tol_rel_obj = tol_rel_obj, tol_grad = tol_grad,
                                         tol_rel_grad = tol_rel_grad, tol_param = tol_param, max_lbfgs_iters = ctl$max_iter),
                        "pathfinder", seed, proc.time()[["elapsed"]] - t0)
  fit$pareto_k <- k
  fit$elbo <- vapply(paths, `[[`, 0, "elbo")
  fit$path <- pid[pick]
  if (isTRUE(save_single_paths)) {
    fit$single_paths <- lapply(seq_along(paths), function(i) {
      p <- paths[[i]]
      r <- run_chains(model, 1L, NULL, 1L, function(id) pathfinder_result(model, p$q, p$log_p, p$log_q))
      f <- make_fit_extra(model, r, fit$config, "pathfinder", seed, NA_real_)
      f$path <- rep(i, nrow(p$q))
      f$elbo_trace <- p$trace
      f
    })
  }
  if (is.finite(k) && k > 0.7) warning(sprintf("Pareto k = %.2f: the importance weights are unreliable", k), call. = FALSE)
  fit
}

pathfinder_result <- function(model, qs, lp, lq) {
  gq <- chain_write_draws(model, qs, FALSE)
  dg <- matrix(0, nrow(qs), 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
  list(q = qs, gq = cbind(log_p__ = lp, log_q__ = lq, gq), diag = dg, step_size = NA_real_,
       inv_metric = NULL, n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
}

pathfinder_control <- function(...) {
  utils::modifyList(list(history_size = 5L, max_iter = 1000L, num_elbo_draws = 25L, num_draws = 1L, calculate_lp = TRUE,
                         init_alpha = 1e-3, tol_obj = 1e-12, tol_rel_obj = 1e4, tol_grad = 1e-8, tol_rel_grad = 1e7,
                         tol_param = 1e-8), list(...))
}

# the best normal approximation along one Pathfinder path
pathfinder_best_normal <- function(lpg, u0, ctl = pathfinder_control()) {
  r <- pathfinder_path(lpg, u0, ctl)
  list(mu = r$mu, L = t(chol(lbfgs_normal_cov(r$approx))))
}

