# Optimization (posterior mode or maximum likelihood) and the Laplace
# approximation around the mode.

# CmdStan reports log densities with all constants (log_prob<false, true>)
# in ADVI's log_p__ and Pathfinder's lp__. Dropped constants do not depend on
# the parameters, so one evaluation gives the shift for every draw.
propto_shift <- function(model, qs, lp) {
  i <- which(is.finite(lp))[1L]
  if (is.na(i)) return(0)
  full <- attr(safe_lp(model, jacobian = TRUE, propto = FALSE), "value")(qs[i, ])
  if (is.finite(full)) full - lp[i] else 0
}

safe_lp <- function(model, jacobian, propto = TRUE) {
  lpg <- model_lpg(model, jacobian, propto)
  lp <- model_lp(model, jacobian, propto)
  f <- function(u) {
    r <- tryCatch(lpg(u), bayescomp_fatal_error = function(e) stop(e),
       bayescomp_domain_error = function(e) NULL,
       bayescomp_reject_error = function(e) NULL)
    if (is.null(r) || !is.finite(r$value)) return(list(value = -Inf, grad = rep(NA_real_, length(u))))
    r
  }
  # value only version, for draws that need no gradient
  attr(f, "value") <- function(u) {
    v <- tryCatch(lp(u), bayescomp_fatal_error = function(e) stop(e),
       bayescomp_domain_error = function(e) -Inf,
       bayescomp_reject_error = function(e) -Inf)
    if (is.finite(v)) v else -Inf
  }
  f
}

# Hessian of the log density by central differences of the gradient
num_hessian <- function(lpg, u, h = 1e-5) {
  n <- length(u)
  H <- matrix(0, n, n)
  for (i in seq_len(n)) {
    e <- numeric(n); e[i] <- h
    H[, i] <- (lpg(u + e)$grad - lpg(u - e)$grad) / (2 * h)
  }
  0.5 * (H + t(H))
}

#' Find a posterior mode or maximum likelihood estimate
#'
#' @param model A `bayescomp_model`.
#' @param algorithm `"lbfgs"`, `"bfgs"` or `"newton"`.
#' @param jacobian Include the Jacobian adjustment (`TRUE` gives the mode of
#'   the posterior on the unconstrained scale, as for Laplace sampling;
#'   `FALSE`, the default, gives the penalized maximum likelihood estimate).
#' @param init Initial values, as in [bayescomp_sample()].
#' @param seed Random seed for initialization.
#' @param iter Maximum number of iterations.
#' @param tol_obj,tol_rel_obj,tol_grad,tol_rel_grad,tol_param Stopping rules of
#'   Stan's optimizer (L-BFGS and BFGS): the absolute change in the objective,
#'   its relative change (in units of machine epsilon), the gradient norm,
#'   the relative gradient magnitude \eqn{g^T H^{-1} g / \max(|f|, 1)} (in
#'   units of machine epsilon) and the change in the parameters. The first
#'   rule met stops the iterations.
#' @param history_size Number of updates L-BFGS keeps.
#' @param init_alpha Initial step size of the first line search.
#' @param save_iterations Return every iterate as a draw (the initial point
#'   first, the optimum last), as CmdStan's `save_iterations`.
#' @param refresh Report the iteration, log density, step length, gradient
#'   norm and line search step every `refresh` iterations of L-BFGS or BFGS
#'   (0, the default, is quiet).
#' @return A `bayescomp_fit` with one draw holding the optimum, and elements
#'   `par` (named constrained values), `value` (the objective), `code` (0
#'   converged, 1 iteration limit, 2 line search failure), `message` and
#'   `iterations`.
#' @examples
#' m <- bayescomp_model("data { int n; int y; } parameters { real<lower=0,upper=1> p; }
#'   model { y ~ binomial(n, p); }", data = list(n = 10, y = 3))
#' bayescomp_optimize(m)$par
#' @export
bayescomp_optimize <- function(model, algorithm = c("lbfgs", "bfgs", "newton"), jacobian = FALSE,
                               init = 2, seed = NULL, iter = 2000L, tol_obj = 1e-12, tol_rel_obj = 1e4,
                               tol_grad = 1e-8, tol_rel_grad = 1e7, tol_param = 1e-8, history_size = 5L,
                               init_alpha = 1e-3, save_iterations = FALSE, refresh = 0L) {
  check_model(model)
  algorithm <- match.arg(algorithm)
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, 1L, seed, 1L, function(id) {
    u0 <- find_init(model, init, 1L)
    lpg <- safe_lp(model, jacobian)
    r <- optimize_u(lpg, u0, algorithm, list(iter = as_count(iter, "iter", 1L), tol_obj = tol_obj, tol_rel_obj = tol_rel_obj,
                                             tol_grad = tol_grad, tol_rel_grad = tol_rel_grad, tol_param = tol_param,
                                             history_size = as_count(history_size, "history_size", 1L), init_alpha = init_alpha,
                                             save_iterations = isTRUE(save_iterations), refresh = as.integer(refresh)))
    # with save_iterations every iterate is a draw, the optimum last
    r$qs <- if (!is.null(r$path)) r$path[, -1L, drop = FALSE] else matrix(r$u, 1L)
    r$lps <- if (!is.null(r$path)) r$path[, 1L] else r$value
    r$gq <- chain_write_draws(model, r$qs, TRUE)
    r
  })[[1L]]
  k <- nrow(res$qs)
  out <- list(list(q = res$qs, gq = res$gq,
                   diag = cbind(lp__ = res$lps, matrix(0, k, 6L, dimnames = list(NULL, c("accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))),
                   step_size = NA_real_, inv_metric = NULL, n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0),
                   profile = res$profile))
  cfg <- list(algorithm = algorithm, jacobian = jacobian, iter = iter, save_iterations = isTRUE(save_iterations))
  if (algorithm != "newton")
    cfg <- c(cfg, list(init_alpha = init_alpha, tol_obj = tol_obj, tol_rel_obj = tol_rel_obj, tol_grad = tol_grad,
                       tol_rel_grad = tol_rel_grad, tol_param = tol_param, history_size = history_size))
  fit <- make_fit(model, out, cfg, "optimize", seed, proc.time()[["elapsed"]] - t0)
  fit$par <- fit$draws[k, 1L, -1L]
  fit$value <- res$value
  fit$code <- res$code
  fit$message <- res$message %||% ""
  fit$return_code <- optim_return_code(fit$message)
  fit$iterations <- res$iterations
  if (res$code != 0L) warning(sprintf("optimization did not converge: %s", res$message %||% "see code"), call. = FALSE)
  fit
}

optim_defaults <- list(iter = 2000L, tol_obj = 1e-12, tol_rel_obj = 1e4, tol_grad = 1e-8,
                       tol_rel_grad = 1e7, tol_param = 1e-8, history_size = 5L, init_alpha = 1e-3)

optimize_u <- function(lpg, u0, algorithm, ctl = optim_defaults) {
  ctl <- utils::modifyList(optim_defaults, ctl)
  n <- length(u0)
  if (n == 0L) return(list(u = u0, value = lpg(u0)$value, code = 0L, iterations = 0L, message = "no parameters"))
  if (algorithm == "newton") return(newton_opt(lpg, u0, ctl$iter, isTRUE(ctl$save_iterations)))
  qn_opt(lpg, u0, dense = algorithm == "bfgs", ctl)
}

# Quasi Newton minimization of -log density with Stan's stopping rules
# (tol_obj, tol_rel_obj, tol_grad, tol_rel_grad, tol_param; relative
# tolerances are multiples of machine epsilon) and a strong Wolfe line search.
# dense = FALSE keeps the last history_size updates (L-BFGS); TRUE keeps the
# dense inverse Hessian (BFGS).
qn_opt <- function(lpg, x, dense, ctl) {
  eps <- .Machine$double.eps
  fg <- function(u) {
    r <- lpg(u)
    if (!is.finite(r$value) || any(!is.finite(r$grad))) return(list(f = Inf, g = rep(NA_real_, length(u))))
    list(f = -r$value, g = -r$grad)
  }
  cur <- fg(x)
  if (!is.finite(cur$f)) stop("the initial point has a non finite log density or gradient", call. = FALSE)
  n <- length(x)
  S <- list(); Y <- list(); H <- if (dense) diag(n) else NULL
  dir <- function(g) {
    if (dense) return(-as.vector(H %*% g))
    q <- g; k <- length(S); a <- numeric(k)
    for (i in rev(seq_len(k))) { a[i] <- sum(S[[i]] * q) / sum(Y[[i]] * S[[i]]); q <- q - a[i] * Y[[i]] }
    if (k) q <- q * sum(S[[k]] * Y[[k]]) / sum(Y[[k]] * Y[[k]])
    for (i in seq_len(k)) { b <- sum(Y[[i]] * q) / sum(Y[[i]] * S[[i]]); q <- q + S[[i]] * (a[i] - b) }
    -q
  }
  p <- dir(cur$g)
  msg <- "maximum number of iterations reached"; code <- 1L
  path <- if (isTRUE(ctl$save_iterations)) list(c(-cur$f, x))
  for (k in seq_len(ctl$iter)) {
    if (sum(p * cur$g) >= 0) { S <- list(); Y <- list(); if (dense) H <- diag(n); p <- -cur$g }
    a0 <- if (k == 1L) ctl$init_alpha else 1
    ls <- wolfe_search(fg, x, cur, p, a0)
    if (is.null(ls)) {
      if (length(S) || (dense && k > 1L)) { S <- list(); Y <- list(); if (dense) H <- diag(n); p <- -cur$g; ls <- wolfe_search(fg, x, cur, p, ctl$init_alpha) }
      if (is.null(ls)) { msg <- "line search failed to make progress"; code <- 2L; break }
    }
    xn <- x + ls$a * p; nxt <- ls$r
    s <- xn - x; y <- nxt$g - cur$g; sy <- sum(s * y)
    if (sy > 1e-10 * sqrt(sum(s^2) * sum(y^2))) {
      if (dense) {
        if (k == 1L) H <- diag(sy / sum(y * y), n)
        rho <- 1 / sy; Hy <- as.vector(H %*% y)
        H <- H - rho * (outer(s, Hy) + outer(Hy, s)) + (rho^2 * sum(y * Hy) + rho) * outer(s, s)
      } else {
        S[[length(S) + 1L]] <- s; Y[[length(Y) + 1L]] <- y
        if (length(S) > ctl$history_size) { S <- S[-1L]; Y <- Y[-1L] }
      }
    }
    df <- abs(nxt$f - cur$f)
    x <- xn; prev <- cur; cur <- nxt
    if (!is.null(path)) path[[length(path) + 1L]] <- c(-cur$f, x)
    if ((ctl$refresh %||% 0L) > 0L && k %% ctl$refresh == 0L)
      message(sprintf("%6d %14.6g %12.4g %12.4g %10.4g", k, -cur$f, sqrt(sum(s^2)), sqrt(sum(cur$g^2)), ls$a))
    p <- dir(cur$g)
    if (df < ctl$tol_obj) { msg <- "absolute change in objective below tol_obj"; code <- 0L; break }
    if (df / max(abs(cur$f), abs(prev$f), 1) < ctl$tol_rel_obj * eps) { msg <- "relative change in objective below tol_rel_obj"; code <- 0L; break }
    if (sqrt(sum(cur$g^2)) < ctl$tol_grad) { msg <- "gradient norm below tol_grad"; code <- 0L; break }
    if (-sum(cur$g * p) / max(abs(cur$f), 1) < ctl$tol_rel_grad * eps) { msg <- "relative gradient magnitude below tol_rel_grad"; code <- 0L; break }
    if (sqrt(sum(s^2)) < ctl$tol_param) { msg <- "change in parameters below tol_param"; code <- 0L; break }
  }
  list(u = x, value = -cur$f, code = code, iterations = k, message = msg, path = if (!is.null(path)) do.call(rbind, path))
}

# strong Wolfe line search (Nocedal and Wright, algorithms 3.5 and 3.6) for
# minimizing along p; returns NULL when no acceptable step is found
wolfe_search <- function(fg, x, cur, p, a1, c1 = 1e-4, c2 = 0.9) {
  f0 <- cur$f; d0 <- sum(cur$g * p)
  if (!(d0 < 0)) return(NULL)
  ev <- function(a) { r <- fg(x + a * p); list(a = a, r = r, f = r$f, d = if (is.finite(r$f)) sum(r$g * p) else NA_real_) }
  zoom <- function(lo, hi) {
    for (i in seq_len(50L)) {
      # quadratic interpolation from f(lo), f'(lo) and f(hi), kept away from
      # the ends of the bracket; bisection when it fails
      w <- hi$a - lo$a
      a <- if (is.finite(hi$f)) lo$a - lo$d * w^2 / (2 * (hi$f - lo$f - lo$d * w)) else lo$a + 0.5 * w
      lim <- sort(c(lo$a + 0.1 * w, hi$a - 0.1 * w))
      if (!is.finite(a) || a < lim[1L] || a > lim[2L]) a <- lo$a + 0.5 * w
      t <- ev(a)
      if (!is.finite(t$f) || t$f > f0 + c1 * a * d0 || t$f >= lo$f) hi <- t
      else {
        if (abs(t$d) <= -c2 * d0) return(t)
        if (t$d * (hi$a - lo$a) >= 0) hi <- lo
        lo <- t
      }
      if (abs(hi$a - lo$a) < 1e-16 * max(1, abs(lo$a))) break
    }
    if (lo$a > 0 && lo$f < f0) lo else NULL
  }
  prev <- list(a = 0, r = cur, f = f0, d = d0)
  a <- a1
  for (i in seq_len(50L)) {
    t <- ev(a)
    if (!is.finite(t$f) || t$f > f0 + c1 * a * d0 || (i > 1L && t$f >= prev$f)) return(zoom(prev, t))
    if (abs(t$d) <= -c2 * d0) return(t)
    if (t$d >= 0) return(zoom(t, prev))
    prev <- t
    a <- 2 * a
  }
  NULL
}

newton_opt <- function(lpg, u, iter, save_iterations = FALSE) {
  r <- lpg(u)
  path <- if (save_iterations) list(c(r$value, u))
  done_path <- function(res) { if (!is.null(path)) res$path <- do.call(rbind, path); res }
  for (k in seq_len(iter)) {
    H <- num_hessian(lpg, u)
    e <- eigen(H, symmetric = TRUE)
    # make the Hessian negative definite
    ev <- -pmax(abs(e$values), 1e-8)
    step <- -as.vector(e$vectors %*% ((t(e$vectors) %*% r$grad) / ev))
    a <- 1
    repeat {
      un <- u + a * step
      rn <- lpg(un)
      if (is.finite(rn$value) && rn$value >= r$value) break
      a <- a / 2
      # no ascent possible: converged only if the gradient vanishes
      if (a < 1e-10) return(done_path(list(u = u, value = r$value, code = if (max(abs(r$grad)) < 1e-6) 0L else 2L, iterations = k,
                                           message = if (max(abs(r$grad)) < 1e-6) "gradient norm below tol_grad" else "line search failed to make progress")))
    }
    done <- abs(rn$value - r$value) < 1e-10 * (1 + abs(r$value))
    u <- un; r <- rn
    if (!is.null(path)) path[[length(path) + 1L]] <- c(r$value, u)
    if (done) return(done_path(list(u = u, value = r$value, code = 0L, iterations = k, message = "relative change in objective below tol_rel_obj")))
  }
  done_path(list(u = u, value = r$value, code = 1L, iterations = iter, message = "maximum number of iterations reached"))
}

#' Laplace approximation
#'
#' Draws from a normal approximation to the posterior on the unconstrained
#' scale, centered at the mode with covariance equal to the inverse negative
#' Hessian of the log density (as Stan's `laplace` method).
#'
#' @param model A `bayescomp_model`.
#' @param mode Optional result of [bayescomp_optimize()] (with
#'   `jacobian = TRUE`), a named list of constrained values at the mode, or
#'   the name of a file with the mode: a Stan CSV file from an optimization
#'   (by CmdStan or [bayescomp_write_csv()]) or a JSON or R dump file.
#' @param jacobian Use the Jacobian adjusted density (recommended).
#' @param draws Number of draws.
#' @param seed Random seed.
#' @param calculate_lp Evaluate the log density of the draws (`log_p__`);
#'   with `FALSE` it is `NaN`, as in CmdStan.
#' @param refresh Report the stages, the optimizer's iterations and every
#'   `refresh` draws (0, the default, is quiet).
#' @return A `bayescomp_fit` with extra columns `log_p__` and `log_g__`.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_laplace(m, draws = 200, seed = 1)
#' sd(bayescomp_extract(fit, "mu"))
#' @export
bayescomp_laplace <- function(model, mode = NULL, jacobian = TRUE, draws = 1000L, seed = NULL, calculate_lp = TRUE,
                              refresh = 0L) {
  check_model(model)
  refresh <- as.integer(refresh)
  say <- function(...) if (refresh > 0L) message(sprintf(...))
  draws <- as_count(draws, "draws", 1L)
  if (is.character(mode))
    mode <- if (grepl("\\.csv$", mode, ignore.case = TRUE)) bayescomp_read_csv(mode, model = model) else bayescomp_read_data(mode)
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, 1L, seed, 1L, function(id) {
    lpg <- safe_lp(model, jacobian)
    u <- if (is.null(mode)) {
      say("Finding the mode with L-BFGS")
      op <- optimize_u(lpg, find_init(model, 2, 1L), "lbfgs", list(refresh = refresh))
      if (op$code != 0L) warning(sprintf("the optimizer did not converge (code %d); the Laplace approximation may be centered away from the mode", op$code), call. = FALSE)
      op$u
    } else if (inherits(mode, "bayescomp_fit")) {
      u1 <- mode$unconstrained[[1L]]
      as.vector(u1[nrow(u1), ])
    } else {
      bayescomp_unconstrain(model, mode)
    }
    say("Calculating the Hessian")
    H <- num_hessian(lpg, u)
    S <- tryCatch(solve(-H), error = function(e) stop("the Hessian at the mode is singular", call. = FALSE))
    S <- 0.5 * (S + t(S))
    L <- tryCatch(t(chol(S)), error = function(e) stop("the negative Hessian at the mode is not positive definite", call. = FALSE))
    n <- length(u)
    z <- matrix(stats::rnorm(draws * n), draws, n)
    qs <- sweep(z %*% t(L), 2L, u, "+")
    # CmdStan's log_q__: the unnormalized Gaussian density, 0 at the mode
    log_g <- -0.5 * rowSums(z^2)
    say("Generating draws")
    log_p <- if (isTRUE(calculate_lp)) vapply(seq_len(draws), function(i) {
      if (refresh > 0L && i %% refresh == 0L) say("iteration: %d / %d", i, draws)
      attr(lpg, "value")(qs[i, ])
    }, 0) else rep(NaN, draws)
    gq <- chain_write_draws(model, qs, FALSE)
    dg <- matrix(0, draws, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- 0
    # Pareto k of the importance ratios posterior / approximation: above
    # 0.7 the approximation is too far from the posterior to rely on, as is
    # common for regressions with many coefficients (van Leeuwen and van Erp
    # 2025)
    lr <- log_p - log_g
    k <- if (isTRUE(calculate_lp) && draws >= 100L && all(is.finite(lr))) psis_smooth(lr)$k else NA_real_
    list(q = qs, gq = cbind(log_p__ = log_p, log_g__ = log_g, gq), diag = dg, step_size = NA_real_,
         inv_metric = S, n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0), pareto_k = k)
  })
  fit <- make_fit_extra(model, res, list(jacobian = jacobian), "laplace", seed, proc.time()[["elapsed"]] - t0)
  fit$pareto_k <- res[[1L]]$pareto_k
  if (isTRUE(fit$pareto_k > 0.7))
    warning(sprintf("the Laplace approximation is poor here (Pareto k = %.2f of the importance ratios); check it with bayescomp_check_approximation() or sample with bayescomp_sample()", fit$pareto_k), call. = FALSE)
  fit
}

# like make_fit but allows extra leading columns in gq (log_p__, log_g__)
make_fit_extra <- function(model, res, cfg, method, seed, elapsed) {
  r <- res[[1L]]
  nm <- colnames(r$gq)
  n <- nrow(r$gq)
  draws <- array(NA_real_, c(n, 1L, 1L + length(nm)), dimnames = list(iteration = NULL, chain = NULL, variable = c("lp__", nm)))
  draws[, 1L, 1L] <- r$diag[, "lp__"]
  draws[, 1L, -1L] <- r$gq
  structure(list(draws = draws, sampler_diagnostics = NULL, method = method, model = model,
                 n_warmup_saved = 0L, step_size = NA_real_, inv_metric = list(r$inv_metric),
                 unconstrained = list(r$q), time = list(total = elapsed), config = cfg, seed = seed %||% attr(res, "seed"),
                 profile = profile_of(res)),
            class = "bayescomp_fit")
}
