# Runtime of the higher order functions.
#
# The user function receives the call's extra arguments unchanged. When some
# of them (or the initial state) are differentiable, the solvers compute the
# derivatives of their result with respect to them: forward sensitivities for
# ODEs, the implicit function theorem for algebraic equations. The Jacobians
# of the user function come from the AD tape, one backward pass per output.

# The differentiable leaves of a list of arguments: their values, the
# originals, and a function rebuilding the arguments from a (possibly
# differentiable) vector of those values.
hof_flatten <- function(args) {
  parents <- list()
  theta <- numeric(0)
  walk <- function(x) {
    if (inherits(x, "bcad")) {
      off <- length(theta)
      parents[[length(parents) + 1L]] <<- x
      theta <<- c(theta, as.numeric(x$v))
      return(list(ad = TRUE, idx = off + seq_along(x$v), dim = dim(x$v)))
    }
    if (is.list(x)) return(list(ad = FALSE, list = TRUE, items = lapply(unclass(x), walk), cls = oldClass(x)))
    list(ad = FALSE, list = FALSE, v = x)
  }
  specs <- lapply(args, walk)
  rebuild <- function(z) {
    rb <- function(s) {
      if (isTRUE(s$ad)) {
        r <- ad_gather(z, s$idx, s$dim)
        if (!inherits(r, "bcad") && length(s$dim) > 1L) dim(r) <- s$dim
        return(r)
      }
      if (isTRUE(s$list)) {
        r <- lapply(s$items, rb)
        oldClass(r) <- s$cls
        return(r)
      }
      s$v
    }
    lapply(specs, rb)
  }
  list(theta = theta, parents = parents, rebuild = rebuild)
}

# value and Jacobian of a vector function at z0
hof_jac <- function(fn, z0) {
  old <- tape_start()
  on.exit(tape_restore(old), add = TRUE)
  zl <- ad_leaf(z0)
  out <- fn(zl)
  v <- as.numeric(val(out))
  J <- matrix(0, length(v), length(z0))
  if (inherits(out, "bcad")) {
    n0 <- .ad$n
    for (i in seq_along(v)) {
      .ad$n <- n0
      oi <- ad_select(out, i)
      .ad$adj <- vector("list", .ad$cap)
      ad_backward(oi)
      g <- .ad$adj[[zl$id]]
      if (!is.null(g)) J[i, ] <- rep_len(as.numeric(g), length(z0))
    }
  }
  list(value = v, J = J)
}

# the adjoint of a parent from the matching slice of a flat gradient
hof_split <- function(gflat, parents) {
  off <- 0L
  lapply(parents, function(p) {
    n <- length(p$v)
    r <- gflat[off + seq_len(n)]
    off <<- off + n
    if (!is.null(dim(p$v))) dim(r) <- dim(p$v)
    r
  })
}

hof_check_finite <- function(name, what, x) {
  v <- as.numeric(val(x))
  bad <- which(!is.finite(v))
  if (length(bad)) stop_domain(sprintf("%s: %s[%d] is %s, but must be finite!", name, what, bad[1L], format_scalar(v[bad[1L]])))
}

# ---- ODE integration ------------------------------------------------------------------------

ode_norm <- function(e, y, ynew, rtol, atol) {
  sc <- atol + rtol * pmax(abs(y), abs(ynew))
  sqrt(mean((e / sc)^2))
}

# Dormand and Prince 5(4) with step size control; returns the states at ts
ode_dopri <- function(Fn, t0, y0, ts, rtol, atol, maxsteps, name) {
  out <- matrix(0, length(ts), length(y0))
  t <- t0; y <- y0
  k1 <- Fn(t, y)
  tend <- ts[length(ts)]
  d0 <- sqrt(mean((y / (atol + rtol * abs(y)))^2)); d1 <- sqrt(mean((k1 / (atol + rtol * abs(y)))^2))
  h <- if (d0 < 1e-5 || d1 < 1e-5) 1e-6 else 0.01 * d0 / d1
  h <- min(h, abs(tend - t0))
  steps <- 0L
  for (j in seq_along(ts)) {
    while (t < ts[j]) {
      steps <- steps + 1L
      if (steps > maxsteps)
        stop_domain(sprintf("%s: Failed to integrate to next output time (%s) in less than max_num_steps steps", name, format_scalar(ts[j])))
      hh <- min(h, ts[j] - t)
      last <- hh >= ts[j] - t
      k2 <- Fn(t + hh / 5, y + hh * (k1 / 5))
      k3 <- Fn(t + 3 * hh / 10, y + hh * (3 / 40 * k1 + 9 / 40 * k2))
      k4 <- Fn(t + 4 * hh / 5, y + hh * (44 / 45 * k1 - 56 / 15 * k2 + 32 / 9 * k3))
      k5 <- Fn(t + 8 * hh / 9, y + hh * (19372 / 6561 * k1 - 25360 / 2187 * k2 + 64448 / 6561 * k3 - 212 / 729 * k4))
      k6 <- Fn(t + hh, y + hh * (9017 / 3168 * k1 - 355 / 33 * k2 + 46732 / 5247 * k3 + 49 / 176 * k4 - 5103 / 18656 * k5))
      ynew <- y + hh * (35 / 384 * k1 + 500 / 1113 * k3 + 125 / 192 * k4 - 2187 / 6784 * k5 + 11 / 84 * k6)
      k7 <- Fn(t + hh, ynew)
      e <- hh * (71 / 57600 * k1 - 71 / 16695 * k3 + 71 / 1920 * k4 - 17253 / 339200 * k5 + 22 / 525 * k6 - 1 / 40 * k7)
      err <- ode_norm(e, y, ynew, rtol, atol)
      if (!is.finite(err)) err <- Inf
      if (err <= 1) {
        t <- if (last) ts[j] else t + hh
        y <- ynew; k1 <- k7
      }
      fac <- if (err == 0) 5 else min(5, max(0.2, 0.9 * err^(-1 / 5)))
      if (err > 1 && hh < 1e-14 * max(1, abs(t)))
        stop_domain(sprintf("%s: step size became too small at t = %s", name, format_scalar(t)))
      if (!(last && err <= 1)) h <- hh * fac else h <- max(h, hh * fac)
    }
    out[j, ] <- y
  }
  out
}

# Rosenbrock method RODAS3 (Sandu et al. 1997; order 3, stiffly accurate,
# L-stable) for stiff problems. jacfn(t, y, F0) gives the Jacobian of the
# system and its time derivative. Stage K_i solves
# (I / (h g) - J) K_i = F(t + alpha_i h, Y_i) + sum_j c_ij K_j / h + h gamma_i dF/dt.
ros_tab <- list(
  A = rbind(c(0, 0, 0, 0), c(0, 0, 0, 0), c(2, 0, 0, 0), c(2, 0, 1, 0)),
  C = rbind(c(0, 0, 0, 0), c(4, 0, 0, 0), c(1, -1, 0, 0), c(1, -1, -8 / 3, 0)),
  M = c(2, 0, 1, 1), E = c(0, 0, 0, 1), alpha = c(0, 0, 1, 1), gamma = c(0.5, 1.5, 0, 0),
  newf = c(TRUE, FALSE, TRUE, TRUE), g = 0.5, order = 3)

# Solver for the stage matrix I / (h gamma) - J. A forward sensitivity
# system has the block lower triangular Jacobian [Jy 0; B I(x)Jy]: every
# diagonal block is the state Jacobian, so one N x N factorization of
# D = I / (h gamma) - Jy solves D k_y = b_y and D K_S = B_S + B k_y for all
# sensitivity columns at once.
ros_stage_solver <- function(JT, hg) {
  if (is.null(JT$Jy)) {
    W <- qr(diag(nrow(JT$J)) / hg - JT$J)
    return(function(r) qr.coef(W, r))
  }
  N <- nrow(JT$Jy)
  D <- qr(diag(N) / hg - JT$Jy)
  function(r) {
    ky <- qr.coef(D, r[seq_len(N)])
    c(ky, as.vector(qr.coef(D, matrix(r[-seq_len(N)] + as.vector(JT$B %*% ky), N))))
  }
}

ode_ros <- function(Fn, t0, y0, ts, rtol, atol, maxsteps, name, jacfn) {
  tb <- ros_tab
  n <- length(y0)
  out <- matrix(0, length(ts), n)
  t <- t0; y <- y0
  F0 <- Fn(t, y)
  tend <- ts[length(ts)]
  sc <- atol + rtol * abs(y)
  h <- min(abs(tend - t0), 0.01 * max(1e-6, sqrt(mean((y / sc)^2))) / max(1e-6, sqrt(mean((F0 / sc)^2))))
  steps <- 0L
  JT <- jacfn(t, y, F0)
  K <- matrix(0, n, 4L)
  for (j in seq_along(ts)) {
    while (t < ts[j]) {
      steps <- steps + 1L
      if (steps > maxsteps)
        stop_domain(sprintf("%s: Failed to integrate to next output time (%s) in less than max_num_steps steps", name, format_scalar(ts[j])))
      hh <- min(h, ts[j] - t)
      last <- hh >= ts[j] - t
      W <- ros_stage_solver(JT, hh * tb$g)
      Fi <- F0
      for (st in 1:4) {
        if (st > 1L && tb$newf[st]) {
          Yi <- y + as.vector(K[, seq_len(st - 1L), drop = FALSE] %*% tb$A[st, seq_len(st - 1L)])
          Fi <- Fn(t + tb$alpha[st] * hh, Yi)
        }
        rhs <- Fi + hh * tb$gamma[st] * JT$T
        if (st > 1L) rhs <- rhs + as.vector(K[, seq_len(st - 1L), drop = FALSE] %*% tb$C[st, seq_len(st - 1L)]) / hh
        K[, st] <- W(rhs)
      }
      ynew <- y + as.vector(K %*% tb$M)
      err <- ode_norm(as.vector(K %*% tb$E), y, ynew, rtol, atol)
      if (!is.finite(err)) err <- Inf
      if (err <= 1) {
        t <- if (last) ts[j] else t + hh
        y <- ynew
        F0 <- Fn(t, y)
        JT <- jacfn(t, y, F0)
      }
      fac <- if (err == 0) 5 else min(5, max(0.2, 0.9 * err^(-1 / tb$order)))
      if (err > 1 && hh < 1e-14 * max(1, abs(t)))
        stop_domain(sprintf("%s: step size became too small at t = %s", name, format_scalar(t)))
      if (!(last && err <= 1)) h <- hh * fac else h <- max(h, hh * fac)
    }
    out[j, ] <- y
  }
  out
}

# Solve y' = f(t, y, args) and differentiate the states at ts with respect
# to the initial state, the times and the differentiable arguments.
ode_solve <- function(callf, y0, t0, ts, args, rtol, atol, maxsteps, method, name) {
  hof_check_finite(name, "initial state", y0)
  if (!is.finite(val(t0))) stop_domain(sprintf("%s: initial time is %s, but must be finite!", name, format_scalar(val(t0))))
  hof_check_finite(name, "times", ts)
  tv <- as.numeric(val(ts)); t0v <- as.numeric(val(t0))
  if (!length(tv)) stop_shape(sprintf("%s: times has size 0, but must have a non-zero size", name))
  if (is.unsorted(tv)) stop_domain(sprintf("%s: times is not a valid sorted vector", name))
  if (tv[1L] <= t0v) stop_domain(sprintf("%s: initial time is %s, but must be less than %s", name, format_scalar(t0v), format_scalar(tv[1L])))
  if (!(rtol > 0) || !all(atol > 0)) stop_domain(sprintf("%s: relative and absolute tolerances must be positive", name))
  if (!(maxsteps > 0)) stop_domain(sprintf("%s: max_num_steps must be positive", name))
  y0v <- as.numeric(val(y0))
  N <- length(y0v)
  fl <- hof_flatten(args)
  pv <- fl$rebuild(fl$theta)
  plain <- function(t, y) {
    r <- as.numeric(val(callf(t, y, pv)))
    if (length(r) != N) stop_shape(sprintf("%s: the ODE right hand side returned %d values for a state of size %d", name, length(r), N))
    r
  }
  sens_y0 <- inherits(y0, "bcad") || inherits(t0, "bcad")
  P <- (if (sens_y0) N else 0L) + length(fl$theta)
  fd_t <- function(Fn, t, y, F0) {
    dt <- sqrt(.Machine$double.eps) * max(abs(t), 1)
    (Fn(t + dt, y) - F0) / dt
  }
  if (P == 0L) {
    jac <- function(t, y, F0) {
      J <- hof_jac(function(z) callf(t, z, pv), y)$J
      list(J = J, T = fd_t(plain, t, y, F0))
    }
    Y <- if (method == "bdf") ode_ros(plain, t0v, y0v, tv, rtol, atol, maxsteps, name, jac)
         else ode_dopri(plain, t0v, y0v, tv, rtol, atol, maxsteps, name)
    if (!inherits(ts, "bcad")) return(Y)
    S <- NULL
  } else {
    np <- length(fl$theta)
    aug <- function(t, Yv) {
      y <- Yv[seq_len(N)]
      S <- matrix(Yv[-seq_len(N)], N, P)
      jj <- hof_jac(function(z) callf(t, ad_gather(z, seq_len(N)), fl$rebuild(ad_gather(z, N + seq_len(np)))), c(y, fl$theta))
      if (length(jj$value) != N) stop_shape(sprintf("%s: the ODE right hand side returned %d values for a state of size %d", name, length(jj$value), N))
      Jy <- jj$J[, seq_len(N), drop = FALSE]
      dS <- Jy %*% S
      if (np) dS[, P - np + seq_len(np)] <- dS[, P - np + seq_len(np)] + jj$J[, N + seq_len(np)]
      c(jj$value, as.vector(dS))
    }
    S0 <- matrix(0, N, P)
    if (sens_y0) S0[, seq_len(N)] <- diag(N)
    # Jacobian of the augmented system: the state block and the identical
    # sensitivity diagonal blocks exactly, the dependence of the
    # sensitivity equations on the state (the first block column B) by
    # finite differences. Only the blocks are kept (see ros_stage_solver()).
    jac_aug <- function(t, Yv, F0) {
      y <- Yv[seq_len(N)]
      Jy <- hof_jac(function(z) callf(t, z, pv), y)$J
      B <- matrix(0, N * P, N)
      for (k in seq_len(N)) {
        dk <- sqrt(.Machine$double.eps) * max(abs(y[k]), 1)
        Yk <- Yv; Yk[k] <- Yk[k] + dk
        B[, k] <- (aug(t, Yk)[-seq_len(N)] - F0[-seq_len(N)]) / dk
      }
      list(Jy = Jy, B = B, T = fd_t(aug, t, Yv, F0))
    }
    # the tolerances apply to the states and their sensitivities
    Y0 <- c(y0v, as.vector(S0))
    full <- if (method == "bdf") ode_ros(aug, t0v, Y0, tv, rtol, atol, maxsteps, name, jac_aug)
            else ode_dopri(aug, t0v, Y0, tv, rtol, atol, maxsteps, name)
    Y <- full[, seq_len(N), drop = FALSE]
    S <- lapply(seq_along(tv), function(j) matrix(full[j, -seq_len(N)], N, P))
  }
  Fts <- if (inherits(ts, "bcad")) t(vapply(seq_along(tv), function(j) plain(tv[j], Y[j, ]), numeric(N))) else NULL
  if (!is.null(Fts) && N == 1L) Fts <- matrix(Fts, ncol = 1L)
  f0 <- if (inherits(t0, "bcad")) plain(t0v, y0v) else NULL
  ad_node(Y, c(list(y0, t0, ts), fl$parents), function(g) {
    g <- matrix(g, length(tv), N)
    gs <- numeric(P)
    if (!is.null(S)) for (j in seq_along(tv)) gs <- gs + as.vector(crossprod(S[[j]], g[j, ]))
    gy0 <- if (inherits(y0, "bcad")) gs[seq_len(N)] else NULL
    gt0 <- if (!is.null(f0)) -sum(gs[seq_len(N)] * f0) else NULL
    gts <- if (!is.null(Fts)) rowSums(g * Fts) else NULL
    gth <- if (length(fl$theta)) hof_split(gs[P - length(fl$theta) + seq_along(fl$theta)], fl$parents) else list()
    c(list(gy0, gt0, gts), gth)
  })
}

hof_call <- function(f, ctx) function(t, y, a) do.call(f, c(list(t, y), a, list(ctx)))

rt_ode <- function(f, fixed, va, ctx, opts) {
  name <- opts$name
  if (opts$method == "adjoint") {
    # Stan's argument checks; the solution and its gradient come from forward
    # sensitivities with the forward tolerances, which give the same values
    # as the adjoint method up to the tolerances
    N <- length(val(fixed[[1L]]))
    pos <- function(k, what) {
      v <- val(fixed[[k]])
      if (!all(is.finite(v)) || any(v <= 0)) stop_domain(sprintf("%s: %s must be positive and finite", name, what))
    }
    pos(4L, "relative_tolerance_forward"); pos(5L, "absolute_tolerance_forward")
    pos(6L, "relative_tolerance_backward"); pos(7L, "absolute_tolerance_backward")
    pos(8L, "relative_tolerance_quadrature"); pos(9L, "absolute_tolerance_quadrature")
    pos(10L, "max_num_steps"); pos(11L, "num_steps_between_checkpoints")
    for (k in c(5L, 7L)) if (length(val(fixed[[k]])) != N)
      stop_shape(sprintf("%s: the absolute tolerance vectors must have the size of the initial state (%d)", name, N))
    for (k in 12:14) if (!(val(fixed[[k]]) %in% c(1, 2)))
      stop_domain(sprintf("%s: %s must be 1 or 2", name, c("interpolation_polynomial", "solver_forward", "solver_backward")[k - 11L]))
    rtol <- fixed[[4L]]; atol <- val(fixed[[5L]]); maxsteps <- fixed[[10L]]
  } else if (opts$tol) {
    rtol <- fixed[[4L]]; atol <- fixed[[5L]]; maxsteps <- fixed[[6L]]
  } else {
    rtol <- 1e-6; atol <- 1e-6; maxsteps <- 1e6
  }
  ode_solve(hof_call(f, ctx), fixed[[1L]], fixed[[2L]], fixed[[3L]], va, val(rtol), val(atol), val(maxsteps), opts$method, name)
}

# ---- algebraic equations ----------------------------------------------------------------------

# damped Newton iterations on f(x) = 0 with a backtracking line search
alg_newton <- function(fx, jx, x, ftol, maxsteps, name) {
  fv <- fx(x)
  it <- 0L
  while (!(max(abs(fv)) <= ftol)) {
    it <- it + 1L
    if (it > maxsteps) stop_domain(sprintf("%s: maximum number of iterations (%s) was exceeded in the solve.", name, format_scalar(maxsteps)))
    J <- jx(x)
    dx <- tryCatch(solve(J, -fv), error = function(e) NULL)
    if (is.null(dx) || !all(is.finite(dx))) dx <- -as.vector(crossprod(J, fv))
    f0 <- sum(fv^2)
    lam <- 1
    repeat {
      xn <- x + lam * dx
      fn <- fx(xn)
      if (all(is.finite(fn)) && sum(fn^2) <= (1 - 1e-4 * lam) * f0) break
      lam <- lam / 2
      if (lam < 1e-12) break
    }
    if (!all(is.finite(fn))) break
    x <- xn; fv <- fn
  }
  if (!(max(abs(fv)) <= ftol))
    stop_domain(sprintf("%s: the norm of the algebraic function is %s but should be lower than the function tolerance: %s.",
                        name, format_scalar(sqrt(sum(fv^2))), format_scalar(ftol)))
  x
}

# Solve f(x, args) = 0; dx/dtheta = -Jx^-1 Jtheta (implicit function theorem)
alg_solve <- function(callf, guess, args, ftol, maxsteps, name) {
  hof_check_finite(name, "initial guess", guess)
  if (!(ftol > 0)) stop_domain(sprintf("%s: function_tolerance must be positive", name))
  if (!(maxsteps > 0)) stop_domain(sprintf("%s: max_num_steps must be positive", name))
  x0 <- as.numeric(val(guess))
  N <- length(x0)
  if (N == 0L) stop_shape(sprintf("%s: initial guess has size 0, but must have a non-zero size", name))
  fl <- hof_flatten(args)
  pv <- fl$rebuild(fl$theta)
  fx <- function(x) {
    r <- as.numeric(val(callf(x, pv)))
    if (length(r) != N) stop_shape(sprintf("%s: the algebraic system returned %d values for %d unknowns", name, length(r), N))
    r
  }
  jx <- function(x) hof_jac(function(z) callf(z, pv), x)$J
  x <- alg_newton(fx, jx, x0, ftol, maxsteps, name)
  if (!length(fl$theta)) return(x)
  np <- length(fl$theta)
  J <- hof_jac(function(z) callf(ad_gather(z, seq_len(N)), fl$rebuild(ad_gather(z, N + seq_len(np)))), c(x, fl$theta))$J
  Jx <- J[, seq_len(N), drop = FALSE]
  Jt <- J[, N + seq_len(np), drop = FALSE]
  ad_node(x, fl$parents, function(g) {
    w <- solve(t(Jx), as.numeric(g))
    hof_split(-as.vector(crossprod(Jt, w)), fl$parents)
  })
}

rt_solve <- function(f, fixed, va, ctx, opts) {
  ftol <- if (opts$tol) fixed[[3L]] else 1e-6
  maxsteps <- if (opts$tol) fixed[[4L]] else 200
  if (opts$tol && opts$method == "newton" && !(val(fixed[[2L]]) > 0)) stop_domain(sprintf("%s: scaling_step_size must be positive", opts$name))
  if (opts$tol && opts$method == "powell" && !(val(fixed[[2L]]) > 0)) stop_domain(sprintf("%s: relative_tolerance must be positive", opts$name))
  alg_solve(function(x, a) do.call(f, c(list(x), a, list(ctx))), fixed[[1L]], va, val(ftol), val(maxsteps), opts$name)
}

rt_solve_legacy <- function(f, fixed, va, ctx, opts) {
  ftol <- if (length(fixed) == 7L) fixed[[6L]] else 1e-6
  maxsteps <- if (length(fixed) == 7L) fixed[[7L]] else 1000
  alg_solve(function(x, a) do.call(f, c(list(x), a, list(ctx))), fixed[[1L]], fixed[2:4], val(ftol), val(maxsteps), opts$name)
}

# ---- 1D integration -----------------------------------------------------------------------------

# Gradient of the integrand in theta at one integration point. Each gradient
# component has its own adaptive quadrature, but they mostly visit the same
# points, so the whole gradient is kept per point (x, xc) for this integral
# at this theta. Integrands cannot call _rng or _lp functions, so the
# gradient at a point is fixed. Storage is capped at 1e6 doubles; points
# past the cap are recomputed.
quad_grad <- function(one, thv) {
  memo <- new.env(hash = TRUE, parent = emptyenv())
  room <- max(1L, 1e6 %/% max(1L, length(thv)))
  function(x, xc) {
    key <- sprintf("%a %a", x, xc)
    g <- memo[[key]]
    if (!is.null(g)) return(g)
    g <- ad_value_grad(function(z) one(x, xc, z), thv)$grad
    if (room > 0L) { assign(key, g, envir = memo); room <<- room - 1L }
    g
  }
}

# integrate_1d(f, a, b, theta, x_r, x_i[, rel_tol]) with f(x, xc, theta, x_r, x_i);
# xc is the signed distance to the nearest finite limit (NaN when a limit is
# infinite), as in Stan; intervals crossing zero are split there
rt_integrate_1d <- function(f, fixed, va, ctx, opts) {
  name <- opts$name
  a <- fixed[[1L]]; b <- fixed[[2L]]; theta <- fixed[[3L]]
  xr <- fixed[[4L]]; xi <- fixed[[5L]]
  rtol <- if (length(fixed) == 6L) val(fixed[[6L]]) else sqrt(.Machine$double.eps)
  av <- val(a); bv <- val(b)
  if (is.nan(av) || is.nan(bv)) stop_domain(sprintf("%s: integration limits must not be NaN", name))
  if (!(rtol >= 0)) stop_domain(sprintf("%s: relative_tolerance must be non-negative", name))
  if (av > bv) stop_domain(sprintf("%s: lower limit (%s) must be less than or equal to the upper limit (%s)", name, format_scalar(av), format_scalar(bv)))
  if (av == bv) return(if (any_ad(a, b, theta)) ad_node(0, list(a, b, theta), function(g) list(NULL, NULL, NULL)) else 0)
  thv <- val(theta)
  one <- function(x, xc, th) do.call(f, list(x, xc, th, xr, xi, ctx))
  infinite <- is.infinite(av) || is.infinite(bv)
  parts <- if (!infinite && av < 0 && bv > 0) list(c(av, 0), c(0, bv)) else list(c(av, bv))
  quad <- function(g) {
    tot <- 0
    for (p in parts) {
      lo <- p[1L]; hi <- p[2L]
      h <- function(x) vapply(x, function(xx) {
        xc <- if (infinite) NaN else if (xx - lo <= hi - xx) lo - xx else hi - xx
        g(xx, xc)
      }, 0)
      r <- tryCatch(stats::integrate(h, lo, hi, rel.tol = max(rtol, 50 * .Machine$double.eps), subdivisions = 1000L),
                    error = function(e) stop_domain(sprintf("%s: %s", name, conditionMessage(e))))
      tot <- tot + r$value
    }
    tot
  }
  value <- quad(function(x, xc) as.numeric(val(one(x, xc, thv))))
  if (!any_ad(a, b, theta)) return(value)
  ga <- if (inherits(a, "bcad")) -as.numeric(val(one(av, 0, thv))) else NULL
  gb <- if (inherits(b, "bcad")) as.numeric(val(one(bv, 0, thv))) else NULL
  gth <- NULL
  if (inherits(theta, "bcad")) {
    gr <- quad_grad(one, thv)
    gth <- vapply(seq_along(thv), function(k) quad(function(x, xc) gr(x, xc)[k]), 0)
  }
  ad_node(value, list(a, b, theta), function(g) list(if (!is.null(ga)) g * ga, if (!is.null(gb)) g * gb, if (!is.null(gth)) g * gth))
}

# integrate_1d_double_exponential and integrate_1d_gauss_kronrod: the
# integrand takes (x, xc, ...). Both use adaptive Gauss-Kronrod quadrature
# (stats::integrate); for the double exponential form the integrand also
# receives xc, the distance to the nearer finite limit, as Stan's does.
# ponytail: max_refinements is accepted but stats::integrate's own
# subdivision limit is used
rt_integrate_1d_v <- function(f, fixed, va, ctx, opts) {
  name <- opts$name
  a <- fixed[[1L]]; b <- fixed[[2L]]
  rtol <- if (opts$tol) val(fixed[[3L]]) else sqrt(.Machine$double.eps)
  atol <- if (opts$tol) val(fixed[[4L]]) else sqrt(.Machine$double.eps)
  if (opts$tol) {
    if (!(rtol >= 0)) stop_domain(sprintf("%s: relative_tolerance must be non-negative", name))
    if (!(atol >= 0)) stop_domain(sprintf("%s: absolute_tolerance must be non-negative", name))
    if (!(val(fixed[[5L]]) >= 1)) stop_domain(sprintf("%s: max_refinements must be positive", name))
  }
  av <- val(a); bv <- val(b)
  if (is.nan(av) || is.nan(bv)) stop_domain(sprintf("%s: integration limits must not be NaN", name))
  if (av > bv) stop_domain(sprintf("%s: lower limit (%s) must be less than or equal to the upper limit (%s)", name, format_scalar(av), format_scalar(bv)))
  fl <- hof_flatten(va)
  thv <- fl$theta
  if (av == bv) return(if (any_ad(a, b) || length(fl$parents)) ad_node(0, c(list(a, b), fl$parents), function(g) vector("list", 2L + length(fl$parents))) else 0)
  one <- function(x, xc, z) do.call(f, c(list(x, xc), fl$rebuild(z), list(ctx)))
  infinite <- is.infinite(av) || is.infinite(bv)
  gk <- opts$method == "gauss_kronrod"
  quad <- function(g) {
    h <- function(x) vapply(x, function(xx) {
      xc <- if (gk || infinite) NaN else if (xx - av <= bv - xx) av - xx else bv - xx
      g(xx, xc)
    }, 0)
    r <- tryCatch(stats::integrate(h, av, bv, rel.tol = max(rtol, 50 * .Machine$double.eps), abs.tol = atol, subdivisions = 1000L),
                  error = function(e) stop_domain(sprintf("%s: %s", name, conditionMessage(e))))
    r$value
  }
  value <- quad(function(x, xc) as.numeric(val(one(x, xc, thv))))
  if (!any_ad(a, b) && !length(fl$parents)) return(value)
  ga <- if (inherits(a, "bcad")) -as.numeric(val(one(av, 0, thv))) else NULL
  gb <- if (inherits(b, "bcad")) as.numeric(val(one(bv, 0, thv))) else NULL
  gr <- quad_grad(one, thv)
  gth <- vapply(seq_along(thv), function(k) quad(function(x, xc) gr(x, xc)[k]), 0)
  ad_node(value, c(list(a, b), fl$parents), function(g) c(list(if (!is.null(ga)) g * ga, if (!is.null(gb)) g * gb),
                                                           lapply(hof_split(gth, fl$parents), function(d) g * d)))
}

# integrate_ode_*: the deprecated interface; the function takes
# (t, y, theta, x_r, x_i) with y an array, which is stored as a vector here
rt_ode_legacy <- function(f, fixed, va, ctx, opts) {
  ctl <- length(fixed) == 9L
  rtol <- if (ctl) val(fixed[[7L]]) else 1e-6
  atol <- if (ctl) val(fixed[[8L]]) else 1e-6
  maxsteps <- if (ctl) val(fixed[[9L]]) else 1e6
  ode_solve(hof_call(f, ctx), fixed[[1L]], fixed[[2L]], fixed[[3L]], list(fixed[[4L]], fixed[[5L]], fixed[[6L]]),
            rtol, atol, maxsteps, opts$method, opts$name)
}

# ---- map_rect and reduce_sum ------------------------------------------------------------------------

rt_map_rect <- function(f, fixed, va, ctx, opts) {
  phi <- fixed[[1L]]; thetas <- fixed[[2L]]; xr <- fixed[[3L]]; xi <- fixed[[4L]]
  J <- storage_dims(thetas)[1L]
  if (storage_dims(xr)[1L] != J || storage_dims(xi)[1L] != J)
    stop_shape(sprintf("%s: the job parameters, real data and int data must have the same number of jobs (%d, %d, %d)",
                       opts$name, J, storage_dims(xr)[1L], storage_dims(xi)[1L]))
  row <- function(x, j) {
    d <- storage_dims(x)
    if (length(d) > 1L && d[2L] == 0L) return(if (is.integer(val(x))) integer(0) else numeric(0))
    rt_index(x, "single", list(j))
  }
  outs <- lapply(seq_len(J), function(j) f(phi, row(thetas, j), row(xr, j), row(xi, j), ctx))
  ad_c(outs)
}

# a single call over the whole range gives the same sum as any slicing
rt_reduce_sum <- function(f, fixed, va, ctx, opts) {
  x <- fixed[[1L]]
  n <- if (is_tarr(x)) length(x) else storage_dims(x)[1L]
  if (n == 0L) return(0)
  do.call(f, c(list(x, 1L, as.integer(n)), va, list(ctx)))
}

# ---- differential algebraic equations --------------------------------------------------------------

# Radau IIA, 3 stages (order 5), for F(t, y, y') = 0 of index 1. The stage
# derivatives Z_i solve F(t + c_i h, y + h sum_j a_ij Z_j, Z_i) = 0 by Newton
# iterations; the error comes from step doubling.
radau_tab <- local({
  r6 <- sqrt(6)
  list(c = c((4 - r6) / 10, (4 + r6) / 10, 1),
       A = rbind(c((88 - 7 * r6) / 360, (296 - 169 * r6) / 1800, (-2 + 3 * r6) / 225),
                 c((296 + 169 * r6) / 1800, (88 + 7 * r6) / 360, (-2 - 3 * r6) / 225),
                 c((16 - r6) / 36, (16 + r6) / 36, 1 / 9)))
})

dae_solve <- function(callf, yy0, yp0, t0, ts, args, rtol, atol, maxsteps, name) {
  hof_check_finite(name, "initial state", yy0)
  hof_check_finite(name, "initial state derivative", yp0)
  if (!is.finite(val(t0))) stop_domain(sprintf("%s: initial time is %s, but must be finite!", name, format_scalar(val(t0))))
  hof_check_finite(name, "times", ts)
  tv <- as.numeric(val(ts)); t0v <- as.numeric(val(t0))
  if (!length(tv)) stop_shape(sprintf("%s: times has size 0, but must have a non-zero size", name))
  if (is.unsorted(tv)) stop_domain(sprintf("%s: times is not a valid sorted vector", name))
  if (tv[1L] <= t0v) stop_domain(sprintf("%s: initial time is %s, but must be less than %s", name, format_scalar(t0v), format_scalar(tv[1L])))
  if (!(rtol > 0 && is.finite(rtol)) || !(atol > 0 && is.finite(atol))) stop_domain(sprintf("%s: relative and absolute tolerances must be positive and finite", name))
  if (!(maxsteps > 0)) stop_domain(sprintf("%s: max_num_steps must be positive", name))
  y <- as.numeric(val(yy0)); yp <- as.numeric(val(yp0))
  N <- length(y)
  if (length(yp) != N) stop_shape(sprintf("%s: initial state has size %d but initial state derivative has size %d", name, N, length(yp)))
  fl <- hof_flatten(args)
  pv <- fl$rebuild(fl$theta)
  np <- length(fl$theta)
  sens_y0 <- inherits(yy0, "bcad") || inherits(t0, "bcad")
  P <- (if (sens_y0) N else 0L) + np
  tb <- radau_tab
  # residual and Jacobians [F_y, F_yp(, F_theta)] at (t, y, yp)
  resj <- function(t, y, yp, with_theta = FALSE) {
    if (with_theta && np) {
      r <- hof_jac(function(z) callf(t, ad_gather(z, seq_len(N)), ad_gather(z, N + seq_len(N)), fl$rebuild(ad_gather(z, 2L * N + seq_len(np)))), c(y, yp, fl$theta))
    } else {
      r <- hof_jac(function(z) callf(t, ad_gather(z, seq_len(N)), ad_gather(z, N + seq_len(N)), pv), c(y, yp))
    }
    if (length(r$value) != N) stop_shape(sprintf("%s: the residual function returned %d values for a state of size %d", name, length(r$value), N))
    r
  }
  # one Radau step; returns the new state, its derivative and the stages
  rstep <- function(t, y, z0, h) {
    Z <- matrix(z0, N, 3L)
    for (it in 1:12) {
      R <- numeric(3L * N); J <- matrix(0, 3L * N, 3L * N)
      for (i in 1:3) {
        Yi <- y + h * as.vector(Z %*% tb$A[i, ])
        r <- resj(t + tb$c[i] * h, Yi, Z[, i])
        ri <- (i - 1L) * N + seq_len(N)
        R[ri] <- r$value
        Fy <- r$J[, seq_len(N), drop = FALSE]; Fp <- r$J[, N + seq_len(N), drop = FALSE]
        for (j in 1:3) J[ri, (j - 1L) * N + seq_len(N)] <- h * tb$A[i, j] * Fy + (if (i == j) Fp else 0)
      }
      if (!all(is.finite(R))) return(NULL)
      dZ <- tryCatch(solve(J, -R), error = function(e) NULL)
      if (is.null(dZ)) return(NULL)
      Z <- Z + matrix(dZ, N, 3L)
      if (max(abs(h * dZ) / (atol + rtol * abs(y))) < 1e-3) break
      if (it == 12L) return(NULL)
    }
    list(y = y + h * as.vector(Z %*% tb$A[3L, ]), Z = Z, t = t, h = h, y0 = y)
  }
  # exact sensitivity of a Radau step: the stage equations differentiated
  sstep <- function(st, S) {
    h <- st$h
    M <- matrix(0, 3L * N, 3L * N)
    B <- matrix(0, 3L * N, P)
    for (i in 1:3) {
      Yi <- st$y0 + h * as.vector(st$Z %*% tb$A[i, ])
      r <- resj(st$t + tb$c[i] * h, Yi, st$Z[, i], with_theta = TRUE)
      ri <- (i - 1L) * N + seq_len(N)
      Fy <- r$J[, seq_len(N), drop = FALSE]; Fp <- r$J[, N + seq_len(N), drop = FALSE]
      for (j in 1:3) M[ri, (j - 1L) * N + seq_len(N)] <- h * tb$A[i, j] * Fy + (if (i == j) Fp else 0)
      B[ri, ] <- Fy %*% S
      if (np) B[ri, P - np + seq_len(np)] <- B[ri, P - np + seq_len(np)] + r$J[, 2L * N + seq_len(np)]
    }
    dZ <- solve(M, -B)
    S + h * Reduce(`+`, lapply(1:3, function(j) tb$A[3L, j] * dZ[(j - 1L) * N + seq_len(N), , drop = FALSE]))
  }
  out <- matrix(0, length(tv), N); Sout <- vector("list", length(tv)); ypout <- matrix(0, length(tv), N)
  S <- matrix(0, N, P)
  if (sens_y0) S[, seq_len(N)] <- diag(N)
  t <- t0v
  h <- min(tv[length(tv)] - t0v, 0.01 * max(1, abs(t0v)))
  steps <- 0L
  for (k in seq_along(tv)) {
    while (t < tv[k]) {
      steps <- steps + 1L
      if (steps > maxsteps)
        stop_domain(sprintf("%s: Failed to integrate to next output time (%s) in less than max_num_steps steps", name, format_scalar(tv[k])))
      hh <- min(h, tv[k] - t)
      last <- hh >= tv[k] - t
      full <- rstep(t, y, yp, hh)
      s1 <- if (!is.null(full)) rstep(t, y, yp, hh / 2) else NULL
      s2 <- if (!is.null(s1)) rstep(t + hh / 2, s1$y, s1$Z[, 3L], hh / 2) else NULL
      if (is.null(s2)) {
        err <- Inf
      } else {
        err <- sqrt(mean(((s2$y - full$y) / 31 / (atol + rtol * pmax(abs(y), abs(s2$y))))^2))
      }
      if (err <= 1) {
        if (P) S <- sstep(s2, sstep(s1, S))
        t <- if (last) tv[k] else t + hh
        y <- s2$y; yp <- s2$Z[, 3L]
      }
      fac <- if (err == 0) 4 else min(4, max(0.2, 0.9 * err^(-1 / 6)))
      if (err > 1 && hh < 1e-14 * max(1, abs(t)))
        stop_domain(sprintf("%s: step size became too small at t = %s", name, format_scalar(t)))
      if (!(last && err <= 1)) h <- hh * fac else h <- max(h, hh * fac)
    }
    out[k, ] <- y; ypout[k, ] <- yp; Sout[[k]] <- S
  }
  if (P == 0L && !inherits(ts, "bcad")) return(out)
  yp0v <- as.numeric(val(yp0))
  ad_node(out, c(list(yy0, yp0, t0, ts), fl$parents), function(g) {
    g <- matrix(g, length(tv), N)
    gs <- numeric(P)
    if (P) for (k in seq_along(tv)) gs <- gs + as.vector(crossprod(Sout[[k]], g[k, ]))
    gy0 <- if (inherits(yy0, "bcad")) gs[seq_len(N)] else NULL
    gt0 <- if (inherits(t0, "bcad")) -sum(gs[seq_len(N)] * yp0v) else NULL
    gts <- if (inherits(ts, "bcad")) rowSums(g * ypout) else NULL
    gth <- if (np) hof_split(gs[P - np + seq_len(np)], fl$parents) else list()
    c(list(gy0, NULL, gt0, gts), gth)
  })
}

rt_dae <- function(f, fixed, va, ctx, opts) {
  ctl <- if (opts$tol) fixed[5:7] else list(1e-10, 1e-10, 1e8)
  dae_solve(function(t, y, yp, a) do.call(f, c(list(t, y, yp), a, list(ctx))), fixed[[1L]], fixed[[2L]], fixed[[3L]], fixed[[4L]],
            va, val(ctl[[1L]]), val(ctl[[2L]]), val(ctl[[3L]]), opts$name)
}
