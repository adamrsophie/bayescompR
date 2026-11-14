# Embedded Laplace approximation.
#
# The latent vector theta has the prior N(0, K(phi_K)) and the log likelihood
# f(theta, phi_l). At the mode theta_hat of f(theta) - theta' K^-1 theta / 2,
# with W = -d2f/dtheta2 and a = K^-1 theta_hat, the approximate log marginal
# likelihood is
#   L = f(theta_hat) - theta_hat' a / 2 - log|I + K W| / 2.
# Its gradient follows Rasmussen and Williams (2006, section 5.5.1) for a
# general likelihood. With C = (I + K W)^-1 K, T(theta, phi_l) = tr(C H) for
# the Hessian H of f and s2 = dT/dtheta / 2:
#   dL/dphi_l = df/dphi_l + dT/dphi_l / 2 + s2' C d2f/dtheta dphi_l
#   dL/dphi_K = sum(M * dK/dphi_K),
#     M = a a' / 2 - t(W (I + K W)^-1) / 2 + u a',  u = (I + K W)^-T s2.
# The Hessian of f is block diagonal with the declared block size; it comes
# from central differences of the AD gradient, and its derivatives from
# central differences of those Hessians.

# ---- semantic analysis -------------------------------------------------------------------------

laplace_builtin <- list(poisson_log = 0L, bernoulli_logit = 0L, neg_binomial_2_log = 1L)
laplace_generic_names <- c("laplace_marginal", "laplace_marginal_tol", "laplace_latent_rng", "laplace_latent_tol_rng",
                           "laplace_latent_solve", "laplace_latent_solve_tol")
laplace_names <- c(laplace_generic_names, unlist(lapply(names(laplace_builtin), function(l)
  c(paste0("laplace_marginal_", l, "_lpmf"), paste0("laplace_marginal_tol_", l, "_lpmf"),
    paste0("laplace_latent_", l, "_rng"), paste0("laplace_latent_tol_", l, "_rng")))))

laplace_opts_ty <- function() ty_tuple(list(T_vec, T_real, T_int, T_int, T_int, T_int))

sem_laplace_fun <- function(s, a, what, nm) {
  if (a$kind != "var") stop_semantic(sprintf("argument %s of '%s' must be the name of a user defined function", what, nm), a$loc)
  cands <- s$ufuns[[a$name]]
  if (is.null(cands)) stop_semantic(sprintf("'%s' is not a user defined function", a$name), a$loc)
  cands
}

sem_laplace <- function(s, e) {
  nm <- e$name
  generic <- nm %in% laplace_generic_names
  tol <- grepl("_tol", nm)
  kind <- if (grepl("_rng$", nm)) "rng" else if (grepl("solve", nm)) "solve" else "marginal"
  if (kind == "rng") {
    ok <- s$block %in% c("transformed data", "generated quantities") || (!is.null(s$fun) && endsWith(s$fun$name, "_rng"))
    if (!ok) stop_semantic(sprintf("random number generator '%s' is only allowed in transformed data, generated quantities and _rng functions", nm), e$loc)
  }
  A <- e$args
  if (generic) {
    if (length(A) < 5L) stop_semantic(sprintf("'%s' needs a likelihood function, its arguments as a tuple, the Hessian block size, a covariance function and its arguments as a tuple", nm), e$loc)
    lcands <- sem_laplace_fun(s, A[[1L]], 1L, nm)
    largs <- sem_expr(s, A[[2L]])
    if (!is_tuple_ty(largs$type) || largs$type$nd > 0L) stop_semantic(sprintf("argument 2 of '%s' must be a tuple of the likelihood arguments", nm), A[[2L]]$loc)
    rest <- A[-(1:2)]
    lik <- NULL
    for (cd in lcands) if (hof_fun_ok(cd, c(list(T_vec), largs$type$elems)) && ty_eq(cd$ret, T_real)) { lik <- cd; break }
    if (is.null(lik)) stop_semantic(sprintf("function '%s' passed to '%s' must take (vector, %s) and return real", A[[1L]]$name, nm,
                                            paste(vapply(largs$type$elems, ty_str, ""), collapse = ", ")), e$loc)
    pre <- list(largs)
  } else {
    fam <- sub("^laplace_(marginal|latent)(_tol)?_(.*)_(lpmf|rng)$", "\\3", nm)
    nextra <- laplace_builtin[[fam]]
    if (kind == "marginal" && !isTRUE(e$cond_bar)) stop_semantic(sprintf("'%s' needs '|' after its first argument", nm), e$loc)
    nfix <- 3L + nextra
    if (length(A) < nfix + 3L) stop_semantic(sprintf("'%s' needs the outcomes, their group indices, %sthe mean, the Hessian block size, a covariance function and its arguments as a tuple",
                                                     nm, if (nextra) "the dispersion, " else ""), e$loc)
    want <- c(list(T_ia, T_ia), if (nextra) list(T_vec), list(T_vec))
    pre <- lapply(seq_len(nfix), function(i) sem_expr(s, A[[i]]))
    for (i in seq_len(nfix)) if (!ty_assignable(want[[i]], pre[[i]]$type))
      stop_semantic(sprintf("argument %d of '%s' must be %s, found %s", i, nm, ty_str(want[[i]]), ty_str(pre[[i]]$type)), A[[i]]$loc)
    for (i in 1:2) if (!isTRUE(pre[[i]]$data_only)) stop_semantic(sprintf("argument %d of '%s' must be data only", i, nm), A[[i]]$loc)
    rest <- A[-seq_len(nfix)]
    lik <- NULL
  }
  hbs <- sem_expr(s, rest[[1L]])
  if (!ty_eq(hbs$type, T_int) || !isTRUE(hbs$data_only)) stop_semantic(sprintf("the Hessian block size of '%s' must be a data only int", nm), rest[[1L]]$loc)
  kcands <- sem_laplace_fun(s, rest[[2L]], "for the covariance", nm)
  kargs <- sem_expr(s, rest[[3L]])
  if (!is_tuple_ty(kargs$type) || kargs$type$nd > 0L) stop_semantic(sprintf("the covariance function arguments of '%s' must be a tuple", nm), rest[[3L]]$loc)
  kf <- NULL
  for (cd in kcands) if (hof_fun_ok(cd, kargs$type$elems) && ty_eq(cd$ret, T_mat)) { kf <- cd; break }
  if (is.null(kf)) stop_semantic(sprintf("covariance function '%s' passed to '%s' must take (%s) and return matrix", rest[[2L]]$name, nm,
                                         paste(vapply(kargs$type$elems, ty_str, ""), collapse = ", ")), e$loc)
  ctl <- rest[-(1:3)]
  if (!tol && length(ctl)) stop_semantic(sprintf("'%s' takes no control arguments; use the _tol form", nm), e$loc)
  if (tol) {
    ctl <- lapply(ctl, function(a) sem_expr(s, a))
    if (length(ctl) == 1L && ty_eq(ctl[[1L]]$type, laplace_opts_ty())) {
      # a tuple of options, as returned by generate_laplace_options
    } else {
      want <- list(T_vec, T_real, T_int, T_int, T_int, T_int)
      if (length(ctl) != 6L) stop_semantic(sprintf("'%s' needs six control arguments (theta_0, tolerance, max_num_steps, solver, max_steps_line_search, allow_fallthrough) or a tuple of them", nm), e$loc)
      for (i in 1:6) {
        if (!ty_assignable(want[[i]], ctl[[i]]$type))
          stop_semantic(sprintf("control argument %d of '%s' must be %s, found %s", i, nm, ty_str(want[[i]]), ty_str(ctl[[i]]$type)), ctl[[i]]$loc)
        if (i > 1L && !isTRUE(ctl[[i]]$data_only)) stop_semantic(sprintf("control argument %d of '%s' must be data only", i, nm), ctl[[i]]$loc)
      }
    }
  }
  e$args <- c(pre, list(hbs, kargs), ctl)
  e$laplace <- list(name = nm, kind = kind, generic = generic, lik = lik, kf = kf, npre = length(pre),
                    fam = if (!generic) fam, nctl = length(ctl))
  if (kind == "solve" && s$block != "generated quantities" && !all(vapply(e$args, function(a) isTRUE(a$data_only), TRUE)))
    stop_semantic(sprintf("all arguments to '%s' must be data only if used outside the generated quantities block", nm), e$loc)
  e$type <- switch(kind, marginal = T_real, rng = T_vec, solve = ty_tuple(list(T_vec, T_mat)))
  e$data_only <- kind == "solve"
  e
}

cg_laplace <- function(e) {
  L <- e$laplace
  codes <- lapply(e$args, cg_expr)
  np <- L$npre
  ctl <- if (L$nctl) codes[np + 2L + seq_len(L$nctl)] else list()
  if (L$nctl == 6L) ctl[[1L]] <- promote_to(ctl[[1L]], e$args[[np + 3L]]$type, T_vec)
  lik <- if (L$generic) as.name(L$lik$mangled) else NULL
  cl(rt_laplace, L$kind, if (L$generic) "generic" else L$fam, lik, cl_list(as.name("list"), codes[seq_len(np)]),
     codes[[np + 1L]], as.name(L$kf$mangled), codes[[np + 2L]], cl_list(as.name("list"), unname(ctl)), CTX, L$name)
}

# ---- likelihoods of the built in forms -----------------------------------------------------------------

# The group counts loop over the first length(theta) outcomes for the
# Poisson and Bernoulli forms and over all outcomes for the negative
# binomial form, as Stan Math 5.2 does.
laplace_groups <- function(y, idx, n, nobs, name) {
  if (length(idx) < nobs || length(y) < nobs) stop_shape(sprintf("%s: the outcomes and their group indices must have at least %d elements", name, nobs))
  if (any(idx[seq_len(nobs)] < 1L | idx[seq_len(nobs)] > n)) stop_domain(sprintf("%s: group indices must be between 1 and %d", name, n))
  counts <- numeric(n); ng <- numeric(n)
  for (i in seq_len(nobs)) {
    counts[idx[i]] <- counts[idx[i]] + y[i]
    ng[idx[i]] <- ng[idx[i]] + 1
  }
  list(counts = counts, n = ng)
}

laplace_builtin_lik <- function(fam, pre, name) {
  y <- as.integer(val(pre[[1L]])); idx <- as.integer(val(pre[[2L]]))
  if (fam == "neg_binomial_2_log") {
    # CmdStan 2.40 does not compile this function (Stan Math multiplies the
    # group sizes and eta as matrices); eta here is one dispersion for all
    # groups or one per group
    eta <- pre[[3L]]; mu <- pre[[4L]]
    lik <- function(theta, a) {
      n <- length(val(theta))
      gr <- laplace_groups(y, idx, n, length(y), name)
      eta <- a[[1L]]; m <- a[[2L]]
      ne <- length(val(eta))
      if (ne != 1L && ne != n) stop_shape(sprintf("%s: eta must have size 1 or %d, but has size %d", name, n, ne))
      eobs <- if (ne == 1L) eta else ad_gather(eta, idx)
      to <- theta + m
      le <- log(eta)
      lse <- ad_log_sum_exp_pair(to, le)
      sum(ad_lchoose(y + eobs - 1, y)) + sum(gr$counts * (to - lse) + (gr$n * eta) * (le - lse))
    }
    return(list(lik = lik, args = list(eta, mu)))
  }
  m <- pre[[3L]]
  lik <- if (fam == "poisson_log") function(theta, a) {
    n <- length(val(theta))
    gr <- laplace_groups(y, idx, n, n, name)
    to <- theta + a[[1L]]
    -sum(lgamma(gr$counts + 1)) + sum(to * gr$counts) - sum(gr$n * exp(to))
  } else function(theta, a) {
    n <- length(val(theta))
    gr <- laplace_groups(y, idx, n, n, name)
    to <- theta + a[[1L]]
    sum(to * gr$counts - gr$n * ad_log1p_exp(to))
  }
  list(lik = lik, args = list(m))
}

ad_log_sum_exp_pair <- function(a, b) {
  mx <- pmax(val(a), val(b))
  mx + log(exp(a - mx) + exp(b - mx))
}

# ---- the approximation ------------------------------------------------------------------------------------

rt_laplace <- function(kind, form, likfun, pre, hbs, kfun, kargs, ctl, ctx, name) {
  if (form == "generic") {
    largs <- unclass(pre[[1L]])
    likf <- function(theta, a) do.call(likfun, c(list(theta), a, list(ctx)))
  } else {
    bl <- laplace_builtin_lik(form, pre, name)
    largs <- bl$args
    likf <- function(theta, a) bl$lik(theta, a)
  }
  kargs <- unclass(kargs)
  Kf <- function(a) do.call(kfun, c(a, list(ctx)))
  if (length(ctl) == 1L) ctl <- unclass(ctl[[1L]])
  hbs <- as.integer(val(hbs))
  fl <- hof_flatten(largs); fk <- hof_flatten(kargs)
  pl <- fl$rebuild(fl$theta); pk <- fk$rebuild(fk$theta)
  K <- as.matrix(val(Kf(pk)))
  n <- nrow(K)
  if (ncol(K) != n) stop_shape(sprintf("%s: the covariance matrix must be square, but has %d rows and %d columns", name, nrow(K), ncol(K)))
  if (!all(is.finite(K))) stop_domain(sprintf("%s: the covariance matrix must be finite", name))
  theta0 <- if (length(ctl)) as.numeric(val(ctl[[1L]])) else numeric(n)
  tolv <- if (length(ctl)) val(ctl[[2L]]) else 1.49012e-08
  maxit <- if (length(ctl)) val(ctl[[3L]]) else 500
  maxls <- if (length(ctl)) val(ctl[[5L]]) else 1000
  if (length(theta0) != n) stop_shape(sprintf("%s: initial guess has size %d, but the covariance matrix has %d rows", name, length(theta0), n))
  if (!all(is.finite(theta0))) stop_domain(sprintf("%s: initial guess must be finite", name))
  if (hbs < 1L || hbs > n || n %% hbs != 0L) stop_domain(sprintf("%s: the Hessian block size (%d) must divide the size of the latent vector (%d)", name, hbs, n))
  if (!(tolv > 0)) stop_domain(sprintf("%s: tolerance must be positive", name))
  if (!(maxit > 0)) stop_domain(sprintf("%s: max_num_steps must be positive", name))

  fval <- function(th, p = pl) as.numeric(val(likf(th, p)))
  fgrad <- function(th, p = pl) ad_value_grad(function(z) likf(z, p), th)$grad
  block <- (seq_len(n) - 1L) %/% hbs
  hstep <- function(x) (.Machine$double.eps)^(1 / 3) * pmax(1, abs(x))
  # block diagonal Hessian from colored central differences of the gradient
  fhess <- function(th, p = pl) {
    H <- matrix(0, n, n)
    h <- hstep(th)
    for (k in seq_len(hbs)) {
      cols <- which((seq_len(n) - 1L) %% hbs == k - 1L)
      d <- numeric(n); d[cols] <- h[cols]
      D <- fgrad(th + d, p) - fgrad(th - d, p)
      for (j in cols) {
        rows <- which(block == block[j])
        H[rows, j] <- D[rows] / (2 * h[j])
      }
    }
    (H + t(H)) / 2
  }

  # Newton iterations on a with theta = K a, halving the step while the
  # objective decreases
  a <- if (any(theta0 != 0)) solve(K, theta0) else numeric(n)
  th <- as.vector(K %*% a)
  obj <- fval(th) - sum(th * a) / 2
  if (!is.finite(obj)) stop_domain(sprintf("%s: the log likelihood is not finite at the initial guess", name))
  it <- 0L
  repeat {
    it <- it + 1L
    if (it > maxit) stop_domain(sprintf("%s: max number of iterations: %s exceeded.", name, format_scalar(maxit)))
    g <- fgrad(th); W <- -fhess(th)
    anew <- solve(diag(n) + W %*% K, W %*% th + g)
    anew <- as.vector(anew)
    da <- anew - a
    ls <- 0L
    repeat {
      tn <- as.vector(K %*% (a + da))
      on <- fval(tn) - sum(tn * (a + da)) / 2
      if (is.finite(on) && on >= obj - 1e-14 * abs(obj)) break
      ls <- ls + 1L
      if (ls > maxls) break
      da <- da / 2
    }
    a <- a + da; th <- tn
    done <- abs(on - obj) < tolv
    obj <- on
    if (done) break
  }
  # one more Newton step at the converged point keeps a = grad f accurate
  g <- fgrad(th); H <- fhess(th); W <- -H
  a <- as.vector(solve(diag(n) + W %*% K, W %*% th + g)); th <- as.vector(K %*% a)
  g <- fgrad(th); H <- fhess(th); W <- -H
  IKW <- diag(n) + K %*% W
  ld <- determinant(IKW, logarithm = TRUE)
  if (ld$sign <= 0) stop_domain(sprintf("%s: the Laplace approximation needs a log concave likelihood at the mode", name))
  C <- solve(IKW, K); C <- (C + t(C)) / 2

  if (kind != "marginal") {
    L <- tryCatch(chol(C), error = function(e) stop_domain(sprintf("%s: the covariance of the approximation is not positive definite", name)))
    if (kind == "solve") return(rt_tuple(list(th, t(L))))
    return(th + as.vector(crossprod(L, stats::rnorm(n))))
  }
  value <- fval(th) - sum(th * a) / 2 - as.numeric(ld$modulus) / 2
  if (!length(fl$theta) && !length(fk$theta)) return(value)

  # derivative pieces
  npl <- length(fl$theta)
  Tfun <- function(th, p) { Hb <- fhess(th, p); vapply(0:(n %/% hbs - 1L), function(b) { r <- which(block == b); sum(C[r, r] * Hb[r, r]) }, 0) }
  h2 <- function(x) (.Machine$double.eps)^(1 / 6) * pmax(1, abs(x))
  s2 <- numeric(n)
  {
    hh <- h2(th)
    for (k in seq_len(hbs)) {
      cols <- which((seq_len(n) - 1L) %% hbs == k - 1L)
      d <- numeric(n); d[cols] <- hh[cols]
      dT <- Tfun(th + d, pl) - Tfun(th - d, pl)
      s2[cols] <- dT[block[cols] + 1L] / (2 * hh[cols]) / 2
    }
  }
  gl <- NULL
  if (npl) {
    gz <- ad_value_grad(function(z) likf(ad_gather(z, seq_len(n)), fl$rebuild(ad_gather(z, n + seq_len(npl)))), c(th, fl$theta))$grad
    dfp <- gz[n + seq_len(npl)]
    # relative steps keep positive arguments (a dispersion) positive
    hp <- (.Machine$double.eps)^(1 / 6) * ifelse(fl$theta == 0, 1, abs(fl$theta))
    Htp <- matrix(0, n, npl); dTp <- numeric(npl)
    for (l in seq_len(npl)) {
      up <- fl$theta; dn <- fl$theta
      up[l] <- up[l] + hp[l]; dn[l] <- dn[l] - hp[l]
      pu <- fl$rebuild(up); pd <- fl$rebuild(dn)
      Htp[, l] <- (fgrad(th, pu) - fgrad(th, pd)) / (2 * hp[l])
      dTp[l] <- (sum(Tfun(th, pu)) - sum(Tfun(th, pd))) / (2 * hp[l])
    }
    CH <- C %*% Htp
  }
  if (npl) gl <- dfp + dTp / 2 + as.vector(crossprod(CH, s2))
  gk <- NULL
  if (length(fk$theta)) {
    u <- solve(t(IKW), s2)
    M <- outer(a, a) / 2 - t(W %*% solve(IKW)) / 2 + outer(u, a)
    gk <- ad_value_grad(function(z) sum(M * Kf(fk$rebuild(z))), fk$theta)$grad
  }
  ad_node(value, c(fl$parents, fk$parents), function(g) {
    c(if (npl) hof_split(g * gl, fl$parents), if (length(fk$theta)) hof_split(g * gk, fk$parents))
  })
}

def("generate_laplace_options", "int", function(types) laplace_opts_ty(), function(n) rt_tuple(list(numeric(n), 1.49012e-08, 500L, 1L, 1000L, 1L)))
def("generate_laplace_options", "vector", function(types) laplace_opts_ty(), function(v) rt_tuple(list(as.numeric(val(v)), 1.49012e-08, 500L, 1L, 1000L, 1L)))
