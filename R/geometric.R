# Samplers that adapt to local geometry or use non standard dynamics:
# Riemannian manifold HMC and simplified manifold MALA (Girolami and
# Calderhead 2011) with the SoftAbs metric (Betancourt 2013), and
# microcanonical Langevin Monte Carlo (Robnik, De Luca, Silverstein and
# Seljak 2023).

# Hessian of the potential V = -log p by central differences of the
# gradient, symmetrized
num_hessian_V <- function(gradV, q, h = 1e-5) {
  n <- length(q)
  H <- matrix(0, n, n)
  for (i in seq_len(n)) {
    e <- replace(numeric(n), i, h)
    H[, i] <- (gradV(q + e) - gradV(q - e)) / (2 * h)
  }
  0.5 * (H + t(H))
}

# SoftAbs metric: eigenvalues lambda of the Hessian mapped to
# lambda * coth(alpha * lambda), a smooth positive version of |lambda|
softabs <- function(H, alpha) {
  e <- eigen(H, symmetric = TRUE)
  lam <- e$values
  sl <- ifelse(abs(alpha * lam) < 1e-8, 1 / alpha, lam / tanh(alpha * lam))
  list(G = e$vectors %*% (sl * t(e$vectors)), Q = e$vectors, l = sl)
}

#' Riemannian manifold samplers
#'
#' Hamiltonian Monte Carlo and the Metropolis adjusted Langevin algorithm
#' with a position dependent metric, the SoftAbs transform of the Hessian of
#' the negative log density (Betancourt 2013), so proposals follow the local
#' curvature of funnels and strongly varying scales.
#'
#' `"smmala"` (simplified manifold MALA) proposes from
#' \eqn{N(q + \epsilon^2 G^{-1} \nabla \log p / 2, \epsilon^2 G^{-1})} with the
#' metric at the current point and corrects with the exact reverse proposal
#' density, so it targets the posterior exactly. `"rmhmc"` integrates the
#' non separable Hamiltonian with the implicit generalized leapfrog; the
#' metric derivatives come from finite differences of Hessians, so the
#' integrator is only approximately volume preserving and a small bias is
#' possible. Hessians are computed by finite differences of gradients, so
#' each iteration costs on the order of `dim` (smMALA) or `dim^2` (RMHMC)
#' gradients: these methods suit models with few parameters.
#'
#' @param model A `bayescomp_model`.
#' @param method `"smmala"` or `"rmhmc"`.
#' @param chains,iter_warmup,iter_sampling,init,seed,cores,refresh As in
#'   [bayescomp_sample()].
#' @param step_size Initial step size, adapted by dual averaging in warmup.
#' @param int_time Integration time for RMHMC.
#' @param alpha SoftAbs sharpness.
#' @param adapt_delta Target acceptance (default 0.574 for smMALA, 0.8 for
#'   RMHMC).
#' @param fixed_point_steps Fixed point iterations of the implicit
#'   generalized leapfrog.
#' @return A `bayescomp_fit` with method `"rmhmc"`.
#' @references Girolami, M. and Calderhead, B. (2011). Riemann manifold
#'   Langevin and Hamiltonian Monte Carlo methods. *Journal of the Royal
#'   Statistical Society B* 73, 123 to 214. Betancourt, M. (2013). A general
#'   metric for Riemannian manifold Hamiltonian Monte Carlo. *Geometric
#'   Science of Information*, 327 to 334.
#' @export
bayescomp_riemannian <- function(model, method = c("smmala", "rmhmc"), chains = 4L, iter_warmup = 500L,
                                 iter_sampling = 500L, step_size = 0.5, int_time = 1.5, alpha = 1e6,
                                 adapt_delta = NULL, fixed_point_steps = 20L, init = 2, seed = NULL,
                                 cores = getOption("mc.cores", 1L), refresh = 0L) {
  check_model(model)
  method <- match.arg(method)
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  if (n > 30L) warning("Riemannian samplers cost O(dim^2) gradients per step; this model has ", n, " parameters", call. = FALSE)
  chains <- as_count(chains, "chains", 1L)
  nw <- as_count(iter_warmup, "iter_warmup")
  ns <- as_count(iter_sampling, "iter_sampling", 1L)
  if (is.null(adapt_delta)) adapt_delta <- if (method == "smmala") 0.574 else 0.8
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, chains, seed, cores, function(id) {
    pot <- hmc_potential(model)
    gradV <- function(q) pot(q)$g
    geom <- function(q) {
      r <- pot(q)
      if (!is.finite(r$V)) return(list(V = Inf))
      H <- num_hessian_V(gradV, q)
      if (any(!is.finite(H))) return(list(V = Inf))
      sa <- softabs(H, alpha)
      Ginv <- sa$Q %*% ((1 / sa$l) * t(sa$Q))
      list(V = r$V, g = r$g, G = sa$G, Ginv = 0.5 * (Ginv + t(Ginv)), logdet = sum(log(sa$l)),
           Lc = sa$Q %*% diag(1 / sqrt(sa$l), n))
    }
    q <- find_init(model, init, id)
    gq <- geom(q)
    eps <- step_size
    da <- dual_avg_new(adapt_delta, mu = log(10 * step_size))
    qs <- matrix(NA_real_, ns, n)
    dg <- matrix(0, ns, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    for (it in seq_len(nw + ns)) {
      if (it == nw + 1L && nw > 0L) eps <- exp(da$x_bar)
      if (method == "smmala") {
        r <- smmala_step(q, gq, eps, geom)
        nl <- 1L
      } else {
        L <- max(1L, min(100L, as.integer(int_time / eps)))
        r <- rmhmc_step(q, gq, eps, L, geom, gradV, fixed_point_steps)
        nl <- L
      }
      q <- r$q; gq <- r$gq
      if (it <= nw) eps <- dual_avg_learn(da, r$acc)
      else {
        i <- it - nw
        qs[i, ] <- q
        dg[i, ] <- c(-gq$V, r$acc, eps, 0, nl, as.numeric(isTRUE(r$div)), r$H %||% NA_real_)
      }
      if (refresh > 0L && it %% refresh == 0L) chain_message(sprintf("Chain %d: Iteration %d / %d", id, it, nw + ns))
    }
    list(q = qs, gq = chain_write_draws(model, qs, FALSE), diag = dg, step_size = eps, inv_metric = gq$Ginv,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  make_fit(model, res, list(algorithm = method), "rmhmc", seed, proc.time()[["elapsed"]] - t0)
}

smmala_step <- function(q, gq, eps, geom) {
  n <- length(q)
  mean_of <- function(x, gx) x - 0.5 * eps^2 * as.vector(gx$Ginv %*% gx$g)
  ldens <- function(y, x, gx) {
    d <- y - mean_of(x, gx)
    -0.5 * sum(d * as.vector(gx$G %*% d)) / eps^2 + 0.5 * gx$logdet
  }
  y <- mean_of(q, gq) + eps * as.vector(gq$Lc %*% stats::rnorm(n))
  gy <- geom(y)
  acc <- 0
  if (is.finite(gy$V)) {
    la <- (gq$V - gy$V) + ldens(q, y, gy) - ldens(y, q, gq)
    acc <- if (is.nan(la)) 0 else min(1, exp(la))
  }
  if (stats::runif(1) < acc) return(list(q = y, gq = gy, acc = acc))
  list(q = q, gq = gq, acc = acc)
}

# dH/dq for H(q, p) = V + logdet(G)/2 + p' G^-1 p / 2, as a function of p
# at a fixed position. The metric derivatives dG/dq_i come from central
# differences of the SoftAbs metric and are computed once per position, so
# the fixed point iterations over p are cheap.
rmhmc_dHdq_fun <- function(q, gq, geom, h = 1e-4) {
  n <- length(q)
  dG <- vector("list", n)
  tr <- numeric(n)
  for (i in seq_len(n)) {
    e <- replace(numeric(n), i, h)
    gp <- geom(q + e); gm <- geom(q - e)
    if (!is.finite(gp$V) || !is.finite(gm$V)) return(NULL)
    dG[[i]] <- (gp$G - gm$G) / (2 * h)
    tr[i] <- sum(gq$Ginv * dG[[i]])
  }
  base <- gq$g + 0.5 * tr
  function(p) {
    a <- as.vector(gq$Ginv %*% p)
    base - 0.5 * vapply(dG, function(D) sum(a * as.vector(D %*% a)), 0)
  }
}

# implicit generalized leapfrog (Girolami and Calderhead 2011, equations 16
# to 18), solved by fixed point iteration. The map is reversible and volume
# preserving only when the implicit equations hold, so a step whose fixed
# point iterations have not converged to `tol` within `fp` iterations is
# rejected (and counted as divergent) instead of entering the Metropolis test.
rmhmc_step <- function(q, gq, eps, L, geom, gradV, fp, tol = 1e-9) {
  conv <- function(new, old) max(abs(new - old)) <= tol * max(1, max(abs(new)))
  n <- length(q)
  Hfun <- function(gx, p) gx$V + 0.5 * gx$logdet + 0.5 * sum(p * as.vector(gx$Ginv %*% p))
  p <- as.vector(t(chol(gq$G)) %*% stats::rnorm(n))
  H0 <- Hfun(gq, p)
  fail <- list(q = q, gq = gq, acc = 0, div = TRUE, H = H0)
  x <- q; gx <- gq
  for (l in seq_len(L)) {
    dH <- rmhmc_dHdq_fun(x, gx, geom)
    if (is.null(dH)) return(fail)
    ph <- p
    ok <- FALSE
    for (k in seq_len(fp)) {
      pn <- p - 0.5 * eps * dH(ph)
      if (any(!is.finite(pn))) return(fail)
      done <- conv(pn, ph)
      ph <- pn
      if (done) { ok <- TRUE; break }
    }
    if (!ok) return(fail)
    v0 <- as.vector(gx$Ginv %*% ph)
    xn <- x + eps * v0
    ok <- FALSE
    for (k in seq_len(fp)) {
      gn <- geom(xn)
      if (!is.finite(gn$V)) return(fail)
      xk <- x + 0.5 * eps * (v0 + as.vector(gn$Ginv %*% ph))
      done <- conv(xk, xn)
      xn <- xk
      if (done) { ok <- TRUE; break }
    }
    if (!ok) return(fail)
    gn <- geom(xn)
    if (!is.finite(gn$V)) return(fail)
    dH2 <- rmhmc_dHdq_fun(xn, gn, geom)
    if (is.null(dH2)) return(fail)
    p <- ph - 0.5 * eps * dH2(ph)
    x <- xn; gx <- gn
    if (any(!is.finite(p))) return(fail)
  }
  H1 <- Hfun(gx, p)
  acc <- if (is.finite(H1)) min(1, exp(H0 - H1)) else 0
  div <- !is.finite(H1) || (H1 - H0) > 1000
  if (stats::runif(1) < acc) return(list(q = x, gq = gx, acc = acc, div = div, H = H1))
  list(q = q, gq = gq, acc = acc, div = div, H = H0)
}

# ---- microcanonical Langevin Monte Carlo ----------------------------------------------

#' Microcanonical Langevin Monte Carlo
#'
#' Unadjusted microcanonical Langevin dynamics (Robnik et al. 2023, 2024):
#' energy sampling Hamiltonian dynamics on a unit speed velocity with partial
#' velocity refreshment, preconditioned by a diagonal scale. During warmup
#' the step size is tuned so the energy error per dimension has variance
#' `energy_var`. At the middle of warmup the diagonal scales are set from
#' the draws of the second quarter; at the end, the momentum decoherence
#' length `L` is set to 0.4 times the step size times the mean number of steps
#' per effective sample of the second half of warmup. Without a Metropolis correction the draws have
#' a small asymptotic bias controlled by the step size.
#'
#' @param model A `bayescomp_model`.
#' @param chains,iter_warmup,iter_sampling,init,seed,cores,refresh As in
#'   [bayescomp_sample()].
#' @param energy_var Target variance of the energy error per dimension.
#' @return A `bayescomp_fit` with method `"mclmc"`; `step_size` and element
#'   `L` per chain.
#' @references Robnik, J., De Luca, G. B., Silverstein, E. and Seljak, U.
#'   (2023). Microcanonical Hamiltonian Monte Carlo. *Journal of Machine
#'   Learning Research* 24, 1 to 34.
#' @export
bayescomp_mclmc <- function(model, chains = 4L, iter_warmup = 1000L, iter_sampling = 1000L,
                            energy_var = 5e-4, init = 2, seed = NULL, cores = getOption("mc.cores", 1L), refresh = 0L) {
  check_model(model)
  n <- model$dim
  if (n < 2L) stop("MCLMC needs at least two parameters", call. = FALSE)
  chains <- as_count(chains, "chains", 1L)
  nw <- as_count(iter_warmup, "iter_warmup", 10L)
  ns <- as_count(iter_sampling, "iter_sampling", 1L)
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, chains, seed, cores, function(id) {
    pot <- hmc_potential(model)
    x <- find_init(model, init, id)
    r <- pot(x)
    sig <- rep(1, n)
    u <- stats::rnorm(n); u <- u / sqrt(sum(u^2))
    eps <- sqrt(n) * 0.25
    L <- sqrt(n)
    est <- est_new(n, FALSE)
    w_sum <- 0; xi_sum <- 0
    # draws of the last warmup stage, in preconditioned coordinates, for L
    s3 <- seq.int(floor(nw / 2) + 1L, nw)
    X3 <- matrix(NA_real_, length(s3), n)
    upd_u <- function(u, g, e) {
      # g: gradient of log p in the preconditioned coordinates
      gn <- sqrt(sum(g^2))
      if (gn == 0) return(list(u = u, dK = 0))
      ev <- g / gn
      ue <- sum(u * ev)
      delta <- e * gn / (n - 1)
      zeta <- exp(-delta)
      nu <- ev * (1 - zeta) * (1 + zeta + ue * (1 - zeta)) + 2 * zeta * u
      list(u = nu / sqrt(sum(nu^2)), dK = (n - 1) * (delta - log(2) + log(1 + ue + (1 - ue) * zeta^2)))
    }
    qs <- matrix(NA_real_, ns, n)
    dg <- matrix(0, ns, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    for (it in seq_len(nw + ns)) {
      # leapfrog in coordinates y = x / sig
      a <- upd_u(u, -r$g * sig, eps / 2)
      xn <- x + eps * sig * a$u
      rn <- pot(xn)
      if (!is.finite(rn$V)) {
        # a step outside the support: resample velocity, and during warmup
        # also shrink the step (the step is fixed while sampling)
        if (it <= nw) eps <- eps / 2
        u <- stats::rnorm(n); u <- u / sqrt(sum(u^2))
        if (it > nw) { i <- it - nw; qs[i, ] <- x; dg[i, ] <- c(-r$V, 0, eps, 0, 1, 1, NA_real_) }
        next
      }
      b <- upd_u(a$u, -rn$g * sig, eps / 2)
      dE <- (rn$V - r$V) + a$dK + b$dK
      x <- xn; r <- rn; u <- b$u
      # partial refreshment
      c1 <- exp(-eps / L)
      u <- c1 * u + sqrt(1 - c1^2) * stats::rnorm(n) / sqrt(n)
      u <- u / sqrt(sum(u^2))
      if (it <= nw) {
        xi <- dE^2 / (n * energy_var) + 1e-8
        w <- exp(-0.5 * (log(xi) / 6 / 1.5)^2)
        xi_sum <- 0.95 * xi_sum + w * xi / eps^6
        w_sum <- 0.95 * w_sum + w
        eps <- (xi_sum / w_sum)^(-1 / 6)
        # precondition at the middle of warmup from the draws of the second
        # quarter, then keep tuning the step size in the new coordinates
        if (it > nw / 4 && it <= nw / 2) est_add(est, x)
        if (it == floor(nw / 2) && est$n > 10) {
          v <- est_value(est)
          sig <- sqrt(pmax(v, 1e-12))
          L <- sqrt(n)
          eps <- eps / sqrt(mean(sig^2))
          xi_sum <- 0; w_sum <- 0
        }
        # third stage (Robnik et al. 2024): L = 0.4 eps times the mean number
        # of steps per effective sample, from the draws of the second half
        if (it > floor(nw / 2)) X3[it - floor(nw / 2), ] <- x / sig
        if (it == nw && nrow(X3) >= 20L) {
          ess <- apply(X3, 2L, function(v) bayescomp_ess_basic(matrix(v, ncol = 1L)))
          ess <- ess[is.finite(ess) & ess > 0]
          if (length(ess)) L <- 0.4 * eps * mean(nrow(X3) / ess)
        }
      } else {
        i <- it - nw
        qs[i, ] <- x
        dg[i, ] <- c(-r$V, NA_real_, eps, 0, 1, 0, dE)
      }
      if (refresh > 0L && it %% refresh == 0L) chain_message(sprintf("Chain %d: Iteration %d / %d", id, it, nw + ns))
    }
    list(q = qs, gq = chain_write_draws(model, qs, FALSE), diag = dg, step_size = eps, inv_metric = sig^2,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0), L = L)
  })
  fit <- make_fit(model, res, list(algorithm = "mclmc"), "mclmc", seed, proc.time()[["elapsed"]] - t0)
  fit$L <- vapply(res, `[[`, 0, "L")
  fit
}
