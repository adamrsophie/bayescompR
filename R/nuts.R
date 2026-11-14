# Hamiltonian Monte Carlo: the No-U-Turn sampler with multinomial trajectory
# sampling and the generalized no-U-turn criterion (Stan's base_nuts), static
# HMC, Euclidean metrics, step size initialization, dual averaging and
# windowed metric adaptation.
#
# The potential energy is V(q) = -log p(q) on the unconstrained space,
# including the Jacobian of the constraining transforms and dropping
# constants (propto). Domain errors and rejections while evaluating V give
# V = Inf, so the transition is rejected and, inside a trajectory, flagged as
# divergent.

hmc_potential <- function(model, print = FALSE) {
  lpg <- model_lpg(model, jacobian = TRUE, propto = TRUE, print = print)
  function(q) {
    r <- tryCatch(lpg(q), bayescomp_fatal_error = function(e) stop(e),
       bayescomp_domain_error = function(e) NULL,
       bayescomp_reject_error = function(e) NULL)
    if (is.null(r) || !is.finite(r$value) || any(!is.finite(r$grad)))
      return(list(V = Inf, g = rep(0, length(q))))
    list(V = -r$value, g = -r$grad)
  }
}

# ---- metrics --------------------------------------------------------------------

make_metric <- function(type, inv) {
  if (type == "diag_e") return(list(type = type, inv = inv))
  U <- chol(inv)
  list(type = type, inv = inv, U = U)
}
metric_sample_p <- function(m, n) {
  u <- stats::rnorm(n)
  if (m$type == "diag_e") u / sqrt(m$inv) else backsolve(m$U, u)
}
metric_dtau <- function(m, p) if (m$type == "diag_e") m$inv * p else as.vector(m$inv %*% p)
metric_tau <- function(m, p) 0.5 * sum(p * metric_dtau(m, p))

# phase space point: q, p, V, g
ps_H <- function(m, z) z$V + metric_tau(m, z$p)

leapfrog <- function(pot, m, z, eps) {
  p <- z$p - 0.5 * eps * z$g
  q <- z$q + eps * metric_dtau(m, p)
  r <- pot(q)
  p <- p - 0.5 * eps * r$g
  list(q = q, p = p, V = r$V, g = r$g)
}

lse2 <- function(a, b) {
  if (a == -Inf) return(b)
  if (b == -Inf) return(a)
  mx <- max(a, b)
  mx + log(exp(a - mx) + exp(b - mx))
}

nuts_criterion <- function(ps_minus, ps_plus, rho) sum(ps_plus * rho) > 0 && sum(ps_minus * rho) > 0

# Recursive tree building. `st` (environment) holds the moving endpoint z,
# the leapfrog count, the summed acceptance probability and the divergence
# flag. Returns the proposal, the sharp momenta and momenta at both ends,
# the summed momentum rho and the log sum of weights of the subtree.
nuts_build_tree <- function(st, depth, sign, H0) {
  if (depth == 0L) {
    st$z <- leapfrog(st$pot, st$metric, st$z, sign * st$eps)
    st$n_leapfrog <- st$n_leapfrog + 1L
    h <- ps_H(st$metric, st$z)
    if (is.nan(h)) h <- Inf
    if ((h - H0) > st$max_deltaH) st$divergent <- TRUE
    st$sum_metro <- st$sum_metro + if (H0 - h > 0) 1 else exp(H0 - h)
    ps <- metric_dtau(st$metric, st$z$p)
    return(list(valid = !st$divergent, propose = st$z, psb = ps, pse = ps, rho = st$z$p,
                pb = st$z$p, pe = st$z$p, lsw = H0 - h))
  }
  init <- nuts_build_tree(st, depth - 1L, sign, H0)
  if (!init$valid) return(list(valid = FALSE))
  fin <- nuts_build_tree(st, depth - 1L, sign, H0)
  if (!fin$valid) return(list(valid = FALSE))
  lsw <- lse2(init$lsw, fin$lsw)
  propose <- init$propose
  if (fin$lsw > lsw || stats::runif(1) < exp(fin$lsw - lsw)) propose <- fin$propose
  rho <- init$rho + fin$rho
  persist <- nuts_criterion(init$psb, fin$pse, rho) &&
    nuts_criterion(init$psb, fin$psb, init$rho + fin$pb) &&
    nuts_criterion(init$pse, fin$pse, fin$rho + init$pe)
  list(valid = persist, propose = propose, psb = init$psb, pse = fin$pse, rho = rho,
       pb = init$pb, pe = fin$pe, lsw = lsw)
}

nuts_transition <- function(st, q) {
  m <- st$metric
  r <- st$pot(q)
  z <- list(q = q, p = metric_sample_p(m, length(q)), V = r$V, g = r$g)
  H0 <- ps_H(m, z)
  z_fwd <- z; z_bck <- z; z_sample <- z
  ps0 <- metric_dtau(m, z$p)
  p_fwd_fwd <- p_fwd_bck <- p_bck_fwd <- p_bck_bck <- z$p
  ps_fwd_fwd <- ps_fwd_bck <- ps_bck_fwd <- ps_bck_bck <- ps0
  rho <- z$p
  lsw <- 0
  st$n_leapfrog <- 0L
  st$sum_metro <- 0
  st$divergent <- FALSE
  depth <- 0L
  while (depth < st$max_depth) {
    if (stats::runif(1) > 0.5) {
      st$z <- z_fwd
      rho_bck <- rho
      p_bck_fwd <- p_fwd_fwd
      ps_bck_fwd <- ps_fwd_fwd
      tr <- nuts_build_tree(st, depth, 1, H0)
      z_fwd <- st$z
      if (!tr$valid) break
      rho_fwd <- tr$rho
      ps_fwd_bck <- tr$psb; ps_fwd_fwd <- tr$pse
      p_fwd_bck <- tr$pb; p_fwd_fwd <- tr$pe
    } else {
      st$z <- z_bck
      rho_fwd <- rho
      p_fwd_bck <- p_bck_bck
      ps_fwd_bck <- ps_bck_bck
      tr <- nuts_build_tree(st, depth, -1, H0)
      z_bck <- st$z
      if (!tr$valid) break
      rho_bck <- tr$rho
      ps_bck_fwd <- tr$psb; ps_bck_bck <- tr$pse
      p_bck_fwd <- tr$pb; p_bck_bck <- tr$pe
    }
    depth <- depth + 1L
    if (tr$lsw > lsw) {
      z_sample <- tr$propose
    } else if (stats::runif(1) < exp(tr$lsw - lsw)) {
      z_sample <- tr$propose
    }
    lsw <- lse2(lsw, tr$lsw)
    rho <- rho_bck + rho_fwd
    persist <- nuts_criterion(ps_bck_bck, ps_fwd_fwd, rho) &&
      nuts_criterion(ps_bck_bck, ps_fwd_bck, rho_bck + p_fwd_bck) &&
      nuts_criterion(ps_bck_fwd, ps_fwd_fwd, rho_fwd + p_bck_fwd)
    if (!persist) break
  }
  list(q = z_sample$q, lp = -z_sample$V, accept_stat = if (st$n_leapfrog > 0L) st$sum_metro / st$n_leapfrog else 0,
       treedepth = depth, n_leapfrog = st$n_leapfrog, divergent = st$divergent,
       energy = ps_H(m, z_sample))
}

# static HMC with integration time int_time (Stan's static_hmc)
hmc_transition <- function(st, q) {
  m <- st$metric
  r <- st$pot(q)
  z <- list(q = q, p = metric_sample_p(m, length(q)), V = r$V, g = r$g)
  H0 <- ps_H(m, z)
  L <- max(1L, as.integer(st$int_time / st$eps))
  z0 <- z
  for (i in seq_len(L)) z <- leapfrog(st$pot, m, z, st$eps)
  h <- ps_H(m, z)
  if (is.nan(h)) h <- Inf
  acc <- min(1, exp(H0 - h))
  if (stats::runif(1) > acc) z <- z0
  list(q = z$q, lp = -z$V, accept_stat = acc, treedepth = 0L, n_leapfrog = L,
       divergent = (h - H0) > st$max_deltaH, energy = ps_H(m, z))
}

# Stan's step size initialization heuristic (base_hmc::init_stepsize)
init_stepsize <- function(st, q) {
  eps <- st$eps
  if (eps == 0 || eps > 1e7 || is.nan(eps)) return(eps)
  m <- st$metric
  r <- st$pot(q)
  one <- function(eps) {
    z <- list(q = q, p = metric_sample_p(m, length(q)), V = r$V, g = r$g)
    H0 <- ps_H(m, z)
    z <- leapfrog(st$pot, m, z, eps)
    h <- ps_H(m, z)
    if (is.nan(h)) h <- Inf
    H0 - h
  }
  direction <- if (one(eps) > log(0.8)) 1 else -1
  repeat {
    dH <- one(eps)
    if (direction == 1 && !(dH > log(0.8))) break
    if (direction == -1 && !(dH < log(0.8))) break
    eps <- if (direction == 1) 2 * eps else 0.5 * eps
    if (eps > 1e7) stop("Posterior is improper. Please check your model.", call. = FALSE)
    if (eps == 0) stop("No acceptably small step size could be found. Perhaps the posterior is not continuous?", call. = FALSE)
  }
  eps
}

# ---- adaptation --------------------------------------------------------------------

dual_avg_new <- function(delta, gamma = 0.05, kappa = 0.75, t0 = 10, mu = log(10)) {
  da <- new.env(parent = emptyenv())
  da$delta <- delta; da$gamma <- gamma; da$kappa <- kappa; da$t0 <- t0; da$mu <- mu
  da$counter <- 0; da$s_bar <- 0; da$x_bar <- 0
  da
}
dual_avg_restart <- function(da) { da$counter <- 0; da$s_bar <- 0; da$x_bar <- 0 }
dual_avg_learn <- function(da, adapt_stat) {
  da$counter <- da$counter + 1
  adapt_stat <- min(1, adapt_stat)
  eta <- 1 / (da$counter + da$t0)
  da$s_bar <- (1 - eta) * da$s_bar + eta * (da$delta - adapt_stat)
  x <- da$mu - da$s_bar * sqrt(da$counter) / da$gamma
  x_eta <- da$counter^(-da$kappa)
  da$x_bar <- (1 - x_eta) * da$x_bar + x_eta * x
  exp(x)
}

windows_new <- function(num_warmup, init_buffer = 75L, term_buffer = 50L, base_window = 25L) {
  w <- new.env(parent = emptyenv())
  w$num_warmup <- num_warmup
  if (num_warmup < 20L) {
    w$enabled <- FALSE
  } else {
    w$enabled <- TRUE
    if (init_buffer + base_window + term_buffer > num_warmup) {
      init_buffer <- as.integer(0.15 * num_warmup)
      term_buffer <- as.integer(0.1 * num_warmup)
      base_window <- num_warmup - (init_buffer + term_buffer)
    }
  }
  w$init_buffer <- init_buffer; w$term_buffer <- term_buffer; w$base_window <- base_window
  w$counter <- 0L
  w$window_size <- base_window
  w$next_window <- init_buffer + base_window - 1L
  w
}
win_in_window <- function(w) {
  w$enabled && w$counter >= w$init_buffer && w$counter < w$num_warmup - w$term_buffer && w$counter != w$num_warmup
}
win_end_window <- function(w) w$enabled && w$counter == w$next_window && w$counter != w$num_warmup
win_compute_next <- function(w) {
  last <- w$num_warmup - w$term_buffer - 1L
  if (w$next_window == last) return(invisible())
  w$window_size <- 2L * w$window_size
  w$next_window <- w$counter + w$window_size
  if (w$next_window != last) {
    if (w$next_window + 2L * w$window_size >= w$num_warmup - w$term_buffer) w$next_window <- last
  }
}

# Welford estimators of the variance or covariance
est_new <- function(n, dense) {
  e <- new.env(parent = emptyenv())
  e$n <- 0; e$mean <- numeric(n)
  e$m2 <- if (dense) matrix(0, n, n) else numeric(n)
  e$dense <- dense
  e
}
est_add <- function(e, q) {
  e$n <- e$n + 1
  d <- q - e$mean
  e$mean <- e$mean + d / e$n
  if (e$dense) e$m2 <- e$m2 + outer(q - e$mean, d) else e$m2 <- e$m2 + (q - e$mean) * d
}
est_value <- function(e) {
  n <- e$n
  v <- e$m2 / (n - 1)
  if (e$dense) (n / (n + 5)) * v + 1e-3 * (5 / (n + 5)) * diag(nrow(v))
  else (n / (n + 5)) * v + 1e-3 * (5 / (n + 5))
}

# ---- one chain ----------------------------------------------------------------------

run_hmc_chain <- function(model, q0, cfg, chain_id) {
  n <- length(q0)
  lp_full <- lp_value_fun(model, print = cfg$print)
  # Metropolis and the slice samplers need only log density values
  pot <- if (cfg$algorithm %in% c("metropolis", "slice", "elliptical_slice")) {
    function(q) list(V = -lp_full(q), g = numeric(n))
  } else hmc_potential(model, print = cfg$print)
  dense <- cfg$metric == "dense_e"
  inv <- cfg$inv_metric %||% (if (dense) diag(n) else rep(1, n))
  st <- new.env(parent = emptyenv())
  st$pot <- pot
  st$metric <- make_metric(cfg$metric, inv)
  st$eps <- cfg$step_size
  st$max_depth <- cfg$max_treedepth
  st$max_deltaH <- 1000
  st$int_time <- cfg$int_time
  lp_off <- lp_full(q0) + pot(q0)$V
  if (!is.finite(lp_off)) lp_off <- 0
  st$lp <- function(q) lp_full(q) - lp_off
  st$lp_cache <- lp_cache_fun(st)
  st$block <- isTRUE(cfg$block)
  st$max_steps <- 32L
  trans <- switch(cfg$algorithm, nuts = nuts_transition, hmc = hmc_transition, mala = mala_transition,
                  metropolis = metropolis_transition, slice = slice_transition,
                  elliptical_slice = elliptical_transition)
  st$center <- q0
  hamiltonian <- cfg$algorithm %in% c("nuts", "hmc")
  nw <- cfg$iter_warmup
  ns <- cfg$iter_sampling
  adapt <- cfg$adapt_engaged && nw > 0L
  tune_step <- adapt && !(cfg$algorithm %in% c("slice", "elliptical_slice"))
  da <- dual_avg_new(cfg$adapt_delta, cfg$adapt_gamma, cfg$adapt_kappa, cfg$adapt_t0, mu = log(10 * cfg$step_size))
  w <- windows_new(nw, cfg$adapt_init_buffer, cfg$adapt_term_buffer, cfg$adapt_window)
  est <- est_new(n, dense)
  if (adapt && hamiltonian) st$eps <- init_stepsize(st, q0)
  keep_w <- if (cfg$save_warmup) seq_len(nw) else integer(0)
  n_keep <- length(keep_w) + ns
  thin_idx <- c(keep_w[(seq_along(keep_w) - 1L) %% cfg$thin == 0L],
                nw + seq_len(ns)[(seq_len(ns) - 1L) %% cfg$thin == 0L])
  qs <- matrix(NA_real_, length(thin_idx), n)
  diag_out <- matrix(NA_real_, length(thin_idx), 7L,
                     dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
  q <- q0
  row <- 0L
  total <- nw + ns
  t_start <- proc.time()[["elapsed"]]
  t_warm <- NA_real_
  for (it in seq_len(total)) {
    if (it == nw + 1L) {
      t_warm <- proc.time()[["elapsed"]] - t_start
      if (tune_step) st$eps <- exp(da$x_bar)
    }
    eps_used <- st$eps
    # elliptical slice sampling relies on the Gaussian fitted during warmup;
    # warmup itself uses scale free slice sampling so that the fit does not
    # shrink along with the moves
    s <- if (cfg$algorithm == "elliptical_slice" && it <= nw) slice_transition(st, q) else trans(st, q)
    q <- s$q
    if (it <= nw && adapt) {
      if (tune_step) st$eps <- dual_avg_learn(da, s$accept_stat)
      if (cfg$adapt_metric) {
        if (win_in_window(w)) est_add(est, q)
        if (win_end_window(w)) {
          win_compute_next(w)
          st$metric <- make_metric(cfg$metric, est_value(est))
          st$center <- est$mean
          est <- est_new(n, dense)
          w$counter <- w$counter + 1L
          if (hamiltonian) st$eps <- init_stepsize(st, q)
          if (tune_step) {
            da$mu <- log(10 * st$eps)
            dual_avg_restart(da)
          }
        } else {
          w$counter <- w$counter + 1L
        }
      }
    }
    if (it %in% thin_idx) {
      row <- row + 1L
      qs[row, ] <- q
      diag_out[row, ] <- c(s$lp, s$accept_stat, eps_used, s$treedepth, s$n_leapfrog, as.numeric(s$divergent), s$energy)
    }
    if (cfg$refresh > 0L && (it %% cfg$refresh == 0L || it == 1L || it == total)) {
      message(sprintf("Chain %d: Iteration: %*d / %d [%3.0f%%]  (%s)", chain_id, nchar(total), it, total,
                      100 * it / total, if (it <= nw) "Warmup" else "Sampling"))
    }
  }
  if (nw == 0L) t_warm <- 0
  t_total <- proc.time()[["elapsed"]] - t_start
  list(q = qs, diag = diag_out, step_size = st$eps, inv_metric = st$metric$inv,
       n_warmup_saved = sum(thin_idx <= nw), time = c(warmup = t_warm, sampling = t_total - t_warm))
}
