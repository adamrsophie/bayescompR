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
  # rejected proposals counted by reason (CmdStan prints each one);
  # attr(pot, "rejections")() returns them
  rej <- new.env(parent = emptyenv())
  note <- function(msg) assign(msg, (rej[[msg]] %||% 0L) + 1L, envir = rej)
  pot <- function(q) {
    r <- tryCatch(lpg(q), bayescomp_fatal_error = function(e) stop(e),
       bayescomp_domain_error = function(e) { note(conditionMessage(e)); NULL },
       bayescomp_reject_error = function(e) { note(conditionMessage(e)); NULL })
    if (is.null(r)) return(list(V = Inf, g = rep(0, length(q))))
    if (!is.finite(r$value) || any(!is.finite(r$grad))) {
      note(if (!is.finite(r$value)) "the log density is not finite" else "the gradient is not finite")
      return(list(V = Inf, g = rep(0, length(q))))
    }
    list(V = -r$value, g = -r$grad)
  }
  attr(pot, "rejections") <- function() {
    nm <- ls(rej, all.names = TRUE)
    stats::setNames(vapply(nm, function(n) rej[[n]], 0L), nm)
  }
  pot
}

# ---- metrics --------------------------------------------------------------------

make_metric <- function(type, inv) {
  if (type == "diag_e") return(list(type = type, inv = inv))
  U <- chol(inv)
  list(type = type, inv = inv, U = U)
}
metric_sample_p <- function(m, n) {
  u <- stats::rnorm(n)
  # the second chain of an antithetic pair draws the negated momenta
  if (isTRUE(m$neg)) u <- -u
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
    q_prev <- st$z$q
    st$z <- if (isTRUE(st$walnuts)) walnuts_step(st, st$z, sign * st$eps)
            else leapfrog(st$pot, st$metric, st$z, sign * st$eps)
    if (!isTRUE(st$walnuts)) st$n_leapfrog <- st$n_leapfrog + 1L
    # the inverse metric times the momentum serves both the kinetic energy
    # and the turning criterion; a dense metric makes it O(dim^2)
    ps <- metric_dtau(st$metric, st$z$p)
    h <- st$z$V + 0.5 * sum(st$z$p * ps)
    if (is.nan(h)) h <- Inf
    if ((h - H0) > st$max_deltaH) {
      # the last point of the trajectory before the energy blew up
      st$div_q <- q_prev
      st$divergent <- TRUE
    }
    st$sum_metro <- st$sum_metro + if (H0 - h > 0) 1 else exp(H0 - h)
    st$n_leaf <- st$n_leaf + 1L
    return(list(valid = !st$divergent, propose = st$z, psb = ps, pse = ps, rho = st$z$p,
                pb = st$z$p, pe = st$z$p, lsw = H0 - h + (st$z$lr %||% 0)))
  }
  init <- nuts_build_tree(st, depth - 1L, sign, H0)
  if (!init$valid) return(list(valid = FALSE))
  fin <- nuts_build_tree(st, depth - 1L, sign, H0)
  if (!fin$valid) return(list(valid = FALSE))
  lsw <- lse2(init$lsw, fin$lsw)
  propose <- init$propose
  # a WALNUTS subtree can have zero weight throughout
  if (fin$lsw > -Inf && (fin$lsw > lsw || stats::runif(1) < exp(fin$lsw - lsw))) propose <- fin$propose
  rho <- init$rho + fin$rho
  persist <- nuts_criterion(init$psb, fin$pse, rho) &&
    nuts_criterion(init$psb, fin$psb, init$rho + fin$pb) &&
    nuts_criterion(init$pse, fin$pse, fin$rho + init$pe)
  list(valid = persist, propose = propose, psb = init$psb, pse = fin$pse, rho = rho,
       pb = init$pb, pe = fin$pe, lsw = lsw)
}

nuts_transition <- function(st, q) {
  m <- st$metric
  # the state a transition ends in is where the next one starts: reuse its
  # potential and gradient instead of evaluating them again
  r <- if (identical(st$last_q, q)) st$last_r else st$pot(q)
  z <- list(q = q, p = metric_sample_p(m, length(q)), V = r$V, g = r$g)
  H0 <- ps_H(m, z)
  z_fwd <- z; z_bck <- z; z_sample <- z
  ps0 <- metric_dtau(m, z$p)
  p_fwd_fwd <- p_fwd_bck <- p_bck_fwd <- p_bck_bck <- z$p
  ps_fwd_fwd <- ps_fwd_bck <- ps_bck_fwd <- ps_bck_bck <- ps0
  rho <- z$p
  lsw <- 0
  st$n_leapfrog <- 0L
  st$n_leaf <- 0L
  st$n_coarse <- 0L
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
  st$last_q <- z_sample$q
  st$last_r <- list(V = z_sample$V, g = z_sample$g)
  # WALNUTS tunes the macro step so that a share adapt_delta of macro steps
  # needs no refinement (Bou-Rabee et al. 2025 aim at about 80%); the energy
  # error is controlled by the micro steps instead
  list(q = z_sample$q, lp = -z_sample$V,
       accept_stat = if (st$n_leaf == 0L) 0 else if (isTRUE(st$walnuts)) st$n_coarse / st$n_leaf else st$sum_metro / st$n_leaf,
       treedepth = depth, n_leapfrog = st$n_leapfrog, divergent = st$divergent,
       energy = ps_H(m, z_sample), p = z_sample$p, g = z_sample$g, div_q = if (st$divergent) st$div_q)
}

# static HMC with integration time int_time (Stan's static_hmc)
hmc_transition <- function(st, q) {
  m <- st$metric
  r <- if (identical(st$last_q, q)) st$last_r else st$pot(q)
  z <- list(q = q, p = metric_sample_p(m, length(q)), V = r$V, g = r$g)
  H0 <- ps_H(m, z)
  L <- max(1L, as.integer(st$int_time / st$eps))
  z0 <- z
  for (i in seq_len(L)) z <- leapfrog(st$pot, m, z, st$eps)
  h <- ps_H(m, z)
  if (is.nan(h)) h <- Inf
  acc <- min(1, exp(H0 - h))
  if (stats::runif(1) > acc) z <- z0
  st$last_q <- z$q
  st$last_r <- list(V = z$V, g = z$g)
  list(q = z$q, lp = -z$V, accept_stat = acc, treedepth = 0L, n_leapfrog = L,
       divergent = (h - H0) > st$max_deltaH, energy = ps_H(m, z), p = z$p, g = z$g)
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
# diagonal inverse metric sqrt(var(q) / var(grad)), the minimizer of the
# Fisher divergence between the draws and a Gaussian with that scale,
# regularized as est_value()
est_fisher <- function(e, eg) {
  k <- e$n
  vq <- e$m2 / (k - 1); vg <- eg$m2 / (k - 1)
  v <- ifelse(vg > 0 & is.finite(vg), sqrt(vq / vg), vq)
  (k / (k + 5)) * v + 1e-3 * (5 / (k + 5))
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
  # Metropolis and the slice samplers need only log density values, of the
  # same (proportional) target as HMC
  lp_full <- lp_value_fun(model, print = cfg$print, propto = TRUE)
  pot <- if (cfg$algorithm %in% c("metropolis", "slice", "elliptical_slice")) {
    function(q) list(V = -lp_full(q), g = numeric(n))
  } else hmc_potential(model, print = cfg$print)
  dense <- cfg$metric == "dense_e"
  inv <- cfg$inv_metric %||% (if (dense) diag(n) else rep(1, n))
  st <- new.env(parent = emptyenv())
  st$pot <- pot
  st$metric <- make_metric(cfg$metric, inv)
  st$metric$neg <- isTRUE(cfg$negate_momentum)
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
  st$walnuts <- cfg$algorithm == "walnuts"
  trans <- switch(cfg$algorithm, nuts = nuts_transition, walnuts = nuts_transition, hmc = hmc_transition, gist = gist_transition, automala = automala_transition, malt = malt_transition, aaps = aaps_transition, mala = mala_transition, barker = barker_transition,
                  metropolis = metropolis_transition, slice = slice_transition,
                  elliptical_slice = elliptical_transition)
  st$center <- cfg$center %||% q0
  hamiltonian <- cfg$algorithm %in% c("nuts", "walnuts", "hmc", "gist", "malt", "aaps")
  nw <- cfg$iter_warmup
  ns <- cfg$iter_sampling
  adapt <- cfg$adapt_engaged && nw > 0L
  tune_step <- adapt && !(cfg$algorithm %in% c("slice", "elliptical_slice", "automala"))
  da <- dual_avg_new(cfg$adapt_delta, cfg$adapt_gamma, cfg$adapt_kappa, cfg$adapt_t0, mu = log(10 * cfg$step_size))
  w <- windows_new(nw, cfg$adapt_init_buffer, cfg$adapt_term_buffer, cfg$adapt_window)
  est <- est_new(n, dense)
  fisher <- isTRUE(cfg$fisher)
  radial <- isTRUE(cfg$radial)
  rad_sd <- cfg$radial_sd %||% (1 / sqrt(n))
  if (fisher) estg <- est_new(n, FALSE)
  if (adapt && hamiltonian) st$eps <- init_stepsize(st, q0)
  off <- as.integer(cfg$sampling_offset %||% 0L)
  keep_w <- if (cfg$save_warmup) seq_len(nw) else integer(0)
  n_keep <- length(keep_w) + ns
  thin_idx <- c(keep_w[(seq_along(keep_w) - 1L) %% cfg$thin == 0L],
                nw + seq_len(ns)[(off + seq_len(ns) - 1L) %% cfg$thin == 0L])
  keep_it <- logical(nw + ns); keep_it[thin_idx] <- TRUE
  qs <- matrix(NA_real_, length(thin_idx), n)
  diag_out <- matrix(NA_real_, length(thin_idx), 7L,
                     dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
  keep_pg <- isTRUE(cfg$diagnostic) && hamiltonian
  if (keep_pg) dp <- dg <- matrix(NA_real_, length(thin_idx), n)
  q <- q0
  row <- 0L
  divs <- list()
  total <- nw + ns
  t_start <- proc.time()[["elapsed"]]
  t_warm <- NA_real_
  # an interrupt (Ctrl-C) ends the chain with the draws kept so far
  n_done <- 0L
  stopped <- tryCatch({ for (it in seq_len(total)) {
    if (it == nw + 1L) {
      t_warm <- proc.time()[["elapsed"]] - t_start
      if (tune_step) st$eps <- exp(da$x_bar)
    }
    eps_nominal <- st$eps
    if (cfg$step_size_jitter > 0) st$eps <- eps_nominal * (1 + cfg$step_size_jitter * (2 * stats::runif(1) - 1))
    eps_used <- st$eps
    # elliptical slice sampling relies on the Gaussian fitted during warmup;
    # warmup itself uses scale free slice sampling so that the fit does not
    # shrink along with the moves
    s <- if (cfg$algorithm == "elliptical_slice" && it <= nw) slice_transition(st, q) else trans(st, q)
    q <- s$q
    if (radial) {
      # radial move (Ostmeyer 2024): scale the distance to the center by
      # exp(g); the Jacobian of the scaling in n dimensions is exp(n g)
      g <- stats::rnorm(1, 0, rad_sd)
      qn <- st$center + (q - st$center) * exp(g)
      lpn <- st$lp(qn)
      la <- lpn - s$lp + n * g
      a <- if (is.finite(la)) min(1, exp(la)) else 0
      if (stats::runif(1) < a) { q <- qn; s$q <- qn; s$lp <- lpn }
      if (it <= nw) rad_sd <- rad_sd * exp((a - 0.5) / sqrt(it))
    }
    st$eps <- eps_nominal
    # autoMALA: the initial step size of the selector follows the geometric
    # mean of the selected step sizes during warmup
    if (!is.null(st$automala_eps) && it <= nw && adapt) {
      st$eps <- exp(log(st$eps) + (log(st$automala_eps) - log(st$eps)) / min(it, 50))
    }
    if (it <= nw && adapt) {
      if (tune_step) st$eps <- dual_avg_learn(da, s$accept_stat)
      if (cfg$adapt_metric) {
        if (win_in_window(w)) { est_add(est, q); if (fisher) est_add(estg, s$g) }
        if (win_end_window(w)) {
          win_compute_next(w)
          # a window of one draw has no variance estimate (as Stan, which
          # needs more than one sample): the metric stays as it is
          if (est$n > 1) st$metric <- make_metric(cfg$metric, if (fisher) est_fisher(est, estg) else est_value(est))
          st$metric$neg <- isTRUE(cfg$negate_momentum)
          st$center <- est$mean
          est <- est_new(n, dense)
          if (fisher) estg <- est_new(n, FALSE)
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
    if (it > nw && !is.null(s$div_q)) divs[[length(divs) + 1L]] <- c(it - nw, s$div_q)
    if (keep_it[it]) {
      row <- row + 1L
      qs[row, ] <- q
      diag_out[row, ] <- c(s$lp, s$accept_stat, eps_used, s$treedepth, s$n_leapfrog, as.numeric(s$divergent), s$energy)
      if (keep_pg) { dp[row, ] <- s$p; dg[row, ] <- s$g }
    }
    if (cfg$refresh > 0L && (it %% cfg$refresh == 0L || it == 1L || it == total)) {
      chain_message(sprintf("Chain %d: Iteration: %*d / %d [%3.0f%%]  (%s)", chain_id, nchar(total), it, total,
                      100 * it / total, if (it <= nw) "Warmup" else "Sampling"))
    }
    n_done <- it
  }; FALSE }, interrupt = function(e) TRUE)
  if (!stopped && ns == 0L && nw > 0L) {
    # a completed warmup without sampling ends as it would before the first
    # sampling iteration, so a run resumed from it continues exactly
    t_warm <- proc.time()[["elapsed"]] - t_start
    if (tune_step) st$eps <- exp(da$x_bar)
  }
  if (stopped) {
    # a forked chain tells the parent it already stopped, so that it is not
    # interrupted a second time while its draws are processed
    if (!is.null(.progress$file)) file.create(paste0(.progress$file, ".stopped"))
    keep <- seq_len(row)
    qs <- qs[keep, , drop = FALSE]
    diag_out <- diag_out[keep, , drop = FALSE]
    if (keep_pg) { dp <- dp[keep, , drop = FALSE]; dg <- dg[keep, , drop = FALSE] }
    thin_idx <- thin_idx[keep]
    if (is.na(t_warm)) t_warm <- proc.time()[["elapsed"]] - t_start
  }
  if (nw == 0L) t_warm <- 0
  t_total <- proc.time()[["elapsed"]] - t_start
  list(q = qs, diag = diag_out, step_size = st$eps, inv_metric = st$metric$inv,
       rejections = if (is.function(attr(pot, "rejections"))) attr(pot, "rejections")(),
       n_warmup_saved = sum(thin_idx <= nw), time = c(warmup = t_warm, sampling = t_total - t_warm), interrupted = stopped,
       divergent_q = if (length(divs)) do.call(rbind, divs),
       diagnostics = if (keep_pg) list(p = dp, g = dg),
       # the state after the last transition, for an exact continuation: an
       # interrupt can stop a transition halfway through its random numbers
       final = list(q = q, rng = get0(".Random.seed", envir = globalenv(), inherits = FALSE), exact = !stopped,
                    center = st$center, radial_sd = if (radial) rad_sd, n_sampling = off + max(0L, n_done - nw)))
}

# GIST (Bou-Rabee, Carpenter and Marsden 2024): the number of leapfrog
# steps L is drawn uniformly from [max(1, floor(psi U)), U], where U is the
# number of steps until the trajectory makes a U-turn; the proposal's own
# U-turn length in the reversed direction enters the acceptance ratio. The
# backward U-turn search reuses the forward states, so only the steps past
# the start are new.
gist_psi <- 0.5
gist_transition <- function(st, q) {
  m <- st$metric
  cap <- 2L^st$max_depth
  r <- if (identical(st$last_q, q)) st$last_r else st$pot(q)
  z <- list(q = q, p = metric_sample_p(m, length(q)), V = r$V, g = r$g)
  H0 <- ps_H(m, z)
  # forward to the U-turn, storing every state
  path <- list(z); H <- H0
  repeat {
    z <- leapfrog(st$pot, m, z, st$eps)
    h <- ps_H(m, z); if (is.nan(h)) h <- Inf
    path[[length(path) + 1L]] <- z; H <- c(H, h)
    k <- length(path) - 1L
    if (h - H0 > st$max_deltaH || sum((z$q - q) * metric_dtau(m, z$p)) < 0 || k >= cap) break
  }
  U1 <- k
  lo1 <- max(1L, floor(gist_psi * U1))
  L <- if (U1 > lo1) sample(lo1:U1, 1L) else U1
  zL <- path[[L + 1L]]; HL <- H[L + 1L]
  # backward from (q_L, -p_L): the states z_{L-1}, ..., z_0 with negated
  # momentum, then new steps past the start
  U2 <- 0L; n_new <- 0L; zb <- NULL
  repeat {
    U2 <- U2 + 1L
    j <- L - U2
    if (j >= 0L) { zj <- path[[j + 1L]]; hj <- H[j + 1L]; pj <- -zj$p }
    else {
      if (is.null(zb)) zb <- list(q = path[[1L]]$q, p = -path[[1L]]$p, V = path[[1L]]$V, g = path[[1L]]$g)
      zb <- leapfrog(st$pot, m, zb, st$eps); n_new <- n_new + 1L
      zj <- zb; hj <- ps_H(m, zb); if (is.nan(hj)) hj <- Inf; pj <- zb$p
    }
    if (hj - HL > st$max_deltaH || sum((zj$q - zL$q) * metric_dtau(m, pj)) < 0 || U2 >= cap) break
  }
  lo2 <- max(1L, floor(gist_psi * U2))
  la <- H0 - HL + log(U1 - lo1 + 1) - log(U2 - lo2 + 1)
  acc <- if (L >= lo2 && L <= U2 && is.finite(la)) min(1, exp(la)) else 0
  out <- if (stats::runif(1) < acc) zL else path[[1L]]
  st$last_q <- out$q
  st$last_r <- list(V = out$V, g = out$g)
  # the step size adapts to the energy error alone, as for NUTS: the
  # path length part of the acceptance does not shrink with the step
  acc_e <- if (is.finite(HL)) min(1, exp(H0 - HL)) else 0
  list(q = out$q, lp = -out$V, accept_stat = acc_e, treedepth = 0L, n_leapfrog = U1 + n_new,
       divergent = (max(H) - H0) > st$max_deltaH, energy = ps_H(m, out), p = out$p, g = out$g)
}

# MALT (Riou-Durand and Vogrinc 2023): a trajectory of OBABO steps of the
# kinetic Langevin diffusion with one Metropolis test for the whole
# trajectory, on the sum of the energy errors of the leapfrog parts (the
# partial momentum refreshes O are exact). The damping is fixed at 1 on the
# scale of the adapted metric, between HMC (0) and overdamped moves.
malt_gamma <- 1
malt_transition <- function(st, q) {
  m <- st$metric
  n <- length(q)
  r <- if (identical(st$last_q, q)) st$last_r else st$pot(q)
  z <- list(q = q, p = metric_sample_p(m, n), V = r$V, g = r$g)
  z0 <- z
  L <- max(1L, as.integer(st$int_time / st$eps))
  eta <- exp(-malt_gamma * st$eps / 2)
  s <- sqrt(1 - eta^2)
  delta <- 0
  for (i in seq_len(L)) {
    z$p <- eta * z$p + s * metric_sample_p(m, n)
    h0 <- ps_H(m, z)
    z <- leapfrog(st$pot, m, z, st$eps)
    h1 <- ps_H(m, z); if (is.nan(h1)) h1 <- Inf
    delta <- delta + (h1 - h0)
    if (!is.finite(delta)) break
    z$p <- eta * z$p + s * metric_sample_p(m, n)
  }
  acc <- if (is.finite(delta)) min(1, exp(-delta)) else 0
  if (stats::runif(1) > acc) z <- z0
  st$last_q <- z$q
  st$last_r <- list(V = z$V, g = z$g)
  list(q = z$q, lp = -z$V, accept_stat = acc, treedepth = 0L, n_leapfrog = L,
       divergent = delta > st$max_deltaH, energy = ps_H(m, z), p = z$p, g = z$g)
}

# AAPS (Sherlock, Urbas and Ludkin 2023): the leapfrog path is cut into
# segments at apogees, where the potential stops rising; the current
# segment is placed uniformly among K + 1 segments and the proposal is drawn
# from all their points with weight pi(z) |x - x_current|^2 (their scheme
# 3). The acceptance ratio compares the weighted squared distances of the
# path to the current point and to the proposal.
# ponytail: K is fixed; their tuning run for K is the upgrade
aaps_K <- 4L
aaps_transition <- function(st, q) {
  m <- st$metric
  cap <- 2L^st$max_depth
  r <- if (identical(st$last_q, q)) st$last_r else st$pot(q)
  z0 <- list(q = q, p = metric_sample_p(m, length(q)), V = r$V, g = r$g)
  cseg <- sample.int(aaps_K + 1L, 1L) - 1L
  rate <- function(z, dir) dir * sum(z$g * metric_dtau(m, z$p))  # dU/dt in forward time
  # walk in one direction until the given number of apogees is passed;
  # returns the points before the last apogee
  walk <- function(dir, n_apo) {
    z <- z0; if (dir < 0) z$p <- -z$p
    out <- list(); s_prev <- rate(z, dir); k <- 0L; hit <- FALSE
    repeat {
      z <- leapfrog(st$pot, m, z, st$eps); k <- k + 1L
      # a point whose position or energy overflowed ends the walk as a
      # divergence; the same segments are walked from any of their points,
      # so the reverse move sees it too
      if (!is.finite(z$V) || !is.finite(ps_H(m, z)) || !all(is.finite(z$q)) || k >= cap) { hit <- TRUE; break }
      s <- rate(z, dir)
      # forward in time an apogee lies between a rise and a fall; walking
      # backward the order of the two is reversed
      if ((dir > 0 && s_prev > 0 && s <= 0) || (dir < 0 && s_prev <= 0 && s > 0)) {
        n_apo <- n_apo - 1L
        if (n_apo == 0L) break
      }
      out[[length(out) + 1L]] <- z; s_prev <- s
    }
    list(z = out, k = k, hit = hit)
  }
  fw <- walk(1, aaps_K - cseg + 1L)
  bw <- walk(-1, cseg + 1L)
  pts <- c(rev(bw$z), list(z0), fw$z)
  X <- do.call(rbind, lapply(pts, `[[`, "q"))
  H <- vapply(pts, function(z) ps_H(m, z), 0)
  lw <- -(H - min(H))
  d0 <- rowSums(sweep(X, 2L, q)^2)
  w <- exp(lw) * d0
  i <- if (sum(w) > 0) sample.int(length(pts), 1L, prob = w) else length(bw$z) + 1L
  di <- rowSums(sweep(X, 2L, X[i, ])^2)
  acc <- if (fw$hit || bw$hit) 0 else min(1, sum(exp(lw) * d0) / sum(exp(lw) * di))
  # every point at the current one (no segment left the start): stay
  if (!is.finite(acc)) acc <- 0
  out <- if (stats::runif(1) < acc) pts[[i]] else z0
  st$last_q <- out$q
  st$last_r <- list(V = out$V, g = out$g)
  # the step size adapts to the energy error along the path, as for NUTS
  H0 <- ps_H(m, z0)
  acc_e <- mean(pmin(1, exp(H0 - H[is.finite(H)])))
  list(q = out$q, lp = -out$V, accept_stat = acc_e, treedepth = 0L, n_leapfrog = fw$k + bw$k,
       divergent = fw$hit || bw$hit, energy = ps_H(m, out), p = out$p, g = out$g)
}

# WALNUTS (Bou-Rabee, Carpenter, Kleppe and Liu 2025): NUTS on a
# grid of macro steps, each integrated with the fewest leapfrog micro steps,
# a power of two, that keep the energy error of the macro step within
# delta. Reversibility needs the reversed macro step to select the same
# number of micro steps; when a smaller number would do in reverse, the
# point and every later point in that direction get weight zero (the
# deterministic choice of micro steps, which they find fastest in practice).
walnuts_delta <- 0.3
walnuts_max_micro <- 1024L
walnuts_err <- function(st, z, h, ell) {
  zz <- z
  for (j in seq_len(ell)) zz <- leapfrog(st$pot, st$metric, zz, h / ell)
  st$n_leapfrog <- st$n_leapfrog + ell
  e <- abs(ps_H(st$metric, zz) - ps_H(st$metric, z))
  list(z = zz, ok = is.finite(e) && e <= walnuts_delta)
}
walnuts_step <- function(st, z, h) {
  lr <- z$lr %||% 0
  ell <- 1L
  repeat {
    f <- walnuts_err(st, z, h, ell)
    if (f$ok || ell >= walnuts_max_micro) break
    ell <- 2L * ell
  }
  if (ell == 1L) st$n_coarse <- st$n_coarse + 1L
  z1 <- f$z
  # the reversed macro step from z1 must not succeed with fewer micro steps
  if (is.finite(lr) && f$ok) {
    k <- 1L
    while (k < ell) {
      if (walnuts_err(st, z1, -h, k)$ok) { lr <- -Inf; break }
      k <- 2L * k
    }
  }
  z1$lr <- lr
  z1
}
