# State space models: bootstrap particle filtering, particle marginal
# Metropolis-Hastings (Andrieu, Doucet and Holenstein 2010) and particle
# Gibbs with ancestor sampling (Lindsten, Jordan and Schon 2014).
#
# The latent process is supplied as R functions; the static parameters and
# their priors come from a bayescomp_model, so transforms, priors and
# output follow the rest of the package. Functions receive `theta`, a named
# numeric vector of the constrained parameters (flat Stan names, so a
# vector b appears as "b[1]", "b[2]", ...).

#' Define a state space model
#'
#' @param y Observations: a vector (one per time) or a matrix with one row
#'   per time.
#' @param init `function(N, theta)` returning `N` draws of the initial state
#'   (a vector for scalar states, else an `N` by `d` matrix).
#' @param transition `function(x, t, theta)` propagating the particles `x`
#'   from time `t - 1` to `t`.
#' @param log_obs `function(y_t, x, t, theta)` returning the log
#'   observation density of each particle.
#' @param log_init,log_transition Optional log densities
#'   `function(x, theta)` and `function(x_new, x_old, t, theta)` (vectorized
#'   over particles), needed by particle Gibbs.
#' @return A `bayescomp_ssm` object.
#' @export
bayescomp_ssm <- function(y, init, transition, log_obs, log_init = NULL, log_transition = NULL) {
  for (nm in c("init", "transition", "log_obs")) if (!is.function(get(nm))) stop(sprintf("'%s' must be a function", nm), call. = FALSE)
  y <- if (is.matrix(y)) y else matrix(y, ncol = 1L)
  structure(list(y = y, T = nrow(y), init = init, transition = transition, log_obs = log_obs,
                 log_init = log_init, log_transition = log_transition), class = "bayescomp_ssm")
}

as_pmat <- function(x) if (is.matrix(x)) x else matrix(x, ncol = 1L)
y_at <- function(ssm, t) { v <- ssm$y[t, ]; if (length(v) == 1L) v[[1L]] else v }

#' Bootstrap particle filter
#'
#' @param ssm A [bayescomp_ssm()].
#' @param theta Named numeric vector of parameters passed to the model
#'   functions.
#' @param num_particles Number of particles.
#' @param resampling Resampling scheme, as in [bayescomp_smc()].
#' @param ess_threshold Resample when the effective sample size falls below
#'   this fraction of the particles.
#' @return A list with `log_lik` (unbiased on the natural scale),
#'   `filtered_mean` (a `T` by `d` matrix), `ess` (per time) and one
#'   trajectory `path` drawn from the final particle system.
#' @export
bayescomp_particle_filter <- function(ssm, theta, num_particles = 500L,
                                      resampling = c("systematic", "stratified", "residual", "multinomial"),
                                      ess_threshold = 1) {
  resampling <- match.arg(resampling)
  N <- as_count(num_particles, "num_particles", 2L)
  pf_run(ssm, theta, N, resampling, ess_threshold, ref = NULL)
}

resample_idx <- function(w, scheme) switch(scheme, systematic = systematic_resample(w), stratified = stratified_resample(w),
                                           residual = residual_resample(w), multinomial = multinomial_resample(w))

# One particle filter pass. With `ref` (a T by d reference path) this is the
# conditional SMC of particle Gibbs with ancestor sampling: particle N is
# pinned to the reference and its ancestor is resampled.
pf_run <- function(ssm, theta, N, scheme, ess_threshold, ref = NULL) {
  Tn <- ssm$T
  X <- as_pmat(ssm$init(N, theta))
  d <- ncol(X)
  if (!is.null(ref)) X[N, ] <- ref[1L, ]
  paths <- array(NA_real_, c(Tn, N, d))
  anc <- matrix(seq_len(N), Tn, N, byrow = TRUE)
  lw <- ssm$log_obs(y_at(ssm, 1L), drop_d(X), 1L, theta)
  ll <- 0
  fm <- matrix(NA_real_, Tn, d); ess <- numeric(Tn)
  paths[1L, , ] <- X
  for (t in seq_len(Tn)) {
    if (t > 1L) {
      m <- max(lw)
      w <- exp(lw - m); w <- w / sum(w)
      e <- 1 / sum(w^2)
      if (!is.null(ref) || e < ess_threshold * N) {
        # conditional SMC resamples the free particles multinomially, which
        # keeps the conditional update exact; the pinned particle's
        # ancestor is drawn by ancestor sampling
        a <- if (is.null(ref)) resample_idx(w, scheme) else c(sample.int(N, N - 1L, TRUE, prob = w), NA_integer_)
        if (!is.null(ref)) {
          la <- log(w) + ssm$log_transition(matrix(ref[t, ], N, d, byrow = TRUE) |> drop_d(), drop_d(X), t, theta)
          la[!is.finite(la)] <- -Inf
          pa <- exp(la - max(la))
          a[N] <- sample.int(N, 1L, prob = pa)
        }
        X <- X[a, , drop = FALSE]
        lw_prev <- rep(0, N)
      } else {
        a <- seq_len(N)
        lw_prev <- lw - m - log(sum(exp(lw - m)))
      }
      anc[t, ] <- a
      X <- as_pmat(ssm$transition(drop_d(X), t, theta))
      if (!is.null(ref)) X[N, ] <- ref[t, ]
      paths[t, , ] <- X
      lo <- ssm$log_obs(y_at(ssm, t), drop_d(X), t, theta)
      lw <- lw_prev + lo
      # incremental likelihood: log sum of normalized previous weights times obs density
    }
    m <- max(lw)
    if (!is.finite(m)) return(list(log_lik = -Inf, filtered_mean = fm, ess = ess, path = NULL))
    w <- exp(lw - m)
    if (t == 1L) ll <- m + log(mean(w)) else ll <- ll + pf_increment(lw, lw_prev)
    w <- w / sum(w)
    fm[t, ] <- colSums(X * w)
    ess[t] <- 1 / sum(w^2)
  }
  # trace one path back through the ancestors
  k <- sample.int(N, 1L, prob = w)
  path <- matrix(NA_real_, Tn, d)
  for (t in Tn:1) {
    path[t, ] <- paths[t, k, ]
    if (t > 1L) k <- anc[t, k]
  }
  list(log_lik = ll, filtered_mean = fm, ess = ess, path = path)
}

drop_d <- function(X) if (ncol(X) == 1L) X[, 1L] else X

# log of sum_i W_{t-1,i} g_t(x_i), with W the normalized previous weights
# (all equal after resampling)
# computed directly so particles with zero weight cannot produce Inf - Inf
pf_increment <- function(lw, lw_prev) log_sum_exp_v(lw) - log_sum_exp_v(lw_prev)

ssm_theta <- function(model, u) {
  v <- unlist(lapply(bayescomp_constrain(model, u, include_tparams = TRUE, include_gqs = FALSE), as.vector), use.names = FALSE)
  names(v) <- bayescomp_param_names(model, include_gqs = FALSE)
  v
}

#' Particle MCMC for state space models
#'
#' `"pmmh"` (particle marginal Metropolis-Hastings) runs adaptive random walk
#' Metropolis on the static parameters with the particle filter estimate of
#' the likelihood, which leaves the exact posterior invariant. `"pgibbs"`
#' (particle Gibbs with ancestor sampling) alternates a conditional SMC
#' update of the latent path with a Metropolis update of the parameters
#' given the path, and needs `log_init` and `log_transition`.
#'
#' @param model A `bayescomp_model` declaring the static parameters and
#'   their priors (its data block may be empty).
#' @param ssm A [bayescomp_ssm()].
#' @param method `"pmmh"` or `"pgibbs"`.
#' @param num_particles Particles per filter.
#' @param chains,iter_warmup,iter_sampling,init,seed,cores,refresh As in
#'   [bayescomp_sample()].
#' @param theta_steps Metropolis steps on the parameters per particle Gibbs
#'   sweep.
#' @return A `bayescomp_fit` with method `"pmcmc"`; element `states` holds
#'   the sampled latent paths as an array (draw, time, state dimension) for
#'   the first chain, and `acceptance` the parameter acceptance per chain.
#' @references Andrieu, C., Doucet, A. and Holenstein, R. (2010). Particle
#'   Markov chain Monte Carlo methods. *Journal of the Royal Statistical
#'   Society B* 72, 269 to 342. Lindsten, F., Jordan, M. I. and Schon, T. B.
#'   (2014). Particle Gibbs with ancestor sampling. *Journal of Machine
#'   Learning Research* 15, 2145 to 2184.
#' @export
bayescomp_pmcmc <- function(model, ssm, method = c("pmmh", "pgibbs"), num_particles = 200L, chains = 2L,
                            iter_warmup = 500L, iter_sampling = 1000L, theta_steps = 5L, init = 2,
                            seed = NULL, cores = 1L, refresh = 0L) {
  check_model(model)
  if (!inherits(ssm, "bayescomp_ssm")) stop("'ssm' must come from bayescomp_ssm()", call. = FALSE)
  method <- match.arg(method)
  if (method == "pgibbs" && (is.null(ssm$log_init) || is.null(ssm$log_transition)))
    stop("particle Gibbs needs 'log_init' and 'log_transition' in bayescomp_ssm()", call. = FALSE)
  n <- model$dim
  N <- as_count(num_particles, "num_particles", 2L)
  chains <- as_count(chains, "chains", 1L)
  nw <- as_count(iter_warmup, "iter_warmup")
  ns <- as_count(iter_sampling, "iter_sampling", 1L)
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, chains, seed, cores, function(id) {
    lprior <- lp_value_fun(model)
    u <- find_init(model, init, id, gradient = FALSE)
    th <- ssm_theta(model, u)
    pf <- pf_run(ssm, th, N, "systematic", 1)
    path <- pf$path
    joint_ll <- function(th, path) {
      v <- sum(ssm$log_init(drop_d(path[1L, , drop = FALSE]), th))
      for (t in seq_len(ssm$T)) {
        if (t > 1L) v <- v + sum(ssm$log_transition(drop_d(path[t, , drop = FALSE]), drop_d(path[t - 1L, , drop = FALSE]), t, th))
        v <- v + sum(ssm$log_obs(y_at(ssm, t), drop_d(path[t, , drop = FALSE]), t, th))
      }
      v
    }
    cur <- lprior(u) + if (method == "pmmh") pf$log_lik else joint_ll(th, path)
    sc <- 2.38 / sqrt(max(n, 1L)); C <- diag(0.1, max(n, 1L))
    est <- est_new(n, TRUE)
    qs <- matrix(NA_real_, ns, n)
    states <- array(NA_real_, c(ns, ssm$T, ncol(path)))
    dg <- matrix(0, ns, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "energy__"] <- NA_real_
    acc_n <- 0
    for (it in seq_len(nw + ns)) {
      warm <- it <= nw
      if (method == "pgibbs") {
        path <- pf_run(ssm, th, N, "systematic", 1, ref = path)$path
        cur <- lprior(u) + joint_ll(th, path)
      }
      steps <- if (method == "pmmh") 1L else theta_steps
      a_mean <- 0
      if (n > 0L) for (s in seq_len(steps)) {
        v <- u + sc * as.vector(t(chol(C)) %*% stats::rnorm(n))
        thv <- ssm_theta(model, v)
        lpv <- lprior(v)
        if (is.finite(lpv)) {
          if (method == "pmmh") { pv <- pf_run(ssm, thv, N, "systematic", 1); llv <- pv$log_lik }
          else llv <- joint_ll(thv, path)
          new <- lpv + llv
        } else new <- -Inf
        a <- if (is.finite(new)) min(1, exp(new - cur)) else 0
        if (stats::runif(1) < a) {
          u <- v; th <- thv; cur <- new
          if (method == "pmmh") path <- pv$path
        }
        a_mean <- a_mean + a / steps
      }
      if (warm && n > 0L) {
        sc <- sc * exp((a_mean - if (method == "pmmh") 0.15 else 0.3) / sqrt(it))
        est_add(est, u)
        if (it %% 50L == 0L && est$n > 2 * n) {
          Cn <- est_value(est)
          if (!inherits(tryCatch(chol(Cn), error = function(e) e), "error")) C <- Cn
        }
      }
      if (!warm) {
        i <- it - nw
        qs[i, ] <- u
        states[i, , ] <- path
        dg[i, c("lp__", "accept_stat__")] <- c(cur, a_mean)
        acc_n <- acc_n + a_mean
      }
      if (refresh > 0L && it %% refresh == 0L) message(sprintf("Chain %d: Iteration %d / %d", id, it, nw + ns))
    }
    list(q = qs, gq = chain_write_draws(model, qs, FALSE), diag = dg, step_size = sc, inv_metric = C,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0), states = states, acc = acc_n / ns)
  })
  fit <- make_fit(model, res, list(algorithm = method, num_particles = N), "pmcmc", seed, proc.time()[["elapsed"]] - t0)
  fit$states <- res[[1L]]$states
  fit$acceptance <- vapply(res, `[[`, 0, "acc")
  fit
}
