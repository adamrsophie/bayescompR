# Piecewise deterministic Markov process samplers: the Zig-Zag sampler
# (Bierkens, Fearnhead and Roberts 2019) and the bouncy particle sampler
# (Bouchard-Cote, Vollmer and Doucet 2018). The state moves in straight
# lines and the velocity changes at the events of a Poisson process whose
# rate depends on the gradient along the path. Event times are found by
# integrating the rate along the path with a fixed time step (Pakman et
# al. 2017 style numerical inversion), so the samplers are exact only in
# the limit of a small step; the step is set relative to the posterior
# scale estimated during warmup. Draws are the positions of the continuous
# trajectory at equally spaced times.

#' Piecewise deterministic samplers
#'
#' @param model A `bayescomp_model`.
#' @param method `"zigzag"` or `"bps"` (bouncy particle sampler).
#' @param chains,iter_warmup,iter_sampling,init,seed,cores,refresh As in
#'   [bayescomp_sample()]; one iteration is one unit of trajectory time
#'   between saved draws.
#' @param refresh_rate Velocity refreshment rate for `"bps"` (needed for
#'   ergodicity) and an optional small refreshment for `"zigzag"`.
#' @param dt Integration step for the event rate, relative to the
#'   trajectory time between draws.
#' @return A `bayescomp_fit` with method `"pdmp"` and element `events`
#'   (mean events per unit time per chain).
#' @references Bierkens, J., Fearnhead, P. and Roberts, G. (2019). The
#'   Zig-Zag process and super efficient sampling for Bayesian analysis of
#'   big data. *Annals of Statistics* 47, 1288 to 1320. Bouchard-Cote, A.,
#'   Vollmer, S. J. and Doucet, A. (2018). The bouncy particle sampler.
#'   *Journal of the American Statistical Association* 113, 855 to 867.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_pdmp(m, chains = 1, iter_warmup = 100, iter_sampling = 200, seed = 1)
#' mean(bayescomp_extract(fit, "mu"))
#' @export
bayescomp_pdmp <- function(model, method = c("zigzag", "bps"), chains = 4L, iter_warmup = 500L,
                           iter_sampling = 1000L, refresh_rate = NULL, dt = 0.05, init = 2, seed = NULL,
                           cores = getOption("mc.cores", 1L), refresh = 0L) {
  check_model(model)
  method <- match.arg(method)
  if (!is.numeric(dt) || length(dt) != 1L || !is.finite(dt) || dt <= 0)
    stop("'dt' must be one positive finite number", call. = FALSE)
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  chains <- as_count(chains, "chains", 1L)
  nw <- as_count(iter_warmup, "iter_warmup")
  ns <- as_count(iter_sampling, "iter_sampling", 1L)
  if (is.null(refresh_rate)) refresh_rate <- if (method == "bps") 1 else 0
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, chains, seed, cores, function(id) {
    pot <- hmc_potential(model)
    x <- find_init(model, init, id)
    r <- pot(x)
    sc <- rep(1, n)          # coordinate scales, set from warmup
    newv <- function() if (method == "zigzag") sample(c(-1, 1), n, TRUE) * sc else stats::rnorm(n) * sc
    v <- newv()
    est <- est_new(n, FALSE)
    qs <- matrix(NA_real_, ns, n)
    dg <- matrix(0, ns, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "energy__"] <- NA_real_
    n_events <- 0
    rates <- function(g, v) if (method == "zigzag") pmax(0, v * g) else max(0, sum(v * g))
    for (it in seq_len(nw + ns)) {
      h <- dt
      t_left <- 1
      E <- stats::rexp(length(rates(r$g, v)))
      Er <- if (refresh_rate > 0) stats::rexp(1) / refresh_rate else Inf
      evals <- 0L
      lam0 <- rates(r$g, v)
      while (t_left > 1e-12) {
        s <- min(h, t_left)
        xn <- x + s * v
        rn <- pot(xn)
        evals <- evals + 1L
        if (!is.finite(rn$V)) {
          # outside the support: reflect off the boundary by reversing
          v <- -v
          lam0 <- rates(r$g, v)
          n_events <- n_events + 1
          t_left <- t_left - s
          next
        }
        lam1 <- rates(rn$g, v)
        inc <- 0.5 * s * (lam0 + lam1)
        # the earliest of the bounce and refresh clocks within this step
        bounce <- any(inc >= E)
        frac <- if (bounce) min(E / pmax(inc, 1e-300)) else Inf
        if (Er <= s && Er < frac * s) {
          # refresh first: move to the refresh time and draw a new velocity
          x <- x + Er * v; r <- pot(x); evals <- evals + 1L
          E <- E - (Er / s) * inc
          t_left <- t_left - Er
          v <- newv(); lam0 <- rates(r$g, v)
          Er <- stats::rexp(1) / refresh_rate
          next
        }
        if (bounce) {
          # event inside the step: locate it linearly within the step
          j <- which.min(E / pmax(inc, 1e-300))
          s2 <- frac * s
          Er <- Er - s2
          x <- x + s2 * v
          r <- pot(x)
          evals <- evals + 1L
          if (method == "zigzag") {
            v[j] <- -v[j]
            E <- E - frac * inc
            E[j] <- stats::rexp(1)
          } else {
            g <- r$g * sc
            vs <- v / sc
            vs <- vs - 2 * sum(vs * g) / sum(g^2) * g
            v <- vs * sc
            E <- stats::rexp(1)
          }
          lam0 <- rates(r$g, v)
          n_events <- n_events + 1
          t_left <- t_left - s2
        } else {
          E <- E - inc
          Er <- Er - s
          x <- xn; r <- rn; lam0 <- lam1
          t_left <- t_left - s
        }
      }
      if (it <= nw) {
        est_add(est, x)
        if (it %in% c(50L, 100L, 200L, 400L) && est$n > 10 && it < nw) {
          # rescale speeds so the trajectory covers about one posterior
          # standard deviation per unit time
          sc <- sqrt(pmax(est_value(est), 1e-10))
          v <- newv()
          E <- stats::rexp(length(rates(r$g, v)))
        }
      } else {
        i <- it - nw
        qs[i, ] <- x
        dg[i, c("lp__", "accept_stat__", "n_leapfrog__")] <- c(-r$V, 1, evals)
      }
      if (refresh > 0L && it %% refresh == 0L) chain_message(sprintf("Chain %d: Iteration %d / %d", id, it, nw + ns))
    }
    list(q = qs, gq = chain_write_draws(model, qs, FALSE), diag = dg, step_size = dt, inv_metric = sc^2,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0), events = n_events / (nw + ns))
  })
  fit <- make_fit(model, res, list(algorithm = method), "pdmp", seed, proc.time()[["elapsed"]] - t0)
  fit$events <- vapply(res, `[[`, 0, "events")
  fit
}
