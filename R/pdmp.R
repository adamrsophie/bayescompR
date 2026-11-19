# Piecewise deterministic Markov process samplers: the Zig-Zag sampler
# (Bierkens, Fearnhead and Roberts 2019) and the bouncy particle sampler
# (Bouchard-Cote, Vollmer and Doucet 2018). The state moves in straight
# lines and the velocity changes at the events of a Poisson process whose
# rate depends on the gradient along the path. Event times are simulated by
# thinning against a piecewise constant bound of the rate built on a grid of
# gradient evaluations along the path over an adaptive horizon (automatic
# Poisson thinning, Andral and Kamatani 2024). The events are exact whenever
# the bound holds; a proposal above the bound is counted as a violation and
# shortens the horizon. Speeds are set relative to the posterior scale
# estimated during warmup. Draws are the positions of the continuous
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
#' @param grid_size Number of grid intervals of the rate bound over one
#'   horizon; each costs one gradient.
#' @param dt Deprecated and ignored; event times no longer use a fixed
#'   integration step.
#' @return A `bayescomp_fit` with method `"pdmp"` and elements `events`
#'   (mean events per unit time per chain) and `bound_violations` (proposals
#'   whose rate exceeded the grid bound, per chain; each makes that event
#'   inexact, so a large count calls for a larger `grid_size`).
#' @references Bierkens, J., Fearnhead, P. and Roberts, G. (2019). The
#'   Zig-Zag process and super efficient sampling for Bayesian analysis of
#'   big data. *Annals of Statistics* 47, 1288 to 1320. Bouchard-Cote, A.,
#'   Vollmer, S. J. and Doucet, A. (2018). The bouncy particle sampler.
#'   *Journal of the American Statistical Association* 113, 855 to 867.
#'   Andral, C. and Kamatani, K. (2024). Automated techniques for efficient
#'   sampling of piecewise-deterministic Markov processes. arXiv:2408.03682.
#' @examples
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_pdmp(m, chains = 1, iter_warmup = 100, iter_sampling = 200, seed = 1)
#' mean(bayescomp_extract(fit, "mu"))
#' @export
bayescomp_pdmp <- function(model, method = c("zigzag", "bps"), chains = 4L, iter_warmup = 500L,
                           iter_sampling = 1000L, refresh_rate = NULL, grid_size = 10L, dt = NULL, init = 2, seed = NULL,
                           cores = getOption("mc.cores", 1L), refresh = 0L) {
  check_model(model)
  method <- match.arg(method)
  if (!is.null(dt)) warning("'dt' is deprecated and ignored: event times are simulated by thinning", call. = FALSE)
  grid_size <- as_count(grid_size, "grid_size", 2L)
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
    n_viol <- 0
    H <- 1                   # horizon of one bound, adapted as in pdmp_jax
    rates <- function(g, v) if (method == "zigzag") pmax(0, v * g) else max(0, sum(v * g))
    for (it in seq_len(nw + ns)) {
      t_left <- 1
      Er <- if (refresh_rate > 0) stats::rexp(1) / refresh_rate else Inf
      evals <- 0L
      while (t_left > 1e-12) {
        Hs <- min(H, t_left, Er)
        # the rates on a grid along the ray; the bound covers [0, Hs]
        tg <- seq(0, Hs, length.out = grid_size + 1L)
        Lam <- matrix(0, grid_size + 1L, length(r$g) ^ (method == "zigzag"))
        Lam[1L, ] <- rates(r$g, v)
        last <- r; K <- grid_size + 1L
        for (k in 2:K) {
          rk <- pot(x + tg[k] * v); evals <- evals + 1L
          if (!is.finite(rk$V)) { K <- k - 1L; break }
          Lam[k, ] <- rates(rk$g, v); last <- rk
        }
        if (K == 1L) {
          # outside the support right away: reflect off the boundary
          v <- -v; n_events <- n_events + 1; H <- H / 2
          next
        }
        if (K < grid_size + 1L) { tg <- tg[seq_len(K)]; Lam <- Lam[seq_len(K), , drop = FALSE]; Hs <- tg[K] }
        B <- pdmp_grid_bound(tg, Lam)
        tot <- rowSums(B)
        cum <- c(0, cumsum(tot * diff(tg)))
        # thinning against the piecewise constant bound
        c_now <- 0; tau <- Hs; hit <- NULL
        repeat {
          target <- c_now + stats::rexp(1)
          if (target >= cum[K]) break
          k <- sum(cum < target)
          tau <- tg[k] + (target - cum[k]) / tot[k]
          rt <- pot(x + tau * v); evals <- evals + 1L
          lam <- rates(rt$g, v)
          j <- if (method == "zigzag") sample.int(length(lam), 1L, prob = B[k, ]) else 1L
          ratio <- lam[j] / B[k, j]
          if (!is.finite(ratio) || ratio > 1) {
            # the grid bound missed a peak: count it and shorten the horizon
            n_viol <- n_viol + 1; H <- H / 2
          }
          if (!is.finite(ratio) || stats::runif(1) < ratio) { hit <- j; break }
          c_now <- target; H <- H / 1.04
        }
        if (is.null(hit)) {
          tau <- Hs
          x <- x + Hs * v; r <- last
          if (Hs == H) H <- H * 1.01
        } else {
          x <- x + tau * v; r <- rt
          if (method == "zigzag") v[hit] <- -v[hit] else {
            g <- r$g * sc
            vs <- v / sc
            vs <- vs - 2 * sum(vs * g) / sum(g^2) * g
            v <- vs * sc
          }
          n_events <- n_events + 1
        }
        t_left <- t_left - tau
        Er <- Er - tau
        if (Er <= 1e-12) { v <- newv(); Er <- stats::rexp(1) / refresh_rate }
      }
      if (it <= nw) {
        est_add(est, x)
        if (it %in% c(50L, 100L, 200L, 400L) && est$n > 10 && it < nw) {
          # rescale speeds so the trajectory covers about one posterior
          # standard deviation per unit time
          sc <- sqrt(pmax(est_value(est), 1e-10))
          v <- newv()
        }
      } else {
        i <- it - nw
        qs[i, ] <- x
        dg[i, c("lp__", "accept_stat__", "stepsize__", "n_leapfrog__")] <- c(-r$V, 1, H, evals)
      }
      if (refresh > 0L && it %% refresh == 0L) chain_message(sprintf("Chain %d: Iteration %d / %d", id, it, nw + ns))
    }
    list(q = qs, gq = chain_write_draws(model, qs, FALSE), diag = dg, step_size = H, inv_metric = sc^2,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0), events = n_events / (nw + ns),
         violations = n_viol)
  })
  fit <- make_fit(model, res, list(algorithm = method), "pdmp", seed, proc.time()[["elapsed"]] - t0)
  fit$events <- vapply(res, `[[`, 0, "events")
  fit$bound_violations <- vapply(res, `[[`, 0, "violations")
  fit
}

# Piecewise constant upper bound of event rates from their values Lam (one
# row per grid time tg, one column per rate) on each grid interval: the
# larger end value, and where the rate bends down, the meeting point of the
# lines through each end with the slope of the neighboring interval
# (Andral and Kamatani 2024 use tangents; secants of the neighbors bound a
# concave stretch the same way without second derivatives). The first and
# last intervals have one neighbor and use its line alone.
pdmp_grid_bound <- function(tg, Lam) {
  K <- nrow(Lam)
  a <- Lam[-K, , drop = FALSE]; b <- Lam[-1L, , drop = FALSE]; d <- diff(tg)
  B <- pmax(a, b)
  if (K >= 4L) {
    sec <- (b - a) / d
    i <- 2:(K - 2L)
    sL <- sec[i - 1L, , drop = FALSE]; sR <- sec[i + 1L, , drop = FALSE]
    ai <- a[i, , drop = FALSE]
    u <- (b[i, , drop = FALSE] - ai - sR * d[i]) / (sL - sR)
    val <- ai + sL * u
    use <- is.finite(u) & sL > sR & u > 0 & u < d[i]
    Bi <- B[i, , drop = FALSE]
    Bi[use] <- pmax(Bi[use], val[use])
    B[i, ] <- Bi
  }
  if (K >= 3L) {
    k1 <- 1L; kl <- K - 1L
    B[k1, ] <- pmax(B[k1, ], b[k1, ] - ((b - a) / d)[2L, ] * d[k1])
    B[kl, ] <- pmax(B[kl, ], a[kl, ] + ((b - a) / d)[kl - 1L, ] * d[kl])
  }
  pmax(B, 0)
}
