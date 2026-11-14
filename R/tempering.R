# Parallel tempering (replica exchange) between a normalized Gaussian
# reference g and the posterior p, through pi_beta = g^(1 - beta) p^beta.
#
# Swaps are non-reversible by default (Syed, Bouchard-Cote, Deligiannidis
# and Doucet 2022): even pairs swap at even iterations and odd pairs at odd
# ones, which turns the index process into a lifted, persistent walk. The
# temperature schedule is adapted during warmup in rounds of doubling length
# so that every pair has the same rejection rate, from the estimated local
# communication barrier. Each replica explores with adaptive random walk
# Metropolis. Adaptation stops at the end of warmup.

#' Parallel tempering
#'
#' @param model A `bayescomp_model`.
#' @param num_temps Number of temperatures, including the reference
#'   (`beta = 0`) and the posterior (`beta = 1`).
#' @param chains,iter_warmup,iter_sampling,seed,cores,refresh,init As in
#'   [bayescomp_sample()].
#' @param explore_steps Metropolis steps per replica between swaps.
#' @param swap `"nonreversible"` (deterministic even and odd swaps) or
#'   `"reversible"` (a random pair at each iteration).
#' @param adapt_schedule Adapt the temperatures during warmup.
#' @param reference,reference_scale Gaussian reference, as in
#'   [bayescomp_smc()].
#' @return A `bayescomp_fit` with draws from the `beta = 1` replica and
#'   elements `betas` (final schedule per chain), `swap_rate` (acceptance per
#'   adjacent pair, per chain), `barrier` (estimated global communication
#'   barrier) and `log_evidence` (stepping stone estimate per chain).
#' @references Syed, S., Bouchard-Cote, A., Deligiannidis, G. and Doucet, A.
#'   (2022). Non-reversible parallel tempering: a scalable highly parallel
#'   MCMC scheme. *Journal of the Royal Statistical Society B* 84, 321 to
#'   350.
#' @export
bayescomp_parallel_tempering <- function(model, num_temps = 10L, chains = 2L, iter_warmup = 1000L,
                                         iter_sampling = 1000L, explore_steps = 3L,
                                         swap = c("nonreversible", "reversible"), adapt_schedule = TRUE,
                                         reference = c("pathfinder", "standard"), reference_scale = 1.5,
                                         init = 2, seed = NULL, cores = 1L, refresh = 0L) {
  check_model(model)
  swap <- match.arg(swap)
  reference <- match.arg(reference)
  K <- as_count(num_temps, "num_temps", 2L)
  chains <- as_count(chains, "chains", 1L)
  nw <- as_count(iter_warmup, "iter_warmup")
  ns <- as_count(iter_sampling, "iter_sampling", 1L)
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, chains, seed, cores, function(id) {
    lp_full <- lp_value_fun(model)
    ref <- gaussian_reference(model, reference, reference_scale)
    betas <- seq(0, 1, length.out = K)
    X <- sweep(matrix(stats::rnorm(K * n), K, n) %*% t(ref$L), 2L, ref$mu, "+")
    X[K, ] <- find_init(model, init, id, gradient = FALSE)
    lp <- vapply(seq_len(K), function(k) lp_full(X[k, ]), 0)
    lg <- ref$log_g(X)
    Ls <- rep(list(ref$L), K)
    sc <- rep(2.38 / sqrt(n), K)
    est <- lapply(seq_len(K), function(k) est_new(n, TRUE))
    rej <- numeric(K - 1L); tries <- numeric(K - 1L)
    round_end <- 32L
    qs <- matrix(NA_real_, ns, n)
    dg <- matrix(0, ns, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    ss <- matrix(NA_real_, ns, K - 1L)   # stepping stone terms
    acc_cold <- 0
    total <- nw + ns
    for (it in seq_len(total)) {
      warm <- it <= nw
      # exploration
      for (k in seq_len(K)) {
        b <- betas[k]
        a_sum <- 0
        for (j in seq_len(explore_steps)) {
          y <- X[k, ] + sc[k] * as.vector(Ls[[k]] %*% stats::rnorm(n))
          lpy <- lp_full(y)
          lgy <- ref$log_g(y)
          # log ratio of g^(1 - b) p^b; zero weight terms are skipped so that
          # 0 * (-Inf) does not arise, and a move between two zero density
          # states is driven by the reference alone
          dp <- lpy - lp[k]
          if (is.nan(dp)) dp <- 0
          la <- (if (b < 1) (1 - b) * (lgy - lg[k]) else 0) + (if (b > 0) b * dp else 0)
          a <- if (is.nan(la)) 0 else min(1, exp(la))
          if (stats::runif(1) < a) { X[k, ] <- y; lp[k] <- lpy; lg[k] <- lgy }
          a_sum <- a_sum + a
        }
        a_mean <- a_sum / explore_steps
        if (k == K) acc_cold <- a_mean
        if (warm) {
          sc[k] <- sc[k] * exp((a_mean - 0.234) / sqrt(it))
          est_add(est[[k]], X[k, ])
        }
      }
      # swaps
      d <- lp - lg
      if (swap == "nonreversible") {
        first <- if (it %% 2L == 0L) 1L else 2L
        pairs <- if (first <= K - 1L) seq(first, K - 1L, by = 2L) else integer(0)
      } else pairs <- sample.int(K - 1L, 1L)
      for (k in pairs) {
        la <- (betas[k + 1L] - betas[k]) * (d[k] - d[k + 1L])
        a <- if (is.finite(la)) min(1, exp(la)) else if (is.nan(la)) 0 else if (la > 0) 1 else 0
        rej[k] <- rej[k] + (1 - a); tries[k] <- tries[k] + 1
        if (stats::runif(1) < a) {
          tmp <- X[k, ]; X[k, ] <- X[k + 1L, ]; X[k + 1L, ] <- tmp
          tmp <- lp[k]; lp[k] <- lp[k + 1L]; lp[k + 1L] <- tmp
          tmp <- lg[k]; lg[k] <- lg[k + 1L]; lg[k + 1L] <- tmp
          tmp <- d[k]; d[k] <- d[k + 1L]; d[k + 1L] <- tmp
        }
      }
      # schedule and proposal adaptation at the end of each warmup round
      if (warm && (it == round_end || it == nw)) {
        r <- ifelse(tries > 0, rej / pmax(tries, 1), 0)
        if (adapt_schedule && sum(r) > 0) betas <- pt_update_schedule(betas, r)
        for (k in seq_len(K)) if (est[[k]]$n > 2 * n) {
          Lk <- tryCatch(t(chol(est_value(est[[k]]))), error = function(e) NULL)
          if (!is.null(Lk)) Ls[[k]] <- Lk
          est[[k]] <- est_new(n, TRUE)
        }
        rej[] <- 0; tries[] <- 0
        round_end <- 2L * round_end
      }
      if (!warm) {
        i <- it - nw
        qs[i, ] <- X[K, ]
        dg[i, c("lp__", "accept_stat__", "stepsize__")] <- c(lp[K], acc_cold, sc[K])
        ss[i, ] <- (betas[-1L] - betas[-K]) * d[-K]
      }
      if (refresh > 0L && it %% refresh == 0L)
        message(sprintf("Chain %d: Iteration %d / %d (%s)", id, it, total, if (warm) "Warmup" else "Sampling"))
    }
    r <- ifelse(tries > 0, rej / pmax(tries, 1), NA_real_)
    log_z <- sum(apply(ss, 2L, function(v) { m <- max(v); m + log(mean(exp(v - m))) }))
    list(q = qs, diag = dg, step_size = sc[K], inv_metric = tcrossprod(Ls[[K]]), n_warmup_saved = 0L,
         time = c(warmup = 0, sampling = 0), betas = betas, swap_rate = 1 - r, barrier = sum(r, na.rm = TRUE),
         log_evidence = log_z, gq = chain_write_draws(model, qs, FALSE))
  })
  fit <- make_fit(model, res, list(algorithm = "parallel_tempering", num_temps = K, swap = swap),
                  "parallel_tempering", seed, proc.time()[["elapsed"]] - t0)
  fit$betas <- lapply(res, `[[`, "betas")
  fit$swap_rate <- lapply(res, `[[`, "swap_rate")
  fit$barrier <- vapply(res, `[[`, 0, "barrier")
  fit$log_evidence <- vapply(res, `[[`, 0, "log_evidence")
  fit
}

# Syed et al. (2022), Algorithm 4: equalize the rejection rates using the
# cumulative communication barrier as a function of beta.
pt_update_schedule <- function(betas, r) {
  K <- length(betas)
  Lam <- c(0, cumsum(r))
  if (Lam[K] <= 0) return(betas)
  target <- seq(0, Lam[K], length.out = K)
  nb <- stats::approx(Lam, betas, xout = target, ties = mean)$y
  nb[1L] <- 0; nb[K] <- 1
  nb <- cummax(nb)
  # keep the schedule strictly increasing
  for (k in 2:K) if (nb[k] <= nb[k - 1L]) nb[k] <- nb[k - 1L] + 1e-6
  nb / nb[K]
}
