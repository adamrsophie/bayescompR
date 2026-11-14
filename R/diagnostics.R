# Convergence diagnostics: rank normalized split R-hat, bulk and tail
# effective sample sizes and Monte Carlo standard errors (Vehtari, Gelman,
# Simpson, Carpenter and Buerkner 2021). The formulas follow the reference
# implementation in the posterior package.

as_draws_matrix_chains <- function(x) {
  # iterations by chains
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  x
}

split_chains <- function(x) {
  x <- as_draws_matrix_chains(x)
  n <- nrow(x)
  if (n < 2L) return(x)
  half <- n %/% 2L
  a <- x[seq_len(half), , drop = FALSE]
  b <- x[(n - half + 1L):n, , drop = FALSE]
  cbind(a, b)
}

# constant up to rounding, relative to the scale of the draws
na_diag <- function(x) any(!is.finite(x)) || (max(x) - min(x)) <= .Machine$double.eps * max(abs(x))

z_scale <- function(x) {
  d <- dim(x)
  r <- rank(as.vector(x), ties.method = "average")
  z <- stats::qnorm((r - 3 / 8) / (length(r) - 2 * 3 / 8 + 1))
  dim(z) <- d
  z
}

autocovariance <- function(x) {
  N <- length(x)
  Mt2 <- 2L * stats::nextn(N)
  yc <- c(x - mean(x), rep.int(0, Mt2 - N))
  tr <- stats::fft(yc)
  ac <- Re(stats::fft(Conj(tr) * tr, inverse = TRUE))[seq_len(N)]
  ac / (N * Mt2)
}

ess_raw <- function(sims) {
  sims <- as_draws_matrix_chains(sims)
  chains <- ncol(sims)
  n <- nrow(sims)
  if (n < 3L || na_diag(sims)) return(NA_real_)
  acov <- vapply(seq_len(chains), function(i) autocovariance(sims[, i]), numeric(n))
  if (is.null(dim(acov))) acov <- matrix(acov, nrow = n)
  chain_mean <- colMeans(sims)
  mean_var <- mean(acov[1L, ]) * n / (n - 1)
  var_plus <- mean_var * (n - 1) / n
  if (chains > 1L) var_plus <- var_plus + stats::var(chain_mean)
  rho <- rep.int(0, n)
  t <- 0L
  rho_even <- 1
  rho[1L] <- rho_even
  rho_odd <- 1 - (mean_var - mean(acov[2L, ])) / var_plus
  rho[2L] <- rho_odd
  while (t < n - 5L && !is.nan(rho_even + rho_odd) && (rho_even + rho_odd > 0)) {
    t <- t + 2L
    rho_even <- 1 - (mean_var - mean(acov[t + 1L, ])) / var_plus
    rho_odd <- 1 - (mean_var - mean(acov[t + 2L, ])) / var_plus
    if ((rho_even + rho_odd) >= 0) {
      rho[t + 1L] <- rho_even
      rho[t + 2L] <- rho_odd
    }
  }
  max_t <- t
  if (rho_even > 0) rho[max_t + 1L] <- rho_even
  t <- 0L
  while (t <= max_t - 4L) {
    t <- t + 2L
    if (rho[t + 1L] + rho[t + 2L] > rho[t - 1L] + rho[t]) {
      rho[t + 1L] <- (rho[t - 1L] + rho[t]) / 2
      rho[t + 2L] <- rho[t + 1L]
    }
  }
  ess <- chains * n
  # posterior keeps rho[1] when max_t is 0 (R's 1:0)
  tau <- -1 + 2 * sum(rho[seq_len(max(max_t, 1L))]) + rho[max_t + 1L]
  tau <- max(tau, 1 / log10(ess))
  ess / tau
}

rhat_raw <- function(x) {
  x <- as_draws_matrix_chains(x)
  if (na_diag(x)) return(NA_real_)
  n <- nrow(x)
  vb <- n * stats::var(colMeans(x))
  vw <- mean(apply(x, 2L, stats::var))
  sqrt((vb / vw + n - 1) / n)
}

#' Convergence diagnostics for one quantity
#'
#' @param x Matrix of draws with one column per chain (iterations by
#'   chains), or a vector for a single chain.
#' @param probs Probabilities for `ess_quantile()`.
#' @return A number.
#' @name diagnostics
NULL

#' @rdname diagnostics
#' @export
bayescomp_rhat <- function(x) {
  x <- as_draws_matrix_chains(x)
  if (na_diag(x)) return(NA_real_)
  s <- split_chains(x)
  r_bulk <- rhat_raw(z_scale(s))
  r_tail <- rhat_raw(z_scale(split_chains(abs(x - stats::median(x)))))
  max(r_bulk, r_tail)
}

#' @rdname diagnostics
#' @export
bayescomp_ess_bulk <- function(x) {
  x <- as_draws_matrix_chains(x)
  if (na_diag(x)) return(NA_real_)
  ess_raw(z_scale(split_chains(x)))
}

#' @rdname diagnostics
#' @export
bayescomp_ess_tail <- function(x) {
  x <- as_draws_matrix_chains(x)
  if (na_diag(x)) return(NA_real_)
  e <- bayescomp_ess_quantile(x, c(0.05, 0.95))
  if (anyNA(e)) NA_real_ else min(e)
}

#' @rdname diagnostics
#' @export
bayescomp_ess_basic <- function(x) ess_raw(split_chains(as_draws_matrix_chains(x)))

#' @rdname diagnostics
#' @export
bayescomp_ess_quantile <- function(x, probs = c(0.05, 0.95)) {
  x <- as_draws_matrix_chains(x)
  # as posterior 1.7.1: the indicator x <= q for every quantile
  vapply(probs, function(p) {
    q <- stats::quantile(x, p, names = FALSE)
    ind <- x <= q
    if (na_diag(ind + 0)) return(NA_real_)
    ess_raw(split_chains(ind + 0))
  }, 0)
}

#' @rdname diagnostics
#' @export
bayescomp_mcse_mean <- function(x) {
  x <- as_draws_matrix_chains(x)
  stats::sd(as.vector(x)) / sqrt(bayescomp_ess_basic(x))
}

#' @rdname diagnostics
#' @export
bayescomp_mcse_sd <- function(x) {
  # delta method: Var(sd) is about Var((x - m)^2) / (4 E[(x - m)^2] ESS)
  x <- as_draws_matrix_chains(x)
  xc <- x - mean(x)
  x2 <- xc^2
  ess <- bayescomp_ess_basic(x2)
  ev <- mean(x2)
  sqrt(mean((x2 - ev)^2) / 4 / ev / ess)
}
