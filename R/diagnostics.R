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
#' @param probs Probabilities for `ess_quantile()` and `mcse_quantile()`.
#' @return A number; for `ess_quantile()` and `mcse_quantile()`, one per
#'   probability. `mcse_quantile()` is the Monte Carlo standard error of a
#'   quantile estimate, computed as in the posterior package.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(400), 100, 4)
#' bayescomp_rhat(x)
#' bayescomp_ess_bulk(x)
#' bayescomp_ess_tail(x)
#' bayescomp_mcse_mean(x)
#' bayescomp_mcse_quantile(x, c(0.05, 0.5, 0.95))
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
  # as posterior 1.8.0: the indicator of the tail beyond the quantile,
  # x <= q below the median and x >= q above it
  vapply(probs, function(p) {
    if (p == 1) p <- (length(x) - 0.5) / length(x)
    q <- stats::quantile(x, p, names = FALSE)
    ind <- if (p > 0.5) x >= q else x <= q
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

#' @rdname diagnostics
#' @export
bayescomp_mcse_quantile <- function(x, probs = c(0.05, 0.95)) {
  # as posterior 1.7.1: the Beta interval of the quantile's rank, mapped
  # back through the sorted draws
  x <- as_draws_matrix_chains(x)
  s <- sort(as.vector(x))
  n <- length(s)
  ess <- bayescomp_ess_quantile(x, probs)
  vapply(seq_along(probs), function(k) {
    e <- ess[k]; p <- probs[k]
    if (is.na(e)) return(NA_real_)
    a <- stats::qbeta(c(0.1586553, 0.8413447), e * p + 1, e * (1 - p) + 1)
    (s[min(ceiling(a[2L] * n), n)] - s[max(floor(a[1L] * n), 1L)]) / 2
  }, 0)
}

#' Gradient test
#'
#' Compares the automatic differentiation gradient of the log density with
#' central finite differences, as CmdStan's `diagnose` method
#' (`test=gradient`) does.
#'
#' @param model A `bayescomp_model`.
#' @param init Point to test: as for [bayescomp_sample()] (a radius, a list
#'   of constrained values, a file name, ...).
#' @param seed Seed for a random initial point.
#' @param epsilon Finite difference step size.
#' @param error Absolute difference above which an element fails.
#' @param jacobian Include the Jacobian adjustment.
#' @return A data frame with one row per unconstrained parameter: `param`
#'   (index from 0, as CmdStan prints it), `value`, `model` (the gradient),
#'   `finite_diff` and `error`, with attribute `"log_prob"` and attribute
#'   `"passed"` (all errors within `error`).
#' @examples
#' m <- bayescomp_model("parameters { real<lower=0> s; } model { s ~ gamma(2, 1); }")
#' bayescomp_check_gradients(m, init = list(s = 1.5))
#' @export
bayescomp_check_gradients <- function(model, init = 2, seed = NULL, epsilon = 1e-6, error = 1e-6, jacobian = TRUE) {
  check_model(model)
  if (!(epsilon > 0) || !(error > 0)) stop("'epsilon' and 'error' must be positive", call. = FALSE)
  u <- with_seed(seed, find_init(model, init, 1L))
  r <- lp_grad(model, u, jacobian, FALSE, print = FALSE)
  lp <- function(x) bayescomp_log_prob(model, x, jacobian = jacobian)
  fd <- vapply(seq_along(u), function(i) {
    e <- replace(numeric(length(u)), i, epsilon)
    (lp(u + e) - lp(u - e)) / (2 * epsilon)
  }, 0)
  out <- data.frame(param = seq_along(u) - 1L, value = u, model = r$grad, finite_diff = fd, error = r$grad - fd)
  structure(out, log_prob = r$value, passed = all(abs(out$error) <= error))
}

#' Multivariate effective sample size
#'
#' The multivariate effective sample size of Vats, Flegal and Jones (2019),
#' `n (det(Lambda) / det(Sigma))^(1 / p)`, where `Lambda` is the sample
#' covariance of the draws and `Sigma` the replicated batch means estimate of
#' the asymptotic covariance of their mean (Vats and Knudson 2021), together
#' with the number of effective draws that estimates every mean to relative
#' precision `epsilon` with confidence `1 - alpha`.
#'
#' @param x A `bayescomp_fit`, or an array of draws (iterations by chains by
#'   variables).
#' @param variables Variables of a fit to use (default all parameters).
#' @param alpha,epsilon Confidence level and relative precision of the
#'   stopping rule.
#' @return A list with `mess`, `min_ess` (the effective sample size the
#'   stopping rule asks for), `enough` (`mess >= min_ess`) and `batch_size`.
#' @references Vats, D., Flegal, J. M. and Jones, G. L. (2019). Multivariate
#'   output analysis for Markov chain Monte Carlo. *Biometrika* 106, 321 to
#'   337. Vats, D. and Knudson, C. (2021). Revisiting the Gelman-Rubin
#'   diagnostic. *Statistical Science* 36, 518 to 529.
#' @examples
#' set.seed(1)
#' x <- array(rnorm(4000), c(500, 4, 2))
#' bayescomp_ess_multi(x)$mess
#' @export
bayescomp_ess_multi <- function(x, variables = NULL, alpha = 0.05, epsilon = 0.05) {
  if (inherits(x, "bayescomp_fit")) {
    x <- bayescomp_draws(x, variables, format = "array")
    if (is.null(variables)) x <- x[, , !grepl("__$", dimnames(x)[[3L]]), drop = FALSE]
  }
  x <- as.array(x)
  if (length(dim(x)) == 2L) x <- array(x, c(dim(x)[1L], 1L, dim(x)[2L]))
  n <- dim(x)[1L]; m <- dim(x)[2L]; p <- dim(x)[3L]
  keep <- vapply(seq_len(p), function(j) !na_diag(x[, , j]), TRUE)
  x <- x[, , keep, drop = FALSE]; p <- sum(keep)
  if (p == 0L || n < 10L) stop("too few draws or no varying variable", call. = FALSE)
  b <- max(1L, floor(sqrt(n)))
  a <- n %/% b
  if (a * m <= p) stop(sprintf("%d batch means for %d variables: the batch means covariance is singular; use more draws or fewer variables", a * m, p), call. = FALSE)
  X <- matrix(x[seq_len(a * b), , , drop = FALSE], ncol = p)
  mu <- colMeans(X)
  # batch means of every chain, centered at the overall mean
  bm <- do.call(rbind, lapply(seq_len(m), function(ch) {
    xc <- matrix(x[seq_len(a * b), ch, , drop = FALSE], ncol = p)
    rowsum(xc, rep(seq_len(a), each = b)) / b
  }))
  bm <- sweep(bm, 2L, mu)
  Sigma <- b * crossprod(bm) / (a * m - 1)
  Lambda <- stats::cov(X)
  ld <- function(S) as.numeric(determinant(S, logarithm = TRUE)$modulus)
  N <- a * b * m
  mess <- N * exp((ld(Lambda) - ld(Sigma)) / p)
  min_ess <- 2^(2 / p) * pi / (p * gamma(p / 2))^(2 / p) * stats::qchisq(1 - alpha, p) / epsilon^2
  list(mess = mess, min_ess = min_ess, enough = mess >= min_ess, batch_size = b)
}

#' Nested R-hat
#'
#' The nested R-hat of Margossian et al. (2024) for many short chains
#' started in groups (superchains) from the same point: it compares the
#' variance between superchain means with the variance within superchains,
#' and is reliable even with one draw per chain.
#'
#' @param x Matrix of draws, iterations by chains (a vector is one draw per
#'   chain).
#' @param superchain Superchain of each chain (a vector with one entry per
#'   column); every superchain needs the same number of chains.
#' @return A number; values near 1 (Margossian et al. suggest below about
#'   1.01 with many chains) indicate the chains have forgotten their starts.
#' @references Margossian, C. C., Hoffman, M. D., Sountsov, P., Riou-Durand,
#'   L., Vehtari, A. and Gelman, A. (2024). Nested R-hat: assessing the
#'   convergence of Markov chain Monte Carlo when running many short chains.
#'   *Bayesian Analysis*.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(10 * 64), 10, 64)
#' bayescomp_rhat_nested(x, rep(1:8, each = 8))
#' @export
bayescomp_rhat_nested <- function(x, superchain) {
  x <- if (is.null(dim(x))) matrix(x, nrow = 1L) else as.matrix(x)
  if (length(superchain) != ncol(x)) stop("'superchain' needs one entry per chain", call. = FALSE)
  g <- split(seq_len(ncol(x)), superchain)
  M <- unique(lengths(g))
  if (length(g) < 2L || length(M) != 1L) stop("need at least two superchains of equal size", call. = FALSE)
  N <- nrow(x)
  cm <- colMeans(x)
  sm <- vapply(g, function(i) mean(cm[i]), 0)
  B <- stats::var(sm)
  W <- mean(vapply(g, function(i) {
    bk <- if (M > 1L) stats::var(cm[i]) else 0
    wk <- if (N > 1L) mean(apply(x[, i, drop = FALSE], 2L, stats::var)) else 0
    bk + wk
  }, 0))
  sqrt(1 + B / W)
}

#' Local R-hat
#'
#' The local R-hat of Moins et al. (2023): the R-hat of the indicator
#' `x <= z`, maximized over thresholds `z` taken at quantiles of the pooled
#' draws (`R-hat_infinity`). It detects chains that agree in location and
#' scale but differ in the shape or the tails of their distributions, which
#' rank normalized R-hat can miss.
#'
#' @param x Matrix of draws, iterations by chains.
#' @param n_grid Number of thresholds.
#' @return A number, `R-hat_infinity`, with attribute `"local"`, a data frame
#'   of thresholds `z` and their R-hat.
#' @references Moins, T., Arbel, J., Dutfoy, A. and Girard, S. (2023). On the
#'   use of a local R-hat to improve MCMC convergence diagnostic. *Bayesian
#'   Analysis*.
#' @examples
#' set.seed(1)
#' x <- cbind(rnorm(500), runif(500, -sqrt(3), sqrt(3)))
#' bayescomp_rhat_local(x)
#' @export
bayescomp_rhat_local <- function(x, n_grid = 100L) {
  x <- as_draws_matrix_chains(x)
  if (ncol(x) < 2L) x <- split_chains(x)
  z <- unique(stats::quantile(x, seq(0, 1, length.out = n_grid + 2L)[-c(1L, n_grid + 2L)], names = FALSE))
  r <- vapply(z, function(zz) { ind <- (x <= zz) + 0; if (na_diag(ind)) NA_real_ else rhat_raw(ind) }, 0)
  structure(max(r, na.rm = TRUE), local = data.frame(z = z, rhat = r))
}

#' Kernel Stein discrepancy
#'
#' The kernel Stein discrepancy of the draws of a fit with the inverse
#' multiquadric kernel `(1 + |x - y|^2)^(-1/2)` on the unconstrained scale
#' (Gorham and Mackey 2017). It needs only the gradient of the log density at
#' the draws, is zero only for the target, and so measures sample quality
#' including bias that convergence diagnostics cannot see. Compare values
#' between samplers or runs of one model; the scale has no absolute meaning.
#'
#' @param fit A `bayescomp_fit` with unconstrained draws.
#' @param max_draws Draws used (thinned evenly); the cost is quadratic.
#' @return A number.
#' @references Gorham, J. and Mackey, L. (2017). Measuring sample quality
#'   with kernels. *Proceedings of the 34th International Conference on
#'   Machine Learning*, 1292 to 1301.
#' @examples
#' \donttest{
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' fit <- bayescomp_sample(m, chains = 1, iter_sampling = 500, seed = 1, refresh = 0)
#' bayescomp_ksd(fit)
#' }
#' @export
bayescomp_ksd <- function(fit, max_draws = 1000L) {
  if (!inherits(fit, "bayescomp_fit") || is.null(fit$unconstrained)) stop("'fit' must be a fit with unconstrained draws", call. = FALSE)
  X <- do.call(rbind, fit$unconstrained)
  if (nrow(X) > max_draws) X <- X[unique(round(seq(1, nrow(X), length.out = max_draws))), , drop = FALSE]
  lpg <- model_lpg(fit$model, jacobian = TRUE, propto = FALSE)
  B <- t(apply(X, 1L, function(u) lpg(u)$grad))
  if (ncol(X) == 1L) B <- matrix(B, ncol = 1L)
  d <- ncol(X); beta <- -0.5
  sq <- rowSums(X^2)
  R2 <- pmax(outer(sq, sq, "+") - 2 * tcrossprod(X), 0)
  S <- 1 + R2
  bx <- rowSums(B * X)
  # b_j . (x_i - x_j) - b_i . (x_i - x_j)
  cross <- (tcrossprod(X, B) - matrix(bx, nrow(X), nrow(X), byrow = TRUE)) -
    (matrix(bx, nrow(X), nrow(X)) - tcrossprod(B, X))
  K0 <- tcrossprod(B) * S^beta + 2 * beta * S^(beta - 1) * cross -
    2 * beta * (d * S^(beta - 1) + 2 * (beta - 1) * S^(beta - 2) * R2)
  sqrt(max(mean(K0), 0))
}
