# Blocked Gibbs sampling with automatic conjugacy detection.
#
# The unconstrained parameters are split into blocks (by default one per
# declared parameter), and each block is updated in turn with its own
# kernel while the others are held fixed. The kernels are
#
# * "conjugate": an independence Metropolis-Hastings step whose proposal is
#   the exact full conditional when the conditional belongs to a recognized
#   family, so it is a Gibbs draw accepted with probability one. The
#   families are detected numerically from the log density:
#     - Gaussian for unconstrained blocks whose log density is quadratic
#       (normal priors and likelihoods in the block, linear models);
#     - Gamma for positive scalars x = exp(u) whose log density is
#       A u + B exp(k u) for k in {1, -1, 2, -2} (gamma precisions, inverse
#       gamma variances, Poisson and exponential rates);
#     - Beta for unit interval scalars x = inv_logit(u) whose log density is
#       A log x + B log(1 - x) (binomial and Bernoulli probabilities).
#   The proposal is fitted at the current state and the reverse proposal
#   at the proposed state, so the step stays exact even when the detection
#   is wrong; the acceptance rate then drops below one.
# * "slice": univariate slice sampling for each coordinate, widths adapted
#   during warmup.
# * "metropolis": random walk Metropolis on the block, covariance and scale
#   adapted during warmup.
# * "mala": Metropolis adjusted Langevin on the block.
# * "auto" (default): "conjugate" when a family is detected at several
#   random states during initialization, else "slice".

#' Blocked Gibbs sampling
#'
#' @param model A `bayescomp_model`.
#' @param blocks `NULL` (one block per parameter, kernel `"auto"`) or a list
#'   of blocks, each a list with `pars` (parameter names) and optionally
#'   `kernel` (`"auto"`, `"conjugate"`, `"slice"`, `"metropolis"` or
#'   `"mala"`).
#' @param chains,iter_warmup,iter_sampling,thin,init,seed,cores,refresh As in
#'   [bayescomp_sample()].
#' @return A `bayescomp_fit` with method `"gibbs"` and element `blocks`
#'   describing each block's kernel, detected family and acceptance rate.
#' @examples
#' \donttest{
#' m <- bayescomp_model("data { int N; vector[N] y; }
#'   parameters { real mu; real<lower=0> tau; }
#'   model { mu ~ normal(0, 10); tau ~ gamma(2, 1); y ~ normal(mu, inv_sqrt(tau)); }",
#'   data = list(N = 5, y = c(1, 2, 0.5, 1.5, 2.5)))
#' fit <- bayescomp_gibbs(m, chains = 2, iter_warmup = 200, iter_sampling = 500, seed = 1)
#' fit$blocks
#' }
#' @export
bayescomp_gibbs <- function(model, blocks = NULL, chains = 4L, iter_warmup = 1000L, iter_sampling = 1000L,
                            thin = 1L, init = 2, seed = NULL, cores = 1L, refresh = 0L) {
  check_model(model)
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  chains <- as_count(chains, "chains", 1L)
  nw <- as_count(iter_warmup, "iter_warmup")
  ns <- as_count(iter_sampling, "iter_sampling", 1L)
  thin <- as_count(thin, "thin", 1L)
  bl <- gibbs_blocks(model, blocks)
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, chains, seed, cores, function(id) {
    lpf <- lp_value_fun(model)
    lpg <- safe_lp(model, jacobian = TRUE, propto = FALSE)
    q <- find_init(model, init, id)
    bs <- lapply(bl, function(b) gibbs_block_init(b, q, lpf))
    lp <- lpf(q)
    keep <- seq_len(ns)[(seq_len(ns) - 1L) %% thin == 0L]
    qs <- matrix(NA_real_, length(keep), n)
    dg <- matrix(0, length(keep), 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "energy__"] <- NA_real_
    row <- 0L
    for (it in seq_len(nw + ns)) {
      warm <- it <= nw
      acc_it <- numeric(length(bs))
      for (j in seq_along(bs)) {
        r <- gibbs_update(bs[[j]], q, lp, lpf, lpg, warm, it)
        q <- r$q; lp <- r$lp; bs[[j]] <- r$b
        acc_it[j] <- r$acc
        if (!warm) { bs[[j]]$acc_sum <- bs[[j]]$acc_sum + r$acc; bs[[j]]$acc_n <- bs[[j]]$acc_n + 1 }
      }
      if (!warm && ((it - nw - 1L) %% thin == 0L)) {
        row <- row + 1L
        qs[row, ] <- q
        dg[row, c("lp__", "accept_stat__")] <- c(lp, mean(acc_it))
      }
      if (refresh > 0L && it %% refresh == 0L)
        message(sprintf("Chain %d: Iteration %d / %d", id, it, nw + ns))
    }
    info <- data.frame(block = vapply(bs, `[[`, "", "label"), kernel = vapply(bs, `[[`, "", "kernel"),
                       family = vapply(bs, function(b) b$family %||% NA_character_, ""),
                       acceptance = vapply(bs, function(b) b$acc_sum / max(1, b$acc_n), 0),
                       stringsAsFactors = FALSE)
    list(q = qs, gq = chain_write_draws(model, qs, FALSE), diag = dg, step_size = NA_real_, inv_metric = NULL,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0), info = info)
  })
  fit <- make_fit(model, res, list(algorithm = "gibbs"), "gibbs", seed, proc.time()[["elapsed"]] - t0)
  fit$blocks <- res[[1L]]$info
  if (chains > 1L) fit$blocks$acceptance <- rowMeans(matrix(vapply(res, function(r) r$info$acceptance, numeric(nrow(res[[1L]]$info))),
                                                            nrow(res[[1L]]$info)))
  fit
}

gibbs_blocks <- function(model, blocks) {
  idx_of <- function(p) {
    d <- model$decls[[p]]
    if (is.null(d) || !(p %in% model$params$pars)) stop(sprintf("'%s' is not a parameter", p), call. = FALSE)
    d$offset + seq_len(d$free_size)
  }
  if (is.null(blocks)) blocks <- lapply(model$params$pars, function(p) list(pars = p))
  out <- lapply(blocks, function(b) {
    if (is.character(b)) b <- list(pars = b)
    ix <- unlist(lapply(b$pars, idx_of))
    list(label = paste(b$pars, collapse = ", "), idx = ix, kernel = b$kernel %||% "auto")
  })
  out <- out[vapply(out, function(b) length(b$idx) > 0L, TRUE)]
  all_idx <- unlist(lapply(out, `[[`, "idx"))
  if (anyDuplicated(all_idx)) stop("a parameter appears in more than one block", call. = FALSE)
  miss <- setdiff(seq_len(model$dim), all_idx)
  if (length(miss)) stop("every parameter must belong to a block", call. = FALSE)
  for (b in out) if (!(b$kernel %in% c("auto", "conjugate", "slice", "metropolis", "mala")))
    stop(sprintf("unknown kernel '%s'", b$kernel), call. = FALSE)
  out
}

gibbs_block_init <- function(b, q, lpf) {
  k <- length(b$idx)
  b$acc_sum <- 0; b$acc_n <- 0
  b$width <- rep(2, k)
  b$scale <- 2.38 / sqrt(k)
  b$step <- 0.5
  b$est <- est_new(k, FALSE)
  b$sd <- rep(1, k)
  if (b$kernel %in% c("auto", "conjugate")) {
    fam <- gibbs_detect_family(b, q, lpf)
    if (is.null(fam)) {
      if (b$kernel == "conjugate") warning(sprintf("no conjugate family detected for block '%s'; the independence proposal may be inefficient", b$label), call. = FALSE)
      if (b$kernel == "auto") b$kernel <- "slice"
      else b$family <- "gaussian"
    } else {
      b$kernel <- "conjugate"
      b$family <- fam
    }
  }
  b
}

# ---- conditional families -----------------------------------------------------------

set_block <- function(q, idx, v) { q[idx] <- v; q }

# fit a family to the block conditional at state q; returns a proposal
# list(draw(), logd(v)) or NULL when the fit fails
gibbs_fit <- function(family, b, q, lpf) {
  idx <- b$idx
  if (family == "gaussian") {
    k <- length(idx)
    f <- function(v) lpf(set_block(q, idx, v))
    v0 <- q[idx]
    f0 <- f(v0)
    if (!is.finite(f0)) return(NULL)
    # exact quadratic fit from function values: Hessian by second differences
    h <- 1
    H <- matrix(0, k, k)
    g <- numeric(k)
    fp <- numeric(k); fm <- numeric(k)
    for (i in seq_len(k)) {
      e <- replace(numeric(k), i, h)
      fp[i] <- f(v0 + e); fm[i] <- f(v0 - e)
      g[i] <- (fp[i] - fm[i]) / (2 * h)
      H[i, i] <- (fp[i] - 2 * f0 + fm[i]) / h^2
    }
    if (k > 1L) for (i in 1:(k - 1L)) for (j in (i + 1L):k) {
      ei <- replace(numeric(k), i, h); ej <- replace(numeric(k), j, h)
      H[i, j] <- H[j, i] <- (f(v0 + ei + ej) - fp[i] - fp[j] + f0) / h^2
    }
    if (any(!is.finite(H)) || any(!is.finite(g))) return(NULL)
    P <- -H
    R <- tryCatch(chol(P), error = function(e) NULL)
    if (is.null(R)) return(NULL)
    m <- v0 + backsolve(R, forwardsolve(t(R), g))
    return(list(draw = function() m + backsolve(R, stats::rnorm(k)),
                logd = function(v) sum(log(diag(R))) - 0.5 * sum((R %*% (v - m))^2)))
  }
  u0 <- q[idx]
  f <- function(u) lpf(set_block(q, idx, u))
  if (family %in% c("gamma1", "gamma-1", "gamma2", "gamma-2")) {
    kk <- as.numeric(sub("gamma", "", family))
    f0 <- f(u0)
    fs <- c(f(u0 - 1), f(u0 + 1))
    if (any(!is.finite(c(f0, fs)))) return(NULL)
    X <- rbind(c(-1, exp(kk * (u0 - 1)) - exp(kk * u0)), c(1, exp(kk * (u0 + 1)) - exp(kk * u0)))
    co <- tryCatch(solve(X, fs - f0), error = function(e) NULL)
    if (is.null(co)) return(NULL)
    A <- co[1L]; B <- co[2L]
    shape <- A / kk; rate <- -B
    if (!is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0) return(NULL)
    return(list(draw = function() log(stats::rgamma(1, shape, rate)) / kk,
                logd = function(u) { t <- exp(kk * u); stats::dgamma(t, shape, rate, log = TRUE) + log(abs(kk)) + kk * u }))
  }
  if (family == "beta") {
    us <- u0 + c(-1, 1)
    fs <- c(f(us[1L]), f(us[2L]))
    f0 <- f(u0)
    if (any(!is.finite(c(fs, f0)))) return(NULL)
    lx <- function(u) c(stats::plogis(u, log.p = TRUE), stats::plogis(-u, log.p = TRUE))
    X <- rbind(lx(us[1L]) - lx(u0), lx(us[2L]) - lx(u0))
    co <- tryCatch(solve(X, fs - f0), error = function(e) NULL)
    if (is.null(co) || any(!is.finite(co)) || any(co <= 0)) return(NULL)
    a <- co[1L]; bb <- co[2L]
    return(list(draw = function() stats::qlogis(stats::rbeta(1, a, bb)),
                logd = function(u) a * stats::plogis(u, log.p = TRUE) + bb * stats::plogis(-u, log.p = TRUE) - lbeta(a, bb)))
  }
  NULL
}

# a family is accepted when its fit predicts the log density at fresh points
gibbs_detect_family <- function(b, q, lpf) {
  k <- length(b$idx)
  fams <- if (k == 1L) c("gaussian", "gamma1", "gamma-1", "gamma2", "gamma-2", "beta") else "gaussian"
  for (fam in fams) {
    ok <- TRUE
    for (trial in 1:3) {
      qq <- q
      qq[b$idx] <- q[b$idx] + stats::rnorm(k, 0, 0.5)
      fit <- gibbs_fit(fam, b, qq, lpf)
      if (is.null(fit)) { ok <- FALSE; break }
      for (check in 1:2) {
        v1 <- qq[b$idx] + stats::rnorm(k, 0, 0.7)
        v2 <- qq[b$idx] + stats::rnorm(k, 0, 0.7)
        d_true <- lpf(set_block(qq, b$idx, v1)) - lpf(set_block(qq, b$idx, v2))
        d_fit <- fit$logd(v1) - fit$logd(v2)
        if (!is.finite(d_true) || !is.finite(d_fit) || abs(d_true - d_fit) > 1e-6 * (1 + abs(d_true))) { ok <- FALSE; break }
      }
      if (!ok) break
    }
    if (ok) return(fam)
  }
  NULL
}

# ---- one block update -----------------------------------------------------------------

gibbs_update <- function(b, q, lp, lpf, lpg, warm, it) {
  idx <- b$idx
  k <- length(idx)
  if (b$kernel == "conjugate") {
    fx <- gibbs_fit(b$family, b, q, lpf)
    if (is.null(fx)) return(list(q = q, lp = lp, b = b, acc = 0))
    v <- fx$draw()
    qn <- set_block(q, idx, v)
    lpn <- lpf(qn)
    # once warmup has shown that every proposal was accepted with
    # probability one, the fitted family is the exact conditional and the
    # reverse fit is skipped: a plain Gibbs draw
    if (isTRUE(b$exact) && !warm && is.finite(lpn)) return(list(q = qn, lp = lpn, b = b, acc = 1))
    fy <- gibbs_fit(b$family, b, qn, lpf)
    acc <- 0
    if (is.finite(lpn) && !is.null(fy)) {
      la <- lpn - lp + fy$logd(q[idx]) - fx$logd(v)
      acc <- if (is.nan(la)) 0 else min(1, exp(la))
    }
    if (stats::runif(1) < acc) { q <- qn; lp <- lpn }
    if (warm) {
      b$warm_n <- (b$warm_n %||% 0) + 1
      b$warm_min <- min(b$warm_min %||% 1, acc)
      b$exact <- b$warm_n >= 20 && b$warm_min > 1 - 1e-8
    }
    return(list(q = q, lp = lp, b = b, acc = acc))
  }
  if (b$kernel == "slice") {
    for (i in seq_len(k)) {
      j <- idx[i]
      logy <- lp - stats::rexp(1)
      w <- b$width[i]
      L <- q[j] - w * stats::runif(1); R <- L + w
      f <- function(x) lpf(replace(q, j, x))
      J <- 31L
      while (J > 0L && f(L) > logy) { L <- L - w; J <- J - 1L }
      J <- 31L
      while (J > 0L && f(R) > logy) { R <- R + w; J <- J - 1L }
      repeat {
        x <- stats::runif(1, L, R)
        fx <- f(x)
        if (fx > logy) { q[j] <- x; lp <- fx; break }
        if (x < q[j]) L <- x else R <- x
        if (R - L < 1e-12) break
      }
    }
    if (warm) {
      est_add(b$est, q[idx])
      if (b$est$n > 10) b$width <- 2 * sqrt(pmax(est_value(b$est), 1e-8))
    }
    return(list(q = q, lp = lp, b = b, acc = 1))
  }
  if (b$kernel == "metropolis") {
    v <- q[idx] + b$scale * b$sd * stats::rnorm(k)
    qn <- set_block(q, idx, v)
    lpn <- lpf(qn)
    acc <- if (is.finite(lpn)) min(1, exp(lpn - lp)) else 0
    if (stats::runif(1) < acc) { q <- qn; lp <- lpn }
    if (warm) {
      b$scale <- b$scale * exp((acc - 0.3) / sqrt(it))
      est_add(b$est, q[idx])
      if (b$est$n > 20) b$sd <- sqrt(pmax(est_value(b$est), 1e-10))
    }
    return(list(q = q, lp = lp, b = b, acc = acc))
  }
  # mala on the block
  e <- b$step
  r0 <- lpg(q)
  m0 <- q[idx] + 0.5 * e^2 * b$sd^2 * r0$grad[idx]
  v <- m0 + e * b$sd * stats::rnorm(k)
  qn <- set_block(q, idx, v)
  r1 <- lpg(qn)
  acc <- 0
  if (is.finite(r1$value) && is.finite(r0$value)) {
    m1 <- v + 0.5 * e^2 * b$sd^2 * r1$grad[idx]
    la <- r1$value - r0$value - 0.5 * sum(((q[idx] - m1) / (e * b$sd))^2) + 0.5 * sum(((v - m0) / (e * b$sd))^2)
    acc <- if (is.nan(la)) 0 else min(1, exp(la))
  }
  if (stats::runif(1) < acc) { q <- qn; lp <- lpf(qn) }
  if (warm) {
    b$step <- b$step * exp((acc - 0.574) / sqrt(it))
    est_add(b$est, q[idx])
    if (b$est$n > 20) b$sd <- sqrt(pmax(est_value(b$est), 1e-10))
  }
  list(q = q, lp = lp, b = b, acc = acc)
}
