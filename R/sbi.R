# Neural posterior estimation (Papamakarios and Murray 2016; Greenberg,
# Nonnenmacher and Macke 2019) for simulation based inference: when the
# likelihood is only available as a simulator, a mixture density network
# q(theta | x) is trained on prior predictive pairs (theta, x) and then
# evaluated at the observed summaries. The prior comes from a
# bayescomp_model (its parameters and their priors; the model block holds
# only the priors), and the network models the unconstrained parameters.
# The network has one tanh hidden layer and is trained with Adam on
# hand-derived gradients; it is amortized, so new observations need no
# retraining (see the returned `posterior` function).

mdn_init <- function(dx, p, K, h) {
  list(W1 = matrix(stats::rnorm(dx * h, 0, 1 / sqrt(dx)), dx, h), b1 = numeric(h),
       Wa = matrix(stats::rnorm(h * K, 0, 0.1 / sqrt(h)), h, K), ba = numeric(K),
       Wm = matrix(stats::rnorm(h * K * p, 0, 1 / sqrt(h)), h, K * p), bm = stats::rnorm(K * p, 0, 0.5),
       Ws = matrix(stats::rnorm(h * K * p, 0, 0.1 / sqrt(h)), h, K * p), bs = numeric(K * p))
}

mdn_forward <- function(w, X) {
  Hd <- tanh(sweep(X %*% w$W1, 2L, w$b1, "+"))
  list(H = Hd, A = sweep(Hd %*% w$Wa, 2L, w$ba, "+"), Mu = sweep(Hd %*% w$Wm, 2L, w$bm, "+"),
       Ls = pmax(pmin(sweep(Hd %*% w$Ws, 2L, w$bs, "+"), 5), -7))
}

row_lse <- function(M) { m <- apply(M, 1L, max); m + log(rowSums(exp(M - m))) }

# mean log density of T given X and its gradient with respect to the weights
mdn_loglik <- function(w, X, T, grad = TRUE) {
  B <- nrow(X); p <- ncol(T); K <- ncol(w$Wa)
  f <- mdn_forward(w, X)
  logpi <- f$A - row_lse(f$A)
  logN <- matrix(0, B, K)
  Zs <- vector("list", K)
  for (k in seq_len(K)) {
    cols <- (k - 1L) * p + seq_len(p)
    Z <- (T - f$Mu[, cols, drop = FALSE]) / exp(f$Ls[, cols, drop = FALSE])
    Zs[[k]] <- Z
    logN[, k] <- rowSums(-0.5 * Z^2 - f$Ls[, cols, drop = FALSE]) - 0.5 * p * log(2 * pi)
  }
  L <- row_lse(logpi + logN)
  if (!grad) return(list(value = mean(L)))
  R <- exp(logpi + logN - L)
  dA <- R - exp(logpi)
  dMu <- matrix(0, B, K * p); dLs <- matrix(0, B, K * p)
  for (k in seq_len(K)) {
    cols <- (k - 1L) * p + seq_len(p)
    sd <- exp(f$Ls[, cols, drop = FALSE])
    dMu[, cols] <- R[, k] * Zs[[k]] / sd
    dLs[, cols] <- R[, k] * (Zs[[k]]^2 - 1)
  }
  clip <- f$Ls > -7 & f$Ls < 5
  dLs <- dLs * clip
  dH <- dA %*% t(w$Wa) + dMu %*% t(w$Wm) + dLs %*% t(w$Ws)
  dZ1 <- dH * (1 - f$H^2)
  g <- list(W1 = crossprod(X, dZ1), b1 = colSums(dZ1), Wa = crossprod(f$H, dA), ba = colSums(dA),
            Wm = crossprod(f$H, dMu), bm = colSums(dMu), Ws = crossprod(f$H, dLs), bs = colSums(dLs))
  list(value = mean(L), grad = lapply(g, function(v) v / B))
}

mdn_sample <- function(w, x, S) {
  p <- length(w$bm) / length(w$ba); K <- length(w$ba)
  f <- mdn_forward(w, matrix(x, 1L))
  pr <- exp(f$A[1L, ] - max(f$A[1L, ])); pr <- pr / sum(pr)
  k <- sample.int(K, S, replace = TRUE, prob = pr)
  out <- matrix(NA_real_, S, p)
  for (j in seq_len(K)) {
    s <- which(k == j)
    if (!length(s)) next
    cols <- (j - 1L) * p + seq_len(p)
    out[s, ] <- sweep(sweep(matrix(stats::rnorm(length(s) * p), length(s), p), 2L, exp(f$Ls[1L, cols]), "*"), 2L, f$Mu[1L, cols], "+")
  }
  out
}

#' Neural posterior estimation
#'
#' @param model A `bayescomp_model` holding the prior: the parameters to
#'   infer and their priors, without the likelihood.
#' @param simulator `function(theta)` returning a numeric vector of summary
#'   statistics simulated at `theta`, a named vector of constrained
#'   parameters (flat Stan names).
#' @param observed The observed summary statistics.
#' @param num_simulations Number of prior predictive simulations.
#' @param components Mixture components of the density network.
#' @param hidden Hidden units.
#' @param epochs Maximum training epochs; training stops early when the
#'   held out loss stops improving.
#' @param batch_size,learning_rate Adam settings.
#' @param output_samples Posterior draws at `observed`.
#' @param seed Random seed.
#' @return A `bayescomp_fit` with method `"npe"`; element `posterior` is a
#'   function of new summaries returning unconstrained draws, and `loss`
#'   the training and validation trace.
#' @references Papamakarios, G. and Murray, I. (2016). Fast epsilon free
#'   inference of simulation models with Bayesian conditional density
#'   estimation. *NeurIPS*. Cranmer, K., Brehmer, J. and Louppe, G. (2020).
#'   The frontier of simulation based inference. *PNAS* 117, 30055 to 30062.
#' @examples
#' \donttest{
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
#' sim <- function(th) mean(rnorm(10, th[["mu"]], 1))
#' fit <- bayescomp_npe(m, sim, observed = 1, num_simulations = 500, epochs = 50, seed = 1)
#' mean(bayescomp_extract(fit, "mu"))
#' }
#' @export
bayescomp_npe <- function(model, simulator, observed, num_simulations = 2000L, components = 3L, hidden = 32L,
                          epochs = 300L, batch_size = 100L, learning_rate = 3e-3, output_samples = 1000L,
                          seed = NULL) {
  check_model(model)
  if (!is.function(simulator)) stop("'simulator' must be a function", call. = FALSE)
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  Ns <- as_count(num_simulations, "num_simulations", 20L)
  K <- as_count(components, "components", 1L); h <- as_count(hidden, "hidden", 1L)
  S <- as_count(output_samples, "output_samples", 1L)
  t0 <- proc.time()[["elapsed"]]
  # separate streams for the prior draws and for the simulator
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max - 1L, 1L)
  pr <- bayescomp_sample(model, chains = 4L, iter_warmup = 500L, iter_sampling = ceiling(Ns / 4) * 2L, thin = 2L,
                         seed = seed, refresh = 0L)
  sim_seed <- as.integer(seed) %% (.Machine$integer.max - 1L) + 1L
  U <- do.call(rbind, pr$unconstrained)[seq_len(Ns), , drop = FALSE]
  info <- NULL
  res <- run_chains(model, 1L, sim_seed, 1L, function(id) {
    X <- t(vapply(seq_len(Ns), function(i) as.numeric(simulator(ssm_theta(model, U[i, ]))), as.numeric(simulator(ssm_theta(model, U[1L, ])))))
    if (length(observed) == 1L) X <- matrix(X, ncol = 1L)
    ok <- apply(is.finite(X), 1L, all)
    X <- X[ok, , drop = FALSE]; Uk <- U[ok, , drop = FALSE]
    if (ncol(X) != length(observed)) stop("the simulator output and 'observed' differ in length", call. = FALSE)
    xm <- colMeans(X); xs <- pmax(apply(X, 2L, stats::sd), 1e-12)
    um <- colMeans(Uk); us <- pmax(apply(Uk, 2L, stats::sd), 1e-12)
    Xs <- sweep(sweep(X, 2L, xm), 2L, xs, "/"); Ts <- sweep(sweep(Uk, 2L, um), 2L, us, "/")
    if (nrow(Xs) < 2L) stop("fewer than two simulations are finite; check the simulator", call. = FALSE)
    nv <- max(1L, floor(0.1 * nrow(Xs)))
    vi <- sample.int(nrow(Xs), nv); tr <- setdiff(seq_len(nrow(Xs)), vi)
    w <- mdn_init(ncol(Xs), n, K, h)
    m1 <- lapply(w, function(v) v * 0); m2 <- m1
    best <- w; best_v <- -Inf; patience <- 0L; step <- 0L
    trace <- matrix(NA_real_, 0, 2, dimnames = list(NULL, c("train", "validation")))
    for (ep in seq_len(epochs)) {
      ord <- tr[sample.int(length(tr))]
      bs <- split(ord, ceiling(seq_along(ord) / batch_size))
      tl <- 0
      for (b in bs) {
        step <- step + 1L
        r <- mdn_loglik(w, Xs[b, , drop = FALSE], Ts[b, , drop = FALSE])
        tl <- tl + r$value * length(b)
        for (nm in names(w)) {
          g <- r$grad[[nm]]
          m1[[nm]] <- 0.9 * m1[[nm]] + 0.1 * g
          m2[[nm]] <- 0.999 * m2[[nm]] + 0.001 * g^2
          w[[nm]] <- w[[nm]] + learning_rate * (m1[[nm]] / (1 - 0.9^step)) / (sqrt(m2[[nm]] / (1 - 0.999^step)) + 1e-8)
        }
      }
      vl <- mdn_loglik(w, Xs[vi, , drop = FALSE], Ts[vi, , drop = FALSE], grad = FALSE)$value
      trace <- rbind(trace, c(tl / length(tr), vl))
      if (vl > best_v) { best_v <- vl; best <- w; patience <- 0L } else patience <- patience + 1L
      if (patience >= 20L) break
    }
    post <- function(x, S) {
      z <- (as.numeric(x) - xm) / xs
      sweep(sweep(mdn_sample(best, z, S), 2L, us, "*"), 2L, um, "+")
    }
    qs <- post(observed, S)
    info <<- list(post = post, loss = trace)
    lp <- lp_value_fun(model)
    dg <- matrix(0, S, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- apply(qs, 1L, lp)
    list(q = qs, gq = chain_write_draws(model, qs, FALSE), diag = dg, step_size = NA_real_, inv_metric = NULL,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  fit <- make_fit_extra(model, res, list(num_simulations = Ns), "npe", seed, proc.time()[["elapsed"]] - t0)
  fit$posterior <- function(x, n = 1000L) info$post(x, n)
  fit$loss <- info$loss
  fit
}
