# Methods that split the likelihood across observations: stochastic
# gradient MCMC (SGLD, Welling and Teh 2011; SGHMC, Chen, Fox and Guestrin
# 2014) and expectation propagation (Vehtari et al. 2020). The model is
# re-instantiated on subsets of its data; by default every data element
# whose length (or number of rows) equals the observation count `N` is
# subset, which covers the usual `int N; vector[N] y; matrix[N, K] X;`
# layout. Supply `subset` for anything else.

model_with_data <- function(model, data) {
  if (!identical(model$backend, "r")) {
    # recompile through the same backend so implementations never mix
    out <- bayescomp_model(model$code, data = data, backend = model$backend,
                           allow_fallback = any(unlist(model$provenance) == "r"))
  } else {
    inst <- instantiate(model$compiled, data)
    out <- model
    out$env <- inst$env; out$decls <- inst$decls; out$data <- data; out$params <- inst$params; out$dim <- inst$dim
  }
  if (out$dim != model$dim)
    stop("the parameter dimension changes with the data subset; supply 'subset' or use a model whose parameter sizes do not depend on N",
         call. = FALSE)
  out
}

default_subset <- function(n_name) {
  function(data, idx) {
    N <- data[[n_name]]
    for (nm in names(data)) {
      v <- data[[nm]]
      if (identical(nm, n_name)) next
      if (is.matrix(v) || is.array(v)) {
        if (dim(v)[1L] == N) {
          d <- dim(v)
          v <- array(v, c(d[1L], prod(d[-1L])))[idx, , drop = FALSE]
          dim(v) <- c(length(idx), d[-1L])
          data[[nm]] <- v
        }
      } else if (length(v) == N && N != 1L) data[[nm]] <- v[idx]
    }
    data[[n_name]] <- length(idx)
    data
  }
}

split_setup <- function(model, n_name, subset) {
  N <- model$data[[n_name]]
  if (is.null(N) || length(N) != 1L) stop(sprintf("the data has no observation count '%s'", n_name), call. = FALSE)
  if (is.null(subset)) subset <- default_subset(n_name)
  prior <- model_with_data(model, subset(model$data, integer(0)))
  list(N = as.integer(N), subset = subset, prior = prior)
}

#' Stochastic gradient MCMC
#'
#' Stochastic gradient Langevin dynamics (`"sgld"`) and stochastic gradient
#' Hamiltonian Monte Carlo (`"sghmc"`) use minibatch estimates of the
#' log likelihood gradient. They have no Metropolis correction, so the draws
#' are approximate with a bias that shrinks with the step size; they pay off
#' only when the data are large and the per observation cost dominates.
#'
#' @param model A `bayescomp_model` with an observation count in its data.
#' @param method `"sgld"` or `"sghmc"`.
#' @param batch_size Observations per minibatch.
#' @param step_size Learning rate (the squared step of the discretized
#'   dynamics); by default `0.01` divided by the largest eigenvalue of the
#'   negative Hessian at the mode.
#' @param friction Friction for SGHMC.
#' @param n_name Name of the observation count in the data.
#' @param subset Optional `function(data, idx)` returning the data for the
#'   observations `idx` (with `idx = integer(0)` for the prior).
#' @param chains,iter_warmup,iter_sampling,thin,init,seed,cores As in
#'   [bayescomp_sample()].
#' @return A `bayescomp_fit` with method `"sgmcmc"`.
#' @references Welling, M. and Teh, Y. W. (2011). Bayesian learning via
#'   stochastic gradient Langevin dynamics. *ICML*. Chen, T., Fox, E. and
#'   Guestrin, C. (2014). Stochastic gradient Hamiltonian Monte Carlo.
#'   *ICML*.
#' @export
bayescomp_sgmcmc <- function(model, method = c("sgld", "sghmc"), batch_size = 32L, step_size = NULL,
                             friction = 0.1, n_name = "N", subset = NULL, chains = 4L, iter_warmup = 500L,
                             iter_sampling = 2000L, thin = 1L, init = 2, seed = NULL, cores = 1L) {
  check_model(model)
  method <- match.arg(method)
  n <- model$dim
  sp <- split_setup(model, n_name, subset)
  b <- min(as_count(batch_size, "batch_size", 1L), sp$N)
  chains <- as_count(chains, "chains", 1L)
  nw <- as_count(iter_warmup, "iter_warmup"); ns <- as_count(iter_sampling, "iter_sampling", 1L)
  thin <- as_count(thin, "thin", 1L)
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, chains, seed, cores, function(id) {
    # ponytail: fixed random partition into minibatches, visited in random
    # order each epoch; one model instance per batch is built once
    perm <- sample.int(sp$N)
    batches <- split(perm, ceiling(seq_along(perm) / b))
    pg <- safe_lp(sp$prior, jacobian = TRUE)
    bg <- lapply(batches, function(ix) safe_lp(model_with_data(model, sp$subset(model$data, ix)), jacobian = TRUE))
    grad_est <- function(u, k) {
      g0 <- pg(u)$grad
      g0 + (sp$N / length(batches[[k]])) * (bg[[k]](u)$grad - g0)
    }
    full <- safe_lp(model, jacobian = TRUE)
    u <- find_init(model, init, id)
    eps <- step_size
    if (is.null(eps)) {
      # scale from the curvature at a nearby mode
      o <- tryCatch(stats::optim(u, function(x) -full(x)$value, function(x) -full(x)$grad, method = "BFGS",
                                 control = list(maxit = 200)), error = function(e) NULL)
      if (!is.null(o)) u <- o$par
      H <- num_hessian_V(function(x) -full(x)$grad, u)
      lam <- tryCatch(max(eigen(H, symmetric = TRUE, only.values = TRUE)$values), error = function(e) NA)
      # learning rate eta = h^2 for a leapfrog style step h = 0.1 / sqrt(lam)
      eps <- if (is.finite(lam) && lam > 0) 0.01 / lam else 1e-6
    }
    keep <- seq_len(ns)[(seq_len(ns) - 1L) %% thin == 0L]
    qs <- matrix(NA_real_, length(keep), n)
    dg <- matrix(0, length(keep), 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "energy__"] <- NA_real_
    v <- numeric(n)
    order_k <- sample.int(length(batches)); pos <- 0L
    row <- 0L
    for (it in seq_len(nw + ns)) {
      pos <- pos + 1L
      if (pos > length(order_k)) { order_k <- sample.int(length(batches)); pos <- 1L }
      g <- grad_est(u, order_k[pos])
      # a non finite minibatch gradient skips the move but keeps the draw
      if (all(is.finite(g))) {
        if (method == "sgld") {
          u <- u + 0.5 * eps * g + sqrt(eps) * stats::rnorm(n)
        } else {
          # Chen et al. (2014), equation 15, with learning rate eta = eps
          v <- (1 - friction) * v + eps * g + sqrt(2 * friction * eps) * stats::rnorm(n)
          u <- u + v
        }
      }
      if (it > nw && ((it - nw - 1L) %% thin == 0L)) {
        row <- row + 1L
        qs[row, ] <- u
        dg[row, c("lp__", "stepsize__")] <- c(full(u)$value, eps)
      }
    }
    list(q = qs, gq = chain_write_draws(model, qs, FALSE), diag = dg, step_size = eps, inv_metric = NULL,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  make_fit(model, res, list(algorithm = method, batch_size = b), "sgmcmc", seed, proc.time()[["elapsed"]] - t0)
}

#' Expectation propagation
#'
#' Expectation propagation with the data split into `sites` partitions
#' (Vehtari et al. 2020). The global approximation is a Gaussian on the
#' unconstrained space, the product of one Gaussian site per partition.
#' Each tilted distribution, the cavity times the partition likelihood times
#' the prior raised to `1 / sites`, is summarized by its Laplace
#' approximation, and the site is updated with damping. Partitions can be
#' processed in parallel with `cores`.
#'
#' @param model A `bayescomp_model` with an observation count in its data.
#' @param sites Number of data partitions.
#' @param iter Maximum number of EP sweeps.
#' @param damping Damping of the site updates, in (0, 1].
#' @param tol Stop when the global mean changes by less than this.
#' @param n_name,subset As in [bayescomp_sgmcmc()].
#' @param output_samples Draws from the final Gaussian approximation.
#' @param seed,cores Random seed and number of parallel workers.
#' @return A `bayescomp_fit` with method `"ep"`; elements `mean` and `cov`
#'   (unconstrained space) and `iterations`.
#' @references Vehtari, A., Gelman, A., Sivula, T. and others (2020).
#'   Expectation propagation as a way of life. *Journal of Machine Learning
#'   Research* 21, 1 to 53.
#' @export
bayescomp_ep <- function(model, sites = 4L, iter = 50L, damping = 0.5, tol = 1e-3, n_name = "N", subset = NULL,
                         output_samples = 1000L, seed = NULL, cores = 1L) {
  check_model(model)
  n <- model$dim
  sp <- split_setup(model, n_name, subset)
  K <- min(as_count(sites, "sites", 1L), sp$N)
  S <- as_count(output_samples, "output_samples", 1L)
  t0 <- proc.time()[["elapsed"]]
  info <- NULL
  res <- run_chains(model, 1L, seed, 1L, function(id) {
    parts <- split(seq_len(sp$N), ceiling(seq_len(sp$N) * K / sp$N))
    pr <- safe_lp(sp$prior, jacobian = TRUE)
    lk <- lapply(parts, function(ix) safe_lp(model_with_data(model, sp$subset(model$data, ix)), jacobian = TRUE))
    # natural parameters: precision Q and shift r = Q m
    Qs <- rep(list(diag(1e-2 / K, n)), K); rs <- rep(list(numeric(n)), K)
    Q <- Reduce(`+`, Qs); r <- Reduce(`+`, rs)
    m_old <- solve(Q, r)
    Q_old <- Q
    it_done <- 0L
    tilted <- function(k, Qc, rc) {
      f <- function(u) {
        a <- lk[[k]](u); p <- pr(u)
        # partition likelihood plus a 1/K share of the prior
        v <- a$value - p$value + p$value / K - 0.5 * sum(u * (Qc %*% u)) + sum(rc * u)
        g <- a$grad - p$grad + p$grad / K - as.vector(Qc %*% u) + rc
        list(value = v, grad = g)
      }
      st <- tryCatch(solve(Qc, rc), error = function(e) m_old)
      o <- stats::optim(st, function(x) { v <- f(x)$value; if (is.finite(v)) -v else 1e100 },
                        function(x) -f(x)$grad, method = "BFGS", control = list(maxit = 500, reltol = 1e-12))
      H <- num_hessian_V(function(x) -f(x)$grad, o$par)
      list(m = o$par, P = H)
    }
    for (sweep in seq_len(iter)) {
      it_done <- sweep
      upd <- par_lapply(seq_len(K), function(k) {
        Qc <- Q - Qs[[k]]; rc <- r - rs[[k]]
        tl <- tryCatch(tilted(k, Qc, rc), error = function(e) NULL)
        if (is.null(tl) || inherits(tryCatch(chol(tl$P), error = function(e) e), "error")) return(NULL)
        list(Q = tl$P - Qc, r = as.vector(tl$P %*% tl$m) - rc)
      }, cores)
      for (k in seq_len(K)) if (!is.null(upd[[k]])) {
        Qn <- (1 - damping) * Qs[[k]] + damping * upd[[k]]$Q
        rn <- (1 - damping) * rs[[k]] + damping * upd[[k]]$r
        if (!inherits(tryCatch(chol(Q - Qs[[k]] + Qn), error = function(e) e), "error")) { Qs[[k]] <- Qn; rs[[k]] <- rn }
      }
      Q <- Reduce(`+`, Qs); r <- Reduce(`+`, rs)
      m <- solve(Q, r)
      dQ <- max(abs(Q - Q_old)) / max(abs(Q))
      Q_old <- Q
      if (max(abs(m - m_old)) < tol && dQ < tol) { m_old <- m; break }
      m_old <- m
    }
    Sig <- solve(Q); Sig <- 0.5 * (Sig + t(Sig))
    L <- t(chol(Sig))
    qs <- sweep(matrix(stats::rnorm(S * n), S, n) %*% t(L), 2L, m_old, "+")
    full <- lp_value_fun(model)
    lp <- apply(qs, 1L, full)
    info <<- list(mean = m_old, cov = Sig, iterations = it_done)
    dg <- matrix(0, S, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- lp
    list(q = qs, gq = chain_write_draws(model, qs, FALSE), diag = dg, step_size = NA_real_, inv_metric = Sig,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  fit <- make_fit_extra(model, res, list(sites = K), "ep", seed, proc.time()[["elapsed"]] - t0)
  fit$mean <- info$mean; fit$cov <- info$cov; fit$iterations <- info$iterations
  fit
}

par_lapply <- function(X, f, cores = 1L) {
  if (cores > 1L && length(X) > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(X, f, mc.cores = min(cores, length(X)))
  else lapply(X, f)
}
