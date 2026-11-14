# Optimization (posterior mode or maximum likelihood) and the Laplace
# approximation around the mode.

safe_lp <- function(model, jacobian, propto = TRUE) {
  lpg <- model_lpg(model, jacobian, propto)
  function(u) {
    r <- tryCatch(lpg(u), bayescomp_fatal_error = function(e) stop(e),
       bayescomp_domain_error = function(e) NULL,
       bayescomp_reject_error = function(e) NULL)
    if (is.null(r) || !is.finite(r$value)) return(list(value = -Inf, grad = rep(NA_real_, length(u))))
    r
  }
}

# Hessian of the log density by central differences of the gradient
num_hessian <- function(lpg, u, h = 1e-5) {
  n <- length(u)
  H <- matrix(0, n, n)
  for (i in seq_len(n)) {
    e <- numeric(n); e[i] <- h
    H[, i] <- (lpg(u + e)$grad - lpg(u - e)$grad) / (2 * h)
  }
  0.5 * (H + t(H))
}

#' Find a posterior mode or maximum likelihood estimate
#'
#' @param model A `bayescomp_model`.
#' @param algorithm `"lbfgs"`, `"bfgs"` or `"newton"`.
#' @param jacobian Include the Jacobian adjustment (`TRUE` gives the mode of
#'   the posterior on the unconstrained scale, as for Laplace sampling;
#'   `FALSE`, the default, gives the penalized maximum likelihood estimate).
#' @param init Initial values, as in [bayescomp_sample()].
#' @param seed Random seed for initialization.
#' @param iter Maximum number of iterations.
#' @param tol_rel_grad,tol_obj Convergence tolerances passed to the optimizer.
#' @return A `bayescomp_fit` with one draw holding the optimum, and elements
#'   `par` (named constrained values), `value` (the objective) and `code`.
#' @examples
#' m <- bayescomp_model("data { int n; int y; } parameters { real<lower=0,upper=1> p; }
#'   model { y ~ binomial(n, p); }", data = list(n = 10, y = 3))
#' bayescomp_optimize(m)$par
#' @export
bayescomp_optimize <- function(model, algorithm = c("lbfgs", "bfgs", "newton"), jacobian = FALSE,
                               init = 2, seed = NULL, iter = 2000L, tol_rel_grad = 1e7, tol_obj = 1e-12) {
  check_model(model)
  algorithm <- match.arg(algorithm)
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, 1L, seed, 1L, function(id) {
    u0 <- find_init(model, init, 1L)
    lpg <- safe_lp(model, jacobian)
    r <- optimize_u(lpg, u0, algorithm, iter, tol_rel_grad, tol_obj)
    r$gq <- chain_write_draws(model, matrix(r$u, 1L), TRUE)
    r
  })[[1L]]
  out <- list(list(q = matrix(res$u, 1L), gq = res$gq,
                   diag = matrix(c(res$value, rep(0, 6L)), 1L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__"))),
                   step_size = NA_real_, inv_metric = NULL, n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0)))
  fit <- make_fit(model, out, list(algorithm = algorithm, jacobian = jacobian), "optimize", seed, proc.time()[["elapsed"]] - t0)
  fit$par <- fit$draws[1L, 1L, -1L]
  fit$value <- res$value
  fit$code <- res$code
  fit$iterations <- res$iterations
  fit
}

optimize_u <- function(lpg, u0, algorithm, iter, tol_rel_grad, tol_obj) {
  n <- length(u0)
  if (n == 0L) return(list(u = u0, value = lpg(u0)$value, code = 0L, iterations = 0L))
  cache <- new.env(parent = emptyenv())
  fn <- function(u) {
    r <- lpg(u)
    cache$u <- u; cache$r <- r
    if (!is.finite(r$value)) return(1e300)
    -r$value
  }
  gr <- function(u) {
    r <- if (identical(cache$u, u)) cache$r else lpg(u)
    if (!is.finite(r$value)) return(rep(0, length(u)))
    -r$grad
  }
  if (algorithm == "newton") return(newton_opt(lpg, u0, iter))
  method <- if (algorithm == "lbfgs") "L-BFGS-B" else "BFGS"
  ctl <- if (method == "L-BFGS-B") list(maxit = iter, factr = tol_rel_grad, pgtol = 1e-8, lmm = 5L)
         else list(maxit = iter, reltol = tol_obj)
  o <- stats::optim(u0, fn, gr, method = method, control = ctl)
  if (o$convergence == 0L && method == "L-BFGS-B") {
    # polish with BFGS to reach a tight gradient tolerance
    o2 <- stats::optim(o$par, fn, gr, method = "BFGS", control = list(maxit = iter, reltol = tol_obj))
    if (o2$value <= o$value) o <- o2
  }
  list(u = o$par, value = -o$value, code = o$convergence, iterations = o$counts[["function"]])
}

newton_opt <- function(lpg, u, iter) {
  r <- lpg(u)
  for (k in seq_len(iter)) {
    H <- num_hessian(lpg, u)
    e <- eigen(H, symmetric = TRUE)
    # make the Hessian negative definite
    ev <- -pmax(abs(e$values), 1e-8)
    step <- -as.vector(e$vectors %*% ((t(e$vectors) %*% r$grad) / ev))
    a <- 1
    repeat {
      un <- u + a * step
      rn <- lpg(un)
      if (is.finite(rn$value) && rn$value >= r$value) break
      a <- a / 2
      # no ascent possible: converged only if the gradient vanishes
      if (a < 1e-10) return(list(u = u, value = r$value, code = if (max(abs(r$grad)) < 1e-6) 0L else 2L, iterations = k))
    }
    done <- abs(rn$value - r$value) < 1e-10 * (1 + abs(r$value))
    u <- un; r <- rn
    if (done) return(list(u = u, value = r$value, code = 0L, iterations = k))
  }
  list(u = u, value = r$value, code = 1L, iterations = iter)
}

#' Laplace approximation
#'
#' Draws from a normal approximation to the posterior on the unconstrained
#' scale, centered at the mode with covariance equal to the inverse negative
#' Hessian of the log density (as Stan's `laplace` method).
#'
#' @param model A `bayescomp_model`.
#' @param mode Optional result of [bayescomp_optimize()] (with
#'   `jacobian = TRUE`) or a named list of constrained values at the mode.
#' @param jacobian Use the Jacobian adjusted density (recommended).
#' @param draws Number of draws.
#' @param seed Random seed.
#' @return A `bayescomp_fit` with extra columns `log_p__` and `log_g__`.
#' @export
bayescomp_laplace <- function(model, mode = NULL, jacobian = TRUE, draws = 1000L, seed = NULL) {
  check_model(model)
  draws <- as_count(draws, "draws", 1L)
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, 1L, seed, 1L, function(id) {
    lpg <- safe_lp(model, jacobian)
    u <- if (is.null(mode)) {
      op <- optimize_u(lpg, find_init(model, 2, 1L), "lbfgs", 2000L, 1e7, 1e-12)
      if (op$code != 0L) warning(sprintf("the optimizer did not converge (code %d); the Laplace approximation may be centered away from the mode", op$code), call. = FALSE)
      op$u
    } else if (inherits(mode, "bayescomp_fit")) {
      as.vector(mode$unconstrained[[1L]][1L, ])
    } else {
      bayescomp_unconstrain(model, mode)
    }
    H <- num_hessian(lpg, u)
    S <- tryCatch(solve(-H), error = function(e) stop("the Hessian at the mode is singular", call. = FALSE))
    S <- 0.5 * (S + t(S))
    L <- tryCatch(t(chol(S)), error = function(e) stop("the negative Hessian at the mode is not positive definite", call. = FALSE))
    n <- length(u)
    z <- matrix(stats::rnorm(draws * n), draws, n)
    qs <- sweep(z %*% t(L), 2L, u, "+")
    log_g <- -0.5 * rowSums(z^2) - 0.5 * n * log(2 * pi) - sum(log(diag(L)))
    log_p <- apply(qs, 1L, function(q) lpg(q)$value)
    gq <- chain_write_draws(model, qs, FALSE)
    dg <- matrix(0, draws, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- 0
    list(q = qs, gq = cbind(log_p__ = log_p, log_g__ = log_g, gq), diag = dg, step_size = NA_real_,
         inv_metric = S, n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  fit <- make_fit_extra(model, res, list(jacobian = jacobian), "laplace", seed, proc.time()[["elapsed"]] - t0)
  fit
}

# like make_fit but allows extra leading columns in gq (log_p__, log_g__)
make_fit_extra <- function(model, res, cfg, method, seed, elapsed) {
  r <- res[[1L]]
  nm <- colnames(r$gq)
  n <- nrow(r$gq)
  draws <- array(NA_real_, c(n, 1L, 1L + length(nm)), dimnames = list(iteration = NULL, chain = NULL, variable = c("lp__", nm)))
  draws[, 1L, 1L] <- r$diag[, "lp__"]
  draws[, 1L, -1L] <- r$gq
  structure(list(draws = draws, sampler_diagnostics = NULL, method = method, model = model,
                 n_warmup_saved = 0L, step_size = NA_real_, inv_metric = list(r$inv_metric),
                 unconstrained = list(r$q), time = list(total = elapsed), config = cfg, seed = seed),
            class = "bayescomp_fit")
}
