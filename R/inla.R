# Integrated nested Laplace approximation (Rue, Martino and Chopin 2009) for
# models with a latent block that is close to Gaussian given the remaining
# hyperparameters. For each hyperparameter value theta the latent block is
# replaced by its Laplace approximation, which gives
#   log p(theta | y) ~ log p(x*, theta, y) - log det(H(theta)) / 2 + const,
# with x* the conditional mode and H the negative Hessian there. The
# hyperparameter posterior is explored on a grid in the eigenbasis of its
# curvature at the mode (the INLA grid strategy), and the latent posterior is
# the mixture of the conditional Gaussians over the grid. All of this happens
# on the unconstrained scale.

#' Integrated nested Laplace approximation
#'
#' @param model A `bayescomp_model`.
#' @param latent Names of the parameters forming the latent block; the
#'   others are hyperparameters (at most 4).
#' @param dz Grid step in standardized hyperparameter units; by default
#'   0.25, 0.5, 0.75 or 1 for one to four hyperparameters.
#' @param threshold Keep grid points whose log density is within this of
#'   the mode.
#' @param output_samples Number of joint draws returned.
#' @param seed Random seed.
#' @return A `bayescomp_fit` with method `"inla"` and element `grid`: the
#'   hyperparameter grid points (unconstrained), their normalized weights
#'   and the log density.
#' @references Rue, H., Martino, S. and Chopin, N. (2009). Approximate
#'   Bayesian inference for latent Gaussian models by using integrated nested
#'   Laplace approximations. *Journal of the Royal Statistical Society B* 71,
#'   319 to 392.
#' @export
bayescomp_inla <- function(model, latent, dz = NULL, threshold = 5, output_samples = 2000L, seed = NULL) {
  check_model(model)
  n <- model$dim
  bad <- setdiff(latent, model$params$pars)
  if (length(bad)) stop(sprintf("'%s' is not a parameter", bad[1L]), call. = FALSE)
  ix <- unlist(lapply(latent, function(p) { d <- model$decls[[p]]; d$offset + seq_len(d$free_size) }))
  it <- setdiff(seq_len(n), ix)
  if (length(it) > 4L) stop("INLA supports at most 4 hyperparameters; move more parameters into 'latent'", call. = FALSE)
  S <- as_count(output_samples, "output_samples", 1L)
  t0 <- proc.time()[["elapsed"]]
  info <- NULL
  res <- run_chains(model, 1L, seed, 1L, function(id) {
    lpg <- safe_lp(model, jacobian = TRUE, propto = FALSE)
    u0 <- find_init(model, 2, id)
    x_prev <- u0[ix]
    cond <- function(th, x0 = x_prev) {
      u <- u0; u[it] <- th
      f <- function(x) { u[ix] <- x; lpg(u) }
      x <- x0
      for (k in 1:50) {
        r <- f(x)
        if (!is.finite(r$value)) return(NULL)
        H <- num_hessian_V(function(z) -f(z)$grad[ix], x)
        R <- tryCatch(chol(H), error = function(e) NULL)
        step <- if (is.null(R)) r$grad[ix] * 0.1 else backsolve(R, forwardsolve(t(R), r$grad[ix]))
        # damped Newton: halve until the log density does not decrease
        a <- 1
        repeat {
          xn <- x + a * step
          vn <- f(xn)$value
          if (is.finite(vn) && vn >= r$value - 1e-10) break
          a <- a / 2
          if (a < 1e-6) { xn <- x; break }
        }
        conv <- max(abs(xn - x)) < 1e-8
        x <- xn
        if (conv) break
      }
      r <- f(x)
      H <- num_hessian_V(function(z) -f(z)$grad[ix], x)
      R <- tryCatch(chol(H), error = function(e) NULL)
      if (is.null(R)) return(NULL)
      list(x = x, R = R, lp = r$value + 0.5 * length(ix) * log(2 * pi) - sum(log(diag(R))))
    }
    if (!length(it)) {
      c0 <- cond(numeric(0))
      grid <- list(list(th = numeric(0), c = c0, lp = c0$lp))
    } else {
      obj <- function(th) { c <- cond(th); if (is.null(c)) return(1e100); x_prev <<- c$x; -c$lp }
      o <- stats::optim(u0[it], obj, method = if (length(it) == 1L) "BFGS" else "Nelder-Mead",
                        control = list(maxit = 2000, reltol = 1e-10))
      th0 <- o$par
      m <- length(it)
      if (is.null(dz)) dz <- c(0.25, 0.5, 0.75, 1)[m]
      Hh <- matrix(0, m, m); h <- 1e-3
      f0 <- -o$value
      fv <- function(d) { c <- cond(th0 + d); if (is.null(c)) NA_real_ else c$lp }
      for (i in seq_len(m)) for (j in i:m) {
        ei <- replace(numeric(m), i, h); ej <- replace(numeric(m), j, h)
        Hh[i, j] <- Hh[j, i] <- -(fv(ei + ej) - fv(ei - ej) - fv(-ei + ej) + fv(-ei - ej)) / (4 * h^2)
      }
      e <- eigen(Hh, symmetric = TRUE)
      lam <- pmax(e$values, 1e-8)
      B <- e$vectors %*% diag(1 / sqrt(lam), m)
      kmax <- ceiling(1.5 * sqrt(2 * threshold) / dz)
      zs <- as.matrix(do.call(expand.grid, rep(list(seq(-kmax, kmax) * dz), m)))
      zs <- zs[rowSums(zs^2) <= (1.5^2) * 2 * threshold * 1.0001, , drop = FALSE]
      grid <- list()
      x_prev <- cond(th0)$x
      for (r in seq_len(nrow(zs))) {
        th <- th0 + as.vector(B %*% zs[r, ])
        c <- cond(th)
        if (is.null(c) || f0 - c$lp > threshold) next
        grid[[length(grid) + 1L]] <- list(th = th, c = c, lp = c$lp)
      }
    }
    lpv <- vapply(grid, `[[`, 0, "lp")
    w <- exp(lpv - max(lpv)); w <- w / sum(w)
    pick <- sample.int(length(grid), S, replace = TRUE, prob = w)
    qs <- matrix(NA_real_, S, n)
    for (s in seq_len(S)) {
      g <- grid[[pick[s]]]
      qs[s, it] <- g$th
      qs[s, ix] <- g$c$x + backsolve(g$c$R, stats::rnorm(length(ix)))
    }
    info <<- data.frame(do.call(rbind, lapply(grid, function(g) g$th)) %||% matrix(numeric(0), length(grid), 0),
                        weight = w, log_density = lpv)
    full <- lp_value_fun(model)
    dg <- matrix(0, S, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- apply(qs, 1L, full)
    list(q = qs, gq = chain_write_draws(model, qs, FALSE), diag = dg, step_size = NA_real_, inv_metric = NULL,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  fit <- make_fit_extra(model, res, list(latent = latent), "inla", seed, proc.time()[["elapsed"]] - t0)
  fit$grid <- info
  fit
}
