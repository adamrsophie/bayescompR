# Reverse diffusion Monte Carlo (Huang, Dong, Chen, Zhang, Lin and Zhang
# 2024), experimental. Draws start from a standard normal and follow the
# time reversal of an Ornstein-Uhlenbeck process that carries the target to
# that normal. The score of the diffused density needs no training: it is
# E[z0 | z_t] under p(z0) N(z_t; e^-t z0, 1 - e^-2t), estimated by self
# normalized importance sampling with the proposal N(e^-t z_t, 1 - e^-2t),
# whose weights reduce to p(z0) / N(z0; 0, I). Coordinates are standardized
# with a Gaussian reference first, so the weights are well behaved when the
# target is close to it. The method can cross between modes that trap
# Markov chains, but its draws are approximate: they carry discretization
# error and importance sampling error in the score. Each draw costs
# steps * inner_samples log density evaluations.

#' Reverse diffusion Monte Carlo
#'
#' @param model A `bayescomp_model`.
#' @param num_samples Number of independent draws.
#' @param steps Time steps of the reverse process.
#' @param horizon Diffusion time `T`; the process is started from its
#'   normal limit at `T`.
#' @param inner_samples Importance samples per score estimate.
#' @param reference,reference_scale Standardizing Gaussian, as in
#'   [bayescomp_smc()].
#' @param seed,cores Random seed and number of parallel workers.
#' @return A `bayescomp_fit` with method `"diffusion"` and element `ess`
#'   (mean importance sampling effective size per step).
#' @references Huang, X., Dong, H., Hao, Y., Ma, Y. and Zhang, T. (2024).
#'   Reverse diffusion Monte Carlo. *International Conference on Learning
#'   Representations*.
#' @export
bayescomp_diffusion <- function(model, num_samples = 300L, steps = 40L, horizon = 4, inner_samples = 32L,
                                reference = c("pathfinder", "standard"), reference_scale = 1.5, seed = NULL,
                                cores = getOption("mc.cores", 1L)) {
  check_model(model)
  reference <- match.arg(reference)
  n <- model$dim
  if (n == 0L) stop("the model has no parameters", call. = FALSE)
  S <- as_count(num_samples, "num_samples", 1L)
  Kst <- as_count(steps, "steps", 2L)
  M <- as_count(inner_samples, "inner_samples", 2L)
  if (!is.numeric(horizon) || length(horizon) != 1L || !is.finite(horizon) || horizon <= 0)
    stop("horizon must be a positive number", call. = FALSE)
  t0 <- proc.time()[["elapsed"]]
  info <- NULL
  res <- run_chains(model, 1L, seed, 1L, function(id) {
    lp <- lp_value_fun(model)
    ref <- gaussian_reference(model, reference, reference_scale)
    # log density of the standardized target, up to a constant
    lz <- function(z) lp(ref$mu + as.vector(ref$L %*% z))
    # times from horizon down to a small floor, denser near zero
    tmin <- min(1e-3, horizon / 10)
    ts <- tmin + (horizon - tmin) * (seq(1, 0, length.out = Kst + 1L))^2
    chunks <- split(seq_len(S), ceiling(seq_len(S) * max(1L, min(cores, S)) / S))
    run_chunk <- function(ix) {
      Z <- matrix(stats::rnorm(length(ix) * n), length(ix), n)
      ess_sum <- 0
      for (k in seq_len(Kst)) {
        t <- ts[k]; dt <- ts[k] - ts[k + 1L]
        s2 <- -expm1(-2 * t); a <- exp(-t)
        for (i in seq_len(nrow(Z))) {
          z0 <- sweep(matrix(stats::rnorm(M * n), M, n) * sqrt(s2), 2L, a * Z[i, ], "+")
          lw <- apply(z0, 1L, lz) + 0.5 * rowSums(z0^2)
          lw[!is.finite(lw)] <- -Inf
          if (all(!is.finite(lw))) { e0 <- a * Z[i, ]; ess <- 0 } else {
            w <- exp(lw - max(lw)); w <- w / sum(w)
            e0 <- colSums(z0 * w); ess <- 1 / sum(w^2)
          }
          score <- (a * e0 - Z[i, ]) / s2
          # Euler Maruyama step of dz = (z + 2 score) dt + sqrt(2) dW backwards in time
          Z[i, ] <- Z[i, ] + (Z[i, ] + 2 * score) * dt + sqrt(2 * dt) * stats::rnorm(n)
          ess_sum <- ess_sum + ess
        }
      }
      list(Z = Z, ess = ess_sum / (Kst * nrow(Z)))
    }
    outs <- par_lapply(chunks, run_chunk, cores)
    Z <- do.call(rbind, lapply(outs, `[[`, "Z"))
    qs <- sweep(Z %*% t(ref$L), 2L, ref$mu, "+")
    info <<- list(ess = mean(vapply(outs, `[[`, 0, "ess")))
    dg <- matrix(0, S, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
    dg[, "lp__"] <- apply(qs, 1L, lp)
    bad <- !is.finite(dg[, "lp__"])
    if (any(bad)) stop(sprintf("%d of %d diffusion draws have zero posterior density; increase inner_samples, steps or reference_scale",
                               sum(bad), S), call. = FALSE)
    list(q = qs, gq = chain_write_draws(model, qs, FALSE), diag = dg, step_size = NA_real_, inv_metric = NULL,
         n_warmup_saved = 0L, time = c(warmup = 0, sampling = 0))
  })
  fit <- make_fit_extra(model, res, list(steps = Kst), "diffusion", seed, proc.time()[["elapsed"]] - t0)
  fit$ess <- info$ess
  fit
}
