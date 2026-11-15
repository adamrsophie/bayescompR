# Simulation based calibration (Talts et al. 2018): if the sampler is
# correct, the rank of a true parameter value among posterior draws fitted
# to data simulated from it is uniform.

#' Simulation based calibration
#'
#' Repeats: draw parameters and data from the prior predictive distribution
#' with `simulate`, fit the model to the data, and record the rank of each
#' true parameter value among `max_rank` thinned posterior draws. With a
#' correct sampler and model the ranks are uniform on `0, ..., max_rank`;
#' a smile shape means the posteriors are too narrow, a frown too wide, and
#' a slope a bias.
#'
#' @param code Stan code or a `.stan` file name of the model to fit.
#' @param simulate Function of the simulation number returning a list with
#'   elements `parameters` (named list of true values of the parameters,
#'   with the shapes of their declarations) and `data` (the data list).
#' @param n_sims Number of simulations.
#' @param max_rank Number of posterior draws each rank is computed among;
#'   the draws are thinned evenly from the fit.
#' @param seed Random seed: simulation `i` uses `seed + i`.
#' @param ... Passed to [bayescomp_sample()] (for example `chains`,
#'   `iter_warmup` and `iter_sampling`).
#' @return An object of class `bayescomp_sbc` with elements `ranks` (data
#'   frame with columns `sim`, `variable` and `rank`), `max_rank`, and
#'   `summary` (one row per variable with the p-value of a chi-squared test
#'   of uniform ranks over up to 20 equal bins, and the number of
#'   divergences summed over simulations). `plot()` draws a rank histogram
#'   per variable.
#' @references Talts, S., Betancourt, M., Simpson, D., Vehtari, A. and
#'   Gelman, A. (2018). Validating Bayesian inference algorithms with
#'   simulation-based calibration. arXiv:1804.06788.
#' @examples
#' \donttest{
#' code <- "data { int N; vector[N] y; } parameters { real mu; }
#'   model { mu ~ normal(0, 1); y ~ normal(mu, 1); }"
#' sim <- function(i) {
#'   mu <- rnorm(1)
#'   list(parameters = list(mu = mu), data = list(N = 5, y = rnorm(5, mu)))
#' }
#' s <- bayescomp_sbc(code, sim, n_sims = 20, max_rank = 49, seed = 1,
#'                    chains = 1, iter_warmup = 200, iter_sampling = 490, refresh = 0)
#' s
#' }
#' @export
bayescomp_sbc <- function(code, simulate, n_sims = 100L, max_rank = 99L, seed = 1L, ...) {
  if (!is.function(simulate)) stop("'simulate' must be a function of the simulation number", call. = FALSE)
  n_sims <- as_count(n_sims, "n_sims", 1L)
  max_rank <- as_count(max_rank, "max_rank", 1L)
  rows <- vector("list", n_sims)
  divs <- numeric(n_sims)
  for (i in seq_len(n_sims)) {
    s <- with_seed(seed + i, simulate(i))
    if (!is.list(s) || !is.list(s$parameters) || !is.list(s$data))
      stop("'simulate' must return a list with list elements 'parameters' and 'data'", call. = FALSE)
    m <- bayescomp_model(code, data = s$data, seed = seed + i)
    fit <- bayescomp_sample(m, seed = seed + i, ...)
    divs[i] <- sum(post_warmup(fit$sampler_diagnostics, fit)[, , "divergent__"])
    a <- unclass(fit_draws(fit))
    nm <- bayescomp_param_names(m, FALSE, FALSE)
    truth <- if (all(m$params$pars %in% names(s$parameters))) flatten_draw(s$parameters[m$params$pars])
    if (length(truth) != length(nm)) stop(sprintf("simulation %d: 'parameters' do not match the declared parameters", i), call. = FALSE)
    d <- matrix(a[, , nm, drop = FALSE], ncol = length(nm))
    if (nrow(d) < max_rank) stop(sprintf("the fits have %d draws, fewer than max_rank = %d", nrow(d), max_rank), call. = FALSE)
    d <- d[round(seq(1, nrow(d), length.out = max_rank)), , drop = FALSE]
    rows[[i]] <- data.frame(sim = i, variable = nm, rank = colSums(sweep(d, 2L, truth, `<`)), row.names = NULL)
  }
  ranks <- do.call(rbind, rows)
  structure(list(ranks = ranks, max_rank = max_rank, summary = sbc_summary(ranks, max_rank, sum(divs))),
            class = "bayescomp_sbc")
}

sbc_summary <- function(ranks, max_rank, divergences, bins = min(20L, max_rank + 1L)) {
  vars <- unique(ranks$variable)
  p <- vapply(vars, function(v) {
    r <- ranks$rank[ranks$variable == v]
    counts <- tabulate(floor(r * bins / (max_rank + 1)) + 1L, bins)
    if (length(r) < 2L) return(NA_real_)
    suppressWarnings(stats::chisq.test(counts)$p.value)
  }, 0)
  data.frame(variable = vars, p_value = unname(p), divergences = divergences, row.names = NULL)
}

#' @export
print.bayescomp_sbc <- function(x, ...) {
  cat(sprintf("Simulation based calibration: %d simulations, ranks among %d draws\n",
              length(unique(x$ranks$sim)), x$max_rank))
  print(x$summary, row.names = FALSE, digits = 3)
  invisible(x)
}

#' @export
plot.bayescomp_sbc <- function(x, variable = NULL, bins = min(20L, x$max_rank + 1L), ...) {
  vars <- variable %||% unique(x$ranks$variable)
  old <- graphics::par(mfrow = grDevices::n2mfrow(length(vars)))
  on.exit(graphics::par(old))
  for (v in vars) {
    r <- x$ranks$rank[x$ranks$variable == v]
    graphics::hist(r, breaks = seq(0, x$max_rank + 1, length.out = bins + 1L), right = FALSE, main = v, xlab = "rank", ...)
    graphics::abline(h = length(r) / bins, lty = 2)
  }
  invisible(x)
}
