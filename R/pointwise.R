# Pointwise log likelihood and posterior predictive draws taken from the
# sampling statements of the model block, so loo and predictive checks work
# without a generated quantities block.

#' Pointwise log likelihood and posterior predictive draws
#'
#' Evaluates the model block at each draw and records, for every sampling
#' statement (`y ~ dist(...)`) and every `target += dist_lpdf(y | ...)` whose
#' left side is data, the normalized log density of each element of `y` and,
#' for `bayescomp_posterior_predict()`, a draw from the matching `_rng`
#' function. Univariate statements contribute one column per element
#' (`y[1]`, `y[2]`, ...); multivariate ones one column per statement.
#' Truncated statements are not supported.
#'
#' @param fit A `bayescomp_fit` with unconstrained draws.
#' @param format `"matrix"` for a draws by observations matrix (chains
#'   stacked), `"array"` for an iterations by chains by observations array.
#' @param seed Seed for the predictive draws.
#' @return A matrix or array with one column (slice) per observation.
#'   `bayescomp_posterior_predict()` gives `NA` for distributions without an
#'   `_rng` function.
#' @examples
#' m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; }
#'   model { y ~ normal(mu, 1); }", data = list(N = 5, y = c(-1, 0, 0.5, 1, 2)))
#' fit <- bayescomp_sample(m, chains = 1, iter_warmup = 100, iter_sampling = 100,
#'                         seed = 1, refresh = 0)
#' dim(bayescomp_log_lik(fit))
#' dim(bayescomp_posterior_predict(fit, seed = 2))
#' @export
bayescomp_log_lik <- function(fit, format = c("matrix", "array")) {
  pointwise_out(pointwise_eval(fit, FALSE), match.arg(format))
}

#' @rdname bayescomp_log_lik
#' @export
bayescomp_posterior_predict <- function(fit, format = c("matrix", "array"), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pointwise_out(pointwise_eval(fit, TRUE), match.arg(format))
}

pointwise_out <- function(a, format) {
  if (format == "array") return(a)
  m <- matrix(a, ncol = dim(a)[3L])
  colnames(m) <- dimnames(a)[[3L]]
  m
}

pointwise_eval <- function(fit, rng) {
  if (!inherits(fit, "bayescomp_fit") || is.null(fit$unconstrained))
    stop("'fit' must be a bayescomp_fit with unconstrained draws", call. = FALSE)
  model <- fit$model
  f <- model_fun(model, "log_prob")
  w <- fit$n_warmup_saved %||% 0L
  qs <- lapply(fit$unconstrained, function(q) q[setdiff(seq_len(nrow(q)), seq_len(w)), , drop = FALSE])
  old <- .ad$tape
  .ad$tape <- NULL
  on.exit(.ad$tape <- old)
  one <- function(u) {
    ctx <- new_ctx(jacobian = FALSE, print = FALSE, propto = FALSE)
    pw <- new.env(parent = emptyenv())
    pw$ll <- list(); pw$yrep <- list(); pw$names <- list(); pw$count <- list(); pw$rng <- rng
    ctx$pw <- pw
    f(as.numeric(u), ctx)
    list(ll = unlist(pw$ll), yrep = unlist(pw$yrep), names = unlist(pw$names))
  }
  first <- one(qs[[1L]][1L, ])
  if (!length(first$names))
    stop("the model block has no sampling statement on data to take pointwise terms from", call. = FALSE)
  nm <- make.unique(first$names)
  n <- nrow(qs[[1L]])
  # only the requested output is kept: the log likelihood, or with `rng` the
  # predictive draws
  key <- if (rng) "yrep" else "ll"
  out <- array(NA_real_, c(n, length(qs), length(nm)), dimnames = list(NULL, NULL, nm))
  for (c in seq_along(qs)) for (i in seq_len(n)) {
    r <- if (c == 1L && i == 1L) first else one(qs[[c]][i, ])
    out[i, c, ] <- r[[key]]
  }
  out
}
