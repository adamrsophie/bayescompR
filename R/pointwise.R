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
#' (`y[1]`, `y[2]`, ...), with arguments such as cutpoints or a vector of
#' log odds shared by all elements, and GLMs one column per row of `x`, as
#' the equivalent unfused statement does. A multivariate statement is one
#' joint observation, so it gives one log likelihood column, while its
#' predictive draws have one column per element of `y`. Truncated
#' statements are renormalized for each observation, and their predictive
#' draws are taken by rejection. A user defined density is one joint
#' observation per statement (one per iteration inside a loop) and has no
#' predictive draws. `group` sums the log likelihood of observations that
#' belong together, for leave one group out cross validation.
#'
#' @param fit A `bayescomp_fit` with unconstrained draws.
#' @param format `"matrix"` for a draws by observations matrix (chains
#'   stacked), `"array"` for an iterations by chains by observations array.
#' @param group Optional vector with one entry per observation column; the
#'   columns of each group are summed and named by the group.
#' @param seed Seed for the predictive draws.
#' @return A matrix or array with one column (slice) per observation or
#'   group. `bayescomp_posterior_predict()` gives `NA` for distributions
#'   without an `_rng` function.
#' @examples
#' m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; }
#'   model { y ~ normal(mu, 1); }", data = list(N = 5, y = c(-1, 0, 0.5, 1, 2)))
#' fit <- bayescomp_sample(m, chains = 1, iter_warmup = 100, iter_sampling = 100,
#'                         seed = 1, refresh = 0)
#' dim(bayescomp_log_lik(fit))
#' dim(bayescomp_posterior_predict(fit, seed = 2))
#' @export
bayescomp_log_lik <- function(fit, format = c("matrix", "array"), group = NULL) {
  a <- pointwise_eval(fit, FALSE)
  if (!is.null(group)) {
    if (length(group) != dim(a)[3L]) stop(sprintf("'group' must have one entry per observation (%d)", dim(a)[3L]), call. = FALSE)
    g <- factor(group, levels = unique(group))
    m <- t(rowsum(t(matrix(a, ncol = dim(a)[3L])), g, reorder = FALSE))
    a <- array(m, c(dim(a)[1:2], nlevels(g)), dimnames = list(NULL, NULL, levels(g)))
  }
  pointwise_out(a, match.arg(format))
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
    pw$out <- list(); pw$names <- list(); pw$count <- list(); pw$rng <- rng
    ctx$pw <- pw
    f(as.numeric(u), ctx)
    list(out = unlist(pw$out), names = unlist(pw$names))
  }
  first <- one(qs[[1L]][1L, ])
  if (!length(first$names))
    stop("the model block has no sampling statement on data to take pointwise terms from", call. = FALSE)
  nm <- make.unique(first$names)
  n <- nrow(qs[[1L]])
  out <- array(NA_real_, c(n, length(qs), length(nm)), dimnames = list(NULL, NULL, nm))
  for (c in seq_along(qs)) for (i in seq_len(n)) {
    r <- if (c == 1L && i == 1L) first else one(qs[[c]][i, ])
    out[i, c, ] <- r$out
  }
  out
}

# Pointwise terms of one density statement on data (see cg_pointwise()).
# `spec$el` marks the elementwise arguments (NULL: the statement is one joint
# observation). Only the requested output is computed: the normalized log
# density of each observation, or with `ctx$pw$rng` a draw from the _rng
# function. The statement itself runs after this and checks the arguments.
rt_pointwise <- function(ctx, label, spec, uf, tr, y, ...) {
  if (!is.null(spec$unsupported)) stop_unsupported(sprintf("pointwise terms of the %s statement on '%s' are not available", spec$unsupported, label))
  vals <- lapply(list(y, ...), val)
  if (!is.null(spec$glm)) vals <- spec$glm(vals)
  full <- uf %||% spec$full
  pw <- ctx$pw
  rng <- isTRUE(pw$rng)
  el <- spec$el
  if (is.null(el)) {
    vy <- vals[[1L]]
    need <- length(unlist(vy))
    if (rng) {
      # one draw with the shape of the left side: a shared mean gives one
      # vector per draw, so array left sides take one draw per row
      draw <- function() as.numeric(unlist(do.call(spec$rng, vals[-1L])))
      r <- if (is.null(spec$rng)) rep(NA_real_, need) else draw()
      if (length(r) < need && need %% length(r) == 0L) r <- as.vector(do.call(rbind, c(list(r), lapply(seq_len(need / length(r) - 1L), function(i) draw()))))
      dims <- dim(vy)
      pw$out[[length(pw$out) + 1L]] <- r
      pw$names[[length(pw$names) + 1L]] <- if (length(r) == 1L) pw_index(pw, label, 1L)
        else if (length(dims) > 1L && prod(dims) == length(r)) sprintf("%s[%s]", label, apply(arrayInd(seq_along(r), dims), 1L, paste, collapse = ","))
        else sprintf("%s[%d]", label, seq_along(r))
    } else {
      pw$out[[length(pw$out) + 1L]] <- as.numeric(do.call(full, vals))
      pw$names[[length(pw$names) + 1L]] <- if (need == 1L) pw_index(pw, label, 1L) else label
    }
    return(invisible())
  }
  n <- max(lengths(vals[el]))
  pw$names[[length(pw$names) + 1L]] <- pw_index(pw, label, n)
  rows <- spec$rows
  at <- function(v, i, e, r) {
    if (r) return(if (is.matrix(v)) v[if (nrow(v) == 1L) 1L else i, ] else v)
    if (!e) v else if (length(v) == 1L) v[[1L]] else v[[i]]
  }
  args_i <- function(args, els, rs, i) Map(at, args, i, els, rs)
  each <- function(f, args, els, rs) vapply(seq_len(n), function(i) as.numeric(do.call(f, args_i(args, els, rs, i))), 0)
  pw$out[[length(pw$out) + 1L]] <- if (rng) {
    if (is.null(spec$rng)) rep(NA_real_, n)
    else if (!is.null(tr)) pw_trunc_rng(spec$rng, tr, n, function(i) args_i(vals[-1L], el[-1L], rows[-1L], i))
    else if (spec$vrng && length(vals) > 1L) {
      # the _rng functions vectorize like R's generators: one draw per element
      r <- as.numeric(do.call(spec$rng, lapply(vals[-1L], function(v) rep_len(as.vector(v), n))))
      if (length(r) == n) r else each(spec$rng, vals[-1L], el[-1L], rows[-1L])
    } else each(spec$rng, vals[-1L], el[-1L], rows[-1L])
  } else {
    ll <- if (!is.null(spec$kern)) as.numeric(rep_len(do.call(spec$kern, lapply(vals, as.vector)), n))
          else each(full, vals, el, rows)
    # a truncated density is renormalized over [lo, hi] for each observation
    # (data outside the bounds make the log density -Inf, so no draw has them)
    if (!is.null(tr)) ll <- ll - vapply(seq_len(n), function(i)
      as.numeric(val(trunc_norm(args_i(vals[-1L], el[-1L], rows[-1L], i), tr$lo, tr$hi, tr$lcdf, tr$lccdf, tr$discrete, tr$fam))), 0)
    ll
  }
  invisible()
}

# names label[k0 + 1], ..., label[k0 + n], counting the terms on `label`
pw_index <- function(pw, label, n) {
  k0 <- pw$count[[label]] %||% 0L
  pw$count[[label]] <- k0 + n
  sprintf("%s[%d]", label, k0 + seq_len(n))
}

# predictive draws of a truncated statement by rejection: draws outside
# [lo, hi] are redrawn up to 1000 times, then NA
pw_trunc_rng <- function(rng, tr, n, args) {
  lo <- if (is.null(tr$lo)) -Inf else val(tr$lo)
  hi <- if (is.null(tr$hi)) Inf else val(tr$hi)
  vapply(seq_len(n), function(i) {
    a <- args(i)
    for (k in 1:1000) {
      r <- as.numeric(do.call(rng, a))
      if (r >= lo && r <= hi) return(r)
    }
    NA_real_
  }, 0)
}

# Fused GLMs as their univariate family on the linear predictor: each row of
# x is one observation. `el` overrides the elementwise arguments of the family.
pw_eta <- function(x, a, b) as.vector(val(glm_linpred("pointwise", x, a, b)))
pw_glm <- list(
  normal_id_glm = list(family = "normal", convert = function(v) list(v[[1L]], pw_eta(v[[2L]], v[[3L]], v[[4L]]), v[[5L]])),
  bernoulli_logit_glm = list(family = "bernoulli_logit", convert = function(v) list(v[[1L]], pw_eta(v[[2L]], v[[3L]], v[[4L]]))),
  poisson_log_glm = list(family = "poisson_log", convert = function(v) list(v[[1L]], pw_eta(v[[2L]], v[[3L]], v[[4L]]))),
  neg_binomial_2_log_glm = list(family = "neg_binomial_2_log", convert = function(v) list(v[[1L]], pw_eta(v[[2L]], v[[3L]], v[[4L]]), v[[5L]])),
  binomial_logit_glm = list(family = "binomial_logit", convert = function(v) list(v[[1L]], v[[2L]], pw_eta(v[[3L]], v[[4L]], v[[5L]])),
                            rng = function(N, eta) .lib$binomial_rng[[1L]]$impl(N, stats::plogis(eta))),
  ordered_logistic_glm = list(family = "ordered_logistic", convert = function(v) list(v[[1L]], pw_eta(v[[2L]], 0, v[[3L]]), v[[4L]])),
  categorical_logit_glm = list(family = "categorical_logit", el = c(TRUE, TRUE), convert = function(v) {
    x <- v[[2L]]
    if (is.null(dim(x))) x <- matrix(x, 1L)
    eta <- sweep(x %*% v[[4L]], 2L, as.vector(v[[3L]]), "+")
    list(v[[1L]], lapply(seq_len(nrow(eta)), function(i) eta[i, ]))
  }))

# Vectorized normalized log densities for the common univariate families,
# used instead of one density call per observation.
pw_kernels <- list(
  normal = function(y, mu, s) stats::dnorm(y, mu, s, log = TRUE),
  student_t = function(y, nu, mu, s) stats::dt((y - mu) / s, nu, log = TRUE) - log(s),
  cauchy = function(y, mu, s) stats::dcauchy(y, mu, s, log = TRUE),
  logistic = function(y, mu, s) stats::dlogis(y, mu, s, log = TRUE),
  lognormal = function(y, mu, s) stats::dlnorm(y, mu, s, log = TRUE),
  double_exponential = function(y, mu, s) -log(2 * s) - abs(y - mu) / s,
  exponential = function(y, b) stats::dexp(y, b, log = TRUE),
  gamma = function(y, a, b) stats::dgamma(y, a, rate = b, log = TRUE),
  weibull = function(y, a, s) stats::dweibull(y, a, s, log = TRUE),
  beta = function(y, a, b) stats::dbeta(y, a, b, log = TRUE),
  poisson = function(y, l) stats::dpois(y, l, log = TRUE),
  poisson_log = function(y, eta) stats::dpois(y, exp(eta), log = TRUE),
  bernoulli = function(y, p) stats::dbinom(y, 1L, p, log = TRUE),
  bernoulli_logit = function(y, eta) stats::plogis((2 * y - 1) * eta, log.p = TRUE),
  binomial = function(y, N, p) stats::dbinom(y, N, p, log = TRUE),
  binomial_logit = function(y, N, eta) lchoose(N, y) + y * stats::plogis(eta, log.p = TRUE) + (N - y) * stats::plogis(-eta, log.p = TRUE),
  neg_binomial_2 = function(y, mu, phi) stats::dnbinom(y, size = phi, mu = mu, log = TRUE),
  neg_binomial_2_log = function(y, eta, phi) stats::dnbinom(y, size = phi, mu = exp(eta), log = TRUE))
