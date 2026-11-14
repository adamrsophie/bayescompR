# Probability distributions with Stan's parameterizations.
#
# Each distribution registers <name>_lpdf (or _lpmf), the unnormalized
# <name>_lupdf (or _lupmf), and where available _cdf, _lcdf, _lccdf and _rng.
# Densities are written with differentiable operations, so gradients come
# from the AD tape. With propto = TRUE (sampling statements and _lupdf),
# terms that do not depend on any differentiable argument are dropped, which
# matches Stan's treatment of constants.
#
# Cumulative distribution functions without closed form partial derivatives
# (incomplete gamma and beta shapes, Student t degrees of freedom, discrete
# distributions) use central finite differences for those partials only; the
# function support matrix lists them as "numerical partials".

# include a term?
inc <- function(propto, ...) !propto || any_ad(...)

# sum of a broadcast term over N elements
bsum <- function(term, N) {
  if (length(val(term)) == 1L) N * term else ad_sum(term)
}

# One tape node for a whole density with analytic partials (the analogue of
# Stan's operands_and_partials). `parts()` returns one partial derivative
# array per operand, of length 1 or N; contributions are summed down to the
# operand's shape. Only differentiable operands receive adjoints.
fused <- function(value, ops, parts, N) {
  keep <- which(vapply(ops, is_ad, TRUE))
  if (!length(keep)) return(value)
  ad_node(value, ops[keep], function(g) {
    pp <- parts()
    lapply(keep, function(k) {
      d <- pp[[k]]
      vx <- ops[[k]]$v
      if (length(vx) == 1L) return(g * (if (length(d) == 1L) N * d else sum(d)))
      if (length(d) == 1L) d <- rep(d, length(vx))
      d <- g * d
      dim(d) <- dim(vx)
      d
    })
  })
}

# c * x with the convention 0 * x = 0 (also for infinite x), recycled
zmul <- function(c, x) {
  r <- c * x
  if (length(c) == 1L) { if (c == 0) r[] <- 0 } else r[c == 0] <- 0
  r
}

# register the density family for a distribution
defdist <- function(name, args, lpdf = NULL, discrete = FALSE, cdf = NULL, lcdf = NULL,
                    lccdf = NULL, rng = NULL, rng_ret = NULL, ret = "real", densities = TRUE) {
  suf <- if (discrete) c("_lpmf", "_lupmf") else c("_lpdf", "_lupdf")
  # the cumulative functions share the density's parameter checks (values
  # only, no AD); the variable itself is left to the cumulative function,
  # whose support handling differs from the density's
  if (!is.null(lpdf)) {
    dens_name <- paste0(name, suf[1L])
    checked <- function(fun, fname) {
      force(fun)
      structure(function(...) {
        a <- lapply(list(...), val)
        ok <- function(y) tryCatch({ do.call(lpdf, c(list(y), a[-1L], list(propto = TRUE))); TRUE },
                                   bayescomp_domain_error = function(e) e)
        r <- ok(a[[1L]])
        if (!isTRUE(r)) {
          # the parameters are valid if the density accepts some variable
          # in its support (the variable itself is the cumulative function's
          # business: it may lie outside the support)
          pv <- unlist(lapply(a[-1L], function(x) as.vector(x)[1L]))
          cand <- c(0.5, 1, 0, 2, -0.5, pv, if (length(pv) >= 2L) mean(pv[1:2]))
          if (discrete) cand <- unique(as.integer(round(c(0, 1, pv))))
          if (!any(vapply(cand, function(y) isTRUE(ok(y)), TRUE)))
            stop_domain(sub("^Domain error: ", "", gsub(dens_name, fname, conditionMessage(r), fixed = TRUE)))
        }
        fun(...)
      }, unchecked = fun)
    }
    # random number generators of univariate families check their parameters
    # the same way, with a scalar variable in the support
    if (!is.null(rng) && args[[1L]] %in% c("reals", "ints")) {
      rng_fun <- rng
      rng <- function(...) {
        a <- lapply(list(...), val)
        pv <- unlist(lapply(a, function(x) as.vector(x)[1L]))
        cand <- c(0.5, 1, 0, 2, -0.5, pv, if (length(pv) >= 2L) mean(pv[1:2]))
        if (discrete) cand <- unique(as.integer(round(c(0, 1, pv))))
        err <- NULL
        for (y in cand) {
          r <- tryCatch({ do.call(lpdf, c(list(y), a, list(propto = TRUE))); TRUE },
                        bayescomp_domain_error = function(e) e, bayescomp_shape_error = function(e) TRUE)
          if (isTRUE(r)) return(rng_fun(...))
          err <- r
        }
        stop_domain(sub("^Domain error: ", "", gsub(dens_name, paste0(name, "_rng"), conditionMessage(err), fixed = TRUE)))
      }
    }
    if (!is.null(lcdf)) lcdf <- checked(lcdf, paste0(name, "_lcdf"))
    if (!is.null(cdf)) cdf <- checked(cdf, paste0(name, "_cdf"))
    if (!is.null(lccdf)) lccdf <- checked(lccdf, paste0(name, "_lccdf"))
  }
  if (!is.null(lpdf) && densities) {
    f1 <- lpdf
    def(paste0(name, suf[1L]), args, ret, function(...) f1(..., propto = FALSE), density = TRUE, dist = name)
    def(paste0(name, suf[2L]), args, ret, function(...) f1(..., propto = TRUE), density = TRUE, dist = name)
  }
  if (!is.null(lcdf)) {
    f2 <- lcdf
    def(paste0(name, "_lcdf"), args, "real", f2, dist = name)
    def(paste0(name, "_cdf"), args, "real", function(...) exp(f2(...)), dist = name)
  } else if (!is.null(cdf)) {
    f3 <- cdf
    def(paste0(name, "_cdf"), args, "real", f3, dist = name)
    def(paste0(name, "_lcdf"), args, "real", function(...) log(f3(...)), dist = name)
  }
  if (!is.null(lccdf)) {
    f4 <- lccdf
    def(paste0(name, "_lccdf"), args, "real", f4, dist = name)
  }
  if (!is.null(rng)) {
    f5 <- rng
    rr <- rng_ret %||% (if (discrete) "int" else "real")
    decl <- args[-1L]
    def(paste0(name, "_rng"), decl, function(types) {
      # only vectorized argument families batch draws: a container passed for
      # "reals" or "ints", or an array passed for "vectors"; a fixed vector
      # argument (categorical, ordered_logistic) does not
      batch <- vapply(seq_along(types), function(i) {
        d <- decl[[i]]; t <- types[[i]]
        (d %in% c("reals", "ints") && !is_scalar_ty(t)) || (d %in% c("vectors", "row_vectors") && t$nd > 0L)
      }, TRUE)
      if (!any(batch)) parse_ty_string(rr) else ty(rr, 1L)
    }, function(...) f5(...), rng = TRUE, dist = name)
  }
  invisible(NULL)
}

# vectorized rng helper: draws N values, N the common length of the arguments
vrng <- function(f, ...) {
  args <- lapply(list(...), function(a) as.vector(val(a)))
  N <- do.call(vlen_check, c(list("rng"), args))
  all_scalar <- all(lengths(args) == 1L)
  r <- do.call(f, c(list(N), lapply(args, rep_len, length.out = N)))
  r
}

# ---- continuous ------------------------------------------------------------------

lp_normal <- function(y, mu, sigma, propto) {
  f <- "normal_lpdf"
  chk_not_nan(f, "Random variable", y); chk_finite(f, "Location parameter", mu); chk_positive_finite(f, "Scale parameter", sigma)
  N <- vlen_check(f, y, mu, sigma)
  if (N == 0L) return(0)
  vy <- as.vector(val(y)); vm <- as.vector(val(mu)); vs <- as.vector(val(sigma))
  inv_s <- 1 / vs
  z <- (vy - vm) * inv_s
  lp <- 0
  if (!propto) lp <- lp - N * 0.9189385332046727
  if (inc(propto, sigma)) lp <- lp - (if (length(vs) == 1L) N * log(vs) else sum(log(vs)))
  if (inc(propto, y, mu, sigma)) lp <- lp - 0.5 * sum(z * z) * (if (length(z) == 1L) N else 1)
  if (!any_ad(y, mu, sigma)) return(lp)
  dy <- -z * inv_s
  fused(lp, list(y, mu, sigma), function() list(dy, -dy, z * z * inv_s - inv_s), N)
}
defdist("normal", c("reals", "reals", "reals"), lp_normal,
  lcdf = function(y, mu, s) ad_sum_v(ad_log_Phi((as_dbl(y) - as_dbl(mu)) / as_dbl(s)), y, mu, s),
  lccdf = function(y, mu, s) ad_sum_v(ad_log1m_Phi((as_dbl(y) - as_dbl(mu)) / as_dbl(s)), y, mu, s),
  rng = function(mu, s) {
    chk_finite("normal_rng", "Location parameter", mu); chk_positive_finite("normal_rng", "Scale parameter", s)
    vrng(function(n, m, s) stats::rnorm(n, m, s), mu, s)
  })

# sum of an elementwise term; scalar terms are counted N times
# y raised to the lower end of the support, so cumulative functions of values
# below the support give cdf 0 and ccdf 1 (the gradient there is zero)
clamp_lo <- function(y, lo = 0) ad_ifelse(val(y) >= lo | is.na(val(y)), y, lo)

ad_sum_v <- function(term, ...) {
  N <- vlen_check("cdf", ...)
  bsum(flat(term), N)
}

lp_std_normal <- function(y, propto) {
  chk_not_nan("std_normal_lpdf", "Random variable", y)
  vy <- as.vector(val(y))
  N <- length(vy)
  lp <- 0
  if (!propto) lp <- lp - N * 0.9189385332046727
  if (inc(propto, y)) lp <- lp - 0.5 * sum(vy * vy)
  if (!is_ad(y)) return(lp)
  fused(lp, list(y), function() list(-vy), N)
}
defdist("std_normal", "reals", lp_std_normal,
  lcdf = function(y) ad_sum(flat(ad_log_Phi(as_dbl(y)))),
  lccdf = function(y) ad_sum(flat(ad_log1m_Phi(as_dbl(y)))))
def("std_normal_rng", character(0), "real", function() stats::rnorm(1), rng = TRUE)

lp_normal_id_glm <- function(y, x, alpha, beta, sigma, propto) {
  mu <- glm_linpred("normal_id_glm_lpdf", x, alpha, beta)
  lp_normal(y, mu, sigma, propto)
}
glm_linpred <- function(f, x, alpha, beta) {
  vx <- val(x)
  if (is.null(dim(vx))) {
    # x is a row_vector or vector: a single observation
    eta <- op_dot(x, beta) + as_dbl(alpha)
    return(eta)
  }
  if (ncol(vx) != nelem(beta)) stop_shape(sprintf("%s: columns of x (%d) must match the size of beta (%d)", f, ncol(vx), nelem(beta)))
  eta <- op_mat_vec(x, beta)
  eta + flat(as_dbl(alpha))
}
for (xt in c("matrix", "row_vector", "vector")) {
  def("normal_id_glm_lpdf", c("reals", xt, "reals", "vector", "reals"), "real", function(...) lp_normal_id_glm(..., propto = FALSE), density = TRUE, dist = "normal_id_glm")
  def("normal_id_glm_lupdf", c("reals", xt, "reals", "vector", "reals"), "real", function(...) lp_normal_id_glm(..., propto = TRUE), density = TRUE, dist = "normal_id_glm")
}

lp_student_t <- function(y, nu, mu, sigma, propto) {
  f <- "student_t_lpdf"
  chk_not_nan(f, "Random variable", y); chk_positive_finite(f, "Degrees of freedom parameter", nu)
  chk_finite(f, "Location parameter", mu); chk_positive_finite(f, "Scale parameter", sigma)
  N <- vlen_check(f, y, nu, mu, sigma)
  if (N == 0L) return(0)
  vy <- as.vector(val(y)); vn <- as.vector(val(nu)); vm <- as.vector(val(mu)); vs <- as.vector(val(sigma))
  rep1 <- function(x) if (length(x) == 1L) N * x else sum(x)
  z <- (vy - vm) / vs
  r <- z * z / vn
  lp <- 0
  if (inc(propto, nu)) lp <- lp + rep1(lgamma((vn + 1) / 2) - lgamma(vn / 2) - 0.5 * log(vn))
  if (!propto) lp <- lp - N * 0.5 * log(pi)
  if (inc(propto, sigma)) lp <- lp - rep1(log(vs))
  if (inc(propto, y, nu, mu, sigma)) lp <- lp - rep1((vn + 1) / 2 * log1p(r))
  if (!any_ad(y, nu, mu, sigma)) return(lp)
  fused(lp, list(y, nu, mu, sigma), function() {
    dy <- -(vn + 1) * z / (vs * vn * (1 + r))
    list(dy,
         0.5 * digamma((vn + 1) / 2) - 0.5 * digamma(vn / 2) - 0.5 / vn - 0.5 * log1p(r) + (vn + 1) / 2 * r / (vn * (1 + r)),
         -dy,
         -1 / vs + (vn + 1) * r / (vs * (1 + r)))
  }, N)
}
pt_ls <- function(y, nu, mu, s, lower, logp) {
  stats::pt((y - mu) / s, nu, lower.tail = lower, log.p = logp)
}
st_cdf_ad <- function(y, nu, mu, s, lower) {
  y <- flat(as_dbl(y)); nu <- flat(as_dbl(nu)); mu <- flat(as_dbl(mu)); s <- flat(as_dbl(s))
  N <- vlen_check("student_t_lcdf", y, nu, mu, s)
  vy <- rep_len(val(y), N); vn <- rep_len(val(nu), N); vm <- rep_len(val(mu), N); vs <- rep_len(val(s), N)
  lp <- pt_ls(vy, vn, vm, vs, lower, TRUE)
  out <- sum(lp)
  if (!any_ad(y, nu, mu, s)) return(out)
  dens <- exp(stats::dt((vy - vm) / vs, vn, log = TRUE) - lp) / vs
  sg <- if (lower) 1 else -1
  ad_node(out, list(y, nu, mu, s), function(g) list(
    if (is_ad(y)) unbroadcast(g * sg * dens, val(y)),
    if (is_ad(nu)) unbroadcast(g * num_partial(function(t) pt_ls(vy, t, vm, vs, lower, TRUE), vn), val(nu)),
    if (is_ad(mu)) unbroadcast(-g * sg * dens, val(mu)),
    if (is_ad(s)) unbroadcast(-g * sg * dens * (vy - vm) / vs, val(s))))
}
defdist("student_t", c("reals", "reals", "reals", "reals"), lp_student_t,
  lcdf = function(y, nu, mu, s) st_cdf_ad(y, nu, mu, s, TRUE),
  lccdf = function(y, nu, mu, s) st_cdf_ad(y, nu, mu, s, FALSE),
  rng = function(nu, mu, s) vrng(function(n, nu, m, s) m + s * stats::rt(n, nu), nu, mu, s))
# vectorized: any container argument gives a vector
def("student_t_qf", c("reals", "reals", "reals", "reals"), function(types) if (all(vapply(types, is_scalar_ty, TRUE))) T_real else ty("vector", 0L),
    function(p, nu, mu, s) {
  if (any_ad(nu)) stop_unsupported("student_t_qf: gradient with respect to the degrees of freedom is not available")
  N <- vlen_check("student_t_qf", p, nu, mu, s)
  scalar <- all(vapply(list(p, nu, mu, s), function(a) length(val(a)) == 1L, TRUE))
  p <- flat(as_dbl(p)); mu <- flat(as_dbl(mu)); s <- flat(as_dbl(s))
  if (length(val(p)) == 1L && N > 1L) p <- p + numeric(N)
  nv <- rep_len(as.vector(val(nu)), N)
  q <- ad_unary(p, function(pp) stats::qt(pp, nv), function(pp, y) 1 / stats::dt(y, nv))
  out <- mu + s * q
  if (scalar) out else flat(out)
})

lp_cauchy <- function(y, mu, sigma, propto) {
  f <- "cauchy_lpdf"
  chk_not_nan(f, "Random variable", y); chk_finite(f, "Location parameter", mu); chk_positive_finite(f, "Scale parameter", sigma)
  N <- vlen_check(f, y, mu, sigma)
  if (N == 0L) return(0)
  vy <- as.vector(val(y)); vm <- as.vector(val(mu)); vs <- as.vector(val(sigma))
  rep1 <- function(x) if (length(x) == 1L) N * x else sum(x)
  z <- (vy - vm) / vs
  lp <- 0
  if (!propto) lp <- lp - N * log(pi)
  if (inc(propto, sigma)) lp <- lp - rep1(log(vs))
  if (inc(propto, y, mu, sigma)) lp <- lp - rep1(log1p(z * z))
  if (!any_ad(y, mu, sigma)) return(lp)
  fused(lp, list(y, mu, sigma), function() {
    dy <- -2 * z / (vs * (1 + z * z))
    list(dy, -dy, -1 / vs + 2 * z * z / (vs * (1 + z * z)))
  }, N)
}
# log of the Cauchy cdf (lower) or ccdf at z; in the far tail 0.5 - atan(z) / pi
# cancels, so it is computed as atan(1 / z) / pi there
cauchy_tail <- function(z, lower) {
  zz <- if (lower) -z else z
  far <- val(zz) > 1
  if (!any(far)) return(log(0.5 - atan(zz) / pi))
  safe <- ad_ifelse(far, zz, 2)
  ad_ifelse(far, log(atan(1 / safe) / pi), log(0.5 - atan(ad_ifelse(far, 0, zz)) / pi))
}
defdist("cauchy", c("reals", "reals", "reals"), lp_cauchy,
  lcdf = function(y, mu, s) ad_sum_v(cauchy_tail((as_dbl(y) - as_dbl(mu)) / as_dbl(s), TRUE), y, mu, s),
  lccdf = function(y, mu, s) ad_sum_v(cauchy_tail((as_dbl(y) - as_dbl(mu)) / as_dbl(s), FALSE), y, mu, s),
  rng = function(mu, s) vrng(function(n, m, s) stats::rcauchy(n, m, s), mu, s))

lp_double_exponential <- function(y, mu, sigma, propto) {
  f <- "double_exponential_lpdf"
  chk_finite(f, "Random variable", y); chk_finite(f, "Location parameter", mu); chk_positive_finite(f, "Scale parameter", sigma)
  y <- flat(as_dbl(y)); mu <- flat(as_dbl(mu)); sigma <- flat(as_dbl(sigma))
  N <- vlen_check(f, y, mu, sigma)
  lp <- 0
  if (!propto) lp <- lp - N * log(2)
  if (inc(propto, sigma)) lp <- lp - bsum(log(sigma), N)
  if (inc(propto, y, mu, sigma)) lp <- lp - bsum(abs(y - mu) / sigma, N)
  lp
}
de_lcdf <- function(y, mu, s) {
  z <- (as_dbl(y) - as_dbl(mu)) / as_dbl(s)
  ad_sum_v(ad_ifelse(val(z) < 0, log(0.5) + z, ad_log1m(0.5 * exp(-z))), y, mu, s)
}
de_lccdf <- function(y, mu, s) {
  z <- (as_dbl(y) - as_dbl(mu)) / as_dbl(s)
  ad_sum_v(ad_ifelse(val(z) < 0, ad_log1m(0.5 * exp(z)), log(0.5) - z), y, mu, s)
}
defdist("double_exponential", c("reals", "reals", "reals"), lp_double_exponential,
  lcdf = de_lcdf, lccdf = de_lccdf,
  rng = function(mu, s) vrng(function(n, m, s) {
    u <- stats::runif(n, -0.5, 0.5)
    m - s * sign(u) * log1p(-2 * abs(u))
  }, mu, s))

lp_logistic <- function(y, mu, sigma, propto) {
  f <- "logistic_lpdf"
  chk_finite(f, "Random variable", y); chk_finite(f, "Location parameter", mu); chk_positive_finite(f, "Scale parameter", sigma)
  y <- flat(as_dbl(y)); mu <- flat(as_dbl(mu)); sigma <- flat(as_dbl(sigma))
  N <- vlen_check(f, y, mu, sigma)
  lp <- 0
  if (inc(propto, sigma)) lp <- lp - bsum(log(sigma), N)
  if (inc(propto, y, mu, sigma)) {
    z <- (y - mu) / sigma
    lp <- lp - bsum(z + 2 * ad_log1p_exp(-z), N)
  }
  lp
}
defdist("logistic", c("reals", "reals", "reals"), lp_logistic,
  lcdf = function(y, mu, s) ad_sum_v(ad_log_inv_logit((as_dbl(y) - as_dbl(mu)) / as_dbl(s)), y, mu, s),
  lccdf = function(y, mu, s) ad_sum_v(ad_log1m_inv_logit((as_dbl(y) - as_dbl(mu)) / as_dbl(s)), y, mu, s),
  rng = function(mu, s) vrng(function(n, m, s) stats::rlogis(n, m, s), mu, s))

lp_gumbel <- function(y, mu, beta, propto) {
  f <- "gumbel_lpdf"
  chk_not_nan(f, "Random variable", y); chk_finite(f, "Location parameter", mu); chk_positive(f, "Scale parameter", beta)
  y <- flat(as_dbl(y)); mu <- flat(as_dbl(mu)); beta <- flat(as_dbl(beta))
  N <- vlen_check(f, y, mu, beta)
  lp <- 0
  if (inc(propto, beta)) lp <- lp - bsum(log(beta), N)
  if (inc(propto, y, mu, beta)) {
    z <- (y - mu) / beta
    lp <- lp - bsum(z + exp(-z), N)
  }
  lp
}
defdist("gumbel", c("reals", "reals", "reals"), lp_gumbel,
  lcdf = function(y, mu, b) ad_sum_v(-exp(-(as_dbl(y) - as_dbl(mu)) / as_dbl(b)), y, mu, b),
  lccdf = function(y, mu, b) ad_sum_v(ad_log1m_exp(-exp(-(as_dbl(y) - as_dbl(mu)) / as_dbl(b))), y, mu, b),
  rng = function(mu, b) vrng(function(n, m, b) m - b * log(-log(stats::runif(n))), mu, b))

lp_skew_normal <- function(y, xi, omega, alpha, propto) {
  f <- "skew_normal_lpdf"
  chk_not_nan(f, "Random variable", y); chk_finite(f, "Location parameter", xi)
  chk_positive(f, "Scale parameter", omega); chk_finite(f, "Shape parameter", alpha)
  y <- flat(as_dbl(y)); xi <- flat(as_dbl(xi)); omega <- flat(as_dbl(omega)); alpha <- flat(as_dbl(alpha))
  N <- vlen_check(f, y, xi, omega, alpha)
  lp <- 0
  if (!propto) lp <- lp + N * (log(2) - 0.5 * log(2 * pi))
  if (inc(propto, omega)) lp <- lp - bsum(log(omega), N)
  if (inc(propto, y, xi, omega, alpha)) {
    z <- (y - xi) / omega
    lp <- lp + bsum(-0.5 * z * z + ad_log_Phi(alpha * z), N)
  }
  lp
}
defdist("skew_normal", c("reals", "reals", "reals", "reals"), lp_skew_normal,
  rng = function(xi, om, al) vrng(function(n, xi, om, al) {
    d <- al / sqrt(1 + al^2)
    u0 <- stats::rnorm(n)
    v <- stats::rnorm(n)
    u1 <- d * u0 + sqrt(1 - d^2) * v
    xi + om * ifelse(u0 >= 0, u1, -u1)
  }, xi, om, al))

lp_lognormal <- function(y, mu, sigma, propto) {
  f <- "lognormal_lpdf"
  chk_nonnegative(f, "Random variable", y); chk_finite(f, "Location parameter", mu); chk_positive_finite(f, "Scale parameter", sigma)
  N <- vlen_check(f, y, mu, sigma)
  if (N == 0L) return(0)
  vy <- as.vector(val(y)); vm <- as.vector(val(mu)); vs <- as.vector(val(sigma))
  if (any(vy == 0)) return(-Inf)
  rep1 <- function(x) if (length(x) == 1L) N * x else sum(x)
  ly <- log(vy)
  z <- (ly - vm) / vs
  lp <- 0
  if (!propto) lp <- lp - N * 0.9189385332046727
  if (inc(propto, sigma)) lp <- lp - rep1(log(vs))
  if (inc(propto, y)) lp <- lp - rep1(ly)
  if (inc(propto, y, mu, sigma)) lp <- lp - 0.5 * rep1(z * z)
  if (!any_ad(y, mu, sigma)) return(lp)
  fused(lp, list(y, mu, sigma), function() list(-1 / vy - z / (vs * vy), z / vs, -1 / vs + z * z / vs), N)
}
ln_z <- function(y, mu, s) (log(as_dbl(y)) - as_dbl(mu)) / as_dbl(s)
defdist("lognormal", c("reals", "reals", "reals"), lp_lognormal,
  lcdf = function(y, mu, s) ad_sum_v(ad_log_Phi(ln_z(y, mu, s)), y, mu, s),
  lccdf = function(y, mu, s) ad_sum_v(ad_log1m_Phi(ln_z(y, mu, s)), y, mu, s),
  rng = function(mu, s) vrng(function(n, m, s) stats::rlnorm(n, m, s), mu, s))

lp_chi_square <- function(y, nu, propto) {
  f <- "chi_square_lpdf"
  chk_nonnegative(f, "Random variable", y); chk_positive_finite(f, "Degrees of freedom parameter", nu)
  y <- flat(as_dbl(y)); nu <- flat(as_dbl(nu))
  N <- vlen_check(f, y, nu)
  lp <- 0
  if (inc(propto, nu)) lp <- lp - bsum(lgamma(nu / 2) + nu / 2 * log(2), N)
  if (inc(propto, y, nu)) lp <- lp + bsum(ad_lmultiply(nu / 2 - 1, y), N)
  if (inc(propto, y)) lp <- lp - bsum(y / 2, N)
  lp
}
# lower incomplete gamma based cdfs (numerical partial for the shape)
gamma_cdf_ad <- function(y, shape, rate, lower) {
  y <- flat(as_dbl(y)); shape <- flat(as_dbl(shape)); rate <- flat(as_dbl(rate))
  N <- vlen_check("gamma_lcdf", y, shape, rate)
  vy <- rep_len(val(y), N); va <- rep_len(val(shape), N); vb <- rep_len(val(rate), N)
  lp <- stats::pgamma(vy, va, vb, lower.tail = lower, log.p = TRUE)
  out <- sum(lp)
  if (!any_ad(y, shape, rate)) return(out)
  dens <- exp(stats::dgamma(vy, va, vb, log = TRUE) - lp)
  sg <- if (lower) 1 else -1
  ad_node(out, list(y, shape, rate), function(g) list(
    if (is_ad(y)) unbroadcast(g * sg * dens * 1, val(y)),
    if (is_ad(shape)) unbroadcast(g * num_partial(function(t) stats::pgamma(vy, t, vb, lower.tail = lower, log.p = TRUE), va), val(shape)),
    if (is_ad(rate)) unbroadcast(g * sg * dens * vy / vb, val(rate))))
}
defdist("chi_square", c("reals", "reals"), lp_chi_square,
  lcdf = function(y, nu) gamma_cdf_ad(y, as_dbl(nu) / 2, 0.5, TRUE),
  lccdf = function(y, nu) gamma_cdf_ad(y, as_dbl(nu) / 2, 0.5, FALSE),
  rng = function(nu) vrng(function(n, nu) stats::rchisq(n, nu), nu))

lp_inv_chi_square <- function(y, nu, propto) {
  f <- "inv_chi_square_lpdf"
  chk_not_nan(f, "Random variable", y); chk_positive_finite(f, "Degrees of freedom parameter", nu)
  y <- flat(as_dbl(y)); nu <- flat(as_dbl(nu))
  N <- vlen_check(f, y, nu)
  if (any(val(y) <= 0)) return(-Inf)
  lp <- 0
  if (inc(propto, nu)) lp <- lp - bsum(nu / 2 * log(2) + lgamma(nu / 2), N)
  if (inc(propto, y, nu)) lp <- lp - bsum((nu / 2 + 1) * log(y), N)
  if (inc(propto, y)) lp <- lp - bsum(0.5 / y, N)
  lp
}
defdist("inv_chi_square", c("reals", "reals"), lp_inv_chi_square,
  lcdf = function(y, nu) gamma_cdf_ad(1 / as_dbl(y), as_dbl(nu) / 2, 0.5, FALSE),
  lccdf = function(y, nu) gamma_cdf_ad(1 / as_dbl(y), as_dbl(nu) / 2, 0.5, TRUE),
  rng = function(nu) vrng(function(n, nu) 1 / stats::rchisq(n, nu), nu))

lp_scaled_inv_chi_square <- function(y, nu, s, propto) {
  f <- "scaled_inv_chi_square_lpdf"
  chk_not_nan(f, "Random variable", y); chk_positive_finite(f, "Degrees of freedom parameter", nu); chk_positive_finite(f, "Scale parameter", s)
  y <- flat(as_dbl(y)); nu <- flat(as_dbl(nu)); s <- flat(as_dbl(s))
  N <- vlen_check(f, y, nu, s)
  if (any(val(y) <= 0)) return(-Inf)
  lp <- 0
  if (inc(propto, nu)) lp <- lp + bsum(nu / 2 * log(nu / 2) - lgamma(nu / 2), N)
  if (inc(propto, nu, s)) lp <- lp + bsum(nu * log(s), N)
  if (inc(propto, y, nu)) lp <- lp - bsum((nu / 2 + 1) * log(y), N)
  if (inc(propto, y, nu, s)) lp <- lp - bsum(nu / 2 * s * s / y, N)
  lp
}
defdist("scaled_inv_chi_square", c("reals", "reals", "reals"), lp_scaled_inv_chi_square,
  lcdf = function(y, nu, s) gamma_cdf_ad(as_dbl(nu) * as_dbl(s) * as_dbl(s) / (2 * as_dbl(y)), as_dbl(nu) / 2, 1, FALSE),
  lccdf = function(y, nu, s) gamma_cdf_ad(as_dbl(nu) * as_dbl(s) * as_dbl(s) / (2 * as_dbl(y)), as_dbl(nu) / 2, 1, TRUE),
  rng = function(nu, s) vrng(function(n, nu, s) nu * s^2 / stats::rchisq(n, nu), nu, s))

lp_exponential <- function(y, beta, propto) {
  f <- "exponential_lpdf"
  chk_nonnegative(f, "Random variable", y); chk_positive_finite(f, "Inverse scale parameter", beta)
  N <- vlen_check(f, y, beta)
  if (N == 0L) return(0)
  vy <- as.vector(val(y)); vb <- as.vector(val(beta))
  rep1 <- function(x) if (length(x) == 1L) N * x else sum(x)
  lp <- 0
  if (inc(propto, beta)) lp <- lp + rep1(log(vb))
  if (inc(propto, y, beta)) lp <- lp - rep1(vb * vy)
  if (!any_ad(y, beta)) return(lp)
  fused(lp, list(y, beta), function() list(-vb + 0 * vy, 1 / vb - vy), N)
}
defdist("exponential", c("reals", "reals"), lp_exponential,
  lcdf = function(y, b) ad_sum_v(ad_log1m_exp(-as_dbl(b) * clamp_lo(as_dbl(y))), y, b),
  lccdf = function(y, b) ad_sum_v(-as_dbl(b) * clamp_lo(as_dbl(y)), y, b),
  rng = function(b) vrng(function(n, b) stats::rexp(n, b), b))

lp_gamma <- function(y, alpha, beta, propto) {
  f <- "gamma_lpdf"
  chk_not_nan(f, "Random variable", y); chk_positive_finite(f, "Shape parameter", alpha); chk_positive_finite(f, "Inverse scale parameter", beta)
  N <- vlen_check(f, y, alpha, beta)
  if (N == 0L) return(0)
  vy <- as.vector(val(y)); va <- as.vector(val(alpha)); vb <- as.vector(val(beta))
  if (any(vy < 0)) return(-Inf)
  rep1 <- function(x) if (length(x) == 1L) N * x else sum(x)
  ly <- log(vy)
  lp <- 0
  if (inc(propto, alpha)) lp <- lp - rep1(lgamma(va))
  if (inc(propto, alpha, beta)) lp <- lp + rep1(va * log(vb))
  if (inc(propto, y, alpha)) lp <- lp + rep1(zmul(va - 1, ly))
  if (inc(propto, y, beta)) lp <- lp - rep1(vb * vy)
  if (!any_ad(y, alpha, beta)) return(lp)
  fused(lp, list(y, alpha, beta), function() list(zmul(va - 1, 1 / vy) - vb, -digamma(va) + log(vb) + ly, va / vb - vy), N)
}
defdist("gamma", c("reals", "reals", "reals"), lp_gamma,
  lcdf = function(y, a, b) gamma_cdf_ad(y, a, b, TRUE),
  lccdf = function(y, a, b) gamma_cdf_ad(y, a, b, FALSE),
  rng = function(a, b) vrng(function(n, a, b) stats::rgamma(n, a, b), a, b))

lp_inv_gamma <- function(y, alpha, beta, propto) {
  f <- "inv_gamma_lpdf"
  chk_not_nan(f, "Random variable", y); chk_positive_finite(f, "Shape parameter", alpha); chk_positive_finite(f, "Scale parameter", beta)
  y <- flat(as_dbl(y)); alpha <- flat(as_dbl(alpha)); beta <- flat(as_dbl(beta))
  N <- vlen_check(f, y, alpha, beta)
  if (any(val(y) <= 0)) return(-Inf)
  lp <- 0
  if (inc(propto, alpha)) lp <- lp - bsum(lgamma(alpha), N)
  if (inc(propto, alpha, beta)) lp <- lp + bsum(alpha * log(beta), N)
  if (inc(propto, y, alpha)) lp <- lp - bsum((alpha + 1) * log(y), N)
  if (inc(propto, y, beta)) lp <- lp - bsum(beta / y, N)
  lp
}
defdist("inv_gamma", c("reals", "reals", "reals"), lp_inv_gamma,
  lcdf = function(y, a, b) gamma_cdf_ad(1 / clamp_lo(as_dbl(y)), a, b, FALSE),
  lccdf = function(y, a, b) gamma_cdf_ad(1 / clamp_lo(as_dbl(y)), a, b, TRUE),
  rng = function(a, b) vrng(function(n, a, b) 1 / stats::rgamma(n, a, b), a, b))

lp_weibull <- function(y, alpha, sigma, propto) {
  f <- "weibull_lpdf"
  chk_finite(f, "Random variable", y); chk_positive_finite(f, "Shape parameter", alpha); chk_positive_finite(f, "Scale parameter", sigma)
  y <- flat(as_dbl(y)); alpha <- flat(as_dbl(alpha)); sigma <- flat(as_dbl(sigma))
  N <- vlen_check(f, y, alpha, sigma)
  if (any(val(y) < 0)) return(-Inf)
  lp <- 0
  if (inc(propto, alpha)) lp <- lp + bsum(log(alpha), N)
  if (inc(propto, y, alpha)) lp <- lp + bsum(ad_lmultiply(alpha - 1, y), N)
  if (inc(propto, alpha, sigma)) lp <- lp - bsum(alpha * log(sigma), N)
  if (inc(propto, y, alpha, sigma)) lp <- lp - bsum((y / sigma)^alpha, N)
  lp
}
defdist("weibull", c("reals", "reals", "reals"), lp_weibull,
  lcdf = function(y, a, s) ad_sum_v(ad_log1m_exp(-(clamp_lo(as_dbl(y)) / as_dbl(s))^as_dbl(a)), y, a, s),
  lccdf = function(y, a, s) ad_sum_v(-(clamp_lo(as_dbl(y)) / as_dbl(s))^as_dbl(a), y, a, s),
  rng = function(a, s) vrng(function(n, a, s) stats::rweibull(n, a, s), a, s))

lp_frechet <- function(y, alpha, sigma, propto) {
  f <- "frechet_lpdf"
  chk_positive(f, "Random variable", y); chk_positive_finite(f, "Shape parameter", alpha); chk_positive_finite(f, "Scale parameter", sigma)
  y <- flat(as_dbl(y)); alpha <- flat(as_dbl(alpha)); sigma <- flat(as_dbl(sigma))
  N <- vlen_check(f, y, alpha, sigma)
  lp <- 0
  if (inc(propto, alpha)) lp <- lp + bsum(log(alpha), N)
  if (inc(propto, y, alpha)) lp <- lp - bsum((alpha + 1) * log(y), N)
  if (inc(propto, alpha, sigma)) lp <- lp + bsum(alpha * log(sigma), N)
  if (inc(propto, y, alpha, sigma)) lp <- lp - bsum((y / sigma)^(-alpha), N)
  lp
}
defdist("frechet", c("reals", "reals", "reals"), lp_frechet,
  lcdf = function(y, a, s) ad_sum_v(-(as_dbl(y) / as_dbl(s))^(-as_dbl(a)), y, a, s),
  lccdf = function(y, a, s) ad_sum_v(ad_log1m_exp(-(as_dbl(y) / as_dbl(s))^(-as_dbl(a))), y, a, s),
  rng = function(a, s) vrng(function(n, a, s) s * (-log(stats::runif(n)))^(-1 / a), a, s))

lp_rayleigh <- function(y, sigma, propto) {
  f <- "rayleigh_lpdf"
  chk_nonnegative(f, "Random variable", y); chk_positive(f, "Scale parameter", sigma)
  y <- flat(as_dbl(y)); sigma <- flat(as_dbl(sigma))
  N <- vlen_check(f, y, sigma)
  lp <- 0
  if (inc(propto, sigma)) lp <- lp - 2 * bsum(log(sigma), N)
  if (inc(propto, y)) lp <- lp + bsum(log(y), N)
  if (inc(propto, y, sigma)) lp <- lp - bsum(y * y / (2 * sigma * sigma), N)
  lp
}
defdist("rayleigh", c("reals", "reals"), lp_rayleigh,
  lcdf = function(y, s) ad_sum_v(ad_log1m_exp(-as_dbl(y) * as_dbl(y) / (2 * as_dbl(s) * as_dbl(s))), y, s),
  lccdf = function(y, s) ad_sum_v(-as_dbl(y) * as_dbl(y) / (2 * as_dbl(s) * as_dbl(s)), y, s),
  rng = function(s) vrng(function(n, s) s * sqrt(-2 * log(stats::runif(n))), s))

lp_pareto <- function(y, ymin, alpha, propto) {
  f <- "pareto_lpdf"
  chk_not_nan(f, "Random variable", y); chk_positive_finite(f, "Scale parameter", ymin); chk_positive_finite(f, "Shape parameter", alpha)
  y <- flat(as_dbl(y)); ymin <- flat(as_dbl(ymin)); alpha <- flat(as_dbl(alpha))
  N <- vlen_check(f, y, ymin, alpha)
  if (any(rep_len(val(y), N) < rep_len(val(ymin), N))) return(-Inf)
  lp <- 0
  if (inc(propto, alpha)) lp <- lp + bsum(log(alpha), N)
  if (inc(propto, ymin, alpha)) lp <- lp + bsum(alpha * log(ymin), N)
  if (inc(propto, y, alpha)) lp <- lp - bsum((alpha + 1) * log(y), N)
  lp
}
defdist("pareto", c("reals", "reals", "reals"), lp_pareto,
  lcdf = function(y, m, a) ad_sum_v(ad_log1m_exp(as_dbl(a) * (log(as_dbl(m)) - log(as_dbl(y)))), y, m, a),
  lccdf = function(y, m, a) ad_sum_v(as_dbl(a) * (log(as_dbl(m)) - log(as_dbl(y))), y, m, a),
  rng = function(m, a) vrng(function(n, m, a) m * stats::runif(n)^(-1 / a), m, a))

lp_pareto_type_2 <- function(y, mu, lambda, alpha, propto) {
  f <- "pareto_type_2_lpdf"
  chk_not_nan(f, "Random variable", y); chk_finite(f, "Location parameter", mu)
  chk_positive_finite(f, "Scale parameter", lambda); chk_positive_finite(f, "Shape parameter", alpha)
  y <- flat(as_dbl(y)); mu <- flat(as_dbl(mu)); lambda <- flat(as_dbl(lambda)); alpha <- flat(as_dbl(alpha))
  N <- vlen_check(f, y, mu, lambda, alpha)
  if (any(rep_len(val(y), N) < rep_len(val(mu), N))) stop_domain(sprintf("%s: Random variable must be greater than or equal to the location parameter", f))
  lp <- 0
  if (inc(propto, alpha)) lp <- lp + bsum(log(alpha), N)
  if (inc(propto, lambda)) lp <- lp - bsum(log(lambda), N)
  if (inc(propto, y, mu, lambda, alpha)) lp <- lp - bsum((alpha + 1) * log1p((y - mu) / lambda), N)
  lp
}
defdist("pareto_type_2", c("reals", "reals", "reals", "reals"), lp_pareto_type_2,
  lcdf = function(y, m, l, a) ad_sum_v(ad_log1m_exp(-as_dbl(a) * log1p((as_dbl(y) - as_dbl(m)) / as_dbl(l))), y, m, l, a),
  lccdf = function(y, m, l, a) ad_sum_v(-as_dbl(a) * log1p((as_dbl(y) - as_dbl(m)) / as_dbl(l)), y, m, l, a),
  rng = function(m, l, a) vrng(function(n, m, l, a) m + l * (stats::runif(n)^(-1 / a) - 1), m, l, a))

lp_beta <- function(theta, a, b, propto) {
  f <- "beta_lpdf"
  chk_positive_finite(f, "First shape parameter", a); chk_positive_finite(f, "Second shape parameter", b)
  chk_bounded(f, "Random variable", theta, 0, 1)
  N <- vlen_check(f, theta, a, b)
  if (N == 0L) return(0)
  vt <- as.vector(val(theta)); va <- as.vector(val(a)); vb <- as.vector(val(b))
  rep1 <- function(x) if (length(x) == 1L) N * x else sum(x)
  lt <- log(vt); l1t <- log1p(-vt)
  lp <- 0
  if (inc(propto, a, b)) lp <- lp - rep1(lbeta(va, vb))
  if (inc(propto, theta, a)) lp <- lp + rep1(zmul(va - 1, lt))
  if (inc(propto, theta, b)) lp <- lp + rep1(zmul(vb - 1, l1t))
  if (!any_ad(theta, a, b)) return(lp)
  fused(lp, list(theta, a, b), function() {
    dab <- digamma(va + vb)
    list(zmul(va - 1, 1 / vt) - zmul(vb - 1, 1 / (1 - vt)),
         -digamma(va) + dab + lt, -digamma(vb) + dab + l1t)
  }, N)
}
# x * log(y) with 0 * log(0) = 0
ad_mult_log <- function(x, y) {
  ad_binary(x, y, function(a, b) zmul(a, log(b)),
            function(a, b, r) { d <- log(b) + 0 * a; d[(a == 0 & b == 0) | is.infinite(d)] <- 0; d },
            function(a, b, r) zmul(a, 1 / b))
}
beta_cdf_ad <- function(y, a, b, lower) {
  y <- flat(as_dbl(y)); a <- flat(as_dbl(a)); b <- flat(as_dbl(b))
  N <- vlen_check("beta_lcdf", y, a, b)
  vy <- rep_len(val(y), N); va <- rep_len(val(a), N); vb <- rep_len(val(b), N)
  lp <- stats::pbeta(vy, va, vb, lower.tail = lower, log.p = TRUE)
  out <- sum(lp)
  if (!any_ad(y, a, b)) return(out)
  dens <- exp(stats::dbeta(vy, va, vb, log = TRUE) - lp)
  sg <- if (lower) 1 else -1
  ad_node(out, list(y, a, b), function(g) list(
    if (is_ad(y)) unbroadcast(g * sg * dens, val(y)),
    if (is_ad(a)) unbroadcast(g * num_partial(function(t) stats::pbeta(vy, t, vb, lower.tail = lower, log.p = TRUE), va), val(a)),
    if (is_ad(b)) unbroadcast(g * num_partial(function(t) stats::pbeta(vy, va, t, lower.tail = lower, log.p = TRUE), vb), val(b))))
}
defdist("beta", c("reals", "reals", "reals"), lp_beta,
  lcdf = function(y, a, b) beta_cdf_ad(y, a, b, TRUE),
  lccdf = function(y, a, b) beta_cdf_ad(y, a, b, FALSE),
  rng = function(a, b) vrng(function(n, a, b) stats::rbeta(n, a, b), a, b))

lp_beta_proportion <- function(theta, mu, kappa, propto) {
  f <- "beta_proportion_lpdf"
  chk_bounded(f, "Location parameter", mu, 0, 1); chk_positive_finite(f, "Precision parameter", kappa)
  mu <- as_dbl(mu); kappa <- as_dbl(kappa)
  lp_beta(theta, mu * kappa, (1 - mu) * kappa, propto)
}
defdist("beta_proportion", c("reals", "reals", "reals"), lp_beta_proportion,
  lcdf = function(y, mu, k) beta_cdf_ad(y, as_dbl(mu) * as_dbl(k), (1 - as_dbl(mu)) * as_dbl(k), TRUE),
  lccdf = function(y, mu, k) beta_cdf_ad(y, as_dbl(mu) * as_dbl(k), (1 - as_dbl(mu)) * as_dbl(k), FALSE),
  rng = function(mu, k) vrng(function(n, mu, k) stats::rbeta(n, mu * k, (1 - mu) * k), mu, k))

lp_von_mises <- function(y, mu, kappa, propto) {
  f <- "von_mises_lpdf"
  chk_finite(f, "Random variable", y); chk_finite(f, "Location parameter", mu); chk_nonnegative(f, "Scale parameter", kappa)
  y <- flat(as_dbl(y)); mu <- flat(as_dbl(mu)); kappa <- flat(as_dbl(kappa))
  N <- vlen_check(f, y, mu, kappa)
  lp <- 0
  if (!propto) lp <- lp - N * log(2 * pi)
  if (inc(propto, kappa)) lp <- lp - bsum(ad_unary(kappa, function(k) log(besselI(k, 0, expon.scaled = TRUE)) + k,
                                                  function(k, y) besselI(k, 1, expon.scaled = TRUE) / besselI(k, 0, expon.scaled = TRUE)), N)
  if (inc(propto, y, mu, kappa)) lp <- lp + bsum(kappa * cos(y - mu), N)
  lp
}
defdist("von_mises", c("reals", "reals", "reals"), lp_von_mises,
  rng = function(mu, k) vrng(function(n, mu, k) {
    # Best and Fisher (1979) rejection sampler
    vapply(seq_len(n), function(i) {
      kk <- k[i]
      if (kk < 1e-8) return(stats::runif(1, -pi, pi) + mu[i])
      tau <- 1 + sqrt(1 + 4 * kk^2)
      rho <- (tau - sqrt(2 * tau)) / (2 * kk)
      r <- (1 + rho^2) / (2 * rho)
      repeat {
        u <- stats::runif(3)
        z <- cos(pi * u[1])
        fz <- (1 + r * z) / (r + z)
        cc <- kk * (r - fz)
        if (cc * (2 - cc) - u[2] > 0 || log(cc / u[2]) + 1 - cc >= 0) break
      }
      th <- sign(u[3] - 0.5) * acos(fz)
      ((th + mu[i] + pi) %% (2 * pi)) - pi
    }, 0)
  }, mu, k))

lp_uniform <- function(y, a, b, propto) {
  f <- "uniform_lpdf"
  chk_not_nan(f, "Random variable", y); chk_finite(f, "Lower bound parameter", a); chk_finite(f, "Upper bound parameter", b)
  y <- flat(as_dbl(y)); a <- flat(as_dbl(a)); b <- flat(as_dbl(b))
  N <- vlen_check(f, y, a, b)
  vy <- rep_len(val(y), N)
  if (any(rep_len(val(b), N) <= rep_len(val(a), N))) stop_domain(sprintf("%s: Upper bound parameter must be greater than the lower bound", f))
  if (any(vy < rep_len(val(a), N) | vy > rep_len(val(b), N))) return(-Inf)
  if (inc(propto, a, b)) return(-bsum(log(b - a), N))
  0
}
defdist("uniform", c("reals", "reals", "reals"), lp_uniform,
  lcdf = function(y, a, b) {
    vy <- val(y)
    if (any(vy < val(a))) return(-Inf)
    ad_sum_v(log(ad_ifelse(vy >= val(b), 1 + 0 * as_dbl(y), (as_dbl(y) - as_dbl(a)) / (as_dbl(b) - as_dbl(a)))), y, a, b)
  },
  lccdf = function(y, a, b) {
    vy <- val(y)
    if (any(vy > val(b))) return(-Inf)
    ad_sum_v(log(ad_ifelse(vy <= val(a), 1 + 0 * as_dbl(y), (as_dbl(b) - as_dbl(y)) / (as_dbl(b) - as_dbl(a)))), y, a, b)
  },
  rng = function(a, b) vrng(function(n, a, b) stats::runif(n, a, b), a, b))

lp_exp_mod_normal <- function(y, mu, sigma, lambda, propto) {
  f <- "exp_mod_normal_lpdf"
  chk_not_nan(f, "Random variable", y); chk_finite(f, "Location parameter", mu)
  chk_positive_finite(f, "Scale parameter", sigma); chk_positive_finite(f, "Inv_scale parameter", lambda)
  y <- flat(as_dbl(y)); mu <- flat(as_dbl(mu)); sigma <- flat(as_dbl(sigma)); lambda <- flat(as_dbl(lambda))
  N <- vlen_check(f, y, mu, sigma, lambda)
  lp <- 0
  if (inc(propto, lambda)) lp <- lp + bsum(log(lambda), N)
  if (inc(propto, y, mu, sigma, lambda)) {
    lp <- lp + bsum(lambda * (mu + 0.5 * lambda * sigma * sigma - y) +
                      ad_log_Phi((y - mu) / sigma - lambda * sigma), N)
  }
  lp
}
defdist("exp_mod_normal", c("reals", "reals", "reals", "reals"), lp_exp_mod_normal,
  rng = function(mu, s, l) vrng(function(n, m, s, l) stats::rnorm(n, m, s) + stats::rexp(n, l), mu, s, l))

lp_skew_double_exponential <- function(y, mu, sigma, tau, propto) {
  f <- "skew_double_exponential_lpdf"
  chk_finite(f, "Random variable", y); chk_finite(f, "Location parameter", mu)
  chk_positive_finite(f, "Scale parameter", sigma); chk_bounded(f, "Skewness parameter", tau, 0, 1)
  y <- flat(as_dbl(y)); mu <- flat(as_dbl(mu)); sigma <- flat(as_dbl(sigma)); tau <- flat(as_dbl(tau))
  N <- vlen_check(f, y, mu, sigma, tau)
  lp <- 0
  if (!propto) lp <- lp + N * log(2)
  if (inc(propto, tau)) lp <- lp + bsum(log(tau) + ad_log1m(tau), N)
  if (inc(propto, sigma)) lp <- lp - bsum(log(sigma), N)
  if (inc(propto, y, mu, sigma, tau)) {
    d <- y - mu
    below <- rep_len(val(y), N) < rep_len(val(mu), N)
    w <- ad_ifelse(below, (1 - tau) * (0 - d), tau * d)
    lp <- lp - bsum(2 / sigma * w, N)
  }
  lp
}
defdist("skew_double_exponential", c("reals", "reals", "reals", "reals"), lp_skew_double_exponential,
  rng = function(mu, s, t) vrng(function(n, m, s, t) {
    u <- stats::runif(n)
    ifelse(u < t, m + s / (2 * (1 - t)) * log(u / t), m - s / (2 * t) * log((1 - u) / (1 - t)))
  }, mu, s, t))

lp_loglogistic <- function(y, alpha, beta, propto) {
  f <- "loglogistic_lpdf"
  chk_nonnegative(f, "Random variable", y); chk_positive_finite(f, "Scale parameter", alpha); chk_positive_finite(f, "Shape parameter", beta)
  y <- flat(as_dbl(y)); alpha <- flat(as_dbl(alpha)); beta <- flat(as_dbl(beta))
  N <- vlen_check(f, y, alpha, beta)
  lp <- 0
  if (inc(propto, beta)) lp <- lp + bsum(log(beta), N)
  if (inc(propto, alpha, beta)) lp <- lp - bsum(beta * log(alpha), N)
  if (inc(propto, y, beta)) lp <- lp + bsum((beta - 1) * log(y), N)
  if (inc(propto, y, alpha, beta)) lp <- lp - 2 * bsum(log1p((y / alpha)^beta), N)
  lp
}
defdist("loglogistic", c("reals", "reals", "reals"), lp_loglogistic,
  cdf = function(y, a, b) exp(ad_sum_v(-log1p((as_dbl(y) / as_dbl(a))^(-as_dbl(b))), y, a, b)),
  rng = function(a, b) vrng(function(n, a, b) { u <- stats::runif(n); a * (u / (1 - u))^(1 / b) }, a, b))

# ---- discrete -----------------------------------------------------------------------

chk_int_range <- function(f, name, n, lo = 0L, hi = Inf) {
  if (any(n < lo | n > hi | is.na(n))) stop_domain(sprintf("%s: %s is %s, but must be in [%s, %s]", f, name, n[which(n < lo | n > hi | is.na(n))[1L]], lo, hi))
}

lp_bernoulli <- function(n, theta, propto) {
  f <- "bernoulli_lpmf"
  chk_int_range(f, "n", n, 0L, 1L); chk_prob(f, "Probability parameter", theta)
  theta <- flat(as_dbl(theta))
  N <- vlen_check(f, n, theta)
  if (!inc(propto, theta)) return(0)
  nn <- rep_len(as.numeric(n), N)
  if (length(val(theta)) == 1L) {
    s <- sum(nn)
    return(ad_mult_log(s, theta) + ad_mult_log(N - s, 1 - theta))
  }
  ad_sum(ad_mult_log(nn, theta) + ad_mult_log(1 - nn, 1 - theta))
}
defdist("bernoulli", c("ints", "reals"), lp_bernoulli, discrete = TRUE,
  cdf = function(n, th) {
    th <- flat(as_dbl(th))
    N <- vlen_check("bernoulli_cdf", n, th)
    nn <- rep_len(n, N)
    if (any(nn < 0)) return(0)
    terms <- ad_ifelse(nn >= 1, 1 + 0 * th, 1 - th)
    exp(bsum(log(terms), if (length(val(terms)) == 1L) N else 1L))
  },
  rng = function(th) vrng(function(n, p) as.integer(stats::rbinom(n, 1L, p)), th))

lp_bernoulli_logit <- function(n, alpha, propto) {
  f <- "bernoulli_logit_lpmf"
  chk_int_range(f, "n", n, 0L, 1L); chk_not_nan(f, "Logit transformed probability parameter", alpha)
  N <- vlen_check(f, n, alpha)
  if (N == 0L || !inc(propto, alpha)) return(0)
  va <- as.vector(val(alpha))
  sgn <- 2 * as.vector(n) - 1
  x <- sgn * va
  lp <- sum(stats::plogis(x, log.p = TRUE)) * (if (length(x) == 1L) N else 1)
  if (!is_ad(alpha)) return(lp)
  fused(lp, list(alpha), function() list(sgn * stats::plogis(-x)), N)
}
defdist("bernoulli_logit", c("ints", "reals"), lp_bernoulli_logit, discrete = TRUE,
  rng = function(a) vrng(function(n, a) as.integer(stats::rbinom(n, 1L, stats::plogis(a))), a))
for (xt in c("matrix", "row_vector", "vector")) {
  def("bernoulli_logit_glm_lpmf", c("ints", xt, "reals", "vector"), "real", function(y, x, a, b) lp_bernoulli_logit(y, glm_linpred("bernoulli_logit_glm_lpmf", x, a, b), FALSE), density = TRUE, dist = "bernoulli_logit_glm")
  def("bernoulli_logit_glm_lupmf", c("ints", xt, "reals", "vector"), "real", function(y, x, a, b) lp_bernoulli_logit(y, glm_linpred("bernoulli_logit_glm_lpmf", x, a, b), TRUE), density = TRUE, dist = "bernoulli_logit_glm")
}
for (xt in c("matrix", "row_vector")) def("bernoulli_logit_glm_rng", c(xt, "reals", "vector"), "int[]", function(x, a, b) {
  eta <- val(glm_linpred("bernoulli_logit_glm_rng", x, a, b))
  as.integer(stats::rbinom(length(eta), 1L, stats::plogis(eta)))
}, rng = TRUE)

lp_binomial <- function(n, Ntr, theta, propto) {
  f <- "binomial_lpmf"
  chk_nonnegative(f, "Population size parameter", Ntr); chk_prob(f, "Probability parameter", theta)
  theta <- flat(as_dbl(theta))
  N <- vlen_check(f, n, Ntr, theta)
  nn <- rep_len(as.numeric(n), N); NN <- rep_len(as.numeric(Ntr), N)
  if (any(nn < 0 | nn > NN)) stop_domain(sprintf("%s: Successes variable must be in [0, N]", f))
  lp <- 0
  if (!propto) lp <- lp + sum(lchoose(NN, nn))
  if (inc(propto, theta)) {
    if (length(val(theta)) == 1L) {
      lp <- lp + ad_mult_log(sum(nn), theta) + ad_mult_log(sum(NN - nn), 1 - theta)
    } else {
      lp <- lp + ad_sum(ad_mult_log(nn, theta) + ad_mult_log(NN - nn, 1 - theta))
    }
  }
  lp
}
binom_cdf_ad <- function(n, Ntr, theta, lower) {
  theta <- flat(as_dbl(theta))
  N <- vlen_check("binomial_lcdf", n, Ntr, theta)
  nn <- rep_len(n, N); NN <- rep_len(Ntr, N); vt <- rep_len(val(theta), N)
  lp <- stats::pbinom(nn, NN, vt, lower.tail = lower, log.p = TRUE)
  out <- sum(lp)
  if (!is_ad(theta)) return(out)
  # d/dtheta P(X <= n) = -N * dbinom(n, N - 1, theta)
  # the ratio on the log scale, so extreme tails do not underflow to 0 / 0
  d <- -NN * exp(suppressWarnings(stats::dbinom(nn, NN - 1, vt, log = TRUE)) - lp) * (if (lower) 1 else -1)
  d[!is.finite(d)] <- 0
  ad_node(out, list(theta), function(g) list(unbroadcast(g * d, val(theta))))
}
defdist("binomial", c("ints", "ints", "reals"), lp_binomial, discrete = TRUE,
  lcdf = function(n, N, th) binom_cdf_ad(n, N, th, TRUE),
  lccdf = function(n, N, th) binom_cdf_ad(n, N, th, FALSE),
  rng = function(N, th) vrng(function(n, N, p) as.integer(stats::rbinom(n, N, p)), N, th))

lp_binomial_logit <- function(n, Ntr, alpha, propto) {
  f <- "binomial_logit_lpmf"
  chk_finite(f, "Probability parameter", alpha)
  N <- vlen_check(f, n, Ntr, alpha)
  if (N == 0L) return(0)
  nn <- rep_len(as.numeric(n), N); NN <- rep_len(as.numeric(Ntr), N)
  if (any(nn < 0 | nn > NN)) stop_domain(sprintf("%s: Successes variable must be in [0, N]", f))
  va <- as.vector(val(alpha))
  lp <- 0
  if (!propto) lp <- lp + sum(lchoose(NN, nn))
  if (inc(propto, alpha)) lp <- lp + sum(nn * stats::plogis(va, log.p = TRUE) + (NN - nn) * stats::plogis(-va, log.p = TRUE))
  if (!is_ad(alpha)) return(lp)
  fused(lp, list(alpha), function() list(nn - NN * stats::plogis(va)), N)
}
defdist("binomial_logit", c("ints", "ints", "reals"), lp_binomial_logit, discrete = TRUE,
  rng = function(N, a) vrng(function(n, N, a) as.integer(stats::rbinom(n, N, stats::plogis(a))), N, a))
for (xt in c("matrix", "row_vector")) {
  def("binomial_logit_glm_lpmf", c("ints", "ints", xt, "reals", "vector"), "real", function(n, N, x, a, b) lp_binomial_logit(n, N, glm_linpred("binomial_logit_glm_lpmf", x, a, b), FALSE), density = TRUE, dist = "binomial_logit_glm")
  def("binomial_logit_glm_lupmf", c("ints", "ints", xt, "reals", "vector"), "real", function(n, N, x, a, b) lp_binomial_logit(n, N, glm_linpred("binomial_logit_glm_lpmf", x, a, b), TRUE), density = TRUE, dist = "binomial_logit_glm")
}

lp_beta_binomial <- function(n, Ntr, a, b, propto) {
  f <- "beta_binomial_lpmf"
  chk_positive_finite(f, "First prior sample size parameter", a); chk_positive_finite(f, "Second prior sample size parameter", b)
  a <- flat(as_dbl(a)); b <- flat(as_dbl(b))
  N <- vlen_check(f, n, Ntr, a, b)
  nn <- rep_len(as.numeric(n), N); NN <- rep_len(as.numeric(Ntr), N)
  if (any(nn < 0 | nn > NN)) return(-Inf)
  lp <- 0
  if (!propto) lp <- lp + sum(lchoose(NN, nn))
  if (inc(propto, a, b)) lp <- lp + ad_sum(flat(ad_lbeta(nn + a, NN - nn + b) - ad_lbeta(a, b)))
  lp
}
bb_cdf <- function(n, Ntr, a, b, lower) {
  N <- vlen_check("beta_binomial_lcdf", n, Ntr, a, b)
  nn <- rep_len(n, N); NN <- rep_len(Ntr, N); va <- rep_len(val(a), N); vb <- rep_len(val(b), N)
  f <- function(va, vb) {
    vapply(seq_len(N), function(i) {
      k <- 0:NN[i]
      lpk <- lchoose(NN[i], k) + lbeta(k + va[i], NN[i] - k + vb[i]) - lbeta(va[i], vb[i])
      sel <- if (lower) k <= nn[i] else k > nn[i]
      if (!any(sel)) return(-Inf)
      m <- max(lpk[sel])
      m + log(sum(exp(lpk[sel] - m)))
    }, 0)
  }
  lp <- f(va, vb)
  out <- sum(lp)
  if (!any_ad(a, b)) return(out)
  ad_node(out, list(a, b), function(g) list(
    if (is_ad(a)) unbroadcast(g * num_partial(function(t) f(t, vb), va), val(a)),
    if (is_ad(b)) unbroadcast(g * num_partial(function(t) f(va, t), vb), val(b))))
}
defdist("beta_binomial", c("ints", "ints", "reals", "reals"), lp_beta_binomial, discrete = TRUE,
  lcdf = function(n, N, a, b) bb_cdf(n, N, as_dbl(a), as_dbl(b), TRUE),
  lccdf = function(n, N, a, b) bb_cdf(n, N, as_dbl(a), as_dbl(b), FALSE),
  rng = function(N, a, b) vrng(function(n, N, a, b) as.integer(stats::rbinom(n, N, stats::rbeta(n, a, b))), N, a, b))

lp_hypergeometric <- function(n, Ntot, a, b, propto) {
  f <- "hypergeometric_lpmf"
  if (n < max(0L, Ntot - b) || n > min(a, Ntot)) stop_domain(sprintf("%s: Successes variable is out of range", f))
  lchoose(a, n) + lchoose(b, Ntot - n) - lchoose(a + b, Ntot)
}
defdist("hypergeometric", c("int", "int", "int", "int"), lp_hypergeometric, discrete = TRUE,
  rng = function(N, a, b) as.integer(stats::rhyper(1L, a, b, N)))

lp_categorical <- function(y, theta, propto) {
  f <- "categorical_lpmf"
  K <- nelem(theta)
  chk_int_range(f, "Number of categories", y, 1L, K)
  m <- simplex_check(val(theta))
  if (!is.null(m)) stop_domain(paste(f, ": Probabilities parameter", m))
  if (!inc(propto, theta)) return(0)
  ad_sum(log(ad_gather(as_dbl(theta), as.integer(y))))
}
defdist("categorical", c("ints", "vector"), lp_categorical, discrete = TRUE,
  rng = function(th) {
    m <- simplex_check(val(th))
    if (!is.null(m)) stop_domain(paste("categorical_rng: Probabilities parameter", m))
    sample.int(length(val(th)), 1L, prob = val(th))
  })
lp_categorical_logit <- function(y, beta, propto) {
  f <- "categorical_logit_lpmf"
  K <- nelem(beta)
  chk_int_range(f, "categorical outcome", y, 1L, K)
  chk_finite(f, "log odds parameter", beta)
  if (!inc(propto, beta)) return(0)
  beta <- as_dbl(beta)
  ad_sum(ad_gather(beta, as.integer(y))) - length(y) * ad_log_sum_exp(beta)
}
defdist("categorical_logit", c("ints", "vector"), lp_categorical_logit, discrete = TRUE,
  rng = function(b) sample.int(length(val(b)), 1L, prob = exp(val(b) - max(val(b)))))
cat_logit_glm <- function(y, x, alpha, beta, propto) {
  vx <- val(x)
  X <- if (is.null(dim(vx))) ad_reshape(as_dbl(x), c(1L, length(vx))) else as_dbl(x)
  eta <- op_mat_mat(X, beta)
  n <- nrow(val(eta)); K <- ncol(val(eta))
  eta <- eta + ad_gather(as_dbl(alpha), rep(seq_len(K), each = n), c(n, K))
  yy <- rep_len(as.integer(y), n)
  chk_int_range("categorical_logit_glm_lpmf", "categorical outcome", yy, 1L, K)
  lp <- ad_sum(ad_gather(eta, (yy - 1L) * n + seq_len(n)))
  for (i in seq_len(n)) lp <- lp - ad_log_sum_exp(ad_gather(eta, i + (seq_len(K) - 1L) * n))
  lp
}
for (xt in c("matrix", "row_vector")) {
  def("categorical_logit_glm_lpmf", c("ints", xt, "vector", "matrix"), "real", function(...) cat_logit_glm(..., propto = FALSE), density = TRUE, dist = "categorical_logit_glm")
  def("categorical_logit_glm_lupmf", c("ints", xt, "vector", "matrix"), "real", function(...) cat_logit_glm(..., propto = TRUE), density = TRUE, dist = "categorical_logit_glm")
}

lp_discrete_range <- function(y, lo, hi, propto) {
  f <- "discrete_range_lpmf"
  N <- vlen_check(f, y, lo, hi)
  yy <- rep_len(y, N); ll <- rep_len(lo, N); hh <- rep_len(hi, N)
  if (any(ll > hh)) stop_domain(sprintf("%s: Lower bound must be less than or equal to the upper bound", f))
  if (any(yy < ll | yy > hh)) return(-Inf)
  if (propto) return(0)
  -sum(log(hh - ll + 1))
}
defdist("discrete_range", c("ints", "ints", "ints"), lp_discrete_range, discrete = TRUE,
  cdf = function(y, lo, hi) {
    N <- vlen_check("discrete_range_cdf", y, lo, hi)
    yy <- rep_len(y, N); ll <- rep_len(lo, N); hh <- rep_len(hi, N)
    prod(pmin(pmax((yy - ll + 1) / (hh - ll + 1), 0), 1))
  },
  lccdf = function(y, lo, hi) {
    N <- vlen_check("discrete_range_lccdf", y, lo, hi)
    yy <- rep_len(y, N); ll <- rep_len(lo, N); hh <- rep_len(hi, N)
    sum(log(pmin(pmax((hh - yy) / (hh - ll + 1), 0), 1)))
  },
  rng = function(lo, hi) vrng(function(n, l, h) as.integer(l + floor(stats::runif(n) * (h - l + 1))), lo, hi))

ordered_core <- function(f, k, eta, cuts, link_lcdf, link_lccdf) {
  cv <- val(cuts)
  # an array of cut point vectors: one per observation
  if (!is.null(dim(cv))) {
    if (nrow(cv) != length(k)) stop_shape(sprintf("%s: the number of cut point vectors must equal the number of observations", f))
    eta <- flat(as_dbl(eta))
    # as Stan Math: one location per cut point vector, also for a real location
    if (length(val(eta)) != nrow(cv)) stop_shape(sprintf("%s: the number of locations must equal the number of cut point vectors", f))
    cuts <- as_dbl(cuts)
    lp <- 0
    for (i in seq_along(k)) {
      ci <- ad_gather(cuts, i + (seq_len(ncol(cv)) - 1L) * nrow(cv))
      e <- if (length(val(eta)) == 1L) eta else ad_select(eta, i)
      lp <- lp + ordered_core(f, k[i], e, ci, link_lcdf, link_lccdf)
    }
    return(lp)
  }
  K <- nelem(cuts) + 1L
  chk_int_range(f, "Random variable", k, 1L, K)
  if (K > 2L && any(!(diff(cv) > 0))) stop_domain(sprintf("%s: Cut-points is not a valid ordered vector", f))
  eta <- flat(as_dbl(eta))
  N <- vlen_check(f, k, eta)
  kk <- rep_len(as.integer(k), N)
  cuts <- as_dbl(cuts)
  lp <- 0
  for (i in seq_len(N)) {
    e <- if (length(val(eta)) == 1L) eta else ad_select(eta, i)
    ki <- kk[i]
    term <- if (ki == 1L) {
      link_lccdf(e - ad_select(cuts, 1L))
    } else if (ki == K) {
      link_lcdf(e - ad_select(cuts, K - 1L))
    } else {
      ad_log_diff_exp(link_lcdf(e - ad_select(cuts, ki - 1L)), link_lcdf(e - ad_select(cuts, ki)))
    }
    lp <- lp + term
  }
  lp
}
lp_ordered_logistic <- function(k, eta, cuts, propto) {
  ordered_core("ordered_logistic_lpmf", k, eta, cuts, ad_log_inv_logit, function(x) ad_log1m_inv_logit(x))
}
lp_ordered_probit <- function(k, eta, cuts, propto) {
  ordered_core("ordered_probit_lpmf", k, eta, cuts, ad_log_Phi, function(x) ad_log1m_Phi(x))
}
defdist("ordered_logistic", c("ints", "reals", "vecarr"), lp_ordered_logistic, discrete = TRUE,
  rng = function(eta, c) {
    p <- -diff(c(1, stats::plogis(val(eta) - val(c)), 0))
    sample.int(length(p), 1L, prob = pmax(p, 0))
  })
defdist("ordered_probit", c("ints", "reals", "vecarr"), lp_ordered_probit, discrete = TRUE,
  rng = function(eta, c) {
    p <- -diff(c(1, stats::pnorm(val(eta) - val(c)), 0))
    sample.int(length(p), 1L, prob = pmax(p, 0))
  })

lp_neg_binomial <- function(n, alpha, beta, propto) {
  f <- "neg_binomial_lpmf"
  chk_int_range(f, "Failures variable", n, 0L); chk_positive_finite(f, "Shape parameter", alpha); chk_positive_finite(f, "Inverse scale parameter", beta)
  alpha <- flat(as_dbl(alpha)); beta <- flat(as_dbl(beta))
  N <- vlen_check(f, n, alpha, beta)
  nn <- rep_len(as.numeric(n), N)
  lp <- 0
  if (inc(propto, alpha)) lp <- lp + ad_sum(flat(lgamma(nn + alpha) - lgamma(alpha)))
  if (!propto) lp <- lp - sum(lgamma(nn + 1))
  if (inc(propto, alpha, beta)) {
    lp <- lp + bsum(alpha * (log(beta) - log1p(beta)), N) - ad_sum(flat(nn * log1p(beta)))
  }
  lp
}
defdist("neg_binomial", c("ints", "reals", "reals"), lp_neg_binomial, discrete = TRUE,
  rng = function(a, b) vrng(function(n, a, b) as.integer(stats::rpois(n, stats::rgamma(n, a, b))), a, b))

# lgamma(n + phi) - lgamma(phi) for counts n, accurate for large phi
# (the direct difference cancels); the phi partial sums 1 / (phi + j)
ad_lgamma_ratio <- function(nn, phi) {
  vp <- rep_len(val(phi), length(nn))
  v <- ifelse(nn == 0, 0, lgamma(pmax(nn, 1)) - lbeta(pmax(nn, 1), vp))
  if (!is_ad(phi)) return(v)
  d <- ifelse(nn <= 100, vapply(seq_along(nn), function(i) sum(1 / (vp[i] + seq_len(nn[i]) - 1)), 0),
              digamma(nn + vp) - digamma(vp))
  ad_node(v, list(phi), function(g) list(unbroadcast(g * d, val(phi))))
}

lp_neg_binomial_2 <- function(n, mu, phi, propto) {
  f <- "neg_binomial_2_lpmf"
  chk_int_range(f, "Failures variable", n, 0L); chk_positive_finite(f, "Location parameter", mu); chk_positive_finite(f, "Precision parameter", phi)
  mu <- flat(as_dbl(mu)); phi <- flat(as_dbl(phi))
  N <- vlen_check(f, n, mu, phi)
  nn <- rep_len(as.numeric(n), N)
  lp <- 0
  if (inc(propto, phi)) lp <- lp + ad_sum(ad_lgamma_ratio(nn, phi))
  if (!propto) lp <- lp - sum(lgamma(nn + 1))
  if (inc(propto, mu, phi)) {
    lmp <- log(mu + phi)
    lp <- lp + ad_sum(flat(ad_mult_log(nn, mu) - nn * lmp - phi * log1p(mu / phi) + 0 * nn))
  }
  lp
}
nb2_cdf_ad <- function(n, mu, phi, lower) {
  mu <- flat(as_dbl(mu)); phi <- flat(as_dbl(phi))
  N <- vlen_check("neg_binomial_2_lcdf", n, mu, phi)
  nn <- rep_len(n, N); vm <- rep_len(val(mu), N); vp <- rep_len(val(phi), N)
  fn <- function(vm, vp) stats::pnbinom(nn, size = vp, mu = vm, lower.tail = lower, log.p = TRUE)
  lp <- fn(vm, vp)
  out <- sum(lp)
  if (!any_ad(mu, phi)) return(out)
  ad_node(out, list(mu, phi), function(g) list(
    if (is_ad(mu)) unbroadcast(g * num_partial(function(t) fn(t, vp), vm), val(mu)),
    if (is_ad(phi)) unbroadcast(g * num_partial(function(t) fn(vm, t), vp), val(phi))))
}
defdist("neg_binomial_2", c("ints", "reals", "reals"), lp_neg_binomial_2, discrete = TRUE,
  lcdf = function(n, mu, phi) nb2_cdf_ad(n, mu, phi, TRUE),
  lccdf = function(n, mu, phi) nb2_cdf_ad(n, mu, phi, FALSE),
  rng = function(mu, phi) vrng(function(n, m, p) as.integer(stats::rnbinom(n, size = p, mu = m)), mu, phi))

lp_neg_binomial_2_log <- function(n, eta, phi, propto) {
  f <- "neg_binomial_2_log_lpmf"
  chk_int_range(f, "Failures variable", n, 0L); chk_finite(f, "Log location parameter", eta); chk_positive_finite(f, "Precision parameter", phi)
  eta <- flat(as_dbl(eta)); phi <- flat(as_dbl(phi))
  N <- vlen_check(f, n, eta, phi)
  nn <- rep_len(as.numeric(n), N)
  lp <- 0
  if (inc(propto, phi)) lp <- lp + ad_sum(ad_lgamma_ratio(nn, phi))
  if (!propto) lp <- lp - sum(lgamma(nn + 1))
  if (inc(propto, eta, phi)) {
    lphi <- log(phi)
    lse <- ad_log_sum_exp2(eta, lphi)
    lp <- lp + ad_sum(flat(nn * (eta - lse) - phi * log1p(exp(eta - lphi)) + 0 * nn))
  }
  lp
}
defdist("neg_binomial_2_log", c("ints", "reals", "reals"), lp_neg_binomial_2_log, discrete = TRUE,
  rng = function(eta, phi) vrng(function(n, e, p) as.integer(stats::rnbinom(n, size = p, mu = exp(e))), eta, phi))
for (xt in c("matrix", "row_vector")) {
  def("neg_binomial_2_log_glm_lpmf", c("ints", xt, "reals", "vector", "reals"), "real", function(y, x, a, b, phi) lp_neg_binomial_2_log(y, glm_linpred("neg_binomial_2_log_glm_lpmf", x, a, b), phi, FALSE), density = TRUE, dist = "neg_binomial_2_log_glm")
  def("neg_binomial_2_log_glm_lupmf", c("ints", xt, "reals", "vector", "reals"), "real", function(y, x, a, b, phi) lp_neg_binomial_2_log(y, glm_linpred("neg_binomial_2_log_glm_lpmf", x, a, b), phi, TRUE), density = TRUE, dist = "neg_binomial_2_log_glm")
}

lp_poisson <- function(n, lambda, propto) {
  f <- "poisson_lpmf"
  chk_int_range(f, "Random variable", n, 0L); chk_nonnegative(f, "Rate parameter", lambda)
  if (any(is.infinite(val(lambda)))) stop_domain(sprintf("%s: Rate parameter is inf, but must be finite", f))
  N <- vlen_check(f, n, lambda)
  if (N == 0L) return(0)
  nn <- as.numeric(n); vl <- as.vector(val(lambda))
  if (any(rep_len(vl, N) == 0 & rep_len(nn, N) > 0)) return(-Inf)
  lp <- 0
  if (!propto) lp <- lp - sum(lgamma(nn + 1)) * (if (length(nn) == 1L) N else 1)
  if (inc(propto, lambda)) {
    t <- zmul(nn, log(vl)) - vl
    lp <- lp + sum(t) * (if (length(t) == 1L) N else 1)
  }
  if (!is_ad(lambda)) return(lp)
  fused(lp, list(lambda), function() list(zmul(nn, 1 / vl) - 1), N)
}
pois_cdf_ad <- function(n, lambda, lower) {
  lambda <- flat(as_dbl(lambda))
  N <- vlen_check("poisson_lcdf", n, lambda)
  nn <- rep_len(n, N); vl <- rep_len(val(lambda), N)
  lp <- stats::ppois(nn, vl, lower.tail = lower, log.p = TRUE)
  out <- sum(lp)
  if (!is_ad(lambda)) return(out)
  d <- -exp(stats::dpois(nn, vl, log = TRUE) - lp) * (if (lower) 1 else -1)
  d[!is.finite(d)] <- 0
  ad_node(out, list(lambda), function(g) list(unbroadcast(g * d, val(lambda))))
}
defdist("poisson", c("ints", "reals"), lp_poisson, discrete = TRUE,
  lcdf = function(n, l) pois_cdf_ad(n, l, TRUE),
  lccdf = function(n, l) pois_cdf_ad(n, l, FALSE),
  rng = function(l) {
    if (any(val(l) >= 2^30)) stop_domain("poisson_rng: Rate parameter must be less than 2^30")
    vrng(function(n, l) as.integer(stats::rpois(n, l)), l)
  })

lp_poisson_log <- function(n, alpha, propto) {
  f <- "poisson_log_lpmf"
  chk_int_range(f, "Random variable", n, 0L); chk_not_nan(f, "Log rate parameter", alpha)
  N <- vlen_check(f, n, alpha)
  if (N == 0L) return(0)
  nn <- as.numeric(n); va <- as.vector(val(alpha))
  if (any(va == Inf)) return(-Inf)
  if (any(rep_len(va, N) == -Inf & rep_len(nn, N) != 0)) return(-Inf)
  ea <- exp(va)
  lp <- 0
  if (!propto) lp <- lp - sum(lgamma(nn + 1)) * (if (length(nn) == 1L) N else 1)
  if (inc(propto, alpha)) {
    t <- zmul(nn, va) - ea
    lp <- lp + sum(t) * (if (length(t) == 1L) N else 1)
  }
  if (!is_ad(alpha)) return(lp)
  fused(lp, list(alpha), function() list(nn - ea), N)
}
defdist("poisson_log", c("ints", "reals"), lp_poisson_log, discrete = TRUE,
  rng = function(a) vrng(function(n, a) as.integer(stats::rpois(n, exp(a))), a))
for (xt in c("matrix", "row_vector")) {
  def("poisson_log_glm_lpmf", c("ints", xt, "reals", "vector"), "real", function(y, x, a, b) lp_poisson_log(y, glm_linpred("poisson_log_glm_lpmf", x, a, b), FALSE), density = TRUE, dist = "poisson_log_glm")
  def("poisson_log_glm_lupmf", c("ints", xt, "reals", "vector"), "real", function(y, x, a, b) lp_poisson_log(y, glm_linpred("poisson_log_glm_lpmf", x, a, b), TRUE), density = TRUE, dist = "poisson_log_glm")
}

# ordered logistic regression: ordered_logistic(y | x * beta, cuts), no intercept
for (xt in c("matrix", "row_vector")) local({
  for (suf in c("_lpmf", "_lupmf")) {
    def(paste0("ordered_logistic_glm", suf), c("ints", xt, "vector", "vector"), "real", function(y, x, b, cuts) {
      eta <- glm_linpred("ordered_logistic_glm_lpmf", x, 0, b)
      if (length(y) == 1L && length(val(eta)) > 1L) y <- rep.int(y, length(val(eta)))
      lp_ordered_logistic(y, eta, cuts, FALSE)
    }, density = TRUE, dist = "ordered_logistic_glm")
  }
})

# Poisson binomial (Stan 2.40): number of successes among independent
# Bernoulli trials with probabilities theta, by the log space recursion.
pb_logprobs <- function(theta) {
  theta <- flat(as_dbl(theta))
  K <- nelem(theta)
  lp <- 0  # log P(S = 0) with zero trials
  lt <- log(theta)
  l1t <- ad_log1m(theta)
  for (i in seq_len(K)) {
    a <- ad_c(list(lp + ad_select(l1t, i), -Inf))
    b <- ad_c(list(-Inf, lp + ad_select(lt, i)))
    lp <- ad_log_sum_exp2(a, b)
  }
  lp
}
lp_poisson_binomial <- function(y, theta, propto) {
  f <- "poisson_binomial_lpmf"
  chk_prob(f, "Probability parameters", theta)
  K <- nelem(theta)
  chk_int_range(f, "Successes variable", y, 0L, K)
  lps <- pb_logprobs(theta)
  ad_sum(ad_gather(lps, as.integer(y) + 1L))
}
defdist("poisson_binomial", c("ints", "vecs1"), lp_poisson_binomial, discrete = TRUE,
  lcdf = function(y, theta) {
    lps <- pb_logprobs(theta)
    out <- 0
    K <- nelem(theta)
    for (yy in y) {
      if (yy < 0L) return(-Inf)
      if (yy < K) out <- out + ad_log_sum_exp(ad_gather(lps, seq_len(yy + 1L)))
    }
    out
  },
  lccdf = function(y, theta) {
    lps <- pb_logprobs(theta)
    K <- nelem(theta)
    out <- 0
    for (yy in y) {
      if (yy >= K) return(-Inf)
      if (yy >= 0L) out <- out + ad_log_sum_exp(ad_gather(lps, (yy + 2L):(K + 1L)))
    }
    out
  },
  rng = function(theta) as.integer(sum(stats::runif(nelem(theta)) < val(theta))))

# Yule Simon (Stan 2.39)
lp_yule_simon <- function(y, alpha, propto) {
  f <- "yule_simon_lpmf"
  chk_int_range(f, "Random variable", y, 1L); chk_positive_finite(f, "Shape parameter", alpha)
  alpha <- flat(as_dbl(alpha))
  N <- vlen_check(f, y, alpha)
  yy <- rep_len(as.numeric(y), N)
  bsum(log(alpha), N) + ad_sum(flat(ad_lbeta(yy, alpha + 1) + 0 * yy))
}
ys_lccdf <- function(y, alpha) {
  # P(Y > y) = y * B(y, alpha + 1)
  alpha <- flat(as_dbl(alpha))
  N <- vlen_check("yule_simon_lccdf", y, alpha)
  yy <- rep_len(as.numeric(y), N)
  # the support starts at 1: P(Y > y) = 1 below it
  if (all(yy < 1)) return(0)
  keep <- yy >= 1
  ad_sum(ad_gather(flat(log(pmax(yy, 1)) + ad_lbeta(pmax(yy, 1), alpha + 1) + 0 * yy), which(keep)))
}
defdist("yule_simon", c("ints", "reals"), lp_yule_simon, discrete = TRUE,
  lccdf = ys_lccdf,
  lcdf = function(y, alpha) {
    alpha <- flat(as_dbl(alpha))
    N <- vlen_check("yule_simon_lcdf", y, alpha)
    yy <- rep_len(as.numeric(y), N)
    if (any(yy < 1)) return(-Inf)
    ad_sum(flat(ad_log1m_exp(log(yy) + ad_lbeta(yy, alpha + 1))))
  },
  rng = function(alpha) vrng(function(n, a) {
    w <- stats::rexp(n, a)
    as.integer(stats::rgeom(n, exp(-w)) + 1L)
  }, alpha))

lp_multinomial <- function(ns, theta, propto) {
  f <- "multinomial_lpmf"
  if (any(ns < 0L)) stop_domain(sprintf("%s: Number of trials variable must be nonnegative", f))
  m <- simplex_check(val(theta))
  if (!is.null(m)) stop_domain(paste(f, ": Probabilities parameter", m))
  if (length(ns) != nelem(theta)) stop_shape(sprintf("%s: size of the counts must match the probabilities", f))
  lp <- 0
  if (!propto) lp <- lp + lgamma(sum(ns) + 1) - sum(lgamma(ns + 1))
  if (inc(propto, theta)) lp <- lp + ad_sum(ad_mult_log(as.numeric(ns), as_dbl(theta)))
  lp
}
defdist("multinomial", c("int[]", "vector"), lp_multinomial, discrete = TRUE,
  rng = NULL)
def("multinomial_rng", c("vector", "int"), "int[]", function(theta, N) as.integer(stats::rmultinom(1L, N, val(theta))[, 1L]), rng = TRUE)
lp_multinomial_logit <- function(ns, gamma, propto) {
  lp_multinomial(ns, softmax_impl(gamma), propto)
}
defdist("multinomial_logit", c("int[]", "vector"), lp_multinomial_logit, discrete = TRUE)
def("multinomial_logit_rng", c("vector", "int"), "int[]", function(g, N) as.integer(stats::rmultinom(1L, N, val(softmax_impl(g)))[, 1L]), rng = TRUE)

lp_dirichlet_multinomial <- function(ns, alpha, propto) {
  f <- "dirichlet_multinomial_lpmf"
  chk_positive_finite(f, "Prior size parameter", alpha)
  if (any(ns < 0L)) stop_domain(sprintf("%s: Number of trials variable must be nonnegative", f))
  if (length(ns) != nelem(alpha)) stop_shape(sprintf("%s: size mismatch", f))
  alpha <- as_dbl(alpha)
  Ntot <- sum(ns)
  if (Ntot == 0L) return(0)
  a0 <- ad_sum(alpha)
  lp <- 0
  if (!propto) lp <- lp + lgamma(Ntot + 1) - sum(lgamma(ns + 1))
  if (inc(propto, alpha)) lp <- lp + lgamma(a0) - lgamma(Ntot + a0) + ad_sum(lgamma(ns + alpha) - lgamma(alpha))
  lp
}
defdist("dirichlet_multinomial", c("int[]", "vector"), lp_dirichlet_multinomial, discrete = TRUE)
def("dirichlet_multinomial_rng", c("vector", "int"), "int[]", function(a, N) {
  g <- stats::rgamma(length(val(a)), val(a))
  as.integer(stats::rmultinom(1L, N, g / sum(g))[, 1L])
}, rng = TRUE)

# ---- multivariate continuous -------------------------------------------------------------

# rows of a "vectors" argument as a list of differentiable vectors
vec_rows <- function(x) {
  v <- val(x)
  if (is.null(dim(v))) return(list(as_dbl(x)))
  n <- nrow(v)
  K <- ncol(v)
  lapply(seq_len(n), function(i) ad_gather(as_dbl(x), i + (seq_len(K) - 1L) * n))
}

lp_dirichlet <- function(theta, alpha, propto) {
  f <- "dirichlet_lpdf"
  th <- vec_rows(theta)
  al <- vec_rows(alpha)
  N <- max(length(th), length(al))
  if (length(th) != 1L && length(al) != 1L && length(th) != length(al)) stop_shape(sprintf("%s: size mismatch", f))
  lp <- 0
  for (i in seq_len(N)) {
    t <- th[[min(i, length(th))]]
    a <- al[[min(i, length(al))]]
    chk_positive_finite(f, "prior sample sizes", a)
    m <- simplex_check(val(t))
    if (!is.null(m)) stop_domain(paste(f, ": probabilities", m))
    if (inc(propto, a)) lp <- lp + lgamma(ad_sum(a)) - ad_sum(lgamma(a))
    if (inc(propto, t, a)) lp <- lp + ad_sum(ad_mult_log(a - 1, t))
  }
  lp
}
defdist("dirichlet", c("vectors", "vectors"), lp_dirichlet,
  rng = function(a) {
    g <- stats::rgamma(nelem(a), val(a))
    g / sum(g)
  }, rng_ret = "vector")
.lib$dirichlet_rng[[1L]]$ret <- "vector"

mvn_core <- function(f, y, mu, logdet_half, solve_fn, propto, depends) {
  ys <- vec_rows(y)
  ms <- vec_rows(mu)
  N <- max(length(ys), length(ms))
  if (length(ys) != 1L && length(ms) != 1L && length(ys) != length(ms)) stop_shape(sprintf("%s: size mismatch between the random variable and the location", f))
  K <- nelem(ys[[1L]])
  lp <- 0
  if (!propto) lp <- lp - N * K * 0.5 * log(2 * pi)
  if (inc(propto, depends)) lp <- lp - N * logdet_half
  for (i in seq_len(N)) {
    yi <- ys[[min(i, length(ys))]]
    mi <- ms[[min(i, length(ms))]]
    if (nelem(yi) != K || nelem(mi) != K) stop_shape(sprintf("%s: size mismatch", f))
    chk_not_nan(f, "Random variable", yi); chk_finite(f, "Location parameter", mi)
    if (inc(propto, yi, mi, depends)) {
      d <- as_dbl(yi) - as_dbl(mi)
      lp <- lp - 0.5 * solve_fn(d)
    }
  }
  lp
}
lp_multi_normal <- function(y, mu, Sigma, propto) {
  sym_check("multi_normal_lpdf", Sigma)
  L <- ad_cholesky(Sigma)
  ldh <- ad_sum(log(.lib$diagonal[[1L]]$impl(L)))
  mvn_core("multi_normal_lpdf", y, mu, ldh, function(d) {
    z <- ad_tri_solve(L, d, "multi_normal_lpdf")
    ad_sum(z * z)
  }, propto, Sigma)
}
lp_multi_normal_cholesky <- function(y, mu, L, propto) {
  ldh <- ad_sum(log(.lib$diagonal[[1L]]$impl(L)))
  mvn_core("multi_normal_cholesky_lpdf", y, mu, ldh, function(d) {
    z <- ad_tri_solve(L, d, "multi_normal_cholesky_lpdf")
    ad_sum(z * z)
  }, propto, L)
}
lp_multi_normal_prec <- function(y, mu, Omega, propto) {
  sym_check("multi_normal_prec_lpdf", Omega)
  L <- ad_cholesky(Omega)
  ldh <- -ad_sum(log(.lib$diagonal[[1L]]$impl(L)))
  mvn_core("multi_normal_prec_lpdf", y, mu, ldh, function(d) {
    ad_sum(d * op_mat_vec(Omega, d))
  }, propto, Omega)
}
mvn_rng <- function(mu, L) {
  ms <- vec_rows(mu)
  K <- nrow(val(L))
  draws <- lapply(ms, function(m) val(m) + as.vector(val(L) %*% stats::rnorm(K)))
  if (length(draws) == 1L && is.null(dim(val(mu)))) return(draws[[1L]])
  do.call(rbind, draws)
}
mv_rng_ret <- function(types) {
  if (types[[1L]]$nd > 0L) ty("vector", 1L) else T_vec
}
for (pair in list(c("multi_normal", "lp_multi_normal"), c("multi_normal_cholesky", "lp_multi_normal_cholesky"),
                  c("multi_normal_prec", "lp_multi_normal_prec"))) {
  local({
    nm <- pair[1L]
    fn <- get(pair[2L])
    def(paste0(nm, "_lpdf"), c("vectors", "vectors", "matrix"), "real", function(...) fn(..., propto = FALSE), density = TRUE, dist = nm)
    def(paste0(nm, "_lupdf"), c("vectors", "vectors", "matrix"), "real", function(...) fn(..., propto = TRUE), density = TRUE, dist = nm)
  })
}
def("multi_normal_rng", c("vectors", "matrix"), mv_rng_ret, function(mu, S) mvn_rng(mu, t(chol(val(S)))), rng = TRUE)
def("multi_normal_cholesky_rng", c("vectors", "matrix"), mv_rng_ret, function(mu, L) mvn_rng(mu, val(L)), rng = TRUE)

mvt_core <- function(f, y, nu, mu, L, propto) {
  ys <- vec_rows(y)
  ms <- vec_rows(mu)
  N <- max(length(ys), length(ms))
  if (!(length(ys) %in% c(1L, N) && length(ms) %in% c(1L, N)))
    stop_shape(sprintf("%s: %d random variables and %d locations", f, length(ys), length(ms)))
  K <- nelem(ys[[1L]])
  chk_positive(f, "Degrees of freedom parameter", nu)
  nu <- as_dbl(nu)
  lp <- 0
  if (inc(propto, nu)) lp <- lp + N * (lgamma((nu + K) / 2) - lgamma(nu / 2) - K / 2 * log(nu))
  if (!propto) lp <- lp - N * K / 2 * log(pi)
  if (inc(propto, L)) lp <- lp - N * ad_sum(log(.lib$diagonal[[1L]]$impl(L)))
  for (i in seq_len(N)) {
    d <- as_dbl(ys[[min(i, length(ys))]]) - as_dbl(ms[[min(i, length(ms))]])
    z <- ad_tri_solve(L, d, f)
    lp <- lp - (nu + K) / 2 * log1p(ad_sum(z * z) / nu)
  }
  lp
}
lp_multi_student_t <- function(y, nu, mu, Sigma, propto) {
  sym_check("multi_student_t_lpdf", Sigma)
  mvt_core("multi_student_t_lpdf", y, nu, mu, ad_cholesky(Sigma), propto)
}
lp_multi_student_t_cholesky <- function(y, nu, mu, L, propto) mvt_core("multi_student_t_cholesky_lpdf", y, nu, mu, as_dbl(L), propto)
for (pair in list(c("multi_student_t", "lp_multi_student_t"), c("multi_student_t_cholesky", "lp_multi_student_t_cholesky"))) {
  local({
    nm <- pair[1L]
    fn <- get(pair[2L])
    def(paste0(nm, "_lpdf"), c("vectors", "real", "vectors", "matrix"), "real", function(...) fn(..., propto = FALSE), density = TRUE, dist = nm)
    def(paste0(nm, "_lupdf"), c("vectors", "real", "vectors", "matrix"), "real", function(...) fn(..., propto = TRUE), density = TRUE, dist = nm)
  })
}
mvt_rng <- function(nu, mu, L) {
  ms <- vec_rows(mu)
  K <- nrow(val(L))
  draws <- lapply(ms, function(m) {
    w <- stats::rchisq(1L, val(nu)) / val(nu)
    val(m) + as.vector(val(L) %*% stats::rnorm(K)) / sqrt(w)
  })
  if (length(draws) == 1L && is.null(dim(val(mu)))) return(draws[[1L]])
  do.call(rbind, draws)
}
mvt_rng_ret <- function(types) if (types[[2L]]$nd > 0L) ty("vector", 1L) else T_vec
def("multi_student_t_rng", c("real", "vectors", "matrix"), mvt_rng_ret, function(nu, mu, S) mvt_rng(nu, mu, t(chol(val(S)))), rng = TRUE)
def("multi_student_t_cholesky_rng", c("real", "vectors", "matrix"), mvt_rng_ret, function(nu, mu, L) mvt_rng(nu, mu, val(L)), rng = TRUE)

# LKJ normalizing constant (Lewandowski, Kurowicka and Joe 2009, section 3.2)
lkj_log_const <- function(eta, K) {
  out <- 0
  for (k in seq_len(K - 1L)) {
    b <- eta + (K - k - 1) / 2
    out <- out + (2 * eta - 2 + K - k) * (K - k) * log(2) + (K - k) * ad_lbeta(b, b)
  }
  -out
}
lp_lkj_corr_cholesky <- function(L, eta, propto) {
  f <- "lkj_corr_cholesky_lpdf"
  chk_positive(f, "Shape parameter", eta)
  m <- tf_cholesky_corr$check(val(L), NULL, NULL)
  if (!is.null(m)) stop_domain(paste(f, ": Random variable", m))
  K <- nrow(val(L))
  eta <- as_dbl(eta)
  lp <- 0
  if (!propto || is_ad(eta)) lp <- lp + lkj_log_const(eta, K)
  if (K > 1L && inc(propto, L, eta)) {
    ld <- log(ad_gather(as_dbl(L), (seq_len(K) - 1L) * K + seq_len(K)))
    w <- K - seq_len(K) + 2 * eta - 2
    lp <- lp + ad_sum(ad_gather(w * ld, 2:K))
  }
  lp
}
defdist("lkj_corr_cholesky", c("matrix", "real"), lp_lkj_corr_cholesky)
def("lkj_corr_cholesky_rng", c("int", "real"), "matrix", function(K, eta) lkj_chol_rng(K, val(eta)), rng = TRUE)
lp_lkj_corr <- function(S, eta, propto) {
  f <- "lkj_corr_lpdf"
  chk_positive(f, "Shape parameter", eta)
  m <- tf_corr_matrix$check(val(S), NULL, NULL)
  if (!is.null(m)) stop_domain(paste(f, ": Correlation matrix", m))
  K <- nrow(val(S))
  eta <- as_dbl(eta)
  lp <- 0
  if (!propto || is_ad(eta)) lp <- lp + lkj_log_const(eta, K)
  if (inc(propto, S, eta)) lp <- lp + (eta - 1) * ad_log_det(S)
  lp
}
defdist("lkj_corr", c("matrix", "real"), lp_lkj_corr)
def("lkj_corr_rng", c("int", "real"), "matrix", function(K, eta) {
  L <- lkj_chol_rng(K, val(eta))
  L %*% t(L)
}, rng = TRUE)
# onion method draw of a correlation Cholesky factor
lkj_chol_rng <- function(K, eta) {
  if (K == 1L) return(matrix(1, 1, 1))
  alpha <- eta + (K - 2) / 2
  cpcs <- numeric(K * (K - 1) / 2)
  pos <- 0L
  for (i in 1:(K - 1L)) {
    alpha <- eta + (K - 1 - i) / 2
    for (j in (i + 1L):K) {
      pos <- pos + 1L
      cpcs[pos] <- 2 * stats::rbeta(1L, alpha, alpha) - 1
    }
  }
  val(corr_L_from_cpcs(cpcs, K))
}

lp_wishart <- function(W, nu, S, propto) {
  f <- "wishart_lpdf"
  sym_check(f, W); sym_check(f, S)
  K <- nrow(val(W))
  if (val(nu) <= K - 1) stop_domain(sprintf("%s: Degrees of freedom must be greater than %d", f, K - 1L))
  nu <- as_dbl(nu)
  LS <- ad_cholesky(S)
  LW <- ad_cholesky(W)
  ldS <- 2 * ad_sum(log(.lib$diagonal[[1L]]$impl(LS)))
  ldW <- 2 * ad_sum(log(.lib$diagonal[[1L]]$impl(LW)))
  lp <- 0
  if (inc(propto, nu)) lp <- lp - nu * K / 2 * log(2) - .lib$lmgamma[[1L]]$impl(K, nu / 2)
  if (inc(propto, nu, S)) lp <- lp - nu / 2 * ldS
  if (inc(propto, nu, W)) lp <- lp + (nu - K - 1) / 2 * ldW
  if (inc(propto, W, S)) lp <- lp - 0.5 * ad_sum(ad_solve(S, W, f) * diag(1, K))
  lp
}
defdist("wishart", c("matrix", "real", "matrix"), lp_wishart)
def("wishart_rng", c("real", "matrix"), "matrix", function(nu, S) stats::rWishart(1L, val(nu), val(S))[, , 1L], rng = TRUE)
lp_inv_wishart <- function(W, nu, S, propto) {
  f <- "inv_wishart_lpdf"
  sym_check(f, W); sym_check(f, S)
  K <- nrow(val(W))
  if (val(nu) <= K - 1) stop_domain(sprintf("%s: Degrees of freedom must be greater than %d", f, K - 1L))
  nu <- as_dbl(nu)
  ldS <- 2 * ad_sum(log(.lib$diagonal[[1L]]$impl(ad_cholesky(S))))
  ldW <- 2 * ad_sum(log(.lib$diagonal[[1L]]$impl(ad_cholesky(W))))
  lp <- 0
  if (inc(propto, nu)) lp <- lp - nu * K / 2 * log(2) - .lib$lmgamma[[1L]]$impl(K, nu / 2)
  if (inc(propto, nu, S)) lp <- lp + nu / 2 * ldS
  if (inc(propto, nu, W)) lp <- lp - (nu + K + 1) / 2 * ldW
  if (inc(propto, W, S)) lp <- lp - 0.5 * ad_sum(ad_solve(W, S, f) * diag(1, K))
  lp
}
defdist("inv_wishart", c("matrix", "real", "matrix"), lp_inv_wishart)
def("inv_wishart_rng", c("real", "matrix"), "matrix", function(nu, S) solve(stats::rWishart(1L, val(nu), solve(val(S)))[, , 1L]), rng = TRUE)
# Cholesky parameterizations: density of the factor, which includes the
# Jacobian of W = L L^T, |J| = 2^K prod_k L_kk^(K - k + 1)
chol_jac <- function(L) {
  K <- nrow(val(L))
  ld <- log(ad_gather(as_dbl(L), (seq_len(K) - 1L) * K + seq_len(K)))
  K * log(2) + ad_sum((K - seq_len(K) + 1) * ld)
}
lp_wishart_cholesky <- function(LW, nu, LS, propto) {
  LW <- as_dbl(LW); LS <- as_dbl(LS)
  W <- ad_matmul(LW, ad_transpose(LW))
  S <- ad_matmul(LS, ad_transpose(LS))
  lp_wishart(0.5 * (W + ad_transpose(W)), nu, 0.5 * (S + ad_transpose(S)), propto) + (if (inc(propto, LW)) chol_jac(LW) else 0)
}
defdist("wishart_cholesky", c("matrix", "real", "matrix"), lp_wishart_cholesky)
def("wishart_cholesky_rng", c("real", "matrix"), "matrix", function(nu, LS) {
  t(chol(stats::rWishart(1L, val(nu), val(LS) %*% t(val(LS)))[, , 1L]))
}, rng = TRUE)
lp_inv_wishart_cholesky <- function(LW, nu, LS, propto) {
  LW <- as_dbl(LW); LS <- as_dbl(LS)
  W <- ad_matmul(LW, ad_transpose(LW))
  S <- ad_matmul(LS, ad_transpose(LS))
  lp_inv_wishart(0.5 * (W + ad_transpose(W)), nu, 0.5 * (S + ad_transpose(S)), propto) + (if (inc(propto, LW)) chol_jac(LW) else 0)
}
defdist("inv_wishart_cholesky", c("matrix", "real", "matrix"), lp_inv_wishart_cholesky)
def("inv_wishart_cholesky_rng", c("real", "matrix"), "matrix", function(nu, LS) {
  S <- val(LS) %*% t(val(LS))
  t(chol(solve(stats::rWishart(1L, val(nu), solve(S))[, , 1L])))
}, rng = TRUE)

# ---- other random number generators --------------------------------------------------


# names that defdist() registers but Stan 2.40 does not define
for (nm in c("beta_proportion_cdf", "loglogistic_lcdf", "binomial_logit_rng"))
  if (exists(nm, envir = .lib, inherits = FALSE)) rm(list = nm, envir = .lib)

# ---- cumulative functions of families defined above ------------------------------

defdist("neg_binomial", c("ints", "reals", "reals"), lp_neg_binomial, discrete = TRUE, densities = FALSE,
  lcdf = function(n, a, b) nb2_cdf_ad(n, as_dbl(a) / as_dbl(b), a, TRUE),
  lccdf = function(n, a, b) nb2_cdf_ad(n, as_dbl(a) / as_dbl(b), a, FALSE))

skew_normal_cdf_terms <- function(y, xi, om, al) {
  z <- (as_dbl(y) - as_dbl(xi)) / as_dbl(om)
  ad_Phi(z) - 2 * owens_t_ad(z, as_dbl(al) + 0 * z)
}
defdist("skew_normal", c("reals", "reals", "reals", "reals"), lp_skew_normal, densities = FALSE,
  lcdf = function(y, xi, om, al) ad_sum_v(log(skew_normal_cdf_terms(y, xi, om, al)), y, xi, om, al),
  lccdf = function(y, xi, om, al) ad_sum_v(ad_log1m(skew_normal_cdf_terms(y, xi, om, al)), y, xi, om, al))

# Phi(u) - exp(-v u + v^2 / 2 + log Phi(u - v)) with u = (y - mu) / sigma, v = lambda sigma
exp_mod_normal_cdf_terms <- function(y, mu, s, l) {
  s <- as_dbl(s)
  u <- (as_dbl(y) - as_dbl(mu)) / s
  v <- as_dbl(l) * s
  ad_Phi(u) - exp(-v * u + 0.5 * v * v + ad_log_Phi(u - v))
}
defdist("exp_mod_normal", c("reals", "reals", "reals", "reals"), lp_exp_mod_normal, densities = FALSE,
  lcdf = function(y, mu, s, l) ad_sum_v(log(exp_mod_normal_cdf_terms(y, mu, s, l)), y, mu, s, l),
  lccdf = function(y, mu, s, l) ad_sum_v(ad_log1m(exp_mod_normal_cdf_terms(y, mu, s, l)), y, mu, s, l))

# below mu: log(tau) - 2 (1 - tau) (mu - y) / sigma; above: log1m((1 - tau) exp(-2 tau (y - mu) / sigma))
sde_terms <- function(y, mu, s, t, lower) {
  y <- as_dbl(y); mu <- as_dbl(mu); s <- as_dbl(s); t <- as_dbl(t)
  below <- val(y) <= val(mu)
  d1 <- ad_ifelse(below, mu - y, 0 * y)
  d2 <- ad_ifelse(below, 0 * y, y - mu)
  if (lower) ad_ifelse(below, log(t) - 2 * (1 - t) * d1 / s, ad_log1m((1 - t) * exp(-2 * t * d2 / s)))
  else ad_ifelse(below, ad_log1m(t * exp(-2 * (1 - t) * d1 / s)), ad_log1m(t) - 2 * t * d2 / s)
}
defdist("skew_double_exponential", c("reals", "reals", "reals", "reals"), lp_skew_double_exponential, densities = FALSE,
  lcdf = function(y, mu, s, t) ad_sum_v(sde_terms(y, mu, s, t, TRUE), y, mu, s, t),
  lccdf = function(y, mu, s, t) ad_sum_v(sde_terms(y, mu, s, t, FALSE), y, mu, s, t))

# von Mises cumulative distribution on [-pi, pi], as in Stan Math: a
# continued fraction series for kappa < 49 and a normal approximation above 50
# (the cut and boundary follow Stan Math, including their dependence on mu)
vm_series <- function(x, k) {
  p <- as.integer(28 + 0.5 * val(k) - 100 / (val(k) + 5) + 1)
  s <- sin(x); c <- cos(x); sn <- sin(p * x); cn <- cos(p * x)
  R <- 0; V <- 0
  for (n in seq_len(p - 1L)) {
    sn_tmp <- sn * c - cn * s
    cn <- cn * c + sn * s
    sn <- sn_tmp
    R <- 1 / (2 * (p - n) / k + R)
    V <- R * (sn / (p - n) + V)
  }
  0.5 + x / (2 * pi) + V / pi
}
vm_normal <- function(x, k) {
  b <- sqrt(2 / pi) / ad_unary(k, function(k) besselI(k, 0, expon.scaled = TRUE),
                               function(k, y) besselI(k, 1, expon.scaled = TRUE) - y)
  ad_Phi(b * sin(x / 2))
}
vm_centered <- function(x, k) {
  kv <- val(k)
  if (kv < 49) {
    f <- vm_series(x, k)
    if (val(f) < 0) return(0)
    if (val(f) > 1) return(1)
    return(f)
  }
  if (kv < 50) return((50 - k) * vm_series(x, 49) + (k - 49) * vm_normal(x, 50))
  vm_normal(x, k)
}
vm_cdf_terms <- function(f, y, mu, k) {
  chk_bounded(f, "Random variable", y, -pi, pi); chk_finite(f, "Location parameter", mu); chk_positive(f, "Scale parameter", k)
  N <- vlen_check(f, y, mu, k)
  y <- flat(as_dbl(y)); mu <- flat(as_dbl(mu)); k <- flat(as_dbl(k))
  pick <- function(a, i) if (length(val(a)) == 1L) a else ad_select(a, i)
  lapply(seq_len(N), function(i) {
    yi <- pick(y, i); mi <- pick(mu, i); ki <- pick(k, i)
    if (val(yi) == -pi) return(0)
    if (val(yi) == pi) return(1)
    # center at mu and wrap both into [-pi, pi)
    x2 <- yi - mi
    x2 <- x2 - 2 * pi * floor((val(x2) + pi) / (2 * pi))
    mu2 <- mi - floor((val(mi) + pi) / (2 * pi)) * 2 * pi
    if (val(mu2) < 0) { cut <- val(mu2) + pi; bound <- -pi - mu2 } else { cut <- val(mu2) - pi; bound <- pi - mu2 }
    fb <- vm_centered(bound, ki)
    cdf <- if (val(yi) <= cut) vm_centered(x2, ki) - fb else vm_centered(x2, ki) + 1 - fb
    if (val(cdf) < 0) 0 else cdf
  })
}
defdist("von_mises", c("reals", "reals", "reals"), lp_von_mises, densities = FALSE,
  cdf = function(y, mu, k) {
    out <- 1
    for (t in vm_cdf_terms("von_mises_cdf", y, mu, k)) out <- out * t
    out
  },
  lccdf = function(y, mu, k) {
    out <- 0
    for (t in vm_cdf_terms("von_mises_lccdf", y, mu, k)) out <- out + ad_log1m(t)
    out
  })
# von_mises_lcdf sums the log terms instead of taking the log of the product
.lib$von_mises_lcdf <- NULL
def("von_mises_lcdf", c("reals", "reals", "reals"), "real", function(y, mu, k) {
  out <- 0
  for (t in vm_cdf_terms("von_mises_lcdf", y, mu, k)) out <- out + log(t)
  out
}, dist = "von_mises")

def("bernoulli_lccdf", c("ints", "reals"), "real", function(n, th) {
  chk_prob("bernoulli_lccdf", "Probability parameter", th)
  N <- vlen_check("bernoulli_lccdf", n, th)
  nn <- rep_len(as.integer(n), N)
  if (any(nn < 0L)) return(0)
  if (any(nn >= 1L)) return(-Inf)
  ad_sum_v(log(flat(as_dbl(th))), n, th)
}, dist = "bernoulli")

# ---- beta negative binomial ----------------------------------------------------

lp_beta_neg_binomial <- function(n, r, a, b, propto) {
  f <- "beta_neg_binomial_lpmf"
  chk_int_range(f, "Failures variable", n, 0L)
  chk_positive_finite(f, "Number of successes parameter", r)
  chk_positive_finite(f, "Prior success parameter", a)
  chk_positive_finite(f, "Prior failure parameter", b)
  N <- vlen_check(f, n, r, a, b)
  nn <- rep_len(as.numeric(n), N)
  r <- flat(as_dbl(r)); a <- flat(as_dbl(a)); b <- flat(as_dbl(b))
  lp <- 0
  if (!propto) lp <- lp - sum(lgamma(nn + 1))
  if (inc(propto, r, a, b))
    lp <- lp + bsum(flat(ad_lbeta(nn + r, a + b) + lgamma(nn + b) - ad_lbeta(r, a) - lgamma(b) + 0 * nn), N)
  lp
}
# log P(N <= n) per element by summing the mass function
# ponytail: O(n) per element; Stan uses a 3F2 closed form, switch if large counts matter
bnb_lcdf_terms <- function(n, r, a, b) {
  N <- vlen_check("beta_neg_binomial_lcdf", n, r, a, b)
  nn <- rep_len(as.integer(n), N)
  r <- flat(as_dbl(r)); a <- flat(as_dbl(a)); b <- flat(as_dbl(b))
  pick <- function(x, i) if (length(val(x)) == 1L) x else ad_select(x, i)
  lapply(seq_len(N), function(i) {
    if (nn[i] < 0L) return(-Inf)
    k <- 0:nn[i]
    ri <- pick(r, i); ai <- pick(a, i); bi <- pick(b, i)
    ad_log_sum_exp(ad_lbeta(k + ri, ai + bi) + lgamma(k + bi) - ad_lbeta(ri, ai) - lgamma(bi) - lgamma(k + 1))
  })
}
defdist("beta_neg_binomial", c("ints", "reals", "reals", "reals"), lp_beta_neg_binomial, discrete = TRUE,
  lcdf = function(n, r, a, b) ad_sum_list(bnb_lcdf_terms(n, r, a, b)),
  lccdf = function(n, r, a, b) ad_sum_list(lapply(bnb_lcdf_terms(n, r, a, b), ad_log1m_exp)),
  rng = function(r, a, b) vrng(function(n, r, a, b) as.integer(stats::rnbinom(n, size = r, prob = stats::rbeta(n, a, b))), r, a, b))

# ---- multivariate Gaussian process and dynamic linear models ---------------------

# rows of y (K x N) are independent normal with covariance Sigma / w[k]
multi_gp_core <- function(f, y, w, logdet, quad, propto, N) {
  chk_positive_finite(f, "Kernel scales", w)
  chk_finite(f, "Random variable", y)
  K <- nrow(val(y))
  if (nelem(w) != K) stop_shape(sprintf("%s: the size of the kernel scales must equal the number of rows of y", f))
  w <- flat(as_dbl(w))
  lp <- 0
  if (!propto) lp <- lp - 0.5 * K * N * log(2 * pi)
  if (inc(propto, logdet)) lp <- lp - 0.5 * K * logdet
  if (inc(propto, w)) lp <- lp + 0.5 * N * ad_sum(log(w))
  lp - 0.5 * ad_sum(w * quad)
}
lp_multi_gp <- function(y, Sigma, w, propto) {
  f <- "multi_gp_lpdf"
  sym_check(f, Sigma)
  N <- ncol(val(y))
  if (any(dim(val(Sigma)) != N)) stop_shape(sprintf("%s: Sigma must be N x N with N the number of columns of y", f))
  L <- ad_cholesky(as_dbl(Sigma))
  Z <- ad_tri_solve(L, ad_transpose(as_dbl(y)), f)
  quad <- ad_matmul(matrix(1, 1L, N), Z * Z)
  multi_gp_core(f, y, w, 2 * ad_sum(log(.lib$diagonal[[1L]]$impl(L))), quad, propto, N)
}
lp_multi_gp_cholesky <- function(y, L, w, propto) {
  f <- "multi_gp_cholesky_lpdf"
  N <- ncol(val(y))
  if (any(dim(val(L)) != N)) stop_shape(sprintf("%s: L must be N x N with N the number of columns of y", f))
  Z <- ad_tri_solve(as_dbl(L), ad_transpose(as_dbl(y)), f)
  quad <- ad_matmul(matrix(1, 1L, N), Z * Z)
  multi_gp_core(f, y, w, 2 * ad_sum(log(.lib$diagonal[[1L]]$impl(L))), quad, propto, N)
}
defdist("multi_gp", c("matrix", "matrix", "vector"), lp_multi_gp)
defdist("multi_gp_cholesky", c("matrix", "matrix", "vector"), lp_multi_gp_cholesky)

# Kalman filter likelihood of y (r x T) for the state space model
#   y_t ~ N(F' theta_t, V), theta_t ~ N(G theta_{t-1}, W), theta_0 ~ N(m0, C0)
lp_gaussian_dlm_obs <- function(y, Fm, G, V, W, m0, C0, propto) {
  f <- "gaussian_dlm_obs_lpdf"
  chk_finite(f, "y", y); chk_finite(f, "F", Fm); chk_finite(f, "G", G); chk_finite(f, "m0", m0)
  sym_check(f, W); sym_check(f, C0)
  r <- nrow(val(y)); Tn <- ncol(val(y)); n <- nrow(val(G))
  if (any(dim(val(Fm)) != c(n, r)) || ncol(val(G)) != n || any(dim(val(W)) != n) || nelem(m0) != n || any(dim(val(C0)) != n))
    stop_shape(sprintf("%s: size mismatch", f))
  if (is.null(dim(val(V)))) {
    chk_nonnegative(f, "V", V)
    if (nelem(V) != r) stop_shape(sprintf("%s: size mismatch", f))
    Vm <- ad_scatter(matrix(0, r, r), (seq_len(r) - 1L) * (r + 1L) + 1L, as_dbl(V))
  } else {
    sym_check(f, V)
    if (any(dim(val(V)) != r)) stop_shape(sprintf("%s: size mismatch", f))
    Vm <- as_dbl(V)
  }
  y <- as_dbl(y); Fm <- as_dbl(Fm); G <- as_dbl(G); W <- as_dbl(W)
  Ft <- ad_transpose(Fm); Gt <- ad_transpose(G)
  m <- ad_reshape(as_dbl(m0), c(n, 1L)); C <- as_dbl(C0)
  lp <- if (propto) 0 else -0.5 * log(2 * pi) * r * Tn
  for (t in seq_len(Tn)) {
    a <- ad_matmul(G, m)
    R <- ad_matmul(ad_matmul(G, C), Gt) + W
    fc <- ad_matmul(Ft, a)
    FR <- ad_matmul(Ft, R)
    Q <- ad_matmul(FR, Fm) + Vm
    Q <- 0.5 * (Q + ad_transpose(Q))
    e <- ad_gather(y, (t - 1L) * r + seq_len(r), c(r, 1L)) - fc
    Qe <- ad_solve(Q, e)
    QFR <- ad_solve(Q, FR)
    m <- a + ad_matmul(ad_transpose(FR), Qe)
    C <- R - ad_matmul(ad_transpose(FR), QFR)
    C <- 0.5 * (C + ad_transpose(C))
    lp <- lp - 0.5 * (ad_log_det(Q) + ad_sum(e * Qe))
  }
  lp
}
for (vt in c("vector", "matrix")) local({
  v_t <- vt
  def("gaussian_dlm_obs_lpdf", c("matrix", "matrix", "matrix", v_t, "matrix", "vector", "matrix"), "real",
      function(...) lp_gaussian_dlm_obs(..., propto = FALSE), density = TRUE, dist = "gaussian_dlm_obs")
  def("gaussian_dlm_obs_lupdf", c("matrix", "matrix", "matrix", v_t, "matrix", "vector", "matrix"), "real",
      function(...) lp_gaussian_dlm_obs(..., propto = TRUE), density = TRUE, dist = "gaussian_dlm_obs")
})

# lognormal standard deviations and an LKJ correlation (deprecated in Stan)
lp_lkj_cov <- function(y, mu, sigma, eta, propto) {
  f <- "lkj_cov_lpdf"
  sym_check(f, y)
  K <- nrow(val(y))
  if (nelem(mu) != K || nelem(sigma) != K) stop_shape(sprintf("%s: size mismatch", f))
  y <- as_dbl(y)
  d <- .lib$diagonal[[1L]]$impl(y)
  sds <- sqrt(d)
  lp <- lp_lognormal(sds, mu, sigma, propto)
  isd <- 1 / sds
  D <- ad_gather(isd, rep(seq_len(K), K), c(K, K)) * ad_gather(isd, rep(seq_len(K), each = K), c(K, K))
  lp + lp_lkj_corr(D * y, eta, propto)
}
defdist("lkj_cov", c("matrix", "vector", "vector", "real"), lp_lkj_cov)

# ---- hidden Markov models -------------------------------------------------------
# log_omega[k, n]: log density of observation n in state k; Gamma[i, j]:
# probability of moving from state i to state j; rho: initial distribution.

hmm_check <- function(f, lo, G, rho) {
  K <- nrow(val(lo))
  if (any(dim(val(G)) != K) || nelem(rho) != K) stop_shape(sprintf("%s: Gamma must be K x K and rho of size K, with K the rows of log_omega", f))
  vG <- val(G)
  if (any(vG < 0) || any(abs(rowSums(vG) - 1) > 1e-8)) stop_domain(sprintf("%s: the rows of Gamma must be simplexes", f))
  m <- simplex_check(val(rho))
  if (!is.null(m)) stop_domain(sprintf("%s: rho %s", f, m))
  chk_not_nan(f, "log_omega", lo)
}

# scaled forward pass; alphas[, n] are the filtered state probabilities
hmm_forward <- function(lo, G, rho) {
  K <- nrow(val(lo)); N <- ncol(val(lo))
  lo <- as_dbl(lo); Gt <- ad_transpose(as_dbl(G))
  col <- function(n) ad_gather(lo, (n - 1L) * K + seq_len(K))
  lp <- 0
  a <- flat(as_dbl(rho))
  alphas <- vector("list", N)
  for (n in seq_len(N)) {
    if (n > 1L) a <- flat(ad_matmul(Gt, ad_reshape(a, c(K, 1L))))
    ln <- col(n)
    mx <- max(val(ln))
    w <- a * exp(ln - mx)
    s <- ad_sum(w)
    lp <- lp + log(s) + mx
    a <- w / s
    alphas[[n]] <- a
  }
  list(lp = lp, alphas = alphas)
}

def("hmm_marginal", c("matrix", "matrix", "vector"), "real", function(lo, G, rho) {
  hmm_check("hmm_marginal", lo, G, rho)
  if (ncol(val(lo)) == 0L) return(0)
  hmm_forward(lo, G, rho)$lp
})

def("hmm_hidden_state_prob", c("matrix", "matrix", "vector"), "matrix", function(lo, G, rho) {
  hmm_check("hmm_hidden_state_prob", lo, G, rho)
  vlo <- val(lo); vG <- val(G)
  K <- nrow(vlo); N <- ncol(vlo)
  if (!N) return(matrix(0, K, 0L))
  al <- vapply(hmm_forward(vlo, vG, val(rho))$alphas, val, numeric(K))
  if (K == 1L) al <- matrix(al, 1L)
  out <- matrix(0, K, N)
  beta <- rep(1, K)
  out[, N] <- al[, N]
  if (N > 1L) for (n in (N - 1L):1L) {
    om <- exp(vlo[, n + 1L] - max(vlo[, n + 1L]))
    beta <- as.vector(vG %*% (om * beta))
    beta <- beta / sum(beta)
    p <- al[, n] * beta
    out[, n] <- p / sum(p)
  }
  out
})

def("hmm_latent_rng", c("matrix", "matrix", "vector"), "int[]", function(lo, G, rho) {
  hmm_check("hmm_latent_rng", lo, G, rho)
  vlo <- val(lo); vG <- val(G)
  K <- nrow(vlo); N <- ncol(vlo)
  if (!N) return(integer(0))
  al <- vapply(hmm_forward(vlo, vG, val(rho))$alphas, val, numeric(K))
  if (K == 1L) al <- matrix(al, 1L)
  z <- integer(N)
  z[N] <- sample.int(K, 1L, prob = al[, N])
  if (N > 1L) for (n in (N - 1L):1L) {
    p <- al[, n] * vG[, z[n + 1L]]
    z[n] <- sample.int(K, 1L, prob = p)
  }
  z
}, rng = TRUE)

# matrix normal with precision matrices Sigma (rows) and D (columns); the
# log determinant weights follow Stan Math
lp_matrix_normal_prec <- function(y, Mu, Sigma, D, propto) {
  f <- "matrix_normal_prec_lpdf"
  sym_check(f, Sigma); sym_check(f, D)
  chk_finite(f, "Random variable", y); chk_finite(f, "Location parameter", Mu)
  N <- nrow(val(y)); P <- ncol(val(y))
  if (any(dim(val(Mu)) != c(N, P)) || any(dim(val(Sigma)) != N) || any(dim(val(D)) != P)) stop_shape(sprintf("%s: size mismatch", f))
  lp <- if (propto) 0 else -0.5 * log(2 * pi) * N * P
  if (inc(propto, Sigma)) lp <- lp + 0.5 * N * ad_log_det(as_dbl(Sigma))
  if (inc(propto, D)) lp <- lp + 0.5 * P * ad_log_det(as_dbl(D))
  if (inc(propto, y, Mu, Sigma, D)) {
    E <- as_dbl(y) - as_dbl(Mu)
    lp <- lp - 0.5 * ad_sum(as_dbl(D) * ad_matmul(ad_matmul(ad_transpose(E), as_dbl(Sigma)), E))
  }
  lp
}
defdist("matrix_normal_prec", c("matrix", "matrix", "matrix", "matrix"), lp_matrix_normal_prec)
