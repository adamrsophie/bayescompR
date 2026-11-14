# Wiener first passage time distributions (drift diffusion models).
#
# wiener_lpdf(y | a, t0, w, v) is the density of hitting the upper boundary
# at time y, computed as Stan Math does: the series of Navarro and Fuss
# (2009) for the standardized zero drift density, with Stan's choice of the
# small or large time series and its number of terms. The five parameter
# form adds the drift variability sv in closed form; the seven parameter
# form averages over a uniform starting point (width sw, Gauss-Legendre
# quadrature) and a uniform non decision time (width st0, through the
# difference of the distribution function). The defective distribution
# functions (wiener_lcdf_unnorm, wiener_lccdf_unnorm) use the large time
# series, with Gauss-Hermite quadrature over the drift when sv > 0. All
# derivatives come from automatic differentiation of these expressions.

gauss_rule <- function(n, kind) {
  i <- seq_len(n - 1L)
  if (kind == "legendre") {
    off <- i / sqrt(4 * i^2 - 1)
    J <- diag(0, n); J[cbind(i, i + 1L)] <- off; J[cbind(i + 1L, i)] <- off
    e <- eigen(J, symmetric = TRUE)
    list(x = e$values, w = 2 * e$vectors[1L, ]^2)
  } else {
    off <- sqrt(i / 2)
    J <- diag(0, n); J[cbind(i, i + 1L)] <- off; J[cbind(i + 1L, i)] <- off
    e <- eigen(J, symmetric = TRUE)
    list(x = e$values, w = sqrt(pi) * e$vectors[1L, ]^2)
  }
}
wiener_gl <- gauss_rule(16L, "legendre")
wiener_gh <- gauss_rule(24L, "hermite")

# log of the standardized zero drift density f(x | w) at the lower boundary
# for x = t / a^2, with the truncation of Stan Math's wiener_lpdf
wiener_log_f0 <- function(x, w) {
  xv <- val(x)
  err <- 1e-6
  sx <- sqrt(xv)
  kl <- ifelse(pi * err * xv < 1, pmax(sqrt(-2 * sqrt(pi) * (log(pi) - 6 * log(10) + log(xv))) / sx, sx / pi), sx / pi)
  e0 <- 2 * sqrt(2 * pi) * err * sx
  ks <- ifelse(e0 < 1, pmax(2 + sx * sqrt(-2 * log(pmin(e0, 1))), sx + 1), 2)
  small <- ks < kl
  lts <- NULL; ltl <- NULL
  if (any(small)) {
    K <- ceiling(ks)
    lo <- -floor((K - 1) / 2); hi <- ceiling((K - 1) / 2)
    tmp <- 0
    for (k in min(lo[small]):max(hi[small])) {
      m <- as.numeric(small & k >= lo & k <= hi)
      if (!any(m > 0)) next
      u <- w + 2 * k
      tmp <- tmp + m * (u * exp(-u * u * 0.5 / x))
    }
    lts <- log(tmp + as.numeric(!small)) - (log(2) / 2 + log(sqrt(pi))) - 1.5 * log(x)
  }
  if (any(!small)) {
    K <- ceiling(kl)
    tmp <- 0
    for (k in seq_len(max(K[!small]))) {
      m <- as.numeric(!small & k <= K)
      tmp <- tmp + m * (k * exp(-(k * k) * (pi * pi / 2) * x) * sin(k * pi * w))
    }
    ltl <- log(tmp + as.numeric(small)) + 2 * log(sqrt(pi))
  }
  if (is.null(lts)) return(ltl)
  if (is.null(ltl)) return(lts)
  ad_ifelse(small, lts, ltl)
}

# log density of the upper boundary with drift variability sv
wiener_log_p5 <- function(t, a, w, v, sv) {
  wl <- 1 - w; vl <- 0 - v
  sv2 <- sv * sv
  wiener_log_f0(t / (a * a), wl) - 2 * log(a) +
    (sv2 * a * a * wl * wl - 2 * a * vl * wl - vl * vl * t) / (2 * (1 + sv2 * t)) - 0.5 * log1p(sv2 * t)
}

# probability of absorption at the lower boundary
wiener_p_lower <- function(a, w, v) {
  z <- val(v) * val(a)
  tiny <- abs(z) < 1e-10
  vs <- v + as.numeric(tiny) * 1e-3
  pos <- val(vs) > 0
  P <- ad_ifelse(pos, (exp(-2 * vs * a * w) - exp(-2 * vs * a)) / (0 - expm1(-2 * vs * a)),
                 expm1(2 * vs * a * (1 - w)) / expm1(2 * vs * a))
  if (any(tiny)) P <- ad_ifelse(tiny, 1 - w + 0 * a * v, P)
  P
}

# sum of x over groups g (1..n), differentiable
ad_group_sum <- function(x, g, n) {
  vx <- as.vector(val(x))
  y <- numeric(n)
  s <- rowsum(vx, g, reorder = TRUE)
  y[as.integer(rownames(s))] <- s[, 1L]
  if (!inherits(x, "bcad")) return(y)
  ad_node(y, list(x), function(gr) list(as.vector(gr)[g]))
}
ex <- function(x, idx) if (inherits(x, "bcad")) ad_gather(x, idx) else as.vector(x)[idx]

# defective distribution function at the lower boundary: the small time
# series (Blurton, Kesselmeier and Gondan 2012) for t < a^2, the large time
# series otherwise, each with enough terms for double precision
wiener_F_lower <- function(t, a, w, v) {
  tv <- as.vector(val(t)); av <- as.vector(val(a))
  pos <- tv > 0
  ts <- t + as.numeric(!pos)
  r <- as.vector(val(ts)) / (av * av)
  small <- r < 1
  Fs <- NULL; Fl <- NULL
  pre <- exp(-v * a * w - v * v * ts / 2)
  if (any(small)) {
    st <- sqrt(ts)
    S <- 0
    for (k in 0:13) {
      rk <- if (k %% 2L == 0L) k * a + a * w else (k + 1) * a - a * w
      u1 <- (rk - v * ts) / st; u2 <- (rk + v * ts) / st
      e0 <- -rk * rk / (2 * ts)
      term <- exp(e0 + u1 * u1 / 2 + ad_log_Phi(0 - u1)) + exp(e0 + u2 * u2 / 2 + ad_log_Phi(0 - u2))
      S <- if (k %% 2L == 0L) S + term else S - term
    }
    Fs <- pre * S
  }
  if (any(!small)) {
    S <- 0
    for (k in 1:8) {
      lam <- k * k * pi * pi / (a * a)
      S <- S + k * sin(k * pi * w) * exp(-lam * ts / 2) / (v * v + lam)
    }
    Fl <- wiener_p_lower(a, w, v) - 2 * pi / (a * a) * pre * S
  }
  F <- if (is.null(Fs)) Fl else if (is.null(Fl)) Fs else ad_ifelse(small, Fs, Fl)
  F * as.numeric(pos)
}

# the upper boundary with drift variability: Gauss-Hermite over the drift
wiener_gh_expand <- function(fun, a, w, v, sv, ...) {
  if (all(val(sv) == 0)) return(fun(a, 1 - w, 0 - v, ...))
  M <- length(val(v)); G <- length(wiener_gh$x)
  idx <- rep(seq_len(M), times = G)
  vj <- ex(v, idx) + sqrt(2) * ex(sv, idx) * rep(wiener_gh$x, each = M)
  extra <- lapply(list(...), ex, idx = idx)
  vals <- do.call(fun, c(list(ex(a, idx), 1 - ex(w, idx), 0 - vj), extra))
  ad_group_sum(vals * rep(wiener_gh$w / sqrt(pi), each = M), idx, M)
}
wiener_F_upper5 <- function(t, a, w, v, sv) wiener_gh_expand(function(a, w, v, t) wiener_F_lower(t, a, w, v), a, w, v, sv, t)
wiener_P_upper5 <- function(a, w, v, sv) wiener_gh_expand(function(a, w, v) wiener_p_lower(a, w, v), a, w, v, sv)

# index, node positions in [-1, 1] and weights of Gauss-Legendre nodes for
# M values (node major)
gl_expand <- function(M) {
  G <- length(wiener_gl$x)
  list(idx = rep(seq_len(M), times = G), x = rep(wiener_gl$x, each = M), w = rep(wiener_gl$w / 2, each = M))
}

# expand the arguments over the starting point interval when sw > 0
wiener_over_w <- function(A) {
  N <- length(val(A$y))
  if (!any(val(A$sw) > 0)) return(list(B = A, q = list(idx = seq_len(N), w = rep(1, N))))
  q <- gl_expand(N)
  B <- lapply(A, ex, idx = q$idx)
  B$w <- B$w + B$sw / 2 * q$x
  list(B = B, q = q)
}
wiener_check <- function(f, y, a, t0, w, v, sv = 0, sw = 0, st0 = 0) {
  chk_positive_finite(f, "Boundary separation", a)
  chk_finite(f, "Nondecision time", t0); chk_nonnegative(f, "Nondecision time", t0)
  chk_bounded(f, "A-priori bias", w, 0, 1)
  chk_finite(f, "Drift rate", v)
  chk_finite(f, "Inter-trial variability in drift rate", sv); chk_nonnegative(f, "Inter-trial variability in drift rate", sv)
  chk_bounded(f, "Inter-trial variability in A-priori bias", sw, 0, 1)
  chk_finite(f, "Inter-trial variability in Nondecision time", st0); chk_nonnegative(f, "Inter-trial variability in Nondecision time", st0)
  wv <- val(w); swv <- val(sw)
  if (any(rep_len(as.vector(wv), max(length(wv), length(swv))) - rep_len(as.vector(swv), max(length(wv), length(swv))) / 2 <= 0) ||
      any(rep_len(as.vector(wv), max(length(wv), length(swv))) + rep_len(as.vector(swv), max(length(wv), length(swv))) / 2 >= 1))
    stop_domain(sprintf("%s: A-priori bias +/- Inter-trial variability in A-priori bias must lie in (0, 1)", f))
  yv <- as.vector(val(y)); t0v <- as.vector(val(t0))
  n <- max(length(yv), length(t0v))
  if (any(rep_len(yv, n) <= rep_len(t0v, n)))
    stop_domain(sprintf("%s: Random variable = %s, but must be greater than nondecision time = %s", f,
                        format_scalar(rep_len(yv, n)[which(rep_len(yv, n) <= rep_len(t0v, n))[1L]]),
                        format_scalar(rep_len(t0v, n)[which(rep_len(yv, n) <= rep_len(t0v, n))[1L]])))
}

lp_wiener4 <- function(y, alpha, tau, beta, delta, propto) {
  f <- "wiener_lpdf"
  chk_positive(f, "Random variable", y)
  chk_positive_finite(f, "Boundary separation", alpha); chk_positive_finite(f, "Nondecision time", tau)
  chk_bounded(f, "A-priori bias", beta, 0, 1); chk_finite(f, "Drift rate", delta)
  y <- flat(as_dbl(y)); alpha <- flat(as_dbl(alpha)); tau <- flat(as_dbl(tau)); beta <- flat(as_dbl(beta)); delta <- flat(as_dbl(delta))
  N <- vlen_check(f, y, alpha, tau, beta, delta)
  if (N == 0L) return(0)
  yv <- rep_len(as.vector(val(y)), N); tv <- rep_len(as.vector(val(tau)), N)
  if (any(yv <= tv)) {
    i <- which(yv <= tv)[1L]
    stop_domain(sprintf("%s: Random variable = %s, but must be greater than nondecision time = %s", f, format_scalar(yv[i]), format_scalar(tv[i])))
  }
  if (!inc(propto, y, alpha, tau, beta, delta)) return(0)
  omb <- 1 - beta
  a2 <- alpha * alpha
  x <- (y - tau) / a2
  bsum(delta * alpha * omb - delta * delta * x * a2 / 2 - log(a2) + wiener_log_f0(bvec(x, N), bvec(omb, N)), N)
}

# broadcast a scalar or vector (possibly AD) to length N
bvec <- function(x, N) {
  if (length(val(x)) == N) return(x)
  if (inherits(x, "bcad")) return(ad_gather(x, rep_len(seq_along(val(x)), N)))
  rep_len(as.vector(x), N)
}

wiener_args <- function(f, args) {
  args <- lapply(args, function(a) flat(as_dbl(a)))
  N <- do.call(vlen_check, c(list(f), args))
  list(args = lapply(args, bvec, N = N), N = N)
}

wiener_named_args <- function(f, y, a, t0, w, v, sv, sw, st0) {
  r <- wiener_args(f, list(y, a, t0, w, v, sv, sw, st0))
  names(r$args) <- c("y", "a", "t0", "w", "v", "sv", "sw", "st0")
  r
}

lp_wiener_full <- function(y, a, t0, w, v, sv, sw = 0, st0 = 0, precision = NULL, propto = FALSE) {
  f <- "wiener_lpdf"
  wiener_check(f, y, a, t0, w, v, sv, sw, st0)
  if (!inc(propto, y, a, t0, w, v, sv, sw, st0)) return(0)
  r <- wiener_named_args(f, y, a, t0, w, v, sv, sw, st0)
  N <- r$N
  if (N == 0L) return(0)
  A <- r$args
  if (all(val(A$sw) == 0) && all(val(A$st0) == 0)) return(sum(wiener_log_p5(A$y - A$t0, A$a, A$w, A$v, A$sv)))
  e <- wiener_over_w(A); B <- e$B; q <- e$q
  p5 <- exp(wiener_log_p5(B$y - B$t0, B$a, B$w, B$v, B$sv))
  stz <- as.vector(val(B$st0)) == 0
  inner <- if (all(stz)) p5 else {
    sts <- B$st0 + as.numeric(stz)
    pst <- (wiener_F_upper5(B$y - B$t0, B$a, B$w, B$v, B$sv) - wiener_F_upper5(B$y - B$t0 - sts, B$a, B$w, B$v, B$sv)) / sts
    if (any(stz)) ad_ifelse(stz, p5, pst) else pst
  }
  sum(log(ad_group_sum(inner * q$w, q$idx, N)))
}

wiener_cdf_full <- function(lower_tail, y, a, t0, w, v, sv = 0, sw = 0, st0 = 0, precision = NULL) {
  f <- if (lower_tail) "wiener_lcdf_unnorm" else "wiener_lccdf_unnorm"
  wiener_check(f, y, a, t0, w, v, sv, sw, st0)
  r <- wiener_named_args(f, y, a, t0, w, v, sv, sw, st0)
  N <- r$N
  if (N == 0L) return(0)
  e <- wiener_over_w(r$args); B <- e$B; q <- e$q
  M <- length(q$idx)
  if (any(val(B$st0) > 0)) {
    q2 <- gl_expand(M)
    C <- lapply(B, ex, idx = q2$idx)
    Fp <- wiener_F_upper5(C$y - C$t0 - C$st0 / 2 * (1 + q2$x), C$a, C$w, C$v, C$sv)
    Fp <- ad_group_sum(Fp * q2$w, q2$idx, M)
  } else {
    Fp <- wiener_F_upper5(B$y - B$t0, B$a, B$w, B$v, B$sv)
  }
  vi <- if (lower_tail) Fp else wiener_P_upper5(B$a, B$w, B$v, B$sv) - Fp
  sum(log(ad_group_sum(vi * q$w, q$idx, N)))
}
defdist("wiener", c("reals", "reals", "reals", "reals", "reals"), lp_wiener4)
# five and seven parameter forms, each with an optional precision argument
# (accepted; the series here are summed well below that precision)
for (nargs in c(6L, 8L)) local({
  na <- nargs
  full <- function(propto) function(y, a, t0, w, v, sv, ...) {
    extra <- list(...)
    if (na == 6L) return(lp_wiener_full(y, a, t0, w, v, sv, propto = propto))
    lp_wiener_full(y, a, t0, w, v, sv, extra[[1L]], extra[[2L]], propto = propto)
  }
  for (k in 0:1) {
    sig <- c(rep("reals", na), if (k) "real")
    def("wiener_lpdf", sig, "real", full(FALSE), density = TRUE, dist = "wiener")
    def("wiener_lupdf", sig, "real", full(TRUE), density = TRUE, dist = "wiener")
  }
})
for (nargs in c(5L, 8L)) local({
  na <- nargs
  wrap <- function(lt) { force(lt); function(y, a, t0, w, v, ...) {
    extra <- list(...)
    if (na == 5L) return(wiener_cdf_full(lt, y, a, t0, w, v, 0, 0, 0))
    wiener_cdf_full(lt, y, a, t0, w, v, extra[[1L]], extra[[2L]], extra[[3L]])
  } }
  for (lt in c(TRUE, FALSE)) {
    force(lt)
    nm <- if (lt) "wiener_lcdf_unnorm" else "wiener_lccdf_unnorm"
    def(nm, rep("reals", na), "real", wrap(lt), dist = "wiener")
    def(nm, c(rep("reals", na), "real"), "real", wrap(lt), dist = "wiener")
  }
})
