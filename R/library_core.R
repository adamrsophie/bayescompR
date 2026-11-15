# The Stan function library: registry, signature patterns and the core
# scalar and container functions.
#
# Every Stan function is registered with one or more signatures. A signature
# lists argument patterns, a return type rule and an R implementation. The
# semantic analyzer resolves calls against these signatures; the code
# generator then calls the implementation directly. Implementations receive
# plain numeric values or differentiable values ("bcad") and must preserve
# Stan semantics, including shape checks.
#
# Argument patterns
#   "int", "real", "vector", "row_vector", "matrix"   exact base types (real accepts int)
#   "int[]", "real[]", "vector[]", ... "real[,]"    arrays with the given dimensions
#   "reals"    real, int, vector, row_vector, or a 1 dimensional array of int or real
#   "ints"     int or a 1 dimensional array of int
#   "vectors"  vector, row_vector or a 1 dimensional array of either
#   "any"      any int or real valued type, including arrays
#   "anyarr"   any array
#   "container" vector, row_vector, matrix or array of reals
#   "string"   a string literal (print, reject)

.lib <- new.env(parent = emptyenv())

def <- function(name, args, ret, impl, ...) {
  sig <- list(name = name, args = args, ret = ret, impl = impl, ...)
  .lib[[name]] <- c(.lib[[name]], list(sig))
  invisible(NULL)
}

parse_ty_string <- function(s) {
  m <- regmatches(s, regexec("^([a-z_]+)(\\[(,*)\\])?$", s))[[1L]]
  if (!length(m)) stop("internal error: bad type string ", s)
  nd <- if (nzchar(m[3L])) nchar(m[4L]) + 1L else 0L
  ty(m[2L], nd)
}

pat_ok <- function(p, t) {
  if (is.null(t)) return(FALSE)
  numeric_bases <- c("int", "real", "vector", "row_vector", "matrix")
  switch(p,
    any = t$base %in% numeric_bases,
    anyarr = t$nd > 0L && t$base %in% numeric_bases,
    container = t$base %in% numeric_bases && !(t$nd == 0L && t$base %in% c("int", "real")),
    reals = (t$nd == 0L && t$base %in% c("int", "real", "vector", "row_vector")) ||
      (t$nd == 1L && t$base %in% c("int", "real")),
    ints = t$base == "int" && t$nd <= 1L,
    vectors = t$nd <= 1L && t$base %in% c("vector", "row_vector"),
    # a vector, row vector or array of reals (Stan's 1 dimensional real containers)
    vecs1 = (t$nd == 0L && t$base %in% c("vector", "row_vector")) || (t$nd == 1L && t$base == "real"),
    intany = t$base == "int",
    vecarr = t$base == "vector" && t$nd <= 1L,
    string = t$base == "string",
    cxall = is_cx_base(t$base),
    cxnum = t$base %in% numeric_bases || is_cx_base(t$base),
    ty_assignable(parse_ty_string(p), t))
}

# a library argument pattern that takes type t without promotion (the family
# patterns such as "reals" stand for one signature per member type)
pat_exact <- function(p, t) {
  if (!pat_ok(p, t)) return(FALSE)
  if (p %in% c("any", "anyarr", "container", "reals", "ints", "vectors", "string")) return(TRUE)
  ty_eq(parse_ty_string(p), t)
}

ret_type <- function(sig, types) {
  r <- sig$ret
  if (is.function(r)) return(r(types))
  parse_ty_string(r)
}

# return rules
r_promote1 <- function(types) {
  t <- types[[1L]]
  if (t$base == "int") ty("real", t$nd) else t
}
r_same1 <- function(types) types[[1L]]
# binary vectorized: the container argument decides
# dimensions of a type counting vectors as one and matrices as two
total_dims <- function(t) t$nd + switch(t$base, vector = , row_vector = 1L, matrix = 2L, 0L)
r_bcast2 <- function(types) {
  a <- types[[1L]]
  b <- types[[2L]]
  if (!is_scalar_ty(a)) {
    if (!is_scalar_ty(b) && !(ty_assignable(r_promote1(list(a)), b) || ty_assignable(r_promote1(list(b)), a))) {
      # an integer array pairs with a real container of as many dimensions
      # (binary_log_loss(array[] int, vector), lmgamma(array[,] int, matrix))
      if (a$base == "int" && b$base != "int" && total_dims(a) == total_dims(b)) return(b)
      if (b$base == "int" && a$base != "int" && total_dims(a) == total_dims(b)) return(a)
      return(NULL)
    }
    return(r_promote1(list(a)))
  }
  if (!is_scalar_ty(b)) return(r_promote1(list(b)))
  T_real
}
r_bcast_n <- function(types) {
  out <- T_real
  for (t in types) if (!is_scalar_ty(t)) {
    tt <- r_promote1(list(t))
    if (!is_scalar_ty(out) && !ty_eq(out, tt)) return(NULL)
    out <- tt
  }
  out
}
r_elem1 <- function(types) {
  t <- types[[1L]]
  if (t$nd > 0L) ty(t$base, t$nd - 1L) else if (t$base == "int") T_int else T_real
}

# ---------------------------------------------------------------------------
# value helpers

as_dbl <- function(x) {
  if (inherits(x, "bcad")) return(x)
  if (is.integer(x) || is.logical(x)) storage.mode(x) <- "double"
  x
}

flat <- function(x) if (inherits(x, "bcad")) ad_reshape(x, NULL) else as.vector(x)

nelem <- function(x) length(val(x))

vlen_check <- function(fname, ...) {
  N <- 1L
  for (k in seq_len(...length())) {
    a <- ...elt(k)
    l <- length(if (inherits(a, "bcad")) a$v else a)
    if (l != 1L) {
      if (N == 1L) N <- l
      else if (l != N) {
        ls <- vapply(list(...), function(a) length(val(a)), 1L)
        stop_shape(sprintf("%s: vectorized arguments have different sizes (%s)", fname,
                           paste(ls, collapse = ", ")))
      }
    }
  }
  N
}

# common shape check for elementwise operations
same_shape <- function(a, b, what) {
  va <- val(a)
  vb <- val(b)
  if (length(va) == 1L && is.null(dim(va))) return(invisible())
  if (length(vb) == 1L && is.null(dim(vb))) return(invisible())
  da <- if (is.null(dim(va))) length(va) else dim(va)
  db <- if (is.null(dim(vb))) length(vb) else dim(vb)
  if (!identical(as.integer(da), as.integer(db)))
    stop_shape(sprintf("%s: size mismatch (%s versus %s)", what,
                       paste(da, collapse = " x "), paste(db, collapse = " x ")))
  invisible()
}

# ---------------------------------------------------------------------------
# Domain checks (Stan throws std::domain_error; we raise bayescomp_domain_error)

chk_finite <- function(f, name, x) {
  v <- val(x)
  if (any(!is.finite(v))) stop_domain(sprintf("%s: %s is %s, but must be finite", f, name, fmt_num(v[!is.finite(v)][1L])))
}
chk_not_nan <- function(f, name, x) {
  v <- val(x)
  if (any(is.nan(v) | is.na(v))) stop_domain(sprintf("%s: %s is nan, but must not be nan", f, name))
}
chk_positive <- function(f, name, x) {
  v <- val(x)
  bad <- !(v > 0)
  if (any(bad | is.na(bad))) stop_domain(sprintf("%s: %s is %s, but must be positive", f, name, fmt_num(v[which(bad | is.na(bad))[1L]])))
}
chk_positive_finite <- function(f, name, x) {
  v <- val(x)
  bad <- !(v > 0 & is.finite(v))
  if (any(bad | is.na(bad))) stop_domain(sprintf("%s: %s is %s, but must be positive and finite", f, name, fmt_num(v[which(bad | is.na(bad))[1L]])))
}
chk_nonnegative <- function(f, name, x) {
  v <- val(x)
  bad <- !(v >= 0)
  if (any(bad | is.na(bad))) stop_domain(sprintf("%s: %s is %s, but must be nonnegative", f, name, fmt_num(v[which(bad | is.na(bad))[1L]])))
}
chk_bounded <- function(f, name, x, lo, hi) {
  v <- val(x)
  bad <- !(v >= lo & v <= hi)
  if (any(bad | is.na(bad))) stop_domain(sprintf("%s: %s is %s, but must be in the interval [%s, %s]", f, name, fmt_num(v[which(bad | is.na(bad))[1L]]), fmt_num(lo), fmt_num(hi)))
}
chk_prob <- function(f, name, x) chk_bounded(f, name, x, 0, 1)

# ---------------------------------------------------------------------------
# Operators. The code generator chooses the variant from static types.

op_add <- function(a, b) {
  same_shape(a, b, "addition")
  as_dbl(a) + as_dbl(b)
}
op_sub <- function(a, b) {
  same_shape(a, b, "subtraction")
  as_dbl(a) - as_dbl(b)
}
op_neg <- function(a) if (is.list(a) && !is_ad(a)) lapply(a, op_neg) else if (is.integer(a)) -a else -as_dbl(a)
op_elt_mul <- function(a, b) {
  same_shape(a, b, "elementwise multiplication")
  as_dbl(a) * as_dbl(b)
}
op_elt_div <- function(a, b) {
  same_shape(a, b, "elementwise division")
  as_dbl(a) / as_dbl(b)
}
op_scale <- function(a, b) as_dbl(a) * as_dbl(b)
op_div_scalar <- function(a, b) as_dbl(a) / as_dbl(b)

int_check <- function(x, what) {
  if (any(is.na(x))) stop_domain(sprintf("integer overflow in %s", what))
  x
}
op_int_add <- function(a, b) int_check(suppressWarnings(a + b), "addition")
op_int_sub <- function(a, b) int_check(suppressWarnings(a - b), "subtraction")
op_int_mul <- function(a, b) int_check(suppressWarnings(a * b), "multiplication")
op_int_neg <- function(a) if (is.list(a) && !is_ad(a)) lapply(a, op_int_neg) else -a
op_int_div <- function(a, b) {
  if (b == 0L) stop_domain("integer division by zero")
  q <- abs(a) %/% abs(b)
  as.integer(sign(a) * sign(b) * q)
}
op_int_mod <- function(a, b) {
  if (b == 0L) stop_domain("modulus by zero")
  as.integer(a - b * op_int_div(a, b))
}
op_pow <- function(a, b) {
  same_shape(a, b, "power")
  ad_pow(as_dbl(a), as_dbl(b))
}

# matrix products: vectors are stored without dims
as_col <- function(x) ad_reshape(x, c(nelem(x), 1L))
as_row <- function(x) ad_reshape(x, c(1L, nelem(x)))
mm_check <- function(ca, rb, what) {
  if (ca != rb) stop_shape(sprintf("%s: columns of the left operand (%d) must match rows of the right operand (%d)", what, ca, rb))
}
op_mat_mat <- function(a, b) {
  mm_check(ncol(val(a)), nrow(val(b)), "matrix multiplication")
  ad_matmul(as_dbl(a), as_dbl(b))
}
op_mat_vec <- function(a, b) {
  mm_check(ncol(val(a)), nelem(b), "matrix times vector")
  ad_reshape(ad_matmul(as_dbl(a), as_col(as_dbl(b))), NULL)
}
op_rvec_mat <- function(a, b) {
  mm_check(nelem(a), nrow(val(b)), "row_vector times matrix")
  ad_reshape(ad_matmul(as_row(as_dbl(a)), as_dbl(b)), NULL)
}
op_dot <- function(a, b) {
  mm_check(nelem(a), nelem(b), "row_vector times vector")
  ia <- inherits(a, "bcad"); ib <- inherits(b, "bcad")
  va <- if (ia) a$v else a; vb <- if (ib) b$v else b
  y <- sum(va * vb)
  if (!ia && !ib) return(y)
  # one node with both partials, in the shapes of the operands
  ad_node(y, list(a, b), function(g) list(
    if (ia) { r <- g * as.vector(vb); dim(r) <- dim(va); r },
    if (ib) { r <- g * as.vector(va); dim(r) <- dim(vb); r }))
}

# arithmetic on two scalars of static type int or real (not both int):
# no shape checks or S3 dispatch
op_sadd <- function(a, b) {
  ia <- inherits(a, "bcad"); ib <- inherits(b, "bcad")
  if (!ia && !ib) return(as.double(a) + b)
  va <- if (ia) a$v else a; vb <- if (ib) b$v else b
  ad_node(va + vb, list(a, b), function(g) list(if (ia) g, if (ib) g))
}
op_ssub <- function(a, b) {
  ia <- inherits(a, "bcad"); ib <- inherits(b, "bcad")
  if (!ia && !ib) return(as.double(a) - b)
  va <- if (ia) a$v else a; vb <- if (ib) b$v else b
  ad_node(va - vb, list(a, b), function(g) list(if (ia) g, if (ib) -g))
}
op_smul <- function(a, b) {
  ia <- inherits(a, "bcad"); ib <- inherits(b, "bcad")
  if (!ia && !ib) return(as.double(a) * b)
  va <- if (ia) a$v else a; vb <- if (ib) b$v else b
  ad_node(va * vb, list(a, b), function(g) list(if (ia) g * vb, if (ib) g * va))
}
op_sdiv <- function(a, b) {
  ia <- inherits(a, "bcad"); ib <- inherits(b, "bcad")
  if (!ia && !ib) return(as.double(a) / b)
  va <- if (ia) a$v else a; vb <- if (ib) b$v else b
  y <- va / vb
  ad_node(y, list(a, b), function(g) list(if (ia) g / vb, if (ib) -g * y / vb))
}
op_outer <- function(a, b) ad_matmul(as_col(as_dbl(a)), as_row(as_dbl(b)))
op_transpose_mat <- function(a) ad_transpose(as_dbl(a))

# comparisons and logic return int
op_cmp <- function(op, a, b) {
  va <- val(a)
  vb <- val(b)
  r <- switch(op, "==" = va == vb, "!=" = va != vb, "<" = va < vb, "<=" = va <= vb,
              ">" = va > vb, ">=" = va >= vb)
  # comparisons with NaN are false, except != which is true
  if (anyNA(r)) return(as.integer(op == "!="))
  as.integer(isTRUE(r))
}
op_not <- function(a) as.integer(val(a) == 0)
truthy <- function(x) {
  v <- val(x)
  if (length(v) != 1L) stop_shape("condition must be a scalar")
  if (is.na(v)) stop_domain("condition is NaN")
  v != 0
}

# ---------------------------------------------------------------------------
# Elementwise unary functions

unary_defs <- list(
  sqrt = function(x) sqrt(x),
  cbrt = function(x) ad_unary(x, function(v) sign(v) * abs(v)^(1 / 3), function(v, y) 1 / (3 * y^2)),
  square = function(x) ad_square(x),
  exp = function(x) exp(x),
  exp2 = function(x) ad_unary(x, function(v) 2^v, function(v, y) y * log(2)),
  expm1 = function(x) expm1(x),
  log = function(x) log(x),
  log2 = function(x) log2(x),
  log10 = function(x) log10(x),
  log1p = function(x) log1p(x),
  log1m = function(x) ad_log1m(x),
  log1p_exp = function(x) ad_log1p_exp(x),
  log1m_exp = function(x) ad_log1m_exp(x),
  log_inv_logit = function(x) ad_log_inv_logit(x),
  log1m_inv_logit = function(x) ad_log1m_inv_logit(x),
  inv_logit = function(x) ad_inv_logit(x),
  logit = function(x) ad_logit(x),
  inv = function(x) 1 / x,
  inv_sqrt = function(x) 1 / sqrt(x),
  inv_square = function(x) 1 / ad_square(x),
  inv_cloglog = function(x) ad_unary(x, function(v) 1 - exp(-exp(v)), function(v, y) exp(v - exp(v))),
  sin = function(x) sin(x), cos = function(x) cos(x), tan = function(x) tan(x),
  asin = function(x) asin(x), acos = function(x) acos(x), atan = function(x) atan(x),
  sinh = function(x) sinh(x), cosh = function(x) cosh(x), tanh = function(x) tanh(x),
  asinh = function(x) asinh(x), acosh = function(x) acosh(x), atanh = function(x) atanh(x),
  lgamma = function(x) lgamma(x),
  tgamma = function(x) gamma(x),
  digamma = function(x) digamma(x),
  trigamma = function(x) trigamma(x),
  erf = function(x) ad_erf(x),
  erfc = function(x) ad_erfc(x),
  inv_erfc = function(x) {
    if (any(val(x) < 0 | val(x) > 2, na.rm = TRUE)) stop_domain("inv_erfc: argument must be in [0, 2]")
    ad_unary(x, function(v) -stats::qnorm(v / 2) / sqrt(2), function(v, y) -sqrt(pi) / 2 * exp(y * y))
  },
  Phi = function(x) ad_Phi(x),
  inv_Phi = function(x) {
    chk_prob("inv_Phi", "Probability variable", x)
    ad_inv_Phi(x)
  },
  Phi_approx = function(x) ad_inv_logit(0.07056 * x * x * x + 1.5976 * x),
  # quantile from the log probability directly, so it does not underflow
  std_normal_log_qf = function(x) ad_unary(x, function(v) stats::qnorm(v, log.p = TRUE),
                                           function(v, y) exp(v - stats::dnorm(y, log = TRUE))),
  std_normal_qf = function(x) ad_inv_Phi(x),
  lambert_w0 = function(x) ad_unary(x, lambert_w_v, function(v, y) ifelse(v == 0, 1, y / (v * (1 + y)))),
  lambert_wm1 = function(x) ad_unary(x, function(v) lambert_w_v(v, branch = -1), function(v, y) y / (v * (1 + y))),
  floor = function(x) floor(val(x)),
  ceil = function(x) ceiling(val(x)),
  round = function(x) round_half_away(val(x)),
  trunc = function(x) trunc(val(x)),
  step = function(x) as.numeric(val(x) >= 0)
)

round_half_away <- function(x) sign(x) * floor(abs(x) + 0.5)

# principal (branch 0) and lower (branch -1) Lambert W by Halley iteration
lambert_w_v <- function(x, branch = 0) {
  vapply(x, function(z) {
    if (is.na(z)) return(NaN)
    if (z < -exp(-1)) return(NaN)
    if (branch == -1 && (z >= 0)) return(if (z == 0) -Inf else NaN)
    if (z == 0) return(0)
    w <- if (branch == 0) {
      if (z < 1) log1p(z) else log(z) - log(log(z) + 1)
    } else {
      if (z > -0.25) log(-z) - log(-log(-z)) else -1 - sqrt(2 * (1 + exp(1) * z))
    }
    for (it in 1:100) {
      ew <- exp(w)
      f <- w * ew - z
      wn <- w - f / (ew * (w + 1) - (w + 2) * f / (2 * w + 2))
      if (!is.finite(wn)) break
      if (abs(wn - w) <= 1e-15 * (1 + abs(wn))) {
        w <- wn
        break
      }
      w <- wn
    }
    w
  }, 0)
}

for (nm in names(unary_defs)) {
  local({
    f <- unary_defs[[nm]]
    def(nm, "any", r_promote1, function(x) f(as_dbl(x)))
  })
}
def("abs", "int", "int", function(x) abs(x))
def("abs", "int[]", "int[]", function(x) abs(x))
def("abs", "any", r_promote1, function(x) abs(as_dbl(x)))
def("int_step", "real", "int", function(x) as.integer(val(x) > 0))
def("is_inf", "real", "int", function(x) as.integer(is.infinite(val(x))))
def("is_nan", "real", "int", function(x) as.integer(is.nan(val(x))))
def("logical_negation", "real", "int", function(x) as.integer(val(x) == 0))
def("to_int", "real", "int", function(x) {
  v <- trunc(val(x))
  if (!is.finite(v) || abs(v) > .Machine$integer.max) stop_domain("to_int: argument is out of the range of a 32 bit integer")
  as.integer(v)
})
for (nd in 2:8) local({
  at <- function(b) paste0(b, "[", strrep(",", nd - 1L), "]")
  def("to_int", at("int"), at("int"), function(x) x)
  def("to_int", at("real"), at("int"), function(x) {
    v <- trunc(val(x))
    if (any(!is.finite(v) | abs(v) > .Machine$integer.max)) stop_domain("to_int: argument is out of the range of a 32 bit integer")
    storage.mode(v) <- "integer"
    v
  })
})
def("to_int", "int[]", "int[]", function(x) x)
def("to_int", "real[]", "int[]", function(x) {
  v <- trunc(val(x))
  if (any(!is.finite(v) | abs(v) > .Machine$integer.max)) stop_domain("to_int: argument is out of the range of a 32 bit integer")
  as.integer(v)
})

# ---------------------------------------------------------------------------
# Elementwise binary functions (vectorized with broadcasting)

bin2 <- function(nm, f) {
  def(nm, c("any", "any"), r_bcast2, function(a, b) {
    same_shape(a, b, nm)
    f(as_dbl(a), as_dbl(b))
  })
}
bin2("pow", function(a, b) ad_pow(a, b))
# a NaN argument is ignored, as in C's fmin and fmax
bin2("fmin", function(a, b) ad_ifelse(is.na(val(b)) | (!is.na(val(a)) & val(a) <= val(b)), a, b))
bin2("fmax", function(a, b) ad_ifelse(is.na(val(b)) | (!is.na(val(a)) & val(a) >= val(b)), a, b))
bin2("fdim", function(a, b) ad_ifelse(val(a) > val(b), a - b, 0 * a))
bin2("fmod", function(a, b) ad_binary(a, b, function(x, y) x - trunc(x / y) * y,
                                      function(x, y, r) 1, function(x, y, r) -trunc(x / y)))
bin2("hypot", function(a, b) sqrt(a * a + b * b))
bin2("atan2", function(a, b) ad_binary(a, b, atan2,
                                       function(y, x, r) x / (x^2 + y^2),
                                       function(y, x, r) -y / (x^2 + y^2)))
bin2("lbeta", function(a, b) ad_lbeta(a, b))
bin2("beta", function(a, b) exp(ad_lbeta(a, b)))
bin2("lchoose", function(a, b) ad_lchoose(a, b))
bin2("log_diff_exp", function(a, b) ad_log_diff_exp(a, b))
bin2("log_inv_logit_diff", function(a, b) ad_log_diff_exp(ad_log_inv_logit(a), ad_log_inv_logit(b)))
# x * log(y) with 0 * log(0) = 0
ad_lmultiply <- function(a, b) ad_binary(a, b, function(x, y) ifelse(x == 0 & y == 0, 0, x * log(y)),
                                         function(x, y, r) ifelse(x == 0 & y == 0, 0, log(y)),
                                         function(x, y, r) ifelse(x == 0 & y == 0, 0, x / y))
bin2("lmultiply", ad_lmultiply)
bin2("falling_factorial", function(a, b) exp(lgamma(a + 1) - lgamma(a - b + 1)))
bin2("log_falling_factorial", function(a, b) lgamma(a + 1) - lgamma(a - b + 1))
bin2("rising_factorial", function(a, b) exp(lgamma(a + b) - lgamma(a)))
bin2("log_rising_factorial", function(a, b) lgamma(a + b) - lgamma(a))
bin2("owens_t", function(h, a) owens_t_ad(h, a))
bin2("gamma_p", function(a, z) ad_binary(a, z, function(a, z) stats::pgamma(z, a),
                                         function(a, z, y) num_partial(function(aa) stats::pgamma(z, aa), a),
                                         function(a, z, y) stats::dgamma(z, a)))
bin2("gamma_q", function(a, z) ad_binary(a, z, function(a, z) stats::pgamma(z, a, lower.tail = FALSE),
                                         function(a, z, y) num_partial(function(aa) stats::pgamma(z, aa, lower.tail = FALSE), a),
                                         function(a, z, y) -stats::dgamma(z, a)))
bin2("binary_log_loss", function(y, p) ad_ifelse(val(y) != 0, -log(p), -ad_log1m(p)))
bin2("log_modified_bessel_first_kind", function(v, z) ad_binary(v, z,
  function(v, z) log(besselI(z, v, expon.scaled = TRUE)) + z,
  function(v, z, y) num_partial(function(vv) log(besselI(z, vv, expon.scaled = TRUE)), v),
  function(v, z, y) besselI(z, v + 1, expon.scaled = TRUE) / besselI(z, v, expon.scaled = TRUE) + v / z))
bin2("log_sum_exp", ad_log_sum_exp2)
# integer first argument, vectorized: f(n, x) with n data
bin_int <- function(nm, f) {
  def(nm, c("intany", "any"), r_bcast2, function(n, x) {
    same_shape(n, x, nm)
    x <- as_dbl(x)
    vn <- val(n)
    # spread a scalar x over an integer array
    if (length(val(x)) == 1L && length(vn) > 1L) x <- x + array(0, if (is.null(dim(vn))) length(vn) else dim(vn))
    f(rep_len(as.integer(vn), length(val(x))), x)
  })
}
# J and I of negative arguments: J_v(-z) = (-1)^v J_v(z), likewise I_v
odd_sign <- function(z, v) ifelse(z < 0 & v %% 2L == 1L, -1, 1)
bessel_j <- function(z, v) odd_sign(z, v) * besselJ(abs(z), v)
bessel_i <- function(z, v) odd_sign(z, v) * besselI(abs(z), v)
bin_int("bessel_first_kind", function(v, z) ad_unary(z, function(z) bessel_j(z, v), function(z, y) ifelse(z == 0, ifelse(v == 1L, 0.5, 0), bessel_j(z, v - 1L) - v / z * y)))
bin_int("bessel_second_kind", function(v, z) ad_unary(z, function(z) besselY(z, v), function(z, y) besselY(z, v - 1) - v / z * y))
bin_int("modified_bessel_first_kind", function(v, z) ad_unary(z, function(z) bessel_i(z, v), function(z, y) ifelse(z == 0, ifelse(v == 1L, 0.5, 0), bessel_i(z, v + 1L) + v / z * y)))
bin_int("modified_bessel_second_kind", function(v, z) ad_unary(z, function(z) besselK(z, v), function(z, y) -besselK(z, v - 1) - v / z * y))
bin_int("lmgamma", function(k, x) {
  out <- k * (k - 1) / 4 * log(pi) + 0 * x
  for (j in seq_len(max(c(k, 0L)))) {
    # terms past an element's own k use lgamma(1) = 0
    out <- out + lgamma(ad_ifelse(j <= k, x + (1 - j) / 2, 1 + 0 * x))
  }
  out
})
def("ldexp", c("any", "intany"), r_bcast2, function(x, n) {
  same_shape(x, n, "ldexp")
  as_dbl(x) * 2^as.numeric(val(n))
})
def("choose", c("intany", "intany"), function(types) { r <- r_bcast2(types); if (is.null(r)) NULL else ty("int", r$nd) },
    function(n, k) {
  same_shape(n, k, "choose")
  n <- val(n); k <- val(k)
  if (any(n < 0L | k < 0L)) stop_domain("choose: arguments must be nonnegative")
  r <- ifelse(k > n, 0, round(choose(n, k)))
  if (any(r > .Machine$integer.max)) stop_domain("choose: result overflows a 32 bit integer")
  d <- if (length(n) >= length(k)) dim(n) else dim(k)
  out <- as.integer(r)
  if (!is.null(d)) dim(out) <- d
  out
})
def("inc_beta", c("real", "real", "real"), "real", function(a, b, x) {
  inc_beta_ad(as_dbl(a), as_dbl(b), as_dbl(x))
})
def("inv_inc_beta", c("real", "real", "real"), "real", function(a, b, p) {
  a <- as_dbl(a); b <- as_dbl(b); p <- as_dbl(p)
  x <- stats::qbeta(val(p), val(a), val(b))
  if (!any_ad(a, b, p)) return(x)
  # implicit function theorem: I(x; a, b) = p
  dens <- stats::dbeta(x, val(a), val(b))
  ad_node(x, list(a, b, p), function(g) list(
    if (is_ad(a)) -g * num_partial(function(aa) stats::pbeta(x, aa, val(b)), val(a)) / dens,
    if (is_ad(b)) -g * num_partial(function(bb) stats::pbeta(x, val(a), bb), val(b)) / dens,
    if (is_ad(p)) g / dens))
})
def("fma", c("any", "any", "any"), r_bcast_n, function(x, y, z) {
  same_shape(x, y, "fma"); same_shape(x, z, "fma"); same_shape(y, z, "fma")
  as_dbl(x) * as_dbl(y) + as_dbl(z)
})
# log(sum(theta * exp(lambda))) with analytic partials (as Stan Math), which
# stay finite where a weight is 0 or 1
ad_log_mix <- function(th, lp) {
  tv <- val(th); lv <- val(lp)
  ft <- as.vector(tv); fl <- as.vector(lv)
  m <- max(fl)
  v <- if (!is.finite(m)) m else m + log(sum(ft * exp(fl - m)))
  ad_node(v, list(th, lp), function(g) {
    e <- if (is.finite(v)) exp(fl - v) else numeric(length(fl))
    gt <- g * e; gl <- g * ft * e
    if (!is.null(dim(tv))) dim(gt) <- dim(tv)
    if (!is.null(dim(lv))) dim(gl) <- dim(lv)
    list(gt, gl)
  })
}
def("log_mix", c("real", "real", "real"), "real", function(th, a, b) {
  th <- as_dbl(th)
  chk_prob("log_mix", "theta", th)
  a <- as_dbl(a); b <- as_dbl(b)
  t <- val(th); va <- val(a); vb <- val(b)
  m <- max(va, vb)
  v <- if (!is.finite(m)) m else m + log(t * exp(va - m) + (1 - t) * exp(vb - m))
  ad_node(v, list(th, a, b), function(g) {
    ea <- if (is.finite(v)) exp(va - v) else 0
    eb <- if (is.finite(v)) exp(vb - v) else 0
    list(g * (ea - eb), g * t * ea, g * (1 - t) * eb)
  })
})
def("log_mix", c("vecs1", "vecs1"), "real", function(th, lp) {
  if (nelem(th) != nelem(lp)) stop_shape("log_mix: theta and lambda must have the same size")
  ad_log_mix(as_dbl(th), as_dbl(lp))
})
# an array of densities: the sum of the mixtures of its elements
for (lt in c("vector[]", "row_vector[]")) def("log_mix", c("vecs1", lt), "real", function(th, lp) {
  K <- nelem(th)
  vl <- val(lp)
  if (ncol(vl) != K) stop_shape("log_mix: theta and each lambda must have the same size")
  th <- as_dbl(th)
  lp <- as_dbl(lp)
  out <- 0
  for (n in seq_len(nrow(vl))) out <- out + ad_log_mix(th, ad_gather(lp, n + (seq_len(K) - 1L) * nrow(vl)))
  out
})
def("hypergeometric_2F1", c("real", "real", "real", "real"), "real", function(a, b, c, z) {
  # the series converges for |z| < 1; for z < -1/2 the Pfaff transformation
  # 2F1(a, b; c; z) = (1 - z)^(-a) 2F1(a, c - b; c; z / (z - 1)) maps z into
  # (0, 1/2] and extends the function to all z < 1
  f <- function(a, b, c, z) {
    if (!is.finite(z) || z >= 1) stop_domain(sprintf("hypergeometric_2F1: z is %s, but must be less than 1", format(z)))
    if (z < -0.5) return((1 - z)^(-a) * hyper_pfq(c(a, c - b), c, z / (z - 1)))
    hyper_pfq(c(a, b), c, z)
  }
  num_ad4(f, as_dbl(a), as_dbl(b), as_dbl(c), as_dbl(z))
})
for (ta in c("vector", "row_vector", "real[]")) for (tb in c("vector", "row_vector", "real[]")) {
  def("hypergeometric_3F2", c(ta, tb, "real"), "real", function(a, b, z) hyper_pfq_ad(flat(as_dbl(a)), flat(as_dbl(b)), as_dbl(z)))
  def("hypergeometric_pFq", c(ta, tb, "real"), "real", function(a, b, z) hyper_pfq_ad(flat(as_dbl(a)), flat(as_dbl(b)), as_dbl(z)))
}
def("hypergeometric_1F0", c("real", "real"), "real", function(a, z) {
  # 1F0(a; ; z) = (1 - z)^(-a)
  (1 - as_dbl(z))^(-as_dbl(a))
})

# int versions of min/max and two argument real versions
def("min", c("int", "int"), "int", function(a, b) min(a, b))
def("max", c("int", "int"), "int", function(a, b) max(a, b))
def("min", c("real", "real"), "real", function(a, b) ad_ifelse(val(a) <= val(b), as_dbl(a), as_dbl(b)))
def("max", c("real", "real"), "real", function(a, b) ad_ifelse(val(a) >= val(b), as_dbl(a), as_dbl(b)))

# Owen's T function by numerical quadrature, with analytic partials
owens_t_v <- function(h, a) {
  mapply(function(h, a) {
    if (is.na(h) || is.na(a)) return(NaN)
    f <- function(x) exp(-0.5 * h^2 * (1 + x^2)) / (1 + x^2)
    stats::integrate(f, 0, a, rel.tol = 1e-12, abs.tol = 0)$value / (2 * pi)
  }, h, a)
}
owens_t_ad <- function(h, a) {
  ad_binary(h, a, owens_t_v,
    function(h, a, y) -exp(-0.5 * h^2) * (stats::pnorm(a * h) - 0.5) / sqrt(2 * pi) *
      ifelse(h == 0, 0, 1) - ifelse(h == 0, 0, 0),
    function(h, a, y) exp(-0.5 * h^2 * (1 + a^2)) / (2 * pi * (1 + a^2)))
}

inc_beta_ad <- function(a, b, x) {
  va <- val(a); vb <- val(b); vx <- val(x)
  y <- stats::pbeta(vx, va, vb)
  if (!any_ad(a, b, x)) return(y)
  ad_node(y, list(a, b, x), function(g) list(
    if (is_ad(a)) g * num_partial(function(aa) stats::pbeta(vx, aa, vb), va),
    if (is_ad(b)) g * num_partial(function(bb) stats::pbeta(vx, va, bb), vb),
    if (is_ad(x)) g * stats::dbeta(vx, va, vb)))
}

# generalized hypergeometric series pFq(a; b; z), converging for |z| < 1
hyper_pfq <- function(a, b, z) {
  term <- 1
  s <- 1
  for (n in 0:100000) {
    term <- term * prod(a + n) / prod(b + n) * z / (n + 1)
    s <- s + term
    if (abs(term) <= 1e-15 * abs(s) || term == 0) return(s)
    if (!is.finite(s)) break
  }
  # a partial sum is not an answer: refuse rather than return it
  stop_domain(sprintf("hypergeometric series did not converge for z = %.10g", z))
}
hyper_pfq_ad <- function(a, b, z) {
  va <- val(a); vb <- val(b); vz <- val(z)
  y <- hyper_pfq(va, vb, vz)
  if (!any_ad(a, b, z)) return(y)
  ad_node(y, list(a, b, z), function(g) list(
    if (is_ad(a)) g * vapply(seq_along(va), function(i) num_partial(function(t) { aa <- va; aa[i] <- t; hyper_pfq(aa, vb, vz) }, va[i]), 0),
    if (is_ad(b)) g * vapply(seq_along(vb), function(i) num_partial(function(t) { bb <- vb; bb[i] <- t; hyper_pfq(va, bb, vz) }, vb[i]), 0),
    if (is_ad(z)) g * prod(va) / prod(vb) * hyper_pfq(va + 1, vb + 1, vz)))
}
num_ad4 <- function(f, a, b, c, z) {
  v <- lapply(list(a, b, c, z), val)
  y <- do.call(f, v)
  if (!any_ad(a, b, c, z)) return(y)
  args <- list(a, b, c, z)
  ad_node(y, args, function(g) lapply(1:4, function(i) {
    if (!is_ad(args[[i]])) return(NULL)
    g * num_partial(function(t) { vv <- v; vv[[i]] <- t; do.call(f, vv) }, v[[i]])
  }))
}

# ---------------------------------------------------------------------------
# Constants

def("pi", character(0), "real", function() pi)
def("e", character(0), "real", function() exp(1))
def("sqrt2", character(0), "real", function() sqrt(2))
def("log2", character(0), "real", function() log(2))
def("log10", character(0), "real", function() log(10))
def("not_a_number", character(0), "real", function() NaN)
def("positive_infinity", character(0), "real", function() Inf)
def("negative_infinity", character(0), "real", function() -Inf)
def("machine_precision", character(0), "real", function() .Machine$double.eps)
