# Complex numbers.
#
# A complex value is stored as a real array with one extra trailing dimension
# of size 2: the real parts come first, then the imaginary parts. A complex
# scalar is a length 2 vector, a complex_vector[n] an n x 2 array, a
# complex_matrix[m, n] an m x n x 2 array and array[k] complex a k x 2
# array. Indexing, assignment and declarations therefore work unchanged, and
# derivatives flow through the ordinary real operations on the two parts.
# Whether a value is complex is decided statically; the code generator
# passes the information the runtime needs.

# ---- storage -------------------------------------------------------------------

# element dimensions (without the trailing real/imaginary dimension)
cx_edims <- function(x) {
  d <- dim(val(x))
  if (is.null(d)) integer(0) else d[-length(d)]
}
cx_n <- function(x) length(val(x)) %/% 2L
part_shape <- function(ed) if (length(ed) > 1L) ed else NULL
cx_re <- function(x) ad_gather(x, seq_len(cx_n(x)), part_shape(cx_edims(x)))
cx_im <- function(x) {
  n <- cx_n(x)
  ad_gather(x, n + seq_len(n), part_shape(cx_edims(x)))
}

# combine parts; ed gives the element dimensions (integer(0) for a scalar)
cx_make <- function(re, im, ed) {
  y <- ad_c(list(re, im))
  sh <- c(ed, 2L)
  if (length(sh) > 1L) ad_reshape(y, sh) else y
}
part_dims <- function(x) {
  v <- val(x)
  if (is.null(dim(v))) length(v) else dim(v)
}
# combine parts of a container result (vector parts have no dim attribute)
cx_make_c <- function(re, im, scalar = FALSE) cx_make(re, im, if (scalar) integer(0) else part_dims(re))

zeros_like <- function(x) {
  v <- val(x)
  z <- numeric(length(v))
  if (!is.null(dim(v))) dim(z) <- dim(v)
  z
}

# real (or int) value promoted to complex
cx_promote <- function(x, scalar) {
  x <- as_dbl(x)
  cx_make(x, zeros_like(x), if (scalar) integer(0) else part_dims(x))
}

# R complex value (dims are the element dimensions when there are two or more)
cx_value <- function(x) {
  v <- as.vector(val(x))
  n <- length(v) %/% 2L
  z <- complex(real = v[seq_len(n)], imaginary = v[n + seq_len(n)])
  ed <- cx_edims(x)
  if (length(ed) > 1L) dim(z) <- ed
  z
}
cx_from_z <- function(z, ed) {
  y <- c(Re(as.vector(z)), Im(as.vector(z)))
  sh <- c(ed, 2L)
  if (length(sh) > 1L) dim(y) <- sh
  y
}
z_dims <- function(z) if (is.null(dim(z))) length(z) else dim(z)

# value written to the output: an R complex scalar, vector or array
cx_out <- function(x) cx_value(x)

# output names in cmdstanr's form: z[real], z[imag], v[1,real], ...
cx_flat_names <- function(name, dims) {
  ed <- dims[-length(dims)]
  base <- if (length(ed)) {
    if (prod(ed) == 0L) return(character(0))
    g <- do.call(expand.grid, lapply(ed, seq_len))
    paste0(name, "[", do.call(paste, c(unname(as.list(g)), sep = ",")), ",")
  } else paste0(name, "[")
  as.vector(rbind(paste0(base, "real]"), paste0(base, "imag]")))
}

# data supplied as R complex (or real) values
read_cx_data <- function(x, info) {
  if (isTRUE(attr(x, "stan_json"))) x <- json_complex(x)
  if (is.list(x)) x <- simplify_data_list(x)
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x))
    stop_data(sprintf("variable '%s' must be complex or numeric", info$name))
  sub <- utils::modifyList(info, list(base = "real", stan_type = "real", eshape = info$eshape[-length(info$eshape)]))
  part <- function(p) {
    d <- list(p)
    names(d) <- info$name
    read_data_value(d, sub)
  }
  re <- part(if (is.complex(x)) array(Re(x), dim = dim(x) %||% length(x)) else x)
  im <- part(if (is.complex(x)) array(Im(x), dim = dim(x) %||% length(x)) else zeros_like(x))
  cx_make(re, im, c(sub$adims, sub$eshape))
}

# ---- printing ------------------------------------------------------------------

cx_fmt <- function(x, nd) {
  z <- cx_value(x)
  one <- function(w) paste0("(", format_scalar(Re(w)), ",", format_scalar(Im(w)), ")")
  ed <- cx_edims(x)
  if (!length(ed)) return(one(z))
  if (length(ed) == 2L) {
    dim(z) <- ed
    rows <- vapply(seq_len(ed[1L]), function(i) paste0("[", paste(vapply(z[i, ], one, ""), collapse = ","), "]"), "")
    return(paste0("[", paste(rows, collapse = ","), "]"))
  }
  paste0("[", paste(vapply(as.vector(z), one, ""), collapse = ","), "]")
}

# ---- operators -----------------------------------------------------------------

scalar_base <- function(b) b %in% c("int", "real", "complex")

# elementwise complex product and quotient on parts (scalars broadcast)
cx_mul_parts <- function(ar, ai, br, bi) list(ar * br - ai * bi, ar * bi + ai * br)
cx_div_parts <- function(ar, ai, br, bi) {
  den <- br * br + bi * bi
  list((ar * br + ai * bi) / den, (ai * br - ar * bi) / den)
}

# a op b where at least one operand is complex; ba, bb are the static bases
cx_binop <- function(op, a, b, ba, bb) {
  sa <- scalar_base(ba)
  sb <- scalar_base(bb)
  if (!is_cx_base(ba)) a <- cx_promote(a, sa)
  if (!is_cx_base(bb)) b <- cx_promote(b, sb)
  ar <- cx_re(a); ai <- cx_im(a); br <- cx_re(b); bi <- cx_im(b)
  ed <- if (sa) cx_edims(b) else cx_edims(a)
  ew_check <- function(what) if (!sa && !sb) same_shape(ar, br, what)
  rb <- function(x) if (is_cx_base(x)) cx_to_real[[x]] else x
  switch(op,
    "==" = , "!=" = {
      eq <- isTRUE(val(ar) == val(br)) && isTRUE(val(ai) == val(bi))
      as.integer(if (op == "==") eq else !eq)
    },
    "+" = { ew_check("addition"); cx_make(ar + br, ai + bi, ed) },
    "-" = { ew_check("subtraction"); cx_make(ar - br, ai - bi, ed) },
    ".*" = { ew_check("elementwise multiplication"); p <- cx_mul_parts(ar, ai, br, bi); cx_make(p[[1L]], p[[2L]], ed) },
    "./" = { ew_check("elementwise division"); p <- cx_div_parts(ar, ai, br, bi); cx_make(p[[1L]], p[[2L]], ed) },
    "^" = , ".^" = { ew_check("power"); cx_pow_ew(a, b) },
    "*" = {
      if (sa || sb) { p <- cx_mul_parts(ar, ai, br, bi); return(cx_make(p[[1L]], p[[2L]], ed)) }
      key <- paste(rb(ba), rb(bb))
      f <- switch(key, "matrix matrix" = op_mat_mat, "matrix vector" = op_mat_vec,
                  "row_vector matrix" = op_rvec_mat, "row_vector vector" = op_dot, "vector row_vector" = op_outer)
      # row_vector * vector is Stan's dot_product, which conjugates the left operand
      if (key == "row_vector vector") return(cx_make(f(ar, br) + f(ai, bi), f(ar, bi) - f(ai, br), integer(0)))
      cx_make_c(f(ar, br) - f(ai, bi), f(ar, bi) + f(ai, br))
    },
    "/" = {
      if (sb || sa) { p <- cx_div_parts(ar, ai, br, bi); return(cx_make(p[[1L]], p[[2L]], ed)) }
      cx_mdivide_right(a, b, rb(ba) == "row_vector")
    },
    stop("internal error: complex operator ", op))
}

# [a, b, ...] with complex elements (items are already complex)
cx_rowvec <- function(items, rows) {
  cx_make_c(rt_rowvec(lapply(items, cx_re), rows), rt_rowvec(lapply(items, cx_im), rows))
}

cx_transpose <- function(x) {
  ar <- cx_re(x); ai <- cx_im(x)
  cx_make_c(ad_transpose(ar), ad_transpose(ai))
}

# ---- elementwise functions -----------------------------------------------------

# holomorphic f with derivative df(z, w); x is complex storage
cx_holo <- function(x, f, df) {
  z <- cx_value(x)
  w <- f(as.vector(z))
  vx <- val(x)
  y <- c(Re(w), Im(w))
  if (!is.null(dim(vx))) dim(y) <- dim(vx)
  if (!inherits(x, "bcad")) return(y)
  d <- df(as.vector(z), w)
  a <- Re(d); b <- Im(d)
  n <- length(w)
  ad_node(y, list(x), function(g) {
    gu <- g[seq_len(n)]; gv <- g[n + seq_len(n)]
    r <- c(gu * a + gv * b, gv * a - gu * b)
    if (!is.null(dim(vx))) dim(r) <- dim(vx)
    list(r)
  })
}

# real valued function of the parts: value f(re, im), partials dre, dim
cx_realfun <- function(x, f, dre, dim_) {
  v <- as.vector(val(x))
  n <- length(v) %/% 2L
  re <- v[seq_len(n)]; im <- v[n + seq_len(n)]
  y <- f(re, im)
  ed <- cx_edims(x)
  if (length(ed) > 1L) dim(y) <- ed
  if (!inherits(x, "bcad")) return(y)
  vx <- val(x)
  gre <- dre(re, im, y); gim <- dim_(re, im, y)
  ad_node(y, list(x), function(g) {
    g <- as.vector(g)
    r <- c(g * gre, g * gim)
    if (!is.null(dim(vx))) dim(r) <- dim(vx)
    list(r)
  })
}

cx_abs <- function(x) cx_realfun(x, function(a, b) Mod(complex(real = a, imaginary = b)),
                                 function(a, b, y) ifelse(y == 0, 0, a / y), function(a, b, y) ifelse(y == 0, 0, b / y))
cx_arg <- function(x) cx_realfun(x, function(a, b) atan2(b, a),
                                 function(a, b, y) -b / (a * a + b * b), function(a, b, y) a / (a * a + b * b))
cx_norm <- function(x) cx_realfun(x, function(a, b) a * a + b * b, function(a, b, y) 2 * a, function(a, b, y) 2 * b)
cx_conj <- function(x) cx_make(cx_re(x), -cx_im(x), cx_edims(x))

cx_exp <- function(x) cx_holo(x, exp, function(z, w) w)
cx_log <- function(x) cx_holo(x, log, function(z, w) 1 / z)
cx_mul_ew <- function(a, b) {
  p <- cx_mul_parts(cx_re(a), cx_im(a), cx_re(b), cx_im(b))
  cx_make(p[[1L]], p[[2L]], if (length(cx_edims(a))) cx_edims(a) else cx_edims(b))
}
# a^b elementwise with scalar broadcasting, as exp(b log a); 0^b is 0
cx_pow_ew <- function(a, b) {
  out <- cx_exp(cx_mul_ew(b, cx_log(a)))
  za <- cx_value(a)
  if (any(za == 0)) {
    # ponytail: the zero base is patched in value only (no derivative through it)
    zv <- cx_value(out)
    zb <- cx_value(b)
    zero <- rep_len(as.vector(za) == 0, length(zv)) & Re(rep_len(as.vector(zb), length(zv))) > 0
    zv[zero] <- 0
    out <- cx_from_z(zv, cx_edims(out))
  }
  out
}

holo_funs <- list(
  exp = list(exp, function(z, w) w),
  log = list(log, function(z, w) 1 / z),
  log10 = list(function(z) log(z) / log(10), function(z, w) 1 / (z * log(10))),
  sqrt = list(sqrt, function(z, w) 1 / (2 * w)),
  sin = list(sin, function(z, w) cos(z)),
  cos = list(cos, function(z, w) -sin(z)),
  tan = list(tan, function(z, w) 1 / cos(z)^2),
  asin = list(asin, function(z, w) 1 / sqrt(1 - z * z)),
  acos = list(acos, function(z, w) -1 / sqrt(1 - z * z)),
  atan = list(atan, function(z, w) 1 / (1 + z * z)),
  sinh = list(sinh, function(z, w) cosh(z)),
  cosh = list(cosh, function(z, w) sinh(z)),
  tanh = list(tanh, function(z, w) 1 / cosh(z)^2),
  asinh = list(asinh, function(z, w) 1 / sqrt(z * z + 1)),
  acosh = list(acosh, function(z, w) 1 / (sqrt(z - 1) * sqrt(z + 1))),
  atanh = list(atanh, function(z, w) 1 / (1 - z * z)))
for (nm in names(holo_funs)) local({
  f <- holo_funs[[nm]]
  def(nm, "complex", "complex", function(x) cx_holo(x, f[[1L]], f[[2L]]))
})
def("proj", "complex", "complex", function(x) {
  z <- cx_value(x)
  if (is.infinite(Re(z)) || is.infinite(Im(z))) return(c(Inf, if (Im(z) < 0 || (Im(z) == 0 && 1 / Im(z) < 0)) -0 else 0))
  x
})
def("polar", c("real", "real"), "complex", function(r, theta) {
  r <- as_dbl(r); theta <- as_dbl(theta)
  if (val(r) < 0 || is.nan(val(r))) stop_domain(sprintf("polar: magnitude is %s, but must be nonnegative", fmt_num(val(r))))
  cx_make(r * cos(theta), r * sin(theta), integer(0))
})
def("arg", "complex", "real", cx_arg)
def("norm", "complex", "real", cx_norm)

# all complex types (any array dimensions)
r_realpart <- function(types) ty_realpart(types[[1L]])
def("abs", "cxall", r_realpart, cx_abs)
def("get_real", "cxall", r_realpart, function(x) cx_re(x))
def("get_imag", "cxall", r_realpart, function(x) cx_im(x))
def("conj", "cxall", r_same1, cx_conj)

# to_complex(re, im), vectorized with scalar broadcasting
def("to_complex", character(0), "complex", function() c(0, 0))
def("to_complex", "real", "complex", function(re) cx_make(as_dbl(re), 0, integer(0)))
def("to_complex", c("real", "real"), "complex", function(re, im) cx_make(as_dbl(re), as_dbl(im), integer(0)))
def("to_complex", c("any", "any"), function(types) {
  rt <- r_bcast2(types)
  if (is.null(rt)) NULL else ty_cx(rt)
}, function(re, im) {
  re <- as_dbl(re); im <- as_dbl(im)
  lr <- length(val(re)); li <- length(val(im))
  # one argument may be a scalar; the other decides the shape
  if (lr == 1L && li != 1L) re <- ad_gather(re, rep(1L, li), dim(val(im)))
  else if (li == 1L && lr != 1L) im <- ad_gather(im, rep(1L, lr), dim(val(re)))
  else same_shape(re, im, "to_complex")
  cx_make_c(re, im)
})

# pow with complex arguments (scalars broadcast)
def("pow", c("cxnum", "cxnum"), function(types) {
  if (!any(vapply(types, is_cx_ty, TRUE))) return(NULL)
  rt <- r_bcast2(lapply(types, function(t) { r <- ty_realpart(t); if (r$base == "int") ty("real", r$nd) else r }))
  if (is.null(rt)) NULL else ty_cx(rt)
}, function(a, b) {
  la <- length(cx_edims(a)); lb <- length(cx_edims(b))
  if (la && lb) same_shape(cx_re(a), cx_re(b), "pow")
  cx_pow_ew(a, b)
}, cxargs = c(TRUE, TRUE))

# ---- reductions and products ----------------------------------------------------

cx_prod <- function(x) {
  re <- flat(cx_re(x)); im <- flat(cx_im(x))
  n <- nelem(re)
  pr <- 1; pi <- 0
  for (k in seq_len(n)) {
    a <- ad_select(re, k); b <- ad_select(im, k)
    p <- cx_mul_parts(pr, pi, a, b)
    pr <- p[[1L]]; pi <- p[[2L]]
  }
  cx_make(as_dbl(pr), as_dbl(pi), integer(0))
}
for (b in c("complex_vector", "complex_row_vector", "complex_matrix", "complex[]")) def("prod", b, "complex", cx_prod)

# ---- matrix division -------------------------------------------------------------

# real embedding [Re -Im; Im Re] of a complex matrix
cx_embed <- function(m) {
  ar <- cx_re(m); ai <- cx_im(m)
  ac <- lib_impl("append_col", list(T_mat, T_mat))
  ar_ <- lib_impl("append_row", list(T_mat, T_mat))
  ar_(ac(ar, -ai), ac(ai, ar))
}
# B / A = B A^-1, solved as A^T X^T = B^T
cx_mdivide_right <- function(b, a, rowvec) {
  d <- cx_edims(a)
  if (d[1L] != d[2L]) stop_shape(sprintf("mdivide_right: matrix must be square, but has %d rows and %d columns", d[1L], d[2L]))
  at <- cx_transpose(a)
  bt <- if (rowvec) b else cx_transpose(b)
  rhs <- if (rowvec) ad_c(list(cx_re(bt), cx_im(bt)))
         else lib_impl("append_row", list(T_mat, T_mat))(cx_re(bt), cx_im(bt))
  if (!rowvec && nrow(val(cx_re(bt))) != d[1L] || rowvec && nelem(cx_re(bt)) != d[1L])
    stop_shape("mdivide_right: size mismatch between the operands")
  sol <- if (rowvec) lib_impl("mdivide_left", list(T_mat, T_vec))(cx_embed(at), rhs)
         else lib_impl("mdivide_left", list(T_mat, T_mat))(cx_embed(at), rhs)
  n <- d[1L]
  if (rowvec) return(cx_make(ad_gather(sol, seq_len(n)), ad_gather(sol, n + seq_len(n)), n))
  k <- ncol(val(sol))
  idx_r <- as.vector(outer(seq_len(n), (seq_len(k) - 1L) * 2L * n, `+`))
  xr <- ad_gather(sol, idx_r, c(n, k)); xi <- ad_gather(sol, idx_r + n, c(n, k))
  cx_transpose(cx_make(xr, xi, c(n, k)))
}
def("mdivide_right", c("complex_row_vector", "complex_matrix"), "complex_row_vector",
    function(b, a) cx_mdivide_right(b, a, TRUE), cxargs = c(TRUE, TRUE))
def("mdivide_right", c("complex_matrix", "complex_matrix"), "complex_matrix",
    function(b, a) cx_mdivide_right(b, a, FALSE), cxargs = c(TRUE, TRUE))

# ---- FFT ---------------------------------------------------------------------------

# complex linear map with forward f and adjoint fa on R complex values
cx_linear <- function(x, f, fa) {
  z <- cx_value(x)
  ed <- cx_edims(x)
  w <- f(z)
  y <- cx_from_z(w, ed)
  if (!inherits(x, "bcad")) return(y)
  ad_node(y, list(x), function(g) {
    gz <- complex(real = g[seq_len(length(w))], imaginary = g[length(w) + seq_len(length(w))])
    if (length(ed) > 1L) dim(gz) <- ed
    list(cx_from_z(fa(gz), ed))
  })
}
def("fft", "complex_vector", "complex_vector", function(x)
  cx_linear(x, function(z) stats::fft(z), function(g) stats::fft(g, inverse = TRUE)), cxargs = TRUE)
def("inv_fft", "complex_vector", "complex_vector", function(x) {
  n <- cx_n(x)
  cx_linear(x, function(z) stats::fft(z, inverse = TRUE) / n, function(g) stats::fft(g) / n)
}, cxargs = TRUE)
def("fft2", "complex_matrix", "complex_matrix", function(x)
  cx_linear(x, function(z) stats::fft(z), function(g) stats::fft(g, inverse = TRUE)), cxargs = TRUE)
def("inv_fft2", "complex_matrix", "complex_matrix", function(x) {
  n <- cx_n(x)
  cx_linear(x, function(z) stats::fft(z, inverse = TRUE) / n, function(g) stats::fft(g) / n)
}, cxargs = TRUE)

# ---- decompositions -----------------------------------------------------------------

# value only decompositions get derivatives by central differences of the
# flattened input (they are rarely applied to parameters)
cx_numeric <- function(x, f) {
  vx <- val(x)
  y <- f(vx)
  if (!inherits(x, "bcad")) return(y)
  n <- length(vx)
  J <- vapply(seq_len(n), function(k) {
    h <- 1e-6 * max(abs(vx[k]), 1)
    up <- vx; dn <- vx
    up[k] <- up[k] + h; dn[k] <- dn[k] - h
    (as.vector(f(up)) - as.vector(f(dn))) / (2 * h)
  }, numeric(length(y)))
  J <- matrix(J, nrow = length(y))
  ad_node(y, list(x), function(g) {
    r <- as.vector(crossprod(J, as.vector(g)))
    if (!is.null(dim(vx))) dim(r) <- dim(vx)
    list(r)
  })
}

as_zmat <- function(v, complex_in) {
  if (complex_in) {
    d <- dim(v)
    n <- length(v) %/% 2L
    z <- complex(real = v[seq_len(n)], imaginary = v[n + seq_len(n)])
    dim(z) <- d[1:2]
    return(z)
  }
  v
}
square_check <- function(f, z) {
  if (nrow(z) != ncol(z)) stop_shape(sprintf("%s: matrix must be square, but has %d rows and %d columns", f, nrow(z), ncol(z)))
  if (length(z) == 0L) stop_shape(sprintf("%s: matrix must not be empty", f))
}

# eigenvectors scaled to unit norm, as Eigen returns them
zeig <- function(z, f) {
  square_check(f, z)
  if (any(!is.finite(z))) stop_domain(sprintf("%s: matrix has non finite elements", f))
  e <- eigen(z, symmetric = FALSE)
  # Eigen's complex solver sorts by increasing magnitude (stable for ties)
  if (is.complex(z)) {
    o <- order(Mod(e$values))
    e$values <- e$values[o]
    e$vectors <- e$vectors[, o, drop = FALSE]
  }
  v <- e$vectors
  v <- sweep(v, 2L, sqrt(colSums(Mod(v)^2)), "/")
  list(values = e$values, vectors = v)
}
zeig_sym <- function(z, f) {
  square_check(f, z)
  if (any(z != t(z))) stop_domain(sprintf("%s: m is not symmetric", f))
  # Stan reads the lower triangle as Hermitian
  h <- z
  h[upper.tri(h)] <- Conj(t(z))[upper.tri(h)]
  e <- eigen(h, symmetric = TRUE)
  o <- order(e$values)
  list(values = complex(real = e$values[o], imaginary = 0), vectors = e$vectors[, o, drop = FALSE])
}
zschur <- function(z, f) {
  square_check(f, z)
  n <- nrow(z)
  # complex Schur form from the eigendecomposition: QR of the eigenvectors
  # gives a unitary U with U^H A U upper triangular (distinct eigenvalues)
  e <- eigen(z, symmetric = FALSE)
  # diagonal in increasing magnitude, as Eigen's ComplexSchur usually ends
  q <- qr.Q(qr(e$vectors[, order(Mod(e$values)), drop = FALSE]))
  t_ <- Conj(t(q)) %*% z %*% q
  low <- t_[lower.tri(t_)]
  # ponytail: Schur form from eigenvectors; a defective matrix (repeated
  # eigenvalues without a full eigenbasis) needs a complex QR algorithm
  if (length(low) && max(Mod(low)) > 1e-8 * max(1, max(Mod(z))))
    stop_unsupported(sprintf("%s: matrix is defective (no complete set of eigenvectors)", f))
  t_[lower.tri(t_)] <- 0
  list(t = t_, u = q)
}
zsvd <- function(z, f) {
  s <- svd(z)
  list(d = s$d, u = s$u, v = s$v)
}

dec_def <- function(name, ret, fun) {
  for (inb in c("matrix", "complex_matrix")) local({
    cin <- inb == "complex_matrix"
    def(name, inb, ret, function(m) {
      vx <- val(m)
      res <- function(v) {
        w <- fun(as_zmat(v, cin), name)
        cx_from_z(as.complex(w), z_dims(w))
      }
      cx_numeric(m, res)
    })
  })
}
dec_def("eigenvalues", "complex_vector", function(z, f) zeig(z, f)$values)
dec_def("eigenvectors", "complex_matrix", function(z, f) zeig(z, f)$vectors)
dec_def("complex_schur_decompose_t", "complex_matrix", function(z, f) zschur(z, f)$t)
dec_def("complex_schur_decompose_u", "complex_matrix", function(z, f) zschur(z, f)$u)
for (nm in c("eigenvalues_sym", "eigenvectors_sym")) local({
  which <- if (nm == "eigenvalues_sym") "values" else "vectors"
  def(nm, "complex_matrix", if (which == "values") "complex_vector" else "complex_matrix", function(m) {
    cx_numeric(m, function(v) { w <- zeig_sym(as_zmat(v, TRUE), nm)[[which]]; cx_from_z(w, z_dims(w)) })
  })
})
def("singular_values", "complex_matrix", "vector", function(m) cx_numeric(m, function(v) zsvd(as_zmat(v, TRUE))$d))
def("svd_U", "complex_matrix", "complex_matrix", function(m) cx_numeric(m, function(v) { w <- zsvd(as_zmat(v, TRUE))$u; cx_from_z(w, dim(w)) }))
def("svd_V", "complex_matrix", "complex_matrix", function(m) cx_numeric(m, function(v) { w <- zsvd(as_zmat(v, TRUE))$v; cx_from_z(w, dim(w)) }))

# ---- container functions lifted from their real versions ----------------------------

# functions linear over the reals: apply the real version to both parts
cx_linear_funs <- c("append_row", "append_col", "append_array", "block", "col", "row", "head", "tail",
                    "segment", "sub_col", "sub_row", "cumulative_sum", "diag_matrix", "diagonal", "add_diag",
                    "rep_vector", "rep_row_vector", "rep_matrix", "rep_array", "reverse", "sum", "trace",
                    "to_matrix", "to_vector", "to_row_vector", "to_array_1d", "to_array_2d", "to_vector_array", "to_row_vector_array",
                    "symmetrize_from_lower_tri", "dims", "size", "num_elements", "rows", "cols")
# bilinear products; conj_first marks dot products that conjugate the left operand
cx_bilinear_funs <- c(diag_pre_multiply = FALSE, diag_post_multiply = FALSE, rows_dot_product = FALSE,
                      columns_dot_product = FALSE, dot_product = TRUE)
# self products |x|^2
cx_self_funs <- c("dot_self", "rows_dot_self", "columns_dot_self")

int_pat <- function(p) is.character(p) && p %in% c("int", "ints", "int[]", "int[,]", "intany")

resolve_cx_lift <- function(s, name, types) {
  if (!(name %in% c(cx_linear_funs, names(cx_bilinear_funs), cx_self_funs))) return(NULL)
  # stanc has dims() for complex scalars and arrays, not complex vectors or matrices
  if (name == "dims" && types[[1L]]$nd == 0L && types[[1L]]$base != "complex") return(NULL)
  rtypes <- lapply(types, function(t) {
    r <- ty_realpart(t)
    if (is_cx_ty(t) && r$base == "real") r else r
  })
  r <- resolve_call(s, name, rtypes)
  if (is.null(r) || isTRUE(r$sig$user)) return(NULL)
  f <- r$sig$impl
  ret_int <- r$ret$base == "int"
  ret <- if (ret_int) r$ret else ty_cx(r$ret)
  scalar <- ret$nd == 0L && ret$base == "complex"
  valargs <- !vapply(r$sig$args, int_pat, TRUE)
  impl <- if (name %in% cx_linear_funs) {
    if (ret_int) function(...) { a <- list(...); a[valargs] <- lapply(a[valargs], cx_re); do.call(f, a) }
    else function(...) {
      a <- list(...)
      re <- a; im <- a
      re[valargs] <- lapply(a[valargs], cx_re)
      im[valargs] <- lapply(a[valargs], cx_im)
      cx_make_c(as_dbl(do.call(f, re)), as_dbl(do.call(f, im)), scalar)
    }
  } else if (name %in% cx_self_funs) {
    function(x) {
      rr <- as_dbl(f(cx_re(x))) + as_dbl(f(cx_im(x)))
      cx_make_c(rr, zeros_like(rr), scalar)
    }
  } else {
    cj <- cx_bilinear_funs[[name]]
    function(a, b) {
      ar <- cx_re(a); ai <- cx_im(a); br <- cx_re(b); bi <- cx_im(b)
      if (cj) cx_make_c(as_dbl(f(ar, br)) + as_dbl(f(ai, bi)), as_dbl(f(ar, bi)) - as_dbl(f(ai, br)), scalar)
      else cx_make_c(as_dbl(f(ar, br)) - as_dbl(f(ai, bi)), as_dbl(f(ar, bi)) + as_dbl(f(ai, br)), scalar)
    }
  }
  sig <- list(name = name, args = r$sig$args, ret = ret, impl = impl, cxargs = valargs)
  list(sig = sig, ret = ret)
}
