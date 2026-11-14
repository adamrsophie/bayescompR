# Container functions: reductions, construction, slicing, sorting.
#
# Storage conventions (see ARCHITECTURE.md): vector and row_vector values are
# plain numeric vectors without dims, a matrix has dim c(rows, cols), and an
# array with dimensions adims of an element with shape s is an R array with
# dim c(adims, s).

shape_of <- function(x) {
  v <- val(x)
  if (is.null(dim(v))) length(v) else dim(v)
}

# ---- reductions --------------------------------------------------------------

def("sum", "int[]", "int", function(x) {
  r <- sum(as.numeric(x))
  if (abs(r) > .Machine$integer.max) stop_domain("sum: integer overflow")
  as.integer(r)
})
def("sum", "container", "real", function(x) ad_sum(as_dbl(x)))
def("sum", "real", "real", function(x) as_dbl(x))
def("prod", "int[]", "int", function(x) as.integer(prod(x)))
def("prod", "container", "real", function(x) ad_prod(as_dbl(x)))
def("mean", "container", "real", function(x) {
  n <- nelem(x)
  if (n == 0L) stop_domain("mean: argument has size zero")
  ad_sum(as_dbl(x)) / n
})
def("variance", "container", "real", function(x) {
  n <- nelem(x)
  if (n < 2L) return(if (n == 1L) 0 else stop_domain("variance: argument has size zero"))
  x <- flat(as_dbl(x))
  m <- ad_sum(x) / n
  d <- x - m
  ad_sum(d * d) / (n - 1)
})
def("sd", "container", "real", function(x) {
  n <- nelem(x)
  if (n < 2L) return(if (n == 1L) 0 else stop_domain("sd: argument has size zero"))
  x <- flat(as_dbl(x))
  vx <- val(x)
  d <- vx - mean(vx)
  sdv <- sqrt(sum(d * d) / (n - 1))
  if (!is_ad(x)) return(sdv)
  # partial (x - mean) / (sd (n - 1)), taken as zero at a constant vector
  ad_node(sdv, list(x), function(g) list(if (sdv > 0) g * d / (sdv * (n - 1)) else numeric(n)))
})
def("log_sum_exp", "container", "real", function(x) ad_log_sum_exp(flat(as_dbl(x))))
def("max", "int[]", "int", function(x) {
  if (!length(x)) stop_domain("max: integer array of size zero")
  max(x)
})
def("min", "int[]", "int", function(x) {
  if (!length(x)) stop_domain("min: integer array of size zero")
  min(x)
})
def("max", "container", "real", function(x) {
  if (!nelem(x)) return(-Inf)
  ad_select(flat(as_dbl(x)), which.max(as.vector(val(x))))
})
def("min", "container", "real", function(x) {
  if (!nelem(x)) return(Inf)
  ad_select(flat(as_dbl(x)), which.min(as.vector(val(x))))
})
def("norm1", "vecs1", "real", function(x) ad_sum(abs(as_dbl(x))))
def("norm2", "vecs1", "real", function(x) sqrt(ad_sum(as_dbl(x) * as_dbl(x))))
def("dot_self", "vectors", "real", function(x) ad_sum(as_dbl(x) * as_dbl(x)))
def("dot_product", c("vecs1", "vecs1"), "real", function(a, b) {
  if (nelem(a) != nelem(b)) stop_shape("dot_product: arguments must have the same size")
  ad_sum(flat(as_dbl(a)) * flat(as_dbl(b)))
})
def("squared_distance", c("vectors", "vectors"), "real", function(a, b) {
  if (nelem(a) != nelem(b)) stop_shape("squared_distance: arguments must have the same size")
  d <- flat(as_dbl(a)) - flat(as_dbl(b))
  ad_sum(d * d)
})
def("squared_distance", c("real", "real"), "real", function(a, b) ad_square(as_dbl(a) - as_dbl(b)))
def("distance", c("vectors", "vectors"), "real", function(a, b) {
  if (nelem(a) != nelem(b)) stop_shape("distance: arguments must have the same size")
  d <- flat(as_dbl(a)) - flat(as_dbl(b))
  sqrt(ad_sum(d * d))
})
def("distance", c("real", "real"), "real", function(a, b) abs(as_dbl(a) - as_dbl(b)))
def("columns_dot_product", c("matrix", "matrix"), "row_vector", function(a, b) {
  same_shape(a, b, "columns_dot_product")
  col_sums(as_dbl(a) * as_dbl(b))
})
def("columns_dot_product", c("vector", "vector"), "row_vector", function(a, b) {
  same_shape(a, b, "columns_dot_product")
  ad_c(list(ad_sum(as_dbl(a) * as_dbl(b))))
})
def("columns_dot_product", c("row_vector", "row_vector"), "row_vector", function(a, b) {
  same_shape(a, b, "columns_dot_product")
  as_dbl(a) * as_dbl(b)
})
def("rows_dot_product", c("matrix", "matrix"), "vector", function(a, b) {
  same_shape(a, b, "rows_dot_product")
  row_sums(as_dbl(a) * as_dbl(b))
})
def("rows_dot_product", c("vector", "vector"), "vector", function(a, b) {
  same_shape(a, b, "rows_dot_product")
  as_dbl(a) * as_dbl(b)
})
def("rows_dot_product", c("row_vector", "row_vector"), "vector", function(a, b) {
  same_shape(a, b, "rows_dot_product")
  ad_c(list(ad_sum(as_dbl(a) * as_dbl(b))))
})
def("columns_dot_self", "matrix", "row_vector", function(x) col_sums(as_dbl(x) * as_dbl(x)))
def("columns_dot_self", "vector", "row_vector", function(x) ad_c(list(ad_sum(as_dbl(x) * as_dbl(x)))))
def("columns_dot_self", "row_vector", "row_vector", function(x) as_dbl(x) * as_dbl(x))
def("rows_dot_self", "matrix", "vector", function(x) row_sums(as_dbl(x) * as_dbl(x)))
def("rows_dot_self", "vector", "vector", function(x) as_dbl(x) * as_dbl(x))
def("rows_dot_self", "row_vector", "vector", function(x) ad_c(list(ad_sum(as_dbl(x) * as_dbl(x)))))

col_sums <- function(X) {
  v <- val(X)
  ad_reshape(ad_matmul(matrix(1, 1L, nrow(v)), X), NULL)
}
row_sums <- function(X) {
  v <- val(X)
  ad_reshape(ad_matmul(X, matrix(1, ncol(v), 1L)), NULL)
}

def("cumulative_sum", "int[]", "int[]", function(x) as.integer(cumsum(x)))
def("cumulative_sum", "real[]", "real[]", function(x) cumsum(as_dbl(x)))
def("cumulative_sum", "vector", "vector", function(x) cumsum(as_dbl(x)))
def("cumulative_sum", "row_vector", "row_vector", function(x) cumsum(as_dbl(x)))

# ---- sizes -------------------------------------------------------------------

def("num_elements", "any", "int", function(x) length(val(x)))
def("size", "anyarr", "int", function(x) as.integer(shape_of(x)[1L]))
def("size", "vector", "int", function(x) nelem(x))
def("size", "row_vector", "int", function(x) nelem(x))
def("size", "matrix", "int", function(x) nelem(x))
def("size", "int", "int", function(x) 1L)
def("size", "real", "int", function(x) 1L)
def("rows", "vector", "int", function(x) nelem(x))
def("rows", "row_vector", "int", function(x) 1L)
def("rows", "matrix", "int", function(x) nrow(val(x)))
def("cols", "vector", "int", function(x) 1L)
def("cols", "row_vector", "int", function(x) nelem(x))
def("cols", "matrix", "int", function(x) ncol(val(x)))
def("dims", "any", "int[]", function(x) {
  v <- val(x)
  if (length(v) == 1L && is.null(dim(v))) return(integer(0))
  as.integer(shape_of(x))
})
# note: dims of a scalar is the empty array; dims of a row_vector is c(1, n)
# and of a vector c(n, 1) in Stan. The code generator passes the static type
# through dims_typed so these cases are exact.
dims_typed <- function(x, base, nd) {
  sh <- as.integer(shape_of(x))
  if (base %in% c("int", "real") && nd == 0L) return(integer(0))
  if (base == "vector") return(c(sh[seq_len(nd)], if (nd) sh[nd + 1L] else sh[1L], 1L))
  if (base == "row_vector") return(c(sh[seq_len(nd)], 1L, if (nd) sh[nd + 1L] else sh[1L]))
  sh
}

# ---- construction --------------------------------------------------------------

def("rep_vector", c("real", "int"), "vector", function(x, n) {
  chk_size("rep_vector", n)
  ad_gather(as_dbl(x), rep(1L, n))
})
def("rep_row_vector", c("real", "int"), "row_vector", function(x, n) {
  chk_size("rep_row_vector", n)
  ad_gather(as_dbl(x), rep(1L, n))
})
def("rep_matrix", c("real", "int", "int"), "matrix", function(x, m, n) {
  chk_size("rep_matrix", m); chk_size("rep_matrix", n)
  ad_gather(as_dbl(x), rep(1L, m * n), c(m, n))
})
def("rep_matrix", c("vector", "int"), "matrix", function(v, n) {
  chk_size("rep_matrix", n)
  k <- nelem(v)
  ad_gather(as_dbl(v), rep(seq_len(k), n), c(k, n))
})
def("rep_matrix", c("row_vector", "int"), "matrix", function(v, m) {
  chk_size("rep_matrix", m)
  k <- nelem(v)
  ad_gather(as_dbl(v), rep(seq_len(k), each = m), c(m, k))
})
rep_array_impl <- function(x, ds, base) {
  for (d in ds) chk_size("rep_array", d)
  v <- val(x)
  es <- if (base %in% c("int", "real")) integer(0) else if (is.null(dim(v))) length(v) else dim(v)
  ne <- length(v)
  idx <- rep(seq_len(ne), each = prod(ds))
  shp <- c(ds, es)
  if (is.integer(x) && !is_ad(x)) {
    out <- x[idx]
    if (length(shp) > 1L) dim(out) <- shp
    return(out)
  }
  ad_gather(as_dbl(x), idx, if (length(shp) > 1L) shp else NULL)
}
for (b in c("int", "real", "vector", "row_vector", "matrix")) local({
  bb <- b
  def("rep_array", c(bb, "int"), paste0(bb, "[]"), function(x, n) rep_array_impl(x, n, bb))
  def("rep_array", c(bb, "int", "int"), paste0(bb, "[,]"), function(x, n, m) rep_array_impl(x, c(n, m), bb))
  def("rep_array", c(bb, "int", "int", "int"), paste0(bb, "[,,]"), function(x, n, m, k) rep_array_impl(x, c(n, m, k), bb))
  # arrays of arrays: the new dimensions come first
  for (nd in 1:2) local({
    at <- paste0(bb, "[", strrep(",", nd - 1L), "]")
    rt <- function(k) paste0(bb, "[", strrep(",", nd + k - 1L), "]")
    def("rep_array", c(at, "int"), rt(1L), function(x, n) rep_array_impl(x, n, "array"))
    def("rep_array", c(at, "int", "int"), rt(2L), function(x, n, m) rep_array_impl(x, c(n, m), "array"))
    def("rep_array", c(at, "int", "int", "int"), rt(3L), function(x, n, m, k) rep_array_impl(x, c(n, m, k), "array"))
  })
})

chk_size <- function(f, n) {
  if (n < 0L) stop_domain(sprintf("%s: size is %d, but must be nonnegative", f, n))
}

zeros_like <- function(n, shape = NULL) {
  r <- numeric(n)
  if (!is.null(shape)) dim(r) <- shape
  r
}
def("zeros_vector", "int", "vector", function(n) { chk_size("zeros_vector", n); numeric(n) })
def("zeros_row_vector", "int", "row_vector", function(n) { chk_size("zeros_row_vector", n); numeric(n) })
def("zeros_array", "int", "real[]", function(n) { chk_size("zeros_array", n); numeric(n) })
def("zeros_int_array", "int", "int[]", function(n) { chk_size("zeros_int_array", n); integer(n) })
def("ones_vector", "int", "vector", function(n) { chk_size("ones_vector", n); rep(1, n) })
def("ones_row_vector", "int", "row_vector", function(n) { chk_size("ones_row_vector", n); rep(1, n) })
def("ones_array", "int", "real[]", function(n) { chk_size("ones_array", n); rep(1, n) })
def("ones_int_array", "int", "int[]", function(n) { chk_size("ones_int_array", n); rep(1L, n) })
def("uniform_simplex", "int", "vector", function(n) {
  if (n < 1L) stop_domain("uniform_simplex: size must be positive")
  rep(1 / n, n)
})
def("identity_matrix", "int", "matrix", function(n) { chk_size("identity_matrix", n); diag(1, n) })
linsp <- function(n, lo, hi) {
  chk_size("linspaced", n)
  if (n == 0L) return(numeric(0))
  if (n == 1L) return(as.numeric(hi))
  seq(val(lo), val(hi), length.out = n)
}
def("linspaced_vector", c("int", "real", "real"), "vector", function(n, a, b) linsp(n, a, b))
def("linspaced_row_vector", c("int", "real", "real"), "row_vector", function(n, a, b) linsp(n, a, b))
def("linspaced_array", c("int", "real", "real"), "real[]", function(n, a, b) linsp(n, a, b))
def("linspaced_int_array", c("int", "int", "int"), "int[]", function(n, a, b) {
  chk_size("linspaced_int_array", n)
  if (b < a) stop_domain("linspaced_int_array: high must be at least low")
  if (n == 0L) return(integer(0))
  by <- (b - a) %/% max(n - 1L, 1L)
  as.integer(a + by * (seq_len(n) - 1L))
})
one_hot <- function(n, k) {
  if (n < 1L) stop_domain("one_hot: size must be positive")
  if (k < 1L || k > n) stop_domain(sprintf("one_hot: index %d out of range 1..%d", k, n))
  r <- numeric(n)
  r[k] <- 1
  r
}
def("one_hot_vector", c("int", "int"), "vector", one_hot)
def("one_hot_row_vector", c("int", "int"), "row_vector", one_hot)
def("one_hot_array", c("int", "int"), "real[]", one_hot)
def("one_hot_int_array", c("int", "int"), "int[]", function(n, k) as.integer(one_hot(n, k)))

# ---- conversions ---------------------------------------------------------------

# to_vector: matrices column major, arrays of reals in Stan's order. For a
# 2 dimensional array Stan's to_vector/to_array_1d use row major order.
def("to_vector", "vector", "vector", function(x) as_dbl(x))
def("to_vector", "row_vector", "vector", function(x) as_dbl(x))
def("to_vector", "matrix", "vector", function(x) flat(as_dbl(x)))
def("to_vector", "real[]", "vector", function(x) as_dbl(flat(x)))
def("to_vector", "int[]", "vector", function(x) as_dbl(flat(x)))
def("to_row_vector", "vector", "row_vector", function(x) as_dbl(x))
def("to_row_vector", "row_vector", "row_vector", function(x) as_dbl(x))
def("to_row_vector", "matrix", "row_vector", function(x) flat(as_dbl(x)))
def("to_row_vector", "real[]", "row_vector", function(x) as_dbl(flat(x)))
def("to_row_vector", "int[]", "row_vector", function(x) as_dbl(flat(x)))
def("to_matrix", "matrix", "matrix", function(x) as_dbl(x))
def("to_matrix", "vector", "matrix", function(x) ad_reshape(as_dbl(x), c(nelem(x), 1L)))
def("to_matrix", "row_vector", "matrix", function(x) ad_reshape(as_dbl(x), c(1L, nelem(x))))
def("to_matrix", "real[,]", "matrix", function(x) {
  v <- val(x)
  if (is.null(dim(v))) return(ad_reshape(as_dbl(x), c(length(v), 0L)))
  as_dbl(x)
})
def("to_matrix", "int[,]", "matrix", function(x) as_dbl(x))
def("to_matrix", "row_vector[]", "matrix", function(x) as_dbl(x))
def("to_matrix", "vector[]", "matrix", function(x) ad_transpose(as_dbl(x)))
to_matrix_mn <- function(x, m, n, col_major = 1L) {
  if (m * n != nelem(x)) stop_shape(sprintf("to_matrix: %d elements cannot fill a %d x %d matrix", nelem(x), m, n))
  if (col_major) return(ad_reshape(flat(as_dbl(x)), c(m, n)))
  ad_transpose(ad_reshape(flat(as_dbl(x)), c(n, m)))
}
for (b in c("vector", "row_vector", "matrix", "real[]", "int[]")) local({
  bb <- b
  def("to_matrix", c(bb, "int", "int"), "matrix", function(x, m, n) to_matrix_mn(x, m, n))
  def("to_matrix", c(bb, "int", "int", "int"), "matrix", function(x, m, n, cm) to_matrix_mn(x, m, n, cm))
})
row_major_flat <- function(x) {
  v <- val(x)
  d <- dim(v)
  if (is.null(d)) return(flat(x))
  p <- as.vector(aperm(array(seq_along(v), d), rev(seq_along(d))))
  ad_gather(x, p)
}
def("to_array_1d", "vector", "real[]", function(x) as_dbl(x))
def("to_array_1d", "row_vector", "real[]", function(x) as_dbl(x))
def("to_array_1d", "matrix", "real[]", function(x) flat(as_dbl(x)))
def("to_array_1d", "int[]", "int[]", function(x) x)
def("to_array_1d", "int[,]", "int[]", function(x) as.integer(val(row_major_flat(x))))
def("to_array_1d", "int[,,]", "int[]", function(x) as.integer(val(row_major_flat(x))))
def("to_array_1d", "real[]", "real[]", function(x) as_dbl(x))
def("to_array_1d", "real[,]", "real[]", function(x) row_major_flat(as_dbl(x)))
def("to_array_1d", "real[,,]", "real[]", function(x) row_major_flat(as_dbl(x)))
for (nd in 4:8) local({
  cm <- strrep(",", nd - 1L)
  def("to_array_1d", sprintf("int[%s]", cm), "int[]", function(x) as.integer(val(row_major_flat(x))))
  def("to_array_1d", sprintf("real[%s]", cm), "real[]", function(x) row_major_flat(as_dbl(x)))
})
def("to_array_2d", "matrix", "real[,]", function(x) as_dbl(x))
def("to_vector_array", "matrix", "vector[]", function(x) ad_transpose(as_dbl(x)))
def("to_row_vector_array", "matrix", "row_vector[]", function(x) as_dbl(x))

# ---- slicing -------------------------------------------------------------------

slice_check <- function(f, start, len, n) {
  if (len < 0L) stop_domain(sprintf("%s: length is %d, but must be nonnegative", f, len))
  if (len > 0L && (start < 1L || start + len - 1L > n))
    stop_domain(sprintf("%s: requested elements %d to %d of a container of size %d", f, start, start + len - 1L, n))
}
take_first <- function(x, idx) {
  v <- val(x)
  d <- dim(v)
  if (is.null(d)) {
    if (is.integer(x)) return(x[idx])
    return(ad_gather(as_dbl(x), idx))
  }
  # array whose first dimension is sliced
  rest <- d[-1L]
  all_idx <- as.vector(outer(idx, (seq_len(prod(rest)) - 1L) * d[1L], `+`))
  shp <- c(length(idx), rest)
  if (is.integer(x)) {
    out <- x[all_idx]
    dim(out) <- shp
    return(out)
  }
  ad_gather(as_dbl(x), all_idx, shp)
}
for (b in c("vector", "row_vector", "anyarr")) local({
  bb <- b
  rr <- if (bb == "anyarr") r_same1 else bb
  def("head", c(bb, "int"), rr, function(x, n) {
    sz <- shape_of(x)[1L]
    slice_check("head", 1L, n, sz)
    take_first(x, seq_len(n))
  })
  def("tail", c(bb, "int"), rr, function(x, n) {
    sz <- shape_of(x)[1L]
    slice_check("tail", sz - n + 1L, n, sz)
    take_first(x, if (n > 0L) (sz - n + 1L):sz else integer(0))
  })
  def("segment", c(bb, "int", "int"), rr, function(x, i, n) {
    sz <- shape_of(x)[1L]
    slice_check("segment", i, n, sz)
    take_first(x, if (n > 0L) i:(i + n - 1L) else integer(0))
  })
})
def("block", c("matrix", "int", "int", "int", "int"), "matrix", function(x, i, j, nr, nc) {
  d <- dim(val(x))
  slice_check("block", i, nr, d[1L]); slice_check("block", j, nc, d[2L])
  rows <- if (nr) i:(i + nr - 1L) else integer(0)
  cols <- if (nc) j:(j + nc - 1L) else integer(0)
  ad_gather(as_dbl(x), as.vector(outer(rows, (cols - 1L) * d[1L], `+`)), c(nr, nc))
})
def("col", c("matrix", "int"), "vector", function(x, j) {
  d <- dim(val(x))
  if (j < 1L || j > d[2L]) stop_domain(sprintf("col: column %d out of range 1..%d", j, d[2L]))
  ad_gather(as_dbl(x), (j - 1L) * d[1L] + seq_len(d[1L]))
})
def("row", c("matrix", "int"), "row_vector", function(x, i) {
  d <- dim(val(x))
  if (i < 1L || i > d[1L]) stop_domain(sprintf("row: row %d out of range 1..%d", i, d[1L]))
  ad_gather(as_dbl(x), i + (seq_len(d[2L]) - 1L) * d[1L])
})
def("sub_col", c("matrix", "int", "int", "int"), "vector", function(x, i, j, n) {
  d <- dim(val(x))
  slice_check("sub_col", i, n, d[1L])
  if (j < 1L || j > d[2L]) stop_domain("sub_col: column out of range")
  ad_gather(as_dbl(x), (j - 1L) * d[1L] + if (n) i:(i + n - 1L) else integer(0))
})
def("sub_row", c("matrix", "int", "int", "int"), "row_vector", function(x, i, j, n) {
  d <- dim(val(x))
  slice_check("sub_row", j, n, d[2L])
  if (i < 1L || i > d[1L]) stop_domain("sub_row: row out of range")
  ad_gather(as_dbl(x), i + ((if (n) j:(j + n - 1L) else integer(0)) - 1L) * d[1L])
})
def("diagonal", "matrix", "vector", function(x) {
  d <- dim(val(x))
  k <- min(d)
  ad_gather(as_dbl(x), (seq_len(k) - 1L) * d[1L] + seq_len(k))
})
def("diag_matrix", "vector", "matrix", function(v) {
  k <- nelem(v)
  ad_scatter(matrix(0, k, k), (seq_len(k) - 1L) * k + seq_len(k), as_dbl(v))
})
def("add_diag", c("matrix", "real"), "matrix", function(m, d) add_diag_impl(m, d))
def("add_diag", c("matrix", "vector"), "matrix", function(m, d) add_diag_impl(m, d))
def("add_diag", c("matrix", "row_vector"), "matrix", function(m, d) add_diag_impl(m, d))
add_diag_impl <- function(m, d) {
  dm <- dim(val(m))
  k <- min(dm)
  if (nelem(d) != 1L && nelem(d) != k) stop_shape("add_diag: size of the diagonal must match the matrix")
  idx <- (seq_len(k) - 1L) * dm[1L] + seq_len(k)
  m <- as_dbl(m)
  ad_scatter(m, idx, ad_gather(m, idx) + as_dbl(d))
}

def("append_row", c("matrix", "matrix"), "matrix", function(a, b) {
  if (ncol(val(a)) != ncol(val(b))) stop_shape("append_row: matrices must have the same number of columns")
  ad_transpose(ad_reshape(ad_c(list(ad_transpose(as_dbl(a)), ad_transpose(as_dbl(b)))), c(ncol(val(a)), nrow(val(a)) + nrow(val(b)))))
})
def("append_row", c("matrix", "row_vector"), "matrix", function(a, b) .lib$append_row[[1L]]$impl(a, as_row(b)))
def("append_row", c("row_vector", "matrix"), "matrix", function(a, b) .lib$append_row[[1L]]$impl(as_row(a), b))
def("append_row", c("row_vector", "row_vector"), "matrix", function(a, b) .lib$append_row[[1L]]$impl(as_row(a), as_row(b)))
def("append_row", c("vector", "vector"), "vector", function(a, b) ad_c(list(as_dbl(a), as_dbl(b))))
def("append_row", c("real", "vector"), "vector", function(a, b) ad_c(list(as_dbl(a), as_dbl(b))))
def("append_row", c("vector", "real"), "vector", function(a, b) ad_c(list(as_dbl(a), as_dbl(b))))
def("append_col", c("matrix", "matrix"), "matrix", function(a, b) {
  if (nrow(val(a)) != nrow(val(b))) stop_shape("append_col: matrices must have the same number of rows")
  ad_reshape(ad_c(list(as_dbl(a), as_dbl(b))), c(nrow(val(a)), ncol(val(a)) + ncol(val(b))))
})
def("append_col", c("matrix", "vector"), "matrix", function(a, b) .lib$append_col[[1L]]$impl(a, as_col(b)))
def("append_col", c("vector", "matrix"), "matrix", function(a, b) .lib$append_col[[1L]]$impl(as_col(a), b))
def("append_col", c("vector", "vector"), "matrix", function(a, b) .lib$append_col[[1L]]$impl(as_col(a), as_col(b)))
def("append_col", c("row_vector", "row_vector"), "row_vector", function(a, b) ad_c(list(as_dbl(a), as_dbl(b))))
def("append_col", c("real", "row_vector"), "row_vector", function(a, b) ad_c(list(as_dbl(a), as_dbl(b))))
def("append_col", c("row_vector", "real"), "row_vector", function(a, b) ad_c(list(as_dbl(a), as_dbl(b))))
def("append_array", c("anyarr", "anyarr"), function(types) {
  a <- types[[1L]]; b <- types[[2L]]
  if (a$nd != b$nd) return(NULL)
  if (a$base == b$base) return(a)
  if (a$base %in% c("int", "real") && b$base %in% c("int", "real")) return(ty("real", a$nd))
  NULL
}, function(a, b) {
  da <- shape_of(a); db <- shape_of(b)
  if (length(da) != length(db) || any(da[-1L] != db[-1L])) stop_shape("append_array: arrays have incompatible sizes")
  if (is.integer(a) && is.integer(b)) {
    if (length(da) == 1L) return(c(a, b))
  }
  n1 <- da[1L]; n2 <- db[1L]
  rest <- prod(da[-1L])
  # put the first dimension last, concatenate, then move it back
  pa <- if (length(da) > 1L) aperm_flat(a, da) else flat(as_dbl(a))
  pb <- if (length(db) > 1L) aperm_flat(b, db) else flat(as_dbl(b))
  out <- ad_c(list(as_dbl(pa), as_dbl(pb)))
  if (length(da) == 1L) return(out)
  nshape <- c(da[-1L], n1 + n2)
  idx <- array(seq_len(prod(nshape)), nshape)
  p <- as.vector(aperm(idx, c(length(nshape), seq_len(length(nshape) - 1L))))
  r <- ad_gather(out, p, c(n1 + n2, da[-1L]))
  if (is.integer(a) && is.integer(b)) { r <- as.integer(r); dim(r) <- c(n1 + n2, da[-1L]) }
  r
})
aperm_flat <- function(x, d) {
  idx <- array(seq_len(prod(d)), d)
  p <- as.vector(aperm(idx, c(seq_along(d)[-1L], 1L)))
  ad_gather(as_dbl(x), p)
}

def("reverse", "vector", "vector", function(x) take_first(x, rev(seq_len(nelem(x)))))
def("reverse", "row_vector", "row_vector", function(x) take_first(x, rev(seq_len(nelem(x)))))
def("reverse", "anyarr", r_same1, function(x) take_first(x, rev(seq_len(shape_of(x)[1L]))))

# ---- sorting --------------------------------------------------------------------

sort_impl <- function(x, decreasing) {
  v <- as.vector(val(x))
  if (any(is.na(v))) stop_domain("sort: argument contains NaN")
  o <- order(v, decreasing = decreasing)
  if (is.integer(x)) return(x[o])
  ad_gather(as_dbl(x), o)
}
for (b in c("vector", "row_vector", "real[]", "int[]")) local({
  bb <- b
  def("sort_asc", bb, bb, function(x) sort_impl(x, FALSE))
  def("sort_desc", bb, bb, function(x) sort_impl(x, TRUE))
  def("sort_indices_asc", bb, "int[]", function(x) order(as.vector(val(x))))
  def("sort_indices_desc", bb, "int[]", function(x) order(as.vector(val(x)), decreasing = TRUE))
  def("rank", c(bb, "int"), "int", function(x, s) {
    v <- as.vector(val(x))
    if (s < 1L || s > length(v)) stop_domain("rank: index out of range")
    as.integer(sum(v < v[s]))
  })
})
def("quantile", c("real[]", "real"), "real", function(x, p) quantile_stan(val(x), val(p)))
def("quantile", c("vector", "real"), "real", function(x, p) quantile_stan(val(x), val(p)))
def("quantile", c("row_vector", "real"), "real", function(x, p) quantile_stan(val(x), val(p)))
def("quantile", c("real[]", "real[]"), "real[]", function(x, p) quantile_stan(val(x), val(p)))
def("quantile", c("vector", "vector"), "vector", function(x, p) quantile_stan(val(x), val(p)))
def("quantile", c("vector", "real[]"), "real[]", function(x, p) quantile_stan(val(x), val(p)))
def("quantile", c("row_vector", "real[]"), "real[]", function(x, p) quantile_stan(val(x), val(p)))
quantile_stan <- function(x, p) {
  if (!length(x)) stop_domain("quantile: container has size zero")
  if (any(p < 0 | p > 1)) stop_domain("quantile: probabilities must be in [0, 1]")
  # Stan uses the sample quantile of type 7
  as.numeric(stats::quantile(x, p, type = 7, names = FALSE))
}

# ---- softmax -------------------------------------------------------------------

softmax_impl <- function(x) {
  x <- as_dbl(x)
  if (!nelem(x)) stop_domain("softmax: argument has size zero")
  exp(x - ad_log_sum_exp(x))
}
log_softmax_impl <- function(x) {
  x <- as_dbl(x)
  if (!nelem(x)) return(numeric(0))
  x - ad_log_sum_exp(x)
}
def("softmax", "vector", "vector", softmax_impl)
def("softmax", "row_vector", "row_vector", softmax_impl)
def("log_softmax", "vector", "vector", log_softmax_impl)
def("log_softmax", "row_vector", "row_vector", log_softmax_impl)
rowwise_apply <- function(x, f) {
  v <- val(x)
  if (is.null(dim(v))) return(f(x))
  n <- nrow(v)
  k <- ncol(v)
  rows <- lapply(seq_len(n), function(i) f(ad_gather(as_dbl(x), i + (seq_len(k) - 1L) * n)))
  ad_transpose(ad_reshape(ad_c(rows), c(k, n)))
}
def("softmax", "vector[]", "vector[]", function(x) rowwise_apply(x, softmax_impl))
def("softmax", "row_vector[]", "row_vector[]", function(x) rowwise_apply(x, softmax_impl))
def("log_softmax", "vector[]", "vector[]", function(x) rowwise_apply(x, log_softmax_impl))
def("log_softmax", "row_vector[]", "row_vector[]", function(x) rowwise_apply(x, log_softmax_impl))
