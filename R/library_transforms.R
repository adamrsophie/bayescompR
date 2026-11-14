# The transform functions of Stan 2.37 and later: for each constrained type
# T, T_constrain(y, ...) maps unconstrained values to the constrained space,
# T_jacobian(y, ...) does the same and adds the log Jacobian to the target
# (only in transformed parameters and _jacobian functions), and
# T_unconstrain(x, ...) is the inverse. They reuse the parameter transforms.

# return type: the first argument promoted to real; bound arguments are a
# real or of the same type as the first
r_bound_args <- function(types) {
  y <- r_promote1(types[1L])
  for (t in types[-1L]) if (!is_scalar_ty(t) && !ty_assignable(y, t)) return(NULL)
  y
}

bound_arg <- function(f, b, y) {
  vb <- val(b)
  if (length(vb) != 1L && length(vb) != length(val(y))) stop_shape(sprintf("%s: the bound must be a real or have the size of the argument", f))
  as_dbl(b)
}

keep_dim <- function(x, like) {
  d <- dim(val(like))
  if (is.null(d) || identical(dim(val(x)), d)) return(x)
  ad_reshape(x, d)
}

bound_constrain <- function(f, y, a, ctx = NULL) {
  r <- tf_bounds$constrain(as_dbl(y), NULL, a)
  if (!is.null(ctx) && isTRUE(ctx$jacobian)) rt_target_add(ctx, r$lj)
  keep_dim(r$x, y)
}

# exact differentiable inverses of the bound transforms
bound_unconstrain <- function(f, x, lo = NULL, up = NULL) {
  x <- as_dbl(x)
  vx <- val(x)
  n <- length(vx)
  z <- 0 * x
  flo <- if (is.null(lo)) rep(FALSE, n) else rep_len(is.finite(val(lo)), n)
  fup <- if (is.null(up)) rep(FALSE, n) else rep_len(is.finite(val(up)), n)
  safe <- function(b, fin) if (is.null(b)) z else ad_ifelse(fin, b + z, z)
  l <- safe(lo, flo); u <- safe(up, fup)
  if (any(flo & vx < val(l))) stop_domain(sprintf("%s: the argument is below the lower bound", f))
  if (any(fup & vx > val(u))) stop_domain(sprintf("%s: the argument is above the upper bound", f))
  one <- 1 + z
  out <- x
  if (any(flo & !fup)) out <- ad_ifelse(flo & !fup, log(ad_ifelse(flo & !fup, x - l, one)), out)
  if (any(fup & !flo)) out <- ad_ifelse(fup & !flo, log(ad_ifelse(fup & !flo, u - x, one)), out)
  if (any(flo & fup)) {
    both <- flo & fup
    w <- ad_ifelse(both, u - l, one)
    p <- ad_ifelse(both, (x - l) / w, 0.5 + z)
    out <- ad_ifelse(both, ad_logit(p), out)
  }
  keep_dim(out, x)
}

def_bounds <- function(stem, nb, constrain, unconstrain) {
  for (suf in c("constrain", "jacobian", "unconstrain")) local({
    s <- suf
    nm <- paste0(stem, "_", s)
    impl <- switch(s,
      constrain = function(y, ...) constrain(nm, y, list(...), NULL),
      jacobian = function(y, ..., .ctx) constrain(nm, y, utils::head(list(...), nb), .ctx),
      unconstrain = function(x, ...) unconstrain(nm, x, list(...)))
    def(nm, c("any", rep("any", nb)), r_bound_args, impl, ctx = s == "jacobian")
  })
}
def_bounds("lower_bound", 1L,
  function(f, y, b, ctx) bound_constrain(f, y, list(lower = bound_arg(f, b[[1L]], y)), ctx),
  function(f, x, b) bound_unconstrain(f, x, lo = bound_arg(f, b[[1L]], x)))
def_bounds("upper_bound", 1L,
  function(f, y, b, ctx) bound_constrain(f, y, list(upper = bound_arg(f, b[[1L]], y)), ctx),
  function(f, x, b) bound_unconstrain(f, x, up = bound_arg(f, b[[1L]], x)))
def_bounds("lower_upper_bound", 2L,
  function(f, y, b, ctx) {
    lo <- bound_arg(f, b[[1L]], y); up <- bound_arg(f, b[[2L]], y)
    if (any(val(lo) >= val(up))) stop_domain(sprintf("%s: the lower bound must be less than the upper bound", f))
    bound_constrain(f, y, list(lower = lo, upper = up), ctx)
  },
  function(f, x, b) bound_unconstrain(f, x, lo = bound_arg(f, b[[1L]], x), up = bound_arg(f, b[[2L]], x)))
def_bounds("offset_multiplier", 2L,
  function(f, y, b, ctx) bound_constrain(f, y, list(offset = bound_arg(f, b[[1L]], y), multiplier = bound_arg(f, b[[2L]], y)), ctx),
  function(f, x, b) {
    o <- bound_arg(f, b[[1L]], x); m <- bound_arg(f, b[[2L]], x)
    if (any(!is.finite(val(m)) | val(m) <= 0)) stop_domain(sprintf("%s: multiplier must be positive and finite", f))
    keep_dim((as_dbl(x) - o) / m, x)
  })

# ---- structured transforms ----------------------------------------------------------

# inverse with gradients: on the constrained manifold du/dx is the
# pseudoinverse of the Jacobian of the constraining map
unconstrain_ad <- function(tf, x, shape) {
  u <- tf$unconstrain(val(x), shape, NULL)
  if (!is_ad(x) || !length(u)) return(u)
  vx <- val(x)
  nx <- length(vx)
  J <- matrix(0, nx, length(u))
  for (j in seq_len(nx)) J[j, ] <- ad_value_grad(function(uu) ad_select(flat(tf$constrain(uu, shape, NULL)$x), j), u)$grad
  P <- J %*% solve(crossprod(J))
  ad_node(u, list(x), function(g) {
    gx <- P %*% g
    dim(gx) <- dim(vx)
    list(if (is.null(dim(vx))) as.vector(gx) else gx)
  })
}

# spec: free shape from the unconstrained argument (and the int arguments),
# constrained shape, and the constrained shape from a constrained value
tf_specs <- list(
  ordered = list(tf = tf_ordered, y = "vector", x = "vector", shape = function(v, n) length(v), xshape = function(v) length(v)),
  positive_ordered = list(tf = tf_positive_ordered, y = "vector", x = "vector", shape = function(v, n) length(v), xshape = function(v) length(v)),
  simplex = list(tf = tf_simplex, y = "vector", x = "vector", shape = function(v, n) length(v) + 1L, xshape = function(v) length(v)),
  unit_vector = list(tf = tf_unit_vector, y = "vector", x = "vector", shape = function(v, n) length(v), xshape = function(v) length(v)),
  sum_to_zero = list(tf = tf_sum_to_zero, y = "vector", x = "vector", shape = function(v, n) length(v) + 1L, xshape = function(v) length(v)),
  sum_to_zero_matrix = list(name = "sum_to_zero", tf = tf_sum_to_zero_matrix, y = "matrix", x = "matrix",
                            shape = function(v, n) dim(v) + 1L, xshape = function(v) dim(v)),
  stochastic_column = list(tf = tf_column_stochastic, y = "matrix", x = "matrix", shape = function(v, n) dim(v) + c(1L, 0L), xshape = function(v) dim(v)),
  stochastic_row = list(tf = tf_row_stochastic, y = "matrix", x = "matrix", shape = function(v, n) dim(v) + c(0L, 1L), xshape = function(v) dim(v)),
  cholesky_factor_corr = list(tf = tf_cholesky_corr, y = "vector", x = "matrix", ints = 1L, shape = function(v, n) c(n[[1L]], n[[1L]]), xshape = function(v) dim(v)),
  cholesky_factor_cov = list(tf = tf_cholesky_cov, y = "vector", x = "matrix", ints = 2L, shape = function(v, n) c(n[[1L]], n[[2L]]), xshape = function(v) dim(v)),
  corr_matrix = list(tf = tf_corr_matrix, y = "vector", x = "matrix", ints = 1L, shape = function(v, n) c(n[[1L]], n[[1L]]), xshape = function(v) dim(v)),
  cov_matrix = list(tf = tf_cov_matrix, y = "vector", x = "matrix", ints = 1L, shape = function(v, n) c(n[[1L]], n[[1L]]), xshape = function(v) dim(v)))

tf_free_check <- function(f, sp, shape, nfree) {
  if (any(shape < 0L)) stop_domain(sprintf("%s: sizes must be nonnegative", f))
  if (sp$tf$free_size(shape) != nfree) stop_shape(sprintf("%s: the unconstrained argument has %d elements, but %d are needed", f, nfree, sp$tf$free_size(shape)))
}

# apply fun to each element of an array (value dims: array dims, then element dims)
map_array_elems <- function(x, nd, edim, fun) {
  vx <- val(x)
  d <- if (is.null(dim(vx))) length(vx) else dim(vx)
  ad_dims <- d[seq_len(nd)]
  nA <- prod(ad_dims)
  ne <- prod(edim)
  outs <- lapply(seq_len(nA), function(i) {
    xi <- ad_gather(as_dbl(x), i + (seq_len(ne) - 1L) * nA, if (length(edim) > 1L) edim else NULL)
    fun(xi)
  })
  od <- if (length(outs)) { v1 <- val(outs[[1L]]); if (is.null(dim(v1))) length(v1) else dim(v1) } else 0L
  no <- prod(od)
  t <- seq_len(nA * no)
  src <- ((t - 1L) %% nA) * no + (t - 1L) %/% nA + 1L
  if (!nA) return(array(0, c(ad_dims, od)))
  ad_gather(ad_c(lapply(outs, function(o) flat(as_dbl(o)))), src, c(ad_dims, od))
}

for (key in names(tf_specs)) local({
  sp <- tf_specs[[key]]
  stem <- sp$name %||% key
  ni <- sp$ints %||% 0L
  ints <- rep("int", ni)
  elem_dim <- function(v, nd) { d <- if (is.null(dim(v))) length(v) else dim(v); d[-seq_len(nd)] }
  for (nd in 0:7) local({
    k <- nd
    arr <- function(b) if (k) paste0(b, "[", strrep(",", k - 1L), "]") else b
    one_constrain <- function(f, y, n, ctx) {
      vy <- val(y)
      shape <- as.integer(sp$shape(vy, n))
      tf_free_check(f, sp, shape, length(vy))
      r <- sp$tf$constrain(flat(as_dbl(y)), shape, NULL)
      if (!is.null(ctx) && isTRUE(ctx$jacobian)) rt_target_add(ctx, r$lj)
      x <- r$x
      if (length(shape) > 1L && is.null(dim(val(x)))) x <- ad_reshape(as_dbl(x), shape)
      x
    }
    do_constrain <- function(f, y, n, ctx) {
      if (!k) return(one_constrain(f, y, n, ctx))
      map_array_elems(y, k, elem_dim(val(y), k), function(yi) one_constrain(f, yi, n, ctx))
    }
    one_unconstrain <- function(f, x) {
      vx <- val(x)
      shape <- as.integer(sp$xshape(vx))
      msg <- sp$tf$check(vx, shape, NULL)
      if (!is.null(msg)) stop_domain(sprintf("%s: the argument %s", f, msg))
      u <- unconstrain_ad(sp$tf, x, shape)
      if (sp$y == "matrix") {
        fd <- switch(stem, sum_to_zero = shape - 1L, stochastic_column = shape - c(1L, 0L), stochastic_row = shape - c(0L, 1L))
        u <- if (is_ad(u)) ad_reshape(u, fd) else array(u, fd)
      }
      u
    }
    do_unconstrain <- function(f, x) {
      if (!k) return(one_unconstrain(f, x))
      map_array_elems(x, k, elem_dim(val(x), k), function(xi) one_unconstrain(f, xi))
    }
    nm_c <- paste0(stem, "_constrain"); nm_j <- paste0(stem, "_jacobian"); nm_u <- paste0(stem, "_unconstrain")
    def(nm_c, c(arr(sp$y), ints), arr(sp$x), function(y, ...) do_constrain(nm_c, y, list(...), NULL))
    def(nm_j, c(arr(sp$y), ints), arr(sp$x), function(y, ..., .ctx) do_constrain(nm_j, y, list(...), .ctx), ctx = TRUE)
    def(nm_u, arr(sp$x), arr(sp$y), function(x) do_unconstrain(nm_u, x))
  })
})
