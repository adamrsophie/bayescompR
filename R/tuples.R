# Tuples.
#
# A tuple value is a list of its element values with class "bctuple". An
# array of tuples is a list of elements with class "bctarr", nested once per
# array dimension. Declarations of tuple type carry a descriptor: a list with
# base "tuple", the array sizes and the descriptors of the elements, and for
# the other elements their base type, array sizes and element shape.

rt_tuple <- function(items) structure(items, class = "bctuple")
rt_tarr <- function(items) structure(items, class = "bctarr")
is_tuple_val <- function(x) inherits(x, "bctuple")
is_tarr <- function(x) inherits(x, "bctarr")

rt_proj <- function(x, k) unclass(x)[[k]]

rt_set_proj <- function(x, k, v, name) {
  y <- unclass(x)
  y[[k]] <- rt_assign(y[[k]], v, name)
  rt_tuple(y)
}

# value of a declaration descriptor, with Stan's NaN filled elements
rt_decl_desc <- function(d) {
  one <- if (identical(d$base, "tuple")) function() rt_tuple(lapply(d$elems, rt_decl_desc))
         else return(rt_decl(if (decl_base_type(d$base) == "int") "int" else "real", d$adims, d$eshape))
  build <- function(dims) {
    if (!length(dims)) return(one())
    if (dims[1L] < 0L) stop_domain(sprintf("declared size is negative (%d)", dims[1L]))
    rt_tarr(lapply(seq_len(dims[1L]), function(i) build(dims[-1L])))
  }
  build(as.integer(d$adims))
}

# ---- indexing and assignment of arrays of tuples ---------------------------------------

tarr_sel <- function(n, kind, v, what) {
  s <- switch(kind,
    single = v, intarray = v,
    range = if (v[2L] < v[1L]) integer(0) else v[1L]:v[2L],
    lower = if (n < v) integer(0) else v:n,
    upper = if (v < 1L) integer(0) else seq_len(v),
    omni = seq_len(n))
  s <- as.integer(s)
  if (length(s) && (anyNA(s) || min(s) < 1L || max(s) > n)) {
    bad <- s[is.na(s) | s < 1L | s > n][1L]
    stop_domain(sprintf("%s out of range: index %s for a dimension of size %d", what, bad, n))
  }
  s
}

tarr_index <- function(x, kinds, vals) {
  if (!length(kinds)) return(x)
  s <- tarr_sel(length(x), kinds[1L], vals[[1L]], "index")
  rest_k <- kinds[-1L]; rest_v <- vals[-1L]
  if (kinds[1L] == "single") return(tarr_index(unclass(x)[[s]], rest_k, rest_v))
  rt_tarr(lapply(s, function(i) tarr_index(unclass(x)[[i]], rest_k, rest_v)))
}

tarr_assign <- function(x, kinds, vals, rhs, name) {
  if (!length(kinds)) return(rt_assign(x, rhs, name))
  s <- tarr_sel(length(x), kinds[1L], vals[[1L]], sprintf("assignment to '%s': index", name))
  y <- unclass(x)
  rest_k <- kinds[-1L]; rest_v <- vals[-1L]
  if (kinds[1L] == "single") {
    y[[s]] <- tarr_assign(y[[s]], rest_k, rest_v, rhs, name)
  } else {
    r <- unclass(rhs)
    if (length(r) != length(s))
      stop_shape(sprintf("assignment to '%s': size mismatch, the left side selects %d elements but the right side has %d", name, length(s), length(r)))
    for (j in seq_along(s)) y[[s[j]]] <- tarr_assign(y[[s[j]]], rest_k, rest_v, r[[j]], name)
  }
  rt_tarr(y)
}

# assignment of a whole tuple or array of tuples (sizes must agree)
tuple_assign <- function(old, rhs, name) {
  o <- unclass(old); r <- unclass(rhs)
  if (length(o) != length(r))
    stop_shape(sprintf("assignment to '%s': size mismatch, the variable has size %d but the value has size %d", name, length(o), length(r)))
  out <- Map(function(a, b) rt_assign(a, b, name), o, r)
  if (is_tarr(old)) rt_tarr(out) else rt_tuple(out)
}

# array expression {t1, t2, ...} of tuples or arrays of tuples
rt_tuple_array <- function(items) {
  n1 <- length(unclass(items[[1L]]))
  if (is_tarr(items[[1L]]) && any(vapply(items, function(it) length(unclass(it)) != n1, TRUE)))
    stop_shape("array expression elements must all have the same size")
  rt_tarr(items)
}

# promotion of int to real and real to complex inside tuples
rt_promote <- function(x, from, to) {
  if (ty_eq(from, to)) return(x)
  if (to$nd > 0L) {
    if (!is_tuple_ty(to)) return(if (is_cx_ty(to)) cx_promote(x, FALSE) else as_dbl(x))
    return(rt_tarr(lapply(unclass(x), function(e) rt_promote(e, ty_nd(from, from$nd - 1L), ty_nd(to, to$nd - 1L)))))
  }
  if (is_tuple_ty(to)) return(rt_tuple(Map(rt_promote, unclass(x), from$elems, to$elems)))
  if (is_cx_ty(to) && !is_cx_ty(from)) return(cx_promote(x, is_scalar_ty(from)))
  if (to$base == "real" && from$base == "int") return(as_dbl(x))
  x
}

# ---- printing ---------------------------------------------------------------------------

fmt_typed <- function(x, t) {
  if (t$nd > 0L && is_tuple_ty(t))
    return(paste0("[", paste(vapply(unclass(x), fmt_typed, "", t = ty_nd(t, t$nd - 1L)), collapse = ","), "]"))
  if (is_tuple_ty(t)) return(paste0("(", paste(unlist(Map(fmt_typed, unclass(x), t$elems)), collapse = ","), ")"))
  if (is_cx_ty(t)) return(cx_fmt(x, t$nd))
  if (!is_scalar_ty(t)) return(fmt_container(x))
  fmt_value(x)
}

# ---- declarations: metadata, data, checks, transforms -----------------------------------

# descriptor of a sized tuple type with sizes evaluated in env
tuple_info <- function(t, env, name) {
  eval_desc <- function(t) {
    adims <- as.integer(eval(cl_list(as.name("c"), lapply(t$adims, cg_expr)), env))
    if (identical(t$base, "tuple"))
      return(list(name = name, kind = "tuple", base = "tuple", stan_type = "tuple", adims = adims,
                  eshape = integer(0), elems = lapply(t$elems, eval_desc)))
    sh <- cg_decl_shape(t)
    list(name = name, kind = decl_kind(t), base = decl_base_type(t$base), stan_type = t$base,
         adims = adims, eshape = as.integer(eval(sh$eshape, env)))
  }
  info <- eval_desc(t)
  set_free <- function(d) {
    if (identical(d$kind, "tuple")) {
      d$elems <- lapply(d$elems, set_free)
      d$free_size <- as.integer(prod(d$adims) * sum(vapply(d$elems, `[[`, 1L, "free_size")))
    } else d$free_size <- free_size_variable(d)
    d
  }
  set_free(info)
}

# code building a descriptor for a local declaration
cg_decl_shape_tuple <- function(t) {
  if (!identical(t$base, "tuple")) {
    sh <- cg_decl_shape(t)
    return(cl(as.name("list"), base = t$base, adims = sh$adims, eshape = sh$eshape))
  }
  cl(as.name("list"), base = "tuple", adims = cl_list(as.name("c"), lapply(t$adims, cg_expr)),
     elems = cl_list(as.name("list"), lapply(t$elems, cg_decl_shape_tuple)))
}

tuple_read_data <- function(x, d) {
  nm <- d$name
  one <- function(v) {
    v <- unclass(v)
    if (!is.list(v) || length(v) != length(d$elems))
      stop_data(sprintf("variable '%s': a tuple of %d elements needs a list of %d values", nm, length(d$elems), length(d$elems)))
    if (!is.null(names(v)) && all(names(v) %in% as.character(seq_along(d$elems)))) v <- v[as.character(seq_along(d$elems))]
    rt_tuple(Map(function(e, ed) {
      if (identical(ed$kind, "tuple")) return(tuple_read_data(e, utils::modifyList(ed, list(name = nm))))
      dl <- list(e); names(dl) <- nm
      read_data_value(dl, utils::modifyList(ed, list(name = nm)))
    }, v, d$elems))
  }
  build <- function(v, dims) {
    if (!length(dims)) return(one(v))
    if (!is.list(v) || length(v) != dims[1L])
      stop_data(sprintf("variable '%s': expected a list of %d tuples", nm, dims[1L]))
    rt_tarr(lapply(v, build, dims = dims[-1L]))
  }
  # a multidimensional array may also come as one flat list of its tuples in
  # column major order, as bayescomp_constrain() returns it
  if (length(d$adims) > 1L && is.list(x) && length(x) == prod(d$adims) && length(x) != d$adims[1L]) {
    flat <- unclass(x)
    sub <- arrayInd(seq_along(flat), d$adims)
    nest <- function(prefix) {
      k <- length(prefix) + 1L
      if (k > length(d$adims)) {
        hit <- which(apply(sub, 1L, function(r) all(r == prefix)))
        return(flat[[hit]])
      }
      lapply(seq_len(d$adims[k]), function(i) nest(c(prefix, i)))
    }
    x <- nest(integer(0))
  }
  build(x, d$adims)
}

tuple_check <- function(x, d, args) {
  for (el in tuple_elems_of(x, d)) {
    for (k in seq_along(d$elems)) {
      ed <- d$elems[[k]]
      m <- if (identical(ed$kind, "tuple")) tuple_check(el[[k]], ed, args[[k]])
           else check_variable(val(el[[k]]), ed, lapply(args[[k]], val))
      if (!is.null(m)) return(sprintf("element %d of the tuple %s", k, m))
    }
  }
  NULL
}
tuple_elems_of <- function(x, d) {
  if (!length(d$adims)) return(list(unclass(x)))
  unlist(lapply(unclass(x), function(e) tuple_elems_of(e, utils::modifyList(d, list(adims = d$adims[-1L])))), recursive = FALSE)
}

tuple_constrain <- function(u, d, args) {
  off <- 0L
  lj <- 0
  take <- function(n) {
    r <- if (n > 0L) ad_gather(u, off + seq_len(n)) else numeric(0)
    off <<- off + n
    r
  }
  one <- function(d) {
    out <- vector("list", length(d$elems))
    for (k in seq_along(d$elems)) {
      ed <- d$elems[[k]]
      if (identical(ed$kind, "tuple")) {
        r <- tuple_constrain(take(ed$free_size), ed, args[[k]])
      } else {
        r <- constrain_variable(take(ed$free_size), ed, args[[k]])
      }
      out[[k]] <- r$x
      lj <<- lj + r$lj
    }
    rt_tuple(out)
  }
  build <- function(dims) {
    if (!length(dims)) return(one(d))
    rt_tarr(lapply(seq_len(dims[1L]), function(i) build(dims[-1L])))
  }
  x <- build(d$adims)
  list(x = x, lj = lj)
}

tuple_unconstrain <- function(x, d, args) {
  unlist(lapply(tuple_elems_of(x, d), function(el) unlist(lapply(seq_along(d$elems), function(k) {
    ed <- d$elems[[k]]
    if (identical(ed$kind, "tuple")) tuple_unconstrain(el[[k]], ed, args[[k]])
    else unconstrain_variable(val(el[[k]]), ed, lapply(args[[k]], val))
  }))))
}

# ---- output -------------------------------------------------------------------------------

# value written to the output: plain lists, complex elements as R complex
# (arrays of two or more dimensions become list arrays, in column major order)
tuple_out <- function(x, d) {
  if (length(d$adims) == 1L) return(lapply(unclass(x), tuple_out, d = utils::modifyList(d, list(adims = integer(0)))))
  if (length(d$adims) > 1L) {
    inner <- utils::modifyList(d, list(adims = integer(0)))
    out <- vector("list", prod(d$adims))
    dim(out) <- d$adims
    g <- as.matrix(expand.grid(lapply(d$adims, seq_len)))
    for (r in seq_len(nrow(g))) {
      e <- x
      for (k in g[r, ]) e <- unclass(e)[[k]]
      out[[r]] <- tuple_out(e, inner)
    }
    return(out)
  }
  Map(function(e, ed) {
    if (identical(ed$kind, "tuple")) tuple_out(e, ed)
    else if (is_cx_base(ed$stan_type)) cx_value(e)
    else val(e)
  }, unclass(x), d$elems)
}

flatten_value <- function(v) {
  if (is.list(v)) return(unlist(lapply(v, flatten_value), use.names = FALSE))
  if (is.complex(v)) return(as.vector(rbind(Re(as.vector(v)), Im(as.vector(v)))))
  as.numeric(v)
}

# CmdStan's dotted names (a.1:2.3) in cmdstanr's form (a[1:2,3])
tuple_flat_names <- function(name, d) {
  dotted <- tuple_dotted(name, d)
  vapply(dotted, function(s) {
    p <- regexpr(".", s, fixed = TRUE)
    if (p < 0L) return(s)
    paste0(substr(s, 1L, p - 1L), "[", gsub(".", ",", substr(s, p + 1L, nchar(s)), fixed = TRUE), "]")
  }, "", USE.NAMES = FALSE)
}
tuple_dotted <- function(prefix, d) {
  if (identical(d$kind, "tuple")) {
    if (length(d$adims)) {
      if (prod(d$adims) == 0L) return(character(0))
      # CmdStan writes arrays of tuples in column major order
      g <- as.matrix(expand.grid(lapply(d$adims, seq_len)))
      pre <- paste0(prefix, ".", apply(g, 1L, paste, collapse = "."))
    } else pre <- prefix
    return(unlist(lapply(pre, function(p) unlist(lapply(seq_along(d$elems), function(k)
      tuple_dotted(paste0(p, ":", k), d$elems[[k]]))))))
  }
  dims <- c(d$adims, d$eshape)
  if (is_cx_base(d$stan_type)) {
    ed <- dims[-length(dims)]
    base <- if (length(ed)) {
      if (prod(ed) == 0L) return(character(0))
      g <- do.call(expand.grid, lapply(ed, seq_len))
      paste0(prefix, ".", do.call(paste, c(unname(as.list(g)), sep = ".")))
    } else prefix
    return(as.vector(rbind(paste0(base, ".real"), paste0(base, ".imag"))))
  }
  if (!length(dims)) return(prefix)
  if (prod(dims) == 0L) return(character(0))
  g <- do.call(expand.grid, lapply(dims, seq_len))
  paste0(prefix, ".", do.call(paste, c(unname(as.list(g)), sep = ".")))
}

# ---- library functions returning tuples ----------------------------------------------------

tuple_ret <- function(...) { el <- list(...); function(types) ty_tuple(el) }
CX_MAT <- ty("complex_matrix"); CX_VEC <- ty("complex_vector")
def("eigendecompose_sym", "matrix", tuple_ret(T_mat, T_vec), function(m)
  rt_tuple(list(lib_impl("eigenvectors_sym", list(T_mat))(m), lib_impl("eigenvalues_sym", list(T_mat))(m))))
def("eigendecompose_sym", "complex_matrix", tuple_ret(CX_MAT, CX_VEC), function(m)
  rt_tuple(list(lib_impl("eigenvectors_sym", list(CX_MAT))(m), lib_impl("eigenvalues_sym", list(CX_MAT))(m))))
for (b in c("matrix", "complex_matrix")) local({
  tb <- ty(b)
  def("eigendecompose", b, tuple_ret(CX_MAT, CX_VEC), function(m)
    rt_tuple(list(lib_impl("eigenvectors", list(tb))(m), lib_impl("eigenvalues", list(tb))(m))))
  def("complex_schur_decompose", b, tuple_ret(CX_MAT, CX_MAT), function(m)
    rt_tuple(list(lib_impl("complex_schur_decompose_u", list(tb))(m), lib_impl("complex_schur_decompose_t", list(tb))(m))))
  um <- if (b == "matrix") T_mat else CX_MAT
  def("svd", b, tuple_ret(um, T_vec, um), function(m)
    rt_tuple(list(lib_impl("svd_U", list(tb))(m), lib_impl("singular_values", list(tb))(m), lib_impl("svd_V", list(tb))(m))))
})
def("qr", "matrix", tuple_ret(T_mat, T_mat), function(m)
  rt_tuple(list(lib_impl("qr_Q", list(T_mat))(m), lib_impl("qr_R", list(T_mat))(m))))
def("qr_thin", "matrix", tuple_ret(T_mat, T_mat), function(m)
  rt_tuple(list(lib_impl("qr_thin_Q", list(T_mat))(m), lib_impl("qr_thin_R", list(T_mat))(m))))
def("csr_extract", "matrix", tuple_ret(T_vec, ty("int", 1L), ty("int", 1L)), function(m)
  rt_tuple(list(lib_impl("csr_extract_w", list(T_mat))(m), lib_impl("csr_extract_v", list(T_mat))(m),
                lib_impl("csr_extract_u", list(T_mat))(m))))
