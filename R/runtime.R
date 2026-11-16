# Runtime support for compiled Stan programs: declarations, indexing,
# assignment, loops, the target accumulator, printing and data conversion.

# shape of a declared variable's storage
decl_shape <- function(adims, eshape) {
  d <- c(adims, eshape)
  as.integer(d)
}

# new uninitialized variable (Stan fills reals with NaN)
rt_decl <- function(base, adims, eshape) {
  d <- decl_shape(adims, eshape)
  if (any(d < 0L)) stop_domain(sprintf("declared size is negative (%s)", paste(d, collapse = ", ")))
  n <- prod(d)
  v <- if (base == "int") rep(NA_integer_, n) else rep(NaN, n)
  if (length(d) > 1L) dim(v) <- d
  v
}

storage_dims <- function(x) {
  v <- val(x)
  d <- dim(v)
  if (is.null(d)) length(v) else d
}

# resolve index specifications to per dimension selections
# kinds: "single", "intarray", "range", "lower", "upper", "omni"
index_selection <- function(D, kinds, vals, what = "index") {
  k <- length(D)
  sel <- vector("list", k)
  keep <- logical(k)
  for (j in seq_len(k)) {
    if (j > length(kinds)) {
      sel[[j]] <- seq_len(D[j])
      keep[j] <- TRUE
      next
    }
    kd <- kinds[j]
    v <- vals[[j]]
    s <- switch(kd,
      single = v,
      intarray = v,
      range = if (v[2L] < v[1L]) integer(0) else v[1L]:v[2L],
      lower = if (D[j] < v) integer(0) else v:D[j],
      upper = if (v < 1L) integer(0) else seq_len(v),
      omni = seq_len(D[j]))
    s <- as.integer(s)
    if (length(s) && (any(is.na(s)) || min(s) < 1L || max(s) > D[j])) {
      bad <- s[is.na(s) | s < 1L | s > D[j]][1L]
      stop_domain(sprintf("%s out of range: index %s for a dimension of size %d", what, bad, D[j]))
    }
    sel[[j]] <- s
    keep[j] <- kd != "single"
  }
  list(sel = sel, keep = keep)
}

flat_index <- function(D, sel) {
  flat <- sel[[1L]]
  stride <- 1L
  if (length(D) > 1L) for (j in 2:length(D)) {
    stride <- stride * D[j - 1L]
    flat <- as.vector(outer(flat, (sel[[j]] - 1L) * stride, `+`))
  }
  as.integer(flat)
}

rt_index <- function(x, kinds, vals) {
  if (is_tarr(x)) return(tarr_index(x, kinds, vals))
  D <- storage_dims(x)
  is <- index_selection(D, kinds, vals)
  flat <- flat_index(D, is$sel)
  shape <- lengths(is$sel)[is$keep]
  v <- val(x)
  if (is.integer(v) && !is_ad(x)) {
    out <- v[flat]
    if (length(shape) > 1L) dim(out) <- shape
    return(out)
  }
  ad_gather(x, flat, if (length(shape) > 1L) shape else NULL)
}

# fast paths for x[i] on a one dimensional value and x[i, j] on a two
# dimensional one; anything else goes through rt_index
rt_index_s1 <- function(x, i) {
  if (is.object(x)) {
    if (!inherits(x, "bcad")) return(rt_index(x, "single", list(i)))
    v <- x$v
  } else v <- x
  d <- dim(v)
  if (length(d) > 1L) return(rt_index(x, "single", list(i)))
  n <- length(v)
  if (is.na(i) || i < 1L || i > n)
    stop_domain(sprintf("index out of range: index %s for a dimension of size %d", i, n))
  i <- as.integer(i)
  if (!is.object(x)) return(v[[i]])
  ad_node(v[[i]], list(x), sparse_grad(i, n, d))
}
# backward rule of a single element read, in an environment that does not
# keep the full container alive
sparse_grad <- function(i, n, d) {
  force(i); force(n); force(d)
  function(g) list(ad_sparse(i, g, n, d))
}
rt_index_s2 <- function(x, i, j) {
  if (is.object(x)) {
    if (!inherits(x, "bcad")) return(rt_index(x, c("single", "single"), list(i, j)))
    v <- x$v
  } else v <- x
  d <- dim(v)
  if (length(d) != 2L) return(rt_index(x, c("single", "single"), list(i, j)))
  if (is.na(i) || i < 1L || i > d[1L])
    stop_domain(sprintf("index out of range: index %s for a dimension of size %d", i, d[1L]))
  if (is.na(j) || j < 1L || j > d[2L])
    stop_domain(sprintf("index out of range: index %s for a dimension of size %d", j, d[2L]))
  k <- as.integer(i + (j - 1L) * d[1L])
  if (!is.object(x)) return(v[[k]])
  n <- length(v)
  ad_node(v[[k]], list(x), sparse_grad(k, n, d))
}

rt_assign_index <- function(x, kinds, vals, rhs, name) {
  if (is_tarr(x)) return(tarr_assign(x, kinds, vals, rhs, name))
  D <- storage_dims(x)
  is <- index_selection(D, kinds, vals, what = sprintf("assignment to '%s': index", name))
  flat <- flat_index(D, is$sel)
  shape <- lengths(is$sel)[is$keep]
  vr <- val(rhs)
  rd <- if (is.null(dim(vr))) length(vr) else dim(vr)
  if (!length(shape)) {
    if (length(vr) != 1L) stop_shape(sprintf("assignment to '%s': a single element needs a scalar value", name))
  } else if (!identical(as.integer(rd), as.integer(shape)) && !(length(shape) == 1L && length(vr) == shape)) {
    stop_shape(sprintf("assignment to '%s': size mismatch, the left side selects %s but the right side has size %s",
                       name, paste(shape, collapse = " x "), paste(rd, collapse = " x ")))
  }
  vx <- val(x)
  if (is.integer(vx) && !is_ad(x)) {
    vx[flat] <- as.integer(vr)
    return(vx)
  }
  ad_scatter(x, flat, as_dbl(rhs))
}

rt_assign <- function(old, rhs, name) {
  if (is_tarr(old) || is_tuple_val(old)) return(tuple_assign(old, rhs, name))
  vo <- val(old)
  vr <- val(rhs)
  do <- if (is.null(dim(vo))) length(vo) else dim(vo)
  dr <- if (is.null(dim(vr))) length(vr) else dim(vr)
  if (!identical(as.integer(do), as.integer(dr)))
    stop_shape(sprintf("assignment to '%s': size mismatch, the variable has size %s but the value has size %s",
                       name, paste(do, collapse = " x "), paste(dr, collapse = " x ")))
  if (is.integer(vo) && !is_ad(old)) { storage.mode(vr) <- "integer"; return(vr) }
  as_dbl(rhs)
}

rt_range <- function(lo, hi) if (hi < lo) integer(0) else lo:hi

rt_foreach <- function(x, nd, base) {
  if (is_tarr(x)) return(unclass(x))
  v <- val(x)
  if (is_cx_base(base) && nd == 0L) {
    # the elements of a complex vector or matrix, in column major order
    n <- length(v) %/% 2L
    return(lapply(seq_len(n), function(k) ad_gather(x, c(k, k + n))))
  }
  if (nd > 0L) {
    D <- storage_dims(x)
    n <- D[1L]
    return(lapply(seq_len(n), function(i) rt_index(x, "single", list(i))))
  }
  n <- length(v)
  if (is.integer(v) && !is_ad(x)) return(as.list(v))
  lapply(seq_len(n), function(i) ad_select(x, i))
}

# array expression {a, b, ...}
# scalar_elem: the static element type is int or real (a vector[1] element
# keeps its dimension)
rt_array <- function(items, is_int, scalar_elem = TRUE) {
  shapes <- lapply(items, function(it) {
    v <- val(it)
    if (is.null(dim(v))) (if (scalar_elem) integer(0) else length(v)) else dim(v)
  })
  s1 <- shapes[[1L]]
  for (sh in shapes[-1L]) if (!identical(as.integer(sh), as.integer(s1)))
    stop_shape("array expression elements must all have the same size")
  n <- length(items)
  if (is_int) {
    out <- unlist(lapply(items, val), use.names = FALSE)
    out <- as.integer(out)
    if (length(s1)) {
      p <- layout_perm(n, s1)
      out <- out[p]
      dim(out) <- c(n, s1)
    }
    return(out)
  }
  cc <- ad_c(lapply(items, as_dbl))
  if (!length(s1)) return(cc)
  ad_gather(cc, layout_perm(n, s1), c(n, s1))
}

# elementwise operation on two containers with a strict size check
rt_elementwise <- function(f, a, b, what) {
  da <- dim(val(a)) %||% length(val(a))
  db <- dim(val(b)) %||% length(val(b))
  if (!identical(as.integer(da), as.integer(db)))
    stop_shape(sprintf("%s: size mismatch (%s versus %s)", what, paste(da, collapse = " x "), paste(db, collapse = " x ")))
  f(a, b)
}

# row vector expression [a, b, ...]
rt_rowvec <- function(items, rows) {
  if (!rows) return(ad_c(lapply(items, as_dbl)))
  lens <- vapply(items, nelem, 1L)
  if (any(lens != lens[1L])) stop_shape("matrix expression rows must have the same size")
  ad_transpose(ad_reshape(ad_c(lapply(items, as_dbl)), c(lens[1L], length(items))))
}

# ---- context --------------------------------------------------------------------------

new_ctx <- function(jacobian = TRUE, print = TRUE, propto = FALSE) {
  ctx <- new.env(parent = emptyenv())
  ctx$propto <- propto
  ctx$target <- vector("list", 16L)
  ctx$n_target <- 0L
  ctx$jacobian <- jacobian
  ctx$print <- print
  ctx
}

# a library distribution call whose container arguments (mask) must agree in size
rt_dist_call <- function(fname, f, mask, ...) {
  a <- list(...)
  ls <- vapply(a[mask], function(x) length(val(x)), 1L)
  if (any(ls != ls[1L]))
    stop_shape(sprintf("%s: container arguments have different sizes (%s)", fname, paste(ls, collapse = ", ")))
  f(...)
}

rt_target_add <- function(ctx, v) {
  if (is.null(ctx)) stop_semantic("target can only be incremented in a log density context")
  vv <- val(v)
  if (any(is.nan(vv))) stop_domain("target increment is NaN")
  n <- ctx$n_target + 1L
  # detach the list before writing so R modifies it in place instead of
  # copying it on every increment, and grow it by doubling
  tl <- ctx$target
  ctx$target <- NULL
  if (n > length(tl)) length(tl) <- 2L * length(tl)
  tl[[n]] <- if (length(vv) == 1L) v else ad_sum(as_dbl(v))
  ctx$target <- tl
  ctx$n_target <- n
  invisible(NULL)
}

rt_target_get <- function(ctx) {
  if (!ctx$n_target) return(0)
  ad_sum_list(ctx$target[seq_len(ctx$n_target)])
}

fmt_value <- function(x) {
  v <- val(x)
  if (length(v) == 1L && is.null(dim(v))) return(format_scalar(v))
  d <- dim(v)
  if (is.null(d)) return(paste0("[", paste(vapply(v, format_scalar, ""), collapse = ","), "]"))
  if (length(d) == 2L) {
    rows <- vapply(seq_len(d[1L]), function(i) paste0("[", paste(vapply(v[i, ], format_scalar, ""), collapse = ","), "]"), "")
    return(paste0("[", paste(rows, collapse = ","), "]"))
  }
  paste0("[", paste(vapply(as.vector(v), format_scalar, ""), collapse = ","), "]")
}
# a container value: a single element still prints in brackets
fmt_container <- function(x) {
  v <- val(x)
  if (length(v) == 1L && is.null(dim(v))) return(paste0("[", format_scalar(v), "]"))
  fmt_value(x)
}
format_scalar <- function(v) {
  if (is.integer(v)) return(as.character(v))
  if (is.nan(v)) return("nan")
  if (v == 0 && 1 / v < 0) return("-0")
  if (is.infinite(v)) return(if (v > 0) "inf" else "-inf")
  # C++ stream output: six significant digits (%g)
  trimws(formatC(v, digits = 6L, format = "g"))
}

rt_print <- function(ctx, items) {
  if (!is.null(ctx) && isFALSE(ctx$print)) return(invisible(NULL))
  msg <- paste(vapply(items, function(it) if (is.character(it)) it else fmt_value(it), ""), collapse = "")
  cat(msg, "\n", sep = "")
  invisible(NULL)
}

rt_reject <- function(items, fatal = FALSE) {
  msg <- paste(vapply(items, function(it) if (is.character(it)) it else fmt_value(it), ""), collapse = "")
  stop_reject(msg, fatal = fatal)
}

# a parameter: slice the unconstrained vector, constrain, add the Jacobian
rt_param <- function(u, decl, args, ctx) {
  n <- decl$free_size
  ui <- if (n > 0L) ad_gather(u, decl$offset + seq_len(n)) else numeric(0)
  for (nm in names(args)) {
    a <- args[[nm]]
    if (any(is.na(val(a)))) stop_domain(sprintf("%s bound of parameter '%s' is NaN", nm, decl$name))
  }
  r <- constrain_variable(ui, decl, args)
  if (isTRUE(ctx$jacobian) && !identical(r$lj, 0)) rt_target_add(ctx, r$lj)
  r$x
}

# validate a top level variable at the end of its block
rt_check <- function(x, decl, args, block) {
  # as in stanc3, an unassigned (NaN) transformed parameter is not an error;
  # only the declared constraints are checked
  v <- val(x)
  m <- check_variable(v, decl, lapply(args, val))
  if (!is.null(m)) stop_domain(sprintf("%s variable '%s' %s", block, decl$name, m))
  invisible(NULL)
}

# run a user density with normalization on (called as foo_lpdf)
rt_normalized <- function(expr, ctx) {
  old <- ctx$propto
  ctx$propto <- FALSE
  on.exit(ctx$propto <- old)
  expr
}

# truncated density normalization
rt_truncate <- function(ctx, y, lo, hi, lcdf, lccdf, args, discrete, fam) {
  vy <- as.vector(val(y))
  N <- length(vy)
  vlo <- if (is.null(lo)) NULL else val(lo)
  vhi <- if (is.null(hi)) NULL else val(hi)
  if ((!is.null(vlo) && any(vy < vlo)) || (!is.null(vhi) && any(vy > vhi))) {
    rt_target_add(ctx, -Inf)
    return(invisible(NULL))
  }
  lo_adj <- if (!is.null(lo) && discrete) lo - 1L else lo
  scalar_args <- all(vapply(args, function(a) length(val(a)) == 1L, TRUE))
  if (!scalar_args && (isTRUE(attr(lcdf, "user")) || isTRUE(attr(lccdf, "user"))))
    stop_unsupported(sprintf("truncation of %s with a user defined cumulative function needs scalar parameters", fam))
  one <- function(a) {
    if (!is.null(lo) && !is.null(hi)) {
      if (!is.null(vlo) && !is.null(vhi) && vlo > vhi) stop_domain(sprintf("truncation of %s: lower bound is greater than upper bound", fam))
      l_lo <- do.call(lcdf, c(list(lo_adj), a))
      l_hi <- do.call(lcdf, c(list(hi), a))
      # a lower bound at the edge of the support: F(lo) = 0, and its
      # gradient (often 0 / 0 there) must not enter
      if (all(val(l_lo) == -Inf)) return(l_hi)
      return(ad_log_diff_exp(l_hi, l_lo))
    }
    if (!is.null(lo)) return(do.call(lccdf, c(list(lo_adj), a)))
    do.call(lcdf, c(list(hi), a))
  }
  if (scalar_args) {
    rt_target_add(ctx, -N * one(args))
  } else {
    for (i in seq_len(N)) {
      ai <- lapply(args, function(a) if (length(val(a)) == 1L) a else ad_select(flat(a), i))
      rt_target_add(ctx, -one(ai))
    }
  }
  invisible(NULL)
}

# ---- profiling -----------------------------------------------------------------
# profile("name") { } blocks accumulate their run time and number of runs,
# split by whether automatic differentiation was active (CmdStan's
# autodiff_calls and no_autodiff_calls)
.prof <- new.env(parent = emptyenv())
.prof_active <- new.env(parent = emptyenv())
prof_new <- function() c(total = 0, ad_time = 0, reverse = 0, chain_stack = 0, autodiff_calls = 0, no_autodiff_calls = 0)
rt_profile_start <- function(name) {
  # as Stan Math, a profile cannot be entered again while it is running
  if (isTRUE(.prof_active[[name]])) stop_domain(sprintf("Profile '%s' already started!", name))
  assign(name, TRUE, envir = .prof_active)
  c(t = proc.time()[["elapsed"]], n = .ad$n)
}
rt_profile_end <- function(name, t0) {
  assign(name, FALSE, envir = .prof_active)
  dt <- proc.time()[["elapsed"]] - t0[["t"]]
  r <- .prof[[name]] %||% prof_new()
  ad <- !is.null(.ad$tape)
  r[["total"]] <- r[["total"]] + dt
  if (ad) {
    r[["ad_time"]] <- r[["ad_time"]] + dt; r[["autodiff_calls"]] <- r[["autodiff_calls"]] + 1
    # nodes made in the block (Stan's chain_stack); their reverse pass is
    # timed by ad_backward()
    k <- .ad$n - t0[["n"]]
    r[["chain_stack"]] <- r[["chain_stack"]] + k
    if (k > 0L) .ad$marks[[length(.ad$marks) + 1L]] <- list(name = name, from = as.integer(t0[["n"]]) + 1L, to = .ad$n)
  } else r[["no_autodiff_calls"]] <- r[["no_autodiff_calls"]] + 1
  assign(name, r, envir = .prof)
  invisible()
}
profile_add_reverse <- function(name, dt) {
  r <- .prof[[name]] %||% prof_new()
  r[["reverse"]] <- r[["reverse"]] + dt
  r[["total"]] <- r[["total"]] + dt
  assign(name, r, envir = .prof)
}
profile_reset <- function() {
  rm(list = ls(.prof, all.names = TRUE), envir = .prof)
  rm(list = ls(.prof_active, all.names = TRUE), envir = .prof_active)
}
profile_table <- function() {
  nm <- sort(ls(.prof, all.names = TRUE))
  if (!length(nm)) return(NULL)
  v <- do.call(rbind, lapply(nm, function(n) .prof[[n]]))
  data.frame(name = nm, thread_id = 1L, total_time = signif(v[, "total"], 6), forward_time = signif(v[, "ad_time"], 6),
             reverse_time = signif(v[, "reverse"], 6), chain_stack = as.integer(v[, "chain_stack"]), no_chain_stack = 0L, autodiff_calls = as.integer(v[, "autodiff_calls"]),
             no_autodiff_calls = as.integer(v[, "no_autodiff_calls"]), row.names = NULL)
}

# a parameter dependent value that lost its derivative node (floor(), an
# assigned constant) as a constant node, so propto keeps its density terms
rt_keep <- function(x) {
  if (inherits(x, "bcad") || !is.numeric(x) || is.null(.ad$tape)) return(x)
  ad_leaf(x)
}

