# Reverse mode automatic differentiation.
#
# A tape records every operation applied to differentiable values while a
# function runs. Each recorded node holds an array value. Its backward rule
# receives the adjoint of the node (an array of the same shape) and adds
# contributions to the adjoints of its parents. Running the rules in reverse
# order of creation gives the gradient of a scalar output with respect to the
# leaves.
#
# Values are ordinary R numeric vectors or arrays. A differentiable value is a
# small list of class "bcad" holding the value `v` and the node index `id` on
# the active tape. Arithmetic, the Math and Summary group generics and a set of
# explicit functions below are defined for "bcad" objects so that numerical
# code can be written in plain R and differentiated without change.

.ad <- new.env(parent = emptyenv())
# The tape lives directly in fields of `.ad` (not in an environment passed
# as an argument): R modifies list elements of a package level environment
# in place, while writing through an argument copies the whole list.
#   tape  TRUE while recording, NULL otherwise
#   n     number of nodes
#   gr    backward rule of each node: function(g) returning parent adjoints
#   ids   parent node ids of each node
#   adj   adjoint of each node
.ad$tape <- NULL
.ad$n <- 0L

tape_start <- function() {
  # an environment, not a list: tape_restore() moves the outer lists back and
  # drops these bindings, so the outer tape is not shared and not copied by
  # its next node
  old <- new.env(parent = emptyenv())
  old$tape <- .ad$tape; old$n <- .ad$n; old$cap <- .ad$cap
  .ad$check <- isTRUE(getOption("bayescompR.check_adjoints"))
  old$gr <- .ad$gr; .ad$gr <- NULL
  old$ids <- .ad$ids; .ad$ids <- NULL
  old$adj <- .ad$adj; .ad$adj <- NULL
  old$marks <- .ad$marks; .ad$marks <- NULL
  .ad$tape <- TRUE
  .ad$n <- 0L
  .ad$gr <- vector("list", 64L)
  .ad$ids <- vector("list", 64L)
  .ad$adj <- vector("list", 64L)
  .ad$cap <- 64L
  old
}
tape_restore <- function(old) {
  .ad$tape <- old$tape; .ad$n <- old$n; .ad$cap <- old$cap
  .ad$gr <- old$gr; old$gr <- NULL
  .ad$ids <- old$ids; old$ids <- NULL
  .ad$adj <- old$adj; old$adj <- NULL
  .ad$marks <- old$marks
}

# never call length(.ad$gr): reading the list marks it shared and the next
# element assignment would copy the whole tape
tape_grow <- function() {
  len <- 2L * .ad$cap
  .ad$cap <- len
  length(.ad$gr) <- len
  length(.ad$ids) <- len
  length(.ad$adj) <- len
}

is_ad <- function(x) inherits(x, "bcad")

any_ad <- function(...) {
  for (a in list(...)) if (inherits(a, "bcad")) return(TRUE)
  FALSE
}

# value of a possibly differentiable object
val <- function(x) if (inherits(x, "bcad")) x$v else x

ad_leaf <- function(v) {
  if (is.null(.ad$tape)) stop("internal error: no active tape", call. = FALSE)
  storage.mode(v) <- "double"
  n <- .ad$n + 1L
  if (n > .ad$cap) tape_grow()
  .ad$n <- n
  x <- list(v = v, id = n)
  class(x) <- "bcad"
  x
}

# Create a node with value `v`. `parents` is a list of operands (bcad or not);
# `grads` is a function(g) returning a list, one adjoint contribution per
# parent (NULL for non differentiable parents).
ad_node <- function(v, parents, grads) {
  if (is.null(.ad$tape)) return(v)
  np <- length(parents)
  ids <- integer(np)
  for (k in seq_len(np)) {
    p <- parents[[k]]
    if (inherits(p, "bcad")) ids[k] <- p$id
  }
  if (all(ids == 0L)) return(v)
  if (isTRUE(.ad$check)) grads <- adj_checked(grads, parents, ids)
  n <- .ad$n + 1L
  if (n > .ad$cap) tape_grow()
  .ad$n <- n
  # detach, modify, reattach: keeps each list unshared so the element
  # assignment happens in place instead of copying the whole tape
  l <- .ad$gr; .ad$gr <- NULL; l[[n]] <- grads; .ad$gr <- l
  l <- .ad$ids; .ad$ids <- NULL; l[[n]] <- ids; .ad$ids <- l
  x <- list(v = v, id = n)
  class(x) <- "bcad"
  x
}

# options(bayescompR.check_adjoints = TRUE), read when a tape starts: every
# backward rule is checked to return, for each differentiable parent, an
# numeric adjoint of its length (or a scalar); NaN is allowed, since a
# trajectory far in the tails can overflow. Used by the test suite to catch shape bugs in new rules.
adj_checked <- function(grads, parents, ids) {
  force(grads); force(ids)
  lens <- vapply(parents, function(p) length(val(p)), 0L)
  parents <- NULL
  function(g) {
    cs <- grads(g)
    for (k in which(ids > 0L)) {
      ck <- cs[[k]]
      if (is.null(ck)) next
      n <- if (inherits(ck, "bcsp")) ck$n else length(ck)
      if (n != lens[k] && !(n == 1L && !inherits(ck, "bcsp")))
        stop(sprintf("internal error: adjoint of length %d for an operand of length %d", n, lens[k]), call. = FALSE)
      v <- if (inherits(ck, "bcsp")) ck$g else ck
      if (!is.numeric(v)) stop("internal error: non-numeric adjoint", call. = FALSE)
    }
    cs
  }
}

ad_backward <- function(out) {
  if (!inherits(out, "bcad")) return(invisible(NULL))
  # the adjoint list is a local variable here, so updates are in place
  adj <- .ad$adj
  .ad$adj <- NULL
  grl <- .ad$gr
  idl <- .ad$ids
  adj[[out$id]] <- 1
  # profile blocks on this tape: node ranges whose reverse pass is timed
  marks <- .ad$marks
  if (length(marks)) {
    from <- vapply(marks, `[[`, 0L, "from"); to <- vapply(marks, `[[`, 0L, "to")
    t_in <- rep(NA_real_, length(marks)); rev_t <- numeric(length(marks))
  }
  for (i in seq.int(out$id, 1L)) {
    if (length(marks)) {
      now <- NULL
      for (k in which(to == i)) { now <- now %||% proc.time()[["elapsed"]]; t_in[k] <- now }
      for (k in which(from == i + 1L & !is.na(t_in))) { now <- now %||% proc.time()[["elapsed"]]; rev_t[k] <- rev_t[k] + now - t_in[k]; t_in[k] <- NA }
    }
    g <- adj[[i]]
    # a zero adjoint contributes nothing; skipping it also keeps 0 * Inf
    # partials (a log at a boundary) from turning the gradient into NaN
    if (is.null(g) || (!anyNA(g) && all(g == 0))) next
    gr <- grl[[i]]
    if (is.null(gr)) next
    cs <- gr(g)
    # only leaf adjoints are read after the sweep; dropping the others keeps
    # memory linear when many nodes hold full size adjoints
    adj[i] <- list(NULL)
    ids <- idl[[i]]
    for (k in seq_along(ids)) {
      id <- ids[k]
      if (id > 0L) {
        ck <- cs[[k]]
        if (!is.null(ck)) {
          cur <- adj[[id]]
          if (inherits(ck, "bcsp")) {
            # sparse contribution (a selection): add into the parent's
            # adjoint in place instead of building a dense array
            if (is.null(cur)) {
              cur <- numeric(ck$n)
              if (!is.null(ck$d)) dim(cur) <- ck$d
            } else {
              adj[id] <- list(NULL)
            }
            cur[ck$i] <- cur[ck$i] + ck$g
            adj[[id]] <- cur
          } else {
            adj[[id]] <- if (is.null(cur)) ck else cur + ck
          }
        }
      }
    }
  }
  .ad$adj <- adj
  if (length(marks)) {
    now <- proc.time()[["elapsed"]]
    for (k in which(!is.na(t_in))) rev_t[k] <- rev_t[k] + now - t_in[k]
    for (k in seq_along(marks)) profile_add_reverse(marks[[k]]$name, rev_t[k])
  }
  invisible(NULL)
}

# Run f on a differentiable copy of x and return value and gradient.
ad_value_grad <- function(f, x, ...) {
  old <- tape_start()
  on.exit(tape_restore(old), add = TRUE)
  xa <- ad_leaf(as.numeric(x))
  out <- f(xa, ...)
  value <- as.numeric(val(out))
  if (length(value) != 1L)
    stop("the function must return a single value", call. = FALSE)
  if (!inherits(out, "bcad")) return(list(value = value, grad = numeric(length(x))))
  ad_backward(out)
  g <- .ad$adj[[xa$id]]
  g <- if (is.null(g)) numeric(length(x)) else as.numeric(g)
  if (length(g) == 1L && length(x) > 1L) g <- rep(g, length(x))
  list(value = value, grad = g)
}

# ---------------------------------------------------------------------------
# Broadcasting helpers

# Sum an adjoint down to the shape of an operand that was broadcast.
unbroadcast <- function(g, x) {
  vx <- val(x)
  lx <- length(vx)
  if (lx == 1L) return(sum(g))
  if (length(g) == 1L) g <- rep(g, lx)
  d <- dim(vx)
  if (is.null(d)) {
    g <- as.vector(g)
  } else {
    dim(g) <- d
  }
  g
}

# result shape for elementwise operations: the larger operand decides
bshape <- function(a, b) {
  va <- val(a)
  vb <- val(b)
  la <- length(va)
  lb <- length(vb)
  if (la != lb && la != 1L && lb != 1L)
    stop(sprintf("internal shape mismatch in elementwise operation: %d versus %d", la, lb),
         call. = FALSE)
  invisible(NULL)
}

# g * d where a zero adjoint gives zero even when d is not finite (an
# element of a container that does not reach the result, such as the off
# diagonal of sqrt(D) at 0, must not turn the gradient into NaN)
gmul <- function(g, d) {
  r <- g * d
  if (anyNA(r)) r[rep_len(g == 0, length(r))] <- 0
  r
}

# elementwise unary op with derivative function df(x, y) giving dy/dx
ad_unary <- function(x, f, df) {
  vx <- val(x)
  y <- f(vx)
  if (!inherits(x, "bcad")) return(y)
  ad_node(y, list(x), function(g) list(gmul(g, df(vx, y))))
}

# elementwise binary op; da(a, b, y), db(a, b, y) are partial derivatives
ad_binary <- function(a, b, f, da, db) {
  va <- val(a)
  vb <- val(b)
  bshape(va, vb)
  y <- f(va, vb)
  if (!inherits(a, "bcad") && !inherits(b, "bcad")) return(y)
  ad_node(y, list(a, b), function(g) {
    list(
      if (inherits(a, "bcad")) unbroadcast(gmul(g, da(va, vb, y)), va),
      if (inherits(b, "bcad")) unbroadcast(gmul(g, db(va, vb, y)), vb)
    )
  })
}

# ---------------------------------------------------------------------------
# S3 methods so that numerical code reads like ordinary R

#' @export
Ops.bcad <- function(e1, e2) {
  if (missing(e2)) {
    switch(.Generic,
      "-" = return(ad_unary(e1, function(x) -x, function(x, y) -1)),
      "+" = return(e1),
      "!" = return(!val(e1)),
      stop("unsupported unary operator ", .Generic, call. = FALSE))
  }
  # scalar arithmetic without the generic machinery (loops over scalars)
  if (.Generic %in% c("+", "-", "*", "/")) {
    ia <- inherits(e1, "bcad"); ib <- inherits(e2, "bcad")
    va <- if (ia) e1$v else e1
    vb <- if (ib) e2$v else e2
    if (length(va) == 1L && length(vb) == 1L && is.null(dim(va)) && is.null(dim(vb)))
      return(switch(.Generic,
        "+" = ad_node(va + vb, list(e1, e2), function(g) list(if (ia) g, if (ib) g)),
        "-" = ad_node(va - vb, list(e1, e2), function(g) list(if (ia) g, if (ib) -g)),
        "*" = ad_node(va * vb, list(e1, e2), function(g) list(if (ia) g * vb, if (ib) g * va)),
        "/" = { y <- va / vb; ad_node(y, list(e1, e2), function(g) list(if (ia) g / vb, if (ib) -g * y / vb)) }))
  }
  switch(.Generic,
    "+" = ad_binary(e1, e2, `+`, function(a, b, y) 1, function(a, b, y) 1),
    "-" = ad_binary(e1, e2, `-`, function(a, b, y) 1, function(a, b, y) -1),
    "*" = ad_binary(e1, e2, `*`, function(a, b, y) b, function(a, b, y) a),
    "/" = ad_binary(e1, e2, `/`, function(a, b, y) 1 / b, function(a, b, y) -y / b),
    "^" = ad_pow(e1, e2),
    "==" = , "!=" = , "<" = , "<=" = , ">" = , ">=" = ,
    "&" = , "|" = get(.Generic)(val(e1), val(e2)),
    stop("unsupported operator ", .Generic, " for differentiable values", call. = FALSE))
}

ad_pow <- function(a, b) {
  ad_binary(a, b, function(x, p) x^p,
    function(x, p, y) { d <- p * x^(p - 1); d[(p + 0 * x) == 0] <- 0; d },
    function(x, p, y) { d <- y * log(x); d[(x + 0 * p) == 0] <- 0; d })
}

#' @export
Math.bcad <- function(x, ...) {
  switch(.Generic,
    exp = ad_unary(x, exp, function(x, y) y),
    log = ad_unary(x, log, function(x, y) 1 / x),
    sqrt = ad_unary(x, sqrt, function(x, y) 0.5 / y),
    abs = ad_unary(x, abs, function(x, y) sign(x)),
    expm1 = ad_unary(x, expm1, function(x, y) y + 1),
    log1p = ad_unary(x, log1p, function(x, y) 1 / (1 + x)),
    log2 = ad_unary(x, log2, function(x, y) 1 / (x * log(2))),
    log10 = ad_unary(x, log10, function(x, y) 1 / (x * log(10))),
    sin = ad_unary(x, sin, function(x, y) cos(x)),
    cos = ad_unary(x, cos, function(x, y) -sin(x)),
    tan = ad_unary(x, tan, function(x, y) 1 + y^2),
    asin = ad_unary(x, asin, function(x, y) 1 / sqrt(1 - x^2)),
    acos = ad_unary(x, acos, function(x, y) -1 / sqrt(1 - x^2)),
    atan = ad_unary(x, atan, function(x, y) 1 / (1 + x^2)),
    sinh = ad_unary(x, sinh, function(x, y) cosh(x)),
    cosh = ad_unary(x, cosh, function(x, y) sinh(x)),
    tanh = ad_unary(x, tanh, function(x, y) 1 - y^2),
    asinh = ad_unary(x, asinh, function(x, y) 1 / sqrt(x^2 + 1)),
    acosh = ad_unary(x, acosh, function(x, y) 1 / sqrt(x^2 - 1)),
    atanh = ad_unary(x, atanh, function(x, y) 1 / (1 - x^2)),
    lgamma = ad_unary(x, lgamma, function(x, y) digamma(x)),
    gamma = ad_unary(x, gamma, function(x, y) y * digamma(x)),
    digamma = ad_unary(x, digamma, function(x, y) trigamma(x)),
    trigamma = ad_unary(x, trigamma, function(x, y) psigamma(x, 2L)),
    floor = , ceiling = , round = , trunc = , sign = get(.Generic)(val(x)),
    cumsum = ad_cumsum(x),
    stop("unsupported function ", .Generic, " for differentiable values", call. = FALSE))
}

#' @export
Summary.bcad <- function(..., na.rm = FALSE) {
  args <- list(...)
  if (length(args) != 1L) {
    if (.Generic == "sum") return(Reduce(function(a, b) ad_sum(a) + ad_sum(b), args))
    args <- list(ad_c(args))
  }
  x <- args[[1L]]
  switch(.Generic,
    sum = ad_sum(x),
    prod = ad_prod(x),
    max = ad_select(x, which_ext(val(x), which.max)),
    min = ad_select(x, which_ext(val(x), which.min)),
    stop("unsupported summary ", .Generic, " for differentiable values", call. = FALSE))
}

#' @export
length.bcad <- function(x) length(x$v)

#' @export
dim.bcad <- function(x) dim(x$v)

#' @export
print.bcad <- function(x, ...) {
  cat("<differentiable value>\n")
  print(x$v)
  invisible(x)
}

#' @export
as.numeric.bcad <- function(x, ...) as.numeric(x$v)

#' @export
as.double.bcad <- function(x, ...) as.double(x$v)

# ---------------------------------------------------------------------------
# Reductions and structural operations

ad_sum <- function(x) {
  if (!inherits(x, "bcad")) return(sum(x))
  vx <- x$v
  ad_node(sum(vx), list(x), function(g) {
    r <- array(g, dim = if (is.null(dim(vx))) length(vx) else dim(vx))
    if (is.null(dim(vx))) r <- as.vector(r)
    list(r)
  })
}

ad_prod <- function(x) {
  vx <- val(x)
  y <- prod(vx)
  if (!inherits(x, "bcad")) return(y)
  ad_node(y, list(x), function(g) {
    # product of the other elements from prefix and suffix products: O(n), exact with zeros
    n <- length(vx)
    d <- c(1, cumprod(vx[-n])) * rev(c(1, cumprod(rev(vx)[-n])))
    if (n == 0L) d <- numeric(0)
    r <- g * d
    if (!is.null(dim(vx))) dim(r) <- dim(vx)
    list(r)
  })
}

ad_cumsum <- function(x) {
  vx <- val(x)
  y <- cumsum(vx)
  if (!inherits(x, "bcad")) return(y)
  ad_node(y, list(x), function(g) list(rev(cumsum(rev(g)))))
}

# select element i (flat index) of x
ad_select <- function(x, i) {
  vx <- val(x)
  y <- vx[i]
  if (!inherits(x, "bcad")) return(y)
  n <- length(vx)
  d <- dim(vx)
  ad_node(y, list(x), function(g) list(ad_sparse(i, g, n, d)))
}

# sparse adjoint contribution g at the distinct flat indices i of a parent
# of length n and dimensions d, accumulated in place by ad_backward
ad_sparse <- function(i, g, n, d) {
  x <- list(i = i, g = as.vector(g), n = n, d = d)
  class(x) <- "bcsp"
  x
}

# Gather: y = x[idx] for a flat integer index vector, reshaped to `shape`.
ad_gather <- function(x, idx, shape = NULL) {
  vx <- val(x)
  y <- vx[idx]
  if (!is.null(shape) && length(shape) > 1L) dim(y) <- shape
  if (!inherits(x, "bcad")) return(y)
  n <- length(vx)
  d <- dim(vx)
  if (4L * length(idx) < n && !anyDuplicated(idx))
    return(ad_node(y, list(x), function(g) list(ad_sparse(idx, g, n, d))))
  ad_node(y, list(x), function(g) {
    r <- numeric(n)
    if (anyDuplicated(idx)) {
      r <- r + tabulate_weighted(idx, as.vector(g), n)
    } else {
      r[idx] <- as.vector(g)
    }
    if (!is.null(d)) dim(r) <- d
    list(r)
  })
}

tabulate_weighted <- function(idx, w, n) {
  r <- numeric(n)
  s <- rowsum(w, idx, reorder = FALSE)
  r[as.integer(rownames(s))] <- s[, 1L]
  r
}

# Scatter: replace x[idx] by y (flat indices). Functional update.
ad_scatter <- function(x, idx, y) {
  vx <- val(x)
  vy <- val(y)
  out <- vx
  storage.mode(out) <- "double"
  out[idx] <- vy
  ix <- inherits(x, "bcad"); iy <- inherits(y, "bcad")
  if (!ix && !iy) return(out)
  ad_node(out, list(x, y), scatter_grad(ix, iy, length(vx), idx, length(vy), dim(vy)))
}
# the backward rule lives in its own small environment: it must not keep the
# full vectors of x and of the result alive, or a loop of element writes
# holds one copy per iteration
scatter_grad <- function(ix, iy, n, idx, ly, dy) {
  # forced now: an unevaluated argument would keep the caller's frame alive
  force(ix); force(iy); force(n); force(idx); force(ly); force(dy)
  # with repeated indices only the last write survives, so only it gets
  # the adjoint
  overwritten <- if (anyDuplicated(idx) > 0L) duplicated(idx, fromLast = TRUE)
  function(g) {
    gx <- NULL
    gy <- NULL
    if (ix) {
      gx <- g
      if (length(gx) == 1L && n > 1L) gx <- rep(gx, n)
      gx[idx] <- 0
    }
    if (iy) {
      gy <- as.vector(g)[idx]
      if (!is.null(overwritten)) gy[overwritten] <- 0
      if (ly == 1L) gy <- sum(gy) else if (!is.null(dy)) dim(gy) <- dy
    }
    list(gx, gy)
  }
}

# Concatenate a list of operands into one flat vector.
ad_c <- function(xs) {
  vals <- lapply(xs, function(a) as.vector(val(a)))
  y <- unlist(vals, use.names = FALSE)
  storage.mode(y) <- "double"
  if (!any(vapply(xs, inherits, TRUE, what = "bcad"))) return(y)
  lens <- lengths(vals)
  ends <- cumsum(lens)
  starts <- ends - lens + 1L
  ad_node(y, xs, function(g) {
    lapply(seq_along(xs), function(k) {
      if (!inherits(xs[[k]], "bcad")) return(NULL)
      r <- g[starts[k]:ends[k]]
      if (lens[k] == 0L) r <- numeric(0)
      d <- dim(val(xs[[k]]))
      if (!is.null(d)) dim(r) <- d
      r
    })
  })
}

# Reshape without changing the column major element order.
ad_reshape <- function(x, shape) {
  vx <- val(x)
  if (is.null(dim(vx)) && (is.null(shape) || length(shape) <= 1L)) return(x)
  if (!is.null(shape) && length(shape) > 1L && identical(as.integer(dim(vx)), as.integer(shape))) return(x)
  y <- vx
  if (is.null(shape) || length(shape) <= 1L) {
    y <- as.vector(y)
  } else {
    dim(y) <- shape
  }
  if (!inherits(x, "bcad")) return(y)
  d <- dim(vx)
  ad_node(y, list(x), function(g) {
    if (is.null(d)) g <- as.vector(g) else dim(g) <- d
    list(g)
  })
}

ad_transpose <- function(x) {
  vx <- val(x)
  if (is.null(dim(vx)) || length(dim(vx)) != 2L) stop("internal error: transpose needs a matrix")
  y <- t(vx)
  if (!inherits(x, "bcad")) return(y)
  ad_node(y, list(x), function(g) list(t(g)))
}

# Matrix product of two R matrices (values carry dims).
ad_matmul <- function(a, b) {
  va <- val(a)
  vb <- val(b)
  ka <- if (is.null(dim(va))) length(va) else ncol(va)
  kb <- if (is.null(dim(vb))) length(vb) else nrow(vb)
  if (!is.null(dim(va)) && !is.null(dim(vb)) && ka != kb)
    stop_shape(sprintf("matrix product: %d columns times %d rows", ka, kb))
  y <- va %*% vb
  if (!inherits(a, "bcad") && !inherits(b, "bcad")) return(y)
  ad_node(y, list(a, b), function(g) {
    list(
      if (inherits(a, "bcad")) tcrossprod(g, vb),
      if (inherits(b, "bcad")) crossprod(va, g)
    )
  })
}

# Elementwise ifelse on a fixed logical condition (no derivative for cond).
ad_ifelse <- function(cond, a, b) {
  va <- val(a)
  vb <- val(b)
  n <- max(length(cond), length(va), length(vb))
  y <- ifelse(rep_len(cond, n), rep_len(va, n), rep_len(vb, n))
  ad_node(y, list(a, b), function(g) {
    cc <- rep_len(cond, n)
    list(
      if (inherits(a, "bcad")) unbroadcast(g * cc, va),
      if (inherits(b, "bcad")) unbroadcast(g * !cc, vb)
    )
  })
}

# ---------------------------------------------------------------------------
# Numerically stable scalar functions (vectorized), used across the library.

log1p_exp_v <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
log1m_exp_v <- function(x) {
  # log(1 - exp(x)) for x < 0
  ifelse(x > -0.693147, log(-expm1(x)), log1p(-exp(x)))
}
inv_logit_v <- function(x) stats::plogis(x)
log_inv_logit_v <- function(x) stats::plogis(x, log.p = TRUE)
log1m_inv_logit_v <- function(x) stats::plogis(x, lower.tail = FALSE, log.p = TRUE)

ad_log1p_exp <- function(x) ad_unary(x, log1p_exp_v, function(x, y) inv_logit_v(x))
ad_log1m_exp <- function(x) ad_unary(x, log1m_exp_v, function(x, y) -1 / expm1(-x))
ad_inv_logit <- function(x) ad_unary(x, inv_logit_v, function(x, y) y * (1 - y))
ad_log_inv_logit <- function(x) ad_unary(x, log_inv_logit_v, function(x, y) inv_logit_v(-x))
ad_log1m_inv_logit <- function(x) ad_unary(x, log1m_inv_logit_v, function(x, y) -inv_logit_v(x))
ad_logit <- function(x) ad_unary(x, stats::qlogis, function(x, y) 1 / (x * (1 - x)))
ad_square <- function(x) ad_unary(x, function(x) x * x, function(x, y) 2 * x)
ad_log1m <- function(x) ad_unary(x, function(x) log1p(-x), function(x, y) -1 / (1 - x))

ad_log_sum_exp <- function(x) {
  vx <- val(x)
  if (length(vx) == 0L) return(-Inf)
  m <- max(vx)
  y <- if (!is.finite(m)) m else m + log(sum(exp(vx - m)))
  if (!inherits(x, "bcad")) return(y)
  ad_node(y, list(x), function(g) {
    w <- if (is.finite(y)) exp(vx - y) else as.numeric(vx == m) / sum(vx == m)
    if (!is.null(dim(vx))) dim(w) <- dim(vx)
    list(g * w)
  })
}

ad_log_sum_exp2 <- function(a, b) {
  ad_binary(a, b,
    function(x, z) {
      m <- pmax(x, z)
      ifelse(is.infinite(m) & m < 0, -Inf, m + log1p(exp(-abs(x - z))))
    },
    function(x, z, y) ifelse(is.finite(y), exp(x - y), 0.5),
    function(x, z, y) ifelse(is.finite(y), exp(z - y), 0.5))
}

ad_log_diff_exp <- function(a, b) {
  ad_binary(a, b,
    function(x, z) {
      r <- x + log1m_exp_v(z - x)
      # log(exp(x) - 0) = x; broadcast without 0 * Inf
      zi <- rep_len(z, length(r)) == -Inf
      r[zi] <- rep_len(x, length(r))[zi]
      r
    },
    function(x, z, y) -1 / expm1(z - x),
    function(x, z, y) -1 / expm1(x - z))
}

ad_Phi <- function(x) ad_unary(x, stats::pnorm, function(x, y) stats::dnorm(x))
ad_log_Phi <- function(x) {
  ad_unary(x, function(x) stats::pnorm(x, log.p = TRUE),
           function(x, y) exp(stats::dnorm(x, log = TRUE) - y))
}
ad_log1m_Phi <- function(x) {
  ad_unary(x, function(x) stats::pnorm(x, lower.tail = FALSE, log.p = TRUE),
           function(x, y) -exp(stats::dnorm(x, log = TRUE) - y))
}
ad_inv_Phi <- function(p) {
  ad_unary(p, stats::qnorm, function(x, y) 1 / stats::dnorm(y))
}
# erf(x) = sign(x) P(1/2, x^2): accurate near 0, where 2 Phi(x sqrt 2) - 1 cancels
ad_erf <- function(x) ad_unary(x, function(x) sign(x) * stats::pgamma(x * x, 0.5),
                               function(x, y) 2 / sqrt(pi) * exp(-x^2))
ad_erfc <- function(x) ad_unary(x, function(x) 2 * stats::pnorm(-x * sqrt(2)),
                                function(x, y) -2 / sqrt(pi) * exp(-x^2))
ad_lbeta <- function(a, b) {
  ad_binary(a, b, lbeta,
    function(a, b, y) digamma(a) - digamma(a + b),
    function(a, b, y) digamma(b) - digamma(a + b))
}
# lchoose rounds k, so real arguments use Stan's -log1p(n) - lbeta(n - k + 1, k + 1)
lchoose_real <- function(n, k) ifelse(k == round(k), lchoose(n, round(k)), -log1p(n) - lbeta(n - k + 1, k + 1))
ad_lchoose <- function(n, k) {
  ad_binary(n, k, lchoose_real,
    function(n, k, y) digamma(n + 1) - digamma(n - k + 1),
    function(n, k, y) digamma(n - k + 1) - digamma(k + 1))
}

# Numerical derivative fallback for functions without closed form partials.
# Central differences with a step scaled to the argument; used only for
# special function arguments (for example the shape of an incomplete gamma
# function) and documented as such in the function support matrix.
num_partial <- function(f, x, ...) {
  h <- 1e-6 * pmax(abs(x), 1)
  # a positive value within one step of zero (a shape near 0) is stepped on
  # the log scale, so the lower point stays in the domain
  small <- x > 0 & x <= h
  lo <- ifelse(small, x * exp(-1e-6), x - h)
  hi <- ifelse(small, x * exp(1e-6), x + h)
  (f(hi, ...) - f(lo, ...)) / (hi - lo)
}

# sum of a list of scalars as a single node
ad_sum_list <- function(xs) {
  y <- 0
  isad <- logical(length(xs))
  for (k in seq_along(xs)) {
    a <- xs[[k]]
    if (inherits(a, "bcad")) { isad[k] <- TRUE; y <- y + sum(a$v) } else y <- y + sum(a)
  }
  if (!any(isad)) return(y)
  ops <- xs[isad]
  ad_node(y, ops, function(g) lapply(ops, function(o) {
    vx <- o$v
    if (length(vx) == 1L) g else array(g, dim = if (is.null(dim(vx))) length(vx) else dim(vx))[seq_along(vx)] + 0 * vx
  }))
}
