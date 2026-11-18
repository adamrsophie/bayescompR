# Loop vectorization.
#
# A loop such as
#
#   for (n in 1:N) y[n] ~ normal(alpha + X[n] * beta + u[g[n]], sigma);
#
# builds one tape node per iteration and statement. When every iteration
# only reads and writes its own element, the loop is the same as the
# statement on slices,
#
#   y[1:N] ~ normal(alpha + X[1:N] * beta + to_vector(u[g[1:N]]), sigma);
#
# which the vectorized densities evaluate with a handful of nodes. The pass
# works on the parse tree, using the types from a first check to decide what
# is safe, and then checks the rewritten program again. A loop that cannot be
# rewritten, or whose rewrite does not check, stays a loop.
#
# A loop qualifies when its bounds do not depend on the loop, its body is a
# sequence of sampling statements, `target +=` of one density call, and
# assignments `v[n] = e`, `v[n] += e` or `v[n] -= e` to a vector `v`, and
# every expression is built from:
#   * expressions that do not use the loop variable,
#   * `x[n]`, `x[n + c]`, `x[n - c]` (c an integer that does not depend on
#     the loop) and `x[idx[n]]` with `x` a vector, row vector or one
#     dimensional array, `X[n] * b` with `X` a matrix and `b` a vector, and
#     `X[n, k]`,
#   * arithmetic (`+ - * / ^`, unary minus) on scalars,
#   * elementwise functions of one scalar (exp, log, inv_logit, ...),
#   * `dot_product(X[n], b)`.
# Scalar locals declared with a value in the body are replaced by that
# value; `target +=` of any scalar adds the sum of the vector.
# A variable assigned in the loop may only be read as `v[n]` in the loop.

vec_unary <- c("exp", "log", "log1p", "expm1", "log2", "log10", "exp2", "inv_logit", "logit", "log_inv_logit",
               "log1m_inv_logit", "inv_cloglog", "log1m", "log1m_exp", "log1p_exp", "square", "sqrt", "cbrt", "inv",
               "inv_sqrt", "inv_square", "Phi", "Phi_approx", "inv_Phi", "fabs", "abs", "lgamma", "digamma", "tgamma",
               "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh", "erf", "erfc",
               "std_normal_log_qf", "std_normal_qf")

vectorize_enabled <- function() !isFALSE(getOption("bayescompR.vectorize"))

# rewrite qualifying loops of the model and transformed parameters blocks;
# returns the checked program
vectorize_program <- function(ast, ir, external = character()) {
  rewrite <- function(allow) {
    st <- new.env(parent = emptyenv())
    st$k <- 0L; st$done <- 0L; st$allow <- allow
    out <- ast
    for (bn in c("transformed parameters", "model")) {
      if (is.null(ast$blocks[[bn]])) next
      out$blocks[[bn]]$body <- vec_stmts(ast$blocks[[bn]]$body, ir$blocks[[bn]]$body, st)
    }
    list(ast = out, found = st$k, done = st$done)
  }
  all <- rewrite(NULL)
  if (!all$done) return(ir)
  ok <- tryCatch(stan_check(all$ast, external = external), bayescomp_error = function(e) NULL)
  if (!is.null(ok)) return(vec_mark(ok, all$done))
  # one of the rewrites does not check: add them one at a time
  allow <- integer(0); best <- ir
  for (k in seq_len(all$found)) {
    r <- rewrite(c(allow, k))
    chk <- tryCatch(stan_check(r$ast, external = external), bayescomp_error = function(e) NULL)
    if (!is.null(chk)) { allow <- c(allow, k); best <- chk }
  }
  vec_mark(best, length(allow))
}

vec_mark <- function(ir, n) { attr(ir, "vectorized_loops") <- n; ir }

vec_as_list <- function(st) if (identical(st$kind, "block_stmt")) st$body else list(st)
vec_block <- function(stmts, loc) list(kind = "block_stmt", body = stmts, loc = loc)

# statements with qualifying loops rewritten (a loop may become several
# statements); descends into blocks, conditionals and other loops. Loops are
# numbered in order; st$allow (NULL for all) selects which are rewritten.
vec_stmts <- function(as, ts, st) {
  if (length(as) != length(ts)) return(as)
  out <- list()
  for (i in seq_along(as)) {
    a <- as[[i]]; t <- ts[[i]]
    if (identical(a$kind, "for") && identical(t$kind, "for")) {
      new <- tryCatch(vec_for(a, t), error = function(e) NULL)
      if (!is.null(new)) {
        st$k <- st$k + 1L
        if (is.null(st$allow) || st$k %in% st$allow) {
          st$done <- st$done + 1L
          out <- c(out, new)
          next
        }
      } else {
        a$body <- vec_block(vec_stmts(vec_as_list(a$body), vec_as_list(t$body), st), a$body$loc %||% a$loc)
      }
    } else if (identical(a$kind, "block_stmt") && identical(t$kind, "block_stmt")) {
      # also the body of a profile block, which stays a profile block
      a$body <- vec_stmts(a$body, t$body, st)
    } else if (identical(a$kind, "if") && identical(t$kind, "if")) {
      for (k in c("then", "els")) if (!is.null(a[[k]]) && !is.null(t[[k]]))
        a[[k]] <- vec_block(vec_stmts(vec_as_list(a[[k]]), vec_as_list(t[[k]]), st), a$loc)
    }
    out[[length(out) + 1L]] <- a
  }
  out
}

# does an AST mention variable `v` anywhere?
vec_mentions <- function(e, v) {
  if (!is.list(e)) return(FALSE)
  if (identical(e$kind, "var") && identical(e$name, v)) return(TRUE)
  for (k in names(e) %||% seq_along(e)) {
    if (k %in% c("loc", "type", "sig")) next
    if (is.list(e[[k]]) && vec_mentions(e[[k]], v)) return(TRUE)
  }
  FALSE
}

# does an AST call target(), a random number generator, an _lp function
# (which can change the target) or a _jacobian function?
vec_effectful <- function(e) {
  if (!is.list(e)) return(FALSE)
  if (identical(e$kind, "call") && is.character(e$name) &&
      (identical(e$name, "target") || grepl("_rng$|_lp$|_jacobian$", e$name))) return(TRUE)
  for (k in names(e) %||% seq_along(e)) {
    if (k %in% c("loc", "type", "sig")) next
    if (is.list(e[[k]]) && vec_effectful(e[[k]])) return(TRUE)
  }
  FALSE
}

vec_is_n <- function(e, n) identical(e$kind, "var") && identical(e$name, n)

# names of the variables an AST reads
ir_var_names <- function(e) {
  if (!is.list(e)) return(character())
  out <- if (identical(e$kind, "var") && is.character(e$name)) e$name else character()
  for (k in names(e) %||% seq_along(e)) if (!k %in% c("loc", "type", "sig") && is.list(e[[k]])) out <- c(out, ir_var_names(e[[k]]))
  unique(out)
}

# does a statement assign to any of the variables `vs`?
ir_assigns_any <- function(st, vs) {
  if (!is.list(st)) return(FALSE)
  if (identical(st$kind, "assign") && any(vapply(vs, function(v) vec_mentions(st$lhs, v), TRUE))) return(TRUE)
  for (k in names(st) %||% seq_along(st)) if (!k %in% c("loc", "type", "sig") && is.list(st[[k]]) && ir_assigns_any(st[[k]], vs)) return(TRUE)
  FALSE
}

# Stan tests the upper bound of a for loop before every iteration; it can
# change when it reads target() or a variable the body assigns
loop_bound_varies <- function(hi, body) vec_effectful(hi) || ir_assigns_any(body, ir_var_names(hi))

vec_scalar <- function(t) is.list(t) && identical(t$nd %||% 0L, 0L) && isTRUE(t$base %in% c("int", "real"))

vec_for <- function(a, t) {
  n <- a$var
  if (!identical(t$iter, NULL) || is.null(a$lo) || is.null(a$hi)) return(NULL)
  if (vec_mentions(a$lo, n) || vec_mentions(a$hi, n)) return(NULL)
  # the bounds are evaluated once (and again by the guard below), which is
  # only the same as the loop when they cannot change between evaluations
  if (vec_effectful(a$lo) || vec_effectful(a$hi) || loop_bound_varies(a$hi, a$body) ||
      ir_assigns_any(a$body, ir_var_names(a$lo))) return(NULL)
  abody <- vec_as_list(a$body); tbody <- vec_as_list(t$body)
  if (!length(abody) || length(abody) != length(tbody)) return(NULL)
  # scalar locals of the body are replaced by their values
  inl <- vec_inline_locals(abody, tbody, n)
  if (is.null(inl)) return(NULL)
  abody <- inl$a; tbody <- inl$t
  if (!length(abody)) return(NULL)
  # variables written in the loop
  written <- character()
  for (k in seq_along(abody)) {
    st <- abody[[k]]
    if (identical(st$kind, "assign")) {
      if (!identical(st$lhs$kind, "index") || !identical(st$lhs$x$kind, "var")) return(NULL)
      written <- c(written, st$lhs$x$name)
    } else if (!st$kind %in% c("tilde", "target_add")) return(NULL)
  }
  if (n %in% written) return(NULL)
  # a written variable may only appear as v[n]
  ok_use <- function(e) {
    if (!is.list(e)) return(TRUE)
    if (identical(e$kind, "index") && identical(e$x$kind, "var") && e$x$name %in% written)
      return(length(e$idx) == 1L && identical(e$idx[[1L]]$kind, "single") && vec_is_n(e$idx[[1L]]$e, n))
    if (identical(e$kind, "var") && e$name %in% written) return(FALSE)
    for (k in names(e) %||% seq_along(e)) if (!k %in% c("loc") && is.list(e[[k]]) && !ok_use(e[[k]])) return(FALSE)
    TRUE
  }
  for (st in abody) if (!ok_use(st)) return(NULL)
  ctx <- list(n = n, lo = a$lo, hi = a$hi, loc = a$loc)
  out <- list()
  for (k in seq_along(abody)) {
    s <- vec_stmt(abody[[k]], tbody[[k]], ctx)
    if (is.null(s)) return(NULL)
    out[[k]] <- s
  }
  # a loop with no iterations evaluates nothing, not even expressions that do
  # not use the loop variable
  cond <- list(kind = "binop", op = "<=", a = a$lo, b = a$hi, loc = a$loc)
  list(list(kind = "if", cond = cond, then = vec_block(out, a$loc), els = NULL, loc = a$loc))
}

# `real mu = e;` in a loop body: later statements read e instead, so the
# body is a sequence of statements the rewrite handles. Only for a scalar
# with a value without effects that the body does not assign again.
vec_inline_locals <- function(as, ts, n) {
  i <- 1L
  while (i <= length(as)) {
    a <- as[[i]]; t <- ts[[i]]
    if (!identical(a$kind, "decl")) { i <- i + 1L; next }
    if (is.null(a$init) || is.null(t$init) || !vec_scalar(t$vtype) || vec_effectful(a$init)) return(NULL)
    rest <- seq_along(as)[-seq_len(i)]
    if (any(vapply(as[rest], ir_assigns_any, TRUE, vs = a$name))) return(NULL)
    ra <- list(kind = "paren", a = a$init, loc = a$init$loc)
    for (k in rest) { as[[k]] <- vec_subst(as[[k]], a$name, ra); ts[[k]] <- vec_subst(ts[[k]], a$name, t$init) }
    as <- as[-i]; ts <- ts[-i]
  }
  list(a = as, t = ts)
}
vec_subst <- function(e, v, r) {
  if (!is.list(e)) return(e)
  if (identical(e$kind, "var") && identical(e$name, v)) return(r)
  for (k in names(e) %||% seq_along(e)) if (!k %in% c("loc", "type", "sig") && is.list(e[[k]])) e[[k]] <- vec_subst(e[[k]], v, r)
  e
}

vec_range <- function(ctx) list(kind = "range", lo = ctx$lo, hi = ctx$hi, loc = ctx$loc)
vec_call <- function(name, args, loc) list(kind = "call", name = name, args = args, cond_bar = FALSE, loc = loc)
vec_size <- function(ctx) {
  one <- list(kind = "int", value = 1L, text = "1", loc = ctx$loc)
  zero <- list(kind = "int", value = 0L, text = "0", loc = ctx$loc)
  d <- list(kind = "binop", op = "+", a = list(kind = "binop", op = "-", a = ctx$hi, b = ctx$lo, loc = ctx$loc), b = one, loc = ctx$loc)
  vec_call("max", list(d, zero), ctx$loc)
}

vec_stmt <- function(a, t, ctx) {
  switch(a$kind,
    tilde = {
      # a user defined density may have unrelated scalar and vector overloads
      if (isTRUE(t$sig$user) || !vec_scalar(t$lhs$type)) return(NULL)
      # truncation bounds are scalars evaluated once
      if (!is.null(a$trunc) && any(vapply(list(a$trunc$lo, a$trunc$hi), function(b) vec_mentions(b, ctx$n) || vec_effectful(b), TRUE)))
        return(NULL)
      if (!all(vec_args_ok(a$args, t$args, t$sig$args[-1L], ctx$n))) return(NULL)
      lhs <- vec_expr(a$lhs, t$lhs, ctx)
      if (is.null(lhs)) return(NULL)
      args <- Map(function(x, y) vec_expr(x, y, ctx), a$args, t$args)
      if (any(vapply(args, is.null, TRUE))) return(NULL)
      if (!lhs$dep && !any(vapply(args, `[[`, TRUE, "dep"))) return(NULL)
      a$lhs <- vec_out(lhs, ctx, arith = FALSE)
      a$args <- lapply(args, vec_out, ctx = ctx, arith = FALSE)
      a
    },
    target_add = {
      v <- a$value; tv <- t$value
      if (!(identical(v$kind, "call") && grepl("_(lpdf|lpmf|lupdf|lupmf|lcdf|lccdf)$", v$name))) {
        # any other scalar: target += of the vector adds its sum
        if (!vec_scalar(tv$type)) return(NULL)
        r <- vec_expr(v, tv, ctx)
        if (is.null(r) || !r$dep || identical(r$kind, "mat")) return(NULL)
        a$value <- vec_out(r, ctx, arith = TRUE)
        return(a)
      }
      if (!identical(tv$kind, "call") || isTRUE(tv$sig$user) || length(v$args) != length(tv$args)) return(NULL)
      if (!all(vec_args_ok(v$args, tv$args, tv$sig$args, ctx$n))) return(NULL)
      args <- Map(function(x, y) vec_expr(x, y, ctx), v$args, tv$args)
      if (any(vapply(args, is.null, TRUE)) || !any(vapply(args, `[[`, TRUE, "dep"))) return(NULL)
      v$args <- lapply(args, vec_out, ctx = ctx, arith = FALSE)
      a$value <- v
      a
    },
    assign = {
      if (!a$op %in% c("=", "+=", "-=")) return(NULL)
      lt <- t$lhs$x$type
      if (!identical(lt$base, "vector") || !identical(lt$nd %||% 0L, 0L)) return(NULL)
      if (length(a$lhs$idx) != 1L || !identical(a$lhs$idx[[1L]]$kind, "single") || !vec_is_n(a$lhs$idx[[1L]]$e, ctx$n)) return(NULL)
      if (!vec_scalar(t$rhs$type)) return(NULL)
      rhs <- vec_expr(a$rhs, t$rhs, ctx)
      if (is.null(rhs)) return(NULL)
      a$lhs$idx <- list(vec_range(ctx))
      a$rhs <- if (rhs$dep) vec_out(rhs, ctx, arith = TRUE) else vec_call("rep_vector", list(rhs$ast, vec_size(ctx)), ctx$loc)
      a
    },
    NULL)
}

# arguments a vectorized density accepts: scalars, and in a position the
# signature shares between elements (cutpoints, a simplex) a vector that
# does not depend on the loop. In an elementwise position ("reals") a vector
# would broadcast the outcome instead.
vec_args_ok <- function(as, ts, pos, n) {
  Map(function(a, t, p) vec_scalar(t$type) ||
        (!is.null(p) && !p %in% c("reals", "ints", "real", "int") && identical(t$type$base, "vector") &&
         (t$type$nd %||% 0L) == 0L && !vec_mentions(a, n) && !vec_effectful(a)),
      as, ts, as.list(pos)[seq_along(ts)]) |> unlist()
}

# the final form of a vectorized expression; in arithmetic and as an
# assigned value, arrays become vectors
vec_out <- function(r, ctx, arith) {
  if (!r$dep) return(r$ast)
  if (arith && r$kind %in% c("arr", "row")) return(vec_call("to_vector", list(r$ast), ctx$loc))
  r$ast
}

# vectorize one expression: list(ast, dep, kind) where kind is "vec"
# (vector), "row" (row vector), "arr" (array of reals or ints) or "mat"
# (rows of a matrix, only valid on the left of a matrix product)
vec_expr <- function(a, t, ctx) {
  n <- ctx$n
  if (!vec_mentions(a, n)) {
    # an expression that does not use the loop variable is evaluated once
    # instead of once per iteration, which is only the same when it has no
    # effects and does not read state the loop changes
    if (vec_effectful(a)) return(NULL)
    return(list(ast = a, dep = FALSE, kind = "scalar"))
  }
  # the typed program has no parentheses
  if (identical(a$kind, "paren") && !identical(t$kind, "paren")) {
    r <- vec_expr(a$a, t, ctx)
    if (is.null(r) || identical(r$kind, "mat")) return(NULL)
    a$a <- r$ast
    return(list(ast = a, dep = r$dep, kind = r$kind))
  }
  if (!identical(a$kind, t$kind)) return(NULL)
  arith <- function(r) if (r$dep) vec_out(r, ctx, arith = TRUE) else r$ast
  switch(a$kind,
    paren = {
      r <- vec_expr(a$x %||% a$a, t$x %||% t$a, ctx)
      if (is.null(r) || identical(r$kind, "mat")) return(NULL)
      if (!is.null(a$x)) a$x <- r$ast else a$a <- r$ast
      list(ast = a, dep = r$dep, kind = r$kind)
    },
    index = vec_index(a, t, ctx),
    unop = {
      if (!identical(a$op, "-") || !vec_scalar(t$type)) return(NULL)
      r <- vec_expr(a$a, t$a, ctx)
      if (is.null(r) || identical(r$kind, "mat")) return(NULL)
      a$a <- arith(r)
      list(ast = a, dep = TRUE, kind = "vec")
    },
    binop = {
      if (!a$op %in% c("+", "-", "*", "/", "^")) return(NULL)
      ra <- vec_expr(a$a, t$a, ctx); rb <- vec_expr(a$b, t$b, ctx)
      if (is.null(ra) || is.null(rb)) return(NULL)
      # X[n] * b: a row of a matrix times a vector
      if (identical(ra$kind, "mat")) {
        if (a$op != "*" || rb$dep || !identical(t$b$type$base, "vector") || (t$b$type$nd %||% 0L) != 0L) return(NULL)
        a$a <- ra$ast
        return(list(ast = a, dep = TRUE, kind = "vec"))
      }
      if (identical(rb$kind, "mat") || !vec_scalar(t$type)) return(NULL)
      if (!vec_scalar(t$a$type) || !vec_scalar(t$b$type)) return(NULL)
      # integer division rounds; on vectors it would not
      if (a$op == "/" && identical(t$a$type$base, "int") && identical(t$b$type$base, "int")) return(NULL)
      both <- ra$dep && rb$dep
      a$op <- switch(a$op, "*" = if (both) ".*" else "*", "/" = if (rb$dep) "./" else "/", "^" = ".^", a$op)
      a$a <- arith(ra); a$b <- arith(rb)
      list(ast = a, dep = TRUE, kind = "vec")
    },
    call = {
      # dot_product(X[n], b): a row of a matrix times a vector
      if (identical(a$name, "dot_product") && length(a$args) == 2L && !isTRUE(t$sig$user)) {
        ra <- vec_expr(a$args[[1L]], t$args[[1L]], ctx)
        tb <- t$args[[2L]]$type
        if (is.null(ra) || !identical(ra$kind, "mat") || vec_mentions(a$args[[2L]], n) || vec_effectful(a$args[[2L]]) ||
            !identical(tb$base, "vector") || (tb$nd %||% 0L) != 0L) return(NULL)
        return(list(ast = list(kind = "binop", op = "*", a = ra$ast, b = a$args[[2L]], loc = a$loc), dep = TRUE, kind = "vec"))
      }
      if (!a$name %in% vec_unary || length(a$args) != 1L || isTRUE(t$sig$user) || !vec_scalar(t$type)) return(NULL)
      r <- vec_expr(a$args[[1L]], t$args[[1L]], ctx)
      if (is.null(r) || identical(r$kind, "mat")) return(NULL)
      a$args <- list(arith(r))
      list(ast = a, dep = TRUE, kind = "vec")
    },
    NULL)
}

vec_index <- function(a, t, ctx) {
  n <- ctx$n
  if (vec_mentions(a$x, n)) return(NULL)
  xt <- t$x$type
  nd <- xt$nd %||% 0L
  # the vectorized position: the loop variable itself, or idx[n] for an
  # integer array idx that does not depend on the loop
  pos <- function(ia, it) {
    if (!identical(ia$kind, "single")) return(NULL)
    e <- ia$e
    if (vec_is_n(e, n)) return(vec_range(ctx))
    # x[n + c] or x[n - c] with an integer c that does not depend on the loop
    if (identical(e$kind, "binop") && e$op %in% c("+", "-") && identical(it$e$kind, "binop") && vec_is_n(e$a, n) &&
        !vec_mentions(e$b, n) && identical(it$e$b$type$base, "int") && (it$e$b$type$nd %||% 0L) == 0L) {
      sh <- function(x) list(kind = "binop", op = e$op, a = list(kind = "paren", a = x, loc = e$loc),
                             b = list(kind = "paren", a = e$b, loc = e$loc), loc = e$loc)
      return(list(kind = "range", lo = sh(ctx$lo), hi = sh(ctx$hi), loc = ctx$loc))
    }
    if (identical(e$kind, "index") && identical(it$e$kind, "index") && !vec_mentions(e$x, n) &&
        identical(it$e$x$type$base, "int") && (it$e$x$type$nd %||% 0L) == 1L && length(e$idx) == 1L &&
        identical(e$idx[[1L]]$kind, "single") && vec_is_n(e$idx[[1L]]$e, n))
      return(list(kind = "single", e = list(kind = "index", x = e$x, idx = list(vec_range(ctx)), loc = e$loc), loc = ia$loc))
    NULL
  }
  if (length(a$idx) == 1L) {
    p <- pos(a$idx[[1L]], t$idx[[1L]])
    if (is.null(p)) return(NULL)
    a$idx <- list(p)
    if (nd == 1L && xt$base %in% c("real", "int")) return(list(ast = a, dep = TRUE, kind = "arr"))
    if (nd == 0L && identical(xt$base, "vector")) return(list(ast = a, dep = TRUE, kind = "vec"))
    if (nd == 0L && identical(xt$base, "row_vector")) return(list(ast = a, dep = TRUE, kind = "row"))
    if (nd == 0L && identical(xt$base, "matrix")) return(list(ast = a, dep = TRUE, kind = "mat"))
    return(NULL)
  }
  # X[n, k] with k fixed: a column slice
  if (length(a$idx) == 2L && nd == 0L && identical(xt$base, "matrix")) {
    p <- pos(a$idx[[1L]], t$idx[[1L]])
    k <- a$idx[[2L]]
    if (is.null(p) || !identical(k$kind, "single") || vec_mentions(k, n) || !identical(t$idx[[2L]]$e$type$base, "int") ||
        (t$idx[[2L]]$e$type$nd %||% 0L) != 0L) return(NULL)
    a$idx <- list(p, k)
    return(list(ast = a, dep = TRUE, kind = "vec"))
  }
  NULL
}
