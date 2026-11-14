# R backend code generation: typed IR to R closures.
#
# Code is built as R language objects with as.call(), never by parsing text,
# so identifiers and strings from the Stan program cannot inject R code.
# Stan variables become symbols "V.<name>" and user functions "F.<name>.<k>";
# neither form can be produced by a Stan identifier, which contains no dots,
# so generated names never collide with runtime objects. Library functions
# are embedded as function objects in the call position.

vsym <- function(name) as.name(paste0("V.", name))
cl <- function(f, ...) as.call(c(list(f), list(...)))
cl_list <- function(f, args) as.call(c(list(f), args))
blk <- function(stmts) as.call(c(list(as.name("{")), stmts))
CTX <- quote(.ctx)

promote_to <- function(code, from, to) {
  if (!is.null(to) && !is.null(from) && to$base == "real" && from$base == "int") return(cl(as_dbl, code))
  code
}

lib_impl <- function(name, types) {
  for (sg in .lib[[name]]) {
    if (length(sg$args) == length(types) && all(mapply(pat_ok, sg$args, types))) return(sg$impl)
  }
  stop("internal error: no library implementation for ", name)
}

# ---- expressions -------------------------------------------------------------------

cg_expr <- function(e) {
  switch(e$kind,
    int = e$value,
    real = e$value,
    var = vsym(e$name),
    paren = cg_expr(e$a),
    unop = {
      a <- cg_expr(e$a)
      switch(e$op,
        "-" = if (e$type$base == "int") cl(op_int_neg, a) else cl(op_neg, a),
        "!" = cl(op_not, a),
        "+" = a)
    },
    binop = cg_binop(e),
    cond = {
      a <- promote_to(cg_expr(e$a), e$a$type, e$type)
      b <- promote_to(cg_expr(e$b), e$b$type, e$type)
      cl(as.name("if"), cl(truthy, cg_expr(e$cond)), a, b)
    },
    transpose = {
      a <- cg_expr(e$a)
      if (e$a$type$base == "matrix") cl(op_transpose_mat, a) else cl(as_dbl, a)
    },
    index = cg_index(e),
    array_expr = {
      items <- lapply(e$items, function(it) promote_to(cg_expr(it), it$type, e$elem_type))
      cl(rt_array, cl_list(as.name("list"), items), e$elem_type$base == "int",
         e$elem_type$base %in% c("int", "real") && (e$elem_type$nd %||% 0L) == 0L)
    },
    rowvec_expr = {
      items <- lapply(e$items, cg_expr)
      cl(rt_rowvec, cl_list(as.name("list"), items), e$type$base == "matrix")
    },
    call = cg_call(e),
    stop("internal error: cannot generate code for expression ", e$kind))
}

cg_index_parts <- function(e) {
  kinds <- vapply(e$idx, `[[`, "", "kind")
  vals <- lapply(e$idx, function(ix) switch(ix$kind,
    single = cg_expr(ix$e),
    intarray = cg_expr(ix$e),
    range = cl(as.name("c"), cg_expr(ix$lo), cg_expr(ix$hi)),
    lower = cg_expr(ix$lo),
    upper = cg_expr(ix$hi),
    omni = 0L))
  list(kinds = kinds, vals = cl_list(as.name("list"), vals))
}

cg_index <- function(e) {
  p <- cg_index_parts(e)
  cl(rt_index, cg_expr(e$x), p$kinds, p$vals)
}

cg_binop <- function(e) {
  a <- cg_expr(e$a)
  b <- cg_expr(e$b)
  ta <- e$a$type
  tb <- e$b$type
  both_int <- is_int_ty(ta) && is_int_ty(tb)
  sca <- is_scalar_ty(ta)
  scb <- is_scalar_ty(tb)
  # two containers: sizes must agree even when one has a single element
  ew <- function(f, what) if (sca || scb) cl(f, a, b) else cl(rt_elementwise, f, a, b, what)
  switch(e$op,
    "+" = if (both_int) cl(op_int_add, a, b) else ew(op_add, "addition"),
    "-" = if (both_int) cl(op_int_sub, a, b) else ew(op_sub, "subtraction"),
    "*" = {
      if (both_int) return(cl(op_int_mul, a, b))
      if (sca || scb) return(cl(op_scale, a, b))
      switch(paste(ta$base, tb$base),
        "matrix matrix" = cl(op_mat_mat, a, b),
        "matrix vector" = cl(op_mat_vec, a, b),
        "row_vector matrix" = cl(op_rvec_mat, a, b),
        "row_vector vector" = cl(op_dot, a, b),
        "vector row_vector" = cl(op_outer, a, b))
    },
    "/" = {
      if (both_int) return(cl(op_int_div, a, b))
      if (sca || scb) return(cl(op_div_scalar, a, b))
      cl(lib_impl("mdivide_right", list(ta, tb)), a, b)
    },
    "%" = cl(op_int_mod, a, b),
    "%/%" = cl(op_int_div, a, b),
    "\\" = cl(lib_impl("mdivide_left", list(ta, tb)), a, b),
    ".*" = if (sca || scb) cl(op_scale, a, b) else ew(op_elt_mul, "elementwise multiplication"),
    "./" = if (sca || scb) cl(op_div_scalar, a, b) else ew(op_elt_div, "elementwise division"),
    "^" = , ".^" = cl(op_pow, a, b),
    "==" = , "!=" = , "<" = , "<=" = , ">" = , ">=" = cl(op_cmp, e$op, a, b),
    "&&" = cl(as.name("as.integer"), cl(as.name("&&"), cl(truthy, a), cl(truthy, b))),
    "||" = cl(as.name("as.integer"), cl(as.name("||"), cl(truthy, a), cl(truthy, b))),
    stop("internal error: operator ", e$op))
}

# containers passed to the vectorized ("reals"/"ints") arguments of a library
# distribution must have equal sizes: a vector[1] is a container, not a
# scalar, so it may not broadcast against a vector[3]
is_container_ty <- function(t) t$nd > 0L || t$base %in% c("vector", "row_vector", "matrix")
cg_size_guard <- function(fname, sg, types, call) {
  if (is.null(sg$dist) || length(sg$args) != length(types)) return(call)
  mask <- vapply(seq_along(types), function(i) sg$args[[i]] %in% c("reals", "ints") && is_container_ty(types[[i]]), TRUE)
  if (sum(mask) < 2L) return(call)
  args <- as.list(call)[-1L]
  cl_list(rt_dist_call, c(list(fname, call[[1L]], mask), unname(args)))
}

cg_call <- function(e) {
  if (identical(e$builtin, "target")) return(cl(rt_target_get, CTX))
  sg <- e$sig
  args <- lapply(e$args, cg_expr)
  if (isTRUE(sg$user)) {
    args <- Map(function(code, a, t) promote_to(code, a$type, t), args, e$args, sg$args)
    call <- cl_list(as.name(sg$mangled), c(unname(args), list(CTX)))
    # a user density called by its normalized name runs with propto off
    if (grepl("_(lpdf|lpmf)$", e$name)) call <- cl(rt_normalized, call, CTX)
    return(call)
  }
  if (e$name == "dims") {
    t <- e$args[[1L]]$type
    return(cl(dims_typed, args[[1L]], t$base, t$nd))
  }
  types <- lapply(e$args, `[[`, "type")
  if (grepl("_(lupdf|lupmf)$", e$name)) {
    full <- lib_impl(sub("_lup", "_lp", e$name), types)
    return(cl(as.name("if"), cl(as.name("isTRUE"), cl(as.name("$"), CTX, as.name("propto"))),
              cg_size_guard(e$name, sg, types, cl_list(sg$impl, unname(args))),
              cg_size_guard(e$name, sg, types, cl_list(full, unname(args)))))
  }
  cg_size_guard(e$name, sg, types, cl_list(sg$impl, unname(args)))
}

# ---- statements ----------------------------------------------------------------------

cg_stmts <- function(stmts, ctx) {
  out <- list()
  for (st in stmts) out <- c(out, cg_stmt(st, ctx))
  out
}

decl_kind <- function(t) {
  if (t$base %in% constrained_bases) return(t$base)
  if (length(t$cons)) return("bounds")
  "none"
}

# code computing c(adims) and eshape for a declaration
cg_decl_shape <- function(t) {
  ad <- cl_list(as.name("c"), lapply(t$adims, cg_expr))
  sz <- lapply(t$sizes, cg_expr)
  es <- switch(t$base,
    int = , real = quote(integer(0)),
    vector = , row_vector = , simplex = , ordered = , positive_ordered = , unit_vector = ,
    sum_to_zero_vector = cl(as.name("c"), sz[[1L]]),
    cholesky_factor_corr = , corr_matrix = , cov_matrix = cl(as.name("c"), sz[[1L]], sz[[1L]]),
    cholesky_factor_cov = if (length(sz) == 1L) cl(as.name("c"), sz[[1L]], sz[[1L]]) else cl(as.name("c"), sz[[1L]], sz[[2L]]),
    cl(as.name("c"), sz[[1L]], sz[[2L]]))
  list(adims = ad, eshape = es)
}

cg_cons_args <- function(t) {
  cl_list(as.name("list"), lapply(t$cons, cg_expr))
}

cg_stmt <- function(st, ctx) {
  switch(st$kind,
    empty = list(),
    decl = {
      t <- st$type
      sh <- cg_decl_shape(t)
      base <- decl_base_type(t$base)
      if (st$top && ctx$block %in% c("transformed parameters", "generated quantities")) {
        # sizes were evaluated when the model was instantiated
        mk <- cl(rt_decl, if (base == "int") "int" else "real",
                 cl(as.name("$"), cl(as.name("[["), quote(.decls), st$name), as.name("adims")),
                 cl(as.name("$"), cl(as.name("[["), quote(.decls), st$name), as.name("eshape")))
      } else {
        mk <- cl(rt_decl, if (base == "int") "int" else "real", sh$adims, sh$eshape)
      }
      out <- list(cl(as.name("<-"), vsym(st$name), mk))
      if (!is.null(st$init)) {
        rhs <- promote_to(cg_expr(st$init), st$init$type, st$vtype)
        out <- c(out, list(cl(as.name("<-"), vsym(st$name), cl(rt_assign, vsym(st$name), rhs, st$name))))
      }
      out
    },
    block_stmt = list(blk(cg_stmts(st$body, ctx))),
    assign = {
      lhs <- st$lhs
      rhs <- st$rhs
      rcode <- cg_expr(rhs)
      if (st$op != "=") {
        bop <- substr(st$op, 1L, nchar(st$op) - 1L)
        fake <- list(kind = "binop", op = bop, a = lhs, b = rhs, type = binop_type(bop, lhs$type, rhs$type))
        rcode <- cg_binop(fake)
        rtype <- fake$type
      } else {
        rtype <- rhs$type
      }
      rcode <- promote_to(rcode, rtype, lhs$type)
      base <- lhs
      while (base$kind == "index") base <- base$x
      name <- base$name
      if (lhs$kind == "var") {
        return(list(cl(as.name("<-"), vsym(name), cl(rt_assign, vsym(name), rcode, name))))
      }
      # nested indexing a[i][j] = v is rewritten as a[i, j] when possible
      chain <- list()
      cur <- lhs
      while (cur$kind == "index") {
        chain <- c(list(cur), chain)
        cur <- cur$x
      }
      if (length(chain) == 1L) {
        p <- cg_index_parts(lhs)
        return(list(cl(as.name("<-"), vsym(name), cl(rt_assign_index, vsym(name), p$kinds, p$vals, rcode, name))))
      }
      # general nested case: read the inner container, assign into it, write it back
      tmp <- rcode
      for (lvl in rev(seq_along(chain))) {
        node <- chain[[lvl]]
        container <- if (lvl == 1L) vsym(name) else cg_expr(chain[[lvl - 1L]])
        p <- cg_index_parts(node)
        tmp <- cl(rt_assign_index, container, p$kinds, p$vals, tmp, name)
      }
      list(cl(as.name("<-"), vsym(name), tmp))
    },
    tilde = {
      lhs <- cg_expr(st$lhs)
      args <- lapply(st$args, cg_expr)
      sg <- st$sig
      call <- if (isTRUE(sg$user)) {
        aa <- Map(function(code, a, t) promote_to(code, a$type, t), c(list(lhs), args), c(list(st$lhs), st$args), sg$args)
        cl_list(as.name(sg$mangled), c(unname(aa), list(CTX)))
      } else {
        types <- c(list(st$lhs$type), lapply(st$args, `[[`, "type"))
        full <- lib_impl(sub("_lup", "_lp", st$fname), types)
        cl(as.name("if"), cl(as.name("isTRUE"), cl(as.name("$"), CTX, as.name("propto"))),
           cg_size_guard(st$fname, sg, types, cl_list(sg$impl, c(list(lhs), unname(args)))),
           cg_size_guard(st$fname, sg, types, cl_list(full, c(list(lhs), unname(args)))))
      }
      out <- list(cl(rt_target_add, CTX, call))
      if (!is.null(st$trunc)) {
        tr <- st$trunc
        mkfun <- function(sig) {
          if (is.null(sig)) return(NULL)
          if (isTRUE(sig$user)) return(as.name(sig$mangled))
          # the sampling statement has already checked the parameters
          attr(sig$impl, "unchecked") %||% sig$impl
        }
        lcdf_f <- mkfun(tr$lcdf_sig)
        lccdf_f <- mkfun(tr$lccdf_sig)
        if (isTRUE(tr$lcdf_sig$user) || isTRUE(tr$lccdf_sig$user))
          stop_unsupported("truncation with user defined cumulative distribution functions is not supported", st$loc)
        out <- c(out, list(cl(rt_truncate, CTX, lhs,
                              if (is.null(tr$lo)) NULL else cg_expr(tr$lo),
                              if (is.null(tr$hi)) NULL else cg_expr(tr$hi),
                              lcdf_f, lccdf_f, cl_list(as.name("list"), unname(args)),
                              isTRUE(tr$discrete), tr$fam)))
      }
      out
    },
    target_add = list(cl(rt_target_add, CTX, cg_expr(st$value))),
    jacobian_add = list(cl(as.name("if"), cl(as.name("isTRUE"), cl(as.name("$"), CTX, as.name("jacobian"))),
                           cl(rt_target_add, CTX, cg_expr(st$value)))),
    "if" = {
      th <- blk(cg_stmt(st$then, ctx))
      if (is.null(st$els)) list(cl(as.name("if"), cl(truthy, cg_expr(st$cond)), th))
      else list(cl(as.name("if"), cl(truthy, cg_expr(st$cond)), th, blk(cg_stmt(st$els, ctx))))
    },
    "while" = list(cl(as.name("while"), cl(truthy, cg_expr(st$cond)), blk(cg_stmt(st$body, ctx)))),
    "for" = list(cl(as.name("for"), vsym(st$var), cl(rt_range, cg_expr(st$lo), cg_expr(st$hi)),
                    blk(cg_stmt(st$body, ctx)))),
    foreach = list(cl(as.name("for"), vsym(st$var), cl(rt_foreach, cg_expr(st$iter), st$iter$type$nd, st$iter$type$base),
                      blk(cg_stmt(st$body, ctx)))),
    "break" = list(quote(break)),
    "continue" = list(quote(next)),
    "return" = {
      if (is.null(st$value)) list(quote(return(invisible(NULL))))
      else list(cl(as.name("return"), promote_to(cg_expr(st$value), st$value$type, st$rtype)))
    },
    print = list(cl(rt_print, CTX, cl_list(as.name("list"), lapply(st$items, function(it) if (it$kind == "string") it$value else cg_expr(it))))),
    reject = list(cl(rt_reject, cl_list(as.name("list"), lapply(st$items, function(it) if (it$kind == "string") it$value else cg_expr(it))), FALSE)),
    fatal_error = list(cl(rt_reject, cl_list(as.name("list"), lapply(st$items, function(it) if (it$kind == "string") it$value else cg_expr(it))), TRUE)),
    call_stmt = list(cg_expr(st$call)),
    stop("internal error: cannot generate code for statement ", st$kind))
}

# ---- user functions -------------------------------------------------------------------

cg_fundef <- function(fd) {
  argnames <- c(paste0("V.", vapply(fd$args, `[[`, "", "name")), ".ctx")
  formals <- rep(list(quote(expr = )), length(argnames))
  names(formals) <- argnames
  body <- blk(c(cg_stmt(fd$body, list(block = "functions")), list(quote(invisible(NULL)))))
  list(name = fd$sig$mangled, formals = formals, body = body)
}

make_function <- function(formals, body, env) {
  f <- as.function(c(formals, list(body)), envir = env)
  f
}
