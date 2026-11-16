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
prof_state <- new.env(parent = emptyenv())
prof_counter <- function() { prof_state$n <- (prof_state$n %||% 0L) + 1L; prof_state$n }

blk <- function(stmts) as.call(c(list(as.name("{")), stmts))
CTX <- quote(.ctx)

promote_to <- function(code, from, to) {
  if (!is.null(to) && !is.null(from) && to$base == "real" && from$base == "int") return(cl(as_dbl, code))
  if (!is.null(to) && !is.null(from) && (is_tuple_ty(to) || is_tuple_ty(from))) return(if (ty_eq(from, to)) code else cl(rt_promote, code, from, to))
  if (!is.null(to) && !is.null(from) && is_cx_ty(to) && !is_cx_ty(from)) return(cl(cx_promote, code, is_scalar_ty(from)))
  code
}

lib_impl <- function(name, types) lib_sig(name, types)$impl

# ---- expressions -------------------------------------------------------------------

cg_expr <- function(e) {
  switch(e$kind,
    int = e$value,
    real = e$value,
    imag = c(0, e$value),
    code = e$code,
    tuple_expr = cl(rt_tuple, cl_list(as.name("list"), lapply(e$items, cg_expr))),
    proj = cl(rt_proj, cg_expr(e$x), e$k),
    var = vsym(e$name),
    paren = cg_expr(e$a),
    unop = {
      a <- cg_expr(e$a)
      switch(e$op,
        "-" = if (e$type$base == "int") cl(op_int_neg, a) else cl(op_neg, a),
        "!" = cl(op_not, a),
        "+" = a)
    },
    binop = if (is_cx_ty(e$type) || is_cx_ty(e$a$type) || is_cx_ty(e$b$type)) {
      cl(cx_binop, e$op, cg_expr(e$a), cg_expr(e$b), e$a$type$base, e$b$type$base)
    } else cg_binop(e),
    cond = {
      a <- promote_to(cg_expr(e$a), e$a$type, e$type)
      b <- promote_to(cg_expr(e$b), e$b$type, e$type)
      cl(as.name("if"), cl(truthy, cg_expr(e$cond)), a, b)
    },
    transpose = {
      a <- cg_expr(e$a)
      if (e$a$type$base == "matrix") cl(op_transpose_mat, a)
      else if (e$a$type$base == "complex_matrix") cl(cx_transpose, a)
      else if (is_cx_ty(e$a$type)) a
      else cl(as_dbl, a)
    },
    index = cg_index(e),
    array_expr = {
      items <- lapply(e$items, function(it) promote_to(cg_expr(it), it$type, e$elem_type))
      if (is_tuple_ty(e$elem_type)) return(cl(rt_tuple_array, cl_list(as.name("list"), items)))
      cl(rt_array, cl_list(as.name("list"), items), e$elem_type$base == "int",
         e$elem_type$base %in% c("int", "real") && (e$elem_type$nd %||% 0L) == 0L)
    },
    rowvec_expr = if (is_cx_ty(e$type)) {
      et <- if (e$type$base == "complex_matrix") ty("complex_row_vector") else ty("complex")
      items <- lapply(e$items, function(it) promote_to(cg_expr(it), it$type, et))
      cl(cx_rowvec, cl_list(as.name("list"), items), e$type$base == "complex_matrix")
    } else {
      items <- lapply(e$items, cg_expr)
      cl(rt_rowvec, cl_list(as.name("list"), items), e$type$base == "matrix")
    },
    call = if (!is.null(e$hof)) cg_hof(e) else if (!is.null(e$laplace)) cg_laplace(e) else cg_call(e),
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
  kinds <- vapply(e$idx, `[[`, "", "kind")
  if (length(kinds) == 1L && kinds == "single")
    return(cl(rt_index_s1, cg_expr(e$x), cg_expr(e$idx[[1L]]$e)))
  if (length(kinds) == 2L && all(kinds == "single"))
    return(cl(rt_index_s2, cg_expr(e$x), cg_expr(e$idx[[1L]]$e), cg_expr(e$idx[[2L]]$e)))
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
  # two real or int scalars, not both int
  rs <- !both_int && sca && scb
  if (rs && e$op %in% c("+", "-", "*", "/", ".*", "./"))
    return(cl(switch(e$op, "+" = op_sadd, "-" = op_ssub, "*" = , ".*" = op_smul, "/" = , "./" = op_sdiv), a, b))
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
  # arguments a complex signature takes as complex
  if (!is.null(sg$cxargs))
    args <- Map(function(code, a, p) if (p && !is_cx_ty(a$type)) promote_to(code, a$type, ty_cx(a$type)) else code,
                args, e$args, rep_len(sg$cxargs, length(args)))
  if (grepl("_(lupdf|lupmf)$", e$name)) {
    full <- lib_impl(sub("_lup", "_lp", e$name), types)
    return(cl(as.name("if"), cl(as.name("isTRUE"), cl(as.name("$"), CTX, as.name("propto"))),
              cg_size_guard(e$name, sg, types, cl_list(sg$impl, cg_keep_args(unname(args), e$args))),
              cg_size_guard(e$name, sg, types, cl_list(full, unname(args)))))
  }
  # library functions that touch the target (the _jacobian transforms)
  if (isTRUE(sg$ctx)) return(cl_list(sg$impl, c(unname(args), list(.ctx = CTX))))
  cg_size_guard(e$name, sg, types, cl_list(sg$impl, unname(args)))
}

# ---- statements ----------------------------------------------------------------------

# top level statements of a block, each preceded by recording its location
# (one assignment per statement, so runtime errors can say where they arose)
cg_block <- function(stmts, ctx) {
  out <- list(quote(.stmt$loc <- NULL))
  for (st in stmts) {
    if (!is.null(st$loc$line)) out <- c(out, list(cl(as.name("<-"), quote(.stmt$loc), st$loc[intersect(c("line", "col", "file"), names(st$loc))])))
    out <- c(out, cg_stmt(st, ctx))
  }
  out
}

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
    # complex values carry a trailing dimension holding the real and imaginary parts
    complex = 2L,
    complex_vector = , complex_row_vector = cl(as.name("c"), sz[[1L]], 2L),
    complex_matrix = cl(as.name("c"), sz[[1L]], sz[[2L]], 2L),
    cl(as.name("c"), sz[[1L]], sz[[2L]]))
  list(adims = ad, eshape = es)
}

cg_cons_args <- function(t) {
  if (identical(t$base, "tuple")) return(cl_list(as.name("list"), lapply(t$elems, cg_cons_args)))
  cl_list(as.name("list"), lapply(t$cons, cg_expr))
}

# counter for the temporaries of unpacking assignments
.unpack_n <- new.env(parent = emptyenv())
.unpack_n$n <- 0L

# assignment to a variable, possibly indexed or projected (a[i].2[j] = v):
# read the inner containers, assign into them, write the variable back
cg_lvalue_set <- function(lhs, rcode) {
  base <- lhs
  while (base$kind %in% c("index", "proj")) base <- base$x
  name <- base$name
  if (lhs$kind == "var") return(list(cl(as.name("<-"), vsym(name), cl(rt_assign, vsym(name), rcode, name))))
  chain <- list()
  cur <- lhs
  while (cur$kind %in% c("index", "proj")) {
    chain <- c(list(cur), chain)
    cur <- cur$x
  }
  tmp <- rcode
  for (lvl in rev(seq_along(chain))) {
    node <- chain[[lvl]]
    container <- if (lvl == 1L) vsym(name) else cg_expr(chain[[lvl - 1L]])
    tmp <- if (node$kind == "proj") cl(rt_set_proj, container, node$k, tmp, name) else {
      p <- cg_index_parts(node)
      cl(rt_assign_index, container, p$kinds, p$vals, tmp, name)
    }
  }
  list(cl(as.name("<-"), vsym(name), tmp))
}

cg_print_items <- function(items) {
  cl_list(as.name("list"), lapply(items, function(it)
    if (it$kind == "string") it$value
    else if (is_tuple_ty(it$type)) cl(fmt_typed, cg_expr(it), it$type)
    else if (is_cx_ty(it$type)) cl(cx_fmt, cg_expr(it), it$type$nd)
    else if (!is_scalar_ty(it$type)) cl(fmt_container, cg_expr(it))
    else cg_expr(it)))
}

# A density statement on data in the model block also records its pointwise
# terms when the context asks for them (bayescomp_log_lik() and
# bayescomp_posterior_predict()): the normalized density of each element
# and, when the distribution has one, a draw from its _rng.
#
# The observation unit comes from the signature: arguments declared `reals`
# or `ints` are elementwise and every element of the left side is one
# observation; other arguments (simplexes, cutpoints, covariance matrices)
# are shared by all of them. A left side that is not elementwise (a vector
# of a multivariate density, a matrix) is one joint observation. GLMs are
# rewritten to their univariate family with the linear predictor, so each
# row of x is one observation as in the equivalent unfused model.
cg_pointwise <- function(fname, lhs_e, arg_es, lhs, args, sg, truncated) {
  if (!identical(.cg$block, "model") || !isTRUE(lhs_e$data_only) || isTRUE(sg$user)) return(NULL)
  root <- lhs_e
  while (identical(root$kind, "index") || identical(root$kind, "paren")) root <- root$x %||% root$a
  label <- if (identical(root$kind, "var")) root$name else "y"
  types <- c(list(lhs_e$type), lapply(arg_es, `[[`, "type"))
  base <- sub("_(lpdf|lpmf|lupdf|lupmf)$", "", fname)
  suffix <- if (grepl("pmf$", fname)) "_lpmf" else "_lpdf"
  glm <- pw_glm[[base]]
  fam <- glm$family %||% base
  sig <- if (is.null(glm)) tryCatch(lib_sig(paste0(base, suffix), types), error = function(e) NULL) else .lib[[paste0(fam, suffix)]][[1L]]
  if (is.null(sig)) return(NULL)
  rng <- if (is.null(glm)) tryCatch(lib_impl(paste0(base, "_rng"), types[-1L]), error = function(e) NULL) else glm$rng %||% .lib[[paste0(fam, "_rng")]][[1L]]$impl
  vec <- glm$el %||% (sig$args %in% c("reals", "ints", "real", "int"))
  # elementwise when the signature vectorizes the left side and every
  # vectorized argument is a scalar or a one dimensional container (a GLM's
  # linear predictor always is). Other arguments are shared, or with one
  # vector per observation (cutpoints as array[N] vector) give row i to
  # observation i.
  if (is.null(glm)) {
    shape_ok <- vapply(types, function(t) (t$nd == 0L && t$base %in% c("int", "real", "vector", "row_vector")) ||
                                            (t$nd == 1L && t$base %in% c("int", "real")), TRUE)
    rows <- !vec & vapply(types, function(t) t$nd == 1L && t$base %in% c("vector", "row_vector"), TRUE)
    shared_ok <- vapply(types, function(t) t$nd == 0L, TRUE) | rows
    split <- vec[1L] && all(shape_ok[vec]) && all(shared_ok[!vec])
  } else {
    rows <- logical(length(vec))
    split <- TRUE
  }
  allvec <- all(sig$args %in% c("reals", "ints", "real", "int"))
  spec <- list(full = sig$impl, rng = rng, glm = glm$convert, el = if (split) vec, rows = rows,
               kern = if (split && allvec) pw_kernels[[fam]], vrng = split && allvec,
               truncated = truncated)
  cl(as.name("if"), cl(as.name("!"), cl(as.name("is.null"), cl(as.name("$"), CTX, as.name("pw")))),
     cl_list(rt_pointwise, c(list(CTX, label, spec, lhs), unname(args))))
}

lib_sig <- function(name, types) {
  for (sg in .lib[[name]]) {
    if (length(sg$args) == length(types) && all(mapply(pat_ok, sg$args, types))) return(sg)
  }
  stop("internal error: no library implementation for ", name)
}

# block being generated, read by cg_call() for _lupdf calls
.cg <- new.env(parent = emptyenv())

# Under propto Stan keeps a density term when an argument is an autodiff
# variable, which in the model block is every real valued expression that is
# not data only, also after floor() or an assignment of a constant. Such
# arguments are wrapped so the density sees them as differentiable.
cg_keep_args <- function(codes, exprs) {
  if (!identical(.cg$block, "model")) return(codes)
  Map(function(code, e) {
    t <- e$type
    if (isTRUE(e$data_only) || identical(t$base, "int") || is.null(t$base)) code else cl(rt_keep, code)
  }, codes, exprs)
}

cg_stmt <- function(st, ctx) {
  old <- .cg$block
  .cg$block <- ctx$block
  on.exit(.cg$block <- old)
  switch(st$kind,
    empty = list(),
    decl = if (identical(st$type$base, "tuple")) {
      desc <- if (st$top && ctx$block %in% c("transformed parameters", "generated quantities"))
        cl(as.name("[["), quote(.decls), st$name) else cg_decl_shape_tuple(st$type)
      out <- list(cl(as.name("<-"), vsym(st$name), cl(rt_decl_desc, desc)))
      if (!is.null(st$init)) {
        rhs <- promote_to(cg_expr(st$init), st$init$type, st$vtype)
        out <- c(out, list(cl(as.name("<-"), vsym(st$name), cl(rt_assign, vsym(st$name), rhs, st$name))))
      }
      out
    } else {
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
    block_stmt = {
      body <- cg_stmts(st$body, ctx)
      if (!is.null(st$profile)) {
        # profile("name") { }: time each execution
        tv <- as.name(sprintf(".prof_%d", prof_counter()))
        # the record is written however the block is left (return, break,
        # continue or an error): tryCatch evaluates the body in this frame
        body <- list(cl(as.name("<-"), tv, cl(rt_profile_start, st$profile)),
                     cl(quote(tryCatch), blk(body), finally = cl(rt_profile_end, st$profile, tv)))
      }
      list(blk(body))
    },
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
      cg_lvalue_set(lhs, rcode)
    },
    unpack = {
      .unpack_n$n <- .unpack_n$n + 1L
      tmp <- as.name(paste0(".unpack", .unpack_n$n))
      out <- list(cl(as.name("<-"), tmp, cg_expr(st$rhs)))
      set_items <- function(lhs, src, rtype) {
        unlist(lapply(seq_along(lhs$items), function(k) {
          it <- lhs$items[[k]]
          code <- cl(rt_proj, src, k)
          if (it$kind == "tuple_expr") return(set_items(it, code, rtype$elems[[k]]))
          cg_lvalue_set(it, promote_to(code, rtype$elems[[k]], it$type))
        }), recursive = FALSE)
      }
      c(out, set_items(st$lhs, tmp, st$rhs$type))
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
        kept <- cg_keep_args(c(list(lhs), unname(args)), c(list(st$lhs), st$args))
        cl(as.name("if"), cl(as.name("isTRUE"), cl(as.name("$"), CTX, as.name("propto"))),
           cg_size_guard(st$fname, sg, types, cl_list(sg$impl, kept)),
           cg_size_guard(st$fname, sg, types, cl_list(full, c(list(lhs), unname(args)))))
      }
      out <- list(cl(rt_target_add, CTX, call))
      pw <- cg_pointwise(st$fname, st$lhs, st$args, lhs, args, sg, truncated = !is.null(st$trunc))
      # after the statement, which checks the arguments first
      if (!is.null(pw)) out <- c(out, list(pw))
      if (!is.null(st$trunc)) {
        tr <- st$trunc
        mkfun <- function(sig) {
          if (is.null(sig)) return(NULL)
          # a user function also takes the evaluation context; the closure is
          # made where .ctx is in scope
          if (isTRUE(sig$user))
            return(cl(structure, cl_list(as.name("function"), list(as.pairlist(alist(... = )),
                                                                  cl_list(as.name(sig$mangled), list(as.name("..."), CTX)))),
                      user = TRUE))
          # the sampling statement has already checked the parameters
          attr(sig$impl, "unchecked") %||% sig$impl
        }
        lcdf_f <- mkfun(tr$lcdf_sig)
        lccdf_f <- mkfun(tr$lccdf_sig)
        out <- c(out, list(cl(rt_truncate, CTX, lhs,
                              if (is.null(tr$lo)) NULL else cg_expr(tr$lo),
                              if (is.null(tr$hi)) NULL else cg_expr(tr$hi),
                              lcdf_f, lccdf_f, cl_list(as.name("list"), unname(args)),
                              isTRUE(tr$discrete), tr$fam)))
      }
      out
    },
    target_add = {
      v <- st$value
      pw <- if (identical(v$kind, "call") && length(v$args) && grepl("_(lpdf|lpmf|lupdf|lupmf)$", v$name))
        cg_pointwise(v$name, v$args[[1L]], v$args[-1L], cg_expr(v$args[[1L]]), lapply(v$args[-1L], cg_expr), v$sig, FALSE)
      c(list(cl(rt_target_add, CTX, cg_expr(st$value))), if (!is.null(pw)) list(pw))
    },
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
    print = list(cl(rt_print, CTX, cg_print_items(st$items))),
    reject = list(cl(rt_reject, cg_print_items(st$items), FALSE)),
    fatal_error = list(cl(rt_reject, cg_print_items(st$items), TRUE)),
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

# ---- hoisting of data only subexpressions ----------------------------------------
# Expressions in the transformed parameters, model and generated quantities
# blocks that depend only on data are computed once when the data are bound
# (as stanc's optimizer moves them to transformed data) and read from a
# symbol `.H<k>` in each evaluation.

hoist_kinds <- c("call", "binop", "unop", "cond", "transpose", "index")

hoistable <- function(e) {
  if (!isTRUE(e$data_only) || !(e$kind %in% hoist_kinds)) return(FALSE)
  ok <- TRUE
  walk <- function(n) {
    if (!ok || !is.list(n)) return()
    k <- n$kind %||% ""
    if (k == "var" && !((n$origin %||% "") %in% c("data", "transformed data"))) { ok <<- FALSE; return() }
    if (k == "call") {
      s <- n$sig
      # calls that read the evaluation context (unnormalized densities drop
      # constants under propto, _jacobian calls depend on jacobian) stay
      if (isTRUE(s$user) || isTRUE(s$rng) || isTRUE(s$ctx) || grepl("_rng$|_lp$|_lupdf$|_lupmf$|_jacobian$", n$name) ||
          isTRUE(n$cond_bar && grepl("_lupdf$|_lupmf$", n$resolved %||% "")) || n$name %in% names(hof_specs) ||
          n$name %in% c("target", "get_lp")) { ok <<- FALSE; return() }
    }
    keys <- if (is.null(names(n))) seq_along(n) else setdiff(names(n), c("sig", "type", "loc"))
    for (nm in keys) if (is.list(n[[nm]])) walk(n[[nm]])
  }
  walk(e)
  ok
}

hoist_ir <- function(ir) {
  hoisted <- list()
  expr <- function(e) {
    if (hoistable(e)) {
      nm <- sprintf(".H%d", length(hoisted) + 1L)
      hoisted[[length(hoisted) + 1L]] <<- list(name = nm, expr = e)
      return(list(kind = "code", code = as.name(nm), type = e$type, data_only = TRUE, loc = e$loc))
    }
    node(e)
  }
  node <- function(n) {
    if (!is.list(n)) return(n)
    k <- n$kind %||% ""
    skip <- c("sig", "type", "vtype", "loc", "cons", if (k %in% c("assign", "unpack")) "lhs")
    keys <- if (is.null(names(n))) seq_along(n) else setdiff(names(n), skip)
    for (nm in keys) {
      v <- n[[nm]]
      if (!is.list(v) || is.function(v)) next
      n[[nm]] <- if (!is.null(v$kind) && !is.null(v$data_only)) expr(v) else node(v)
    }
    n
  }
  for (b in c("transformed parameters", "model", "generated quantities"))
    if (!is.null(ir$blocks[[b]])) ir$blocks[[b]]$body <- node(ir$blocks[[b]]$body)
  list(ir = ir, hoisted = hoisted)
}
