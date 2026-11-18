# Semantic analysis.
#
# Walks the untyped AST, resolves every identifier and function call,
# assigns a Stan type to every expression, and enforces Stan's rules about
# where statements and functions may appear. The result is the typed
# intermediate representation (IR) used by every execution backend: the AST
# with `type`, `data_only` and resolved `sig` fields added, plus symbol and
# block metadata. See ARCHITECTURE.md, "Intermediate representation".

sem_new <- function() {
  s <- new.env(parent = emptyenv())
  s$scopes <- list(list())
  s$ufuns <- list()          # name -> list of signatures
  s$block <- NULL
  s$fun <- NULL              # current function definition (when in functions block)
  s$loop_depth <- 0L
  s$decls <- list()          # top level declarations by block
  s
}

sem_push <- function(s) s$scopes <- c(list(list()), s$scopes)
sem_pop <- function(s) s$scopes <- s$scopes[-1L]
sem_lookup <- function(s, name) {
  for (sc in s$scopes) if (!is.null(sc[[name]])) return(sc[[name]])
  NULL
}
sem_declare <- function(s, name, info, loc) {
  if (!is.null(sem_lookup(s, name)))
    stop_semantic(sprintf("identifier '%s' is already in use", name), loc)
  if (name %in% stan_reserved_words)
    stop_semantic(sprintf("'%s' is a reserved word", name), loc)
  if (endsWith(name, "__")) stop_semantic(sprintf("identifier '%s' cannot end in two underscores", name), loc)
  s$scopes[[1L]][[name]] <- info
}

origin_data <- c("data", "transformed data")

stan_check <- function(prog, external = character()) {
  s <- sem_new()
  ext <- list()
  b <- prog$blocks
  # functions: register signatures first (allows recursion and forward use)
  if (!is.null(b$functions)) {
    for (fd in b$functions$body) sem_register_fun(s, fd)
    for (fd in b$functions$body) if (!is.null(fd$body)) sem_fundef(s, fd)
    for (nm in names(s$ufuns)) for (sg in s$ufuns[[nm]]) if (!isTRUE(sg$defined)) {
      # a declaration without a body may be implemented outside Stan (as
      # with stanc's --allow-undefined), here by an R function
      if (!nm %in% external)
        stop_semantic(sprintf("function '%s' is declared but never defined; supply an R implementation with 'functions'", nm), sg$loc)
      ext[[length(ext) + 1L]] <- list(name = nm, mangled = sg$mangled)
    }
  }
  out <- prog
  out$externals <- ext
  if (!is.null(b$functions)) out$blocks$functions <- b$functions
  for (bn in setdiff(stan_block_names, "functions")) {
    if (is.null(b[[bn]])) next
    s$block <- bn
    body <- b[[bn]]$body
    newbody <- vector("list", length(body))
    # model block variables are local: not visible in generated quantities
    if (bn == "model") sem_push(s)
    for (i in seq_along(body)) newbody[[i]] <- sem_stmt(s, body[[i]], top = TRUE)
    if (bn == "model") sem_pop(s)
    out$blocks[[bn]]$body <- newbody
  }
  if (!is.null(b$functions)) out$blocks$functions$body <- s$typed_funs
  out$ufuns <- s$ufuns
  out
}

# ---- functions ---------------------------------------------------------------------

utype <- function(ut) if (identical(ut$base, "tuple")) ty_tuple(lapply(ut$elems, utype), ut$nd) else ty(ut$base, ut$nd)

sem_register_fun <- function(s, fd) {
  argt <- lapply(fd$args, function(a) utype(a$type))
  rt <- utype(fd$rtype)
  nm <- fd$name
  # a library function may be overloaded, but only with a new signature
  for (lg in .lib[[nm]]) if (length(lg$args) == length(argt) && all(mapply(pat_exact, lg$args, argt)))
    stop_semantic(sprintf("function '%s' with these argument types is already defined in the Stan library", nm), fd$loc)
  if (grepl("_(lpdf|lpmf)$", nm)) {
    if (length(argt) == 0L) stop_semantic(sprintf("probability function '%s' needs at least one argument", nm), fd$loc)
    if (!ty_eq(rt, T_real)) stop_semantic(sprintf("probability function '%s' must return real", nm), fd$loc)
    if (endsWith(nm, "_lpmf") && argt[[1L]]$base != "int") stop_semantic(sprintf("the first argument of '%s' must be an integer type (_lpmf)", nm), fd$loc)
    if (endsWith(nm, "_lpdf") && argt[[1L]]$base == "int") stop_semantic(sprintf("the first argument of '%s' must not be an integer type (_lpdf); use _lpmf", nm), fd$loc)
  }
  existing <- s$ufuns[[nm]]
  for (k in seq_along(existing)) {
    eg <- existing[[k]]
    if (length(eg$args) == length(argt) && all(mapply(ty_eq, eg$args, argt))) {
      if (!ty_eq(eg$ret, rt)) stop_semantic(sprintf("function '%s' redeclared with a different return type", nm), fd$loc)
      if (!is.null(fd$body)) {
        if (isTRUE(eg$defined)) stop_semantic(sprintf("function '%s' is defined twice", nm), fd$loc)
        s$ufuns[[nm]][[k]]$defined <- TRUE
      }
      return(invisible())
    }
  }
  sg <- list(name = nm, args = argt, ret = rt, data_only = vapply(fd$args, function(a) a$data_only, TRUE),
             defined = !is.null(fd$body), loc = fd$loc, user = TRUE,
             index = length(existing) + 1L, rng = endsWith(nm, "_rng"), lp = endsWith(nm, "_lp"),
             jacobian = endsWith(nm, "_jacobian"))
  sg$mangled <- sprintf("F.%s.%d", nm, sg$index)
  s$ufuns[[nm]] <- c(existing, list(sg))
}

sem_fundef <- function(s, fd) {
  argt <- lapply(fd$args, function(a) utype(a$type))
  idx <- NULL
  for (k in seq_along(s$ufuns[[fd$name]])) {
    eg <- s$ufuns[[fd$name]][[k]]
    if (length(eg$args) == length(argt) && all(mapply(ty_eq, eg$args, argt))) idx <- k
  }
  sg <- s$ufuns[[fd$name]][[idx]]
  s$block <- "functions"
  s$fun <- sg
  sem_push(s)
  for (i in seq_along(fd$args)) {
    a <- fd$args[[i]]
    sem_declare(s, a$name, list(type = argt[[i]], origin = "argument", readonly = TRUE,
                                data_only = a$data_only), a$loc)
  }
  body <- sem_stmt(s, fd$body)
  sem_pop(s)
  if (!ty_eq(sg$ret, T_void) && !returns_all_paths(body))
    stop_semantic(sprintf("function '%s' may reach the end of its body without returning a value", fd$name), fd$loc)
  s$fun <- NULL
  s$block <- NULL
  td <- fd
  td$body <- body
  td$sig <- sg
  s$typed_funs <- c(s$typed_funs, list(td))
}

returns_all_paths <- function(st) {
  if (is.null(st)) return(FALSE)
  switch(st$kind,
    "return" = TRUE,
    "reject" = TRUE, "fatal_error" = TRUE,
    "block_stmt" = length(st$body) > 0L && any(vapply(st$body, returns_all_paths, TRUE)),
    "if" = !is.null(st$els) && returns_all_paths(st$then) && returns_all_paths(st$els),
    "while" = FALSE, "for" = FALSE, "foreach" = FALSE,
    FALSE)
}

# ---- declarations --------------------------------------------------------------------

constrained_bases <- c("simplex", "ordered", "positive_ordered", "unit_vector", "sum_to_zero_vector",
                       "sum_to_zero_matrix", "cholesky_factor_corr", "cholesky_factor_cov",
                       "corr_matrix", "cov_matrix", "column_stochastic_matrix", "row_stochastic_matrix")

sem_decl <- function(s, d, top) {
  t <- d$type
  bn <- s$block
  is_local <- !top || bn %in% c("model", "functions")
  if (is_local && (length(t$cons) > 0L || t$base %in% constrained_bases))
    stop_semantic(sprintf("local variable '%s' cannot have a constrained type", d$name), d$loc)
  if (bn %in% c("data", "parameters") && !is.null(d$init))
    stop_semantic(sprintf("variables in the %s block cannot be assigned in their declaration", bn), d$loc)
  if (bn == "parameters" && t$base == "int") stop_semantic("parameters cannot be integers", d$loc)
  if (bn %in% c("transformed parameters") && top && t$base == "int") stop_semantic("transformed parameters cannot be integers", d$loc)
  if (t$base == "int" && (!is.null(t$cons$offset) || !is.null(t$cons$multiplier)))
    stop_semantic("integer variables cannot have offset or multiplier", d$loc)
  if (identical(t$base, "tuple")) {
    if (is_local && tuple_constrained(t))
      stop_semantic(sprintf("local variable '%s' cannot have a constrained type", d$name), d$loc)
    r <- sem_decl_type(s, t, top && !is_local, bn, d$loc)
    vt <- r$vt
    d$type <- r$t
    init <- NULL
    if (!is.null(d$init)) {
      init <- sem_expr(s, d$init)
      if (!ty_assignable(vt, init$type))
        stop_semantic(sprintf("cannot assign a value of type %s to variable '%s' of type %s", ty_str(init$type), d$name, ty_str(vt)), init$loc)
    }
    origin <- if (is_local) "local" else bn
    sem_declare(s, d$name, list(type = vt, origin = origin, readonly = bn %in% c("data", "parameters"),
                                data_only = bn %in% origin_data || (is_local && bn %in% origin_data),
                                block = bn, decl_block = bn), d$loc)
    d$vtype <- vt
    d$init <- init
    d$top <- top && !is_local
    return(d)
  }
  r <- sem_decl_type(s, t, top && !is_local, bn, d$loc)
  t <- r$t
  vt <- r$vt
  init <- NULL
  if (!is.null(d$init)) {
    init <- sem_expr(s, d$init)
    if (!ty_assignable(vt, init$type))
      stop_semantic(sprintf("cannot assign a value of type %s to variable '%s' of type %s", ty_str(init$type), d$name, ty_str(vt)), init$loc)
  }
  origin <- if (is_local) "local" else bn
  sem_declare(s, d$name, list(type = vt, origin = origin, readonly = bn %in% c("data", "parameters"),
                              data_only = bn %in% origin_data || (is_local && bn %in% origin_data),
                              block = bn, decl_block = bn), d$loc)
  d$type <- t
  d$vtype <- vt
  d$init <- init
  d$top <- top && !is_local
  d
}

tuple_constrained <- function(t) {
  if (!identical(t$base, "tuple")) return(length(t$cons) > 0L || t$base %in% constrained_bases)
  any(vapply(t$elems, tuple_constrained, TRUE))
}

# sizes and bounds of a sized type; returns the checked type and its value type
sem_decl_type <- function(s, t, need_data, bn, loc) {
  if (bn == "parameters" && t$base == "int") stop_semantic("parameters cannot be integers", loc)
  if (identical(t$base, "tuple")) {
    t$adims <- lapply(t$adims, function(e) sem_size_expr(s, e, need_data))
    rs <- lapply(t$elems, function(et) sem_decl_type(s, et, need_data, bn, loc))
    t$elems <- lapply(rs, `[[`, "t")
    return(list(t = t, vt = ty_tuple(lapply(rs, `[[`, "vt"), length(t$adims))))
  }
  if (t$base == "int" && (!is.null(t$cons$offset) || !is.null(t$cons$multiplier)))
    stop_semantic("integer variables cannot have offset or multiplier", loc)
  if (bn == "transformed parameters" && need_data && t$base == "int") stop_semantic("transformed parameters cannot be integers", loc)
  # sizes must be integer expressions; at top level they must be data only
  t$adims <- lapply(t$adims, function(e) sem_size_expr(s, e, need_data))
  t$sizes <- lapply(t$sizes, function(e) sem_size_expr(s, e, need_data))
  for (cn in names(t$cons)) {
    ce <- sem_expr(s, t$cons[[cn]])
    db <- decl_base_type(t$base)
    shape_ok <- is_scalar_ty(ce$type) || (ce$type$base == db && ce$type$nd %in% c(0L, length(t$adims)))
    if (!shape_ok) stop_semantic(sprintf("bound '%s' has type %s, which does not match the declared type", cn, ty_str(ce$type)), ce$loc)
    if (t$base == "int" && ce$type$base != "int") stop_semantic(sprintf("bound '%s' of an integer variable must be an integer", cn), ce$loc)
    if (bn == "parameters" || bn == "transformed parameters" || bn == "generated quantities" || bn == "data" || bn == "transformed data") {
      # bounds of data may only depend on data
      if (bn %in% origin_data && !ce$data_only) stop_semantic(sprintf("bound '%s' must only depend on data", cn), ce$loc)
    }
    t$cons[[cn]] <- ce
  }
  list(t = t, vt = ty(decl_base_type(t$base), length(t$adims)))
}

sem_size_expr <- function(s, e, need_data) {
  te <- sem_expr(s, e)
  if (!is_int_ty(te$type)) stop_semantic(sprintf("size must be an integer, found %s", ty_str(te$type)), te$loc)
  # integer expressions count as data elsewhere, so the variables read are
  # checked too (to_int(floor(theta)) is an integer that is not data)
  if (need_data && (!te$data_only || reads_nondata(te))) stop_semantic("sizes of top level variables must only depend on data", te$loc)
  te
}
reads_nondata <- function(e) {
  if (!is.list(e)) return(FALSE)
  if (identical(e$kind, "var") && isTRUE(e$origin %in% c("parameters", "transformed parameters", "generated quantities", "model"))) return(TRUE)
  for (k in names(e) %||% seq_along(e)) if (!k %in% c("loc", "type", "sig") && is.list(e[[k]]) && reads_nondata(e[[k]])) return(TRUE)
  FALSE
}

# ---- statements ---------------------------------------------------------------------

in_lp_context <- function(s) {
  s$block == "model" || (!is.null(s$fun) && (isTRUE(s$fun$lp)))
}

sem_stmt <- function(s, st, top = FALSE) {
  k <- st$kind
  bn <- s$block
  if (k == "decl") return(sem_decl(s, st, top))
  if (bn == "data" || bn == "parameters") stop_semantic(sprintf("only declarations are allowed in the %s block", bn), st$loc)
  switch(k,
    empty = st,
    block_stmt = {
      sem_push(s)
      st$body <- lapply(st$body, function(x) sem_stmt(s, x))
      sem_pop(s)
      st
    },
    assign = {
      lhs <- sem_lvalue(s, st$lhs)
      rhs <- sem_expr(s, st$rhs)
      if (st$op != "=") {
        bop <- substr(st$op, 1L, nchar(st$op) - 1L)
        rt <- binop_type(bop, lhs$type, rhs$type)
        if (is.null(rt)) stop_semantic(sprintf("operator %s is not defined for %s and %s", st$op, ty_str(lhs$type), ty_str(rhs$type)), st$loc)
        if (!ty_assignable(lhs$type, rt)) stop_semantic(sprintf("result of %s has type %s, which cannot be assigned to %s", st$op, ty_str(rt), ty_str(lhs$type)), st$loc)
      } else if (!ty_assignable(lhs$type, rhs$type)) {
        stop_semantic(sprintf("cannot assign a value of type %s to a variable of type %s", ty_str(rhs$type), ty_str(lhs$type)), st$loc)
      }
      st$lhs <- lhs
      st$rhs <- rhs
      st
    },
    unpack = {
      rhs <- sem_expr(s, st$rhs)
      lhs <- sem_unpack_lhs(s, st$lhs)
      if (!ty_assignable(lhs$type, rhs$type))
        stop_semantic(sprintf("cannot unpack a value of type %s into %s", ty_str(rhs$type), ty_str(lhs$type)), st$loc)
      st$lhs <- lhs
      st$rhs <- rhs
      st
    },
    tilde = {
      if (!in_lp_context(s)) stop_semantic("sampling statements (~) are only allowed in the model block or in _lp functions", st$loc)
      if (paste0(st$dist, "_lpmf") %in% laplace_names && is.null(st$trunc))
        return(sem_stmt(s, list(kind = "target_add", loc = st$loc,
                                value = list(kind = "call", name = paste0(st$dist, "_lpmf"), args = c(list(st$lhs), st$args),
                                             cond_bar = TRUE, loc = st$loc))))
      lhs <- sem_expr(s, st$lhs)
      args <- lapply(st$args, function(a) sem_expr(s, a))
      types <- c(list(lhs$type), lapply(args, `[[`, "type"))
      res <- resolve_dist(s, st$dist, types, lhs$type$base == "int", st$loc)
      st$lhs <- lhs
      st$args <- args
      st$fname <- res$fname
      st$sig <- res$sig
      # T[ , ] has no bounds and means no truncation
      if (!is.null(st$trunc) && is.null(st$trunc$lo) && is.null(st$trunc$hi)) st$trunc <- NULL
      if (!is.null(st$trunc)) {
        if (!is_scalar_ty(lhs$type) && !(lhs$type$base %in% c("int", "real") && lhs$type$nd <= 1L) && !(lhs$type$base %in% c("vector", "row_vector")))
          stop_semantic("truncation needs a univariate distribution", st$loc)
        fam <- st$dist
        for (bnd in c("lo", "hi")) {
          e <- st$trunc[[bnd]]
          if (is.null(e)) next
          te <- sem_expr(s, e)
          if (!(te$type$base %in% c("int", "real")) || te$type$nd != 0L)
            stop_semantic(sprintf("truncation bound must be of type int or real, found %s", ty_str(te$type)), e$loc %||% st$loc)
          st$trunc[[bnd]] <- te
        }
        need <- c(if (!is.null(st$trunc$lo)) (if (!is.null(st$trunc$hi)) "lcdf" else "lccdf"), if (!is.null(st$trunc$hi)) "lcdf")
        for (nd in unique(need)) {
          fn <- paste0(fam, "_", nd)
          bt <- if (nd == "lccdf") st$trunc$lo$type else st$trunc$hi$type
          rr <- resolve_call(s, fn, c(list(bt), lapply(args, `[[`, "type")), quiet = TRUE)
          if (is.null(rr))
            stop_semantic(sprintf("truncation requires %s, which is not available for these arguments", fn), st$loc)
          st$trunc[[paste0(nd, "_sig")]] <- rr$sig
        }
        if (!is.null(st$trunc$lo) && !is.null(st$trunc$hi)) {
          # optional: two sided truncation in the upper tail uses the
          # difference of upper tail probabilities
          rr <- resolve_call(s, paste0(fam, "_lccdf"), c(list(st$trunc$lo$type), lapply(args, `[[`, "type")), quiet = TRUE)
          if (!is.null(rr)) st$trunc$lccdf_sig <- rr$sig
        }
        st$trunc$fam <- fam
        st$trunc$discrete <- lhs$type$base == "int"
      }
      st
    },
    target_add = {
      if (!in_lp_context(s)) stop_semantic("'target +=' is only allowed in the model block or in _lp functions", st$loc)
      st$value <- sem_expr(s, st$value)
      if (!(st$value$type$base %in% c("int", "real", "vector", "row_vector", "matrix")))
        stop_semantic("'target +=' needs a real valued expression", st$loc)
      st
    },
    jacobian_add = {
      ok <- bn == "transformed parameters" || (!is.null(s$fun) && isTRUE(s$fun$jacobian))
      if (!ok) stop_semantic("'jacobian +=' is only allowed in the transformed parameters block or in _jacobian functions", st$loc)
      st$value <- sem_expr(s, st$value)
      st
    },
    "if" = {
      st$cond <- sem_cond(s, st$cond)
      sem_push(s); st$then <- sem_stmt(s, st$then); sem_pop(s)
      if (!is.null(st$els)) { sem_push(s); st$els <- sem_stmt(s, st$els); sem_pop(s) }
      st
    },
    "while" = {
      st$cond <- sem_cond(s, st$cond)
      s$loop_depth <- s$loop_depth + 1L
      sem_push(s); st$body <- sem_stmt(s, st$body); sem_pop(s)
      s$loop_depth <- s$loop_depth - 1L
      st
    },
    "for" = {
      st$lo <- sem_expr(s, st$lo)
      st$hi <- sem_expr(s, st$hi)
      if (!is_int_ty(st$lo$type) || !is_int_ty(st$hi$type)) stop_semantic("for loop bounds must be integers", st$loc)
      s$loop_depth <- s$loop_depth + 1L
      sem_push(s)
      sem_declare(s, st$var, list(type = T_int, origin = "loop", readonly = TRUE, data_only = TRUE), st$loc)
      st$body <- sem_stmt(s, st$body)
      sem_pop(s)
      s$loop_depth <- s$loop_depth - 1L
      st
    },
    foreach = {
      st$iter <- sem_expr(s, st$iter)
      it <- st$iter$type
      et <- if (it$nd > 0L) ty_nd(it, it$nd - 1L) else switch(it$base, vector = T_real, row_vector = T_real, matrix = T_real,
        complex_vector = , complex_row_vector = , complex_matrix = ty("complex"), NULL)
      if (is.null(et)) stop_semantic("foreach loops need an array, vector, row_vector or matrix", st$loc)
      st$elem_type <- et
      s$loop_depth <- s$loop_depth + 1L
      sem_push(s)
      sem_declare(s, st$var, list(type = et, origin = "loop", readonly = TRUE, data_only = st$iter$data_only), st$loc)
      st$body <- sem_stmt(s, st$body)
      sem_pop(s)
      s$loop_depth <- s$loop_depth - 1L
      st
    },
    "break" = , "continue" = {
      if (s$loop_depth == 0L) stop_semantic(sprintf("'%s' outside of a loop", k), st$loc)
      st
    },
    "return" = {
      if (is.null(s$fun)) stop_semantic("'return' is only allowed inside function bodies", st$loc)
      rt <- s$fun$ret
      if (is.null(st$value)) {
        if (!ty_eq(rt, T_void)) stop_semantic("missing return value in a non void function", st$loc)
      } else {
        if (ty_eq(rt, T_void)) stop_semantic("void functions cannot return a value", st$loc)
        st$value <- sem_expr(s, st$value)
        if (!ty_assignable(rt, st$value$type)) stop_semantic(sprintf("returned type %s does not match the declared return type %s", ty_str(st$value$type), ty_str(rt)), st$loc)
      }
      st$rtype <- rt
      st
    },
    print = , reject = , fatal_error = {
      st$items <- lapply(st$items, function(it) if (it$kind == "string") it else sem_expr(s, it))
      st
    },
    call_stmt = {
      ce <- sem_expr(s, st$call, allow_void = TRUE)
      if (!ty_eq(ce$type, T_void)) stop_semantic(sprintf("non void function '%s' used as a statement; its value would be discarded", st$call$name), st$loc)
      st$call <- ce
      st
    },
    stop_semantic(sprintf("unexpected statement '%s'", k), st$loc))
}

sem_cond <- function(s, e) {
  te <- sem_expr(s, e)
  if (!is_scalar_ty(te$type)) stop_semantic(sprintf("condition must be int or real, found %s", ty_str(te$type)), te$loc)
  te
}

sem_unpack_lhs <- function(s, e) {
  if (e$kind != "tuple_expr") return(sem_lvalue(s, e))
  e$items <- lapply(e$items, function(it) sem_unpack_lhs(s, it))
  e$type <- ty_tuple(lapply(e$items, `[[`, "type"))
  e
}

sem_lvalue <- function(s, e) {
  base <- e
  while (base$kind %in% c("index", "proj")) base <- base$x
  if (base$kind != "var") stop_semantic("invalid assignment target", e$loc)
  info <- sem_lookup(s, base$name)
  if (is.null(info)) stop_semantic(sprintf("identifier '%s' not in scope", base$name), base$loc)
  bn <- s$block
  assignable <- !isTRUE(info$readonly) &&
    (info$origin == "local" || (info$origin == bn) || (bn == "functions" && info$origin == "local"))
  if (info$origin %in% c("loop", "argument")) assignable <- FALSE
  if (!assignable) stop_semantic(sprintf("cannot assign to variable '%s' (declared in %s) in this block", base$name, info$origin), e$loc)
  if (e$kind == "index") {
    te <- sem_expr(s, e)
    for (ix in te$idx) if (ix$kind != "single") te$multi_assign <- TRUE
    return(te)
  }
  sem_expr(s, e)
}

# ---- expressions ---------------------------------------------------------------------

sem_expr <- function(s, e, allow_void = FALSE) {
  k <- e$kind
  out <- switch(k,
    int = { e$type <- T_int; e$data_only <- TRUE; e },
    real = { e$type <- T_real; e$data_only <- TRUE; e },
    imag = { e$type <- ty("complex"); e$data_only <- TRUE; e },
    tuple_expr = {
      items <- lapply(e$items, function(x) sem_expr(s, x))
      e$items <- items
      e$type <- ty_tuple(lapply(items, `[[`, "type"))
      e$data_only <- all(vapply(items, function(x) isTRUE(x$data_only), TRUE))
      e
    },
    proj = {
      x <- sem_expr(s, e$x)
      t <- x$type
      if (!is_tuple_ty(t) || t$nd > 0L) stop_semantic(sprintf("'.%d' needs a tuple, found %s", e$k, ty_str(t)), e$loc)
      if (e$k < 1L || e$k > length(t$elems)) stop_semantic(sprintf("tuple of %d elements has no element %d", length(t$elems), e$k), e$loc)
      e$x <- x
      e$type <- t$elems[[e$k]]
      e$data_only <- x$data_only
      e
    },
    var = {
      info <- sem_lookup(s, e$name)
      if (is.null(info)) {
        if (!is.null(.lib[[e$name]]) || !is.null(s$ufuns[[e$name]]))
          stop_semantic(sprintf("'%s' is a function; call it with parentheses", e$name), e$loc)
        stop_semantic(sprintf("identifier '%s' not in scope", e$name), e$loc)
      }
      e$type <- info$type
      e$data_only <- isTRUE(info$data_only)
      e$origin <- info$origin
      e
    },
    paren = {
      a <- sem_expr(s, e$a)
      a
    },
    unop = {
      a <- sem_expr(s, e$a)
      t <- a$type
      if (e$op == "!") {
        if (!is_scalar_ty(t)) stop_semantic("'!' needs an int or real operand", e$loc)
        e$type <- T_int
      } else {
        if (!(t$base %in% c("int", "real", "vector", "row_vector", "matrix", names(cx_to_real))))
          stop_semantic(sprintf("unary %s is not defined for %s", e$op, ty_str(t)), e$loc)
        e$type <- t
      }
      e$a <- a
      e$data_only <- a$data_only
      e
    },
    binop = {
      a <- sem_expr(s, e$a)
      b <- sem_expr(s, e$b)
      rt <- binop_type(e$op, a$type, b$type)
      if (is.null(rt)) stop_semantic(sprintf("operator %s is not defined for operands of type %s and %s", e$op, ty_str(a$type), ty_str(b$type)), e$loc)
      e$a <- a
      e$b <- b
      e$type <- rt
      e$data_only <- a$data_only && b$data_only
      e
    },
    cond = {
      c0 <- sem_cond(s, e$cond)
      a <- sem_expr(s, e$a)
      b <- sem_expr(s, e$b)
      rt <- if (ty_eq(a$type, b$type)) a$type else if (ty_assignable(a$type, b$type)) a$type else if (ty_assignable(b$type, a$type)) b$type else NULL
      if (is.null(rt)) stop_semantic(sprintf("branches of a conditional expression have incompatible types %s and %s", ty_str(a$type), ty_str(b$type)), e$loc)
      e$cond <- c0; e$a <- a; e$b <- b
      e$type <- rt
      e$data_only <- c0$data_only && a$data_only && b$data_only
      e
    },
    transpose = {
      a <- sem_expr(s, e$a)
      t <- a$type
      if (t$nd > 0L || !(t$base %in% c("vector", "row_vector", "matrix", "complex_vector", "complex_row_vector", "complex_matrix")))
        stop_semantic(sprintf("transpose is not defined for %s", ty_str(t)), e$loc)
      e$type <- ty(switch(t$base, vector = "row_vector", row_vector = "vector", matrix = "matrix",
                          complex_vector = "complex_row_vector", complex_row_vector = "complex_vector",
                          complex_matrix = "complex_matrix"))
      e$a <- a
      e$data_only <- a$data_only
      e
    },
    index = sem_index(s, e),
    array_expr = {
      items <- lapply(e$items, function(x) sem_expr(s, x))
      ts <- lapply(items, `[[`, "type")
      et <- ts[[1L]]
      for (t in ts[-1L]) {
        if (ty_eq(et, t)) next
        if (ty_assignable(t, et)) et <- t else if (!ty_assignable(et, t))
          stop_semantic("array expression elements have incompatible types", e$loc)
      }
      e$items <- items
      e$elem_type <- et
      e$type <- ty_nd(et, et$nd + 1L)
      e$data_only <- all(vapply(items, `[[`, TRUE, "data_only"))
      e
    },
    rowvec_expr = {
      items <- lapply(e$items, function(x) sem_expr(s, x))
      ts <- lapply(items, `[[`, "type")
      is_cxs <- function(t) t$nd == 0L && t$base == "complex"
      is_rv <- function(t) t$nd == 0L && t$base %in% c("row_vector", "complex_row_vector")
      if (length(items) == 0L) {
        e$type <- T_rvec
      } else if (all(vapply(ts, is_scalar_ty, TRUE))) {
        e$type <- T_rvec
      } else if (all(vapply(ts, function(t) is_scalar_ty(t) || is_cxs(t), TRUE))) {
        e$type <- ty("complex_row_vector")
      } else if (all(vapply(ts, function(t) ty_eq(t, T_rvec), TRUE))) {
        e$type <- T_mat
      } else if (all(vapply(ts, is_rv, TRUE))) {
        e$type <- ty("complex_matrix")
      } else {
        stop_semantic("row vector expressions need all real or all row_vector elements", e$loc)
      }
      e$items <- items
      e$data_only <- all(vapply(items, function(x) isTRUE(x$data_only), TRUE))
      e
    },
    call = sem_call(s, e),
    string = stop_semantic("string literals are only allowed in print, reject and fatal_error", e$loc),
    stop_semantic(sprintf("unexpected expression '%s'", k), e$loc))
  if (!allow_void && !is.null(out$type) && ty_eq(out$type, T_void))
    stop_semantic("void function used in an expression", out$loc)
  # integers are never parameters, so integer expressions are data
  if (identical(out$type$base, "int")) out$data_only <- TRUE
  out
}

sem_index <- function(s, e) {
  x <- sem_expr(s, e$x)
  t <- x$type
  idx <- e$idx
  ndim_total <- ty_total_dims(t)
  if (length(idx) > ndim_total) stop_semantic(sprintf("too many indexes (%d) for a value of type %s", length(idx), ty_str(t)), e$loc)
  data_only <- x$data_only
  kinds <- character(length(idx))
  for (j in seq_along(idx)) {
    ix <- idx[[j]]
    if (ix$kind == "single") {
      te <- sem_expr(s, ix$e)
      if (is_int_ty(te$type)) {
        kinds[j] <- "single"
      } else if (te$type$base == "int" && te$type$nd == 1L) {
        kinds[j] <- "multi"
        ix$kind <- "intarray"
      } else {
        stop_semantic(sprintf("index must be an int or int array, found %s", ty_str(te$type)), te$loc)
      }
      ix$e <- te
      data_only <- data_only && te$data_only
    } else {
      kinds[j] <- "multi"
      for (fld in c("lo", "hi")) if (!is.null(ix[[fld]])) {
        te <- sem_expr(s, ix[[fld]])
        if (!is_int_ty(te$type)) stop_semantic("range bounds must be integers", te$loc)
        ix[[fld]] <- te
        data_only <- data_only && te$data_only
      }
    }
    idx[[j]] <- ix
  }
  # result type
  nd <- t$nd
  ka <- kinds[seq_len(min(nd, length(kinds)))]
  res_nd <- sum(ka == "multi") + max(0L, nd - length(kinds))
  ke <- if (length(kinds) > nd) kinds[(nd + 1L):length(kinds)] else character(0)
  base <- t$base
  cx <- is_cx_base(base)
  if (cx) base <- cx_to_real[[base]]
  res_base <- switch(base,
    int = "int", real = "real",
    vector = , row_vector = if (length(ke) >= 1L && ke[1L] == "single") "real" else base,
    matrix = {
      kr <- if (length(ke) >= 1L) ke[1L] else "multi"
      kc <- if (length(ke) >= 2L) ke[2L] else "multi"
      if (kr == "single" && kc == "single") "real"
      else if (kr == "single") "row_vector"
      else if (kc == "single") "vector"
      else "matrix"
    })
  if (cx) res_base <- real_to_cx[[res_base]]
  e$x <- x
  e$idx <- idx
  e$kinds <- kinds
  e$type <- if (is_tuple_ty(t)) ty_nd(t, res_nd) else ty(res_base, res_nd)
  e$data_only <- data_only
  e
}

# operator typing
binop_type <- function(op, a, b) {
  if (is_cx_ty(a) || is_cx_ty(b)) return(binop_type_cx(op, a, b))
  binop_type_real(op, a, b)
}

# complex operands: the real typing rules on the real parts, then complex
binop_type_cx <- function(op, a, b) {
  ok <- function(t) t$nd == 0L && t$base %in% c("int", "real", "vector", "row_vector", "matrix", names(cx_to_real))
  if (!ok(a) || !ok(b)) return(NULL)
  ra <- ty_realpart(a)
  rb <- ty_realpart(b)
  if (op %in% c("==", "!=")) return(if (ra$base %in% c("int", "real") && rb$base %in% c("int", "real")) T_int else NULL)
  if (!(op %in% c("+", "-", "*", "/", ".*", "./", "^", ".^"))) return(NULL)
  if (ra$base == "int") ra <- T_real
  if (rb$base == "int") rb <- T_real
  rt <- binop_type_real(op, ra, rb)
  if (is.null(rt)) NULL else ty_cx(rt)
}

binop_type_real <- function(op, a, b) {
  sc <- function(t) is_scalar_ty(t)
  arith_base <- c("vector", "row_vector", "matrix")
  isC <- function(t) t$nd == 0L && t$base %in% arith_base
  promote <- function(a, b) if (a$base == "int" && b$base == "int") T_int else T_real
  switch(op,
    "+" = , "-" = {
      if (sc(a) && sc(b)) return(promote(a, b))
      if (sc(a) && isC(b)) return(b)
      if (isC(a) && sc(b)) return(a)
      if (isC(a) && isC(b) && ty_eq(a, b)) return(a)
      NULL
    },
    "*" = {
      if (sc(a) && sc(b)) return(promote(a, b))
      if (sc(a) && isC(b)) return(b)
      if (isC(a) && sc(b)) return(a)
      if (!isC(a) || !isC(b)) return(NULL)
      key <- paste(a$base, b$base)
      switch(key,
        "matrix matrix" = T_mat, "matrix vector" = T_vec, "row_vector matrix" = T_rvec,
        "row_vector vector" = T_real, "vector row_vector" = T_mat, NULL)
    },
    "/" = {
      if (sc(a) && sc(b)) return(promote(a, b))
      if (isC(a) && sc(b)) return(a)
      if (sc(a) && isC(b)) return(b)
      if (isC(a) && isC(b)) {
        key <- paste(a$base, b$base)
        return(switch(key, "row_vector matrix" = T_rvec, "matrix matrix" = T_mat, NULL))
      }
      NULL
    },
    "%" = , "%/%" = if (is_int_ty(a) && is_int_ty(b)) T_int else NULL,
    "\\" = {
      if (!isC(a) || !isC(b)) return(NULL)
      switch(paste(a$base, b$base), "matrix matrix" = T_mat, "matrix vector" = T_vec, NULL)
    },
    ".*" = , "./" = {
      if (isC(a) && isC(b) && ty_eq(a, b)) return(a)
      if (isC(a) && sc(b)) return(a)
      if (sc(a) && isC(b)) return(b)
      if (sc(a) && sc(b)) return(T_real)
      NULL
    },
    # since Stan 2.30 ^ is elementwise like .^
    "^" = ,
    ".^" = {
      if (sc(a) && sc(b)) return(T_real)
      if (isC(a) && sc(b)) return(a)
      if (sc(a) && isC(b)) return(b)
      if (isC(a) && isC(b) && ty_eq(a, b)) return(a)
      if (a$nd > 0L || b$nd > 0L) {
        # elementwise power on arrays of reals
        if (sc(b) && a$base %in% c("int", "real")) return(ty("real", a$nd))
        if (ty_assignable(ty("real", b$nd), a) && ty_assignable(ty("real", a$nd), b)) return(ty("real", max(a$nd, b$nd)))
      }
      NULL
    },
    "==" = , "!=" = , "<" = , "<=" = , ">" = , ">=" = if (sc(a) && sc(b)) T_int else NULL,
    "&&" = , "||" = if (sc(a) && sc(b)) T_int else NULL,
    NULL)
}

# resolve a call against user functions and the library
# library calls must also match a signature of stanc (R/sysdata.rda, built
# by data-raw/stan_signatures.R), with Stan's int to real and real to
# complex promotions; the R implementations are often more general
sig_cache <- new.env(parent = emptyenv())
stan_sig_key <- function(name) sub("_lup([dm])f$", "_lp\\1f", name)
stan_sig_parsed <- function(key) {
  p <- sig_cache[[key]]
  if (is.null(p)) {
    p <- lapply(stan_signatures[[key]], function(a) {
      if (!nzchar(a)) return(list())
      lapply(strsplit(a, ",", fixed = TRUE)[[1L]], function(x) {
        b <- strsplit(x, ":", fixed = TRUE)[[1L]]
        ty(b[1L], as.integer(b[2L]))
      })
    })
    assign(key, p, envir = sig_cache)
  }
  p
}
stan_sig_ok <- function(name, types) {
  key <- stan_sig_key(name)
  # a library name stanc does not know is rejected like an unknown function
  if (is.null(stan_signatures[[key]])) return(FALSE)
  if (any(vapply(types, is_tuple_ty, TRUE))) return(TRUE)
  for (p in stan_sig_parsed(key))
    if (length(p) == length(types) && all(mapply(ty_assignable, p, types))) return(TRUE)
  FALSE
}
stan_sig_strings <- function(name) {
  vapply(stan_sig_parsed(stan_sig_key(name)), function(p) sprintf("  %s(%s)", name, paste(vapply(p, ty_str, ""), collapse = ", ")), "")
}

resolve_call <- function(s, name, types, quiet = FALSE, loc = NULL) {
  cands <- c(s$ufuns[[name]], .lib[[name]])
  if (!length(cands)) return(NULL)
  # user overloads and library signatures form one overload set; the unique
  # viable signature with the fewest promotions wins (Stan Reference Manual,
  # "Resolving overloads")
  best <- NULL; best_cost <- Inf; tie <- FALSE
  for (sg in cands) {
    if (!isTRUE(sg$user) || length(sg$args) != length(types)) next
    if (!all(mapply(ty_assignable, sg$args, types))) next
    cost <- sum(unlist(mapply(promo_cost, sg$args, types)))
    if (cost < best_cost) { best <- sg; best_cost <- cost; tie <- FALSE }
    else if (cost == best_cost) tie <- TRUE
  }
  lib <- NULL
  for (sg in cands) {
    if (isTRUE(sg$user) || length(sg$args) != length(types)) next
    if (all(mapply(pat_ok, sg$args, types))) {
      rt <- ret_type(sg, types)
      if (!is.null(rt)) { lib <- list(sig = sg, ret = rt); break }
    }
  }
  if (!is.null(best)) {
    if (!is.null(lib)) {
      # promotions of the closest signature in stanc's table
      lc <- vapply(Filter(function(p) length(p) == length(types) && all(mapply(ty_assignable, p, types)),
                          stan_sig_parsed(stan_sig_key(name))), function(p) sum(unlist(mapply(promo_cost, p, types))), 0)
      lc <- if (length(lc)) min(lc) else 0
      if (lc < best_cost) return(lib)
      if (lc == best_cost) tie <- TRUE
    }
    if (tie) stop_semantic(sprintf("ambiguous call to overloaded function '%s'", name), loc)
    return(list(sig = best, ret = best$ret))
  }
  if (!is.null(lib)) return(lib)
  # container functions of complex values reuse their real versions
  if (any(vapply(types, is_cx_ty, TRUE))) return(resolve_cx_lift(s, name, types))
  NULL
}

describe_sigs <- function(s, name) {
  cands <- c(s$ufuns[[name]], .lib[[name]])
  vapply(cands, function(sg) {
    a <- if (isTRUE(sg$user)) vapply(sg$args, ty_str, "") else sg$args
    sprintf("  %s(%s)", name, paste(a, collapse = ", "))
  }, "")
}

# Stan's functional names of the operators (add(a, b) is a + b, ...)
op_functions <- c(add = "+", subtract = "-", multiply = "*", divide = "/", elt_multiply = ".*",
                  elt_divide = "./", modulus = "%", logical_and = "&&", logical_or = "||",
                  logical_eq = "==", logical_neq = "!=", logical_lt = "<", logical_lte = "<=",
                  logical_gt = ">", logical_gte = ">=")
op_unary_functions <- c(minus = "-", plus = "+", logical_negation = "!")

sem_call <- function(s, e) {
  nm <- e$name
  if (is.null(s$ufuns[[nm]]) && !isTRUE(e$cond_bar)) {
    na <- length(e$args)
    if ((na == 2L && nm %in% names(op_functions)) || (na == 1L && nm %in% names(op_unary_functions))) {
      # the functional forms have fewer signatures than the operators
      at <- lapply(e$args, function(a) sem_expr(s, a)$type)
      if (!stan_sig_ok(nm, at))
        stop_semantic(sprintf("no matching signature for %s(%s); Stan %s defines:\n%s", nm, paste(vapply(at, ty_str, ""), collapse = ", "),
                              stan_signatures_version, paste(utils::head(stan_sig_strings(nm), 20L), collapse = "\n")), e$loc)
    }
    # fun_name records the call for bayescomp_info()
    if (na == 2L && nm %in% names(op_functions))
      return(sem_expr(s, list(kind = "binop", op = op_functions[[nm]], a = e$args[[1L]], b = e$args[[2L]], loc = e$loc, fun_name = nm)))
    if (na == 1L && nm %in% names(op_unary_functions))
      return(sem_expr(s, list(kind = "unop", op = op_unary_functions[[nm]], a = e$args[[1L]], loc = e$loc, fun_name = nm)))
    if (na == 1L && nm == "transpose")
      return(sem_expr(s, list(kind = "transpose", a = e$args[[1L]], loc = e$loc)))
  }
  # higher order functions (ODE solvers, reduce_sum, ...) take a function name
  # as an argument; report them before that argument fails to resolve
  if (nm %in% unsupported_functions && is.null(s$ufuns[[nm]]) && is.null(.lib[[nm]]))
    stop_unsupported(sprintf("function '%s' is not implemented in this version of bayescompR", nm), e$loc)
  if (nm %in% names(hof_specs) && is.null(s$ufuns[[nm]])) return(sem_hof(s, e))
  if (nm %in% laplace_names && is.null(s$ufuns[[nm]])) return(sem_laplace(s, e))
  args <- lapply(e$args, function(a) sem_expr(s, a))
  types <- lapply(args, `[[`, "type")
  if (nm == "target") {
    if (length(args)) stop_semantic("target() takes no arguments", e$loc)
    # reading target is also allowed in transformed parameters (stanc3),
    # incrementing it is not
    if (!in_lp_context(s) && !identical(s$block, "transformed parameters"))
      stop_semantic("target() is only allowed in the model block, transformed parameters or _lp functions", e$loc)
    e$type <- T_real
    e$data_only <- FALSE
    e$args <- args
    e$builtin <- "target"
    return(e)
  }
  if (nm %in% c("get_lp", "increment_log_prob", "integrate_ode", "integrate_ode_rk45", "integrate_ode_bdf", "integrate_ode_adams", "if_else", "multiply_log", "binomial_coefficient_log", "cov_exp_quad", "fabs")) {
    stop_semantic(sprintf("'%s' was removed from the Stan language", nm), e$loc)
  }
  if (e$cond_bar && !grepl("_(lpdf|lupdf|lpmf|lupmf|lcdf|lccdf|cdf)$", nm))
    stop_semantic(sprintf("'|' is only allowed in calls to probability functions, not in '%s'", nm), e$loc)
  if (!e$cond_bar && grepl("_(lpdf|lupdf|lpmf|lupmf)$", nm) && length(args) > 1L)
    stop_semantic(sprintf("probability density and mass functions need '|' after the first argument: %s(y | ...)", nm), e$loc)
  # unnormalized user densities: foo_lupdf resolves to the user foo_lpdf
  uname <- nm
  if (grepl("_lupdf$", nm) && !is.null(s$ufuns[[sub("_lupdf$", "_lpdf", nm)]])) uname <- sub("_lupdf$", "_lpdf", nm)
  if (grepl("_lupmf$", nm) && !is.null(s$ufuns[[sub("_lupmf$", "_lpmf", nm)]])) uname <- sub("_lupmf$", "_lpmf", nm)
  if (grepl("_(lupdf|lupmf)$", nm)) {
    ok <- s$block == "model" || (!is.null(s$fun) && grepl("_(lpdf|lpmf|lp)$", s$fun$name))
    if (!ok) stop_semantic(sprintf("unnormalized function '%s' is only allowed in the model block or in the body of _lpdf, _lpmf or _lp functions", nm), e$loc)
  }
  if (is.null(s$ufuns[[uname]]) && is.null(.lib[[uname]])) {
    # the transform functions of Stan 2.36 (lower_bound_jacobian, simplex_constrain, ...)
    if (nm %in% unsupported_functions) stop_unsupported(sprintf("function '%s' is not implemented in this version of bayescompR", nm), e$loc)
    stop_semantic(sprintf("unknown function '%s'", nm), e$loc)
  }
  res <- resolve_call(s, uname, types)
  if (is.null(res)) {
    stop_semantic(sprintf("no matching signature for %s(%s); available signatures:\n%s", nm,
                          paste(vapply(types, ty_str, ""), collapse = ", "),
                          paste(describe_sigs(s, uname), collapse = "\n")), e$loc)
  }
  sg <- res$sig
  if (!isTRUE(sg$user) && is.null(stan_signatures[[stan_sig_key(uname)]])) stop_semantic(sprintf("unknown function '%s'", nm), e$loc)
  if (!isTRUE(sg$user) && !stan_sig_ok(uname, types)) {
    stop_semantic(sprintf("no matching signature for %s(%s); Stan %s defines:\n%s", nm,
                          paste(vapply(types, ty_str, ""), collapse = ", "), stan_signatures_version,
                          paste(utils::head(stan_sig_strings(uname), 20L), collapse = "\n")), e$loc)
  }
  is_rng <- isTRUE(sg$rng) || endsWith(nm, "_rng")
  if (is_rng) {
    ok <- s$block %in% c("transformed data", "generated quantities") || (!is.null(s$fun) && endsWith(s$fun$name, "_rng"))
    if (!ok) stop_semantic(sprintf("random number generator '%s' is only allowed in transformed data, generated quantities and _rng functions", nm), e$loc)
  }
  if (endsWith(nm, "_lp")) {
    ok <- s$block %in% c("model", "transformed parameters") || (!is.null(s$fun) && endsWith(s$fun$name, "_lp"))
    if (!ok) stop_semantic(sprintf("'%s' (an _lp function) is only allowed in the model block, the transformed parameters block or other _lp functions", nm), e$loc)
  }
  if (endsWith(nm, "_jacobian")) {
    ok <- s$block == "transformed parameters" || (!is.null(s$fun) && endsWith(s$fun$name, "_jacobian"))
    if (!ok) stop_semantic(sprintf("'%s' (a _jacobian function) is only allowed in transformed parameters or other _jacobian functions", nm), e$loc)
  }
  if (isTRUE(sg$user)) {
    for (i in seq_along(args)) if (isTRUE(sg$data_only[i]) && !args[[i]]$data_only)
      stop_semantic(sprintf("argument %d of '%s' must be data only", i, nm), args[[i]]$loc)
  } else if (nm %in% data_only_functions && (s$block %in% c("transformed parameters", "model") || !is.null(s$fun))) {
    # in generated quantities every value is a double, so parameters qualify;
    # in a function body only arguments declared data (and values built from
    # them) do, since the function may be called with parameters
    for (i in seq_along(args)) if (!isTRUE(args[[i]]$data_only))
      stop_semantic(sprintf("argument %d of '%s' must be data only", i, nm), args[[i]]$loc)
  }
  e$args <- args
  e$sig <- sg
  e$type <- res$ret
  e$data_only <- all(vapply(args, function(a) isTRUE(a$data_only), TRUE)) && !is_rng
  e$resolved <- uname
  e
}

resolve_dist <- function(s, dist, types, int_lhs, loc) {
  tries <- if (int_lhs) c("_lupmf", "_lpmf", "_lupdf", "_lpdf") else c("_lupdf", "_lpdf", "_lupmf", "_lpmf")
  for (suf in tries) {
    fn <- paste0(dist, suf)
    ufn <- sub("_lup", "_lp", fn)
    if (!is.null(s$ufuns[[ufn]])) {
      r <- resolve_call(s, ufn, types)
      if (!is.null(r)) return(list(fname = ufn, sig = r$sig, propto = TRUE))
    }
    if (!is.null(.lib[[fn]])) {
      r <- resolve_call(s, fn, types)
      if (!is.null(r)) {
        if (is.null(stan_signatures[[stan_sig_key(fn)]])) stop_semantic(sprintf("unknown distribution '%s'", dist), loc)
        if (!stan_sig_ok(fn, types))
          stop_semantic(sprintf("no matching signature for %s(%s); Stan %s defines:\n%s", stan_sig_key(fn),
                                paste(vapply(types, ty_str, ""), collapse = ", "), stan_signatures_version,
                                paste(utils::head(stan_sig_strings(fn), 20L), collapse = "\n")), loc)
        return(list(fname = fn, sig = r$sig))
      }
    }
  }
  any_exists <- any(vapply(paste0(dist, c("_lpdf", "_lpmf")), function(f) !is.null(.lib[[f]]) || !is.null(s$ufuns[[f]]), TRUE))
  if (!any_exists) {
    if (dist %in% unsupported_distributions) stop_unsupported(sprintf("distribution '%s' is not implemented in this version of bayescompR", dist), loc)
    stop_semantic(sprintf("unknown distribution '%s'", dist), loc)
  }
  fam <- if (!is.null(.lib[[paste0(dist, "_lpdf")]]) || !is.null(s$ufuns[[paste0(dist, "_lpdf")]])) paste0(dist, "_lpdf") else paste0(dist, "_lpmf")
  stop_semantic(sprintf("no matching signature for %s(%s); available signatures:\n%s", fam,
                        paste(vapply(types, ty_str, ""), collapse = ", "),
                        paste(describe_sigs(s, fam), collapse = "\n")), loc)
}

# Stan functions recognized but not implemented in this version. Using them
# raises bayescomp_unsupported_error instead of an unknown function error.
# library functions whose arguments are all data only in Stan 2.40
data_only_functions <- c("linspaced_array", "linspaced_int_array", "linspaced_vector", "linspaced_row_vector",
                         "one_hot_array", "one_hot_int_array", "one_hot_vector", "one_hot_row_vector",
                         "zeros_array", "zeros_int_array", "zeros_vector", "zeros_row_vector",
                         "ones_array", "ones_int_array", "ones_vector", "ones_row_vector",
                         "uniform_simplex", "quantile", "to_int")

unsupported_functions <- character(0)
unsupported_distributions <- character(0)
