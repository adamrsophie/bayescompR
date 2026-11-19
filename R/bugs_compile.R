# Compiling a BUGS program into a graph of nodes.
#
# Every statement, repeated over its loop indices, defines a family of
# nodes: a stochastic family `x[i] ~ d(...)` or a logical family
# `x[i] <- e`. The compiler
#   1. unrolls the loops (their bounds must be constants: data or numbers),
#   2. evaluates the indices on the left sides, sizes the arrays and checks
#      that every element is defined once,
#   3. marks stochastic nodes whose elements all have data as observed,
#   4. finds the parents of every node; an index that depends on a node
#      (`mu[T[i]]`, a mixture) makes the node depend on that node and on
#      every element the index can reach,
#   5. sorts the nodes so that parents come first, and
#   6. turns each statement into an R function of its instance rows. When
#      every index is a single position and every function acts elementwise
#      the function evaluates all requested instances at once with
#      vectorized R; otherwise it loops over them.
#
# The values of data and nodes live in the environment `S` as arrays, so
# `x[i, j]` for vectors i and j is `S$x[cbind(i, j)]`.

bugs_compile <- function(prog, data = list(), dialect = "jags") {
  M <- new.env(parent = emptyenv())
  M$dialect <- dialect
  M$prog <- prog
  data <- bugs_check_values(data, "data")
  # JAGS data block: deterministic transformations of the data, computed once
  if (length(prog$data)) data <- bugs_data_block(prog$data, data, dialect)
  M$data <- data
  stmts <- list()
  flatten <- function(ss, loops) for (s in ss) {
    if (identical(s$kind, "for")) flatten(s$body, c(loops, list(s)))
    else {
      if (identical(s$kind, "derivative")) stop_unsupported("differential equations (D(...) statements) are not supported", s$loc)
      s$loops <- loops
      # one canonical (BUGS) parameter list, used by code generation and the
      # choice of samplers alike
      if (identical(s$kind, "stochastic") && !is.null(s$argnames) && any(nzchar(s$argnames))) s <- bugs_named_args(s)
      # NIMBLE's I() is truncation
      if (identical(s$kind, "stochastic") && identical(s$bound$type, "I") && dialect == "nimble") s$bound$type <- "T"
      if (identical(s$kind, "logical")) s$cut <- identical(s$expr$kind, "call") && identical(s$expr$name, "cut")
      cut_bad <- if (identical(s$kind, "logical")) bugs_has_cut(s$expr, top = TRUE)
                 else any(vapply(s$args, bugs_has_cut, TRUE))
      if (cut_bad) stop_unsupported("cut() must be the whole right side of a logical node, as in y <- cut(x)", s$loc)
      stmts[[length(stmts) + 1L]] <<- s
    }
  }
  flatten(prog$model, list())
  M$stmts <- stmts
  # 1. instances of each statement
  for (k in seq_along(stmts)) M$stmts[[k]]$inst <- bugs_unroll(M$stmts[[k]], data)
  # 2. left side cells and array sizes
  bugs_size_arrays(M, prog$var)
  # 3. nodes and code
  bugs_make_nodes(M)
  for (k in seq_along(M$stmts)) M$stmts[[k]]$fn <- bugs_codegen(M, k)
  # logical nodes computable from the data are constants (this covers
  # OpenBUGS data transformations), then observation status
  bugs_fold(M)
  bugs_observed(M)
  # 4, 5. parents and order
  bugs_find_parents(M)
  bugs_order(M)
  for (k in seq_along(M$stmts)) if (M$stmts[[k]]$kind == "stochastic")
    environment(M$stmts[[k]]$fn$lp)$cobs <- M$nodes$obs[M$stmts[[k]]$nodes]
  M
}

# ---- constant evaluation ----------------------------------------------------------------
# Loop bounds and left side indices are evaluated from the data and the
# loop variables; anything else is an error.

bugs_const <- function(e, env, data, what = "this expression") {
  ev <- function(e) {
    switch(e$kind,
      num = e$value,
      var = {
        if (is.null(e$index) && !is.null(env[[e$name]])) return(env[[e$name]])
        v <- data[[e$name]]
        if (is.null(v)) stop_semantic(sprintf("%s must be a constant, but '%s' is not data", what, e$name), e$loc)
        if (is.null(e$index)) {
          if (anyNA(v)) stop_semantic(sprintf("%s uses a missing value of '%s'", what, e$name), e$loc)
          return(v)
        }
        idx <- lapply(e$index, function(ix) switch(ix$kind, single = ev(ix$e), range = ev(ix$lo):ev(ix$hi), all = NULL))
        d <- dim(v) %||% length(v)
        if (length(idx) != length(d)) stop_semantic(sprintf("'%s' has %d dimensions but %d indices", e$name, length(d), length(idx)), e$loc)
        idx <- Map(function(ix, n) if (is.null(ix)) seq_len(n) else ix, idx, d)
        r <- do.call(`[`, c(list(if (is.null(dim(v))) array(v, length(v)) else v), idx, list(drop = TRUE)))
        if (anyNA(r)) stop_semantic(sprintf("%s uses a missing value of '%s'", what, e$name), e$loc)
        r
      },
      unop = if (e$op == "-") -ev(e$a) else as.numeric(!ev(e$a)),
      binop = bugs_apply_op(e$op, ev(e$a), ev(e$b)),
      call = {
        f <- bugs_fun_get(e$name, e$loc)
        do.call(f$r, lapply(e$args, ev))
      },
      stop_semantic(sprintf("%s must be a constant", what), e$loc))
  }
  ev(e)
}

bugs_apply_op <- function(op, a, b) {
  switch(op,
    "+" = a + b, "-" = a - b, "*" = a * b, "/" = a / b, "^" = a^b,
    "<" = as.numeric(a < b), "<=" = as.numeric(a <= b), ">" = as.numeric(a > b), ">=" = as.numeric(a >= b),
    "==" = as.numeric(a == b), "!=" = as.numeric(a != b),
    "&&" = , "&" = as.numeric(a & b), "||" = , "|" = as.numeric(a | b),
    "%*%" = a %*% b,
    stop_unsupported(sprintf("the operator '%s' is not supported", op)))
}

# the loop index values of every instance of a statement: a matrix with one
# column per enclosing loop (zero columns outside loops)
bugs_unroll <- function(s, data) {
  loops <- s$loops
  vars <- vapply(loops, `[[`, "", "var")
  rows <- list()
  rec <- function(k, env) {
    if (k > length(loops)) { rows[[length(rows) + 1L]] <<- unlist(env[vars]); return(invisible()) }
    lp <- loops[[k]]
    lo <- bugs_const(lp$lo, env, data, "a loop bound")
    hi <- bugs_const(lp$hi, env, data, "a loop bound")
    if (length(lo) != 1L || length(hi) != 1L) stop_semantic("loop bounds must be scalars", lp$loc)
    if (hi < lo) return(invisible())
    for (v in lo:hi) { env[[lp$var]] <- v; rec(k + 1L, env) }
  }
  rec(1L, list())
  if (!length(vars)) return(matrix(numeric(0), 1L, 0L))
  if (!length(rows)) return(matrix(numeric(0), 0L, length(vars), dimnames = list(NULL, vars)))
  m <- matrix(unlist(rows), ncol = length(vars), byrow = TRUE)
  colnames(m) <- vars
  m
}

bugs_inst_env <- function(s, r) {
  if (!ncol(s$inst)) return(list())
  as.list(s$inst[r, ])
}

# ---- arrays and nodes -------------------------------------------------------------------

# the index ranges of the left side of every instance: a list (one element
# per instance) of lists of index vectors
bugs_lhs_index <- function(s, data, dims = NULL) {
  idx <- s$lhs$index
  lapply(seq_len(nrow(s$inst)), function(r) {
    env <- bugs_inst_env(s, r)
    if (is.null(idx)) return(NULL)
    lapply(seq_along(idx), function(k) {
      ix <- idx[[k]]
      v <- switch(ix$kind,
        single = bugs_const(ix$e, env, data, "an index on the left side"),
        range = bugs_const(ix$lo, env, data, "an index on the left side"):bugs_const(ix$hi, env, data, "an index on the left side"),
        all = if (is.null(dims)) NA else seq_len(dims[k]))
      if (any(is.na(v)) && !identical(ix$kind, "all")) stop_semantic("an index on the left side is missing", s$loc)
      if (!identical(ix$kind, "all") && (any(v < 1) || any(v != round(v)))) stop_semantic(sprintf("index %s of '%s' on the left side is not a positive integer", paste(v, collapse = ","), s$lhs$name), s$loc)
      v
    })
  })
}

bugs_size_arrays <- function(M, decls) {
  data <- M$data
  dims <- list()
  for (d in decls) {
    env <- list()
    dims[[d$name]] <- if (!length(d$dims)) 1L else vapply(d$dims, function(ix) as.integer(bugs_const(ix$e, env, data, "a declared size")), 1L)
  }
  for (nm in names(data)) if (is.null(dims[[nm]])) dims[[nm]] <- dim(data[[nm]]) %||% length(data[[nm]])
  # sizes of variables defined only on the left side: the largest index
  for (k in seq_along(M$stmts)) {
    s <- M$stmts[[k]]
    nm <- s$lhs$name
    li <- bugs_lhs_index(s, data, dims[[nm]])
    M$stmts[[k]]$lhs_idx <- li
    if (!is.null(dims[[nm]]) && (nm %in% names(data) || nm %in% vapply(decls, `[[`, "", "name"))) next
    for (ix in li) {
      ext <- if (is.null(ix)) 1L else vapply(ix, function(v) as.integer(max(v)), 1L)
      if (any(is.na(ext))) stop_semantic(sprintf("the size of '%s' is unknown; give it in the data or a var declaration", nm), s$loc)
      cur <- dims[[nm]]
      if (is.null(cur)) dims[[nm]] <- ext
      else {
        if (length(cur) != length(ext)) stop_semantic(sprintf("'%s' is used with %d and %d indices", nm, length(cur), length(ext)), s$loc)
        dims[[nm]] <- pmax(cur, ext)
      }
    }
  }
  # left sides with empty indices need the size, now known
  for (k in seq_along(M$stmts)) {
    s <- M$stmts[[k]]
    if (any(vapply(s$lhs$index %||% list(), function(ix) identical(ix$kind, "all"), TRUE)))
      M$stmts[[k]]$lhs_idx <- bugs_lhs_index(s, data, dims[[s$lhs$name]])
  }
  M$dims <- dims
  S <- new.env(parent = emptyenv())
  for (nm in names(dims)) {
    a <- array(NA_real_, dims[[nm]])
    if (!is.null(data[[nm]])) {
      v <- data[[nm]]
      if (length(v) != length(a)) stop_data(sprintf("'%s' has %d values but the model needs %d", nm, length(v), length(a)))
      a[] <- as.numeric(v)
    }
    assign(nm, a, envir = S)
  }
  M$S <- S
}

# linear positions of index vectors in an array of dimensions d
bugs_cells <- function(ix, d) {
  if (is.null(ix)) return(1L)
  for (k in seq_along(ix)) if (any(ix[[k]] > d[k])) stop_semantic("an index is beyond the size of the array")
  g <- as.matrix(expand.grid(ix))
  as.integer(1L + colSums((t(g) - 1L) * cumprod(c(1L, d[-length(d)]))))
}

bugs_make_nodes <- function(M) {
  owner <- lapply(M$dims, function(d) array(0L, d))
  nodes <- list(stmt = integer(), row = integer(), var = character(), kind = character())
  cells <- list()
  n <- 0L
  transforms <- integer()
  for (k in seq_along(M$stmts)) {
    s <- M$stmts[[k]]
    nm <- s$lhs$name
    d <- M$dims[[nm]]
    ids <- integer(nrow(s$inst))
    for (r in seq_len(nrow(s$inst))) {
      cl <- bugs_cells(s$lhs_idx[[r]], d)
      prev <- unique(owner[[nm]][cl])
      prev <- prev[prev != 0L]
      if (length(prev)) {
        # OpenBUGS: a logical node that is also stochastic transforms data
        pk <- nodes$kind[prev[1L]]
        if (length(prev) == 1L && pk != s$kind && setequal(cells[[prev]], cl)) {
          transforms <- c(transforms, if (s$kind == "logical") n + 1L else prev)
        } else stop_semantic(sprintf("%s[%s] is defined more than once", nm, paste(arrayInd(cl[1L], d), collapse = ",")), s$loc)
      }
      n <- n + 1L
      nodes$stmt[n] <- k; nodes$row[n] <- r; nodes$var[n] <- nm; nodes$kind[n] <- s$kind
      cells[[n]] <- cl
      if (!length(prev)) owner[[nm]][cl] <- n
      ids[r] <- n
    }
    M$stmts[[k]]$nodes <- ids
  }
  nodes$obs <- logical(n)
  nodes$folded <- logical(n)
  M$owner <- owner
  M$nodes <- nodes
  M$cells <- cells
  M$transforms <- transforms
  M$fixed <- list()
}

# Logical nodes whose values follow from the data alone are computed once
# and become constants. Unknown values are NA at this point, so a node is
# folded when its value has no NA (R's NA rules make such a value hold for
# every value of the unknowns).
bugs_fold <- function(M) {
  L <- which(M$nodes$kind == "logical")
  repeat {
    changed <- FALSE
    for (i in L[!M$nodes$folded[L]]) {
      s <- M$stmts[[M$nodes$stmt[i]]]
      v <- tryCatch(suppressWarnings(s$fn$value(M$nodes$row[i])), error = function(e) NULL)
      if (is.list(v)) v <- v[[1L]]
      if (is.null(v) || anyNA(v) || length(v) != length(M$cells[[i]]) && length(v) != 1L) next
      cur <- get(M$nodes$var[i], envir = M$S)[M$cells[[i]]]
      if (!i %in% M$transforms && any(!is.na(cur)))
        stop_semantic(sprintf("%s is a logical node and cannot be given as data", bugs_node_label(M, i)), s$loc)
      bugs_set(M, M$nodes$var[i], M$cells[[i]], v)
      M$nodes$folded[i] <- TRUE
      changed <- TRUE
    }
    if (!changed) break
  }
}

# observed stochastic nodes: every element has a value. A multivariate
# normal may be partly observed: its missing elements are sampled.
bugs_observed <- function(M) {
  n <- length(M$nodes$kind)
  for (i in seq_len(n)) {
    v <- get(M$nodes$var[i], envir = M$S)[M$cells[[i]]]
    s <- M$stmts[[M$nodes$stmt[i]]]
    if (M$nodes$kind[i] == "stochastic") {
      if (all(!is.na(v))) M$nodes$obs[i] <- TRUE
      else if (any(!is.na(v))) {
        if (!s$dist %in% c("dmnorm", "dmnorm.vcov", "dmt"))
          stop_semantic(sprintf("%s is partly observed; only multivariate normal and t nodes may be", bugs_node_label(M, i)), s$loc)
        M$fixed[[as.character(i)]] <- !is.na(v)
      }
    } else if (!M$nodes$folded[i] && any(!is.na(v)) && !i %in% M$transforms) {
      stop_semantic(sprintf("%s is a logical node and cannot be given as data", bugs_node_label(M, i)), s$loc)
    }
  }
  # values that never change: observed stochastic and folded logical nodes
  M$const <- M$nodes$obs | M$nodes$folded
}

bugs_node_label <- function(M, i) {
  nm <- M$nodes$var[i]
  d <- M$dims[[nm]]
  cl <- M$cells[[i]]
  if (length(d) == 1L && d == 1L && !nm %in% names(M$data) && is.null(M$stmts[[M$nodes$stmt[i]]]$lhs$index)) return(nm)
  ai <- arrayInd(cl, d)
  if (length(cl) == 1L) return(sprintf("%s[%s]", nm, paste(ai, collapse = ",")))
  rng <- apply(ai, 2L, function(v) if (min(v) == max(v)) as.character(v[1L]) else sprintf("%d:%d", min(v), max(v)))
  sprintf("%s[%s]", nm, paste(rng, collapse = ","))
}

# ---- parents ----------------------------------------------------------------------------

bugs_stmt_exprs <- function(s) {
  ex <- switch(s$kind, stochastic = s$args, logical = list(s$expr))
  if (!is.null(s$bound)) ex <- c(ex, Filter(Negate(is.null), list(s$bound$lo, s$bound$hi)))
  ex
}

# the cells (as owner node ids) an expression reads in one instance, and
# whether an index depends on nodes
bugs_refs <- function(M, e, env) {
  out <- integer(); poss <- integer()
  walk <- function(e, possible = FALSE) {
    add <- function(ids) if (possible) poss <<- c(poss, ids) else out <<- c(out, ids)
    switch(e$kind,
      num = NULL,
      unop = walk(e$a, possible),
      binop = { walk(e$a, possible); walk(e$b, possible) },
      call = for (a in e$args) walk(a, possible),
      var = {
        if (is.null(e$index) && !is.null(env[[e$name]])) return(NULL)
        own <- M$owner[[e$name]]
        if (is.null(own)) stop_semantic(sprintf("unknown variable '%s'", e$name), e$loc)
        d <- dim(own)
        if (is.null(e$index)) { add(as.vector(own)); return(NULL) }
        if (length(e$index) != length(d) && !(length(d) == 1L && d == 1L && length(e$index) == 1L))
          stop_semantic(sprintf("'%s' has %d dimensions but is used with %d indices", e$name, length(d), length(e$index)), e$loc)
        unresolved <- FALSE
        ix <- lapply(seq_along(e$index), function(k) {
          ii <- e$index[[k]]
          r <- tryCatch(switch(ii$kind,
            single = bugs_const(ii$e, env, M$known),
            range = bugs_const(ii$lo, env, M$known):bugs_const(ii$hi, env, M$known),
            all = seq_len(d[k])), bayescomp_semantic_error = function(err) NULL)
          if (is.null(r)) {
            # an index that depends on nodes: every position it can take
            for (sub in Filter(Negate(is.null), list(ii$e, ii$lo, ii$hi))) walk(sub, possible)
            unresolved <<- TRUE
            r <- seq_len(d[k])
          }
          if (any(r < 1 | r > d[k])) stop_semantic(sprintf("index %s of '%s' is out of range 1:%d", paste(range(r), collapse = ":"), e$name, d[k]), e$loc)
          r
        })
        ids <- as.vector(own[as.matrix(expand.grid(ix))])
        if (unresolved) poss <<- c(poss, ids) else add(ids)
      })
  }
  walk(e)
  list(definite = out, possible = poss)
}

bugs_find_parents <- function(M) {
  n <- length(M$nodes$kind)
  # known values for index evaluation: data, folded and observed values
  M$known <- lapply(stats::setNames(nm = ls(M$S)), function(v) bugs_drop(get(v, envir = M$S)))
  parents <- vector("list", n)
  for (k in seq_along(M$stmts)) {
    s <- M$stmts[[k]]
    ex <- bugs_stmt_exprs(s)
    for (r in seq_len(nrow(s$inst))) {
      id <- s$nodes[r]
      if (M$const[id] && M$nodes$kind[id] == "logical") next
      env <- bugs_inst_env(s, r)
      rf <- lapply(ex, function(e) bugs_refs(M, e, env))
      def <- unique(unlist(lapply(rf, `[[`, "definite"), use.names = FALSE))
      pos <- unique(unlist(lapply(rf, `[[`, "possible"), use.names = FALSE))
      if (id %in% def) stop_semantic(sprintf("the model has a directed cycle: %s depends on itself", bugs_node_label(M, id)), s$loc)
      p <- unique(c(def, pos))
      p <- p[p != 0L & p != id]
      # constants are not parents
      parents[[id]] <- p[!M$const[p]]
    }
  }
  M$parents <- parents
  children <- vector("list", n)
  for (i in seq_len(n)) for (p in parents[[i]]) children[[p]] <- c(children[[p]], i)
  M$children <- children
}

bugs_order <- function(M) {
  n <- length(M$nodes$kind)
  indeg <- lengths(M$parents)
  ord <- integer(n); k <- 0L
  q <- which(indeg == 0L)
  while (length(q)) {
    i <- q[1L]; q <- q[-1L]
    k <- k + 1L; ord[k] <- i
    for (c in M$children[[i]]) { indeg[c] <- indeg[c] - 1L; if (indeg[c] == 0L) q <- c(q, c) }
  }
  if (k < n) {
    bad <- which(indeg > 0L)[1L]
    stop_semantic(sprintf("the model has a directed cycle through %s", bugs_node_label(M, bad)), M$stmts[[M$nodes$stmt[bad]]]$loc)
  }
  M$order <- ord
  M$rank <- integer(n); M$rank[ord] <- seq_len(n)
}

# JAGS data blocks transform the data before the model is compiled
bugs_data_block <- function(stmts, data, dialect) {
  prog <- list(var = list(), data = list(), model = stmts)
  D <- bugs_compile(prog, data, dialect)
  for (i in seq_along(D$nodes$kind)) {
    if (D$nodes$kind[i] == "stochastic" && !D$nodes$obs[i]) stop_unsupported("stochastic relations in a data block are not supported")
    if (D$nodes$kind[i] == "logical" && !D$nodes$folded[i]) stop_semantic(sprintf("%s in the data block cannot be computed from the data", bugs_node_label(D, i)))
  }
  out <- data
  for (nm in names(D$dims)) out[[nm]] <- bugs_drop(get(nm, envir = D$S))
  out
}

# cut() anywhere but at the top of a logical node's expression
bugs_has_cut <- function(e, top = FALSE) {
  switch(e$kind,
    call = (!top && identical(e$name, "cut")) || any(vapply(e$args, bugs_has_cut, TRUE)),
    unop = bugs_has_cut(e$a),
    binop = bugs_has_cut(e$a) || bugs_has_cut(e$b),
    var = any(vapply(e$index %||% list(), function(ix) any(vapply(Filter(Negate(is.null), list(ix$e, ix$lo, ix$hi)), bugs_has_cut, TRUE)), TRUE)),
    FALSE)
}

bugs_drop <- function(a) if (length(dim(a)) == 1L) as.vector(a) else a
