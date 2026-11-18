# MCMC for compiled BUGS models (see bugs_compile.R).
#
# The unobserved stochastic nodes are updated one family at a time. The
# nodes of a statement whose full conditionals do not share terms (no node
# is a descendant of another and their stochastic children differ) form a
# block and are updated together: one evaluation of the vectorized
# statement code gives every node's full conditional, so the samplers below
# work on all of them at once. Nodes that interact are split into blocks by
# greedy coloring.
#
# Samplers:
#   conjugate gamma     a gamma prior on the precision of normal children
#   conjugate Wishart   a Wishart prior on the precision of dmnorm children
#   enumeration         discrete nodes with a few possible values
#   slice               every other scalar node (Neal 2003, stepping out and
#                       shrinkage), on floor(z) for unbounded integers
#   elementwise slice   multivariate normal and Dirichlet nodes (the latter
#                       on log ratios)

# ---- graph helpers --------------------------------------------------------------------

# the logical descendants of a node and its stochastic children (through
# logical nodes)
bugs_desc <- function(M, i) {
  lg <- integer(); st <- integer()
  q <- M$children[[i]]
  seen <- integer()
  while (length(q)) {
    j <- q[1L]; q <- q[-1L]
    if (j %in% seen) next
    seen <- c(seen, j)
    if (M$nodes$kind[j] == "logical") { lg <- c(lg, j); q <- c(q, M$children[[j]]) }
    else st <- c(st, j)
  }
  list(logical = lg, stoch = st)
}

# logical nodes in groups that can be evaluated together: by level within
# the set, then statement
bugs_logical_groups <- function(M, L) {
  if (!length(L)) return(list())
  L <- L[order(M$rank[L])]
  lev <- integer(length(L)); names(lev) <- L
  for (k in seq_along(L)) {
    p <- intersect(M$parents[[L[k]]], L)
    lev[k] <- if (length(p)) 1L + max(lev[as.character(p)]) else 1L
  }
  st <- M$nodes$stmt[L]
  key <- paste(lev, st)
  ord <- order(lev, st, M$nodes$row[L])
  key <- key[ord]; L <- L[ord]; st <- st[ord]
  lapply(split(seq_along(L), factor(key, unique(key))), function(k) list(stmt = st[k[1L]], rows = M$nodes$row[L[k]]))
}

# stochastic children grouped by statement, with the block member each one
# belongs to
bugs_child_groups <- function(M, C, owner) {
  if (!length(C)) return(list())
  st <- M$nodes$stmt[C]
  lapply(split(seq_along(C), st), function(k) list(stmt = st[k[1L]], rows = M$nodes$row[C[k]], owner = owner[k]))
}

# split the nodes of a family into blocks with separate full conditionals
bugs_color <- function(M, ids, desc) {
  touched <- lapply(seq_along(ids), function(k) c(ids[k], desc[[k]]$logical, desc[[k]]$stoch))
  blocks <- list(); used <- list()
  for (k in seq_along(ids)) {
    placed <- FALSE
    for (b in seq_along(blocks)) {
      if (!any(touched[[k]] %in% used[[b]])) {
        blocks[[b]] <- c(blocks[[b]], k); used[[b]] <- c(used[[b]], touched[[k]]); placed <- TRUE; break
      }
    }
    if (!placed) { blocks[[length(blocks) + 1L]] <- k; used[[length(used) + 1L]] <- touched[[k]] }
  }
  blocks
}

# ---- blocks -----------------------------------------------------------------------------

bugs_make_block <- function(M, ids, desc) {
  k <- M$nodes$stmt[ids[1L]]
  s <- M$stmts[[k]]
  d <- bugs_dist_get(s$dist, M$dialect)
  L <- unique(unlist(lapply(desc, `[[`, "logical")))
  C <- unlist(lapply(desc, `[[`, "stoch"))
  owner <- rep(seq_along(ids), vapply(desc, function(x) length(x$stoch), 1L))
  b <- list(ids = ids, stmt = k, rows = M$nodes$row[ids], var = s$lhs$name, dist = d, s = s,
            lgroups = bugs_logical_groups(M, L), cgroups = bugs_child_groups(M, C, owner), children = C, n = length(ids))
  b$sampler <- bugs_choose_sampler(M, b)
  b
}

# full conditional log densities of the block members, after their values
# were written
bugs_cond <- function(M, b) {
  for (g in b$lgroups) M$stmts[[g$stmt]]$fn$write(g$rows)
  lp <- M$stmts[[b$stmt]]$fn$lp(b$rows)
  for (g in b$cgroups) {
    l <- M$stmts[[g$stmt]]$fn$lp(g$rows)
    lp <- lp + if (b$n == 1L) sum(l) else as.vector(rowsum(l, factor(g$owner, seq_len(b$n)), reorder = TRUE))
  }
  lp[is.nan(lp)] <- -Inf
  lp
}

bugs_set <- function(M, var, cells, x) {
  a <- get(var, envir = M$S)
  a[cells] <- x
  assign(var, a, envir = M$S)
}

# ---- sampler choice ---------------------------------------------------------------------

bugs_is_ref <- function(e, var) identical(e$kind, "var") && identical(e$name, var) &&
  (is.null(e$index) || all(vapply(e$index, function(ix) identical(ix$kind, "all"), TRUE)))

bugs_choose_sampler <- function(M, b) {
  d <- b$dist; s <- b$s
  bounded <- !is.null(s$bound)
  # conjugate samplers: every stochastic child is linked directly
  # (logical descendants are allowed only when no child depends on them)
  L <- unlist(lapply(b$lgroups, function(g) M$stmts[[g$stmt]]$nodes[g$rows]))
  direct <- !any(unlist(M$parents[b$children]) %in% L)
  if (!bounded && direct && length(b$children) && b$n == 1L) {
    ch <- lapply(b$cgroups, function(g) M$stmts[[g$stmt]])
    if (d$name == "dgamma" && all(vapply(ch, function(c) identical(c$dist, "dnorm") && is.null(c$bound) &&
                                         bugs_is_ref(c$args[[2L]], b$var) && !bugs_mentions(c$args[[1L]], b$var), TRUE)))
      return("gamma_normal")
    if (d$name == "dwish" && all(vapply(ch, function(c) identical(c$dist, "dmnorm") && is.null(c$bound) &&
                                        bugs_is_ref(c$args[[2L]], b$var) && !bugs_mentions(c$args[[1L]], b$var), TRUE)))
      return("wishart_mnorm")
  }
  if (d$name == "dwish") stop_unsupported(sprintf("a Wishart node (%s) whose children are not multivariate normals with it as precision", b$var), s$loc)
  if (d$support == "multi") stop_unsupported(sprintf("an unobserved multinomial node (%s)", b$var), s$loc)
  if (d$support == "simplex") return("dirichlet")
  if (d$support %in% c("mvreal")) return("mv_slice")
  if (isTRUE(d$discrete)) {
    if (d$support == "finite") return("enumerate")
    return("discrete_slice")
  }
  "slice"
}

bugs_mentions <- function(e, var) {
  switch(e$kind,
    num = FALSE,
    var = identical(e$name, var) || any(vapply(e$index %||% list(), function(ix)
      any(vapply(Filter(Negate(is.null), list(ix$e, ix$lo, ix$hi)), bugs_mentions, TRUE, var = var)), TRUE)),
    unop = bugs_mentions(e$a, var),
    binop = bugs_mentions(e$a, var) || bugs_mentions(e$b, var),
    call = any(vapply(e$args, bugs_mentions, TRUE, var = var)),
    FALSE)
}

# ---- samplers ---------------------------------------------------------------------------

# vectorized slice sampling of the scalar values at `cells` (one per block
# member); `w` are the widths. With `discrete`, the slice variable is
# continuous and the value its floor.
bugs_slice <- function(M, b, cells, w, discrete = FALSE, m = 50L) {
  n <- length(cells)
  a <- get(b$var, envir = M$S)
  x0 <- a[cells]
  z0 <- if (discrete) x0 + stats::runif(n) else x0
  val <- function(z) if (discrete) floor(z) else z
  f <- function(z) { bugs_set(M, b$var, cells, val(z)); bugs_cond(M, b) }
  y <- f(z0) - stats::rexp(n)
  if (any(!is.finite(y))) {
    bad <- which(!is.finite(y))[1L]
    stop_sampler(sprintf("the full conditional of %s is not finite at its current value", bugs_node_label(M, b$ids[bad])))
  }
  L <- z0 - w * stats::runif(n); R <- L + w
  act <- rep(TRUE, n); steps <- 0L
  while (any(act) && steps < m) {
    z <- z0; z[act] <- L[act]
    fl <- f(z)
    act <- act & fl > y
    L[act] <- L[act] - w[act]
    steps <- steps + 1L
  }
  act <- rep(TRUE, n); steps <- 0L
  while (any(act) && steps < m) {
    z <- z0; z[act] <- R[act]
    fr <- f(z)
    act <- act & fr > y
    R[act] <- R[act] + w[act]
    steps <- steps + 1L
  }
  z1 <- z0
  todo <- rep(TRUE, n); it <- 0L
  while (any(todo)) {
    it <- it + 1L
    if (it > 200L) stop_sampler(sprintf("the slice sampler for %s did not find a new value", b$var))
    z <- z1
    z[todo] <- L[todo] + stats::runif(sum(todo)) * (R[todo] - L[todo])
    fz <- f(z)
    ok <- todo & is.finite(fz) & fz >= y
    z1[ok] <- z[ok]
    todo <- todo & !ok
    shrink <- todo & z < z0
    L[shrink] <- z[shrink]
    grow <- todo & z >= z0
    R[grow] <- z[grow]
  }
  bugs_set(M, b$var, cells, val(z1))
  for (g in b$lgroups) M$stmts[[g$stmt]]$fn$write(g$rows)
  abs(z1 - z0)
}

stop_sampler <- function(message) stop(bc_condition("bayescomp_sampler_error", paste("Sampler error:", message)))

bugs_enumerate <- function(M, b, cells) {
  fn <- M$stmts[[b$stmt]]$fn
  a <- fn$args(b$rows)
  bd <- if (fn$mode == "vec") { r <- do.call(b$dist$bounds, a); list(rep_len(r[[1L]], b$n), rep_len(r[[2L]], b$n)) }
        else { r <- lapply(a, function(x) do.call(b$dist$bounds, x)); list(vapply(r, `[[`, 0, 1L), vapply(r, `[[`, 0, 2L)) }
  if (!is.null(b$s$bound)) {
    bb <- fn$bounds(b$rows)
    if (fn$mode == "vec") { bd[[1L]] <- pmax(bd[[1L]], ceiling(bb[[1L]])); bd[[2L]] <- pmin(bd[[2L]], floor(bb[[2L]])) }
  }
  lo <- min(bd[[1L]]); hi <- max(bd[[2L]])
  if (!is.finite(hi - lo) || hi - lo > 100) {
    return(bugs_slice(M, b, cells, rep(2, b$n), discrete = TRUE))
  }
  vals <- lo:hi
  lp <- matrix(-Inf, b$n, length(vals))
  for (k in seq_along(vals)) {
    inr <- vals[k] >= bd[[1L]] & vals[k] <= bd[[2L]]
    if (!any(inr)) next
    bugs_set(M, b$var, cells, ifelse(inr, vals[k], bd[[1L]]))
    l <- bugs_cond(M, b)
    lp[inr, k] <- l[inr]
  }
  mx <- apply(lp, 1L, max)
  if (any(!is.finite(mx))) stop_sampler(sprintf("every value of %s has zero probability", bugs_node_label(M, b$ids[which(!is.finite(mx))[1L]])))
  p <- exp(lp - mx)
  cp <- p / rowSums(p)
  u <- stats::runif(b$n)
  pick <- 1L + rowSums(u > t(apply(cp, 1L, cumsum)))
  bugs_set(M, b$var, cells, vals[pick])
  for (g in b$lgroups) M$stmts[[g$stmt]]$fn$write(g$rows)
  NULL
}

bugs_gamma_normal <- function(M, b) {
  fn <- M$stmts[[b$stmt]]$fn
  a <- fn$args(b$rows)
  if (fn$mode != "vec") a <- a[[1L]]
  shape <- a[[1L]]; rate <- a[[2L]]
  for (g in b$cgroups) {
    cf <- M$stmts[[g$stmt]]$fn
    ca <- cf$args(g$rows)
    x <- cf$value(g$rows)
    if (cf$mode != "vec") { mu <- unlist(lapply(ca, `[[`, 1L)); x <- unlist(x) } else mu <- ca[[1L]]
    r <- x - mu
    shape <- shape + length(r) / 2
    rate <- rate + sum(r * r) / 2
  }
  bugs_set(M, b$var, M$cells[[b$ids]], stats::rgamma(1L, shape, rate))
}

bugs_wishart_mnorm <- function(M, b) {
  fn <- M$stmts[[b$stmt]]$fn
  a <- fn$args(b$rows)[[1L]]
  R <- a[[1L]]; k <- a[[2L]]
  for (g in b$cgroups) {
    cf <- M$stmts[[g$stmt]]$fn
    ca <- cf$args(g$rows)
    x <- cf$value(g$rows)
    for (j in seq_along(g$rows)) {
      r <- x[[j]] - ca[[j]][[1L]]
      R <- R + tcrossprod(r)
      k <- k + 1
    }
  }
  W <- stats::rWishart(1L, k, solve(R))[, , 1L]
  bugs_set(M, b$var, M$cells[[b$ids]], as.vector(W))
}

# elementwise slice sampling of multivariate normal nodes
bugs_mv_slice <- function(M, b, w) {
  K <- length(M$cells[[b$ids[1L]]])
  for (j in seq_len(K)) {
    cells <- vapply(b$ids, function(i) M$cells[[i]][j], 1L)
    jump <- bugs_slice(M, b, cells, w[, j])
    w[, j] <- bugs_adapt_w(w[, j], jump, b$adapting)
  }
  w
}

# Dirichlet nodes: slice sampling of z_j = log(x_j / x_K) one coordinate at
# a time, with the Jacobian prod(x)
bugs_dirichlet <- function(M, b) {
  K <- length(M$cells[[b$ids[1L]]])
  cells <- lapply(b$ids, function(i) M$cells[[i]])
  get_x <- function() { a <- get(b$var, envir = M$S); t(vapply(cells, function(c) a[c], numeric(K))) }
  for (j in seq_len(K - 1L)) {
    X <- matrix(get_x(), b$n)
    z0 <- log(X[, j] / X[, K])
    set_z <- function(z) {
      Z <- log(X / X[, K]); Z[, j] <- z
      E <- exp(Z - apply(Z, 1L, max))
      Y <- E / rowSums(E)
      a <- get(b$var, envir = M$S)
      for (i in seq_len(b$n)) a[cells[[i]]] <- Y[i, ]
      assign(b$var, a, envir = M$S)
      rowSums(log(Y))
    }
    f <- function(z) { jac <- set_z(z); bugs_cond(M, b) + jac }
    y <- f(z0) - stats::rexp(b$n)
    w <- rep(1, b$n)
    L <- z0 - w * stats::runif(b$n); R <- L + w
    for (s in 1:20) { z <- z0; act <- f(replace(z, TRUE, L)) > y; if (!any(act)) break; L[act] <- L[act] - w[act] }
    for (s in 1:20) { act <- f(replace(z0, TRUE, R)) > y; if (!any(act)) break; R[act] <- R[act] + w[act] }
    z1 <- z0; todo <- rep(TRUE, b$n)
    for (it in 1:200) {
      z <- z1; z[todo] <- L[todo] + stats::runif(sum(todo)) * (R[todo] - L[todo])
      fz <- f(z)
      ok <- todo & is.finite(fz) & fz >= y
      z1[ok] <- z[ok]; todo <- todo & !ok
      if (!any(todo)) break
      L[todo & z < z0] <- z[todo & z < z0]; R[todo & z >= z0] <- z[todo & z >= z0]
    }
    set_z(z1)
  }
  for (g in b$lgroups) M$stmts[[g$stmt]]$fn$write(g$rows)
}

# slice widths adapt to twice the mean jump during adaptation
bugs_adapt_w <- function(w, jump, adapting) {
  if (!adapting) return(w)
  pmax(0.95 * w + 0.05 * 2 * pmax(jump, 1e-10), 1e-8)
}

bugs_update <- function(M, b) {
  cells <- if (b$sampler %in% c("slice", "discrete_slice", "enumerate")) vapply(b$ids, function(i) M$cells[[i]], 1L)
  switch(b$sampler,
    slice = { jump <- bugs_slice(M, b, cells, b$w); b$w <- bugs_adapt_w(b$w, jump, b$adapting) },
    discrete_slice = { jump <- bugs_slice(M, b, cells, b$w, discrete = TRUE); b$w <- pmax(1, bugs_adapt_w(b$w, jump, b$adapting)) },
    enumerate = bugs_enumerate(M, b, cells),
    gamma_normal = bugs_gamma_normal(M, b),
    wishart_mnorm = bugs_wishart_mnorm(M, b),
    mv_slice = { b$w <- bugs_mv_slice(M, b, b$w) },
    dirichlet = bugs_dirichlet(M, b))
  if (b$sampler %in% c("gamma_normal", "wishart_mnorm")) for (g in b$lgroups) M$stmts[[g$stmt]]$fn$write(g$rows)
  b
}

# ---- setting up ---------------------------------------------------------------------------

bugs_blocks <- function(M) {
  U <- which(M$nodes$kind == "stochastic" & !M$nodes$obs)
  U <- U[order(M$rank[U])]
  desc <- lapply(U, function(i) bugs_desc(M, i))
  names(desc) <- U
  blocks <- list()
  for (k in unique(M$nodes$stmt[U])) {
    ids <- U[M$nodes$stmt[U] == k]
    dd <- desc[as.character(ids)]
    for (cl in bugs_color(M, ids, dd)) {
      b <- bugs_make_block(M, ids[cl], dd[cl])
      K <- length(M$cells[[b$ids[1L]]])
      b$w <- if (b$sampler == "mv_slice") matrix(1, b$n, K) else rep(if (b$sampler == "discrete_slice") 2 else 1, b$n)
      b$adapting <- TRUE
      blocks[[length(blocks) + 1L]] <- b
    }
  }
  # blocks in the order of their first node, as in a sweep through the graph
  blocks[order(vapply(blocks, function(b) min(M$rank[b$ids]), 1L))]
}

# initial values: given ones, then the typical value of the prior (JAGS) or
# a draw from it (OpenBUGS gen inits), in graph order
bugs_initialize <- function(M, inits, how = c("typical", "prior")) {
  how <- match.arg(how)
  inits <- bugs_check_values(inits %||% list(), "inits")
  U <- which(M$nodes$kind == "stochastic" & !M$nodes$obs)
  uvars <- unique(M$nodes$var[U])
  for (nm in names(inits)) {
    if (!nm %in% names(M$dims)) { warning(sprintf("the initial value of '%s' is ignored: not a variable of the model", nm), call. = FALSE); next }
    v <- inits[[nm]]
    a <- get(nm, envir = M$S)
    if (length(v) != length(a)) stop_data(sprintf("the initial value of '%s' has %d values but the variable has %d", nm, length(v), length(a)))
    ok <- !is.na(v)
    # only unobserved stochastic elements take initial values
    own <- M$owner[[nm]]
    settable <- own != 0L & M$nodes$kind[pmax(own, 1L)] == "stochastic" & !M$nodes$obs[pmax(own, 1L)]
    bad <- ok & !settable & is.na(a)
    if (any(bad)) warning(sprintf("initial values of '%s' for logical nodes are ignored", nm), call. = FALSE)
    if (any(ok & !is.na(a))) warning(sprintf("initial values of '%s' for data are ignored", nm), call. = FALSE)
    a[ok & settable] <- as.numeric(v)[ok & settable]
    assign(nm, a, envir = M$S)
  }
  given <- lapply(stats::setNames(nm = unique(M$nodes$var)), function(v) get(v, envir = M$S))
  # draws from vague priors can have zero probability (a precision drawn
  # from dgamma(0.001, 0.001) is often 0): try again, then use typical values
  for (attempt in seq_len(if (how == "prior") 21L else 1L)) {
    use <- if (attempt == 21L) "typical" else how
    for (v in names(given)) assign(v, given[[v]], envir = M$S)
    bugs_init_fill(M, use)
    if (bugs_init_ok(M)) return(invisible())
  }
  for (i in M$order) if (M$nodes$kind[i] == "stochastic") {
    lp <- M$stmts[[M$nodes$stmt[i]]]$fn$lp(M$nodes$row[i])
    if (!is.finite(sum(lp))) stop_sampler(sprintf("%s has zero probability at the initial values%s", bugs_node_label(M, i),
                                                 if (M$nodes$obs[i]) "; check the data and initial values" else ""))
  }
  invisible()
}

bugs_init_ok <- function(M) {
  for (k in seq_along(M$stmts)) if (M$stmts[[k]]$kind == "stochastic" && nrow(M$stmts[[k]]$inst))
    if (!is.finite(sum(M$stmts[[k]]$fn$lp(seq_len(nrow(M$stmts[[k]]$inst)))))) return(FALSE)
  TRUE
}

bugs_init_fill <- function(M, how) {
  for (i in M$order) {
    s <- M$stmts[[M$nodes$stmt[i]]]
    r <- M$nodes$row[i]
    if (M$nodes$kind[i] == "logical") { s$fn$write(r); next }
    if (M$nodes$obs[i]) next
    cur <- get(M$nodes$var[i], envir = M$S)[M$cells[[i]]]
    if (all(!is.na(cur))) next
    v <- if (how == "prior") s$fn$rng(r) else s$fn$typical(r)
    v <- unlist(v)
    if (!is.null(s$bound) && how == "typical") {
      bb <- s$fn$bounds(r); bb <- if (s$fn$mode == "vec") bb else bb[[1L]]
      lo <- bb[[1L]]; hi <- bb[[2L]]
      if (any(v < lo | v > hi)) v <- if (is.finite(lo) && is.finite(hi)) (lo + hi) / 2 else if (is.finite(lo)) lo + 1 else hi - 1
      if (isTRUE(s$fn$discrete)) v <- round(v)
    }
    bugs_set(M, M$nodes$var[i], M$cells[[i]], v)
  }
}

# the deviance: -2 times the log density of the observed stochastic nodes
bugs_deviance_fun <- function(M) {
  O <- which(M$nodes$kind == "stochastic" & M$nodes$obs)
  g <- lapply(split(M$nodes$row[O], M$nodes$stmt[O]), identity)
  ks <- as.integer(names(g))
  function() {
    s <- 0
    for (j in seq_along(g)) s <- s + sum(M$stmts[[ks[j]]]$fn$lp(g[[j]]))
    -2 * s
  }
}

# names and cells of the monitored variables, elements in column major order
bugs_monitor <- function(M, vars) {
  out <- list()
  for (v in vars) {
    if (identical(v, "deviance")) next
    d <- M$dims[[v]]
    if (is.null(d)) stop_semantic(sprintf("cannot monitor '%s': not a variable of the model", v))
    n <- prod(d)
    nm <- if (n == 1L && is.null(M$data[[v]]) && length(d) == 1L) v
          else paste0(v, "[", apply(arrayInd(seq_len(n), d), 1L, paste, collapse = ","), "]")
    out[[v]] <- list(var = v, names = nm)
  }
  out
}

# one chain: n_adapt adaptation iterations, then n_iter iterations of which
# those after n_burnin are kept every n_thin
bugs_chain <- function(M, inits, n_adapt, n_burnin, n_iter, n_thin, monitor, deviance, init_how, refresh = 0L, chain = 1L) {
  bugs_initialize(M, inits, init_how)
  blocks <- bugs_blocks(M)
  mon <- bugs_monitor(M, monitor)
  devf <- if (deviance) bugs_deviance_fun(M)
  vnames <- c(unlist(lapply(mon, `[[`, "names"), use.names = FALSE), if (deviance) "deviance")
  keep <- n_burnin + n_thin * seq_len(max(0L, (n_iter - n_burnin) %/% n_thin))
  out <- matrix(NA_real_, length(keep), length(vnames), dimnames = list(NULL, vnames))
  k <- 0L
  total <- n_adapt + n_iter
  # running means of the unobserved stochastic variables, for OpenBUGS's Dhat
  uvars <- unique(M$nodes$var[M$nodes$kind == "stochastic" & !M$nodes$obs])
  sums <- lapply(stats::setNames(nm = uvars), function(v) 0)
  for (t in seq_len(total)) {
    adapting <- t <= n_adapt
    for (j in seq_along(blocks)) {
      blocks[[j]]$adapting <- adapting
      blocks[[j]] <- bugs_update(M, blocks[[j]])
    }
    it <- t - n_adapt
    if (it >= 1L && it %in% keep) {
      k <- k + 1L
      row <- unlist(lapply(mon, function(m) as.vector(get(m$var, envir = M$S))), use.names = FALSE)
      if (deviance) row <- c(row, devf())
      out[k, ] <- row
      for (v in uvars) sums[[v]] <- sums[[v]] + get(v, envir = M$S)
    }
    if (refresh > 0L && t %% refresh == 0L) message(sprintf("chain %d: iteration %d / %d", chain, t, total))
  }
  list(draws = out, means = lapply(sums, function(x) x / max(k, 1L)), samplers = table(vapply(blocks, `[[`, "", "sampler")))
}
