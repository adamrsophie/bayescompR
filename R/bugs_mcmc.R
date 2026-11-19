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
# Samplers, chosen in this order:
#   conjugate        gamma (normal and lognormal precision, Poisson mean,
#                    exponential and gamma rate), beta (Bernoulli and
#                    binomial probability), Dirichlet (categorical and
#                    multinomial probabilities), Wishart (multivariate
#                    normal precision) and normal (a normal mean or a
#                    coefficient entering normal means linearly) priors
#   enumeration      discrete nodes with at most 101 possible values
#   slice            every other scalar node (Neal 2003, stepping out with
#                    a randomly split step budget, then shrinkage); for
#                    unbounded integers the slice variable is continuous
#                    and the value its floor
#   elementwise      multivariate normal and t nodes (only their missing
#   slice            elements when partly observed) and Dirichlet nodes
#                    (on log ratios of their nonzero elements)
#   pairwise         multinomial nodes: the counts of two categories are
#   enumeration      redrawn with their sum fixed
#
# cut(): a node y <- cut(x) passes values of x down but no likelihood up,
# so updates of x ignore the stochastic descendants of y (OpenBUGS).

# ---- graph helpers --------------------------------------------------------------------

# the logical descendants of a node and its stochastic children (through
# logical nodes); below a cut() node logical descendants are recomputed but
# their stochastic children are not part of the full conditional
bugs_desc <- function(M, i) {
  lg <- integer(); st <- integer()
  q <- M$children[[i]]; qcut <- rep(FALSE, length(q))
  seen <- integer()
  while (length(q)) {
    j <- q[1L]; cut <- qcut[1L]; q <- q[-1L]; qcut <- qcut[-1L]
    if (j %in% seen) next
    seen <- c(seen, j)
    if (M$nodes$kind[j] == "logical") {
      lg <- c(lg, j)
      ch <- M$children[[j]]
      q <- c(q, ch); qcut <- c(qcut, rep(cut || isTRUE(M$stmts[[M$nodes$stmt[j]]]$cut), length(ch)))
    } else if (!cut) st <- c(st, j)
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
    p <- M$parents[[L[k]]]
    p <- p[as.character(p) %in% names(lev)[seq_len(k - 1L)]]
    lev[k] <- if (length(p)) 1L + max(lev[as.character(p)]) else 1L
  }
  st <- M$nodes$stmt[L]
  key <- paste(lev, st)
  ord <- order(lev, st, M$nodes$row[L])
  key <- key[ord]; L <- L[ord]; st <- st[ord]
  lapply(split(seq_along(L), factor(key, unique(key))), function(k) list(stmt = st[k[1L]], rows = M$nodes$row[L[k]], ids = L[k]))
}

# stochastic children grouped by statement, with the block member each one
# belongs to
bugs_child_groups <- function(M, C, owner, n) {
  if (!length(C)) return(list())
  st <- M$nodes$stmt[C]
  lapply(split(seq_along(C), st), function(k) {
    o <- owner[k]
    list(stmt = st[k[1L]], rows = M$nodes$row[C[k]], ids = C[k], owner = o, uowner = sort(unique(o)),
         identity = length(o) == n && all(o == seq_len(n)))
  })
}

# split the nodes of a family into blocks with separate full conditionals
bugs_color <- function(M, ids, desc) {
  N <- length(M$nodes$kind)
  blocks <- list(); used <- list()
  for (k in seq_along(ids)) {
    tk <- c(ids[k], desc[[k]]$logical, desc[[k]]$stoch)
    placed <- FALSE
    for (b in seq_along(blocks)) {
      if (!any(used[[b]][tk])) {
        blocks[[b]] <- c(blocks[[b]], k); used[[b]][tk] <- TRUE; placed <- TRUE; break
      }
    }
    if (!placed) {
      u <- logical(N); u[tk] <- TRUE
      blocks[[length(blocks) + 1L]] <- k; used[[length(used) + 1L]] <- u
    }
  }
  blocks
}

# ---- blocks -----------------------------------------------------------------------------

bugs_make_block <- function(M, ids, desc) {
  k <- M$nodes$stmt[ids[1L]]
  s <- M$stmts[[k]]
  L <- unique(unlist(lapply(desc, `[[`, "logical")))
  # logical nodes that no stochastic node depends on are computed once per
  # recorded iteration instead (bugs_feeds)
  L <- L[M$feeds[L]]
  C <- unlist(lapply(desc, `[[`, "stoch"))
  owner <- rep(seq_along(ids), vapply(desc, function(x) length(x$stoch), 1L))
  K <- length(M$cells[[ids[1L]]])
  b <- list(ids = ids, stmt = k, rows = M$nodes$row[ids], var = s$lhs$name, dist = bugs_dist_get(s$dist, M$dialect), s = s,
            lgroups = bugs_logical_groups(M, L), lids = L, cgroups = bugs_child_groups(M, C, owner, length(ids)),
            children = C, child_owner = owner, n = length(ids), K = K,
            cells = matrix(unlist(M$cells[ids]), length(ids), K, byrow = TRUE))
  b$fixed <- matrix(FALSE, b$n, K)
  for (j in seq_along(ids)) { m <- M$fixed[[as.character(ids[j])]]; if (!is.null(m)) b$fixed[j, ] <- m }
  b$sampler <- bugs_choose_sampler(M, b)
  b$dyn <- startsWith(b$sampler, "conj_dyn_")
  if (b$dyn) {
    b$sampler <- sub("conj_dyn_", "conj_", b$sampler, fixed = TRUE)
    b$dyn_owner <- lapply(b$cgroups, bugs_dyn_owner_fun, M = M, b = b)
  }
  b$w <- if (b$sampler == "mv_slice") matrix(1, b$n, K) else rep(if (b$sampler == "discrete_slice") 2 else 1, b$n)
  b$adapting <- TRUE
  b
}

bugs_write_desc <- function(M, b) for (g in b$lgroups) M$stmts[[g$stmt]]$fn$write(g$rows)

# full conditional log densities of the block members at the values now
# stored. `cells` and `safe` (the coordinate being changed and its current
# values) let a value outside the support be rejected before the logical
# descendants see it (an index out of range would fail there).
bugs_cond <- function(M, b, cells = NULL, safe = NULL) {
  lp <- M$stmts[[b$stmt]]$fn$lp(b$rows)
  lp[is.nan(lp)] <- -Inf
  bad <- !is.finite(lp)
  if (any(bad) && !is.null(cells)) bugs_set(M, b$var, cells[bad], safe[bad])
  bugs_write_desc(M, b)
  for (g in b$cgroups) {
    l <- M$stmts[[g$stmt]]$fn$lp(g$rows)
    if (b$n == 1L) lp <- lp + sum(l)
    else if (g$identity) lp <- lp + l
    else lp[g$uowner] <- lp[g$uowner] + as.vector(rowsum(l, g$owner, reorder = TRUE))
  }
  lp[is.nan(lp) | bad] <- -Inf
  lp
}

bugs_set <- function(M, var, cells, x) {
  a <- get(var, envir = M$S)
  a[cells] <- x
  assign(var, a, envir = M$S)
}

# ---- sampler choice ---------------------------------------------------------------------

# is the parameter `slot` of the child exactly the block member (all of it)?
bugs_child_ref <- function(M, child_id, slot, target_id) {
  s <- M$stmts[[M$nodes$stmt[child_id]]]
  e <- s$args[[slot]]
  if (!identical(e$kind, "var")) return(FALSE)
  env <- bugs_inst_env(s, M$nodes$row[child_id])
  r <- tryCatch(bugs_refs(M, e, env), error = function(err) NULL)
  !is.null(r) && !length(r$possible) && all(r$definite == target_id) &&
    length(r$definite) == length(M$cells[[target_id]])
}

# conjugate families: prior -> child distribution -> the slot holding the node
.bugs_conj <- list(
  dgamma = list(dnorm = 2L, dlnorm = 2L, dpois = 1L, dexp = 1L, dgamma = 2L),
  dbeta = list(dbern = 1L, dbin = 1L, dbinom = 1L),
  ddirch = list(dcat = 1L, dmulti = 1L),
  dwish = list(dmnorm = 2L))

bugs_choose_sampler <- function(M, b) {
  d <- b$dist; s <- b$s
  bounded <- !is.null(s$bound)
  conj <- .bugs_conj[[d$name]]
  dyn <- FALSE
  if (!bounded && !any(b$fixed) && !is.null(conj) && length(b$children)) {
    # children may not depend on the node through logical nodes
    direct <- !any(unlist(M$parents[b$children]) %in% b$lids)
    ok <- direct && all(vapply(seq_along(b$children), function(k) {
      c <- b$children[k]
      cs <- M$stmts[[M$nodes$stmt[c]]]
      slot <- conj[[cs$dist]]
      if (is.null(slot) || !is.null(cs$bound)) return(FALSE)
      tgt <- b$ids[b$child_owner[k]]
      env <- bugs_inst_env(cs, M$nodes$row[c])
      if (!bugs_child_ref(M, c, slot, tgt)) {
        # a mixture: the parameter is a node of the family chosen by a
        # latent index, as in p[z[i], ]; which node it is is found at each
        # update
        e <- cs$args[[slot]]
        if (!identical(e$kind, "var") || e$name != b$var || !length(e$index)) return(FALSE)
        r <- bugs_refs(M, e, env)
        if (!length(r$possible)) return(FALSE)
        if (any(M$nodes$var[c(r$definite, r$possible)] == b$var & M$nodes$kind[c(r$definite, r$possible)] != "stochastic")) return(FALSE)
        dyn <<- TRUE
        tgt <- b$ids
      }
      # the other parameters of the child do not involve the node
      !any(vapply(setdiff(seq_along(cs$args), slot), function(o) {
        r <- bugs_refs(M, cs$args[[o]], env); any(c(r$definite, r$possible) %in% tgt)
      }, TRUE))
    }, TRUE))
    if (ok) return(paste0(if (dyn) "conj_dyn_" else "conj_", d$name))
  }
  if (d$name == "dwish") {
    if (!length(b$children) && !bounded) return("conj_dwish")
    stop_unsupported(sprintf("a Wishart node (%s) whose children are not multivariate normals with it as precision", b$var), s$loc)
  }
  if (d$name == "dnorm" && !bounded && length(b$children) && bugs_linear_normal(M, b)) return("conj_dnorm")
  if (d$support == "multi") return("multinomial")
  if (d$support == "simplex") return("dirichlet")
  if (d$support == "mvreal") return("mv_slice")
  if (d$dim %in% c("vector", "matrix")) stop_unsupported(sprintf("sampling %s nodes (%s)", d$name, b$var), s$loc)
  if (isTRUE(d$discrete)) return(if (d$support == "finite") "enumerate" else "discrete_slice")
  "slice"
}

# Is every stochastic child a normal whose mean is linear in the node and
# whose precision does not involve it? The mean is followed through logical
# nodes: sums, differences, products with terms free of the node, and
# division by such terms keep it linear.
bugs_linear_normal <- function(M, b) {
  memo <- new.env(parent = emptyenv())
  for (k in seq_along(b$children)) {
    c <- b$children[k]
    cs <- M$stmts[[M$nodes$stmt[c]]]
    if (!identical(cs$dist, "dnorm") || !is.null(cs$bound)) return(FALSE)
    tgt <- b$ids[b$child_owner[k]]
    env <- bugs_inst_env(cs, M$nodes$row[c])
    if (bugs_lin(M, cs$args[[1L]], env, tgt, b$lids, memo) != "l") return(FALSE)
    if (bugs_lin(M, cs$args[[2L]], env, tgt, b$lids, memo) != "c") return(FALSE)
  }
  TRUE
}

# "c" free of the target, "l" linear in it, "n" otherwise
bugs_lin <- function(M, e, env, tgt, L, memo) {
  lin <- function(e) bugs_lin(M, e, env, tgt, L, memo)
  switch(e$kind,
    num = "c",
    var = {
      r <- bugs_refs(M, e, env)
      ids <- unique(c(r$definite, r$possible))
      if (!any(ids %in% c(tgt, L))) return("c")
      if (any(r$possible %in% c(tgt, L)) || length(ids) != 1L) return("n")
      if (ids == tgt) return("l")
      key <- paste(ids, tgt)
      if (!is.null(memo[[key]])) return(memo[[key]])
      ls <- M$stmts[[M$nodes$stmt[ids]]]
      res <- if (!is.null(ls$link) || isTRUE(ls$cut) || length(M$cells[[ids]]) != 1L) "n"
             else bugs_lin(M, ls$expr, bugs_inst_env(ls, M$nodes$row[ids]), tgt, L, memo)
      assign(key, res, envir = memo)
      res
    },
    unop = { a <- lin(e$a); if (e$op == "-") a else if (a == "c") "c" else "n" },
    binop = {
      a <- lin(e$a); bb <- lin(e$b)
      if (a == "n" || bb == "n") return("n")
      switch(e$op,
        "+" = , "-" = if (a == "l" || bb == "l") "l" else "c",
        "*" = if (a == "l" && bb == "l") "n" else if (a == "l" || bb == "l") "l" else "c",
        "/" = if (bb == "l") "n" else a,
        if (a == "c" && bb == "c") "c" else "n")
    },
    call = if (all(vapply(e$args, lin, "") == "c")) "c" else "n",
    "n")
}

# ---- the slice sampler ------------------------------------------------------------------

# One slice update of n independent coordinates at once. f maps a vector of
# n values to their n log densities; w are the widths and m the step budget,
# split at random between the left and right stepping out (Neal 2003,
# figure 3), which keeps the target invariant.
bugs_slice_core <- function(f, z0, w, m = 1000L, what = "a node") {
  n <- length(z0)
  y <- f(z0) - stats::rexp(n)
  if (any(!is.finite(y))) stop_sampler(sprintf("the full conditional of %s is not finite at its current value", what))
  L <- z0 - w * stats::runif(n); R <- L + w
  J <- floor(m * stats::runif(n)); K <- (m - 1L) - J
  act <- J > 0
  while (any(act)) {
    z <- z0; z[act] <- L[act]
    act <- act & f(z) > y
    L[act] <- L[act] - w[act]; J[act] <- J[act] - 1L
    act <- act & J > 0
  }
  act <- K > 0
  while (any(act)) {
    z <- z0; z[act] <- R[act]
    act <- act & f(z) > y
    R[act] <- R[act] + w[act]; K[act] <- K[act] - 1L
    act <- act & K > 0
  }
  z1 <- z0
  todo <- rep(TRUE, n); it <- 0L
  while (any(todo)) {
    it <- it + 1L
    if (it > 500L) stop_sampler(sprintf("the slice sampler for %s did not find a new value", what))
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
  z1
}

stop_sampler <- function(message) stop(bc_condition("bayescomp_sampler_error", paste("Sampler error:", message)))

# slice sampling of one scalar coordinate (cells, one per member) of the
# members `sel` of a block
bugs_slice_coord <- function(M, b, cells, sel, w, discrete = FALSE) {
  x0 <- get(b$var, envir = M$S)[cells]
  val <- if (discrete) floor else identity
  f <- function(z) {
    v <- x0; v[sel] <- val(z)
    bugs_set(M, b$var, cells, v)
    bugs_cond(M, b, cells, x0)[sel]
  }
  z0 <- if (discrete) x0[sel] + stats::runif(length(sel)) else x0[sel]
  z1 <- bugs_slice_core(f, z0, w, what = b$var)
  v <- x0; v[sel] <- val(z1)
  bugs_set(M, b$var, cells, v)
  bugs_write_desc(M, b)
  abs(z1 - z0)
}

# Slice widths adapt during adaptation. For a uniform draw on a slice of
# length l the mean jump is l / 3, and stepping out costs least at a width of
# about 1.36 l (Campbell 2026; Yallup et al. 2026), so the width tracks
# 3.6 times the mean jump.
bugs_adapt_w <- function(w, jump, adapting) {
  if (!adapting) return(w)
  pmax(0.95 * w + 0.05 * 3.6 * pmax(jump, 1e-10), 1e-8)
}

# ---- other samplers ---------------------------------------------------------------------

# the parameters of each row, as a list of lists
bugs_args_list <- function(fn, rows) {
  a <- fn$args(rows)
  if (fn$mode == "vec") lapply(seq_along(rows), function(j) lapply(a, function(x) x[(j - 1L) %% length(x) + 1L])) else a
}

bugs_enumerate <- function(M, b, cells) {
  fn <- M$stmts[[b$stmt]]$fn
  a <- bugs_args_list(fn, b$rows)
  bd <- vapply(a, function(x) unlist(do.call(b$dist$bounds, x)), numeric(2L))
  lo <- bd[1L, ]; hi <- bd[2L, ]
  if (!is.null(b$s$bound)) {
    bb <- fn$bounds(b$rows)
    if (fn$mode == "vec") { blo <- rep_len(bb[[1L]], b$n); bhi <- rep_len(bb[[2L]], b$n) }
    else { blo <- vapply(bb, function(x) x[[1L]], 0); bhi <- vapply(bb, function(x) x[[2L]], 0) }
    lo <- pmax(lo, ceiling(blo)); hi <- pmin(hi, floor(bhi))
  }
  if (!all(is.finite(c(lo, hi))) || max(hi - lo) > 100) {
    bugs_slice_coord(M, b, cells, seq_len(b$n), rep(2, b$n), discrete = TRUE)
    return(invisible())
  }
  x0 <- get(b$var, envir = M$S)[cells]
  vals <- min(lo):max(hi)
  lp <- matrix(-Inf, b$n, length(vals))
  for (k in seq_along(vals)) {
    inr <- vals[k] >= lo & vals[k] <= hi
    if (!any(inr)) next
    bugs_set(M, b$var, cells, ifelse(inr, vals[k], x0))
    l <- bugs_cond(M, b, cells, x0)
    lp[inr, k] <- l[inr]
  }
  mx <- apply(lp, 1L, max)
  if (any(!is.finite(mx))) stop_sampler(sprintf("every value of %s has zero probability", bugs_node_label(M, b$ids[which(!is.finite(mx))[1L]])))
  p <- exp(lp - mx)
  cp <- p / rowSums(p)
  # Metropolized Gibbs (Liu 1996): propose a value other than the current one
  # with probability proportional to its conditional probability and accept
  # with min(1, (1 - p_current) / (1 - p_proposed)). It leaves the conditional
  # invariant and never wastes an update on staying put when another value is
  # likely, so it is never worse than plain Gibbs (Neal 2024).
  k0 <- match(x0, vals)
  q <- cp
  cur <- !is.na(k0)
  q[cbind(which(cur), k0[cur])] <- 0
  rs <- rowSums(q)
  pick <- 1L + rowSums(stats::runif(b$n) > t(apply(q / pmax(rs, 1e-300), 1L, cumsum)))
  pick <- pmin(pick, length(vals))
  pc <- ifelse(cur, cp[cbind(seq_len(b$n), ifelse(cur, k0, 1L))], 0)
  pp <- cp[cbind(seq_len(b$n), pick)]
  acc <- rs > 1e-12 & stats::runif(b$n) < pmin(1, (1 - pc) / pmax(1 - pp, 1e-300))
  # a current value outside the support is always replaced
  new <- ifelse(acc | !cur, vals[pick], x0)
  bugs_set(M, b$var, cells, new)
  bugs_write_desc(M, b)
}

# sums of child terms by block member
bugs_by_owner <- function(b, g, v) {
  out <- numeric(b$n)
  if (b$n == 1L) out[1L] <- sum(v)
  else out[g$uowner] <- as.vector(rowsum(v, g$owner, reorder = TRUE))
  out
}

bugs_child_data <- function(M, g) {
  cf <- M$stmts[[g$stmt]]$fn
  list(args = bugs_args_list(cf, g$rows), x = if (cf$mode == "vec") as.list(cf$value(g$rows)) else cf$value(g$rows))
}

# conjugate updates: the posterior is in the prior's family
# the block member each child of a mixture refers to now (NA for none):
# the index expressions are evaluated with the current values
bugs_dyn_owner_fun <- function(M, b, g) {
  cs <- M$stmts[[g$stmt]]
  e <- cs$args[[.bugs_conj[[b$dist$name]][[cs$dist]]]]
  envs <- lapply(g$rows, bugs_inst_env, s = cs)
  own <- M$owner[[e$name]]
  d <- dim(own) %||% length(own)
  ids <- b$ids
  # Returns, per child, the index of the node in b$ids that the slot reads
  # (NA when none) and the slot positions holding its cells. A node may fill
  # part of the slot when the other cells are fixed at zero, as in
  # error[2, 1:2] ~ ddirich() read through error[true[i], ].
  none <- list(NA_integer_, NULL)
  function() {
    out <- lapply(envs, function(env) {
      lin <- 1; stride <- 1
      for (k in seq_along(e$index)) {
        ii <- e$index[[k]]
        i <- switch(ii$kind, single = bugs_const(ii$e, env, M$S), all = seq_len(d[k]),
                    range = bugs_const(ii$lo, env, M$S):bugs_const(ii$hi, env, M$S))
        if (any(i < 1 | i > d[k])) return(none)
        lin <- rep.int(lin, length(i)) + rep((i - 1) * stride, each = length(lin))
        stride <- stride * d[k]
      }
      o <- own[lin]; id <- max(o)
      if (id <= 0L || is.na(j <- match(id, ids))) return(none)
      pos <- which(o == id)
      if (length(pos) != length(M$cells[[id]])) return(none)
      if (length(pos) < length(o) && (any(o[-pos] != 0L) || any(get(e$name, envir = M$S)[lin[-pos]] != 0))) return(none)
      list(j, if (length(pos) < length(o)) pos)
    })
    list(owner = vapply(out, `[[`, 1L, 1L), pos = lapply(out, `[[`, 2L))
  }
}

bugs_dyn_groups <- function(M, b) {
  out <- lapply(seq_along(b$cgroups), function(k) {
    g <- b$cgroups[[k]]
    r <- b$dyn_owner[[k]](); o <- r$owner
    keep <- !is.na(o)
    g$rows <- g$rows[keep]; g$ids <- g$ids[keep]; g$owner <- o[keep]; g$uowner <- sort(unique(o[keep]))
    g$pos <- if (!all(vapply(r$pos[keep], is.null, NA))) r$pos[keep]
    g
  })
  Filter(function(g) length(g$rows) > 0L, out)
}

bugs_conjugate <- function(M, b) {
  if (isTRUE(b$dyn)) b$cgroups <- bugs_dyn_groups(M, b)
  fn <- M$stmts[[b$stmt]]$fn
  a <- bugs_args_list(fn, b$rows)
  nm <- b$dist$name
  if (nm == "dgamma") {
    shape <- vapply(a, `[[`, 0, 1L); rate <- vapply(a, `[[`, 0, 2L)
    for (g in b$cgroups) {
      cd <- bugs_child_data(M, g); cs <- M$stmts[[g$stmt]]
      x <- unlist(cd$x)
      add <- switch(cs$dist,
        dnorm = list(0.5, 0.5 * (x - vapply(cd$args, `[[`, 0, 1L))^2),
        dlnorm = list(0.5, 0.5 * (log(x) - vapply(cd$args, `[[`, 0, 1L))^2),
        dpois = list(x, 1),
        dexp = list(1, x),
        dgamma = list(vapply(cd$args, `[[`, 0, 1L), x))
      shape <- shape + bugs_by_owner(b, g, rep_len(add[[1L]], length(x)))
      rate <- rate + bugs_by_owner(b, g, rep_len(add[[2L]], length(x)))
    }
    v <- stats::rgamma(b$n, shape, rate)
  } else if (nm == "dbeta") {
    s1 <- vapply(a, `[[`, 0, 1L); s2 <- vapply(a, `[[`, 0, 2L)
    for (g in b$cgroups) {
      cd <- bugs_child_data(M, g); cs <- M$stmts[[g$stmt]]
      x <- unlist(cd$x)
      n <- if (cs$dist == "dbern") 1 else vapply(cd$args, `[[`, 0, 2L)
      s1 <- s1 + bugs_by_owner(b, g, x); s2 <- s2 + bugs_by_owner(b, g, rep_len(n - x, length(x)))
    }
    v <- stats::rbeta(b$n, s1, s2)
  } else if (nm == "ddirch") {
    A <- matrix(unlist(lapply(a, function(x) as.numeric(x[[1L]]))), b$n, b$K, byrow = TRUE)
    for (g in b$cgroups) {
      cd <- bugs_child_data(M, g); cs <- M$stmts[[g$stmt]]
      for (j in seq_along(g$rows)) {
        o <- g$owner[j]; x <- cd$x[[j]]; pos <- g$pos[[j]]
        if (cs$dist == "dcat") {
          if (!is.null(pos)) x <- match(x, pos)
          A[o, x] <- A[o, x] + 1
        } else A[o, ] <- A[o, ] + (if (is.null(pos)) x else x[pos])
      }
    }
    G <- matrix(stats::rgamma(length(A), A), nrow(A))
    G[A == 0] <- 0
    v <- t(G / rowSums(G))
  } else if (nm == "dwish") {
    R <- a[[1L]][[1L]]; k <- a[[1L]][[2L]]
    for (g in b$cgroups) {
      cd <- bugs_child_data(M, g)
      for (j in seq_along(g$rows)) { r <- cd$x[[j]] - cd$args[[j]][[1L]]; R <- R + tcrossprod(r); k <- k + 1 }
    }
    v <- stats::rWishart(1L, k, solve(R))[, , 1L]
  } else if (nm == "dnorm") {
    # a normal mean or a coefficient: child means are m0 + c x
    cells <- b$cells[, 1L]
    x0 <- get(b$var, envir = M$S)[cells]
    mu <- vapply(a, `[[`, 0, 1L); tau <- vapply(a, `[[`, 0, 2L)
    means_at <- function(v) {
      bugs_set(M, b$var, cells, rep_len(v, b$n)); bugs_write_desc(M, b)
      lapply(b$cgroups, function(g) vapply(bugs_args_list(M$stmts[[g$stmt]]$fn, g$rows), `[[`, 0, 1L))
    }
    m0 <- means_at(0); m1 <- means_at(1)
    prec <- tau; num <- tau * mu
    for (k in seq_along(b$cgroups)) {
      g <- b$cgroups[[k]]
      cd <- bugs_child_data(M, g)
      ct <- vapply(cd$args, `[[`, 0, 2L)
      cc <- m1[[k]] - m0[[k]]
      prec <- prec + bugs_by_owner(b, g, ct * cc^2)
      num <- num + bugs_by_owner(b, g, ct * cc * (unlist(cd$x) - m0[[k]]))
    }
    v <- stats::rnorm(b$n, num / prec, 1 / sqrt(prec))
    if (any(!is.finite(v))) stop_sampler(sprintf("the conjugate normal update of %s failed", b$var))
  }
  bugs_set(M, b$var, as.vector(t(b$cells)), as.vector(v))
  bugs_write_desc(M, b)
}

# elementwise slice sampling of multivariate normal and t nodes; observed
# elements of partly observed nodes stay fixed
bugs_mv_slice <- function(M, b) {
  for (j in seq_len(b$K)) {
    sel <- which(!b$fixed[, j])
    if (!length(sel)) next
    jump <- bugs_slice_coord(M, b, b$cells[, j], sel, b$w[sel, j])
    b$w[sel, j] <- bugs_adapt_w(b$w[sel, j], jump, b$adapting)
  }
  b$w
}

# Dirichlet nodes: slice sampling of z_j = log(x_j / x_r), r the last
# element with a positive parameter, with the Jacobian prod(x); elements
# whose parameter is zero stay at zero
bugs_dirichlet <- function(M, b) {
  fn <- M$stmts[[b$stmt]]$fn
  A <- matrix(unlist(lapply(bugs_args_list(fn, b$rows), function(x) as.numeric(x[[1L]]) > 0)), b$n, b$K, byrow = TRUE)
  ref <- apply(A, 1L, function(r) max(which(r)))
  cells <- as.vector(t(b$cells))
  for (j in seq_len(b$K)) {
    sel <- which(A[, j] & ref != j)
    if (!length(sel)) next
    X <- matrix(get(b$var, envir = M$S)[cells], b$n, b$K, byrow = TRUE)
    Z <- log(X / X[cbind(seq_len(b$n), ref)])
    Z[!A] <- -Inf
    setz <- function(z) {
      Zt <- Z; Zt[sel, j] <- z
      E <- exp(Zt - apply(Zt, 1L, max)); E[!A] <- 0
      Y <- E / rowSums(E)
      bugs_set(M, b$var, cells, as.vector(t(Y)))
      rowSums(ifelse(A, log(Y), 0))
    }
    f <- function(z) { jac <- setz(z); (bugs_cond(M, b) + jac)[sel] }
    z1 <- bugs_slice_core(f, Z[sel, j], rep(1, length(sel)), what = b$var)
    setz(z1)
  }
  bugs_write_desc(M, b)
}

# multinomial nodes: for a random pair of categories the two counts are
# redrawn from their full conditional with their sum fixed
bugs_multinomial <- function(M, b) {
  if (b$K < 2L) return(invisible())
  for (step in seq_len(b$K - 1L)) {
    jk <- sample.int(b$K, 2L)
    cj <- b$cells[, jk[1L]]; ck <- b$cells[, jk[2L]]
    a <- get(b$var, envir = M$S)
    xj <- a[cj]; xk <- a[ck]; tot <- xj + xk
    if (max(tot) <= 100) {
      vals <- 0:max(tot)
      lp <- matrix(-Inf, b$n, length(vals))
      for (t in seq_along(vals)) {
        inr <- vals[t] <= tot
        bugs_set(M, b$var, cj, ifelse(inr, vals[t], xj)); bugs_set(M, b$var, ck, ifelse(inr, tot - vals[t], xk))
        lp[inr, t] <- bugs_cond(M, b)[inr]
      }
      mx <- apply(lp, 1L, max)
      p <- exp(lp - mx); cp <- p / rowSums(p)
      pick <- vals[pmin(1L + rowSums(stats::runif(b$n) > t(apply(cp, 1L, cumsum))), length(vals))]
    } else {
      f <- function(z) {
        t <- pmin(pmax(floor(z), 0), tot)
        bugs_set(M, b$var, cj, t); bugs_set(M, b$var, ck, tot - t)
        l <- bugs_cond(M, b); l[z < 0 | z >= tot + 1] <- -Inf; l
      }
      pick <- floor(bugs_slice_core(f, xj + stats::runif(b$n), rep(2, b$n), what = b$var))
    }
    bugs_set(M, b$var, cj, pick); bugs_set(M, b$var, ck, tot - pick)
  }
  bugs_write_desc(M, b)
}

bugs_update <- function(M, b) {
  switch(b$sampler,
    slice = { jump <- bugs_slice_coord(M, b, b$cells[, 1L], seq_len(b$n), b$w); b$w <- bugs_adapt_w(b$w, jump, b$adapting) },
    discrete_slice = { jump <- bugs_slice_coord(M, b, b$cells[, 1L], seq_len(b$n), b$w, discrete = TRUE); b$w <- pmax(1, bugs_adapt_w(b$w, jump, b$adapting)) },
    enumerate = bugs_enumerate(M, b, b$cells[, 1L]),
    mv_slice = { b$w <- bugs_mv_slice(M, b) },
    dirichlet = bugs_dirichlet(M, b),
    multinomial = bugs_multinomial(M, b),
    bugs_conjugate(M, b))
  b
}

# ---- setting up ---------------------------------------------------------------------------

# logical nodes with a stochastic descendant
bugs_feeds <- function(M) {
  n <- length(M$nodes$kind)
  f <- logical(n)
  for (i in rev(M$order)) if (M$nodes$kind[i] == "logical") {
    ch <- M$children[[i]]
    f[i] <- any(M$nodes$kind[ch] == "stochastic") || any(f[ch])
  }
  f
}

bugs_blocks <- function(M) {
  M$feeds <- bugs_feeds(M)
  U <- which(M$nodes$kind == "stochastic" & !M$nodes$obs)
  U <- U[order(M$rank[U])]
  desc <- lapply(U, function(i) bugs_desc(M, i))
  names(desc) <- U
  blocks <- list()
  for (k in unique(M$nodes$stmt[U])) {
    ids <- U[M$nodes$stmt[U] == k]
    dd <- desc[as.character(ids)]
    for (cl in bugs_color(M, ids, dd)) blocks[[length(blocks) + 1L]] <- bugs_make_block(M, ids[cl], dd[cl])
  }
  # OpenBUGS updates nodes in the order the model declares them; JAGS in
  # the order of the graph
  key <- if (M$dialect == "openbugs") vapply(blocks, function(b) b$stmt * 1e9 + min(M$rank[b$ids]), 0)
         else vapply(blocks, function(b) min(M$rank[b$ids]), 1L)
  blocks[order(key)]
}

# initial values: given ones, then the typical value of the prior (JAGS) or
# a draw from it (OpenBUGS gen inits), in graph order
bugs_initialize <- function(M, inits, how = c("typical", "prior")) {
  how <- match.arg(how)
  inits <- inits %||% list()
  # JAGS's .RNG.name and .RNG.seed choose generators; the chain streams of
  # bayescompR are used instead
  inits <- inits[!startsWith(names(inits) %||% character(), ".")]
  inits <- bugs_check_values(inits, "inits")
  for (nm in names(inits)) {
    if (!nm %in% names(M$dims)) { warning(sprintf("the initial value of '%s' is ignored: not a variable of the model", nm), call. = FALSE); next }
    v <- inits[[nm]]
    a <- get(nm, envir = M$S)
    if (length(v) != length(a)) stop_data(sprintf("the initial value of '%s' has %d values but the variable has %d", nm, length(v), length(a)))
    ok <- !is.na(v)
    # only missing elements of unobserved stochastic nodes take initial values
    own <- M$owner[[nm]]
    settable <- own != 0L & M$nodes$kind[pmax(own, 1L)] == "stochastic" & !M$nodes$obs[pmax(own, 1L)] & is.na(a)
    if (any(ok & !settable & is.na(a))) warning(sprintf("initial values of '%s' for logical nodes are ignored", nm), call. = FALSE)
    if (any(ok & !is.na(a))) warning(sprintf("initial values of '%s' for data are ignored", nm), call. = FALSE)
    a[ok & settable] <- as.numeric(v)[ok & settable]
    assign(nm, a, envir = M$S)
  }
  given <- lapply(stats::setNames(nm = unique(M$nodes$var)), function(v) get(v, envir = M$S))
  # draws from vague priors can have zero probability (a precision drawn
  # from dgamma(0.001, 0.001) is often 0) and typical values can be
  # impossible: try both before giving up
  plan <- if (how == "prior") c(rep("prior", 20L), "typical") else c("typical", rep("prior", 20L))
  repair <- NULL
  for (use in plan) {
    for (v in names(given)) assign(v, given[[v]], envir = M$S)
    bugs_init_fill(M, use)
    if (bugs_init_ok(M)) return(invisible())
    # discrete latent nodes drawn from their priors can contradict the data
    # (a category whose error probability is zero): one Gibbs sweep over
    # them picks values the data allow, as the first update of OpenBUGS does
    repair <- repair %||% Filter(function(b) b$sampler == "enumerate", bugs_blocks(M))
    if (length(repair) && !is.null(tryCatch({ for (b in repair) bugs_update(M, b); TRUE }, error = function(e) NULL)) && bugs_init_ok(M))
      return(invisible())
  }
  for (v in names(given)) assign(v, given[[v]], envir = M$S)
  bugs_init_fill(M, how)
  for (i in M$order) if (M$nodes$kind[i] == "stochastic") {
    lp <- M$stmts[[M$nodes$stmt[i]]]$fn$lp(M$nodes$row[i])
    if (!is.finite(sum(lp))) stop_sampler(sprintf("%s has zero probability at the initial values; give initial values", bugs_node_label(M, i)))
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
    if (M$nodes$kind[i] == "logical") { if (!M$const[i]) s$fn$write(r); next }
    if (M$nodes$obs[i]) next
    cur <- get(M$nodes$var[i], envir = M$S)[M$cells[[i]]]
    if (all(!is.na(cur))) next
    d <- bugs_dist_get(s$dist, M$dialect)
    v <- tryCatch(unlist(if (how == "prior") s$fn$rng(r) else s$fn$typical(r)), error = function(e) rep(NA_real_, length(cur)))
    if (isTRUE(d$discrete)) v <- round(v)
    if (!is.null(s$bound) && length(v) == 1L) {
      bb <- s$fn$bounds(r); if (s$fn$mode != "vec") bb <- bb[[1L]]
      lo <- bb[[1L]]; hi <- bb[[2L]]
      if (isTRUE(d$discrete)) { lo <- ceiling(lo); hi <- floor(hi) }
      if (is.na(v) || v < lo || v > hi) v <- if (is.finite(lo) && is.finite(hi)) { if (isTRUE(d$discrete)) lo else (lo + hi) / 2 }
                                         else if (is.finite(lo)) { if (isTRUE(d$discrete)) lo else lo + 1 } else hi - 1
    }
    v <- rep_len(v, length(cur))
    v[!is.na(cur)] <- cur[!is.na(cur)]
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

# The monitored values: a variable name ("b"), or a part of one ("b[2]",
# "Y[1:3, 2]"). Elements that the model never defines are left out.
bugs_monitor <- function(M, vars) {
  out <- list()
  for (v in setdiff(vars, "deviance")) {
    e <- tryCatch(bugs_parse(sprintf("model { monitor <- %s }", v))$model[[1L]]$expr, error = function(err) NULL)
    if (is.null(e) || !identical(e$kind, "var")) stop_semantic(sprintf("cannot monitor '%s'", v))
    d <- M$dims[[e$name]]
    if (is.null(d)) stop_semantic(sprintf("cannot monitor '%s': not a variable of the model", v))
    ix <- if (is.null(e$index)) lapply(d, seq_len) else {
      if (length(e$index) != length(d)) stop_semantic(sprintf("cannot monitor '%s': '%s' has %d dimensions", v, e$name, length(d)))
      lapply(seq_along(e$index), function(k) {
        ii <- e$index[[k]]
        r <- switch(ii$kind, all = seq_len(d[k]), single = bugs_const(ii$e, list(), M$known),
                    range = bugs_const(ii$lo, list(), M$known):bugs_const(ii$hi, list(), M$known))
        if (any(r < 1 | r > d[k])) stop_semantic(sprintf("cannot monitor '%s': index out of range", v))
        r
      })
    }
    g <- as.matrix(expand.grid(ix))
    cells <- as.integer(1L + colSums((t(g) - 1L) * cumprod(c(1L, d[-length(d)]))))
    defined <- M$owner[[e$name]][cells] != 0L | !is.na(get(e$name, envir = M$S)[cells])
    whole <- is.null(e$index) && all(defined)
    cells <- cells[defined]; g <- g[defined, , drop = FALSE]
    if (!length(cells)) stop_semantic(sprintf("cannot monitor '%s': no element is defined", v))
    scalar <- prod(d) == 1L && is.null(e$index) && length(d) == 1L && is.null(M$data[[e$name]])
    nm <- if (scalar) e$name else paste0(e$name, "[", apply(g, 1L, paste, collapse = ","), "]")
    key <- e$name
    if (!is.null(out[[key]])) {
      keep <- !cells %in% out[[key]]$cells
      out[[key]]$cells <- c(out[[key]]$cells, cells[keep]); out[[key]]$names <- c(out[[key]]$names, nm[keep])
      out[[key]]$index <- rbind(out[[key]]$index, g[keep, , drop = FALSE]); out[[key]]$dim <- NULL
    } else out[[key]] <- list(var = e$name, cells = cells, names = nm, index = g, dim = if (whole && !scalar) d)
  }
  out
}

# the state of one chain: the compiled model with its values and samplers
bugs_chain_init <- function(M, inits, init_how) {
  bugs_initialize(M, inits, init_how)
  blocks <- bugs_blocks(M)
  rest <- which(M$nodes$kind == "logical" & !M$const & !M$feeds)
  list(M = M, blocks = blocks, rest = bugs_logical_groups(M, rest))
}
bugs_write_rest <- function(st) for (g in st$rest) st$M$stmts[[g$stmt]]$fn$write(g$rows)

# run n iterations; every thin-th one after `skip` is recorded, and the
# running sums of the unobserved variables are kept in `means` if given
bugs_chain_run <- function(st, n, adapting, mon, deviance, skip = 0L, thin = 1L, means = NULL, refresh = 0L, chain = 1L) {
  M <- st$M
  devf <- if (deviance) bugs_deviance_fun(M)
  vnames <- c(unlist(lapply(mon, `[[`, "names"), use.names = FALSE), if (deviance) "deviance")
  keep <- skip + thin * seq_len(max(0L, (n - skip) %/% thin))
  out <- matrix(NA_real_, length(keep), length(vnames), dimnames = list(NULL, vnames))
  k <- 0L
  for (t in seq_len(n)) {
    for (j in seq_along(st$blocks)) {
      st$blocks[[j]]$adapting <- adapting
      st$blocks[[j]] <- bugs_update(M, st$blocks[[j]])
    }
    if (k < length(keep) && t == keep[k + 1L]) {
      k <- k + 1L
      bugs_write_rest(st)
      row <- unlist(lapply(mon, function(m) get(m$var, envir = M$S)[m$cells]), use.names = FALSE)
      if (deviance) row <- c(row, devf())
      out[k, ] <- row
      if (!is.null(means)) {
        for (v in names(means$sums)) means$sums[[v]] <- means$sums[[v]] + get(v, envir = M$S)
        means$k <- means$k + 1L
      }
    }
    if (refresh > 0L && t %% refresh == 0L) message(sprintf("chain %d: iteration %d / %d", chain, t, n))
  }
  bugs_write_rest(st)
  list(state = st, draws = out)
}

bugs_means_env <- function(M) {
  means <- new.env(parent = emptyenv())
  uvars <- unique(M$nodes$var[M$nodes$kind == "stochastic" & !M$nodes$obs])
  means$sums <- lapply(stats::setNames(nm = uvars), function(v) 0)
  means$k <- 0L
  means
}

# one chain: n_adapt adaptation iterations, then n_iter iterations of which
# those after n_burnin are kept every n_thin
bugs_chain <- function(M, inits, n_adapt, n_burnin, n_iter, n_thin, monitor, deviance, init_how, refresh = 0L, chain = 1L) {
  st <- bugs_chain_init(M, inits, init_how)
  mon <- bugs_monitor(M, monitor)
  if (n_adapt > 0L) st <- bugs_chain_run(st, n_adapt, TRUE, list(), FALSE, n_adapt, 1L, refresh = refresh, chain = chain)$state
  means <- bugs_means_env(M)
  r <- bugs_chain_run(st, n_iter, FALSE, mon, deviance, n_burnin, n_thin, means, refresh, chain)
  list(draws = r$draws, means = lapply(means$sums, function(x) x / max(means$k, 1L)),
       samplers = table(vapply(r$state$blocks, `[[`, "", "sampler")))
}
