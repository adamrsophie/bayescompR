# R code for the statements of a compiled BUGS model (see bugs_compile.R).
#
# A statement becomes functions of `rows`, the instances to evaluate:
#   logical:    value(rows), write(rows)
#   stochastic: lp(rows)      log densities, with -Inf outside the support
#                             and the normalizer of a T() truncation
#               args(rows)    the parameters
#               rng(rows), typical(rows), support(rows)
# In "vec" mode the loop variables are vectors over the rows and every
# operation is elementwise; in "loop" mode the rows are evaluated one at a
# time (vector valued functions, ranges, multivariate distributions).

# ---- functions --------------------------------------------------------------------------
# elementwise: whether the function acts on each element (so that it can be
# applied to vectors over instances); dialect differences are noted

.bugs_fun <- new.env(parent = emptyenv())
bfun <- function(names, r, elementwise = TRUE) for (nm in names) assign(nm, list(r = r, elementwise = elementwise), envir = .bugs_fun)
icloglog <- function(x) -expm1(-exp(x))
bfun("abs", abs); bfun("sqrt", sqrt); bfun("exp", exp); bfun("log", log); bfun("log10", log10)
bfun("cos", cos); bfun("sin", sin); bfun("tan", tan); bfun("cosh", cosh); bfun("sinh", sinh); bfun("tanh", tanh)
bfun(c("arccos", "acos"), acos); bfun(c("arcsin", "asin"), asin); bfun(c("arctan", "atan"), atan)
bfun(c("arccosh", "acosh"), acosh); bfun(c("arcsinh", "asinh"), asinh); bfun(c("arctanh", "atanh"), atanh)
bfun("logit", stats::qlogis); bfun("ilogit", stats::plogis); bfun("cloglog", function(x) log(-log1p(-x))); bfun("icloglog", icloglog)
bfun("probit", stats::qnorm); bfun(c("phi", "Phi"), stats::pnorm)
bfun("pow", function(a, b) a^b); bfun("step", function(x) as.numeric(x >= 0)); bfun("equals", function(a, b) as.numeric(a == b))
bfun("round", function(x) sign(x) * floor(abs(x) + 0.5)); bfun("loggam", lgamma); bfun("logfact", lfactorial); bfun("lgamma", lgamma); bfun("gamma", gamma)
bfun("ifelse", function(c, a, b) { n <- max(length(c), length(a), length(b)); ifelse(rep_len(c != 0, n), rep_len(a, n), rep_len(b, n)) })
bfun("trunc", trunc)   # OpenBUGS: floor, set per dialect in bugs_fun_get()
bfun("max", function(...) max(...), elementwise = FALSE)
bfun("min", function(...) min(...), elementwise = FALSE)
bfun("sum", function(...) sum(...), elementwise = FALSE); bfun("mean", function(x) mean(x), elementwise = FALSE)
bfun("prod", function(x) prod(x), elementwise = FALSE); bfun("sd", function(x) stats::sd(x), elementwise = FALSE)
bfun("inprod", function(a, b) sum(a * b), elementwise = FALSE)
bfun("inverse", function(a) solve(a), elementwise = FALSE)
bfun("logdet", function(a) as.numeric(determinant(a, logarithm = TRUE)$modulus), elementwise = FALSE)
bfun("sort", function(x) sort(x), elementwise = FALSE)
# rank(v[]) gives the ranks (JAGS, ties in order); rank(v[], s) the number
# of elements of v that are at most v[s] (OpenBUGS)
bfun("rank", function(x, s) { if (missing(s)) rank(x, ties.method = "first") else sum(x <= x[s]) }, elementwise = FALSE)
bfun("ranked", function(x, s) sort(x)[s], elementwise = FALSE)
bfun("order", function(x) order(x), elementwise = FALSE)
bfun("t", function(x) t(x), elementwise = FALSE)
bfun("length", function(x) length(x), elementwise = FALSE)
bfun("rep", function(x, n) rep(x, n), elementwise = FALSE)
bfun("interp.lin", function(e, x, y) stats::approx(x, y, e, rule = 2)$y, elementwise = FALSE)
bfun("cut", function(x) x)   # cut() only changes how parents are updated (see the samplers)

bugs_fun_get <- function(name, loc = NULL, dialect = "jags") {
  f <- .bugs_fun[[name]]
  if (is.null(f)) stop_unsupported(sprintf("the function '%s' is not supported", name), loc)
  if (name == "trunc" && dialect == "openbugs") f$r <- floor
  f
}

# probabilities from link functions are inside (0, 1); keep them there when
# they round to 0 or 1, as at extreme initial values
open01 <- function(p) pmin(pmax(p, .Machine$double.xmin), 1 - .Machine$double.eps / 2)
bugs_inv_link <- list(log = exp, logit = function(x) open01(stats::plogis(x)), cloglog = function(x) open01(icloglog(x)),
                      probit = function(x) open01(stats::pnorm(x)))

# log cumulative distribution functions for T() truncation: log P(X <= q)
# with lower = TRUE, log P(X > q) otherwise
bugs_cdf <- list(
  dnorm = function(q, lower, mu, tau) stats::pnorm(q, mu, sd_of(tau), lower.tail = lower, log.p = TRUE),
  dlnorm = function(q, lower, mu, tau) stats::plnorm(q, mu, sd_of(tau), lower.tail = lower, log.p = TRUE),
  dt = function(q, lower, mu, tau, k) stats::pt((q - mu) * sqrt(tau), k, lower.tail = lower, log.p = TRUE),
  dlogis = function(q, lower, mu, tau) stats::plogis(q, mu, 1 / tau, lower.tail = lower, log.p = TRUE),
  ddexp = function(q, lower, mu, tau) {
    near <- log(0.5) - tau * abs(q - mu)
    far <- log1p(-0.5 * exp(-tau * abs(q - mu)))
    if (lower) ifelse(q < mu, near, far) else ifelse(q < mu, far, near)
  },
  dgamma = function(q, lower, r, mu) stats::pgamma(q, r, mu, lower.tail = lower, log.p = TRUE),
  dexp = function(q, lower, lambda) stats::pexp(q, lambda, lower.tail = lower, log.p = TRUE),
  dchisqr = function(q, lower, k) stats::pchisq(q, k, lower.tail = lower, log.p = TRUE),
  dweib = function(q, lower, v, lambda) {
    h <- lambda * pmax(q, 0)^v
    if (lower) ifelse(q > 0, log(-expm1(-h)), -Inf) else ifelse(q > 0, -h, 0)
  },
  dbeta = function(q, lower, a, b) stats::pbeta(q, a, b, lower.tail = lower, log.p = TRUE),
  dunif = function(q, lower, a, b) stats::punif(q, a, b, lower.tail = lower, log.p = TRUE),
  dpar = function(q, lower, alpha, c) {
    r <- pmin(c / pmax(q, c), 1)
    if (lower) ifelse(q > c, log1p(-r^alpha), -Inf) else ifelse(q > c, alpha * log(r), 0)
  },
  dpois = function(q, lower, lambda) stats::ppois(q, lambda, lower.tail = lower, log.p = TRUE),
  dbin = function(q, lower, p, n) stats::pbinom(q, n, p, lower.tail = lower, log.p = TRUE),
  dnegbin = function(q, lower, p, r) stats::pnbinom(q, r, p, lower.tail = lower, log.p = TRUE),
  dbern = function(q, lower, p) stats::pbinom(q, 1, p, lower.tail = lower, log.p = TRUE))
bugs_cdf$dchisq <- bugs_cdf$dchisqr
bugs_cdf$dbinom <- bugs_cdf$dbin

# log(1 - exp(x)) for x <= 0
log1mexp <- function(x) ifelse(x > -log(2), log(-expm1(x)), log1p(-exp(x)))

# log P(lo <= X <= hi), computed from the tail where it is accurate; for a
# discrete X the lower limit is the smallest admitted integer
bugs_log_mass <- function(cdf, lo, hi, a, discrete) {
  if (discrete) lo <- ceiling(lo) - 1
  n <- max(length(lo), length(hi))
  lo <- rep_len(lo, n); hi <- rep_len(hi, n)
  out <- numeric(n)
  up <- is.infinite(hi) & hi > 0
  dn <- is.infinite(lo) & lo < 0
  if (any(up & !dn)) out[up & !dn] <- do.call(cdf, c(list(lo, FALSE), a))[up & !dn]
  if (any(dn & !up)) out[dn & !up] <- do.call(cdf, c(list(hi, TRUE), a))[dn & !up]
  two <- !up & !dn
  if (any(two)) {
    flo <- do.call(cdf, c(list(lo, TRUE), a)); fhi <- do.call(cdf, c(list(hi, TRUE), a))
    slo <- do.call(cdf, c(list(lo, FALSE), a)); shi <- do.call(cdf, c(list(hi, FALSE), a))
    upper <- flo > log(0.5)
    m <- ifelse(upper, slo + log1mexp(pmin(shi - slo, 0)), fhi + log1mexp(pmin(flo - fhi, 0)))
    out[two] <- rep_len(m, n)[two]
  }
  out
}

# ---- translation ------------------------------------------------------------------------

bugs_loop_name <- function(v) as.name(paste0(".L_", v))

# can the expression be evaluated for many instances at once?
bugs_is_elementwise <- function(e, M) {
  switch(e$kind,
    num = TRUE,
    var = is.null(e$index) && (is.null(M$dims[[e$name]]) || prod(M$dims[[e$name]]) == 1L) ||
          (!is.null(e$index) && all(vapply(e$index, function(ix) identical(ix$kind, "single") && bugs_is_elementwise(ix$e, M), TRUE))),
    unop = bugs_is_elementwise(e$a, M),
    binop = !grepl("^%", e$op) && bugs_is_elementwise(e$a, M) && bugs_is_elementwise(e$b, M),
    call = isTRUE(.bugs_fun[[e$name]]$elementwise) && all(vapply(e$args, bugs_is_elementwise, TRUE, M = M)),
    FALSE)
}

bugs_tr <- function(e, M, loopvars, vec) {
  tr <- function(e) bugs_tr(e, M, loopvars, vec)
  switch(e$kind,
    num = e$value,
    var = {
      if (is.null(e$index) && e$name %in% loopvars) return(bugs_loop_name(e$name))
      sv <- call("$", as.name("S"), as.name(e$name))
      d <- M$dims[[e$name]]
      if (is.null(e$index)) return(if (prod(d) == 1L) call("[", sv, 1L) else call(".bugs_whole", sv))
      if (vec) {
        ix <- lapply(e$index, function(i) tr(i$e))
        if (length(ix) == 1L) return(call("[", sv, ix[[1L]]))
        return(call("[", sv, as.call(c(list(as.name("cbind")), ix))))
      }
      ix <- lapply(seq_along(e$index), function(k) {
        i <- e$index[[k]]
        switch(i$kind, single = tr(i$e), range = call(":", tr(i$lo), tr(i$hi)), all = NULL)
      })
      as.call(c(list(as.name(".bugs_sub"), sv), lapply(ix, function(x) if (is.null(x)) quote(NULL) else x)))
    },
    unop = if (e$op == "-") call("-", tr(e$a)) else call("as.numeric", call("!", tr(e$a))),
    binop = {
      op <- switch(e$op, "&&" = "&", "||" = "|", e$op)
      if (!op %in% c("+", "-", "*", "/", "^", "<", "<=", ">", ">=", "==", "!=", "&", "|", "%*%"))
        stop_unsupported(sprintf("the operator '%s' is not supported", e$op), e$loc)
      r <- call(op, tr(e$a), tr(e$b))
      if (op %in% c("<", "<=", ">", ">=", "==", "!=", "&", "|")) call("as.numeric", r) else r
    },
    call = {
      f <- bugs_fun_get(e$name, e$loc, M$dialect)
      as.call(c(list(call("$", as.name(".bf"), as.name(e$name))), lapply(e$args, tr)))
    },
    stop_semantic("unsupported expression", e$loc))
}

# x[i, 1:K] with drop, as an R array; a missing index is NULL
# (linear indices are computed directly: this runs once per node and row)
.bugs_sub <- function(a, ...) {
  ix <- list(...)
  d <- dim(a)
  if (is.null(d)) d <- length(a)
  lin <- 1; stride <- 1; keep <- integer()
  for (k in seq_along(ix)) {
    i <- ix[[k]]
    i <- if (is.null(i)) seq_len(d[k]) else as.integer(i)
    if (any(i < 1L | i > d[k], na.rm = TRUE)) stop("subscript out of bounds", call. = FALSE)
    if (length(i) == 1L) lin <- lin + (i - 1L) * stride
    else { keep <- c(keep, length(i)); lin <- rep.int(lin, length(i)) + rep((i - 1L) * stride, each = length(lin)) }
    stride <- stride * d[k]
  }
  r <- a[lin]
  if (length(keep) <= 1L) r else array(r, keep)
}
.bugs_whole <- function(a) if (length(dim(a)) == 1L) as.vector(a) else a

# ---- statements -------------------------------------------------------------------------

bugs_codegen <- function(M, k) {
  s <- M$stmts[[k]]
  loopvars <- colnames(s$inst) %||% character()
  ex <- bugs_stmt_exprs(s)
  ids <- s$nodes
  single_lhs <- all(lengths(M$cells[ids]) == 1L)
  d <- if (s$kind == "stochastic") bugs_dist_get(s$dist, M$dialect, s$loc)
  if (s$kind == "stochastic" && length(s$args) != d$npar)
    stop_semantic(sprintf("%s takes %d parameters but %d were given", s$dist, d$npar, length(s$args)), s$loc)
  vec <- single_lhs && all(vapply(ex, bugs_is_elementwise, TRUE, M = M)) && (s$kind == "logical" || identical(d$dim, "scalar"))
  fenv <- new.env(parent = asNamespace("bayescompR"))
  fenv$S <- M$S
  fenv$.bf <- lapply(as.list(.bugs_fun), function(f) if (identical(f, .bugs_fun$trunc) && M$dialect == "openbugs") floor else f$r)
  fenv$inst <- s$inst
  fenv$cells <- M$cells[ids]
  fenv$cellv <- if (single_lhs) unlist(M$cells[ids]) else NULL
  fenv$var <- s$lhs$name
  fenv$d <- d
  fenv$shape <- if (!single_lhs) lapply(s$lhs_idx, function(ix) { l <- lengths(ix); l[l > 1L] })
  bind <- lapply(seq_along(loopvars), function(j) call("<-", bugs_loop_name(loopvars[j]), call("[", as.name("inst"), as.name("rows"), j)))
  mk <- function(body, env = fenv) {
    f <- function(rows) NULL
    body(f) <- as.call(c(list(as.name("{")), bind, body))
    environment(f) <- env
    f
  }
  mk_loop <- function(per_row) {
    # per_row: a body evaluated with scalar loop variables for row `r`
    f <- function(rows) lapply(rows, function(r) NULL)
    inner <- function(r) NULL
    body(inner) <- as.call(c(list(as.name("{")),
                             lapply(seq_along(loopvars), function(j) call("<-", bugs_loop_name(loopvars[j]), call("[", as.name("inst"), as.name("r"), j))),
                             per_row))
    environment(inner) <- fenv
    # each closure carries its own row evaluator
    body(f) <- bquote(lapply(rows, .(inner)))
    environment(f) <- fenv
    f
  }
  out <- list(mode = if (vec) "vec" else "loop", kind = s$kind)
  if (s$kind == "logical") {
    rhs <- bugs_tr(s$expr, M, loopvars, vec)
    if (!is.null(s$link)) rhs <- call(".il", rhs)
    fenv$.il <- bugs_inv_link[[s$link %||% "log"]]
    if (vec) {
      out$value <- mk(list(rhs))
      fenv$value <- out$value
      out$write <- function(rows) { v <- value(rows); S[[var]][cellv[rows]] <- rep_len(v, length(rows)); invisible() }
      environment(out$write) <- fenv
    } else {
      out$value <- mk_loop(list(rhs))
      fenv$value <- out$value
      out$write <- function(rows) {
        vals <- value(rows)
        a <- S[[var]]
        for (j in seq_along(rows)) {
          v <- vals[[j]]
          if (length(v) != length(cells[[rows[j]]])) stop_shape(sprintf("the value of %s has %d elements for %d places", var, length(v), length(cells[[rows[j]]])))
          a[cells[[rows[j]]]] <- v
        }
        assign(var, a, envir = S)
        invisible()
      }
      environment(out$write) <- fenv
    }
    return(out)
  }
  # stochastic
  args <- lapply(s$args, bugs_tr, M = M, loopvars = loopvars, vec = vec)
  lo <- if (!is.null(s$bound$lo)) bugs_tr(s$bound$lo, M, loopvars, vec)
  hi <- if (!is.null(s$bound$hi)) bugs_tr(s$bound$hi, M, loopvars, vec)
  fenv$btype <- s$bound$type %||% ""
  # a flat prior restricted by T() has no normalizer
  if (identical(fenv$btype, "T") && s$dist == "dflat") fenv$btype <- "I"
  # OpenBUGS: C() only restricts missing values; an observed value is used
  # as is (cobs is set once the observed nodes are known)
  fenv$cobs <- logical(length(ids))
  if (identical(fenv$btype, "T")) {
    fenv$cdf <- bugs_cdf[[s$dist]]
    if (is.null(fenv$cdf)) stop_unsupported(sprintf("truncation of %s", s$dist), s$loc)
  }
  fenv$discrete <- isTRUE(d$discrete)
  arg_call <- as.call(c(list(as.name("list")), args))
  lo_call <- if (is.null(lo)) -Inf else lo
  hi_call <- if (is.null(hi)) Inf else hi
  if (vec) {
    out$args <- mk(list(arg_call))
    out$bounds <- mk(list(call("list", lo_call, hi_call)))
    fenv$args_f <- out$args; fenv$bounds_f <- out$bounds
    out$value <- function(rows) S[[var]][cellv[rows]]
    out$lp <- function(rows, x = S[[var]][cellv[rows]]) {
      a <- args_f(rows)
      lp <- suppressWarnings(do.call(d$lpdf, c(list(x), a)))
      if (nzchar(btype)) {
        b <- bounds_f(rows)
        out <- x < b[[1L]] | x > b[[2L]]
        if (btype == "C") out <- out & !cobs[rows]
        lp[out] <- -Inf
        if (btype == "T") lp <- lp - bugs_log_mass(cdf, b[[1L]], b[[2L]], a, discrete)
      }
      rep_len(lp, length(rows))
    }
    out$rng <- function(rows) { a <- args_f(rows); rep_len(do.call(d$rng, c(list(length(rows)), a)), length(rows)) }
    out$typical <- function(rows) { a <- args_f(rows); rep_len(do.call(d$typical, a), length(rows)) }
  } else {
    out$args <- mk_loop(list(arg_call))
    out$bounds <- mk_loop(list(call("list", lo_call, hi_call)))
    fenv$args_f <- out$args; fenv$bounds_f <- out$bounds
    out$value <- function(rows) lapply(rows, function(r) bugs_shape_value(S[[var]][cells[[r]]], shape[[r]]))
    out$lp <- function(rows, x = NULL) {
      a <- args_f(rows)
      b <- if (nzchar(btype)) bounds_f(rows)
      sv <- if (is.null(x)) S[[var]]
      suppressWarnings(vapply(seq_along(rows), function(j) {
        xv <- if (is.null(x)) bugs_shape_value(sv[cells[[rows[j]]]], shape[[rows[j]]]) else x[[j]]
        lp <- do.call(d$lpdf, c(list(xv), a[[j]]))
        if (!is.null(b)) {
          bj <- b[[j]]
          if (!(btype == "C" && cobs[rows[j]]) && any(xv < bj[[1L]] | xv > bj[[2L]])) return(-Inf)
          if (btype == "T") lp <- lp - bugs_log_mass(cdf, bj[[1L]], bj[[2L]], a[[j]], discrete)
        }
        sum(lp)
      }, 0))
    }
    out$rng <- function(rows) { a <- args_f(rows); lapply(seq_along(rows), function(j) do.call(d$rng, c(list(1L), a[[j]]))) }
    out$typical <- function(rows) { a <- args_f(rows); lapply(seq_along(rows), function(j) do.call(d$typical, a[[j]])) }
  }
  for (nm in c("value", "lp", "rng", "typical")) environment(out[[nm]]) <- fenv
  out
}

bugs_shape_value <- function(v, shape) if (length(shape) >= 2L) array(v, shape) else v

# NIMBLE's named parameters: dnorm(mu, sd = s), dgamma(shape, scale = r).
# Positional arguments fill the BUGS parameters in order; a name either
# names a BUGS parameter or an alternative that is converted to one.
# Unknown names are errors.
.bugs_named <- list(
  dnorm = list(par = c("mean", "tau"), default = list(0, 1), alt = list(sd = c("tau", "sd2prec"), var = c("tau", "inv"))),
  dlnorm = list(par = c("meanlog", "taulog"), default = list(0, 1), alt = list(sdlog = c("taulog", "sd2prec"), varlog = c("taulog", "inv"))),
  dt = list(par = c("mu", "tau", "df"), default = list(0, 1, 1), alt = list(sigma = c("tau", "sd2prec"))),
  dlogis = list(par = c("location", "rate"), default = list(0, 1), alt = list(scale = c("rate", "inv"))),
  ddexp = list(par = c("location", "rate"), default = list(0, 1), alt = list(scale = c("rate", "inv"), var = c("rate", "var2rate"))),
  dgamma = list(par = c("shape", "rate"), default = list(NULL, 1), alt = list(scale = c("rate", "inv"))),
  dexp = list(par = "rate", default = list(1), alt = list(scale = c("rate", "inv"))),
  dbeta = list(par = c("shape1", "shape2"), default = list(NULL, NULL), alt = list()),
  dunif = list(par = c("min", "max"), default = list(0, 1), alt = list()),
  dweib = list(par = c("shape", "lambda"), default = list(NULL, 1), alt = list(scale = c("lambda", "scale2lambda"), rate = c("lambda", "rate2lambda"))),
  dbern = list(par = "prob", default = list(NULL), alt = list()),
  dbin = list(par = c("prob", "size"), default = list(NULL, NULL), alt = list()),
  dbinom = list(par = c("prob", "size"), default = list(NULL, NULL), alt = list()),
  dnegbin = list(par = c("prob", "size"), default = list(NULL, NULL), alt = list()),
  dpois = list(par = "lambda", default = list(NULL), alt = list()),
  dchisqr = list(par = "df", default = list(NULL), alt = list()),
  dchisq = list(par = "df", default = list(NULL), alt = list()),
  dcat = list(par = "prob", default = list(NULL), alt = list()),
  dmulti = list(par = c("prob", "size"), default = list(NULL, NULL), alt = list()),
  ddirch = list(par = "alpha", default = list(NULL), alt = list()),
  dmnorm = list(par = c("mean", "prec"), default = list(NULL, NULL), alt = list(cov = c("prec", "inverse"))),
  dwish = list(par = c("R", "df"), default = list(NULL, NULL), alt = list()))

bugs_named_args <- function(s) {
  spec <- .bugs_named[[s$dist]]
  if (is.null(spec)) stop_unsupported(sprintf("named parameters of %s", s$dist), s$loc)
  nm <- s$argnames
  a <- s$args
  num <- function(v) list(kind = "num", value = v, loc = NULL)
  bin <- function(op, x, y) list(kind = "binop", op = op, a = x, b = y, loc = NULL)
  conv <- function(how, e, others) switch(how,
    sd2prec = bin("/", num(1), bin("^", e, num(2))),
    inv = bin("/", num(1), e),
    var2rate = list(kind = "call", name = "sqrt", args = list(bin("/", num(2), e)), loc = NULL),
    scale2lambda = bin("^", e, list(kind = "unop", op = "-", a = others$shape, loc = NULL)),
    rate2lambda = bin("^", e, others$shape),
    inverse = list(kind = "call", name = "inverse", args = list(e), loc = NULL))
  out <- stats::setNames(vector("list", length(spec$par)), spec$par)
  unnamed <- which(!nzchar(nm))
  if (length(unnamed) > length(spec$par)) stop_semantic(sprintf("%s takes %d parameters", s$dist, length(spec$par)), s$loc)
  for (k in seq_along(unnamed)) out[[k]] <- a[[unnamed[k]]]
  named <- which(nzchar(nm))
  alts <- list()
  for (k in named) {
    if (nm[k] %in% spec$par) {
      if (!is.null(out[[nm[k]]])) stop_semantic(sprintf("parameter '%s' of %s is given twice", nm[k], s$dist), s$loc)
      out[[nm[k]]] <- a[[k]]
    } else if (nm[k] %in% names(spec$alt)) alts[[nm[k]]] <- a[[k]]
    else if (s$dist == "dgamma" && nm[k] %in% c("mean", "sd")) alts[[nm[k]]] <- a[[k]]
    else stop_semantic(sprintf("%s has no parameter '%s'", s$dist, nm[k]), s$loc)
  }
  if (s$dist == "dgamma" && xor(!is.null(alts$mean), !is.null(alts$sd)))
    stop_semantic("dgamma needs both mean and sd", s$loc)
  if (s$dist == "dgamma" && !is.null(alts$mean)) {
    v <- bin("^", alts$sd, num(2))
    out$shape <- bin("/", bin("^", alts$mean, num(2)), v)
    out$rate <- bin("/", alts$mean, v)
    alts$mean <- alts$sd <- NULL
  }
  for (an in names(alts)) {
    target <- spec$alt[[an]][1L]
    if (!is.null(out[[target]])) stop_semantic(sprintf("%s: give only one of '%s' and '%s'", s$dist, target, an), s$loc)
    out[[target]] <- conv(spec$alt[[an]][2L], alts[[an]], out)
  }
  for (k in seq_along(out)) if (is.null(out[[k]])) {
    if (is.null(spec$default[[k]])) stop_semantic(sprintf("parameter '%s' of %s is missing", spec$par[k], s$dist), s$loc)
    out[[k]] <- num(spec$default[[k]])
  }
  s$args <- unname(out)
  s$argnames <- NULL
  s
}

# evaluate logical nodes (ids, in order), grouping consecutive nodes of a
# statement that is evaluated elementwise
bugs_eval_nodes <- function(M, ids) {
  if (!length(ids)) return(invisible())
  st <- M$nodes$stmt[ids]
  run <- rle(st)
  pos <- 0L
  for (k in seq_along(run$lengths)) {
    sel <- ids[pos + seq_len(run$lengths[k])]
    pos <- pos + run$lengths[k]
    M$stmts[[run$values[k]]]$fn$write(M$nodes$row[sel])
  }
  invisible()
}
