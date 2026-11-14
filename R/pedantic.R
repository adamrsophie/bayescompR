# Pedantic mode and uninitialized variable warnings (stanc's --warn-pedantic
# and --warn-uninitialized). The analyses follow stanc3's: a dependency graph
# of every variable on the parameters and data, the factor graph of the
# model's density terms, and a definite assignment pass.

#' Pedantic warnings
#'
#' Runs the checks of stanc's pedantic mode (`--warn-pedantic`) on a Stan
#' program: distribution arguments that do not match their domain,
#' parameters without a prior or with several, unused parameters, control
#' flow that depends on parameters, parameters on the left of several
#' distribution statements, very large or small constants in distributions,
#' hard bounds, uniform, lkj_corr and gamma(epsilon, epsilon) priors,
#' nonlinear transforms on the left of `~`, and variables used before they
#' are assigned. `bayescomp_uninitialized()` runs only the last check
#' (`--warn-uninitialized`).
#'
#' @param x A `bayescomp_model`, Stan code, or a `.stan` file name.
#' @param include_paths Directories searched for `#include` files.
#' @return A data frame of warnings (`line`, `column`, `message`) sorted by
#'   position, with class `bayescomp_warnings`; it prints like stanc.
#' @examples
#' bayescomp_pedantic("data { int N; array[N] real x; } parameters { real sigma; }
#'                     model { real mu; x ~ normal(mu, sigma); }")
#' @export
bayescomp_pedantic <- function(x, include_paths = NULL) {
  ir <- ped_ir(x, include_paths)
  ped_result(c(ped_checks(ir), ped_uninit(ir)))
}

#' @rdname bayescomp_pedantic
#' @export
bayescomp_uninitialized <- function(x, include_paths = NULL) ped_result(ped_uninit(ped_ir(x, include_paths)))

#' @export
print.bayescomp_warnings <- function(x, ...) {
  if (!nrow(x)) cat("No warnings.\n")
  for (i in seq_len(nrow(x))) {
    cat(sprintf("Warning at line %d, column %d:\n", x$line[i], x$column[i]))
    cat(paste0("    ", strwrap(x$message[i], 74L)), sep = "\n")
  }
  invisible(x)
}

ped_ir <- function(x, include_paths) {
  if (inherits(x, "bayescomp_model")) return(x$ir)
  f <- if (is.character(x) && length(x) == 1L && grepl("\\.stan$", x) && file.exists(x)) x
  stan_check(stan_parse(read_stan_source(if (is.null(f)) x, f, include_paths)))
}

ped_result <- function(w) {
  df <- if (length(w)) data.frame(line = vapply(w, function(z) as.integer(z$loc$line %||% NA_integer_), 0L),
                                  column = vapply(w, function(z) as.integer(z$loc$col %||% NA_integer_), 0L),
                                  message = vapply(w, `[[`, "", "msg"), stringsAsFactors = FALSE)
        else data.frame(line = integer(), column = integer(), message = character())
  df <- unique(df[order(df$line, df$column), , drop = FALSE])
  rownames(df) <- NULL
  class(df) <- c("bayescomp_warnings", "data.frame")
  df
}

pw <- function(loc, msg) list(loc = loc, msg = msg)

# ---- expression helpers --------------------------------------------------------

# variable nodes of an expression, in source order
ped_vars <- function(e) {
  out <- list()
  walk <- function(n) {
    if (!is.list(n)) return()
    if (identical(n$kind, "var")) { out[[length(out) + 1L]] <<- n; return() }
    for (k in setdiff(names(n) %||% seq_along(n), c("sig", "type", "loc", "vtype"))) {
      v <- n[[k]]
      if (is.list(v) && !is.function(v)) walk(v)
    }
  }
  walk(e)
  out
}

# value of a constant expression, or NULL
ped_fold <- function(e) {
  if (!is.list(e)) return(NULL)
  switch(e$kind %||% "",
    int = , real = as.numeric(e$value),
    paren = ped_fold(e$x %||% e$a),
    unop = { a <- ped_fold(e$a); if (is.null(a)) NULL else if (e$op == "-") -a else if (e$op == "+") a else NULL },
    binop = {
      a <- ped_fold(e$a); b <- ped_fold(e$b)
      if (is.null(a) || is.null(b)) NULL
      else switch(e$op, "+" = a + b, "-" = a - b, "*" = a * b, "/" = if (b != 0) a / b, NULL)
    },
    NULL)
}

ped_num <- function(v) format(v, digits = 15)

# ---- the checks ----------------------------------------------------------------

ped_checks <- function(ir) {
  w <- list()
  pars <- top_decls(ir, "parameters")
  pinfo <- stats::setNames(lapply(pars, function(d) {
    cons <- d$type$cons %||% list()
    list(lower = ped_fold(cons$lower), upper = ped_fold(cons$upper), has_lower = !is.null(cons$lower),
         has_upper = !is.null(cons$upper), base = d$type$base, loc = d$loc)
  }), vapply(pars, `[[`, "", "name"))

  # hard and nonsensical bounds
  for (nm in names(pinfo)) {
    p <- pinfo[[nm]]
    if (!is.null(p$lower) && !is.null(p$upper) && p$lower >= p$upper)
      w[[length(w) + 1L]] <- pw(p$loc, sprintf("Parameter %s has constraints that don't make sense. The lower bound should be strictly less than the upper bound.", nm))
    else if (p$has_lower && p$has_upper && !(identical(p$lower, 0) && identical(p$upper, 1)) && !(identical(p$lower, -1) && identical(p$upper, 1)))
      w[[length(w) + 1L]] <- pw(p$loc, sprintf("Your Stan program has a parameter %s with a lower and upper bound in its declaration. These hard constraints are not recommended, for two reasons: (a) Except when there are logical or physical constraints, it is very unusual for you to be sure that a parameter will fall inside a specified range, and (b) The infinite gradient induced by a hard constraint can cause difficulties for Stan's sampling algorithm. As a consequence, we recommend soft constraints rather than hard constraints; for example, instead of constraining an elasticity parameter to fall between 0, and 1, leave it unconstrained and give it a normal(0.5,0.5) prior distribution.", nm))
  }

  # statements of the transformed parameters and model blocks
  tp_body <- ir$blocks[["transformed parameters"]]$body %||% list()
  mb <- ir$blocks$model$body %||% list()
  stmts <- c(tp_body, mb)

  # dependency of every variable on parameters and data (two passes reach the
  # fixed point of simple loops), the density factors, and control flow
  dep <- new.env(parent = emptyenv())
  roots <- function(e) {
    r <- character()
    for (v in ped_vars(e)) {
      o <- v$origin %||% ""
      if (o == "parameters") r <- c(r, paste0("p:", v$name))
      else if (o %in% c("data", "transformed data")) {
        if (!(identical(v$type$base, "int") && identical(as.integer(v$type$nd %||% 0L), 0L))) r <- c(r, paste0("d:", v$name))
      } else r <- c(r, dep[[v$name]])
    }
    unique(r)
  }
  factors <- list(); cflow <- list()
  add_dep <- function(nm, r) assign(nm, unique(c(dep[[nm]], r)), envir = dep)
  walk <- function(ss, ctrl, record) {
    for (s in ss) {
      k <- s$kind %||% ""
      switch(k,
        decl = if (!is.null(s$init)) add_dep(s$name, c(roots(s$init), ctrl)) else if (length(ctrl)) add_dep(s$name, ctrl),
        assign = {
          b <- s$lhs; while (b$kind %in% c("index", "proj")) b <- b$x
          add_dep(b$name, c(roots(s$rhs), roots(s$lhs), ctrl))
        },
        unpack = for (nm in ped_unpack_names(s$lhs)) add_dep(nm, c(roots(s$rhs), ctrl)),
        "if" = {
          r <- roots(s$cond)
          if (record && any(startsWith(r, "p:"))) cflow[[length(cflow) + 1L]] <<- list(loc = s$loc, r = r)
          walk(list(s$then), c(ctrl, r), record); if (!is.null(s$els)) walk(list(s$els), c(ctrl, r), record)
        },
        "while" = {
          r <- roots(s$cond)
          if (record && any(startsWith(r, "p:"))) cflow[[length(cflow) + 1L]] <<- list(loc = s$loc, r = r)
          walk(list(s$body), c(ctrl, r), record)
        },
        "for" = {
          r <- unique(c(roots(s$lo), roots(s$hi)))
          if (record && any(startsWith(r, "p:"))) cflow[[length(cflow) + 1L]] <<- list(loc = s$loc, r = r)
          add_dep(s$var, c(r, ctrl)); walk(list(s$body), c(ctrl, r), record)
        },
        foreach = { r <- roots(s$iter); add_dep(s$var, c(r, ctrl)); walk(list(s$body), c(ctrl, r), record) },
        block_stmt = walk(s$body, ctrl, record),
        tilde = if (record) factors[[length(factors) + 1L]] <<- list(r = unique(c(roots(s[c("lhs", "args", "trunc")]), ctrl))),
        target_add = , jacobian_add = if (record) factors[[length(factors) + 1L]] <<- list(r = unique(c(roots(s$value), ctrl))),
        call_stmt = if (record && isTRUE(s$call$sig$lp %||% s$sig$lp)) factors[[length(factors) + 1L]] <<- list(r = unique(c(roots(s), ctrl))),
        NULL)
    }
  }
  walk(stmts, character(), FALSE)
  walk(stmts, character(), FALSE)
  # control flow is reported in the transformed parameters and model blocks
  walk(stmts, character(), TRUE)

  for (c in cflow) {
    p <- sub("^p:", "", c$r[startsWith(c$r, "p:")])
    w[[length(w) + 1L]] <- pw(c$loc, sprintf("A control flow statement depends on parameter(s): %s.", paste(p, collapse = ", ")))
  }

  # unused parameters and the number of priors (factor graph)
  used <- unique(unlist(lapply(factors, `[[`, "r")))
  for (nm in names(pinfo)) {
    pn <- paste0("p:", nm)
    if (!(pn %in% used)) {
      w[[length(w) + 1L]] <- pw(pinfo[[nm]]$loc, sprintf("The parameter %s was declared but was not used in the density calculation.", nm))
      next
    }
    n_prior <- 0L
    for (f in factors) {
      if (!(pn %in% f$r)) next
      # a prior factor reaches no data except through the parameter itself
      seen_f <- vapply(factors, function(g) identical(g, f), TRUE)
      nodes <- setdiff(f$r, pn)
      repeat {
        grow <- which(!seen_f & vapply(factors, function(g) any(g$r %in% nodes), TRUE))
        if (!length(grow)) break
        seen_f[grow] <- TRUE
        nodes <- unique(c(nodes, setdiff(unlist(lapply(factors[grow], `[[`, "r")), pn)))
      }
      if (!any(startsWith(nodes, "d:"))) n_prior <- n_prior + 1L
    }
    if (n_prior == 0L)
      w[[length(w) + 1L]] <- pw(pinfo[[nm]]$loc, sprintf("The parameter %s has no priors. This means either no prior is provided, or the prior(s) depend on data variables. In the later case, this may be a false positive.", nm))
    else if (n_prior > 1L)
      w[[length(w) + 1L]] <- pw(pinfo[[nm]]$loc, sprintf("The parameter %s has %d priors.", nm, n_prior))
  }

  # statement level checks on distribution statements
  lhs_seen <- list()
  sat <- function(p, cons) {
    if (is.null(p)) return(TRUE)
    lo <- p$lower; hi <- p$upper
    if (p$base %in% c("simplex", "stochastic_row_matrix", "stochastic_column_matrix", "column_stochastic_matrix", "row_stochastic_matrix"))
      return(cons != "strictly positive" || TRUE)
    if (p$base %in% c("positive_ordered", "cov_matrix", "cholesky_factor_cov")) return(cons %in% c("strictly positive", "non-negative"))
    switch(cons,
      "strictly positive" = , "non-negative" = !is.null(lo) && lo >= 0,
      "[0,1]" = , "(0,1)" = !is.null(lo) && lo >= 0 && !is.null(hi) && hi <= 1,
      TRUE)
  }
  in_dom <- function(v, cons) switch(cons, "strictly positive" = v > 0, "non-negative" = v >= 0,
                                     "[0,1]" = v >= 0 && v <= 1, "(0,1)" = v > 0 && v < 1, TRUE)
  lin <- function(e) {
    # TRUE when the expression is a linear function of its variables
    k <- e$kind %||% ""
    if (k %in% c("var", "index", "int", "real", "paren")) return(if (k == "paren") lin(e$x %||% e$a) else TRUE)
    if (k == "unop") return(lin(e$a))
    if (k == "binop" && e$op %in% c("+", "-")) return(lin(e$a) && lin(e$b))
    if (k == "binop" && e$op %in% c("*", ".*")) {
      ca <- !ped_param_or_local(e$a); cb <- !ped_param_or_local(e$b)
      return((ca && lin(e$b)) || (cb && lin(e$a)))
    }
    if (k == "binop" && e$op %in% c("/", "./")) return(!ped_param_or_local(e$b) && lin(e$a))
    !ped_param_or_local(e)
  }
  st_walk <- function(ss) {
    for (s in ss) {
      k <- s$kind %||% ""
      if (k %in% c("if")) { st_walk(list(s$then)); if (!is.null(s$els)) st_walk(list(s$els)) }
      else if (k %in% c("while", "for", "foreach")) st_walk(list(s$body))
      else if (k == "block_stmt") st_walk(s$body)
      else if (k == "target_add" && identical(s$value$kind, "call") && grepl("_(lpdf|lpmf|lupdf|lupmf)$", s$value$name %||% "")) {
        a1 <- s$value$args[[1L]]
        if (identical(a1$kind, "var")) lhs_seen[[length(lhs_seen) + 1L]] <<- list(name = a1$name, loc = s$loc)
        if (!lin(a1)) w[[length(w) + 1L]] <<- pw(s$loc, ped_jacobian_msg)
      } else if (k == "tilde") {
        dist <- s$dist
        user <- isTRUE(s$sig$user)
        if (identical(s$lhs$kind, "var")) lhs_seen[[length(lhs_seen) + 1L]] <<- list(name = s$lhs$name, loc = s$loc)
        if (!lin(s$lhs)) w[[length(w) + 1L]] <<- pw(s$loc, ped_jacobian_msg)
        if (user) next
        lhs_p <- if (identical(s$lhs$kind, "var") && identical(s$lhs$origin, "parameters")) pinfo[[s$lhs$name]]
        # variate support
        sup <- ped_variate_support[dist]
        if (!is.na(sup)) {
          lv <- ped_fold(s$lhs)
          if (!is.null(lhs_p) && !sat(lhs_p, sup))
            w[[length(w) + 1L]] <<- pw(s$lhs$loc, sprintf("Parameter %s is given a %s distribution, which has %s support, but %s was not constrained to be %s.", s$lhs$name, dist, sup, s$lhs$name, sup))
          else if (!is.null(lv) && !in_dom(lv, sup))
            w[[length(w) + 1L]] <<- pw(s$lhs$loc, sprintf("Value %s is given a %s distribution, which has %s support, but %s is not %s.", ped_num(lv), dist, sup, ped_num(lv), sup))
        }
        # arguments
        for (t in ped_arg_table) {
          if (t[[1L]] != dist || t[[2L]] > length(s$args)) next
          a <- s$args[[t[[2L]]]]
          av <- ped_fold(a)
          ap <- if (identical(a$kind, "var") && identical(a$origin, "parameters")) pinfo[[a$name]]
          if (!is.null(ap) && !sat(ap, t[[4L]]))
            w[[length(w) + 1L]] <<- pw(a$loc, sprintf("A %s distribution is given parameter %s as %s (argument %d), but %s was not constrained to be %s.", dist, a$name, t[[3L]], as.integer(t[[2L]]), a$name, t[[4L]]))
          else if (!is.null(av) && !in_dom(av, t[[4L]]))
            w[[length(w) + 1L]] <<- pw(a$loc, sprintf("A %s distribution is given value %s as %s (argument %d), but %s is not %s.", dist, ped_num(av), t[[3L]], as.integer(t[[2L]]), t[[3L]], t[[4L]]))
        }
        # constants that suggest a parameter is not on unit scale
        for (a in s$args) {
          av <- ped_fold(a)
          if (!is.null(av) && av != 0 && (abs(av) < 0.1 || abs(av) > 10))
            w[[length(w) + 1L]] <<- pw(a$loc, sprintf("Argument %s suggests there may be parameters that are not unit scale; consider rescaling with a multiplier, see: https://mc-stan.org/docs/stan-users-guide/efficiency-tuning.html#standardizing-predictors", ped_num(av)))
        }
        if (dist %in% c("gamma", "inv_gamma") && length(s$args) == 2L) {
          a1 <- ped_fold(s$args[[1L]]); a2 <- ped_fold(s$args[[2L]])
          if (!is.null(a1) && !is.null(a2) && a1 == a2 && a1 < 1)
            w[[length(w) + 1L]] <<- pw(s$args[[1L]]$loc, "There is a gamma or inverse-gamma distribution with parameters that are equal to each other and set to values less than 1. This is mathematically acceptable and can make sense in some problems, but typically we see this model used as an attempt to assign a noninformative prior distribution. In fact, priors such as inverse-gamma(.001,.001) can be very strong, as explained by Gelman (2006). Instead we recommend something like a normal(0,1) or student_t(4,0,1), with parameter constrained to be positive.")
        }
        if (dist == "lkj_corr")
          w[[length(w) + 1L]] <<- pw(s$loc, "It is suggested to reparameterize your model to replace lkj_corr with lkj_corr_cholesky, the Cholesky factor variant. lkj_corr tends to run slower, consume more memory, and has higher risk of numerical errors.")
        if (dist == "uniform" && !is.null(lhs_p) && length(s$args) >= 1L) {
          lo <- ped_fold(s$args[[1L]])
          if (is.null(lo) || is.null(lhs_p$lower) || lhs_p$lower != lo)
            w[[length(w) + 1L]] <<- pw(s$loc, sprintf("Parameter %s is given a uniform distribution. The uniform distribution is not recommended, for two reasons: (a) Except when there are logical or physical constraints, it is very unusual for you to be sure that a parameter will fall inside a specified range, and (b) The infinite gradient induced by a uniform density can cause difficulties for Stan's sampling algorithm. As a consequence, we recommend soft constraints rather than hard constraints; for example, instead of giving an elasticity parameter a uniform(0,1) distribution, try normal(0.5,0.5).", s$lhs$name))
        }
      }
    }
  }
  st_walk(mb)
  if (length(lhs_seen)) {
    nms <- vapply(lhs_seen, `[[`, "", "name")
    for (nm in unique(nms[duplicated(nms)])) {
      first <- lhs_seen[[match(nm, nms)]]
      w[[length(w) + 1L]] <- pw(first$loc, sprintf("The parameter %s is on the left-hand side of more than one tilde statement.", nm))
    }
  }
  w
}

# variables assigned by an unpacking statement (nested tuples included)
ped_unpack_names <- function(lhs) {
  if (identical(lhs$kind, "tuple_expr")) return(unlist(lapply(lhs$items, ped_unpack_names)))
  while (lhs$kind %in% c("index", "proj")) lhs <- lhs$x
  lhs$name
}

ped_param_or_local <- function(e) any(vapply(ped_vars(e), function(v) (v$origin %||% "") %in% c("parameters", "transformed parameters", "local"), TRUE))

ped_jacobian_msg <- "Left-hand side of distribution statement (~) may contain a non-linear transform of a parameter or local variable. If it does, you need to include a target += statement with the log absolute determinant of the Jacobian of the transform. You could also consider defining a transformed parameter and using jacobian += in the transformed parameters block."

# ---- definite assignment -------------------------------------------------------

ped_uninit <- function(ir) {
  w <- list()
  warned <- character()
  use <- function(e, st) {
    for (v in ped_vars(e)) {
      if (v$name %in% st$declared && !(v$name %in% st$assigned) && !(v$name %in% warned)) {
        warned <<- c(warned, v$name)
        w[[length(w) + 1L]] <<- pw(v$loc, sprintf("The variable %s may not have been assigned a value before its first use.", v$name))
      }
    }
  }
  run <- function(ss, st) {
    for (s in ss) {
      k <- s$kind %||% ""
      if (k == "decl") {
        if (!is.null(s$init)) { use(s$init, st); st$assigned <- c(st$assigned, s$name) }
        else { st$declared <- c(st$declared, s$name); st$assigned <- setdiff(st$assigned, s$name) }
      } else if (k == "assign") {
        use(s$rhs, st)
        b <- s$lhs; idx <- list()
        while (b$kind %in% c("index", "proj")) { idx <- c(idx, list(b[setdiff(names(b), c("x", "kind", "loc", "type"))])); b <- b$x }
        use(idx, st)
        if (!identical(s$op, "=")) use(b, st)
        st$assigned <- c(st$assigned, b$name)
      } else if (k == "unpack") {
        use(s$rhs, st)
        st$assigned <- c(st$assigned, ped_unpack_names(s$lhs))
      } else if (k == "if") {
        use(s$cond, st)
        a <- run(list(s$then), st)
        b <- if (!is.null(s$els)) run(list(s$els), st) else st
        st$assigned <- intersect(a$assigned, b$assigned)
      } else if (k %in% c("while")) { use(s$cond, st); run(list(s$body), st) }
      else if (k == "for") { use(list(s$lo, s$hi), st); run(list(s$body), st) }
      else if (k == "foreach") { use(s$iter, st); run(list(s$body), st) }
      else if (k == "block_stmt") st <- run(s$body, st)
      else use(s[setdiff(names(s), c("kind", "loc", "sig"))], st)
    }
    st
  }
  st <- list(declared = character(), assigned = character())
  for (fd in ir$blocks$functions$body %||% list()) if (!is.null(fd$body)) run(list(fd$body), st)
  for (b in c("transformed data", "transformed parameters", "model", "generated quantities")) {
    body <- ir$blocks[[b]]$body %||% list()
    st <- run(body, st)
    # top level variables of a block are checked when the block ends
    tops <- vapply(Filter(function(s) identical(s$kind, "decl"), body), `[[`, "", "name")
    st$assigned <- c(st$assigned, tops)
  }
  w
}

# argument domains of the built-in distributions, as stanc checks them
ped_arg_table <- list(
    list("bernoulli", 1, "chance of success", "[0,1]"),
    list("beta_binomial", 2, "a prior success count", "strictly positive"),
    list("beta_binomial", 3, "a prior failure count", "strictly positive"),
    list("beta", 1, "a count parameter", "strictly positive"),
    list("beta", 2, "a count parameter", "strictly positive"),
    list("beta_proportion", 1, "a unit mean parameter", "(0,1)"),
    list("beta_proportion", 2, "a precision parameter", "strictly positive"),
    list("binomial", 2, "chance of success", "[0,1]"),
    list("cauchy", 2, "a scale parameter", "strictly positive"),
    list("chi_square", 1, "degrees of freedom", "strictly positive"),
    list("double_exponential", 2, "a scale parameter", "strictly positive"),
    list("exp_mod_normal", 2, "a scale parameter", "strictly positive"),
    list("exp_mod_normal", 3, "a shape parameter", "strictly positive"),
    list("exponential", 1, "a scale parameter", "strictly positive"),
    list("frechet", 1, "a shape parameter", "strictly positive"),
    list("frechet", 2, "a scale parameter", "strictly positive"),
    list("gamma", 1, "a shape parameter", "strictly positive"),
    list("gamma", 2, "an inverse scale parameter", "strictly positive"),
    list("gumbel", 2, "a scale parameter", "strictly positive"),
    list("inv_chi_square", 1, "degrees of freedom", "strictly positive"),
    list("inv_gamma", 1, "a shape parameter", "strictly positive"),
    list("inv_gamma", 2, "a scale parameter", "strictly positive"),
    list("logistic", 2, "a scale parameter", "strictly positive"),
    list("lognormal", 2, "a scale parameter", "strictly positive"),
    list("neg_binomial_2_log", 2, "an inverse overdispersion control parameter", "strictly positive"),
    list("neg_binomial_2", 1, "a shape parameter", "strictly positive"),
    list("neg_binomial_2", 2, "a precision parameter", "strictly positive"),
    list("neg_binomial", 1, "a shape parameter", "strictly positive"),
    list("neg_binomial", 2, "an inverse scale parameter", "strictly positive"),
    list("normal", 2, "a scale parameter", "strictly positive"),
    list("pareto", 1, "a positive minimum parameter", "strictly positive"),
    list("pareto", 2, "a shape parameter", "strictly positive"),
    list("pareto_type_2", 2, "a scale parameter", "strictly positive"),
    list("pareto_type_2", 3, "a shape parameter", "strictly positive"),
    list("poisson", 1, "a rate parameter", "strictly positive"),
    list("rayleigh", 1, "a scale parameter", "strictly positive"),
    list("scaled_inv_chi_square", 1, "degrees of freedom", "strictly positive"),
    list("scaled_inv_chi_square", 2, "a scale parameter", "strictly positive"),
    list("skew_normal", 2, "a scale parameter", "strictly positive"),
    list("student_t", 1, "degrees of freedom", "strictly positive"),
    list("student_t", 3, "a scale parameter", "strictly positive"),
    list("von_mises", 2, "a scale parameter", "strictly positive"),
    list("weibull", 1, "a shape parameter", "strictly positive"),
    list("weibull", 2, "a scale parameter", "strictly positive"),
    list("wiener", 1, "a boundary separation parameter", "strictly positive"),
    list("wiener", 2, "a non-decision time parameter", "strictly positive"),
    list("wiener", 3, "an a-priori bias parameter", "[0,1]")
)

ped_variate_support <- c(
  beta = "(0,1)",
  beta_proportion = "(0,1)",
  chi_square = "strictly positive",
  exponential = "strictly positive",
  frechet = "strictly positive",
  gamma = "strictly positive",
  inv_chi_square = "strictly positive",
  inv_gamma = "strictly positive",
  lognormal = "strictly positive",
  pareto = "strictly positive",
  rayleigh = "non-negative",
  scaled_inv_chi_square = "strictly positive",
  weibull = "non-negative",
  wiener = "strictly positive"
)
