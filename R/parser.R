# Recursive descent parser for Stan programs (Stan 2.40 grammar), with a
# precedence climbing expression parser.
#
# The parser produces an untyped abstract syntax tree made of plain R lists.
# Every node has a `kind` and a `loc` (line and column). The semantic
# analyzer later annotates expression nodes with their Stan type.

stan_block_names <- c("functions", "data", "transformed data", "parameters",
                      "transformed parameters", "model", "generated quantities")

basic_types <- c("int", "real", "complex", "vector", "row_vector", "matrix",
                 "simplex", "ordered", "positive_ordered", "unit_vector",
                 "sum_to_zero_vector", "cholesky_factor_corr", "cholesky_factor_cov",
                 "corr_matrix", "cov_matrix", "column_stochastic_matrix",
                 "row_stochastic_matrix", "complex_vector", "complex_row_vector",
                 "complex_matrix", "tuple", "sum_to_zero_matrix")

stan_parse <- function(src) {
  toks <- stan_lex(src)
  p <- new.env(parent = emptyenv())
  p$toks <- toks
  p$pos <- 1L
  prog <- parse_program(p)
  prog
}

# --- token helpers -----------------------------------------------------------

tk <- function(p, k = 0L) p$toks[[min(p$pos + k, length(p$toks))]]
tk_loc <- function(t) list(line = t$line, col = t$col)
tk_is <- function(p, value, type = NULL, k = 0L) {
  t <- tk(p, k)
  t$value == value && (is.null(type) || t$type == type) && t$type %in% c("op", "ident")
}
tk_next <- function(p) {
  t <- tk(p)
  p$pos <- p$pos + 1L
  t
}
expect <- function(p, value, what = NULL) {
  t <- tk(p)
  if (!(t$type %in% c("op", "ident") && t$value == value)) {
    found <- if (t$type == "eof") "end of input" else sprintf("'%s'", t$value)
    stop_syntax(sprintf("expected '%s'%s but found %s", value,
                        if (is.null(what)) "" else paste0(" ", what), found), tk_loc(t))
  }
  tk_next(p)
}
expect_ident <- function(p, what = "an identifier") {
  t <- tk(p)
  if (t$type != "ident") {
    found <- if (t$type == "eof") "end of input" else sprintf("'%s'", t$value)
    stop_syntax(sprintf("expected %s but found %s", what, found), tk_loc(t))
  }
  tk_next(p)
}

# --- program -----------------------------------------------------------------

parse_program <- function(p) {
  blocks <- list()
  order_seen <- 0L
  while (tk(p)$type != "eof") {
    t <- tk(p)
    if (t$type != "ident") stop_syntax(sprintf("expected a program block but found '%s'", t$value), tk_loc(t))
    name <- t$value
    if (name %in% c("transformed", "generated")) {
      tk_next(p)
      t2 <- expect_ident(p, "a block name")
      name <- paste(name, t2$value)
    } else {
      tk_next(p)
    }
    idx <- match(name, stan_block_names)
    if (is.na(idx)) stop_syntax(sprintf("unknown program block '%s'", name), tk_loc(t))
    if (idx <= order_seen) {
      if (!is.null(blocks[[name]])) stop_syntax(sprintf("block '%s' appears more than once", name), tk_loc(t))
      stop_syntax(sprintf("block '%s' is out of order; Stan blocks must appear in the order: %s",
                          name, paste(stan_block_names, collapse = ", ")), tk_loc(t))
    }
    order_seen <- idx
    expect(p, "{", sprintf("to open block '%s'", name))
    if (name == "functions") {
      body <- list()
      while (!tk_is(p, "}")) {
        if (tk(p)$type == "eof") stop_syntax("unterminated functions block", tk_loc(tk(p)))
        body[[length(body) + 1L]] <- parse_fundef(p)
      }
    } else {
      body <- parse_stmts_until_close(p, name)
    }
    expect(p, "}", sprintf("to close block '%s'", name))
    blocks[[name]] <- list(kind = "block", name = name, body = body, loc = tk_loc(t))
  }
  list(kind = "program", blocks = blocks)
}

parse_stmts_until_close <- function(p, block) {
  body <- list()
  while (!tk_is(p, "}")) {
    if (tk(p)$type == "eof") stop_syntax(sprintf("missing '}' to close block '%s'", block), tk_loc(tk(p)))
    st <- parse_statement(p, block)
    if (st$kind == "decls") body <- c(body, st$items) else body[[length(body) + 1L]] <- st
  }
  body
}

# --- declarations -------------------------------------------------------------

starts_decl <- function(p) {
  t <- tk(p)
  if (t$type != "ident") return(FALSE)
  if (t$value == "array") return(tk_is(p, "[", k = 1L))
  if (t$value %in% c("data")) return(FALSE)
  if (!(t$value %in% basic_types)) return(FALSE)
  # a type keyword followed by an identifier, '<', '[' or '('
  nx <- tk(p, 1L)
  nx$type == "ident" || (nx$type == "op" && nx$value %in% c("<", "[", "("))
}

parse_expr_list_in <- function(p, open, close) {
  expect(p, open)
  out <- list()
  if (!tk_is(p, close)) {
    repeat {
      out[[length(out) + 1L]] <- parse_expr(p)
      if (tk_is(p, ",")) {
        tk_next(p)
        next
      }
      break
    }
  }
  expect(p, close)
  out
}

# parse `<lower = e, upper = e>` or `<offset = e, multiplier = e>`
parse_constraint <- function(p) {
  out <- list()
  if (!tk_is(p, "<")) return(out)
  expect(p, "<")
  repeat {
    t <- expect_ident(p, "'lower', 'upper', 'offset' or 'multiplier'")
    if (!(t$value %in% c("lower", "upper", "offset", "multiplier")))
      stop_syntax(sprintf("unknown constraint '%s'", t$value), tk_loc(t))
    if (!is.null(out[[t$value]])) stop_syntax(sprintf("constraint '%s' given twice", t$value), tk_loc(t))
    expect(p, "=")
    out[[t$value]] <- parse_expr(p, no_gt = TRUE)
    if (tk_is(p, ",")) {
      tk_next(p)
      next
    }
    break
  }
  expect(p, ">", "to close the constraint")
  if (!is.null(out$offset) || !is.null(out$multiplier)) {
    if (!is.null(out$lower) || !is.null(out$upper))
      stop_syntax("offset/multiplier cannot be combined with lower/upper bounds", tk_loc(tk(p)))
  }
  out
}

parse_type <- function(p) {
  t0 <- tk(p)
  adims <- list()
  if (tk_is(p, "array")) {
    tk_next(p)
    adims <- parse_expr_list_in(p, "[", "]")
    if (length(adims) == 0L) stop_syntax("array declaration needs at least one dimension", tk_loc(t0))
  }
  t <- expect_ident(p, "a type")
  base <- t$value
  if (!(base %in% basic_types)) stop_syntax(sprintf("unknown type '%s'", base), tk_loc(t))
  if (base %in% c("complex", "complex_vector", "complex_row_vector", "complex_matrix"))
    stop_unsupported("complex numbers are not supported in this version", tk_loc(t))
  if (base == "tuple") stop_unsupported("tuple types are not supported in this version", tk_loc(t))
  cons <- list()
  if (base %in% c("int", "real", "vector", "row_vector", "matrix")) cons <- parse_constraint(p)
  sizes <- list()
  needs <- switch(base,
    int = 0L, real = 0L,
    vector = 1L, row_vector = 1L, simplex = 1L, ordered = 1L, positive_ordered = 1L,
    unit_vector = 1L, sum_to_zero_vector = 1L, cholesky_factor_corr = 1L, corr_matrix = 1L,
    cov_matrix = 1L, matrix = 2L, column_stochastic_matrix = 2L, row_stochastic_matrix = 2L,
    cholesky_factor_cov = -1L, sum_to_zero_matrix = 2L)
  if (needs != 0L) {
    if (!tk_is(p, "[")) stop_syntax(sprintf("type '%s' needs a size in brackets", base), tk_loc(tk(p)))
    sizes <- parse_expr_list_in(p, "[", "]")
    if (needs > 0L && length(sizes) != needs)
      stop_syntax(sprintf("type '%s' takes %d size%s", base, needs, if (needs > 1L) "s" else ""), tk_loc(t))
    if (needs < 0L && !(length(sizes) %in% 1:2))
      stop_syntax("cholesky_factor_cov takes one or two sizes", tk_loc(t))
  }
  if (tk_is(p, "[") && length(adims) == 0L)
    stop_syntax("old array syntax is not allowed; write 'array[...] type name' (removed in Stan 2.33)", tk_loc(tk(p)))
  list(base = base, adims = adims, sizes = sizes, cons = cons, loc = tk_loc(t0))
}

parse_decl <- function(p, block) {
  ty <- parse_type(p)
  items <- list()
  repeat {
    t <- expect_ident(p, "a variable name")
    if (t$value %in% stan_reserved_words)
      stop_syntax(sprintf("'%s' is a reserved word and cannot be used as a variable name", t$value), tk_loc(t))
    if (tk_is(p, "[")) stop_syntax("old array syntax is not allowed; write 'array[...] type name' (removed in Stan 2.33)", tk_loc(tk(p)))
    init <- NULL
    if (tk_is(p, "=")) {
      tk_next(p)
      init <- parse_expr(p)
    }
    items[[length(items) + 1L]] <- list(kind = "decl", name = t$value, type = ty, init = init,
                                         block = block, loc = tk_loc(t))
    if (tk_is(p, ",")) {
      tk_next(p)
      next
    }
    break
  }
  expect(p, ";", "after declaration")
  if (length(items) == 1L) return(items[[1L]])
  list(kind = "decls", items = items, loc = items[[1L]]$loc)
}

# --- functions ----------------------------------------------------------------

# unsized type for function arguments and return values
parse_unsized_type <- function(p) {
  t0 <- tk(p)
  nd <- 0L
  if (tk_is(p, "array")) {
    tk_next(p)
    expect(p, "[")
    nd <- 1L
    while (tk_is(p, ",")) {
      tk_next(p)
      nd <- nd + 1L
    }
    expect(p, "]")
  }
  t <- expect_ident(p, "a type")
  if (!(t$value %in% c("int", "real", "vector", "row_vector", "matrix", "void", "complex",
                       "complex_vector", "complex_row_vector", "complex_matrix", "tuple")))
    stop_syntax(sprintf("'%s' is not a valid function argument or return type", t$value), tk_loc(t))
  if (t$value %in% c("complex", "complex_vector", "complex_row_vector", "complex_matrix", "tuple"))
    stop_unsupported(sprintf("type '%s' is not supported in this version", t$value), tk_loc(t))
  if (tk_is(p, "[")) stop_syntax("function argument types are unsized; old array syntax is not allowed", tk_loc(tk(p)))
  list(base = t$value, nd = nd, loc = tk_loc(t0))
}

parse_fundef <- function(p) {
  rt <- parse_unsized_type(p)
  name <- expect_ident(p, "a function name")
  expect(p, "(")
  args <- list()
  if (!tk_is(p, ")")) {
    repeat {
      data_only <- FALSE
      if (tk_is(p, "data")) {
        tk_next(p)
        data_only <- TRUE
      }
      at <- parse_unsized_type(p)
      if (at$base == "void") stop_syntax("function arguments cannot be void", at$loc)
      an <- expect_ident(p, "an argument name")
      args[[length(args) + 1L]] <- list(name = an$value, type = at, data_only = data_only, loc = tk_loc(an))
      if (tk_is(p, ",")) {
        tk_next(p)
        next
      }
      break
    }
  }
  expect(p, ")")
  body <- NULL
  if (tk_is(p, ";")) {
    tk_next(p)
  } else {
    body <- parse_statement(p, "functions")
    if (body$kind != "block_stmt") stop_syntax("function body must be a block in braces", body$loc)
  }
  list(kind = "fundef", name = name$value, rtype = rt, args = args, body = body, loc = tk_loc(name))
}

# --- statements -----------------------------------------------------------------

assign_ops <- c("=", "+=", "-=", "*=", "/=", ".*=", "./=")

parse_statement <- function(p, block) {
  t <- tk(p)
  loc <- tk_loc(t)
  if (t$type == "op" && t$value == ";") {
    tk_next(p)
    return(list(kind = "empty", loc = loc))
  }
  if (t$type == "op" && t$value == "{") {
    tk_next(p)
    body <- parse_stmts_until_close(p, block)
    expect(p, "}")
    return(list(kind = "block_stmt", body = body, loc = loc))
  }
  if (t$type == "ident") {
    if (starts_decl(p)) return(parse_decl(p, block))
    switch(t$value,
      "if" = {
        tk_next(p)
        expect(p, "(")
        cond <- parse_expr(p)
        expect(p, ")")
        then <- parse_statement(p, block)
        els <- NULL
        if (tk_is(p, "else")) {
          tk_next(p)
          els <- parse_statement(p, block)
        }
        return(list(kind = "if", cond = cond, then = then, els = els, loc = loc))
      },
      "for" = {
        tk_next(p)
        expect(p, "(")
        v <- expect_ident(p, "a loop variable")
        expect(p, "in")
        e1 <- parse_expr(p)
        if (tk_is(p, ":")) {
          tk_next(p)
          e2 <- parse_expr(p)
          expect(p, ")")
          body <- parse_statement(p, block)
          return(list(kind = "for", var = v$value, lo = e1, hi = e2, body = body, loc = loc))
        }
        expect(p, ")")
        body <- parse_statement(p, block)
        return(list(kind = "foreach", var = v$value, iter = e1, body = body, loc = loc))
      },
      "while" = {
        tk_next(p)
        expect(p, "(")
        cond <- parse_expr(p)
        expect(p, ")")
        body <- parse_statement(p, block)
        return(list(kind = "while", cond = cond, body = body, loc = loc))
      },
      "break" = , "continue" = {
        tk_next(p)
        expect(p, ";")
        return(list(kind = t$value, loc = loc))
      },
      "return" = {
        tk_next(p)
        if (tk_is(p, ";")) {
          tk_next(p)
          return(list(kind = "return", value = NULL, loc = loc))
        }
        e <- parse_expr(p)
        expect(p, ";")
        return(list(kind = "return", value = e, loc = loc))
      },
      "print" = , "reject" = , "fatal_error" = {
        tk_next(p)
        expect(p, "(")
        items <- list()
        repeat {
          if (tk(p)$type == "string") {
            s <- tk_next(p)
            items[[length(items) + 1L]] <- list(kind = "string", value = s$value, loc = tk_loc(s))
          } else {
            items[[length(items) + 1L]] <- parse_expr(p)
          }
          if (tk_is(p, ",")) {
            tk_next(p)
            next
          }
          break
        }
        expect(p, ")")
        expect(p, ";")
        return(list(kind = t$value, items = items, loc = loc))
      },
      "profile" = {
        tk_next(p)
        expect(p, "(")
        s <- tk_next(p)
        if (s$type != "string") stop_syntax("profile() needs a string name", tk_loc(s))
        expect(p, ")")
        body <- parse_statement(p, block)
        if (body$kind != "block_stmt") stop_syntax("profile needs a block in braces", loc)
        return(list(kind = "block_stmt", body = body$body, loc = loc))
      },
      "target" = {
        if (tk_is(p, "+=", k = 1L)) {
          tk_next(p)
          tk_next(p)
          e <- parse_expr(p)
          expect(p, ";")
          return(list(kind = "target_add", value = e, loc = loc))
        }
      },
      "jacobian" = {
        if (tk_is(p, "+=", k = 1L)) {
          tk_next(p)
          tk_next(p)
          e <- parse_expr(p)
          expect(p, ";")
          return(list(kind = "jacobian_add", value = e, loc = loc))
        }
      },
      "increment_log_prob" = stop_syntax("increment_log_prob() was removed from Stan; use 'target +='", loc),
      NULL)
  }
  # expression based statements
  lhs <- parse_expr(p)
  nt <- tk(p)
  if (nt$type == "op" && nt$value %in% assign_ops) {
    tk_next(p)
    rhs <- parse_expr(p)
    expect(p, ";", "after assignment")
    if (!(lhs$kind %in% c("var", "index")))
      stop_syntax("the left side of an assignment must be a variable, optionally indexed", lhs$loc)
    return(list(kind = "assign", op = nt$value, lhs = lhs, rhs = rhs, loc = loc))
  }
  if (nt$type == "op" && nt$value == "~") {
    tk_next(p)
    dt <- expect_ident(p, "a distribution name")
    args <- parse_expr_list_in(p, "(", ")")
    trunc <- NULL
    if (tk(p)$type == "ident" && tk(p)$value == "T" && tk_is(p, "[", k = 1L)) {
      tk_next(p)
      expect(p, "[")
      lo <- NULL
      hi <- NULL
      if (!tk_is(p, ",")) lo <- parse_expr(p)
      expect(p, ",", "in truncation bounds")
      if (!tk_is(p, "]")) hi <- parse_expr(p)
      expect(p, "]")
      trunc <- list(lo = lo, hi = hi)
    }
    expect(p, ";", "after sampling statement")
    return(list(kind = "tilde", lhs = lhs, dist = dt$value, args = args, trunc = trunc, loc = loc))
  }
  if (lhs$kind == "call") {
    expect(p, ";", "after function call")
    return(list(kind = "call_stmt", call = lhs, loc = loc))
  }
  found <- if (nt$type == "eof") "end of input" else sprintf("'%s'", nt$value)
  stop_syntax(sprintf("expected an assignment, sampling statement or function call, found %s", found), tk_loc(nt))
}

# --- expressions ------------------------------------------------------------------

binary_prec <- c(
  "||" = 1L, "&&" = 2L, "==" = 3L, "!=" = 3L,
  "<" = 4L, "<=" = 4L, ">" = 4L, ">=" = 4L,
  "+" = 5L, "-" = 5L,
  # as in stanc3: .* and ./ share the level of * and /, and \ binds tighter
  "*" = 6L, "/" = 6L, "%" = 6L, "%/%" = 6L, ".*" = 6L, "./" = 6L,
  "\\" = 7L
)

parse_expr <- function(p, no_gt = FALSE) {
  # no_gt: a constraint expression, which may only use arithmetic operators
  cond <- parse_binary(p, if (no_gt) 5L else 1L, no_gt)
  if (!no_gt && tk_is(p, "?")) {
    q <- tk_next(p)
    a <- parse_expr(p, no_gt)
    expect(p, ":", "in conditional expression")
    b <- parse_expr(p, no_gt)
    return(list(kind = "cond", cond = cond, a = a, b = b, loc = tk_loc(q)))
  }
  cond
}

parse_binary <- function(p, min_prec, no_gt) {
  lhs <- parse_unary(p, no_gt)
  repeat {
    t <- tk(p)
    if (t$type != "op") break
    pr <- binary_prec[t$value]
    if (is.na(pr) || pr < min_prec) break
    if (no_gt && t$value %in% c(">", ">=")) break
    tk_next(p)
    rhs <- parse_binary(p, pr + 1L, no_gt)
    lhs <- list(kind = "binop", op = t$value, a = lhs, b = rhs, loc = tk_loc(t))
  }
  lhs
}

parse_unary <- function(p, no_gt) {
  t <- tk(p)
  if (t$type == "op" && t$value %in% c("-", "!", "+")) {
    tk_next(p)
    e <- parse_unary(p, no_gt)
    return(list(kind = "unop", op = t$value, a = e, loc = tk_loc(t)))
  }
  parse_power(p, no_gt)
}

parse_power <- function(p, no_gt) {
  base <- parse_postfix(p, no_gt)
  t <- tk(p)
  if (t$type == "op" && t$value %in% c("^", ".^")) {
    tk_next(p)
    ex <- parse_unary(p, no_gt)
    return(list(kind = "binop", op = t$value, a = base, b = ex, loc = tk_loc(t)))
  }
  base
}

parse_postfix <- function(p, no_gt) {
  e <- parse_primary(p, no_gt)
  repeat {
    t <- tk(p)
    if (t$type == "op" && t$value == "[") {
      e <- list(kind = "index", x = e, idx = parse_indices(p), loc = tk_loc(t))
      next
    }
    if (t$type == "op" && t$value == "'") {
      tk_next(p)
      e <- list(kind = "transpose", a = e, loc = tk_loc(t))
      next
    }
    break
  }
  e
}

parse_indices <- function(p) {
  expect(p, "[")
  out <- list()
  repeat {
    t <- tk(p)
    if (t$type == "op" && t$value %in% c(",", "]")) {
      out[[length(out) + 1L]] <- list(kind = "omni", loc = tk_loc(t))
    } else if (t$type == "op" && t$value == ":") {
      tk_next(p)
      if (tk_is(p, ",") || tk_is(p, "]")) {
        out[[length(out) + 1L]] <- list(kind = "omni", loc = tk_loc(t))
      } else {
        out[[length(out) + 1L]] <- list(kind = "upper", hi = parse_expr(p), loc = tk_loc(t))
      }
    } else {
      e <- parse_expr(p)
      if (tk_is(p, ":")) {
        tk_next(p)
        if (tk_is(p, ",") || tk_is(p, "]")) {
          out[[length(out) + 1L]] <- list(kind = "lower", lo = e, loc = e$loc)
        } else {
          out[[length(out) + 1L]] <- list(kind = "range", lo = e, hi = parse_expr(p), loc = e$loc)
        }
      } else {
        out[[length(out) + 1L]] <- list(kind = "single", e = e, loc = e$loc)
      }
    }
    if (tk_is(p, ",")) {
      tk_next(p)
      next
    }
    break
  }
  expect(p, "]")
  out
}

parse_primary <- function(p, no_gt) {
  t <- tk(p)
  loc <- tk_loc(t)
  if (t$type == "int") {
    tk_next(p)
    v <- suppressWarnings(as.numeric(t$value))
    if (is.na(v) || v > .Machine$integer.max)
      stop_syntax(sprintf("integer literal %s is too large; Stan integers are 32 bit", t$value), loc)
    return(list(kind = "int", value = as.integer(v), loc = loc))
  }
  if (t$type == "real") {
    tk_next(p)
    return(list(kind = "real", value = as.numeric(t$value), loc = loc))
  }
  if (t$type == "imag") stop_unsupported("imaginary literals (complex numbers) are not supported", loc)
  if (t$type == "string") stop_syntax("string literals are only allowed in print, reject and fatal_error", loc)
  if (t$type == "op" && t$value == "(") {
    tk_next(p)
    e <- parse_expr(p)
    if (tk_is(p, ",")) stop_unsupported("tuple expressions are not supported in this version", loc)
    expect(p, ")")
    return(list(kind = "paren", a = e, loc = loc))
  }
  if (t$type == "op" && t$value == "{") {
    items <- parse_expr_list_in(p, "{", "}")
    if (length(items) == 0L) stop_syntax("array expressions cannot be empty", loc)
    return(list(kind = "array_expr", items = items, loc = loc))
  }
  if (t$type == "op" && t$value == "[") {
    items <- parse_expr_list_in(p, "[", "]")
    return(list(kind = "rowvec_expr", items = items, loc = loc))
  }
  if (t$type == "ident") {
    if (t$value %in% c("true", "false")) stop_syntax("Stan has no boolean literals; use 1 and 0", loc)
    tk_next(p)
    if (tk_is(p, "(")) {
      tk_next(p)
      args <- list()
      cond_bar <- FALSE
      if (!tk_is(p, ")")) {
        repeat {
          args[[length(args) + 1L]] <- parse_expr(p)
          if (tk_is(p, "|")) {
            if (length(args) != 1L) stop_syntax("'|' may only follow the first argument", tk_loc(tk(p)))
            tk_next(p)
            cond_bar <- TRUE
            # nothing after the bar: std_normal_lcdf(x|)
            if (tk_is(p, ")")) break
            next
          }
          if (tk_is(p, ",")) {
            tk_next(p)
            next
          }
          break
        }
      }
      expect(p, ")", "to close the function call")
      return(list(kind = "call", name = t$value, args = args, cond_bar = cond_bar, loc = loc))
    }
    return(list(kind = "var", name = t$value, loc = loc))
  }
  found <- if (t$type == "eof") "end of input" else sprintf("'%s'", t$value)
  stop_syntax(sprintf("expected an expression but found %s", found), loc)
}
