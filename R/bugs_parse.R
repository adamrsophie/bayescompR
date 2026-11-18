# Parser for the BUGS language, in the dialects of JAGS and OpenBUGS.
#
# The language is described in the OpenBUGS manual ("The BUGS Model
# Specification Language") and the JAGS user manual. A program is
#
#   [var declarations;]  [data { statements }]  model { statements }
#
# where a statement is a stochastic relation `lhs ~ dist(args)`, optionally
# followed by a bound `T(lo, hi)` (truncation), `I(lo, hi)` (OpenBUGS: the
# same as T for a node without observed parents) or `C(lo, hi)` (OpenBUGS
# censoring), a logical relation `lhs <- expr` or `link(lhs) <- expr` with
# link in log, logit, cloglog and probit, or a loop `for (i in a:b) { ... }`.
# Expressions follow the precedence of R, from which the language comes.
#
# The result is a list of statements; every node carries `loc` (line and
# column) for error messages.

bugs_parse <- function(code, file = NULL) {
  if (length(code) > 1L) code <- paste(code, collapse = "\n")
  toks <- bugs_lex(code, file)
  st <- new.env(parent = emptyenv())
  st$t <- toks; st$i <- 1L; st$file <- file
  prog <- list(var = list(), data = list(), model = NULL)
  while (!bugs_at(st, "eof")) {
    tk <- bugs_peek(st)
    if (bugs_at(st, "name", "var")) { bugs_next(st); prog$var <- c(prog$var, bugs_var_decls(st)); next }
    if (bugs_at(st, "name", "data")) {
      bugs_next(st)
      prog$data <- c(prog$data, bugs_block(st))
      next
    }
    if (bugs_at(st, "name", "model")) {
      bugs_next(st)
      if (!is.null(prog$model)) bugs_err(st, "the program has two model blocks")
      prog$model <- bugs_block(st)
      next
    }
    bugs_err(st, sprintf("expected 'model', 'data' or 'var' but found '%s'", tk$text))
  }
  if (is.null(prog$model)) bugs_err(st, "the program has no model block")
  structure(prog, class = "bugs_program")
}

# ---- lexer ------------------------------------------------------------------------------

bugs_lex <- function(code, file = NULL) {
  chars <- strsplit(code, "")[[1L]]
  n <- length(chars)
  out <- vector("list", 256L); k <- 0L
  i <- 1L; line <- 1L; col <- 1L
  push <- function(type, text, l, c) {
    k <<- k + 1L
    if (k > length(out)) length(out) <<- 2L * length(out)
    out[[k]] <<- list(type = type, text = text, loc = list(line = l, col = c, file = file))
  }
  adv <- function(m = 1L) {
    for (j in seq_len(m)) {
      if (chars[i] == "\n") { line <<- line + 1L; col <<- 1L } else col <<- col + 1L
      i <<- i + 1L
    }
  }
  is_digit <- function(ch) ch >= "0" && ch <= "9"
  is_alpha <- function(ch) grepl("^[A-Za-z]$", ch)
  while (i <= n) {
    ch <- chars[i]
    if (ch %in% c(" ", "\t", "\r", "\n", "\f")) { adv(); next }
    if (ch == "#") { while (i <= n && chars[i] != "\n") adv(); next }
    l0 <- line; c0 <- col
    nx <- if (i < n) chars[i + 1L] else ""
    if (is_alpha(ch) || (ch == "." && !is_digit(nx)) || ch == "_") {
      j <- i
      while (j <= n && grepl("^[A-Za-z0-9._]$", chars[j])) j <- j + 1L
      txt <- paste(chars[i:(j - 1L)], collapse = "")
      adv(j - i)
      push("name", txt, l0, c0)
      next
    }
    if (is_digit(ch) || (ch == "." && is_digit(nx))) {
      j <- i
      while (j <= n && (is_digit(chars[j]) || chars[j] == ".")) j <- j + 1L
      if (j <= n && chars[j] %in% c("e", "E")) {
        j2 <- j + 1L
        if (j2 <= n && chars[j2] %in% c("+", "-")) j2 <- j2 + 1L
        if (j2 <= n && is_digit(chars[j2])) {
          j <- j2
          while (j <= n && is_digit(chars[j])) j <- j + 1L
        }
      }
      txt <- paste(chars[i:(j - 1L)], collapse = "")
      if (is.na(suppressWarnings(as.numeric(txt)))) stop_syntax(sprintf("malformed number '%s'", txt), list(line = l0, col = c0, file = file))
      adv(j - i)
      push("num", txt, l0, c0)
      next
    }
    if (ch == "%") {
      j <- i + 1L
      while (j <= n && chars[j] != "%" && chars[j] != "\n") j <- j + 1L
      if (j > n || chars[j] != "%") stop_syntax("unterminated %operator%", list(line = l0, col = c0, file = file))
      txt <- paste(chars[i:j], collapse = "")
      adv(j - i + 1L)
      push("op", txt, l0, c0)
      next
    }
    two <- paste0(ch, nx)
    if (two %in% c("<-", "<=", ">=", "==", "!=", "&&", "||")) { adv(2L); push("op", two, l0, c0); next }
    if (ch %in% c("~", "(", ")", "[", "]", "{", "}", ",", ":", ";", "+", "-", "*", "/", "^", "<", ">", "!", "=", "&", "|")) {
      adv(); push("op", ch, l0, c0); next
    }
    stop_syntax(sprintf("unexpected character '%s'", ch), list(line = l0, col = c0, file = file))
  }
  push("eof", "end of input", line, col)
  out[seq_len(k)]
}

# ---- token stream -----------------------------------------------------------------------

bugs_peek <- function(st, ahead = 0L) st$t[[min(st$i + ahead, length(st$t))]]
bugs_next <- function(st) { tk <- st$t[[st$i]]; if (st$i < length(st$t)) st$i <- st$i + 1L; tk }
bugs_at <- function(st, type, text = NULL, ahead = 0L) {
  tk <- bugs_peek(st, ahead)
  identical(tk$type, type) && (is.null(text) || tk$text %in% text)
}
bugs_accept <- function(st, text) if (bugs_at(st, "op", text)) { bugs_next(st); TRUE } else FALSE
bugs_expect <- function(st, text) {
  if (!bugs_at(st, "op", text)) bugs_err(st, sprintf("expected '%s' but found '%s'", text, bugs_peek(st)$text))
  bugs_next(st)
}
bugs_err <- function(st, msg) stop_syntax(msg, bugs_peek(st)$loc)

# ---- statements -------------------------------------------------------------------------

# JAGS: var x[N], y[N, M], z;
bugs_var_decls <- function(st) {
  out <- list()
  repeat {
    if (!bugs_at(st, "name")) bugs_err(st, "expected a variable name in the var declarations")
    nm <- bugs_next(st)
    dims <- list()
    if (bugs_accept(st, "[")) {
      dims <- bugs_index_list(st)
      bugs_expect(st, "]")
    }
    out[[length(out) + 1L]] <- list(name = nm$text, dims = dims, loc = nm$loc)
    if (!bugs_accept(st, ",")) break
  }
  bugs_accept(st, ";")
  out
}

bugs_block <- function(st) {
  bugs_expect(st, "{")
  stmts <- list()
  while (!bugs_at(st, "op", "}")) {
    if (bugs_at(st, "eof")) bugs_err(st, "missing '}'")
    if (bugs_accept(st, ";")) next
    stmts[[length(stmts) + 1L]] <- bugs_statement(st)
  }
  bugs_next(st)
  stmts
}

bugs_links <- c("log", "logit", "cloglog", "probit")

bugs_statement <- function(st) {
  tk <- bugs_peek(st)
  if (bugs_at(st, "name", "for") && bugs_at(st, "op", "(", 1L)) {
    bugs_next(st); bugs_expect(st, "(")
    if (!bugs_at(st, "name")) bugs_err(st, "expected the loop variable")
    v <- bugs_next(st)$text
    if (!bugs_at(st, "name", "in")) bugs_err(st, "expected 'in'")
    bugs_next(st)
    lo <- bugs_expr(st)
    bugs_expect(st, ":")
    hi <- bugs_expr(st)
    bugs_expect(st, ")")
    body <- if (bugs_at(st, "op", "{")) bugs_block(st) else list(bugs_statement(st))
    return(list(kind = "for", var = v, lo = lo, hi = hi, body = body, loc = tk$loc))
  }
  if (!bugs_at(st, "name")) bugs_err(st, sprintf("expected a statement but found '%s'", tk$text))
  # OpenBUGS ordinary differential equations: D(x[1:n], t) <- expr
  if (tk$text == "D" && bugs_at(st, "op", "(", 1L)) {
    bugs_next(st); bugs_next(st)
    lhs <- bugs_lhs(st)
    bugs_expect(st, ",")
    tv <- bugs_lhs(st)
    bugs_expect(st, ")")
    if (!bugs_accept(st, "<-")) bugs_err(st, "expected '<-' after D(...)")
    return(list(kind = "derivative", lhs = lhs, time = tv, expr = bugs_expr(st), loc = tk$loc))
  }
  # link(lhs) <- expr
  link <- NULL
  if (tk$text %in% bugs_links && bugs_at(st, "op", "(", 1L)) {
    bugs_next(st); bugs_expect(st, "(")
    link <- tk$text
    lhs <- bugs_lhs(st)
    bugs_expect(st, ")")
    if (!bugs_accept(st, "<-") && !bugs_accept(st, "=")) bugs_err(st, sprintf("a link function on the left side needs '<-' after %s(...)", link))
    return(list(kind = "logical", lhs = lhs, link = link, expr = bugs_expr(st), loc = tk$loc))
  }
  lhs <- bugs_lhs(st)
  if (bugs_accept(st, "~")) {
    if (!bugs_at(st, "name")) bugs_err(st, "expected a distribution after '~'")
    d <- bugs_next(st)
    bugs_expect(st, "(")
    args <- bugs_args(st, ")")
    bugs_expect(st, ")")
    bound <- NULL
    if (bugs_at(st, "name", c("T", "I", "C")) && bugs_at(st, "op", "(", 1L)) {
      b <- bugs_next(st)$text
      bugs_next(st)
      lo <- if (bugs_at(st, "op", ",")) NULL else bugs_expr(st)
      bugs_expect(st, ",")
      hi <- if (bugs_at(st, "op", ")")) NULL else bugs_expr(st)
      bugs_expect(st, ")")
      bound <- list(type = b, lo = lo, hi = hi)
    }
    return(list(kind = "stochastic", lhs = lhs, dist = d$text, args = args, bound = bound, loc = tk$loc))
  }
  if (bugs_accept(st, "<-") || bugs_accept(st, "=")) return(list(kind = "logical", lhs = lhs, link = NULL, expr = bugs_expr(st), loc = tk$loc))
  bugs_err(st, sprintf("expected '~' or '<-' after %s but found '%s'", lhs$name, bugs_peek(st)$text))
}

bugs_lhs <- function(st) {
  if (!bugs_at(st, "name")) bugs_err(st, "expected a variable name")
  nm <- bugs_next(st)
  idx <- NULL
  if (bugs_accept(st, "[")) {
    idx <- bugs_index_list(st)
    bugs_expect(st, "]")
  }
  list(kind = "var", name = nm$text, index = idx, loc = nm$loc)
}

# indices: expr, a:b, or empty (the whole extent)
bugs_index_list <- function(st) {
  out <- list()
  repeat {
    if (bugs_at(st, "op", c(",", "]"))) out[[length(out) + 1L]] <- list(kind = "all")
    else {
      lo <- bugs_expr(st)
      if (bugs_accept(st, ":")) out[[length(out) + 1L]] <- list(kind = "range", lo = lo, hi = bugs_expr(st))
      else out[[length(out) + 1L]] <- list(kind = "single", e = lo)
    }
    if (!bugs_accept(st, ",")) break
  }
  out
}

bugs_args <- function(st, close) {
  out <- list()
  if (bugs_at(st, "op", close)) return(out)
  repeat {
    out[[length(out) + 1L]] <- bugs_expr(st)
    if (!bugs_accept(st, ",")) break
  }
  out
}

# ---- expressions ------------------------------------------------------------------------
# precedence, lowest first: || , && , ! , comparisons , + - , * / , %op% ,
# unary minus , ^ (right associative)

bugs_expr <- function(st) bugs_or(st)
bugs_binary <- function(st, ops, next_fn) {
  a <- next_fn(st)
  while (bugs_at(st, "op", ops)) {
    op <- bugs_next(st)
    b <- next_fn(st)
    a <- list(kind = "binop", op = op$text, a = a, b = b, loc = op$loc)
  }
  a
}
bugs_or <- function(st) bugs_binary(st, c("||", "|"), bugs_and)
bugs_and <- function(st) bugs_binary(st, c("&&", "&"), bugs_not)
bugs_not <- function(st) {
  if (bugs_at(st, "op", "!")) { op <- bugs_next(st); return(list(kind = "unop", op = "!", a = bugs_not(st), loc = op$loc)) }
  bugs_cmp(st)
}
bugs_cmp <- function(st) bugs_binary(st, c("<", "<=", ">", ">=", "==", "!="), bugs_add)
bugs_add <- function(st) bugs_binary(st, c("+", "-"), bugs_mul)
bugs_mul <- function(st) bugs_binary(st, c("*", "/"), bugs_special)
bugs_special <- function(st) {
  a <- bugs_unary(st)
  while (identical(bugs_peek(st)$type, "op") && grepl("^%.*%$", bugs_peek(st)$text)) {
    op <- bugs_next(st)
    a <- list(kind = "binop", op = op$text, a = a, b = bugs_unary(st), loc = op$loc)
  }
  a
}
bugs_unary <- function(st) {
  if (bugs_at(st, "op", c("-", "+"))) {
    op <- bugs_next(st)
    a <- bugs_unary(st)
    if (op$text == "+") return(a)
    if (identical(a$kind, "num")) { a$value <- -a$value; return(a) }
    return(list(kind = "unop", op = "-", a = a, loc = op$loc))
  }
  bugs_pow(st)
}
bugs_pow <- function(st) {
  a <- bugs_primary(st)
  if (bugs_at(st, "op", "^")) {
    op <- bugs_next(st)
    # the exponent binds tighter than a unary minus in front of the base,
    # and may itself carry one: 2^-1
    return(list(kind = "binop", op = "^", a = a, b = bugs_unary(st), loc = op$loc))
  }
  a
}
bugs_primary <- function(st) {
  tk <- bugs_peek(st)
  if (identical(tk$type, "num")) { bugs_next(st); return(list(kind = "num", value = as.numeric(tk$text), loc = tk$loc)) }
  if (bugs_accept(st, "(")) {
    e <- bugs_expr(st)
    bugs_expect(st, ")")
    return(e)
  }
  if (identical(tk$type, "name")) {
    bugs_next(st)
    if (bugs_accept(st, "(")) {
      args <- bugs_args(st, ")")
      bugs_expect(st, ")")
      return(list(kind = "call", name = tk$text, args = args, loc = tk$loc))
    }
    idx <- NULL
    if (bugs_accept(st, "[")) {
      idx <- bugs_index_list(st)
      bugs_expect(st, "]")
    }
    return(list(kind = "var", name = tk$text, index = idx, loc = tk$loc))
  }
  bugs_err(st, sprintf("unexpected '%s' in an expression", tk$text))
}
