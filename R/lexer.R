# Tokenizer for the Stan language (Stan 2.40 grammar).
#
# Produces a list of tokens. Each token is a list with fields
#   type   one of "ident", "int", "real", "imag", "string", "op", "eof"
#   value  the token text (for literals the source text)
#   line, col  1 based source position of the first character

stan_reserved_words <- c(
  # Stan keywords
  "for", "in", "while", "repeat", "until", "if", "then", "else", "true", "false",
  "target", "functions", "model", "data", "parameters", "quantities",
  "transformed", "generated", "profile", "return", "break", "continue",
  "print", "reject", "fatal_error", "void", "int", "real", "complex", "vector",
  "row_vector", "matrix", "array", "tuple", "ordered", "positive_ordered",
  "simplex", "unit_vector", "sum_to_zero_vector", "cholesky_factor_corr",
  "cholesky_factor_cov", "corr_matrix", "cov_matrix", "complex_vector",
  "complex_row_vector", "complex_matrix", "column_stochastic_matrix",
  "row_stochastic_matrix", "sum_to_zero_matrix", "lower", "upper", "offset", "multiplier",
  "struct", "typedef", "export", "auto", "extern", "var", "static",
  "jacobian"
)

# multi character operators, longest first
stan_operators <- c(
  ".*=", "./=", "%/%",
  "<-", "+=", "-=", "*=", "/=", "==", "!=", "<=", ">=", "&&", "||", ".*", "./", ".^",
  "+", "-", "*", "/", "%", "\\", "^", "'", "!", "<", ">", "=", "~", "?", ":",
  ",", ";", "(", ")", "[", "]", "{", "}", "|"
)

stan_lex <- function(src) {
  if (length(src) > 1L) src <- paste(src, collapse = "\n")
  chars <- strsplit(src, "", fixed = TRUE)[[1L]]
  n <- length(chars)
  toks <- vector("list", 256L)
  nt <- 0L
  i <- 1L
  line <- 1L
  col <- 1L
  add <- function(type, value, l, c) {
    nt <<- nt + 1L
    if (nt > length(toks)) length(toks) <<- 2L * length(toks)
    toks[[nt]] <<- list(type = type, value = value, line = l, col = c)
  }
  advance <- function(k = 1L) {
    for (j in seq_len(k)) {
      if (i <= n && chars[i] == "\n") {
        line <<- line + 1L
        col <<- 1L
      } else {
        col <<- col + 1L
      }
      i <<- i + 1L
    }
  }
  peek <- function(k = 0L) if (i + k <= n) chars[i + k] else ""
  is_digit <- function(ch) ch %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
  is_alpha <- function(ch) grepl("^[A-Za-z_]$", ch)
  ops3 <- stan_operators[nchar(stan_operators) == 3L]
  ops2 <- stan_operators[nchar(stan_operators) == 2L]
  ops1 <- stan_operators[nchar(stan_operators) == 1L]

  while (i <= n) {
    ch <- chars[i]
    if (ch %in% c(" ", "\t", "\n", "\r", "\f")) {
      advance()
      next
    }
    l0 <- line
    c0 <- col
    # comments
    if (ch == "/" && peek(1L) == "/") {
      while (i <= n && chars[i] != "\n") advance()
      next
    }
    if (ch == "/" && peek(1L) == "*") {
      advance(2L)
      closed <- FALSE
      while (i <= n) {
        if (chars[i] == "*" && peek(1L) == "/") {
          advance(2L)
          closed <- TRUE
          break
        }
        advance()
      }
      if (!closed) stop_syntax("unterminated block comment", list(line = l0, col = c0))
      next
    }
    if (ch == "#") {
      rest <- paste(chars[i:min(n, i + 7L)], collapse = "")
      if (startsWith(rest, "#include"))
        stop_unsupported("#include directives are not supported; paste the included file into the program", list(line = l0, col = c0))
      stop_syntax("'#' comments were removed in Stan 2.33; use '//' instead", list(line = l0, col = c0))
    }
    # tuple projection: x.1, f(x).2, t.2.1
    if (ch == "." && is_digit(peek(1L)) && nt > 0L) {
      pt <- toks[[nt]]
      proj_ok <- pt$type == "ident" || (pt$type == "op" && pt$value %in% c(")", "]")) ||
        (pt$type == "int" && nt > 1L && identical(toks[[nt - 1L]]$value, ".") && toks[[nt - 1L]]$type == "op")
      if (proj_ok) {
        advance()
        add("op", ".", l0, c0)
        start <- i
        while (i <= n && is_digit(chars[i])) advance()
        add("int", paste(chars[start:(i - 1L)], collapse = ""), l0, c0 + 1L)
        next
      }
    }
    # numbers
    if (is_digit(ch) || (ch == "." && is_digit(peek(1L)))) {
      start <- i
      is_real <- FALSE
      # Stan allows _ as a digit separator, between two digits only
      digit_or_sep <- function() is_digit(chars[i]) || (chars[i] == "_" && i > start && is_digit(chars[i - 1L]) && is_digit(peek(1L)))
      while (i <= n && digit_or_sep()) advance()
      if (i <= n && chars[i] == "." && !(peek(1L) %in% c("*", "/", "^"))) {
        is_real <- TRUE
        advance()
        while (i <= n && digit_or_sep()) advance()
      }
      if (i <= n && chars[i] %in% c("e", "E")) {
        j <- 1L
        if (peek(1L) %in% c("+", "-")) j <- 2L
        if (is_digit(peek(j))) {
          is_real <- TRUE
          advance(j)
          while (i <= n && digit_or_sep()) advance()
        } else {
          stop_syntax("malformed exponent in numeric literal", list(line = l0, col = c0))
        }
      }
      text <- gsub("_", "", paste(chars[start:(i - 1L)], collapse = ""), fixed = TRUE)
      if (i <= n && chars[i] == "i" && !(i < n && (is_alpha(chars[i + 1L]) || is_digit(chars[i + 1L])))) {
        advance()
        add("imag", text, l0, c0)
        next
      }
      if (i <= n && (is_alpha(chars[i]))) stop_syntax(sprintf("unexpected character '%s' after numeric literal", chars[i]), list(line = line, col = col))
      add(if (is_real) "real" else "int", text, l0, c0)
      next
    }
    # identifiers and keywords
    if (is_alpha(ch)) {
      start <- i
      while (i <= n && (is_alpha(chars[i]) || is_digit(chars[i]))) advance()
      text <- paste(chars[start:(i - 1L)], collapse = "")
      add("ident", text, l0, c0)
      next
    }
    # strings (print, reject, profile)
    if (ch == "\"") {
      advance()
      start <- i
      while (i <= n && chars[i] != "\"") {
        if (chars[i] == "\n") stop_syntax("unterminated string literal", list(line = l0, col = c0))
        advance()
      }
      if (i > n) stop_syntax("unterminated string literal", list(line = l0, col = c0))
      text <- if (i > start) paste(chars[start:(i - 1L)], collapse = "") else ""
      advance()
      add("string", text, l0, c0)
      next
    }
    # operators
    s3 <- paste(chars[i:min(n, i + 2L)], collapse = "")
    s2 <- substr(s3, 1L, 2L)
    if (nchar(s3) == 3L && s3 %in% ops3) {
      advance(3L)
      add("op", s3, l0, c0)
      next
    }
    if (nchar(s2) == 2L && s2 %in% ops2) {
      if (s2 == "<-") stop_syntax("the '<-' assignment operator was removed from Stan; use '='", list(line = l0, col = c0))
      advance(2L)
      add("op", s2, l0, c0)
      next
    }
    if (ch %in% ops1) {
      advance()
      add("op", ch, l0, c0)
      next
    }
    stop_syntax(sprintf("unexpected character '%s'", ch), list(line = l0, col = c0))
  }
  add("eof", "", line, col)
  toks[seq_len(nt)]
}
