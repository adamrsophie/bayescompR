# Formatter for Stan programs (stanc's --auto-format, --print-canonical,
# --canonicalize and --max-line-length).

#' Format a Stan program
#'
#' Prints a Stan program in the layout of stanc's `--auto-format`: two space
#' indentation, one statement per line, spaced operators, and long
#' expressions broken before binary operators to fit `max_line_length`.
#' Comments and `#include` lines are kept.
#'
#' @param code Stan program as a character string, or `NULL` with `file`.
#' @param file Path to a `.stan` file.
#' @param canonicalize Changes made before printing, as stanc's
#'   `--canonicalize`: any of `"deprecations"` (nothing to do, since removed
#'   syntax is rejected), `"parentheses"` (drop redundant parentheses),
#'   `"braces"` (put the bodies of `if`, `else`, `for` and `while` in
#'   braces), `"includes"` (replace `#include` lines by the files) and
#'   `"strip-comments"`. `TRUE` applies all but `"strip-comments"`, as
#'   stanc's `--print-canonical`.
#' @param max_line_length Column at which long lines are broken.
#' @param include_paths Directories searched for `#include` files (used with
#'   `canonicalize = "includes"`).
#' @return The formatted program, a character string of class
#'   `bayescomp_code` that prints as plain text.
#' @examples
#' bayescomp_format("parameters{real mu;}model{mu~normal(0,1);}")
#' @export
bayescomp_format <- function(code = NULL, file = NULL, canonicalize = character(), max_line_length = 78L,
                             include_paths = NULL) {
  if (isTRUE(canonicalize)) canonicalize <- c("deprecations", "parentheses", "braces", "includes")
  if (isFALSE(canonicalize)) canonicalize <- character()
  bad <- setdiff(canonicalize, c("deprecations", "parentheses", "braces", "includes", "strip-comments"))
  if (length(bad)) stop(sprintf("unknown canonicalization '%s'", bad[1L]), call. = FALSE)
  if ("includes" %in% canonicalize) {
    src <- read_stan_source(code, file, include_paths)
  } else {
    if (is.null(code)) {
      if (is.null(file)) stop("supply either 'code' or 'file'", call. = FALSE)
      code <- readLines(file, warn = FALSE)
    }
    lines <- unlist(strsplit(paste(code, collapse = "\n"), "\n", fixed = TRUE))
    # keep #include lines as comments that are printed back unchanged
    inc <- grepl("^[[:space:]]*#include", lines)
    lines[inc] <- paste0("//#include", sub("^[[:space:]]*#include", "", lines[inc]))
    src <- paste(lines, collapse = "\n")
  }
  prog <- stan_parse(src, remap = FALSE)
  stan_check(prog)
  st <- new.env(parent = emptyenv())
  st$comments <- if ("strip-comments" %in% canonicalize) list() else prog$comments %||% list()
  st$width <- as.integer(max_line_length)
  st$parens <- "parentheses" %in% canonicalize
  st$braces <- "braces" %in% canonicalize
  out <- character()
  bl <- prog$blocks
  for (k in seq_along(bl)) {
    b <- bl[[k]]
    out <- c(out, fmt_comments_before(st, b$loc$line, ""))
    body <- if (b$name == "functions") lapply(b$body, fmt_fundef, st = st, ind = "  ") else fmt_stmts(b$body, st, "  ")
    nxt <- if (k < length(bl)) bl[[k + 1L]]$loc$line else .Machine$integer.max
    tail <- fmt_comments_before(st, nxt, "  ")
    body <- unlist(body)
    if (!length(body) && !length(tail)) body <- "  "
    out <- c(out, paste0(b$name, " {"), body, tail, "}")
  }
  out <- c(out, fmt_comments_before(st, .Machine$integer.max, ""))
  structure(paste0(paste(out, collapse = "\n"), "\n"), class = "bayescomp_code")
}

#' @export
print.bayescomp_code <- function(x, ...) {
  cat(x)
  invisible(x)
}

# comments that start before a line, printed on their own lines
fmt_comments_before <- function(st, line, ind) {
  out <- character()
  while (length(st$comments) && st$comments[[1L]]$line < line) {
    out <- c(out, paste0(ind, fmt_comment_text(st$comments[[1L]]$text)))
    st$comments <- st$comments[-1L]
  }
  out
}

fmt_comment_text <- function(t) if (startsWith(t, "//#include")) paste0("#include", substring(t, 11L)) else t

# a comment on the same line as the end of a statement goes after it
fmt_trailing <- function(st, line, lines) {
  if (length(st$comments) && st$comments[[1L]]$line == line && length(lines) && !(line %in% st$hold)) {
    lines[length(lines)] <- paste(lines[length(lines)], fmt_comment_text(st$comments[[1L]]$text))
    st$comments <- st$comments[-1L]
  }
  lines
}

# evaluate expr while comments on `line` are reserved for the enclosing
# statement
fmt_hold <- function(st, line, expr) {
  old <- st$hold
  st$hold <- c(old, line)
  on.exit(st$hold <- old)
  force(expr)
}

# consecutive declarations split from `real a, b;` print together again
fmt_merge_decls <- function(ss) {
  out <- list()
  for (s in ss) {
    n <- length(out)
    if (n && identical(s$kind, "decl") && identical(out[[n]]$kind, "decl") && identical(out[[n]]$type, s$type) && identical(out[[n]]$type$loc, s$type$loc)) {
      out[[n]]$names <- c(out[[n]]$names %||% list(list(name = out[[n]]$name, init = out[[n]]$init)), list(list(name = s$name, init = s$init)))
    } else out[[n + 1L]] <- s
  }
  out
}

fmt_last_line <- function(s) {
  # the last source line a node mentions
  m <- 0L
  walk <- function(n) {
    if (!is.list(n)) return()
    if (!is.null(n$line) && is.numeric(n$line)) m <<- max(m, n$line)
    for (e in n) if (is.list(e)) walk(e)
  }
  walk(s)
  m
}

fmt_stmts <- function(ss, st, ind) {
  # stanc prints an empty body as one indented blank line
  if (!length(ss)) return(ind)
  out <- character()
  for (s in fmt_merge_decls(ss)) {
    out <- c(out, fmt_comments_before(st, s$loc$line %||% 0L, ind))
    last <- fmt_last_line(s)
    # a compound statement written on one line keeps its comment after it
    lines <- fmt_hold(st, last, fmt_stmt(s, st, ind))
    out <- c(out, fmt_trailing(st, last, lines))
  }
  out
}

# ---- types -----------------------------------------------------------------------

fmt_sized_type <- function(t, st) {
  if (identical(t$base, "tuple")) {
    core <- paste0("tuple(", paste(vapply(t$elems, fmt_sized_type, "", st = st), collapse = ", "), ")")
  } else {
    cons <- t$cons %||% list()
    cs <- if (length(cons)) paste0("<", paste(sprintf("%s=%s", names(cons), vapply(cons, fmt_e, "", st = st)), collapse = ", "), ">") else ""
    sz <- if (length(t$sizes)) paste0("[", paste(vapply(t$sizes, fmt_e, "", st = st), collapse = ", "), "]") else ""
    core <- paste0(t$base, cs, sz)
  }
  if (length(t$adims)) paste0("array[", paste(vapply(t$adims, fmt_e, "", st = st), collapse = ", "), "] ", core) else core
}

fmt_unsized_type <- function(t) {
  core <- if (identical(t$base, "tuple")) paste0("tuple(", paste(vapply(t$elems, fmt_unsized_type, ""), collapse = ", "), ")") else t$base
  if ((t$nd %||% 0L) > 0L) paste0("array[", strrep(",", t$nd - 1L), "] ", core) else core
}

fmt_fundef <- function(fd, st, ind) {
  args <- vapply(fd$args, function(a) paste0(if (isTRUE(a$data_only)) "data " else "", fmt_unsized_type(a$type), " ", a$name), "")
  head <- paste0(ind, fmt_unsized_type(fd$rtype), " ", fd$name, "(", paste(args, collapse = ", "), ")")
  pre <- fmt_comments_before(st, fd$loc$line, ind)
  if (is.null(fd$body)) return(c(pre, fmt_trailing(st, fd$loc$line, paste0(head, ";"))))
  last <- fmt_last_line(fd)
  body <- fmt_hold(st, last, fmt_stmts(fd$body$body, st, paste0(ind, "  ")))
  c(pre, fmt_trailing(st, last, c(paste0(head, " {"), body, paste0(ind, "}"))))
}

# ---- statements ------------------------------------------------------------------

# a statement as lines; bodies of control statements go on the next line
# unless they are blocks
fmt_body <- function(head, body, st, ind) {
  if (st$braces && !identical(body$kind, "block_stmt")) body <- list(kind = "block_stmt", body = list(body), loc = body$loc)
  if (identical(body$kind, "block_stmt") && is.null(body$profile)) {
    inner <- fmt_stmts(body$body, st, paste0(ind, "  "))
    return(c(paste0(ind, head, " {"), inner, paste0(ind, "}")))
  }
  c(paste0(ind, head, " "), fmt_stmt(body, st, paste0(ind, "  ")))
}

fmt_stmt <- function(s, st, ind) {
  e <- function(x, col = nchar(ind)) fmt_e(x, st, col)
  line <- function(txt) paste0(ind, txt)
  switch(s$kind,
    decl = {
      ty <- fmt_sized_type(s$type, st)
      items <- s$names %||% list(list(name = s$name, init = s$init))
      parts <- vapply(items, function(it) if (is.null(it$init)) it$name else paste0(it$name, " = ", fmt_e(it$init, st, nchar(ind) + nchar(ty) + nchar(it$name) + 4L)), "")
      line(paste0(ty, " ", paste(parts, collapse = ", "), ";"))
    },
    assign = { l <- e(s$lhs); line(paste0(l, " ", s$op, " ", e(s$rhs, nchar(ind) + nchar(l) + nchar(s$op) + 2L), ";")) },
    unpack = line(paste0(e(s$lhs), " = ", e(s$rhs), ";")),
    tilde = {
      lhs <- e(s$lhs)
      args <- paste(vapply(s$args, e, ""), collapse = ", ")
      tr <- if (!is.null(s$trunc)) paste0(" T[", if (!is.null(s$trunc$lo)) e(s$trunc$lo) else "", ", ",
                                          if (!is.null(s$trunc$hi)) e(s$trunc$hi) else "", "]") else ""
      line(paste0(lhs, " ~ ", s$dist, "(", args, ")", tr, ";"))
    },
    target_add = line(paste0("target += ", e(s$value, nchar(ind) + 10L), ";")),
    jacobian_add = line(paste0("jacobian += ", e(s$value, nchar(ind) + 12L), ";")),
    call_stmt = line(paste0(e(s$call), ";")),
    print = , reject = , fatal_error = line(paste0(s$kind, "(", paste(vapply(s$items, e, ""), collapse = ", "), ");")),
    "return" = line(if (is.null(s$value)) "return;" else paste0("return ", e(s$value, nchar(ind) + 7L), ";")),
    "break" = line("break;"), "continue" = line("continue;"), empty = line(";"),
    block_stmt = {
      head <- if (!is.null(s$profile)) paste0("profile(\"", s$profile, "\") {") else "{"
      c(line(head), fmt_stmts(s$body, st, paste0(ind, "  ")), line("}"))
    },
    "if" = {
      out <- fmt_body(paste0("if (", e(s$cond, nchar(ind) + 4L), ")"), s$then, st, ind)
      if (!is.null(s$els)) {
        closed <- endsWith(out[length(out)], "}")
        if (identical(s$els$kind, "if") && !st$braces) {
          rest <- fmt_stmt(s$els, st, ind)
          rest[1L] <- sub(paste0("^", ind), "", rest[1L])
          if (closed) out[length(out)] <- paste0(out[length(out)], " else ", rest[1L]) else out <- c(out, paste0(ind, "else ", rest[1L]))
          out <- c(out, rest[-1L])
        } else {
          rest <- fmt_body("else", s$els, st, ind)
          rest[1L] <- sub(paste0("^", ind), "", rest[1L])
          if (closed) out[length(out)] <- paste0(out[length(out)], " ", rest[1L]) else out <- c(out, paste0(ind, rest[1L]))
          out <- c(out, rest[-1L])
        }
      }
      out
    },
    "for" = fmt_body(paste0("for (", s$var, " in ", e(s$lo), " : ", e(s$hi), ")"), s$body, st, ind),
    foreach = fmt_body(paste0("for (", s$var, " in ", e(s$iter), ")"), s$body, st, ind),
    "while" = fmt_body(paste0("while (", e(s$cond), ")"), s$body, st, ind),
    stop(sprintf("cannot format statement '%s'", s$kind), call. = FALSE))
}

# ---- expressions -----------------------------------------------------------------

fmt_atomic <- c("var", "int", "real", "imag", "call", "index", "paren", "proj", "array_expr", "rowvec_expr", "tuple_expr", "string", "transpose")

# expression text starting at column col; long chains of binary operators
# are broken before an operator and aligned with the start of the chain
fmt_e <- function(x, st, col = 0L) {
  k <- x$kind
  r <- function(y, c = col) fmt_e(y, st, c)
  if (k == "binop") {
    # flatten a left associative chain of operators of one precedence
    ops <- character(); xs <- list(); cur <- x
    p <- fmt_prec(x$op)
    while (identical(cur$kind, "binop") && fmt_prec(cur$op) == p && x$op != "^") {
      ops <- c(cur$op, ops); xs <- c(list(cur$b), xs); cur <- cur$a
    }
    if (!length(ops)) { ops <- x$op; xs <- list(x$b); cur <- x$a }
    first <- r(cur)
    txt <- first
    linew <- col + nchar(fmt_lastline(first))
    for (i in seq_along(ops)) {
      opnd <- r(xs[[i]], linew + nchar(ops[i]) + 2L)
      piece <- paste0(" ", ops[i], " ", opnd)
      if (!grepl("\n", piece) && linew + nchar(piece) <= st$width) {
        txt <- paste0(txt, piece)
        linew <- if (grepl("\n", piece)) nchar(fmt_lastline(piece)) else linew + nchar(piece)
      } else {
        opnd <- r(xs[[i]], col + nchar(ops[i]) + 1L)
        piece <- paste0(ops[i], " ", opnd)
        txt <- paste0(txt, "\n", strrep(" ", col), piece)
        linew <- if (grepl("\n", piece)) nchar(fmt_lastline(piece)) else col + nchar(piece)
      }
    }
    return(txt)
  }
  switch(k,
    var = x$name,
    int = , real = x$text %||% format(x$value),
    imag = paste0(x$text %||% format(x$value), "i"),
    string = paste0("\"", x$value, "\""),
    paren = if (st$parens && (x$a$kind %in% fmt_atomic)) r(x$a) else paste0("(", r(x$a, col + 1L), ")"),
    unop = paste0(x$op, r(x$a, col + 1L)),
    transpose = paste0(r(x$a), "'"),
    proj = paste0(r(x$x), ".", x$k),
    cond = paste0(r(x$cond), " ? ", r(x$a), " : ", r(x$b)),
    array_expr = paste0("{", paste(vapply(x$items, r, ""), collapse = ", "), "}"),
    rowvec_expr = paste0("[", paste(vapply(x$items, r, ""), collapse = ", "), "]"),
    tuple_expr = paste0("(", paste(vapply(x$items, r, ""), collapse = ", "), if (length(x$items) == 1L) ",", ")"),
    index = paste0(r(x$x), "[", paste(vapply(x$idx, fmt_index, "", st = st), collapse = ", "), "]"),
    call = {
      a <- vapply(x$args, r, "")
      if (isTRUE(x$cond_bar) && length(a)) paste0(x$name, "(", a[1L], if (length(a) > 1L) paste0(" | ", paste(a[-1L], collapse = ", ")) else " |", ")")
      else paste0(x$name, "(", paste(a, collapse = ", "), ")")
    },
    stop(sprintf("cannot format expression '%s'", k), call. = FALSE))
}

fmt_firstline <- function(s) sub("\n.*$", "", s)
fmt_lastline <- function(s) sub("^.*\n", "", s)

fmt_index <- function(i, st) {
  switch(i$kind,
    omni = " : ",
    single = , intarray = fmt_e(i$e, st),
    lower = paste0(fmt_e(i$lo, st), " : "),
    upper = paste0(" : ", fmt_e(i$hi, st)),
    range = paste0(fmt_e(i$lo, st), " : ", fmt_e(i$hi, st)))
}

fmt_prec <- function(op) {
  switch(op, "||" = 1L, "&&" = 2L, "==" = , "!=" = 3L, "<" = , "<=" = , ">" = , ">=" = 4L,
         "+" = , "-" = 5L, "*" = , "/" = , "%" = , "%/%" = 6L, "\\" = 7L, ".*" = , "./" = 8L, "^" = , ".^" = 10L, 0L)
}
