# Replay the signatures of stanc --dump-stan-math-signatures through
# bayescomp_check(): Rscript validation/signatures/sigdiff.R <dump.txt> [out.csv]
# Signatures with complex, tuple or function arguments, operators and the
# transform functions are skipped (they are counted separately).
suppressMessages(pkgload::load_all(".", quiet = TRUE))
a <- commandArgs(TRUE)
lines <- readLines(a[1])
split_top <- function(s) {
  if (!nzchar(s)) return(character())
  ch <- strsplit(s, "")[[1]]; depth <- 0L; out <- character(); cur <- ""
  for (c in ch) {
    if (c %in% c("[", "(")) depth <- depth + 1L
    if (c %in% c("]", ")")) depth <- depth - 1L
    if (c == "," && depth == 0L) { out <- c(out, trimws(cur)); cur <- "" } else cur <- paste0(cur, c)
  }
  c(out, trimws(cur))
}
parse_sig <- function(l) {
  m <- regmatches(l, regexec("^([A-Za-z_0-9]+)\\((.*)\\) => (.*)$", l))[[1]]
  if (!length(m)) return(NULL)
  list(name = m[2], args = split_top(m[3]), ret = m[4], line = l)
}
decl <- function(t, nm) {
  t <- sub("^data ", "", t)
  arr <- regmatches(t, regexec("^array\\[([,]*)\\] (.*)$", t))[[1]]
  nd <- if (length(arr)) nchar(arr[2]) + 1L else 0L
  base <- if (length(arr)) arr[3] else t
  el <- switch(base, int = "int", real = "real", vector = "vector[2]", row_vector = "row_vector[2]", matrix = "matrix[2, 2]", NA)
  if (is.na(el)) return(NA_character_)
  init <- switch(base, int = "1", real = "0.5", vector = "[0.2, 0.3]'", row_vector = "[0.2, 0.3]", matrix = "[[1, 0.1], [0.1, 1]]")
  if (nd) {
    for (k in seq_len(nd)) init <- sprintf("{%s, %s}", init, init)
    sprintf("array[%s] %s %s = %s;", paste(rep("2", nd), collapse = ", "), el, nm, init)
  } else sprintf("%s %s = %s;", el, nm, init)
}
ops <- c("add", "subtract", "multiply", "divide", "minus", "plus", "elt_multiply", "elt_divide", "modulus", "transpose",
         "logical_negation", "logical_and", "logical_or", "logical_eq", "logical_neq", "logical_lt", "logical_lte", "logical_gt", "logical_gte", "pow")
res <- list()
for (l in lines) {
  s <- parse_sig(l)
  if (is.null(s)) next
  if (s$name %in% ops || grepl("_(constrain|unconstrain|jacobian)$", s$name)) next
  if (any(grepl("complex|tuple|\\(|=>", s$args)) || grepl("complex|tuple", s$ret)) next
  ds <- vapply(seq_along(s$args), function(i) decl(s$args[i], paste0("a", i)), "")
  if (anyNA(ds)) next
  argn <- if (length(s$args)) paste0("a", seq_along(s$args)) else character()
  density <- grepl("_(lpdf|lpmf|lupdf|lupmf)$", s$name)
  cdf <- grepl("_(lcdf|lccdf|cdf)$", s$name)
  call <- if ((density || cdf) && length(argn)) sprintf("%s(%s | %s)", s$name, argn[1], paste(argn[-1], collapse = ", ")) else sprintf("%s(%s)", s$name, paste(argn, collapse = ", "))
  call <- sub("\\|  ?\\)", "|)", call)
  rng <- grepl("_rng$", s$name)
  lp <- grepl("_lp$", s$name) || grepl("_lupdf$|_lupmf$", s$name)
  stmt <- if (s$ret == "void") paste0(call, ";") else sprintf("print(%s);", call)
  body <- paste(c(ds, stmt), collapse = " ")
  prog <- if (rng) sprintf("generated quantities { %s }", body)
          else if (lp) sprintf("parameters { real p__; } model { %s }", body)
          else sprintf("transformed data { %s }", body)
  err <- tryCatch({ bayescomp_check(prog); "" }, error = function(e) paste(class(e)[1], conditionMessage(e)))
  res[[length(res) + 1L]] <- data.frame(name = s$name, sig = l, ok = !nzchar(err), error = substr(err, 1, 200))
}
r <- do.call(rbind, res)
cat(sprintf("tested %d signatures over %d names; accepted %d (%.1f%%); rejected %d in %d names\n", nrow(r), length(unique(r$name)),
            sum(r$ok), 100 * mean(r$ok), sum(!r$ok), length(unique(r$name[!r$ok]))))
if (length(a) > 1) utils::write.csv(r[!r$ok, ], a[2], row.names = FALSE)
