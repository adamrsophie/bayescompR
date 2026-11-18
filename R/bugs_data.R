# Data and initial values for BUGS models.
#
# From R, data and initial values are lists of numbers, vectors and arrays
# with R's own layout (`x[i, j]` is `x[i, j]`), as R2jags, rjags and
# R2OpenBUGS take them. From text they come in one of three formats:
#
# * OpenBUGS lists, `list(N = 3, Y = structure(.Data = c(...), .Dim = c(2, 3)))`.
#   OpenBUGS fills an array from `.Data` with the last index changing
#   fastest (row major), the reverse of R (OpenBUGS manual, "Formatting of
#   data").
# * R dump files as written by `dump()` and read by JAGS: `N <- 3` and
#   `Y <- structure(c(...), .Dim = c(2L, 3L))`, in R's column major order.
# * OpenBUGS rectangular files: a header of `name[]` or `Y[, 2]` columns,
#   rows of numbers and a final `END`.
#
# The text is parsed with R's parser but evaluated by a small interpreter
# that allows only constants, c(), list(), structure(), rep(), seq(), `:`
# and arithmetic, so a data file cannot run code.

bugs_read_values <- function(x, what = "data") {
  if (is.null(x)) return(list())
  if (is.list(x)) return(bugs_check_values(x, what))
  if (is.environment(x)) return(bugs_check_values(as.list(x), what))
  if (!is.character(x)) stop_data(sprintf("%s must be a list or the name of a file", what))
  # a character vector of variable names, as R2jags and R2OpenBUGS accept,
  # names objects in the calling environment; this is resolved by the caller
  txt <- if (length(x) == 1L && file.exists(x)) readLines(x, warn = FALSE) else x
  bugs_parse_values(paste(txt, collapse = "\n"), what)
}

bugs_parse_values <- function(text, what = "data") {
  if (grepl("^\\s*$", text)) return(list())
  if (bugs_is_rectangular(text)) return(bugs_read_rectangular(text, what))
  # OpenBUGS reads `.Dim` split over lines (Jaws data does this)
  text <- gsub("\\.\\s+(Dim|Data)\\b", ".\\1", text, perl = TRUE)
  ex <- tryCatch(parse(text = text, keep.source = FALSE), error = function(e) stop_data(sprintf("cannot parse the %s: %s", what, conditionMessage(e))))
  out <- list()
  for (e in as.list(ex)) {
    if (is.call(e) && identical(e[[1L]], as.name("list"))) {
      # OpenBUGS list: structures are row major
      v <- bugs_eval_value(e, row_major = TRUE)
      for (nm in names(v)) out[[nm]] <- v[[nm]]
    } else if (is.call(e) && as.character(e[[1L]]) %in% c("<-", "=") && length(e) == 3L) {
      nm <- if (is.character(e[[2L]])) e[[2L]] else as.character(e[[2L]])
      out[[nm]] <- bugs_eval_value(e[[3L]], row_major = FALSE)
    } else stop_data(sprintf("the %s must be a list(...) or assignments 'name <- value'", what))
  }
  bugs_check_values(out, what)
}

bugs_eval_value <- function(e, row_major) {
  ev <- function(e) {
    if (is.numeric(e) || is.logical(e)) return(e)
    if (is.character(e)) stop_data(sprintf("text value '%s' in data", e))
    if (is.name(e)) {
      nm <- as.character(e)
      if (nm %in% c("NA", "TRUE", "FALSE", "Inf")) return(switch(nm, `NA` = NA, `TRUE` = TRUE, `FALSE` = FALSE, `Inf` = Inf))
      stop_data(sprintf("unknown name '%s' in data", nm))
    }
    if (!is.call(e)) stop_data("unsupported value in data")
    f <- as.character(e[[1L]])
    args <- as.list(e)[-1L]
    switch(f,
      "(" = ev(args[[1L]]),
      "-" = if (length(args) == 1L) -ev(args[[1L]]) else ev(args[[1L]]) - ev(args[[2L]]),
      "+" = if (length(args) == 1L) ev(args[[1L]]) else ev(args[[1L]]) + ev(args[[2L]]),
      "*" = ev(args[[1L]]) * ev(args[[2L]]),
      "/" = ev(args[[1L]]) / ev(args[[2L]]),
      ":" = ev(args[[1L]]):ev(args[[2L]]),
      "c" = unlist(lapply(args, ev), use.names = FALSE),
      "rep" = do.call(rep, lapply(args, ev)),
      "seq" = do.call(seq, lapply(args, ev)),
      "as.integer" = as.integer(ev(args[[1L]])),
      "as.double" = as.double(ev(args[[1L]])),
      "list" = {
        vals <- lapply(args, ev)
        if (is.null(names(args)) || any(names(args) == "")) stop_data("every element of a data list needs a name")
        vals
      },
      "structure" = {
        nms <- names(args) %||% rep("", length(args))
        dat <- ev(args[[if (".Data" %in% nms) which(nms == ".Data") else 1L]])
        dimi <- which(nms %in% c(".Dim", "dim"))
        if (!length(dimi)) return(dat)
        d <- as.integer(ev(args[[dimi[1L]]]))
        if (prod(d) != length(dat)) stop_data(sprintf("structure with %d values and dimensions %s", length(dat), paste(d, collapse = " x ")))
        if (row_major && length(d) > 1L) aperm(array(dat, rev(d)), rev(seq_along(d))) else array(dat, d)
      },
      stop_data(sprintf("function '%s' is not allowed in data", f)))
  }
  ev(e)
}

bugs_is_rectangular <- function(text) {
  first <- sub("^\\s+", "", strsplit(text, "\n")[[1L]])
  first <- first[nzchar(first) & !startsWith(first, "#")][1L]
  !is.na(first) && grepl("^[A-Za-z.][A-Za-z0-9._]*\\[[^]]*\\]", first) && !grepl("<-|=|list\\(", first)
}

bugs_read_rectangular <- function(text, what) {
  lines <- trimws(strsplit(text, "\n")[[1L]])
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  end <- which(toupper(lines) == "END")
  if (!length(end)) stop_data(sprintf("rectangular %s must end with END", what))
  head <- strsplit(lines[1L], "\\s+(?![^[]*\\])", perl = TRUE)[[1L]]
  head <- head[nzchar(head)]
  rows <- lines[seq.int(2L, end[1L] - 1L)]
  vals <- lapply(rows, function(r) {
    tk <- strsplit(r, "[[:space:],]+")[[1L]]
    tk <- tk[nzchar(tk)]
    if (length(tk) != length(head)) stop_data(sprintf("a row of the rectangular %s has %d values for %d columns", what, length(tk), length(head)))
    suppressWarnings(as.numeric(ifelse(tk == "NA", NA, tk)))
  })
  M <- do.call(rbind, vals)
  out <- list()
  for (k in seq_along(head)) {
    h <- gsub("\\s", "", head[k])
    nm <- sub("\\[.*$", "", h)
    inner <- sub("^[^[]*\\[(.*)\\]$", "\\1", h)
    parts <- strsplit(inner, ",", fixed = TRUE)[[1L]]
    if (length(parts) && nzchar(parts[1L])) stop_data(sprintf("the first index of column '%s' must be empty", h))
    if (length(parts) <= 1L) { out[[nm]] <- M[, k]; next }
    pos <- as.integer(parts[-1L])
    cur <- out[[nm]]
    if (is.null(cur)) cur <- list()
    cur[[length(cur) + 1L]] <- list(pos = pos, v = M[, k])
    out[[nm]] <- cur
  }
  # columns of a multidimensional array: assemble them, NA where absent
  for (nm in names(out)) if (is.list(out[[nm]])) {
    cols <- out[[nm]]
    ext <- apply(do.call(rbind, lapply(cols, `[[`, "pos")), 2L, max)
    a <- array(NA_real_, c(nrow(M), ext))
    for (cl in cols) a[as.matrix(cbind(seq_len(nrow(M)), matrix(cl$pos, nrow(M), length(cl$pos), byrow = TRUE)))] <- cl$v
    out[[nm]] <- a
  }
  out
}

bugs_check_values <- function(x, what) {
  if (length(x) && (is.null(names(x)) || any(!nzchar(names(x))))) stop_data(sprintf("every element of the %s needs a name", what))
  for (nm in names(x)) {
    v <- x[[nm]]
    if (is.data.frame(v)) v <- as.matrix(v)
    if (is.factor(v)) v <- as.integer(v)
    if (is.logical(v)) v <- as.integer(v)
    if (!is.numeric(v)) stop_data(sprintf("%s '%s' is not numeric", what, nm))
    x[[nm]] <- v
  }
  x
}
