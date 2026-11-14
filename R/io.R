# Data and initial value files in CmdStan's JSON and R dump formats

#' Read and write data files
#'
#' `bayescomp_read_data()` reads a data or initial values file in either of
#' CmdStan's formats: JSON (the `.json` extension, or content starting with
#' `{`) or R dump (`name <- value` assignments, as written by
#' [base::dump()]). `bayescomp_write_data()` writes a named list as CmdStan
#' JSON. Any function that takes `data` or `init` also accepts a file name.
#'
#' JSON arrays are row major: `[[1, 2], [3, 4]]` is a matrix whose first row
#' is `(1, 2)`, and `array[N] vector[K]` is an `N` by `K` matrix in R.
#' `Infinity`, `-Infinity` and `NaN` are accepted bare or as strings (also
#' `"inf"`, `"+inf"`, `"-inf"`), complex values are `[re, im]` pairs, and
#' tuples are objects with keys `"1"`, `"2"`, ... . In R dump files, arrays
#' are column major with a `.Dim` attribute, as R stores them.
#'
#' @param file Path to the file.
#' @param data Named list of values: numbers, vectors, matrices, arrays and
#'   complex values. A tuple is a list named `"1"`, `"2"`, ... (written as a
#'   JSON object); other lists are written as JSON arrays.
#' @return `bayescomp_read_data()` returns a named list;
#'   `bayescomp_write_data()` returns `file`, invisibly.
#' @examples
#' f <- tempfile(fileext = ".json")
#' bayescomp_write_data(list(N = 2L, y = c(1.5, Inf), X = diag(2)), f)
#' readLines(f)
#' str(bayescomp_read_data(f))
#' @export
bayescomp_read_data <- function(file) {
  if (!is.character(file) || length(file) != 1L || !file.exists(file))
    stop(sprintf("data file '%s' does not exist", paste(file, collapse = ", ")), call. = FALSE)
  txt <- paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (grepl("\\.json$", file, ignore.case = TRUE) || grepl("^\\s*\\{", txt)) read_stan_json(txt, file) else read_rdump(txt, file)
}

#' @rdname bayescomp_read_data
#' @export
bayescomp_write_data <- function(data, file) {
  if (!is.list(data) || (length(data) && (is.null(names(data)) || any(!nzchar(names(data))))))
    stop("'data' must be a named list", call. = FALSE)
  body <- vapply(names(data), function(nm) paste0("  ", json_string(nm), ": ", json_value(data[[nm]])), "")
  writeLines(c("{", paste(body, collapse = ",\n"), "}"), file)
  invisible(file)
}

# ---- JSON ------------------------------------------------------------------

read_stan_json <- function(txt, file) {
  # bare Infinity and NaN are not standard JSON: quote them outside strings
  m <- gregexpr('"(?:[^"\\\\]|\\\\.)*"', txt, perl = TRUE)
  outside <- regmatches(txt, m, invert = TRUE)[[1L]]
  outside <- gsub("(-?)\\b(Infinity|NaN)\\b", "\"\\1\\2\"", outside, perl = TRUE)
  regmatches(txt, m, invert = TRUE) <- list(outside)
  x <- tryCatch(jsonlite::parse_json(txt, simplifyVector = FALSE),
                error = function(e) stop_data(sprintf("cannot parse JSON file '%s': %s", file, conditionMessage(e))))
  if (!is.list(x) || (length(x) && is.null(names(x))))
    stop_data(sprintf("JSON file '%s' must contain one object of name and value pairs", file))
  if (anyDuplicated(names(x))) stop_data(sprintf("JSON file '%s' defines '%s' more than once", file, names(x)[anyDuplicated(names(x))]))
  lapply(stats::setNames(names(x), names(x)), function(nm) json_convert(x[[nm]], nm))
}

json_special <- function(s, nm) {
  t <- tolower(s)
  if (t %in% c("inf", "+inf", "infinity", "+infinity")) return(Inf)
  if (t %in% c("-inf", "-infinity")) return(-Inf)
  if (t %in% c("nan", "+nan", "-nan")) return(NaN)
  s  # other strings (such as metric_type in metric files) stay strings
}

# leaves become numbers; arrays stay lists marked as coming from JSON, so
# that complex declarations read [re, im] pairs; objects become tuples
json_convert <- function(v, nm) {
  if (is.null(v)) stop_data(sprintf("variable '%s' is null", nm))
  if (is.character(v)) return(json_special(v, nm))
  if (is.logical(v)) return(as.integer(v))
  if (!is.list(v)) return(v)
  out <- lapply(v, json_convert, nm = nm)
  if (!is.null(names(v))) {
    k <- suppressWarnings(as.integer(names(v)))
    if (anyNA(k) || !setequal(k, seq_along(k))) stop_data(sprintf("variable '%s': tuple keys must be \"1\", \"2\", ...", nm))
    out <- out[order(k)]
    names(out) <- seq_along(out)
    return(out)
  }
  attr(out, "stan_json") <- TRUE
  out
}

# [re, im] pairs to complex values, for complex declarations
json_complex <- function(x) {
  if (!is.list(x)) return(x)
  if (length(x) == 2L && !is.list(x[[1L]]) && !is.list(x[[2L]])) return(complex(real = x[[1L]], imaginary = x[[2L]]))
  lapply(x, json_complex)
}

json_string <- function(s) paste0('"', gsub('(["\\\\])', "\\\\\\1", s), '"')

json_num <- function(x) {
  out <- formatC(x, digits = 15, format = "g")
  redo <- is.finite(x) & as.numeric(out) != x
  out[redo] <- formatC(x[redo], digits = 17, format = "g")
  out[is.nan(x)] <- '"NaN"'
  out[is.infinite(x)] <- ifelse(x[is.infinite(x)] > 0, '"Infinity"', '"-Infinity"')
  out[is.na(x) & !is.nan(x)] <- "null"
  gsub(" ", "", out)
}

# nested row major arrays; the last dimension varies fastest
json_nested <- function(x, dims, leaf) {
  if (length(dims) == 1L) return(paste0("[", paste(leaf(x), collapse = ", "), "]"))
  inner <- vapply(seq_len(dims[1L]), function(i) {
    sl <- x[(seq_len(prod(dims[-1L])) - 1L) * dims[1L] + i]
    json_nested(sl, dims[-1L], leaf)
  }, "")
  paste0("[", paste(inner, collapse = ", "), "]")
}

json_value <- function(v) {
  if (is.data.frame(v)) v <- as.matrix(v)
  if (is.list(v)) {
    tup <- inherits(v, "bctuple") || (!is.null(names(v)) && setequal(names(v), as.character(seq_along(v))))
    if (tup) {
      if (!is.null(names(v))) v <- v[as.character(seq_along(v))]
      return(paste0("{", paste(sprintf('"%d": %s', seq_along(v), vapply(unclass(v), json_value, "")), collapse = ", "), "}"))
    }
    return(paste0("[", paste(vapply(v, json_value, ""), collapse = ", "), "]"))
  }
  if (is.logical(v)) v <- as.integer(v)
  if (is.character(v)) return(if (length(v) == 1L) json_string(v) else paste0("[", paste(json_string(v), collapse = ", "), "]"))
  leaf <- if (is.complex(v)) function(z) sprintf("[%s, %s]", json_num(Re(z)), json_num(Im(z)))
          else if (is.integer(v)) function(z) ifelse(is.na(z), "null", as.character(z))
          else json_num
  d <- dim(v)
  if (is.null(d)) {
    if (length(v) == 1L && !inherits(v, "AsIs")) return(leaf(v))
    return(paste0("[", paste(leaf(v), collapse = ", "), "]"))
  }
  if (prod(d) == 0L) return("[]")
  json_nested(as.vector(v), d, leaf)
}

# ---- R dump -----------------------------------------------------------------

read_rdump <- function(txt, file) {
  ex <- tryCatch(parse(text = txt, keep.source = FALSE),
                 error = function(e) stop_data(sprintf("cannot parse R dump file '%s': %s", file, conditionMessage(e))))
  # evaluated in an environment that holds only the constructors R dump
  # files use (its parent is the empty environment), so no other code runs
  safe <- new.env(parent = emptyenv())
  # R dump files write .Dim, which structure() no longer accepts silently
  assign("structure", function(.Data, ...) {
    a <- list(...)
    names(a)[names(a) == ".Dim"] <- "dim"
    names(a)[names(a) == ".Dimnames"] <- "dimnames"
    do.call(base::structure, c(list(.Data), a))
  }, envir = safe)
  for (f in c("c", "list", ":", "-", "+", "integer", "double", "numeric", "seq", "rep", "array", "matrix", "as.integer", "as.double"))
    assign(f, get(f, envir = baseenv()), envir = safe)
  assign("Inf", Inf, envir = safe); assign("NaN", NaN, envir = safe)
  assign("T", TRUE, envir = safe); assign("F", FALSE, envir = safe)
  out <- list()
  for (e in ex) {
    if (!is.call(e) || !(as.character(e[[1L]]) %in% c("<-", "=")) || length(e) != 3L)
      stop_data(sprintf("R dump file '%s' may only contain assignments", file))
    nm <- as.character(e[[2L]])
    v <- tryCatch(eval(e[[3L]], safe), error = function(err) stop_data(sprintf("R dump file '%s', variable '%s': %s", file, nm, conditionMessage(err))))
    d <- attr(v, ".Dim") %||% dim(v)
    attributes(v) <- NULL
    if (!is.null(d)) dim(v) <- d
    out[[nm]] <- v
  }
  out
}

# data or init given as a file name
data_arg <- function(x, what = "data") {
  if (is.character(x)) {
    if (length(x) != 1L) stop(sprintf("'%s' must be one file name", what), call. = FALSE)
    return(bayescomp_read_data(x))
  }
  x
}
