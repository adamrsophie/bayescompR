# Builds R/sysdata.rda: the library signatures of stanc 2.40.0, used by the
# semantic analyzer so that calls stanc rejects are rejected here too.
# Run from the package root with STANC pointing to the stanc binary:
#   STANC=~/.cmdstan/cmdstan-2.40.0/bin/stanc Rscript data-raw/stan_signatures.R
stanc <- Sys.getenv("STANC", "stanc")
lines <- system2(path.expand(stanc), "--dump-stan-math-signatures", stdout = TRUE)
lines <- lines[grepl("^[A-Za-z_][A-Za-z0-9_]*\\(", lines)]
split_top <- function(s) {
  if (!nzchar(s)) return(character(0))
  ch <- strsplit(s, "")[[1L]]
  depth <- cumsum((ch %in% c("(", "[")) - (ch %in% c(")", "]")))
  cut <- which(ch == "," & depth == 0L)
  trimws(substring(s, c(1L, cut + 1L), c(cut - 1L, nchar(s))))
}
enc <- function(a) {
  a <- sub("^data ", "", a)
  if (grepl("[(]|tuple", a)) return(NA_character_)
  nd <- 0L
  if (grepl("^array\\[", a)) {
    nd <- nchar(gsub("[^,]", "", sub("\\].*", "", a))) + 1L
    a <- sub("^array\\[,*\\] ", "", a)
  }
  paste0(a, ":", nd)
}
nm <- sub("\\(.*", "", lines)
args <- vapply(lines, function(l) {
  inner <- sub(") => .*$", "", sub("^[^(]*\\(", "", l))
  e <- vapply(split_top(inner), enc, "")
  if (anyNA(e)) NA_character_ else paste(e, collapse = ",")
}, "", USE.NAMES = FALSE)
keep <- !startsWith(nm, "operator")
nm <- nm[keep]; args <- args[keep]
# functions with function or tuple arguments are checked by their own rules
skip <- unique(nm[is.na(args)])
stan_signatures <- lapply(split(args, nm)[setdiff(unique(nm), skip)], unique)
stan_signatures_version <- "2.40.0"
save(stan_signatures, stan_signatures_version, file = "R/sysdata.rda", compress = "xz")
cat(length(stan_signatures), "functions,", sum(lengths(stan_signatures)), "signatures\n")
