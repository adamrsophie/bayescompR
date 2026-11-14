# Condition classes. Every error raised by the package carries one of these
# classes so that callers can tell syntax, semantic, shape and runtime domain
# errors apart.
#
#   bayescomp_syntax_error       malformed Stan source (lexer or parser)
#   bayescomp_semantic_error     well formed source that breaks Stan rules
#   bayescomp_unsupported_error  valid Stan that this version does not implement
#   bayescomp_shape_error        container sizes do not agree at run time
#   bayescomp_domain_error       an argument is outside a function's domain
#   bayescomp_data_error         data supplied by the user do not match the program
#   bayescomp_reject_error       reject() or fatal_error() was called

bc_condition <- function(class, message, loc = NULL, call = NULL) {
  if (!is.null(loc) && !is.null(loc$line)) {
    message <- sprintf("%s (line %d, column %d)", message, loc$line, loc$col)
  }
  structure(
    class = c(class, "bayescomp_error", "error", "condition"),
    list(message = message, call = call, loc = loc)
  )
}

stop_syntax <- function(message, loc = NULL) stop(bc_condition("bayescomp_syntax_error", paste("Syntax error:", message), loc))
stop_semantic <- function(message, loc = NULL) stop(bc_condition("bayescomp_semantic_error", paste("Semantic error:", message), loc))
stop_unsupported <- function(message, loc = NULL) stop(bc_condition("bayescomp_unsupported_error", paste("Unsupported:", message), loc))
stop_shape <- function(message, loc = NULL) stop(bc_condition("bayescomp_shape_error", paste("Shape error:", message), loc))
stop_domain <- function(message, loc = NULL) stop(bc_condition("bayescomp_domain_error", paste("Domain error:", message), loc))
stop_data <- function(message) stop(bc_condition("bayescomp_data_error", paste("Data error:", message)))
stop_reject <- function(message, fatal = FALSE) {
  cls <- if (fatal) c("bayescomp_fatal_error", "bayescomp_reject_error") else "bayescomp_reject_error"
  stop(structure(class = c(cls, "bayescomp_error", "error", "condition"),
                 list(message = message, call = NULL)))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
