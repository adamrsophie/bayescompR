# Execution backends.
#
# The parser and semantic analyzer produce a typed IR that does not depend on
# any backend. A backend turns the IR into the functions the inference
# algorithms call. The "r" backend (this package) is always available. Native
# backends live in separate packages and register themselves when loaded.

.backends <- new.env(parent = emptyenv())

#' Register an execution backend
#'
#' Used by backend packages such as bayescompRCranelift and bayescompRLLVM.
#'
#' @param name Backend name.
#' @param capabilities Character vector of supported capabilities, from
#'   `"log_density"`, `"gradient"`, `"constrain"`, `"generated_quantities"`,
#'   `"rng"`.
#' @param ir_versions IR schema versions the backend accepts.
#' @param compile Function taking an IR (see [bayescomp_ir()]) and the data
#'   bound model, and returning a list with function `log_density_gradient(u,
#'   jacobian, propto)` returning `list(value, grad)`.
#' @param description One line description.
#' @return Invisibly, the backend record.
#' @export
bayescomp_register_backend <- function(name, capabilities, ir_versions, compile, description = "") {
  stopifnot(is.character(name), length(name) == 1L, is.function(compile))
  rec <- list(name = name, capabilities = capabilities, ir_versions = ir_versions,
              compile = compile, description = description)
  assign(name, rec, envir = .backends)
  invisible(rec)
}

#' List the available execution backends
#' @return A data frame with one row per backend.
#' @export
bayescomp_backends <- function() {
  nms <- sort(ls(.backends))
  data.frame(
    name = nms,
    capabilities = vapply(nms, function(n) paste(.backends[[n]]$capabilities, collapse = ", "), ""),
    ir_versions = vapply(nms, function(n) paste(.backends[[n]]$ir_versions, collapse = ", "), ""),
    description = vapply(nms, function(n) .backends[[n]]$description, ""),
    row.names = NULL, stringsAsFactors = FALSE)
}

backend_get <- function(name) {
  if (!is.character(name) || length(name) != 1L) stop("'backend' must be a single string", call. = FALSE)
  rec <- .backends[[name]]
  if (is.null(rec)) {
    stop(sprintf("backend '%s' is not available; installed backends: %s", name,
                 paste(ls(.backends), collapse = ", ")), call. = FALSE)
  }
  if (!(ir_version %in% rec$ir_versions))
    stop(sprintf("backend '%s' does not support IR version %s", name, ir_version), call. = FALSE)
  rec
}

bayescomp_register_backend("r", c("log_density", "gradient", "constrain", "generated_quantities", "rng"),
                           ir_versions = "1.0.0", compile = function(ir, model) NULL,
                           description = "Pure R closures on a reverse mode AD tape")

ir_version <- "1.0.0"
