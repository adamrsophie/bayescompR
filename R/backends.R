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

# capabilities this version takes from a native backend
native_capabilities <- c("log_density", "gradient")

#' Inputs for a native backend
#'
#' Collects what a native backend needs to compile a model: the IR as JSON
#' and the values of every data and transformed data variable, already
#' computed by the R implementation. Values are flattened in Stan order
#' (array elements in row major order, each element column major).
#'
#' @param model A `bayescomp_model` (any backend) or the data bound model
#'   passed to a backend's `compile` function.
#' @return A list with `ir` (JSON string), `data` (JSON string mapping each
#'   name to `dims` and `values`) and `dim` (number of unconstrained
#'   parameters).
#' @export
bayescomp_native_inputs <- function(model) {
  check_model(model)
  ir <- bayescomp_ir(model, json = TRUE)
  names_d <- c(vapply(model$compiled$data_decls, `[[`, "", "name"),
               vapply(model$compiled$td_decls, `[[`, "", "name"))
  parts <- vapply(names_d, function(nm) {
    info <- model$decls[[nm]]
    v <- val(get(paste0("V.", nm), envir = model$env))
    dims <- c(info$adims, info$eshape)
    x <- to_element_order(as.numeric(v), info$adims, info$eshape)
    paste0(json_str(nm), ":{\"dims\":", json_array(as.numeric(dims)), ",\"values\":", json_array(x), "}")
  }, "")
  list(ir = ir, data = paste0("{", paste(parts, collapse = ","), "}"), dim = model$dim)
}

#' Signal an error raised by a native backend
#'
#' Native runtimes report errors as strings with a code prefix. This turns
#' them into the condition classes the R backend uses, so samplers treat a
#' rejection the same way whichever backend raised it.
#'
#' @param message Error string such as `"E_REJECT: ..."`.
#' @return Does not return.
#' @export
bayescomp_native_error <- function(message) {
  code <- sub(":.*$", "", message)
  body <- sub("^[A-Z_]+:[[:space:]]*", "", message)
  switch(code,
    E_REJECT = stop_reject(body),
    E_FATAL = stop_reject(body, fatal = TRUE),
    E_UNSUPPORTED = stop_unsupported(body),
    E_DOMAIN = stop_domain(body),
    E_SHAPE = stop_shape(body),
    stop(sprintf("native backend error: %s", body), call. = FALSE))
}

#' Compatibility matrix
#'
#' A machine readable list of Stan 2.40 features and their status in this
#' version: `"implemented"`, `"partial"`, `"not_implemented"` or
#' `"not_applicable"`, with the test files that cover each feature.
#'
#' @return A data frame with columns `area`, `feature`, `status`, `tests`
#'   and `notes`.
#' @export
#' @examples
#' cm <- bayescomp_compatibility()
#' table(cm$status)
bayescomp_compatibility <- function() {
  utils::read.csv(system.file("compatibility.csv", package = "bayescompR"), stringsAsFactors = FALSE)
}
