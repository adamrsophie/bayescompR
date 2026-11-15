# Native backends downloaded as prebuilt shared libraries.
#
# The libraries are built from github.com/adamrsophie/bayescomp-native-rs by
# its release workflow, one per backend, operating system and architecture.
# They do not link R: their entry points follow the calling convention of
# .C(), so this package calls them without compiled code of its own. Each
# download is checked against its checksum before it is loaded: the one
# recorded below for the assets of the release this version was tested with,
# otherwise the one in the release's SHA256SUMS or MD5SUMS file.

native_release <- "v0.1.0"
native_abi <- 2L
native_engines <- c(interpreter = 0L, cranelift = 1L, llvm = 2L)
native_state <- new.env(parent = emptyenv())

# sha256 and md5 checksums of the release assets (md5 for R < 4.5, which has
# no sha256sum)
native_assets <- list()

native_root <- function() file.path(tools::R_user_dir("bayescompR", "data"), "backends")
native_dir <- function() file.path(native_root(), native_release)

native_sha256 <- function() {
  ns <- asNamespace("tools")
  if (exists("sha256sum", envir = ns)) get("sha256sum", envir = ns) else NULL
}

native_platform <- function() {
  os <- tolower(Sys.info()[["sysname"]])
  arch <- R.version$arch
  if (arch %in% c("arm64", "aarch64")) arch <- "aarch64"
  if (arch %in% c("x86_64", "x64", "amd64")) arch <- "x86_64"
  paste(os, arch, sep = "-")
}

native_asset <- function(backend) sprintf("bcn-%s-%s%s", backend, native_platform(), .Platform$dynlib.ext)

native_file <- function(backend) file.path(native_dir(), paste0("bcn-", backend, .Platform$dynlib.ext))

native_hash_ok <- function(file, sums) {
  sha <- native_sha256()
  if (!is.null(sha)) {
    identical(unname(sha(file)), sums[["sha256"]])
  } else {
    identical(unname(tools::md5sum(file)), sums[["md5"]])
  }
}

#' Install a native backend
#'
#' Downloads the prebuilt shared library of a native backend from the
#' releases of <https://github.com/adamrsophie/bayescomp-native-rs>, checks
#' its checksum, stores it in the user data directory
#' (`tools::R_user_dir("bayescompR", "data")`) and registers the backend, so
#' models can then be built with `backend = "cranelift"` or `"llvm"`. The
#' backend is registered again in every later session. Nothing is
#' downloaded unless this function is called.
#'
#' The native backends evaluate the log density and its gradient in compiled
#' code; constraining, generated quantities and random numbers stay in R, so
#' sampling needs `allow_fallback = TRUE`. Programs with features the native
#' lowering does not cover (user defined functions, truncation and a few
#' others) are refused with `bayescomp_unsupported_error` when the model is
#' built. `options(bayescompR.native_engine = "interpreter")` uses the
#' reference interpreter of the native runtime instead of the JIT compiler.
#' An evaluation in native code cannot be interrupted, so a Stan `while`
#' loop that never ends blocks the session.
#'
#' The checksum is SHA-256 on R 4.5 and newer and MD5 on older versions of
#' R; MD5 detects a damaged download but is not a security check. The
#' Cranelift libraries need glibc 2.28 or newer on Linux and macOS 11 or
#' newer. The LLVM libraries link LLVM statically; on Linux they need glibc
#' 2.35 or newer, libstdc++, zlib, zstd and libxml2 (`libxml2.so.2`), and on
#' macOS they are built and tested on macOS 15 (Apple silicon). Installing
#' removes the libraries of other releases from the user data directory.
#'
#' @param backend `"cranelift"` or `"llvm"`.
#' @param path Install a library built locally (for example with `cargo
#'   build -p bcn-ffi --features cranelift` in bayescomp-native-rs) instead
#'   of downloading one; its checksum is not checked.
#' @param force Download again even when the backend is installed.
#' @param quiet Suppress download progress.
#' @return Invisibly, the path of the installed library.
#' @seealso [bayescomp_remove_backend()], [bayescomp_validate_backend()]
#' @examples
#' \dontrun{
#' bayescomp_install_backend("cranelift")
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }",
#'                      backend = "cranelift", allow_fallback = TRUE)
#' bayescomp_validate_backend(m)
#' }
#' @export
bayescomp_install_backend <- function(backend = c("cranelift", "llvm"), path = NULL, force = FALSE, quiet = FALSE) {
  backend <- match.arg(backend)
  dest <- native_file(backend)
  if (file.exists(dest) && !force && is.null(path)) {
    native_register(backend)
    return(invisible(dest))
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(fileext = .Platform$dynlib.ext)
  on.exit(unlink(tmp), add = TRUE)
  if (!is.null(path)) {
    if (!file.exists(path)) stop(sprintf("file '%s' does not exist", path), call. = FALSE)
    file.copy(path, tmp, overwrite = TRUE)
  } else {
    asset <- native_asset(backend)
    sums <- native_assets[[asset]] %||% native_release_sums(asset)
    if (is.null(sums))
      stop(sprintf("no prebuilt %s backend for %s in release %s of bayescomp-native-rs; build bcn-ffi and pass 'path'",
                   backend, native_platform(), native_release), call. = FALSE)
    url <- native_url(asset)
    ok <- native_download(url, tmp, quiet)
    if (!isTRUE(ok)) stop(sprintf("could not download %s: %s", url, ok), call. = FALSE)
    if (!native_hash_ok(tmp, sums))
      stop(sprintf("the checksum of %s does not match; the file was not installed", asset), call. = FALSE)
  }
  # check the new library before it replaces a working one
  # under a unique name in the temporary directory: a library is never
  # unloaded, since the Rust runtime registers thread local destructors and
  # unloading it can crash R (dyn.load() of a loaded path unloads it first).
  # So a library loaded in this session stays in use until R restarts.
  staged <- tempfile("bcn-staged-", fileext = .Platform$dynlib.ext)
  if (!file.copy(tmp, staged, overwrite = TRUE)) stop(sprintf("could not write '%s'", staged), call. = FALSE)
  native_open(staged, backend)
  if (!is.null(native_state[[backend]]))
    message(sprintf("The %s backend was already loaded in this session; restart R to use the new library.", backend))
  # copied next to the destination and renamed onto it: the rename is atomic
  # and gives the file a new inode, so a copy mapped by another session is
  # never rewritten
  part <- tempfile("bcn-part-", tmpdir = dirname(dest), fileext = .Platform$dynlib.ext)
  if (!file.copy(tmp, part)) stop(sprintf("could not write '%s'", part), call. = FALSE)
  if (!file.rename(part, dest)) {
    unlink(part)
    stop(sprintf("could not replace '%s' (it may be in use); restart R and try again", dest), call. = FALSE)
  }
  native_prune()
  native_load(backend)
  native_register(backend)
  invisible(dest)
}

native_url <- function(file)
  sprintf("https://github.com/adamrsophie/bayescomp-native-rs/releases/download/%s/%s", native_release, file)

# TRUE, or the reason the download failed
native_download <- function(url, dest, quiet = TRUE) {
  old <- options(timeout = max(600, getOption("timeout")))
  on.exit(options(old))
  tryCatch(if (utils::download.file(url, dest, mode = "wb", quiet = quiet) == 0L) TRUE else "download.file() failed",
           error = function(e) conditionMessage(e), warning = function(w) conditionMessage(w))
}

# the checksums published with the release; NULL when the release has no such
# asset, an error when the checksums cannot be downloaded
native_release_sums <- function(asset) {
  kind <- if (is.null(native_sha256())) "md5" else "sha256"
  f <- tempfile()
  on.exit(unlink(f))
  url <- native_url(if (kind == "sha256") "SHA256SUMS" else "MD5SUMS")
  ok <- native_download(url, f)
  if (!isTRUE(ok)) stop(sprintf("could not download %s: %s", url, ok), call. = FALSE)
  x <- strsplit(trimws(readLines(f, warn = FALSE)), "[[:space:]]+")
  hit <- Filter(function(l) length(l) == 2L && sub("^\\*", "", l[2L]) == asset, x)
  if (!length(hit)) return(NULL)
  stats::setNames(list(hit[[1L]][1L]), kind)
}

#' Remove an installed native backend
#'
#' Unregisters a backend installed with [bayescomp_install_backend()] and
#' deletes its library. A library already loaded in the session stays in
#' memory until R restarts (on Windows the file can only be deleted then).
#' @param backend `"cranelift"` or `"llvm"`.
#' @return Invisibly, `TRUE` when a library was removed.
#' @examples
#' \dontrun{
#' bayescomp_remove_backend("cranelift")
#' }
#' @export
bayescomp_remove_backend <- function(backend = c("cranelift", "llvm")) {
  backend <- match.arg(backend)
  rec <- .backends[[backend]]
  if (!is.null(rec) && isTRUE(rec$downloaded)) rm(list = backend, envir = .backends)
  f <- native_file(backend)
  done <- file.exists(f) && unlink(f) == 0L
  # empty directories go too
  for (d in c(native_dir(), native_root())) if (dir.exists(d) && !length(list.files(d, all.files = TRUE, no.. = TRUE))) unlink(d, recursive = TRUE)
  invisible(done)
}

# the libraries of other releases
native_prune <- function() {
  old <- setdiff(list.dirs(native_root(), full.names = FALSE, recursive = FALSE), native_release)
  if (length(old)) unlink(file.path(native_root(), old), recursive = TRUE)
}

# load a library and check its interface
native_open <- function(f, backend) {
  dll <- tryCatch(dyn.load(f, local = TRUE, now = TRUE), error = function(e)
    stop(sprintf("'%s' is not a loadable %s backend library: %s", f, backend, conditionMessage(e)), call. = FALSE))
  sym <- function(n) tryCatch(getNativeSymbolInfo(n, dll), error = function(e)
    stop(sprintf("'%s' is not a %s backend library (no %s)", f, backend, n), call. = FALSE))
  lib <- list(path = f, info = sym("bcn_rc_info"), compile = sym("bcn_rc_compile"),
              eval = sym("bcn_rc_eval"), free = sym("bcn_rc_free"))
  info_sym <- lib$info
  inf <- .C(info_sym, abi = integer(1), engines = integer(1), generation = integer(1))
  if (inf$abi != native_abi)
    stop(sprintf("the %s backend library has interface version %d, this bayescompR needs %d; run bayescomp_install_backend(\"%s\", force = TRUE)",
                 backend, inf$abi, native_abi, backend), call. = FALSE)
  need <- c(cranelift = 2L, llvm = 4L)[[backend]]
  if (bitwAnd(inf$engines, need) == 0L)
    stop(sprintf("'%s' was built without the %s engine", f, backend), call. = FALSE)
  lib$engines <- inf$engines
  lib$generation <- inf$generation
  lib
}

native_load <- function(backend) {
  lib <- native_state[[backend]]
  if (!is.null(lib)) return(lib)
  f <- native_file(backend)
  if (!file.exists(f))
    stop(sprintf("the %s backend is not installed; run bayescomp_install_backend(\"%s\")", backend, backend), call. = FALSE)
  lib <- native_open(f, backend)
  assign(backend, lib, envir = native_state)
  lib
}

native_register <- function(backend) {
  rec <- bayescomp_register_backend(
    backend, c("log_density", "gradient"), ir_versions = ir_version, compile = native_compile(backend),
    description = sprintf("%s machine code from the downloaded native library (log density and gradient)",
                          c(cranelift = "Cranelift", llvm = "LLVM")[[backend]]))
  .backends[[backend]]$downloaded <- TRUE
  invisible(rec)
}

native_compile <- function(backend) {
  force(backend)
  function(ir, model) {
    inp <- bayescomp_native_inputs(model)
    engine <- if (identical(getOption("bayescompR.native_engine", "jit"), "interpreter")) "interpreter" else backend
    code <- native_engines[[engine]]
    errlen <- 4096L
    # a model is named by an id in the library's table and the generation of
    # the loaded library; it is valid only in the process that compiled it, so
    # another process (a socket cluster worker, a forked chain, a saved model)
    # or a library loaded again compiles its own
    st <- new.env(parent = emptyenv())
    st$id <- 0L
    st$gen <- 0L
    st$pid <- -1L
    reg.finalizer(st, function(e) {
      lib <- native_state[[backend]]
      if (e$id > 0L && e$pid == Sys.getpid() && !is.null(lib) && lib$generation == e$gen) {
        free_sym <- lib$free
        .C(free_sym, as.integer(e$gen), as.integer(e$id))
      }
    }, onexit = TRUE)
    build <- function() {
      lib <- native_load(backend)
      compile_sym <- lib$compile
      r <- .C(compile_sym, enc2utf8(inp$ir), enc2utf8(inp$data), as.integer(inp$dim), code, id = integer(1),
              err = strrep(" ", errlen - 1L), errlen)
      if (r$id == 0L) bayescomp_native_error(native_utf8(r$err))
      st$id <- r$id
      st$gen <- lib$generation
      st$pid <- Sys.getpid()
    }
    build()
    eval1 <- function(u, jacobian, propto) {
      lib <- native_load(backend)
      eval_sym <- lib$eval
      .C(eval_sym, as.integer(st$gen), as.integer(st$id), as.double(u), length(u), as.integer(isTRUE(jacobian)),
         as.integer(isTRUE(propto)), lp = double(1), grad = double(length(u)), status = integer(1),
         err = strrep(" ", errlen - 1L), errlen, NAOK = TRUE)
    }
    ldg <- function(u, jacobian, propto) {
      lib <- native_load(backend)
      if (st$pid != Sys.getpid() || st$gen != lib$generation) build()
      r <- eval1(u, jacobian, propto)
      if (r$status == 2L) {
        build()
        r <- eval1(u, jacobian, propto)
      }
      if (r$status != 0L) bayescomp_native_error(native_utf8(r$err))
      list(value = r$lp, grad = r$grad)
    }
    list(log_density_gradient = ldg)
  }
}

native_utf8 <- function(x) {
  Encoding(x) <- "UTF-8"
  x
}

.onLoad <- function(libname, pkgname) {
  for (b in c("cranelift", "llvm")) {
    if (is.null(.backends[[b]]) && file.exists(native_file(b))) try(native_register(b), silent = TRUE)
  }
}
