test_that("native backend assets are named by platform and checked by checksum", {
  expect_match(native_platform(), "^[a-z]+-[a-z0-9_]+$")
  expect_match(native_asset("cranelift"), paste0("^bcn-cranelift-", native_platform(), "\\", .Platform$dynlib.ext, "$"))
  f <- tempfile()
  writeBin(charToRaw("abc\n"), f)
  skip_if(is.null(native_sha256()))
  expect_true(native_hash_ok(f, "edeaaff3f1774ad2888673770c6d64097e391bc362d7d6fb34982ddf0efd18cb"))
  expect_false(native_hash_ok(f, "0"))
})

test_that("a missing native backend gives an informative error", {
  old <- Sys.getenv("R_USER_DATA_DIR", NA)
  Sys.setenv(R_USER_DATA_DIR = tempfile())
  on.exit(if (is.na(old)) Sys.unsetenv("R_USER_DATA_DIR") else Sys.setenv(R_USER_DATA_DIR = old), add = TRUE)
  expect_error(native_load("llvm"), "bayescomp_install_backend")
  expect_error(bayescomp_install_backend("cranelift", path = tempfile()), "does not exist")
  local_mocked_bindings(native_platform = function() "plan9-mips", native_release_sums = function(asset) NULL)
  if (is.null(native_sha256())) expect_error(bayescomp_install_backend("cranelift"), "R 4.5")
  else expect_error(bayescomp_install_backend("cranelift"), "no prebuilt cranelift backend for plan9-mips")
})

test_that("the checksums published with a release are read", {
  local_mocked_bindings(native_download = function(url, dest, quiet = TRUE) {
    writeLines(c("abc123  bcn-cranelift-linux-x86_64.so", "def456 *bcn-llvm-linux-x86_64.so"), dest)
    TRUE
  })
  s <- native_release_sums("bcn-llvm-linux-x86_64.so")
  expect_equal(s, "def456")
  expect_null(native_release_sums("bcn-llvm-windows-x86_64.dll"))
  local_mocked_bindings(native_download = function(url, dest, quiet = TRUE) "cannot open URL")
  expect_error(native_release_sums("bcn-llvm-linux-x86_64.so"), "could not download .*cannot open URL")
})

test_that("installing removes the libraries of other releases", {
  old <- Sys.getenv("R_USER_DATA_DIR", NA)
  Sys.setenv(R_USER_DATA_DIR = tempfile())
  on.exit(if (is.na(old)) Sys.unsetenv("R_USER_DATA_DIR") else Sys.setenv(R_USER_DATA_DIR = old), add = TRUE)
  dir.create(file.path(native_root(), "v0.0.9"), recursive = TRUE)
  dir.create(native_dir(), recursive = TRUE)
  native_prune()
  expect_equal(list.files(native_root()), native_release)
  expect_false(bayescomp_remove_backend("llvm"))
  expect_false(dir.exists(native_root()))
})

test_that("a locally built native library matches the R backend", {
  lib <- Sys.getenv("BAYESCOMPR_NATIVE_LIB")
  skip_if(!nzchar(lib) || !file.exists(lib))
  old <- Sys.getenv("R_USER_DATA_DIR", NA)
  Sys.setenv(R_USER_DATA_DIR = tempfile())
  on.exit(if (is.na(old)) Sys.unsetenv("R_USER_DATA_DIR") else Sys.setenv(R_USER_DATA_DIR = old), add = TRUE)
  on.exit(bayescomp_remove_backend("cranelift"), add = TRUE)
  bayescomp_install_backend("cranelift", path = lib)
  expect_true("cranelift" %in% bayescomp_backends()$name)
  code <- "data { int N; vector[N] y; } parameters { real mu; real<lower=0> s; } model { y ~ normal(mu, s); }"
  d <- list(N = 5, y = c(0.1, 0.5, -0.3, 1.2, 0.8))
  m <- bayescomp_model(code, data = d, backend = "cranelift", allow_fallback = TRUE)
  r <- bayescomp_model(code, data = d)
  expect_equal(bayescomp_grad_log_prob(m, c(0.2, 0.1)), bayescomp_grad_log_prob(r, c(0.2, 0.1)), tolerance = 1e-10)
  expect_error(bayescomp_log_prob(bayescomp_model("parameters { real<lower=0> s; } model { target += normal_lpdf(1 | 0, -s); }",
                                                  backend = "cranelift", allow_fallback = TRUE), 0),
               class = "bayescomp_domain_error")
  # a bad file does not replace the working library
  bad <- tempfile(fileext = .Platform$dynlib.ext)
  writeLines("not a library", bad)
  expect_error(bayescomp_install_backend("cranelift", path = bad), "not a loadable")
  expect_equal(list.files(dirname(native_file("cranelift"))), basename(native_file("cranelift")))
  expect_equal(bayescomp_log_prob(m, c(0.2, 0.1)), bayescomp_log_prob(r, c(0.2, 0.1)), tolerance = 1e-10)
  # forked chains compile their own copies
  skip_on_os("windows")
  fit <- bayescomp_sample(m, chains = 2, cores = 2, iter_warmup = 50, iter_sampling = 50, seed = 1, refresh = 0)
  expect_equal(dim(fit$draws)[2], 2L)
})

test_that("strings in the IR are valid JSON with control characters", {
  x <- json_str(c(paste0("a", intToUtf8(1), "b\tc\r"), "q\"\\"))
  expect_equal(x, c("\"a\\u0001b\\tc\\r\"", "\"q\\\"\\\\\""))
})
