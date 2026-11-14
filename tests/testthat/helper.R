# shared helpers for the test suite

num_grad <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(i) {
    e <- numeric(length(x)); e[i] <- h
    (f(x + e) - f(x - e)) / (2 * h)
  }, 0)
}

# gradient of a scalar function written with differentiable operations
expect_grad_ok <- function(f, x, tol = 1e-5, label = "") {
  r <- ad_value_grad(f, x)
  ng <- num_grad(function(z) val(f(z)), x)
  expect_lt(max(abs(r$grad - ng) / pmax(1, abs(ng))), tol, label = label)
}

long_tests <- function() identical(Sys.getenv("BAYESCOMPR_LONG_TESTS"), "true")

lib_fun <- function(name, i = 1L) .lib[[name]][[i]]$impl
