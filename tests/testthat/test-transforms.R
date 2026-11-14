jac_logdet <- function(fx, u) {
  J <- vapply(seq_along(u), function(i) {
    e <- numeric(length(u)); e[i] <- 1e-6
    (fx(u + e) - fx(u - e)) / 2e-6
  }, numeric(length(fx(u))))
  J <- matrix(J, ncol = length(u))
  as.numeric(determinant(J, logarithm = TRUE)$modulus)
}

test_that("constraining transforms round trip and have Stan's Jacobians", {
  set.seed(42)
  # kind, eshape, args, function selecting the free coordinates of x
  cases <- list(
    list("bounds", 3L, list(lower = 1), function(x) x),
    list("bounds", 3L, list(upper = 2), function(x) x),
    list("bounds", 3L, list(lower = -1, upper = 4), function(x) x),
    list("bounds", 3L, list(offset = 1, multiplier = 3), function(x) x),
    list("ordered", 4L, list(), function(x) x),
    list("positive_ordered", 4L, list(), function(x) x),
    list("simplex", 4L, list(), function(x) x[1:3]),
    list("cholesky_factor_corr", c(3L, 3L), list(), function(x) x[lower.tri(x)]),
    list("corr_matrix", c(3L, 3L), list(), function(x) x[lower.tri(x)]),
    list("cov_matrix", c(3L, 3L), list(), function(x) x[lower.tri(x, diag = TRUE)]),
    list("cholesky_factor_cov", c(4L, 3L), list(), function(x) x[lower.tri(x, diag = TRUE)]),
    list("column_stochastic_matrix", c(3L, 2L), list(), function(x) as.vector(x[1:2, ])),
    list("row_stochastic_matrix", c(2L, 3L), list(), function(x) as.vector(t(x[, 1:2]))))
  for (cs in cases) {
    decl <- list(name = "x", kind = cs[[1L]], adims = integer(0), eshape = cs[[2L]])
    n <- free_size_variable(decl)
    u <- stats::rnorm(n, 0, 0.7)
    r <- constrain_variable(u, decl, cs[[3L]])
    x <- val(r$x)
    expect_null(check_variable(x, decl, cs[[3L]]), label = cs[[1L]])
    expect_equal(unconstrain_variable(x, decl, cs[[3L]]), u, tolerance = 1e-8, label = cs[[1L]])
    fx <- function(uu) cs[[4L]](val(constrain_variable(uu, decl, cs[[3L]])$x))
    expect_equal(val(r$lj), jac_logdet(fx, u), tolerance = 1e-5, label = paste("Jacobian", cs[[1L]]))
  }
})

test_that("unit vector, sum to zero and arrays of constrained types work", {
  decl <- list(name = "x", kind = "sum_to_zero_vector", adims = integer(0), eshape = 5L)
  u <- c(0.3, -1, 2, 0.5)
  x <- val(constrain_variable(u, decl, list())$x)
  expect_equal(sum(x), 0, tolerance = 1e-12)
  expect_equal(unconstrain_variable(x, decl, list()), u, tolerance = 1e-10)
  decl <- list(name = "x", kind = "sum_to_zero_matrix", adims = integer(0), eshape = c(3L, 4L))
  u <- stats::rnorm(6)
  x <- val(constrain_variable(u, decl, list())$x)
  expect_equal(rowSums(x), rep(0, 3), tolerance = 1e-12)
  expect_equal(colSums(x), rep(0, 4), tolerance = 1e-12)
  decl <- list(name = "x", kind = "unit_vector", adims = integer(0), eshape = 3L)
  x <- val(constrain_variable(c(1, 2, 2), decl, list())$x)
  expect_equal(sqrt(sum(x^2)), 1)
  decl <- list(name = "x", kind = "simplex", adims = c(2L, 3L), eshape = 4L)
  u <- stats::rnorm(free_size_variable(decl))
  x <- val(constrain_variable(u, decl, list())$x)
  expect_equal(dim(x), c(2L, 3L, 4L))
  expect_equal(apply(x, c(1, 2), sum), matrix(1, 2, 3))
  expect_equal(unconstrain_variable(x, decl, list()), u, tolerance = 1e-10)
})
