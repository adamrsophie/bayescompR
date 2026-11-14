test_that("elementwise math functions have correct gradients", {
  pts <- list(
    default = c(0.3, 1.7, -0.8),
    pos = c(0.3, 1.7, 2.5),
    unit = c(0.15, 0.5, 0.85),
    ge1 = c(1.3, 2, 4))
  domain <- list(sqrt = "pos", log = "pos", log2 = "pos", log10 = "pos", log1m = "unit", logit = "unit",
                 inv_sqrt = "pos", acosh = "ge1", asin = "unit", acos = "unit", atanh = "unit",
                 lgamma = "pos", tgamma = "pos", digamma = "pos", trigamma = "pos", inv_Phi = "unit",
                 std_normal_qf = "unit", log1m_exp = "default", lambert_w0 = "pos", cbrt = "pos")
  skip_fns <- c("floor", "ceil", "round", "trunc", "step", "fabs", "lambert_wm1", "std_normal_log_qf")
  for (nm in setdiff(names(unary_defs), skip_fns)) {
    x <- pts[[domain[[nm]] %||% "default"]]
    if (nm == "log1m_exp") x <- -abs(x) - 0.1
    if (nm == "trigamma") next
    f <- unary_defs[[nm]]
    expect_grad_ok(function(z) ad_sum(f(z)), x, tol = 1e-4, label = nm)
  }
})

test_that("arithmetic, broadcasting and structural operations differentiate", {
  expect_grad_ok(function(x) ad_sum(ad_gather(x, 1:2) * ad_gather(x, 3:4) / (1 + ad_gather(x, 5L))), c(1, 2, 3, 4, 0.5))
  expect_grad_ok(function(x) ad_sum(x^2.5 + 3^x), c(0.5, 1.2))
  expect_grad_ok(function(x) ad_log_sum_exp(x) + ad_log_sum_exp2(ad_select(x, 1L), ad_select(x, 3L)), c(0.1, -2, 3))
  expect_grad_ok(function(x) ad_sum(ad_cumsum(x) * c(1, 2, 3)), c(0.1, -2, 3))
  expect_grad_ok(function(x) ad_prod(x), c(0.5, -2, 3))
  expect_grad_ok(function(x) {
    A <- ad_reshape(x, c(2L, 2L))
    ad_sum(ad_matmul(A, ad_transpose(A)) * matrix(1:4, 2))
  }, c(1, 0.5, -0.3, 2))
  expect_grad_ok(function(x) ad_sum(ad_scatter(rep(0, 4), c(2L, 4L), x) * (1:4)), c(1.5, -2))
})

test_that("linear algebra gradients match finite differences", {
  spd <- function(x) {
    L <- ad_reshape(x, c(3L, 3L))
    S <- ad_matmul(L, ad_transpose(L))
    S + diag(3)
  }
  x0 <- c(1, 0.2, -0.3, 0.1, 1.1, 0.4, -0.2, 0.3, 0.9)
  expect_grad_ok(function(x) ad_log_det(spd(x)), x0)
  expect_grad_ok(function(x) ad_sum(ad_cholesky(spd(x)) * matrix(1:9, 3)), x0)
  expect_grad_ok(function(x) ad_sum(lib_fun("inverse")(spd(x)) * matrix(1:9, 3)), x0)
  expect_grad_ok(function(x) ad_sum(ad_solve(spd(x), c(1, 2, 3), "mdivide_left")), x0)
  expect_grad_ok(function(x) ad_sum(ad_eigen_sym(spd(x), "values") * (1:3)), x0, tol = 1e-4)
  expect_grad_ok(function(x) ad_sum(ad_matrix_exp(0.3 * ad_reshape(x, c(3L, 3L))) * matrix(1:9, 3)), x0, tol = 1e-4)
})
