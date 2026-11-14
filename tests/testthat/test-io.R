test_that("CmdStan JSON data files are read as CmdStan reads them", {
  code <- "
data {
  int N; array[N] int y; matrix[2, 3] A; array[2, 3, 4] int ar;
  complex z; complex_vector[2] cv; complex_matrix[2, 2] cm;
  tuple(real, array[2] int) t; array[2] tuple(int, vector[2]) at;
  real pinf; real ninf; real nn; array[0] real e; array[2, 0] real e2; vector[1] v1;
}
generated quantities {
  real a12 = A[1, 2]; int ar234 = ar[2, 3, 4]; real cm21 = get_real(cm[2, 1]); real cm12 = get_imag(cm[1, 2]);
  real zi = get_imag(z); real cv2 = get_imag(cv[2]); real t1 = t.1; int at21 = at[2].1; real at22 = at[2].2[2];
  real p = pinf; real n = ninf; int isn = is_nan(nn); real v = v1[1];
}"
  f <- tempfile(fileext = ".json")
  writeLines('{ "N": 3, "y": [0, 1, 0], "A": [[1, 2.7, -9.8], [4.2, 1.8, -7.3]],
    "ar": [[[0,1,2,3], [4,5,6,7], [8,9,10,11]], [[12,13,14,15], [16,17,18,19], [20,21,22,23]]],
    "z": [1, -2.3], "cv": [[1, 2], [3, 4]], "cm": [[[1, 2], [3, -4.2]], [[3.1, 0.5], [0, -5]]],
    "t": {"1": 1.4, "2": [1, 2]}, "at": [{"1": 1, "2": [1, 2]}, {"2": [3, 4], "1": 5}],
    "pinf": Infinity, "ninf": "-Infinity", "nn": NaN, "e": [], "e2": [[], []], "v1": [7] }', f)
  m <- bayescomp_model(code, data = f)
  v <- unlist(bayescomp_constrain(m, numeric(0)))
  expect_equal(unname(v), c(2.7, 23, 3.1, -4.2, -2.3, 4, 1.4, 5, 4, Inf, -Inf, 1, 7))
  # writing and reading back gives the same data
  g <- tempfile(fileext = ".json")
  bayescomp_write_data(m$data, g)
  expect_identical(bayescomp_constrain(bayescomp_model(code, data = g), numeric(0)), bayescomp_constrain(m, numeric(0)))
  expect_match(paste(readLines(g), collapse = ""), '"t": \\{"1": 1.4, "2": \\[1, 2\\]\\}')
  expect_error(bayescomp_read_data(tempfile()), "does not exist")
})

test_that("R dump files and init files are read", {
  f <- tempfile()
  writeLines(c("N <- 2L", "y <- structure(c(1, 2, 3, 4, 5, 6), .Dim = c(2L, 3L))", "s <- Inf"), f)
  d <- bayescomp_read_data(f)
  expect_equal(d$y, matrix(1:6, 2))
  expect_identical(d$s, Inf)
  writeLines("x <- system('echo hi')", f)
  expect_error(bayescomp_read_data(f), class = "bayescomp_data_error")
  m <- bayescomp_model("parameters { real<lower=0> s; real mu; } model { s ~ exponential(1); mu ~ normal(0, 1); }")
  i <- tempfile(fileext = ".json")
  bayescomp_write_data(list(s = 0.25, mu = 0.5), i)
  q <- find_init(m, i, 1L)
  expect_equal(q, c(log(0.25), 0.5))
})

test_that("library calls follow stanc's signatures", {
  expect_error(bayescomp_check("data { array[2, 2] int a; } generated quantities { int s = sum(a); }"),
               "no matching signature", class = "bayescomp_semantic_error")
  expect_error(bayescomp_check("data { array[2] int a; vector[2] v; } generated quantities { vector[2] p = pow(a, v); }"),
               class = "bayescomp_semantic_error")
  expect_error(bayescomp_check("data { real a; } generated quantities { int b = logical_negation(a); }"),
               class = "bayescomp_semantic_error")
  expect_silent(bayescomp_check("data { array[2] int a; vector[2] v; } generated quantities { real s = sum(a) + max(a) + sum(v); real l = binomial_logit_lpmf(1 | 2, 0.1); }"))
})

test_that("log_prob takes constrained values, files and matrices of points", {
  m <- bayescomp_model("parameters { real<lower=0> s; } model { s ~ exponential(2); }")
  f <- tempfile(fileext = ".json")
  bayescomp_write_data(list(s = 0.5), f)
  lp <- bayescomp_log_prob(m, log(0.5))
  expect_equal(bayescomp_log_prob(m, list(s = 0.5)), lp)
  expect_equal(bayescomp_log_prob(m, f), lp)
  U <- matrix(log(c(0.5, 1, 2)), 3L)
  expect_equal(bayescomp_log_prob(m, U), vapply(U, function(u) bayescomp_log_prob(m, u), 0))
  g <- bayescomp_grad_log_prob(m, U)
  expect_identical(dim(g), c(3L, 1L))
  expect_equal(g[, 1L], 1 - 2 * exp(U[, 1L]))
})

test_that("the gradient test agrees with finite differences", {
  m <- bayescomp_model("parameters { real<lower=0> s; vector[2] b; } model { s ~ gamma(2, 1); b ~ normal(s, 2); }")
  g <- bayescomp_check_gradients(m, init = list(s = 1.5, b = c(0.2, -0.3)))
  expect_true(attr(g, "passed"))
  expect_identical(g$param, 0:2)
  expect_equal(g$value[1L], log(1.5))
  expect_equal(attr(g, "log_prob"), bayescomp_log_prob(m, g$value))
})

test_that("print formats numbers as CmdStan does", {
  # reference output from CmdStan 2.40
  expect_output(bayescomp_model("transformed data { print(1, \" \", 1.5, \" \", 123456789.123, \" \", 1.0/3, \" \", 100000.0, \" \", 1234567.0, \" \", [1, 2.5], \" \", {1.5, 2.5}, \" \", to_complex(1.5, -2), \" \", not_a_number(), \" \", negative_infinity()); }"),
                "1 1.5 1.23457e+08 0.333333 100000 1.23457e+06 [1,2.5] [1.5,2.5] (1.5,-2) nan -inf", fixed = TRUE)
})
