# Reference output from CmdStan 2.40 for the same statements.

cx_print <- function(decls, code) {
  utils::capture.output(invisible(bayescomp_model(paste("transformed data {", decls, code, "}"))))
}
cx_decls <- "complex z = to_complex(1, 2);
  complex_vector[2] u = [to_complex(1, 2), 3i]';
  complex_vector[2] v = [to_complex(0.5, -1), to_complex(2, 1)]';
  complex_matrix[2,2] A = [[1, 2i], [to_complex(0, 1), 3]];
  complex_matrix[2,2] B = [[to_complex(1, 1), 2], [to_complex(0, -1), 1]];"

test_that("complex arithmetic and functions print as in CmdStan 2.40", {
  ref <- c(
    'print(z, " ", u, " ", A);' = "(1,2) [(1,2),(0,3)] [[(1,0),(0,2)],[(0,1),(3,0)]]",
    'print(dot_product(u, v), " ", dot_self(u), " ", u\' * v, " ", sum(u .* v));' = "(1.5,-8) (14,0) (1.5,-8) (-0.5,6)",
    'print(z * z, " ", z / to_complex(1, 1), " ", pow(z, 2.5), " ", abs(z), " ", arg(z), " ", norm(z), " ", conj(z), " ", sqrt(z), " ", exp(z), " ", log(z));' =
      "(-3,4) (1.5,0.5) (-6.96066,2.72962) 2.23607 1.10715 5 (1,-2) (1.27202,0.786151) (-1.1312,2.47173) (0.804719,1.10715)",
    'print(A * u, " ", A\', " ", transpose(A), " ", rows_dot_product(A, A));' =
      "[(-5,2),(-2,10)] [[(1,0),(0,1)],[(0,2),(3,0)]] [[(1,0),(0,1)],[(0,2),(3,0)]] [(-3,0),(8,0)]",
    'print(fft([1, 2, 3, 4]\'), " ", inv_fft(fft([1, 2, 3, 4]\')), " ", fft2(A), " ", inv_fft2(A));' =
      "[(10,0),(-2,2),(-2,0),(-2,-2)] [(1,0),(2,0),(3,0),(4,0)] [[(4,3),(-2,-1)],[(-2,1),(4,-3)]] [[(1,0.75),(-0.5,-0.25)],[(-0.5,0.25),(1,-0.75)]]",
    'print(prod(u), " ", trace(A), " ", polar(2, 0.5), " ", proj(z), " ", get_real(u), " ", get_imag(A));' =
      "(-6,3) (4,0) (1.75517,0.958851) (1,2) [1,0] [[0,2],[1,0]]",
    'print(z == to_complex(1, 2), " ", z != 1, " ", 1 + z, " ", -z, " ", z^2, " ", u ./ v, " ", 2 * u, " ", u / 2);' =
      "1 1 (2,2) (-1,-2) (-3,4) [(-1.2,1.6),(0.6,1.2)] [(2,4),(0,6)] [(0.5,1),(0,1.5)]",
    'print(columns_dot_product(A, B), " ", rows_dot_self(A), " ", rows_dot_product(u, v), " ", columns_dot_product(u, v));' =
      "[(2,1),(3,4)] [(5,0),(10,0)] [(2.5,0),(-3,6)] [(-0.5,6)]",
    'print(u * v\', " ", A * B, " ", diag_pre_multiply(u, A), " ", diag_post_multiply(A, v), " ", dot_product(u\', v\'));' =
      "[[(2.5,0),(0,5)],[(3,1.5),(-3,6)]] [[(3,1),(2,2)],[(-1,-2),(3,2)]] [[(1,2),(-4,2)],[(-3,0),(0,9)]] [[(0.5,-1),(-2,4)],[(1,0.5),(6,3)]] (1.5,-8)",
    'print(dims({A, A}), " ", size(u), " ", num_elements(A), " ", rows(A), " ", cols(u\'), " ", dims(z));' = "[2,2,2] 2 4 2 2 []",
    'print(mdivide_right(u\', A), " ", A / B, " ", u / to_complex(1, 1), " ", to_complex(2, 1) / u\', " ", pow(u, v));' =
      "[(1.2,1.2),(0.8,0.2)] [[(-0.1,0.3),(0.2,1.4)],[(1.2,0.4),(0.6,-0.8)]] [(1.5,0.5),(1.5,1.5)] [(0.8,-0.6),(0.333333,-0.666667)] [(4.38257,-1.1244),(-0.850953,-1.6662)]",
    'print(cumulative_sum(u), " ", sum(A), " ", prod(A), " ", add_diag(A, 1i), " ", append_row(u, 1), " ", to_matrix(u, 1, 2));' =
      "[(1,2),(1,5)] (4,3) (-6,0) [[(1,1),(0,2)],[(0,1),(3,1)]] [(1,2),(0,3),(1,0)] [[(1,2),(0,3)]]",
    'print(sin(z), " ", cosh(z), " ", atanh(to_complex(0.3, 0.2)), " ", log10(z), " ", acos(z), " ", asinh(z));' =
      "(3.16578,1.9596) (-0.642148,1.06861) (0.29575,0.215474) (0.349485,0.480829) (1.14372,-1.52857) (1.46935,1.06344)",
    'print(to_complex([1, 2]\', [3, 4]\'), " ", to_complex(1, {2.0, 3.0}), " ", get_imag({to_complex(1, 2)}), " ", abs(u), " ", conj(A), " ", to_complex(), " ", to_complex(3));' =
      "[(1,3),(2,4)] [(1,2),(1,3)] [2] [2.23607,3] [[(1,-0),(0,-2)],[(0,-1),(3,-0)]] (0,0) (3,0)",
    'print(eigenvalues_sym([[2, 1i], [1i, 3]]), " ", singular_values([[1, 2], [3, 4]]));' = "[(1.38197,0),(3.61803,0)] [5.46499,0.365966]")
  for (code in names(ref)) expect_identical(cx_print(cx_decls, code), ref[[code]], label = code)
})

test_that("complex decompositions agree with CmdStan up to phase", {
  vals <- function(code) {
    out <- cx_print(cx_decls, code)
    num <- as.numeric(strsplit(trimws(gsub("[][(),]", " ", out)), " +")[[1L]])
    complex(real = num[c(TRUE, FALSE)], imaginary = num[c(FALSE, TRUE)])
  }
  ev <- vals("print(eigenvalues([[1, 2], [-3, 1]]));")
  expect_equal(sort(Im(ev)), c(-sqrt(6), sqrt(6)), tolerance = 1e-5)
  expect_equal(Re(ev), c(1, 1), tolerance = 1e-5)
  # A V = V diag(lambda) for the complex eigenvectors
  m <- matrix(c(1, -3, 2, 1), 2)
  V <- matrix(vals("print(eigenvectors([[1, 2], [-3, 1]]));"), 2, byrow = TRUE)
  expect_equal(m %*% V, V %*% diag(ev), tolerance = 1e-5)
  expect_equal(colSums(Mod(V)^2), c(1, 1), tolerance = 1e-5)
  # Schur: U T U^H reproduces the matrix, T upper triangular
  Tm <- matrix(vals("print(complex_schur_decompose_t([[1, 2], [3, 4]]));"), 2, byrow = TRUE)
  U <- matrix(vals("print(complex_schur_decompose_u([[1, 2], [3, 4]]));"), 2, byrow = TRUE)
  expect_equal(U %*% Tm %*% Conj(t(U)), matrix(c(1, 3, 2, 4), 2) + 0i, tolerance = 1e-5)
  expect_equal(Re(diag(Tm)), c(-0.372281, 5.37228), tolerance = 1e-5)
  # SVD of a complex matrix
  A <- matrix(c(1, 1i, 2i, 3), 2)
  Us <- matrix(vals("print(svd_U(A));"), 2, byrow = TRUE)
  Vs <- matrix(vals("print(svd_V(A));"), 2, byrow = TRUE)
  d <- c(3.61803, 1.38197)
  expect_equal(Us %*% diag(d) %*% Conj(t(Vs)), A, tolerance = 1e-5)
})

test_that("complex values in the model block have CmdStan gradients", {
  code <- "data { complex_vector[3] d; }
  parameters { vector[4] x; real<lower=0> s; }
  model {
    complex_vector[4] f = fft(x);
    complex z = to_complex(x[1], x[2] * s);
    complex_matrix[2,2] M = [[z, 1], [x[3] * 1i, to_complex(s, x[4])]];
    complex_vector[2] w = M * [x[1], z]';
    target += -0.5 * sum(abs(f)) - 0.1 * norm(z * z + exp(z / 3));
    target += get_real(sum(w)) + get_imag(prod(w)) * 0.01 + arg(z) + abs(pow(z, 1.5));
    target += abs(dot_product(w, d[1:2])) * 0.1 + get_real(sum(inv_fft(f) .* d[1] ./ (1 + x[1]^2)));
    target += abs(sin(z) + log(z) + sqrt(z) + tanh(z) / 4) * 0.1 + get_real(conj(z) / (z + 2));
    target += get_imag(sum((M / to_complex(2, 1))[1])) + abs(trace(M * M)) * 0.01;
    target += sum(abs(mdivide_right(w', M))) * 0.1;
    x ~ normal(0, 1);
    s ~ lognormal(0, 1);
  }"
  m <- bayescomp_model(code, data = list(d = complex(real = c(1, -0.3, 0.7), imaginary = c(0.5, 2, -1))))
  u <- bayescomp_unconstrain(m, list(x = c(0.3, -0.7, 1.1, 0.4), s = 1.3))
  g <- bayescomp_grad_log_prob(m, u, jacobian = TRUE, propto = TRUE)
  # CmdStan 2.40 log_prob method
  expect_equal(attr(g, "log_prob"), -2.7687107, tolerance = 1e-7)
  expect_equal(as.numeric(g), c(3.9453588, 0.57173266, -1.8120959, 0.76657787, 1.6538513), tolerance = 1e-7)
})

test_that("complex output follows cmdstanr names and order", {
  m <- bayescomp_model("parameters { real y; } model { y ~ normal(0, 1); }
    generated quantities { complex g = to_complex(y, 1); complex_vector[2] gv = [g, 2 * g]';
      array[2] complex ga = {g, conj(g)}; complex_matrix[2,2] gm = [[g, 1], [2, 3i]]; }")
  nm <- bayescomp_param_names(m)
  expect_identical(nm[1:7], c("y", "g[real]", "g[imag]", "gv[1,real]", "gv[1,imag]", "gv[2,real]", "gv[2,imag]"))
  expect_identical(nm[12:15], c("gm[1,1,real]", "gm[1,1,imag]", "gm[2,1,real]", "gm[2,1,imag]"))
  v <- bayescomp_constrain(m, 0.5)
  expect_identical(v$gm, matrix(c(0.5 + 1i, 2, 1, 3i), 2))
  expect_identical(v$ga, c(0.5 + 1i, 0.5 - 1i))
  expect_identical(flatten_draw(v)[1:5], c(0.5, 0.5, 1, 0.5, 1))
})

test_that("complex data, functions and loops", {
  code <- "functions { complex twice(complex z) { return 2 * z; } }
    data { array[2] complex zs; complex_matrix[2, 2] Z; }
    transformed data {
      complex_vector[2] v = [1, 2]';
      v[2] = 3i;
      v[1] += 1;
      for (z in zs) print(twice(z));
      for (w in Z) print(w);
      print(v, \" \", Z[2], \" \", Z[:, 1], \" \", Z[1, 2]);
    }"
  out <- utils::capture.output(invisible(bayescomp_model(code, data = list(zs = c(1 + 2i, complex(real = 0, imaginary = -1)),
                                                                       Z = matrix(c(1, 2i, 3, 4 + 1i), 2)))))
  expect_identical(out, c("(2,4)", "(0,-2)", "(1,0)", "(0,2)", "(3,0)", "(4,1)",
                          "[(2,0),(0,3)] [(0,2),(4,1)] [(1,0),(0,2)] (3,0)"))
  expect_error(bayescomp_model("data { complex_vector[2] v; }", data = list(v = 1:3)), class = "bayescomp_data_error")
  expect_error(bayescomp_check("parameters { real y; } model { complex z = y; target += z; }"), class = "bayescomp_semantic_error")
  expect_error(bayescomp_check("transformed data { complex_vector[2] u; print(dims(u)); }"), class = "bayescomp_semantic_error")
})
