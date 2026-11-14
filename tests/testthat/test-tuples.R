# Reference output from CmdStan 2.40 for the same programs.

test_that("tuple expressions, projection and unpacking print as in CmdStan", {
  code <- "
  functions {
    tuple(real, int) f(tuple(real, int) x) { return (x.1 * 2, x.2 + 1); }
  }
  data { tuple(real, array[2] int) td; array[2] tuple(int, vector[2]) atd; }
  transformed data {
    tuple(real, int) t = (1.5, 2);
    tuple(real, tuple(int, vector[2])) nt = (0.5, (3, [1, 2]'));
    real a; int b;
    (a, b) = f(t);
    t.1 = 4;
    print(t, \" \", nt, \" \", a, \" \", b, \" \", nt.2.2[1], \" \", td, \" \", atd[2].2);
    tuple(matrix[2,2], vector[2]) es = eigendecompose_sym([[2, 1], [1, 3]]);
    print(es.2, \" \", qr([[1,2],[3,4]]).1, \" \", csr_extract([[1, 0], [0, 2]]));
  }"
  out <- utils::capture.output(invisible(bayescomp_model(code, data = list(
    td = list(2.5, c(3L, 4L)), atd = list(list(1L, c(0.1, 0.2)), list(2L, c(0.3, 0.4)))))))
  expect_identical(out, c("(4,2) (0.5,(3,[1,2])) 3 3 1 (2.5,[3,4]) [0.3,0.4]",
                          "[1.38197,3.61803] [[0.316228,0.948683],[0.948683,-0.316228]] ([1,2],[1,2],[1,2,3])"))
})

test_that("tuple parameters match CmdStan log_prob", {
  code <- "
  functions { real f(tuple(real, vector) t) { return t.1 * sum(t.2); } }
  parameters {
    tuple(real<lower=0>, simplex[3]) p;
    array[2] tuple(real, tuple(real<lower=-1, upper=1>, vector<lower=0>[2])) q;
  }
  model {
    target += f(p) + normal_lpdf(q[1].1 | 0, 1) + q[2].2.1 * 2 + sum(log(q[1].2.2)) - dot_self(q[2].2.2);
    p.1 ~ exponential(1);
    q[2].1 ~ normal(1, 2);
  }
  generated quantities {
    array[2, 2] tuple(int, real) g2;
    for (i in 1:2) for (j in 1:2) g2[i, j] = (i * 10 + j, i + 0.5 * j);
    tuple(real, vector[2]) h = (1.5, [1, 2]');
  }"
  m <- bayescomp_model(code)
  pars <- list(p = list(1.3, c(0.2, 0.3, 0.5)), q = list(list(0.4, list(0.3, c(1.2, 0.7))), list(-0.6, list(-0.5, c(0.4, 2.1)))))
  u <- bayescomp_unconstrain(m, pars)
  g <- bayescomp_grad_log_prob(m, u, jacobian = TRUE, propto = TRUE)
  expect_equal(attr(g, "log_prob"), -11.875173, tolerance = 1e-7)
  expect_equal(as.numeric(g), c(1, 0.21213203, 0.61237244, -0.4, -0.3, 2, 2, 0.4, 1.25, 0.68, -7.82), tolerance = 1e-7)
  v <- bayescomp_constrain(m, u)
  expect_equal(v$q, pars$q)
  expect_identical(v$g2[[2, 1]], list(21L, 2.5))
  # cmdstanr names; arrays of tuples are written in column major order
  nm <- bayescomp_param_names(m)
  expect_identical(nm[1:6], c("p:1", "p:2[1]", "p:2[2]", "p:2[3]", "q[1:1]", "q[1:2:1]"))
  expect_identical(nm[13:16], c("g2[1,1:1]", "g2[1,1:2]", "g2[2,1:1]", "g2[2,1:2]"))
  expect_identical(flatten_draw(v)[13:20], c(11, 1.5, 21, 2.5, 12, 2, 22, 3))
})

test_that("tuple semantics follow stanc", {
  expect_error(bayescomp_check("transformed data { tuple(real, int) t = (1.5, 2.5); }"), class = "bayescomp_semantic_error")
  expect_error(bayescomp_check("transformed data { tuple(real, int) t = (1.5, 2); print(t.3); }"), class = "bayescomp_semantic_error")
  expect_error(bayescomp_check("model { tuple(real<lower=0>, int) t; }"), class = "bayescomp_semantic_error")
  # promotion inside tuples, nested unpacking, loops over arrays of tuples
  out <- utils::capture.output(invisible(bayescomp_model("transformed data {
    tuple(real, complex) t = (1, 2);
    real a; real b; complex c;
    ((a, b), c) = ((1, 2.5), 3i);
    array[2] tuple(int, real) at = {(1, 0.5), (2, 1.5)};
    for (x in at) print(x.1 + x.2);
    print(t, \" \", a + b, \" \", c);
  }")))
  expect_identical(out, c("1.5", "3.5", "(1,(2,0)) 3.5 (0,3)"))
  # as in stanc, size() and other array functions do not take arrays of tuples
  expect_error(bayescomp_check("transformed data { array[2] tuple(int, real) xs = {(1, 0.5), (2, 1.5)}; print(size(xs)); }"),
               class = "bayescomp_semantic_error")
})
