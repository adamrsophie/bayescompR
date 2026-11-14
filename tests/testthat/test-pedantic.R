# expected warnings come from stanc 2.40 --warn-pedantic on the same programs

test_that("pedantic mode finds what stanc finds", {
  w <- bayescomp_pedantic("data { int N; vector[N] y; real x; }
parameters { real a; real b; real<lower=-2, upper=1012> d; real sig; real g; cov_matrix[2] S; real<lower=0> e; real<lower=3, upper=1> h; real u; }
model {
  real m;
  int k;
  y ~ normal(m, sig);
  a ~ normal(-100, 100);
  a ~ normal(x, 0.01);
  if (a > 0) k = 1; else k = 2;
  for (i in 1:k) target += 1;
  while (b > 10) { }
  1 ~ poisson(g);
  e ~ gamma(0.001, 0.001);
  u ~ uniform(0, 1);
  h ~ normal(0, 1);
  d ~ normal(0, 1);
  y ~ student_t(g, 0, 1);
  y ~ normal(0, 1 - 2);
  log(e) ~ normal(0, 1);
  target += normal_lpdf(b | 0, 1);
}")
  has <- function(p) expect_true(any(grepl(p, w$message, fixed = TRUE)), info = p)
  has("The parameter a has 2 priors.")
  has("parameter d with a lower and upper bound")
  has("The parameter S was declared but was not used in the density calculation.")
  has("Parameter h has constraints that don't make sense.")
  has("The parameter y is on the left-hand side of more than one tilde statement.")
  has("The variable m may not have been assigned a value before its first use.")
  has("A normal distribution is given parameter sig as a scale parameter (argument 2), but sig was not constrained to be strictly positive.")
  has("Argument -100 suggests")
  has("Argument 0.01 suggests")
  has("A control flow statement depends on parameter(s): a.")
  has("A control flow statement depends on parameter(s): b.")
  has("A poisson distribution is given parameter g as a rate parameter (argument 1)")
  has("There is a gamma or inverse-gamma distribution")
  has("Parameter u is given a uniform distribution.")
  has("A student_t distribution is given parameter g as degrees of freedom (argument 1)")
  has("A normal distribution is given value -1 as a scale parameter (argument 2), but a scale parameter is not strictly positive.")
  has("Left-hand side of distribution statement (~) may contain a non-linear transform")
  expect_false(any(grepl("parameter b is on the left", w$message)))
  expect_identical(sum(grepl("control flow", w$message)), 3L)
  expect_output(print(w), "Warning at line")
})

test_that("priors follow the factor graph and int scalars are not data", {
  w <- bayescomp_pedantic("data { real x; array[3] int k; int j; }
    parameters { real a; real b; real c; real d; real<lower=0> l; real<lower=0> h; }
    model { a ~ normal(0, 1); x ~ normal(a, 1); b ~ normal(x, 1); c ~ normal(b, 1);
            d ~ normal(0, 1); 1 ~ normal(d, 1); k ~ poisson(l); j ~ poisson(h); }")
  expect_identical(sort(w$message[grepl("prior", w$message)]),
                   sort(c(sprintf("The parameter %s has no priors. This means either no prior is provided, or the prior(s) depend on data variables. In the later case, this may be a false positive.", c("c", "l")),
                          "The parameter d has 2 priors.")))
})

test_that("uninitialized variables are found by definite assignment", {
  code <- "transformed data { real t; if (1 > 2) t = 1; else t = 2; print(t); real v; if (1 > 2) v = 1; print(v); }
           parameters { real a; } model { a ~ normal(0, 1); }
           generated quantities { real u; u = u * 2; real z = a; }"
  w <- bayescomp_uninitialized(code)
  expect_identical(w$message, c("The variable v may not have been assigned a value before its first use.",
                                "The variable u may not have been assigned a value before its first use."))
  expect_output(bayescomp_check(code, warn_uninitialized = TRUE), "variable v")
  expect_identical(nrow(bayescomp_pedantic("parameters { real<lower=0> s; } model { s ~ normal(0, 1); }")), 0L)
})
