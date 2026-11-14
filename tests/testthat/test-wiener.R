# Reference values from CmdStan 2.40 for the same programs.

test_that("wiener densities and distribution functions print as in CmdStan", {
  code <- "
transformed data {
  real a = 1.3; real t0 = 0.12; real w = 0.45; real v = 0.7; real sv = 0.4; real sw = 0.15; real st0 = 0.06;
  vector[3] y = [0.45, 0.8, 1.6]';
  print(wiener_lpdf(y | a, t0, w, v));
  print(wiener_lpdf(y[1] | a, t0, w, v, sv));
  print(wiener_lpdf(y[1] | a, t0, w, v, 0));
  print(wiener_lpdf(y[2] | a, t0, w, v, sv, sw, st0));
  print(wiener_lpdf(y[3] | a, t0, w, v, 0, sw, 0));
  print(wiener_lpdf(y[1] | a, t0, w, v, 0, 0, st0));
  print(wiener_lcdf_unnorm(y[2], a, t0, w, v));
  print(wiener_lccdf_unnorm(y[3], a, t0, w, v));
  print(wiener_lcdf_unnorm(y[1], a, t0, w, v, sv, sw, st0));
  print(wiener_lccdf_unnorm(y[2], a, t0, w, v, sv, 0, st0));
  print(wiener_lcdf_unnorm(5.0, a, t0, w, v), \" \", wiener_lccdf_unnorm(0.13, a, t0, w, v));
}
parameters { real z; } model { z ~ normal(0, 1); }
"
  out <- utils::capture.output(invisible(bayescomp_model(code)))
  expect_identical(out, c("-4.59476", "0.0191125", "0.0270385", "-1.00335", "-3.58544", "0.110086", "-0.586871", "-4.72824", "-1.18215", "-2.16824", "-0.40461 -0.404609"))
})

test_that("wiener gradients match CmdStan log_prob", {
  code <- "
data { int N; vector[N] y; }
parameters { real<lower=0> a; real<lower=0, upper=0.15> t0; real<lower=0.2, upper=0.8> w; real v; real<lower=0> sv; real<lower=0, upper=0.3> sw; real<lower=0, upper=0.1> st0; }
model {
  a ~ normal(1.5, 0.5); v ~ normal(0, 1); sv ~ normal(0.5, 0.3); sw ~ normal(0.1, 0.05); st0 ~ normal(0.05, 0.02);
  target += wiener_lpdf(y | a, t0, w, v);
  target += wiener_lpdf(y[1] | a, t0, w, v, sv);
  target += wiener_lpdf(y | rep_vector(a, N), rep_vector(t0, N), rep_vector(w, N), rep_vector(v, N), rep_vector(sv, N), 1e-6);
  target += wiener_lpdf(y[2] | a, t0, w, v, sv, sw, st0);
  target += wiener_lpdf(y[3] | a, t0, w, v, 0, sw, 0);
  target += wiener_lpdf(y[1] | a, t0, w, v, 0, 0, st0);
  target += wiener_lcdf_unnorm(y[2], a, t0, w, v);
  target += wiener_lccdf_unnorm(y[3], a, t0, w, v);
  target += wiener_lcdf_unnorm(y[1], a, t0, w, v, sv, sw, st0);
  target += wiener_lccdf_unnorm(y[2], a, t0, w, v, sv, 0, st0);
}
"
  m <- bayescomp_model(code, data = list(N = 3L, y = c(0.45, 0.8, 1.6)))
  u <- bayescomp_unconstrain(m, list(a = 1.3, t0 = 0.12, w = 0.45, v = 0.7, sv = 0.4, sw = 0.15, st0 = 0.06))
  g <- bayescomp_grad_log_prob(m, u, jacobian = TRUE, propto = TRUE)
  expect_equal(attr(g, "log_prob"), -36.1041, tolerance = 1e-6)
  expect_equal(as.numeric(g), c(42.181027, 0.18586072, 0.62076518, 1.5677987, 0.95333643, -1.5129501, -0.73175035),
               tolerance = 1e-6)
})

test_that("wiener arguments are checked", {
  expect_error(bayescomp_model("transformed data { print(wiener_lpdf(0.1 | 1, 0.2, 0.5, 0)); }"),
               "nondecision time", class = "bayescomp_domain_error")
  expect_error(bayescomp_model("transformed data { print(wiener_lpdf(1 | 1, 0.2, 0.3, 0, 0, 0.8, 0)); }"),
               class = "bayescomp_domain_error")
})
