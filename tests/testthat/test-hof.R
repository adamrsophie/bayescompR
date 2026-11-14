test_that("algebraic solvers, integrate_1d, map_rect and reduce_sum match CmdStan log_prob", {
  code <- "
functions {
  vector sys(vector x, vector th, real c) {
    vector[2] r;
    r[1] = x[1] + th[1] * x[2]^3 - c;
    r[2] = exp(x[2]) - th[2] * x[1];
    return r;
  }
  vector sysl(vector x, vector th, array[] real x_r, array[] int x_i) { return sys(x, th, x_r[1]); }
  real ig(real x, real xc, array[] real th, array[] real x_r, array[] int x_i) {
    return exp(-th[1] * x^2) * (1 + th[2] * x);
  }
  vector mr(vector phi, vector th, data array[] real x_r, data array[] int x_i) {
    return [normal_lpdf(x_r[1] | th[1], phi[1]) + x_i[1] * th[2]]';
  }
  real ps(array[] real ys, int st, int en, real mu, real s) { return normal_lpdf(ys | mu, s); }
}
data { array[3] real y; array[2, 1] real xr; array[2, 1] int xi; }
parameters { vector<lower=0>[2] th; real<lower=0> sig; array[2] vector[2] mth; real mu; }
transformed parameters {
  vector[2] xs = solve_newton(sys, [1, 0.5]', th, 2.0);
  vector[2] xp = solve_powell_tol(sys, [1, 0.5]', 1e-10, 1e-8, 500, th, 2.0);
  vector[2] xl = algebra_solver(sysl, [1, 0.5]', th, {2.0}, {0});
  real I1 = integrate_1d(ig, -1.0, 2.0, to_array_1d(th), {0.0}, {0});
  real I2 = integrate_1d(ig, negative_infinity(), positive_infinity(), {th[1], 0}, {0.0}, {0}, 1e-8);
}
model {
  th ~ lognormal(0, 0.5); sig ~ exponential(1); mu ~ normal(0, 1);
  for (j in 1:2) mth[j] ~ normal(0, 1);
  target += sum(xs) + 2 * xp[2] - xl[1] + I1 + 0.3 * I2;
  target += sum(map_rect(mr, [sig]', mth, xr, xi));
  target += reduce_sum(ps, y, 1, mu, sig);
}
"
  m <- bayescomp_model(code, data = list(y = c(0.1, -0.4, 1.2), xr = matrix(c(0.3, -0.2), 2), xi = matrix(1:2, 2)))
  u <- bayescomp_unconstrain(m, list(th = c(0.8, 1.5), sig = 1.3, mth = rbind(c(0.2, -0.1), c(0.5, 0.4)), mu = 0.25))
  g <- bayescomp_grad_log_prob(m, u, jacobian = TRUE, propto = TRUE)
  expect_equal(attr(g, "log_prob"), -4.0715339, tolerance = 1e-7)
  expect_equal(as.numeric(g), c(-1.0548739, 0.19799047, -4.2068047, -0.1408284, 1.1, -0.91420118, 1.6, -0.1612426),
               tolerance = 1e-6)
})

test_that("solver failures raise domain errors", {
  f <- "functions { vector s(vector x, real c) { return [x[1]^2 + c]'; } }\n"
  expect_error(bayescomp_model(paste0(f, "transformed data { vector[1] x = solve_newton(s, [1]', 1.0); }")),
               class = "bayescomp_domain_error")
  expect_error(bayescomp_model(paste0(f, "transformed data { vector[1] x = solve_powell_tol(s, [1]', 1e-10, 1e-6, 3, 1.0); }")),
               class = "bayescomp_domain_error")
  out <- utils::capture.output(invisible(bayescomp_model(paste0(f, "transformed data { print(solve_newton(s, [1]', -4.0)); }"))))
  expect_identical(out, "[2]")
})
