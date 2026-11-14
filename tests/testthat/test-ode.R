rhs_fun <- "functions { vector rhs(real t, vector y, real k) { return -k * y; } }\n"

test_that("ODE solvers match the exact solution", {
  out <- utils::capture.output(invisible(bayescomp_model(paste0(rhs_fun, "transformed data {
    print(ode_rk45(rhs, [1]', 0, {1, 2, 3}, 0.5));
    print(ode_adams(rhs, [1]', 0, {1, 2, 3}, 0.5));
    print(ode_ckrk_tol(rhs, [1]', 0, {1, 2, 3}, 1e-8, 1e-8, 10000, 0.5));
    print(ode_bdf(rhs, [1]', 0, {1, 2, 3}, 0.5));
    print(ode_adjoint_tol_ctl(rhs, [2]', 0.0, {1.0}, 1e-8, [1e-8]', 1e-8, [1e-8]', 1e-8, 1e-8, 10000, 150, 1, 1, 1, 0.5));
  }"))))
  nums <- lapply(regmatches(out, gregexpr("[0-9.]+", out)), as.numeric)
  for (k in 1:4) expect_equal(nums[[k]], exp(-0.5 * 1:3), tolerance = 1e-5)
  expect_equal(nums[[5]], 2 * exp(-0.5), tolerance = 1e-5)
})

test_that("ODE sensitivities match CmdStan log_prob", {
  skip_on_cran()
  code <- "
  functions {
    vector sho(real t, vector y, real theta, vector k) {
      vector[2] dydt;
      dydt[1] = y[2];
      dydt[2] = -y[1] - theta * y[2] + k[1] * sin(t);
      return dydt;
    }
  }
  data { int T; array[T] real ts; }
  parameters { vector[2] y0; real<lower=0> theta; vector[1] k; real<lower=0> t1; }
  transformed parameters {
    array[T] vector[2] y = ode_rk45(sho, y0, 0.0, ts, theta, k);
    array[T] vector[2] yb = ode_bdf_tol(sho, y0, 0.0, ts, 1e-8, 1e-8, 100000, theta, k);
    array[1] vector[2] yt = ode_ckrk(sho, y0, 0.0, {t1}, theta, k);
  }
  model {
    y0 ~ normal(1, 1); theta ~ normal(0.5, 0.2); k ~ normal(0, 1); t1 ~ normal(2, 0.5);
    target += sum(y[3]) + 2 * y[T, 2] - square(yb[2, 1]) + yt[1, 1] * 0.5;
  }"
  m <- bayescomp_model(code, data = list(T = 4L, ts = c(0.5, 1, 2, 3.5)))
  u <- bayescomp_unconstrain(m, list(y0 = c(1.2, -0.3), theta = 0.4, k = 0.7, t1 = 1.7))
  g <- bayescomp_grad_log_prob(m, u, jacobian = TRUE, propto = TRUE)
  # CmdStan 2.40 log_prob; both solve to a relative tolerance of 1e-6
  expect_equal(attr(g, "log_prob"), -2.1549362, tolerance = 1e-5)
  expect_equal(as.numeric(g), c(-1.3500363, 0.1844482, 1.9335833, -0.22949528, 2.7557964), tolerance = 1e-5)
})

test_that("ODE calls are checked as in stanc", {
  chk <- function(code) bayescomp_check(paste0(rhs_fun, code))
  expect_error(chk("parameters { real k; } transformed parameters { array[2] vector[1] y = ode_bdf_tol(rhs, [1]', 0, {1, 2}, 1e-6, k, 100, k); }"),
               "must be data only", class = "bayescomp_semantic_error")
  expect_error(chk("transformed data { array[2] vector[1] y = ode_rk45(rhs, [1]', 0, {1, 2}, 1, 2); }"),
               "must take", class = "bayescomp_semantic_error")
  expect_error(chk("transformed data { array[2] vector[1] y = ode_rk45(nofun, [1]', 0, {1, 2}, 1); }"),
               class = "bayescomp_semantic_error")
  expect_error(bayescomp_model(paste0(rhs_fun, "transformed data { array[1] vector[1] y = ode_rk45(rhs, [1]', 1, {0.5}, 0.5); }")),
               "initial time", class = "bayescomp_domain_error")
  expect_error(bayescomp_model(paste0(rhs_fun, "transformed data { array[1] vector[1] y = ode_rk45_tol(rhs, [1]', 0, {50.0}, 1e-6, 1e-6, 5, 0.5); }")),
               "max_num_steps", class = "bayescomp_domain_error")
})

test_that("DAE solutions and sensitivities match CmdStan log_prob", {
  skip_on_cran()
  code <- "
functions {
  vector chem(real t, vector yy, vector yp, real p1, real p2, real p3) {
    vector[3] res;
    res[1] = yp[1] + p1 * yy[1] - p2 * yy[2] * yy[3];
    res[2] = yp[2] - p1 * yy[1] + p2 * yy[2] * yy[3] + p3 * yy[2] * yy[2];
    res[3] = yy[1] + yy[2] + yy[3] - 1.0;
    return res;
  }
}
data { array[3] real ts; }
parameters { real<lower=0> p1; real<lower=0> p2; real<lower=0> p3; real<lower=0, upper=1> a; }
transformed parameters {
  vector[3] yy0 = [a, 0, 1 - a]';
  vector[3] yp0 = [-p1 * a, p1 * a, 0]';
  array[3] vector[3] y = dae(chem, yy0, yp0, 0.0, ts, p1, p2, p3);
  array[3] vector[3] yt = dae_tol(chem, yy0, yp0, 0.0, ts, 1e-8, 1e-8, 10000, p1, p2, p3);
}
model {
  p1 ~ lognormal(-1, 1); p2 ~ lognormal(0, 1); p3 ~ lognormal(1, 1);
  target += 3 * y[3, 1] + 5 * y[2, 2] - y[1, 3] + 2 * yt[3, 2];
}
"
  m <- bayescomp_model(code, data = list(ts = c(0.4, 1, 2.5)))
  u <- bayescomp_unconstrain(m, list(p1 = 0.04, p2 = 1.5, p3 = 3, a = 0.9))
  g <- bayescomp_grad_log_prob(m, u, jacobian = TRUE, propto = TRUE)
  expect_equal(attr(g, "log_prob"), -5.0506134, tolerance = 1e-7)
  expect_equal(as.numeric(g), c(2.2594183, -0.39934006, -0.12348525, -0.44099808), tolerance = 1e-6)
})
