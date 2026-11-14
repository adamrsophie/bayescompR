# Library functions checked against CmdStan 2.40 values and finite differences

stan_print <- function(code) {
  out <- utils::capture.output(invisible(bayescomp_model(paste("transformed data {", code, "}"))))
  txt <- gsub("][", " ", paste(out, collapse = " "), fixed = TRUE)
  as.numeric(strsplit(trimws(gsub("[][]", " ", txt)), "[ ,]+")[[1L]])
}

test_that("new cumulative functions match CmdStan 2.40", {
  ref <- list(
    "print(neg_binomial_lcdf(3 | 2.5, 0.7), \" \", neg_binomial_lccdf({1, 4} | 2.5, 0.7), \" \", neg_binomial_cdf(2 | 1.5, 0.4));" = c(-0.546753, -1.50501, 0.462448),
    "print(skew_normal_lcdf(0.3 | 0.1, 1.2, 2.0), \" \", skew_normal_lccdf([0.3, -1]' | 0.1, 1.2, -1.5), \" \", skew_normal_cdf(0.5 | 0, 1, 0.5));" = c(-1.50297, -2.48855, 0.562485),
    "print(exp_mod_normal_lcdf(0.3 | 0.1, 1.2, 2.0), \" \", exp_mod_normal_lccdf([0.3, 2]' | 0.1, 1.2, 0.5), \" \", exp_mod_normal_cdf(1 | 0, 1, 1));" = c(-0.882441, -1.04197, 0.538079),
    "print(skew_double_exponential_lcdf([0.3, -1]' | 0.1, 1.2, 0.3), \" \", skew_double_exponential_lccdf([0.3, -1]' | 0.1, 1.2, 0.3), \" \", skew_double_exponential_cdf(0.5 | 0, 1, 0.7));" = c(-3.49075, -0.543468, 0.851024),
    "print(von_mises_lcdf(0.3 | 0.1, 2.0), \" \", von_mises_lccdf([0.3, -2]' | 2.5, 0.7), \" \", von_mises_cdf(1 | -3, 60), \" \", von_mises_cdf(-1 | 0.4, 49.5));" = c(-0.506222, -0.742143, 0.862925, 6.43929e-14),
    "print(bernoulli_lccdf(0 | 0.3), \" \", bernoulli_lccdf({0, 0} | [0.3, 0.6]'));" = c(-1.20397, -1.7148),
    "print(beta_neg_binomial_lpmf(3 | 2.5, 3.0, 1.5), \" \", beta_neg_binomial_lcdf({0, 4} | 2.5, 3.0, 1.5), \" \", beta_neg_binomial_lccdf(2 | 1.5, 2.0, 2.5), \" \", beta_neg_binomial_cdf(5 | 1.5, 2.0, 2.5));" = c(-2.62753, -0.977226, -1.04784, 0.82925),
    "print(poisson_binomial_lpmf(1 | [0.2, 0.5]), \" \", poisson_binomial_lcdf({1, 2} | {0.2, 0.5, 0.7}));" = c(-0.693147, -0.670408),
    "print(ordered_logistic_lpmf({1, 3} | [0.1, 0.5]', {[-1, 1]', [-0.5, 0.7]'}), \" \", ordered_probit_lpmf({2, 3} | [0.2, -0.1]', {[-1, 1]', [-0.5, 0.7]'}));" = c(-2.18547, -1.94775),
    "print(log_mix([0.3, 0.7]', {[1, 2]', [3, 4]'}), \" \", log_mix({0.3, 0.7}, {1.0, 2.0}), \" \", lmgamma({1, 3}, [2.5, 3]'), \" \", trace_dot([[1, 0.2], [0.1, 1]], [[0.9, 0.1], [0, 0.8]]));" = c(5.57946, 1.78973, 0.284683, 2.69492, 1.71))
  for (code in names(ref)) expect_equal(stan_print(code), ref[[code]], tolerance = 1e-5, label = code)
})

test_that("multivariate additions match CmdStan 2.40", {
  td <- "matrix[2,3] Y = [[0.1, -0.3, 0.5], [0.2, 0.4, -0.1]]; matrix[3,3] S = [[2, 0.3, 0.1], [0.3, 1.5, 0.2], [0.1, 0.2, 1]];
    matrix[2,2] Fm = [[1, 0.2], [0.1, 1]]; matrix[2,2] G = [[0.9, 0.1], [0, 0.8]]; matrix[2,2] W = [[0.5, 0.1], [0.1, 0.4]];
    matrix[2,2] C0 = [[1, 0.2], [0.2, 1]]; matrix[2,2] V = [[0.3, 0.05], [0.05, 0.2]];"
  expect_equal(stan_print(paste(td, "print(multi_gp_lpdf(Y | S, [1.5, 0.7]'), \" \", multi_gp_cholesky_lpdf(Y | cholesky_decompose(S), [1.5, 0.7]'));")),
               c(-6.80439, -6.80439), tolerance = 1e-5)
  expect_equal(stan_print(paste(td, "print(gaussian_dlm_obs_lpdf(Y | Fm, G, V, W, [0.1, -0.2]', C0), \" \", gaussian_dlm_obs_lpdf(Y | Fm, G, [0.3, 0.2]', W, [0.1, -0.2]', C0));")),
               c(-6.09595, -6.08138), tolerance = 1e-5)
  expect_equal(stan_print(paste(td, "print(lkj_cov_lpdf(S | [0.1, 0.2, 0]', [1, 0.5, 2]', 2.5));")), -3.71861, tolerance = 1e-5)
})

test_that("gradients of the new functions agree with finite differences", {
  progs <- c(
    "parameters { real<lower=0> a; real<lower=0> b; } model { target += neg_binomial_lcdf(3 | a, b) + neg_binomial_lccdf({1, 4} | a, b); }",
    "parameters { real x; real<lower=0> s; real al; } model { target += skew_normal_lcdf([0.3, -1]' | x, s, al) + skew_normal_lccdf(0.2 | x, s, al); }",
    "parameters { real x; real<lower=0> s; real<lower=0> l; } model { target += exp_mod_normal_lcdf([0.3, 2]' | x, s, l) + exp_mod_normal_lccdf(0.1 | x, s, l); }",
    "parameters { real x; real<lower=0> s; real<lower=0, upper=1> t; } model { target += skew_double_exponential_lcdf([0.3, -1]' | x, s, t) + skew_double_exponential_lccdf([0.5, -2]' | x, s, t); }",
    "parameters { real<lower=-3, upper=3> m; real<lower=0> k; } model { target += von_mises_lcdf([0.3, -2]' | m, k) + von_mises_lccdf(1 | m, k); }",
    "parameters { real<lower=0> r; real<lower=0> a; real<lower=0> b; } model { target += beta_neg_binomial_lpmf({0, 3} | r, a, b) + beta_neg_binomial_lcdf(4 | r, a, b) + beta_neg_binomial_lccdf(2 | r, a, b); }",
    "parameters { vector<lower=0>[2] w; cholesky_factor_corr[3] L; } model { matrix[2,3] Y = [[0.1, -0.3, 0.5], [0.2, 0.4, -0.1]]; target += multi_gp_cholesky_lpdf(Y | L, w) + multi_gp_lpdf(Y | multiply_lower_tri_self_transpose(L), w); }",
    "parameters { real<lower=0> v; real g; } model { matrix[2,3] Y = [[0.1, -0.3, 0.5], [0.2, 0.4, -0.1]]; target += gaussian_dlm_obs_lpdf(Y | [[1, 0.2], [0.1, 1]], [[g, 0.1], [0, 0.8]], [v, 0.2]', [[0.5, 0.1], [0.1, 0.4]], [0.1, -0.2]', [[1, 0.2], [0.2, 1]]); }",
    "parameters { cov_matrix[3] S; real<lower=0> eta; } model { target += lkj_cov_lpdf(S | [0.1, 0.2, 0]', [1, 0.5, 2]', eta); }",
    "parameters { vector[2] x; } model { target += sum(bessel_first_kind({1, 2}, x)) + sum(lmgamma(2, x + 3)) + log_sum_exp(1, x)[1] + sum(fma(x, x, 1)) + sum(ldexp(x, 2)); }",
    "parameters { vector[3] x; real<lower=0> l; } model { array[3] vector[2] p = {[0, 0]', [1, x[1]]', [x[2], 2]'}; target += sum(gp_exp_quad_cov(p, 2.0, {l, 2.0})) + sum(gp_dot_prod_cov(p, l)) + sum(gp_periodic_cov(p, 1.0, l, 3.0)) + trace_dot(to_matrix(x[1:2], 2, 1), to_matrix(x[2:3], 1, 2)); }")
  for (p in progs) {
    m <- bayescomp_model(p)
    set.seed(3)
    u <- stats::runif(m$dim, -0.5, 0.5)
    g <- as.vector(bayescomp_grad_log_prob(m, u))
    ng <- num_grad(function(z) bayescomp_log_prob(m, z), u)
    expect_lt(max(abs(g - ng) / pmax(1, abs(ng))), 1e-5, label = p)
  }
})

test_that("ordered regressions check the number of cut point vectors", {
  expect_error(stan_print("print(ordered_probit_lpmf({2, 3} | 0.2, {[-1, 1]', [-0.5, 0.7]'}));"), class = "bayescomp_shape_error")
})

test_that("transform functions match CmdStan 2.40", {
  y <- "vector[3] y = [0.3, -0.7, 1.1]';"
  chk <- function(code, ref) expect_equal(stan_print(paste(y, code)), ref, tolerance = 1e-5, label = code)
  chk("print(lower_bound_constrain(y, 1), lower_bound_unconstrain([2, 3]', [1, 0.5]'));", c(2.34986, 1.49659, 4.00417, 0, 0.916291))
  chk("print(lower_upper_bound_constrain(y, -1, [1, 2, 3]'), lower_upper_bound_unconstrain(0.5, 0, 2));", c(0.148885, -0.00456332, 2.00104, -1.09861))
  chk("print(simplex_constrain(y), simplex_unconstrain([0.2, 0.3, 0.5]'));", c(0.258877, 0.16937, 0.493509, 0.0782433, -0.286707, -0.582618))
  chk("print(sum_to_zero_constrain([[0.1, 0.2], [0.3, -0.4]]));",
      c(-0.161004, 0.248803, -0.0877992, 0.306538, -0.266667, -0.0398717, -0.145534, 0.0178633, 0.127671))
  chk("print(stochastic_row_constrain([[0.1, 0.2], [0.3, -0.4]]));", c(0.384996, 0.334224, 0.28078, 0.336199, 0.219959, 0.443842))
  chk("print(cholesky_factor_corr_unconstrain(cholesky_factor_corr_constrain(y, 3)), cov_matrix_unconstrain([[2, 0.3], [0.3, 1]]));",
      c(0.3, -0.7, 1.1, 0.346574, 0.212132, -0.023022))
  chk("print(corr_matrix_constrain(y, 3));", c(1, 0.291313, -0.604368, 0.291313, 1, 0.434041, -0.604368, 0.434041, 1))
  chk("print(simplex_constrain({y, -y}));", c(0.258877, 0.16937, 0.493509, 0.0782433, 0.157192, 0.240263, 0.0824572, 0.520088))
})

test_that("_jacobian transforms add the log Jacobian and differentiate", {
  m <- bayescomp_model("parameters { vector[3] z; } transformed parameters { simplex[4] s = simplex_jacobian(z);
    real l = lower_upper_bound_jacobian(z[1], -1, 2); cholesky_factor_corr[3] L = cholesky_factor_corr_jacobian(z, 3); } model { z ~ normal(0, 1); }")
  g <- bayescomp_grad_log_prob(m, c(0.2, -0.4, 0.9), propto = TRUE)
  # CmdStan 2.40 log_prob method
  expect_equal(as.vector(g), c(-0.89169978, 2.0392445, -3.0145085), tolerance = 1e-7)
  expect_equal(attr(g, "log_prob"), -7.1020609, tolerance = 1e-7)
  expect_error(bayescomp_model("parameters { vector[3] z; } model { vector[4] s = simplex_jacobian(z); }"), class = "bayescomp_semantic_error")
  # unconstrain of a parameter dependent value differentiates on the manifold
  u <- bayescomp_model("parameters { simplex[3] s; cov_matrix[2] S; } model { target += sum(simplex_unconstrain(s)) + sum(cov_matrix_unconstrain(S)); }")
  set.seed(1)
  p <- stats::runif(u$dim, -0.5, 0.5)
  expect_equal(as.vector(bayescomp_grad_log_prob(u, p)), num_grad(function(z) bayescomp_log_prob(u, z), p), tolerance = 1e-6)
})

test_that("hidden Markov model functions match CmdStan 2.40", {
  td <- "matrix[2,4] lo = [[-1, -2.5, -0.3, -4], [-2, -0.4, -1.7, -0.2]]; matrix[2,2] G = [[0.8, 0.2], [0.3, 0.7]]; vector[2] rho = [0.6, 0.4]';"
  expect_equal(stan_print(paste(td, "print(hmm_marginal(lo, G, rho), hmm_hidden_state_prob(lo, G, rho));")),
               c(-4.79664, 0.629946, 0.243183, 0.466711, 0.043397, 0.370054, 0.756817, 0.533289, 0.956603), tolerance = 1e-5)
  m <- bayescomp_model("parameters { simplex[2] r1; simplex[2] r2; simplex[2] rho; vector[2] mu; } model {
    matrix[2,4] lo; for (n in 1:4) for (k in 1:2) lo[k, n] = normal_lpdf(n * 0.3 | mu[k], 1);
    target += hmm_marginal(lo, [r1', r2'], rho); }
    generated quantities { array[4] int z; { matrix[2,4] lo; for (n in 1:4) for (k in 1:2) lo[k, n] = normal_lpdf(n * 0.3 | mu[k], 1); z = hmm_latent_rng(lo, [r1', r2'], rho); } }")
  set.seed(2)
  u <- stats::runif(m$dim, -0.5, 0.5)
  expect_equal(as.vector(bayescomp_grad_log_prob(m, u)), num_grad(function(z) bayescomp_log_prob(m, z), u), tolerance = 1e-6)
  z <- bayescomp_constrain(m, u, seed = 1)
  expect_true(all(unlist(z[grep("^z", names(z))]) %in% 1:2))
})
