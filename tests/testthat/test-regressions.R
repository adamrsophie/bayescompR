# regressions for numerical and semantic edge cases

test_that("log_diff_exp, boundary truncation and tail CDF gradients are finite", {
  m <- bayescomp_model("generated quantities { real a = log_diff_exp(0, negative_infinity()); }")
  expect_equal(bayescomp_constrain(m, numeric(0))$a, 0)
  tr <- bayescomp_model("data { real y; } parameters { real<lower=0> lam; }
                         model { y ~ exponential(lam) T[0, 1]; lam ~ gamma(2, 1); }", data = list(y = 0.4))
  fd <- (bayescomp_log_prob(tr, 0.2 + 1e-6) - bayescomp_log_prob(tr, 0.2 - 1e-6)) / 2e-6
  expect_equal(as.vector(bayescomp_grad_log_prob(tr, 0.2)), fd, tolerance = 1e-5)
  tl <- bayescomp_model("parameters { real<lower=0> l; real<lower=0, upper=1> p; }
                         model { target += poisson_lcdf(0 | l) + binomial_lcdf(0 | 1000, p); }")
  g <- as.vector(bayescomp_grad_log_prob(tl, c(log(1000), stats::qlogis(0.9))))
  # chain rule through the transforms: d/du = -lambda + 1 and -1000 theta (1 - theta) / (1 - theta) + Jacobian
  expect_equal(g, c(-999, -900.8), tolerance = 1e-6)
})

test_that("repeated assignment indices give the adjoint to the last write only", {
  m <- bayescomp_model("parameters { vector[2] x; }
    model { vector[2] z = rep_vector(0, 2); z[{1, 1}] = x; target += z[1] - 0.5 * dot_self(x); }")
  expect_equal(as.vector(bayescomp_grad_log_prob(m, c(0, 0))), c(0, 1))
})

test_that("shapes: transformed data sizes, singleton arrays and size mismatches", {
  m <- bayescomp_model("transformed data { int K = 2; vector[K] v = rep_vector(1, K); }
                        parameters { real a; } model { a ~ normal(sum(v), 1); }")
  expect_equal(bayescomp_log_prob(m, 2), stats::dnorm(0, log = TRUE))
  a <- bayescomp_model("transformed data { array[2] vector[1] a = {rep_vector(1, 1), rep_vector(2, 1)}; }
                        generated quantities { real s = a[2][1]; }")
  expect_equal(bayescomp_constrain(a, numeric(0))$s, 2)
  expect_error(bayescomp_model("transformed data { vector[1] x = rep_vector(1, 1); vector[2] z = x + rep_vector(1, 2); }"),
               class = "bayescomp_shape_error")
})

test_that("CDFs check their parameters and impossible bounds are rejected", {
  m <- bayescomp_model("generated quantities { real a = normal_lcdf(0 | 0, -1); }")
  expect_error(bayescomp_constrain(m, numeric(0)), "Scale parameter", class = "bayescomp_domain_error")
  ok <- bayescomp_model("generated quantities { real b = poisson_lccdf(-1 | 3); }")
  expect_equal(bayescomp_constrain(ok, numeric(0))$b, 0)
  ib <- bayescomp_model("parameters { real<lower=positive_infinity()> a; } model { }")
  expect_error(bayescomp_log_prob(ib, 0), class = "bayescomp_domain_error")
})

test_that("diagnostics match the posterior package", {
  expect_equal(bayescomp_ess_basic(rep(c(-1, 1), 10)), 10)
  set.seed(3)
  y <- matrix(stats::rnorm(2000), 500)
  expect_equal(bayescomp_rhat(y * 1e-20), bayescomp_rhat(y))
  skip_if_not_installed("posterior")
  x <- y + matrix(rep(cumsum(stats::rnorm(500)) * 0.05, 4), 500)
  z <- matrix(stats::rbinom(2000, 1, 0.3), 500)
  a <- matrix(rep(c(-1, 1), 1000), 500)
  for (f in c("rhat", "ess_bulk", "ess_tail", "ess_basic", "mcse_mean", "mcse_sd"))
    for (m in list(x, z, a))
      expect_equal(get(paste0("bayescomp_", f))(m), get(f, asNamespace("posterior"))(m), label = f)
})

test_that("unequal minibatches give an unbiased stochastic gradient", {
  skip_on_cran()
  # N = 3, batch size 2: batches of sizes 2 and 1 used equally often
  m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; }
                        model { mu ~ normal(0, 1); y ~ normal(mu, 1); }", data = list(N = 3, y = c(0, 0, 3)))
  f <- bayescomp_sgmcmc(m, "sgld", batch_size = 2, step_size = 0.01, chains = 1, iter_warmup = 500,
                        iter_sampling = 8000, seed = 1)
  expect_equal(mean(bayescomp_extract(f, "mu")), 0.75, tolerance = 0.15)
})

test_that("particle filter increments ignore zero weight particles", {
  expect_equal(pf_increment(c(-Inf, 0), c(-Inf, 0)), 0)
})

test_that("warmup only fits and impossible diffusion horizons", {
  skip_on_cran()
  m <- bayescomp_model("parameters { real x; } model { x ~ normal(0, 1); }")
  f <- bayescomp_sample(m, chains = 1, iter_warmup = 20, iter_sampling = 20, seed = 1, refresh = 0)
  f$n_warmup_saved <- dim(f$draws)[1L]
  expect_length(kept(f), 0L)
  expect_error(bayescomp_diffusion(m, horizon = -1), "horizon")
  d <- bayescomp_diffusion(m, num_samples = 20, steps = 100, horizon = 4, inner_samples = 8, seed = 1)
  expect_true(all(is.finite(bayescomp_extract(d, "x"))))
})

test_that("generated quantities refuse a fit with different parameters", {
  skip_on_cran()
  m <- bayescomp_model("parameters { real a; real b; } model { a ~ normal(0, 1); b ~ normal(0, 1); }")
  f <- bayescomp_sample(m, chains = 1, iter_warmup = 20, iter_sampling = 20, seed = 1, refresh = 0)
  swapped <- bayescomp_model("parameters { real b; real a; } model { } generated quantities { real s = a; }")
  expect_error(bayescomp_generate_quantities(swapped, f), "parameter declarations")
})

test_that("Gibbs slice kernel leaves a bounded uniform invariant", {
  skip_on_cran()
  m <- bayescomp_model("parameters { real<lower=0, upper=1> x; } model { }")
  f <- bayescomp_gibbs(m, blocks = list(list(pars = "x", kernel = "slice")), chains = 2,
                       iter_warmup = 100, iter_sampling = 1500, seed = 1)
  x <- bayescomp_extract(f, "x")
  expect_equal(mean(x), 0.5, tolerance = 0.06)
  expect_equal(stats::var(x), 1 / 12, tolerance = 0.15)
})

test_that("gradient free samplers initialize without gradients", {
  m <- bayescomp_model("parameters { real x; } model { matrix[2, 2] A = [[x, 1], [0, 2]];
                        target += -0.5 * square(x) + 0 * svd_U(A)[1, 1]; }")
  local_mocked_bindings(ad_value_grad = function(...) stop("gradient called"))
  f <- bayescomp_sample(m, algorithm = "slice", chains = 1, iter_warmup = 20, iter_sampling = 20, seed = 1, refresh = 0)
  expect_equal(dim(f$draws)[1:2], c(20L, 1L))
})

test_that("state space parameter names line up with the values", {
  m <- bayescomp_model("parameters { matrix[2, 2] M; real s; } model { to_vector(M) ~ normal(0, 1); s ~ normal(0, 1); }")
  u <- c(1, 2, 3, 4, 5)
  v <- ssm_theta(m, u)
  expect_equal(unname(v[c("M[2,1]", "M[1,2]", "s")]), c(2, 3, 5))
})

test_that("the compatibility matrix is well formed and its test files exist", {
  cm <- bayescomp_compatibility()
  expect_true(all(cm$status %in% c("implemented", "partial", "not_implemented", "not_applicable")))
  files <- unique(unlist(strsplit(cm$tests[nzchar(cm$tests)], ";")))
  expect_true(all(file.exists(test_path(files))))
})

test_that("random number generators have Stan's return types and check parameters", {
  gq <- function(code) bayescomp_constrain(bayescomp_model(paste("generated quantities {", code, "}")), numeric(0))
  expect_true(gq("int y = categorical_rng([0.5, 0.5]');")$y %in% 1:2)
  expect_true(gq("int y = ordered_logistic_rng(0, [-1, 1]');")$y %in% 1:3)
  expect_length(gq("array[2] real y = normal_rng({0, 1}, 1);")$y, 2L)
  expect_error(gq("real a = student_t_rng(3, 0, -1);"), "Scale parameter", class = "bayescomp_domain_error")
})

test_that("cumulative functions and densities at the edge of the support", {
  v <- function(code) bayescomp_constrain(bayescomp_model(paste("generated quantities { real a =", code, "; }")), numeric(0))$a
  expect_equal(v("exponential_lccdf(-1 | 1)"), 0)
  expect_equal(v("exponential_lcdf(-1 | 1)"), -Inf)
  expect_equal(v("inv_gamma_lcdf(-1 | 2, 1)"), -Inf)
  expect_equal(v("inv_gamma_lccdf(-1 | 2, 1)"), 0)
  expect_equal(v("weibull_lccdf(-1 | 2, 1)"), 0)
  expect_equal(v("chi_square_lpdf(0 | 2)"), log(0.5))
  expect_equal(v("weibull_lpdf(0 | 1, 1)"), 0)
  expect_equal(v("std_normal_log_qf(-1000)"), stats::qnorm(-1000, log.p = TRUE))
  expect_equal(v("cauchy_lccdf(1e20 | 0, 1)"), log(stats::pcauchy(1e20, lower.tail = FALSE)), tolerance = 1e-12)
  expect_equal(v("cauchy_lcdf(-3 | 0, 1)"), stats::pcauchy(-3, log.p = TRUE))
  expect_equal(v("yule_simon_lcdf(0 | 2)"), -Inf)
  expect_equal(v("yule_simon_lccdf(0 | 2)"), 0)
  expect_equal(v("poisson_binomial_lcdf(5 | [0.5, 0.5]')"), 0)
  expect_equal(v("poisson_binomial_lccdf(-1 | [0.5, 0.5]')"), 0)
  expect_error(v("dirichlet_multinomial_lpmf({-1, 1} | [1, 1]')"), class = "bayescomp_domain_error")
  expect_error(v("hypergeometric_2F1(1, 1, 2, 0.9999999)"), "did not converge", class = "bayescomp_domain_error")
  expect_error(v("multi_student_t_lpdf({[0, 0]', [1, 1]'} | 3, {[0, 0]', [0, 0]', [0, 0]'}, diag_matrix(rep_vector(1, 2)))"),
               class = "bayescomp_shape_error")
})

test_that("gradients of the rewritten tails and densities match finite differences", {
  m <- bayescomp_model("parameters { real<lower=0> s; } model { target += cauchy_lccdf(30 | 0, s) + cauchy_lcdf(-30 | 0, s) + std_normal_log_qf(-s - 5) + chi_square_lpdf(1 | s + 2) + weibull_lpdf(0.5 | s + 1, 1); }")
  for (u in c(-1, 0.5)) {
    fd <- (bayescomp_log_prob(m, u + 1e-6) - bayescomp_log_prob(m, u - 1e-6)) / 2e-6
    expect_equal(as.vector(bayescomp_grad_log_prob(m, u)), fd, tolerance = 1e-5)
  }
})

test_that("functions bodies enforce data only arguments of library functions", {
  expect_error(bayescomp_model("functions { int f(real x) { return to_int(x); } } parameters { real x; } model { x ~ normal(f(x), 1); }"),
               "data only", class = "bayescomp_semantic_error")
  expect_s3_class(bayescomp_model("functions { int f(data real x) { return to_int(x); } } data { real d; } parameters { real x; } model { x ~ normal(f(d), 1); }",
                                  data = list(d = 1)), "bayescomp_model")
})

test_that("invalid sampler controls are refused", {
  m <- bayescomp_model("parameters { vector[2] x; } model { x ~ normal(0, 1); }")
  expect_error(bayescomp_sample(m, metric = "dense_e", inv_metric = matrix(c(1, 0.5, 0, 1), 2), chains = 1, seed = 1, refresh = 0),
               "inv_metric")
  expect_error(bayescomp_pdmp(m, dt = 0), "dt")
})

test_that("generated quantities of tempering, ensembles and RJMCMC are reproducible", {
  skip_on_cran()
  m <- bayescomp_model("parameters { real x; } model { x ~ normal(0, 1); } generated quantities { real z = normal_rng(0, 1); }")
  run <- function() bayescomp_extract(bayescomp_parallel_tempering(m, chains = 1, iter_warmup = 20, iter_sampling = 5, seed = 1), "z")
  expect_identical(run(), run())
  ens <- function() bayescomp_extract(bayescomp_ensemble(m, iter_warmup = 20, iter_sampling = 5, seed = 2), "z")
  expect_identical(ens(), ens())
})

test_that("large precision negative binomials, real lchoose and NaN comparisons match Stan", {
  v <- function(code) bayescomp_constrain(bayescomp_model(paste("generated quantities { real a =", code, "; }")), numeric(0))$a
  expect_equal(v("neg_binomial_2_lpmf(3 | 2, 1e15)"), stats::dpois(3, 2, log = TRUE))
  expect_equal(v("neg_binomial_2_log_lpmf(3 | log(2), 1e15)"), stats::dpois(3, 2, log = TRUE))
  expect_equal(v("neg_binomial_2_lpmf(3 | 2, 5)"), stats::dnbinom(3, size = 5, mu = 2, log = TRUE))
  expect_equal(v("lchoose(5.5, 2.5)"), lgamma(6.5) - lgamma(3.5) - lgamma(4), tolerance = 1e-12)
  expect_equal(v("not_a_number() != 1"), 1L)
  expect_equal(v("not_a_number() < 1"), 0L)
  expect_equal(v("fmax(not_a_number(), 2)"), 2)
  expect_equal(v("fmin(not_a_number(), 2)"), 2)
  m <- bayescomp_model("parameters { real<lower=0> p; } model { target += neg_binomial_2_lpmf(3 | 2, p); }")
  fd <- (bayescomp_log_prob(m, 1 + 1e-6) - bayescomp_log_prob(m, 1 - 1e-6)) / 2e-6
  expect_equal(as.vector(bayescomp_grad_log_prob(m, 1)), fd, tolerance = 1e-6)
  b <- bayescomp_model("parameters { real x; } model { target += binary_log_loss(1, inv_logit(x)); }")
  expect_true(is.finite(bayescomp_grad_log_prob(b, 40)))
})

test_that("elementwise ^, unary minus on arrays and unsupported higher order calls", {
  m <- bayescomp_model("parameters { vector[3] v; array[2] real a; } model { vector[3] w = v ^ 2; array[2] real b = -a; target += sum(w) + sum(b); }
                        generated quantities { array[2] vector[2] q = -{[1, 2]', [3, 4]'}; array[2] int k = -{1, 2}; }")
  expect_equal(as.vector(bayescomp_grad_log_prob(m, 1:5)), c(2, 4, 6, -1, -1))
  r <- bayescomp_constrain(m, 1:5)
  expect_equal(r$k, c(-1L, -2L))
  expect_equal(as.vector(t(r$q)), -(1:4))
  m2 <- bayescomp_model("functions { real ps(array[] real s, int a, int b, real x) { return x * sum(s); } }
                         data { array[2] real y; } parameters { real x; } model { target += reduce_sum(ps, y, 1, x); }",
                         data = list(y = c(1, 2)))
  expect_equal(as.numeric(bayescomp_grad_log_prob(m2, 0.5)), 3)
})

test_that("parallel chains print progress and match sequential draws", {
  skip_on_cran()
  skip_on_os("windows")
  m <- bayescomp_model("parameters { real y; } model { y ~ normal(0, 1); }")
  a <- bayescomp_sample(m, chains = 2, iter_warmup = 20, iter_sampling = 20, seed = 1, refresh = 0)
  expect_message(b <- bayescomp_sample(m, chains = 2, iter_warmup = 20, iter_sampling = 20, seed = 1,
                                       cores = 2, refresh = 20), "Chain 2: Iteration")
  expect_identical(bayescomp_extract(a, "y"), bayescomp_extract(b, "y"))
  expect_identical(formals(bayescomp_sample)$cores, quote(getOption("mc.cores", 1L)))
})

test_that("a failing generated quantities block keeps the row's parameters", {
  m <- bayescomp_model("parameters { real y; } model { y ~ normal(0, 1); }
                        generated quantities { real z = y; if (y > 0) reject(\"positive\"); }")
  f <- suppressWarnings(bayescomp_sample(m, chains = 1, iter_warmup = 50, iter_sampling = 50, seed = 2, refresh = 0))
  y <- bayescomp_extract(f, "y"); z <- bayescomp_extract(f, "z")
  expect_false(anyNA(y))
  expect_true(all(is.na(z[y > 0])))
  expect_equal(z[y <= 0], y[y <= 0])
})

test_that("ADVI reports convergence status and ensembles run as independent chains", {
  m <- bayescomp_model("parameters { real y; } model { y ~ normal(1, 2); }")
  v <- suppressWarnings(bayescomp_variational(m, seed = 1, iter = 2000))
  expect_true(is.logical(v$vi$converged))
  expect_type(v$vi$status, "character")
  e <- bayescomp_ensemble(m, chains = 2, iter_warmup = 100, iter_sampling = 100, seed = 1)
  expect_equal(e$config$ensembles, 2L)
  expect_equal(dim(e$draws)[2], 2L * e$config$walkers)
})

test_that("container arguments of a distribution must agree in size", {
  m <- bayescomp_model("data { vector[3] y; vector[1] mu; } parameters { real<lower=0> s; } model { y ~ normal(mu, s); }",
                       data = list(y = c(1, 2, 3), mu = 0))
  expect_error(bayescomp_log_prob(m, 0), class = "bayescomp_shape_error")
})

test_that("sd has a finite gradient at a constant vector", {
  m <- bayescomp_model("parameters { real a; } model { vector[3] v = rep_vector(a, 3); target += sd(v) - 0.5 * a^2; }")
  expect_true(all(is.finite(bayescomp_grad_log_prob(m, 0.3))))
})

test_that("functions removed from Stan are rejected", {
  for (f in c("fabs(1.0)", "multiply_log(1.0, 2.0)", "binomial_coefficient_log(3.0, 1.0)"))
    expect_error(bayescomp_model(sprintf("transformed data { real x = %s; }", f)), "removed from the Stan language")
})

test_that("operator precedence, empty bar and operator functions match stanc", {
  v <- function(e, ty = "real") bayescomp_constrain(bayescomp_model(sprintf("generated quantities { %s r = %s; }", ty, e)), numeric(0))$r
  # .* binds like *, so a' * a .* a is (a' * a) .* a, a vector
  expect_equal(v("[1, 2]' * 2 .* [3, 4]'", "vector[2]"), c(6, 16))
  expect_equal(v("[1, 2] * [3, 4]' ./ 2"), 5.5)
  expect_equal(v("std_normal_lcdf(0 |)"), log(0.5))
  expect_equal(v("add(1.0, multiply(2, 3))"), 7)
  expect_equal(v("minus(3.0)"), -3)
  expect_equal(v("logical_negation(0)", "int"), 1L)
  expect_equal(v("modulus(7, 3)", "int"), 1L)
  expect_equal(v("transpose([1, 2]')", "row_vector[2]"), c(1, 2))
})

test_that("multidimensional int arrays keep their shape when assigned", {
  m <- bayescomp_model("functions { array[,] int mk(int n, int t) { array[t, n] int V; for (j in 1:n) for (i in 1:t) V[i, j] = j; return V; } }
                        transformed data { array[1, 3] int V = mk(3, 1); }
                        parameters { vector[3] b; } model { b ~ normal(0, 1); }
                        generated quantities { vector[3] z = b[V[1]]; array[2] int d = dims(V); }")
  r <- bayescomp_constrain(m, c(1, 2, 3))
  expect_equal(r$d, c(1L, 3L))
  expect_equal(r$z, c(1, 2, 3))
})

test_that("CSR functions match dense algebra, including empty rows", {
  A <- matrix(c(0, 0, 2, 0, 0, 3, 4, 0, 0, 0, 5, 6), 3)
  m <- bayescomp_model("data { matrix[3, 4] A; }
    transformed data { vector[size(csr_extract_w(A))] w = csr_extract_w(A); array[size(w)] int v = csr_extract_v(A); array[4] int u = csr_extract_u(A); }
    parameters { vector[4] b; vector[size(w)] e; }
    model { target += -0.5 * dot_self(csr_matrix_times_vector(3, 4, w .* exp(e), v, u, b)) - 0.5 * dot_self(b) - 0.5 * dot_self(e); }
    generated quantities { matrix[3, 4] D = csr_to_dense_matrix(3, 4, w, v, u); vector[3] y = csr_matrix_times_vector(3, 4, w, v, u, b); }",
    data = list(A = A))
  u <- c(0.3, -0.2, 0.5, 1, 0.1, -0.1, 0.2, 0.05, -0.3)
  r <- bayescomp_constrain(m, u)
  expect_equal(r$D, A)
  expect_equal(as.vector(r$y), as.vector(A %*% u[1:4]))
  fd <- vapply(seq_along(u), function(i) { h <- replace(numeric(length(u)), i, 1e-6); (bayescomp_log_prob(m, u + h) - bayescomp_log_prob(m, u - h)) / 2e-6 }, 0)
  expect_equal(as.vector(bayescomp_grad_log_prob(m, u)), fd, tolerance = 1e-6)
})

test_that("ordered_logistic_glm equals ordered_logistic of the linear predictor", {
  set.seed(1)
  X <- matrix(rnorm(20), 10)
  m <- bayescomp_model("data { matrix[10, 2] X; array[10] int y; } parameters { vector[2] b; ordered[2] c; }
                        model { target += ordered_logistic_glm_lpmf(y | X, b, c); }
                        generated quantities { real d = ordered_logistic_glm_lpmf(y | X, b, c) - ordered_logistic_lpmf(y | X * b, c); }",
                       data = list(X = X, y = sample(1:3, 10, TRUE)))
  expect_equal(bayescomp_constrain(m, c(0.3, -0.2, -0.5, 0.1))$d, 0)
})

test_that("functions declared without a body can be supplied from R", {
  code <- "functions { vector ols(matrix X, vector y); }
           data { matrix[5, 2] X; vector[5] y; } transformed data { vector[2] b0 = ols(X, y); }
           parameters { real a; } model { a ~ normal(0, 1); } generated quantities { vector[2] b = b0; }"
  X <- cbind(1, 1:5); y <- c(2, 4, 5, 8, 9)
  expect_error(bayescomp_model(code, data = list(X = X, y = y)), "declared but never defined", class = "bayescomp_semantic_error")
  m <- bayescomp_model(code, data = list(X = X, y = y), functions = list(ols = function(X, y) as.vector(qr.coef(qr(X), y))))
  expect_equal(bayescomp_constrain(m, 0)$b, as.vector(qr.coef(qr(X), y)))
})

test_that("library functions can be overloaded with new signatures only", {
  m <- bayescomp_model("functions { real log1p(vector x, real k) { return sum(log1p(x)) * k; } }
                        generated quantities { real a = log1p([1, 2]', 2); real b = log1p(1.0); }")
  r <- bayescomp_constrain(m, numeric(0))
  expect_equal(r$a, 2 * sum(log1p(1:2)))
  expect_equal(r$b, log(2))
  expect_error(bayescomp_model("functions { real log1p(real x) { return x; } }"), "already defined in the Stan library")
})

test_that("digit separators are allowed between digits anywhere in a literal", {
  r <- bayescomp_constrain(bayescomp_model("generated quantities { real a = 1_000.000_5; real b = 1.5e1_0; int c = 1_000_000; }"), numeric(0))
  expect_equal(r$a, 1000.0005)
  expect_equal(r$b, 1.5e10)
  expect_identical(r$c, 1000000L)
  expect_error(bayescomp_model("generated quantities { int c = 1_; }"), class = "bayescomp_syntax_error")
})

test_that("an unassigned unconstrained transformed parameter is not an error, as in stanc3", {
  m <- bayescomp_model("parameters { real y; } transformed parameters { real w; real<lower=0> z; } model { y ~ normal(0, 1); }")
  expect_error(bayescomp_log_prob(m, 0), "must be greater than or equal to 0")
  m2 <- bayescomp_model("parameters { real y; } transformed parameters { real w; } model { y ~ normal(0, 1); }")
  expect_equal(bayescomp_log_prob(m2, 0), dnorm(0, log = TRUE))
})

test_that("#include directives are expanded from the include paths", {
  d <- tempfile("inc"); dir.create(d); sub <- file.path(d, "sub"); dir.create(sub)
  writeLines("real twice(real x) { return 2 * x; }", file.path(sub, "twice.stanfunctions"))
  writeLines("functions {\n  #include twice.stanfunctions\n}", file.path(sub, "funs.stan"))
  writeLines("#include <sub/funs.stan>\nparameters { real y; } model { y ~ normal(0, 1); }\ngenerated quantities { real z = twice(y); }", file.path(d, "main.stan"))
  m <- bayescomp_model(file = file.path(d, "main.stan"))
  expect_equal(bayescomp_constrain(m, 1.5)$z, 3)
  code <- readLines(file.path(d, "main.stan"))
  expect_error(bayescomp_model(code), "not found", class = "bayescomp_syntax_error")
  expect_equal(bayescomp_constrain(bayescomp_model(code, include_paths = d), 1)$z, 2)
  writeLines("#include self.stan", file.path(d, "self.stan"))
  expect_error(bayescomp_model(file = file.path(d, "self.stan")), "includes itself")
})

test_that("target() may be read in transformed parameters but not incremented there", {
  m <- bayescomp_model("parameters { real y; } transformed parameters { real t = target(); } model { y ~ normal(0, 1); }")
  expect_equal(bayescomp_log_prob(m, 1), dnorm(1, log = TRUE))
  expect_error(bayescomp_model("parameters { real y; } transformed parameters { real t = 1; target += 1; } model { y ~ normal(0, 1); }"),
               class = "bayescomp_semantic_error")
})

test_that("an unsupplied CHAIN_ID is the chain number, as in rstan", {
  m <- bayescomp_model("data { int<lower=1> CHAIN_ID; } transformed data { real c = CHAIN_ID; }
                        parameters { real y; } model { y ~ normal(0, 1); } generated quantities { real cid = c; }")
  f <- bayescomp_sample(m, chains = 3, iter_warmup = 20, iter_sampling = 5, seed = 1, refresh = 0)
  expect_equal(apply(f$draws[, , "cid"], 2, unique), c(1, 2, 3))
  expect_identical(m$env$V.CHAIN_ID, 1L)
  m2 <- bayescomp_model("data { int CHAIN_ID; } parameters { real y; } model { y ~ normal(0, 1); } generated quantities { int c = CHAIN_ID; }",
                        data = list(CHAIN_ID = 7L))
  f2 <- bayescomp_sample(m2, chains = 2, iter_warmup = 20, iter_sampling = 5, seed = 1, refresh = 0)
  expect_true(all(f2$draws[, , "c"] == 7))
})

test_that("posterior, bayesplot and loo generics accept a fit", {
  skip_if_not_installed("posterior")
  m <- bayescomp_model("data { vector[4] y; } parameters { real mu; } model { y ~ normal(mu, 1); }
                        generated quantities { vector[4] log_lik; for (i in 1:4) log_lik[i] = normal_lpdf(y[i] | mu, 1); }",
                       data = list(y = c(-1, 0, 0.5, 2)))
  f <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100, seed = 1, refresh = 0)
  d <- posterior::as_draws_df(f)
  expect_equal(posterior::ndraws(d), 200)
  expect_true(all(c("lp__", "mu", "log_lik[4]") %in% posterior::variables(d)))
  expect_equal(posterior::variables(posterior::as_draws(f, variable = "mu")), "mu")
  skip_if_not_installed("bayesplot")
  np <- bayesplot::nuts_params(f)
  expect_named(np, c("Chain", "Iteration", "Parameter", "Value"))
  expect_true("divergent__" %in% levels(np$Parameter))
  expect_equal(nrow(bayesplot::log_posterior(f)), 200)
  skip_if_not_installed("loo")
  l <- loo::loo(f)
  expect_s3_class(l, "psis_loo")
  expect_equal(dim(l$pointwise)[1], 4)
  skip_if_not_installed("rstantools")
  ll <- rstantools::log_lik(f)
  expect_equal(dim(ll), c(200L, 4L))
  expect_equal(unname(ll[101, 2]), unname(f$draws[1, 2, "log_lik[2]"]))
  # no y_rep generated quantity: draws from the sampling statements instead
  expect_equal(dim(rstantools::posterior_predict(f)), dim(ll))
  expect_error(rstantools::posterior_predict(f, variable = "yy"), "unknown variable")
})

test_that("truncation works with user defined cumulative functions", {
  code <- "functions { real my_lpdf(real y, real mu) { return normal_lpdf(y | mu, 1); }
             real my_lcdf(real y, real mu) { return normal_lcdf(y | mu, 1); }
             real my_lccdf(real y, real mu) { return normal_lccdf(y | mu, 1); } }
           data { real y; } parameters { real mu; } model { mu ~ normal(0, 1); y ~ my(mu) TR; }"
  ref <- "data { real y; } parameters { real mu; } model { mu ~ normal(0, 1); y ~ normal(mu, 1) TR; }"
  for (tr in c("T[0, 3]", "T[0, ]", "T[, 3]")) {
    a <- bayescomp_model(sub("TR", tr, code, fixed = TRUE), data = list(y = 1.2))
    b <- bayescomp_model(sub("TR", tr, ref, fixed = TRUE), data = list(y = 1.2))
    expect_equal(bayescomp_log_prob(a, 0.4), bayescomp_log_prob(b, 0.4))
    expect_equal(bayescomp_grad_log_prob(a, 0.4), bayescomp_grad_log_prob(b, 0.4))
  }
})

test_that("the optimizer uses Stan's stopping rules and reports why it stopped", {
  m <- bayescomp_model("parameters { real x; real y; } model { target += -(1 - x)^2 - 100 * (y - x^2)^2; }")
  for (alg in c("lbfgs", "bfgs")) {
    o <- bayescomp_optimize(m, algorithm = alg, seed = 1)
    expect_identical(o$code, 0L)
    expect_equal(unname(o$par), c(1, 1), tolerance = 1e-3)
    expect_type(o$message, "character")
  }
  o <- bayescomp_optimize(m, seed = 1, history_size = 2, tol_rel_grad = 0, tol_rel_obj = 0, tol_obj = 0, tol_param = 0)
  expect_match(o$message, "tol_grad")
  expect_warning(o <- bayescomp_optimize(m, seed = 1, iter = 3), "did not converge")
  expect_identical(o$code, 1L)
})

test_that("unit_e metric, step size jitter, fixed_param and adapt_metric", {
  m <- bayescomp_model("parameters { vector[2] y; } model { y ~ normal(0, [1, 3]'); } generated quantities { real s = sum(y); }")
  f <- bayescomp_sample(m, chains = 1, iter_warmup = 150, iter_sampling = 100, seed = 1, refresh = 0, metric = "unit_e")
  expect_equal(as.vector(f$inv_metric[[1]]), c(1, 1))
  j <- bayescomp_sample(m, chains = 1, iter_warmup = 150, iter_sampling = 100, seed = 1, refresh = 0, step_size_jitter = 0.5)
  expect_gt(length(unique(j$sampler_diagnostics[, 1, "stepsize__"])), 1)
  a <- bayescomp_sample(m, chains = 1, iter_warmup = 150, iter_sampling = 50, seed = 1, refresh = 0,
                        adapt_metric = FALSE, inv_metric = c(2, 5))
  expect_equal(as.vector(a$inv_metric[[1]]), c(2, 5))
  fp <- bayescomp_sample(m, chains = 2, iter_sampling = 20, seed = 1, refresh = 0, algorithm = "fixed_param",
                         init = list(y = c(0.5, -1)))
  expect_true(all(fp$draws[, , "y[1]"] == 0.5))
  expect_true(all(fp$draws[, , "s"] == -0.5))
  expect_error(bayescomp_sample(m, step_size_jitter = 2), "step_size_jitter")
})

test_that("socket cluster chains match sequential draws", {
  skip_on_cran()
  # the workers load the installed package, as under R CMD check
  skip_if(testthat::is_checking() == FALSE && !identical(Sys.getenv("BAYESCOMPR_TEST_PSOCK"), "true"), "needs the installed package")
  m <- bayescomp_model("parameters { real y; } model { y ~ normal(0, 1); }")
  a <- bayescomp_sample(m, chains = 2, iter_warmup = 30, iter_sampling = 20, seed = 1, refresh = 0)
  old <- options(bayescompR.parallel = "psock"); on.exit(options(old))
  b <- bayescomp_sample(m, chains = 2, iter_warmup = 30, iter_sampling = 20, seed = 1, refresh = 0, cores = 2)
  expect_identical(a$draws, b$draws)
})

test_that("runtime errors report the Stan statement they arose in", {
  m <- bayescomp_model("parameters { real y; }\nmodel {\n  y ~ normal(0, 1);\n  target += normal_lpdf(y | 0, -1);\n}")
  expect_error(bayescomp_log_prob(m, 0), "Scale parameter is -1.*\\(line 4, column 3\\)", class = "bayescomp_domain_error")
  expect_error(bayescomp_model("transformed data {\n  vector[2] a = [1, 2]';\n  vector[3] b = a;\n}"), "\\(line 3, column 13\\)", class = "bayescomp_shape_error")
  r <- bayescomp_model("parameters { real y; } model { y ~ normal(0, 1); }\ngenerated quantities {\n  if (y > -100) reject(\"stop\");\n}")
  expect_error(bayescomp_constrain(r, 0), "stop \\(line 3, column 3\\)", class = "bayescomp_reject_error")
  # parameter transforms and checks after the statements carry no stale location
  k <- bayescomp_model("transformed data { real a = 1; } parameters { real y; } transformed parameters { real<lower=0> z = y; } model { y ~ normal(0, 1); }")
  e <- tryCatch(bayescomp_log_prob(k, -1), error = identity)
  expect_null(e$loc)
})

test_that("the negative binomial precision partial handles mixed small and huge counts", {
  m <- bayescomp_model("data { int N; array[N] int y; } parameters { real<lower=0> mu; real<lower=0> phi; } model { y ~ neg_binomial_2(mu, phi); }",
                       data = list(N = 3, y = c(0L, 2L, 1000000L)))
  u <- c(log(5), log(2))
  f <- function(v) bayescomp_log_prob(m, v)
  fd <- vapply(1:2, function(i) { e <- replace(numeric(2), i, 1e-6); (f(u + e) - f(u - e)) / 2e-6 }, 0)
  expect_equal(as.vector(bayescomp_grad_log_prob(m, u)), fd, tolerance = 1e-6)
  expect_lt(system.time(bayescomp_grad_log_prob(m, u))[["elapsed"]], 1)
})

test_that("the right triangular solve matches the dense solve and its gradient", {
  L <- matrix(c(2, 0.3, -0.4, 0, 1.5, 0.2, 0, 0, 1.1), 3)
  B <- matrix(c(0.2, -1, 0.5, 0.3, 1.2, -0.7), 2)
  expect_equal(ad_tri_solve_right(B, L, FALSE), B %*% solve(L))
  expect_equal(ad_tri_solve_right(B[1, ], L, TRUE), as.vector(B[1, ] %*% solve(L)))
  m <- bayescomp_model("parameters { matrix[2, 3] B; cholesky_factor_cov[3] L; } model { target += -0.5 * sum(square(mdivide_right_tri_low(B, L))); to_vector(B) ~ normal(0, 1); }")
  u <- seq(-0.5, 0.5, length.out = bayescomp_dim(m))
  f <- function(v) bayescomp_log_prob(m, v)
  fd <- vapply(seq_along(u), function(i) { e <- replace(numeric(length(u)), i, 1e-6); (f(u + e) - f(u - e)) / 2e-6 }, 0)
  expect_equal(as.vector(bayescomp_grad_log_prob(m, u)), fd, tolerance = 1e-6)
})

test_that("several fixed covariance matrices are each factored once", {
  .mvn_cache$entries <- NULL
  A <- diag(3) + 0.5; B <- diag(3) * 2
  for (i in 1:5) { a <- mvn_factor(A); b <- mvn_factor(B) }
  expect_length(.mvn_cache$entries, 2L)
  expect_equal(a$L, t(chol(A)))
  expect_equal(b$L, t(chol(B)))
})

test_that("MALA and static HMC reuse the density of the retained state", {
  calls <- 0L
  bayescomp_register_backend("count_ldg", c("log_density", "gradient"), ir_versions = "1.0.0",
    compile = function(ir, model) list(log_density_gradient = function(u, jacobian, propto) {
      calls <<- calls + 1L
      list(value = -0.5 * sum(u^2), grad = -u)
    }))
  on.exit(rm("count_ldg", envir = bayescompR:::.backends))
  m <- bayescomp_model("parameters { vector[2] x; } model { x ~ normal(0, 1); }", backend = "count_ldg",
                       allow_fallback = TRUE)
  run <- function(alg, n, eps) {
    calls <<- 0L
    f <- bayescomp_sample(m, algorithm = alg, chains = 1, iter_warmup = 0, iter_sampling = n, step_size = eps,
                          seed = 1, refresh = 0)
    d <- as.data.frame(f$sampler_diagnostics)
    list(calls = calls, nl = d[[grep("n_leapfrog__", names(d))]], acc = d[[grep("accept_stat__", names(d))]])
  }
  for (alg in c("mala", "hmc")) {
    one <- run(alg, 1, 1.5)
    many <- run(alg, 40, 1.5)
    # both accepted and rejected proposals occur
    expect_true(any(many$acc < 0.5) && any(many$acc > 0.5))
    # every transition after the first starts from the retained state, so
    # the density evaluations beyond the reported counts stay those of the
    # first transition and the initialization
    expect_equal(many$calls - sum(many$nl), one$calls - sum(one$nl))
  }
  # one density evaluation per MALA transition, the proposal's
  expect_equal(run("mala", 40, 1.5)$nl[-1L], rep(1, 39L))
})

test_that("multivariate normal and Student t on matrices match the row by row sum", {
  vg <- bayescompR:::ad_value_grad
  set.seed(3)
  K <- 3; N <- 4
  Y <- matrix(rnorm(N * K), N, K); M <- matrix(rnorm(N * K), N, K); mv <- rnorm(K)
  L <- t(chol(crossprod(matrix(rnorm(K * K), K)) + diag(K)))
  mvn <- bayescompR:::lp_multi_normal_cholesky
  mvt <- bayescompR:::lp_multi_student_t_cholesky
  as_mat <- function(z, d) { z$v <- array(z$v, d); z }
  pick <- function(z, i, d) as_mat(bayescompR:::ad_gather(z, i), d)
  for (shared in c(FALSE, TRUE)) {
    mu <- if (shared) mv else M
    x0 <- c(as.vector(Y), as.vector(mu))
    whole <- function(z, st) {
      y <- pick(z, seq_len(N * K), c(N, K))
      m <- if (shared) bayescompR:::ad_gather(z, N * K + seq_len(K)) else pick(z, N * K + seq_len(N * K), c(N, K))
      if (st) mvt(y, 5, m, L, FALSE) else mvn(y, m, L, FALSE)
    }
    rows <- function(z, st) {
      out <- 0
      for (i in seq_len(N)) {
        y <- bayescompR:::ad_gather(z, i + (seq_len(K) - 1L) * N)
        m <- bayescompR:::ad_gather(z, N * K + if (shared) seq_len(K) else i + (seq_len(K) - 1L) * N)
        out <- out + if (st) mvt(y, 5, m, L, FALSE) else mvn(y, m, L, FALSE)
      }
      out
    }
    for (st in c(FALSE, TRUE)) {
      a <- vg(whole, x0, st = st); b <- vg(rows, x0, st = st)
      expect_equal(a$value, b$value, tolerance = 1e-12)
      expect_equal(a$grad, b$grad, tolerance = 1e-12)
    }
  }
  # one random vector against a matrix of locations
  yv <- rnorm(K)
  one <- function(z) mvn(bayescompR:::ad_gather(z, seq_len(K)), pick(z, K + seq_len(N * K), c(N, K)), L, FALSE)
  each <- function(z) {
    out <- 0
    for (i in seq_len(N)) out <- out + mvn(bayescompR:::ad_gather(z, seq_len(K)),
                                           bayescompR:::ad_gather(z, K + i + (seq_len(K) - 1L) * N), L, FALSE)
    out
  }
  a <- vg(one, c(yv, as.vector(M))); b <- vg(each, c(yv, as.vector(M)))
  expect_equal(a$value, b$value, tolerance = 1e-12)
  expect_equal(a$grad, b$grad, tolerance = 1e-12)
  # invalid elements are reported in the order of the rows
  bm <- M; bm[3, 1] <- Inf; bm[2, 2] <- -Inf
  expect_error(mvn(Y, bm, L, FALSE), "Location parameter is -Inf", class = "bayescomp_domain_error")
  expect_error(mvn(Y, M[1:3, ], L, FALSE), class = "bayescomp_shape_error")
})

test_that("the Cholesky correlation transform has a finite Jacobian far out and a correct gradient", {
  tf <- bayescompR:::tf_cholesky_corr
  vg <- bayescompR:::ad_value_grad
  r <- vg(function(u) tf$constrain(u, 4, NULL)$lj, c(30, -40, 0.5, 25, 1, -1))
  expect_true(is.finite(r$value) && all(is.finite(r$grad)))
  set.seed(4)
  u <- rnorm(10); W <- matrix(rnorm(25), 5)
  f <- function(u) { c <- tf$constrain(u, 5, NULL); bayescompR:::ad_sum(W * c$x) + c$lj }
  fv <- function(u) { c <- tf$constrain(u, 5, NULL); sum(W * c$x) + c$lj }
  fd <- vapply(seq_along(u), function(i) {
    e <- replace(numeric(10), i, 1e-6); (fv(u + e) - fv(u - e)) / 2e-6
  }, 0)
  expect_equal(vg(f, u)$grad, fd, tolerance = 1e-6)
  expect_equal(tf$unconstrain(tf$constrain(u, 5, NULL)$x, 5, NULL), u, tolerance = 1e-10)
})

test_that("normal_id_glm with fixed data matches the observation level density", {
  vg <- bayescompR:::ad_value_grad
  set.seed(5)
  X <- matrix(rnorm(60), 20); y <- rnorm(20)
  glm <- bayescompR:::lp_normal_id_glm
  for (propto in c(FALSE, TRUE)) {
    f <- function(z) glm(y, X, bayescompR:::ad_select(z, 1), bayescompR:::ad_gather(z, 2:4), bayescompR:::ad_select(z, 5), propto)
    g <- function(z) bayescompR:::lp_normal(y, bayescompR:::ad_select(z, 1) + bayescompR:::op_mat_vec(X, bayescompR:::ad_gather(z, 2:4)),
                                            bayescompR:::ad_select(z, 5), propto)
    z0 <- c(0.2, 0.5, -1, 0.3, 1.4)
    vg(f, z0)
    a <- vg(f, z0)
    expect_true(!is.null(bayescompR:::normal_glm_prep(y, X)))
    b <- vg(g, z0)
    expect_equal(a$value, b$value, tolerance = 1e-12)
    expect_equal(a$grad, b$grad, tolerance = 1e-10)
  }
  expect_error(bayescompR:::glm_linpred("bernoulli_logit_glm_lpmf", X, rnorm(3), rnorm(3)),
               "size of alpha", class = "bayescomp_shape_error")
})

test_that("joint objectives evaluate once per point", {
  n <- 0L
  ob <- bayescompR:::joint_objective(function(x) { n <<- n + 1L; list(value = -sum(x^2), grad = -2 * x) })
  ob$fn(c(1, 2)); ob$gr(c(1, 2)); ob$gr(c(1, 3))
  expect_equal(n, 2L)
})

test_that("SG-MCMC computes lp__ only where asked and keeps the draws", {
  d <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; } model { mu ~ normal(0, 10); y ~ normal(mu, 1); }",
                       data = list(N = 40, y = seq(-1, 3, length.out = 40)))
  a <- bayescomp_sgmcmc(d, batch_size = 5, chains = 1, iter_warmup = 20, iter_sampling = 10, seed = 1)
  b <- bayescomp_sgmcmc(d, batch_size = 5, chains = 1, iter_warmup = 20, iter_sampling = 10, seed = 1, lp_every = 3)
  col <- function(f, nm) { d <- as.data.frame(f$draws); d[[grep(nm, names(d), fixed = TRUE)]] }
  la <- col(a, "lp__"); lb <- col(b, "lp__")
  expect_equal(col(a, "mu"), col(b, "mu"))
  expect_equal(which(!is.na(lb)), c(1L, 4L, 7L, 10L))
  expect_equal(lb[c(1, 4, 7, 10)], la[c(1, 4, 7, 10)])
})

test_that("matrix Dirichlet densities equal the sum over rows", {
  vg <- bayescompR:::ad_value_grad
  set.seed(6); K <- 3; n <- 6
  TH <- matrix(rgamma(n * K, 2), n); TH <- TH / rowSums(TH); AL <- matrix(rexp(n * K) + 0.5, n)
  lp <- bayescompR:::lp_dirichlet
  sel <- function(z, i, d) { r <- bayescompR:::ad_gather(z, i); dim(r$v) <- d; r }
  for (shared in c(FALSE, TRUE)) for (propto in c(FALSE, TRUE)) {
    al <- if (shared) AL[1, ] else AL
    z0 <- c(as.vector(TH), as.vector(al))
    f <- function(z) lp(sel(z, seq_len(n * K), c(n, K)), sel(z, n * K + seq_along(al), dim(al)), propto)
    g <- function(z) Reduce(`+`, lapply(seq_len(n), function(i)
      lp(bayescompR:::ad_gather(z, i + (seq_len(K) - 1L) * n),
         if (shared) bayescompR:::ad_gather(z, n * K + seq_len(K)) else bayescompR:::ad_gather(z, n * K + i + (seq_len(K) - 1L) * n), propto)))
    a <- vg(f, z0); b <- vg(g, z0)
    expect_equal(a$value, b$value, tolerance = 1e-12)
    expect_equal(a$grad, b$grad, tolerance = 1e-10)
  }
  bad <- TH; bad[4, 2] <- -0.1
  expect_error(lp(bad, AL, FALSE), "elements must be nonnegative", class = "bayescomp_domain_error")
})

test_that("integrate_1d gradients match finite differences", {
  m <- bayescomp_model("functions { real f(real x, real xc, array[] real th, array[] real xr, array[] int xi) { return exp(-th[1] * x * x) * (1 + th[2] * x + th[3] * x * x); } }
    parameters { array[3] real<lower=0.5, upper=2> th; } model { target += log(integrate_1d(f, -1.0, 2.0, th, {0.0}, {0}, 1e-10)); }")
  u <- c(0.1, -0.3, 0.4)
  g <- bayescomp_grad_log_prob(m, u)
  lp <- function(v) attr(bayescomp_grad_log_prob(m, v), "log_prob")
  fd <- vapply(1:3, function(k) { e <- replace(numeric(3), k, 1e-5); (lp(u + e) - lp(u - e)) / 2e-5 }, 0)
  expect_equal(as.numeric(g), fd, tolerance = 1e-6)
})

test_that("each model instance keeps its own Gaussian GLM factorizations", {
  code <- "data { int N; matrix[N, 2] X; vector[N] y; } parameters { real a; vector[2] b; real<lower=0> s; } model { y ~ normal_id_glm(X, a, b, s); }"
  set.seed(7)
  ms <- lapply(1:3, function(i) bayescomp_model(code, data = list(N = 20, X = matrix(rnorm(40), 20), y = rnorm(20))))
  u <- c(0.1, 0.2, -0.3, 0.4)
  r <- lapply(ms, bayescomp_grad_log_prob, upars = u)
  for (m in ms) expect_length(m$env$.glm_store$entries, 1L)
  expect_equal(lapply(ms, bayescomp_grad_log_prob, upars = u), r)
})

test_that("squared exponential ARD covariance of fixed inputs matches the general path", {
  vg <- bayescompR:::ad_value_grad
  set.seed(8); n <- 7; D <- 2
  X <- matrix(rnorm(n * D), n); X[3, ] <- X[1, ]; X2 <- matrix(rnorm(4 * D), 4)
  W <- matrix(rnorm(n * n), n); W2 <- matrix(rnorm(n * 4), n)
  sel <- function(z, i, d) { r <- bayescompR:::ad_gather(z, i); dim(r$v) <- d; r }
  L <- bayescompR:::.lib$gp_exp_quad_cov
  f1 <- L[[which(vapply(L, function(s) identical(s$args, c("vector[]", "real", "real[]")), TRUE))]]$impl
  f2 <- L[[which(vapply(L, function(s) identical(s$args, c("vector[]", "vector[]", "real", "real[]")), TRUE))]]$impl
  z0 <- c(1.2, 0.7, 1.9, as.vector(X))
  # the inputs as a parameter take the general path
  for (cross in c(FALSE, TRUE)) {
    k <- function(x, s, l) if (cross) bayescompR:::ad_sum(f2(x, X2, s, l) * W2) else bayescompR:::ad_sum(f1(x, s, l) * W)
    a <- vg(function(z) k(X, bayescompR:::ad_select(z, 1L), bayescompR:::ad_gather(z, 2:3)), z0[1:3])
    b <- vg(function(z) k(sel(z, 3 + seq_len(n * D), c(n, D)), bayescompR:::ad_select(z, 1L), bayescompR:::ad_gather(z, 2:3)), z0)
    expect_equal(a$value, b$value, tolerance = 1e-12)
    expect_equal(a$grad, b$grad[1:3], tolerance = 1e-10)
  }
})

test_that("Cholesky Wishart densities match the covariance densities", {
  vg <- bayescompR:::ad_value_grad
  set.seed(9); K <- 3
  mk <- function() t(chol(crossprod(matrix(rnorm(K * K), K)) + diag(K)))
  LW <- mk(); LS <- mk(); lo <- which(lower.tri(LW, diag = TRUE))
  sc <- function(z, off) bayescompR:::ad_scatter(matrix(0, K, K), lo, bayescompR:::ad_gather(z, off + seq_along(lo)))
  for (inv in c(FALSE, TRUE)) {
    chol_f <- if (inv) bayescompR:::lp_inv_wishart_cholesky else bayescompR:::lp_wishart_cholesky
    full_f <- if (inv) bayescompR:::lp_inv_wishart else bayescompR:::lp_wishart
    z0 <- c(LW[lo], LS[lo], 6.5)
    nu <- function(z) bayescompR:::ad_select(z, 2L * length(lo) + 1L)
    a <- vg(function(z) chol_f(sc(z, 0), nu(z), sc(z, length(lo)), FALSE), z0)
    b <- vg(function(z) {
      W <- bayescompR:::ad_matmul(sc(z, 0), bayescompR:::ad_transpose(sc(z, 0)))
      S <- bayescompR:::ad_matmul(sc(z, length(lo)), bayescompR:::ad_transpose(sc(z, length(lo))))
      full_f(0.5 * (W + bayescompR:::ad_transpose(W)), nu(z), 0.5 * (S + bayescompR:::ad_transpose(S)), FALSE) + bayescompR:::chol_jac(sc(z, 0))
    }, z0)
    expect_equal(a$value, b$value, tolerance = 1e-12)
    expect_equal(a$grad, b$grad, tolerance = 1e-10)
  }
})

test_that("the block stage solver of sensitivity systems matches a dense solve", {
  set.seed(10); N <- 3; P <- 2
  Jy <- matrix(rnorm(N * N), N); B <- matrix(rnorm(N * P * N), N * P)
  J <- matrix(0, N * (P + 1), N * (P + 1))
  for (b in 0:P) J[b * N + seq_len(N), b * N + seq_len(N)] <- Jy
  J[-seq_len(N), seq_len(N)] <- B
  r <- rnorm(N * (P + 1))
  dense <- bayescompR:::ros_stage_solver(list(J = J), 0.1)(r)
  expect_equal(bayescompR:::ros_stage_solver(list(Jy = Jy, B = B), 0.1)(r), dense, tolerance = 1e-12)
})

test_that("Cholesky Wishart densities read only the lower triangles", {
  set.seed(11); K <- 3
  LW <- t(chol(crossprod(matrix(rnorm(K * K), K)) + diag(K))); LS <- t(chol(diag(K) * 2))
  up <- which(upper.tri(LW))
  for (f in list(bayescompR:::lp_wishart_cholesky, bayescompR:::lp_inv_wishart_cholesky)) {
    r <- bayescompR:::ad_value_grad(function(z) { L <- bayescompR:::ad_gather(z, seq_len(K * K)); dim(L$v) <- c(K, K); f(L, 5, LS, FALSE) }, as.vector(LW))
    expect_equal(r$grad[up], rep(0, length(up)))
  }
})

test_that("Hessians and Hessian vector products match the analytic values", {
  m <- bayescomp_model("parameters { vector[2] x; } model { x ~ multi_normal([0, 0]', [[2, 0.5], [0.5, 1]]); }")
  H <- bayescomp_hessian(m, c(0.1, 0.2))
  expect_equal(unclass(H)[1:2, 1:2], -solve(matrix(c(2, 0.5, 0.5, 1), 2)), tolerance = 1e-7, ignore_attr = TRUE)
  expect_equal(as.numeric(bayescomp_hessian_vector(m, c(0.1, 0.2), c(1, -2))), as.numeric(H %*% c(1, -2)), tolerance = 1e-7)
  expect_error(bayescomp_hessian_vector(m, c(0.1, 0.2), 1), "length 2")
})

test_that("ordered logistic draws take one category per linear predictor", {
  f <- bayescompR:::.lib$ordered_logistic_rng[[1L]]$impl
  set.seed(1)
  r <- f(c(-20, 0, 20), c(-1, 1))
  expect_equal(r[c(1, 3)], c(1L, 3L))
  expect_equal(f(c(0, 0), rbind(c(-20, -19), c(19, 20))), c(3L, 1L))
})

test_that("the capability report lists every solver interface", {
  cap <- bayescomp_capabilities()
  expect_true(all(c("functions", "algorithm", "controls", "ignored", "derivatives", "parallel") %in% names(cap)))
  expect_true(any(grepl("ode_bdf", cap$functions)) && any(grepl("reduce_sum", cap$functions)))
})
