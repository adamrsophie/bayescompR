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
  # posterior 1.8.0 takes the upper tail indicator above the median
  skip_if_not_installed("posterior", "1.8.0")
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

test_that("fused density partials match finite differences", {
  ns <- asNamespace("bayescompR")
  fd <- function(f, x, h = 1e-6) vapply(seq_along(x), function(k) { e <- replace(numeric(length(x)), k, h); (f(x + e) - f(x - e)) / (2 * h) }, 0)
  cases <- list(
    list(ns$lp_student_t, c(0.3, -0.4), c(4, 0.1, 1.3)),
    list(ns$lp_cauchy, c(0.3, -0.4), c(0.1, 1.3)),
    list(ns$lp_gamma, c(0.7, 1.9), c(2.2, 1.4)),
    list(ns$lp_beta, c(0.3, 0.8), c(2.2, 1.4)))
  for (cs in cases) for (k in seq_along(cs[[3]]) + 1L) {
    y <- cs[[2]]; p <- cs[[3]]
    g <- function(z) { a <- as.list(p); a[[k - 1L]] <- ns$ad_select(z, 1L); do.call(cs[[1]], c(list(y), a, list(FALSE))) }
    r <- ns$ad_value_grad(g, p[k - 1L])
    expect_equal(r$grad, fd(function(v) ns$val(g(v)), p[k - 1L]), tolerance = 1e-6)
    r <- ns$ad_value_grad(function(z) do.call(cs[[1]], c(list(z), as.list(p), list(FALSE))), y)
    expect_equal(r$grad, fd(function(v) ns$val(do.call(cs[[1]], c(list(v), as.list(p), list(FALSE)))), y), tolerance = 1e-6)
  }
  expect_error(ns$lp_normal(c(0, NaN), 0, 1, FALSE), "is nan", class = "bayescomp_domain_error")
})

test_that("matrix multivariate and Dirichlet densities check their shapes", {
  ns <- asNamespace("bayescompR")
  th <- matrix(1 / 3, 2, 3)
  expect_error(ns$lp_dirichlet(th, matrix(1, 3, 3), FALSE), "size mismatch", class = "bayescomp_shape_error")
  expect_error(ns$lp_dirichlet(th, c(1, 1), FALSE), "size mismatch", class = "bayescomp_shape_error")
  expect_equal(ns$lp_dirichlet(matrix(0, 0, 3), c(1, 1, 1), FALSE), 0)
  expect_equal(ns$lp_dirichlet(th[1, , drop = FALSE], c(2, 2, 2), FALSE), ns$lp_dirichlet(th[1, ], c(2, 2, 2), FALSE))
  mvn <- ns$.lib$multi_normal_lpdf[[1L]]$impl
  expect_equal(mvn(matrix(0, 0, 2), c(0, 0), diag(2)), 0)
  expect_error(mvn(matrix(0, 3, 2), c(0, 0, 0), diag(2)), "size mismatch", class = "bayescomp_shape_error")
})

test_that("Cholesky Wishart densities reject invalid factors and check the degrees of freedom", {
  ns <- asNamespace("bayescompR")
  L <- t(chol(matrix(c(2, 0.5, 0.5, 1), 2)))
  full <- L; full[1, 2] <- 0.3
  for (f in list(ns$lp_wishart_cholesky, ns$lp_inv_wishart_cholesky)) {
    expect_error(f(full, 4, L, FALSE), "not lower triangular", class = "bayescomp_domain_error")
    expect_error(f(L, 4, full, FALSE), "not lower triangular", class = "bayescomp_domain_error")
    expect_error(f(L %*% diag(c(1, -1)), 4, L, FALSE), "diagonal must be positive", class = "bayescomp_domain_error")
    expect_error(f(L, 4, diag(3), FALSE), "size mismatch", class = "bayescomp_shape_error")
    expect_error(f(L, 0.5, L, FALSE), "Degrees of freedom", class = "bayescomp_domain_error")
  }
  m <- bayescomp_model("parameters { real z; } transformed parameters { matrix[2,2] L = [[1,z],[0.3,1.2]]; }
    model { target += wishart_cholesky_lpdf(L | 5, diag_matrix(rep_vector(1,2))); }")
  expect_error(bayescomp_log_prob(m, 0.1, jacobian = FALSE), "not lower triangular")
})

test_that("fixed input GP covariances agree with and without gradients at extreme scales", {
  code <- "data { array[2] vector[1] x; array[1] real ell; } parameters { real s; } model { target += sum(gp_exp_quad_cov(x, s, ell)); }"
  for (off in c(0, 1e8)) {
    m <- bayescomp_model(code, data = list(x = matrix(c(off, off + 1), 2, 1), ell = 1))
    g <- bayescomp_grad_log_prob(m, 1, jacobian = FALSE)
    expect_identical(bayescomp_log_prob(m, 1, jacobian = FALSE), attr(g, "log_prob"))
    expect_equal(attr(g, "log_prob"), 2 + 2 * exp(-0.5), tolerance = 1e-14)
    expect_equal(as.numeric(g), 2 * (2 + 2 * exp(-0.5)), tolerance = 1e-14)
  }
  m <- bayescomp_model(code, data = list(x = matrix(c(0, 1), 2, 1), ell = 1))
  expect_equal(as.numeric(bayescomp_grad_log_prob(m, 1e-170, jacobian = FALSE)), 2e-170 * (2 + 2 * exp(-0.5)), tolerance = 1e-12)
  m <- bayescomp_model("data { array[3] vector[2] x; } parameters { array[2] real<lower=0> ell; } model { target += sum(gp_exp_quad_cov(x, 1.5, ell)); }",
                       data = list(x = matrix(c(0, 1, 1, 0, 2, 2), 3, 2)))
  g <- bayescomp_grad_log_prob(m, log(c(1e-160, 0.7)), jacobian = FALSE)
  expect_true(all(is.finite(g)))
  u <- log(c(0.8, 0.7)); h <- 1e-6
  fd <- vapply(1:2, function(k) { e <- replace(c(0, 0), k, h); (bayescomp_log_prob(m, u + e, jacobian = FALSE) - bayescomp_log_prob(m, u - e, jacobian = FALSE)) / (2 * h) }, 0)
  g <- bayescomp_grad_log_prob(m, u, jacobian = FALSE)
  expect_equal(as.numeric(g), fd, tolerance = 1e-7)
  expect_identical(bayescomp_log_prob(m, u, jacobian = FALSE), attr(g, "log_prob"))
})

test_that("triangular right division and fixed ARD covariances check their arguments", {
  ns <- asNamespace("bayescompR")
  L <- matrix(c(2, 1, 0, 3), 2)
  expect_error(ns$ad_tri_solve_right(matrix(1, 2, 3), L, FALSE), "columns of the right side", class = "bayescomp_shape_error")
  expect_error(ns$ad_tri_solve_right(c(1, 2), matrix(c(0, 1, 0, 3), 2), TRUE), "singular", class = "bayescomp_domain_error")
  r <- ns$ad_value_grad(function(z) ns$ad_sum(ns$ad_tri_solve_right(ns$ad_reshape(z, c(1L, 2L)), L, TRUE)), c(1, 2))
  expect_equal(r$grad, solve(L, c(1, 1)))
  expect_error(ns$gp_sqexp_ard_fixed("gp_exp_quad_cov", matrix(0, 2, 2), matrix(0, 2, 2), 1, c(1, 1, 1)),
               "number of length scales", class = "bayescomp_shape_error")
})

test_that("pointwise output reports unsupported statements", {
  m2 <- bayescomp_model("functions { real my_lpdf(real y, real mu) { return normal_lpdf(y | mu, 1); } real my_lcdf(real y, real mu) { return normal_lcdf(y | mu, 1); }
    real my_lccdf(real y, real mu) { return normal_lccdf(y | mu, 1); } }
    data { real y; } parameters { real mu; } model { y ~ my(mu) T[-5, 5]; }", data = list(y = 0.5))
  fit2 <- bayescomp_sample(m2, chains = 1, iter_warmup = 20, iter_sampling = 2, seed = 1, refresh = 0)
  # (expect_error() fails to take a backtrace through this generated code)
  e <- tryCatch(bayescomp_log_lik(fit2), error = function(e) e)
  expect_s3_class(e, "bayescomp_unsupported_error")
  expect_match(conditionMessage(e), "truncated user defined")
  expect_error(bayescompR:::trunc_norm(list(0, 1), 2, 1, NULL, NULL, FALSE, "normal"), "lower bound is greater", class = "bayescomp_domain_error")
  m3 <- bayescomp_model("parameters { vector[2] x; } model { x ~ normal(0, 1); }")
  expect_equal(as.numeric(bayescomp_hessian_vector(m3, c(0, 0), c(0, 0))), c(0, 0))
})

pw_fixed_fit <- function(code, data, init = list(a = 0), draws = 2L) {
  m <- bayescomp_model(code, data = data)
  bayescomp_sample(m, algorithm = "fixed_param", chains = 1, iter_warmup = 0, iter_sampling = draws,
                   init = init, seed = 42, refresh = 0)
}

test_that("pointwise terms use the argument values the statement used", {
  fit <- pw_fixed_fit("data { real y; } parameters { real a; } model { y ~ normal(target(), 1); }", list(y = 1))
  expect_equal(unname(bayescomp_log_lik(fit)[1, ]), dnorm(1, 0, 1, log = TRUE))
  fit <- pw_fixed_fit("data { real y; } parameters { real a; }
    model { target += normal_lpdf(y | target(), 1); target += normal_lpdf(y | target(), 1); }", list(y = 1))
  l1 <- dnorm(1, 0, 1, log = TRUE)
  expect_equal(unname(bayescomp_log_lik(fit)[1, ]), c(l1, dnorm(1, l1, 1, log = TRUE)))
  fit <- pw_fixed_fit("functions { real my_lpdf(real y, real mu) { return normal_lpdf(y | mu, 1); } }
    data { vector[2] y; } parameters { real a; }
    model { for (i in 1:2) y[i] ~ my(a + target()); y ~ normal(a, 1) T[target() - 10, 10]; }", list(y = c(1, 2)))
  l1 <- dnorm(1, 0, 1, log = TRUE); l2 <- dnorm(2, l1, 1, log = TRUE)
  lo <- l1 + l2 - 10
  tn <- function(y) dnorm(y, log = TRUE) - log(pnorm(10) - pnorm(lo))
  expect_equal(unname(bayescomp_log_lik(fit)[1, ]), c(l1, l2, tn(1), tn(2)))
  # the model density is unchanged
  expect_equal(bayescomp_log_prob(fit$model, 0, jacobian = FALSE, propto = FALSE), sum(c(l1, l2, tn(1), tn(2))))
})

test_that("log link pointwise terms stay finite when the mean underflows", {
  for (d in c("poisson_log(eta)", "neg_binomial_2_log(eta, 2)")) {
    fit <- pw_fixed_fit(paste0("data { int y; real eta; } parameters { real a; } model { y ~ ", d, "; }"), list(y = 1L, eta = -1000))
    expect_equal(unname(bayescomp_log_lik(fit)[1, ]), -1000)
  }
  k <- bayescompR:::pw_kernels
  for (e in c(-1000, -3, 0, 5, 30)) for (y in c(0, 1, 7)) {
    expect_equal(k$poisson_log(y, e), bayescompR:::lp_poisson_log(y, e, FALSE), tolerance = 1e-12)
    expect_equal(k$neg_binomial_2_log(y, e, 2), bayescompR:::lp_neg_binomial_2_log(y, e, 2, FALSE), tolerance = 1e-12)
  }
  expect_equal(k$poisson_log(c(0, 1, 0), c(-Inf, -Inf, Inf)), c(0, -Inf, -Inf))
  fit <- pw_fixed_fit("data { array[2] int y; matrix[2, 1] x; } parameters { real a; } model { y ~ poisson_log_glm(x, -1000, [1]'); }",
                      list(y = c(0L, 1L), x = matrix(c(0, 1), 2)))
  expect_equal(unname(bayescomp_log_lik(fit)[1, ]), c(-exp(-1000), -999))
})

test_that("grouped log likelihoods reject missing labels and sum by group", {
  fit <- pw_fixed_fit("data { vector[4] y; } parameters { real a; } model { y ~ normal(a, 1); }", list(y = 1:4))
  ll <- bayescomp_log_lik(fit)
  expect_error(bayescomp_log_lik(fit, group = c("a", NA, "b", "b")), "missing values")
  g <- bayescomp_log_lik(fit, group = factor(c("b", "a", "b", "a"), levels = c("z", "a", "b")))
  expect_equal(colnames(g), c("b", "a"))
  expect_equal(unname(g), unname(cbind(ll[, 1] + ll[, 3], ll[, 2] + ll[, 4])))
  expect_equal(unname(bayescomp_log_lik(fit, group = c(2, 2, 2, 2))), unname(matrix(rowSums(ll))))
  expect_equal(dim(bayescomp_log_lik(fit, "array", group = c(1, 1, 2, 2))), c(2L, 1L, 2L))
})

test_that("two sided truncation is normalized in the far tails", {
  code <- "data { real y; real lo; real hi; } parameters { real a; } model { y ~ normal(a, 1) T[lo, hi]; }"
  ref <- function(y, a, lo, hi) {
    if (lo > a) { l1 <- pnorm(lo, a, 1, FALSE, TRUE); l2 <- pnorm(hi, a, 1, FALSE, TRUE) }
    else { l1 <- pnorm(hi, a, 1, TRUE, TRUE); l2 <- pnorm(lo, a, 1, TRUE, TRUE) }
    dnorm(y, a, log = TRUE) - (l1 + log(-expm1(l2 - l1)))
  }
  for (cs in list(c(40.5, 40, 41), c(-40.5, -41, -40), c(0.5, -1, 2), c(1, 0.9, 1.1))) {
    m <- bayescomp_model(code, data = list(y = cs[1], lo = cs[2], hi = cs[3]))
    expect_equal(bayescomp_log_prob(m, 0.3, jacobian = FALSE, propto = FALSE), ref(cs[1], 0.3, cs[2], cs[3]), tolerance = 1e-10)
    h <- 1e-6
    fd <- (bayescomp_log_prob(m, 0.3 + h, jacobian = FALSE) - bayescomp_log_prob(m, 0.3 - h, jacobian = FALSE)) / (2 * h)
    expect_equal(as.numeric(bayescomp_grad_log_prob(m, 0.3, jacobian = FALSE)), fd, tolerance = 1e-6)
  }
  m <- bayescomp_model("data { int y; } parameters { real<lower=0> l; } model { y ~ poisson(l) T[30, 40]; }", data = list(y = 32L))
  expect_equal(bayescomp_log_prob(m, log(2), jacobian = FALSE, propto = FALSE), dpois(32, 2, log = TRUE) - log(sum(dpois(30:40, 2))), tolerance = 1e-10)
})

test_that("truncated normal predictive draws are taken by inversion", {
  set.seed(1)
  for (b in list(c(0, Inf), c(8, Inf), c(40, 41), c(-41, -40), c(-1, 2))) {
    z <- bayescompR:::rtrunc_normal(rep(0, 2e4), 1, b[1], b[2])
    expect_true(all(z >= b[1] & z <= b[2]))
    q <- stats::median(z)
    cc <- if (b[1] > 0) { la <- pnorm(b[1], lower.tail = FALSE, log.p = TRUE); lb <- pnorm(b[2], lower.tail = FALSE, log.p = TRUE)
                          (1 - exp(pnorm(q, lower.tail = FALSE, log.p = TRUE) - la)) / (1 - exp(lb - la)) }
          else { la <- pnorm(b[2], log.p = TRUE); lb <- pnorm(b[1], log.p = TRUE); (exp(pnorm(q, log.p = TRUE) - la) - exp(lb - la)) / (1 - exp(lb - la)) }
    expect_equal(cc, 0.5, tolerance = 0.02)
  }
  fit <- pw_fixed_fit("data { vector[3] y; } parameters { real a; } model { y ~ normal(a, 1) T[8, ]; }", list(y = rep(8.1, 3)))
  p <- bayescomp_posterior_predict(fit, seed = 1)
  expect_false(anyNA(p)); expect_true(all(p >= 8))
})

test_that("resume keeps the controls, a finished warmup and the elliptical center", {
  m <- bayescomp_model("parameters { vector[2] x; } model { x ~ normal(0, 1); }")
  same <- function(alg, extra = list(), nw = 20L, segs = c(6L, 6L), thin = 1L) {
    common <- c(list(model = m, chains = 1L, algorithm = alg, iter_warmup = nw, init = list(x = c(.3, -.2)),
                     seed = 531L, refresh = 0L, thin = thin), extra)
    full <- do.call(bayescomp_sample, c(common, list(iter_sampling = sum(segs))))
    prev <- do.call(bayescomp_sample, c(common, list(iter_sampling = segs[1])))
    u <- prev$unconstrained[[1]]
    for (k in segs[-1]) {
      prev <- bayescomp_sample(m, resume = prev, iter_sampling = k, seed = 819L, refresh = 0)
      u <- rbind(u, prev$unconstrained[[1]])
    }
    expect_identical(u, full$unconstrained[[1]])
    expect_true(prev$config$resume$exact)
  }
  same("nuts", list(step_size_jitter = 0.3)); same("nuts", list(max_treedepth = 1L))
  same("hmc", list(int_time = 0.7)); same("metropolis", list(block = FALSE))
  for (a in c("nuts", "mala", "elliptical_slice")) same(a, segs = c(0L, 8L))
  same("elliptical_slice", nw = 0L)
  same("nuts", segs = c(4L, 5L, 6L), thin = 3L)
  f <- bayescomp_sample(m, chains = 1, iter_warmup = 20, iter_sampling = 3, seed = 1, refresh = 0)
  r <- bayescomp_sample(m, resume = f, iter_sampling = 3, refresh = 0, max_treedepth = 3)
  expect_false(r$config$resume$exact); expect_identical(r$config$resume$changed, "max_treedepth")
  m2 <- bayescomp_model("parameters { vector[2] x; } model { x ~ normal(1, 1); }")
  expect_warning(bayescomp_sample(m2, resume = f, iter_sampling = 3, refresh = 0), "not an exact continuation")
})

test_that("finite difference curvature handles scale and checks its arguments", {
  m <- bayescomp_model("parameters { vector[2] x; } model { x ~ normal(0, 1); }")
  expect_equal(unclass(bayescomp_hessian(m, c(1e12, -1e12))), -diag(2), ignore_attr = TRUE, tolerance = 1e-6)
  for (v in list(c(1e200, 2e200), c(1e-200, 2e-200)))
    expect_equal(as.numeric(bayescomp_hessian_vector(m, c(.2, -.3), v)), -v, tolerance = 1e-8)
  expect_error(bayescomp_hessian(m, c(.2, -.3), h = 0), "finite positive")
  expect_error(bayescomp_hessian_vector(m, c(.2, -.3), c(NA, 1)), "finite numeric")
  expect_equal(attr(bayescomp_hessian_vector(m, c(.2, -.3), c(1, 1)), "log_prob"), bayescomp_log_prob(m, c(.2, -.3)))
})

test_that("Gaussian GLM factorizations are kept for many statements of one model", {
  set.seed(1); K <- 70; N <- 8; P <- 2
  m <- bayescomp_model("data { int K; int N; int P; array[K] matrix[N, P] X; array[K] vector[N] Y; } parameters { real a; vector[P] b; real<lower=0> s; }
    model { for (k in 1:K) Y[k] ~ normal_id_glm(X[k], a, b, s); }",
    data = list(K = K, N = N, P = P, X = array(rnorm(K * N * P), c(K, N, P)), Y = matrix(rnorm(K * N), K)))
  u <- c(0.1, 0.2, -0.1, 0)
  g1 <- bayescomp_grad_log_prob(m, u)
  expect_length(m$env$.glm_store$entries, K)
  expect_identical(bayescomp_grad_log_prob(m, u), g1)
  expect_length(m$env$.glm_store$entries, K)
})

lp_and_grad <- function(m, u) c(bayescomp_log_prob(m, u, propto = FALSE), as.numeric(bayescomp_grad_log_prob(m, u)))

test_that("loop vectorization keeps user overloads, empty loops and effectful bounds", {
  ov <- "functions { real obs_lpdf(real y, real m) { return normal_lpdf(y | m, 1); }
    real obs_lpdf(vector y, real m) { return normal_lpdf(y | m, 2); } }
    data { vector[3] y; } parameters { real x; } model { for (i in 1:3) y[i] ~ obs(x); }"
  m <- bayescomp_model(ov, data = list(y = c(1, 2, 3)))
  expect_equal(lp_and_grad(m, 0.1), c(sum(dnorm(1:3, 0.1, 1, log = TRUE)), sum(1:3 - 0.1)))
  m <- bayescomp_model("data { int N; array[N] real y; } parameters { real x; }
    model { x ~ std_normal(); for (i in 1:N) y[i] ~ normal(0, -1); }", data = list(N = 0, y = numeric(0)))
  expect_equal(bayescomp_log_prob(m, -1, propto = FALSE), dnorm(-1, log = TRUE))
})

test_that("the prepared Gaussian GLM is accurate far from the origin and at extreme scales", {
  gl <- "data { int N; matrix[N,1] X; vector[N] y; } parameters { real a; real b; real ls; } model { y ~ normal_id_glm(X, a, [b]', exp(ls)); }"
  gd <- "data { int N; matrix[N,1] X; vector[N] y; } parameters { real a; real b; real ls; } model { y ~ normal(a + X[,1] * b, exp(ls)); }"
  both <- function(d, u) lapply(list(gl, gd), function(cd) lp_and_grad(bayescomp_model(cd, data = d), u))
  X <- matrix(0:8)
  r <- both(list(N = 9, X = X, y = 1e16 + 2 * c(X)), c(1e16, 2, 0))
  expect_equal(r[[1]], r[[2]], tolerance = 1e-10)
  expect_equal(r[[1]][1], -9 * 0.5 * log(2 * pi))
  X4 <- matrix(0:3)
  r <- both(list(N = 4, X = X4, y = 1 + c(X4)), c(1, 1, -400))
  expect_equal(r[[1]], r[[2]])
  set.seed(3); X <- matrix(rnorm(40)); y <- 1 + 2 * c(X) + rnorm(40)
  r <- both(list(N = 40, X = X, y = y), c(0.7, 1.6, -0.2))
  expect_equal(r[[1]], r[[2]], tolerance = 1e-10)
})

test_that("truncation normalizes every term of a broadcast density", {
  m <- bayescomp_model("data { real y; vector[3] off; } parameters { real x; } model { y ~ normal(off + x, 1) T[0, ]; }",
                       data = list(y = 1, off = c(-1, 0, 1)))
  mu <- c(-1, 0, 1) + 0.2
  expect_equal(bayescomp_log_prob(m, 0.2, propto = FALSE),
               sum(dnorm(1, mu, 1, log = TRUE) - pnorm(0, mu, lower.tail = FALSE, log.p = TRUE)))
})

test_that("value only samplers use the proportional target of programs that read target()", {
  m <- bayescomp_model("parameters { real x; } model { x ~ std_normal(); target += 8 * log1p_exp(target()); }")
  f <- bayescompR:::lp_value_fun(m, propto = TRUE)
  expect_equal(f(2) - f(0), bayescomp_log_prob(m, 2, propto = TRUE) - bayescomp_log_prob(m, 0, propto = TRUE))
})

test_that("ARD covariances of fixed inputs do not depend on the origin", {
  for (k in c("gp_exp_quad_cov", "gp_matern32_cov")) {
    code <- sprintf("data { array[2] vector[1] x; } parameters { real<lower=0> l; } model { target += sum(%s(x, 1.0, {l})); }", k)
    r <- lapply(list(c(0, 2), c(1e16, 1e16 + 2)), function(xx)
      lp_and_grad(bayescomp_model(code, data = list(x = matrix(xx, 2, 1))), log(3)))
    expect_equal(r[[1]], r[[2]])
  }
})

test_that("cholesky_factor_corr stays a valid factor far in the tails", {
  m <- bayescomp_model("parameters { cholesky_factor_corr[2] L; } model { target += log(L[2,2]); }")
  for (u in c(20, 300)) expect_equal(lp_and_grad(m, u), c(3 * log(1 / cosh(u)), -3 * tanh(u)))
  m3 <- bayescomp_model("parameters { cholesky_factor_corr[3] L; corr_matrix[3] R; } model { L ~ lkj_corr_cholesky(2); R ~ lkj_corr(2); }")
  u <- c(.3, -.4, .5, .2, -.1, .7)
  expect_equal(as.numeric(bayescomp_grad_log_prob(m3, u)), num_grad(function(v) bayescomp_log_prob(m3, v), u), tolerance = 1e-5)
})

test_that("negative binomial values and gradients at extreme precision and mean", {
  for (fam in c("neg_binomial_2_log(e, exp(lp))", "neg_binomial_2(exp(e), exp(lp))")) {
    m <- bayescomp_model(sprintf("data { array[3] int y; } parameters { real e; real lp; } model { y ~ %s; }", fam),
                         data = list(y = c(0L, 2L, 101L)))
    for (u in list(c(0.3, 1.2), c(0, -40), c(4, 2), c(700, 0))) {
      expect_equal(bayescomp_log_prob(m, u, propto = FALSE), sum(dnbinom(c(0, 2, 101), size = exp(u[2]), mu = exp(u[1]), log = TRUE)))
      expect_true(all(is.finite(bayescomp_grad_log_prob(m, u))))
    }
    expect_equal(as.numeric(bayescomp_grad_log_prob(m, c(0.3, 1.2))), num_grad(function(v) bayescomp_log_prob(m, v), c(0.3, 1.2)), tolerance = 1e-5)
  }
  m <- bayescomp_model("parameters { real x; } model { 1 ~ neg_binomial_2_log(0, exp(x)); }")
  expect_equal(as.numeric(bayescomp_grad_log_prob(m, -40)), 1, tolerance = 1e-8)
  m <- bayescomp_model("parameters { real x; } model { 2 ~ neg_binomial_2_log(x, 1); }")
  expect_equal(lp_and_grad(m, 710), c(-710, -1))
})

test_that("data subsets keep the transformed data draws of the full model", {
  code <- "data { int N; vector[N] y; } transformed data { real bias = normal_rng(0, 1); } parameters { real q; } model { y ~ normal(q + bias, 1); }"
  d <- list(N = 4, y = c(-1, 0, .25, 2))
  m <- bayescomp_model(code, data = d, seed = 7)
  sp <- bayescompR:::default_subset("N")
  gb <- sum(vapply(list(1:2, 3:4), function(ix) as.numeric(bayescomp_grad_log_prob(bayescompR:::model_with_data(m, sp(d, ix)), 0)), 0))
  expect_equal(gb, as.numeric(bayescomp_grad_log_prob(m, 0)))
  expect_error(bayescompR:::model_with_data(bayescomp_model(code, data = d), sp(d, 1:2)), "give the model a seed")
})

test_that("pointwise count terms broadcast a scalar outcome", {
  k <- bayescompR:::pw_kernels
  expect_equal(k$poisson_log(2, c(0, 1, 2)), dpois(2, exp(0:2), log = TRUE))
  expect_equal(k$neg_binomial_2_log(2, c(0, 1, 2), 2), dnbinom(2, size = 2, mu = exp(0:2), log = TRUE))
  expect_equal(k$neg_binomial_2_log(2, 0, 1e16), dnbinom(2, size = 1e16, mu = 1, log = TRUE))
})

test_that("pointwise replay uses each chain's CHAIN_ID and keeps indexed names", {
  m <- bayescomp_model("data { int CHAIN_ID; real y; } parameters { real a; } model { y ~ normal(CHAIN_ID + a, 1); }", data = list(y = 0))
  fit <- bayescomp_sample(m, algorithm = "fixed_param", chains = 2, iter_warmup = 0, iter_sampling = 2,
                          init = list(a = 0), seed = 42, refresh = 0)
  ll <- bayescomp_log_lik(fit, format = "array")
  expect_equal(ll[1, , 1], dnorm(0, 1:2, 1, log = TRUE))
  fit <- pw_fixed_fit("data { int N; vector[N] obs; } parameters { real a; } model { for (j in 1:N) obs[j] ~ normal(a + 0, 1); }",
                      list(N = 2, obs = c(1, 2)))
  expect_equal(colnames(bayescomp_log_lik(fit)), c("obs[1]", "obs[2]"))
})

test_that("truncated pointwise terms are -Inf outside parameter bounds", {
  fit <- pw_fixed_fit("data { real y; } parameters { real a; } model { y ~ normal(a, 1) T[a, ]; }", list(y = 0), init = list(a = -1))
  expect_true(is.finite(bayescomp_log_lik(fit)[1, 1]))
  # a draw with the bound above the observation (as from an approximation)
  fit$unconstrained[[1]][] <- 1
  expect_equal(unname(bayescomp_log_lik(fit)[, 1]), c(-Inf, -Inf))
})

test_that("resume of fixed_param keeps the thinning phase and exactness is checked", {
  m <- bayescomp_model("parameters { real a; } model { a ~ normal(0, 1); }")
  f1 <- bayescomp_sample(m, algorithm = "fixed_param", chains = 1, iter_warmup = 0, iter_sampling = 1, thin = 3,
                         init = list(a = 0), seed = 1, refresh = 0)
  f2 <- bayescomp_sample(m, algorithm = "fixed_param", iter_sampling = 1, resume = f1, refresh = 0)
  expect_equal(nrow(f1$unconstrained[[1]]) + nrow(f2$unconstrained[[1]]), 1L)
  expect_true(f2$config$resume$exact)
  f1$final[[1]]$rng <- NULL
  expect_false(bayescomp_sample(m, algorithm = "fixed_param", iter_sampling = 1, resume = f1, refresh = 0)$config$resume$exact)
})

test_that("multivariate densities check sizes and Cholesky factors", {
  expect_error(lib_fun("multi_student_t_lpdf")(c(.1, .2), 3, 0, diag(2)), "size mismatch")
  expect_error(lib_fun("multi_normal_cholesky_lpdf")(rbind(c(1, 2), c(0, 1)), c(0, 0), matrix(c(1, 0, 8, 1), 2)), "lower triangular")
})

test_that("an adaptation window of one draw keeps the metric", {
  m <- bayescomp_model("parameters { vector[2] x; } model { x ~ normal(0, 1); }")
  expect_no_error(bayescomp_sample(m, chains = 1, iter_warmup = 12, iter_sampling = 5, adapt_init_buffer = 5,
                                   adapt_term_buffer = 5, adapt_window = 1, seed = 19, refresh = 0))
})

test_that("Hessian vector products scale the step by the moved coordinates", {
  m <- bayescomp_model("parameters { vector[2] x; } model { x[1] ~ normal(0, 1); target += -log1p_exp(-x[2]) - log1p_exp(x[2]); }")
  expect_equal(as.numeric(bayescomp_hessian_vector(m, c(1e12, 0), c(0, 1))), c(0, -0.5), tolerance = 1e-6)
})

test_that("the ODE step limit applies between output times", {
  m <- bayescomp_model("functions { vector f(real t, vector y) { return rep_vector(0, 1); } }
    data { array[24] real ts; } parameters { real x; }
    model { array[24] vector[1] y = ode_rk45_tol(f, [1]', 0, ts, 1e-6, 1e-6, 12); target += x * y[24, 1]; }",
    data = list(ts = (1:24) / 24))
  expect_equal(bayescomp_log_prob(m, 0.4), 0.4)
})

test_that("SG-MCMC rejects a nonpositive step size", {
  m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; } model { y ~ normal(mu, 1); }", data = list(N = 4, y = 1:4))
  expect_error(bayescomp_sgmcmc(m, step_size = -0.01, chains = 1, iter_warmup = 2, iter_sampling = 2), "positive number")
})

test_that("evidence estimators keep the normalized density of programs that read target()", {
  m <- bayescomp_model("parameters { real x; } model { x ~ normal(0, 1); if (0) print(target()); }")
  lp <- bayescompR:::lp_value_fun(m)
  expect_equal(lp(0.5), dnorm(0.5, log = TRUE))
  expect_equal(bayescompR:::lp_value_fun(m, propto = TRUE)(0.5), -0.125)
})

test_that("value only samplers drop data only terms in parameter dependent branches", {
  m <- bayescomp_model("parameters { real x; } model { x ~ normal(0, 1); if (x > 0) 0 ~ normal(0, 1); else 0 ~ normal(0, 5); }")
  lp <- bayescompR:::lp_value_fun(m, propto = TRUE)
  expect_equal(lp(1) - lp(-1), 0)
})

test_that("for loops test their upper bound before every iteration", {
  m <- bayescomp_model("parameters { real x; } model { int K = 3; for (i in 1:K) { K -= 1; x ~ normal(0, 1); } }")
  expect_equal(bayescomp_log_prob(m, 1, propto = FALSE), 2 * dnorm(1, log = TRUE))
  m <- bayescomp_model("parameters { real x; } model { x ~ normal(0, 1); for (i in 1:2) if (target() > -2) x ~ normal(0, 1); }")
  expect_equal(bayescomp_log_prob(m, 0, propto = FALSE), 3 * dnorm(0, log = TRUE))
})

test_that("resumed chains continue generated quantities without repeating draws", {
  m <- bayescomp_model("parameters { real x; } model { x ~ normal(0, 1); } generated quantities { real z = normal_rng(0, 1); }")
  d <- function(f) as.numeric(posterior::as_draws_matrix(f$draws)[, "z"])
  r1 <- bayescomp_sample(m, chains = 1, iter_warmup = 30, iter_sampling = 10, seed = 3, refresh = 0)
  r2 <- bayescomp_sample(m, resume = r1, iter_sampling = 10, refresh = 0)
  r0 <- bayescomp_sample(m, chains = 1, iter_warmup = 30, iter_sampling = 20, seed = 3, refresh = 0)
  expect_equal(c(d(r1), d(r2)), d(r0))
  g <- bayescomp_model("generated quantities { real z = normal_rng(0, 1); }")
  f1 <- bayescomp_sample(g, algorithm = "fixed_param", chains = 1, iter_sampling = 3, seed = 2, refresh = 0)
  f2 <- bayescomp_sample(g, resume = f1, iter_sampling = 3, refresh = 0)
  expect_false(any(d(f2) %in% d(f1)))
  # a model without a seed whose transformed data draws random numbers
  td <- "transformed data { real c = normal_rng(0, 1); } parameters { real b; } model { b ~ normal(c, 1); }"
  a <- bayescomp_sample(bayescomp_model(td), chains = 1, iter_warmup = 5, iter_sampling = 2, seed = 1, refresh = 0)
  expect_warning(r <- bayescomp_sample(bayescomp_model(td), resume = a, iter_sampling = 2, refresh = 0), "differ")
  expect_false(r$config$resume$exact)
})

test_that("ARD distances, negative binomial and Student t stay finite at extremes", {
  m <- bayescomp_model("parameters { real x; } model { target += sum(gp_exp_quad_cov({[0.0]', [1.0]'}, 1.0, {exp(-400 + x)})); }")
  expect_equal(lp_and_grad(m, 0), c(2, 0))
  m <- bayescomp_model("parameters { real x; } model { 1e160 ~ student_t(3, x, 1); }")
  v <- lp_and_grad(m, 0)
  expect_equal(v[1], lgamma(2) - lgamma(1.5) - 0.5 * log(3 * pi) - 2 * (2 * log(1e160) - log(3)))
  expect_true(is.finite(v[2]))
  nb <- bayescompR:::lgamma_ratio_dphi(c(0, 3, 500, 5), c(2, 1e-300, 1e12, 0.7))
  expect_equal(nb[c(1, 2, 4)], c(0, digamma(3 + 1e-300) - digamma(1e-300), digamma(5.7) - digamma(0.7)))
  expect_equal(nb[3], 500 / 1e12, tolerance = 1e-8)
})

test_that("determinants of singular matrices have a finite gradient", {
  m <- bayescomp_model("parameters { real x; } model { target += determinant(diag_matrix([x, 1]')); }")
  expect_equal(lp_and_grad(m, 0), c(0, 1))
  m <- bayescomp_model("parameters { real x; } model { target += determinant([[x, 2], [3, x * x]]); }")
  expect_equal(as.numeric(bayescomp_grad_log_prob(m, 1.5)), 3 * 1.5^2)
  expect_error(lib_fun("mdivide_left_spd")(matrix(c(1, 2, 2, 1), 2), c(1, 1)), "positive definite")
})

test_that("complex compound assignment, powers at zero and real promotion", {
  m <- bayescomp_model("generated quantities { complex z = to_complex(1, 2); z += 1; complex w = to_complex(1, 2); w *= to_complex(2, 3);
    real n = norm(3.0); complex q = pow(to_complex(0, 0), to_complex(0, 0)); }")
  v <- unclass(posterior::as_draws_matrix(bayescomp_sample(m, algorithm = "fixed_param", chains = 1, iter_sampling = 1, refresh = 0, seed = 1)$draws))[1, ]
  expect_equal(unname(v[c("z[real]", "z[imag]", "w[real]", "w[imag]", "n", "q[real]", "q[imag]")]), c(2, 2, -4, 7, 9, 1, 0))
  m <- bayescomp_model("parameters { real t; } model { complex_vector[2] z = [to_complex(0, 0), to_complex(t, 0)]'; target += get_real((z ^ 2)[2]); }")
  expect_equal(lp_and_grad(m, 3), c(9, 6))
})

test_that("GLM densities with a shared row of x, a vector x, and short intercepts", {
  m <- bayescomp_model("data { array[4] int y; row_vector[2] x; } parameters { vector[3] a; matrix[2, 3] B; } model { y ~ categorical_logit_glm(x, a, B); }",
                       data = list(y = c(1L, 2L, 3L, 2L), x = c(0.5, -1)))
  expect_equal(bayescomp_log_prob(m, rep(0, 9)), -4 * log(3))
  m <- bayescomp_model("data { int y; row_vector[2] x; } parameters { vector[1] a; matrix[2, 3] B; } model { y ~ categorical_logit_glm(x, a, B); }",
                       data = list(y = 2L, x = c(0.5, -1)))
  expect_error(bayescomp_log_prob(m, rep(0, 7)), "size of alpha")
  m <- bayescomp_model("data { vector[4] x; vector[4] y; } parameters { real a; vector[1] b; real<lower=0> s; } model { y ~ normal_id_glm(x, a, b, s); }",
                       data = list(x = 1:4, y = 2 * (1:4)))
  expect_equal(bayescomp_log_prob(m, c(0, 2, 0), propto = FALSE), 4 * dnorm(0, log = TRUE))
})

test_that("checks and edge cases match Stan", {
  expect_error(lib_fun("rayleigh_lpdf")(0, 1), "must be positive")
  expect_error(lib_fun("loglogistic_lpdf")(0, 1, 2), "must be positive")
  m <- bayescomp_model("parameters { real x; } model { vector[2] v; real m = max(v); x ~ normal(0, 1); }")
  expect_equal(bayescomp_log_prob(m, 0, propto = FALSE), dnorm(0, log = TRUE))
  expect_equal(bayescomp_log_prob(bayescomp_model("functions { real h() { return .5; } } parameters { real x; } model { x ~ normal(h(), 1); }"), 0.5, propto = FALSE),
               dnorm(0, log = TRUE))
  expect_error(bayescomp_check("parameters { real<lower=1> s; } generated quantities { vector[to_int(floor(s))] v; }"), "only depend on data")
  expect_warning(bayescomp_model("data { array[3] int x; array[2] int idx; } model { real z = x[idx][1] / 2; }",
                                 data = list(x = 1:3, idx = 1:2)), "x\\[idx\\]\\[1\\] / 2")
  m <- bayescomp_model("functions { real my_lpdf(real y, real mu) { return normal_lupdf(y | mu, 1); } }
    parameters { real x; } model { 1.5 ~ my(floor(x)); }")
  expect_equal(bayescomp_log_prob(m, 1.2, propto = TRUE), -0.125)
  m <- bayescomp_model("parameters { real<lower=0> t; } model { matrix[2, 2] K = sqrt([[0, 1], [1, 0]] / square(t)); target += K[1, 2]; }")
  expect_true(is.finite(bayescomp_grad_log_prob(m, 0)))
  m <- bayescomp_model("data { int<lower=1> N; vector[N] y; } parameters { real mu; } model { y ~ normal(mu, 1); }", data = list(N = 4, y = 1:4))
  expect_no_error(bayescomp_sgmcmc(m, batch_size = 2, chains = 1, iter_warmup = 3, iter_sampling = 3, seed = 1))
})

test_that("proportional values without a tape equal the taped value", {
  ns <- asNamespace("bayescompR")
  progs <- list(
    list("parameters { real x; } model { x ~ normal(0, 1); if (x > 0) 0 ~ normal(0, 1); else 0 ~ normal(0, 5); }", list()),
    list("parameters { real x; } model { x ~ std_normal(); target += 8 * log1p_exp(target()); }", list()),
    list("data { int N; array[N] int y; vector[N] x; } parameters { real a; real b; real<lower=0> phi; } model { y ~ neg_binomial_2_log(a + b * x, phi); }",
         list(N = 3, y = c(0L, 3L, 1L), x = c(-1, 0.5, 2))),
    list("functions { real my_lpdf(real y, real mu) { return normal_lupdf(y | mu, 1); } } parameters { real x; } model { 1.5 ~ my(floor(x)); }", list()),
    list("data { int N; vector[N] y; matrix[N, 2] X; } parameters { real a; vector[2] b; real<lower=0> s; } transformed parameters { real c = 2; }
          model { y ~ normal_id_glm(X, a, b, s); for (n in 1:N) y[n] ~ normal(a + c, s); c ~ normal(a, 1); }",
         list(N = 3, y = c(0.1, -0.4, 1.2), X = matrix(c(1, 0, 2, -1, 0.5, 1), 3))),
    list("data { real y; } parameters { real mu; real<lower=0> s; } model { real k = 3; y ~ normal(mu, s) T[0, ]; k ~ normal(floor(mu), 1);
          target += normal_lupdf(1 | mu, s) + normal_lupdf(1 | 0, 1); }", list(y = 0.4)),
    list("functions { vector f(real t, vector y, real k) { return -k * y; } } data { array[2] real ts; } parameters { real<lower=0> k; }
          model { array[2] vector[1] y = ode_rk45(f, [1]', 0, ts, k); target += normal_lupdf(0.5 | y[2, 1], 0.1); }", list(ts = c(0.5, 1))),
    # a GLM computes its linear predictor inside the library
    list("data { int N; matrix[N, 2] X; array[N] int k; } parameters { real a; vector[2] b; } model { a ~ normal(0, 1); k ~ poisson_log_glm(X, a, b); }",
         list(N = 3, X = matrix(c(1, 0, -1, 0.5, 2, 1), 3), k = c(0L, 2L, 5L))))
  for (pr in progs) {
    m <- bayescomp_model(pr[[1]], data = pr[[2]])
    u <- seq(-0.3, 0.4, length.out = m$dim)
    ref <- ns$model_lpg(m, TRUE, TRUE)(u)$value
    expect_equal(ns$model_lp(m, TRUE, TRUE)(u), ref)
    expect_equal(ns$lp_value_fun(m, propto = TRUE)(u), ref)
  }
})

test_that("one sided truncation with vector parameters uses one normalizer call", {
  set.seed(3); N <- 25
  d <- list(N = N, y = abs(rnorm(N)) + 0.1, x = rnorm(N), k = rpois(N, 3) + 1L)
  hdr <- "data { int N; vector[N] y; vector[N] x; array[N] int k; } parameters { real a; real<lower=0> s; } model { %s }"
  u <- c(0.2, -0.1); mu <- 0.2 + d$x; s <- exp(-0.1); l <- exp(0.2 + 0.1 * d$x)
  cases <- list(
    list("y ~ normal(a + x, s) T[0, ];", sum(dnorm(d$y, mu, s, log = TRUE) - pnorm(0, mu, s, lower.tail = FALSE, log.p = TRUE))),
    list("y ~ normal(a + x, s) T[, 5];", sum(dnorm(d$y, mu, s, log = TRUE) - pnorm(5, mu, s, log.p = TRUE))),
    list("k ~ poisson(exp(a + 0.1 * x)) T[1, ];", sum(dpois(d$k, l, log = TRUE) - ppois(0, l, lower.tail = FALSE, log.p = TRUE))))
  for (cs in cases) {
    m <- bayescomp_model(sprintf(hdr, cs[[1]]), data = d)
    expect_equal(bayescomp_log_prob(m, u, propto = FALSE), cs[[2]] - 0.1, tolerance = 1e-10)
    g <- bayescomp_grad_log_prob(m, u)
    fd <- vapply(1:2, function(j) { e <- c(0, 0); e[j] <- 1e-6
      (bayescomp_log_prob(m, u + e, propto = FALSE) - bayescomp_log_prob(m, u - e, propto = FALSE)) / 2e-6 }, 0)
    expect_equal(as.numeric(g), fd, tolerance = 1e-5)
  }
  # a parameter bound gets the gradient of the normalizer
  m <- bayescomp_model("data { int N; vector[N] y; vector[N] x; } parameters { real a; real<lower=0> s; real<upper=0.4> L; } model { y ~ normal(a + x, s) T[L, ]; }",
                       data = list(N = 4, y = c(0.6, 1, 1.5, 0.9), x = c(-1, 0, 1, 0.5)))
  u <- c(0.1, -0.2, 0.3)
  fd <- vapply(1:3, function(j) { e <- 0 * u; e[j] <- 1e-6
    (bayescomp_log_prob(m, u + e, propto = FALSE) - bayescomp_log_prob(m, u - e, propto = FALSE)) / 2e-6 }, 0)
  expect_equal(as.numeric(bayescomp_grad_log_prob(m, u)), fd, tolerance = 1e-5)
})

test_that("multi_normal with a parameter covariance has finite difference gradients", {
  set.seed(2); N <- 6; x <- sort(runif(N, 0, 5)); Y <- matrix(rnorm(3 * N), 3)
  progs <- list(
    list("data { int N; array[N] real x; vector[N] y; } parameters { real<lower=0> rho; real<lower=0> a; real<lower=0> s; real m; }
          model { y ~ multi_normal(rep_vector(m, N), add_diag(gp_exp_quad_cov(x, a, rho), square(s))); }", list(N = N, x = x, y = Y[1, ])),
    list("data { int N; array[N] real x; array[3] vector[N] y; } parameters { real<lower=0> rho; real<lower=0> a; real<lower=0> s; vector[N] m; }
          model { y ~ multi_normal(m, add_diag(gp_exp_quad_cov(x, a, rho), square(s))); }", list(N = N, x = x, y = Y)),
    list("data { int N; array[N] real x; vector[N] y; } parameters { real<lower=0> rho; real<lower=0> a; real<lower=0> s; }
          model { y ~ multi_normal(rep_vector(0, N), add_diag(gp_exp_quad_cov(x, a, rho), square(s))); }", list(N = N, x = x, y = Y[1, ])),
    list("parameters { vector[3] y; cov_matrix[3] S; } model { S ~ wishart(5, diag_matrix(rep_vector(1, 3))); y ~ multi_normal([1, 0, -1]', S); }", list()))
  for (pr in progs) {
    m <- bayescomp_model(pr[[1]], data = pr[[2]])
    u <- seq(-0.3, 0.4, length.out = m$dim)
    for (pt in c(FALSE, TRUE)) {
      g <- bayescomp_grad_log_prob(m, u, propto = pt)
      fd <- vapply(seq_along(u), function(j) { e <- 0 * u; e[j] <- 1e-5
        (bayescomp_log_prob(m, u + e, propto = pt) - bayescomp_log_prob(m, u - e, propto = pt)) / 2e-5 }, 0)
      expect_equal(as.numeric(g), fd, tolerance = 1e-6)
    }
  }
})
