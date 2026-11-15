test_that("context dependent calls are not hoisted", {
  m <- bayescomp_model("parameters { real x; }
    transformed parameters { real z = 0; if (x > 0) z = lower_bound_jacobian(2.0, 0.0); }
    model { x ~ normal(0, 1); }")
  expect_equal(bayescomp_log_prob(m, 1) - bayescomp_log_prob(m, 1, jacobian = FALSE), 2)
  expect_equal(bayescomp_log_prob(m, -1) - bayescomp_log_prob(m, -1, jacobian = FALSE), 0)
  m <- bayescomp_model("data { real y; } parameters { real x; } model { x ~ std_normal(); target += x * normal_lupdf(y | 0, 1); }",
                       data = list(y = 2))
  # propto drops a density whose arguments are all data, as Stan does
  expect_equal(as.vector(bayescomp_grad_log_prob(m, 0, propto = TRUE)), 0)
  expect_equal(as.vector(bayescomp_grad_log_prob(m, 0)), -2 - 0.5 * log(2 * pi))
})

test_that("a failing hoisted expression errors on every use without warnings", {
  m <- bayescomp_model("data { int N; vector[N] y; } parameters { real mu; }
    model { if (mu > 0) target += y[N + 1] * mu; mu ~ normal(0, 1); }", data = list(N = 3L, y = c(1, 2, 3)))
  for (i in 1:3) expect_no_warning(expect_error(bayescomp_log_prob(m, 0.5), class = "bayescomp_domain_error"))
  expect_true(is.finite(bayescomp_log_prob(m, -0.5)))
})

test_that("transformed data with CHAIN_ID is the same in every chain and run", {
  skip_on_cran()
  code <- "data { int CHAIN_ID; } transformed data { real z = normal_rng(0, 1); }
    parameters { real mu; } model { mu ~ normal(z, 1); } generated quantities { real zz = z; }"
  m <- bayescomp_model(code)
  run <- function(cores = 1L) bayescomp_sample(m, chains = 3, iter_warmup = 20, iter_sampling = 5, seed = 1, refresh = 0, cores = cores)
  a <- run()
  expect_length(unique(a$draws[1, , "zz"]), 1L)
  expect_identical(a$draws, run()$draws)
  expect_identical(a$draws, run(2L)$draws)
})

test_that("seeded calls leave a fresh session without a seed", {
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, globalenv()))
  if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv())
  bayescomp_sample(m, chains = 1, iter_warmup = 10, iter_sampling = 5, seed = 1, refresh = 0)
  expect_false(exists(".Random.seed", globalenv()))
  expect_error(bayescomp_sample(m, chains = 3, init = list(list(mu = 0.1), list(mu = 0.2)), iter_warmup = 5, iter_sampling = 5),
               "one per chain")
  expect_s3_class(bayescomp_pathfinder(m, num_paths = 1, num_draws = 2, num_psis_draws = 2, max_lbfgs_iters = 3,
                                       num_elbo_draws = 2, seed = .Machine$integer.max), "bayescomp_fit")
})

test_that("sampling statements and truncations follow stanc's rules", {
  expect_error(bayescomp_check("data { int a1; matrix[2,2] a2; row_vector[2] a3; vector[2] a4; } parameters { real mu; }
    model { a1 ~ bernoulli_logit_glm(a2, a3, a4); mu ~ normal(0, 1); }"), class = "bayescomp_semantic_error")
  expect_error(bayescomp_check("data { matrix[2,2] y; matrix[2,2] S; matrix[2,2] T; } parameters { matrix[2,2] mu; }
    model { y ~ matrix_normal_prec(mu, S, T); }"), "unknown distribution")
  m <- bayescomp_model("parameters { real y; } model { y ~ normal(0, 1) T[ , ]; }")
  expect_equal(bayescomp_log_prob(m, 0.3), stats::dnorm(0.3, log = TRUE))
  expect_error(bayescomp_check("data { vector[2] L; vector[2] y; } parameters { real mu; } model { y ~ normal(mu, 1) T[L, ]; }"),
               "truncation bound")
})

test_that("#include lines in comments are ignored and errors name the file and line", {
  d <- tempfile(); dir.create(d)
  f <- file.path(d, "a.stan")
  writeLines(c("/*", "#include nothere.stan", "*/", "parameters { real y; } model { y ~ normal(0,1); }"), f)
  expect_s3_class(bayescomp_model(file = f), "bayescomp_model")
  writeLines(c("real g(real x) {", "  return x;", "}"), file.path(d, "inc.stan"))
  writeLines(c("real h(real x) {", "  return x + q;", "}"), file.path(d, "bad.stan"))
  writeLines(c("functions {", "#include inc.stan", "}", "parameters {", "  real y;", "}", "model { y ~ normal(0, w); }"), f)
  expect_error(bayescomp_model(file = f), "line 7, column 23", fixed = TRUE)
  writeLines(c("functions {", "#include bad.stan", "}", "parameters { real y; }", "model { y ~ normal(0, 1); }"), f)
  expect_error(bayescomp_model(file = f), "in 'bad.stan', line 2", fixed = TRUE)
})

test_that("overloads rank user and library signatures by promotions", {
  m <- bayescomp_model("functions { real log1p(complex z) { return 99; } } generated quantities { real y = log1p(1.0); }")
  expect_equal(bayescomp_constrain(m, numeric(0))$y, log(2))
  m <- bayescomp_model("functions { real f(tuple(real, real) x) { return 1; } real f(tuple(real, int) x) { return 2; } }
    generated quantities { real z = f((1, 1)); }")
  expect_equal(bayescomp_constrain(m, numeric(0))$z, 2)
  m <- bayescomp_model("functions { real f(real x) { return 1; } real f(complex x) { return 2; } } generated quantities { real z = f(1); }")
  expect_equal(bayescomp_constrain(m, numeric(0))$z, 1)
  f <- bayescomp_expose_functions("functions { complex f(complex z) { return z * z; } real g(real x) { return 1; } real g(int x) { return 2; } }")
  expect_equal(f$f(1 + 2i), -3 + 4i)
  expect_equal(c(f$g(1L), f$g(1.5)), c(2, 1))
})

test_that("log_mix and erf have the right derivatives at the boundary", {
  m <- bayescomp_model("parameters { real p; } model { target += log_mix(p, log(4.0), 0.0); }")
  expect_equal(as.vector(bayescomp_grad_log_prob(m, 0)), 3)
  m <- bayescomp_model("parameters { real x; } model { target += log(erf(1e-20 * exp(x))); }")
  expect_equal(bayescomp_log_prob(m, 0), log(2 / sqrt(pi) * 1e-20))
  expect_equal(as.vector(bayescomp_grad_log_prob(m, 0)), 1)
})

test_that("generated quantities carry parameter values over to a model with changed bounds", {
  code <- "data { real s; } parameters { real<multiplier=s> x; } model {}"
  f1 <- bayescomp_generate_quantities(bayescomp_model(code, list(s = 1)), matrix(c(0, 1), ncol = 1, dimnames = list(NULL, "x")))
  f2 <- bayescomp_generate_quantities(bayescomp_model(code, list(s = 2)), f1)
  expect_equal(as.vector(bayescomp_draws(f2, "x", format = "matrix")), c(0, 1))
})

test_that("complex parameters and tuple arrays round trip", {
  m <- bayescomp_model("parameters { complex z; complex_vector[2] w; array[2, 3] tuple(real, real) a; }
    model { target += -0.5 * (norm(z) + sum(abs(w)^2)); }")
  u <- seq(0.1, 1.8, length.out = 18)
  v <- bayescomp_constrain(m, u)
  expect_true(is.complex(v$z))
  expect_equal(bayescomp_unconstrain(m, v), u)
  expect_equal(bayescomp_unconstrain(m, unflatten_params(m, flatten_draw(v))), u)
})

test_that("data generation keeps shapes and makes tuples", {
  code <- "data { int N; vector[1] y; tuple(real<lower=0>, array[2] int<lower=3, upper=4>) t; matrix[N, 1] X; } parameters { real mu; } model {}"
  f <- tempfile(fileext = ".json")
  d <- bayescomp_generate_data(code, seed = 1, file = f)
  expect_match(json_value(d$y), "^\\[")
  expect_true(d$t[["1"]] > 0 && all(d$t[["2"]] %in% 3:4))
  expect_s3_class(bayescomp_model(code, data = f), "bayescomp_model")
})

test_that("diagnose, stansummary and the density columns follow CmdStan", {
  r <- bayescomp_diagnose(test_path("fixtures", "cmdstan_sample.csv"), quiet = TRUE)
  expect_false(r$problems)
  set.seed(1)
  x <- matrix(stats::arima.sim(list(ar = 0.7), 400), 200, 2)
  d <- cmdstan_mcmc_summary(x)
  expect_equal(d[["mcse"]], stats::sd(as.vector(x)) / sqrt(ess_raw(x)))
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
  l <- bayescomp_laplace(m, draws = 5, seed = 1)
  expect_equal(l$draws[, 1, "log_g__"], -0.5 * l$draws[, 1, "mu"]^2)
  p <- bayescomp_pathfinder(m, num_paths = 1, num_draws = 20, num_psis_draws = 20, seed = 1)
  expect_equal(p$draws[, 1, "log_p__"], stats::dnorm(p$draws[, 1, "mu"], log = TRUE))
})

test_that("ADVI output files do not depend on when they are written", {
  m <- bayescomp_model("parameters { real x; } model { x ~ normal(0, 1); } generated quantities { real y = normal_rng(0, 1); }")
  v <- suppressWarnings(bayescomp_variational(m, iter = 3, eta = 0.1, output_samples = 2, elbo_samples = 2, eval_elbo = 1, seed = 1))
  expect_identical(csv_lines(v, 1, 8), csv_lines(v, 1, 8))
})

test_that("profile records survive early returns and optimization", {
  m <- bayescomp_model('functions { real f(real x) { profile("f") { return x; } } }
    parameters { real x; } model { profile("lik") { target += -0.5 * x^2; } target += f(x) - x; }')
  profile_reset()
  bayescomp_grad_log_prob(m, 1)
  expect_setequal(profile_table()$name, c("f", "lik"))
  expect_setequal(bayescomp_optimize(m, seed = 1)$profile$name, c("f", "lik"))
})

test_that("info and the formatter match stanc", {
  code <- "transformed data { print((1.5,)); real q = add(1.0, 2.0); } parameters { real x; } model { if (x > 10) reject(\"no\"); x ~ normal(0, 1); } generated quantities { }"
  expect_equal(bayescomp_info(code)$functions, c("add", "print", "reject"))
  out <- strsplit(as.character(bayescomp_format(code)), "\n")[[1L]]
  expect_true("  print((1.5,));" %in% out)
  expect_equal(tail(out, 3L), c("generated quantities {", "  ", "}"))
})

test_that("no two internal functions share a name", {
  dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(dir))
  nm <- unlist(lapply(list.files(dir, "[.]R$", full.names = TRUE), function(f) {
    l <- readLines(f)
    sub(" <- function.*", "", grep("^[A-Za-z_.][A-Za-z0-9_.]* <- function", l, value = TRUE))
  }))
  expect_length(unique(nm[duplicated(nm)]), 0L)
})

test_that("refresh reports progress without changing the results", {
  m <- bayescomp_model("parameters { vector[2] mu; } model { mu ~ normal(1, 1); }")
  expect_message(o <- bayescomp_optimize(m, seed = 1, refresh = 1))
  expect_message(p <- bayescomp_pathfinder(m, num_paths = 1, num_draws = 20, num_psis_draws = 20, seed = 1, refresh = 1), "Best ELBO")
  expect_message(l <- bayescomp_laplace(m, draws = 4, seed = 1, refresh = 2), "iteration: 4 / 4")
  expect_identical(p$draws, bayescomp_pathfinder(m, num_paths = 1, num_draws = 20, num_psis_draws = 20, seed = 1)$draws)
  expect_identical(l$draws, bayescomp_laplace(m, draws = 4, seed = 1)$draws)
  expect_silent(bayescomp_laplace(m, draws = 4, seed = 1))
})

test_that("propto keeps density terms of parameter dependent values", {
  d <- function(code, a, b) {
    m <- bayescomp_model(code)
    f <- function(u) bayescomp_log_prob(m, u) - bayescomp_log_prob(m, u, propto = TRUE)
    f(a) - f(b)
  }
  expect_equal(d("parameters { real x; } model { 0 ~ normal(floor(x), 1); }", 0.2, 1.2), 0)
  expect_equal(d("parameters { real x; } model { real z = 0; if (x > 0) z = 1; 0 ~ normal(z, 1); }", -1, 1), 0)
  expect_equal(d("parameters { real x; } model { target += normal_lupdf(0 | floor(x), 1); }", 0.2, 1.2), 0)
})

test_that("integer division binds and warns as in stanc", {
  code <- "transformed data { int a = 8 * 3 %/% 2; int h = 8 * 3 %/% 2 * 5; int n = 7; int b = n / 2; }"
  expect_warning(bayescomp_model(code), "Found int division")
  m <- suppressWarnings(bayescomp_model(code))
  expect_equal(c(get("V.a", m$env), get("V.h", m$env), get("V.b", m$env)), c(8L, 40L, 3L))
  expect_match(bayescomp_pedantic("transformed data { int a = 1 / 2; }")$message, "Found int division: 1 / 2")
  expect_error(bayescomp_sample(bayescomp_model("parameters { real x; } model { x ~ std_normal(); }"), iter_warmup = Inf),
               "must be an integer")
})

test_that("profiles time the reverse pass and cannot nest under one name", {
  m <- bayescomp_model('parameters { vector[50] mu; } model { profile("lp") { for (n in 1:50) target += -0.5 * mu[n]^2; } }')
  profile_reset()
  bayescomp_grad_log_prob(m, rep(0.1, 50))
  p <- profile_table()
  expect_gt(p$chain_stack, 0)
  expect_true(p$reverse_time >= 0)
  expect_error(bayescomp_model('transformed data { real s = 0; profile("p") { profile("p") { s += 1; } } }'), "already started")
  profile_reset()
})

test_that("rejections are counted by reason", {
  m <- bayescomp_model('parameters { real x; } model { if (x > 1.5) reject("too big"); x ~ normal(0, 1); }')
  f <- bayescomp_sample(m, chains = 1, iter_warmup = 100, iter_sampling = 100, seed = 1, refresh = 0)
  expect_match(f$rejections$reason, "too big")
  expect_gt(sum(f$rejections$count), 0)
  expect_equal(bayescomp_diagnose(f, quiet = TRUE)$rejections, f$rejections)
})

test_that("a nested tape leaves the outer tape intact", {
  m <- bayescomp_model("functions { vector f(vector y, real k) { return y - k; } }
    parameters { real<lower=0> k; } model { vector[1] s = solve_newton(f, [1]', k); target += -0.5 * square(s[1] - 2); }")
  g <- bayescomp_grad_log_prob(m, log(2))
  h <- 1e-6
  fd <- (bayescomp_log_prob(m, log(2) + h) - bayescomp_log_prob(m, log(2) - h)) / (2 * h)
  expect_equal(as.vector(g), fd, tolerance = 1e-5)
})

test_that("diagnose runs only the checks that apply to the method", {
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
  f <- bayescomp_sample(m, algorithm = "metropolis", chains = 2, iter_warmup = 50, iter_sampling = 50, seed = 1, refresh = 0)
  out <- capture.output(r <- bayescomp_diagnose(f))
  expect_true(any(grepl("Checks of tree depth and divergences and E-BFMI do not apply to the metropolis sampler", out)))
  expect_true(is.na(r$num_divergent) && is.na(r$num_max_treedepth))
  expect_false(any(grepl("Checking sampler transitions", out)))
  p <- suppressWarnings(bayescomp_laplace(m, refresh = 0))
  out <- capture.output(r <- bayescomp_diagnose(p))
  expect_true(any(grepl("does not draw Markov chains", out)))
  expect_false(r$problems)
  f <- bayescomp_sample(m, chains = 1, iter_warmup = 50, iter_sampling = 50, seed = 1, refresh = 0)
  out <- capture.output(bayescomp_diagnose(f))
  expect_true(any(grepl("Checking sampler transitions treedepth", out)))
  expect_false(any(grepl("do not apply", out)))
})
