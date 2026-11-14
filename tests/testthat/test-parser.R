test_that("programs parse into the expected blocks", {
  p <- stan_parse("
    functions { real f(real x) { return 2 * x; } }
    data { int<lower=0> N; vector[N] y; }
    transformed data { real m = mean(y); }
    parameters { real mu; real<lower=0> sigma; }
    transformed parameters { real z = mu / sigma; }
    model { y ~ normal(mu, sigma); }
    generated quantities { real yrep = normal_rng(mu, sigma); }")
  expect_equal(names(p$blocks), stan_block_names)
  expect_equal(p$blocks$data$body[[2L]]$name, "y")
})

test_that("operator precedence and associativity follow Stan", {
  e <- stan_parse("transformed data { real x = 1 + 2 * 3 ^ 2 ^ 2; }")$blocks[["transformed data"]]$body[[1L]]$init
  expect_equal(e$op, "+")
  expect_equal(e$b$op, "*")
  expect_equal(e$b$b$op, "^")
  expect_equal(e$b$b$b$op, "^")  # right associative
  m <- bayescomp_model("transformed data { real x = -2 ^ 2; int a = 7 %/% 2; int b = -7 %/% 2; int c = 7 % 3; }
                        generated quantities { real xx = x; int aa = a; int bb = b; int cc = c; }")
  v <- bayescomp_constrain(m, numeric(0))
  expect_equal(v$xx, -4)
  expect_equal(c(v$aa, v$bb, v$cc), c(3L, -3L, 1L))
  # as in stanc3: left division binds tighter than *, and .* shares the
  # level of * (left associative), so v * 2 .* v is (v * 2) .* v
  e <- stan_parse("transformed data { matrix[2, 2] A; vector[2] v; vector[2] w = A * A \\ v; vector[2] u = v * 2 .* v; }")$blocks[["transformed data"]]$body
  expect_equal(e[[3L]]$init$op, "*")
  expect_equal(e[[3L]]$init$b$op, "\\")
  expect_equal(e[[4L]]$init$op, ".*")
  expect_equal(e[[4L]]$init$a$op, "*")
})

test_that("syntax errors carry the right class and location", {
  expect_error(stan_parse("parameters { real x }"), class = "bayescomp_syntax_error")
  expect_error(stan_parse("model { } data { }"), "out of order", class = "bayescomp_syntax_error")
  expect_error(stan_parse("data { real y[3]; }"), "old array syntax", class = "bayescomp_syntax_error")
  expect_error(stan_parse("model { # comment\n }"), class = "bayescomp_syntax_error")
  err <- tryCatch(stan_parse("parameters {\n  real x\n}"), error = identity)
  expect_match(conditionMessage(err), "line")
})

test_that("newer language features check", {
  expect_silent(bayescomp_check("parameters { real<lower=0> y; } model { 0.5 ~ wiener(1, y, 0.5, 0.2); }"))
  # stanc rejects complex parameters
  expect_silent(bayescomp_check("parameters { complex z; }"))
  expect_silent(bayescomp_check("data { tuple(real, int) x; }"))
})

test_that("comments, scientific literals and truncation parse", {
  p <- stan_parse("/* block */ data { real a; } // line\nparameters { real<lower=0> s; }
                   model { s ~ normal(0, 1e1) T[0, ]; a ~ normal(.5, 2.5E-1); }")
  st <- p$blocks$model$body[[1L]]
  expect_equal(st$kind, "tilde")
  expect_false(is.null(st$trunc$lo))
  expect_null(st$trunc$hi)
})
