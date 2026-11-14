test_that("type errors are reported as semantic errors", {
  se <- function(code) expect_error(bayescomp_check(code), class = "bayescomp_semantic_error")
  se("parameters { real x; } model { x = 1; }")                      # assign to parameter
  se("data { int N; } parameters { vector[N] x; } model { x ~ normal(0, 1); real y = x; }")
  se("model { real x; x ~ foo(1); }")
  se("parameters { real x; } model { x ~ normal(0); }")
  se("parameters { real x; } transformed parameters { real y = normal_rng(0, 1); } model { }")
  se("parameters { real x; } generated quantities { target += 1; }")
  se("data { real y; } parameters { real x; } model { y ~ normal(x, 1) T[0, 1]; int k = 1.5; }")
  se("functions { real f(real x) { if (x > 0) return 1; } } model { }")
  se("parameters { real x; } model { break; }")
  se("parameters { real x; real x; } model { }")
  se("parameters { int k; } model { }")
  se("parameters { real x; } model { target += normal_lpdf(x, 0, 1); }")
  se("data { vector[3] v; } model { real s = v[1, 2]; }")
})

test_that("valid programs pass checking", {
  expect_silent(bayescomp_check("
    functions {
      real my_lpdf(real y, real mu) { return normal_lupdf(y | mu, 1); }
      void hello_lp(real x) { target += -x^2; }
      real twice(data real x) { return 2 * x; }
    }
    data { int N; array[N] int<lower=0, upper=1> y; matrix[N, 2] X; }
    transformed data { real t = twice(1.0); }
    parameters { vector[2] b; real<lower=0> s; }
    model {
      s ~ my(0);
      hello_lp(s);
      y ~ bernoulli_logit(X * b);
      for (i in 1:N) { if (y[i] == 1) continue; }
      while (0) break;
    }
    generated quantities { array[N] int yrep = bernoulli_logit_rng(X * b); }"))
})

test_that("rng and lp functions respect their contexts", {
  expect_error(bayescomp_check("parameters { real x; } model { real z = normal_rng(0, 1); }"), class = "bayescomp_semantic_error")
  expect_error(bayescomp_check("functions { void f_lp() { target += 1; } } generated quantities { f_lp(); }"), class = "bayescomp_semantic_error")
  expect_error(bayescomp_check("parameters { real x; } generated quantities { real y = x; y ~ normal(0, 1); }"), class = "bayescomp_semantic_error")
})
