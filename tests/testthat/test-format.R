test_that("the formatter reproduces stanc --auto-format", {
  code <- "// lead
functions {
  real foo(real x, data array[] int k) { return x*sum(k); }  // trailing
}
data { int<lower=0> N; vector[N] y; array[N, 2] int<lower=0, upper=1> k; real a, b; }
parameters { real mu; real<lower=0> sigma; vector<offset=mu, multiplier=sigma>[N] eta; }
model {
  real m = (mu + 1) * 2;
  real s2;
  if (mu > 0) m = 1; else if (mu < -1) { m = 2; } else m = 3;
  for (i in 1:N) { y[i] ~ normal(mu, sigma) T[0, ]; }
  y[:3] ~ normal(mu, sigma);
  target += normal_lupdf(y | mu, sigma) + foo(mu, k[:, 1]) + log_sum_exp(mu, sigma) * fma(mu, sigma, 1.0e-3) * log_mix(0.5, mu, sigma);
  profile(\"p\") { mu ~ std_normal(); }
  (m, s2) = (1, 2.0);
  m = -mu^2 + get_real(3i) * 0;
}"
  out <- bayescomp_format(code)
  expect_identical(strsplit(out, "\n")[[1]], c(
    "// lead", "functions {", "  real foo(real x, data array[] int k) {", "    return x * sum(k);", "  } // trailing", "}",
    "data {", "  int<lower=0> N;", "  vector[N] y;", "  array[N, 2] int<lower=0, upper=1> k;", "  real a, b;", "}",
    "parameters {", "  real mu;", "  real<lower=0> sigma;", "  vector<offset=mu, multiplier=sigma>[N] eta;", "}",
    "model {", "  real m = (mu + 1) * 2;", "  real s2;", "  if (mu > 0) ", "    m = 1;", "  else if (mu < -1) {", "    m = 2;",
    "  } else ", "    m = 3;", "  for (i in 1 : N) {", "    y[i] ~ normal(mu, sigma) T[0, ];", "  }", "  y[ : 3] ~ normal(mu, sigma);",
    "  target += normal_lupdf(y | mu, sigma) + foo(mu, k[ : , 1])",
    "            + log_sum_exp(mu, sigma) * fma(mu, sigma, 1.0e-3)",
    "              * log_mix(0.5, mu, sigma);",
    "  profile(\"p\") {", "    mu ~ std_normal();", "  }", "  (m, s2) = (1, 2.0);", "  m = -mu ^ 2 + get_real(3i) * 0;", "}"))
  # the output parses to the same program
  expect_identical(bayescomp_format(out), out)
  expect_output(print(out), "functions \\{")
})

test_that("canonicalize options and include lines", {
  b <- bayescomp_format("parameters { real a; } model { if (a > 0) a ~ normal((0), 1); // c\n }", canonicalize = c("braces", "parentheses", "strip-comments"))
  expect_identical(as.character(b), "parameters {\n  real a;\n}\nmodel {\n  if (a > 0) {\n    a ~ normal(0, 1);\n  }\n}\n")
  d <- tempfile(); dir.create(d)
  writeLines("real twice(real x) { return 2 * x; }", file.path(d, "tw.stan"))
  src <- "functions {\n#include tw.stan\n}"
  expect_match(bayescomp_format(src), "#include tw.stan", fixed = TRUE)
  expect_match(bayescomp_format(src, canonicalize = TRUE, include_paths = d), "real twice(real x) {", fixed = TRUE)
  expect_error(bayescomp_format(src, canonicalize = "bogus"), "unknown canonicalization")
})
