# Contributing

* Open an issue before large changes.
* Tests run on GitHub Actions; `BAYESCOMPR_LONG_TESTS=true` enables the
  posterior validation and simulation based calibration tests.
* New library functions need value tests against base R or a reference
  formula and gradient tests against finite differences
  (`expect_grad_ok()` in `tests/testthat/helper.R`).
* New inference algorithms need moment tests on analytic targets and, for
  approximate methods, accuracy tests against analytic posteriors.
* Style: American English; comments and documentation do not use dashes as
  sentence punctuation.
* Commits on `main` for 0.1.0; feature branches and pull requests afterward.
