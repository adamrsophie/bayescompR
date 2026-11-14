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

## Releasing to CRAN

* Statistical tests longer than about two seconds call `skip_on_cran()`, and
  vignette code runs only when `NOT_CRAN=true`. The CI job named `cran` checks
  the package the way CRAN does.
* Build the submission tarball with every vignette evaluated:
  `NOT_CRAN=true R CMD build .`
* Run `R CMD check --as-cran` on the tarball, the R-hub workflow
  (`gh workflow run rhub.yaml -f config=linux,macos,windows`) and
  `devtools::check_win_devel()`, then update `cran-comments.md`.
