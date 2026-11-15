## Submission

This is the first submission of bayescompR.

## Test environments

* GitHub Actions: Ubuntu (R release and R devel), macOS (R release),
  Windows (R release)
* GitHub Actions, Ubuntu R devel with `NOT_CRAN` unset, as CRAN runs the checks
* R-hub: Linux, macOS, Windows
* win-builder: R devel

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

## Notes for the reviewers

* The package is an independent implementation of the Stan
  modeling language in R. It does not call or bundle Stan. The Description states that
  it is not affiliated with or endorsed by the Stan Development Team.
* Statistical validation tests (long Markov chain runs) and the vignette code
  are skipped on CRAN to keep the check time short; they run on our continuous
  integration with `NOT_CRAN=true`. The vignettes in the tarball were built
  with all code evaluated.
* Examples that run samplers for more than a few seconds are in `\donttest{}`.
