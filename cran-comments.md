## Resubmission

This is a resubmission. In this version I have:

* Replaced the relative link to CODE_OF_CONDUCT.md in README.md, a file not
  included in the package, with its full URL.

## Submission

This is the first submission of bayescompR.

## Test environments

* GitHub Actions: Ubuntu (R release and R devel), macOS (R release),
  Windows (R release)
* GitHub Actions, Ubuntu R devel with `NOT_CRAN` unset, as CRAN runs the checks
* R-hub: Linux, macOS, Windows
* win-builder: R devel (Windows, 1 note: new submission and the words below)
* local: macOS (R release), `R CMD check --as-cran --run-donttest`

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.
* Possibly misspelled words in DESCRIPTION: Betancourt, Gelman, Vehtari and
  Zhang are author names of the cited references; Langevin and Riemannian are
  the names of the sampling methods (Langevin dynamics, Riemannian manifold
  Hamiltonian Monte Carlo).

## Notes for the reviewers

* The package is an independent implementation of the Stan
  modeling language in R. It does not call or bundle Stan. The Description states that
  it is not affiliated with or endorsed by the Stan Development Team.
* Statistical validation tests (long Markov chain runs) and the vignette code
  are skipped on CRAN to keep the check time short; they run on our continuous
  integration with `NOT_CRAN=true`. The vignettes in the tarball were built
  with all code evaluated.
* Examples that run samplers for more than a few seconds are in `\donttest{}`.
* `bayescomp_install_backend()` downloads an optional prebuilt shared
  library (a native backend) from a fixed release of the package's GitHub
  repository over https, but only when the user calls it. We are aware that
  the CRAN policy allows downloading precompiled code only as a last resort
  and with the agreement of the CRAN team, and we ask for that agreement
  here: the libraries are compiled Rust (Cranelift or LLVM code generators)
  that would otherwise require a Rust toolchain and, for LLVM, an LLVM
  installation on the user's machine. The package works fully without them.
  The file is checked against its SHA-256 checksum, recorded in the package
  source (R/native.R), before it is loaded (on R
  older than 4.5, which has no `tools::sha256sum()`, the function refuses to
  download), written only to
  `tools::R_user_dir("bayescompR", "data")`, which the package manages
  (installing removes the files of other releases), and loaded with
  `dyn.load()`. Nothing is downloaded at installation, loading, or in
  examples, tests or vignettes, and the package has no compiled code of its
  own. If this is not acceptable, we will remove the installer.
