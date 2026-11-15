# bayescompR <img src="man/figures/logo.png" align="right" height="139" alt="bayescompR logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/adamrsophie/bayescompR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/adamrsophie/bayescompR/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/adamrsophie/bayescompR/actions/workflows/pkgdown.yaml/badge.svg)](https://adamrsophie.github.io/bayescompR/)
[![test-coverage](https://github.com/adamrsophie/bayescompR/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/adamrsophie/bayescompR/actions/workflows/test-coverage.yaml)
[![R-hub](https://github.com/adamrsophie/bayescompR/actions/workflows/rhub.yaml/badge.svg)](https://github.com/adamrsophie/bayescompR/actions/workflows/rhub.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R >= 4.3.0](https://img.shields.io/badge/R-%E2%89%A5%204.3.0-276DC3.svg)](https://www.r-project.org)
<!-- CRAN badge: uncomment once the package is on CRAN. -->
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/bayescompR)](https://CRAN.R-project.org/package=bayescompR) -->
<!-- DOI badge: Zenodo mints a concept DOI on the first GitHub release; replace XXXXXXXX and uncomment. -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXXX) -->
<!-- Software Heritage archives public repositories only; uncomment once the repository is public. -->
<!-- [![Archived on Software Heritage](https://archive.softwareheritage.org/badge/origin/https://github.com/adamrsophie/bayescompR/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/adamrsophie/bayescompR) -->
<!-- badges: end -->

**bayescompR** runs Stan programs in pure R. It implements the Stan
modeling language (version 2.40): a parser and type checker, a compiler to
R closures over a reverse mode automatic differentiation tape, the full
math and distribution library with gradients, and Stan's inference
algorithms with CmdStan's adaptation, output files and diagnostics.

There is nothing to compile. Models start in a fraction of a second, run
wherever R runs, and need no C++ toolchain, CmdStan or rstan.

bayescompR is an independent implementation. It is not affiliated with
or endorsed by the Stan Development Team or NumFOCUS; "Stan" refers to the
language specification it implements.

## Installation

```r
# install.packages("remotes")
remotes::install_github("adamrsophie/bayescompR")
```

## Quick start

```r
library(bayescompR)
options(mc.cores = parallel::detectCores())  # run chains in parallel

code <- "
data {
  int<lower=0> J;
  array[J] real y;
  array[J] real<lower=0> sigma;
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[J] eta;
}
transformed parameters {
  vector[J] theta = mu + tau * eta;
}
model {
  mu ~ normal(0, 5);
  tau ~ cauchy(0, 5);
  eta ~ std_normal();
  y ~ normal(theta, sigma);
}
"
schools <- list(J = 8, y = c(28, 8, -3, 7, -1, 1, 18, 12),
                sigma = c(15, 10, 16, 11, 9, 11, 10, 18))

model <- bayescomp_model(code, data = schools)
fit <- bayescomp_sample(model, chains = 4, seed = 1)

summary(fit, c("mu", "tau", "theta"))
bayescomp_diagnose(fit)
loo::loo(fit)                       # pointwise terms come from the model block
yrep <- bayescomp_posterior_predict(fit)
```

## Features

**The Stan language.** All program blocks; user defined functions,
including recursion and `_lpdf`, `_lp`, `_rng` and `_jacobian` functions;
every constrained type (simplexes, `sum_to_zero_vector`, stochastic
matrices, correlation and covariance matrices and their Cholesky factors,
`offset` and `multiplier`); tuples and complex numbers; truncation;
`profile` blocks; ODE, DAE, algebraic solver, 1D integration, `map_rect`
and `reduce_sum`; and every function and distribution of the Stan 2.40
library, with gradients.

**CmdStan parity.** NUTS with Stan's windowed adaptation of the step size
and a diagonal or dense metric, static HMC, optimization (L-BFGS, BFGS,
Newton), Laplace approximation, Pathfinder, ADVI and standalone generated
quantities. Stan CSV, metric, diagnostic, profile and configuration files
are written and read in CmdStan's formats, so `bayescomp_read_csv()` also
reads CmdStan output. `bayescomp_stansummary()` and `bayescomp_diagnose()`
reproduce `stansummary` and `diagnose`, and the stanc tools are there too:
`bayescomp_format()`, `bayescomp_info()`, `bayescomp_pedantic()`,
`bayescomp_generate_data()` and `bayescomp_expose_functions()`.

**Beyond CmdStan.**

* Pointwise log likelihood and posterior predictive draws taken from the
  model block (`bayescomp_log_lik()`, `bayescomp_posterior_predict()`), so
  `loo` and predictive checks need no generated quantities.
* Exact refits for observations with large Pareto k (`bayescomp_reloo()`)
  and simulation based calibration (`bayescomp_sbc()`).
* Where divergent trajectories failed, rejected proposals counted by reason,
  and diagnostics that report only the checks that apply to the method.
* Interrupt with Ctrl-C and keep the draws; continue chains with
  `bayescomp_sample(resume = fit)`.
* More samplers: MALA, Metropolis, slice and elliptical slice sampling,
  parallel tempering, tempered SMC, nested sampling, ensemble samplers,
  Gibbs sampling with conjugacy detection, Riemannian HMC, MCLMC,
  piecewise deterministic samplers, SG-MCMC, reversible jump, particle
  MCMC for state space models, and approximations including normalizing
  flow VI, EP, INLA and neural posterior estimation.
* Direct access to the log density, its gradient and the transforms
  (`bayescomp_log_prob()`, `bayescomp_grad_log_prob()`,
  `bayescomp_constrain()`), for use with other software.

**The ecosystem.** Fits work with posterior, bayesplot, loo and
rstantools generics. Programs written by brms (`make_stancode()` and
`make_standata()`), rstanarm and multinma run unchanged; see
`vignette("ecosystem")`.

## Performance

Everything runs in the R interpreter, so bayescompR is slower than
compiled Stan. Its compiler works to keep the gap small:

* Loops whose iterations are independent are rewritten into vectorized
  statements, so models written observation by observation (as brms
  writes them) run 10 to 30 times faster than the loops would.
* Densities, GLM likelihoods and constraining transforms are single tape
  nodes with hand written partial derivatives, and data only expressions
  are computed once.
* Chains run in parallel with `cores` or `options(mc.cores = ...)`, with
  results identical to sequential runs for the same seed.

A gradient of a small hierarchical model costs well under a millisecond.
The typed, versioned IR (`bayescomp_ir()`) and the backend registry let
native backends compile the same model to machine code.
`bayescomp_install_backend("cranelift")` (or `"llvm"`) downloads a
prebuilt backend library, so no compiler is needed; see
`vignette("backends")` and `vignette("performance")`.

## Documentation

The website <https://adamrsophie.github.io/bayescompR/> has the reference
and these vignettes:

| Vignette | Contents |
|---|---|
| `vignette("bayescompR")` | Getting started: the full workflow on one model |
| `vignette("stan-language")` | What the implementation supports, and how it differs from stanc |
| `vignette("inference")` | NUTS and every argument of `bayescomp_sample()`, HMC, MALA, Metropolis, slice, SMC, optimization, Laplace, ADVI, Pathfinder |
| `vignette("algorithms")` | The other samplers and approximations |
| `vignette("diagnostics")` | R-hat, ESS, divergences and where they occur, rejections |
| `vignette("cmdstan-workflow")` | The stanc and CmdStan tools: files, CSV, stansummary, diagnose |
| `vignette("ecosystem")` | posterior, bayesplot, loo, rstantools, brms, rstanarm, multinma |
| `vignette("performance")` | Writing fast models; parallel chains |
| `vignette("autodiff-and-transforms")` | The AD tape and the constraining transforms |
| `vignette("backends")` | The IR and the interface for native backends |
| `vignette("validation")` | How correctness is checked, including simulation based calibration |

## Citation

Run `citation("bayescompR")` for the package citation and the Stan
reference. Machine-readable metadata is in `CITATION.cff`, `codemeta.json`,
`.zenodo.json` and `ro-crate-metadata.json`.

## Contributing

See `CONTRIBUTING.md`. Report bugs and request features in the
[issue tracker](https://github.com/adamrsophie/bayescompR/issues), and
security problems as described in `SECURITY.md`. The project follows a
[Code of Conduct](CODE_OF_CONDUCT.md).

## License

GPL-3. Data and texts derived from stanc3 and CmdStan are under the BSD
3-Clause License; see `inst/COPYRIGHTS`.
