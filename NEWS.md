# bayescompR (development version)

## Bug fixes

* Loops whose loop invariant expressions read `target()` or call `_lp`,
  `_rng` or `_jacobian` functions are no longer vectorized; the rewrite
  evaluated them once instead of once per iteration and changed the density.
* `cov_matrix` and `cholesky_factor_cov` parameters no longer give a NaN
  gradient for large off diagonal values.
* Metropolis and slice sampling work with the native Cranelift and LLVM
  backends; they evaluated every point as -Inf.

## Performance

* MALA and static HMC reuse the density and gradient of the retained state
  (one evaluation per MALA transition instead of two), and NUTS computes the
  inverse metric product once per leaf.
* ADVI, flow VI and Laplace evaluate value only log densities without a
  reverse sweep.
* Fused densities compute partials only for differentiable arguments.
* `mdivide_right_tri_low()` uses triangular solves.
* The negative binomial precision partial no longer sums over large counts
  whose result is discarded.
* Cholesky factors are kept for up to four fixed covariance matrices.
* Covariance transforms build their index layout once per shape.

# bayescompR 0.1.0

First release.

* A pure R implementation of the Stan modeling language (version 2.40):
  parser, type checker and compiler for every program block, user defined
  functions and constrained type, with the Stan math and distribution
  library and reverse mode automatic differentiation. No C++ toolchain,
  CmdStan or rstan is needed.
* CmdStan's inference methods with its adaptation and defaults: NUTS,
  static HMC, optimization, Laplace approximation, ADVI, Pathfinder and
  standalone generated quantities. CmdStan's CSV, metric, data and
  configuration files are read and written, and the stanc tools
  (formatting, pedantic warnings, model information, generated data) and
  `stansummary` and `diagnose` have R equivalents.
* Further samplers and approximations: MALA, Metropolis, slice and
  elliptical slice sampling, parallel tempering, tempered SMC, nested
  sampling, ensemble samplers, Gibbs sampling, Riemannian HMC, MCLMC,
  piecewise deterministic samplers, stochastic gradient MCMC, reversible
  jump, particle MCMC for state space models, normalizing flow VI,
  expectation propagation, INLA and neural posterior estimation.
* Model checking: convergence diagnostics, divergence locations,
  pointwise log likelihood and posterior predictive draws from the model
  block, exact refits for LOO, and simulation based calibration. Runs can
  be interrupted and resumed.
* Works with posterior, bayesplot, loo and rstantools, and runs the Stan
  programs written by brms, rstanarm and multinma.
* A typed, versioned intermediate representation and a backend registry
  for native backends; `bayescomp_install_backend()` downloads prebuilt
  Cranelift and LLVM backends.
