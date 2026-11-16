# bayescompR (development version)

## Bug fixes

* Loops whose loop invariant expressions read `target()` or call `_lp`,
  `_rng` or `_jacobian` functions are no longer vectorized; the rewrite
  evaluated them once instead of once per iteration and changed the density.
* `cov_matrix` and `cholesky_factor_cov` parameters no longer give a NaN
  gradient for large off diagonal values.
* Metropolis and slice sampling work with the native Cranelift and LLVM
  backends; they evaluated every point as -Inf.
* `cholesky_factor_corr` and `corr_matrix` parameters have a finite log
  Jacobian where `tanh()` of the unconstrained value rounds to 1; it was
  -Inf with a NaN gradient.
* GLM densities reject an intercept vector whose length is neither 1 nor the
  number of rows of `x`; it was recycled.

## New features

* `bayescomp_sgmcmc()` gains `lp_every`, so the exact full data `lp__` can
  be computed at only some saved draws or not at all.

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
* Multivariate normal and Student t densities transpose arrays of vectors
  whole instead of splitting them into rows; a shared mean vector is
  broadcast by one operation.
* `cholesky_factor_corr` uses one node for the factor and one for the log
  Jacobian instead of a scalar graph.
* GLM linear predictors are one node, and `normal_id_glm` with fixed data, a
  scalar intercept and a scalar scale uses a QR factorization prepared once,
  so each evaluation costs O(P^2) instead of O(N P).
* SG-MCMC initialization and expectation propagation reuse the joint value
  and gradient that the optimizer asks for in two calls.
* RMHMC carries the metric derivatives at the end of one leapfrog step into
  the next.
* NaN checks no longer allocate, and pointwise log likelihoods and
  predictive draws keep only the requested array.

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
