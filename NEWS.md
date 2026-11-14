# bayescompR 0.1.0

First release.

* Lexer, parser and semantic analyzer for Stan 2.40 programs, with
  classed errors that report line and column.
* Code generation into R closures over a reverse mode automatic
  differentiation tape with fused density kernels.
* All constrained types of Stan 2.40 with Stan's transforms and Jacobians.
* Math, container, linear algebra and distribution library.
* NUTS (multinomial, generalized U-turn) and static HMC with dual averaging
  and windowed diagonal or dense metric adaptation.
* MALA, block and componentwise random walk Metropolis with adaptive
  covariance during warmup, and slice sampling with stepping out, all
  sharing Stan's warmup schedule.
* Adaptive tempered sequential Monte Carlo with a Gaussian reference and a
  log evidence estimate.
* Pathfinder and ADVI fits can initialize MCMC chains.
* SMC rejuvenation by Metropolis, MALA or HMC, with systematic, stratified,
  residual or multinomial resampling.
* Generalized elliptical slice sampling, non reversible parallel tempering
  with schedule adaptation, differential evolution and stretch move
  ensembles, blocked Gibbs with numerical conjugacy detection.
* Simplified manifold MALA and RMHMC with the SoftAbs metric, microcanonical
  Langevin Monte Carlo, Zig-Zag and bouncy particle samplers.
* State space models with particle filters, PMMH and particle Gibbs with
  ancestor sampling.
* SGLD and SGHMC, expectation propagation, nested sampling, reversible jump
  MCMC, INLA, planar normalizing flow VI, neural posterior estimation and
  reverse diffusion Monte Carlo.
* Rank histograms and local, quantile and folded ESS plots.
* Optimization (L-BFGS, BFGS, Newton), Laplace approximation, ADVI (mean
  field and full rank) and Pathfinder with Pareto smoothed importance
  resampling.
* Standalone generated quantities.
* Rank normalized R-hat, bulk and tail ESS and MCSE matching the posterior
  package.
* Typed, versioned intermediate representation (IR version 1.0.0) and a
  backend registry for native code backends. Missing backend capabilities
  require `allow_fallback = TRUE`; `bayescomp_validate_backend()` compares
  a native backend with the R implementation.
* `bayescomp_native_inputs()` gives a native backend the IR and the data and
  transformed data values in Stan order; `bayescomp_native_error()` turns
  native error codes into the package's condition classes. Only the log
  density and gradient are taken from a native backend in this version, and
  `provenance` says so.
* Operator precedence follows Stan: `\` binds tighter than `*`, and `.*`,
  `./` tighter still.
* Without `seed`, a seed is drawn from the global stream (so repeated calls
  differ and `set.seed()` reproduces them) and recorded in `fit$seed`.
* Bounds may mix finite and infinite values elementwise.
* `inv_erfc()`; `hypergeometric_2F1()` for all `z < 1`. Stan functions that
  are not implemented yet raise `bayescomp_unsupported_error`.
