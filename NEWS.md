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
  for native backends.
