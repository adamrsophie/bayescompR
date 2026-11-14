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
