# bayescompR 0.1.0

First release.

* Lexer, parser and semantic analyzer for Stan 2.40 programs, with
  classed errors that report line and column.
* Code generation into R closures over a reverse mode automatic
  differentiation tape with fused density kernels.
* All constrained types of Stan 2.40 with Stan's transforms and Jacobians.
* Math, container, linear algebra and distribution library.
* The `wiener` distribution with four, five and seven parameters and its
  defective distribution functions.
* `rstantools::log_lik()` and `rstantools::posterior_predict()` methods.
* Exact gradients of `svd_U()`, `svd_V()`, `qr_Q()` and `qr_R()` of tall
  matrices, and of `student_t_qf()` in the degrees of freedom.
* Embedded Laplace approximation (`laplace_marginal` and related
  functions, including the built in Poisson, Bernoulli and negative
  binomial forms).
* ODE, DAE and algebraic solvers, `integrate_1d`, `map_rect` and
  `reduce_sum`, differentiated with respect to all their parameters.
* Tuples in every block, including constrained tuple parameters, with
  projection, unpacking and the library functions that return tuples.
* Complex numbers: complex scalars, vectors, row vectors and matrices with
  differentiable arithmetic, elementwise functions, FFT, eigenvalues, SVD and
  Schur decompositions, and the transform, HMM and further distribution
  functions of Stan 2.40.
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
* `inv_erfc()`; `hypergeometric_2F1()` for all `z < 1`.
* Gibbs conjugate blocks are always Metropolis corrected, and the slice
  kernel splits its stepping out budget at random, so both are exact.
* SG-MCMC scales minibatch gradients by the number of batches, which is
  unbiased when batch sizes differ.
* Overloaded user functions resolve to the signature with the fewest
  promotions; ties are errors. Model block variables are not visible in
  generated quantities, and data only library arguments are enforced.
* Transformed data sizes may depend on earlier transformed data; array
  expressions keep singleton vector elements, and elementwise operations
  on containers of different sizes are shape errors.
* Cumulative distribution functions check their parameters like the
  densities; truncation at the edge of the support and extreme tail CDF
  gradients are finite; `log_diff_exp(a, negative_infinity())` is `a`.
* Repeated indices in an assignment give the adjoint to the last write.
* R-hat, bulk and tail ESS, basic ESS and MCSE of the mean and SD equal
  those of posterior 1.7.1 (tested when posterior is installed), including
  ESS at lag zero; the constant check is relative to the draws' scale.
* Nested sampling reports `converged` and computes the information from the
  dead and live points; incomplete SMC sets `complete = FALSE` and
  `log_evidence = NA`; a failed Newton line search is not reported as
  converged, and Laplace warns when the optimizer did not converge.
* Diffusion time grids end at a positive time and diffusion fails rather
  than return draws outside the support; MCLMC keeps its step size fixed
  while sampling; ADVI drops draws with zero density from the ELBO as Stan does.
* Zig-Zag and BPS compare refresh and bounce times; particle filter
  increments ignore zero weight particles; vector parameters work in PMCMC
  and neural posterior estimation.
* `bayescomp_generate_quantities()` checks that the parameter declarations
  match the fit; impossible infinite bounds and nonpositive multipliers are
  domain errors; slice, Metropolis, ensemble, tempering and PMCMC start and
  run without a gradient.
* `bayescomp_compatibility()` returns the machine readable feature matrix
  (implemented, partial, not implemented, not applicable) with the test
  files covering each feature.
