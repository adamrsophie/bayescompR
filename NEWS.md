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
* Parallel chains on a socket cluster (Windows) run sequentially, with a
  message, when the worker processes cannot load an installed bayescompR,
  and see arguments such as `init = fit` correctly (they were evaluated in
  the worker, where `pf` meant `stats::pf()`).
* `bayescomp_diagnose()` runs only the checks that apply: tree depth for
  NUTS, divergences and E-BFMI for NUTS and static HMC, and Pareto k for
  approximate methods (it failed on Pathfinder and Laplace fits).
* `bayescomp_sbc()` runs simulation based calibration, with rank
  uniformity tests and rank histograms (`plot()`).
* `bayescomp_reloo()` refits the model without each observation whose
  PSIS-LOO Pareto k is too large and updates the LOO estimates with its
  exact log predictive density.
* An interrupt (Ctrl-C or Esc) while chains run sequentially returns the
  finished chains, or the draws of the interrupted chain, instead of
  losing them; `bayescomp_sample(resume = fit)` continues each chain from
  its last draw with its adapted step size and metric.
* Fits from `bayescomp_sample()` record where divergent trajectories left
  the typical set (element `divergences`: chain, iteration and parameter
  values).
* `bayescomp_log_lik()` and `bayescomp_posterior_predict()` take the
  pointwise log likelihood and predictive draws from the sampling
  statements of the model block, so `loo::loo()`, `rstantools::log_lik()`
  and `rstantools::posterior_predict()` work on models without a
  generated quantities block.
* Loops whose iterations only read and write their own element (sampling
  statements, `target +=` of a density and assignments to `v[n]`) run as
  one statement on slices: 13 to 34 times faster gradients on brms models
  and up to 300 times on scalar regression loops, with the same density
  and gradient. Loops that do not qualify, or whose rewrite does not check,
  stay loops; `options(bayescompR.vectorize = FALSE)` turns the pass off.
* Scalar arithmetic, `row_vector * vector` and matrix product adjoints
  skip generic dispatch and transposes (1.4 to 1.8 times faster gradients
  on small models).
* Correlation matrix and Cholesky correlation transforms, and lower and
  upper bounds that are data, run as single nodes with hand written
  adjoints (a model with a 3 x 3 correlation matrix: 1.9 to 1.1 ms per
  gradient).
* Element writes in loops that cannot be vectorized (`mu[n] = r * mu[n - 1]`)
  use linear rather than quadratic memory; multivariate normal densities
  factor a data covariance or precision matrix once.
* Faster gradients: scalar arithmetic and single indexes skip the generic
  machinery, selections accumulate adjoints in place, sum to zero and
  simplex transforms run in linear time, multivariate densities and ordered
  and categorical likelihoods evaluate all observations at once, Pathfinder
  uses the compact L-BFGS form, and `bayescomp_log_prob()` skips the
  gradient.
* `rstantools::log_lik()` and `rstantools::posterior_predict()` methods.
* Data, initial values and metrics from CmdStan JSON and R dump files
  (`bayescomp_read_data()`, `bayescomp_write_data()`; any `data` or `init`
  argument also takes a file name).
* Stan CSV output and input for every method (`bayescomp_write_csv()`,
  `bayescomp_read_csv()`), readable by CmdStan's `stansummary` and by
  cmdstanr, including the diagnostic file (`save_diagnostics = TRUE`), the
  adapted metric (`save_metric = TRUE`) and single Pathfinder paths.
* `bayescomp_check_gradients()`, the gradient test of CmdStan's `diagnose`
  method. `bayescomp_log_prob()` and `bayescomp_grad_log_prob()` take
  constrained values, files and matrices of points.
* Pathfinder takes all of CmdStan's arguments (`num_psis_draws`,
  `calculate_lp`, `save_single_paths`, `init_alpha` and the L-BFGS
  tolerances) and records the path of each draw; Laplace takes its mode from
  a file and has `calculate_lp`.
* stanc's compiler options have R equivalents: `bayescomp_format()`
  (`--auto-format`, `--print-canonical`, `--canonicalize`,
  `--max-line-length`; comments and `#include` lines are kept),
  `bayescomp_info()` (`--info`), `bayescomp_pedantic()` and
  `bayescomp_check(pedantic = TRUE)` (`--warn-pedantic`, with stanc's
  messages), `bayescomp_uninitialized()` (`--warn-uninitialized`),
  `bayescomp_expose_functions()` (`--standalone-functions`),
  `bayescomp_generate_data()` and `bayescomp_generate_inits()`
  (`--debug-generate-data`, `--debug-generate-inits`) and
  `bayescomp_stan_version()`.
* Subexpressions that depend only on data are computed once when the data
  are bound instead of in every evaluation of the log density.
* `print()` writes numbers as CmdStan does (six significant digits), and
  `bayescomp_write_csv(save_config = TRUE)` writes CmdStan's configuration
  JSON.
* `bayescomp_optimize(save_iterations = TRUE)` keeps every iterate.
  `bayescomp_generate_quantities()` also takes CSV files, posterior draws
  objects, arrays and matrices of constrained draws.
* `profile()` blocks are timed (element `profile` of a fit, written as
  CmdStan's profile file). `bayescomp_stansummary()` reproduces CmdStan's
  `stansummary` (percentiles, significant figures, autocorrelations, CSV
  output) and `bayescomp_diagnose()` CmdStan's `diagnose` report; both take
  CSV files.
* Library calls are checked against the signatures of stanc 2.40, so calls
  stanc rejects (such as `sum()` of a two dimensional array) are rejected.
* Exact gradients of `svd_U()`, `svd_V()`, `qr_Q()` and `qr_R()` of tall
  matrices, and of `student_t_qf()` in the degrees of freedom.
* Embedded Laplace approximation (`laplace_marginal` and related
  functions, including the built in Poisson, Bernoulli and negative
  binomial forms).
* ODE, DAE and algebraic solvers (including the deprecated
  `integrate_ode_*` interface), the 1D integrators (`integrate_1d`,
  `integrate_1d_double_exponential`, `integrate_1d_gauss_kronrod` and their
  `_tol` forms), `map_rect` and `reduce_sum`, differentiated with respect to
  all their parameters.
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
* Complex parameters (scalars, vectors, row vectors and matrices) as in
  stanc 2.40; tuple data can be generated and multidimensional arrays of
  tuples are read back from CSV files.
* `#include` directives inside block comments are ignored, and errors in
  included files name the file and line. Sampling statements must name a
  distribution stanc knows; `T[ , ]` is no truncation and truncation bounds
  must be scalars.
* Transformed data with `CHAIN_ID` is the same in every chain, run and
  number of cores; seeded calls restore the random number generator.
* `bayescomp_stansummary()` and `bayescomp_diagnose()` compute MCSE, ESS and
  R-hat as `bin/stansummary` and `bin/diagnose` do. The density columns of
  ADVI, Pathfinder and Laplace output (`log_p__`, `log_g__`, `log_q__`)
  follow CmdStan's conventions, and the CSV header uses CmdStan's layout.
* `bayescomp_optimize()`, `bayescomp_pathfinder()` and `bayescomp_laplace()` gain
  `refresh` for progress output; `bayescomp_format()` keeps the
  comma of a one element tuple and prints empty blocks as stanc does;
  `bayescomp_info()` lists `print`, `reject` and operators called by name.
* `bayescomp_write_csv()` needs an explicit `dir`.
* `%/%` binds tighter than `*`, as in stanc; dividing two integers with `/`
  gives stanc's "Found int division" warning, also in `bayescomp_pedantic()`.
* Under `propto`, densities keep the terms of parameter dependent arguments
  that are plain numbers (after `floor()` or an assigned constant), as Stan does.
* Profiles record the reverse pass time and node count, and a profile cannot
  start again under a name that is still running. Sampling fits record
  rejected proposals by reason (`fit$rejections`), which `bayescomp_diagnose()`
  reports; it also checks E-BFMI per chain.
* LOO relative efficiencies no longer underflow; non finite counts are
  rejected; a nested tape no longer copies the outer tape.
* `matrix_normal_prec` is removed, as stanc 2.40 rejects it.
* `bayescomp_compatibility()` returns the machine readable feature matrix
  (implemented, partial, not implemented, not applicable) with the test
  files covering each feature.
