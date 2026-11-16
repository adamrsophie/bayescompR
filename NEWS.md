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
* `bayescomp_log_lik()` gives fused GLMs (`normal_id_glm`,
  `bernoulli_logit_glm`, `poisson_log_glm`, `neg_binomial_2_log_glm`,
  `binomial_logit_glm`, `ordered_logistic_glm`, `categorical_logit_glm`) one
  column per row of `x`, as the equivalent unfused statement; it gave one
  joint column. `bayescomp_posterior_predict()` draws from them.
* `bayescomp_posterior_predict()` works for multivariate statements such as
  `multi_normal`, with one column per element of the left side; it failed
  with "number of items to replace is not a multiple of replacement length".
* Pointwise terms of statements with a shared vector argument, such as
  `categorical_logit` log odds or `ordered_logistic` cutpoints, no longer
  split that argument across observations; they failed with a domain error.
  Cutpoints given as one vector per observation now give one term per
  observation instead of one joint term.
* `ordered_logistic_rng()` and `ordered_probit_rng()` with a vector of
  linear predictors draw one category per element; they mixed the
  cutpoints of all elements into one draw.

* The squared distances of GP covariance functions with fixed inputs are
  computed from coordinate differences, so points far from the origin keep
  their separation, and values with and without gradients agree.
* Two sided truncation of a distribution with an upper tail function
  (`T[40, 41]` on a standard normal) uses upper tail probabilities for
  intervals above the median; the density was infinite.
* `wishart_cholesky_lpdf()` and `inv_wishart_cholesky_lpdf()` reject
  matrices that are not Cholesky factors (lower triangular with a positive
  diagonal), as Stan does, instead of evaluating them through the
  covariance matrices.

## New features

* `bayescomp_sgmcmc()` gains `lp_every`, so the exact full data `lp__` can
  be computed at only some saved draws or not at all.
* `bayescomp_log_lik()` and `bayescomp_posterior_predict()` support
  truncated statements (renormalized per observation; predictive draws by
  inversion for the normal, which also works far in the tails, and by
  rejection otherwise) and user defined densities (one joint term per
  statement), and `bayescomp_log_lik(group = )` sums observations by group
  for leave one group out cross validation. Pointwise terms use the
  argument values the statement used, also when an argument reads
  `target()`.
* `bayescomp_hessian()` and `bayescomp_hessian_vector()` give the Hessian
  of the log density and its product with a vector, by central differences
  of the gradient with steps scaled to the parameter values.
* `bayescomp_capabilities()` lists the algorithm behind each ODE, DAE,
  algebraic solver, integrator, `map_rect` and `reduce_sum` interface, the
  controls it honors and ignores, and how it is differentiated.
* `bayescomp_sample(resume = )` continues a finished fit exactly: from the
  state after the last transition (not the last kept draw) with the chain's
  random number stream, transition controls and thinning phase, also after
  a warmup only run, so the parameter draws equal those of one longer run.
  `config$resume` of the result records whether a continuation was exact.

## Performance

* The factorizations of Gaussian GLMs with fixed data are kept within a
  memory budget instead of 64 per model, so models with many such
  statements do not factor them again on every evaluation.

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
  predictive draws compute and keep only the requested array.
* Pointwise log likelihoods of the common univariate families are computed
  for all observations at once, and predictive draws call the vectorized
  `_rng` function once per statement.
* Dirichlet densities of arrays of simplexes work on whole matrices instead
  of splitting them into rows.
* `integrate_1d()` and its variants evaluate the gradient of the integrand
  once per integration point instead of once per point and parameter.
* The `normal_id_glm` factorizations belong to the model instance, so
  SG-MCMC and expectation propagation with many data subsets no longer
  refactor at every evaluation.
* `gp_exp_quad_cov()` with one length scale per dimension and fixed inputs
  is one node with a dedicated reverse pass.
* `wishart_cholesky` and `inv_wishart_cholesky` densities use the factors
  directly (one triangular solve) instead of forming and refactoring the
  covariance matrices.
* The stiff ODE solver solves the forward sensitivity stages with one
  factorization of the state block instead of the full augmented system.

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
