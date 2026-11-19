# bayescompR (development version)

## New features

* Models written in the BUGS language run without JAGS, OpenBUGS, WinBUGS,
  MultiBUGS or NIMBLE installed. `jags()` takes the arguments of
  `R2jags::jags()` and returns its result; `bugs()` does the same for
  `R2OpenBUGS::bugs()`; `jags.model()`, `update()`, `coda.samples()` and
  `jags.samples()` follow rjags; and `nimbleCode()`, `nimbleModel()`,
  `configureMCMC()`, `buildMCMC()`, `compileNimble()`, `runMCMC()` and
  `nimbleMCMC()` follow NIMBLE. The models are compiled and sampled in R
  with Gibbs updates (conjugate gamma, beta, Dirichlet, Wishart and normal
  updates for linear models, enumeration of discrete nodes) and slice
  sampling. Models can be files, text, or R functions; data and initial
  values can be R lists, R dump files, or OpenBUGS list and rectangular
  files. `vignette("bugs")` lists what is supported.

* `bayescomp_ensemble()` gains `method = "side"` (a Gaussian step along the
  difference of two walkers, tuned toward acceptance 0.5) and
  `method = "slice"` (ensemble slice sampling with a scale tuned during
  warmup). The stretch move warns above 50 parameters, where it mixes poorly.

* `bayescomp_parallel_tempering()` gains `adapt_reference`, which refits the
  Gaussian reference to the posterior replica after each full warmup round,
  and reports `suggested_num_temps` from the estimated barrier.

* Finite discrete BUGS nodes use metropolized Gibbs, which proposes a value
  other than the current one and so moves more often than plain Gibbs while
  keeping the same conditional. Slice widths adapt to four times the mean
  jump, near the width that minimizes stepping out cost.

## Bug fixes

* `bayescomp_ess_quantile()` and `bayescomp_ess_tail()` use the indicator
  of the upper tail (`x >= q`) for probabilities above one half, as
  posterior 1.8.0 does.

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
* Loop vectorization no longer rewrites a loop over a user defined density
  (whose scalar and vector overloads can differ), evaluates nothing for a
  loop with no iterations, and leaves loops whose bounds read `target()`,
  have effects or change in the body as loops.
* A `for` loop tests its upper bound before every iteration, as Stan does,
  when the bound reads `target()` or a variable the body assigns; it was
  evaluated once.
* The prepared `normal_id_glm()` factorization centers the data and
  standardizes the residuals, so data far from the origin and very small
  or large scales give the same value and gradient as the unfused
  density; it falls back to the observation level sum when its rounding
  error is not negligible against the scale.
* Truncation of a density with a scalar outcome and vector parameters
  subtracts one normalizing term per element.
* Metropolis, slice, ensemble, Gibbs and reverse diffusion samplers use the
  density with constants dropped, as HMC does. The two targets differ by more than a
  constant when a program reads `target()` or drops a data only statement
  in a parameter dependent branch. Evidence estimators and model
  probabilities keep the normalized density.
* ARD covariance functions of fixed inputs (`gp_exp_quad_cov`, Matern and
  exponential with one length scale per dimension) take differences
  before scaling, so the covariance does not depend on the origin.
* `cholesky_factor_corr` and `corr_matrix` parameters compute
  `1 - tanh(u)^2` from `u`, so the factor stays valid (with a positive
  diagonal) far in the tails.
* `neg_binomial_2` and `neg_binomial_2_log` are one node with stable
  partials: finite gradients for very small and very large precision and
  log means up to the largest finite values.
* SG-MCMC and expectation propagation build the data subsets with the seed
  of the full model, so random numbers drawn in transformed data are the
  same in every subset (a model without a seed whose transformed data
  draws random numbers is an error).
* Pointwise `poisson_log` and `neg_binomial_2_log` terms broadcast a scalar
  outcome over vector parameters; pointwise replay binds each chain's
  `CHAIN_ID`; observations keep the name of an indexed variable; draws that
  record different observations are an error; truncated terms of
  observations outside their bounds are `-Inf`.
* `bayescomp_sample(resume = )` reports an exact continuation only when
  every chain finished and kept its state and random number stream, and
  transformed data random numbers are the same (with a seed);
  `fixed_param` chains keep their thinning phase. Generated quantities
  draw from their own substream of each chain's stream, so a continuation
  also continues them instead of repeating or overlapping draws.
  `step_size = NULL` and `inv_metric = NULL` keep the adapted values, and
  a different `algorithm` warns that it is ignored.
* Expectation propagation warns about sites whose update never succeeded
  and does not report convergence while one is missing.
* Adaptation windows of one draw keep the metric; `bayescomp_hessian_vector()`
  scales its step by the coordinates the direction moves; the ODE
  `max_num_steps` applies between output times, as in Stan;
  `multi_normal_cholesky` with arrays of vectors checks its Cholesky
  factor (as Stan does) and `multi_student_t` the sizes of its arguments; SG-MCMC rejects a nonpositive step size;
  Laplace fits record the seed they drew.
* Complex compound assignment (`z += 1`, `w *= z`) uses complex
  arithmetic; `z ^ b` keeps the derivatives of the other elements when one
  base is zero and gives `0 ^ 0 = 1`; `norm()`, `arg()`, `proj()` and the
  complex elementary functions accept real arguments.
* `categorical_logit_glm` with a `row_vector` of predictors uses every
  outcome and checks the size of the intercepts, and its `_lupmf` is 0 for
  data only arguments; `normal_id_glm` accepts a `vector` of predictors as
  one column; a target increment that is `NA` is an error.
* `gamma`, `rayleigh` and `loglogistic` densities require a positive
  variate, as in Stan; `mdivide_left_spd` rejects a matrix that is not
  positive definite.
* `determinant()` has the cofactor gradient also at a singular matrix, and
  `log_determinant()` rejects the point instead of failing.
* `student_t` and `cauchy` densities are finite for residuals whose square
  overflows.
* `max()` and `min()` of a container with NaN return NaN.
* Arguments of user defined functions called in the model block keep the
  terms of `_lupdf` densities they reach, also after `floor()` of a
  parameter.
* A zero adjoint no longer turns an infinite elementwise derivative into
  NaN, so `sqrt()` of a matrix with zeros has a finite gradient.
* `return .5;` parses, an integer division of a multi indexed value
  (`x[idx][1] / 2`) compiles, and sizes of top level variables that depend
  on parameters through an integer are a semantic error.
* SG-MCMC, expectation propagation and nested sampling accept a model that
  declares `int<lower=1> N`; parallel SMC and EP workers draw from seeds
  per stage instead of overlapping streams, and a worker error is
  reported.
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

* `multi_normal()` with a covariance that depends on parameters (a
  marginal Gaussian process) computes its partials in the covariance
  directly, as Stan does, instead of differentiating through the Cholesky
  factor: about 4 times faster gradients with 200 observations. The
  adjoint of `cholesky_decompose()` uses triangular solves instead of an
  explicit inverse (about 1.6 times faster for latent Gaussian processes).

* A truncated density with one bound and vector parameters computes its
  normalizer with one vectorized `_lccdf` or `_lcdf` call instead of one
  call per element (about 200 times faster over 1000 observations), and
  loop vectorization now also rewrites truncated statements whose bounds
  do not depend on the loop.

* Loop vectorization also rewrites loops whose body declares scalar
  locals (`real mu = ...;`, replaced by their values), uses
  `dot_product(X[n], b)`, adds any scalar to `target` (including `_lcdf`
  and `_lccdf` terms), or sits in a `profile` block: 10 to 28 times faster
  gradients for such loops over 200 observations. A loop such as
  `for (n in 1:N) y[n] ~ ordered_logistic(eta[n], c)` with a shared vector
  argument (cutpoints, a simplex, logits) that does not depend on the loop
  is rewritten as well (about 50 times faster over 300 observations).

* Proportional log density values (Metropolis, slice, ensemble and Gibbs
  sampling, ADVI, Pathfinder, optimization) no longer record a tape. The
  program runs on plain numbers and only the parameter dependent arguments
  of densities are marked, so each density drops the terms Stan drops
  (about 1.4 times faster on a hierarchical regression). Programs with
  user defined functions or solvers keep the class of parameter dependent
  values without recording them.

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
* The negative binomial precision partial is one vectorized digamma
  difference (with a series where the precision is large) instead of a sum
  per observation.
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
