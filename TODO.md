# Roadmap

bayescompR aims at a capability aware set of inference engines that fail in
different ways, not at the largest number of algorithms. Each method states
what it needs from the compiled model (log density only, gradient, Gaussian
blocks, conjugate structure, sequential structure) and is offered only where
those needs are met.

## Shipped in 0.1.0

* NUTS and static HMC (diagonal or dense metric, Stan warmup).
* MALA, block and componentwise random walk Metropolis with adaptive
  covariance during warmup, coordinate slice sampling (stepping out and
  shrinkage).
* L-BFGS, BFGS and Newton optimization; Laplace approximation.
* Mean field and full rank ADVI; single and multipath Pathfinder with Pareto
  k diagnostics; Pathfinder or ADVI fits as MCMC initializers.
* Adaptive tempered SMC with a Gaussian reference and log evidence estimate.
* Rank normalized and folded split R-hat, bulk and tail ESS, MCSE,
  divergences, tree depth, E-BFMI; simulation based calibration in the long
  CI tests.

## Next (0.2.0)

* GPU, TPU and NPU execution through a backend package that lowers the IR to
  jaxR (StableHLO, XLA through PJRT: CUDA, ROCm, TPU, Metal and PJRT NPU
  plugins). The IR already carries the shapes and types needed. Apple's
  Neural Engine has no public general compute API.
* Parallel tempering (replica exchange) with temperature ladder adaptation,
  then non reversible swap schedules.
* Elliptical slice sampling for Gaussian prior blocks detected in the IR.
* Differential evolution MCMC and affine invariant ensemble samplers.
* User composable blocked samplers (different kernels for parameter
  blocks).
* HMC or MALA rejuvenation and stratified or residual resampling in SMC.
* Folded and localized ESS plots; rank plots.

## Later

* Compiler pass that detects a whitelist of conjugate structures, then
  automatic Gibbs and block Gibbs.
* State space abstractions, particle filters and particle MCMC.
* Riemannian manifold HMC, microcanonical Langevin Monte Carlo, Zig-Zag and
  bouncy particle samplers, normalizing flow VI: realistic once the native
  backends are mature.

## Not planned for the core package

INLA, generic expectation propagation, reversible jump MCMC, neural
simulation based inference, diffusion posterior sampling, stochastic
gradient MCMC and nested sampling. They target narrower or structurally
different problems; tempered SMC already covers evidence estimation.

## Language features not yet supported

Complex numbers, tuples, ODE and DAE solvers, algebraic solvers, 1D
integrators, HMM functions, the embedded Laplace functions, `map_rect` and
`reduce_sum`, sparse CSR functions, FFT, SVD and general eigendecomposition,
and the Wiener, Gaussian DLM, multi GP and beta negative binomial
distributions. Programs that use them fail at checking with
`bayescomp_unsupported_error`.

## Developer bibliography

* Haario, Saksman and Tamminen (2001). An adaptive Metropolis algorithm.
* Neal (2003). Slice sampling.
* Del Moral, Doucet and Jasra (2006). Sequential Monte Carlo samplers.
* Murray, Adams and MacKay (2010). Elliptical slice sampling.
* Andrieu, Doucet and Holenstein (2010). Particle Markov chain Monte Carlo
  methods.
* Girolami and Calderhead (2011). Riemann manifold Langevin and Hamiltonian
  Monte Carlo methods.
* Hoffman and Gelman (2014). The No-U-Turn sampler.
* Betancourt (2017). A conceptual introduction to Hamiltonian Monte Carlo.
* Kucukelbir, Tran, Ranganath, Gelman and Blei (2017). Automatic
  differentiation variational inference.
* Talts, Betancourt, Simpson, Vehtari and Gelman (2018). Validating Bayesian
  inference algorithms with simulation based calibration.
* Vehtari, Gelman, Simpson, Carpenter and Buerkner (2021). Rank
  normalization, folding, and localization: an improved R-hat.
* Zhang, Carpenter, Gelman and Vehtari (2022). Pathfinder: parallel quasi
  Newton variational inference.
* Vehtari, Simpson, Gelman, Yao and Gabry (2024). Pareto smoothed importance
  sampling.
