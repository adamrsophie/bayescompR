# Roadmap

bayescompR aims at a capability aware set of inference engines that fail in
different ways, not at the largest number of algorithms. Each method states
what it needs from the compiled model (log density only, gradient, Gaussian
blocks, conjugate structure, sequential structure) and is offered only where
those needs are met.

## Shipped in 0.1.0

Gradient based MCMC

* NUTS and static HMC (diagonal or dense metric, Stan warmup), MALA.
* Simplified manifold MALA and Riemannian manifold HMC with the SoftAbs
  metric (`bayescomp_riemannian()`), for small models.
* Microcanonical Langevin Monte Carlo, unadjusted (`bayescomp_mclmc()`).
* Zig-Zag and bouncy particle samplers with numerical event times
  (`bayescomp_pdmp()`).
* Stochastic gradient Langevin dynamics and stochastic gradient HMC on
  minibatches (`bayescomp_sgmcmc()`).

Gradient free MCMC

* Block and componentwise random walk Metropolis with adaptive covariance,
  coordinate slice sampling, generalized elliptical slice sampling.
* Differential evolution MCMC and the affine invariant stretch move
  (`bayescomp_ensemble()`).
* Non reversible and reversible parallel tempering with schedule adaptation
  and a stepping stone evidence (`bayescomp_parallel_tempering()`).
* Blocked Gibbs with per block kernels and numerical detection of Gaussian,
  Gamma and Beta full conditionals (`bayescomp_gibbs()`).
* Reversible jump MCMC across a list of models (`bayescomp_rjmcmc()`).

Sequential and evidence methods

* Adaptive tempered SMC with Metropolis, MALA or HMC rejuvenation and four
  resampling schemes; nested sampling.
* State space models: bootstrap particle filter, particle marginal
  Metropolis-Hastings, particle Gibbs with ancestor sampling.

Optimization and approximations

* L-BFGS, BFGS and Newton optimization; Laplace approximation.
* Mean field and full rank ADVI; multipath Pathfinder with Pareto k; planar
  normalizing flow VI.
* Expectation propagation over data partitions; INLA for a user chosen
  latent block.
* Neural posterior estimation for simulators; reverse diffusion Monte Carlo
  (experimental).

Diagnostics

* Rank normalized and folded split R-hat, bulk, tail and quantile ESS, MCSE,
  divergences, tree depth, E-BFMI; trace, density, rank, local ESS and
  quantile ESS plots; simulation based calibration in the long CI tests.

## Next (0.2.0)

* GPU, TPU and NPU execution through a backend package that lowers the IR to
  jaxR (StableHLO, XLA through PJRT: CUDA, ROCm, TPU, Metal and PJRT NPU
  plugins). The IR already carries the shapes and types needed. Apple's
  Neural Engine has no public general compute API.
* Conjugacy detection on the IR instead of numerically, which would also
  cover Dirichlet, Wishart and Poisson Gamma blocks.
* Exact event times for Zig-Zag and BPS on models whose gradient bounds can
  be derived from the IR.
* Richer flows (coupling layers) for normalizing flow VI and neural
  posterior estimation.

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
* ter Braak (2006). Differential evolution Markov chain.
* Goodman and Weare (2010). Ensemble samplers with affine invariance.
* Welling and Teh (2011). Stochastic gradient Langevin dynamics.
* Chen, Fox and Guestrin (2014). Stochastic gradient HMC.
* Lindsten, Jordan and Schon (2014). Particle Gibbs with ancestor sampling.
* Rezende and Mohamed (2015). Variational inference with normalizing flows.
* Green (1995). Reversible jump MCMC.
* Skilling (2006). Nested sampling.
* Rue, Martino and Chopin (2009). INLA.
* Betancourt (2013). The SoftAbs metric.
* Bouchard-Cote, Vollmer and Doucet (2018). The bouncy particle sampler.
* Bierkens, Fearnhead and Roberts (2019). The Zig-Zag process.
* Vehtari et al. (2020). Expectation propagation as a way of life.
* Syed, Bouchard-Cote, Deligiannidis and Doucet (2022). Non reversible
  parallel tempering.
* Robnik, De Luca, Silverstein and Seljak (2023). Microcanonical HMC.
* Huang et al. (2024). Reverse diffusion Monte Carlo.
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
