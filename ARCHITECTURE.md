# Architecture

```
Stan source
  -> lexer.R        tokens with line and column
  -> parser.R       untyped AST
  -> semantic.R     typed AST (the IR): types, data_only flags, resolved
                    signatures, block placement rules
  -> codegen.R      R language objects, one closure per model function
  -> model.R        bayescomp_model: data reading and validation, transformed
                    data, log density, gradient, constrain and unconstrain
```

The runtime (runtime.R) handles declarations, indexing, assignment, target
accumulation and printing. Values are plain R numbers or `bcad` nodes on the
reverse mode tape in ad.R. Library functions (library_*.R) accept both and
create fused tape nodes with analytic partial derivatives. transforms.R holds
Stan's constraining transforms and their log Jacobians.

## Inference

Every algorithm uses the same model interface: the unconstrained log density
with or without its gradient, constraining for output, and generated
quantities.

| Method | Needs | Exact | File |
|--------|-------|-------|------|
| NUTS, static HMC | gradient | yes | nuts.R |
| MALA | gradient | yes | mcmc.R |
| Random walk Metropolis (block, componentwise) | log density | yes | mcmc.R |
| Slice sampling | log density | yes | mcmc.R |
| Tempered SMC | log density | consistent as particles grow | smc.R |
| Optimization, Laplace | gradient, Hessian by differences | approximate | optimize.R |
| ADVI, Pathfinder | gradient | approximate | vi.R |
| Elliptical slice sampling | log density | yes | mcmc.R |
| Parallel tempering | log density | yes | tempering.R |
| DE MCMC, stretch move | log density | yes | ensemble.R |
| Blocked Gibbs, conjugate detection | log density (gradient for MALA blocks) | yes | gibbs.R |
| smMALA, RMHMC | gradient, Hessian by differences | yes (RMHMC up to finite differences) | geometric.R |
| MCLMC | gradient | approximate (unadjusted) | geometric.R |
| Zig-Zag, BPS | gradient | approximate (numerical event times) | pdmp.R |
| Planar flow VI | gradient | approximate | flow.R |
| Particle filter, PMMH, particle Gibbs | user state space functions | yes | ssm.R |
| SGLD, SGHMC | gradient, data split | approximate | bigdata.R |
| Expectation propagation | gradient, data split | approximate | bigdata.R |
| Nested sampling | log density, data split | consistent | evidence.R |
| Reversible jump | log density per model | yes | evidence.R |
| INLA | gradient, latent block | approximate | inla.R |
| Neural posterior estimation | prior model, simulator | approximate | sbi.R |
| Reverse diffusion Monte Carlo | log density | approximate | diffusion.R |

All MCMC kernels share the chain driver `run_hmc_chain()`: Stan's warmup
windows, dual averaging toward `adapt_delta`, and Welford estimation of the
inverse metric, which is the proposal covariance for Metropolis and MALA.
Adaptation stops at the end of warmup. `run_chains()` gives each chain its
own L'Ecuyer-CMRG stream, so results do not depend on `cores`.

Methods with their own drivers (parallel tempering, ensembles, Gibbs,
Riemannian samplers, MCLMC, PDMPs, particle MCMC, SG-MCMC) return the same
`bayescomp_fit` structure, so summaries, R-hat, ESS and plots apply. The
data split methods rebind subsets of the data to the compiled model
(`model_with_data()`), with the prior given by the empty data set.

## Backends and the IR

`bayescomp_ir()` serializes the typed AST to JSON with an `ir_version`
(semantic versioning). Backends register through
`bayescomp_register_backend()` and declare capabilities (log density,
gradient, transforms, generated quantities, RNG). The pure R backend is
always available; bayescompRCranelift and bayescompRLLVM compile the IR to
machine code. A backend never falls back to R silently: fallback requires
`allow_fallback = TRUE`.

## Performance notes

* Never call `length()` on the tape lists in `.ad`: reading a list marks it
  shared, and the next element assignment copies the whole tape.
* Tape lists are detached, modified and reattached so assignment is in
  place.
* Densities are single tape nodes (`fused()` in library_dist.R).
* Value only evaluation (Metropolis, slice, SMC) disables the tape; propto
  is off there because every value is a constant, and the chain driver shifts
  lp__ by a constant to match the propto value.
