# Ecosystem compatibility checks

These scripts check that the Stan programs written by brms, rstanarm,
multinma and edstan run in bayescompR with the same log density gradient
as rstan. They are not part of the package build.

1. Capture programs and data into a directory, one `.rds` per model:
   `Rscript validation/ecosystem/brms_programs.R out/brms` (likewise
   `rstanarm_programs.R` and `multinma_programs.R`, which intercept
   `rstan::sampling` before any sampling happens).
2. Compare: `Rscript validation/ecosystem/compare_rstan.R out/brms out/rstanarm ...`.
   Each model prints `MATCH` when the gradients agree to 1e-6 (relative)
   at three random unconstrained points and the log density differs from
   rstan's by a constant.

Programs that declare functions without a body (rstanarm's `CODOLS`)
receive an R implementation through the `functions` argument of
`bayescomp_model()`.
