test_that("Stan CSV files written by CmdStan are read", {
  m <- bayescomp_model("parameters { real mu; vector<lower=0>[2] s; matrix[2,2] A; } model { mu ~ normal(0,1); s ~ lognormal(0,1); to_vector(A) ~ normal(0,1); } generated quantities { array[2] int k = {1, 2}; }")
  fit <- bayescomp_read_csv(test_path("fixtures", "cmdstan_sample.csv"), model = m)
  expect_identical(fit$method, "sample")
  expect_identical(dimnames(fit$draws)[[3L]], c("lp__", bayescomp_param_names(m)))
  expect_identical(fit$n_warmup_saved, 3L)
  expect_equal(fit$step_size, 0.396327)
  expect_equal(fit$inv_metric[[1L]], rep(1, 7))
  expect_equal(unname(fit$draws[5L, 1L, "A[2,1]"]), -0.35642127)
  expect_equal(unname(fit$sampler_diagnostics[4L, 1L, "n_leapfrog__"]), 15)
  # the unconstrained draws reproduce the constrained ones
  expect_equal(bayescomp_constrain(m, fit$unconstrained[[1L]][5L, ])$s, c(4.608055, 0.77509752), tolerance = 1e-7)
  expect_identical(dim(bayescomp_generate_quantities(m, fit, seed = 1)$draws)[1L], 5L)
})

test_that("fits written as Stan CSV read back unchanged", {
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }
                        generated quantities { complex z = to_complex(mu, 2); tuple(real, array[2] int) t = (mu, {1, 2}); }")
  fit <- bayescomp_sample(m, chains = 2, iter_warmup = 30, iter_sampling = 10, save_warmup = TRUE, seed = 1, refresh = 0)
  f <- bayescomp_write_csv(fit, tempfile(), sig_figs = 18)
  expect_length(f, 2L)
  h <- readLines(f[1L])
  expect_true("lp__,accept_stat__,stepsize__,treedepth__,n_leapfrog__,divergent__,energy__,mu,z.real,z.imag,t:1,t:2.1,t:2.2" %in% h)
  back <- bayescomp_read_csv(f, model = m)
  expect_equal(back$draws, fit$draws)
  expect_equal(back$sampler_diagnostics, fit$sampler_diagnostics)
  expect_equal(back$unconstrained, fit$unconstrained)
  expect_identical(back$n_warmup_saved, fit$n_warmup_saved)
  o <- bayescomp_optimize(m, seed = 1)
  h <- readLines(bayescomp_write_csv(o, tempfile()))
  expect_true(any(startsWith(h, "lp__,converged__,mu,")))
  v <- bayescomp_variational(m, seed = 1, refresh = 0, output_samples = 4)
  vf <- bayescomp_write_csv(v, tempfile())
  vr <- bayescomp_read_csv(vf)
  expect_identical(dim(vr$draws), dim(v$draws))
  expect_equal(from_cmdstan_names(c("A.1.2", "z.real", "t:2.1", "at.1:2.real", "mu")), c("A[1,2]", "z[real]", "t:2[1]", "at[1:2,real]", "mu"))
})

test_that("diagnostic and metric files follow CmdStan", {
  m <- bayescomp_model("parameters { simplex[3] th; real mu; } model { mu ~ normal(0, 1); }")
  f <- bayescomp_sample(m, chains = 1, iter_warmup = 30, iter_sampling = 3, seed = 1, refresh = 0, save_diagnostics = TRUE)
  d <- tempfile()
  bayescomp_write_csv(f, d, save_metric = TRUE)
  x <- readLines(file.path(d, "output_diagnostic.csv"))
  hdr <- x[!startsWith(x, "#")][1L]
  expect_match(hdr, "energy__,th.1,th.2,mu,p_th.1,p_th.2,p_mu,g_th.1,g_th.2,g_mu$")
  # the gradient of the potential of mu ~ normal(0, 1) is mu
  expect_equal(f$diagnostics[[1L]]$g[, 3L], f$unconstrained[[1L]][, 3L])
  mj <- bayescomp_read_data(file.path(d, "output_metric.json"))
  expect_identical(mj$metric_type, "diag_e")
  expect_equal(unlist(mj$inv_metric), f$inv_metric[[1L]])
  f2 <- bayescomp_sample(m, chains = 1, iter_warmup = 0, iter_sampling = 2, seed = 1, refresh = 0,
                         inv_metric = file.path(d, "output_metric.json"), adapt_engaged = FALSE)
  expect_equal(f2$inv_metric[[1L]], f$inv_metric[[1L]])
})

test_that("Pathfinder has CmdStan's outputs and options", {
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(1, 2); }")
  p <- bayescomp_pathfinder(m, seed = 1, num_paths = 2, num_draws = 40, num_psis_draws = 60, save_single_paths = TRUE)
  expect_identical(dim(p$draws), c(60L, 1L, 4L))
  expect_true(all(p$path %in% 1:2))
  expect_length(p$single_paths, 2L)
  d <- tempfile()
  bayescomp_write_csv(p, d)
  expect_setequal(list.files(d), c("output.csv", "output_path_1.csv", "output_path_1.json", "output_path_2.csv", "output_path_2.json"))
  x <- readLines(file.path(d, "output.csv"))
  expect_true("lp_approx__,lp__,path__,mu" %in% x)
  one <- bayescomp_pathfinder(m, seed = 1, num_paths = 1, num_draws = 30, calculate_lp = FALSE)
  expect_identical(dim(one$draws)[1L], 30L)
  expect_true(all(is.nan(one$draws[, 1L, "log_p__"])))
})

test_that("Laplace takes its mode from files", {
  m <- bayescomp_model("parameters { real<lower=0> s; } model { s ~ gamma(3, 2); }")
  o <- bayescomp_optimize(m, jacobian = TRUE, seed = 1)
  f <- bayescomp_write_csv(o, tempfile(), sig_figs = 18)
  a <- bayescomp_laplace(m, mode = o, draws = 20, seed = 3)
  b <- bayescomp_laplace(m, mode = f, draws = 20, seed = 3)
  expect_equal(a$draws, b$draws)
  j <- tempfile(fileext = ".json")
  bayescomp_write_data(list(s = unname(o$par[["s"]])), j)
  expect_equal(bayescomp_laplace(m, mode = j, draws = 20, seed = 3)$draws, a$draws)
  expect_true(all(is.nan(bayescomp_laplace(m, mode = o, draws = 5, seed = 1, calculate_lp = FALSE)$draws[, 1L, "log_p__"])))
})

test_that("profile blocks are timed and written as CmdStan's profile file", {
  m <- bayescomp_model("parameters { real mu; } model { profile(\"prior\") { mu ~ normal(0, 1); } }
                        generated quantities { profile(\"gq\") { real z = mu; } }")
  f <- bayescomp_sample(m, chains = 1, iter_warmup = 20, iter_sampling = 10, seed = 1, refresh = 0)
  p <- f$profile
  expect_setequal(p$name, c("prior", "gq"))
  expect_identical(p$no_autodiff_calls[p$name == "gq"], 10L)
  expect_gt(p$autodiff_calls[p$name == "prior"], 0L)
  d <- tempfile()
  bayescomp_write_csv(f, d)
  x <- readLines(file.path(d, "output_profile.csv"))
  expect_identical(x[1L], "name,thread_id,total_time,forward_time,reverse_time,chain_stack,no_chain_stack,autodiff_calls,no_autodiff_calls")
})

test_that("stansummary and diagnose reports follow CmdStan", {
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
  f <- bayescomp_sample(m, chains = 2, iter_warmup = 100, iter_sampling = 100, seed = 1, refresh = 0)
  out <- utils::capture.output(tab <- bayescomp_stansummary(f, percentiles = c(10, 90), csv_file = csv <- tempfile()))
  expect_identical(rownames(tab), c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__", "mu"))
  expect_identical(names(tab), c("Mean", "MCSE", "StdDev", "MAD", "10%", "90%", "ESS_bulk", "ESS_tail", "ESS_bulk/s", "R_hat"))
  expect_true(is.nan(tab["stepsize__", "ESS_bulk"]))
  expect_equal(tab["mu", "Mean"], mean(bayescomp_draws(f, "mu")))
  expect_identical(nrow(utils::read.csv(csv)), 8L)
  rep <- utils::capture.output(d <- bayescomp_diagnose(bayescomp_write_csv(f, tempfile())))
  expect_identical(rep[1L], "Checking sampler transitions treedepth.")
  expect_true("Rank-normalized split effective sample size satisfactory for all parameters." %in% rep)
  expect_equal(d$ebfmi_all, {e <- as.vector(f$sampler_diagnostics[, , "energy__"]); mean(diff(e)^2) / stats::var(e)}, tolerance = 1e-6)
})

test_that("generated quantities run on draws from files, matrices and posterior objects", {
  m <- bayescomp_model("parameters { real<lower=0> s; vector[2] b; } model { s ~ exponential(1); b ~ normal(0, s); }
                        generated quantities { real t = s + sum(b); }")
  f <- bayescomp_sample(m, chains = 2, iter_warmup = 30, iter_sampling = 5, seed = 1, refresh = 0)
  ref <- bayescomp_draws(bayescomp_generate_quantities(m, f, seed = 1), "t")
  mat <- bayescomp_draws(f, format = "matrix")
  g <- bayescomp_generate_quantities(m, mat, seed = 1)
  expect_equal(as.vector(bayescomp_draws(g, "t")), as.vector(mat[, "s"] + mat[, "b[1]"] + mat[, "b[2]"]))
  csv <- bayescomp_write_csv(f, tempfile(), sig_figs = 18)
  expect_equal(bayescomp_draws(bayescomp_generate_quantities(m, csv, seed = 1), "t"), ref)
  dotted <- mat[, c("s", "b[1]", "b[2]")]
  colnames(dotted) <- c("s", "b.1", "b.2")
  expect_equal(as.vector(bayescomp_draws(bayescomp_generate_quantities(m, dotted), "t")), as.vector(bayescomp_draws(g, "t")))
  skip_if_not_installed("posterior")
  pd <- posterior::as_draws_array(bayescomp_draws(f))
  expect_equal(bayescomp_draws(bayescomp_generate_quantities(m, pd, seed = 1), "t"), ref)
  expect_error(bayescomp_generate_quantities(m, mat[, "s", drop = FALSE]), "no values for b")
})

test_that("optimization can save its iterations", {
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(3, 1); }")
  for (alg in c("lbfgs", "bfgs", "newton")) {
    o <- bayescomp_optimize(m, algorithm = alg, seed = 1, save_iterations = TRUE)
    k <- dim(o$draws)[1L]
    expect_gt(k, 1L)
    expect_equal(unname(o$par[["mu"]]), 3, tolerance = 1e-5)
    expect_equal(unname(o$draws[k, 1L, "mu"]), unname(o$par[["mu"]]))
  }
  x <- readLines(bayescomp_write_csv(o, tempfile()))
  rows <- x[!startsWith(x, "#")][-1L]
  expect_identical(length(rows), k)
  expect_true(all(sub(",.*", "", sub("^[^,]*,", "", rows[-k])) == "0"))
  expect_equal(bayescomp_laplace(m, mode = o, draws = 2, seed = 1)$draws[, 1L, "lp__"], c(0, 0))
})

test_that("the run configuration is saved as JSON", {
  m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(0, 1); }")
  f <- bayescomp_sample(m, chains = 1, iter_warmup = 10, iter_sampling = 5, seed = 3, refresh = 0)
  d <- tempfile()
  bayescomp_write_csv(f, d, save_config = TRUE)
  j <- jsonlite::read_json(file.path(d, "output_config.json"))
  expect_identical(j$method$value, "sample")
  expect_identical(j$method$sample$num_chains, 1L)
  expect_identical(j$stan_minor_version, "40")
  expect_identical(j$random$seed, 3L)
})
