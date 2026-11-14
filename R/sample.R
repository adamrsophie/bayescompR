# Sampling service: initialization, per chain random number streams,
# parallel chains and assembly of the fit object.

#' Draw posterior samples with Markov chain Monte Carlo
#'
#' Runs the No-U-Turn sampler (default), static HMC, the Metropolis adjusted
#' Langevin algorithm, random walk Metropolis or slice sampling. All of them
#' use Stan's warmup schedule: dual averaging of the step size toward
#' `adapt_delta` and windowed estimation of a diagonal or dense inverse
#' metric, which for Metropolis and MALA is the proposal covariance and for
#' slice sampling sets the bracket directions and widths. Adaptation stops at
#' the end of warmup, so the sampling phase uses a fixed Markov kernel.
#'
#' NUTS is the recommended default for differentiable models. MALA is a
#' cheaper gradient based kernel that can mix poorly for heavy tails.
#' Metropolis and slice sampling need no gradients: they are useful as
#' reference implementations, for checking results, and for models whose
#' log density is not smooth.
#'
#' @param model A `bayescomp_model`.
#' @param chains Number of chains.
#' @param iter_warmup,iter_sampling Warmup and sampling iterations per chain.
#' @param seed Integer seed. Each chain uses its own L'Ecuyer-CMRG stream, so
#'   results do not depend on `cores`. With a seed, the global random number
#'   state is left unchanged. Without one, a seed is drawn from the global
#'   stream (so `set.seed()` makes the call reproducible and repeated calls
#'   differ) and recorded in the fit's `seed` element.
#' @param init Initial values: a radius `r` (uniform on `(-r, r)` in the
#'   unconstrained space, default 2), `0`, a named list of constrained values,
#'   a list of such lists (one per chain), a function returning a list, or a
#'   `bayescomp_fit` (for example from [bayescomp_pathfinder()]) whose draws
#'   are used as starting points, or the names of JSON or R dump files (one
#'   per chain, recycled). Parameters missing from a list or file are
#'   initialized at random.
#' @param algorithm `"nuts"`, `"hmc"` (static integration time), `"mala"`,
#'   `"metropolis"`, `"slice"`, `"elliptical_slice"` (generalized
#'   elliptical slice sampling around a Gaussian whose mean and covariance
#'   are estimated during warmup) or `"fixed_param"` (parameters stay at
#'   their initial values; only generated quantities change, as Stan's
#'   fixed_param sampler).
#' @param metric `"diag_e"`, `"dense_e"` or `"unit_e"` (identity metric,
#'   never adapted).
#' @param inv_metric Optional initial inverse metric (vector or matrix), or
#'   the name of a JSON or R dump file with element `inv_metric` (CmdStan's
#'   `metric_file`, as written by `bayescomp_write_csv(save_metric = TRUE)`).
#' @param save_diagnostics Keep the unconstrained values, momenta and
#'   potential gradients of each saved draw (CmdStan's `diagnostic_file`) in
#'   element `diagnostics` of the fit; [bayescomp_write_csv()] writes them.
#' @param step_size Initial step size: leapfrog step for HMC, proposal scale
#'   for Metropolis and MALA, bracket width in metric units for slice
#'   sampling. `NULL` selects a default for the algorithm.
#' @param adapt_delta Target acceptance statistic. `NULL` selects 0.8 for
#'   NUTS and HMC, 0.574 for MALA, and 0.234 (block) or 0.44 (componentwise)
#'   for Metropolis. Unused for slice sampling.
#' @param block For Metropolis: propose all coordinates jointly (`TRUE`) or
#'   update one coordinate at a time.
#' @param max_treedepth Maximum tree depth for NUTS.
#' @param int_time Integration time for static HMC.
#' @param step_size_jitter Uniform random jitter of the step size in each
#'   transition, as a fraction in \[0, 1\] (Stan's `stepsize_jitter`).
#' @param adapt_engaged Adapt during warmup.
#' @param adapt_metric Adapt the metric during warmup (when
#'   `adapt_engaged`); `FALSE` adapts the step size only.
#' @param adapt_gamma,adapt_kappa,adapt_t0 Dual averaging constants.
#' @param adapt_init_buffer,adapt_term_buffer,adapt_window Warmup windows.
#' @param thin Keep every `thin`th draw.
#' @param save_warmup Keep warmup draws.
#' @param refresh Report progress every `refresh` iterations (0 for quiet).
#' @param cores Number of chains to run in parallel. Defaults to
#'   `getOption("mc.cores", 1L)`, so `options(mc.cores = 4)` makes every
#'   sampler use four cores. Chains run in forked processes, or on Windows
#'   (and with `options(bayescompR.parallel = "psock")`) in a socket cluster
#'   whose workers load the installed package and print their progress
#'   when they finish.
#' @param show_print Show output of Stan `print()` statements.
#' @return A `bayescomp_fit` object.
#' @examples
#' \donttest{
#' m <- bayescomp_model("parameters { real mu; } model { mu ~ normal(1, 2); }")
#' fit <- bayescomp_sample(m, chains = 2, iter_warmup = 200, iter_sampling = 200,
#'                         seed = 1, refresh = 0)
#' fit
#' }
#' @export
bayescomp_sample <- function(model, chains = 4L, iter_warmup = 1000L, iter_sampling = 1000L,
                             seed = NULL, init = 2, algorithm = c("nuts", "hmc", "mala", "metropolis", "slice", "elliptical_slice", "fixed_param"),
                             metric = c("diag_e", "dense_e", "unit_e"), inv_metric = NULL,
                             step_size = NULL, adapt_delta = NULL, block = TRUE, max_treedepth = 10L,
                             int_time = 2 * pi, step_size_jitter = 0, adapt_engaged = TRUE, adapt_metric = TRUE, adapt_gamma = 0.05,
                             adapt_kappa = 0.75, adapt_t0 = 10, adapt_init_buffer = 75L,
                             adapt_term_buffer = 50L, adapt_window = 25L, thin = 1L,
                             save_warmup = FALSE, refresh = 200L, cores = getOption("mc.cores", 1L), show_print = TRUE,
                             save_diagnostics = FALSE) {
  check_model(model)
  algorithm <- match.arg(algorithm)
  metric <- match.arg(metric)
  chains <- as_count(chains, "chains", min = 1L)
  iter_warmup <- as_count(iter_warmup, "iter_warmup")
  iter_sampling <- as_count(iter_sampling, "iter_sampling")
  thin <- as_count(thin, "thin", min = 1L)
  n <- model$dim
  if (!(is.numeric(step_size_jitter) && length(step_size_jitter) == 1L && step_size_jitter >= 0 && step_size_jitter <= 1))
    stop("'step_size_jitter' must be in [0, 1]", call. = FALSE)
  if (metric == "unit_e") {
    if (!is.null(inv_metric)) stop("'inv_metric' cannot be given with metric = \"unit_e\"", call. = FALSE)
    metric <- "diag_e"; inv_metric <- rep(1, n); adapt_metric <- FALSE
  }
  if (algorithm == "fixed_param") adapt_engaged <- FALSE
  if (is.null(adapt_delta)) adapt_delta <- default_adapt_delta(algorithm, block)
  if (!(algorithm %in% c("slice", "elliptical_slice", "fixed_param")) && !(adapt_delta > 0 && adapt_delta < 1)) stop("'adapt_delta' must be in (0, 1)", call. = FALSE)
  if (is.null(step_size)) {
    step_size <- switch(algorithm, nuts = 1, hmc = 1, slice = 2, elliptical_slice = 1, fixed_param = 1, mala = 1 / max(1, n)^(1 / 6),
                        metropolis = if (block) 2.38 / sqrt(max(1, n)) else 2.4)
  }
  if (!(is.numeric(step_size) && length(step_size) == 1L && step_size > 0)) stop("'step_size' must be positive", call. = FALSE)
  if (is.character(inv_metric)) {
    # CmdStan's metric_file: a JSON or R dump file with element inv_metric
    im <- bayescomp_read_data(inv_metric)$inv_metric
    if (is.null(im)) stop("the metric file has no 'inv_metric' element", call. = FALSE)
    if (is.list(im)) im <- simplify_data_list(im)
    inv_metric <- if (metric == "dense_e") as.matrix(im) else as.numeric(im)
  }
  if (!is.null(inv_metric)) {
    ok <- if (metric == "diag_e") is.numeric(inv_metric) && length(inv_metric) == n && all(is.finite(inv_metric)) && all(inv_metric > 0)
          else is.matrix(inv_metric) && all(dim(inv_metric) == n) && all(is.finite(inv_metric)) &&
            isSymmetric(unname(inv_metric)) && is_pos_def(inv_metric)
    if (!ok) stop("'inv_metric' does not match the metric type and model dimension", call. = FALSE)
  }
  cfg <- list(algorithm = algorithm, block = block, metric = metric, inv_metric = inv_metric, step_size = step_size,
              adapt_delta = adapt_delta, max_treedepth = as_count(max_treedepth, "max_treedepth", 1L),
              int_time = int_time, step_size_jitter = step_size_jitter, adapt_engaged = adapt_engaged, adapt_metric = isTRUE(adapt_metric),
              adapt_gamma = adapt_gamma, adapt_kappa = adapt_kappa, adapt_t0 = adapt_t0,
              adapt_init_buffer = as.integer(adapt_init_buffer), adapt_term_buffer = as.integer(adapt_term_buffer),
              adapt_window = as.integer(adapt_window), iter_warmup = iter_warmup, iter_sampling = iter_sampling,
              thin = thin, save_warmup = save_warmup, refresh = as.integer(refresh), print = show_print,
              diagnostic = isTRUE(save_diagnostics))
  t0 <- proc.time()[["elapsed"]]
  res <- run_chains(model, chains, seed, cores, function(id) {
    q0 <- find_init(model, init, id, gradient = !(algorithm %in% c("metropolis", "slice", "elliptical_slice", "fixed_param")))
    out <- if (n == 0L || algorithm == "fixed_param") fixed_param_chain(cfg, q0) else run_hmc_chain(model, q0, cfg, id)
    out$gq <- chain_write_draws(model, out$q, show_print)
    out
  })
  make_fit(model, res, cfg, method = "sample", seed = seed, elapsed = proc.time()[["elapsed"]] - t0)
}

as_count <- function(x, name, min = 0L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < min || x != round(x))
    stop(sprintf("'%s' must be an integer >= %d", name, min), call. = FALSE)
  as.integer(x)
}

check_model <- function(model) {
  if (!inherits(model, "bayescomp_model")) stop("'model' must be a bayescomp_model", call. = FALSE)
}

# fixed_param (and models without parameters): the parameters stay at q0
# and only generated quantities are produced; lp__ is 0 as in Stan
fixed_param_chain <- function(cfg, q0 = numeric(0)) {
  ns <- length(seq_len(cfg$iter_sampling)[(seq_len(cfg$iter_sampling) - 1L) %% cfg$thin == 0L])
  dg <- matrix(0, ns, 7L, dimnames = list(NULL, c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")))
  list(q = matrix(q0, ns, length(q0), byrow = TRUE), diag = dg, step_size = NA_real_, inv_metric = numeric(0), n_warmup_saved = 0L,
       time = c(warmup = 0, sampling = 0))
}

# Run `f(chain_id)` for every chain with independent RNG streams.
# Without a seed, one is drawn from the caller's stream (which therefore
# advances, so repeated calls differ); the drawn seed is returned as the
# "seed" attribute and recorded in the fit.
# progress lines of a chain: forked chains write them to a file that the
# parent prints, since R GUIs such as RStudio do not show a child's output
.progress <- new.env(parent = emptyenv())
chain_message <- function(text) {
  if (is.null(.progress$file)) message(text) else cat(text, "\n", sep = "", file = .progress$file, append = TRUE)
}

# run f(1), ..., f(n) in at most `cores` forked processes, relaying progress
run_forked <- function(n, cores, f) {
  files <- vapply(seq_len(n), function(i) tempfile(sprintf("chain%d_", i)), "")
  on.exit(unlink(files), add = TRUE)
  seen <- integer(n)
  relay <- function() for (i in seq_len(n)) if (file.exists(files[i])) {
    lines <- readLines(files[i], warn = FALSE)
    if (length(lines) > seen[i]) {
      # an unfinished last line is printed on the next pass
      done <- length(lines) - if (grepl("\n$", readChar(files[i], file.size(files[i]), useBytes = TRUE))) 0L else 1L
      if (done > seen[i]) message(paste(lines[(seen[i] + 1L):done], collapse = "\n"))
      seen[i] <<- max(seen[i], done)
    }
  }
  out <- vector("list", n)
  todo <- seq_len(n)
  running <- list()
  while (length(todo) || length(running)) {
    while (length(todo) && length(running) < cores) {
      i <- todo[1L]
      todo <- todo[-1L]
      job <- parallel::mcparallel({ .progress$file <- files[i]; f(i) }, silent = FALSE)
      running[[as.character(job$pid)]] <- list(job = job, i = i)
    }
    res <- parallel::mccollect(lapply(running, `[[`, "job"), wait = FALSE, timeout = 0.25)
    relay()
    for (pid in names(res)) {
      out[[running[[pid]]$i]] <- res[[pid]]
      running[[pid]] <- NULL
    }
  }
  relay()
  out
}

# forked processes where available; socket clusters on Windows (or when
# options(bayescompR.parallel = "psock"))
parallel_backend <- function() {
  b <- getOption("bayescompR.parallel", if (.Platform$OS.type == "windows") "psock" else "fork")
  match.arg(b, c("fork", "psock"))
}

# run f(1), ..., f(n) on a socket cluster of `cores` R processes, which
# load the installed bayescompR; progress is printed when the chains finish
run_psock <- function(n, cores, f) {
  files <- vapply(seq_len(n), function(i) tempfile(sprintf("chain%d_", i)), "")
  on.exit(unlink(files), add = TRUE)
  cl <- parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  out <- parallel::clusterApplyLB(cl, seq_len(n), function(i, chain_fun, files) {
    assign("file", files[i], envir = get(".progress", envir = asNamespace("bayescompR")))
    try(chain_fun(i), silent = TRUE)
  }, chain_fun = f, files = files)
  for (fl in files) if (file.exists(fl)) message(paste(readLines(fl, warn = FALSE), collapse = "\n"))
  out
}

bind_chain_id <- function(model, id) {
  if (identical(model$env$V.CHAIN_ID, as.integer(id))) return(invisible())
  d <- model$data
  d$CHAIN_ID <- as.integer(id)
  inst <- with_seed(model$compiled$seed, instantiate(model$compiled, d))
  for (nm in ls(inst$env, all.names = TRUE)) {
    v <- tryCatch(get(nm, envir = inst$env), error = function(e) e)
    # a hoisted expression that fails stays a promise that fails where used
    if (inherits(v, "error")) eval(call("delayedAssign", nm, call("get", nm, envir = inst$env), model$env, model$env))
    else assign(nm, v, envir = model$env)
  }
  invisible()
}

run_chains <- function(model, chains, seed, cores, f) {
  old_kind <- RNGkind()
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit({
    RNGkind(old_kind[1L], old_kind[2L], old_kind[3L])
    if (had) assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(seed)
  streams <- vector("list", chains)
  s <- get(".Random.seed", envir = globalenv())
  for (i in seq_len(chains)) {
    s <- parallel::nextRNGStream(s)
    streams[[i]] <- s
  }
  rebind <- isTRUE(model$compiled$chain_id_auto) && chains > 1L
  if (rebind) on.exit(bind_chain_id(model, 1L), add = TRUE)
  one <- function(i) {
    # each chain sees its own CHAIN_ID, so transformed data is rebuilt
    if (rebind) bind_chain_id(model, i)
    assign(".Random.seed", streams[[i]], envir = globalenv())
    profile_reset()
    r <- f(i)
    if (is.list(r)) r$profile <- profile_table()
    r
  }
  if (cores > 1L && chains > 1L) {
    out <- if (parallel_backend() == "fork") run_forked(chains, min(cores, chains), one)
           else run_psock(chains, min(cores, chains), one)
    bad <- vapply(out, inherits, TRUE, what = "try-error")
    if (any(bad)) stop(sprintf("chain %d failed: %s", which(bad)[1L], attr(out[[which(bad)[1L]]], "condition")$message), call. = FALSE)
  } else {
    out <- lapply(seq_len(chains), one)
  }
  attr(out, "seed") <- seed
  out
}

# Stan's initialization: up to 100 attempts at a point with finite log
# density and gradient.
# gradient = FALSE: the method needs only log density values, so a model
# whose gradient is unavailable can still be initialized
find_init <- function(model, init, chain_id, gradient = TRUE) {
  n <- model$dim
  if (n == 0L) return(numeric(0))
  if (inherits(init, "bayescomp_fit")) {
    # e.g. a Pathfinder or ADVI fit: start from one of its draws
    u <- do.call(rbind, init$unconstrained)
    if (ncol(u) != n) stop("the initialization fit has a different number of parameters", call. = FALSE)
    pot <- init_check_fun(model, gradient)
    for (attempt in seq_len(100L)) {
      q <- u[sample.int(nrow(u), 1L), ]
      r <- pot(q)
      if (is.finite(r$V) && all(is.finite(r$g))) return(q)
    }
    stop(sprintf("chain %d: no draw of the initialization fit has a finite log density", chain_id), call. = FALSE)
  }
  user <- init
  if (is.character(user)) user <- bayescomp_read_data(user[[(chain_id - 1L) %% length(user) + 1L]])
  if (is.function(user)) user <- user()
  if (is.list(user) && length(user) && is.list(user[[1L]]) && is.null(names(user))) user <- user[[chain_id]]
  radius <- 2
  if (is.numeric(user) && length(user) == 1L) {
    radius <- user
    user <- NULL
  }
  pot <- init_check_fun(model, gradient)
  for (attempt in seq_len(100L)) {
    q <- if (radius == 0) numeric(n) else stats::runif(n, -radius, radius)
    if (is.list(user) && length(user)) {
      cur <- write_array(model, q, FALSE, FALSE, print = FALSE)[[1L]]
      for (nm in names(user)) {
        if (!(nm %in% model$params$pars)) next
        cur[[nm]] <- user[[nm]]
      }
      q <- tryCatch(bayescomp_unconstrain(model, cur), error = function(e) {
        stop(sprintf("chain %d: invalid initial values: %s", chain_id, conditionMessage(e)), call. = FALSE)
      })
    }
    r <- pot(q)
    if (is.finite(r$V) && all(is.finite(r$g))) return(q)
    if (is.list(user) && length(user) && length(intersect(names(user), model$params$pars)) == length(model$params$pars))
      stop(sprintf("chain %d: the log density or its gradient is not finite at the supplied initial values", chain_id), call. = FALSE)
    if (radius == 0) break
  }
  stop(sprintf("chain %d: initialization failed after 100 attempts; try smaller initial values or check the model", chain_id), call. = FALSE)
}

init_check_fun <- function(model, gradient) {
  if (gradient) return(hmc_potential(model))
  lp <- lp_value_fun(model)
  function(q) list(V = -lp(q), g = 0)
}

# constrained parameters, transformed parameters and generated quantities
chain_write_draws <- function(model, qs, print) {
  nm <- bayescomp_param_names(model)
  out <- matrix(NA_real_, nrow(qs), length(nm), dimnames = list(NULL, nm))
  for (i in seq_len(nrow(qs))) {
    r <- tryCatch(write_array(model, qs[i, ], TRUE, TRUE, print = print),
                  bayescomp_fatal_error = function(e) stop(e),
                  bayescomp_error = function(e) {
                    warning(sprintf("generated quantities failed at draw %d: %s", i, conditionMessage(e)), call. = FALSE)
                    NULL
                  })
    if (is.null(r)) {
      # keep the parameters (and transformed parameters, if they evaluate);
      # only the failed block is left missing
      r <- tryCatch(write_array(model, qs[i, ], TRUE, FALSE, print = FALSE),
                    bayescomp_fatal_error = function(e) stop(e), bayescomp_error = function(e) NULL)
      if (is.null(r)) r <- write_array(model, qs[i, ], FALSE, FALSE, print = FALSE)
      v <- flatten_draw(c(r[[1L]], r[[2L]]))
      out[i, seq_along(v)] <- v
      next
    }
    out[i, ] <- flatten_draw(c(r[[1L]], r[[2L]], r[[3L]]))
  }
  out
}
