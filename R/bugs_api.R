# Running BUGS models: interfaces with the arguments and results of R2jags
# (jags()), R2OpenBUGS (bugs()) and rjags (jags.model(), update(),
# coda.samples(), jags.samples()). No external program is needed.

bugs_model_source <- function(model.file) {
  if (is.function(model.file) || is.call(model.file) || is.expression(model.file)) return(bugs_from_r(model.file))
  if (inherits(model.file, "connection")) return(bugs_parse(readLines(model.file, warn = FALSE)))
  if (is.character(model.file)) {
    if (length(model.file) == 1L && file.exists(path.expand(model.file)))
      return(bugs_parse(readLines(path.expand(model.file), warn = FALSE), file = model.file))
    if (any(grepl("model\\s*\\{", model.file))) return(bugs_parse(model.file))
    if (length(model.file) == 1L) stop(sprintf("the model file '%s' does not exist", model.file), call. = FALSE)
  }
  stop("model.file must be a file name, BUGS code, a connection or an R function holding the model", call. = FALSE)
}

# data as a list, a vector of names of objects, or one or more files
bugs_data_arg <- function(data, env) {
  if (is.null(data)) return(list())
  if (is.character(data) && all(file.exists(path.expand(data)))) {
    out <- list()
    for (f in data) { v <- bugs_read_values(path.expand(f), "data"); out[names(v)] <- v }
    return(out)
  }
  if (is.character(data)) {
    miss <- data[!vapply(data, exists, TRUE, envir = env)]
    if (length(miss)) stop(sprintf("data object '%s' not found", miss[1L]), call. = FALSE)
    return(mget(data, envir = env, inherits = TRUE))
  }
  if (is.environment(data)) return(as.list(data))
  if (!is.list(data)) stop("data must be a list, a vector of object names or a file", call. = FALSE)
  if (length(data) && (is.null(names(data)) || any(!nzchar(names(data))))) stop("every element of data must be named", call. = FALSE)
  data
}

# the initial values of chain i from NULL, a list of lists, one list, a
# function or files; functions are called inside the chain's random stream
bugs_inits_for <- function(inits, i) {
  if (is.null(inits)) return(NULL)
  if (is.function(inits)) return(if (length(formals(inits))) inits(i) else inits())
  if (is.character(inits)) inits <- lapply(path.expand(inits), bugs_read_values, what = "inits")
  if (is.list(inits) && length(inits) && all(vapply(inits, function(x) is.list(x) || is.null(x), TRUE)) &&
      (is.null(names(inits)) || all(!nzchar(names(inits))))) {
    # fewer sets than chains are recycled
    return(inits[[(i - 1L) %% length(inits) + 1L]])
  }
  if (is.list(inits)) return(inits)
  stop("inits must be NULL, a function, a list or a list of lists", call. = FALSE)
}

# the variables monitored by default: the unobserved stochastic variables
bugs_default_monitor <- function(M) {
  U <- M$nodes$kind == "stochastic" & !M$nodes$obs
  unique(M$nodes$var[U])
}

# a fresh copy of the values of a compiled model for each chain
bugs_snapshot <- function(M) as.list(M$S)
bugs_restore <- function(M, snap) { for (v in names(snap)) assign(v, snap[[v]], envir = M$S); M }

bugs_run <- function(prog, data, inits, monitor, n.chains, n.adapt, n.burnin, n.iter, n.thin, deviance, dialect,
                     seed, cores, refresh) {
  M <- bugs_compile(prog, data, dialect)
  if (is.null(monitor)) monitor <- bugs_default_monitor(M)
  monitor <- setdiff(unique(monitor), "deviance")
  mon <- bugs_monitor(M, monitor)
  snap <- bugs_snapshot(M)
  # OpenBUGS's gen inits and NIMBLE draw missing initial values from the
  # priors; JAGS starts from typical values
  how <- if (dialect %in% c("openbugs", "nimble")) "prior" else "typical"
  t0 <- Sys.time()
  res <- run_chains(list(), n.chains, seed, cores, function(id) {
    bugs_restore(M, snap)
    bugs_chain(M, bugs_inits_for(inits, id), n.adapt, n.burnin, n.iter, n.thin, monitor, deviance, how, refresh, id)
  })
  list(draws = lapply(res, `[[`, "draws"), means = lapply(res, `[[`, "means"), samplers = res[[1L]]$samplers,
       model = M, mon = mon, prog = prog, data = data, dialect = dialect,
       seed = attr(res, "seed"), elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# the sims.array / sims.list / summary structure of R2WinBUGS and R2jags;
# R2jags sorts the variables by name, R2OpenBUGS keeps the monitored order
bugs_output <- function(run, n.chains, n.iter, n.burnin, n.thin, DIC, model.file, dic, sort_vars) {
  mon <- run$mon
  if ("deviance" %in% colnames(run$draws[[1L]])) mon$deviance <- list(var = "deviance", names = "deviance", index = NULL, dim = NULL)
  roots <- names(mon)
  roots <- if (sort_vars) sort(roots) else c(setdiff(roots, "deviance"), intersect(roots, "deviance"))
  vn <- unlist(lapply(roots, function(r) mon[[r]]$names), use.names = FALSE)
  root <- rep(roots, vapply(roots, function(r) length(mon[[r]]$names), 1L))
  n.keep <- nrow(run$draws[[1L]])
  sa <- array(NA_real_, c(n.keep, n.chains, length(vn)), dimnames = list(NULL, NULL, vn))
  for (c in seq_len(n.chains)) sa[, c, ] <- run$draws[[c]][, vn, drop = FALSE]
  sm <- matrix(sa, n.keep * n.chains, length(vn), dimnames = list(NULL, vn))
  shape_of <- function(r) {
    m <- mon[[r]]
    if (!is.null(m$dim)) m$dim else if (length(m$names) == 1L) NULL else length(m$names)
  }
  sl <- lapply(roots, function(r) {
    x <- sm[, root == r, drop = FALSE]
    d <- shape_of(r)
    if (is.null(d)) as.vector(x) else if (length(d) == 1L) unname(x) else array(x, c(nrow(x), d))
  })
  names(sl) <- roots
  stat <- function(x) {
    xm <- matrix(x, n.keep)
    c(mean = mean(x), sd = stats::sd(x), stats::quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)),
      Rhat = if (n.chains > 1L) bayescomp_rhat(xm) else NA_real_, n.eff = round(bayescomp_ess_bulk(xm)))
  }
  summ <- t(vapply(seq_along(vn), function(j) stat(sa[, , j]), numeric(9L)))
  dimnames(summ) <- list(vn, c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "Rhat", "n.eff"))
  reshape <- function(v) {
    out <- lapply(roots, function(r) {
      x <- unname(v[root == r]); d <- shape_of(r)
      if (is.null(d) || length(d) == 1L) x else array(x, d)
    })
    names(out) <- roots
    out
  }
  last <- lapply(seq_len(n.chains), function(c) reshape(sa[n.keep, c, ])[setdiff(roots, "deviance")])
  out <- list(n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.keep = n.keep,
              n.sims = n.keep * n.chains, sims.array = sa, sims.list = sl, sims.matrix = sm, summary = summ,
              mean = reshape(summ[, "mean"]), sd = reshape(summ[, "sd"]), median = reshape(summ[, "50%"]),
              root.short = roots, long.short = lapply(roots, function(r) which(root == r)),
              dimension.short = vapply(roots, function(r) length(shape_of(r) %||% integer()), 1L),
              indexes.short = lapply(roots, function(r) {
                ix <- mon[[r]]$index
                if (is.null(ix) || is.null(shape_of(r))) NULL else lapply(seq_len(nrow(ix)), function(k) unname(ix[k, ]))
              }),
              last.values = last, program = "bayescompR", model.file = model.file, isDIC = DIC, DICbyR = TRUE)
  if (DIC && "deviance" %in% roots) {
    dev <- matrix(sa[, , "deviance"], n.keep)
    if (dic == "openbugs" && !is.null(run$Dhat)) {
      out$pD <- mean(dev) - run$Dhat
      out$DICbyR <- FALSE
    } else {
      # R2jags: half the mean of the within-chain variances
      out$pD <- mean(apply(dev, 2L, stats::var)) / 2
    }
    out$DIC <- mean(dev) + out$pD
  }
  class(out) <- "bugs"
  out
}

bugs_in_dir <- function(dir, expr) {
  if (is.null(dir)) return(expr)
  old <- setwd(dir)
  on.exit(setwd(old))
  expr
}

#' Run JAGS models without JAGS
#'
#' Runs a model written in the BUGS language as JAGS understands it, with the
#' arguments and the result of `R2jags::jags()`. The model is compiled and
#' sampled in R by bayescompR: neither JAGS nor R2jags is needed. The draws
#' come from the same posterior as JAGS gives, not the same random numbers.
#'
#' @param data A named list of data, a character vector naming objects to
#'   take from the calling environment, or data files in R dump or BUGS list
#'   format. Missing values (`NA`) are unobserved and sampled.
#' @param inits Initial values: `NULL` (JAGS's typical values), a list of
#'   lists (one per chain, recycled), a single list used for every chain, or
#'   a function returning a list (called once per chain).
#' @param parameters.to.save Names of the variables to monitor, or parts of
#'   them such as `"b[2]"`. `NULL` monitors every unobserved stochastic
#'   variable. `"deviance"` monitors the deviance.
#' @param model.file The model: a file name, a character vector of BUGS code,
#'   or an R function whose body is the model (as in R2jags).
#' @param n.chains,n.iter,n.burnin,n.thin Chains, iterations per chain
#'   (including burn-in), burn-in and thinning, with R2jags's defaults.
#' @param DIC Monitor the deviance and compute DIC with R2jags's penalty: half
#'   the mean within-chain variance of the deviance.
#' @param pD Accepted for compatibility; the penalty of `dic.samples()` is
#'   not computed.
#' @param working.directory Directory of the model and data files.
#' @param n.adapt Adaptation iterations run before `n.iter`.
#' @param jags.seed Seed of the random numbers.
#' @param cores Chains run in parallel on this many cores.
#' @param refresh Print progress every `refresh` iterations; 0 is silent.
#' @param digits Accepted for compatibility.
#' @param ... Other arguments of `R2jags::jags()` (`progress.bar`, `quiet`,
#'   `RNGname`, `jags.module`) are accepted and ignored.
#' @return An object of class `rjags` with a `BUGSoutput` element of class
#'   `bugs` holding `sims.list`, `sims.array`, `sims.matrix`, `summary`,
#'   `mean`, `sd`, `median`, `last.values`, `DIC` and `pD`.
#' @seealso [bugs()], [jags.model()], [nimbleMCMC()] and `vignette("bugs")`.
#' @examples
#' code <- "model {
#'   for (i in 1:N) { y[i] ~ dnorm(mu, tau) }
#'   mu ~ dnorm(0, 1.0E-4)
#'   tau ~ dgamma(0.01, 0.01)
#' }"
#' fit <- jags(list(y = c(1.2, 0.8, 1.9, 1.1), N = 4), NULL, c("mu", "tau"), code,
#'             n.chains = 2, n.iter = 1000, jags.seed = 1)
#' fit$BUGSoutput$summary
#' @export
jags <- function(data, inits = NULL, parameters.to.save = NULL, model.file = "model.bug", n.chains = 3, n.iter = 2000,
                 n.burnin = floor(n.iter / 2), n.thin = max(1, floor((n.iter - n.burnin) / 1000)), DIC = TRUE, pD = FALSE,
                 working.directory = NULL, jags.seed = 123, n.adapt = if (n.burnin > 0) 0L else 100L,
                 cores = getOption("mc.cores", 1L), refresh = 0L, digits = 5, ...) {
  env <- parent.frame()
  if (isTRUE(pD)) warning("pD = TRUE is not supported; DIC uses the variance penalty", call. = FALSE)
  prog <- bugs_in_dir(working.directory, bugs_model_source(model.file))
  data <- bugs_in_dir(working.directory, bugs_data_arg(data, env))
  if (is.character(inits)) inits <- bugs_in_dir(working.directory, lapply(path.expand(inits), bugs_read_values, what = "inits"))
  n.chains <- as.integer(n.chains); n.iter <- as.integer(n.iter); n.burnin <- as.integer(n.burnin); n.thin <- as.integer(n.thin)
  if (n.burnin >= n.iter) stop("n.burnin must be smaller than n.iter", call. = FALSE)
  # JAGS adapts during the first part of the burn-in
  na <- min(n.burnin, 1000L)
  dev <- DIC || "deviance" %in% parameters.to.save
  run <- bugs_run(prog, data, inits, parameters.to.save, n.chains, as.integer(n.adapt) + na, n.burnin - na, n.iter - na, n.thin,
                  dev, "jags", jags.seed, cores, refresh)
  bo <- bugs_output(run, n.chains, n.iter, n.burnin, n.thin, DIC, if (is.character(model.file)) model.file else "", "jags", TRUE)
  structure(list(model = run$model, BUGSoutput = bo, parameters.to.save = c(names(run$mon), if (dev) "deviance"),
                 model.file = model.file, n.iter = n.iter, DIC = DIC, seed = run$seed, elapsed = run$elapsed),
            class = "rjags")
}

#' Run OpenBUGS, WinBUGS and MultiBUGS models without BUGS
#'
#' Runs a model written in the BUGS language as OpenBUGS understands it, with
#' the arguments and the result of `R2OpenBUGS::bugs()`. Nothing outside R
#' is needed.
#'
#' @inheritParams jags
#' @param inits Initial values as in [jags()]; `NULL` draws them from the
#'   priors (OpenBUGS's gen inits).
#' @param n.chains Number of chains.
#' @param n.burnin Burn-in iterations per chain, discarded.
#' @param n.iter Iterations per chain, including the burn-in.
#' @param n.thin Thinning: as in OpenBUGS, every iteration after the burn-in
#'   runs `n.thin` updates and keeps the last, so `n.iter - n.burnin` draws
#'   are kept per chain.
#' @param DIC Monitor the deviance and compute DIC with OpenBUGS's
#'   `pD = mean(deviance) - deviance(posterior means)`. When an unobserved
#'   discrete node is a parent of the data, the plug-in deviance is not
#'   defined and half the variance of the deviance is used instead.
#' @param codaPkg Write CODA files to `working.directory` (or a temporary
#'   directory) and return their names instead of the results.
#' @param bugs.seed Seed of the random numbers.
#' @param debug,digits,summary.only,saveExec,restart,clearWD,over.relax,save.history
#'   Accepted for compatibility and ignored.
#' @param ... Other arguments of `R2OpenBUGS::bugs()` (`OpenBUGS.pgm`,
#'   `useWINE`, `WINE`, `newWINE`, `WINEPATH`) are accepted and ignored.
#' @return An object of class `bugs`, as R2OpenBUGS returns.
#' @examples
#' code <- "model {
#'   for (i in 1:N) { y[i] ~ dnorm(mu, tau) }
#'   mu ~ dnorm(0, 1.0E-4)
#'   tau ~ dgamma(0.01, 0.01)
#' }"
#' fit <- bugs(list(y = c(1.2, 0.8, 1.9, 1.1), N = 4), list(list(mu = 0, tau = 1), list(mu = 1, tau = 2)),
#'             c("mu", "tau"), 1000, code, n.chains = 2)
#' fit$summary
#' @export
bugs <- function(data, inits, parameters.to.save, n.iter, model.file = "model.txt", n.chains = 3, n.burnin = floor(n.iter / 2),
                 n.thin = 1, saveExec = FALSE, restart = FALSE, debug = FALSE, DIC = TRUE, digits = 5, codaPkg = FALSE,
                 working.directory = NULL, clearWD = FALSE, bugs.seed = 1, summary.only = FALSE, save.history = FALSE,
                 over.relax = FALSE, cores = getOption("mc.cores", 1L), refresh = 0L, ...) {
  env <- parent.frame()
  if (missing(inits)) inits <- NULL
  if (missing(parameters.to.save)) parameters.to.save <- NULL
  prog <- bugs_in_dir(working.directory, bugs_model_source(model.file))
  data <- bugs_in_dir(working.directory, bugs_data_arg(data, env))
  if (is.character(inits)) inits <- bugs_in_dir(working.directory, lapply(path.expand(inits), bugs_read_values, what = "inits"))
  n.chains <- as.integer(n.chains); n.iter <- as.integer(n.iter); n.burnin <- as.integer(n.burnin); n.thin <- as.integer(n.thin)
  if (n.burnin >= n.iter) stop("n.burnin must be smaller than n.iter", call. = FALSE)
  na <- min(n.burnin, 1000L)
  dev <- DIC || "deviance" %in% parameters.to.save
  # the burn-in, then n.iter - n.burnin kept iterations of n.thin updates each
  total <- (n.burnin - na) + (n.iter - n.burnin) * n.thin
  run <- bugs_run(prog, data, inits, parameters.to.save, n.chains, na, n.burnin - na, total, n.thin, dev, "openbugs",
                  bugs.seed, cores, refresh)
  if (DIC) run$Dhat <- bugs_dhat(run)
  bo <- bugs_output(run, n.chains, n.iter, n.burnin, n.thin, DIC, if (is.character(model.file)) model.file else "", "openbugs", FALSE)
  if (isTRUE(codaPkg)) return(bugs_write_coda(bo, working.directory %||% tempdir()))
  bo
}

# OpenBUGS's CODA files: an index of variables and one file per chain
bugs_write_coda <- function(bo, dir) {
  vn <- dimnames(bo$sims.array)[[3L]]
  n <- bo$n.keep
  writeLines(sprintf("%s\t%d\t%d", vn, (seq_along(vn) - 1L) * n + 1L, seq_along(vn) * n), file.path(dir, "CODAindex.txt"))
  vapply(seq_len(bo$n.chains), function(c) {
    f <- file.path(dir, sprintf("CODAchain%d.txt", c))
    writeLines(sprintf("%d\t%.10g", rep(bo$n.burnin + seq_len(n), length(vn)), as.vector(bo$sims.array[, c, ])), f)
    f
  }, "")
}

# the deviance at the posterior means of the unobserved stochastic nodes
# (OpenBUGS's Dhat); NULL when it is not defined (an unobserved discrete
# node is a parent of the data) or not finite
bugs_dhat <- function(run) {
  M <- bugs_compile(run$prog, run$data, run$dialect)
  U <- which(M$nodes$kind == "stochastic" & !M$nodes$obs)
  if (any(vapply(U, function(i) isTRUE(bugs_dist_get(M$stmts[[M$nodes$stmt[i]]]$dist, M$dialect)$discrete), TRUE))) return(NULL)
  for (v in names(run$means[[1L]])) {
    a <- get(v, envir = M$S)
    miss <- is.na(a)
    a[miss] <- (Reduce(`+`, lapply(run$means, `[[`, v)) / length(run$means))[miss]
    assign(v, a, envir = M$S)
  }
  for (i in M$order) if (M$nodes$kind[i] == "logical" && !M$const[i]) M$stmts[[M$nodes$stmt[i]]]$fn$write(M$nodes$row[i])
  d <- bugs_deviance_fun(M)()
  if (is.finite(d)) d
}

#' @export
print.bugs <- function(x, digits = 3, ...) {
  cat(sprintf("Inference for a BUGS model (bayescompR): %d chains, each with %d iterations (first %d discarded)",
              x$n.chains, x$n.iter, x$n.burnin))
  if (x$n.thin > 1L) cat(sprintf(", n.thin = %d", x$n.thin))
  cat(sprintf("\n n.sims = %d iterations saved\n", x$n.sims))
  print(round(x$summary, digits))
  if (isTRUE(x$isDIC) && !is.null(x$DIC)) cat(sprintf("\npD = %.1f and DIC = %.1f\n", x$pD, x$DIC))
  invisible(x)
}

#' @export
print.rjags <- function(x, digits = 3, ...) {
  print(x$BUGSoutput, digits = digits)
  invisible(x)
}

# ---- rjags -------------------------------------------------------------------------------

#' rjags interface: models without JAGS
#'
#' `jags.model()` compiles a model and initializes its chains, `update()`
#' runs iterations without keeping them, and `coda.samples()` and
#' `jags.samples()` run iterations and return the draws, with the arguments
#' and results of the rjags functions of the same names. The model object
#' keeps the state of its chains between calls.
#'
#' @param file The model: a file name, BUGS code, a connection or an R
#'   function.
#' @param data,inits Data and initial values, as in [jags()].
#' @param n.chains Number of chains.
#' @param n.adapt Adaptation iterations run at once.
#' @param quiet Accepted for compatibility.
#' @param seed Seed of the random numbers; every chain gets its own stream.
#' @param object,model A model from `jags.model()`.
#' @param n.iter Iterations to run.
#' @param variable.names Variables to monitor, or parts of them.
#' @param thin Thinning.
#' @param by,progress.bar Accepted for compatibility.
#' @param x A model from `jags.model()`.
#' @param ... Unused.
#' @return `jags.model()`: an object of class `jags`. `coda.samples()`: an
#'   `mcmc.list` (a list of `mcmc` matrices, as the coda package defines).
#'   `jags.samples()`: a list of `mcarray` arrays with dimensions variable,
#'   iteration and chain.
#' @examples
#' m <- jags.model(textConnection("model { y ~ dnorm(mu, 1)  mu ~ dnorm(0, 0.01) }"),
#'                 data = list(y = 1.5), n.chains = 2, seed = 1)
#' update(m, 500)
#' s <- coda.samples(m, "mu", n.iter = 1000)
#' mean(unlist(lapply(s, as.vector)))
#' @export
jags.model <- function(file, data = NULL, inits = NULL, n.chains = 1, n.adapt = 1000, quiet = FALSE, seed = NULL) {
  env <- parent.frame()
  prog <- bugs_model_source(file)
  data <- bugs_data_arg(data, env)
  n.chains <- as.integer(n.chains)
  M0 <- bugs_compile(prog, data, "jags")
  rngs <- bugs_streams(n.chains, seed)
  # each chain keeps its own model values and random stream
  chains <- lapply(seq_len(n.chains), function(i) {
    M <- if (i == 1L) M0 else bugs_compile(prog, data, "jags")
    bugs_with_stream(rngs, i, bugs_chain_init(M, bugs_inits_for(inits, i), "typical"))
  })
  st <- new.env(parent = emptyenv())
  st$chains <- chains; st$rngs <- rngs; st$iter <- 0L
  out <- structure(list(state = st, model = M0), class = "jags")
  if (n.adapt > 0) bugs_jags_advance(out, as.integer(n.adapt), adapting = TRUE)
  out
}

# random number streams of the chains of a jags model, and evaluating code
# in one of them
bugs_streams <- function(n, seed) {
  e <- new.env(parent = emptyenv())
  e$s <- with_rng_state({
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    set.seed(seed %||% sample.int(.Machine$integer.max, 1L))
    s <- get(".Random.seed", envir = globalenv())
    out <- vector("list", n)
    for (i in seq_len(n)) { s <- parallel::nextRNGStream(s); out[[i]] <- s }
    out
  })
  e
}
bugs_with_stream <- function(rngs, i, expr) with_rng_state({
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  assign(".Random.seed", rngs$s[[i]], envir = globalenv())
  val <- force(expr)
  rngs$s[[i]] <- get(".Random.seed", envir = globalenv())
  val
})

# evaluate code and restore the caller's generator and seed afterwards
with_rng_state <- function(expr) {
  old_kind <- RNGkind()
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit({
    RNGkind(old_kind[1L], old_kind[2L], old_kind[3L])
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  })
  expr
}

bugs_jags_advance <- function(object, n, adapting = FALSE, monitor = NULL, thin = 1L) {
  st <- object$state
  draws <- lapply(seq_along(st$chains), function(i) {
    r <- bugs_with_stream(st$rngs, i, {
      mon <- if (length(monitor)) bugs_monitor(st$chains[[i]]$M, setdiff(monitor, "deviance")) else list()
      bugs_chain_run(st$chains[[i]], n, adapting, mon, "deviance" %in% monitor, 0L, thin)
    })
    st$chains[[i]] <- r$state
    r$draws
  })
  st$iter <- st$iter + n
  draws
}

#' @rdname jags.model
#' @importFrom stats update
#' @export
update.jags <- function(object, n.iter = 1, by = NULL, progress.bar = NULL, ...) {
  bugs_jags_advance(object, as.integer(n.iter))
  invisible(NULL)
}

#' @rdname jags.model
#' @export
coda.samples <- function(model, variable.names, n.iter, thin = 1, ...) {
  start <- model$state$iter
  d <- bugs_jags_advance(model, as.integer(n.iter), monitor = variable.names, thin = as.integer(thin))
  out <- lapply(d, function(x) {
    it <- start + thin * seq_len(nrow(x))
    structure(x, mcpar = c(min(it), max(it), thin), class = "mcmc")
  })
  structure(out, class = "mcmc.list")
}

#' @rdname jags.model
#' @export
jags.samples <- function(model, variable.names, n.iter, thin = 1, ...) {
  mon <- bugs_monitor(model$state$chains[[1L]]$M, setdiff(variable.names, "deviance"))
  d <- bugs_jags_advance(model, as.integer(n.iter), monitor = variable.names, thin = as.integer(thin))
  nk <- nrow(d[[1L]]); nc <- length(d)
  roots <- c(names(mon), intersect(variable.names, "deviance"))
  out <- lapply(roots, function(r) {
    nm <- if (r == "deviance") "deviance" else mon[[r]]$names
    dd <- if (r != "deviance" && !is.null(mon[[r]]$dim)) mon[[r]]$dim else length(nm)
    a <- array(NA_real_, c(length(nm), nk, nc))
    for (c in seq_len(nc)) a[, , c] <- t(d[[c]][, nm, drop = FALSE])
    structure(array(a, c(dd, nk, nc)), class = "mcarray")
  })
  names(out) <- roots
  out
}

#' @rdname jags.model
#' @export
print.jags <- function(x, ...) {
  M <- x$model
  cat(sprintf("JAGS model (bayescompR): %d chains, %d iterations run\n", length(x$state$chains), x$state$iter))
  cat(sprintf("  %d observed and %d unobserved stochastic nodes, %d logical nodes\n",
              sum(M$nodes$kind == "stochastic" & M$nodes$obs), sum(M$nodes$kind == "stochastic" & !M$nodes$obs),
              sum(M$nodes$kind == "logical")))
  invisible(x)
}
