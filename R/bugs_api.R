# Running BUGS models: interfaces with the arguments and results of R2jags
# (jags()), R2OpenBUGS and R2WinBUGS (bugs()), and rjags (jags.model(),
# update(), coda.samples()). No external program is needed.

bugs_model_source <- function(model.file, env = parent.frame(2L)) {
  if (is.function(model.file) || is.call(model.file) || is.expression(model.file)) return(bugs_from_r(model.file))
  if (inherits(model.file, "connection")) return(bugs_parse(readLines(model.file, warn = FALSE)))
  if (is.character(model.file)) {
    if (length(model.file) == 1L && !grepl("[{~]", model.file) && file.exists(model.file))
      return(bugs_parse(readLines(model.file, warn = FALSE), file = model.file))
    if (any(grepl("model\\s*\\{", model.file))) return(bugs_parse(model.file))
    if (length(model.file) == 1L && !file.exists(model.file)) stop(sprintf("the model file '%s' does not exist", model.file), call. = FALSE)
  }
  stop("model.file must be a file name, BUGS code, a connection or an R function holding the model", call. = FALSE)
}

# data as a list, a vector of names of objects, or a file
bugs_data_arg <- function(data, env) {
  if (is.null(data)) return(list())
  if (is.character(data) && length(data) == 1L && file.exists(data)) return(bugs_read_values(data, "data"))
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

# per chain initial values from NULL, a list of lists, one list, a function
# or files
bugs_inits_arg <- function(inits, n.chains) {
  if (is.null(inits)) return(rep(list(NULL), n.chains))
  if (is.function(inits)) {
    return(lapply(seq_len(n.chains), function(i) if (length(formals(inits))) inits(i) else inits()))
  }
  if (is.character(inits)) inits <- lapply(inits, bugs_read_values, what = "inits")
  if (is.list(inits) && length(inits) && all(vapply(inits, function(x) is.list(x) || is.null(x), TRUE)) &&
      (is.null(names(inits)) || all(!nzchar(names(inits))))) {
    # fewer sets than chains are recycled
    return(rep_len(inits, n.chains))
  }
  if (is.list(inits)) return(rep(list(inits), n.chains))
  stop("inits must be NULL, a function, a list or a list of lists", call. = FALSE)
}

# the variables to monitor by default: unobserved stochastic variables
bugs_default_monitor <- function(M) {
  U <- M$nodes$kind == "stochastic" & !M$nodes$obs
  unique(M$nodes$var[U])
}

bugs_run <- function(prog, data, inits, monitor, n.chains, n.adapt, n.burnin, n.iter, n.thin, DIC, dialect,
                     seed, cores, refresh) {
  M0 <- bugs_compile(prog, data, dialect)
  if (is.null(monitor)) monitor <- bugs_default_monitor(M0)
  monitor <- setdiff(unique(monitor), "deviance")
  if (DIC) monitor_dev <- TRUE else monitor_dev <- FALSE
  bugs_monitor(M0, monitor)   # checks the names
  ini <- bugs_inits_arg(inits, n.chains)
  how <- if (dialect == "jags") "typical" else "prior"
  t0 <- Sys.time()
  res <- run_chains(list(), n.chains, seed, cores, function(id) {
    M <- bugs_compile(prog, data, dialect)
    r <- bugs_chain(M, ini[[id]], n.adapt, n.burnin, n.iter, n.thin, monitor, monitor_dev, how, refresh, id)
    r$last <- lapply(stats::setNames(bugs_default_monitor(M), bugs_default_monitor(M)), function(v) bugs_drop(get(v, envir = M$S)))
    r
  })
  draws <- lapply(res, `[[`, "draws")
  list(draws = draws, means = lapply(res, `[[`, "means"), last = lapply(res, `[[`, "last"), samplers = res[[1L]]$samplers, model = M0,
       seed = attr(res, "seed"), elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# the sims.array / sims.list / summary structure of R2WinBUGS and R2jags
bugs_output <- function(run, n.chains, n.iter, n.burnin, n.thin, DIC, model.file, dic_type) {
  vn <- colnames(run$draws[[1L]])
  n.keep <- nrow(run$draws[[1L]])
  root <- sub("\\[.*$", "", vn)
  ord <- order(root, seq_along(vn))
  ord <- unlist(lapply(sort(unique(root)), function(r) which(root == r)))
  vn <- vn[ord]; root <- root[ord]
  sa <- array(NA_real_, c(n.keep, n.chains, length(vn)), dimnames = list(NULL, NULL, vn))
  for (c in seq_len(n.chains)) sa[, c, ] <- run$draws[[c]][, ord, drop = FALSE]
  sm <- matrix(sa, n.keep * n.chains, length(vn), dimnames = list(NULL, vn))
  roots <- unique(root)
  dims <- lapply(roots, function(r) {
    if (r == "deviance") return(NULL)
    d <- run$model$dims[[r]]
    if (prod(d) == 1L && sum(root == r) == 1L && !grepl("\\[", vn[root == r][1L])) NULL else d
  })
  names(dims) <- roots
  sl <- lapply(roots, function(r) {
    x <- sm[, root == r, drop = FALSE]
    d <- dims[[r]]
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
  reshape <- function(col) {
    out <- lapply(roots, function(r) {
      v <- summ[root == r, col]
      d <- dims[[r]]
      if (is.null(d)) unname(v) else if (length(d) == 1L) unname(v) else array(v, d)
    })
    names(out) <- roots
    out
  }
  out <- list(n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.keep = n.keep,
              n.sims = n.keep * n.chains, sims.array = sa, sims.list = sl, sims.matrix = sm, summary = summ,
              mean = reshape("mean"), sd = reshape("sd"), median = reshape("50%"),
              root.short = roots, long.short = lapply(roots, function(r) which(root == r)),
              dimension.short = vapply(roots, function(r) length(dims[[r]]), 1L),
              indexes.short = lapply(roots, function(r) NULL), last.values = run$last,
              isDIC = DIC, DICbyR = TRUE, model.file = model.file, program = "bayescompR")
  if (DIC && "deviance" %in% roots) {
    dev <- sa[, , "deviance"]
    out$pD <- if (dic_type == "openbugs" && !is.null(run$Dhat)) mean(dev) - run$Dhat else stats::var(as.vector(dev)) / 2
    out$DIC <- mean(dev) + out$pD
  }
  class(out) <- "bugs"
  out
}

#' Run JAGS models without JAGS
#'
#' Runs a model written in the BUGS language as JAGS understands it, with the
#' arguments and the result of `R2jags::jags()`. The model is compiled and
#' sampled in R by bayescompR: neither JAGS nor R2jags is needed. The draws
#' come from the same posterior as JAGS gives, not the same random numbers.
#'
#' @param data A named list of data, or a character vector naming objects to
#'   take from the calling environment, or a data file in R dump or BUGS list
#'   format. Missing values (`NA`) are unobserved and sampled.
#' @param inits Initial values: `NULL` (JAGS's typical values), a list of
#'   lists (one per chain), a single list used for every chain, or a
#'   function returning a list.
#' @param parameters.to.save Names of the variables to monitor. `NULL`
#'   monitors every unobserved stochastic variable.
#' @param model.file The model: a file name, a character vector of BUGS code,
#'   or an R function whose body is the model (as in R2jags).
#' @param n.chains,n.iter,n.burnin,n.thin Chains, iterations per chain
#'   (including burn-in), burn-in and thinning, with R2jags's defaults.
#' @param n.adapt Adaptation iterations run before `n.iter`.
#' @param DIC Monitor the deviance and compute DIC with `pD = var(deviance)/2`.
#' @param jags.seed Seed of the random numbers.
#' @param cores Chains run in parallel on this many cores.
#' @param refresh Print progress every `refresh` iterations; 0 is silent.
#' @param ... Other arguments of `R2jags::jags()` (`working.directory`,
#'   `progress.bar`, `quiet`, `RNGname`, `jags.module`) are accepted and
#'   ignored.
#' @return An object of class `rjags` with a `BUGSoutput` element of class
#'   `bugs` holding `sims.list`, `sims.array`, `sims.matrix`, `summary`,
#'   `mean`, `sd`, `median`, `DIC` and `pD`.
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
                 n.burnin = floor(n.iter / 2), n.thin = max(1, floor((n.iter - n.burnin) / 1000)), DIC = TRUE,
                 n.adapt = if (n.burnin > 0) 0L else 100L, jags.seed = 123, cores = getOption("mc.cores", 1L), refresh = 0L, ...) {
  env <- parent.frame()
  prog <- bugs_model_source(model.file)
  data <- bugs_data_arg(data, env)
  n.chains <- as.integer(n.chains); n.iter <- as.integer(n.iter); n.burnin <- as.integer(n.burnin); n.thin <- as.integer(n.thin)
  if (n.burnin >= n.iter) stop("n.burnin must be smaller than n.iter", call. = FALSE)
  # JAGS adapts during the burn-in: the first n.burnin iterations adapt
  run <- bugs_run(prog, data, inits, parameters.to.save, n.chains, as.integer(n.adapt) + min(n.burnin, 1000L),
                  n.burnin - min(n.burnin, 1000L), n.iter - min(n.burnin, 1000L), n.thin, DIC, "jags", jags.seed, cores, refresh)
  bo <- bugs_output(run, n.chains, n.iter, n.burnin, n.thin, DIC, if (is.character(model.file)) model.file else "", "jags")
  structure(list(model = run$model, BUGSoutput = bo, parameters.to.save = c(parameters.to.save, if (DIC) "deviance"),
                 model.file = model.file, n.iter = n.iter, DIC = DIC, seed = run$seed, elapsed = run$elapsed),
            class = "rjags")
}

#' Run OpenBUGS, WinBUGS and MultiBUGS models without BUGS
#'
#' Runs a model written in the BUGS language as OpenBUGS understands it, with
#' the arguments and the result of `R2OpenBUGS::bugs()` and
#' `R2WinBUGS::bugs()`. Nothing outside R is needed.
#'
#' @inheritParams jags
#' @param inits Initial values as in [jags()]; `NULL` draws them from the
#'   priors (OpenBUGS's gen inits).
#' @param n.thin Thinning (default 1, as R2OpenBUGS).
#' @param DIC Monitor the deviance and compute DIC with OpenBUGS's
#'   `pD = mean(deviance) - deviance(posterior means)`.
#' @param seed Seed of the random numbers.
#' @param ... Other arguments of `bugs()` (`debug`, `working.directory`,
#'   `OpenBUGS.pgm`, `bugs.directory`, `codaPkg`, `useWINE`) are accepted and
#'   ignored.
#' @return An object of class `bugs`, as R2OpenBUGS returns.
#' @examples
#' code <- "model {
#'   for (i in 1:N) { y[i] ~ dnorm(mu, tau) }
#'   mu ~ dnorm(0, 1.0E-4)
#'   tau ~ dgamma(0.01, 0.01)
#' }"
#' fit <- bugs(list(y = c(1.2, 0.8, 1.9, 1.1), N = 4), list(list(mu = 0, tau = 1), list(mu = 1, tau = 2)),
#'             c("mu", "tau"), code, n.chains = 2, n.iter = 1000, seed = 1)
#' fit$summary
#' @export
bugs <- function(data, inits = NULL, parameters.to.save = NULL, model.file = "model.txt", n.chains = 3, n.iter = 2000,
                 n.burnin = floor(n.iter / 2), n.thin = 1, DIC = TRUE, seed = NULL, cores = getOption("mc.cores", 1L), refresh = 0L, ...) {
  env <- parent.frame()
  prog <- bugs_model_source(model.file)
  data <- bugs_data_arg(data, env)
  n.chains <- as.integer(n.chains); n.iter <- as.integer(n.iter); n.burnin <- as.integer(n.burnin); n.thin <- as.integer(n.thin)
  if (n.burnin >= n.iter) stop("n.burnin must be smaller than n.iter", call. = FALSE)
  na <- min(n.burnin, 1000L)
  run <- bugs_run(prog, data, inits, parameters.to.save, n.chains, na, n.burnin - na, n.iter - na, n.thin, DIC, "openbugs", seed, cores, refresh)
  if (DIC) run$Dhat <- bugs_dhat(prog, data, run)
  bugs_output(run, n.chains, n.iter, n.burnin, n.thin, DIC, if (is.character(model.file)) model.file else "", "openbugs")
}

# the deviance at the posterior means of the unobserved stochastic nodes
# (OpenBUGS's Dhat); NULL when the means are not monitored or invalid
bugs_dhat <- function(prog, data, run) {
  M <- bugs_compile(prog, data, "openbugs")
  for (v in names(run$means[[1L]])) {
    a <- get(v, envir = M$S)
    miss <- is.na(a)
    a[miss] <- Reduce(`+`, lapply(run$means, `[[`, v))[miss] / length(run$means)
    assign(v, a, envir = M$S)
  }
  for (i in M$order) if (M$nodes$kind[i] == "logical") M$stmts[[M$nodes$stmt[i]]]$fn$write(M$nodes$row[i])
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
