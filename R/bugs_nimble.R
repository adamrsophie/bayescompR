# NIMBLE interface: models written with nimbleCode() run with nimbleMCMC(),
# or nimbleModel(), configureMCMC(), buildMCMC(), compileNimble() and
# runMCMC(), with NIMBLE's arguments and results. No C++ compilation
# happens: compileNimble() returns its argument.

#' NIMBLE interface: models without NIMBLE
#'
#' Write a model with `nimbleCode()` and sample it with `nimbleMCMC()`, or
#' build it step by step with `nimbleModel()`, `configureMCMC()`,
#' `buildMCMC()`, `compileNimble()` and `runMCMC()`, as in NIMBLE. The model
#' code is NIMBLE's dialect of BUGS: named parameters such as
#' `dnorm(mu, sd = s)`, truncation with `T()` and `I()`, and vectorized
#' declarations. The samples come from the same posterior as NIMBLE gives;
#' the samplers are bayescompR's (see [jags()]), and user-defined
#' `nimbleFunction`s are not supported.
#'
#' @param code Model code: an expression in braces.
#' @param constants,data Named lists of constants and data. Missing values
#'   (`NA`) in the data are sampled.
#' @param inits Initial values: a list, a function returning one, or (in
#'   `runMCMC()` and `nimbleMCMC()`) a list of lists, one per chain. Values
#'   not given are drawn from the priors.
#' @param dimensions Accepted for compatibility (sizes come from the data
#'   and the code).
#' @param name,check,calculate,print,project,resetFunctions,showCompilerOutput,progressBar,enableWAIC
#'   Accepted for compatibility.
#' @param model A model from `nimbleModel()`.
#' @param monitors Variables to monitor; the default is the top-level
#'   stochastic variables (those with no stochastic parents), as in NIMBLE.
#' @param thin Thinning.
#' @param conf,mcmc A configuration from `configureMCMC()` or a model, and an
#'   MCMC from `buildMCMC()`.
#' @param niter,nburnin,nchains Iterations per chain (including burn-in),
#'   burn-in and chains.
#' @param setSeed `TRUE` (seed 1), a number, or `FALSE`.
#' @param samples,samplesAsCodaMCMC,summary,WAIC Return the samples, as coda
#'   `mcmc` objects, with summary statistics; `WAIC` is not computed.
#' @param cores Chains run in parallel on this many cores.
#' @param ... Unused, or the objects to "compile".
#' @return `nimbleMCMC()` and `runMCMC()` return NIMBLE's results: a matrix
#'   of samples (one chain) or a list of matrices `chain1`, `chain2`, ...,
#'   and with `summary = TRUE` a list with `samples` and `summary`.
#' @examples
#' code <- nimbleCode({
#'   for (i in 1:N) y[i] ~ dnorm(mu, sd = sigma)
#'   mu ~ dnorm(0, sd = 100)
#'   sigma ~ dunif(0, 10)
#' })
#' out <- nimbleMCMC(code, constants = list(N = 5), data = list(y = c(1.1, 0.4, 1.8, 0.9, 1.3)),
#'                   inits = list(mu = 0, sigma = 1), niter = 2000, nburnin = 1000, setSeed = 1)
#' colMeans(out)
#' @name nimble
NULL

#' @rdname nimble
#' @export
nimbleCode <- function(code) substitute(code)

#' @rdname nimble
#' @export
nimbleModel <- function(code, constants = list(), data = list(), inits = list(), dimensions = list(), name = NULL,
                        check = TRUE, calculate = TRUE, ...) {
  prog <- bugs_from_r(code)
  d <- c(bugs_nimble_values(constants), bugs_nimble_values(data))
  M <- bugs_compile(prog, d, "nimble")
  structure(list(prog = prog, data = d, inits = inits, compiled = M), class = "bayescomp_nimble_model")
}

bugs_nimble_values <- function(x) {
  if (is.null(x)) return(list())
  if (!is.list(x) || (length(x) && (is.null(names(x)) || any(!nzchar(names(x)))))) stop("constants and data must be named lists", call. = FALSE)
  x
}

# the top-level stochastic variables: no stochastic ancestors
bugs_top_level <- function(M) {
  U <- which(M$nodes$kind == "stochastic" & !M$nodes$obs)
  has_stoch_anc <- function(i) {
    seen <- integer(); q <- M$parents[[i]]
    while (length(q)) {
      j <- q[1L]; q <- q[-1L]
      if (j %in% seen) next
      seen <- c(seen, j)
      if (M$nodes$kind[j] == "stochastic") return(TRUE)
      q <- c(q, M$parents[[j]])
    }
    FALSE
  }
  top <- vapply(U, function(i) !has_stoch_anc(i), TRUE)
  vars <- unique(M$nodes$var[U])
  vars[vapply(vars, function(v) all(top[M$nodes$var[U] == v]), TRUE)]
}

#' @rdname nimble
#' @export
configureMCMC <- function(model, monitors = NULL, thin = 1, print = FALSE, enableWAIC = FALSE, ...) {
  if (!inherits(model, "bayescomp_nimble_model")) stop("model must come from nimbleModel()", call. = FALSE)
  conf <- new.env(parent = emptyenv())
  conf$model <- model
  conf$monitors <- monitors %||% bugs_top_level(model$compiled)
  conf$thin <- thin
  conf$addMonitors <- function(...) { conf$monitors <- unique(c(conf$monitors, unlist(list(...)))); invisible(conf$monitors) }
  conf$setMonitors <- function(...) { conf$monitors <- unique(unlist(list(...))); invisible(conf$monitors) }
  conf$getMonitors <- function() conf$monitors
  conf$printSamplers <- function(...) {
    st <- bugs_chain_init(model$compiled, NULL, "typical")
    for (b in st$blocks) cat(sprintf("%s sampler: %s (%d node%s)\n", sub("^conj_", "conjugate_", b$sampler), b$var, b$n, if (b$n > 1L) "s" else ""))
    invisible(NULL)
  }
  class(conf) <- "bayescomp_nimble_conf"
  conf
}

#' @rdname nimble
#' @export
buildMCMC <- function(conf, ...) {
  if (inherits(conf, "bayescomp_nimble_model")) conf <- configureMCMC(conf, ...)
  if (!inherits(conf, "bayescomp_nimble_conf")) stop("buildMCMC() takes a configuration or a model", call. = FALSE)
  structure(list(conf = conf), class = "bayescomp_nimble_mcmc")
}

#' @rdname nimble
#' @export
compileNimble <- function(..., project = NULL, resetFunctions = FALSE, showCompilerOutput = FALSE) {
  x <- list(...)
  if (length(x) == 1L) x[[1L]] else x
}

#' @rdname nimble
#' @export
runMCMC <- function(mcmc, niter = 10000, nburnin = 0, thin = NULL, nchains = 1, inits = NULL, setSeed = FALSE,
                    progressBar = FALSE, samples = TRUE, samplesAsCodaMCMC = FALSE, summary = FALSE, WAIC = FALSE,
                    cores = getOption("mc.cores", 1L), ...) {
  if (inherits(mcmc, "bayescomp_nimble_model")) mcmc <- buildMCMC(mcmc)
  if (!inherits(mcmc, "bayescomp_nimble_mcmc")) stop("mcmc must come from buildMCMC()", call. = FALSE)
  conf <- mcmc$conf; model <- conf$model
  if (isTRUE(WAIC)) warning("WAIC is not computed", call. = FALSE)
  thin <- as.integer(thin %||% conf$thin)
  niter <- as.integer(niter); nburnin <- as.integer(nburnin); nchains <- as.integer(nchains)
  if (nburnin >= niter) stop("nburnin must be smaller than niter", call. = FALSE)
  inits <- inits %||% (if (length(model$inits)) model$inits)
  seed <- if (isTRUE(setSeed)) 1L else if (is.numeric(setSeed)) as.integer(setSeed[1L]) else NULL
  na <- min(nburnin, 1000L)
  run <- bugs_run(model$prog, model$data, inits, conf$monitors, nchains, na, nburnin - na, niter - na, thin, FALSE,
                  "nimble", seed, cores, 0L)
  vn <- unlist(lapply(sort(names(run$mon)), function(r) run$mon[[r]]$names), use.names = FALSE)
  ch <- lapply(run$draws, function(d) d[, vn, drop = FALSE])
  names(ch) <- paste0("chain", seq_len(nchains))
  smp <- if (samplesAsCodaMCMC) {
    m <- lapply(ch, function(x) structure(x, mcpar = c(nburnin + thin, nburnin + thin * nrow(x), thin), class = "mcmc"))
    if (nchains == 1L) m[[1L]] else structure(m, class = "mcmc.list")
  } else if (nchains == 1L) ch[[1L]] else ch
  if (!summary) return(smp)
  sfun <- function(x) {
    out <- t(apply(x, 2L, function(v) c(Mean = mean(v), Median = stats::median(v), St.Dev. = stats::sd(v),
                                       `95%CI_low` = unname(stats::quantile(v, 0.025)), `95%CI_upp` = unname(stats::quantile(v, 0.975)))))
    out
  }
  summ <- if (nchains == 1L) sfun(ch[[1L]]) else c(lapply(ch, sfun), list(all.chains = sfun(do.call(rbind, ch))))
  out <- list(summary = summ)
  if (samples) out <- c(list(samples = smp), out)
  out
}

#' @rdname nimble
#' @export
nimbleMCMC <- function(code = NULL, constants = list(), data = list(), inits = NULL, dimensions = list(), model = NULL,
                       monitors = NULL, thin = 1, niter = 10000, nburnin = 0, nchains = 1, check = TRUE, setSeed = FALSE,
                       progressBar = FALSE, samples = TRUE, samplesAsCodaMCMC = FALSE, summary = FALSE, WAIC = FALSE,
                       cores = getOption("mc.cores", 1L)) {
  if (is.null(model)) model <- nimbleModel(code, constants, data, if (is.list(inits) && !is.null(names(inits))) inits else list())
  conf <- configureMCMC(model, monitors = monitors, thin = thin)
  runMCMC(buildMCMC(conf), niter = niter, nburnin = nburnin, thin = thin, nchains = nchains, inits = inits, setSeed = setSeed,
          samples = samples, samplesAsCodaMCMC = samplesAsCodaMCMC, summary = summary, WAIC = WAIC, cores = cores)
}
