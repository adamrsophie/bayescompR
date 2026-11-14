# Compare log density gradients of bayescompR against rstan for captured
# programs: Rscript validation/ecosystem/compare_rstan.R <dir> ...
# Each <dir> holds .rds files list(code, data) made by the *_programs.R scripts.
suppressMessages({ pkgload::load_all(".", quiet = TRUE); library(rstan) })
rstan_options(auto_write = FALSE)
args <- commandArgs(TRUE)
pkgs <- c("multinma", "rstanarm")
known <- list()
for (p in pkgs) if (requireNamespace(p, quietly = TRUE)) { sm <- get("stanmodels", asNamespace(p)); for (nm in names(sm)) known[[length(known) + 1L]] <- list(pkg = p, name = nm, sm = sm[[nm]]) }
norm <- function(s) gsub("[[:space:]]+", "", s)
codols <- list(CODOLS = function(X, y) as.vector(qr.coef(qr(X), y)))
for (D in args) for (f in setdiff(list.files(D, "\\.rds$", full.names = TRUE), list.files(D, "\\.sm\\.rds$", full.names = TRUE))) {
  x <- readRDS(f); if (grepl("CHAIN_ID", x$code) && is.null(x$data$CHAIN_ID)) x$data$CHAIN_ID <- 1L
  tag <- sprintf("%-10s %-14s", basename(D), sub("\\.rds$", "", basename(f)))
  r <- tryCatch({
    hit <- Filter(function(k) identical(norm(k$sm@model_code), norm(x$code)), known)
    sm <- if (length(hit)) hit[[1]]$sm else {
      cf <- sub("\\.rds$", ".sm.rds", f)
      if (file.exists(cf)) readRDS(cf) else { s <- suppressMessages(stan_model(model_code = x$code)); saveRDS(s, cf); s }
    }
    src <- if (length(hit)) paste0(hit[[1]]$pkg, "::", hit[[1]]$name) else "compiled"
    rd <- x$data
    for (nm in c("weights", "offset_", "offset")) if (is.null(rd[[nm]])) rd[[nm]] <- numeric(0)
    fit <- if (grepl("CHAIN_ID", x$code)) suppressMessages(suppressWarnings(sampling(sm, data = rd, chains = 1, iter = 1, algorithm = "Fixed_param", refresh = 0, init = 0, chain_id = 1L)))
           else suppressMessages(suppressWarnings(sampling(sm, data = rd, chains = 0)))
    m <- bayescomp_model(x$code, data = x$data, functions = if (grepl("CODOLS", x$code)) codols else list())
    if (get_num_upars(fit) != m$dim) stop(sprintf("dims differ: rstan %d, bayescompR %d", get_num_upars(fit), m$dim))
    set.seed(3)
    worst <- 0; dl <- 0; tried <- 0
    pts <- list()
    while (length(pts) < 3 && tried < 30) {
      tried <- tried + 1
      u <- runif(m$dim, -0.3, 0.3)
      a <- tryCatch(log_prob(fit, u, adjust_transform = TRUE, gradient = TRUE), error = function(e) NULL)
      if (is.null(a) || !is.finite(a)) next
      g <- bayescomp_grad_log_prob(m, u)
      worst <- max(worst, max(abs(attr(a, "gradient") - g) / (1 + abs(attr(a, "gradient")))))
      pts[[length(pts) + 1L]] <- c(rs = as.numeric(a), bc = attr(g, "log_prob"))
    }
    if (!length(pts)) stop("no finite point")
    d <- vapply(pts, function(p) p[["rs"]] - p[["bc"]], 0)
    sprintf("%s dim %3d  max rel grad err %.1e  lp offset spread %.1e  [%s]", if (worst < 1e-6 && diff(range(d)) < 1e-6) "MATCH" else "DIFF ", m$dim, worst, diff(range(d)), src)
  }, error = function(e) paste("ERR", gsub("\n", " ", substr(conditionMessage(e), 1, 200))))
  cat(tag, r, "\n")
}
