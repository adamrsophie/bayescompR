# Every sampling algorithm and metric must recover the mean, scales and
# correlation of a correlated Gaussian within Monte Carlo error.
test_that("all samplers target a correlated Gaussian", {
  skip_on_cran()
  mu <- c(1, -2); s <- c(1, 3); r <- 0.8
  S <- diag(s) %*% matrix(c(1, r, r, 1), 2) %*% diag(s)
  m <- bayescomp_model("data { vector[2] mu; cov_matrix[2] S; } parameters { vector[2] x; } model { x ~ multi_normal(mu, S); }",
                       data = list(mu = mu, S = S))
  cases <- list(
    list(algorithm = "nuts", metric = "diag_e"), list(algorithm = "nuts", metric = "dense_e"),
    list(algorithm = "nuts", metric = "unit_e"), list(algorithm = "hmc", metric = "diag_e"),
    list(algorithm = "mala", metric = "dense_e"), list(algorithm = "barker", metric = "dense_e"),
    list(algorithm = "metropolis", metric = "dense_e"),
    list(algorithm = "slice", metric = "diag_e"), list(algorithm = "elliptical_slice", metric = "dense_e"))
  for (cs in cases) {
    label <- paste(cs$algorithm, cs$metric)
    fit <- do.call(bayescomp_sample, c(list(m, chains = 2, iter_warmup = 1000, iter_sampling = 2000, seed = 11, refresh = 0), cs))
    a <- unclass(bayescomp_draws(fit, "x", format = "array"))
    for (k in 1:2) {
      d <- a[, , k]
      ess <- bayescomp_ess_bulk(d)
      expect_gt(ess, 100, label = paste(label, "ESS", k))
      expect_lt(abs(mean(d) - mu[k]) / (s[k] / sqrt(ess)), 5, label = paste(label, "mean", k))
      expect_lt(abs(stats::sd(d) / s[k] - 1), 0.15, label = paste(label, "sd", k))
    }
    expect_lt(abs(stats::cor(as.vector(a[, , 1]), as.vector(a[, , 2])) - r), 0.08, label = paste(label, "correlation"))
  }
})
