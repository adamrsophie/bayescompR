test_that("BUGS statements parse in both dialects", {
  p <- bugs_parse("
    var x[N];   # JAGS declarations
    data { for (i in 1:N) { z[i] <- sqrt(y[i]) } }
    model {
      for (i in 1 : N) {
        r[i] ~ dbin(p[i], n[i])
        logit(p[i]) <- alpha0 + alpha1 * x1[i] + b[i]
        b[i] ~ dnorm(0.0, tau)
      }
      tau ~ dgamma(1.0E-3, 1.0E-3)
      y.cens ~ dnorm(mu, 1)C(0, )
      w ~ dweib(2, 1)T(, 3)
      v ~ dnorm(0, 1)I(-1, 1)
      M[1:2, 1:2] ~ dwish(R[ , ], 3)
      s <- 1 / sqrt(tau)
    }")
  expect_s3_class(p, "bugs_program")
  expect_length(p$var, 1L)
  expect_length(p$data, 1L)
  m <- p$model
  expect_equal(m[[1L]]$kind, "for")
  expect_equal(m[[1L]]$body[[2L]]$link, "logit")
  expect_equal(m[[2L]]$args[[1L]]$value, 1e-3)
  expect_equal(m[[3L]]$bound$type, "C")
  expect_null(m[[3L]]$bound$hi)
  expect_equal(m[[4L]]$bound$type, "T")
  expect_null(m[[4L]]$bound$lo)
  expect_equal(m[[5L]]$bound$type, "I")
  expect_equal(m[[6L]]$args[[1L]]$index[[1L]]$kind, "all")
  expect_equal(m[[6L]]$lhs$index[[1L]]$kind, "range")
})

test_that("BUGS expressions follow R's precedence", {
  e <- function(s) bugs_parse(paste("model { a <-", s, "}"))$model[[1L]]$expr
  x <- e("-2^2")                       # -(2^2)
  expect_equal(x$kind, "unop"); expect_equal(x$a$op, "^")
  x <- e("1 + 2 * 3")
  expect_equal(x$op, "+"); expect_equal(x$b$op, "*")
  x <- e("2^-1")
  expect_equal(x$op, "^"); expect_equal(x$b$value, -1)
  x <- e("a %*% b + c")
  expect_equal(x$op, "+"); expect_equal(x$a$op, "%*%")
  x <- e("x > 0 && y <= 1 || !z")
  expect_equal(x$op, "||"); expect_equal(x$a$op, "&&"); expect_equal(x$b$op, "!")
})

test_that("BUGS syntax errors give the location", {
  expect_error(bugs_parse("model {\n  x ~ dnorm(0 1)\n}"), "line 2", class = "bayescomp_syntax_error")
  expect_error(bugs_parse("model { x ~ dnorm(0, 1)"), "missing '}'", class = "bayescomp_syntax_error")
  expect_error(bugs_parse("x ~ dnorm(0, 1)"), "expected 'model'", class = "bayescomp_syntax_error")
})

test_that("OpenBUGS lists fill arrays row major, R dumps column major", {
  v <- bugs_parse_values("list(N = 2, Y = structure(.Data = c(1, 2, 3, 4, 5, 6), .Dim = c(2, 3)), z = c(1, NA))")
  expect_equal(v$Y, matrix(1:6, 2, byrow = TRUE))
  expect_equal(v$z, c(1, NA))
  a <- bugs_parse_values("list(A = structure(.Data = c(1:24), .Dim = c(3, 2, 4)))")$A
  expect_equal(a[2, 1, 3], 11)
  expect_equal(a[3, 2, 4], 24)
  d <- bugs_parse_values('N <- 2\n"Y" <- structure(c(1, 2, 3, 4, 5, 6), .Dim = c(2L, 3L))')
  expect_equal(d$Y, matrix(1:6, 2))
  # .Dim split over two lines, as in the Jaws data
  expect_equal(dim(bugs_parse_values("list(R = structure(.Data = c(1, 0, 0, 1), .\n Dim = c(2, 2)))")$R), c(2L, 2L))
})

test_that("rectangular data are read by column", {
  v <- bugs_parse_values("age[] sex[] Y[, 1] Y[, 2]\n26 0 1 2\n52 1 3 NA\nEND\n")
  expect_equal(v$age, c(26, 52))
  expect_equal(v$Y, matrix(c(1, 3, 2, NA), 2))
})

test_that("data files cannot run code", {
  expect_error(bugs_parse_values("list(x = system('echo hi'))"), "not allowed", class = "bayescomp_data_error")
  expect_error(bugs_parse_values("list(x = y)"), "unknown name", class = "bayescomp_data_error")
})
