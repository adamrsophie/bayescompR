test_that("package namespace loads", {
  expect_true("bayescompR" %in% loadedNamespaces())
})
