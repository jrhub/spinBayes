test_that("selection for default method", {
  # data(gExp)
  spbayes=BVCfit(X, Y, Z, E, clin)
  selected = BVSelection(spbayes)
  expect_true(grepl("Median Probability Model", selected$method))
  expect_equal(length(selected$indices), 3)
})

test_that("selection for non-structural method", {
  # data(gExp)
  spbayes=BVCfit(X=X, Y=Y, Z=Z, E=E, clin=clin, structural=FALSE)
  selected = BVSelection(spbayes)
  expect_true(grepl("Median Probability Model", selected$method))
  expect_true(is.null(selected$indices$Constant))
  expect_true(selected$summary["Constant effect",]==0)
})

test_that("selection for non-sparse method", {
  # data(gExp)
  spbayes=BVCfit(X=X, Y=Y, Z=Z, clin=clin, sparse=FALSE)
  selected = BVSelection(spbayes)
  expect_match(selected$method, "95% credible interval")
  expect_true(is.null(selected$indices$Linear))
  expect_true(selected$summary["Linear interaction",]==0)
})

test_that("selection for BLasso method", {
  # data(gExp)
  spbayes=BVCfit(X=X, Y=Y, Z=Z, E=E, clin=clin, VC=FALSE)
  selected = BVSelection(spbayes, prob=0.9)
  expect_match(selected$method, "90% credible interval")
})
