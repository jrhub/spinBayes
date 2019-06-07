test_that("correct length of returned coefficient list", {
  # skip_on_cran()
  data(gExp)
  spbayes=BVCfit(X, Y, Z, E, clin)
  expect_equal(length(spbayes$coefficient), 4)
  expect_equal(ncol(spbayes$coefficient$ZX), ncol(X)+1)
})
