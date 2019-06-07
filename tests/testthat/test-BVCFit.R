test_that("correct returned coefficient list for default method", {
  # skip_on_cran()
  spbayes=BVCfit(X, Y, Z, E, clin)
  expect_equal(length(spbayes$coefficient), 4)
  expect_equal(ncol(spbayes$coefficient$ZX), ncol(X)+1)
  expect_equal(nrow(spbayes$coefficient$ZX), spbayes$basis$q)
  expect_true(sum(spbayes$coefficient$ZX == 0)>0)
  expect_true(sum(spbayes$coefficient$EX == 0)>0)
})


test_that("correct returned coefficient list for non-sparse", {
  # skip_on_cran()
  spbayes=BVCfit(X=X, Y=Y, Z=Z, E=E, sparse=FALSE)
  expect_equal(length(spbayes$coefficient), 3)
  expect_equal(ncol(spbayes$coefficient$ZX), ncol(X)+1)
  expect_equal(nrow(spbayes$coefficient$ZX), spbayes$basis$q)
  expect_true(sum(spbayes$coefficient$ZX == 0)==0)
  expect_true(sum(spbayes$coefficient$EX == 0)==0)
})


test_that("correct returned coefficient list for non-structural", {
  # skip_on_cran()
  spbayes=BVCfit(X=X, Y=Y, Z=Z, clin=clin, structural=FALSE)
  expect_equal(length(spbayes$coefficient), 2)
  expect_equal(ncol(spbayes$coefficient$ZX), ncol(X)+1)
  expect_equal(nrow(spbayes$coefficient$ZX), spbayes$basis$q)
  expect_true(all(apply((spbayes$coefficient$ZX !=0), 2, sum)) %in% c(spbayes$basis$q, 0))
})

test_that("correct returned coefficient list for non-sparse and non-structural", {
  # skip_on_cran()
  spbayes=BVCfit(X=X, Y=Y, Z=Z, clin=clin, sparse=FALSE, structural=FALSE, kn = 3, degree = 3)
  expect_output(print(spbayes))
  expect_equal(length(spbayes$coefficient), 2)
  expect_equal(ncol(spbayes$coefficient$ZX), ncol(X)+1)
  expect_equal(nrow(spbayes$coefficient$ZX), spbayes$basis$q)
  expect_true(sum(spbayes$coefficient$ZX == 0)==0)
})

test_that("correct returned coefficient list for BLASSO", {
  # skip_on_cran()
  iterations = 5000
  spbayes=BVCfit(X, Y, Z, E, clin, VC=FALSE, iterations = iterations)
  expect_equal(length(spbayes$coefficient), 4)
  expect_equal(length(spbayes$coefficient$ZX), 4)
  expect_equal(length(spbayes$coefficient$ZX$Main), ncol(X))
  expect_equal(length(spbayes$coefficient$ZX$Interaction), ncol(X))
  expect_equal(length(spbayes$coefficient$ZX$Main), ncol(X))
  expect_equal(nrow(spbayes$posterior$GS.rs), iterations)
})

