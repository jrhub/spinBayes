test_that("Design.matrix", {
  kn=3; degree=2; q=kn+degree+1
  design = Design.matrix(Z, cbind(1, X), kn, degree)
  expect_equal(ncol(design$Xns), (ncol(X)+1)*q)
  expect_true(all(design$X[,(1:q)] == design$Xns[,(1:q)]))
  B = design$Xns[,-(1:q)]
  B = B[,c(TRUE, rep(FALSE, (q-1)))]
  expect_equal(ncol(B), ncol(X))
  expect_equal(B, design$X[,((1:ncol(X))+q)])
})
