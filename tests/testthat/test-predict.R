test_that("prediction for VarLin class", {
  spbayes=BVCfit(X, Y, Z, E, clin)
  pred = predict(spbayes, X.new=X.new, Z.new=Z.new, E.new=E.new, clin.new=clin.new, Y.new=Y.new)
  expect_equal(length(pred$y.pred), length(Y.new))
  expect_true(!is.null(pred$pmse))
  expect_error(predict(spbayes, X.new=X.new, Z.new=Z.new, clin.new=clin.new, Y.new=Y.new))
  expect_error(predict(spbayes, X.new=X.new, Z.new=Z.new, E.new=E.new, Y.new=Y.new))
  expect_error(predict(spbayes, X.new=X.new, Z.new=Z.new, E.new=E.new, clin.new=clin.new[,1], Y.new=Y.new))
})


test_that("prediction for VarOnly method", {
  spbayes=BVCfit(X=X, Y=Y, Z=Z, clin=clin, sparse=FALSE)
  pred = predict(spbayes, X.new=X.new, Z.new=Z.new, clin.new=clin.new, Y.new=Y.new)
  expect_equal(length(pred$y.pred), length(Y.new))
  expect_true(!is.null(pred$pmse))
  expect_silent(predict(spbayes, X.new=X.new, Z.new=Z.new, E.new=E.new, clin.new=clin.new))
  expect_error(predict(spbayes, X.new=X.new, Z.new=Z.new, E.new=E.new, Y.new=Y.new))
})


test_that("prediction for LinOnly method", {
  spbayes=BVCfit(X=X, Y=Y, Z=Z, E=E, clin=clin, VC=FALSE)
  pred = predict(spbayes, X.new=X.new, Z.new=Z.new, E.new=E.new, clin.new=clin.new)
  expect_output(print(pred))
  expect_equal(length(pred$y.pred), length(Y.new))
  expect_true(is.null(pred$pmse))
  expect_error(predict(spbayes, X.new=X.new, Z.new=Z.new, clin.new=clin.new, Y.new=Y.new))
  expect_error(predict(spbayes, X.new=X.new, Z.new=Z.new, E.new=E.new, Y.new=Y.new))
})
