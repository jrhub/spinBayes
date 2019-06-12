create_BVCSparse <- function(p, iterations, burn.in){
  GS.r0 = GS.zeta = GS.phi = matrix(0, iterations, p)
  GS.r0[,1] = GS.phi[,1] = rbinom(iterations, 1, 0.9)
  GS.r0[,2] = GS.phi[,3] = rbinom(iterations, 1, 0.8)
  GS.zeta[,p] = rbinom(iterations, 1, 0.85)
  posterior = list(GS.r0=GS.r0, GS.phi=GS.phi, GS.zeta=GS.zeta)
  coefficient = NULL
  obj = list(posterior=posterior, coefficient=coefficient, iterations=iterations, burn.in=burn.in)
  class(obj) = c("BVCfit", "BVCSparse")
  obj
}

create_BVCNonSparse <- function(p, q, iterations, burn.in){
  basis = list(q=q); L=q-1
  GS.r0 = matrix(0, iterations, p)
  # GS.zeta = NULL
  GS.rs = matrix(0, iterations, p*L)
  GS.r0[,1] = GS.rs[,1] = rnorm(iterations, mean = -15, sd = 1)
  GS.r0[,2] = GS.rs[,(3*L)] = rnorm(iterations, mean = 15, sd = 1)
  GS.r0[,3] = rnorm(iterations, mean = 0, sd = 1)
  posterior = list(GS.r0=GS.r0, GS.rs=GS.rs)
  coefficient = NULL
  obj = list(posterior=posterior, coefficient=coefficient, iterations=iterations,
             burn.in=burn.in, basis=basis)
  class(obj) = c("BVCfit", "BVCNonSparse")
  obj
}

test_that("selection for BVCSparse class", {
  p = 5
  spbayes=create_BVCSparse(p, 1000, 100)
  selected = BVSelection(spbayes)
  expect_true(grepl("Median Probability Model", selected$method))
  expect_equal(selected$indices$Constant, 2)
  expect_equal(selected$indices$Varying, c(1,3))
  expect_equal(selected$indices$Linear, p)
})

test_that("selection for BVCNonSparse class", {
  spbayes=create_BVCNonSparse(5, 3, 1000, 100)
  selected = BVSelection(spbayes)
  expect_match(selected$method, "95% credible interval")
  expect_equal(selected$indices$Constant, 2)
  expect_equal(selected$indices$Varying, c(1,3))
  expect_true(is.null(selected$indices$Linear))
})

test_that("selection for default method", {
  spbayes=BVCfit(X=X, Y=Y, Z=Z, clin=clin, hyper=list(r.v=10))
  selected = BVSelection(spbayes)
  expect_true(grepl("Median Probability Model", selected$method))
  expect_equal(length(selected$indices), 3)
})

test_that("selection for non-structural method", {
  spbayes=BVCfit(X=X, Y=Y, Z=Z, E=E, clin=clin, structural=FALSE)
  selected = BVSelection(spbayes)
  expect_true(grepl("Median Probability Model", selected$method))
  expect_true(is.null(selected$indices$Constant))
  expect_true(selected$summary["Constant effect",]==0)
})

test_that("selection for non-sparse method", {
  spbayes=BVCfit(X=X, Y=Y, Z=Z, clin=clin, sparse=FALSE)
  selected = BVSelection(spbayes, prob=0.9)
  expect_match(selected$method, "90% credible interval")
  expect_output(print(selected))
  expect_true(is.null(selected$indices$Linear))
  expect_true(selected$summary["Linear interaction",]==0)
})

# test_that("selection for BLasso method", {
#   spbayes=BVCfit(X=X, Y=Y, Z=Z, E=E, clin=clin, VC=FALSE)
#   selected = BVSelection(spbayes, prob=0.9)
#   expect_output(print(selected))
#   expect_match(selected$method, "90% credible interval")
# })
