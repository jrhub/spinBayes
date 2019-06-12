#include <testthat.h>
#include<RcppArmadillo.h>
#include<Rmath.h>
#include"BVCUtilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

context("mvrnormCpp") {

  test_that("negative definite covariance matrix") {
	arma::vec mean(3,fill::zeros);
	arma::mat sigmaSq(3,3,fill::eye);
	sigmaSq(0,0) = -1;
    expect_error_as(mvrnormCpp(mean, sigmaSq), std::runtime_error);
  }
  
    test_that("non-positive definite covariance matrix") {
	arma::vec mean(4,fill::zeros);
	arma::mat X(3, 4, fill::randu);
	arma::mat sigmaSq = X.t() * X;
    expect_error_as(mvrnormCpp(mean, sigmaSq), std::runtime_error);
	arma::vec out = mvrnormCpp(mean, sigmaSq, 0.1);
	expect_true(out.n_elem == 4);
  }
  
}