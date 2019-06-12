#include<RcppArmadillo.h>
#include<Rmath.h>
// #include<stdio.h>
#include"BVCUtilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;

double rinvgaussian(double mu, double lambda){
	if(mu>1000){
		mu = 1000;
	}
	double random_sample;
    double z,y,x,u;
	z=R::rnorm(0,1);
	y=z*z;
	x=mu+0.5*mu*mu*y/lambda - 0.5*(mu/lambda)*sqrt(4*mu*lambda*y+mu*mu*y*y);
	u=R::runif(0,1);
	if(u <= mu/(mu+x)){
		random_sample = x;
	}else{
		random_sample = mu*mu/x;
	};
    return(random_sample);
}


arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma, double tol){
	unsigned int p = mu.n_elem;
	arma::vec eigval;
	arma::mat eigvec;
	// double tol = std::pow(10, -6);
	arma::eig_sym(eigval, eigvec, sigma);
	if(arma::any(eigval < (-1*tol*std::abs(eigval(eigval.n_elem -1))))){
		std::string error = std::string("covariance matrix is not positive definite");
		throw std::runtime_error(error);
	}
	arma::vec z = arma::randn(p);
	arma::vec rs = mu + eigvec * diagmat(sqrt(arma::clamp(eigval, 0, eigval.max()))) * z;
	return rs;
}

arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma){
	unsigned int p = mu.n_elem;
	arma::vec eigval;
	arma::mat eigvec;
	arma::eig_sym(eigval, eigvec, sigma);
	double scale = std::pow(10,5);
	if(arma::any((arma::round(eigval * scale)/scale) <= 0)){
		std::string error = std::string("covariance matrix is not positive definite");
		throw std::runtime_error(error);
	}
	arma::vec z = arma::randn(p);
	arma::vec rs = mu + eigvec * diagmat(sqrt(eigval)) * z;
	return rs;
}

