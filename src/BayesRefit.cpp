#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
// #include"BVCTests.h"
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BayesRefit (arma::mat& xx, arma::vec& y, unsigned int s0, unsigned int q, int maxSteps, double hatM, arma::vec& hatRStar, double invSigM0, arma::vec& hatInvTauSq, double hatLambdaSq, double hatSigmaSq, double aStar, double bStar, double alpha, double gamma, int progress, bool debug)
{
	unsigned int n = xx.n_rows, s = hatRStar.n_elem;
	arma::mat gsRStar(maxSteps, s),
			gsInvTauSq(maxSteps, s);
		
	arma::vec gsM(maxSteps),
			gsLambda(maxSteps),
			gsSigmaSq(maxSteps);
		
	arma::mat Br = xx.cols(1, xx.n_cols-1);
	arma::vec tBrBrDiag = sum(square(Br), 0).t();
	
	arma::vec res, muInvTauSq; // mu_m, mu_alpha, 
	double varM, meanM, tempS, meanRs, varRs;
	
	for (int k = 0; k < maxSteps; k++) {

		// m|y, r.star
		varM = 1/(n/hatSigmaSq + invSigM0);
		res = y - Br * hatRStar;
		meanM = varM * arma::accu(res/hatSigmaSq);
		hatM = R::rnorm(meanM, varM);
		res -= hatM;
		gsM(k) = hatM;
		
		for(unsigned int j=0; j<s; j++){			
			tempS = 1/(tBrBrDiag(j) + hatInvTauSq(j));
			varRs = hatSigmaSq * tempS;
			res += Br.col(j) * hatRStar(j);
			meanRs = arma::as_scalar(tempS * Br.col(j).t() * res);
			hatRStar(j) = R::rnorm(meanRs, sqrt(varRs));
			res -= Br.col(j) * hatRStar(j);
		}
		gsRStar.row(k) = hatRStar.t();
		
		// invTAUsq.star|
		muInvTauSq = arma::sqrt(hatLambdaSq * hatSigmaSq / arma::square(hatRStar));
		for(unsigned int j = 0; j<s; j++){
			hatInvTauSq(j) = rinvgaussian(muInvTauSq(j), hatLambdaSq);
		}
		gsInvTauSq.row(k) = hatInvTauSq.t();
		
		// sigma.sq|
		double shapeSig = alpha + (n+s)/2;
		double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) + 
									arma::accu(square(hatRStar) % hatInvTauSq));
		hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
		gsSigmaSq(k) = hatSigmaSq;
		
		
		// lambdaSq
		double shapeS = aStar + s;
		double rateS = bStar + arma::accu(1/hatInvTauSq)/2;
		hatLambdaSq = R::rgamma(shapeS, 1/rateS);
		gsLambda(k) = hatLambdaSq;
		
		if(progress != 0 && k % progress == 0){
			Rcpp::Rcout << "Iteration: " << k << std::endl;
			if(debug){
				Rcpp::Rcout << "hatRStar: " << hatRStar.t() << std::endl;
				Rcpp::Rcout << "hatSigmaSq: " << hatSigmaSq << std::endl;
				Rcpp::Rcout << "hatLambdaSq: " << hatLambdaSq << std::endl;
			}
		}
		
	}
	
	return Rcpp::List::create(Rcpp::Named("GS.m") = gsM,
							Rcpp::Named("GS.rs") = gsRStar,
							Rcpp::Named("GS.invTAUsq") = gsInvTauSq,
							Rcpp::Named("GS.lambda") = gsLambda,
							Rcpp::Named("GS.sigma.sq") = gsSigmaSq);
}

