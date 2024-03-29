#' @useDynLib spinBayes, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @keywords overview
"_PACKAGE"

#' @name spinBayes-package
#' @title Semi-Parametric Gene-Environment Interaction via Bayesian Variable Selection
#' @aliases spinBayes-package
#' @description Within the Bayesian framework, we propose a partially linear varying coefficient model (PLVC) for G×E interactions.  The varying coefficient functions capture the possible non-linear G×E interaction, and the linear component  models the G×E interactions with linear assumptions. The changing of basis with B splines is adopted to separate the coefficient functions with varying, non-zero constant and zero forms, corresponding to cases of nonlinear interaction, main effect only (no interaction) and no genetic interaction at all.
#' @details The user friendly, integrated interface BVCfit() allows users to flexibly choose the fitting methods they prefer. There are three arguments in BVCfit() that control the fitting method
#' \tabular{rl}{
#' sparse: \tab whether to use the spike-and-slab priors to achieve sparsity. \cr\cr
#' VC: \tab whether to separate the coefficient functions with varying effects \cr \tab and non-zero constant (main) effects.\cr\cr
#' structural: \tab whether to use varying coefficient functions for modeling \cr \tab non-linear GxE interactions.
#' }
#' BVCfit() returns a BVCfit object that contains the posterior estimates of each coefficients.
#' S3 generic functions BVSelection(), predict(), plot() and print() are implemented for BVCfit objects.
#' BVSelection() takes a BVCfit object and returns the variable selection results.
#' predict() takes a BVCfit object and returns the predicted values for new observations.
#'
#' @references
#' Ren, J., Zhou, F., Li, X., Chen, Q., Zhang, H., Ma, S., Jiang, Y., Wu, C. (2020). Semiparametric Bayesian variable selection
#' for gene-environment interactions. \emph{Statistics in Medicine}, 39(5): 617–638. \doi{10.1002/sim.8434}.
#'
#' Zhou, F., Ren, J., Lu, X., Ma, S., and Wu, C. (2021). Gene-Environment Interaction: A Variable Selection Perspective.
#' \emph{Methods in Molecular Biology}, 2212:191-223. \doi{10.1007/978-1-0716-0947-7_13}. PMID: 33733358.
#'
#' Wu, C., Jiang, Y., Ren, J., Cui, Y., Ma, S. (2018). Dissecting gene-environment interactions: A penalized robust approach
#' accounting for hierarchical structures. \emph{Statistics in Medicine}, 37:437–456. \doi{10.1002/sim.7518}.
#'
#' Wu, C., Zhong, P.-S., and Cui, Y. (2018). Additive varying–coefficient model for nonlinear gene–environment interactions.
#' \emph{Statistical Applications in Genetics and Molecular Biology}, 17(2). \doi{10.1515/sagmb-2017-0008}.
#'
#' Jiang, Y., Huang, Y., Du, Y., Zhao, Y., Ren, J., Ma, S., Wu, C. (2017). Identification of prognostic genes and pathways in
#' lung adenocarcinoma using a Bayesian Approach. \emph{Cancer Informatics}, 1(7).
#'
#' Wu, C., Shi, X., Cui, Y., and Ma, S. (2015). A penalized robust semiparametric approach for gene-environment interactions.
#' \emph{Statistics in Medicine}, 34(30): 4016–4030. \doi{10.1002/sim.6609}.
#'
#' Wu, C., and Ma, S. (2015). A selective review of robust variable selection with applications in bioinformatics.
#' \emph{Briefings in Bioinformatics}, 16(5), 873–883. \doi{10.1093/bib/bbu046}.
#'
#' Wu, C., Cui, Y., and Ma, S. (2014). Integrative analysis of gene–environment interactions under a multi–response partially
#' linear varying coefficient model. \emph{Statistics in Medicine}, 33(28), 4988–4998. \doi{10.1002/sim.6287}.
#'
#' Wu, C., and Cui, Y. (2013). Boosting signals in gene–based association studies via efficient SNP selection.
#' \emph{Briefings in Bioinformatics}, 15(2):279–291. \doi{10.1093/bib/bbs087}.
#'
#' Wu, C., and Cui, Y. (2013). A novel method for identifying nonlinear gene–environment interactions in case–control
#' association studies. \emph{Human Genetics}, 132(12):1413–1425. \doi{10.1007/s00439-013-1350-z}.
#'
#' Wu, C., Zhong, P.S., and Cui, Y. (2013). High dimensional variable selection for gene-environment interactions.
#' \emph{Technical Report. Michigan State University}.
#'
#' Wu, C., Li, S., and Cui, Y. (2012). Genetic Association Studies: An Information Content Perspective.
#' \emph{Current Genomics}, 13(7), 566–573. \doi{10.2174/138920212803251382}.
#'
#' @seealso \code{\link{BVCfit}}
NULL

