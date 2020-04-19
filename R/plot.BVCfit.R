#' plot a BVCfit object
#'
#' plot the identified varying effects
#'
#' @param x BVCfit object.
#' @param prob probability for credible interval, between 0 and 1. e.g. prob=0.95 leads to 95\% credible interval
#' @param ... other plot arguments
#' @usage \method{plot}{BVCfit}(x, prob=0.95, \dots)
#' @seealso \code{\link{BVCfit}}
#'
#' @examples
#' data(gExp)
#' spbayes=BVCfit(X, Y, Z, E, clin)
#' plot(spbayes)
#'
#'@export
plot.BVCfit=function(x, prob=0.95,...){

  if("LinOnly" %in% class(x)) stop("no varying effect is identified.")

  q = x$basis$q
  kn = x$basis$kn
  degree = x$basis$degree
  Z = x$basis$Z
  # s = x$basis$s

  adj=0.1
  s = ncol(x$coefficient$ZX)-1
  varName=colnames(x$coefficient$ZX)
  xlab = colnames(Z)
  if(is.null(xlab)) xlab="Z"
  lt = (1-prob)/2; ut= 1-lt

  u.plot = seq(min(Z), max(Z), length.out = 100)
  u.star = seq(0, 1, length=kn+2)[-c(1,kn+2)]
  Knots.star = as.numeric(stats::quantile(u.plot, u.star))
  pi.star = splines::bs(u.plot, knots=Knots.star, intercept=TRUE, degree=degree)#[,1:(q)]
  pi.star = cbind(1,pi.star[,-1]); dim(pi.star)

  selected = BVSelection(x)
  ind.V = selected$indices$Varying

  if(length(ind.V)>0){
    temp = matrix(rbind(x$posterior$GS.r0, matrix(x$posterior$GS.rs, ncol = s)), ncol = s*q)
  }

  for(j in c(0, ind.V)){
    if(j==0){
      coeff.mat = pi.star %*% t(x$posterior$GS.m)
    }else{
      last = j*q; first = last-q+1
      coeff.mat = pi.star %*% t(temp[,first:last])
    }
    pe = apply(coeff.mat, 1, stats::median)
    LL = apply(coeff.mat, 1, function(t) stats::quantile(t, lt))
    UL = apply(coeff.mat, 1, function(t) stats::quantile(t, ut))

    # lower = min(LL)-abs(min(LL))*adj; upper = max(UL)+abs(max(UL))*adj
    # graphics::plot(u.plot, pe, ylab =  bquote(~ beta[.(j)] ~ (Z)), xlab=xlab,
    #      lwd=1.5, type="l", ylim=c(lower, upper), main=paste(varName[j+1]))
    # graphics::lines(u.plot, LL, col="blue", lwd=1.5, lty = 2)
    # graphics::lines(u.plot, UL, col="blue", lwd=1.5, lty = 2)
    data <- data.frame(u.plot,pe,LL,UL)
    p <- ggplot2::ggplot(data, ggplot2::aes(x=u.plot)) +
            ggplot2::geom_line(ggplot2::aes(y = pe), color = "gray30", size=1) +
            ggplot2::geom_line(ggplot2::aes(y = LL), color="steelblue", linetype=2, size=1) +
            ggplot2::geom_line(ggplot2::aes(y = UL), color="steelblue", linetype=2, size=1) +
            ggplot2::labs(title=varName[j+1],
                   x= xlab,
                   y= bquote(~ beta[.(j)] ~ (Z)))
    print(p)
  }
}

