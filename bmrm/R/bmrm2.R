


#' Bundle Methods for Regularized Risk Minimization
#' 
#' Implement Bundle Methods for Regularized Risk Minimization as described in Teo et. al 2007.
#' Find w that minimize: LAMBDA*regularization_norm(w) + lossfun(w)
#' where regularization_norm is L1 (for L2, use nrbm()).
#' 
#' @param riskFun the loss function to use in the optimization (e.g.: hingeLoss, softMarginVectorLoss). 
#'   The function must evaluate the loss value and its gradient for a given point vector (w).
#' @param LAMBDA control the regularization strength in the optimization process. This is the value used as coefficient of the regularization term.
#' @param MAX_ITER the maximum number of iteration to perform. The function stop with a warning message if the number of iteration exceed this value
#' @param EPSILON_TOL control optimization stoping criteria: the optimization end when the optimization gap is below this threshold
#' @param verbose a length one logical. Show progression of the convergence on stdout
#' @param w0 initial weight vector where optimization start
#' @return the optimized weight vector, with attribute "log" beging a data.frame storing a trace of important values of the optimization process.
#' @export
#' @import LowRankQP
#' @import lpSolve
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @author Julien Prados
#' @seealso \code{\link{hingeLoss}} \code{\link{softMarginVectorLoss}}
bmrm2 <- function(riskFun,LAMBDA=1,MAX_ITER=100,EPSILON_TOL=0.01,w0=0,verbose=TRUE) {
	regfun <- match.arg(regfun)
	
	loss <- riskFun(w0)
  g <- as.vector(gradient(loss))
  A <- matrix(numeric(0),0L,length(g))
	b <- numeric(0)
  
	w <- rep(w0,length.out=length(g))
	ub <- LAMBDA*sum(abs(w)) + loss
	log <- list(loss=numeric(),ub=numeric(),epsilon=numeric(),nnz=integer())
	for (i in 1:MAX_ITER) {
	  # add the new cutting plane to the working set
	  A <- rbind(A,g)
	  b <- c(b,loss - crossprod(w,g))

    # optimize underestimator
    opt <- lp(direction = "max",
              objective.in = c(-1,rep_len(-LAMBDA,2L*ncol(A))),
              const.mat = cbind(-1,A,-A),
              const.dir = rep("<=",nrow(A)),
              const.rhs = -b
    )
    if (opt$status!=0) warning("issue in the LP solver:",opt$status)
    opt$W <- matrix(opt$solution[-1L],ncol=2L)
    w <- opt$W[,1L] - opt$W[,2L]
    lb <- -opt$objval

	  # log optimization status
	  log$loss[i]<-loss;log$ub[i]<-ub;log$epsilon[i]<-ub-lb;log$nnz[i]<-sum(w!=0)
	  if (verbose) {cat(sprintf("%d:gap=%g obj=%g reg=%g risk=%g w=[%g,%g] nnz=%d\n",i,log$epsilon[i],log$ub[i],log$ub[i]-log$loss[i],log$loss[i],min(w),max(w),log$nnz[i]))}
	  
	  # test for the end of convergence
	  if (ub-lb < EPSILON_TOL) break
	  
    # estimate loss and regularization at new optimum
	  loss <- riskFun(w)
	  ub <- min(ub,LAMBDA*sum(abs(w)) + loss)
	  g <- as.vector(gradient(loss))    
	}
	if (i >= MAX_ITER) warning('max # of itertion exceeded')
  attr(w,"log") <- as.data.frame(log)
	return(w)
}





