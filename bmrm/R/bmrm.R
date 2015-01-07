
#' @export
gradient <- function(x,...) UseMethod("gradient")

#' @export
gradient.default <- function(x) attr(x, "gradient")

#' @export
"gradient<-" <- function(x,...) UseMethod("gradient<-")

#' @export
"gradient<-.default" <- function(x,value) {attr(x, "gradient") <- value;x}




#
# -- Define Solver for L1 regularization
# 
newL1Solver <- function(LAMBDA) {
  lp <- initProbCLP()    
  setLogLevelCLP(lp,0)
  setObjDirCLP(lp,-1)

  within(list(),{
    destroy <- function() {
      delProbCLP(lp)
    }
    addCuttingPlane <- function(a,b) {
      nc <- 2L*length(a)+1L
      if (getNumRowsCLP(lp)<1L) {
        resizeCLP(lp,0L,nc)
        chgObjCoefsCLP(lp,c(-1,rep_len(-LAMBDA,nc-1L)))
        chgColLowerCLP(lp,rep_len(0,nc))
        chgColUpperCLP(lp,rep_len(.Machine$double.xmax,nc))
      }
      addRowsCLP(lp,1L,-.Machine$double.xmax,-b,c(0L,nc),seq_len(nc)-1L,c(-1,a,-a))
    }
    regval <- function(w) {
      sum(abs(w))
    }
    optimize <- function() {
      solveInitialCLP(lp)
      if (getSolStatusCLP(lp)!=0) warning("issue in the LP solver:",status_codeCLP(getSolStatusCLP(lp))) 
      W <- getColPrimCLP(lp)
      W <- matrix(W[-1],,2L)
      w <- W[,1L]-W[,2L]
      return(list(w = w, obj = -getObjValCLP(lp)))
    }
  })
}


#
# -- Define Solver for L2 regularization
#
newL2Solver <- function(LAMBDA) {
  A <- matrix(NA_real_,0L,0L)
  b <- numeric(0L)

  within(list(),{
    destroy <- function() {}
    addCuttingPlane <- function(a,bt) {
      if (ncol(A)!=length(a)) dim(A)[2L] <- length(a)
      A <<- rbind(A,a)
      b <<- c(b,bt)
    }
    regval <- function(w) {
      0.5*crossprod(w)
    }
    optimize <- function() {
      Ale <- matrix(1,1L,nrow(A)+1L)
      H <- matrix(0,1L+nrow(A),1L+nrow(A))
      H[-1,-1] <- tcrossprod(A)
      opt <- LowRankQP(H,c(0,-LAMBDA*b),Ale,1,rep(1,nrow(A)+1L),method="LU")
      alpha <- opt$alpha[-1L]
      w <- as.vector(-crossprod(A,alpha) / LAMBDA)
      R <- max(0,A %*% w + b)
      return(list(w = w, obj = LAMBDA*regval(w)+R))
    }
  })
}



#' Bundle Methods for Regularized Risk Minimization
#' 
#' Implement Bundle Methods for Regularized Risk Minimization as described in Teo et. al 2007.
#' Find w that minimize: LAMBDA*norm(w) + lossfun(w)
#' 
#' @param lossfun the loss function to use in the optimization (e.g.: hingeLoss, softMarginVectorLoss). 
#'   The function must evaluate the loss value and its gradient for a given point vector (w).
#'   The function must be of the form lossfun(w,...,cache=NULL), i.e. accept as first parameter the vector of weight w, and unused arguments to bmrm().
#'   The return value must be a list(value,gardient,cache), where value is the numeric value of the loss for w, and gradient is the gradient vector of the function at point w.
#'   The "cache" parameter and the "cache" element in the return value can be used to store variable from one call to the next call. 
#'   The "cache" parameter is set to NULL at the first call, and is set to the previous returned "cache" value at next calls.
#' @param LAMBDA control the regularization strength in the optimization process. This is the value used as coefficient of the regularization term.
#' @param MAX_ITER the maximum number of iteration to perform. The function stop with a warning message if the number of iteration exceed this value
#' @param EPSILON_TOL control optimization stoping criteria: the optimization end when the optimization gap is below this threshold
#' @param regfun type of regularization to consider in the optimization. It can either be the character string "l1" for L1-norm regularization, 
#'   or "l2" (default) for L2-norm regularization.
#' @param verbose a length one logical. Show progression of the convergence on stdout
#' @param w0 initial weight vector where optimization start
#' @return a list of 2 fileds: "w" the optimized weight vector; "log" a data.frame showing the trace of important values in the optimization process.
#' @export
#' @import clpAPI
#' @import kernlab
#' @import methods
#' @import LowRankQP
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @author Julien Prados
#' @seealso \code{\link{hingeLoss}} \code{\link{softMarginVectorLoss}}
#' @examples
#'   # -- Create a 2D dataset with the first 2 features of iris, with binary labels
#'   x <- data.matrix(iris[1:2])
#'   y <- c(-1,1,1)[iris$Species]
#'   
#'   # -- Add a constant dimension to the dataset to learn the intercept
#'   x <- cbind(x,1)
#'   
#'   # -- train scalar prediction models with maxMarginLoss and fbetaLoss 
#'   models <- list(
#'     svm_L1 = bmrm(hingeLoss(x,y),LAMBDA=0.1,regfun='l1',verbose=TRUE),
#'     svm_L2 = bmrm(hingeLoss(x,y),LAMBDA=0.1,regfun='l2',verbose=TRUE),
#'     f1_L1 = bmrm(fbetaLoss(x,y),LAMBDA=0.01,regfun='l1',verbose=TRUE)
#'   )
#'   
#'   # -- Plot the dataset and the predictions
#'   layout(matrix(1:2,1,2))
#'   plot(x,pch=20+y,main="dataset & hyperplanes")
#'   legend('bottomright',legend=names(models),col=seq_along(models),lty=1,cex=0.75,lwd=3)
#'   for(i in seq_along(models)) {
#'     m <- models[[i]]
#'     if (m$w[2]!=0) abline(-m$w[3]/m$w[2],-m$w[1]/m$w[2],col=i,lwd=3)
#'   }
#'   
#'   rx <- range(na.rm=TRUE,1,unlist(lapply(models,function(e) nrow(e$log))))
#'   ry <- range(na.rm=TRUE,0,unlist(lapply(models,function(e) e$log$epsilon)))
#'   plot(rx,ry,type="n",ylab="epsilon gap",xlab="iteration",main="evolution of the epsilon gap")
#'   for(i in seq_along(models)) {
#'     m <- models[[i]]
#'     lines(m$log$epsilon,type="o",col=i,lwd=3)
#'   }
#'   
#'   
#'   # -- 
#'   set.seed(123)
#'   X <- matrix(rnorm(4000*200), 4000, 200)
#'   beta <- c(rep(1,ncol(X)-4),0,0,0,0)
#'   Y <- X%*%beta + rnorm(nrow(X))
#'   model <- bmrm(ladRegressionLoss(X,Y),regfun="l2",LAMBDA=100,MAX_ITER=150)
#'   layout(1)
#'   barplot(model$w)
#'   
#'   
bmrm <- function(lossfun,LAMBDA=1,MAX_ITER=100,EPSILON_TOL=0.01,regfun=c('l1','l2'),w0=0,verbose=TRUE) {
	regfun <- match.arg(regfun)
  rrm <- switch(regfun,l1=newL1Solver(LAMBDA),l2=newL2Solver(LAMBDA))
  on.exit(rrm$destroy())

	opt <- list(w=w0)
	loss <- lossfun(opt$w)
	regval <- rrm$regval(opt$w)
	ub <- LAMBDA*regval + loss
  
	log <- list(loss=numeric(),regVal=numeric(),lb=numeric(),ub=numeric(),epsilon=numeric(),nnz=integer())
	for (i in 1:MAX_ITER) {
    # reformat current gradient and weight vector
	  g <- as.vector(gradient(loss))
    w <- opt$w <- rep(opt$w,length.out=length(g))
    
    # add the new cutting plane for current weight vector and optimize
	  rrm$addCuttingPlane(g,loss - crossprod(opt$w,g))
	  opt <- rrm$optimize()
    
    # estimate loss and regularization at new optimum
	  loss <- lossfun(opt$w)
    regval <- rrm$regval(opt$w)
	  ub <- min(ub,LAMBDA*regval + loss)
	  lb <- opt$obj
    
    # log optimization status
    log$loss[i]<-loss;log$regVal[i]<-regval;log$lb[i]<-lb;log$ub[i]<-ub;log$epsilon[i]<-ub-lb;log$nnz[i]<-sum(opt$w!=0);log$commonNZ[i]<-sum(w!=0 & opt$w!=0)
	  if (verbose) {cat(sprintf("i=%d,eps=%g (=%g-%g),nnz=%d(%d),loss=%g,reg=%g\n",i,log$epsilon[i],log$ub[i],log$lb[i],log$nnz[i],log$commonNZ[i],log$loss[i],log$regVal[i]))}
    
    # test end of convergence
    if (ub-lb < EPSILON_TOL) break
	}
	if (i >= MAX_ITER) warning('max # of itertion exceeded')
	return(list(w=opt$w,log=as.data.frame(log)))
}





