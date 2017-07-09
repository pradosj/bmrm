
#' Convex and non-convex risk minimization with L2 regularization and limited memory
#' 
#' Use algorithm of Do and Artieres, JMLR 2012 to find w minimizing: 
#' f(w) = 0.5*LAMBDA*l2norm(w) + riskFun(w)
#' where riskFun is either a convex or a non-convex risk function.
#' @param riskFun the risk function to use in the optimization (e.g.: hingeLoss, softMarginVectorLoss). 
#'   The function must evaluate the loss value and its gradient for a given point vector (w).
#' @param LAMBDA control the regularization strength in the optimization process. 
#'   This is the value used as coefficient of the regularization term.
#' @param MAX_ITER the maximum number of iteration to perform. 
#'   The function stop with a warning message if the number of iteration exceed this value
#' @param EPSILON_TOL control optimization stoping criteria: 
#'   the optimization end when the optimization gap is below this threshold
#' @param w0 initial weight vector where optimization start
#' @param maxCP mximal number of cutting plane to use to limit memory footprint
#' @param convexRisk a length 1 logical telling if the risk function riskFun is convex. 
#'    If TRUE, use CRBM algorithm; if FALSE use NRBM algorithm from Do and Artieres, JMLR 2012
#' @param LowRankQP.method a single character value defining the method used by LowRankQP (should be either "LU" or "CHOL")
#' @param regularizer 
#' @return the optimal weight vector (w)
#' @references Do and Artieres
#'   Regularized Bundle Methods for Convex and Non-Convex Risks
#'   JMLR 2012
#' @export
#' @import LowRankQP
#' @importFrom utils head
#' @examples
#'   # -- Create a 2D dataset with the first 2 features of iris, with binary labels
#'   x <- data.matrix(iris[1:2])
#'   y <- c(-1,1,1)[iris$Species]
#'
#'   # -- Add a constant dimension to the dataset to learn the intercept
#'   x <- cbind(intercept=1000,x)
#'
#'   # -- train scalar prediction models with maxMarginLoss and fbetaLoss
#'   models <- list(
#'     svm_L1 = nrbm(hingeLoss(x,y),LAMBDA=0.1,reg='l1'),
#'     svm_L2 = nrbm(hingeLoss(x,y),LAMBDA=0.1,reg='l2'),
#'     f1_L1 = nrbm(fbetaLoss(x,y),LAMBDA=0.01,reg='l1')
#'   )
#'
#'   # -- Plot the dataset and the predictions
#'   plot(x[,-1],pch=20+y,main="dataset & hyperplanes")
#'   legend('bottomright',legend=names(models),col=seq_along(models),lty=1,cex=0.75,lwd=3)
#'   for(i in seq_along(models)) {
#'     w <- models[[i]]
#'     if (w[3]!=0) abline(-w[1]*1000/w[3],-w[2]/w[3],col=i,lwd=3)
#'   }
#'
#'
#'   # -- fit a least absolute deviation linear model on a synthetic dataset
#'   # -- containing 196 meaningful features and 4 noisy features. Then
#'   # -- check if the model has detected the noise
#'   set.seed(123)
#'   X <- matrix(rnorm(4000*200), 4000, 200)
#'   beta <- c(rep(1,ncol(X)-4),0,0,0,0)
#'   Y <- X%*%beta + rnorm(nrow(X))
#'   w <- nrbm(ladRegressionLoss(X/100,Y/100),maxCP=50)
#'   barplot(as.vector(w))
nrbm <- function(riskFun,LAMBDA=1,MAX_ITER=1000L,EPSILON_TOL=0.01,w0=0,maxCP=100L,convexRisk=TRUE,LowRankQP.method="LU",regularizer=c("l2","l1")) {
  # check parameters
  if (maxCP<3) stop("maxCP should be >=3")
  regularizer <- match.arg(regularizer)
  if (regularizer=="l1" && !convexRisk) stop("l1 regularizer only support convex risk")

  # intialize first point estimation
  w <- riskFun(w0)
  at <- as.vector(gradient(w))
  bt <- as.vector(lvalue(w)) - crossprod(w,at)
  switch(regularizer,
         l1 = {f <- LAMBDA*sum(abs(w)) + lvalue(w)},
         l2 = {f <- LAMBDA*0.5*crossprod(w) + lvalue(w)}
  )
  st <- 0
  
  # initialize aggregated working plane and working set
  A <- rbind(at,at);b <- c(bt,bt);s <- c(st,st)
  inactivity.score <- c(NA_real_,NA_real_)
  ub <- f;ub.w <- w;ub.t <- 2
  for (i in 1:MAX_ITER) {
    
    switch(regularizer,
           l1 = {
             # optimize underestimator
             opt <- lp(compute.sens = 1, direction = "max",
                       objective.in = c(-1,rep_len(-LAMBDA,2L*ncol(A))),
                       const.mat = cbind(-1,A,-A), const.dir = rep("<=",nrow(A)),const.rhs = -b
             )
             if (opt$status!=0) warning("issue in the LP solver:",opt$status)
             opt$W <- matrix(opt$solution[-1L],ncol=2L)
             alpha <- opt$duals[1:nrow(A)]
             
             # compute the optimum point and corresponding objective value    
             w <- opt$W[,1L] - opt$W[,2L]
             lb <- -opt$objval
             
             # estimate loss at the new underestimator optimum
             w <- riskFun(w)
             f <- LAMBDA*sum(abs(w)) + lvalue(w)
           },
           l2 = {
             # optimize underestimator
             H <- matrix(0,1L+nrow(A),1L+nrow(A))
             H[-1,-1] <- tcrossprod(A)
             alpha <- LowRankQP(H,c(0,-LAMBDA*b),matrix(1,1L,nrow(H)),1,rep(1,nrow(H)),method=LowRankQP.method)$alpha[-1L]
             
             # compute the optimum point and corresponding objective value
             w <- as.vector(-crossprod(A,alpha) / LAMBDA)
             lb <- LAMBDA*0.5*crossprod(w) + max(0,A %*% w + b)
             
             # estimate loss at the new underestimator optimum
             w <- riskFun(w)
             f <- LAMBDA*0.5*crossprod(w) + lvalue(w)
           }
    )
    # update inactivity score
    inactivity.score <- inactivity.score + pmax(1-alpha,0)

    # deduce parameters of the new cutting plane
    at <- as.vector(gradient(w))
    bt <- lvalue(w) - crossprod(w,at)

    if (!convexRisk) {
      # L2 regularization only
      # solve possible conflicts with the new cutting plane
      if (f<ub) {
        st <- 0
        s <- s + as.vector(0.5*LAMBDA*crossprod(ub.w-w))
        b <- pmin(b,lvalue(w) - (A %*% w) - s)
      } else { # null step
        st <- 0.5*LAMBDA*crossprod(w-ub.w)
        if (lvalue(ub.w) < st + crossprod(at,ub.w) + bt) {
          U <- lvalue(ub.w) - crossprod(at,ub.w) - st
          L <- ub - crossprod(at,w) - 0.5*LAMBDA*crossprod(w)
          if (L<=U) {
            bt <- L
          } else {
            at <- -LAMBDA*ub.w
            bt <- ub - 0.5*LAMBDA*crossprod(w) - crossprod(at,w)
          }
        }
      }
    }

    # update aggregated cutting plane
    A[1,] <- alpha %*% A
    b[1] <- alpha %*% b
    
    # add the new cutting plane to the working set
    if (nrow(A)<maxCP) {
      A <- rbind(A,at);b <- c(b,bt);s <- c(s,st)
      inactivity.score <- c(inactivity.score,0)
      t <- length(b)
    } else {
      t <- which.max(inactivity.score)
      A[t,] <- at;b[t] <- bt;s[t] <- st
      inactivity.score[t] <- 0
    }
    
    # if new best found
    if (f<ub) {
      inactivity.score[ub.t] <- 0
      ub <- f;ub.w <- w;ub.t <- t
      inactivity.score[ub.t] <- NA
    }
    
    # test end of the convergence
    cat(sprintf("%d:gap=%g obj=%g reg=%g risk=%g w=[%g,%g]\n",i,ub-lb,ub,LAMBDA*0.5*crossprod(ub.w),lvalue(ub.w),min(ub.w),max(ub.w)))
    if (ub-lb < EPSILON_TOL) break
  }
  if (i >= MAX_ITER) warning('max # of itertion exceeded')
  return(ub.w)
}




