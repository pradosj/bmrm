
#' Soft Margin Vector Loss function for multiclass SVM
#' 
#' @param x instance matrix, where x(t,) defines the features of instance t
#' @param y target vector where y(t) is an integer encoding target of x(t,)
#' @param l loss matrix. l(t,p(t)) must be the loss for predicting target p(t) instead of y(t) 
#'        for instance t. By default, the parameter is set to character value "0/1" so that the loss is set to a 0/1 loss matrix.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @examples
#'   # -- Load the data
#'   x <- data.matrix(iris[1:2])
#'   y <- iris$Species
#'   
#'   # -- Add a constant dimension to the dataset to learn the intercept
#'   cst <- sqrt(max(rowSums(x*x)))
#'   x <- cbind(x,cst)
#'   
#'   # -- train a multi-class SVM & compute the predictions
#'   train.multiclassSVM <- function(x,y,...) {
#'     w <- bmrm(softMarginVectorLoss(x,y),...)
#'     m <- list(w=as.vector(w),log=attr(w,"log"))
#'     m$w <- matrix(m$w,ncol(x))
#'     m$f <- x %*% m$w
#'     m$y <- max.col(m$f)
#'     m$contingencyTable <- table(y,m$y)
#'     print(m$contingencyTable)
#'     return(m)    
#'   }
#'   # train a L1-regularized multi-class SVM
#'   m <- train.multiclassSVM(x,y,MAX_ITER=50,regfun='l1',LAMBDA=1)
#'   
#'   # -- Plot the dataset and the decision boundaries
#'   gx <- seq(min(x[,1]),max(x[,1]),length=200) # positions of the probes on x-axis
#'   gy <- seq(min(x[,2]),max(x[,2]),length=200) # positions of the probes on y-axis
#'   Y <- outer(gx,gy,function(a,b){
#'      max.col(cbind(a,b,cst) %*% m$w)
#'   }) # matrix of predictions for all probes
#'   layout(matrix(c(1,3,2,3),2,2))
#'   image(gx,gy,Y,asp=1,main="dataset & decision boundaries",xlab=colnames(x)[1],ylab=colnames(x)[2])
#'   points(x,pch=19+as.integer(y))
#'   plot(m$log$epsilon,type="o",ylab="epsilon gap",xlab="iteration")
#'   plot(row(m$f),m$f,pch=19+col(m$f),ylab="prediction values",xlab="sample")
softMarginVectorLoss <- function(x,y,l=1 - table(seq_along(y),y)) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.factor(y)) stop('y must be a factor')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  if (!identical(nrow(x),nrow(l))) stop('dimensions of x and l mismatch')
  if (nlevels(y)>ncol(l)) stop('some values in y are out of range of the loss matrix')
  
  function(w) {
    w <- matrix(w,ncol(x),ncol(l))
    fp <- x %*% w
    fy <- rowSums(x * t(w[,y]))
    lp <- fp - fy + l
    p <- max.col(lp,ties.method='first')
    lp <- lp[cbind(1:length(p),p)]
    
    # compute gradient
    gy <- gp <- matrix(0,length(y),ncol(w))
    gp[cbind(seq_along(y),p)] <- 1
    gy[cbind(seq_along(y),y)] <- 1
    grad <- gp - gy
    
    val <- sum(lp)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}







