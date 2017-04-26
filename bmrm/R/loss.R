


#' The loss function to perform a least mean square regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param loss.weights numeric vector of loss weights to incure for each instance of x. 
#'        Vector length should match length(y), but values are cycled if not of identical size.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
#' @examples
#'   x <- cbind(data.matrix(iris[1:2]),1)
#'   y <- iris[[3]]
#'   w <- bmrm(lmsRegressionLoss(x,y),LAMBDA=0.1,verbose=TRUE)
lmsRegressionLoss <- function(x,y,loss.weights=1) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  
  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    loss <- loss.weights * 0.5*(f-y)^2
    grad <- loss.weights * (f-y)
    val <- colSums(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}



#' The loss function to perform a least absolute deviation regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param loss.weights numeric vector of loss weights to incure for each instance of x. 
#'        Vector length should match length(y), but values are cycled if not of identical size.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
#' @examples
#'   x <- cbind(data.matrix(iris[1:2]),1)
#'   y <- iris[[3]]
#'   w <- bmrm(ladRegressionLoss(x,y),LAMBDA=0.1,verbose=TRUE)
ladRegressionLoss <- function(x,y,loss.weights=1) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))

  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    loss <- loss.weights * abs(f-y)
    grad <- loss.weights * sign(f-y)
    val <- colSums(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}



#' The loss function to perform a logistic regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param loss.weights numeric vector of loss weights to incure for each instance of x. 
#'        Vector length should match length(y), but values are cycled if not of identical size.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
#' @examples
#'   x <- cbind(data.matrix(iris[1:2]),1)
#'   y <- c(-1,1,1)[iris$Species]
#'   w <- bmrm(logisticRegressionLoss(x,y),LAMBDA=1,verbose=TRUE)
logisticRegressionLoss <- function(x,y,loss.weights=1) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  
  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    loss <- loss.weights * log(1+exp(-y*f))
    grad <- loss.weights * -y/(1+exp(-y*f))
    val <- colSums(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}



#' The loss function to perform a quantile regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param q a numeric value in the range [0-1] defining quantile value to consider
#' @param loss.weights numeric vector of loss weights to incure for each instance of x. 
#'        Vector length should match length(y), but values are cycled if not of identical size.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
#' @examples
#'   x <- cbind(data.matrix(iris[1:2]),1)
#'   y <- iris[[3]]
#'   w <- bmrm(quantileRegressionLoss(x,y),LAMBDA=0.1,verbose=TRUE)
quantileRegressionLoss <- function(x,y,q=0.5,loss.weights=1) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')    
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  if (length(q)!=1 || q<0 || q>1) stop('q must be a length one numeric in the range [0-1]')
  loss.weights <- rep(loss.weights,length.out=length(y))
  
  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    loss <- loss.weights * pmax(q*(f-y),(1-q)*(y-f))
    grad <- loss.weights * ifelse(f>y,q,q-1)
    val <- colSums(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}



#' The loss function to perform a epsilon-insensitive regression (Vapnik et al. 1997)
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param epsilon a numeric value setting tolerance of the epsilon-regression
#' @param loss.weights numeric vector of loss weights to incure for each instance of x. 
#'        Vector length should match length(y), but values are cycled if not of identical size.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
#' @examples
#'   x <- cbind(data.matrix(iris[1:2]),1)
#'   y <- iris[[3]]
#'   w <- bmrm(epsilonInsensitiveRegressionLoss(x,y,1),LAMBDA=0.1,verbose=TRUE)
epsilonInsensitiveRegressionLoss <- function(x,y,epsilon,loss.weights=1) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  
  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    loss <- loss.weights * pmax(abs(f-y)-epsilon,0)
    grad <- loss.weights * ifelse(abs(f-y)<epsilon,0,sign(f-y))
    val <- colSums(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}











#' Hinge Loss function for SVM
#'
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values in (-1,+1) representing the training labels for each instance in x
#' @param loss.weights numeric vector of loss weights to incure for each instance of x in case of misprediction. 
#'        Vector length should match length(y), but values are cycled if not of identical size. 
#'        Default to 1 so we define a standard 0/1 loss for SVM classifier. 
#'        The parameter might be useful to adapt SVM learning in case of unbalanced class distribution.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @seealso bmrm
#' @examples
#'   x <- cbind(data.matrix(iris[1:2]),1)
#'   y <- c(-1,1,1)[iris$Species]
#'   w <- bmrm(hingeLoss(x,y),LAMBDA=0.1,verbose=TRUE)
hingeLoss <- function(x,y,loss.weights=1) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (is.numeric(y)) {
    if (!all(y %in% c(-1,1))) stop('y must be a numeric vector of either -1 or +1, or a 2-levels factor')
  } else {
    y <- as.factor(y)
    if (nlevels(y)!=2) stop("y must have exatly 2 levels")
    y <- c(-1,1)[y]
  }
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  
  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    loss <- loss.weights * pmax(1-y*f,0)
    grad <- loss.weights * (loss>0) * (-y)
    val <- colSums(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}



#' The loss function to maximize area under the ROC curve
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y a 2-levels factor representing the training labels for each instance in x
#' @param loss.weights numeric vector of loss weights to incure for each instance of x. 
#'        Vector length should match length(y), but values are cycled if not of identical size.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @import matrixStats
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
#' @examples
#'   x <- cbind(data.matrix(iris[1:2]),1)
#'   y <- c(-1,1,1)[iris$Species]
#'   w <- bmrm(rocLoss(x,y),LAMBDA=0.1,verbose=TRUE)
rocLoss <- function(x,y,loss.weights=1) {
  y <- as.factor(y)
  if (nlevels(y)!=2) stop("y must have exatly 2 levels")
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  y <- c(-1,1)[y]
  
  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    c <- x %*% w - 0.5*y
    o <- matrix(row(c)[order(col(c),c)],nrow(c))
    
    sp <- matrixStats::colCumsums(0+matrix(y[o]==+1,nrow(o)))
    sm <- sum(y==-1) - matrixStats::colCumsums(0+matrix(y[o]==-1,nrow(o)))
    l <- 0*o
    l[cbind(as.vector(o),as.vector(col(o)))] <- ifelse(y[o]==-1,sp,-sm)
    l <- l/(sum(y==-1)*sum(y==+1))
    l <- loss.weights * l
    
    val <- colSums(l*c)
    gradient(val) <- crossprod(x,l)
    return(val)
  }
}



#' F beta score loss function
#'
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values in (-1,+1) representing the training labels for each instance in x
#' @param beta a numeric value setting the beta parameter is the f-beta score
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @seealso bmrm
#' @examples
#'   x <- cbind(data.matrix(iris[1:2]),1)
#'   y <- c(-1,1,1)[iris$Species]
#'   w <- bmrm(fbetaLoss(x,y),LAMBDA=0.01,verbose=TRUE)
fbetaLoss <- function(x,y,beta=1) {
  y <- as.factor(y)
  if (nlevels(y)!=2) stop("y must have exatly 2 levels")
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  y <- c(-1,1)[y]
  
  .fbeta <- function(TP,TN,P,N,beta) {
    beta2 <- beta*beta
    (1+beta2)*TP / (TP+N-TN+beta2*P)
  }
  
  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    
    o <- matrix(row(f)[order(col(f),-f)],nrow(f))
    op <- matrix(o[y[o]==1],ncol=ncol(o))
    on <- matrix(o[y[o]==-1],ncol=ncol(o))
    on <- on[rev(seq(nrow(on))),,drop=FALSE]
    p <- local({
      F <- matrix(f[cbind(as.vector(op),as.vector(col(op)))],nrow(op))
      2*t(colSums(F,na.rm=TRUE) - t(matrixStats::colCumsums(rbind(0,F))))
    })
    n <- local({
      F <- matrix(f[cbind(as.vector(on),as.vector(col(on)))],nrow(on))
      2*t(colSums(F,na.rm=TRUE) - t(matrixStats::colCumsums(rbind(0,F))))
    })
    
    ij <- expand.grid(i=seq(nrow(p)),j=seq(nrow(n)))
    
    # warning: matrix R might be memory consuming
    R <- 1 - .fbeta(ij$i-1,ij$j-1,nrow(op),nrow(on),beta) - p[ij$i,,drop=FALSE] + n[ij$j,,drop=FALSE]
    
    mi <- max.col(t(R),ties.method="first")
    Y <- matrix(-y,length(y),ncol(w))
    
    msk <- t(t(row(Y))<ij$i[mi])
    Y[cbind(op[msk],col(Y)[msk])] <- 1
    
    msk <- t(t(row(Y))<ij$j[mi])
    Y[cbind(on[msk],col(Y)[msk])] <- -1
    
    val <- R[cbind(mi,seq_along(mi))]
    gradient(val) <- crossprod(x,Y-y)
    return(val)
  }
}
