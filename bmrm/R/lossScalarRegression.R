


#' The loss function to perform a least mean square regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
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
lmsRegressionLoss <- function(x,y) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  
  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    loss <- 0.5*(f-y)^2
    grad <- f-y
    val <- colSums(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}



#' The loss function to perform a least absolute deviation regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
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
ladRegressionLoss <- function(x,y) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')

  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    loss <- abs(f-y)
    grad <- sign(f-y)
    val <- colSums(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}



#' The loss function to perform a logistic regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
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
logisticRegressionLoss <- function(x,y) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')

  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    loss <- log(1+exp(-y*f))
    grad <- -y/(1+exp(-y*f))
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
quantileRegressionLoss <- function(x,y,q=0.5) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')    
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  if (length(q)!=1 || q<0 || q>1) stop('q must be a length one numeric in the range [0-1]')

  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    loss <- pmax(q*(f-y),(1-q)*(y-f))
    grad <- ifelse(f>y,q,q-1)
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
epsilonInsensitiveRegressionLoss <- function(x,y,epsilon) {
  
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')

  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    loss <- pmax(abs(f-y)-epsilon,0)
    grad <- ifelse(abs(f-y)<epsilon,0,sign(f-y))
    val <- colSums(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}





