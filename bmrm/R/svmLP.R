

#' SVM Linear Program
#' 
#' Solve a linear program implementing a linear SVM 
#' with L1 regularization and L1 loss. It solves:
#' min_w LAMBDA*|w| + sum_i(e_i);
#' s.t. y_i * <w.x_i> >= 1-e_i; e_i >= 0
#' 
#' @param x a numeric data matrix.
#' @param y a response factor of 2 levels for each row of x.
#' @param LAMBDA control the regularization strength in the optimization process. This is the value used as coefficient of the regularization term.
#' @param instance.weights a numeric vector of weight for each row of x, recycled as necessary.
#' @return the optimized weight vector.
#' @export
#' @import lpSolve
#' @author Julien Prados
#' @examples
#'   x <- cbind(100,data.matrix(iris[1:4]))
#'   y <- iris$Species
#'   levels(y)[3] <- levels(y)[2]
#'   w <- svmLP(x,y)
#'   table(predict(w,x),y)
svmLP <- function(x,y,LAMBDA=1,instance.weights=1) {
  y <- as.factor(y)
  if (nlevels(y)!=2) stop("nlevels(y) must be 2")
  if (nrow(x)!=length(y)) stop("length(y) must match nrow(x)")
  instance.weights <- rep_len(instance.weights,nrow(x))
  y.num <- c(-1,+1)[y]
  opt <- lp(direction = "min",
            objective.in = c(rep(c(LAMBDA,LAMBDA),c(ncol(x),ncol(x))),instance.weights),
            const.mat = cbind(y.num*x,-y.num*x,diag(nrow(x))),
            const.dir = rep(">=",nrow(x)),
            const.rhs = 1
  )
  u <- opt$solution[1:ncol(x)] 
  v <- opt$solution[ncol(x) + 1:ncol(x)]
  w <- u - v
  attr(w,"levels") <- levels(y)
  class(w) <- "svmLP"
  w
}


#' Predict function for svmLP models
#'
#' @param object an object of class svmLP
#' @param x data.matrix to predict
#' @return prediction for row of x, with an attribute "decision.value"
#' @export
#' @author Julien Prados
predict.svmLP <- function(object,x,...) {
  f <- x %*% w
  y <- attr(object,"levels")[1+(f<0)]
  attr(y,"decision.value") <- f
  return(y)
}




