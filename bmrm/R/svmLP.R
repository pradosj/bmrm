

#' SVM Linear Program
#' 
#' Solve a linear program implementing a linear SVM 
#' with L1 regularization and L1 loss. It solves:
#' min_w LAMBDA*|w| + sum_i(e_i);
#' s.t. y_i * <w.x_i> >= 1-e_i; e_i >= 0
#' where |w| is the L1-norm of w
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
            objective.in = c(rep(LAMBDA,2L*ncol(x)),instance.weights),
            const.mat = cbind(y.num*x,-y.num*x,diag(nrow(x))),
            const.dir = ">=",
            const.rhs = 1
  )
  u <- opt$solution[seq(1L,length.out=ncol(x))] 
  v <- opt$solution[seq(ncol(x)+1L,length.out=ncol(x))]
  w <- cbind(v-u,u-v)
  colnames(w) <- levels(y)
  class(w) <- "svmLP"
  w
}


#' SVM Linear Program
#' 
#' Solve a linear program implementing multiclass-SVM 
#' with L1 regularization and L1 loss. It solves:
#' min_w LAMBDA*|w| + sum_i(e_i);
#' s.t. <w.x_i> - <w.x_j> >= 1-e_i; e_i >= 0
#' where |w| is the L1-norm of w
#' 
#' @param x a numeric data matrix.
#' @param y a response factor of >=2 levels for each row of x.
#' @param LAMBDA control the regularization strength in the optimization process. This is the value used as coefficient of the regularization term.
#' @param instance.weights a numeric vector of weight for each row of x, recycled as necessary.
#' @return the optimized weight matrix
#' @export
#' @import lpSolve
#' @author Julien Prados
#' @examples
#'   x <- cbind(100,data.matrix(iris[1:4]))
#'   y <- iris$Species
#'   w <- svmMulticlassLP(x,y)
#'   table(predict(w,x),y)
svmMulticlassLP <- function(x,y,LAMBDA=1,instance.weights=1) {
  y <- as.factor(y)
  if (nlevels(y)<2) stop("nlevels(y) must be >=2")
  if (nrow(x)!=length(y)) stop("length(y) must match nrow(x)")
  instance.weights <- rep_len(instance.weights,nrow(x))

  opt <- lp(direction = "min",
            objective.in = c(rep(LAMBDA,2L*ncol(x)*nlevels(y)),instance.weights),
            const.mat = local({
              L <- merge(data.frame(i=seq_along(y),y=y),data.frame(Y=factor(levels(y),levels(y))))
              L <- L[L$y != L$Y,]
              
              m <- matrix(rep(seq_len(nlevels(y)),each=ncol(x)),nlevels(y),ncol(x)*nlevels(y),byrow=TRUE)
              m <- (m == row(m)) + 0
              m <- m[L$y,] - m[L$Y,]
              m <- m * as.vector(x[L$i,])
              
              n <- matrix(0,nrow(L),nrow(x))
              n[cbind(seq_len(nrow(L)),L$i)] <- 1
              
              cbind(m,-m,n)
            }),
            const.dir = ">=",
            const.rhs = 1
  )
  u <- matrix(opt$solution[seq(1L,length.out=ncol(x)*nlevels(y))],ncol(x)) 
  v <- matrix(opt$solution[seq(ncol(x)*nlevels(y)+1L,length.out=ncol(x)*nlevels(y))],ncol(x))
  w <- u - v
  colnames(w) <- levels(y)
  class(w) <- "svmLP"
  w
}

    

#' Predict function for svmLP models
#'
#' @param object an object of class svmLP
#' @param x data.matrix to predict
#' @param ... unused, present to satisfy the generic predict() prototype
#' @return prediction for row of x, with an attribute "decision.value"
#' @export
#' @author Julien Prados
predict.svmLP <- function(object,x,...) {
  f <- x %*% object
  y <- colnames(f)[max.col(f,ties.method = "first")]
  attr(y,"decision.values") <- f
  return(y)
}













