



#' @describeIn binaryClassificationLoss Hinge Loss for Linear Support Vector Machine (SVM)
#' @export
hingeLoss <- function(x,y,loss.weights=1) {
  if (!is.logical(y)) stop("y must be logical")
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  loss.weights <- loss.weights/sum(loss.weights)
  
  f <- function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    y[is.na(y)] <- f[is.na(y)]>0
    loss <- loss.weights * pmax(f*ifelse(y,-1,+1)+1,0)
    grad <- loss.weights * (loss>0) * ifelse(y,-1,+1)
    lvalue(w) <- colSums(loss)
    gradient(w) <- crossprod(x,grad)
    class(w) <- c("hingeLoss","binaryClassificationLoss")
    return(w)
  }
  is.convex(f) <- all(!is.na(y))
  return(f)
}




#' Ontology Loss Function
#' 
#' Ontology loss function may be used when the class labels are organized has an ontology structure
#' 
#' @param x instance matrix, where x(t,) defines the features of instance t
#' @param y target vector where y(t) is an integer encoding target of x(t,)
#' @param l loss matrix. l(t,p(t)) must be the loss for predicting target p(t) instead of y(t) 
#'        for instance t. By default, the parameter is set to a 0/1 loss matrix.
#' @param dag a numeric matrix defining the path in the Direct Acyclic Graph (DAG) to each class label
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @examples
#'   # -- Load the data
#'   x <- cbind(intercept=100,data.matrix(iris[1:4]))
#'   dag <- matrix(nrow=nlevels(iris$Species),byrow=TRUE,dimnames=list(levels(iris$Species)),c(
#'       1,0,0,0,
#'       0,1,1,0,
#'       0,1,0,1
#'   ))
#'   w <- nrbm(ontologyLoss(x,iris$Species,dag=dag))
#'   table(predict(w,x),iris$Species)
#'   
#'   
#'   
#'   
#'   # -- Build a 2D dataset from iris, and add an intercept
#'   x <- cbind(intercept=100,data.matrix(iris[c(1,2)]))
#'   y <- iris$Species
#'   
#'   # -- build the multiclass SVM model
#'   w <- nrbm(ontologyLoss(x,y))
#'   table(predict(w,x),y)
#'   
#'   # -- Plot the dataset, the decision boundaries, the convergence curve, and the predictions
#'   gx <- seq(min(x[,2]),max(x[,2]),length=200) # positions of the probes on x-axis
#'   gy <- seq(min(x[,3]),max(x[,3]),length=200) # positions of the probes on y-axis
#'   Y <- outer(gx,gy,function(a,b) {predict(w,cbind(100,a,b))})
#'   image(gx,gy,unclass(Y),asp=1,main="dataset & decision boundaries",
#'         xlab=colnames(x)[2],ylab=colnames(x)[3])
#'   points(x[,-1],pch=19+as.integer(y))
ontologyLoss <- function(x,y,l = 1 - table(seq_along(y),y),dag = NULL) {
  if (!is.factor(y)) stop('y must be a factor')
  if (is.null(dag)) {
    dag <- diag(nlevels(y))
    rownames(dag) <- levels(y)
  }
  if (nrow(dag) != nlevels(y)) stop('ncol(dag) should match with nlevels(y)')
  if (nrow(dag) > ncol(dag)) stop('dag matrix must have more row than column (or equal)')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  if (nrow(x) != nrow(l)) stop('dimensions of x and l mismatch')
  if (nlevels(y) != ncol(l)) stop('ncol(l) do not match with nlevels(y)')
  
  f <- function(w) {
    W <- matrix(w,ncol(x),ncol(dag),dimnames = list(colnames(x),colnames(dag)))
    fp <- tcrossprod(x %*% W,dag)
    y[is.na(y)] <- levels(y)[max.col(fp[is.na(y),],ties.method="first")]
    
    z <- fp + l
    Y <- max.col(z,ties.method = "first")
    G <- dag[Y,] - dag[y,]
    
    w <- as.vector(W)
    attr(w,"model.dim") <- dim(W)
    attr(w,"model.dimnames") <- dimnames(W)
    attr(w,"model.dag") <- dag
    lvalue(w) <- sum(z[cbind(seq_along(Y),Y)] - z[cbind(seq_along(y),y)])
    gradient(w) <- as.vector(crossprod(x,G))
    class(w) <- "ontologyLoss"
    return(w)
  }
  is.convex(f) <- all(!is.na(y))
  return(f)
}



#' @export
predict.ontologyLoss <- function(object,x,...) {
  W <- array(object,attr(object,"model.dim"),attr(object,"model.dimnames"))
  W <- tcrossprod(W,attr(object,"model.dag"))
  f <- x %*% W
  y <- max.col(f,ties.method="first")
  if (!is.null(colnames(f))) y <- factor(colnames(f)[y],colnames(f))
  attr(y,"decision.value") <- f
  y
}




