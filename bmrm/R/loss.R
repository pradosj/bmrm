

#' Loss functions to perform a regression
#' 
#' @name regressionLosses
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param loss.weights numeric vector of loss weights to incure for each instance of x. 
#'        Vector length should match length(y), but values are cycled if not of identical size.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
#' @examples
#'   x <- cbind(intercept=100,data.matrix(iris[1:2]))
#'   y <- iris[[3]]
#'   w <- bmrm(lmsRegressionLoss(x,y))
#'   w <- bmrm(ladRegressionLoss(x,y))
#'   w <- bmrm(quantileRegressionLoss(x,y,q=0.5))
#'   w <- bmrm(epsilonInsensitiveRegressionLoss(x,y,epsilon=1))
NULL


#' @describeIn regressionLosses Least Mean Square regression
#' @export
lmsRegressionLoss <- function(x,y,loss.weights=1/length(y)) {
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


#' @describeIn regressionLosses Least Absolute Deviation regression
#' @export
ladRegressionLoss <- function(x,y,loss.weights=1/length(y)) {
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


#' @describeIn regressionLosses Quantile Regression
#' @param q a numeric value in the range [0-1] defining quantile value to consider
#' @export
quantileRegressionLoss <- function(x,y,q=0.5,loss.weights=1/length(y)) {
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


#' @describeIn regressionLosses epsilon-insensitive regression (Vapnik et al. 1997)
#' @param epsilon a numeric value setting tolerance of the epsilon-regression
#' @export
epsilonInsensitiveRegressionLoss <- function(x,y,epsilon,loss.weights=1/length(y)) {
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


















#' Loss functions for binary classification
#' 
#' @name binaryClassificationLosses
#' @param x matrix of training instances (one instance by row)
#' @param y a 2-levels factor representing the training labels for each instance in x
#' @param loss.weights numeric vector of loss weights to incure for each instance of x. 
#'        Vector length should match length(y), but values are cycled if not of identical size.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @seealso bmrm
#' @examples
#'   x <- cbind(intercept=100,data.matrix(iris[1:2]))
#'   y <- ifelse(iris$Species=="setosa","setosa","not_setosa")
#'   w <- bmrm(hingeLoss(x,y))
#'   w <- bmrm(logisticLoss(x,y))
#'   w <- bmrm(rocLoss(x,y))
#'   w <- bmrm(fbetaLoss(x,y))
NULL


#' @describeIn binaryClassificationLosses Hinge Loss for Linear Support Vector Machine (SVM)
#' @export
hingeLoss <- function(x,y,loss.weights=1/length(y)) {
  y <- as.factor(y)
  if (nlevels(y)!=2) stop("y must have exatly 2 levels")
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  y <- c(-1,1)[y]
  
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


#' @describeIn binaryClassificationLosses logistic regression
#' @export
logisticLoss <- function(x,y,loss.weights=1/length(y)) {
  y <- as.factor(y)
  if (nlevels(y)!=2) stop("y must have exatly 2 levels")
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  y <- c(-1,1)[y]
  
  function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    loss <- loss.weights * log(1+exp(-y*f))
    grad <- loss.weights * y*(1/(1+exp(-y*f)) - 1)
    val <- colSums(loss)
    gradient(val) <- crossprod(x,grad)
    return(val)
  }
}


#' @describeIn binaryClassificationLosses Find linear weights maximize area under its ROC curve
#' @export
#' @import matrixStats
rocLoss <- function(x,y) {
  y <- as.factor(y)
  if (nlevels(y)!=2) stop("y must have exatly 2 levels")
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
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

    val <- colSums(l*c)
    gradient(val) <- crossprod(x,l)
    return(val)
  }
}



#' @describeIn binaryClassificationLosses F-beta score loss function
#' @param beta a numeric value setting the beta parameter is the f-beta score
#' @export
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
#'   # -- Build a 2D dataset from iris, and add an intercept
#'   x <- cbind(intercept=100,data.matrix(iris[c(1,2)]))
#'   y <- iris$Species
#'   
#'   # -- build the multiclass SVM model
#'   w <- bmrm(softMarginVectorLoss(x,y))
#'   dim(w) <- c(ncol(x),nlevels(y))
#'   dimnames(w) <- list(colnames(x),levels(y))
#'   F <- x %*% w
#'   pred <- colnames(F)[max.col(F)]
#'   table(pred,y)
#'   
#'   # -- Plot the dataset, the decision boundaries, the convergence curve, and the predictions
#'   gx <- seq(min(x[,2]),max(x[,2]),length=200) # positions of the probes on x-axis
#'   gy <- seq(min(x[,3]),max(x[,3]),length=200) # positions of the probes on y-axis
#'   Y <- outer(gx,gy,function(a,b){
#'      max.col(cbind(100,a,b) %*% w)
#'   })
#'   layout(matrix(c(1,3,2,3),2,2))
#'   image(gx,gy,Y,asp=1,main="dataset & decision boundaries",xlab=colnames(x)[1],ylab=colnames(x)[2])
#'   points(x[,-1],pch=19+as.integer(y))
#'   plot(attr(w,"log")$epsilon,type="o",ylab="epsilon gap",xlab="iteration")
#'   plot(row(F),F,pch=19+col(F),ylab="prediction values",xlab="sample")
softMarginVectorLoss <- function(x,y,l=1 - table(seq_along(y),y)) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.factor(y)) stop('y must be a factor')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  if (!identical(nrow(x),nrow(l))) stop('dimensions of x and l mismatch')
  if (nlevels(y)>ncol(l)) stop('some values in y are out of range of the loss matrix')
  
  function(w) {
    w <- cbind(matrix(numeric(),ncol(x)*ncol(l),0),w)

    fp <- x %*% matrix(w,ncol(x))
    fp <- matrix(fp[,t(matrix(seq_len(ncol(fp)),ncol(l)))],ncol=ncol(l)) # resize fp matrix
    fy <- fp[cbind(seq_len(nrow(fp)),y)]
    lp <- fp - fy + l[arrayInd(seq_len(nrow(fp)),nrow(x)),]
    p <- matrix(max.col(lp,ties.method='first'),nrow(x))
    lp <- matrix(lp[cbind(seq_along(p),as.vector(p))],nrow(p))
    
    # compute gradient
    gy <- gp <- array(0,c(length(y),ncol(l),ncol(w)))
    gp[cbind(as.vector(row(p)),as.vector(p),as.vector(col(p)))] <- 1
    gy[cbind(as.vector(row(p)),y,as.vector(col(p)))] <- 1
    grad <- gp - gy

    val <- colSums(lp)
    gradient(val) <- array(crossprod(x,matrix(grad,nrow(grad))),dim(w))
    return(val)
  }
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
#'   y <- iris$Species
#'   dag <- matrix(c(1,0,0,0,
#'                   0,1,1,0,
#'                   0,1,0,1),3,4,byrow=TRUE)
#'   w <- nrbm(ontologyLoss(x,y,dag=dag))
#'   w <- matrix(w,ncol(x))
#'   f <- x %*% tcrossprod(w,dag)
#'   table(y,max.col(f))
ontologyLoss <- function(x,y,l=1 - table(seq_along(y),y),dag=diag(nlevels(y))) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.factor(y)) stop('y must be a factor')
  if (!is.matrix(dag)) stop('x must be a numeric matrix')
  if (nrow(dag)!=nlevels(y)) stop('ncol(dag) should match with nlevels(y)')
  if (nrow(dag)>ncol(dag)) stop('dag matrix must have more row than column (or equal)')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  if (!identical(nrow(x),nrow(l))) stop('dimensions of x and l mismatch')
  if (nlevels(y)!=ncol(l)) stop('ncol(l) do not match with nlevels(y)')
  
  function(w) {
    w <- cbind(matrix(numeric(),ncol(x)*ncol(dag),0),w)
    fp <- x %*% matrix(w,ncol(x))
    fp <- matrix(fp[,t(matrix(seq_len(ncol(fp)),ncol(dag)))],ncol=ncol(dag)) # resize fp matrix
    z <- tcrossprod(fp,dag) + l[arrayInd(seq_len(nrow(fp)),nrow(x)),]
    Y <- matrix(max.col(z,ties.method="first"),nrow(x))
    val <- colSums(matrix(z[cbind(seq_along(Y),as.vector(Y))] - z[cbind(seq_along(Y),y[row(Y)])],nrow(x)))
    G <- dag[Y,] - dag[y[row(Y)],]
    G <- crossprod(x,matrix(G,nrow(x)))
    G <- array(G[,t(matrix(seq_len(ncol(G)),ncol(w)))],dim(w))
    gradient(val) <- G
    return(val)
  }
}

