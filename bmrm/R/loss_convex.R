


set.convex <- function(f) {
  is.convex(f) <- TRUE
  f
}


#' Loss functions to perform a regression
#' 
#' @name linearRegressionLoss
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param loss.weights numeric vector of loss weights to incure for each instance of x. 
#'        Vector length should match length(y), but values are cycled if not of identical size.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso nrbm
#' @examples
#'   x <- cbind(intercept=100,data.matrix(iris[1:2]))
#'   y <- iris[[3]]
#'   w <- nrbm(lmsRegressionLoss(x,y))
#'   w <- nrbm(ladRegressionLoss(x,y))
#'   w <- nrbm(quantileRegressionLoss(x,y,q=0.5))
#'   w <- nrbm(epsilonInsensitiveRegressionLoss(x,y,epsilon=1))
NULL

#' @export
predict.linearRegressionLoss <- function(object,x,...) {
  x %*% object
}

#' @describeIn linearRegressionLoss Least Mean Square regression
#' @export
lmsRegressionLoss <- function(x,y,loss.weights=1) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  
  set.convex(function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    loss <- loss.weights * 0.5*(f-y)^2
    grad <- loss.weights * (f-y)
    lvalue(w) <- sum(loss)
    gradient(w) <- crossprod(x,grad)
    class(w) <- "linearRegressionLoss"
    return(w)
  })
}


#' @describeIn linearRegressionLoss Least Absolute Deviation regression
#' @export
ladRegressionLoss <- function(x,y,loss.weights=1) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))

  set.convex(function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    loss <- loss.weights * abs(f-y)
    grad <- loss.weights * sign(f-y)
    lvalue(w) <- sum(loss)
    gradient(w) <- crossprod(x,grad)
    class(w) <- "linearRegressionLoss"
    return(w)
  })
}


#' @describeIn linearRegressionLoss Quantile Regression
#' @param q a numeric value in the range [0-1] defining quantile value to consider
#' @export
quantileRegressionLoss <- function(x,y,q=0.5,loss.weights=1) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')    
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  if (length(q)!=1 || q<0 || q>1) stop('q must be a length one numeric in the range [0-1]')
  loss.weights <- rep(loss.weights,length.out=length(y))
  
  set.convex(function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    loss <- loss.weights * pmax(q*(f-y),(1-q)*(y-f))
    grad <- loss.weights * ifelse(f>y,q,q-1)
    lvalue(w) <- sum(loss)
    gradient(w) <- crossprod(x,grad)
    class(w) <- "linearRegressionLoss"
    return(w)
  })
}


#' @describeIn linearRegressionLoss epsilon-insensitive regression (Vapnik et al. 1997)
#' @param epsilon a numeric value setting tolerance of the epsilon-regression
#' @export
epsilonInsensitiveRegressionLoss <- function(x,y,epsilon,loss.weights=1) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.numeric(y)) stop('y must be a numeric vector')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  
  set.convex(function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    loss <- loss.weights * pmax(abs(f-y)-epsilon,0)
    grad <- loss.weights * ifelse(abs(f-y)<epsilon,0,sign(f-y))
    lvalue(w) <- sum(loss)
    gradient(w) <- crossprod(x,grad)
    class(w) <- "linearRegressionLoss"
    return(w)
  })
}


















#' Loss functions for binary classification
#' 
#' @name binaryClassificationLoss
#' @param x matrix of training instances (one instance by row)
#' @param y a logical vector representing the training labels for each instance in x
#' @param loss.weights numeric vector of loss weights to incure for each instance of x. 
#'        Vector length should match length(y), but values are cycled if not of identical size.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @seealso nrbm
#' @examples
#'   x <- cbind(intercept=100,data.matrix(iris[1:2]))
#'   w <- nrbm(hingeLoss(x,iris$Species=="setosa"));predict(w,x)
#'   w <- nrbm(logisticLoss(x,iris$Species=="setosa"));predict(w,x)
#'   w <- nrbm(rocLoss(x,iris$Species=="setosa"));predict(w,x)
#'   w <- nrbm(fbetaLoss(x,iris$Species=="setosa"));predict(w,x)
NULL

#' @export
predict.binaryClassificationLoss <- function(object,x,...) {
  f <- x %*% object
  y <- f>0
  attr(y,"decision.value") <- f
  y
}

#' @describeIn binaryClassificationLoss Hinge Loss for Linear Support Vector Machine (SVM)
#' @export
hingeLoss <- function(x,y,loss.weights=1) {
  if (!is.logical(y)) stop("y must be logical")
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))

  set.convex(function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    loss <- loss.weights * pmax(ifelse(y,-f+1,f+1),0)
    grad <- loss.weights * (loss>0) * ifelse(y,-1,+1)
    lvalue(w) <- sum(loss)
    gradient(w) <- crossprod(x,grad)
    class(w) <- c("hingeLoss","binaryClassificationLoss")
    return(w)
  })
}


#' @describeIn binaryClassificationLoss logistic regression
#' @export
logisticLoss <- function(x,y,loss.weights=1) {
  if (!is.logical(y)) stop("y must be logical")
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))

  set.convex(function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    loss <- loss.weights * log(1+exp(ifelse(y,-f,+f)))
    grad <- loss.weights * ifelse(y, 1/(1+exp(-f))-1, -1/(1+exp(+f))+1)
    lvalue(w) <- sum(loss)
    gradient(w) <- crossprod(x,grad)
    class(w) <- c("logisticLoss","binaryClassificationLoss")
    return(w)
  })
}

#' @export
predict.logisticLoss <- function(object,x,...) {
  y <- NextMethod()
  f <- attr(y,"decision.value")
  attr(y,"probability") <- exp(f) / (1+exp(f))
  y
}


#' @describeIn binaryClassificationLoss Find linear weights maximize area under its ROC curve
#' @export
#' @import matrixStats
rocLoss <- function(x,y) {
  if (!is.logical(y)) stop("y must be logical")
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')

  set.convex(function(w) {
    w <- rep(w,length.out=ncol(x))
    c <- x %*% w - ifelse(y,0.5,-0.5)
    o <- order(c)
    
    sp <- cumsum(y[o])
    sm <- sum(!y) - cumsum(!y[o])
    l <- numeric(length(o))
    l[o] <- ifelse(!y[o],sp,-sm)
    l <- l/(sum(!y)*sum(y))

    lvalue(w) <- as.vector(crossprod(l,c))
    gradient(w) <- crossprod(x,l)
    class(w) <- "rocLoss"
    return(w)
  })
}

#' @export
predict.rocLoss <- function(object,x,...) {
  f <- x %*% object
  f
}

#' @describeIn binaryClassificationLoss F-beta score loss function
#' @param beta a numeric value setting the beta parameter is the f-beta score
#' @export
fbetaLoss <- function(x,y,beta=1) {
  if (!is.logical(y)) stop("y must be logical")
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')

  .fbeta <- function(TP,TN,P,N,beta) {
    beta2 <- beta*beta
    (1+beta2)*TP / (TP+N-TN+beta2*P)
  }
  
  set.convex(function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    
    o <- order(f,decreasing=TRUE)
    op <- o[y[o]]
    on <- rev(o[!y[o]])
    p <- 2*(sum(f[op]) - cumsum(c(0,f[op])))
    n <- 2*(sum(f[on]) - cumsum(c(0,f[on])))
    
    # warning: matrix R might be memory consuming
    R <- outer(seq_along(p),seq_along(n),function(i,j) {
      1 - .fbeta(i-1,j-1,length(op),length(on),beta) - p[i] + n[j]
    })      
    
    ij <- arrayInd(which.max(R),dim(R))
    Y <- -ifelse(y,1,-1)
    Y[op[seq_len(ij[1,1]-1L)]] <- 1L
    Y[on[seq_len(ij[1,2]-1L)]] <- -1L
    
    lvalue(w) <- R[ij]
    gradient(w) <- crossprod(x,Y-ifelse(y,1,-1))
    class(w) <- c("fbetaLoss","binaryClassificationLoss")
    return(w)
  })
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
#'   w <- nrbm(softMarginVectorLoss(x,y))
#'   table(predict(w,x),y)
#'   
#'   # -- Plot the dataset, the decision boundaries, the convergence curve, and the predictions
#'   gx <- seq(min(x[,2]),max(x[,2]),length=200) # positions of the probes on x-axis
#'   gy <- seq(min(x[,3]),max(x[,3]),length=200) # positions of the probes on y-axis
#'   Y <- outer(gx,gy,function(a,b) {predict(w,cbind(100,a,b))})
#'   image(gx,gy,unclass(Y),asp=1,main="dataset & decision boundaries",
#'         xlab=colnames(x)[2],ylab=colnames(x)[3])
#'   points(x[,-1],pch=19+as.integer(y))
softMarginVectorLoss <- function(x,y,l=1 - table(seq_along(y),y)) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.factor(y)) stop('y must be a factor')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  if (!identical(nrow(x),nrow(l))) stop('dimensions of x and l mismatch')
  if (any(levels(y)!=colnames(l))) stop('colnames(l) must match with levels(y)')
  
  set.convex(function(w) {
    w <- matrix(w,ncol(x),ncol(l),dimnames=list(colnames(x),levels(y)))
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

    lvalue(w) <- sum(lp)
    gradient(w) <- crossprod(x,grad)
    class(w) <- "softMarginVectorLoss"
    return(w)
  })
}

#' @export
predict.softMarginVectorLoss <- function(object,x,...) {
  f <- x %*% object
  y <- max.col(f,ties.method="first")
  y <- factor(colnames(object)[y],colnames(object))
  attr(f,"decision.value") <- f
  y
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
#'   table(predict(w,x),y)
ontologyLoss <- function(x,y,l=1 - table(seq_along(y),y),dag=diag(nlevels(y))) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.factor(y)) stop('y must be a factor')
  if (!is.matrix(dag)) stop('x must be a numeric matrix')
  if (nrow(dag)!=nlevels(y)) stop('ncol(dag) should match with nlevels(y)')
  if (nrow(dag)>ncol(dag)) stop('dag matrix must have more row than column (or equal)')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  if (!identical(nrow(x),nrow(l))) stop('dimensions of x and l mismatch')
  if (nlevels(y)!=ncol(l)) stop('ncol(l) do not match with nlevels(y)')
  
  set.convex(function(w) {
    w <- matrix(w,ncol(x),ncol(dag))
    fp <- x %*% w
    z <- tcrossprod(fp,dag) + l
    Y <- max.col(z,ties.method = "first")
    G <- dag[Y,] - dag[y,]
    lvalue(w) <- sum(z[cbind(seq_along(Y),Y)] - z[cbind(seq_along(y),y)])
    gradient(w) <- crossprod(x,G)
    class(w) <- "ontologyLoss"
    attr(w,"dag") <- dag
    return(w)
  })
}

#' @export
predict.ontologyLoss <- function(object,x,...) {
  f <- x %*% tcrossprod(object,attr(object,"dag"))
  y <- max.col(f,ties.method="first")
  attr(y,"decision.value") <- f
  y
}








#' Compute or check the structure of a cost matrix 
#' 
#' @param y a factor representing the labels of the instances
#' @param C either a cost matrix to check for consistency with labels in y, or a character string defining the standard matrix to compute. 
#'        If a character string the accepted values are "0/1" for a 0-1 cost matrix or "linear" for linear cost.
#' @return the cost matrix object
#' @export
#' @seealso nrbm, ordinalRegressionLoss
costMatrix <- function(y,C=c("0/1","linear")) {
  y <- as.factor(y)
  if (is.character(C)) {
    C <- match.arg(C)
    C <- switch(C,
                "0/1" = {
                  C <- matrix(1,nlevels(y),nlevels(y),dimnames = list(levels(y),levels(y)))
                  diag(C) <- 0
                  return(C)
                },
                "linear" = abs(outer(levels(y),levels(y),'-'))
    )
  } else {
    C <- as.matrix(C)
    if (nrow(C)!=ncol(C)) stop("C must be a square matrix")
    if (nlevels(y)!=nrow(C)) stop("dimension of the square matrix C doesn't match with number of levels in y")
  }
  return(C)
}


#' The loss function for ordinal regression
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y integer vector of positive values (>=1) representing the training labels for each instance in x
#' @param C the cost matrix to use, C[i,j] being the cost for predicting label i instead of label j.
#' @param impl either the string "loglin" or "quadratic", that define the implementation to use for the computation of the loss.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso nrbm
#' @import matrixStats
#' @examples
#' # -- Load the data
#' x <- data.matrix(iris[1:4])
#' y <- as.integer(iris$Species)
#' 
#' # -- Train the model
#' w <- nrbm(ordinalRegressionLoss(x,y),LAMBDA=0.001,EPSILON_TOL=0.0001)
#' w2 <- nrbm(ordinalRegressionLoss(x,y,impl="quadratic"),LAMBDA=0.001,EPSILON_TOL=0.0001)
#' 
#' # -- plot predictions
#' f <- x %*% w
#' f2 <- x %*% w2
#' layout(1:2)
#' plot(y,f)
#' plot(f,f2,main="compare predictions of quadratic and loglin implementations")
#' 
#' # -- Compute accuracy
#' ij <- expand.grid(i=seq(nrow(x)),j=seq(nrow(x)))
#' n <- tapply(f[ij$i] - f[ij$j]>0,list(y[ij$i],y[ij$j]),sum)
#' N <- table(y[ij$i],y[ij$j])
#' print(n/N)
ordinalRegressionLoss <- function(x,y,C="0/1",impl=c("loglin","quadratic")) {
  impl <- match.arg(impl)
  
  # check parameters at first call
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')  
  C <- costMatrix(y,C)
  y <- as.integer(y)
  m <- length(y)
  mi <- tabulate(y,nbins=ncol(C))
  M <- (m*m - sum(mi*mi))/2
  C <- C / M
  
  .loglin <- function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    c <- c(f-0.5,f+0.5)
    o <- order(c)
    j <- ((o-1)%%m)+1
    
    l <- matrix(0,2*m,length(mi))
    l[cbind(which(o<=m),y[j[o<=m]])] <- 1
    l <- matrixStats::colCumsums(l)
    
    u <- matrix(0,2*m,length(mi))
    u[cbind(which(o>m),y[j[o>m]])] <- 1
    u <- mi[col(u)] - matrixStats::colCumsums(u)
    
    Gu <- t(C)[y[j],] * u
    Gu[col(Gu)>=y[j]] <- 0
    Gl <- C[y[j],] * l
    Gl[col(Gl)<=y[j]] <- 0
    
    v <- ifelse(o<=m,-rowSums(Gu), rowSums(Gl))
    lvalue(w) <- sum(v*c[o])
    
    g <- matrix(NA,m,2)
    g[cbind(j,1 + (o-1)%/%m)] <- v
    g <- rowSums(g)
    gradient(w) <- crossprod(g,x)
    attr(w,"class") <- "ordinalRegressionLoss"
    return(w)
  }
  
  .quadratic <- function(w) {
    w <- rep(w,length.out=ncol(x))
    f <- x %*% w
    
    # alternative computation in quadratic time for debugging purpose only
    z <- expand.grid(i=factor(1:m),j=factor(1:m))
    z <- z[y[z$i] < y[z$j],]
    z$Cij <- C[cbind(y[z$i],y[z$j])]
    lvalue(w) <- sum(z$Cij * pmax(1+f[z$i,drop=FALSE]-f[z$j,drop=FALSE],0))
    gradient(w) <- crossprod(z$Cij*(x[z$i,]-x[z$j,]),(1+f[z$i]-f[z$j])>0)
    attr(w,"class") <- "ordinalRegressionLoss"
    return(w)
  }
  
  set.convex(switch(impl,loglin=.loglin,quadratic=.quadratic))
}

#' @export
predict.ordinalRegressionLoss <- function(object,x,...) {
  x %*% object
}

