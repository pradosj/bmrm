

#' Compute or check the structure of a cost matrix 
#' 
#' @param y a factor representing the labels of the instances
#' @param C either a cost matrix to check for consistency with labels in y, or a character string defining the standard matrix to compute. 
#'        If a character string the accepted values are "0/1" for a 0-1 cost matrix or "linear" for linear cost.
#' @return the cost matrix object
#' @export
#' @seealso bmrm, ordinalRegressionLoss
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
#' @seealso bmrm
#' @import matrixStats
#' @examples
#' # -- Load the data
#' x <- data.matrix(iris[1:4])
#' y <- as.integer(iris$Species)
#' 
#' # -- Train the model
#' w <- bmrm(ordinalRegressionLoss(x,y),LAMBDA=0.001,EPSILON_TOL=0.0001)
#' w2 <- bmrm(ordinalRegressionLoss(x,y,impl="quadratic"),LAMBDA=0.001,EPSILON_TOL=0.0001)
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
  
  colOrder <- function(m) {
    array(row(m)[order(col(m),m)],dim(m))
  }
  
  .loglin <- function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    c <- rbind(f-0.5,f+0.5)
    o <- colOrder(c)
    j <- array(arrayInd(o,m),dim(o))
    
    l <- array(0,c(2*m,ncol(w),length(mi)))
    l[cbind(which(o<=m,arr.ind=TRUE),y[j[o<=m]])] <- 1
    l <- array(matrixStats::colCumsums(matrix(l,nrow(l))),dim(l))

    ijk <- arrayInd(seq_along(l),dim(l))
    
    u <- array(0,c(2*m,ncol(w),length(mi)))
    u[cbind(which(o>m,arr.ind=TRUE),y[j[o>m]])] <- 1
    u <- mi[ijk[,3]] - array(matrixStats::colCumsums(matrix(u,nrow(u))),dim(u))
    
    Gu <- array(t(C)[y[j],],dim(u)) * u
    Gu[ijk[,3]>=y[j]] <- 0
    Gl <- array(C[y[j],],dim(l)) * l
    Gl[ijk[,3]<=y[j]] <- 0
    
    v <- ifelse(o<=m,-rowSums(Gu,dims=2), rowSums(Gl,dims=2))
    r <- colSums(v*c[cbind(as.vector(o),as.vector(col(o)))])
    
    g <- array(NA,c(m,ncol(w),2))
    g[cbind(as.vector(j),as.vector(col(j)),as.vector(1 + (o-1)%/%m))] <- v
    g <- rowSums(g,dims=2)
    gradient(r) <- crossprod(x,g)
    return(r)
  }
  
  .quadratic <- function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    
    # alternative computation in quadratic time for debugging purpose only
    z <- expand.grid(i=factor(1:m),j=factor(1:m))
    z <- z[y[z$i] < y[z$j],]
    z$Cij <- C[cbind(y[z$i],y[z$j])]
    R <- colSums(z$Cij * pmax(1+f[z$i,,drop=FALSE]-f[z$j,,drop=FALSE],0))
    gradient(R) <- crossprod(z$Cij*(x[z$i,]-x[z$j,]),(1+f[z$i,,drop=FALSE]-f[z$j,,drop=FALSE])>0)
    return(R)
  }
  
  switch(impl,loglin=.loglin,quadratic=.quadratic)
}




