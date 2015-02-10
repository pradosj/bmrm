


#' Find cluster assignment minimizing max-margin clustering risk
#' 
#' @export
#' @param R numeric matrix of risks: element (i,j) is the loss penalty for assigning cluster j to sample i
#' @param lb an integer specifying the minimum number of sample per cluster
#' @return a binary matrix of the same dimension as R with the solution to the assignment problem
mmcMinClusterAssign <- function(R,lb=1L) {
  eq <- cbind(as.vector(row(R)),seq_along(R),1)
  gt <- cbind(as.vector(col(R)) + nrow(R),seq_along(R),1)
  opt <- lp("min",
            objective.in = as.vector(R),
            dense.const = rbind(eq,gt),
            const.dir = rep(c("==",">="),c(nrow(R),ncol(R))),
            const.rhs = rep(c(1,lb),c(nrow(R),ncol(R))),
            binary.vec = seq_along(R)
  )
  if (opt$status!=0) stop("not feasible")
  matrix(opt$solution,nrow(R),ncol(R))
}


#' Compute max-margin-clustering risk associated to a given prediction matrix 
#' 
#' @export
#' @param F numeric prediction matrix of the MMC model: (i,j) being prediction of the model for sample i being part of cluster j
#' @return numeric matrix of risks where element (i,j) is the loss penalty for assigning cluster j to sample i
mmcRisk <- function(F) {
  R <- matrix(NA_real_,nrow(F),ncol(F))
  for(k in  1:ncol(F)) {
    xi <- pmax(1 + F - F[,k],0)
    xi[,k] <- NA
    R[,k] <- rowSums(xi,na.rm=TRUE)
  }
  return(R)
}


#' Loss function for max-margin clustering
#' 
#' @export
#' @param x numeric matrix representing the dataset (one sample per row)
#' @param k an integer specifying number of clusters to find
#' @param ... arguments for mmcMinClusterAssign()
#' @return the loss function to optimize for max margin clustering of the given dataset
mmcLoss <- function(x, k=3L, ...) {
  if (!is.matrix(x)) stop("x must be a numeric matrix")    
  
  function(w) {
    W <- matrix(w, ncol(x),k)
    F <- x %*% W
    R <- mmcRisk(F)
    Y <- mmcMinClusterAssign(R,...)
    
    G <- 1-Y+F-rowSums(F*Y)
    G <- ifelse(G>0,1,0)
    G <- Y*rowSums(G) - G
    
    val <- sum(R*Y)
    gradient(val) <- crossprod(x,-G)
    val
  }
}


#' Convenient wrapper function to solve max-margin clustering problem on a dataset
#' 
#' Solve max-margin clustering problem with multiple random starting points to avoid being trap by local minima.
#' The random starting points are determined by randomly assigning 3 samples to each cluster and solving for multi-class SVM
#' 
#' @export
#' @param x numeric matrix representing the dataset (one sample per row)
#' @param k an integer specifying number of clusters to find
#' @param LAMBDA the complexity parameter for nrbm()
#' @param NUM_RAMDOM_START number of random starting points to start with
#' @param seed the random seed basis to use
#' @param nrbmArgsSmv arguments to nrbm() when solving for multi-class SVM problem
#' @param nrbmArgsMmc arguments to nrbm() when solving for max-margin clustering problem
#' @param ... additional arguments are passed to mmcLoss()
#' @return the MMC model matrix.
#' @examples
#'    # -- Prepare a 2D dataset to cluster with an intercept
#'    x <- data.matrix(iris[c(1,3)])
#'    x <- x - colMeans(x)[col(x)]
#'    colSds <- sqrt(colMeans((x - colMeans(x)[col(x)])^2))
#'    x <- x / colSds[col(x)]
#'    x <- cbind(x,intercept=1)
#' 
#'    # -- Find max-margin clusters
#'    W <- mmc(x,k=3,LAMBDA=0.001,lb=10,NUM_RAMDOM_START=10)
#'    y <- max.col(x %*% W)
# 
#'    # -- Plot the dataset and the MMC decision boundaries
#'    gx <- seq(min(x[,1]),max(x[,1]),length=200)
#'    gy <- seq(min(x[,2]),max(x[,2]),length=200)
#'    Y <- outer(gx,gy,function(a,b){max.col(cbind(a,b,1) %*% W)})
#'    image(gx,gy,Y,asp=1,main="MMC clustering",xlab=colnames(x)[1],ylab=colnames(x)[2])
#'    points(x,pch=19+y)
mmc <- function(x,k=2L,LAMBDA=1,NUM_RAMDOM_START=50L,seed=123,nrbmArgsSmv=list(maxCP=50L,MAX_ITER=300L,LAMBDA=LAMBDA),nrbmArgsMmc=list(),...) {
  nrbmArgsSmv$convexRisk <- TRUE
  nrbmArgsMmc$convexRisk <- FALSE
  nrbmArgsMmc$riskFun <- mmcLoss(x,k=k,...)
  nrbmArgsMmc$LAMBDA <- LAMBDA
  
  models <- lapply(seq_len(NUM_RAMDOM_START),function(i) {
    # select a starting point w0 by randomly selecting 3 samples in each cluster and train a multi-class SVM on them
    set.seed(seed+i)
    n0 <- min(nrow(x),k*3L)
    i0 <- sample(seq_len(nrow(x)),n0)
    y0 <- sample(rep_len(seq_len(k),n0))
    nrbmArgsSmv$riskFun <- softMarginVectorLoss(x[i0,],y0)
    w0 <- do.call(nrbm,nrbmArgsSmv)
    
    # run MMC solver starting at w0
    nrbmArgsMmc$w0 <- w0
    w <- do.call(nrbm,nrbmArgsMmc)
    f <- 0.5*LAMBDA*crossprod(w) + nrbmArgsMmc$riskFun(w)
    list(f=f,W=matrix(w,ncol(x)))
  })
  # find the minimum of all runs
  f <- sapply(models,'[[',"f")
  W <- models[[which.min(f)]]$W
  return(W)
}












