


#' Find cluster assignment minimizing max-margin clustering risk
#' 
#' @export
#' @param R numeric matrix of risks: element (i,j) is the loss penalty for assigning cluster j to sample i
#' @param minClusterSize an integer vector specifying the minimum number of sample per cluster. 
#'   Given values are reclycled if necessary to have one value per cluster.
#' @param groups a data.frame defining groups of instances. The data.frame must contains the 2 fields "groupName" and "sampleId". 
#'   "groupName" is a character vector identifying a group of samples; "sampleId" is an integer identifying the samples that are member of the group
#' @param minGroupOverlap a data.frame to constrain the minimum overlap size between a group and a cluster. The data.frame may contains 
#'   the 3 fields "groupName","clusterId","overlapSize". Each row of the data.frame specify a constraint on the minimum overlap between a group and a cluster.
#'   field "clusterId" is optional and in this case the constraints apply to all clusters
#' @return a binary matrix of the same dimension as R with the solution to the assignment problem
mmcBestClusterAssignment <- function(R,minClusterSize=1L,groups=NULL,minGroupOverlap=NULL) {
  
  # Helper function building LP constraints for the groups
  grpConst <- function(R,groups,minGroupOverlap) {
    if (is.null(minGroupOverlap)) return(list())
    # Simplify arguments
    groups <- unique(groups[intersect(names(groups),c("groupName","sampleId"))])
    minGroupOverlap <- unique(minGroupOverlap[intersect(names(minGroupOverlap),c("groupName","clusterId","overlapSize"))])
    
    # Expand groups constraints to all clusters when "clusterId" is missing or if it contains NAs
    if (!exists("clusterId",minGroupOverlap)) minGroupOverlap$clusterId <- NA
    D <- data.frame(clusterId=c(rep(NA,ncol(R)),seq_len(ncol(R))),extClusterId=c(seq_len(ncol(R)),seq_len(ncol(R))))
    minGroupOverlap <- unique(merge(minGroupOverlap,D))
    minGroupOverlap$clusterId <- minGroupOverlap$extClusterId
    minGroupOverlap$extClusterId <- NULL
  
    # Add an identifier for the groups constraints
    minGroupOverlap$constraintId <- seq_len(nrow(minGroupOverlap))
  
    # Build dense constraint matrix
    dense.const <- data.frame(sampleId=as.vector(row(R)),clusterId=as.vector(col(R)),varId=seq_along(R))
    dense.const <- merge(dense.const,groups,by="sampleId")
    dense.const <- merge(dense.const,minGroupOverlap,by=c("groupName","clusterId"))
    dense.const <- cbind(dense.const$constraintId,dense.const$varId,1)
    
    list(dense.const = dense.const,
         const.dir = rep_len(">=",nrow(minGroupOverlap)),
         const.rhs = minGroupOverlap$overlapSize
    )
  }
  
  # constrain the instances to belong to one and only one cluster
  eq <- cbind(as.vector(row(R)),seq_along(R),1)
  eq.dir <- rep_len("==",nrow(R))
  eq.rhs <- rep_len(1,nrow(R))
  
  # constraint on the minimum size of the clusters
  gt <- cbind(as.vector(col(R)) + nrow(R),seq_along(R),1)
  gt.dir <- rep_len(">=",ncol(R))
  gt.rhs <- rep_len(minClusterSize,ncol(R))
  
  # retreive groups constraints
  grp <- grpConst(R,groups,minGroupOverlap)
  grp$dense.const[,1] <- grp$dense.const[,1] + max(gt[,1])
  
  opt <- lp("min",
            objective.in = as.vector(R),
            dense.const = rbind(eq,gt,grp$dense.const),
            const.dir = c(eq.dir,gt.dir,grp$const.dir),
            const.rhs = c(eq.rhs,gt.rhs,grp$const.rhs),
            binary.vec = seq_along(R)
  )
  if (opt$status!=0) stop("LP problem not feasible")
  matrix(opt$solution,nrow(R),ncol(R))
}


#' Compute max-margin-clustering risk associated to a given prediction matrix 
#' 
#' @export
#' @param F numeric prediction matrix of the MMC model: (i,j) being prediction of the model for sample i being part of cluster j
#' @return numeric matrix of risks where element (i,j) is the loss penalty for assigning cluster j to sample i
mmcClusterAssignmentRisk <- function(F) {
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
    R <- mmcClusterAssignmentRisk(F)
    Y <- mmcBestClusterAssignment(R,...)
    
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
#' The random starting points are determined by randomly assigning N0 samples to each cluster and solving for multi-class SVM
#' 
#' @export
#' @import parallel
#' @param x numeric matrix representing the dataset (one sample per row)
#' @param k an integer specifying number of clusters to find
#' @param N0 number of instance to randomly assign per cluster when determining a random starting point.
#'        The classification dataset it defines is used to train a multi-class SVM whose solution is used 
#'        as the starting point of current MMC iteration.
#' @param LAMBDA the complexity parameter for nrbm()
#' @param NUM_RAMDOM_START number of MMC iteration to perform with a random starting point
#' @param seed the random seed basis to use
#' @param nrbmArgsSmv arguments to nrbm() when solving for multi-class SVM problem
#' @param nrbmArgsMmc arguments to nrbm() when solving for max-margin clustering problem
#' @param mc.cores number of core to use when running the random iterations in parallel
#' @param ... additional arguments are passed to mmcLoss()
#' @return the MMC model matrix
#' @examples
#'    # -- Prepare a 2D dataset to cluster with an intercept
#'    x <- data.matrix(iris[c(1,3)])
#'    x <- x - colMeans(x)[col(x)]
#'    colSds <- sqrt(colMeans((x - colMeans(x)[col(x)])^2))
#'    x <- x / colSds[col(x)]
#'    x <- cbind(x,intercept=1)
#' 
#'    # -- Find max-margin clusters
#'    W <- mmc(x,k=3,LAMBDA=0.001,minClusterSize=10,NUM_RAMDOM_START=10)
#'    y <- max.col(x %*% W)
# 
#'    # -- Plot the dataset and the MMC decision boundaries
#'    gx <- seq(min(x[,1]),max(x[,1]),length=200)
#'    gy <- seq(min(x[,2]),max(x[,2]),length=200)
#'    Y <- outer(gx,gy,function(a,b){max.col(cbind(a,b,1) %*% W)})
#'    image(gx,gy,Y,asp=1,main="MMC clustering",xlab=colnames(x)[1],ylab=colnames(x)[2])
#'    points(x,pch=19+y)
mmc <- function(x,k=2L,N0=3L,LAMBDA=1,NUM_RAMDOM_START=50L,seed=123,
                svmCall=call("nrbm",convexRisk=TRUE,maxCP=50L,MAX_ITER=300L,LAMBDA=LAMBDA),
                mmcCall=call("nrbm",convexRisk=FALSE,LAMBDA=LAMBDA),
                mc.cores=getOption("mc.cores",1L),...) {  
  mmcCall$riskFun <- mmcLoss(x,k=k,...)
  models <- mclapply(mc.cores=mc.cores,seq_len(NUM_RAMDOM_START),function(i) {
    # select a starting point w0 by randomly selecting 3 samples in each cluster and train a multi-class SVM on them
    set.seed(seed+i)
    n0 <- min(nrow(x),k*N0)
    i0 <- sample(seq_len(nrow(x)),n0)
    y0 <- sample(rep_len(seq_len(k),n0))
    
    svmCall$riskFun <- softMarginVectorLoss(x[i0,],y0)
    w0 <- eval(svmCall)
    
    # run MMC solver starting at w0
    mmcCall$w0 <- w0
    w <- eval(mmcCall)
    f <- 0.5*LAMBDA*crossprod(w) + mmcCall$riskFun(w)
    list(f=f,W=matrix(w,ncol(x)))
  })
  # find the minimum of all runs
  f <- sapply(models,'[[',"f")
  W <- models[[which.min(f)]]$W
  return(W)
}












