

# Internal method that build the constraints of the MMC assignment problem for n instances and k clusters
mmcBuildAssignmentConstraints <- function(instanceIds,k,minClusterSize,groups,minGroupOverlap,maxGroupOverlap) {
  n <- length(instanceIds)
  
  # Helper function building LP constraints for the groups
  grpConst <- function(instanceIds,k,groups,minGroupOverlap,maxGroupOverlap) {
    
    # check argument groups
    if (is.null(groups)) {
      if (!is.null(minGroupOverlap) | !is.null(maxGroupOverlap)) stop("min/maxGroupOverlap may be NULL when group is NULL")
      return(list())
    } else if (is.matrix(groups)) {
      if (!is.logical(groups)) stop("groups must be a logical matrix or its sparse representation data.frame")
      groups <- subset(melt(groups,varnames=c("sampleId","groupName")),value)
    } else {
      if (!all(c("sampleId","groupName") %in% names(groups))) stop("when a data.frame, groups must contain fields sampleId and groupName")
    }
    groups$groupName <- as.factor(groups$groupName)
    
    group.overlap.arg <- function(groupOverlap) {
      if (is.null(groupOverlap)) groupOverlap <- data.frame(clusterId=integer(),groupName=character(),value=numeric())
      if (is.matrix(groupOverlap)) {
        if (!is.numeric(groupOverlap)) stop("(min/max)GroupOverlap must be a numeric matrix or its sparse representation data.frame")
        if (is.null(colnames(groupOverlap)))
          colnames(groupOverlap) <- levels(groups$groupName)
        else if (any(colnames(groupOverlap) != levels(groups$groupName)))
          stop("when a matrix colnames(min/maxGroupOverlap) must match levels(groups$groupName)")
        groupOverlap <- subset(melt(groupOverlap,varnames=c("clusterId","groupName")),!is.na(value) & value>0)
      } else {
        if (!all(c("clusterId","groupName","value") %in% names(groupOverlap))) 
          stop("when a data.frame, (min/max)GroupOverlap must contain fields clusterId, groupName and value")
      }
      groupOverlap <- unique(groupOverlap[intersect(names(groupOverlap),c("groupName","clusterId","value"))])
      
      # Expand groups constraints to all clusters when "clusterId" is missing or if it contains NAs
      if (!exists("clusterId",groupOverlap)) groupOverlap$clusterId <- NA
      D <- data.frame(clusterId=c(rep(NA,k),seq_len(k)),extClusterId=c(seq_len(k),seq_len(k)))
      groupOverlap <- unique(merge(groupOverlap,D))
      groupOverlap$clusterId <- groupOverlap$extClusterId
      groupOverlap$extClusterId <- NULL
      groupOverlap$constraintId <- seq_len(nrow(groupOverlap))
      groupOverlap
    }
    
    # Simplify arguments
    groups <- unique(groups[intersect(names(groups),c("groupName","sampleId"))])
    minGroupOverlap <- group.overlap.arg(minGroupOverlap)
    maxGroupOverlap <- group.overlap.arg(maxGroupOverlap)
    
    # Add an identifier for the groups constraints
    maxGroupOverlap$constraintId <- maxGroupOverlap$constraintId + max(minGroupOverlap$constraintId,0)
    minGroupOverlap$dir <- rep_len(">=",nrow(minGroupOverlap))
    maxGroupOverlap$dir <- rep_len("<=",nrow(maxGroupOverlap))
    
    # Build dense constraint matrix
    n <- length(instanceIds)
    dense.const <- data.frame(sampleId=instanceIds,clusterId=rep(1:k,each=n),varId=seq_len(n*k))
    dense.const <- merge(dense.const,groups,by="sampleId")
    dense.min.const <- merge(dense.const,minGroupOverlap,by=c("groupName","clusterId"))
    dense.max.const <- merge(dense.const,maxGroupOverlap,by=c("groupName","clusterId"))
    
    list(dense.const = rbind(
      with(dense.min.const,cbind(constraintId,varId,rep_len(1,length(constraintId)))),
      with(dense.max.const,cbind(constraintId,varId,rep_len(1,length(constraintId))))
    ),
    const.dir = c(minGroupOverlap$dir,maxGroupOverlap$dir),
    const.rhs = c(minGroupOverlap$value,maxGroupOverlap$value)
    )
  }
  
  # constrain the instances to belong to one and only one cluster
  eq <- cbind(seq_len(n),seq_len(n*k),1)
  eq.dir <- rep_len("==",n)
  eq.rhs <- rep_len(1,n)
  
  # constraint on the minimum size of the clusters
  gt <- cbind(as.vector(col(matrix(NA,n,k))) + n,seq_len(n*k),1)
  gt.dir <- rep_len(">=",k)
  gt.rhs <- rep_len(minClusterSize,k)
  
  # retreive groups constraints
  grp <- grpConst(instanceIds,k,groups,minGroupOverlap,maxGroupOverlap)
  grp$dense.const[,1] <- grp$dense.const[,1] + max(gt[,1])
  
  list(
    dense.const = rbind(eq,gt,grp$dense.const),
    const.dir = c(eq.dir,gt.dir,grp$const.dir),
    const.rhs = c(eq.rhs,gt.rhs,grp$const.rhs)
  )
}

# Internal method to find cluster assignment minimizing max-margin-clustering risk
# 
# @param R numeric matrix of risks: element (i,j) is the loss penalty for assigning cluster j to sample i
# @param contraints the constraints list build using mmcBuildAssignmentConstraints
# @return a binary matrix of the same dimension as R with the solution to the assignment problem
#' @import lpSolve
mmcBestClusterAssignment <- function(R,lp.constraints) {
  opt <- lp("min",
            objective.in = as.vector(R),
            dense.const = lp.constraints$dense.const,
            const.dir = lp.constraints$const.dir,
            const.rhs = lp.constraints$const.rhs,
            binary.vec = seq_along(R)
  )
  if (opt$status!=0) stop("LP problem not feasible")
  matrix(opt$solution,nrow(R),ncol(R))
}



#' Loss function for max-margin clustering
#' 
#' @export
#' @param x numeric matrix representing the dataset (one sample per row)
#' @param k an integer specifying number of clusters to find
#' @param minClusterSize an integer vector specifying the minimum number of sample per cluster. 
#'   Given values are reclycled if necessary to have one value per cluster.
#' @param groups a logical matrix where groups[i,j] is TRUE when sample i belong to group j.
#'   Alternatively, a sparse representation of the groups matrix can be given as a data.frame. 
#'   This data.frame must contains the 2 fields "groupName" and "sampleId". 
#'   "groupName" is a character vector identifying a group of samples; 
#'   "sampleId" is an integer identifying the samples that are member of the group.
#' @param minGroupOverlap a data.frame to constrain the minimum overlap size between a group and a cluster. The data.frame may contains 
#'   the 3 fields "groupName","clusterId","overlapSize". Each row of the data.frame specify a constraint on the minimum overlap between a group and a cluster.
#'   field "clusterId" is optional and in this case the constraints apply to all clusters
#' @param weight a weight vector for each instance
#' @return the loss function to optimize for max margin clustering of the given dataset
#' @import reshape2
mmcLoss <- function(x, k=3L, minClusterSize=1L, groups=NULL, minGroupOverlap=NULL, maxGroupOverlap=NULL, weight=1/nrow(x)) {
  if (!is.matrix(x)) stop("x must be a numeric matrix")
  if (is.null(rownames(x))) rownames(x) <- seq_len(nrow(x))
  weight <- rep(weight,length.out=nrow(x))
  lp.constraints <- mmcBuildAssignmentConstraints(rownames(x),k,minClusterSize,groups,minGroupOverlap,maxGroupOverlap)
  
  function(w) {
    W <- matrix(w, ncol(x),k)
    F <- x %*% W
    R <- weight*array(rowSums(pmax(1 + (rep(1,ncol(F)) %x% F) - as.vector(F),0))-1,dim(F))
    Y <- mmcBestClusterAssignment(R,lp.constraints)
    
    G <- 1-Y+F-rowSums(F*Y)
    G <- ifelse(G>0,1,0)
    G <- Y*rowSums(G) - G
    G <- G*weight
    
    val <- sum(R*Y)
    gradient(val) <- crossprod(x,-G)
    attr(val,"R") <- R
    attr(val,"Y") <- Y
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
#' @import bmrm
#' @param x numeric matrix representing the dataset (one sample per row)
#' @param k an integer specifying number of clusters to find
#' @param N0 number of instance to randomly assign per cluster when determining a random starting point.
#'        The classification dataset it defines is used to train a multi-class SVM whose solution is used 
#'        as the starting point of current MMC iteration.
#' @param LAMBDA the complexity parameter for nrbm()
#' @param NUM_RAMDOM_START number of MMC iteration to perform with a random starting point
#' @param seed the random seed basis to use
#' @param nrbmArgsSvm arguments to nrbm() when solving for multi-class SVM problem
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
#'    y <- mmc(x,k=3,LAMBDA=0.001,minClusterSize=10,NUM_RAMDOM_START=10)
#'    table(y,iris$Species)
#'    
#'    # -- Plot the dataset and the MMC decision boundaries
#'    gx <- seq(min(x[,1]),max(x[,1]),length=200)
#'    gy <- seq(min(x[,2]),max(x[,2]),length=200)
#'    Y <- outer(gx,gy,function(a,b){predict(y,cbind(a,b,1))})
#'    image(gx,gy,Y,asp=1,main="MMC clustering",xlab=colnames(x)[1],ylab=colnames(x)[2])
#'    points(x,pch=19+max.col(y))
#'    
#'    # -- show support vectors
#'    L <- attr(y,"loss")
#'    is.sv <- rowSums(attr(L,"Y")*attr(L,"R"))>0
#'    points(x[is.sv,,drop=FALSE],col="blue",pch=8)
mmc <- function(x,k=2L,N0=3L,LAMBDA=1,NUM_RAMDOM_START=50L,seed=123,
                nrbmArgsSvm=list(maxCP=20L,MAX_ITER=100L),
                nrbmArgsMmc=list(maxCP=20L,MAX_ITER=300L),
                mc.cores=getOption("mc.cores",1L),...) {  
  nrbmArgsMmc$riskFun <- mmcLoss(x,k=k,...)
  nrbmArgsMmc$convexRisk <- FALSE
  nrbmArgsSvm$convexRisk <- TRUE
  nrbmArgsSvm$LAMBDA <- nrbmArgsMmc$LAMBDA <- k*LAMBDA
  
  mmcFit <- function(i) {
    # select a starting point w0 by randomly selecting 3 samples in each cluster and train a multi-class SVM on them
    set.seed(seed+i)
    n0 <- min(nrow(x),k*N0)
    i0 <- sample(seq_len(nrow(x)),n0)
    y0 <- as.factor(sample(rep_len(seq_len(k),n0)))
    
    nrbmArgsSvm$riskFun <- softMarginVectorLoss(x[i0,],y0)
    w0 <- do.call(nrbm,nrbmArgsSvm)
    
    # run MMC solver starting at w0
    nrbmArgsMmc$w0 <- w0
    w <- do.call(nrbm,nrbmArgsMmc)
    return(w)
  }
  mmcError <- function(w) {
    as.vector(0.5*nrbmArgsMmc$LAMBDA*crossprod(w) + nrbmArgsMmc$riskFun(w))
  }
  
  err <- mclapply(mc.cores=mc.cores,seq_len(NUM_RAMDOM_START),function(i) mmcError(mmcFit(i)))
  err <- simplify2array(err)
  
  # find the minimum of all runs
  W <- mmcFit(which.min(err))
  W <- matrix(W,ncol(x),k)
  rownames(W) <- colnames(x)
  
  L <- nrbmArgsMmc$riskFun(W)
  y <- max.col(attr(L,"Y"))
  attr(y,"W") <- W
  attr(y,"loss") <- L
  attr(y,"riskFun") <- nrbmArgsMmc$riskFun
  attr(y,"errors") <- err
  class(y) <- "mmc"
  return(y)
}



#' Predict class of new instances according to a mmc model
#' 
#' @export
#' @param object a mmc object
#' @param x a matrix similar to the dataset used, i.e. where rows are instances for
#'          which class must be predicted
#' @return a integer vector whose length match nrow(x) and containing the predicted 
#'         class for each of the given instances.
predict.mmc <- function(object,x) {
  max.col(x %*% attr(object,"W"))
}




