


#' Split a dataset for Cross Validation taking into account class balance
#' 
#' @param y the class labels of each sample of the dataset
#' @param num.cv number of cross validation required
#' @return a factor of num.cv levels that assign to each sample a test fold
#' @export
balanced.cv.fold <- function(y,num.cv=10) {
  fold <- factor(rep_len(sample(seq_len(num.cv)),length.out=length(y)))
  i <- sample(seq_along(y))
  o <- order(y[i])
  fold[i[o]] <- fold
  fold
}

#' Compute loss.weights so that total losses of each class is balanced
#' 
#' @param y a object coerced to factor that represent the class labels of each sample of the dataset
#' @return a numeric vector of the same length as y
#' @export
balanced.loss.weights <- function(y) {
  #if (is.logical(y)) return(ifelse(y,mean(!y),mean(y)))
  y <- as.factor(y)
  cw <- 1/(tabulate(y)/length(y))
  cw <- cw / sum(cw)
  cw[y]
}


#' Rank linear weight of a linear model
#' 
#' @param w a numeric vector of linear weights
#' @return a data.frame with a rank for each feature as well as z-score, p-value, and false discovery rate.
#' @export
rank.linear.weights <- function(w) {
  w <- as.vector(w)
  R <- data.frame(stringsAsFactors = FALSE,
    feature.name = if (is.null(names(w))) seq_along(w) else names(w),
    w = w,
    rk = pmin(rank(+w,ties.method="first"),rank(-w,ties.method="first"))
  )
  R$z <- as.vector((w-mean(w))/sd(w))
  R$pval <- as.vector(2*pnorm(abs(R$z),lower.tail=FALSE))
  R$fdr <- as.vector(p.adjust(R$pval,method="fdr"))
  R
}


#' Compute statistics for ROC curve plotting
#' 
#' @param f decision value for each instance
#' @param y a logical that specify binary labels
#' @return a data.frame() that compute for each threshold value 'f' roc curve statistics: TP, FP, TN, FN, FPR, TPR, sensitivity, specificity, precision, recall, accuracy
#' @author Julien Prados, adapted from Bob Horton code
#' @export
#' @examples
#'   x <- cbind(data.matrix(iris[1:4]))
#'   w <- nrbmL1(rocLoss(x,iris$Species=="versicolor"),LAMBDA=0.01)
#'   with(roc.stat(x %*% w,iris$Species=="versicolor"),plot(FPR,TPR,type="l"))
#'   with(roc.stat(-x[,2],iris$Species=="versicolor"),lines(FPR,TPR,col="blue"))
roc.stat <- function(f,y) {
  if (!is.logical(y)) stop("y must be a logical vector")
  if (length(f)!=length(y)) stop("f and y must have same length")
  o <- order(f, decreasing=TRUE)
  roc <- data.frame(
    f=c(-Inf,f[o]),
    TP = c(0,cumsum(y[o])),
    FP = c(0,cumsum(!y[o]))
  )
  roc <- roc[rev(!duplicated(rev(roc$f))),]
  roc$TN <- sum(!y) - roc$FP
  roc$FN <- sum(y) - roc$TP
  roc$FPR <- roc$FP/sum(!y)
  roc$sensitivity <- roc$recall <- roc$TPR <- roc$TP/sum(y)
  roc$accuracy <- (roc$TP+roc$TN)/length(y)
  roc$specificity <- roc$TNR <- roc$TN / (roc$TN + roc$FP)
  roc$precision <- roc$TP/(roc$TP+roc$FP)
  
  dx <- diff(roc$FPR)
  dy <- diff(roc$TPR)
  attr(roc,"AUC") <- sum(dx*roc$TPR[seq_along(dx)] + dx*dy/2)
  return(roc)
}





#' Perform multiple hierachical clustering on random subsets of a dataset
#' 
#' @param x the numeric matrix containing the data to cluster (one instance per row)
#' @param seeds a vector of random seed to use.
#' @param mc.cores number of core to use for parallelization
#' @param row.rate,col.rate numeric value in [0,1] to specify the proportion of instance 
#'        (resp. feature) to subset at each random iteration.
#' @param max.cluster upper bound on the number of expected cluster (can by +Inf).
#' @param hc.method a clustering method of arity 1, taking as input a random subset of the 
#'        input matrix x and returning an hclust object
#' @param ... additional arguments are passed to the hc.method
#' @return a list of 3 square matrices N,H,K of size nrow(x): N is the number of 
#'         time each pair of instance as been seen in the random subsets; H is the
#'         average heights where the pair of sample as been merged in the tree; K is 
#'         the average number of split possible into the trees still preserving the 
#'         two samples into the same cluster.
#' @author Julien Prados
#' @import stats
#' @export
iterative.hclust <- function(x,seeds=1:100,mc.cores=getOption("mc.cores",1L),
  row.rate=0.3,col.rate=0.1,max.cluster=10,
  hc.method=function(x,PCs=1:6,...) {hclust(dist(prcomp(x,rank.=max(PCs))$x[,PCs,drop=FALSE]),...)},
  ...
) {
  N0 <- matrix(0,nrow(x),nrow(x))
  fun <- function(n0,seed) {
    set.seed(seed)
    i <- sample(nrow(x),nrow(x)*row.rate)
    j <- sample(ncol(x),ncol(x)*col.rate)
    M <- x[i,j]
    hc <- hc.method(M,...)
    A <- outer(seq_along(hc$order),seq_along(hc$order),hclust_fca,hc=hc)
    H <- array(hc$height[A],dim(A))
    
    n0$N[i,i] <- n0$N[i,i] + 1
    n0$K[i,i] <- n0$K[i,i] + pmin(nrow(hc$merge)-A+1,max.cluster)
    n0$H[i,i] <- n0$H[i,i] + H
    n0
  }
  #N <- Reduce(fun,seeds,list(N=N0,K=N0,H=N0))
  N <- mclapply(split(seeds,seq_along(seeds)%%mc.cores),Reduce,f=fun,init=list(N=N0,K=N0,H=N0),mc.cores=mc.cores)
  N <- Reduce(function(n0,x) {n0$N <- n0$N + x$N;n0$K <- n0$K + x$K;n0$H <- n0$H + x$H;n0},N,list(N=N0,K=N0,H=N0))
  
  N$K <- ifelse(N$N>0,N$K/N$N,NA)
  N$H <- ifelse(N$N>0,N$H/N$N,NA)
  return(N)
}







