


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




#' Compute data for ROC curve plotting
#' 
#' @param f decision value for each instance
#' @param y a logical that specify binary labels
#' @return a 3 columns data.frame() with 'FPR' and 'TPR' at each threshold value 'f'
#' @author Julien Prados, inspired by Bob Horton code
#' @export
roc.curve <- function(f,y) {
  y <- as.logical(y)
  if (length(f)!=length(y)) stop("scores and y not of the same length")
  o <- order(f, decreasing=TRUE)
  data.frame(FPR=cumsum(!y[o])/sum(!y),TPR=cumsum(y[o])/sum(y),f=f[o])
}

