


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

