require(bmrm)


test <- function() {
  n = 4000
  d = 200
  set.seed(123)
  X = matrix(rnorm(n*d), n, d)
  beta = c(rep(1,d-4),0,0,0,0)
  eps = rnorm(n)
  Y = X%*%beta + eps
  model <- bmrm(X,Y, lossfun=ladRegressionLoss,regfun="l2",LAMBDA=100,verbose=TRUE,MAX_ITER=3000)
  traceback()
}



