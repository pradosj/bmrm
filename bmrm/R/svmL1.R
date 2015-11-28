
#svmL1 solve:
# min_w LAMBDA*|w| + sum(e_i)
# y_i * w.x_i >= 1-e_i
# e_i >= 0
svmL1 <- function(x,y,LAMBDA=1) {
  opt <- lp(direction = "min",
            objective.in = rep(c(LAMBDA,LAMBDA,1),c(ncol(x),ncol(x),nrow(x))),
            const.mat = cbind(y*x,-y*x,diag(nrow(x))),
            const.dir = rep(">=",nrow(x)),
            const.rhs = 1
  )
  u <- opt$solution[1:ncol(x)] 
  v <- opt$solution[ncol(x) + 1:ncol(x)]
  u - v
}






