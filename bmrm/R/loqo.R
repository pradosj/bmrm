
#' A copy of the method ipop from kernlab package
#' 
#' A copy of the method ipop from kernlab package that sove QP using LOQO method
#' 
#' @param c 
#' @param H
#' @param A
#' @param b
#' @param l
#' @param u
#' @param r
#' @param sigf
#' @param maxiter
#' @param margin
#' @param bound
#' @param
#' @return 
#' @export
#' @seealso \code{\link{ipop}}
loqo <- function (c, H, A, b, l, u, r, sigf = 7, maxiter = 40, margin = 0.05, bound = 10, verb = 0) {
  if (!is.matrix(H)) stop("H must be a matrix")
  if (!is.matrix(A) && !is.vector(A)) stop("A must be a matrix or a vector")
  if (!is.matrix(c) && !is.vector(c)) stop("c must be a matrix or a vector")
  if (!is.matrix(l) && !is.vector(l)) stop("l must be a matrix or a vector")
  if (!is.matrix(u) && !is.vector(u)) stop("u must be a matrix or a vector")
  if (is.vector(A)) A <- matrix(A, 1)
  if (ncol(H) > nrow(H)) H <- t(H)
  is.square <- (ncol(H) == nrow(H))
  primal <- rep(0, nrow(H))
  if (missing(b)) bvec <- rep(0, nrow(A))
  if (nrow(H) != length(c)) stop("H and c are incompatible!")
  if (nrow(H) != ncol(A)) stop("A and c are incompatible!")
  if (nrow(A) != length(b)) stop("A and b are incompatible!")
  if (nrow(H) != length(u)) stop("u is incopatible with H")
  if (nrow(H) != length(l)) stop("l is incopatible with H")
  c <- matrix(c)
  l <- matrix(l)
  u <- matrix(u)
  m <- nrow(A)
  n <- ncol(A)
  b.plus.1 <- max(svd(b)$d) + 1
  c.plus.1 <- max(svd(c)$d) + 1
  H.y <- diag(1, m)
  c.x <- c
  c.y <- b
  if (is.square) {
    H.x <- H
    diag(H.x) <- diag(H) + 1
    AP <- matrix(0, m + n, m + n)
    AP[1:n,1:n] <- -H.x
    AP[-(1:n),1:n] <- A
    AP[1:n,-(1:n)] <- t(A)
    AP[-(1:n),-(1:n)] <- H.y
    s.tmp <- solve(AP, c(c.x, c.y))
    x <- s.tmp[1:n]
    y <- s.tmp[-(1:n)]
  } else {
    H.x <- t(H)
    V <- diag(ncol(H))
    smwinner <- chol(V + crossprod(H))
    smwa1 <- t(A)
    smwa2 <- smwa1 - (H %*% solve(smwinner, solve(t(smwinner), crossprod(H, smwa1))))
    smwc2 <- c.x - (H %*% solve(smwinner, solve(t(smwinner), crossprod(H, c.x))))
    y <- solve(A %*% smwa2 + H.y, c.y + A %*% smwc2)
    x <- smwa2 %*% y - smwc2
  }
  g <- pmax(abs(x - l), bound)
  z <- pmax(abs(x), bound)
  t <- pmax(abs(u - x), bound)
  s <- pmax(abs(x), bound)
  v <- pmax(abs(y), bound)
  w <- pmax(abs(y), bound)
  p <- pmax(abs(r - w), bound)
  q <- pmax(abs(y), bound)
  mu <- as.vector(crossprod(z, g) + crossprod(v, w) + crossprod(s, t) + crossprod(p, q))/(2 * (m + n))
  sigfig <- 0
  alfa <- 1
  if (verb > 0) cat("Iter    PrimalInf  DualInf  SigFigs  Rescale  PrimalObj  DualObj\n")
  for(counter in seq_len(maxiter)) {
    if (is.square) H.dot.x <- H %*% x else H.dot.x <- H %*% crossprod(H, x)
    rho <- b - A %*% x + w
    nu <- l - x + g
    tau <- u - x - t
    alpha <- r - w - p
    sigma <- c - crossprod(A, y) - z + s + H.dot.x
    beta <- y + q - v
    x.dot.H.dot.x <- crossprod(x, H.dot.x)
    primal.infeasibility <- max(svd(rbind(rho, tau, matrix(alpha), nu))$d)/b.plus.1
    dual.infeasibility <- max(svd(rbind(sigma, t(t(beta))))$d)/c.plus.1
    primal.obj <- crossprod(c, x) + 0.5 * x.dot.H.dot.x
    dual.obj <- crossprod(b, y) - 0.5 * x.dot.H.dot.x + crossprod(l, z) - crossprod(u, s) - crossprod(r, q)
    sigfig <- max(-log10(abs(primal.obj - dual.obj)/(abs(primal.obj) + 1)), 0)
    if (sigfig >= sigf) break
    if (verb > 0) cat(counter, "\t", signif(primal.infeasibility, 6), signif(dual.infeasibility, 6), sigfig, alfa, primal.obj, dual.obj, "\n")
    hat.beta <- beta + v
    hat.alpha <- alpha + p
    hat.nu <- nu - g
    hat.tau <- tau + t
    d <- z/g + s/t
    e <- 1/(v/w + q/p)
    if (is.square) diag(H.x) <- diag(H) + d
    diag(H.y) <- e
    c.x <- sigma - z * hat.nu/g - s * hat.tau/t
    c.y <- rho - e * (hat.beta - q * hat.alpha/p)
    if (is.square) {
      AP[1:n,1:n] <- -H.x
      AP[-(1:n),-(1:n)] <- H.y
      s1.tmp <- solve(AP, c(c.x, c.y))
      delta.x <- s1.tmp[1:n]
      delta.y <- s1.tmp[-(1:n)]
    } else {
      V <- diag(ncol(H))
      smwinner <- chol(V + chunkmult(t(H), 2000, d))
      smwa2 <- t(A) - (H %*% solve(smwinner, solve(t(smwinner), crossprod(H, t(A)/d))))
      smwa2 <- smwa2/d
      smwc2 <- (c.x - (H %*% solve(smwinner, solve(t(smwinner), crossprod(H, c.x/d)))))/d
      delta.y <- solve(A %*% smwa2 + H.y, c.y + A %*% smwc2)
      delta.x <- smwa2 %*% delta.y - smwc2
    }
    delta.w <- -e * (hat.beta - q * hat.alpha/p + delta.y)
    delta.s <- s * (delta.x - hat.tau)/t
    delta.z <- z * (hat.nu - delta.x)/g
    delta.q <- q * (delta.w - hat.alpha)/p
    delta.v <- v * (-w - delta.w)/w
    delta.p <- p * (-q - delta.q)/q
    delta.g <- g * (-z - delta.z)/z
    delta.t <- t * (-s - delta.s)/s
    alfa <- -(1 - margin)/min(c(delta.g/g, delta.w/w, delta.t/t, delta.p/p, delta.z/z, delta.v/v, delta.s/s, delta.q/q, -1))
    newmu <- (crossprod(z, g) + crossprod(v, w) + crossprod(s, t) + crossprod(p, q))/(2 * (m + n))
    newmu <- mu * ((alfa - 1)/(alfa + 10))^2
    gamma.z <- mu/g - z - delta.z * delta.g/g
    gamma.w <- mu/v - w - delta.w * delta.v/v
    gamma.s <- mu/t - s - delta.s * delta.t/t
    gamma.q <- mu/p - q - delta.q * delta.p/p
    hat.beta <- beta - v * gamma.w/w
    hat.alpha <- alpha - p * gamma.q/q
    hat.nu <- nu + g * gamma.z/z
    hat.tau <- tau - t * gamma.s/s
    c.x <- sigma - z * hat.nu/g - s * hat.tau/t
    c.y <- rho - e * (hat.beta - q * hat.alpha/p)
    if (is.square) {
      AP[1:n, 1:n] <- -H.x
      AP[-(1:n),-(1:n)] <- H.y
      s1.tmp <- solve(AP, c(c.x, c.y))
      delta.x <- s1.tmp[1:n]
      delta.y <- s1.tmp[-(1:n)]
    } else {
      smwc2 <- (c.x - (H %*% solve(smwinner, solve(t(smwinner), crossprod(H, c.x/d)))))/d
      delta.y <- solve(A %*% smwa2 + H.y, c.y + A %*% smwc2)
      delta.x <- smwa2 %*% delta.y - smwc2
    }
    delta.w <- -e * (hat.beta - q * hat.alpha/p + delta.y)
    delta.s <- s * (delta.x - hat.tau)/t
    delta.z <- z * (hat.nu - delta.x)/g
    delta.q <- q * (delta.w - hat.alpha)/p
    delta.v <- v * (gamma.w - delta.w)/w
    delta.p <- p * (gamma.q - delta.q)/q
    delta.g <- g * (gamma.z - delta.z)/z
    delta.t <- t * (gamma.s - delta.s)/s
    alfa <- -(1 - margin)/min(c(delta.g/g, delta.w/w, delta.t/t, delta.p/p, delta.z/z, delta.v/v, delta.s/s, delta.q/q, -1))
    x <- x + delta.x * alfa
    g <- g + delta.g * alfa
    w <- w + delta.w * alfa
    t <- t + delta.t * alfa
    p <- p + delta.p * alfa
    y <- y + delta.y * alfa
    z <- z + delta.z * alfa
    v <- v + delta.v * alfa
    s <- s + delta.s * alfa
    q <- q + delta.q * alfa
    mu <- newmu
  }
  if (verb > 0) cat(counter, primal.infeasibility, dual.infeasibility, sigfig, alfa, primal.obj, dual.obj)
  
  ret <- list(
    primal = x,
    dual = drop(y),
    how = if (sigfig >= sigf && counter<maxiter) {
              "converged"
            } else {
              msg <- matrix(c("slow convergence, change bound?",
                      "primal infeasible",
                      "dual infeasible",
                      "primal and dual infeasible"),2,2)
              msg[ifelse(primal.infeasibility>1e+06,1,2),ifelse(dual.infeasibility>1e+06,1,2)]
            }
  )
  return(ret)
}