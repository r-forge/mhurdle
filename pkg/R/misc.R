log2 <- function(x) ifelse(x>0,log(x),0)

mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))


if (FALSE){
  mypbivnorm <- function(x1, x2, rho){
    if (is.null(x1) &&  is.null(x2)) result <- list(f = 1, a = 0, b = 0, rho = 0, rho.rho = 0)
    if (is.null(x1) && !is.null(x2)) result <- list(f = pnorm(x2), a = 0, b = dnorm(x2), rho = 0, rho.rho = 0)
    if (is.null(x2) && !is.null(x1)) result <- list(f = pnorm(x1), a = dnorm(x1), b = 0, rho = 0, rho.rho = 0)
    if (!is.null(x1) && !is.null(x2)){
      Phix1 <- pnorm(x1)
      Phix2 <- pnorm(x2)
      phix1 <- dnorm(x1)
      phix2 <- dnorm(x2)
      frho <- rho + 1 / 2 * rho^2 * x1 * x2 + 1/6 * rho^3 * (x1^2 - 1)*(x2^2 - 1)
      f <- Phix1 * Phix2 + phix1 * phix2 * frho
      frho.x1 <- 1/2 * rho^2 * x2 + 1/3 * rho^3 * x1 * (x2^2-1)
      frho.x2 <- 1/2 * rho^2 * x1 + 1/3 * rho^3 * x2 * (x1^2-1)
      frho.rho <- 1 + rho * x1 * x2 + 1/2 * rho^2 * (x1^2-1) * (x2^2-1)
      frho.rho.rho <- x1 * x2 + rho * (x1^2-1) * (x2^2-1)
      result <- list(f = f,
                     a   = phix1 * Phix2-x1 * phix1 * phix2 * frho + phix1 * phix2 * frho.x1,
                     b   = phix2 * Phix1-x2 * phix1 * phix2 * frho + phix1 * phix2 * frho.x2,
                     rho = phix1 * phix2 * frho.rho,
                     rho.rho = phix1 * phix2 * frho.rho.rho)
    }
    result
  }
} else {
  mypbivnorm <- function(x1, x2, rho){
    if (is.null(x1) &&  is.null(x2)) result <- list(f = 1, a = 0, b = 0, rho = 0, rho.rho = 0)
    if (is.null(x1) && !is.null(x2)) result <- list(f = pnorm(x2), a = 0, b = dnorm(x2), rho = 0, rho.rho = 0)
    if (is.null(x2) && !is.null(x1)) result <- list(f = pnorm(x1), a = dnorm(x1), b = 0, rho = 0, rho.rho = 0)
    if (!is.null(x1) && !is.null(x2)){
      f <- pbivnorm(x1, x2, rho)
      b <- dnorm(x2) * pnorm( (x1 - rho * x2) / sqrt(1 - rho^2) )
      a <- dnorm(x1) * pnorm( (x2 - rho * x1) / sqrt(1 - rho^2) )
      eps <- 1E-07
      rho <- (pbivnorm(x1, x2, rho + 1E-07) - pbivnorm(x1, x2, rho)) / 1E-07
      result <- list(f=f, a=a, b=b, rho=rho)
    }
    result
  }
}

bdiag <- function(...){
  if (nargs() == 1)
    x <- as.list(...)
  else
    x <- list(...)
  n <- length(x)
  if(n==0) return(NULL)
  x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
              stop("Zero-length component in x"))
  d <- array(unlist(lapply(x, dim)), c(2, n))
  rr <- d[1,]
  cc <- d[2,]
  rsum <- sum(rr)
  csum <- sum(cc)
  out <- array(0, c(rsum, csum))
  ind <- array(0, c(4, n))
  rcum <- cumsum(rr)
  ccum <- cumsum(cc)
  ind[1,-1] <- rcum[-n]
  ind[2,] <- rcum
  ind[3,-1] <- ccum[-n]
  ind[4,] <- ccum
  imat <- array(1:(rsum * csum), c(rsum, csum))
  iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
                                               (y[3]+1):y[4]], imat=imat)
  iuse <- as.vector(unlist(iuse))
  out[iuse] <- unlist(x)
  return(out)
} 


