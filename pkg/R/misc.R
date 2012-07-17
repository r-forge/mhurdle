log2 <- function(x) ifelse(x > 0,log(x), 0)

mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))

# compute the bivariate normal cdf and its derivates
mypbivnorm <- function(x1, x2, rho = 0){
  if (is.null(x1) &&  is.null(x2)) result <- list(f = 1, a = 0, b = 0, rho = 0, rho.rho = 0)
  if (is.null(x1) && !is.null(x2)) result <- list(f = pnorm(x2), a = 0, b = dnorm(x2), rho = 0, rho.rho = 0)
  if (is.null(x2) && !is.null(x1)) result <- list(f = pnorm(x1), a = dnorm(x1), b = 0, rho = 0, rho.rho = 0)
  if (!is.null(x1) && !is.null(x2)){
    f <- pbivnorm(x1, x2, rho)
    b <- dnorm(x2) * pnorm( (x1 - rho * x2) / sqrt(1 - rho^2) )
    a <- dnorm(x1) * pnorm( (x2 - rho * x1) / sqrt(1 - rho^2) )
    eps <- 1E-07
    rho <- (pbivnorm(x1, x2, rho + 1E-07) - pbivnorm(x1, x2, rho)) / 1E-07
    result <- list(f = f, a = a, b = b, rho = rho)
  }
  result
}

# a function to construct block-diagonal matrix (can't remember from
# which package it is borrowed)
bdiag <- function(...){
  if (nargs() == 1)
    x <- as.list(...)
  else
    x <- list(...)
  n <- length(x)
  if(n == 0) return(NULL)
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

# a function used to compute the expected value of the variable for
# some models
psy <- function(x1 , x2, rho, terms = "tot", degree = 3){
  d0 <- degree >= 0
  d1 <- degree >= 1
  d2 <- degree >= 2
  d3 <- degree >= 3
  
  A0 <- pnorm(x2) * dnorm(x1)
  A1 <- dnorm(x2) * (pnorm(x1) - x1 * dnorm(x1))
  A2 <- -dnorm(x2) * x2 * dnorm(x1)*(1 + x1 ^ 2)
  A3 <- dnorm(x2) * (1 - x2 ^ 2) * x1 ^ 3 * dnorm(x1)

  B0 <- dnorm(x2) * pnorm(x1)
  B1 <- - x2 * dnorm(x2) * dnorm(x1)
  B2 <- - dnorm(x2) * (x1 * dnorm(x1) * (x2^2 - 1) + pnorm(x1))
  B3 <- dnorm(x2) * x2 * dnorm(x1) * (3 * x1 ^ 2 + x2 ^ 2 - x2 ^ 2 * x1 ^ 2)
  
  A <- rho * (d0 * A0 + A1 * d1 * rho + A2 * d2 * rho ^ 2 / 2 + A3 * d3 * rho ^ 3 / 6)
  B <- sqrt(abs(1 - rho ^ 2) )*(d0 * B0 + B1 * rho + B2 * d2 * rho ^ 2 / 2 + B3 * d3 * rho ^ 3 / 6)
  
  result <- switch(terms,
                   "tot" = A + B,
                   "A"   = A,
                   "B"   = B
                   )
  result
}
