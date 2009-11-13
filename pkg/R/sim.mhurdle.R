sim.mhurdle <- function(reg = c(0.5, 0.5, 1), sel = NULL, ifr = NULL, rho = 0,
                        n = 200,
                        dist = c("l","t","n")){
  dist <- match.arg(dist)
  if (rho <= -1 | rho >= 1) stop("rho must be in the -1 +1 range")
  corr <- rho != 0

  if (is.null(sel) && corr) stop("correlated models only implemented for the selection process")

  x <- rnorm(n)
  bx <- reg[1] + reg[2] * x
  sigma <- reg[3]

  if (!is.null(sel)){
    s <- rnorm(n)
    bs <- sel[1] + sel[2] * s
  }
  
  if (!is.null(ifr)){
    p <- rnorm(n)
    bp <- ifr[1] + ifr[2] * p
    epsp <- rnorm(n)
    PR <- pnorm(bp)
    yp <- bp + epsp
  }
  else{
    yp <- 1
    PR <- 1
  }
  
  if (dist != "t"){
    if (corr){
      # if corr, sel is TRUE
      V <- matrix(c(1, rho*sigma, rho*sigma, sigma^2), ncol = 2)
      cV <- chol(V)
      eps <- matrix(rnorm(2*n), ncol = 2)
      eps <- crossprod(t(eps), cV)
      epss <- eps[,1]
      epsx <- eps[,2]
    }
    else{
      epsx <- rnorm(n, sd = sigma)
      if (!is.null(sel)) epss <- rnorm(n)
    }
  }
  else{
    u <- runif(n)
    epsx <- sigma*qnorm(pnorm(bx/sigma)*u+pnorm(-bx/sigma))
    if (corr){
      epss <- sqrt(1-rho^2)*rnorm(n, mean = rho * epsx, sd = sqrt(1-rho^2)) + rho * epsx
    }
    else{
      if (!is.null(sel)) epss <- rnorm(n)
    }
  }
  
  yx <- bx + epsx
  if (!is.null(sel)) ys <- bs + epss else ys <- 1

  if (dist == "l") y <- (ys > 0 & exp(yx) > 0 & yp >0)* exp(yx)/PR
  else y <- (ys > 0 & yx > 0 & yp >0)* yx/PR
  data <- data.frame(y = y, x = x)
  if (!is.null(sel)) data <- cbind(data, s = s)
  if (!is.null(ifr)) data <- cbind(data, p = p)
  data
}


