sim.mhurdle <- function(param = c(.5, .5, .5, .5, .5),
                     rho = 0,
                     n = 200,
                     sel = TRUE, ifr = FALSE,
                     dist = c("l","t","n"),
                     res = FALSE){

  if (rho <= -1 | rho >= 1) stop("rho must be in the -1 +1 range")
  corr <- rho != 0

  if (!is.null(corr) && !is.logical(corr))
    stop("corr should be a logical value")
  if (!is.null(res) && !is.logical(res))
    stop("res should be a logical value")
  if (!is.null(sel) && !is.logical(sel))
    stop("sel should be a logical value")
  if (!is.null(ifr) && !is.logical(ifr))
    stop("ifr should be a logical value")

  dist <- match.arg(dist)
  z <- typo.mhurdle()
  
  model <- z[as.character(res),as.character(sel),
             as.character(ifr),as.character(corr),
             dist]

  if (is.na(model)) stop("not implemented  model")
  if (model == "irrelevant") stop("irrelevant model")

  x <- rnorm(n)
  s <- rnorm(n)
  int1 <- param[1]
  beta1 <- param[2]
  int2 <- param[3]
  beta2 <- param[4]
  sigma <- param[5]
  bs <- int1 + beta1 * s
  bx <- int2 + beta2 * x

  if (dist != "t"){
    if (corr){
      V <- matrix(c(1,rho*sigma,rho*sigma,sigma^2),ncol=2)
      cV <- chol(V)
      eps <- matrix(rnorm(2*n), ncol = 2)
      eps <- crossprod(t(eps), cV)
      eps1 <- eps[,1]
      eps2 <- eps[,2]
    }
    else{
      eps1 <- rnorm(n)
      eps2 <- rnorm(n, sd = sigma)
    }
  }
  else{
    u <- runif(n)
    eps2 <- sigma*qnorm(pnorm(bx/sigma)*u+pnorm(-bx/sigma))
    if (corr){
      eps1 <- sqrt(1-rho^2)*rnorm(n,mean=rho*eps2,sd=sqrt(1-rho^2))+rho*eps2
    }
    else{
      eps1 <- rnorm(n)
    }
  }
  y1s <- bs + eps1
  y2s <- bx + eps2
  if (dist == "l") y2s <- exp(y2s)
  if (res) I <-  (y1s > 0 & y2s > 0)  else I <- (y1s > 0)
  if (ifr) y <- I * y2s / pnorm(bs) else y <- I * y2s
  data.frame(y = y, x = x, s = s)
}

sim.est <- function(param = c(.5, .5, .5, .5, .5),
                     rho = 0,
                     n = 200,
                     sel = TRUE, ifr = FALSE,
                     dist = c("l","t","n"),
                     res = FALSE){
  dd <- sim.mhurdle(param = param,
                     rho = rho,
                     n = n,
                     sel = sel, ifr = ifr,
                     dist = dist,
                     res = res)
  corr <- rho != 0
  result <- mhurdle(y ~ x | s, dd,
                 sel = sel, ifr = ifr,
                 dist = dist,
                 res = res,
                 corr = corr,
                 print.level=2)
  result
}

