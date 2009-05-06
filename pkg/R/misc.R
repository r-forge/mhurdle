log2 <- function(x) ifelse(x>0,log(x),0)

mills <- function(x) dnorm(x)/pnorm(x)

ldnorm <- function(x) -0.5*log(2*pi)-0.5*x^2


print.est.stat <- function(x, ...){
  et <- x$elaps.time[3]
  i <- x$nb.iter[1]
  eps <- as.numeric(x$eps)
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  tstr <- paste(h,"h:",m,"m:",s,"s",sep="")
  cat(paste(x$method,"\n"))
  cat(paste(x$message,"\n"))
  cat(paste(i,"iterations,",tstr,"\n"))
  cat(paste("g'(-H)^-1g =",sprintf("%5.3G",eps),"\n"))
}


der <- function(f, param, eps = 1E-6, ...){
  fp <- do.call(f, list(param, ...))
  d <- rep(NA,length(param))
  for (i in 1:length(param)){
    parami <- param
    params <- param
    parami[i] <- parami[i]-eps
    params[i] <- params[i]+eps
    d[i] <- (do.call(f, list(params, ...)) - do.call(f, list(parami, ...)))/(2*eps)
  }
  d
}


pbivnorm <- function(x1, x2, rho, gradient = TRUE){
  Phix1 <- pnorm(x1)
  Phix2 <- pnorm(x2)
  phix1 <- dnorm(x1)
  phix2 <- dnorm(x2)
  frho <- rho+1/2*rho^2*x1*x2+1/6*rho^3*(x1^2-1)*(x2^2-1)
  
  res <- Phix1*Phix2+phix1*phix2*frho

  if (gradient){
    frho.x1 <- 1/2*rho^2*x2+1/3*rho^3*x1*(x2^2-1)
    frho.x2 <- 1/2*rho^2*x1+1/3*rho^3*x2*(x1^2-1)

    frho.rho <- 1+rho*x1*x2+1/2*rho^2*(x1^2-1)*(x2^2-1)
    frho.rho.rho <- x1*x2+rho*(x1^2-1)*(x2^2-1)
    gradient <- list(x1 = phix1*Phix2-x1*phix1*phix2*frho+phix1*phix2*frho.x1,
                     x2 = phix2*Phix1-x2*phix1*phix2*frho+phix1*phix2*frho.x2,
#                     x1 = pnorm(x1) * dnorm( (x2-rho*x1)/sqrt(1-rho^2)),
#                     x2 = pnorm(x2) * dnorm( (x1-rho*x2)/sqrt(1-rho^2)),
                     rho = phix1*phix2*frho.rho,
                     rho.rho = phix1*phix2*frho.rho.rho)
    attr(res,"gradient") <- gradient
  }
  res
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


pb <- function(x1,x2,mean=c(0,0),sigma=c(1,1),rho){
  x1 <- (x1-mean[1])/sigma[1]
  x2 <- (x2-mean[2])/sigma[2]
  Phix1 <- pnorm(x1)
  Phix2 <- pnorm(x2)
  phix1 <- dnorm(x1)
  phix2 <- dnorm(x2)
  frho <- rho+1/2*rho^2*x1*x2+1/6*rho^3*(x1^2-1)*(x2^2-1)
  res <- Phix1*Phix2+phix1*phix2*frho/sigma[1]/sigma[2]
  res
}




psy <- function(x1 , x2, rho, terms = "tot", degree = 3){
  d0 <- degree>=0
  d1 <- degree>=1
  d2 <- degree>=2
  d3 <- degree>=3
  
  A0 <- pnorm(x2)*dnorm(x1)
  A1 <- dnorm(x2)*(pnorm(x1)-x1*dnorm(x1))
  A2 <- -dnorm(x2)*x2*dnorm(x1)*(1+(x1)^2)
  A3 <- dnorm(x2)*(1-(x2)^2)*x1^3*dnorm(x1)

  B0 <- dnorm(x2)*pnorm(x1)
  B1 <- -x2*dnorm(x2)*dnorm(x1)
  B2 <- -dnorm(x2)*(x1*dnorm(x1)*((x2)^2-1)+pnorm(x1))
  B3 <- dnorm(x2)*x2*dnorm(x1)*(3*(x1)^2+(x2)^2-(x2)^2*(x1)^2)
  
  A <- rho*(d0*A0 + A1*d1*rho + A2*d2*rho^2/2 + A3*d3*rho^3/6)
  B <- sqrt(abs(1-rho^2))*(d0*B0 + B1*rho + B2*d2*rho^2/2 + B3*d3*rho^3/6)
  
  result <- switch(terms,
                   "tot"=A+B,
                   "A"=A,
                   "B"=B
                   )
  result
}

typo.mhurdle <- function(){
  z <- array(NA,dim=c(2,2,2,2,3),dimnames=list(c("FALSE","TRUE"),c("FALSE","TRUE"),c("FALSE","TRUE"),c("FALSE","TRUE"),c("n","t","l")))
  
  z["TRUE","FALSE","FALSE","FALSE","n"] <- "tobit"
  z["FALSE","TRUE","FALSE","TRUE","n"] <- "tobit2"
  z["FALSE","TRUE","FALSE","FALSE","n"] <- "tobit2i"
  z["TRUE","TRUE","FALSE","TRUE","n"] <- "snd"
  z["TRUE","TRUE","FALSE","FALSE","n"] <- "sni"
  z["FALSE","TRUE","FALSE","FALSE","t"] <- "sti"
  z["FALSE","TRUE","FALSE","TRUE","t"] <- "std"
  z["FALSE","TRUE","FALSE","FALSE","l"] <- "sli"
  z["FALSE","TRUE","FALSE","TRUE","l"] <- "sld"
  z["TRUE","FALSE","TRUE","FALSE","n"] <- "pni"
  z["FALSE","FALSE","TRUE","TRUE","t"] <- "ptd"
  z["FALSE","FALSE","TRUE","FALSE","t"] <- "pti"
  z["FALSE","FALSE","TRUE","TRUE","l"] <- "pld"
  z["FALSE","FALSE","TRUE","FALSE","l"] <- "pli"
  z["FALSE","FALSE","TRUE","FALSE","n"] <- "pni1"
  z["FALSE","FALSE","TRUE","TRUE","n"] <- "pnd1"
  
  

  z["FALSE","FALSE","FALSE","FALSE","n"] <- "irrelevant"
  z["TRUE","FALSE","FALSE","TRUE","n"] <- "irrelevant"

  z["FALSE","FALSE","FALSE",,c("t","l")] <- "irrelevant"
  z["TRUE",,,,c("t","l")] <- "irrelevant"
  z
}




namemhurdle <- function(x){
  res <- as.logical(substr(x,1,1))
  sel <- as.logical(substr(x,2,2))
  ifr <- as.logical(substr(x,3,3))
  corr <- as.logical(substr(x,4,4))
  dist <- substr(x,5,5)
  toprint <- c()
  if (res) toprint <- c(toprint,"res")
  if (sel) toprint <- c(toprint,"sel")
  if (ifr) toprint <- c(toprint,"ifr")
  if (corr) toprint <- c(toprint,"corr")
  toprint <- c(toprint,dist)
  paste(toprint,sep=" ",collapse=" ")
}
