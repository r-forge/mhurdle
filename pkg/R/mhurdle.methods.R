sumrho <- function(x){
  ScoreTest <- x['prime']/sqrt(-x['second'])
  names(ScoreTest) <- "z"
  ImpliedRho <- -x['prime']/x['second']
  ScoreTest <- list(statistic = ScoreTest,
                    p.value = 1 - pnorm(ScoreTest),
                    method = "Score Test"
                    )
  class(ScoreTest) <- "htest"
  list(test = ScoreTest, value = ImpliedRho)
}

nm.mhurdle <- function(object,
                    which = c("all", "sel", "ifr", "reg", "other", "sigma", "rho"),
                    ...){

  which <- match.arg(which)
  K <- sapply(object$coef.names,length)
  soi <- ifelse(describe(object, "sel"), "sel", "ifr")
  if (which == "all"){
    soi.names <- paste(soi, object$coef.names[[soi]],sep = ".")
    reg.names <- paste("reg", object$coef.names$reg,sep = ".")
    result <- c(soi.names, reg.names, object$coef.names$other)
  }
  if (which == "sel") result <- object$coef.names[[1]]
  if (which == "ifr") result <- object$coef.names[[1]]
  if (which == "reg") result <- object$coef.names[[2]]
  if (which == "other"){
    if (describe(object, "corr")) result <- c("sigma","rho")
    else result <- "sigma"
  }
  if (which == "rho") result <- "rho"
  if (which == "sigma") result <- "sigma"
  result
}

sub.mhurdle <- function(object,
                      which = c("all", "sel", "ifr", "reg", "other", "sigma", "rho"),
                      ...){

  which <- match.arg(which)
  K <- sapply(object$coef.names,length)
  if (which == "all") sub <- NULL
  if (which == "sel"){
    if (!describe(object, "sel")) stop("no sel equation")
    sub <- 1:K[[1]]
  }
  if (which == "ifr"){
    if (!describe(object, "ifr")) stop("no ifruency equation")
    sub<- 1:K[[1]]
  }
  if (which == "reg"){
    sub <- (K[[1]]+1):(K[[1]]+K[[2]])
  }
  if (which == "other"){
    sub <- -c(1:(K[[1]]+K[[2]]))
  }
  if (which == "rho"){
    if (!describe(object, "corr")) stop("no corr coefficient estimated")
    sub <- K[[1]]+K[[2]]+2
  }
  if (which == "sigma"){
    sub <- K[[1]]+K[[2]]+1
  }    
  sub
}

coef.mhurdle <- function(object,
                      which = c("all", "sel", "ifr", "reg", "other", "sigma", "rho"),
                      ...){
  which <- match.arg(which)
  sub <- sub.mhurdle(object, which)
  nm <- nm.mhurdle(object, which)
  result <- object$coefficients
  if (!is.null(sub)) result <- result[sub]
  names(result) <- nm
  result
}

vcov.mhurdle <- function(object,
                      which = c("all", "sel", "ifr", "reg", "other", "sigma", "rho"),
                      ...){
  which <- match.arg(which)
  sub <- sub.mhurdle(object, which)
  nm <- nm.mhurdle(object, which)
  result <- object$vcov
  if (!is.null(sub)) result <- result[sub, sub, drop = FALSE]
  rownames(result) <- colnames(result) <- nm
  if (length(result) == 1){
    nm <- rownames(result)
    result <- as.numeric(result)
    names(result) <- nm
  }
  result
}



logLik.mhurdle <- function(object, which = c("all", "zero", "positive"), ...){
  which <- match.arg(which)
  x <- object$logLik
  y <- attr(x, "y")
  x <- switch(which,
              all      = sum(x),
              zero     = sum(x[y == 0]),
              positive = sum(x[y > 0])
              )
  attr(x,"df") <- NULL
  attr(x,"y")  <- NULL
  sum(x)
}

print.mhurdle <- function (x, digits = max(3, getOption("digits") - 2),
                        width = getOption("width"), ...){
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

summary.mhurdle <- function (object,...){
  b <- coef(object)
  std.err <- sqrt(diag(vcov(object)))
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  object$CoefTable <- CoefTable
  class(object) <- c("summary.mhurdle","mhurdle")
  return(object)
}

print.summary.mhurdle <- function(x, digits = max(3, getOption("digits") - 2),
                                  width = getOption("width"), ...){

  cat("\nCall:\n")
  print(x$call)

  cat("\n")
  
  y <- x$model[,1]
  zeros <- length(y[y==0])/length(y)
  if (zeros>0) cat(paste("Frequency of 0: ",round(zeros,digits=digits),"\n"))
  
  if (!is.null(x$est.stat)){
    cat("\n")
    print(x$est.stat)
  }
  
  cat("\nCoefficients :\n")
  printCoefmat(x$CoefTable,digits=digits)
  cat("\n")
  df <- attr(x$logLik,"df")

  cat(paste("Log-Likelihood: ",
            signif(logLik(x),digits),
            " on ",df," Df\n",sep=""))

  if (!is.null(x$rho)){
    rhotest <- sumrho(x$rho)
    cat("rho: ")
    cat(paste("implied value :", round(rhotest$value,3)," ",sep=""))
    z <- 
    cat(paste("score test : z = ", round(rhotest$test$statistic,3),
              " (p.value = ",round(rhotest$test$p.value,3),")\n",sep=""))
  }
  invisible(x)
}

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

## a simple copy from mlogit. update with formula doesn't work
## otherwise ????

update.mhurdle <- function (object, new, ...){
  call <- object$call
  if (is.null(call))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(new))
    call$formula <- update(formula(object), new)
  if(length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    ## do these individually to allow NULL to remove entries.
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  eval(call, parent.frame())
}

compute.fitted.mhurdle <- function(param, X, S, dist, res, sel, ifr, corr){
  z <- typo.mhurdle()
  model <- z[as.character(res),
             as.character(sel),
             as.character(ifr),
             as.character(corr),
             dist]
  KX <- ncol(X)
  KS <- ncol(S)
  betaS <- param[1:KS]
  betaX <- param[(KS+1):(KS+KX)]
  sigma <- param[KX+KS+1]
  bX <- as.numeric(crossprod(t(X),betaX))
  bS <- as.numeric(crossprod(t(S),betaS))
  PhiX <- pnorm(bX/sigma)
  phiX <- dnorm(bX/sigma)
  PhiS <- pnorm(bS)
  phiS <- dnorm(bS)
  if (corr){
    rho <- param[KX+KS+2]
    Phib <- pbivnorm(bS, bX/sigma, rho, gradient = FALSE)
    psi <- psy(bS,bX/sigma,rho)
    millsG <- psi/Phib
  }
  else{
    Phib <- pnorm(bS) * pnorm(bX/sigma)
    millsG <- mills(bS)
  }

  prob.null <- 1 - PhiS
  if (dist == "t" && corr) prob.null <- 1 - Phib/PhiX
  if (res) prob.null <- 1 - Phib
  
  if (dist == "l") esp.cond <- exp(bX+sigma^2/2) else{
    esp.cond <- bX + sigma * millsG
  }
  if (ifr) esp.cond <- esp.cond / PhiS
  if (dist == "l" && corr){
    dl <- (PhiS+rho*sigma*phiS+
           (rho*sigma)^2/2*(PhiS-bS*phiS)+
           (rho*sigma)^3/6*(2+bS^2)*phiS+
           (rho*sigma)^4/12*(3*PhiS-bS*(3+bS^2)*phiS))/PhiS
    if (ifr) esp.cond <- exp(bX+sigma^2*(1-rho^2))*dl/PhiS
    if (sel) esp.cond <- exp(bX+sigma^2*(1-rho^2))*dl
  }
  
  result <- cbind(zero = prob.null,
                  positive = esp.cond)
  
  result
}

predict.mhurdle <- function(object, newdata = NULL, ...){
  if (is.null(newdata)){
    result <- fitted(object, ...)
  }
  else{
    cl <- object$call
    sel <- ifelse(is.null(cl$sel), TRUE, cl$sel)
    res <- ifelse(is.null(cl$resd), TRUE, cl$res)
    ifr <- ifelse(is.null(cl$ifr), TRUE, cl$ifr)
    dist <- ifelse(is.null(cl$dist), TRUE, cl$dist)
    corr <- ifelse(is.null(cl$corr), FALSE, cl$corr)
    m <- model.frame(formula(object), newdata)
    X <- model.matrix(formula(object),m)
    S <- model.matrix(formula(object),m, part = "second")
    result <- compute.fitted.mhurdle(coef(object), X, S, dist, res, sel, ifr, corr)
    
  }
  result
}

fitted.mhurdle <- function(object, which = c("zero", "positive"), ...){
  which <- match.arg(which)
  switch(which,
         zero = object$fitted[,1],
         positive = object$fitted[,2])
}
            
r.squared <- function(object,
                      which = c("all", "zero", "positive"),
                      type = c("linear", "mcfadden")){
  type <- match.arg(type)
  which <- match.arg(which)
  
#  if (!is.logical(adjusted)) stop("adjusted should be logical") else

  if (type == "mcfadden"){
    result <- 1 - logLik(object, which)/logLik(object$naive, which)
  }
  if (type == "linear"){
    y <- model.response(model.frame(object))
    p <- mean(y == 0)
    if (which != "positive"){
      # number of null observations
      no <- sum(y == 0)
      # mean probabilities of 0
      pmo <- mean(y == 0)
      # fitted values for 0
      pfo <- fitted(object, "zero")
      # residuals for 0
      eo <- 0 - pfo
      R2n <- sum( (pfo - pmo)^2) / (no * pmo^2)
      ## alternative based on residuals
      ## R2n <- 1- sum(eo^2) / (no * pmo^2)
    }
    if (which != "zero"){
      yp <- y[y>0]
      yf <- fitted(object, "positive")[y > 0]
      ym <- mean(yp)
      ## expression based on fitted values
      ## R2p <- sum( (yf - ym)^2 ) / sum( (yp - ym)^2)
      ## alternative based on residuals
      ## R2p <- 1 - sum( (yp-yf)^2 )/sum( (yp-ym)^2)
      ## alternative base on the coefficient of correlation
      R2p <- cor(yf,yp)^2
    }
    if (which == "all"){
      R2 <- p * R2n + (1-p) * R2p
    }
    result <- switch(which,
                     all = R2,
                     zero = R2n,
                     positive = R2p
                     )
  }
  result
}
  
  
vuongtest <- function(x, y){
  lx <- x$logLik
  ly <- y$logLik
  attr(lx, "y") <- attr(ly, "y") <- attr(lx, "df") <- attr(ly, "df") <- NULL
  n <- length(lx)
  if (length(ly) != n) stop("the number of observations of the two models differ")
  LR <- sum(lx - ly)
  w2 <- 1/n*sum((lx-ly)^2)-(1/n*LR)^2
  value <- LR/sqrt(n*w2)
  names(value) <- "normal"
  pval<- pnorm(abs(value), lower.tail = FALSE)*2
  data <- paste(deparse(x$call$formula))
  result <- list(statistic = value,
                 p.value = pval,
                 data.name = NULL,
                 method = "Vuong Test")
  class(result) <- "htest"
  result
}
