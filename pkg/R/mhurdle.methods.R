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
  if (which == "all"){
    reg.names <- paste("reg", object$coef.names$reg,sep = ".")
    sel.names <- ifr.names <- NULL
    if (describe(object, "sel")) sel.names <- paste("sel", object$coef.names$sel,sep = ".")
    if (describe(object, "ifr")) ifr.names <- paste("ifr", object$coef.names$ifr,sep = ".")

    result <- c(sel.names, reg.names, ifr.names, object$coef.names$other)
  }
  if (which == "sel") result <- object$coef.names[[1]]
  if (which == "reg") result <- object$coef.names[[2]]
  if (which == "ifr") result <- object$coef.names[[3]]
  if (which == "other")
    if (describe(object, "corr")) {result <- c("sigma","rho")}
    else {result <- "sigma"}
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
  if (which == "reg") sub <- 1:K[[1]] 
  if (which == "sel"){
    if (!describe(object, "sel")) stop("no sel equation")
    else sub <- (K[[1]]+1):(K[[1]]+K[[2]])
  }
  if (which == "ifr"){
    if (!describe(object, "ifr")) stop("no ifrequency equation")
    else sub <- (K[[1]]+K[[2]]+1):(K[[1]]+K[[2]]+K[[3]])
  }
  if (which == "other")
    if (describe(object,"corr")) sub <- (K[[1]]+K[[2]]+K[[3]]+1):(K[[1]]+K[[2]]+K[[3]]+2)
    else sub <- (K[[1]]+K[[2]]+K[[3]]+1) 
  if (which == "rho"){
    if (!describe(object, "corr")) stop("no corr coefficient estimated")
    else sub <- K[[1]]+K[[2]]+K[[3]]+2
  }
  if (which == "sigma") sub <- K[[1]]+K[[2]]+K[[3]]+1
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

logLik.mhurdle <- function(object,
                           which = c("all", "zero", "positive"),
                           naive = FALSE, ...){
  which <- match.arg(which)
  if (!naive){
    x <- object$logLik
    y <- attr(x, "y")
    x <- switch(which,
                all      = sum(x),
                zero     = sum(x[y == 0]),
                positive = sum(x[y > 0])
              )
    attr(x,"df") <- NULL
    attr(x,"y")  <- NULL
    result <- sum(x)
  }
  else{
    naive <- object$naive
    result <- switch(which,
                     all      = naive$logLik,
                     zero     = naive$logLik.part['zero'],
                     positive = naive$logLik.part['positive'])
  }
  result
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
  object$r.squared <- c(all      = r.squared(object, "all"),
                        zero     = r.squared(object, "zero"),
                        positive = r.squared(object, "positive"))
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

  cat("\nMc Fadden R^2 :\n")
  rs <- x$r.squared
  cat(paste(" - all      :", signif(rs[1], digits), "\n"))
  cat(paste(" - zero     :", signif(rs[2], digits), "\n"))
  cat(paste(" - positive :", signif(rs[3], digits), "\n"))

  if (FALSE){
  cat(paste("model r-squared: ",
            signif(r.squared(x,which="all",type="linear",dfcor=FALSE,r2pos = "rss"),
                   digits)))
  cat("\n")
  cat(paste("censored mechanism r-squared: ",
            signif(r.squared(x,which="zero",type="linear",dfcor=FALSE,r2pos = "rss"),
                   digits)))
  cat("\n")
  cat(paste("positive r-squared: ",
            signif(r.squared(x,which="positive",type="linear",dfcor=FALSE,r2pos = "rss"),
                   digits)))
  cat("\n")}
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
  for (i in 1:length(attr(call$formula, "rhs"))){
    # update.Formula returns "1 - 1" instead of 0 for empty parts
    zero <- paste(deparse(attr(call$formula, "rhs")[[i]])) == as.character("1 - 1")
    if (zero) attr(call$formula, "rhs")[[i]] <- 0
  }
  eval(call, parent.frame())
}

compute.fitted.mhurdle <- function(param, X, S, P, dist, corr){

  KS <- ifelse(is.null(S), 0, ncol(S))
  KP <- ifelse(is.null(P), 0, ncol(P))  
  KX <- ncol(X)
  sel <- KS > 0
  ifr <- KP > 0
  if (sel) betaS <- param[1:KS] else betaS <- NULL
  betaX <- param[(KS+1):(KS+KX)]
  if (ifr) betaP <- param[(KX+KS+1):(KX+KS+KP)] else betaP <- NULL
  sigma <- param[KX+KS+KP+1]
  if (corr) rho <- param[KX+KS+KP+2] else rho <- 0
  bX <- as.numeric(crossprod(t(X),betaX))
  PhiX <- pnorm(bX/sigma)
  phiX <- dnorm(bX/sigma)

  if (sel){
    bS <- as.numeric(crossprod(t(S),betaS))
    PhiS <- pnorm(bS)
    phiS <- dnorm(bS)
    Phib <- pbivnorm(bS,bX/sigma,rho)
    psi <- psy(bS,bX/sigma,rho)
    millsG <- psi/Phib
    Phib.grad <- attr(Phib,"gradient")
    attr(Phib, "gradient") <- NULL
    Phib.S <- Phib.grad$x1
    Phib.X <- Phib.grad$x2
    Phib.rho <- Phib.grad$rho
    Phib.rho.rho <- Phib.grad$rho.rho
  }
  else{
    bS <- 0
    Phib <- PhiX
    PhiS <- 1
    psi <- psy(bS,bX/sigma,rho)
    millsG <- psi/Phib
  }
  if (ifr){
    bP <- as.numeric(crossprod(t(P),betaP))
    PhiP <- pnorm(bP)
  }
  else{
    PhiP <- 1
  }
 
  if (dist == "l" ) prob.null <- 1 - PhiS*PhiP
  if (dist == "t" ) prob.null <- 1 - Phib*PhiP/PhiX
  if (dist == "n" ) prob.null <- 1 - Phib*PhiP
  if ((dist == "l") && sel){
    phiS <- dnorm(bS)
    esp.cond <-
      exp(bX+sigma^2*(1-rho^2))/PhiP*(PhiS+rho*sigma*phiS+
                                      (rho*sigma)^2/2*(PhiS-bS*phiS)+
                                      (rho*sigma)^3/6*(2+bS^2)*phiS+
                                      (rho*sigma)^4/12*(3*PhiS-bS*(3+bS^2)*phiS))/PhiS}
  if ((dist == "l") && !sel){
    esp.cond <- exp(bX+sigma^2*(1-rho^2))/PhiP}
  
  if (dist !="l")
    esp.cond <- bX/PhiP + sigma * millsG/PhiP^2
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
    dist <- ifelse(is.null(cl$dist), TRUE, cl$dist)
    corr <- ifelse(is.null(cl$corr), FALSE, cl$corr)
    m <- model.frame(formula(object), newdata)
    X <- model.matrix(formula(object),m,rhs=1)
    S <- model.matrix(formula(object),m,rhs=2)
    P <- model.matrix(formula(object),m,rhs=3)
    result <- compute.fitted.mhurdle(coef(object), X, S, P,  dist, corr)
  }
  result
}

fitted.mhurdle <- function(object, which = c("all", "zero", "positive"), ...){
  which <- match.arg(which)
  switch(which,
         all = cbind(
           "zero" = object$fitted.values[,1],
           "positive" = object$fitted.values[,2]
           ),
         zero = object$fitted.values[,1],
         positive = object$fitted.values[,2]
         )
}


r.squared <- function(object,
                      which = c("all", "zero", "positive")){
  lnL <- logLik(object, which = which, naive = FALSE)
  lnLo <- logLik(object, which = which, naive = TRUE)
  as.numeric(1 - lnL/lnLo)
}
