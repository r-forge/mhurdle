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
                       which = c("all", "h1", "h3", "h2", "other", "sigma", "rho"),
                       ...){

  which <- match.arg(which)
  K <- sapply(object$coef.names,length)
  if (which == "all"){
    h2.names <- paste("h2", object$coef.names$h2,sep = ".")
    h1.names <- h3.names <- NULL
    if (describe(object, "h1")) h1.names <- paste("h1", object$coef.names$h1,sep = ".")
    if (describe(object, "h3")) h3.names <- paste("h3", object$coef.names$h3,sep = ".")

    result <- c(h1.names, h2.names, h3.names, object$coef.names$other)
  }
  if (which == "h1") result <- object$coef.names[[1]]
  if (which == "h2") result <- object$coef.names[[2]]
  if (which == "h3") result <- object$coef.names[[3]]
  if (which == "other")
    if (!is.null(describe(object, "corr"))) {result <- c("sigma","rho")}
    else {result <- "sigma"}
  if (which == "rho") result <- "rho"
  if (which == "sigma") result <- "sigma"
  result
}

sub.mhurdle <- function(object,
                        which = c("all", "h1", "h3", "h2", "other", "sigma", "rho"),
                        ...){

  which <- match.arg(which)
  K <- sapply(object$coef.names,length)
  if (which == "all") sub <- NULL
  if (which == "h2") sub <- (K[[1]] + 1):(K[[1]] + K[[2]])
  if (which == "h1"){
    if (!describe(object, "h1")) stop("no sel equation")
    else sub <- 1:K[[1]]
  }
  if (which == "h3"){
    if (!describe(object, "h3")) stop("no ifrequency equation")
    else sub <- (K[[1]]+K[[2]]+1):(K[[1]]+K[[2]]+K[[3]])
  }
  if (which == "other")
    if (!is.null(describe(object,"corr"))) sub <- (K[[1]]+K[[2]]+K[[3]]+1):(K[[1]]+K[[2]]+K[[3]]+2)
    else sub <- (K[[1]]+K[[2]]+K[[3]]+1) 
  if (which == "rho"){
    if (is.null(describe(object, "corr"))) stop("no corr coefficient estimated")
    else sub <- K[[1]]+K[[2]]+K[[3]]+2
  }
  if (which == "sigma") sub <- K[[1]]+K[[2]]+K[[3]]+1
  sub
}


coef.mhurdle <- function(object,
                      which = c("all", "h1", "h2", "h3", "other", "sigma", "rho"),
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
                      which = c("all", "h1", "h2", "h3", "other", "sigma", "rho"),
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

logLik.mhurdle <- function(object, naive = FALSE, ...){
  if (!naive){
    x <- object$logLik
    y <- attr(x, "y")
    x <- sum(x)
    attr(x,"df") <- NULL
    attr(x,"y")  <- NULL
    result <- sum(x)
  }
  else{
    naive <- object$naive
    result <- naive$logLik
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

coef.summary.mhurdle <- function(object,
                                 which = c("all", "h1", "h3", "h2", "other", "sigma", "rho"),
                                 ...){
  which <- match.arg(which)
  sub <- sub.mhurdle(object, which)
  nm <- nm.mhurdle(object, which)
  result <- object$CoefTable
  if (!is.null(sub)) result <- result[sub, , drop=FALSE]
  rownames(result) <- nm
  result
}

summary.mhurdle <- function (object,...){
  b <- coef(object)
  std.err <- sqrt(diag(vcov(object)))
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  object$CoefTable <- CoefTable
  object$r.squared <- c(regression = r.squared(object, type = "regression"),
                        mcfadden = r.squared(object, type = "mcfadden"))
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
#    cat(paste("implied value :", round(rhotest$value,3)," ",sep=""))
    z <- 
    cat(paste("score test : z = ", round(rhotest$test$statistic,3),
              " (p.value = ",round(rhotest$test$p.value,3),")\n",sep=""))
  }

  cat("\nR^2 :\n")
  rs <- x$r.squared
  cat(paste(" McFadden   :", signif(rs['mcfadden'], digits), "\n"))
  cat(paste(" Regression :", signif(rs['regression'], digits), "\n"))
  invisible(x)
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

compute.fitted.mhurdle <- function(param, X1, X2, X3, dist, corr){

  h1 <- !is.null(X1) ;  K1 <- ifelse(is.null(X1), 0, ncol(X1))
  h3 <- !is.null(X3) ;  K3 <- ifelse(is.null(X3), 0, ncol(X3))
  
  K2 <- ncol(X2)
  beta2 <- param[(K1 + 1) : (K1 + K2)]

  if (h1){
    beta1 <- param[1:K1]
    bX1 <- as.numeric(crossprod(t(X1), beta1))
    Phi1 <- pnorm(bX1) ; phi1 <- dnorm(bX1)
  }
  else{
    bX1 <- beta1 <- NULL
    Phi1 <- 1 ; phi1 <- 0;
  }

  if (h3){
    beta3 <- param[(K1 + K2 + 1):(K1 + K2 + K3)]
    bX3 <- as.numeric(crossprod(t(X3), beta3))
    Phi3 <- pnorm(bX3) ; phi3 <- dnorm(bX3)
  }
  else{
    bX3 <- beta3 <- NULL
    Phi3 <- 1 ; phi3 <- 0
  }
  
  sigma <- param[K1 + K2 + K3 + 1]
  bX2 <- as.numeric(crossprod(t(X2),beta2))
  Phi2 <- pnorm(bX2 / sigma)
  phi2 <- dnorm(bX2 / sigma)

  ATAN <- FALSE
  rho1 <- rho3 <- 0
  if (!is.null(corr)){
    if (ATAN){
      if (corr == 'h1') rho1 <- atan(param[K1 + K2 + K3 + 2]) * 2 / pi
      if (corr == 'h3') rho3 <- atan(param[K1 + K2 + K3 + 2]) * 2 / pi
    }
    else{
      if (corr == 'h1') rho1 <- param[K1 + K2 + K3 + 2]
      if (corr == 'h3') rho3 <- param[K1 + K2 + K3 + 2]
      if (rho1 < -1) rho1 <- - 0.99
      if (rho1 >  1) rho1 <-   0.99
      if (rho3 < -1) rho3 <- - 0.99
      if (rho3 >  1) rho3 <-   0.99
    }      
  }
  Phi12 <- mypbivnorm(bX1, bX2 / sigma, rho1)
  Phi23 <- mypbivnorm(bX2 / sigma, bX3, rho3)
  
  prob.null <- switch(dist,
                      "t" = 1 - Phi12$f * Phi23$f / Phi2 ^ 2,
                      "n" = 1 - Phi12$f * Phi23$f / Phi2,
                      "l" = 1 - Phi3 * Phi1
                      )
  
  # (2i) et (2d)
  if ( (dist == "l") && h1 && !h3){
    if (is.null(corr)) esp.cond <- exp(bX2 + 0.5 * sigma^2 * (1 - rho1 ^ 2) )
    else esp.cond <- exp(bX2 + 0.5 * sigma^2 * (1 - rho1 ^ 2) ) * psylog(bX1, rho1 * sigma) / Phi1
  }

  # (2ti) et (2td)
  if ( (dist == "t") && h1 && !h3){
    phi1 <- dnorm(bX1)
    if (is.null(corr)) esp.cond <- bX2 + sigma * mills(bX2 / sigma)
    else esp.cond <- bX2 + sigma * psy(bX1, bX2 / sigma, rho1) / Phi12$f
  }

  # (3)
  if ( (dist == "n") && !h1 && !h3 ) esp.cond <- bX2 + sigma * phi2 / Phi2
  
  # (4)
  if ( (dist == "l") && !h1 && h3){
    if (is.null(corr)) esp.cond <- exp(bX2 + 0.5 * sigma^2 * (1 - rho3 ^ 2) ) / Phi3
    else esp.cond <- exp(bX2 + 0.5 * sigma^2 * (1 - rho3 ^ 2)) *  psylog(bX3, rho3 * sigma) / Phi3
  }
    
  # (4t)
  if ( (dist == "t") && !h1 && h3){
    if (is.null(corr)) esp.cond <- bX2 + sigma * mills(bX2 / sigma)
    else esp.cond <- (bX2 + sigma * psy(bX3, bX2 / sigma, rho3) / Phi23$f) / Phi3
  }
  
  # (5i) et (5d)
  if ( (dist == "n") && h1 && !h3){
    if (is.null(corr)) esp.cond <- bX2 + sigma * mills(bX2 / sigma)
    else esp.cond <- bX2 + sigma * psy(bX1, bX2 / sigma, rho1) / Phi12$f
  }
    
  # (6)
  if ( (dist == "l") && h1 && h3){
    if (is.null(corr)) esp.cond <- exp(bX2 + 0.5 * sigma^2)
    else{
      if (corr == "h1") esp.cond <- exp(bX2 + 0.5 * sigma^2 * (1 - rho1 ^ 2)) * psylog(bX1, rho1 * sigma) / Phi1
      if (corr == "h3") esp.cond <- exp(bX2 + 0.5 * sigma^2 * (1 - rho1 ^ 2)) * psylog(bX3, rho3 * sigma) / Phi3
    }
  }

  # (6t)
  if ( (dist == "t") && h1 && h3){
    esp.cond <- NA
  }

  # (7)
  if ( (dist == "n") && !h1 && h3){
    if (is.null(corr)) esp.cond <- (bX2 + sigma * mills(bX2 / sigma)) / Phi3
    else esp.cond <- (bX2 + sigma * psy(bX3, bX2 / sigma, rho3) / Phi23$f) / Phi3
  }
  
  # (8)
  if ( (dist == "n") && h1 && h3){
    if (is.null(corr)) esp.cond <- (bX2 + sigma * mills(bX2 / sigma)) / Phi3
    else{
      if (corr == "h1") esp.cond <- (bX2 + sigma * psy(bX1, bX2 / sigma, rho1) / Phi12$f) / Phi3
      if (corr == "h3") esp.cond <- (bX2 + sigma * psy(bX2 / sigma, bX3, rho3) / Phi23$f) / Phi3
    }          
  }
  result <- cbind(zero     = prob.null,
                  positive = esp.cond)
  result
}


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

psylog <- function(x, rho){
  phi <- dnorm(x)
  Phi <- pnorm(x)
  (Phi + rho * phi +
   rho ^ 2 / 2 * (Phi - x * phi) +
   rho ^ 3 / 6 * (2 + x ^ 2) * phi +
   rho ^ 4 / 12 * (3 * Phi - x * (3 + x ^ 2) * phi))
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
    X1 <- model.matrix(formula(object), m, rhs = 1)
    X2 <- model.matrix(formula(object), m, rhs = 2)
    X3 <- model.matrix(formula(object), m, rhs = 3)
    if (length(X1) == 0) X1 <- NULL
    if (length(X3) == 0) X3 <- NULL
    result <- compute.fitted.mhurdle(coef(object), X1, X2, X3,  dist, corr)
  }
  result
}

fitted.mhurdle <- function(object, which = c("all", "zero", "positive"), ...){
  which <- match.arg(which)
  switch(which,
         all      = cbind(
           "zero" = object$fitted.values[,1],
           "positive" = object$fitted.values[,2]
           ),
         zero = object$fitted.values[,1],
         positive = object$fitted.values[,2]
         )
}

print.est.stat <- function(x, ...){
  et <- x$elaps.time[3]
  i <- x$nb.iter[1]
  halton <- x$halton
  method <- x$method
  if (!is.null(x$type) && x$type != "simple"){
    R <- x$nb.draws
    cat(paste("Simulated maximum likelihood with", R, "draws\n"))
  }
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  cat(paste(method, "method\n"))
  tstr <- paste(h, "h:", m, "m:", s, "s", sep="")
  cat(paste(i,"iterations,",tstr,"\n"))
  if (!is.null(halton)) cat("Halton's sequences used\n")
  if (!is.null(x$eps)) cat(paste("g'(-H)^-1g =", sprintf("%5.3G", as.numeric(x$eps)),"\n"))
  if (is.numeric(x$code)){
    msg <- switch(x$code,
                  "1" = "gradient close to zero",
                  "2" = "successive fonction values within tolerance limits",
                  "3" = "last step couldn't find higher value",
                  "4" = "iteration limit exceeded"
                  )
    cat(paste(msg, "\n"))
  }
  else cat(paste(x$code, "\n"))
}

