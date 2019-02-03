## nm.mhurdle : extract the names of a relevant subset of the coefficients
## sub.murdle : extract the indexes of a relevant subset of the coefficients
## coef.mhurdle
## vcov.mhurdle
## logLik.mhurdle
## print.mhurdle
## summary.mhurdle
## coef.summary.mhurdle
## print.summary.mhurdle
## fitted.mhurdle
## predict.mhurdle
## update.mhurdle
## rsq
## effects.mhurdle

nm.mhurdle <- function(object,
                       which = c("all", "h1", "h2", "h3", "sd", "h4", "corr", "tr", "pos"),
                       ...){
    coefnames <- object$coef.names
    which <- match.arg(which)
    K <- sapply(coefnames,length)
    if (which == "all"){
        h2.names <- paste("h2", coefnames$h2,sep = ".")
        h1.names <- h3.names <- h4.names <- NULL
        if (! is.null(coefnames$h1)) h1.names <- paste("h1", coefnames$h1,sep = ".")
        if (! is.null(coefnames$h3)) h3.names <- paste("h3", coefnames$h3,sep = ".")
        if (length(coefnames$sd) == 1) sd.names <- "sd"
        if (! is.null(coefnames$h4)) h4.names <- paste("h4", coefnames$h4,sep = ".")
        else sd.names <- paste("sd", coefnames$sd, sep = ".")
        corr.names <- coefnames$corr
        tr.names <- coefnames$tr
        mu.names <- coefnames$pos
        result <- c(h1.names, h2.names, h3.names, sd.names, h4.names, corr.names, tr.names, mu.names)
    }
    else{
        result <- coefnames[[which]]
        if (is.null(result)) stop(paste("no", which, "coefficient\n"))
    }
    result
}

sub.mhurdle <- function(object,
                        which = c("all", "h1", "h2", "h3", "sd", "h4", "corr", "tr", "pos"),
                        ...){
  # there is no need to check if the coefficient is relevant at it has
  # been checked previously by the nm.mhurdle function
    which <- match.arg(which)
    if ("mhurdle" %in% class(object)) K <- lapply(object$coef.names, length)
    else K <- lapply(object, length)
    if (which == "all")  sub <- 1:sum(Reduce("c", K))
    if (which == "h2")   sub <- (K$h1 + 1):(K$h1 + K$h2)
    if (which == "h1")   sub <- 1:K$h1
    if (which == "h3")   sub <- (K$h1 + K$h2 + 1):(K$h1 + K$h2 + K$h3)
    if (which == "sd")   sub <- (K$h1 + K$h2 + K$h3 + 1)
    if (which == "h4")   sub <- (K$h1 + K$h2 + K$h3 + 1 + 1):(K$h1 + K$h2 + K$h3 + 1 + K$h4)
    if (which == "corr") sub <- (K$h1 + K$h2 + K$h3 + 1 + K$h4 + 1) : (K$h1 + K$h2 + K$h3 + 1 + K$h4 + K$corr)
    if (which == "tr")   sub <- (K$h1 + K$h2 + K$h3 + 1 + K$h4 + K$corr + 1)
    if (which == "pos")  sub <- K$h1 + K$h2 + K$h3 + 1 + K$h4 + K$corr + K$tr + 1
    sub
}

coef.mhurdle <- function(object,
                         which = c("all", "h1", "h2", "h3", "h4", "sd", "corr", "tr", "pos"),
                      ...){
  which <- match.arg(which)
  nm <- nm.mhurdle(object, which)
  sub <- sub.mhurdle(object, which)
  result <- object$coefficients[sub]
  names(result) <- nm
  result
}

vcov.mhurdle <- function(object,
                         which = c("all", "h1", "h2", "h3", "h4", "sd", "corr", "tr", "pos"),
                      ...){
  which <- match.arg(which)
  nm <- nm.mhurdle(object, which)
  sub <- sub.mhurdle(object, which)
  result <- object$vcov
  result <- result[sub, sub]
  if (is.matrix(result)) rownames(result) <- colnames(result) <- nm
  else names(result) <- nm
  result
}

logLik.mhurdle <- function(object, naive = FALSE, ...){
    if (naive) result <- object$naive$logLik
    else{
        result <- object$logLik
        attrlogLik <- attributes(result)
        result <- sum(object$logLik)
        attributes(result) <- attrlogLik
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
  coefficients <- cbind(b,std.err,z,p)
  colnames(coefficients) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  object$coefficients <- coefficients
  object$r.squared <- c(coefdet = rsq(object, type = "coefdet"),
                        lratio  = rsq(object, type = "lratio"))
  class(object) <- c("summary.mhurdle","mhurdle")
  return(object)
}


coef.summary.mhurdle <- function(object,
                                 which = c("all", "h1", "h2", "h3", "sd", "corr", "tr", "pos"),
                                 ...){
  which <- match.arg(which)
  sub <- sub.mhurdle(object, which)
  nm <- nm.mhurdle(object, which)
  result <- object$coefficients
  if (!is.null(sub)) result <- result[sub, , drop = FALSE]
  rownames(result) <- nm
  result
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
  printCoefmat(x$coefficients, digits = digits)
  cat("\n")
  df <- attr(x$logLik, "df")

  cat(paste("Log-Likelihood: ",
            signif(logLik(x), digits),
            " on ",df," Df\n",sep=""))

  cat("\nR^2 :\n")
  rs <- x$r.squared
  cat(paste(" Coefficient of determination :", signif(rs['coefdet'], digits), "\n"))
  cat(paste(" Likelihood ratio index       :", signif(rs['lratio'], digits), "\n"))
  invisible(x)
}

fitted.mhurdle <- function(object, which = c("all", "zero", "positive"), mean = FALSE, ...){
  which <- match.arg(which)
  res <- switch(which,
                all      = object$fitted.values,
                zero = object$fitted.values[, 1],
                positive = object$fitted.values[, 2]
                )
  if (mean){
      if (is.matrix(res)) res <- apply(res, 2, mean)
      else res <- mean(res)
  }
  res
}

predict.mhurdle <- function(object, newdata = NULL, ...){
    geomean <- attr(object$model, "geomean")
    if (is.null(newdata)){
        result <- fitted(object, ...)
    }
    else{
        cl <- object$call
        h2 <- ifelse(is.null(cl$h2), FALSE, TRUE)
        dist <- if (is.null(cl$dist)) ifelse(h2, "ln2", "ln")
        corr <- ifelse(is.null(cl$corr), FALSE, cl$corr)
        robust <- FALSE
        m <- model.frame(formula(object), newdata)
        X1 <- model.matrix(formula(object), m, rhs = 1)
        X2 <- model.matrix(formula(object), m, rhs = 2)
        if (length(formula(object))[2] > 2) X3 <- model.matrix(formula(object), m, rhs = 3) else X3 <- NULL
        if (length(formula(object))[2] == 4) X4 <- model.matrix(formula(object), m, rhs = 4) else X4 <- NULL
        y <- model.response(m)
        if (! is.null(geomean)) attr(y, "geomean") <- geomean
        if (length(X1) == 0) X1 <- NULL
        result <- attr(mhurdle.lnl(coef(object), X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y,
                                   gradient = FALSE, fitted = TRUE, robust = FALSE,
                                   dist = dist, corr = corr), "fitted")
    }
    result
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



#' R squared and pseudo R squared
#' 
#' This function computes the R squared for multiple hurdle models. The measure
#' is a pseudo coefficient of determination or may be based on the likelihood.
#' 
#' 
#' @param object an object of class \code{"mhurdle"},
#' @param type one of \code{"coefdet"} or \code{"lratio"} to select a pseudo
#' coefficient of correlation or a Mc Fadden like measure based on the
#' likelihood function,
#' @param adj if \code{TRUE} a correction for the degrees of freedom is
#' performed,
#' @param r2pos only for pseudo coefficient of determination, should the
#' positive part of the R squared be computed using the residual sum of squares
#' (\code{"rss"}), the explained sum of squares (\code{"ess"}) or the
#' coefficient of correlation between the fitted values and the response
#' (\code{cor}).
#' @return a numerical value
#' @references
#' 
#' McFadden D (1974). The Measurement of Urban Travel Demand. Journal of Public
#' Economics, 3, 303-328.
#' @keywords htest
#' @examples
#' 
#' data("Interview", package = "mhurdle")
#' # independent double hurdle model
#' idhm <- mhurdle(vacations ~ car + size | linc + linc2 | 0, Interview,
#'               dist = "ln", h2 = TRUE, method = "bfgs")
#' rsq(idhm, type = "lratio")
#' rsq(idhm, type = "coefdet", r2pos = "rss")
#' 
rsq <- function(object,
                type = c("coefdet", "lratio"),
                adj = FALSE,
                r2pos=c("rss","ess","cor")){
    
    type <- match.arg(type)
    r2pos <- match.arg(r2pos)
    
    K1 <- length(object$coef.names$h1)
    K2 <- length(object$coef.names$h2)
    K3 <- length(object$coef.names$h3)
    
    y <- model.response(model.frame(object))
    n <- length(y)
    no <- sum(y == 0)
    po <- mean(y == 0)
    pp <- 1 - po

    K <- length(coef(object))
    Ko <- length(object$naive$coefficients)
    
    if (type == "lratio"){
        ## print(logLik(object))
        ## print(logLik(object, naive = TRUE))
        ## print(c(K, Ko))
              
        if (!adj) R2 <- 1 - logLik(object) / logLik(object, naive = TRUE)
        else R2 <- 1 - (logLik(object) - K) / (logLik(object, naive = TRUE) - Ko)
        R2 <- as.numeric(R2)
    }
    if (type == "coefdet"){
        ym <- mean(y)
        yf <- fitted(object, "positive") * (1 - fitted(object, "zero"))
        R2 <- switch(r2pos,
                     ess = ifelse(adj,
                         sum( (yf - ym) ^ 2) / sum( (y - ym) ^ 2) * (n - K) / (n - Ko),
                         sum( (yf - ym) ^ 2) / sum( (y - ym) ^ 2)),
                     rss = ifelse(adj,
                         1 - (n - Ko) / (n - K) * sum( (y - yf) ^ 2) / sum( (y - ym) ^ 2),
                         1 - sum( (y - yf) ^ 2) / sum( (y - ym) ^ 2)),
                     cor = ifelse(adj,
                         stop("no adjusted R2 using the correlation formula"),
                         cor(y, yf) ^ 2
                         )
                     )
    }
    R2
}


nobs.mhurdle <- function(object, which = c("all", "null", "positive"), ...){
    y <- model.response(model.frame(object))
    which <- match.arg(which)
    switch(which,
           all = length(y),
           null = sum(y == 0),
           positive = sum(y > 0))
}


effects.mhurdle <- function(object, covariate = NULL, data = NULL, reflevel = NULL, mean = FALSE, ...){
    if (is.null(covariate)) stop("the name of a covariate should be indicated")
    if (is.null(data)) odata <- eval(object$call$data) else odata <- data
    eps <- 1E-04
    ndata <- odata
    thecov <- ndata[[covariate]]
    if (is.numeric(thecov)){
        step <- ifelse(is.integer(thecov), 1, eps)
        ndata[[covariate]] <- ndata[[covariate]] + step
        ofitted <- predict(object, odata)
        nfitted <- predict(object, ndata)
        mfx <- (nfitted - ofitted) / step
        mfx[abs(mfx) < 1E-08] <- 0
        if (mean) mfx <- apply(mfx, 2, mean)
    }
    if (is.factor(thecov)){
        levs <- levels(thecov)
        if (is.null(reflevel)) reflevel <- levs[1]
        else{
            if (! reflevel %in% levs) stop("undefined level")
        }
        nx <- vector(mode = "list", length = length(levs))
        nx <- lapply(levs, function(d){
            ndata[[covariate]] <- factor(d, levels = levels(thecov))
            predict(object, ndata)
        }
        )
        zero <- sapply(nx, function(x) x[, 1])
        pos <- sapply(nx, function(x) x[, 2])
        colnames(zero) <- colnames(pos) <- levs
        zero <- zero - zero[, reflevel]
        pos <- pos - pos[, reflevel]
        zero <- zero[, - which(levs == reflevel)]
        pos <- pos[, - which(levs == reflevel)]
        mfx <- list(zero = zero, pos = pos)
        if (mean){
            mfx <- sapply(mfx, apply, 2, mean)
        }
    }
    mfx
}


getindex <- function(X1, X2, X3, X4, corr, dist, which){
    K1 <- ifelse(is.null(X1), 0, ncol(X1))
    K2 <- ncol(X2)
    K3 <- ifelse(is.null(X3), 0, ncol(X3))
    K4 <- ifelse(is.null(X4), 0, ncol(X4))
    cumul <- 0

    if (which == "h1"){
        if (K1 == 0) return(numeric(0)) else return(1:K1)
    }
    cumul <- cumul + K1

    if (which == "h2") return( (cumul + 1):(cumul + K2) )
    cumul <- cumul + K2

    if (which == "h3"){
        if (K3 == 0) return(numeric(0)) else return( (cumul + 1) : (cumul + K3) )
    }
    cumul <- cumul + K3

    if (which == "sd") return(cumul + 1)
    cumul <- cumul + 1

    if (which == "h4"){
        if (K4 == 0) return(numeric(0)) else return( (cumul + 1) : (cumul + K4) )
    }
    cumul <- cumul + K4

    if (corr){
        h1 <- K1 > 0
        h3 <- K3 > 0
        if (! h1 & ! h3) return(numeric(0))
        else{
            if (h1 + h3 == 2){
                if (which == "corr") return( (cumul + 1) : (cumul + 3) )
                cumul <- cumul + 3
            }
            else{
                if (which == "corr") return( (cumul + 1) )
                cumul <- cumul + 1
            }
        }
    }
    else{
        if (which == "corr") return(numeric(0))
    }

    if (dist %in% c("ihs", "bc", "bc2")){
        if (which == "tr") return(cumul + 1)
        cumul <- cumul + 1
    }
    else{
        if (which == "tr") return(numeric(0))
    }
    
    if (dist %in% c("ln2", "bc2")){
        if (which == "pos") return(cumul + 1)
        cumul <- cumul + 1
    }
    else{
        if (which == "pos") return(numeric(0))
    }
}


## extract.maxLik <- function (model, include.nobs = TRUE, ...){
##     s <- summary(model, ...)
##     names <- rownames(s$estimate)
##     class(names) <- "character"
##     co <- s$estimate[, 1]
##     se <- s$estimate[, 2]
##     pval <- s$estimate[, 4]
##     class(co) <- class(se) <- class(pval) <- "numeric"
##     n <- nrow(model$gradientObs)
##     lik <- logLik(model)
##     gof <- numeric()
##     gof.names <- character()
##     gof.decimal <- logical()
##     gof <- c(gof, n, lik)
##     gof.names <- c(gof.names, "Num. obs.", "Log Likelihood")
##     gof.decimal <- c(gof.decimal, FALSE, TRUE)
##     tr <- createTexreg(coef.names = names, coef = co, se = se, pvalues = pval,
##                        gof.names = gof.names, gof = gof, gof.decimal = gof.decimal)
##     return(tr)
## }

## setMethod("extract", signature = className("maxLik", "maxLik"), definition = extract.maxLik)

extract.mhurdle <- function (model, include.nobs = TRUE, ...){
    s <- summary(model, ...)
    names <- rownames(s$coefficients)
    class(names) <- "character"
    co <- s$coefficients[, 1]
    se <- s$coefficients[, 2]
    pval <- s$coefficients[, 4]
    class(co) <- class(se) <- class(pval) <- "numeric"
    n <- nobs(model)
    lik <- logLik(model)
    coefdet <- s$r.squared['coefdet']
    lratio <- s$r.squared['lratio']
    R2dich <- model$R2[1]
    R2pos <- model$R2[2]
    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    gof <- c(gof, n, lik, coefdet, lratio, R2dich, R2pos)
#    gof.names <- c(gof.names, "Num. obs.", "Log Likelihood", "$R^2$", "McFadden $R^2$")
    gof.names <- c(gof.names, "Num. obs.", "Log Likelihood", "$R^2$", "McFadden $R^2$", "$R^2 (y = 0)$", "$R^2 (y > 0)$")    
    gof.decimal <- c(gof.decimal, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
    tr <- createTexreg(coef.names = names, coef = co, se = se, pvalues = pval,
                       gof.names = gof.names, gof = gof, gof.decimal = gof.decimal)
    return(tr)
}

setMethod("extract", signature = className("mhurdle", "mhurdle"), definition = extract.mhurdle)

