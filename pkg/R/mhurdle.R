mhurdle <- function(formula, data, subset, weights, na.action,
                 start = NULL, dist = c("l","t","n"), corr = NULL, ...){
  dots <- list(...)
  oldoptions <- options(warn=-1)
  on.exit(options(oldoptions))
  cl <- match.call()
  posT <- as.list(cl) == "T"
  posF <- as.list(cl) == "F"
  cl[posT] <- TRUE
  cl[posF] <- FALSE
  cl.save <- cl
  dist <- match.arg(dist)

  ##
  # compute the model.frame and the model.matrix
  ##

  if (!inherits(formula, "Formula")) formula <- Formula(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action","weights"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf$formula <- formula
  mf <- eval(mf, parent.frame())
  X1 <- model.matrix(formula, data=mf , rhs = 1)
  X2 <- model.matrix(formula, data=mf , rhs = 2)
  X3 <- model.matrix(formula, data=mf , rhs = 3)
  y <- model.response(mf)
  n <- length(y)
  if (length(X1) == 0) X1 <- NULL
  if (length(X3) == 0) X3 <- NULL
  if (length(X2) == 0) stop("the second hurdle (consumption equation) is mandatory")
  h1 <- !is.null(X1)
  h3 <- !is.null(X3)

  if (!is.null(corr)){
    if (h1 & h3) if (!(corr %in% c('h1', 'h3'))) stop('irrelevant value for corr')
    if (h1 & !h3) if (corr != 'h1') stop('irrelevant value for corr')
    if (!h1 & h3) if (corr != 'h3') stop('irrelevant value for corr')
  }
  
  ##
  # compute the "naive" model
  ##

  Pnull <- mean(y == 0)
  if (dist != "l"){
    Ec <- mean(y[y > 0])
    Vc <- var(y[y > 0])}
  else{
    Ec <- mean(log(y[y > 0]))
    Vc <- var(log(y[y > 0]))
  }

  start.naive <- c(rep(0.1, 1 + h1 + h3), 1)
  if (!is.null(corr)) start.naive <- c(start.naive, 0)
  moments <- c(Pnull, Ec, Vc)
  naive <- maxLik(lnl.naive, start = start.naive,
                  dist = dist, moments = moments,
                  h1 = h1, h3 = h3);
  coef.naive <- naive$est
  logLik.naive <- naive$max * n
  naive <- list(coefficients = coef.naive, logLik = logLik.naive, code = naive$code)

  #  Special cases where the models can be estimated directly, without
  #  relevant starting values
  
  # model (1l 1t), no censure, estimate the log-linear or truncated
  # model
  if (!h1 && !h3 && dist != "n"){
    if (dist == "l") result <- lm( log(y) ~ X2 - 1)
    if (dist == "t") result <- truncreg(y ~ X2 - 1)
    return(result)
  }
  # models (2li, 2ti) h1 without corr, can be estimated simply in two
  # parts using fit.simple.mhurdle() used as starting values for (3d,
  # 4d)
  if (h1 && !h3 &&  dist != "n" && is.null(corr)){
    result <- fit.simple.mhurdle(X1, X2, y, dist = dist)
    result$naive <- naive
    result$call <- cl.save
    result$model <- mf
    result$formula <- formula
    return(result)
  }

  # Compute the starting values if not provided
  if (is.null(start)) start <- start.mhurdle(X1, X2, X3, y, dist, corr)
  # Fit the model
  result <- mhurdle.fit(start, X1, X2, X3, y,
                        gradient = TRUE, fit = FALSE,
                        dist = dist, corr = corr, ...)
  result$naive <- naive
  result$call <- cl.save
  result$formula <- formula
  names(result$coefficients) <- colnames(result$vcov) <-
    rownames(result$vcov) <- nm.mhurdle(result)
  result$model <- mf
  result
}

mhurdle.fit <- function(start, X1, X2, X3, y, gradient = FALSE, fit = FALSE,
                        dist = c("l","n","t"), corr = NULL, ...){

  start.time <- proc.time()
  if (!is.null(corr)) other.coef <- c("sigma", "rho")
  else other.coef <- c("sigma")
  K1 <- ifelse(is.null(X1), 0, ncol(X1))
  K2 <- ncol(X2)
  K3 <- ifelse(is.null(X3), 0, ncol(X3))  
  h1 <- K1 > 0
  h3 <- K3 > 0

  f <- function(param) mhurdle.lnl(param, X1 = X1, X2 = X2, X3 = X3, y = y,
                              gradient = TRUE, fit = FALSE,
                              dist = dist, corr = corr)

  check.gradient <- FALSE
  if (check.gradient){
    ngrad <- c()
    oparam <- start
    fo <- f(start)
    agrad <- apply(attr(fo, "gradient"), 2, sum)
    eps <- 1E-5
    for (i in 1:length(start)){
      oparam[i] <- oparam[i] + eps
      ngrad <- c(ngrad, sum((as.numeric(f(oparam)) - fo) / eps))
      oparam <- start
    }
    print(cbind(start, agrad, ngrad))
  }
  maxl <- maxLik(f, start = start,...)
  coefficients <- maxl$estimate
  logLik <- f(coefficients)
  gradi <- attr(logLik, "gradi")
  attr(logLik,"df") <- length(coefficients)
  attr(logLik, "y") <- y
  attr(logLik, "gradi") <- NULL
  hessian <- maxl$hessian
  fitted <- compute.fitted.mhurdle(coefficients, X1, X2, X3, dist, corr)
  convergence.OK <- maxl$code <= 2
  grad.conv <- maxl$gradient
  elaps.time <- proc.time() - start.time
  nb.iter <- maxl$iterations
  eps <- grad.conv%*%solve(-maxl$hessian)%*%grad.conv
  est.stat <- list(elaps.time = elaps.time,
                   nb.iter = nb.iter,
                   eps = eps,
                   method = maxl$type,
                   message = maxl$message
                   )
  class(est.stat) <- "est.stat"
  coef.names <- list(h1   = colnames(X1),
                     h2   = colnames(X2),
                     h3   = colnames(X3),
                     other = other.coef)
  result <- list(coefficients  = coefficients,
                 vcov          = - solve(maxl$hessian),
                 fitted.values = fitted,
                 logLik        = logLik,
                 gradient      = gradi,
                 formula       = NULL,
                 model         = NULL,
                 coef.names    = coef.names,
                 call          = NULL,
                 est.stat      = est.stat,
                 naive         = NULL
                 )
  if (ncol(X2) > 1) class(result) <- c("mhurdle","maxLik")
  result
}

describe <- function(x, which){
  if (which == "h1") result <- attr(x$formula, "rhs")[[1]] != 0
  if (which == "h2") result <- attr(x$formula, "rhs")[[2]] != 0
  if (which == "h3") result <- attr(x$formula, "rhs")[[3]] != 0
  if (which == "corr") if(!is.null(x$call$corr)) result <- x$call$corr else result <- NULL
  result
}
 
# Compute the starting values of the model

start.mhurdle <- function(X1, X2, X3, y, dist, corr){
  h1 <- !is.null(X1)
  h3 <- !is.null(X3)
  # for models (2ld, 2td), estimates of models (2li, 2ti) which can be
  # estimated simply in two parts using fit.simple.mhurdle() are used
  # as starting values
  if (h1 && !h3 &&  dist != "n"){
    result <- fit.simple.mhurdle(X1, X2, y, dist = dist)
    start <- c(result$coefficients, 0.1)
  }

  # model (3) tobit : use linear model as starting values
  if (!h1 && !h3 &&  dist == "n"){
    lin <- lm(y ~ X2 - 1)
    start <- c(coef(lin), summary(lin)$sigma)
  }
  
  # model (4, 7) h3 whithout h1
  if (!h1 && h3){
    probit <- glm( (y != 0) ~ X3 - 1, family = binomial(link = 'probit'))
    bX3 <- as.numeric(crossprod(coef(probit), t(X3)))
    Phi3 <- pnorm(bX3)
    yPP <- y * Phi3
    lin <- switch(dist,
                  "l" = lm(log(yPP)  ~ X2 - 1, subset = y != 0),
                  "n" = lm(yPP       ~ X2 - 1, subset = y != 0),
                  "t" = truncreg(yPP ~ X2 - 1, subset = y != 0)
                  )
    if(dist != "t") start <- c(coef(lin), coef(probit), summary(lin)$sigma)
    else start <- c(coef(lin)[- length(coef(lin))], coef(probit), coef(lin)[length(coef(lin))])
    if (!is.null(corr)) start <- c(start, 0.1)
  }

  # model (5), double hurdle use model (3i) as starting values
  if (h1 && !h3 && dist == "n"){
    result <- fit.simple.mhurdle(X1, X2, y, dist = dist)
    start <- result$coefficients
    if (!is.null(corr)) start <- c(start, 0.1)
  }
  
  # model (6 and 8)
  if (h1 && h3){
    probit.h3 <- glm( (y != 0) ~ X3 - 1, family = binomial(link = 'probit'))
    probit.h1 <- glm( (y != 0) ~ X1 - 1, family = binomial(link = 'probit'))
    beta3 <- coef(probit.h3)
    beta1 <- coef(probit.h1)
    bX3 <- as.numeric(crossprod(beta3, t(X3)))
    Phi3 <- pnorm(bX3)
    bX1 <- as.numeric(crossprod(beta1, t(X1)))
    P0 <- mean(y > 0)
    yPP <- y * Phi3
    lin <- switch(dist,
                  "l" = lm(log(yPP)  ~ X2 - 1, subset = y!= 0),
                  "n" = lm(yPP       ~ X2 - 1, subset = y!= 0),
                  "t" = truncreg(yPP ~ X2 - 1, subset = y!= 0)
                  )
    if(dist != "t")
      start <- c(beta1, coef(lin), beta3, summary(lin)$sigma)
    else
      start <- c(beta1, coef(lin)[-length(coef(lin))],
                 beta3, coef(lin)[length(coef(lin))])
    if (!is.null(corr)) start <- c(start, 0.1)
  }
  start
}


  ## V <- 1985.696
  ## P1 <- 0.783
  ## yb <- 19.553

  ## za <- function(b2){
  ##   s <- sqrt( V + (yb - b2)^2 + (yb - b2) )
  ##   (yb-b2)/s - dnorm(b2/s) / pnorm(b2/s)
  ## }
