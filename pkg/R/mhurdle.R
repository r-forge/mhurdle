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
  X <- model.matrix(formula, data=mf , rhs = 2)
  S <- model.matrix(formula, data=mf , rhs = 1)
  P <- model.matrix(formula, data=mf , rhs = 3)
  y <- model.response(mf)
  n <- length(y)
  if (length(S) == 0) S <- NULL
  if (length(P) == 0) P <- NULL
  sel <- !is.null(S)
  ifr <- !is.null(P)

  if (!is.null(corr)){
    if (sel & ifr) if (!(corr %in% c('sel', 'ifr'))) stop('irrelevant value for corr')
    if (sel & !ifr) if (corr != 'sel') stop('irrelevant value for corr')
    if (!sel & ifr) if (corr != 'ifr') stop('irrelevant value for corr')
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

  start.naive <- c(rep(0.1, 1 + sel + ifr), 1)
  moments <- c(Pnull, Ec, Vc)
  naive <- maxLik(lnl.naive, start = start.naive,
                  dist = dist, moments = moments,
                  sel = sel, ifr = ifr)
  coef.naive <- naive$est
  logLik.naive <- naive$max * n
  logLik.naive.zero <- lnl.naive(coef.naive, dist = dist,
                                 moments = moments,
                                 sel = sel, ifr = ifr,
                                 which = "zero")  * n
  logLik.naive.pos <- lnl.naive(coef.naive, dist = dist,
                                moments = moments,
                                sel = sel, ifr = ifr,
                                which = "positive")  * n
  naive <- list(coefficients = coef.naive,
                logLik = logLik.naive,
                logLik.part = c(zero = logLik.naive.zero, positive = logLik.naive.pos)
                )
  ####
  #  Special cases where the models can be estimated directly, without
  #  relevant starting values
  ####
  
  # model (1, 2), no censure, estimate the log-linear or truncated
  # model
  if (!sel && !ifr && dist != "n"){
    if (dist == "l") result <- lm(log(y)~ X - 1)
    if (dist == "t") result <- truncreg(y ~ X - 1)
    return(result)
  }
  # models (3i, 4i) sel without corr, can be estimated simply in two
  # parts using fit.simple.mhurdle()
  # used as starting values for (3d, 4d)
  if (sel && !ifr &&  dist != "n" && is.null(corr)){
    result <- fit.simple.mhurdle(X, S, y, dist = dist)
    result$naive <- naive
    result$call <- cl.save
    result$model <- mf
    result$formula <- formula
    return(result)
  }
  
  ####
  #  Compute the starting values if not provided
  ####

  if (is.null(start)) start <- start.mhurdle(X, S, P, y, dist, corr)

  ####
  # Fit the model
  ####

  result <- mhurdle.fit(start, X, S, P, y,
                        gradient = TRUE, fit = FALSE, score = FALSE,
                        dist = dist, corr = corr, ...)
  result$naive <- naive
  result$call <- cl.save
  result$formula <- formula
  names(result$coefficients) <- colnames(result$vcov) <-
    rownames(result$vcov) <- nm.mhurdle(result)
  result$model <- mf
  result
}

mhurdle.fit <- function(start, X, S, P, y, gradient = FALSE, fit = FALSE,
                        score = FALSE, dist = c("l","n","t"), corr = FALSE, ...){
  start.time <- proc.time()
  if (!is.null(corr)) other.coef <- c("sigma","rho")
  else other.coef <- c("sigma")
  KS <- ifelse(is.null(S), 0, ncol(S))
  KP <- ifelse(is.null(P), 0, ncol(P))  
  KX <- ncol(X)
  sel <- KS > 0
  ifr <- KP > 0

  s <- function(param) attr(mhurdle.lnl(param, X = X, S = S, P = P, y = y,
                                   gradient = FALSE, fit = FALSE, score = TRUE,
                                   dist = dist, corr = corr),
                            "score")
  f <- function(param) mhurdle.lnl(param, X = X, S = S, P = P, y = y,
                              gradient = TRUE, fit = FALSE, score = FALSE,
                              dist = dist, corr = corr)
  check.gradient <- TRUE
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
  maxl <- maxLik(f, start = start, iterlim = 2000, ...)
  coefficients <- maxl$estimate
  logLik <- f(coefficients)
  attr(logLik,"df") <- length(coefficients)
  attr(logLik, "y") <- y
  hessian <- maxl$hessian
  fitted <- compute.fitted.mhurdle(coefficients, X, S, P, dist, corr)
  convergence.OK <- maxl$code <= 2
  grad.conv <- maxl$gradient
  elaps.time <- proc.time() - start.time
  nb.iter <- maxl$iterations
  eps <- grad.conv%*%solve(-maxl$hessian)%*%grad.conv
  if (!is.null(corr) && sel){
    rho <-attr(mhurdle.lnl(
                           coefficients, X = X, S = S, P = P, y = y,
                           gradient = FALSE, fit = FALSE, score = TRUE,
                           dist = dist, corr = corr),
               "score")
  }
  else rho <- NULL
  est.stat <- list(elaps.time = elaps.time,
                   nb.iter = nb.iter,
                   eps = eps,
                   method = maxl$type,
                   message = maxl$message
                   )
  class(est.stat) <- "est.stat"
  coef.names <- list(reg   = colnames(X),
                     sel   = colnames(S),
                     ifr   = colnames(P),
                     other = other.coef)
  result <- list(coefficients  = coefficients,
                 vcov          = - solve(maxl$hessian),
                 fitted.values = fitted,
                 logLik        = logLik,
                 gradient      = grad.conv,
                 formula       = NULL,
                 model         = NULL,
                 coef.names    = coef.names,
                 rho           = rho,
                 call          = NULL,
                 est.stat      = est.stat,
                 naive         = NULL
                 )
  if (ncol(X) > 1) class(result) <- c("mhurdle","maxLik")
  result

}

describe <- function(x, which){
  if (which == "sel") result <- attr(x$formula, "rhs")[[1]] != 0
  if (which == "reg") result <- attr(x$formula, "rhs")[[2]] != 0
  if (which == "ifr") result <- attr(x$formula, "rhs")[[3]] != 0
  if (which == "corr") result <- ifelse(is.null(x$call$corr), NULL, x$call$corr)
  result
}
 
# Compute the starting values of the model

start.mhurdle <- function(X, S, P, y, dist, corr){
  sel <- !is.null(S)
  ifr <- !is.null(P)
  
  # for models (3d, 4d), estimates of models (3i, 4i) which can be
  # estimated simply in two parts using fit.simple.mhurdle() are used as
  # starting values
  if (sel && !ifr &&  dist != "n"){
    result <- fit.simple.mhurdle(X, S, y, dist = dist)
    start <- c(result$coefficients,0.1)
  }

  # model (5) tobit : use linear model as starting values
  if (!sel && !ifr &&  dist == "n"){
    lin <- lm(y ~ X - 1)
    start <- c(coef(lin), summary(lin)$sigma)
  }
  
  # model (6, 7, 11) ifr whithout sel
  if (!sel && ifr){
    probit <- glm( (y!=0) ~ P - 1, family = binomial(link='probit'))
    bP <- as.numeric(crossprod(coef(probit),t(P)))
    PP <- pnorm(bP)
    yPP <- y/PP
    lin <- switch(dist,
                  "l" = lm(log(yPP)~ X - 1, subset=y!=0),
                  "n" = lm(yPP ~ X - 1, subset=y!=0),
                  "t" = truncreg(yPP ~ X - 1, subset = y!=0)
                  )
    if(dist != "t") start <- c(coef(lin), coef(probit), summary(lin)$sigma)
    else start <- c(coef(lin)[-length(coef(lin))], coef(probit), coef(lin)[length(coef(lin))])
    if (!is.null(corr)) start <- c(start,0.1)
  }
  
  # model (8), double hurdle
  # use model (3i) as starting values
  if (sel && !ifr && dist == "n"){
    result <- fit.simple.mhurdle(X, S, y, dist = dist)
    start <- result$coefficients
    if (!is.null(corr)) start <- c(start,0.1)
  }
  
  # model (9, 10, 12)
  if (sel && ifr){
    probit.ifr <- glm( (y!=0) ~ P - 1, family = binomial(link='probit'))
    probit.sel <- glm( (y!=0) ~ S - 1, family = binomial(link='probit'))
    coef.ifr <- coef(probit.ifr)
    coef.sel <- coef(probit.sel)
    bP <- as.numeric(crossprod(coef.ifr,t(P)))
    PP <- pnorm(bP)
    bS <- as.numeric(crossprod(coef.sel,t(S)))
    P0 <- mean(y > 0)
    yPP <- y/PP
    lin <- switch(dist,
                  "l" = lm(log(yPP)~ X - 1, subset=y!=0),
                  "n" = lm(yPP ~ X - 1, subset=y!=0),
                  "t" = truncreg(yPP ~ X - 1, subset = y!=0)
                  )
    if(dist != "t")
      start <- c(coef.sel, coef(lin), coef.ifr, summary(lin)$sigma)
    else
      start <- c(coef.sel, coef(lin)[-length(coef(lin))],
                 coef.ifr, coef(lin)[length(coef(lin))])
    if (!is.null(corr)) start <- c(start,0.1)
  }
  start
}

