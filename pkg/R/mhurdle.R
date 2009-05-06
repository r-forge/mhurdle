mhurdle <- function(formula, data, subset, weights, na.action,
                 start = NULL,
                 dist = c("l","t","n"),
                 res = FALSE,
                 sel = TRUE,
                 ifr = FALSE,
                 corr = FALSE, ...){
  dots <- list(...)
  oldoptions <- options(warn=-1)
  on.exit(options(oldoptions))
  cl <- match.call()
  posT <- as.list(cl) == "T"
  posF <- as.list(cl) == "F"
  cl[posT] <- TRUE
  cl[posF] <- FALSE
  cl.save <- cl
  method <- dots$method

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
  model <- z[as.character(res),
             as.character(sel),
             as.character(ifr),
             as.character(corr),
             dist]

  if (is.na(model)) stop("not implemented  model")
  if (model == "irrelevant") stop("irrelevant model")

  # compute the model.frame and the model.matrix
  formula <- Formula(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action","weights"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf$formula <- formula
  mf <- eval(mf, parent.frame())
  X <- model.matrix(formula, data = mf, part = "first")
  if (length(formula) == 2) S <- model.matrix(formula, data = mf, part = "second")
  else S <- NULL
  y <- model.response(mf)
  fit <- NULL

  # compute the "naive" model
  n <- length(y)
  Xo <- So <- matrix(1, nrow = n)
  colnames(Xo) <- colnames(So) <- "(intercept)"
  alphaX <- ifelse(dist == "l", mean(log(y)[y>0]), mean(y[y>0]))
  sigma <-  ifelse(dist == "l", sd(log(y)[y>0]), sd(y[y>0]))
  # pas correct pour tronqué
  alphaS <- qnorm(mean(y>0))
  starto <- c(alphaS, alphaX, sigma)
  if (corr) starto <- c(starto, 0)
  za <- mhurdle.fit(starto, Xo, So, y, dist, res, sel, ifr, corr, ...)
  logLiko <- za$logLik
  naive <- list(coefficients = coef(za), logLik = za$logLik, est.stat = za$est.stat)
  # sel without corr, can be estimated simply in two parts

  if (sel && !res && !corr){
    result <- sel.simple(X, S, y, dist = dist, ...)
    result$call <- cl.save
    result$model <- mf
    result$formula <- formula
  }
  else{
    if (is.null(start)){
      if (sel){
        if (res) cl$res <- FALSE
        cl$corr <- FALSE
        start <- coef(eval(cl,parent.frame()))
      }
      if (ifr){
        probit <- glm( (y!=0) ~ S - 1, family = binomial(link='probit'))
        bS <- as.numeric(crossprod(coef(probit),t(S)))
        PS <- pnorm(bS)
        yPS <- y/PS
        lin <- switch(dist,
                      "l" = lm(log(yPS)~ X - 1, subset=y!=0),
                      "n" = lm(yPS ~ X - 1, subset=y!=0),
                      "t" = truncreg(yPS ~ X - 1, subset = y!=0)
                      )
        start <- c(coef(probit),coef(lin))
        if(dist != "t") start <- c(start,summary(lin)$sigma)
      }
      if (corr) start <- c(start,0)
      result <- mhurdle.fit(start, X, S, y, dist, res, sel, ifr, corr, ...)
      result$call <- cl.save
      result$formula <- formula
      names(result$coefficients) <- colnames(result$vcov) <-
        rownames(result$vcov) <- nm.mhurdle(result)
      result$model <- mf
    }
  }
  class(naive) <- "mhurdle"
  result$naive <- naive
  result
}

mhurdle.fit <- function(start, X, S, y, dist, res, sel, ifr, corr, ...){

  start.time <- proc.time()
  if (corr) other.coef <- c("sigma","rho")
  else other.coef <- c("sigma")

  f <- function(param) ml.tot(param, X = X, S = S, y = y,
                              gradient = FALSE, fit = FALSE, score = FALSE,
                              dist = dist, res = res,
                              sel = sel, ifr = ifr,
                              corr = corr)
  g <- function(param) attr(ml.tot(param, X = X, S = S, y = y,
                                   gradient = TRUE, fit = FALSE, score = FALSE,
                                   dist = dist, res = res,
                                   sel = sel, ifr = ifr,
                                   corr = corr),
                            "gradient")
  s <- function(param) attr(ml.tot(param, X = X, S = S, y = y,
                                   gradient = FALSE, fit = FALSE, score = TRUE,
                                   dist = dist, res = res,
                                   sel = sel, ifr = ifr,
                                   corr = corr),
                            "score")
  maxl <- maxLik(f, g, start = start, ...)
  coefficients <- maxl$estimate
  grad.conv <- g(coefficients)
  logLik <- f(coefficients)
  attr(logLik,"df") <- length(coefficients)
  attr(logLik, "y") <- y
  hessian <- maxl$hessian
  fitted <- compute.fitted.mhurdle(coefficients, X, S, dist, res, sel, ifr, corr)
  convergence.OK <- maxl$code <= 2
  grad.conv <- maxl$gradient
  elaps.time <- proc.time() - start.time
  nb.iter <- maxl$iterations
  eps <- grad.conv%*%solve(-maxl$hessian)%*%grad.conv
  if (!corr && !ifr) rho <- s(coefficients) else rho <- NULL
  est.stat <- list(elaps.time = elaps.time,
                   nb.iter = nb.iter,
                   eps = eps,
                   method = maxl$type,
                   message = maxl$message
                   )
  
  class(est.stat) <- "est.stat"
  if (sel){
    coef.names <- list(sel = colnames(S),
                       reg = colnames(X),
                       other = other.coef)
  }
  if (ifr){
    coef.names <- list(ifr = colnames(S),
                       reg = colnames(X),
                       other = other.coef)
  }
  result <- list(coefficients = coefficients,
                 vcov = - solve(maxl$hessian),
                 fitted.values = fitted,
                 logLik = logLik,
                 gradient = grad.conv,
                 formula = NULL,
                 model = NULL,
                 coef.names = coef.names,
                 rho = rho,
                 call = NULL,
                 est.stat = est.stat
                 )
  class(result) <- c("mhurdle","maxLik")
  result
}

describe <- function(x, which){
  cl <- x$call
  if (which == "sel")  z <- ifelse(is.null(cl$sel) , TRUE , cl$sel)
  if (which == "corr") z <- ifelse(is.null(cl$corr), FALSE, cl$corr)
  if (which == "ifr")  z <- ifelse(is.null(cl$ifr) , FALSE, cl$ifr)
  if (which == "dist") z <- ifelse(is.null(cl$dist), "l"  , cl$dist)
  z
}
  
