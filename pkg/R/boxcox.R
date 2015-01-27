boxcoxreg <- function(formula, data, subset, weights, na.action,
                   model = TRUE, y = FALSE, x = FALSE, scaled = FALSE, lambda = 0.5, ...){
    formula.type <- TRUE
    if (class(formula[[3]]) == "name"){
        X <- try(eval(formula[[3]],sys.frame(which=-3)),silent = TRUE)
        if (class(X) == "matrix") formula.type <- FALSE else formula.type <- TRUE
    }
    
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    
    if (formula.type){
        m <- match(c("formula", "data", "subset", "na.action","weights"),
                   names(mf), 0L)
        mf <- mf[c(1L, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1L]] <- as.name("model.frame")
        mf <- eval(mf, parent.frame())
        X <- model.matrix(formula, data = mf)
        Y <- model.response(mf)
        mt <- attr(mf, "terms")
    }
    else{
        Y <- eval(formula[[2]], sys.frame(which=-3))
        mt <- terms(formula)
    }
  
    ## process options
    result <- boxcox.fit(X, Y,  scaled, lambda, ...)
    result$call <- cl
    result$terms <- mt
    if(model) result$model <- mf
    if(y) result$y <- Y
    if(x) result$x <- X
    result
}

boxcox.fit <- function(X, y, scaled, lambda, ...)
{
  dots <- list(...)
  if (is.null(dots$method)) method <- "bfgs" else method <- dots$method
  if (is.null(dots$iterlim)) iterlim <- 50 else iterlim <- dots$iterlim
  if (is.null(dots$print.level)) print.level <- 0 else print.level <- dots$print.level
  if (is.null(dots$fixed)) fixed <- NULL else fixed <- dots$fixed
  
  oldoptions <- options(warn=-1)
  on.exit(options(oldoptions))
  start.time <- proc.time()

  f <- function(param)  ml.boxcox(param,  X = X, y = y,
                                  gradient = FALSE, hessian = FALSE,
                                  fit = FALSE, scaled = scaled)
  g <- function(param){
      attr(ml.boxcox(param, X = X, y = y,
                     gradient = TRUE, hessian = FALSE,
                     fit = FALSE, scaled = scaled), "gradient")
  } 
  h <- function(param){
      attr(ml.boxcox(param, X = X, y = y,
                     gradient = FALSE, hessian = TRUE,
                     fit = FALSE, scaled = scaled), "hessian")
  } 
  Tyinit <- (y ^ lambda - 1) / lambda
  linmod <- lm.fit(X, Tyinit)
#  linmod <- lm.fit(X, y)
  start <- c(coef(linmod), sqrt(sum(linmod$residuals ^ 2) / linmod$df.residual), lambda)
  if (scaled) start[1:ncol(X)] <- start[1:ncol(X)] / start[ncol(X) + 1]
  if (FALSE){
      f0 <-  ml.boxcox(start,  X = X, y = y,
                      gradient = TRUE, hessian = FALSE,
                      fit = FALSE, scaled = scaled)
      agrad <- apply(attr(f0, "gradient"), 2, sum)
      ostart <- start
      of <- sum(f(start))
      ngrad <- c()
      eps <- 1E-5
      for (i in 1:length(start)){
          start <- ostart
          start[i] <- start[i] + eps
          ngrad <- c(ngrad, (sum(f(start)) - of) / eps)
      }
      print(cbind(ngrad, agrad))
      start <- ostart
  }
  ## maxl <- maxLik(f, g, h, start = start, method = method,
  ##                iterlim = iterlim, print.level = print.level)
  maxl <- maxLik(f, g,  h, start = start, method = method,
                 iterlim = iterlim, print.level = print.level, fixed = fixed)
  grad.conv <- g(maxl$estimate)
  coefficients <- maxl$estimate
  vcov <- -solve(maxl$hessian)
  ## fit <- attr(ml.boxcox(coefficients, X = X, y = y,
  ##                         gradient = FALSE, hessian = FALSE,
  ##                         fit = TRUE), "fit")
  ## names(fit) <- rownames(X)
  logLik <- maxl$maximum
  attr(logLik,"df") <- length(coefficients)
  hessian <- maxl$hessian
  convergence.OK <- maxl$code<=2
  elaps.time <- proc.time() - start.time
  nb.iter <- maxl$iterations
  eps <- maxl$gradient%*%solve(-maxl$hessian)%*%maxl$gradient
  est.stat <- list(elaps.time = elaps.time,
                   nb.iter = nb.iter,
                   eps = eps,
                   method = maxl$type,
                   message = maxl$message
                   )
  class(est.stat) <- "est.stat"
  coef.names <- c(colnames(X),"sigma", "lambda")
  names(coefficients) <- rownames(vcov) <- colnames(vcov) <- coef.names

  result <- list(coefficients = coefficients,
                 vcov = vcov,
#                 fitted.values = fit,
                 logLik = logLik,
                 gradient = grad.conv,
		 nobs = length(y),
                 call = NULL,
		 terms = NULL,
                 model = NULL,
		 y = NULL,
		 x = NULL,
                 est.stat = est.stat
                 )
  class(result) <- c("truncreg", "maxLik")
  result
}

ml.boxcox <- function(param, X, y, gradient = FALSE, hessian = FALSE, fit = FALSE, scaled = FALSE){
    beta <- param[1:ncol(X)]
    sigma <- param[ncol(X) + 1]
    lambda <- param[ncol(X) + 2]
    sgn <- sign(lambda)
    bX <- as.numeric(crossprod(beta, t(X)))
    Ty <- (y ^ lambda - 1) / lambda
    DTy <- (log(y) * y ^ lambda - Ty)/ lambda
    HTy <- (log(y) ^ 2 * y ^ lambda - 2 * DTy) / lambda

    resid <- Ty - bX
    trunc <- bX + 1 / lambda
    z <- sgn * trunc / sigma
    mills <- dnorm(z) / pnorm(z)
    lnL <-  - log(sigma) + (lambda - 1) * log(y) + dnorm(resid / sigma, log = TRUE) - pnorm(z, log.p = TRUE)
    if (gradient){
        gbX <- resid / sigma ^ 2 -  sgn * mills / sigma
        gsigma <- - 1 / sigma + resid ^ 2 / sigma ^ 3  +  sgn * mills * trunc / sigma ^ 2
        glambda <- log(y) - resid * DTy / sigma ^ 2 + sgn * mills / sigma / lambda ^ 2
        gradi <- cbind(gbX * X, as.numeric(gsigma), as.numeric(glambda))
        attr(lnL, "gradient") <- gradi
    }
    if (fit){
        attr(lnL, "fit") <- bX
    }
    if (hessian){
        bb <- - 1 / sigma ^ 2 + mills * (z + mills) / sigma ^ 2
        ss <- 1 / sigma ^ 2 - 3 * resid ^ 2 / sigma ^ 4 - 2 * mills * sgn * trunc / sigma ^ 3 +
            mills * (z + mills) * trunc / sigma ^ 4
        bs <- - 2 * resid / sigma ^ 3 + sgn * mills / sigma ^ 2 -
            mills * (mills + z) * trunc / sigma ^ 3
        bl <- DTy /sigma ^ 2 - mills * (z + mills) / sigma ^ 2 / lambda ^ 2
        sl <- 2 * resid * DTy / sigma ^ 3 + trunc * mills * (z + mills) / sigma ^ 3 / lambda ^ 2 - sgn * mills / sigma ^ 2 / lambda ^ 2
        ll <- - (resid * HTy + DTy ^ 2) / sigma ^ 2 + mills * (z + mills) / sigma ^ 2 / lambda ^ 4 - 2 * mills / sigma / lambda ^ 3
        bb <- crossprod(bb * X, X)
        bs <- apply(bs * X, 2, sum)
        bl <- apply(bl * X, 2, sum)
        ss <- sum(ss)
        ll <- sum(ll)
        sl <- sum(sl)
        h <- rbind(cbind(bb, bs), c(bs, ss))
        h <- cbind(h, c(bl, sl))
        h <- rbind(h, c(bl, sl, ll))
        attr(lnL,"hessian") <- h
    }
    lnL
}
