# Box Cox with a location parameter

boxcoxreg <- function(formula, data, subset, weights, na.action,
                      model = TRUE, y = FALSE, x = FALSE,
                      lambda = 0.5, alpha = 0.5, start = NULL,
                      truncated = TRUE, check.gradient = FALSE, robust = TRUE, ...){
    
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    
    m <- match(c("formula", "data", "subset", "na.action","weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    X <- model.matrix(formula, data = mf)
    Y <- model.response(mf)
    mt <- attr(mf, "terms")
    
    if (truncated){
        X <- X[Y > 0, ]
        Y <- Y[Y > 0]
    }
    
    ## process options
    result <- boxcox.fit(X = X, y = Y, lambda = lambda, alpha = alpha, start = start, check.gradient = check.gradient, truncated = truncated, robust = robust, ...)
    result$call <- cl
    result$terms <- mt
    if(model) result$model <- mf
    if(y) result$y <- Y
    if(x) result$x <- X
    result
}

boxcox.fit <- function(X, y, lambda, alpha, start, check.gradient, truncated, robust,  ...){

    if (robust){
        fisigma <- function(x) log(x)
        fsigma <- function(x) exp(x)
        gsigma <- function(x) exp(x)
        fialpha <- function(x) log(x)
        falpha <- function(x) exp(x)
        galpha <- function(x) exp(x)
        ## fialpha <- function(x) qnorm(x)
        ## falpha <- function(x) pnorm(x)
        ## galpha <- function(x) dnorm(x)
    }
    else{
        fisigma <- function(x) x
        fsigma <- function(x) x
        gsigma <- function(x) 1
        fialpha <- function(x) x
        falpha <- function(x) x
        galpha <- function(x) 1
    }


    dots <- list(...)
    if (is.null(dots$method)) method <- "bfgs" else method <- dots$method
    if (is.null(dots$iterlim)) iterlim <- 100 else iterlim <- dots$iterlim
    if (is.null(dots$print.level)) print.level <- 0 else print.level <- dots$print.level
    if (is.null(dots$fixed)) fixed <- NULL else fixed <- dots$fixed
  
    oldoptions <- options(warn = -1)
    on.exit(options(oldoptions))
    start.time <- proc.time()

    f <- function(param)      ml.boxcox(param, X = X, y = y, gradient = FALSE, truncated = truncated, robust = robust)
    g <- function(param) attr(ml.boxcox(param, X = X, y = y, gradient = TRUE,  truncated = truncated, robust = robust), "gradient")

    if (is.null(start)){
        if (lambda != 0) Tyinit <- (y ^ lambda - 1) / lambda else Tyinit <- log(y)
        linmod <- lm.fit(X[y > 0, , drop = FALSE], Tyinit[y > 0])
        sigma <- sqrt(sum(linmod$residuals ^ 2) / length(Tyinit[y > 0]))
        sigma <- fisigma(sigma)
        alpha <- fialpha(alpha)
        start <- c(coef(linmod), sigma = sigma, lambda = lambda, alpha = alpha)
    }
    else{
        start["alpha"] <- fialpha(start["alpha"])
        start["sigma"] <- fisigma(start["sigma"])
    }
    
    if (check.gradient){
        f0 <-  ml.boxcox(start,  X = X, y = y, gradient = TRUE, truncated = truncated, robust = robust)
        agrad <- apply(attr(f0, "gradient"), 2, sum)
        ostart <- start
        of <- sum(f(start))
        ngrad <- c()
        eps <- 1E-6
        for (i in 1:length(start)){
            start <- ostart
            start[i] <- start[i] + eps
            ngrad <- c(ngrad, (sum(f(start)) - of) / eps)
        }
        print(ostart)
        print(cbind(ngrad, agrad))
        start <- ostart
        stop()
    }
    
    maxl <- maxLik(f,  g, start = start, method = method, finalHessian = "BHHH",
                   iterlim = iterlim, print.level = print.level, fixed = fixed)
    
    actPars <- activePar(maxl)
    grad.conv <- g(maxl$estimate)

    
    coefficients <- maxl$estimate

    if (robust){
        fsigma <- function(x) exp(x)
        gsigma <- function(x) exp(x)
        falpha <- function(x) exp(x)
        galpha <- function(x) exp(x)
        ## falpha <- function(x) pnorm(x)
        ## galpha <- function(x) dnorm(x)
    }
    else{
        fsigma <- function(x) x
        gsigma <- function(x) 1
        falpha <- function(x) x
        galpha <- function(x) 1
    }
    gtheta <- rep(1, length(coefficients))
    names(gtheta) <- names(coefficients)
    gtheta["sigma"] <- gsigma(coefficients["sigma"])
    gtheta["alpha"] <- galpha(coefficients["alpha"])
    coefficients["sigma"] <- fsigma(coefficients["sigma"])
    coefficients["alpha"] <- falpha(coefficients["alpha"])
    vcov <- diag(gtheta) %*% vcov(maxl) %*% diag(gtheta)
    logLik <- maxl$maximum
    attr(logLik,"df") <- length(coefficients)
    hessian <- maxl$hessian
    convergence.OK <- maxl$code <= 2
    elaps.time <- proc.time() - start.time
    nb.iter <- maxl$iterations
    eps <- maxl$gradient[actPars] %*% vcov[actPars, actPars] %*% maxl$gradient[actPars]
    est.stat <- list(elaps.time = elaps.time,
                     nb.iter = nb.iter,
                     eps = eps,
                     method = maxl$type,
                     message = maxl$message
                     )
    class(est.stat) <- "est.stat"
    
    score <- NULL
    if (truncated & ( (alpha == 0 & ! robust) | (alpha == -Inf & robust))){
        galpha <- apply(grad.conv, 2, sum)["lnL.alpha"]
        valpha <- solve(crossprod(grad.conv))["lnL.alpha", "lnL.alpha"]
        score <- c(gradient = galpha, var = valpha, stat = galpha * sqrt(valpha))
    }

    if (! truncated & lambda == 0){
        glambda <- apply(grad.conv, 2, sum)["lnL.lambda"]
        vlambda <- solve(crossprod(grad.conv))["lnL.lambda", "lnL.lambda"]
        score <- c(gradient = glambda, var = vlambda, stat = glambda * sqrt(vlambda))
    }
    
    result <- list(coefficients = coefficients,
                   vcov = vcov,
                   logLik = logLik,
                   gradient = maxl$gradient,
                   gradientObs = maxl$gradientObs,
                   nobs = length(y),
                   call = NULL,
                   terms = NULL,
                   model = NULL,
                   y = NULL,
                   x = NULL,
                   est.stat = est.stat,
                   score = score
                   )
    class(result) <- c("truncreg", "maxLik")
    result
}

ml.boxcox <- function(param, X, y, gradient, truncated, robust){

    if (robust){
        fsigma <- function(x) exp(x)
        gsigma <- function(x) exp(x)
        falpha <- function(x) exp(x)
        galpha <- function(x) exp(x)
        ## falpha <- function(x) pnorm(x)
        ## galpha <- function(x) dnorm(x)
    }
    else{
        fsigma <- function(x) x
        gsigma <- function(x) 1
        falpha <- function(x) x
        galpha <- function(x) 1
    }


    beta <- param[1:ncol(X)]
    
    sigma <- param[ncol(X) + 1]
    g.sigma <- gsigma(sigma)
    sigma <- fsigma(sigma)
    
    lambda <- param[ncol(X) + 2]
    
    alpha <- param[ncol(X) + 3]
    g.alpha <- galpha(alpha)
    alpha <- falpha(alpha)
    
    sgn <- sign(lambda)
    myInf <- 1000
    bX <- as.numeric(crossprod(beta, t(X)))
    
    if (lambda == 0){
        T0 <- log(alpha)
        Ty <- log(y + alpha)
        T.lambda <- 1 / 2 * log(y + alpha) ^ 2
    }
    else{
        T0 <- (alpha ^ lambda - 1) / lambda
        Ty <- ((y + alpha) ^ lambda - 1) / lambda
        T.lambda <- Ty * log(y + alpha) -  (Ty - log(y + alpha)) / lambda
    }

    T.alpha <- (y + alpha) ^ (lambda - 1)
    resid <- Ty - bX
    mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))

    if (truncated){
        # truncated model
        if (alpha == 0){
            if (lambda > 0){
                B1 <- (- 1 / lambda - bX) / sigma
                B2 <- + myInf
                B1.sigma <- (1 / lambda + bX) / sigma ^ 2
                B1.lambda <- 1 / sigma / lambda ^ 2
                B1.beta <- - 1 / sigma
                B1.alpha <- 0
                B2.sigma <- B2.beta <- B2.alpha <- B2.lambda <- 0
            }
            if (lambda < 0){
                B1 <- - myInf
                B2 <- (- 1 / lambda - bX) / sigma
                B1.sigma <- B1.beta <- B1.alpha <- B1.lambda <- 0
                B2.sigma <- 1 / lambda / sigma ^ 2
                B2.lambda <- 1 / sigma / lambda ^ 2
                B2.beta <- - 1 / sigma
                B2.alpha <- 0
            }
            if (lambda == 0){
                B1 <- - myInf
                B2 <- myInf
                B1.sigma <- B1.beta <- B1.alpha <- B1.lambda <- 0
                B2.sigma <- B2.beta <- B2.alpha <- B2.lambda <- 0
                
            }
        }
        else{
            B1 <-  (T0  - bX) / sigma
            B1.sigma <- - (T0 - bX) / sigma ^ 2
            B1.beta <- - 1 / sigma
            B1.alpha <- alpha ^ (lambda - 1) / sigma
            if (lambda != 0) B1.lambda <- T0 * log(alpha) - (T0 - log(alpha)) / lambda
            else B1.lambda <- 1 / 2 * log(alpha) ^ 2
            if (lambda >= 0){
                B2 <- myInf
                B2.sigma <- B2.beta <- B2.alpha <- B2.lambda <- 0
            }
            if (lambda < 0){
                B2 <- (- 1 / lambda - bX) / sigma
                B2.sigma <- (1 / lambda + bX) / sigma ^ 2
                B2.beta <- - 1 /sigma
                B2.lambda <- 1 / sigma / lambda ^ 2
                B2.alpha <- 0
            }
        }
        PB1 <- pnorm(B1)
        PB2 <- pnorm(B2)
        lnL <-  - log(sigma) + (lambda - 1) * log(y + alpha) + dnorm(resid / sigma, log = TRUE) - log(PB2 - PB1)
        if (gradient){
            MU1 <- mills(B1)
            MU2 <- mills(B2)
            OM1 <- PB1 / (PB2 - PB1)
            OM2 <- PB2 / (PB2 - PB1)
            lnL.sigma <- - 1 / sigma +  resid ^ 2 / sigma ^ 3 - (OM2 * MU2 * B2.sigma - OM1 * MU1 * B1.sigma)
            lnL.beta <-  resid / sigma ^ 2 - (OM2 * MU2 * B2.beta  - OM1 * MU1 * B1.beta)
            lnL.lambda <- log(y + alpha) - resid / sigma ^ 2 * T.lambda - (OM2 * MU2 * B2.lambda - OM1 * MU1 * B1.lambda)
            lnL.alpha <- (lambda - 1) / (y + alpha) - resid / sigma ^ 2 * T.alpha - (OM2 * MU2 * B2.alpha - OM1 * MU1 * B1.alpha)
            gradi <- cbind(lnL.beta * X, lnL.sigma * g.sigma, lnL.lambda, lnL.alpha)
            attr(lnL, "gradient") <- gradi
        }
    }
    else{
        # censored model
        # truncation point
        if (lambda != 0){
            ZT <- - (1 / lambda + bX) / sigma
            sgn <- sign(lambda)
            ZT.lambda <- 1  / lambda ^ 2 / sigma
            ZT.beta <- - 1 / sigma
            ZT.sigma <- (1 / lambda + bX)/ sigma ^ 2
            ZT.alpha <- 0
            LD <- pnorm(- sgn * ZT, log.p = TRUE)
            MUD <- mills(- sgn * ZT)
            LD.beta <-   - MUD * sgn * ZT.beta
            LD.sigma <-  - MUD * sgn * ZT.sigma
            LD.lambda <- - MUD * sgn * ZT.lambda
            LD.alpha <-  - MUD * sgn * ZT.alpha
        }
        else LD <- LD.beta <- LD.sigma <- LD.lambda <- LD.alpha <- 0
        
        # z(y = 0)
        Z0 <- (T0 - bX) / sigma
        Z0.beta <- - 1 / sigma
        Z0.sigma <- - (T0 - bX) / sigma ^ 2 
        Z0.alpha <- alpha ^ (lambda - 1) / sigma
        if (lambda == 0) Z0.lambda <- 1 / 2 * log(alpha) ^ 2 / sigma
        else Z0.lambda <- (T0 * log(alpha) - (T0 - log(alpha)) / lambda) / sigma

        # PI 0
        if (lambda > 0){
            MU0 <- mills(Z0)
            MUT <- mills(ZT)
            P0 <- pnorm(Z0)
            PT <- pnorm(ZT)
            OM0 <- P0 / (P0 - PT)
            OMT <- PT / (P0 - PT)
            LPI0 <- log(P0 - PT)
            PI0.beta <-   OM0 * MU0 * Z0.beta   - OMT * MUT * ZT.beta
            PI0.sigma <-  OM0 * MU0 * Z0.sigma  - OMT * MUT * ZT.sigma
            PI0.lambda <- OM0 * MU0 * Z0.lambda - OMT * MUT * ZT.lambda
            PI0.alpha <-  OM0 * MU0 * Z0.alpha  - OMT * MUT * ZT.alpha
        }
        else{
            LPI0 <- pnorm(Z0, log.p = TRUE)
            MU0 <- mills(Z0)
            PI0.beta <-   MU0 * Z0.beta
            PI0.sigma <-  MU0 * Z0.sigma
            PI0.lambda <- MU0 * Z0.lambda
            PI0.alpha <-  MU0 * Z0.alpha
        }
        lnL <- vector(mode = "numeric", length = length(y))
        lnL[y == 0] <- LPI0[y == 0]
        lnL[y >  0] <-  - log(sigma) + (lambda - 1) * log(y[y > 0] + alpha) + dnorm(resid[y > 0] / sigma, log = TRUE)
        lnLP <- sum(lnL[y > 0])
        lnL0 <- sum(lnL[y == 0])
        lnLD <- sum(- LD)
        lnL <- lnL - LD
        if (any(is.na(lnL))){
            stop("NA values in lnL")
        }
        if (any(abs(lnL) > myInf)){
            lnL[abs(lnL) > myInf] <- sign(lnL[abs(lnL) > myInf]) * myInf
            warning("Huge values of lnL")
        }
        if (gradient){
            lnL.beta <- lnL.sigma <- lnL.lambda <- lnL.alpha <- vector(mode = "numeric", length = length(y))
            
            lnL.beta  [y == 0] <- PI0.beta  [y == 0]
            lnL.sigma [y == 0] <- PI0.sigma [y == 0]
            lnL.lambda[y == 0] <- PI0.lambda[y == 0]
            lnL.alpha [y == 0] <- PI0.alpha [y == 0]

            lnL.beta  [y > 0] <-  (resid / sigma ^ 2)                                      [y > 0]
            lnL.sigma [y > 0] <- (- 1 / sigma +  resid ^ 2 / sigma ^ 3)                    [y > 0]
            lnL.lambda[y > 0] <- (log(y + alpha) - resid / sigma ^ 2 * T.lambda)           [y > 0]
            lnL.alpha [y > 0] <- ((lambda - 1) / (y + alpha) - resid / sigma ^ 2 * T.alpha)[y > 0]

            
            lnL.beta <-   lnL.beta   - LD.beta
            lnL.sigma <-  lnL.sigma  - LD.sigma
            lnL.lambda <- lnL.lambda - LD.lambda
            lnL.alpha <-  lnL.alpha  - LD.alpha

            gradi <- cbind(lnL.beta * X, lnL.sigma * g.sigma, lnL.lambda, lnL.alpha * g.alpha)
            attr(lnL, "gradient") <- gradi
        }
    }
    lnL
}
