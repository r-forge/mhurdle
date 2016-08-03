
log2 <- function(x) ifelse(x > 0,log(x), 0)

mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))

# a function to construct block-diagonal matrix (can't remember from
# which package it is borrowed)
bdiag <- function(...){
  if (nargs() == 1)
    x <- as.list(...)
  else
    x <- list(...)
  n <- length(x)
  if(n == 0) return(NULL)
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

# print a summary of the non-linear opimization
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

# Compute the estimation of one-equation models. Those are only
# relevant for distributions that admit negative values of y^*, namely
# ln2, bc2, ihs and n

onequation.mhurdle <- function(X2, y, dist = NULL){
    if (dist %in% c("ln", "tn", "bc")) stop("the specified model doesn't allow zero observations")
    if (dist == "n"){
        # the tobit model, efficiently estimated using Survreg from
        # the survival package
        tobit <- survreg(Surv(y, y > 0, type = "left") ~ X2 - 1, dist = "gaussian")
        beta <- coef(tobit)
        bX <- as.numeric(crossprod(t(X2), beta))
        resid <- y - bX
        sigma <- tobit$scale
        pmills <- mills(- bX / sigma)
        gbeta <- - (y == 0) * mills(- bX / sigma) / sigma + (y > 0) * (y - bX) / sigma ^ 2
        gsigma <- (y == 0) * mills(- bX / sigma) * bX / sigma ^ 2 +
            (y > 0) * (- 1 / sigma + (y - bX) ^ 2 / sigma ^ 3)
        gradi <- cbind(gbeta * X2, sigma = gsigma)
        result <- list(coefficients = c(beta, sigma = sigma),
                       vcov = solve(crossprod(gradi)),
                       logLik = logLik(tobit),
                       gradient = gradi,
                       fitted = cbind(1 - pnorm(bX / sigma), bX + sigma * mills(bX / sigma)))
    }
    if (dist == "ln2"){
        thelm <- lm(log(y + 1) ~ X2 - 1)
        sigma <- summary(thelm)$sigma
        beta <- coef(thelm)
        start <- c(beta, sigma = sigma, alpha = 1)
        ltobit.lnl <- function(param, fitted = FALSE){
            beta <- param[1 : (length(param) - 2)]
            sigma <- param[length(param) - 1]
            alpha <- param[length(param)]
            bX <- as.numeric(crossprod(t(X2), beta))
            lnL <- (y == 0) * pnorm( (log(    alpha) - bX) / sigma, log.p = TRUE) +
                (y > 0) * ( dnorm(   (log(y + alpha) - bX) / sigma, log   = TRUE) -
                               log(y + alpha) - log(sigma)
                           )
            gbeta <- - (y == 0) * mills( (log(alpha) - bX) / sigma) / sigma +
                (y > 0) * (log(y + alpha) - bX) / sigma ^ 2
            gsigma <- - (y == 0) * mills( (log(alpha) - bX) / sigma) * (log(alpha) - bX) / sigma ^ 2 +
                (y > 0) * ( - 1 / sigma + (log(alpha + y) - bX) ^ 2 / sigma ^ 3)
            galpha <- (y == 0) * mills( (log(alpha) - bX) / sigma) / (alpha * sigma) -
                (y > 0) * ( 1 / (y + alpha) + ( log(y + alpha) - bX) / (sigma ^ 2 * (y + alpha)))
            attr(lnL, "gradient") <- cbind(gbeta * X2, gsigma, galpha)
            if (fitted){
                P0 <- pnorm( (log(alpha) - bX) / sigma)
                Ep <- exp(bX + 0.5 * sigma ^ 2) * pnorm( (bX - log(alpha)) / sigma + sigma) / (1 - P0)
                attr(lnL, "fitted") <- cbind(P0, Ep)
            }
            lnL
        }
        result <- maxLik(ltobit.lnl, start = start, method = "bhhh", print.level = 0)
        fitted <- attr(ltobit.lnl(coef(result), fitted = TRUE), "fitted")
        result <- list(coefficients = coef(result),
                       vcov = vcov(result),
                       logLik = logLik(result),
                       gradient = result$gradientObs,
                       fitted = fitted)

    }
    if (dist == "bc2"){
        result <- boxcox.fit(X2, y, lambda = 0, alpha = 1, robust = TRUE, start = NULL, check.gradient = FALSE, truncated = FALSE)
    }

    coef.names <- list(h1    = NULL,
                       h2    = colnames(X2),
                       h3    = NULL,
                       sd    = "sd",
                       h4    = NULL,
                       corr  = NULL,
                       tr    = NULL,
                       pos   = NULL)
    if (dist == "bc2") coef.names$tr <- "tr"
    if (dist %in% c("bc2", "ln2")) coef.names$pos <- "pos"

    
    result <- list(coefficients = result$coefficients, 
                   vcov = result$vcov,
                   fitted.values = result$fitted,
                   logLik = result$logLik,
                   gradient = result$gradient,
                   model = NULL,
                   formula = NULL,
                   coef.names = coef.names,
                   call = NULL
                   )
    class(result) <- c("mhurdle","maxLik")
    result

    
    return(result)
}
        

# Compute the estimation of hurdle models in the cases where it can be
# done using two independent estimations (a binomial logit model and a
# normal/log-normal/truncated linear model). This is relevant for
# uncorrelated models with selection

seperate.mhurdle <- function(X1, X2, y, dist = NULL){
    probit <- glm(y != 0 ~ X1 - 1, family = binomial(link = "probit"))
    # Computation of the likelihood for zero observations
    beta1 <- coef(probit)
    bX1 <- as.numeric(crossprod(beta1, t(X1)))
    mills1 <- mills(bX1)
    mills1m <- mills(- bX1)
    L.null <- (y == 0) * log(1 - pnorm(bX1))
    gbX1 <- (y == 0) * (- mills1m) + (y != 0) * mills1
    
    # Computation of the likelihood for positive observations for log-normal distribution
    if (dist %in% "ln"){
        lin <- lm(log(y) ~ X2 - 1, subset = y != 0)
        beta2 <- coef(lin)[1:ncol(X2)]
        bX2 <- as.numeric(crossprod(beta2, t(X2)))
        logy <- rep(0, length(y))
        logy[y != 0] <- log(y[y != 0])
        resid <- (logy - bX2)
        df <- df.residual(lin)
        np <- sum(y != 0)
        scr <- sum(resid[y != 0] ^ 2)
        sigma <- sqrt(scr / np)
        L.pos <- (y != 0) * (pnorm(bX1, log.p = TRUE) + dnorm(resid / sigma, log = TRUE) -
                                 log(sigma) - logy)
        gbX2 <- (y != 0) * (resid / sigma ^ 2)
        gsigma <- (y != 0) * (resid ^ 2 / sigma ^ 3 - 1 / sigma)
        gradi <- cbind(gbX1 * X1, gbX2 * X2, as.numeric(gsigma))
        dss <- - 3 * scr / sigma ^ 4 + sum(y != 0) / sigma ^ 2
        vcov <- bdiag(vcov(probit), vcov(lin) / np * df, - 1 / dss)
        coef <- c(coef(probit), coef(lin), sigma)
        P0 <- pnorm(bX1)
        E <- exp(bX2 + 0.5 * sigma ^ 2) / (1 - P0)
        fitted <- cbind(zero = pnorm(bX1), pos = E)
    }
    else{
        if (dist == "tn"){
                          lin <- truncreg(y ~ X2 - 1, subset = y != 0)
                          lin$gradient <- lin$gradientObs
                      }
        if (dist == "bc"){
            lin <- boxcoxreg(y ~ X2 - 1, subset = y != 0, alpha = 0, lambda = 0, fixed = "alpha", method = "bhhh", print.level = 0)
            K <- length(coef(lin)) - 1
            lin$vcov <- lin$vcov[1:K, 1:K]
            lin$coefficients <- lin$coefficients[1:K]
            lin$gradient <- lin$gradientObs[, 1:K]
        }
        bX2 <- as.numeric(crossprod(coef(lin)[1:ncol(X2)], t(X2)))
        L.pos <- as.numeric(logLik(lin)) + sum( (y != 0) * (pnorm(bX1, log.p = TRUE)))
        vcov <- bdiag(vcov(probit), vcov(lin))
        coef <- c(coef(probit), coef(lin))
        fit <- cbind(zero = bX1, pos = fitted(lin))
        g2 <- matrix(0, nrow = length(y), ncol = ncol(lin$gradient))
        g2[y != 0, ] <- lin$gradient
        gradi <- cbind(gbX1 * X1, g2)
        vcov <- bdiag(vcov(probit) , vcov(lin))
        if (dist == "tn"){
            P0 <- pnorm(bX1)
            E <- (bX2 + mills(bX2)) / (1 - P0)
            fitted <- cbind(zero = P0, pos = E)}
        else fitted <- NULL
    }
    coef.names <- list(h1    = colnames(X1),
                       h2    = colnames(X2),
                       h3    = NULL,
                       sd    = "sd",
                       h4    = NULL,
                       corr  = NULL,
                       tr    = NULL,
                       pos   = NULL)
    if (dist == "bc") coef.names$tr <- "tr"
    logLik <- structure(sum(L.null) + sum(L.pos), df = length(coef), nobs = length(y), class = "logLik")
    result <- list(coefficients = coef, 
                   vcov = vcov,
                   fitted.values = fitted,
                   logLik = logLik,
                   gradient = gradi,
                   model = NULL,
                   formula = NULL,
                   coef.names = coef.names,
                   call = NULL
                   )
    class(result) <- c("mhurdle","maxLik")
    result
}

# Compute the "naive" model, i.e. a model with no explanatory
# variables.  Full version with correlation ; not used because
# identification problems

lnl.naive <- function(param, dist = c("ln", "tn", "n", "ln2"), moments,
                      h1 = TRUE, h3 = FALSE){
    dist <- match.arg(dist)
    n <- moments[1]
    ym <- moments[2]
    s2 <- moments[3]
    if (h1){
        alpha1 <- param[1]
        alpha2 <- param[2]
        param <- param[-c(1,2)]
    }
    else{
        alpha2 <- param[1]
        param <- param[-1]
    }
    if (h3){
        alpha3 <- param[1]
        param <- param[-1]
    }
    sigma <- param[1]
    
    if (h1){
        Phi1 <- pnorm(alpha1)
        phi1 <- dnorm(alpha1)
    }
    else Phi1 <- 1
    if (h3){
        Phi3 <- pnorm(alpha3)
        phi3 <- dnorm(alpha3)
    }
    else Phi3 <- 1
    Phi2 <- pnorm(alpha2/sigma)
    phi2 <- dnorm(alpha2/sigma)
    scr <- ifelse(dist == "ln",
                  s2 + (ym + log(Phi3) - alpha2)^2,
                  Phi3^2*(s2 + (ym - alpha2/Phi3)^2)
                  )
    Pbiv <- Phi1 * Phi2
    Phi1bis <- Phi1
    s2term <- 0
    P0 <- switch(dist,
                 "ln" = 1 - Phi1 * Phi3,
                 "ln2" = 1 - Phi1 * Phi3,
                 "tn" = 1 - Pbiv/Phi2 * Phi3,
                 "n" = 1 - Pbiv * Phi3
                 )
    lnPos <-
        -log(sigma) - 0.5 * log(2*pi) -
            scr/(2*sigma^2) +
                log(Phi3) +
                    (log(Phi1bis)+s2term) * h1 -
                        ym * (dist == "ln") +
                            log(Phi3) * (dist != "ln") -
                                log(Phi2) * (dist == "tn")
    
    n * log(P0) + (1 - n) * lnPos
}


## Ostart.mhurdle <- function(X1, X2, X3, y, dist){
##     h1 <- ! is.null(X1)
##     h3 <- ! is.null(X3)
##     # for models (2ld, 2td), estimates of models (2li, 2ti) which can be
##     # estimated simply in two parts using seperate.mhurdle() are used
##     # as starting values
##     if (h1 && ! h3 &&  dist %in% c("ln", "tn")){
##         result <- seperate.mhurdle(X1, X2, y, dist = dist)
##         start <- result$coefficients
##     }
##     # model (3) tobit : use linear model as starting values
##     if (!h1 && !h3 &&  dist == "n"){
##         lin <- lm(y ~ X2 - 1)
##         sdest <- summary(lin)$sigma
##         start <- c(coef(lin), sigma = sdest)
##     }
##     # model (4, 7) h3 whithout h1
##     if (!h1 && h3){
##         probit <- glm( (y != 0) ~ X3 - 1, family = binomial(link = 'probit'))
##         bX3 <- as.numeric(crossprod(coef(probit), t(X3)))
##         Phi3 <- pnorm(bX3)
##         yPP <- y * Phi3
##         lin <- switch(dist,
##                       "ln" = lm(log(yPP)  ~ X2 - 1, subset = y != 0),
##                       "n" = lm(yPP        ~ X2 - 1, subset = y != 0),
##                       "tn" = truncreg(yPP ~ X2 - 1, subset = y != 0, scaled = TRUE)
##                       )
##         if (dist %in% c("n", "ln")) start <- c(coef(lin), coef(probit), sigma = summary(lin)$sigma)
##         if (dist == "tn") start <- c(coef(lin)[- length(coef(lin))], coef(probit), sigma = coef(lin)[length(coef(lin))])
##     }
##     # model (5), double hurdle use model (3i) as starting values
##     if (h1 && !h3 && dist == "n"){
##         result <- seperate.mhurdle(X1, X2, y, dist = dist)
##         start <- result$coefficients
##     }
##     # model (6 and 8)
##     if (h1 && h3){
##         probit.h3 <- glm( (y != 0) ~ X3 - 1, family = binomial(link = 'probit'))
##         probit.h1 <- glm( (y != 0) ~ X1 - 1, family = binomial(link = 'probit'))
##         beta3 <- coef(probit.h3)
##         beta1 <- coef(probit.h1)
##         bX3 <- as.numeric(crossprod(beta3, t(X3)))
##         Phi3 <- pnorm(bX3)
##         bX1 <- as.numeric(crossprod(beta1, t(X1)))
##         P0 <- mean(y > 0)
##         yPP <- y * Phi3
##         lin <- switch(dist,
##                       "ln" = lm(log(yPP)  ~ X2 - 1, subset = y!= 0),
##                       "n" = lm(yPP       ~ X2 - 1, subset = y!= 0),
##                       "tn" = truncreg(yPP ~ X2 - 1, subset = y!= 0, scaled = TRUE)
##                       )
##         if(dist != "tn")
##             start <- c(beta1, coef(lin), beta3, sigma = summary(lin)$sigma)
##         else
##             start <- c(beta1, coef(lin)[- length(coef(lin))],
##                        beta3, sigma = coef(lin)[length(coef(lin))])
##     }
##     start
## }


start.mhurdle <- function(X1, X2, X3, y, dist){
    if (!is.null(X1)) beta1 <- coef(glm( (y != 0) ~ X1 - 1, family = binomial(link = 'probit'))) else beta1 <- NULL
    if (!is.null(X3)) beta3 <- coef(glm( (y != 0) ~ X3 - 1, family = binomial(link = 'probit'))) else beta3 <- NULL
    if (dist == "n") lin <- lm(y ~ X2 - 1)
    if (dist == "tn") lin <- truncreg(y ~ X2 - 1, subset = y > 0)
    if (dist == "ln") lin <- lm(y ~ X2 - 1, subset = y > 0)
    beta2 <- coef(lin[1:ncol(X2)])
    if (dist %in% c("n", "ln")) sigma <- summary(lin)$sigma
    else sigma <- coef(lin)[ncol(X2) + 1]
    c(beta1, beta2, beta3, sigma)
}
