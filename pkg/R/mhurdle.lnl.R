# Version du 10 mai ; pas normalisé
# This function computes the log-likelihood for the hurdles models and
# returns also as attributes, if required, the gradient, the fitted
# values and the score vector.


mhurdle.lnl <- function(param, X1, X2, X3, X4, y, gradient = FALSE,
                        fitted = FALSE, dist = NULL, corr = NULL){
#    print(param)
    frho <- function(x) atan(x) * 2 / pi
    grho <- function(x) 2 / pi / (1 +  x ^ 2)
    fmu <- function(x) exp(x)
    gmu <- function(x) exp(x)
    myInf <- 1000
    N <- length(y)
    #  Extract the elements of the model
    h1 <- !is.null(X1) ;  K1 <- ifelse(is.null(X1), 0, ncol(X1))
    h3 <- !is.null(X3) ;  K3 <- ifelse(is.null(X3), 0, ncol(X3))
    if (is.null(X4)) K4 <- 1 else K4 <- ncol(X4)

    K2 <- ncol(X2)
    if (dist %in% c("bc", "bc2", "ihs")){
        lambda <- param[K1 + K2 + K3 + K4 + (!is.null(corr)) + 1]
        if (dist == "bc2"){
            mu <- fmu(param[K1 + K2 + K3 + K4 + (!is.null(corr)) + 2])
            gradientmu <- gmu(param[K1 + K2 + K3 + K4 + (!is.null(corr)) + 2])
        }
        if (dist == "bc") mu <- 0
    }
    if (dist %in% c("bc", "bc2")) sgn <- sign(lambda) else sgn <- + 1
    if (dist == "ln2") mu <- param[K1 + K2 + K3 + K4 + (!is.null(corr)) + 1]

    # equation 1
    if (h1){
        beta1 <- param[1:K1]
        bX1 <- as.numeric(crossprod(t(X1), beta1))
        Phi1 <- pnorm(bX1) ; phi1 <- dnorm(bX1)
    }
    else{
        bX1 <- rep(myInf, N) ; beta1 <- NULL
        Phi1 <- 1 ; phi1 <- 0;
    }

    # equation 4
    if (is.null(X4)) K4 <- 1 else K4 <- ncol(X4)
    beta4 <- param[(K1 + K2 + K3 + 1):(K1 + K2 + K3 + K4)]
    if (is.null(X4)) sigma <- beta4 else sigma <- as.numeric(exp(crossprod(t(X4), beta4)))
    
    # equation 2
    beta2 <- param[(K1 + 1):(K1 + K2)]
    bX2 <- as.numeric(crossprod(t(X2), beta2))
    # ymin is the minimum value of y^*, which is 0 except for bc2 and
    # ln2. For box-cox, the minimum value of y^* is 0, T(0) = T(ymin) = - 1
    # /lambda if labmbda > 0 or -inf if lambda < 0
    T0 <- switch(dist,
                 "bc" = (- 1 / lambda) * (lambda > 0) - myInf * (lambda < 0),
                 "bc2" = (mu ^ lambda - 1) / lambda,
                 "ln" = - myInf,
                 "ln2" = log(exp(mu)),
                 0)

    Tymin <- switch(dist,
                    "bc" =  (- 1 / lambda) * (lambda > 0) - myInf * (lambda < 0),
                    "bc2" = (- 1 / lambda) * (lambda > 0) - myInf * (lambda < 0),
                    "tn" = 0,
                    - myInf)

    Tymax <- switch(dist,
                    "bc" =  myInf * (lambda > 0) - (1 / lambda) * (lambda < 0),
                    "bc2" = myInf * (lambda > 0) - (1 / lambda) * (lambda < 0),
                    myInf)
                    
    # equation 3
    if (h3){
        beta3 <- param[(K1 + K2 + 1):(K1 + K2 + K3)]
        bX3 <- as.numeric(crossprod(t(X3), beta3))
        Phi3 <- pnorm(bX3) ; phi3 <- dnorm(bX3)
    }
    else{
        bX3 <- rep(myInf, N) ; beta3 <- NULL
        Phi3 <- 1 ; phi3 <- 0
    }
    
    # correlation coefficients
    rho1 <- rho3 <- 0
    if (!is.null(corr)){
        if (corr == 'h1'){
            rho1 <- param[K1 + K2 + K3 + K4 + 1]
            if (rho1 < -1) rho1 <- - 0.99
            if (rho1 >  1) rho1 <-   0.99
        }
        if (corr == 'h3'){
            rho3 <- param[K1 + K2 + K3 + K4 + 1]
            if (rho3 < -1) rho3 <- - 0.99
            if (rho3 >  1) rho3 <-   0.99
        }
    }
    rho <- c(rho1, 0, rho3)

    rho1 <- rho3 <- 0
    if (!is.null(corr)){
        if (corr == 'h1'){
            rho1 <- frho(param[K1 + K2 + K3 + K4 + 1])
        }
        if (corr == 'h3'){
            rho3 <- frho(param[K1 + K2 + K3 + K4 + 1])
        }
        gradientrho <- grho(param[K1 + K2 + K3 + K4 + 1])
    }
    rho <- c(rho1, 0, rho3)

    
    # Compute the trivariate probability with the relevant derivatives
    Phi123 <- function(x, rho){
        Phi1 <- pnorm(x[, 1]) ; phi1 <- dnorm(x[, 1])
        Phi2 <- pnorm(x[, 2]) ; phi2 <- dnorm(x[, 2])
        Phi3 <- pnorm(x[, 3]) ; phi3 <- dnorm(x[, 3])
        if (rho[1] == 0 & rho[3] == 0){
            Pr <- Phi1 * Phi2 * Phi3
            a <-  phi1 * Phi2 * Phi3
            b <-  Phi1 * phi2 * Phi3
            c <-  Phi1 * Phi2 * phi3
            drho <- 0
        }
        else{
            if (rho[1] != 0){
                P2 <-   p2norm(x[, 1], x[, 2], rho[1])
                Pr <-   P2$f * Phi3
                a <-    P2$a * Phi3
                b <-    P2$b * Phi3
                c <-    P2$f * phi3
                drho <- P2$rho * Phi3
            }
            if (rho[3] != 0){
                P2 <-   p2norm(x[, 2], x[, 3], rho[3])
                Pr <-   P2$f * Phi1
                a <-    P2$f * phi1
                b <-    P2$a * Phi1
                c <-    P2$b * Phi1
                drho <- P2$rho * Phi1
            }
        }
    list(f = Pr, a = a, b = b, c = c, rho = drho)
    }
    
    # Transformation of the dependent variable
    Ty <- switch(dist,
                 "ln" = log2(y) + log(Phi3),
                 "ln2" = log2(y * Phi3 + exp(mu)),
                 "bc" = (exp(lambda * log(y * Phi3)) - 1) / lambda,#((y * Phi3) ^ lambda - 1) / lambda,
                 "bc2" = (exp(lambda * log(y * Phi3 + mu)) - 1) / lambda,
                 "ihs" = log(lambda * y * Phi3 + sqrt(1 + (lambda  * y * Phi3) ^ 2)) / lambda,
                 y * Phi3
                 )
    # logarithm of the jacobian
    lnJ <- switch(dist,
                  "ln" = - log2(y),
                  "ln2" = log(Phi3) - log(exp(mu) + Phi3 * y),
                  "bc" = (lambda - 1) * log2(y) + lambda * log(Phi3),
                  "bc2" = (lambda - 1) * log(Phi3 * y + mu) + log(Phi3),
                  "ihs" = - 0.5 * log(1 + (lambda * Phi3 * y) ^ 2) + log(Phi3),
                  log(Phi3)
                  )
    
    # derivative of lnJ respective with lambda
    lnJlb <- switch(dist,
                    "bc" = log2(y) + log(Phi3),
                    "bc2" = log(Phi3 * y + mu),
                    "ihs" = - lambda * y ^ 2 * Phi3 ^ 2 / (1 + (lambda * y * Phi3) ^ 2)
                    )

    # derivative of lnJ respective with mu
    if (dist == "ln2") lnJmu <- - 1 / (exp(mu) + Phi3 * y) * exp(mu)
    if (dist == "bc2") lnJmu <- (lambda - 1) / (Phi3 * y + mu)

    # the  residual of the consumption equation
    resid <- (Ty - bX2)
    # pb with bc and lambda < 0, for y = 0, resid = -inf and log(dnorm(resid)) = -inf
    resid[y == 0] <- 0
    
    z <- function(x, resid, rho){
        f <- (x + rho * resid / sigma) / sqrt(1 - rho ^ 2) 
        g <- x * rho * (1 - rho ^ 2) ^ - 1.5 + ((1 - rho ^ 2) ^ - 0.5 + rho ^ 2 * (1 - rho ^ 2) ^ - 1.5) * resid / sigma
        h <- x * ((1 - rho ^ 2) ^ -1.5 + 3 * rho ^ 2 * (1 - rho ^ 2) ^ -2.5) +
            resid / sigma * (3 * rho * (1 - rho ^ 2) ^ -1.5 + 3 * rho ^ 3 * (1 - rho ^ 2) ^ -2.5)
        result <- list(f = f, g = g, h = h)
    }
    z1 <- z(bX1, resid, rho1)
    z3 <- z(bX3, resid, rho3)
    Phi2 <- pnorm((bX2 - T0) / sigma)
    phi2 <- dnorm((bX2 - T0) / sigma)
    Pr123A <- Phi123(cbind(bX1, (bX2 - T0   ) / sigma, bX3), rho)
    Pr123B <- Phi123(cbind(bX1, (bX2 - Tymax) / sigma, bX3), rho)
    PI <- pnorm( (bX2 - Tymin) / sigma) - pnorm( (bX2 - Tymax) / sigma)

    Denom <- PI - Pr123A$f + Pr123B$f
    lnL.null <- log(Denom) - log(PI)
    lnL.pos <-
        - log(sigma) +
            dnorm(resid / sigma, log = TRUE) +
                    pnorm(z1$f, log.p = TRUE) +
                        pnorm(z3$f, log.p = TRUE) +
                            lnJ - log(PI)

    lnL <- lnL.null * (y == 0) + lnL.pos * (y != 0)
    if (gradient){
        gradi <- c()
        if (h1){
            lnL.beta1 <- (y == 0) * ( (- Pr123A$a + Pr123B$a) / Denom) +
                (y != 0) * (mills(z1$f) / sqrt(1 - rho1 ^ 2))
            gradi <- cbind(gradi, lnL.beta1 * X1)
        }
        
        PI2 <-  (dnorm( (bX2 - Tymin) / sigma) - dnorm( (bX2 - Tymax) / sigma) ) / sigma
        lnL.beta2 <- (y == 0) * ( (PI2 - Pr123A$b / sigma + Pr123B$b / sigma) / Denom - PI2 / PI) +
            (y != 0) * ( resid / sigma ^ 2 - mills(z1$f) * rho1 / sqrt(1 - rho1 ^ 2) / sigma -
                 mills(z3$f) * rho3 / sqrt(1 - rho3 ^ 2) / sigma - PI2 / PI
                 )
        gradi <- cbind(gradi, lnL.beta2 * X2)
        
        if (h3){
            # derivative of T(Phi3 y) with bX3
            Ty3 <- switch(dist,
                          "ln" = mills(bX3),
                          "ln2" = y * phi3 / (exp(mu) + y * Phi3),
                          "bc" = exp(lambda * log(y * Phi3)) * mills(bX3),#(y * Phi3) ^ lambda  * mills(bX3),
                          "bc2" = exp((lambda - 1) * log(y * Phi3 + mu)) * phi3 * y, # a revoir
                          "ihs" = y * phi3 / sqrt( 1 + (lambda * y * Phi3) ^ 2),
                          y * phi3
                          )
            # derivative of lnJ with bX3
            lnJ3 <- switch(dist,
                           "ln" = 0,
                           "ln2" = mills(bX3) - y * phi3 / (exp(mu) + y * Phi3),
                           "bc" = lambda * mills(bX3),
                           "bc2" = (lambda - 1) * phi3 * y / (Phi3 * y + mu) + mills(bX3),
                           "ihs" = - phi3 * Phi3 * lambda ^ 2 * y ^ 2 / (1 + (lambda * y * Phi3) ^ 2) + mills(bX3),
                           mills(bX3)
                          )

            lnL.beta3 <- (y == 0) * (- Pr123A$c + Pr123B$c) / Denom + 
                (y != 0) * (- resid / sigma ^ 2 * Ty3
                     + mills(z1$f) * rho1 / sqrt(1 - rho1 ^ 2) / sigma * Ty3
                     + mills(z3$f) * (1 + rho3 / sigma * Ty3) / sqrt(1 - rho3 ^ 2)
                     + lnJ3
                     )
            gradi <- cbind(gradi, lnL.beta3 * X3)
        }
        
#        PIs <- sgn * dnorm( (bX2 - Tymin) / sigma) * (- (bX2 - Tymin) / sigma ^ 2)
        PIs <- dnorm( (bX2 - Tymin) / sigma) * (- (bX2 - Tymin) / sigma ^ 2) -
            dnorm( (bX2 - Tymax) / sigma) * (- (bX2 - Tymax) / sigma ^ 2)
        
        Pr123As <- Pr123A$b * (- (bX2 - T0) / sigma ^ 2)
        Pr123Bs <- Pr123B$b * (- (bX2 - Tymax) / sigma ^ 2)

        lnL.sigma <- (y == 0) * ( (PIs - Pr123As + Pr123Bs) / Denom - PIs / PI) +
            (y != 0) * (- 1 / sigma + resid ^ 2 / sigma ^ 3 - mills(z1$f) * rho1 / sqrt(1 - rho1 ^ 2) * resid / sigma ^ 2 -
                 mills(z3$f) * rho3 / sqrt(1 - rho3 ^ 2) * resid / sigma ^ 2 - PIs / PI)
                 
        if (!is.null(X4)) lnL.sigma <- lnL.sigma * sigma * X4
        gradi <- cbind(gradi, lnL.sigma)

        if (!is.null(corr)){
            if (h1 & corr == 'h1') lnL.rho <- (y == 0) * ( - Pr123A$rho + Pr123B$rho) / Denom + (y != 0) * mills(z1$f) * z1$g
            if (h3 & corr == 'h3') lnL.rho <- (y == 0) * ( - Pr123A$rho + Pr123B$rho) + (y != 0) * mills(z3$f) * z3$g
            gradi <- cbind(gradi, rho = lnL.rho * gradientrho)
        }
        
        if (dist %in% c("bc", "bc2")){
#            Tylb <- (log(Phi3 * y) * (Phi3 * y) ^ lambda * lambda - Ty) / lambda
            Tylb <- (log(Phi3 * y + mu) * (Phi3 * y + mu) ^ lambda - Ty) / lambda
            if (mu == 0) T0lb <- ( 1 / lambda ^ 2 * (lambda > 0) + 0 * (lambda < 0))
            else T0lb <- (log(mu) * mu ^ lambda - T0) / lambda
            Tymaxlb <- 0 * (lambda > 0) + (1 / lambda ^ 2) * (lambda < 0)
            PIl <- - sign(lambda) * dnorm( (bX2 + 1 / lambda) / sigma ) / (sigma * lambda ^ 2)
            lnL.lambda <- vector(mode = "numeric", length = length(y))
            lnL.lambda[y == 0] <- ( (PIl - Pr123A$b * (- T0lb /  sigma) + Pr123B$b * (- Tymaxlb / sigma)) / Denom - PIl / PI)[y == 0]
            lnL.lambda[y > 0] <- (( - resid / sigma ^ 2 +  mills(z1$f) * rho1 / sigma / sqrt(1 - rho1 ^ 2) +
                                   mills(z3$f) * rho3 / sigma / sqrt(1 - rho3 ^ 2) ) * Tylb + lnJlb - PIl / PI)[y > 0]
            gradi <- cbind(gradi, tr = lnL.lambda)
        }

        if (dist == "bc2"){
            Tymu <- (Phi3 * y + mu) ^ (lambda - 1)
            T0mu <- mu ^ (lambda - 1)
            lnL.mu <- vector(mode = "numeric", length = length(y))
            lnL.mu[y == 0] <- (- Pr123A$b * (- T0mu / sigma) / Denom)[y == 0]
            ## if (lambda > 0) lnL.mu[y == 0] <- (- Pr123$b / (PI - Pr123$f) * (- T0mu / sigma) )[y == 0]
            ## else lnL.mu[y == 0] <- (- Pr123a$b / (PI - Pr123$f) * (- T0mu / sigma) )[y == 0]
            lnL.mu[y > 0] <- (( - resid / sigma ^ 2 +  mills(z1$f) * rho1 / sigma / sqrt(1 - rho1 ^ 2)  +
                               mills(z3$f) * rho3 / sigma / sqrt(1 - rho3 ^ 2) ) * Tymu + lnJmu)[y != 0]
#            lnL.mu[y == 0] <- 0
            #ZA
            gradi <- cbind(gradi, mu = lnL.mu * gradientmu)
        }

        if (dist == "ihs"){
            Tylb <- (y * Phi3) / lambda / sqrt(1 + (lambda * y * Phi3) ^ 2) - Ty / lambda
            lnL.lambda <- vector(mode = "numeric", length = length(y))
            lnL.lambda[y != 0] <- (( -resid / sigma ^ 2 +  mills(z1$f) * rho1 / sigma / sqrt(1 - rho1 ^ 2) +
                           mills(z3$f) * rho3 / sigma / sqrt(1 - rho3 ^ 2) ) * Tylb + lnJlb)[y != 0]
            gradi <- cbind(gradi, lnL.lambda)
        }

        if (dist == "ln2"){
            Tymu <- 1 / (exp(mu) + Phi3 * y) * exp(mu)
            lnL.mu <- vector(mode = "numeric", length = length(y))
            lnL.mu[y == 0] <- ( Pr123$b / Denom * (- 1 / (sigma * exp(mu)) ))[y == 0] * exp(mu)
#            lnL.mu[y == 0] <- ((- sgn * Pr123$b / (1 - sgn * Pr123$f)) * (- 1 / (sigma * exp(mu)) ))[y == 0] * exp(mu)
            lnL.mu[y != 0] <- (( -resid / sigma ^ 2 +  mills(z1$f) * rho1 / sigma / sqrt(1 - rho1 ^ 2) +
                           mills(z3$f) * rho3 / sigma / sqrt(1 - rho3 ^ 2) ) * Tymu + lnJmu)[y != 0]
            gradi <- cbind(gradi, lnL.mu)
#            if (dist == "ln2") print(apply(gradi, 2, sum))
        }
        attr(lnL, "gradient") <- gradi
    }
    if (fitted){
        P0 <- exp(lnL.null)
        if (dist != "ln"){
            if (h3) Psi23 <- rho3 * phi3 * pnorm( (bX2 / sigma - rho3 * bX3) / sqrt(1 - rho3 ^ 2)) +
                phi2 * pnorm( (bX3 - rho3 * bX2 / sigma) / sqrt(1 - rho3 ^ 2))
            else Psi23 <- phi2
            if (h1) Psi21 <- rho1 * phi1 * pnorm( (bX2 / sigma - rho1 * bX1) / sqrt(1 - rho1 ^ 2)) +
                phi2 * pnorm( (bX1 - rho1 * bX2 / sigma) / sqrt(1 - rho1 ^ 2))
            else Psi21 <- phi2
            # PAS MODIFI2 PAR FLEMME !!!!!!!!
            Phi12 <- mypbivnorm(bX1, bX2  / sigma , rho1)
            Phi23 <- mypbivnorm( bX2 / sigma, bX3, rho3)
            Econd <- bX2 / Phi3 + sigma * (Psi23 * Psi21 * Phi2) / (Phi12$f * Phi23$f * Phi3 * phi2)
        }
        else{
            if (h3) Psi23 <- pnorm(bX3 + sigma * rho3) else Psi23 <- 1
            if (h1) Psi21 <- pnorm(bX1 + sigma * rho1) else Psi21 <- 1
            Econd <- exp(bX2 + sigma ^ 2 / 2) * Psi21 * Psi23 / (Phi1 * Phi3 ^ 2)
        }
        attr(lnL, "fitted") <- cbind("P(y=0)" = P0, "E(y|y>0)" = Econd)
    }
    lnL
}

# Compute the estimation of hurdle models in the cases where it can be
# done using two independent estimations (a binomial logit model and a
# normal/log-normal/truncated linear model). This is relevant for
# uncorrelated models with selection

fit.simple.mhurdle <- function(X1, X2, y, dist = NULL){
    probit <- glm(y != 0 ~ X1 - 1, family = binomial(link = "probit"))
    # Computation of the likelihood for zero observations
    beta1 <- coef(probit)
    bX1 <- as.numeric(crossprod(beta1, t(X1)))
    mills1 <- mills(bX1)
    mills1m <- mills(- bX1)
    L.null <- (y == 0) * log(1 - pnorm(bX1))
    gbX1 <- (y == 0) * (- mills1m) + (y != 0) * mills1
    
    # Computation of the likelihood for positive observations for (log-)normal distribution
    if (dist %in% c("ln", "n")){
        if (dist == "ln") lin <- lm(log(y) ~ X2 - 1, subset = y != 0)
        else lin <- lm(y ~ X2 - 1, subset = y != 0)
        beta2 <- coef(lin)[1:ncol(X2)]
        bX2 <- as.numeric(crossprod(beta2, t(X2)))
        logy <- rep(0, length(y))
        logy[y != 0] <- log(y[y != 0])
        if (dist == "ln"){
            resid <- (logy - bX2)
        }
        else resid <- y - bX2
        df <- df.residual(lin)
        np <- sum(y != 0)
        scr <- sum(resid[y != 0] ^ 2)
        sigma <- sqrt(scr / np)
        L.pos <- (y != 0) * (pnorm(bX1, log.p = TRUE) + dnorm(resid / sigma, log = TRUE) - log(sigma) -
                      (dist == "ln") * logy)
        gbX2 <- (y != 0) * (resid / sigma ^ 2)
        gsigma <- (y != 0) * (resid ^ 2 / sigma ^ 3 - 1 / sigma)
        gradi <- cbind(gbX1 * X1, gbX2 * X2, as.numeric(gsigma))
        dss <- - 3 * scr / sigma ^ 4 + sum(y != 0) / sigma ^ 2
        vcov <- bdiag(vcov(probit), vcov(lin) / np * df, - 1 / dss)
        coef <- c(coef(probit), coef(lin), sigma)
        fit <- cbind(zero = bX1, pos = bX2)
    }
    else{
        if (dist == "tn") lin <- truncreg(y ~ X2 - 1, subset = y != 0)
        if (dist == "bc") lin <- boxcoxreg(y ~ X2 - 1, subset = y != 0)
        L.pos <- as.numeric(logLik(lin)) + sum( (y != 0) * (pnorm(bX1, log.p = TRUE)))
        vcov <- bdiag(vcov(probit), vcov(lin))
        coef <- c(coef(probit), coef(lin))
        fit <- cbind(zero = bX1, pos = fitted(lin))
        g2 <- matrix(0, nrow = length(y), ncol = ncol(lin$gradient))
        g2[y != 0, ] <- lin$gradient
        gradi <- cbind(gbX1 * X1, g2)
        vcov <- bdiag(vcov(probit) , vcov(lin))
    }
    
    coef.names <- list(h1    = colnames(X1),
                       h2    = colnames(X2),
                       sd    = "sd")
    if (dist == "bc") coef.names <- c(coef.names, tr = "tr")

    fitted <- attr(mhurdle.lnl(coef, X1 = X1, X2 = X2, X3 = NULL, X4 = NULL, y = y,
                               gradient = FALSE, fitted = TRUE,
                               dist = dist, corr = NULL), "fitted")
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

if (FALSE){
lnl.naive <- function(param, dist = c("ln", "tn", "n"), moments,
                     h1 = TRUE, h3 = FALSE,
                     which = c("all", "zero", "positive")){
  dist <- match.arg(dist)
  which <- match.arg(which)
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
  
  if (length(param) == 2) rho <- param[2] else rho <- 0
  if (rho < - 1) rho <- - 0.99
  if (rho > 1) rho <- 0.99
  rho1 <- rho
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
  if (!rho){
    Pbiv <- Phi1 * Phi2
    Phi1bis <- Phi1
    s2term <- 0
  }
  else{
    Pbiv <- mypbivnorm(alpha1, alpha2/sigma, rho)$f
    zo <- switch(dist,
                 "ln" = ym + log(Phi3) - alpha2 ,
                 Phi3 * (ym - alpha2/Phi3)
                 )
    millso <- dnorm(zo)/pnorm(zo)
    Phi1bis <-(alpha1 + rho/sigma * zo)/sqrt(1 - rho^2)
    s2term <- 0.5 * s2 * (rho / (sigma * sqrt(1 - rho^2) * Phi3^(dist != "ln")))^2 *
      millso * (zo + millso)
  }
  P0 <- switch(dist,
               "ln" = 1 - Phi1 * Phi3,
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
  switch(which,
         "all" = n * log(P0) + (1 - n) * lnPos,
         "zero" = n * log(P0),
         "positive" = (1 - n) * lnPos)
}
}
# Version without correlation

if (TRUE){
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
}
