verbal <- FALSE
mhurdle.lnl <- function(param, X1, X2, X3, X4, y, gradient = FALSE,
                        fitted = FALSE, dist = NULL, corr = FALSE, robust = TRUE){
    if (robust){
        frho <- function(x) atan(x) * 2 / pi
        grho <- function(x) 2 / pi / (1 +  x ^ 2)
        fmu <- function(x) exp(x)
        gmu <- function(x) exp(x)
        fsd <- function(x) exp(x)
        gsd <- function(x) exp(x)
    }
    else{
        frho <- function(x) x
        grho <- function(x) rep(1, length(x))
        fmu <- function(x) x
        gmu <- function(x) 1
        fsd <- function(x) x
        gsd <- function(x) 1
    }

    myeps <- 1E-07
    myInf <- 1000
    N <- length(y)
    #  dummies for the existing equations hi and number of
    #  coefficients Ki
    h1 <- ! is.null(X1) ;  K1 <- ifelse(h1, ncol(X1), 0)
    h3 <- ! is.null(X3) ;  K3 <- ifelse(h3, ncol(X3), 0)
    h4 <- ! is.null(X4) ;  K4 <- ifelse(h4, ncol(X4), 0)
    # KR is either 1 (h1 or h3) or 3 (h1 and h3)
    KR <- ifelse(corr, h1 + h3 + h1 * h3, 0)
    K2 <- ncol(X2)
    
    # shape and scale parameters for box-cox and ihs transformations
    if (dist %in% c("ln2", "bc", "bc2", "ihs")){
        lambda <- param[K1 + K2 + K3 + 1 + K4 + KR + 1]
        if (dist %in% c("bc2", "ln2")){
            if (dist == "bc2") posmu <- K1 + K2 + K3 + 1 + K4 + KR + 2
            else posmu <- K1 + K2 + K3 + 1 + K4 + KR + 1
            mu <- fmu(param[posmu])
            gradientmu <- gmu(param[posmu])
        }
        if (dist == "bc") mu <- 0
    }
    if (dist %in% c("bc", "bc2")) sgn <- sign(lambda) else sgn <- + 1

    # equation 1
    if (h1){
        beta1 <- param[1:K1]
        bX1 <- as.numeric(crossprod(t(X1), beta1))
        Phi1 <- pnorm(bX1) ; phi1 <- dnorm(bX1)
        if (min(Phi1) < myeps){
            nval <- as.numeric(table(Phi1 < myeps)[2])
            if (verbal) cat(paste(nval, " null values of Phi1\n", sep = ""))
            Phi1[Phi1 < myeps] <- myeps
        }
    }
    else{
        bX1 <- rep(myInf, N) ; beta1 <- NULL
        Phi1 <- 1 ; phi1 <- 0;
    }

    # equation 2
    beta2 <- param[(K1 + 1):(K1 + K2)]
    bX2 <- as.numeric(crossprod(t(X2), beta2))
    Phi2 <- pnorm(bX2)

    # equation 3
    if (h3){
        beta3 <- param[(K1 + K2 + 1):(K1 + K2 + K3)]
        bX3 <- as.numeric(crossprod(t(X3), beta3))
        Phi3 <- pnorm(bX3) ; phi3 <- dnorm(bX3)
        if (min(Phi3) < myeps){
            # ROB
            nval <- table(Phi3 < myeps)
            nval <- as.numeric(table(Phi3[[1]] < 1E-07)[2])
            if (verbal) cat(paste(nval, "null values of Phi3\n", sep = ""))
            Phi3[Phi3 < myeps] <- myeps
        }
    }
    else{
        bX3 <- rep(myInf, N) ; beta3 <- NULL
        Phi3 <- 1 ; phi3 <- 0
    }

    # standard deviation
    sd <- param[K1 + K2 + K3 + 1]
    gradientsd <- gsd(sd)
    sd <- fsd(sd)
    # equation 4
    if (h4){
        beta4 <- param[(K1 + K2 + K3 + 2):(K1 + K2 + K3 + K4 + 1)]
        bX4 <- as.numeric(crossprod(t(X4), beta4))
        pbX4 <- pnorm(bX4)
        sigma <- sd * pbX4
    }
    else{
        sigma <- sd
        pbX4 <- 1
    }
    # ymin is the minimum value of y^*, which is 0 except for bc2 and
    # ln2. For box-cox, the minimum value of y^* is 0, T(0) = T(ymin) = - 1
    # /lambda if labmbda > 0 or -inf if lambda < 0
    T0 <- switch(dist,
                 "bc" = (- 1 / lambda) * (lambda > 0) - myInf * (lambda < 0),
                 "bc2" = (mu ^ lambda - 1) / lambda,
                 "ln" = - myInf,
                 "ln2" = log2(mu),
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

    # correlation coefficients
    if (corr){
        posrho <- (K1 + K2 + K3 + K4 + 2):(K1 + K2 + K3 + K4 + 1 + KR)
        rho <- param[posrho]
        # In case of only one correlation coefficient, use the whole
        # vector with only one non-null component
        if (h1 & ! h3) rho <- c(rho, 0, 0)
        if (! h1 & h3) rho <- c(0, 0, rho)
        gradientrho <- grho(rho)
        rho <- frho(rho)
        if (h1 & h3){
            # In case of a trivariate normal distribution, check the
            # joint relation of the three coefficients
            rho3 <- function(rho) rho[1] * rho[2] + c(-1, 1) *
                sqrt(1 + rho[1] ^ 2 * rho[2] ^ 2 - rho[1] ^ 2 - rho[2] ^ 2)
            if (rho[3] < rho3(rho)[1]) rho[3] <- rho3(rho)[1] + 1E-04 
            if (rho[3] > rho3(rho)[2]) rho[3] <- rho3(rho)[2] - 1E-04
        }
    }
    else rho <- rep(0, 3)

    # Transformation of the dependent variable
    Ty <- switch(dist,
                 "ln" = log2(y) + log(Phi3),
                 "ln2" = log2(y * Phi3 + mu),
                 "bc" = (exp(lambda * log(y * Phi3)) - 1) / lambda,
                 "bc2" = (exp(lambda * log(y * Phi3 + mu)) - 1) / lambda,
                 "ihs" = log(lambda * y * Phi3 + sqrt(1 + (lambda  * y * Phi3) ^ 2)) / lambda,
                 y * Phi3
                 )
    
    # logarithm of the jacobian
    lnJ <- switch(dist,
                  "ln" = - log2(y),
                  "ln2" = log(Phi3) - log(mu + Phi3 * y),
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
    if (dist == "ln2") lnJmu <- - 1 / (mu + Phi3 * y)
    if (dist == "bc2") lnJmu <- (lambda - 1) / (Phi3 * y + mu)
    
    # the  residual of the consumption equation
    resid <- (Ty - bX2)
    # problem with bc and lambda < 0, for y = 0, resid = -inf and
    # log(dnorm(resid)) = -inf
    resid[y == 0] <- 0

    # compute the relevant bivariate and trivariate cumulative normals
    Phi2 <- pnorm((bX2 - T0) / sigma)
    phi2 <- dnorm((bX2 - T0) / sigma)
    Pr123A <- PHI3(bX1, (bX2 - T0   ) / sigma, bX3, rho)
    Pr123B <- PHI3(bX1, (bX2 - Tymax) / sigma, bX3, rho)
    Pr13 <- PHI2((bX1 + rho[1] * resid / sigma) / sqrt(1 - rho[1] ^ 2),
                 (bX3 + rho[3] * resid / sigma) / sqrt(1 - rho[3] ^ 2),
                 (rho[2] - rho[1] * rho[3]) / sqrt(1 - rho[1] ^ 2) / sqrt(1 - rho[3] ^ 2)
                 )
    if (min(Pr13$f) < myeps){
        nval <- as.numeric(table(Pr13$f < 1E-07)[2])
        if (verbal) cat(paste(nval, "null values of Phi13\n"))
        Pr13$f[Pr13$f < myeps] <- myeps
    }
    # PI is the correction of the truncature
    PI <- pnorm( (bX2 - Tymin) / sigma) - pnorm( (bX2 - Tymax) / sigma)
    Numerator <- PI - Pr123A$f + Pr123B$f
    lnL.null <- log(Numerator) - log(PI)
    lnL.null[y != 0] <- 0
    lnL.pos <-
        - log(sigma) +
            dnorm(resid / sigma, log = TRUE) +
                log(Pr13$f) +
                    lnJ - log(PI)
    lnL.pos[y == 0] <- 0
    
    lnL <- lnL.null * (y == 0) + lnL.pos * (y != 0)
    if (any(is.na(lnL) | is.infinite(lnL))){
        warnings("infinite or missing values of lnLi")
    }

    if (gradient){
        gradi <- c()
        
        # derivatives respective to beta1
        if (h1){
            lnL.beta1 <- (y == 0) * ( (- Pr123A$d1 + Pr123B$d1) / Numerator) +
                (y != 0) * ( Pr13$d1 / sqrt(1 - rho[1] ^ 2) / Pr13$f)
            gradi <- cbind(gradi, lnL.beta1 * X1)
        }
        
        # derivatives respective to beta2
        PI2 <-  (dnorm( (bX2 - Tymin) / sigma) - dnorm( (bX2 - Tymax) / sigma) ) / sigma
        lnL.beta2 <- (y == 0) * ( (PI2 - Pr123A$d2 / sigma + Pr123B$d2 / sigma) / Numerator -
                                     PI2 / PI) +
            (y != 0) * ( resid / sigma ^ 2 -
                            ( Pr13$d1 * rho[1] / sigma / sqrt(1 - rho[1] ^ 2) +
                                 Pr13$d2 * rho[3] / sigma / sqrt(1 - rho[3] ^ 2) ) / Pr13$f -
                                     PI2 / PI)
        gradi <- cbind(gradi, lnL.beta2 * X2)
        
        # derivatives respective to beta3
        if (h3){
            # derivative of T(Phi3 y) with bX3
            Ty3 <- switch(dist,
                          "ln" = mills(bX3),
                          "ln2" = y * phi3 / (mu + y * Phi3),
                          "bc" = exp(lambda * log(y * Phi3)) * mills(bX3),
                          "bc2" = exp((lambda - 1) * log(y * Phi3 + mu)) * phi3 * y, # a revoir
                          "ihs" = y * phi3 / sqrt( 1 + (lambda * y * Phi3) ^ 2),
                          y * phi3
                          )
            # derivative of lnJ with bX3
            lnJ3 <- switch(dist,
                           "ln" = 0,
                           "ln2" = mills(bX3) - y * phi3 / (mu + y * Phi3),
                           "bc" = lambda * mills(bX3),
                           "bc2" = (lambda - 1) * phi3 * y / (Phi3 * y + mu) + mills(bX3),
                           "ihs" = - phi3 * Phi3 * lambda ^ 2 * y ^ 2 /
                               (1 + (lambda * y * Phi3) ^ 2) + mills(bX3),
                           mills(bX3)
                          )
            
            lnL.beta3 <- (y == 0) * (- Pr123A$d3 + Pr123B$d3) / Numerator + 
                (y != 0) * (- resid / sigma ^ 2 * Ty3 +
                                ( Pr13$d1 * Ty3 * rho[1] / sigma / sqrt(1 - rho[1] ^ 2) +
                                     Pr13$d2 * (1 + Ty3 * rho[3] / sigma) /
                                         sqrt(1 - rho[3] ^ 2) ) / Pr13$f +
                                             lnJ3
                            )
            gradi <- cbind(gradi, lnL.beta3 * X3)
        }
        
        # derivatives respective to sigma
        PIs <- dnorm( (bX2 - Tymin) / sigma) * (- (bX2 - Tymin) / sigma ^ 2) -
            dnorm( (bX2 - Tymax) / sigma) * (- (bX2 - Tymax) / sigma ^ 2)
        Pr123As <- Pr123A$d2 * (- (bX2 - T0) / sigma ^ 2)
        Pr123Bs <- Pr123B$d2 * (- (bX2 - Tymax) / sigma ^ 2)
        lnL.sigma <- (y == 0) * ( (PIs - Pr123As + Pr123Bs) / Numerator - PIs / PI) +
            (y != 0) * (- 1 / sigma + resid ^ 2 / sigma ^ 3 +
                            ( - Pr13$d1 * rho[1] * resid / sigma ^ 2 / sqrt(1 - rho[1] ^ 2) -
                                 Pr13$d2 * rho[3] * resid / sigma ^ 2 /
                                     sqrt(1 - rho[3] ^ 2))  / Pr13$f -
                                         PIs / PI)
        gradi <- cbind(gradi, sigma = lnL.sigma * pbX4  * gradientsd)
        
        # derivatives respective to beta4
        if (!is.null(X4)){
            gradi <- cbind(gradi,  lnL.sigma * sd * dnorm(bX4) * X4)
        }

        # derivatives respective to rho
        if (corr){
            Drho12 <- (rho[1] * rho[2] - rho[3]) / (1 - rho[1] ^ 2) ^ 1.5  / sqrt(1 - rho[3] ^ 2)
            Drho13 <- 1 / sqrt(1 - rho[1] ^ 2) / sqrt(1 - rho[3] ^ 2)
            Drho23 <- (rho[3] * rho[2] - rho[1]) / (1 - rho[3] ^ 2) ^ 1.5  / sqrt(1 - rho[1] ^ 2)
            lnL.rho12 <- (y == 0) * (- Pr123A$dr[, 1] + Pr123B$dr[, 1]) / Numerator  +
                (y != 0) * ( Pr13$d1 * (resid / sigma + rho[1] / (1 - rho[1] ^ 2) *
                                            (bX1 + rho[1] * resid / sigma) ) / (1 - rho[1] ^ 2) ^ 0.5 +
                                                Pr13$dr * Drho12) / Pr13$f
            lnL.rho13 <- (y == 0) * (- Pr123A$dr[, 2] + Pr123B$dr[, 2]) / Numerator  +
                (y != 0) * ( Pr13$dr * Drho13) / Pr13$f
            lnL.rho23 <- (y == 0) * (- Pr123A$dr[, 3] + Pr123B$dr[, 3]) / Numerator  +
                (y != 0) * ( Pr13$d2 * (resid / sigma + rho[3] / (1 - rho[3] ^ 2) *
                                            (bX3 + rho[3] * resid / sigma) ) / (1 - rho[3] ^ 2) ^ 0.5 +
                                                Pr13$dr * Drho23) / Pr13$f
            # return one or three derivates whether there are two or
            # three equations
            lnL.rho <- c()
            if (h1) lnL.rho <- cbind(lnL.rho, rho12 = lnL.rho12 * gradientrho[1])
            if (h1 & h3) lnL.rho <- cbind(lnL.rho, rho13 = lnL.rho13 * gradientrho[2])
            if (h3) lnL.rho <- cbind(lnL.rho, rho23 = lnL.rho23 * gradientrho[3])
            gradi <- cbind(gradi, lnL.rho)
        }

        # derivative respective to lambda
        if (dist %in% c("bc", "bc2")){
            Tylb <- (log(Phi3 * y + mu) * (Phi3 * y + mu) ^ lambda - Ty) / lambda
            if (mu == 0) T0lb <- ( 1 / lambda ^ 2 * (lambda > 0) + 0 * (lambda < 0))
            else T0lb <- (log(mu) * mu ^ lambda - T0) / lambda
            Tymaxlb <- 0 * (lambda > 0) + (1 / lambda ^ 2) * (lambda < 0)
            PIl <- - sign(lambda) * dnorm( (bX2 + 1 / lambda) / sigma ) / (sigma * lambda ^ 2)
            lnL.lambda <- vector(mode = "numeric", length = length(y))
            
            lnL.lambda[y == 0] <- ( (PIl - Pr123A$d2 * (- T0lb /  sigma) + Pr123B$d2 * (- Tymaxlb / sigma)) / Numerator - PIl / PI)[y == 0]

            
            lnL.lambda[y > 0] <- ( (- resid / sigma ^ 2 + (Pr13$d1 * rho[1] / sqrt(1 - rho[1] ^ 2) / sigma +
                                                               Pr13$d2 * rho[3] / sqrt(1 - rho[3] ^ 2) / sigma) / Pr13$f ) * Tylb +
                                                                   lnJlb - PIl / PI)[y > 0]
            gradi <- cbind(gradi, tr = lnL.lambda)
        }

        if (dist == "bc2"){
            Tymu <- (Phi3 * y + mu) ^ (lambda - 1)
            T0mu <- mu ^ (lambda - 1)
            lnL.mu <- vector(mode = "numeric", length = length(y))
            lnL.mu[y == 0] <- (- Pr123A$d2 * (- T0mu / sigma) / Numerator)[y == 0]
            lnL.mu[y > 0] <- (( - resid / sigma ^ 2 +
                                   (Pr13$d1 * rho[1] / sqrt(1 - rho[1] ^ 2) / sigma +
                                        Pr13$d2 * rho[3] / sqrt(1 - rho[3] ^ 2) / sigma) / Pr13$f
                               ) * Tymu + lnJmu)[y != 0]
            gradi <- cbind(gradi, mu = lnL.mu * gradientmu)
        }

        if (dist == "ihs"){
            Tylb <- (y * Phi3) / lambda / sqrt(1 + (lambda * y * Phi3) ^ 2) - Ty / lambda
            lnL.lambda <- vector(mode = "numeric", length = length(y))
            ## lnL.lambda[y != 0] <- (( -resid / sigma ^ 2 +  mills(z1$f) * rho1 / sigma / sqrt(1 - rho1 ^ 2) +
            ##                             mills(z3$f) * rho3 / sigma / sqrt(1 - rho3 ^ 2) ) * Tylb + lnJlb)[y != 0]
            lnL.lambda[y != 0]
            gradi <- cbind(gradi, lnL.lambda)
        }
        
        if (dist == "ln2"){
            Tymu <- 1 / (mu + Phi3 * y)
            T0mu <- 1 / mu
            lnL.mu <- vector(mode = "numeric", length = length(y))
            lnL.mu[y == 0] <- ( - Pr123A$d2 * (- T0mu / sigma ) / Numerator)[y == 0]
            lnL.mu[y != 0] <- (( - resid / sigma ^ 2 +
                                    (Pr13$d1 * rho[1] / sqrt(1 - rho[1] ^ 2) / sigma +
                                         Pr13$d2 * rho[3] / sqrt(1 - rho[3] ^ 2) / sigma) / Pr13$f
                                ) * Tymu + lnJmu)[y != 0]
            gradi <- cbind(gradi, mu = lnL.mu * gradientmu)
        }
        
        if (any(is.na(gradi))){
            warnings("NA values in the gradient\n")
        }
        attr(lnL, "gradient") <- gradi
    }

    fitted <- FALSE
    lnL
}
