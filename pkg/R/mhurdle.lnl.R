# This function computes the log-likelihood for the hurdles models and
# returns also as attributes, if required, the gradient, the fitted
# values and the score vector.

mhurdle.lnl <- function(param, X1, X2, X3, y, gradient = FALSE,
                   fit = FALSE, score = FALSE,
                   dist = NULL, corr = NULL){
  ####
  #  Extract the elements of the model
  ####
  ATAN <- FALSE
  
  h1 <- !is.null(X1) ;  K1 <- ifelse(is.null(X1), 0, ncol(X1))
  h3 <- !is.null(X3) ;  K3 <- ifelse(is.null(X3), 0, ncol(X3))

  K2 <- ncol(X2)
  beta2 <- param[(K1 + 1):(K1 + K2)]
  bX2 <- as.numeric(crossprod(t(X2), beta2))

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
#  if (sigma < 0) sigma <- 0.01  
  Phi2 <- pnorm(bX2 / sigma)
  phi2 <- dnorm(bX2 / sigma)

  rho1 <- rho3 <- 0

  if (!is.null(corr)){
    if (corr == 'h1'){
      if (ATAN){
        eta <- param[K1 + K2 + K3 + 2]
        rho1 <- atan(eta) * 2 / pi
      }
      else{
        rho1 <- param[K1 + K2 + K3 + 2]
        if (rho1 < -1) rho1 <- - 0.99
        if (rho1 >  1) rho1 <-   0.99
      }
    }
    if (corr == 'h3'){
      if (ATAN){
        eta <- param[K1 + K2 + K3 + 2]
        rho3 <- atan(eta) * 2 / pi
      }
      else{
        rho3 <- param[K1 + K2 + K3 + 2]
        if (rho3 < -1) rho3 <- - 0.99
        if (rho3 >  1) rho3 <-   0.99
      }
    }
  }
  Phi12 <- mypbivnorm(bX1, bX2 / sigma, rho1)
  Phi23 <- mypbivnorm(bX2 / sigma, bX3, rho3)
  
  if (dist == "l") resid <- log2(y) + log(Phi3) - bX2
  else resid <- y * Phi3 - bX2
  z <- function(x, rho){
    if (is.null(x)) result <- list(f = 100, g = 0, h = 0)
    else{
      f <- (x + rho / sigma * resid) / sqrt(1 - rho ^ 2) 
      g <- x * rho * (1 - rho ^ 2) ^ - 1.5 + ((1 - rho ^ 2) ^ - 0.5 + rho ^ 2 * (1 - rho ^ 2) ^ - 1.5) * resid / sigma
      h <- x * ((1 - rho ^ 2) ^ -1.5 + 3 * rho ^ 2 * (1 - rho ^ 2) ^ -2.5) +
        resid / sigma * (3 * rho * (1 - rho ^ 2) ^ -1.5 + 3 * rho ^ 3 * (1 - rho ^ 2) ^ -2.5)
      result <- list(f = f, g = g, h = h)
    }
    result
  }
  z1 <- z(bX1, rho1)
  z3 <- z(bX3, rho3)

  lnL.null <- switch(dist,
                     "t" = log(1 - Phi12$f * Phi23$f / Phi2 ^ 2),
                     "n" = log(1 - Phi12$f * Phi23$f / Phi2),
                     "l" = log(1 - Phi1 * Phi3)
                     )
  lnL.pos <-
    - log(sigma) +
      dnorm(resid / sigma, log = TRUE) +
      pnorm(z1$f, log.p = TRUE) +
      pnorm(z3$f, log.p = TRUE) +
      (dist == "l") * (- log2(y)) +
      (dist == "t") * (- log(Phi2) + log(Phi3)) +
      (dist == "n") * log(Phi3)
  lnL <- lnL.null * (y == 0) + lnL.pos * (y != 0)

  if (gradient){
    if (h1){
      lnL.beta1 <- (y == 0)*(switch(dist,
                      "t"  = - (Phi12$a * Phi23$f)/(Phi2 ^ 2 - Phi12$f * Phi23$f),
                      "n"  = - (Phi12$a * Phi23$f)/(Phi2 - Phi12$f * Phi23$f),
                      "l"  = - phi1 * Phi3 / (1 - Phi1 * Phi3) ) ) +
                   (y != 0) * (mills(z1$f) / sqrt(1 - rho1 ^ 2))
      gradi <- lnL.beta1 * X1
    }
    else gradi <- c()
    
    lnL.beta2 <- (y == 0)*(switch(dist,
                    "t" =  (2 * Phi2 * phi2 - Phi12$b * Phi23$f - Phi12$f * Phi23$a) /
                    (Phi2 ^ 2 - Phi12$f * Phi23$f) / sigma - 2 * mills(bX2 / sigma) / sigma,
                    "n" =  (phi2 - Phi12$b * Phi23$f - Phi12$f * Phi23$a) /
                    (Phi2 - Phi12$f * Phi23$f) / sigma - mills(bX2 / sigma) / sigma,
                    "l" = 0 ) ) +
                 (y != 0) * (
                    resid / sigma ^ 2 -  mills(z1$f) * rho1 / sigma / sqrt(1 - rho1 ^ 2) -
                    mills(z3$f) * rho3 / sigma / sqrt(1 - rho3 ^ 2) -
                    (dist == "t") * mills(bX2 / sigma) / sigma
                    )
    gradi <- cbind(gradi, lnL.beta2 * X2)
    
    if (h3){
      deS <- (dist == "l") * mills(bX3) + (dist != "l") * phi3 * y
      lnL.beta3 <- (y == 0) * (switch(dist,
                      "t" = - (Phi12$f * Phi23$b)/(Phi2 ^ 2 - Phi12$f * Phi23$f),
                      "n" = - (Phi12$f * Phi23$b)/(Phi2 - Phi12$f * Phi23$f),
                      "l" = - Phi1 * phi3 / (1 - Phi1 * Phi3) ) ) +
                   (y != 0) * (
                      - resid / sigma^2 * deS
                      + mills(z1$f) * rho1 / sqrt(1 - rho1 ^ 2) / sigma * deS
                      + mills(z3$f) * (1 + rho3 / sigma * deS) / sqrt(1 - rho3 ^ 2)
                      + (dist != 'l') * mills(bX3)
                      )
      gradi <- cbind(gradi, lnL.beta3 * X3)
    }
    
    lnL.sigma <- (y == 0)*(
                    switch(dist,
                           "t" =  - (2 * Phi2 * phi2 - Phi12$b * Phi23$f - Phi12$f * Phi23$a) /
                           (Phi2^2 - Phi12$f * Phi23$f) * bX2 / sigma^2 + 2 * mills(bX2 / sigma) * bX2 / sigma^2,
                           "n" =  - (phi2 - Phi12$b * Phi23$f - Phi12$f * Phi23$a) /
                           (Phi2 - Phi12$f * Phi23$f) * bX2 / sigma^2 + mills(bX2 / sigma) * bX2 / sigma^2,
                           "l" = 0 ) ) +
                 (y != 0) * (
                    - 1 / sigma + resid ^ 2 / sigma ^ 3 - mills(z1$f) * rho1 / sqrt(1 - rho1 ^ 2) * resid / sigma ^ 2 -
                    mills(z3$f) * rho3 / sqrt(1 - rho3 ^ 2) * resid / sigma ^ 2 +
                    (dist == "t") * mills(bX2 / sigma) * bX2 / sigma ^ 2
                    )
    gradi <- cbind(gradi, lnL.sigma)
    
    if (h1 && !is.null(corr) && corr == 'h1'){
      lnL.rho1 <- (y == 0) * (
                     switch(dist,
                            "t" = - Phi12$rho * Phi23$f / (Phi2^2 - Phi12$f * Phi23$f),
                            "n" = - Phi12$rho * Phi23$f / (Phi2   - Phi12$f * Phi23$f),
                            "l" = 0
                            )
                     ) +
                  (y != 0) * (
                     mills(z1$f) * z1$g
                     )
      if (ATAN) gradi <- cbind(gradi, lnL.rho1 * 2 / pi / (1 + eta ^ 2))
      else gradi <- cbind(gradi, lnL.rho1)
    
    }
    if (h3 && !is.null(corr) && corr == 'h3'){
      lnL.rho3 <- (y == 0) * (
                     switch(dist,
                            "t" = - Phi12$f * Phi23$rho / (Phi2^2 - Phi12$f * Phi23$f),
                            "n" = - Phi12$f * Phi23$rho / (Phi2 - Phi12$f * Phi23$f),
                            "l" = 0
                            )
                     ) +
                  (y != 0) * (
                     mills(z3$f) * z3$g
                     )
      if (ATAN) gradi <- cbind(gradi, lnL.rho3 * 2 / pi / (1 + eta ^ 2))
      else gradi <- cbind(gradi, lnL.rho3)
    }
    attr(lnL, "gradient") <- gradi
    attr(lnL, "score") <- c(prime = - 1, second = -1)
  }
  lnL
}

# Compute the estimation of hurdle models in the cases where it can be
# done using two independent estimations (a binomial logit model and a
# normal/log-normal/truncated linear model). This is relevant for
# uncorrelated models with selection

fit.simple.mhurdle <- function(X1, X2, y, dist = NULL){
  probit <- glm(y != 0 ~ X1 - 1, family = binomial(link = "probit"))
  lin <- switch(dist,
                "l" = lm(log(y) ~ X2 - 1, subset = y != 0),
                "n" = lm(y ~ X2 - 1, subset = y != 0),
                "t" = truncreg(y ~ X2 - 1, subset = y != 0)
                )
  df <- df.residual(lin)
  np <- sum(y != 0)
  K1 <- ncol(X1)
  beta1 <- coef(probit)
  bX1 <- as.numeric(crossprod(beta1, t(X1)))
  K2 <- ncol(X2)
  if (dist == "t"){
    sigma <- coef(lin)[ncol(X2)+1]
    beta2 <- coef(lin)[-(ncol(X2)+1)]
  }
  else beta2 <- coef(lin)
  bX2 <- as.numeric(crossprod(beta2, t(X2)))
  L.null <- (y == 0) * log(1 - pnorm(bX1))
  if (dist == "l"){
    logy <- rep(0, length(y))
    logy[y != 0] <- log(y[y != 0])
    resid <- (logy - bX2)
  }
  else resid <- y - bX2
  scr <- sum(resid[y != 0] ^ 2)

  if (dist != "t") sigma <- sqrt(scr / np)
  
  mills1 <- mills(bX1)
  mills2 <- mills(bX2 / sigma) / pnorm(bX2 / sigma)
  mills1m <- mills(- bX1)
  
  L.pos <- switch(dist,
                  "l" = (y != 0) * (- logy + pnorm(bX1, log.p = TRUE) + dnorm(resid / sigma, log = TRUE) - log(sigma)),
                  "n" = (y != 0) * (pnorm(bX1, log.p = TRUE) + dnorm(resid / sigma, log = TRUE) - log(sigma)),
                  "t" = (y != 0) * (pnorm(bX1, log.p = TRUE)
                           + dnorm(resid / sigma, log = TRUE) - log(sigma)
                           - pnorm(bX2 / sigma, log.p = TRUE))
                  )
  gbX1 <- switch(dist,
                "l" = (y == 0) * (- mills1m) + (y != 0) * mills1,
                "n" = (y == 0) * (- mills1m) + (y != 0) * mills1,
                "t" = (y == 0) * (- mills1m) + (y != 0) * mills1
                )

  gbX2 <- switch(dist,
                "l" = (y != 0) * (resid / sigma^2),
                "n" = (y != 0) * (resid / sigma^2),
                "t" = (y != 0) * (resid / sigma^2 - 1 / sigma * mills2)
                ) 

  gsigma <- switch(dist,
                "l" = (y != 0) * (resid^2 / sigma^3 - 1/ sigma),
                "n" = (y != 0) * (resid^2 / sigma^3 - 1/ sigma),
                "t" = (y != 0) * (resid^2 / sigma^3 - 1/ sigma + bX2 / sigma^2 * mills2)
                )

  gradi <- cbind(gbX1 * X1, gbX2 * X2, as.numeric(gsigma))
  dss <- - 3 * scr / sigma ^ 4 + sum(y != 0) / sigma ^ 2

  if (dist == "t"){
    vcov <- bdiag(vcov(probit), vcov(lin))
    coef <- c(coef(probit), coef(lin))
  }
  else{    
    vcov <- bdiag(vcov(probit), vcov(lin) / np * df, - 1 / dss)
    coef <- c(coef(probit), coef(lin), sigma)
  }
  lnL <- sum(L.null + L.pos)
  attr(lnL,"df") <- length(coef)
  fit <- cbind(zero = bX1, pos = bX2)

  other.coef <- c("sigma")

  coef.names <- list(h1    = colnames(X1),
                     h2    = colnames(X2),
                     other = other.coef)
  
  rho <- attr(mhurdle.lnl(c(coef, 0), X1, X2, X3 = NULL, y, score = TRUE,
                     dist = dist),
              "score")
  
  fitted <- compute.fitted.mhurdle(coef, X1, X2, X3 = NULL, dist,
                                   corr = FALSE)

  logLik <- L.null + L.pos
  attr(logLik, "y") <- y
  result <- list(coefficients = coef, 
                 vcov = vcov,
                 fitted.values = fitted,
                 logLik = logLik,
                 gradient = apply(gradi, 2, sum),
                 model = NULL,
                 formula = NULL,
                 coef.names = coef.names,
                 call = NULL,
                 rho = rho
                 )
  
  class(result) <- c("mhurdle","maxLik")
  result
}

# Compute the likelihood model, i.e. a model with no explanatory
# variables.

# Full version with correlation

lnl.naive <- function(param, dist = c("l", "t", "n"), moments,
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
  scr <- ifelse(dist == "l",
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
                 "l" = ym + log(Phi3) - alpha2 ,
                 Phi3 * (ym - alpha2/Phi3)
                 )
    millso <- dnorm(zo)/pnorm(zo)
    Phi1bis <-(alpha1 + rho/sigma * zo)/sqrt(1 - rho^2)
    s2term <- 0.5 * s2 * (rho / (sigma * sqrt(1 - rho^2) * Phi3^(dist != "l")))^2 *
      millso * (zo + millso)
  }
  P0 <- switch(dist,
               "l" = 1 - Phi1 * Phi3,
               "t" = 1 - Pbiv/Phi2 * Phi3,
               "n" = 1 - Pbiv * Phi3
               )
  lnPos <-
    -log(sigma) - 0.5 * log(2*pi) -
      scr/(2*sigma^2) +
        log(Phi3) +
          (log(Phi1bis)+s2term) * h1 -
            ym * (dist == "l") +
              log(Phi3) * (dist != "l") -
                log(Phi2) * (dist == "t")
  switch(which,
         "all" = n * log(P0) + (1 - n) * lnPos,
         "zero" = n * log(P0),
         "positive" = (1 - n) * lnPos)
}

# Version without correlation
lnl.naive <- function(param, dist = c("l", "t", "n"), moments,
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
  scr <- ifelse(dist == "l",
                s2 + (ym + log(Phi3) - alpha2)^2,
                Phi3^2*(s2 + (ym - alpha2/Phi3)^2)
                )
  Pbiv <- Phi1 * Phi2
  Phi1bis <- Phi1
  s2term <- 0
  P0 <- switch(dist,
               "l" = 1 - Phi1 * Phi3,
               "t" = 1 - Pbiv/Phi2 * Phi3,
               "n" = 1 - Pbiv * Phi3
               )
  lnPos <-
    -log(sigma) - 0.5 * log(2*pi) -
      scr/(2*sigma^2) +
        log(Phi3) +
          (log(Phi1bis)+s2term) * h1 -
            ym * (dist == "l") +
              log(Phi3) * (dist != "l") -
                log(Phi2) * (dist == "t")

  n * log(P0) + (1 - n) * lnPos
}
