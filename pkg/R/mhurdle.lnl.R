# This function computes the log-likelihood for the hurdles models and
# returns also as attributes, if required, the gradient, the fitted
# values and the score vector.

mhurdle.lnl <- function(param, X, S, P, y, gradient = FALSE,
                   fit = FALSE, score = FALSE,
                   dist = NULL, corr = NULL){
  ####
  #  Extract the elements of the model
  ####

  sel <- !is.null(S) ;  KS <- ifelse(is.null(S), 0, ncol(S))
  ifr <- !is.null(P) ;  KP <- ifelse(is.null(P), 0, ncol(P))

  KX <- ncol(X)
  betaX <- param[(KS+1):(KS+KX)]
  bX <- as.numeric(crossprod(t(X), betaX))

  if (sel){
    betaS <- param[1:KS]
    bS <- as.numeric(crossprod(t(S), betaS))
    PhiS <- pnorm(bS) ; phiS <- dnorm(bS)
  }
  else{
    bS <- betaS <- NULL
    PhiS <- 1 ; phiS <- 0;
  }

  if (ifr){
    betaP <- param[(KX + KS + 1):(KX + KS + KP)]
    bP <- as.numeric(crossprod(t(P), betaP))
    PhiP <- pnorm(bP) ; phiP <- dnorm(bP)
  }
  else{
    bP <- betaP <- NULL
    PhiP <- 1 ; phiP <- 0
  }
  
  sigma <- param[KX + KS + KP + 1]
  PhiX <- pnorm(bX / sigma)
  phiX <- dnorm(bX / sigma)

  rhoS <- rhoP <- 0
  if (!is.null(corr)){
    if (corr == 'sel') rhoS <- param[KX + KS + KP + 2]
    if (corr == 'ifr') rhoP <- param[KX + KS + KP + 2]
  }

  PhiSX <- pbivnorm(bS, bX / sigma, rhoS)
  PhiXP <- pbivnorm(bX / sigma, bP, rhoP)
  
  if (dist == "l") resid <- log2(y) + log(PhiP) - bX
  else resid <- y*PhiP - bX

  z <- function(x, rho){
    if (is.null(x)) result <- list(f = 100, g = 0, h = 0)
    else{
      f <- (x + rho/sigma*resid)/sqrt(1 - rho^2) 
      g <- x * rho * (1 - rho^2)^ - 1.5 + ((1 - rho^2)^ - 0.5+rho^2*(1 - rho^2)^ - 1.5)*resid/sigma
      h <- x* ((1 - rho^2)^-1.5 + 3 * rho^2 * (1 - rho^2)^-2.5) +
        resid / sigma * (3 * rho * (1 - rho^2)^-1.5 + 3 * rho^3 * (1 - rho^2)^-2.5)
      result <- list(f = f, g = g, h = h)
    }
    result
  }
  zS <- z(bS, rhoS)
  zP <- z(bP, rhoP)
  lnL.null <- switch(dist,
                     "t" = log(1 - PhiSX$f * PhiXP$f / PhiX^2),
                     "n" = log(1 - PhiSX$f * PhiXP$f / PhiX),
                     "l" = log(1 - PhiP * PhiS)
                     )
  lnL.pos <-
    - log(sigma) +
      dnorm(resid / sigma, log = TRUE) +
      pnorm(zS$f, log.p = TRUE) +
      pnorm(zP$f, log.p = TRUE) +
      (dist == "l") * (- log2(y)) +
      (dist == "t") * (- log(PhiX) + log(PhiP)) +
      (dist == "n") * log(PhiP)
  
  lnL <- lnL.null * (y == 0) + lnL.pos * (y != 0)

  if (gradient){
    if (sel){
      lnL.betaS <- (y == 0)*(switch(dist,
                      "t"  = - (PhiSX$a * PhiXP$f)/(PhiX^2 - PhiSX$f * PhiXP$f),
                      "n"  = - (PhiSX$a * PhiXP$f)/(PhiX - PhiSX$f * PhiXP$f),
                      "l"  = - phiS * PhiP / (1 - PhiS * PhiP) ) ) +
                   (y != 0) * (mills(zS$f) / sqrt(1-rhoS^2))
      gradi <- lnL.betaS * S
    }
    else gradi <- c()
    
    lnL.betaX <- (y == 0)*(switch(dist,
                    "t" =  (2 * PhiX * phiX - PhiSX$b * PhiXP$f - PhiSX$f * PhiXP$a) /
                    (PhiX^2 - PhiSX$f * PhiXP$f) / sigma - 2 * mills(bX / sigma) / sigma,
                    "n" =  (phiX - PhiSX$b * PhiXP$f - PhiSX$f * PhiXP$a) /
                    (PhiX - PhiSX$f * PhiXP$f) / sigma - mills(bX / sigma) / sigma,
                    "l" = 0 ) ) +
                 (y != 0) * (
                    resid / sigma^2 -  mills(zS$f) * rhoS / sigma / sqrt(1 - rhoS^2) -
                    mills(zP$f) * rhoP / sigma / sqrt(1 - rhoP^2) -
                    (dist == "t") * mills(bX / sigma) / sigma
                    )
    gradi <- cbind(gradi, lnL.betaX * X)
    
    if (ifr){
      deS <- (dist == "l") * mills(bP) + (dist != "l") * phiP * y
      lnL.betaP <- (y == 0) * (switch(dist,
                      "t" = - (PhiSX$f * PhiXP$b)/(PhiX^2 - PhiSX$f * PhiXP$f),
                      "n" = - (PhiSX$f * PhiXP$b)/(PhiX - PhiSX$f * PhiXP$f),
                      "l" = - PhiS * phiP / (1 - PhiS * PhiP) ) ) +
                   (y != 0) * (
                      - resid / sigma^2 * deS
                      + mills(zS$f) * rhoS / sqrt(1 - rhoS^2) / sigma * deS
                      + mills(zP$f) * (1 + rhoP / sigma * deS) / sqrt(1 - rhoS^2)
                      + (dist != 'l') * mills(bP)
                      )
      gradi <- cbind(gradi, lnL.betaP * P)
    }
    
    lnL.sigma <- (y == 0)*(
                    switch(dist,
                           "t" =  - (2 * PhiX * phiX - PhiSX$b * PhiXP$f - PhiSX$f * PhiXP$a) /
                           (PhiX^2 - PhiSX$f * PhiXP$f) * bX / sigma^2 + 2 * mills(bX / sigma) * bX / sigma^2,
                           "n" =  - (phiX - PhiSX$b * PhiXP$f - PhiSX$f * PhiXP$a) /
                           (PhiX - PhiSX$f * PhiXP$f) * bX / sigma^2 + mills(bX / sigma) * bX / sigma^2,
                           "l" = 0 ) ) +
                 (y != 0) * (
                    - 1 / sigma + resid^2 / sigma^3 - mills(zS$f) *rhoS / sqrt(1 - rhoS^2) * resid / sigma^2 -
                    mills(zP$f) *rhoP / sqrt(1 - rhoP^2) * resid / sigma^2 -
                    (dist == "t") * mills(bX / sigma) * bX / sigma^2
                    )
    gradi <- cbind(gradi, lnL.sigma)
    
    if (sel && !is.null(corr)){
      lnL.rhoS <- (y == 0) * (
                     switch(dist,
                            "t" = - PhiSX$rho * PhiXP$f / (PhiX^2 - PhiSX$f * PhiXP$f),
                            "n" = - PhiSX$rho * PhiXP$f / (PhiX   - PhiSX$f * PhiXP$f),
                            "l" = 0
                            )
                     ) +
                  (y != 0) * (
                     mills(zS$f) * zS$g
                     )
      gradi <- cbind(gradi, lnL.rhoS)
    }
    if (ifr && !is.null(corr)){
      lnL.rhoP <- (y == 0) * (
                     switch(dist,
                            "t" = - PhiSX$f * PhiXP$rho / (PhiX - PhiSX$f * PhiXP$f),
                            "n" = - PhiSX$rho * PhiXP$f / (PhiX^2 - PhiSX$f * PhiXP$f),
                            "l" = 0
                            )
                     ) +
                  (y != 0) * (
                     mills(zP$f) * zP$g
                     )
      gradi <- cbind(gradi, lnL.rhoP)
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

fit.simple.mhurdle <- function(X, S, y, dist = NULL){
  probit <- glm(y != 0 ~ S - 1, family = binomial(link = "probit"))
  lin <- switch(dist,
                "l" = lm(log(y) ~ X - 1, subset = y != 0),
                "n" = lm(y ~ X - 1, subset = y != 0),
                "t" = truncreg(y ~ X - 1, subset = y != 0)
                )
  df <- df.residual(lin)
  np <- sum(y != 0)
  KS <- ncol(S)
  betaS <- coef(probit)
  bS <- as.numeric(crossprod(betaS, t(S)))
  KX <- ncol(X)
  if (dist == "t"){
    sigma <- coef(lin)[ncol(X)+1]
    betaX <- coef(lin)[-(ncol(X)+1)]
  }
  else betaX <- coef(lin)
  bX <- as.numeric(crossprod(betaX, t(X)))
  L.null <- (y == 0) * log(1 - pnorm(bS))
  if (dist == "l"){
    logy <- rep(0, length(y))
    logy[y != 0] <- log(y[y != 0])
    resid <- (logy - bX)
  }
  else resid <- y - bX
  scr <- sum(resid[y != 0]^2)
  
  if (dist != "t") sigma <- sqrt(scr / np)
  
  millsS <- mills(bS)
  millsX <- mills(bX / sigma) / pnorm(bX / sigma)
  millsSm <- mills(- bS)
  
  L.pos <- switch(dist,
                  "l" = (y != 0) * (-logy + pnorm(bS, log.p = TRUE) + dnorm(resid / sigma, log = TRUE) - log(sigma)),
                  "n" = (y != 0) * (pnorm(bS, log.p = TRUE) + dnorm(resid / sigma, log = TRUE) - log(sigma)),
                  "t" = (y != 0) * (pnorm(bS, log.p = TRUE)
                           + dnorm(resid / sigma, log = TRUE) - log(sigma)
                           - pnorm(bX / sigma, log.p = TRUE))
                  )
  gbS <- switch(dist,
                "l" = (y == 0) * (- millsSm) + (y != 0) * millsS,
                "n" = (y == 0) * (- millsSm) + (y != 0) * millsS,
                "t" = (y == 0) * (- millsSm) + (y != 0) * millsS
                )

  gbX <- switch(dist,
                "l" = (y != 0) * (resid / sigma^2),
                "n" = (y != 0) * (resid / sigma^2),
                "t" = (y != 0) * (resid / sigma^2 - 1 / sigma * millsX)
                ) 

  gsigma <- switch(dist,
                "l" = (y != 0) * (resid^2 / sigma^3 - 1/ sigma),
                "n" = (y != 0) * (resid^2 / sigma^3 - 1/ sigma),
                "t" = (y != 0) * (resid^2 / sigma^3 - 1/ sigma + bX / sigma^2 * millsX)
                )

  gradi <- cbind(gbS * S, gbX * X, as.numeric(gsigma))
  dss <- - 3 * scr / sigma^4 + sum(y != 0) / sigma^2

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
  fit <- cbind(zero = bS, pos = bX)

  other.coef <- c("sigma")

  coef.names <- list(sel   = colnames(S),
                     reg   = colnames(X),
                     other = other.coef)
  
  
  rho <- attr(mhurdle.lnl(c(coef, 0), X, S, P = NULL, y, score = TRUE,
                     dist = dist),
              "score")
  
  fitted <- compute.fitted.mhurdle(coef, X, S, P = NULL, dist,
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

lnl.naive <- function(param, dist = c("l", "t", "n"), moments,
                     sel = TRUE, ifr = FALSE,
                     which = c("all", "zero", "positive")){
  dist <- match.arg(dist)
  which <- match.arg(which)
  n <- moments[1]
  ym <- moments[2]
  s2 <- moments[3]
  if (sel){
    alpha1 <- param[1]
    alpha2 <- param[2]
    param <- param[-c(1,2)]
  }
  else{
    alpha2 <- param[1]
    param <- param[-1]
  }
  if (ifr){
    alpha3 <- param[1]
    param <- param[-1]
  }
  sigma <- param[1]
  
  if (length(param) == 2) rho <- param[2] else rho <- 0
  rhoS <- rho
  if (sel){
    Phi1 <- pnorm(alpha1)
    phi1 <- dnorm(alpha1)
  }
  else Phi1 <- 1
  if (ifr){
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
    Pbiv <- pbivnorm(alpha1, alpha2/sigma, rho)
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
          (log(Phi1bis)+s2term) * sel -
            ym * (dist == "l") +
              log(Phi3) * (dist != "l") -
                log(Phi2) * (dist == "t")
  switch(which,
         "all" = n * log(P0) + (1 - n) * lnPos,
         "zero" = n * log(P0),
         "positive" = (1 - n) * lnPos)
}


