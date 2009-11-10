ml.tot <- function(param, X, S, P, y, gradient = FALSE,
                   fit = FALSE, score = FALSE,
                   dist = NULL, corr = FALSE){

  ####
  #  Extract the elements of the model
  ####
  KS <- ifelse(is.null(S), 0, ncol(S))
  KP <- ifelse(is.null(P), 0, ncol(P))  
  KX <- ncol(X)
  sel <- KS > 0
  ifr <- KP > 0
  if (sel) betaS <- param[1:KS] else betaS <- NULL
  betaX <- param[(KS+1):(KS+KX)]
  if (ifr) betaP <- param[(KX+KS+1):(KX+KS+KP)] else betaP <- NULL
  sigma <- param[KX+KS+KP+1]
  if (corr) rho <- param[KX+KS+KP+2] else rho <- 0

  ####
  #  Compute intermediate variables for the log-likelihood
  ####
  bX <- as.numeric(crossprod(t(X),betaX))
  PhiX <- pnorm(bX/sigma)
  phiX <- dnorm(bX/sigma)
  if (!is.null(S)){
    bS <- as.numeric(crossprod(t(S),betaS))
    PhiS <- pnorm(bS)
    phiS <- dnorm(bS)
    Phib <- pbivnorm(bS,bX/sigma,rho)
    Phib.grad <- attr(Phib,"gradient")
    attr(Phib, "gradient") <- NULL
    Phib.S <- Phib.grad$x1
    Phib.X <- Phib.grad$x2
    Phib.rho <- Phib.grad$rho
    Phib.rho.rho <- Phib.grad$rho.rho
  }
  else{
    Phib <- PhiX
    Phib.X <- phiX
    Phib.rho <- 0
    Phib.rho.rho <- 0
    bS <- 0
    PhiS <- 1
  }
  if (!is.null(P)){
    bP <- as.numeric(crossprod(t(P),betaP))
    PhiP <- pnorm(bP)
    phiP <- dnorm(bP)
  }
  else{
    bP <- 0
    PhiP <- 1
    phiP <- 0
  }

  ####
  #  Compute the log-likelihood
  ####
  if (dist == "l") resid <- log2(y) + log(PhiP) - bX
  else resid <- y*PhiP - bX
  # argument de la proba cond 1|2
  zc1 <- (bS + rho/sigma*resid)/sqrt(1-rho^2) 
  zc1.rho <- bS*rho*(1-rho^2)^-1.5+((1-rho^2)^-0.5+rho^2*(1-rho^2)^-1.5)*resid/sigma
  zc1.rho.rho <- bS*((1-rho^2)^-1.5+3*rho^2*(1-rho^2)^-2.5)+
    resid/sigma*(3*rho*(1-rho^2)^-1.5+3*rho^3*(1-rho^2)^-2.5)
  lnL.null <- switch(dist,
                     "t" = log(1-Phib/PhiX*PhiP),
                     "n" = log(1-Phib*PhiP),
                     "l" = log(1-PhiS*PhiP)
                     )
  lnL.pos <-
    -log(sigma)+
      ldnorm(resid/sigma)+
        sel*log(pnorm(zc1))+
          log(PhiP)+
            (dist == "l")*(-log2(y))+
              (dist == "t")*(-log(PhiX)+log(PhiP))+
                (dist == "n")*log(PhiP)
  lnL <- lnL.null*(y==0)+lnL.pos*(y!=0)

  ####
  #  Compute the gradient if required
  ####

  if (gradient){
    if (sel){
      lnL.betaS <-
        (y==0)*(
           switch(dist,
                  "t"  = -Phib.S*PhiP/PhiX/(1-Phib*PhiP/PhiX),
                  "n" = - Phib.S*PhiP/(1-Phib*PhiP),
                  "l"  = - phiS*PhiP/(1-PhiS*PhiP)
                  )
           )+
             (y!=0)*(
                sel*(mills(zc1)/sqrt(1-rho^2)))
      gradi <- lnL.betaS*S
    }
    else gradi <- c()

    lnL.betaX <-
      (y==0)*(
         switch(dist,
                "t" =  (phiX/sigma-Phib.X*PhiP/sigma)/(PhiX-Phib*PhiP) - mills(bX/sigma)/sigma,
                "n" = - Phib.X*PhiP/(1-Phib*PhiP)/sigma,
                "l" = 0
                )
         )+
           (y!=0)*(
              resid/sigma^2-sel*mills(zc1)*rho/sigma/sqrt(1-rho^2)-
              (dist=="t")*1/sigma*mills(bX/sigma)
              )
    gradi <- cbind(gradi, lnL.betaX*X)
    if (ifr){
      lnL.betaP <-
        (y==0)*(
           switch(dist,
                  "t" = - Phib*phiP/PhiX/(1-Phib*PhiP/PhiX),
                  "n" = - Phib*phiP/(1-Phib*PhiP),
                  "l" = - PhiS*phiP/(1-PhiS*PhiP)
                  )
           )+
             (y!=0)*(sel*mills(zc1)/sqrt(1-rho^2)*(rho/sigma*(
                (dist=="l")*mills(bP)+
                (dist!="l")*dnorm(bP)*y
                )
                )-
                resid/sigma^2*(
                               (dist=="l")*mills(bP)+
                               (dist!="l")*dnorm(bP)*y
                               )+
                mills(bP)+
                (dist!="l")*mills(bP)
                )
      gradi <- cbind(gradi, lnL.betaP*P)
    }
    
    lnL.sigma <- (y==0)*(
                    switch(dist,
                           "t" = (-phiX*bX/sigma^2+Phib.X*PhiP*bX/sigma^2)/(PhiX-Phib*PhiP) +
                           mills(bX/sigma)*bX/sigma^2,
                           "n" = Phib.X*PhiP/(1-Phib*PhiP)*bX/sigma^2,
                           "l" = 0
                           )
                    )+
                      (y!=0)*(
                         -1/sigma+resid^2/sigma^3-
                         sel*rho*resid/sqrt(1-rho^2)/sigma^2*mills(zc1)+
                         (dist=="t")*mills(bX/sigma)*bX/sigma^2
                         )
    gradi <- cbind(gradi, lnL.sigma)
  }                           

  ####
  # rho-derivates : relevant if there is a selection process and if
  # either the gradient or the score are required
  ####
  if (sel && (gradient || score) ){
    lnL.rho <- (y==0)*(
                  switch(dist,
                         "t" = -Phib.rho*PhiP/(PhiX-Phib*PhiP),
                         "n" = -Phib.rho*PhiP/(1-Phib*PhiP),
                         "l" = 0
                         )
                  )+
                    (y!=0)*(
                       sel*mills(zc1)*zc1.rho
                       )
    lnL.rho.rho <-
      (y==0)*(
         switch(dist,
                "t" = -Phib.rho.rho*PhiP/(PhiX-Phib*PhiP)-(Phib.rho*PhiP/(PhiX-Phib*PhiP))^2,
                "n" = -Phib.rho.rho*PhiP/(1-Phib*PhiP)-(Phib.rho*PhiP/(1-Phib*PhiP))^2,
                "l" = 0
                )
         )+
           (y!=0)*sel*(
              -mills(zc1)*(zc1+mills(zc1))*zc1.rho^2+mills(zc1)*zc1.rho.rho
              )
    if (score && !corr) attr(lnL, "score") <- c(prime = sum(lnL.rho), second = sum(lnL.rho.rho))
    if (gradient && corr) gradi <- cbind(gradi, lnL.rho)
  }
  if (gradient) attr(lnL, "gradient") <- gradi
  lnL
}

sel.simple <- function(X, S, y, dist = NULL){
  probit <- glm(y!=0 ~ S - 1, family = binomial(link = "probit"))
  lin <- switch(dist,
                "l" = lm(log(y) ~ X - 1, subset = y!=0),
                "n" = lm(y ~ X - 1, subset = y!=0),
                "t" = truncreg(y ~ X - 1, subset = y!=0)
                )
  df <- df.residual(lin)
  np <- sum(y!=0)
  KS <- ncol(S)
  betaS <- coef(probit)
  bS <- as.numeric(crossprod(betaS,t(S)))
  KX <- ncol(X)
  if (dist == "t"){
    sigma <- coef(lin)[ncol(X)+1]
    betaX <- coef(lin)[-(ncol(X)+1)]
  }
  else betaX <- coef(lin)
  bX <- as.numeric(crossprod(betaX,t(X)))
  L.null <- (y==0)*log(1-pnorm(bS))
  if (dist == "l"){
    logy <- rep(0,length(y))
    logy[y!=0] <- log(y[y!=0])
    resid <- (logy - bX)
  }
  else resid <- y - bX
  scr <- sum(resid[y!=0]^2)
  
  if (dist != "t") sigma <- sqrt(scr/np)
  
  millsS <- dnorm(bS)/pnorm(bS)
  millsX <- dnorm(bX/sigma)/pnorm(bX/sigma)
  
  L.pos <- switch(dist,
                  "l" = (y!=0)*(-logy+log(pnorm(bS))+ldnorm(resid/sigma)-log(sigma)),
                  "n" = (y!=0)*(log(pnorm(bS))+ldnorm(resid/sigma)-log(sigma)),
                  "t" = (y!=0)*(log(pnorm(bS))
                               +ldnorm(resid/sigma)-log(sigma)
                               -log(pnorm(bX/sigma)))
                  )
  gbS <- switch(dist,
                "l" = (y==0)*(-dnorm(bS)/pnorm(-bS))+(y!=0)*(dnorm(bS)/pnorm(bS)),
                "n" = (y==0)*(-dnorm(bS)/pnorm(-bS))+(y!=0)*(dnorm(bS)/pnorm(bS)),
                "t" = (y==0)*(-dnorm(bS)/pnorm(-bS))+(y!=0)*millsS
                )

  gbX <- switch(dist,
                "l" = (y!=0)*(resid/sigma^2),
                "n" = (y!=0)*(resid/sigma^2),
                "t" = (y!=0)*(resid/sigma^2 - 1/sigma*millsX)
                )

  gsigma <- switch(dist,
                "l" = (y!=0)*(resid^2/sigma^3-1/sigma),
                "n" = (y!=0)*(resid^2/sigma^3-1/sigma),
                "t" = (y!=0)*(resid^2/sigma^3-1/sigma+bX/sigma^2*millsX)
                )

  gradi <- cbind(gbS*S, gbX*X, as.numeric(gsigma))
  dss <- -3*scr/sigma^4+sum(y!=0)/sigma^2

  if (dist == "t"){
    vcov <- bdiag(vcov(probit),vcov(lin))
    coef <- c(coef(probit),coef(lin))
  }
  else{    
    vcov <- bdiag(vcov(probit),vcov(lin)/np*df,-1/dss)
    coef <- c(coef(probit),coef(lin),sigma)
  }

  lnL <- sum(L.null+L.pos)
  attr(lnL,"df") <- length(coef)
  fit <- cbind(zero = bS, pos = bX)

  other.coef <- c("sigma")

  coef.names <- list(sel = colnames(S),
                     reg = colnames(X),
                     other = other.coef)
  
  
  rho <- attr(ml.tot(c(coef,0), X, S, P = NULL, y, score = TRUE,
                     dist = dist),
              "score")
  
  fitted <- compute.fitted.mhurdle(coef, X, S, P=NULL, dist,
                                   corr = FALSE)

  logLik <- L.null + L.pos
  attr(logLik, "y") <- y
  result <- list(coefficients = coef,
                 vcov = vcov,
                 fitted.values = fitted,
                 logLik = logLik,
                 gradient = apply(gradi,2,sum),
                 model = NULL,
                 formula = NULL,
                 coef.names = coef.names,
                 call = NULL,
                 rho = rho
                 )
  
  class(result) <- c("mhurdle","maxLik")
  result
}

