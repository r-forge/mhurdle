ml.tot <- function(param, X, S, y, gradient = FALSE,
                   fit = FALSE, score = FALSE,
                   dist = NULL, res = FALSE,
                   sel = FALSE, ifr = FALSE, corr = FALSE){
  n <- length(y)
  KS <- ncol(S)
  KX <- ncol(X) 
  betaS <- param[1:KS]
  betaX <- param[(KS+1):(KS+KX)]
  sigma <- param[KS+KX+1]
  if (corr) rho <- param[length(param)] else rho <- 0
  bS <- as.numeric(crossprod(t(S),betaS))
  bX <- as.numeric(crossprod(t(X),betaX))
  PhiX <- pnorm(bX/sigma)
  phiX <- dnorm(bX/sigma)
  PhiS <- pnorm(bS)

  Phib <- pbivnorm(bS,bX/sigma,rho)
  Phib.grad <- attr(Phib,"gradient")
  attr(Phib, "gradient") <- NULL
  Phib.S <- Phib.grad$x1
  Phib.X <- Phib.grad$x2
  Phib.rho <- Phib.grad$rho
  Phib.rho.rho <- Phib.grad$rho.rho

  if (dist == "l") resid <- log2(y) + ifr*log(PhiS) - bX
  else resid <- y*PhiS^ifr - bX

  zc1 <- (bS+rho/sigma*resid)/sqrt(1-rho^2) # argument de la proba cond 1|2
  zc1.rho <- bS*rho*(1-rho^2)^-1.5+((1-rho^2)^-0.5+rho^2*(1-rho^2)^-1.5)*resid/sigma
  zc1.rho.rho <- bS*((1-rho^2)^-1.5+3*rho^2*(1-rho^2)^-2.5)+
    resid/sigma*(3*rho*(1-rho^2)^-1.5+3*rho^3*(1-rho^2)^-2.5)
  
  if (dist == "n") dist <- ifelse(res, "nd", "ns")

  lnL.null <- switch(dist,
                     "t" = log(1-Phib/PhiX),
                     "nd" = log(1-Phib),
                     log(1-PhiS)
                     )

  lnL.pos <-
    -log(sigma)+
      ldnorm(resid/sigma)+
        log(pnorm(zc1))-
          (dist == "l")*log2(y)-
            (dist == "t")*log(PhiX)+
              (dist != "l")*ifr*log(PhiS)

  

  lnL <- lnL.null*(y==0)+lnL.pos*(y!=0)
  # lnL is a vector of lenght n
  
  lnL.beta1 <-
    (y==0)*(
       switch(dist,
              "t" = -Phib.S/PhiX/(1-Phib/PhiX),
              "nd" = - Phib.S/(1-Phib),
              - mills(-bS)
              )
       )+
         (y!=0)*(
            (mills(zc1)/sqrt(1-rho^2)*
             (1+
              ifr*rho/sigma*(
                                (dist == "l")*mills(bS)+
                                (dist != "l")*dnorm(bS)*y
                                )
              )
             )+
             ifr*(-resid/sigma^2*
                     ( (dist != "l")*y*dnorm(bS) + (dist == "l")*mills(bS))+
                     (dist != "l")*mills(bS))
            )
           

  lnL.beta2 <-
    (y==0)*(
       switch(dist,
              "l" = 0,
              "t" =  (phiX/sigma-Phib.X/sigma)/(PhiX-Phib) - mills(bX/sigma)/sigma,
              "ns" = 0,
              "nd" = - Phib.X/(1-Phib)/sigma
              )
       )+
         (y!=0)*(
            resid/sigma^2-mills(zc1)*rho/sigma/sqrt(1-rho^2)-
            (dist=="t")*1/sigma*mills(bX/sigma)
            )
  
  lnL.sigma <-
    (y==0)*(
       switch(dist,
              "t" = (-phiX*bX/sigma^2+Phib.X*bX/sigma^2)/(PhiX-Phib) +
              mills(bX/sigma)*bX/sigma^2,
              "nd" = Phib.X/(1-Phib)*bX/sigma^2,
              0
              )
       )+
         (y!=0)*(
            -1/sigma+resid^2/sigma^3-
            rho*resid/sqrt(1-rho^2)/sigma^2*mills(zc1)+
            (dist=="t")*mills(bX/sigma)*bX/sigma^2
            )

  lnL.rho <-
    (y==0)*(
       switch(dist,
              "t" = -Phib.rho/(PhiX-Phib),
              "nd" = -Phib.rho/(1-Phib),
              0
              )
       )+
         (y!=0)*(
            mills(zc1)*zc1.rho
            )
  
  lnL.rho.rho <-
    (y==0)*(
       switch(dist,
              "l" = 0,
              "t" = -Phib.rho.rho/(PhiX-Phib)-(Phib.rho/(PhiX-Phib))^2,
              "ns" = 0,
              "nd" = -Phib.rho.rho/(1-Phib)-(Phib.rho/(1-Phib))^2,
              )
       )+
         (y!=0)*(
            -mills(zc1)*(zc1+mills(zc1))*zc1.rho^2+mills(zc1)*zc1.rho.rho
            )
  
  if (score) attr(lnL,"score") <- c(prime = sum(lnL.rho), second = sum(lnL.rho.rho))
  if (fit) attr(lnL,"fitted") <- cbind(zero = bS, pos = bX)
  
  if (gradient){
    gradi <- cbind(lnL.beta1*S, lnL.beta2*X, lnL.sigma)
    if (corr) gradi <- cbind(gradi, lnL.rho)
    attr(lnL,"gradient") <- gradi
  }
  lnL
}


sel.simple <- function(X, S, y, dist = NULL, ...){
  probit <- glm(y != 0~S-1, family = binomial(link = "probit"))
  lin <- switch(dist,
                "l" = lm(log(y)~X-1,subset=y!=0),
                "n" = lm(y~X-1,subset=y!=0),
                "t" = truncreg(y~X-1,subset = y!=0)
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
  else{
    betaX <- coef(lin)
  }
  
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

  coef.names <- list(sel = colnames(S),
                     reg = colnames(X),
                     other = "sigma")

  rho <- attr(ml.tot(c(coef,0), X, S, y, score = TRUE,
                          dist = dist, res = FALSE),
             "score")

  fitted <- compute.fitted.mhurdle(coef, X, S, dist,
                                   res = FALSE, sel = TRUE, ifr = FALSE, corr = FALSE)

  result <- list(coefficients = coef,
                 vcov = vcov,
                 fitted.values = fit,
                 logLik = lnL,
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

