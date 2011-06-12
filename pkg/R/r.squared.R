r.squared <- function(object,
                      which = c("all", "zero", "positive"),
                      type = c("regression", "mcfadden"),
                      dfcor = FALSE,
                      r2pos=c("rss","ess","cor")){

  type <- match.arg(type)
  which <- match.arg(which)
  r2pos <- match.arg(r2pos)
  
  sel <- ifelse(!is.null(object$coef.names$sel), TRUE, FALSE)
  reg <- ifelse(!is.null(object$coef.names$reg), TRUE, FALSE)
  ifr <- ifelse(!is.null(object$coef.names$ifr), TRUE, FALSE)
  
  KX <- length(object$coef.names$reg)
  KS <- length(object$coef.names$sel)
  KP <- length(object$coef.names$ifr)
  if (is.null(KS)) KS <- 0
  if (is.null(KP)) KP <- 0

  y <- model.response(model.frame(object))
  n <- length(y)
  no <- sum(y == 0)
  po <- mean(y==0)
  pp <- 1 - po

  ## On attribute KS et KP et po * KX à l'explication des 0 et (1 - p)
  ## * KX à celle des valeurs positives.

  K <- switch(which,
              positive = (1 - po) * KX,
              zero     = KS + KP + po * KX,
              all      = KS + KP + KX
              )
  Ko <- switch(which,
               positive = (1 - po) * 2,
               zero     = (KS > 0) + (KP > 0) + po * 2,
               all      = (KS > 0) + (KP > 0) + 2
               )
  
  if (type == "mcfadden"){
    if (!dfcor) R2 <- 1 - logLik(object, which, naive = TRUE)/logLik(object, which)
    else R2 <- 1 - (logLik(object, which, naive = TRUE)- K)/(logLik(object, which) - Ko)
  }
  if (type == "regression"){
    if (which != "positive"){
      pfo <- fitted(object, "zero")
      eo <- 0 - pfo
      if (po != 0){
        if (!dfcor) R2n <- sum( (pfo - po)^2) / (no * po^2)
        else R2n <- sum( (pfo - po)^2) / (no * po^2)
      }
      else R2n <- 0
    }
    if (which != "zero"){
      yp <- y[y > 0]
      yf <- fitted(object, "positive")[y > 0]
      ym <- mean(yp)
      if (r2pos=="ess") R2p <- sum( (yf - ym)^2 ) / sum( (yp - ym)^2)
      if (r2pos=="rss") R2p <- 1 - sum( (yp - yf)^2 )/sum( (yp - ym)^2)
      if (r2pos=="cor") R2p <- cor(yf, yp)^2
    }
    if (which == "all") R2 <- po * R2n + (1 - po) * R2p
  }
  result <- switch(which,
                   all = R2,
                   zero = R2n,
                   positive = R2p
                   )
  result
}
