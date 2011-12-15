r.squared <- function(object,
                      type = c("regression", "mcfadden"),
                      dfcor = FALSE,
                      r2pos=c("ess","rss","cor")){

  type <- match.arg(type)
  r2pos <- match.arg(r2pos)
  
  K1 <- length(object$coef.names$h1)
  K2 <- length(object$coef.names$h2)
  K3 <- length(object$coef.names$h3)
  
  y <- model.response(model.frame(object))
  n <- length(y)
  no <- sum(y == 0)
  po <- mean(y==0)
  pp <- 1 - po

  K <- length(coef(object))
  Ko <- length(object$naive$coefficients)
  
  if (type == "mcfadden"){
    if (!dfcor) R2 <- 1 - logLik(object) / logLik(object, naive = TRUE)
    else R2 <- 1 - (logLik(object)- K) / (logLik(object, naive = TRUE) - Ko)
  }
  if (type == "regression"){
    pfo <- fitted(object, "zero")
    eo <- 0 - pfo
    if (po != 0){
      if (!dfcor) R2n <- sum( (pfo - po)^2) / (no * po^2)
      else R2n <- sum( (pfo - po)^2) / (no * po^2)
    }
    else R2n <- 0
    yp <- y[y > 0]
    yf <- fitted(object, "positive")[y > 0]
    ym <- mean(yp)
    if (r2pos == "ess") R2p <- sum( (yf - ym) ^ 2 ) / sum( (yp - ym) ^ 2)
    if (r2pos == "rss") R2p <- 1 - sum( (yp - yf) ^ 2 )/sum( (yp - ym) ^ 2)
    if (r2pos == "cor") R2p <- cor(yf, yp) ^ 2
    R2 <- po * R2n + (1 - po) * R2p
  }
  R2
}
