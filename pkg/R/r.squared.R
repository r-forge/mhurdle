rsq <- function(object,
                type = c("coefdet", "lratio"),
                adj = FALSE,
                r2pos=c("rss","ess","cor")){
  
  type <- match.arg(type)
  r2pos <- match.arg(r2pos)
  
  K1 <- length(object$coef.names$h1)
  K2 <- length(object$coef.names$h2)
  K3 <- length(object$coef.names$h3)
  
  y <- model.response(model.frame(object))
  n <- length(y)
  no <- sum(y == 0)
  po <- mean(y == 0)
  pp <- 1 - po

  K <- length(coef(object))
  Ko <- length(object$naive$coefficients)

  if (type == "lratio"){
    if (!adj) R2 <- 1 - logLik(object) / logLik(object, naive = TRUE)
    else R2 <- 1 - (logLik(object) - K) / (logLik(object, naive = TRUE) - Ko)
  }
  if (type == "coefdet"){
    ym <- mean(y)
    yf <- fitted(object, "positive") * (1 - fitted(object, "zero"))
    R2 <- switch(r2pos,
                 ess = ifelse(adj,
                   sum( (yf - ym) ^ 2) / sum( (y - ym) ^ 2) * (n - K) / (n - Ko),
                   sum( (yf - ym) ^ 2) / sum( (y - ym) ^ 2)),
                 rss = ifelse(adj,
                   1 - (n - Ko) / (n - K) * sum( (y - yf) ^ 2) / sum( (y - ym) ^ 2),
                   1 - sum( (y - yf) ^ 2) / sum( (y - ym) ^ 2)),
                 cor = ifelse(adj,
                   stop("no adjusted R2 using the correlation formula"),
                   cor(y, yf) ^ 2
                   )
                 )
  }
  R2
}
