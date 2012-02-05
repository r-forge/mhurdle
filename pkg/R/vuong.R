

############################################################
## les fonctions donnant les quantiles et les proba pour la
## loi des Khi Deux pondérés
############################################################

pwchisq <- function(q, weights, lower.tail = TRUE, n = 1000){
  K <- length(weights)
  e <- matrix(rnorm(n * K) ^ 2, n, K)
  wcs <- apply(e, 1, function(x) sum(x * weights))
  F <- ecdf(wcs)
  ifelse(lower.tail, F(q), 1 - F(q))
}

qwchisq <- function(p, weights, lower.tail = TRUE, n = 1000){
  K <- length(weights)
  e <- matrix(rnorm(n * K) ^ 2, n, K)
  wcs <- apply(e, 1, function(x) sum(x * weights))
  ifelse(lower.tail, quantile(wcs, p), quantile(wcs, 1 - p))
}  

##########################################################
## Fonction Vuong 2 qui permet l'ensemble des tests 
## pour modèles emboîtés, non emboîtés et partiellement
## emboîtés, avec ou sans hypothèse vraie a priori
##########################################################

vuongtest <- function(x, y,
                      type = c("non.nested", "nested", "overlapping"),
                      hyp = FALSE){
  type <- match.arg(type)
  data.name <- c(
                 paste(deparse(substitute(x))),
                 paste(deparse(substitute(y)))
                 )
  set.seed(100)

  ## for convenience, call x the larger model
  ##   if (length(coef(x)) < length(coef(y))){
  ##     oldy <- y
  ##     y <- x
  ##     x <- oldy
  ##   }  
  lx <- as.numeric(x$logLik)
  ly <- as.numeric(y$logLik)
  logLx <- sum(lx)
  logLy <- sum(ly)
  Kx <- length(coef(x))
  Ky <- length(coef(y))
  n <- length(lx)
  
  if (length(ly) != length(lx))
    stop("the number of observations of the two models differ")
  
  LR <- logLx - logLy
  w2 <- 1 / n * sum((lx - ly) ^ 2) - (1 / n * LR) ^ 2
  
  if (type == "non.nested"){
    statistic <- c(z = LR / sqrt(n * w2))
    alternative <- ifelse(statistic > 0,
                          "The first model is better",
                          "The second model is better")
    method <- "Vuong Test (non-nested)"
    pval <- pnorm(abs(statistic), lower.tail = FALSE)
    parameter <- NULL
  }
  else{
    gradx <- attr(x$logLik, "gradi")
    grady <- attr(y$logLik, "gradi")
    BF <- t(gradx) %*% gradx / n
    BG <- t(grady) %*% grady / n
    AF <- -solve(vcov(x)) / n
    AG <- -solve(vcov(y)) / n
  }
  if (type == "nested"){
    if (Kx == Ky) stop("the two models are not nested")
    nestedOK <- prod(names(coef(y)) %in% names(coef(x)))
    if (! nestedOK || (Kx == Ky)) stop("the two models are not nested")
    parameter <- c(df = Kx - Ky)
    statistic <- 2 * LR
    method <- "Vuong Test (nested)"
    alternative <- "The larger model is better"
    if (! hyp){
      names(statistic) <- "wchisq"
      common.coef <- names(coef(x)) %in% names(coef(y))
      trans <- diag(1, Kx)[common.coef, ]
      W <- BF %*% (t(trans) %*% solve(AG) %*% trans - solve(AF))
      Z <- eigen(W)$values
      pval <- pwchisq(statistic, Z, lower.tail = FALSE)
    }
    else{
      names(statistic) <- "chisq"
      pval <- pchisq(statistic, parameter, lower.tail = FALSE)
    }
  }
  if (type == "overlapping"){
    BFG <- t(gradx) %*% grady / n
    BGF <- t(grady) %*% gradx / n
    W <- rbind(cbind(- BF %*% solve(AF), - BFG %*% solve(AG)),
               cbind(- BGF %*% solve(AF), - BG %*% solve(AG))
               )
    method <- "Vuong Test (overlapping)"
    if (! hyp){
      statistic <- c(wchisq = n * w2)
      Z <- eigen(W)$values ^ 2
      pval <- pwchisq(statistic, Z, lower.tail = FALSE)
      alternative <- c("two different models : use non-nested version of the Vuong test")
      parameter <- c(df = length(Z))
    }
    else{
      statistic <- c(wchisq = 2 * LR)
      Z <- eigen(W)$values
      q1 <- qwchisq(0.95, Z)
      q2 <- qwchisq(0.05, Z)
      pval <- pwchisq(statistic, Z)
      alternative <- ""
      parameter <- c(df = length(Z))
    }
  }
  result <- list(statistic = statistic,
                 method = method,
                 p.value = pval,
                 data.name = data.name,
#                 alternative = alternative,
                 parameter = parameter)
  class(result) <- "htest"
  result
}

