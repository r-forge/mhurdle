#' Estimation of limited dependent variable models
#' 
#' mhurdle fits a large set of models relevant when the dependent variable is 0
#' for a part of the sample.
#' 
#' 
#' \code{mhurdle} fits models for which the dependent variable is zero for a
#' part of the sample. Null values of the dependent variable may occurs because
#' of one or several mechanisms : good rejection, lack of ressources and
#' purchase infrequency. The model is described using a three-parts formula :
#' the first part describes the selection process if any, the second part the
#' regression equation and the third part the purchase infrequency process.
#' \code{y ~ 0 | x1 + x2 | z1 + z2} means that there is no selection process.
#' \code{y ~ w1 + w2 | x1 + x2 | 0} and \code{y ~ w1 + w2 | x1 + x2} describe
#' the same model with no purchase infrequency process. The second part is
#' mandatory, it explains the positive values of the dependant variable. The
#' \code{dist} argument indicates the distribution of the error term. If
#' \code{dist = "n"}, the error term is normal and (at least part of) the zero
#' observations are also explained by the second part as the result of a corner
#' solution. Several models described in the litterature are obtained as
#' special cases :
#' 
#' A model with a formula like \code{y~0|x1+x2} and \code{dist="n"} is the
#' Tobit model proposed by Tobin (1958).
#' 
#' \code{y~w1+w2|x1+x2} and \code{dist="l"} or \code{dist="t"} is the single
#' hurdle model proposed by Cragg (1971). With \code{dist="n"}, the double
#' hurdle model also proposed by Cragg (1971) is obtained. With
#' \code{corr="h1"} we get the correlated version of this model described by
#' Blundell (1987).
#' 
#' \code{y~0|x1+x2|z1+z2} is the P-Tobit model of Deaton and Irish (1984),
#' which can be a single hurdle model if \code{dist="t"} or \code{dist="l"} or
#' a double hurdle model if \code{dist="n"}.
#' 
#' @aliases mhurdle coef.mhurdle vcov.mhurdle logLik.mhurdle print.mhurdle
#' summary.mhurdle print.summary.mhurdle predict.mhurdle update.mhurdle
#' fitted.mhurdle effects.mhurdle
#' @param formula a symbolic description of the model to be fitted,
#' @param data a \code{data.frame},
#' @param newdata a \code{data.frame} for which the predictions should be
#' computed,
#' @param subset see \code{\link{lm}},
#' @param weights see \code{\link{lm}},
#' @param na.action see \code{\link{lm}},
#' @param start starting values,
#' @param dist the distribution of the error of the consumption equation: one
#' of \code{"n"} (normal), \code{"ln"} (log-normal) \code{"bc"} (box-cox
#' normal) and \code{"ihs"} (inverse hyperbolic sinus transformation),
#' @param h2 if \code{TRUE} the second hurdle is effective, it is not
#' otherwise,
#' @param scaled if \code{TRUE}, the dependent variable is divided by its
#' geometric mean,
#' @param corr a boolean indicating whether the errors of the different
#' equations are correlated or not,
#' @param robust transformation of the structural parameters in order to avoid
#' numerical problems,
#' @param check.grad if \code{TRUE}, a matrix containing the analytical and the
#' numerical gradient for the starting values are returned,
#' @param naive a boolean, it \code{TRUE}, the likelihood of the naive model is
#' returned,
#' @param object,x an object of class \code{"mhurdle"},
#' @param new an updated formula for the \code{update} method,
#' @param digits see \code{\link{print}},
#' @param width see \code{\link{print}},
#' @param which which coefficients or covariances should be extracted ? Those
#' of the selection (\code{"h1"}), consumption (\code{"h2"}) or purchase
#' (\code{"h3"}) equation, the other coefficients \code{"other"} (the standard
#' error and the coefficient of corr), the standard error (\code{"sigma"}) or
#' the coefficient of correlation (\code{"rho"}),
#' @param covariate the covariate for which the effect has to be computed,
#' @param reflevel for the computation of effects for a factor, the reference
#' level,
#' @param mean if \code{TRUE}, the mean of the effects is returned,
#' @param \dots further arguments.
#' @return
#' 
#' an object of class \code{c("mhurdle", "maxLik")}.
#' 
#' A \code{"mhurdle"} object has the following elements :
#' 
#' \describe{ \item{coefficients}{the vector of coefficients,} \item{vcov}{the
#' covariance matrix of the coefficients,} \item{fitted.values}{a matrix of
#' fitted.values, the first column being the probability of 0 and the second
#' one the mean values for the positive observations,} \item{logLik}{the
#' log-likelihood,} \item{gradient}{the gradient at convergence,}
#' \item{model}{a data.frame containing the variables used for the estimation,}
#' \item{coef.names}{a list containing the names of the coefficients in the
#' selection equation, the regression equation, the infrequency of purchase
#' equation and the other coefficients (the standard deviation of the error
#' term and the coefficient of correlation if \code{corr = TRUE}),}
#' \item{formula}{the model formula, an object of class \code{Formula},}
#' \item{call}{the call,} \item{rho}{the lagrange multiplier test of no
#' correlation.}
#' 
#' }
#' @references
#' 
#' Blundell R, Meghir C (1987). Bivariate Alternatives to the Tobit Model.
#' Journal of Econometrics, 34, 179-200.
#' 
#' Cragg JG (1971). Some Statistical Models for Limited Dependent Variables
#' with Applications for the Demand for Durable Goods. Econometrica, 39(5),
#' 829-44.
#' 
#' Deaton A, Irish M (1984). A Statistical Model for Zero Expenditures in
#' Household Budgets.  Journal of Public Economics, 23, 59-80.
#' 
#' Tobin J (1958). Estimation of Relationships for Limited Dependent Variables.
#' Econometrica, 26(1), 24-36.
#' @keywords regression
#' @examples
#' 
#' 
#' data("Interview", package = "mhurdle")
#' 
#' # independent double hurdle model
#' idhm <- mhurdle(vacations ~ car + size | linc + linc2 | 0, Interview,
#'               dist = "ln", h2 = TRUE, method = "bfgs")
#' 
#' # dependent double hurdle model
#' ddhm <- mhurdle(vacations ~ car + size | linc + linc2  | 0, Interview,
#'               dist = "ln", h2 = TRUE, method = "bfgs", corr = TRUE)
#' 
#' # a double hurdle p-tobit model
#' ptm <- mhurdle(vacations ~ 0 | linc + linc2 | car + size, Interview,
#'               dist = "ln", h2 = TRUE, method = "bfgs", corr = TRUE)
#' 
#' 
mhurdle <- function(formula, data, subset, weights, na.action,
                    start = NULL, dist = c("ln", "n", "bc", "ihs"), h2 = FALSE,
                    scaled = TRUE, corr = FALSE, robust = TRUE,
                    check.grad = FALSE, ...){
    fitted = TRUE
    dots <- list(...)
    oldoptions <- options(warn = -1)
    on.exit(options(oldoptions))
    cl <- match.call()
    posT <- as.list(cl) == "T" ; posF <- as.list(cl) == "F"
    cl[posT] <- TRUE           ; cl[posF] <- FALSE
    dist <- match.arg(dist)

    if (dist == "ln" & h2) dist <- "ln2"
    if (dist == "n" & ! h2) dist <- "tn"
    if (dist == "bc" & h2) dist <- "bc2"
    
    # 1. Compute the model.frame and the model.matrix

    if (! inherits(formula, "Formula")) formula <- Formula(formula)
    if (length(formula)[2] > 4) stop("at most 4 rhs should be provided in the formula")
    mf <- match.call(expand.dots = FALSE)
    mc <- match.call(expand.dots = TRUE)
    
    m <- match(c("formula", "data", "subset", "na.action", "weights"),
               names(mc), 0L)
    mf <- mc[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf$formula <- formula
    mf <- eval(mf, parent.frame())
    X1 <- model.matrix(formula, data = mf, rhs = 1)
    X2 <- model.matrix(formula, data = mf, rhs = 2)
    X3 <- model.matrix(formula, data = mf, rhs = 3)
    X4 <- model.matrix(formula, data = mf, rhs = 4)
    # Remove the intercept if any for the heteroscedastic model
    if (ncol(X4)){
        int4 <- attr(terms(formula, data = mf, rhs = 4), "intercept")
        X4 <- X4[, - 1, drop = FALSE]
    }
    y <- model.response(mf)
    if (scaled){
        geomean <- exp(mean(log(y[y > 0])))
        y <- y / geomean
        attr(y, "geomean") <- attr(mf, "geomean") <- geomean
    }
    N <- length(y)
    if (length(X1) == 0) X1 <- NULL
    if (length(X2) == 0) stop("the second hurdle (consumption equation) is mandatory")
    if (length(X3) == 0) X3 <- NULL
    if (length(X4) == 0) X4 <- NULL
    h1 <- ! is.null(X1)
    h3 <- ! is.null(X3)
    
    #  2. One equation models
    if (! h1 && ! h3){
        result <- onequation.mhurdle(X2, y, dist)
        result$naive <- NULL#naive
        result$call <- cl
        result$model <- mf
        result$formula <- formula
        result$R2 <- c(null = NA, positive = NA)
        return(result)
    }
    
    # 3. Selection single hurdle models without correlation that can
    # be estimated simply in two parts using seperate.mhurdle()
    if (h1 && !h3  && ! corr && !h2){
        result <- seperate.mhurdle(X1, X2, y, dist = dist)
        result$naive <- NULL#naive
        result$call <- cl
        result$model <- mf
        result$formula <- formula
        result$R2 <- c(null = NA, positive = NA)
        return(result)
    }

    # 5. Compute the starting values if not provided (use the linear
    # specification as the starting value for ihs and the log-linear
    # specification for Box-Cox)

    if (is.null(start)){
        dist.start <- dist
        if (dist %in% c("bc", "bc2", "ln2")) dist.start <- "ln"
        if (dist == "ihs") dist.start <- "n"
        start <- start.mhurdle(X1, X2, X3, y, dist.start)
        
        # in case of heteroscedasctic model, add K4 zeros to the start
        # vector and the intercept should be ln(sigma_o) (not sigma_o)
        # because of the exp form
        sd.pos <- getindex(X1, X2, X3, X4, corr, dist, which = "sd")
        if (robust) start[sd.pos] <- log(start[sd.pos])
        if (! is.null(X4)) start <- c(start[1:sd.pos], rep(0, ncol(X4)))
        
        # add shape and/or scale parameters
        if (corr){
            if (robust) rhoinit <- tan(0.1 * pi / 2) else rhoinit <- 0.1
            if (h1 + h3 == 2) start <- c(start, rho12 = rhoinit, rho13 = rhoinit, rho23 = rhoinit)
            else start <- c(start, rhoinit)
        }
        if (dist %in% c("bc", "bc2", "ihs")) start <- c(start, tr = -0.1)
        if (dist %in% c("bc2", "ln2")) start <- c(start, pos = 1)
    }
    else{
        if (robust){
            sd.pos <- getindex(X1, X2, X3, X4, corr, dist, which = "sd")
            mu.pos <- getindex(X1, X2, X3, X4, corr, dist, which = "pos")
            rho.pos <- getindex(X1, X2, X3, X4, corr, dist, which = "corr")
            start[c(sd.pos, mu.pos)] <- log(start[c(sd.pos, mu.pos)])
            start[rho.pos] <- tan(start[rho.pos] * pi / 2)
        }
    }
    result <- mhurdle.fit(start, X1, X2, X3, X4, y, gradient = TRUE,
                          dist = dist, corr = corr, robust = robust,
                          fitted = fitted, check.grad = check.grad,
                          ...)
    if (check.grad) return(result)
    
    # 3. Compute the naive model

    Pnull <- mean(y == 0)
    dist.naive <- dist
    if (dist %in% c("ihs")) dist.naive <- "n"
    if (dist %in% c("bc", "bc2", "ln2")) dist.naive <- "ln"
    if (dist.naive != "ln"){
        Ec <- mean(y[y > 0])
        Vc <- var(y[y > 0])}
    else{
        Ec <- mean(log(y[y > 0]))
        Vc <- var(log(y[y > 0]))
    }

    start.naive <- c(rep(0.1, 1 + h1 + h3), 1)
    moments <- c(Pnull, Ec, Vc)
    naive <- maxLik(lnl.naive, start = start.naive,
                    dist = dist.naive, moments = moments,
                    h1 = h1, h3 = h3)
    coef.naive <- naive$est
    parts <- attr(naive$max, "parts") * c(Pnull, 1 - Pnull, 1 - Pnull) * N
    logLik.naive <- structure(as.numeric(naive$max) * N, nobs = length(y),
                              parts = parts, 
                              df = length(coef.naive), class = "logLik")
    naive <- list(coefficients = coef.naive, logLik = logLik.naive, code = naive$code)

    result$naive <- naive
    result$call <- cl
    result$formula <- formula
    result$model <- mf
    
    lnL1c <- attr( naive$logLik, "parts")[1] + attr( naive$logLik, "parts")[2]
    lnL1u <- attr(result$logLik, "parts")[1] + attr(result$logLik, "parts")[2]

    lnL2c <- attr( naive$logLik, "parts")[3] - attr( naive$logLik, "parts")[2]
    lnL2u <- attr(result$logLik, "parts")[3] - attr(result$logLik, "parts")[2]

    R1 <- 1 - (lnL1c / lnL1u) ^ (2 / N * lnL1c)
    R2 <- 1 - (exp(lnL2c) / exp(lnL2u)) ^ (2 / (N * (1 - Pnull)))
    
    result$R2 <- c(null = R1, positive = R2)
    result
}

mhurdle.fit <- function(start, X1, X2, X3, X4, y, gradient = FALSE, fit = FALSE,
                        dist = c("ln", "n", "tn", "bc", "ihs", "bc2", "ln2"),
                        corr = FALSE, robust = TRUE,  fitted = FALSE,
                        check.grad = FALSE, ...){
    start.time <- proc.time()
    h1 <- ! is.null(X1)
    h3 <- ! is.null(X3)
    h4 <- ! is.null(X4)
    
    # fancy coefficients names
    sd.names <- "sd"
    
    if (corr){
        if (h1 & h3) rho.names <- c("corr12", "corr13", "corr23")
        else rho.names <- ifelse(h1, "corr12", "corr23")
    }
    else rho.names <- NULL
    if (dist %in% c("bc", "bc2", "ihs")) tr.names <- "tr" else tr.names <- NULL
    if (dist %in% c("ln2", "bc2")) mu.names <- "pos" else mu.names <- NULL

    coef.names <- list(h1   = colnames(X1),
                       h2   = colnames(X2),
                       h3   = colnames(X3),
                       sd   = sd.names,
                       h4   = colnames(X4),
                       corr = rho.names,
                       tr   = tr.names,
                       pos   = mu.names)

    start.names <- coef.names
    if (h1) start.names$h1 <- paste("h1", start.names$h1, sep = ".")
    start.names$h2 <- paste("h2", start.names$h2, sep = ".")
    if (h3) start.names$h3 <- paste("h3", start.names$h3, sep = ".")
    if (h4) start.names$h4 <- paste("h4", start.names$h4, sep = ".")
    names(start) <- Reduce("c", start.names)

    f <- function(param) mhurdle.lnl(param, X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y,
                                     gradient = TRUE, fitted = FALSE,
                                     dist = dist, corr = corr,
                                     robust = robust)
    if (check.grad){
        ngrad <- c()
        oparam <- start
        fo <- f(start)
        agrad <- apply(attr(fo, "gradient"), 2, sum)
        eps <- 1E-05
        for (i in 1:length(start)){
            oparam[i] <- oparam[i] + eps
            ngrad <- c(ngrad, sum((as.numeric(f(oparam)) - fo) / eps))
            oparam <- start
        }
        return(cbind(start, agrad, ngrad))
    }
    maxl <- maxLik(f, start = start, control = list(lambdatol = 1E-20),  ...)
    nb.iter <- maxl$iterations
    convergence.OK <- maxl$code <= 2
    coefficients <- maxl$estimate

    if (fitted) fitted.values <- attr(mhurdle.lnl(coefficients, X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y,
                                                  gradient = FALSE, fitted = TRUE, robust = robust,
                                                  dist = dist, corr = corr), "fitted")
    else fitted.values <- NULL

    # contribution of every single observation to the likelihood and
    # its gradient (as an attribute)
    logLik <- f(coefficients)
    gradi <- attr(logLik, "gradi")
    logLik <- structure(as.numeric(logLik), df = length(coefficients),
                        parts = attr(logLik, "parts"),
                        nobs = length(y), class = "logLik")
    hessian <- maxl$hessian
    elaps.time <- proc.time() - start.time
    eps <- with(maxl, gradient %*% solve(- hessian) %*% gradient)
    est.stat <- list(elaps.time = elaps.time,
                     nb.iter = nb.iter,
                     eps = eps,
                     method = maxl$type,
                     message = maxl$message
                     )
    class(est.stat) <- "est.stat"
    gtheta <- rep(1, length(coefficients))
    if (robust){
        if (corr){
            poscor <- sub.mhurdle(coef.names, "corr")
            gtheta[poscor] <- 2 / pi / (1 + coefficients[poscor] ^ 2)
            coefficients[poscor] <- atan(coefficients[poscor]) * 2 / pi
        }
        if (dist %in% c("bc2", "ln2")){
            posmu <- sub.mhurdle(coef.names, "pos")
            gtheta[posmu] <- exp(coefficients[posmu])
            coefficients[posmu] <- exp(coefficients[posmu])
        }
        possd <- sub.mhurdle(coef.names, "sd")
        gtheta[possd] <- exp(coefficients[possd])
        coefficients[possd] <- exp(coefficients[possd])
    }
    result <- list(coefficients  = coefficients,
                   vcov          = diag(gtheta) %*% (- solve(maxl$hessian) ) %*% diag(gtheta),
                   fitted.values = fitted.values,
                   logLik        = logLik,
                   gradient      = gradi,
                   hessian       = hessian,
                   formula       = NULL,
                   model         = NULL,
                   coef.names    = coef.names,
                   call          = NULL,
                   est.stat      = est.stat,
                   naive         = NULL
                   )
    class(result) <- c("mhurdle", class(result))
    result
    
}


sanitize <- function(x, myeps = 1E-07, mymax = 1E02, string = c("", ""), replace = TRUE, verbal = TRUE){
    string <- paste("of", string[1], "in", string[2])
    if (replace){
        if (any(is.na(x))){
            if (verbal) cat(paste(sum(is.na(x)), "NA values", string, "replaced by 0\n"))
            x[is.na(x)] <- 0
        }
        if ( any(x > 0 & x < myeps)){
            if (verbal) cat(paste(sum(x > 0 & x < myeps), "values", string, "lower than",  myeps,"replaced by", myeps, "\n"))
            x[x > 0 & x < myeps] <- myeps
        }
        if (any(x < - mymax)){
            if (verbal) cat(paste(sum(x < - mymax), "values", string, "lower than", - mymax, "replaced by", - mymax, "\n"))
            x[x < - mymax] <- - mymax
        }
        if (any(x > mymax)){
            if (verbal) cat(paste(sum(x >  mymax), "values", string, "greater than",  mymax, "replaced by", mymax, "\n"))
            x[x > mymax] <- mymax
        }
    }
    x
}
    
