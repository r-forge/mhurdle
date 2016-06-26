mhurdle <- function(formula, data, subset, weights, na.action,
                    start = NULL, dist = c("ln", "n", "bc", "ihs"), h2 = FALSE,
                    scaled = TRUE, corr = FALSE, robust = TRUE, ...){
    fitted = TRUE
    check.grad <- FALSE
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
    }
    n <- length(y)
    if (length(X1) == 0) X1 <- NULL
    if (length(X2) == 0) stop("the second hurdle (consumption equation) is mandatory")
    if (length(X3) == 0) X3 <- NULL
    if (length(X4) == 0) X4 <- NULL
    h1 <- ! is.null(X1)
    h3 <- ! is.null(X3)
    
    #  2. One equation models
    if (!h1 && !h3){
        result <- onequation.mhurdle(X2, y, dist)
        result$naive <- NULL#naive
        result$call <- cl
        result$model <- mf
        result$formula <- formula
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
            else start <- c(start, rho = rhoinit)
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
    result <- mhurdle.fit(start, X1, X2, X3, X4, y,
                          gradient = TRUE, dist = dist, corr = corr,
                          robust = robust, fitted = fitted, ...)
    if (fitted & scaled) result$fitted.values[, 2] <- result$fitted.values[, 2] * geomean
    
    # 3. Compute the naive model

#    if (FALSE){
    Pnull <- mean(y == 0)
    if (dist != "ln"){
        Ec <- mean(y[y > 0])
        Vc <- var(y[y > 0])}
    else{
        Ec <- mean(log(y[y > 0]))
        Vc <- var(log(y[y > 0]))
    }
    start.naive <- c(rep(0.1, 1 + h1 + h3), 1)
    moments <- c(Pnull, Ec, Vc)
    dist.naive <- dist
    if (dist %in% c("ihs")) dist.naive <- "n"
    if (dist %in% c("bc", "bc2", "ln2")) dist.naive <- "ln"
    naive <- maxLik(lnl.naive, start = start.naive,
                    dist = dist.naive, moments = moments,
                    h1 = h1, h3 = h3)
    coef.naive <- naive$est
    logLik.naive <- structure(naive$max * n, nobs = length(y),
                              df = length(coef.naive), class = "logLik")
    naive <- list(coefficients = coef.naive, logLik = logLik.naive, code = naive$code)
#}
#    else naive <- NULL

    result$naive <- naive
    result$call <- cl
    result$formula <- formula
    result$model <- mf
    result
}

mhurdle.fit <- function(start, X1, X2, X3, X4, y, gradient = FALSE, fit = FALSE,
                        dist = c("ln", "n", "tn", "bc", "ihs", "bc2", "ln2"),
                        corr = FALSE, robust = TRUE,  fitted = FALSE, ...){
    start.time <- proc.time()
    h1 <- ! is.null(X1)
    h3 <- ! is.null(X3)
    h4 <- ! is.null(X4)
    
    # fancy coefficients names
    sd.names <- "sd"
    
    if (corr){
        if (h1 & h3) rho.names <- c("corr12", "corr13", "corr23")
        else rho.names <- c("corr")
    }
    else rho.names <- NULL
    if (dist %in% c("bc", "bc2", "ihs")) tr.names <- "tr" else tr.names <- NULL
    if (dist %in% c("ln2", "bc2")) mu.names <- "mu" else mu.names <- NULL

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
    check.grad <- FALSE
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
        print(cbind(start, agrad, ngrad))
        print(as.numeric(sum(fo)))
        stop()
    }
    maxl <- maxLik(f, start = start, ...)
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
    
