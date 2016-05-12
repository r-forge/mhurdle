mhurdle <- function(formula, data, subset, weights, na.action,
                    start = NULL, dist = c("ln", "tn", "n", "bc", "ihs", "ln2", "bc2"),
                    corr = FALSE, robust = TRUE, ...){
    fitted = FALSE
    check.grad = FALSE
    dots <- list(...)
    oldoptions <- options(warn = -1)
    on.exit(options(oldoptions))
    cl <- match.call()
    posT <- as.list(cl) == "T"
    posF <- as.list(cl) == "F"
    cl[posT] <- TRUE
    cl[posF] <- FALSE
    cl.save <- cl
    dist <- match.arg(dist)
    isMu <- dist == "bc2"
    if (robust){
        if (! (isMu | corr)){
            robust <- FALSE
            cat("robust irrelevant\n")
        }
    }

    # 1. Compute the model.frame and the model.matrix

    if (!inherits(formula, "Formula")) formula <- Formula(formula)
    if (length(formula)[2] > 4) stop("at most 4 rhs should be provided in the formula")
    mf <- match.call(expand.dots = FALSE) 
    m <- match(c("formula", "data", "subset", "na.action", "weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf$formula <- formula
    mf <- eval(mf, parent.frame())
    X1 <- model.matrix(formula, data = mf, rhs = 1)
    X2 <- model.matrix(formula, data = mf, rhs = 2)
    X3 <- model.matrix(formula, data = mf, rhs = 3)
    X4 <- model.matrix(formula, data = mf, rhs = 4)
    y <- model.response(mf)
    n <- length(y)
    if (length(X1) == 0) X1 <- NULL
    if (length(X2) == 0) stop("the second hurdle (consumption equation) is mandatory")
    if (length(X3) == 0) X3 <- NULL
    if (length(X4) == 0) X4 <- NULL
    h1 <- ! is.null(X1)
    h3 <- ! is.null(X3)
    
    #  2. One equation models

    if (!h1 && !h3 && dist != "n"){
        if (dist == "ln") result <- lm( log(y) ~ X2 - 1, subset = y > 0)
        if (dist == "tn") result <- truncreg(y ~ X2 - 1, scaled = TRUE, subset = y > 0)
        if (dist == "bc") result <- boxcoxreg(y ~ X2 - 1, subset = y > 0)
        if (dist != "ln") names(result$coefficients) <- c(colnames(X2), "sigma")
        else names(result$coefficients) <- colnames(X2)
        return(result)
    }

    # 3. Compute the naive model

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

    # 4. Selection single hurdle models without correlation that can
    # be estimated simply in two parts using fit.simple.mhurdle()

    if (h1 && !h3 && (dist %in% c("ln", "bc", "tn")) && ! corr){
        result <- fit.simple.mhurdle(X1, X2, y, dist = dist)
        result$naive <- naive
        result$call <- cl.save
        result$model <- mf
        result$formula <- formula
        return(result)
    }

    # 5. Compute the starting values if not provided (use the linear
    # specification as the starting value for ihs and the log-linear
    # specification for Box-Cox)
    if (is.null(start)){
        # for (100d), use the results of (100i) as starting values
        if (h1 && !h3 & (dist %in% c("ln", "bc", "tn"))){
            start <- coef(fit.simple.mhurdle(X1, X2, y, dist = dist))
            if (dist == "bc"){
                start.lambda <- start[length(start)]
                start <- start[- length(start)]
            }
        }
        else{
            dist.start <- dist
            if (dist %in% c("bc", "bc2", "ln2")) dist.start <- "ln"
            if (dist == "bc") start.lambda <- 0.01
            if (dist == "ihs") dist.start <- "n"
            start <- start.mhurdle(X1, X2, X3, y, dist.start)
        }
        # in case of heteroscedasctic model, add K4 zeros to the start
        # vector and the intercept should be ln(sigma_o) (not sigma_o)
        # because of the exp form
        
        sd.pos <- ifelse(h1, ncol(X1), 0) + ncol(X2) + ifelse(h3, ncol(X3), 0) + 1
        start[sd.pos] <- log(start[sd.pos])
        
        if (!is.null(X4)){
#            start <- c(start[1:sd.pos], rep(0, ncol(X4)))
            start <- c(start[1:sd.pos], rnorm(ncol(X4), sd = 0.001))
        }
        
        # add shape and/or scale parameters
        if (corr){
            if (robust) rhoinit <- tan(0.0 * pi / 2) else rhoinit <- 0.01
            if (h1 + h3 == 2){
                start <- c(start, rho12 = rhoinit, rho13 = rhoinit, rho23 = rhoinit)
            }
            else start <- c(start, rho = rhoinit)
        }
        if (dist == "bc") start <- c(start, tr = start.lambda)
        if (dist == "bc2") start <- c(start, tr = 0.01, pos = 0.1)
        if (dist == "ihs") start <- c(start, tr = 0.01)
        if (dist == "ln2") start <- c(start, pos = 1)
    }
    result <- mhurdle.fit(start, X1, X2, X3, X4, y,
                          gradient = TRUE, fit = FALSE,
                          dist = dist, corr = corr,
                          robust = robust, fitted = fitted,
                          check.grad = check.grad, ...)
    cat("\n\n")
    result$naive <- naive
    result$call <- cl.save
    result$formula <- formula
    result$model <- mf
    result
}

mhurdle.fit <- function(start, X1, X2, X3, X4, y, gradient = FALSE, fit = FALSE,
                        dist = c("ln", "n", "tn", "bc", "ihs", "bc2", "ln2"),
                        corr = FALSE, robust = TRUE,  fitted = FALSE,
                        check.grad = FALSE, ...){
    start.time <- proc.time()
    
    # fancy coefficients names
    if (corr){
        nbeqs <- 1 + (! is.null(X1)) + (! is.null(X3))
        KR <- ifelse(nbeqs == 2, 1, 3)
    }
    else KR <- 0

    sd.names <- "sd"
    
    if (corr){
        if (KR == 3) rho.names <- c("corr12", "corr13", "corr23")
        else{
            if (! is.null(X1)) rho.names <- c("corr12")
            else rho.names <- c("corr23")
        }
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
    if (! is.null(X1)) start.names$h1 <- paste("h1", start.names$h1, sep = ".")
    start.names$h2 <- paste("h2", start.names$h2, sep = ".")
    if (! is.null(X3)) start.names$h3 <- paste("h3", start.names$h3, sep = ".")
    if (! is.null(X4)) start.names$h4 <- paste("h4", start.names$h4, sep = ".")
    names(start) <- Reduce("c", start.names)

#    cat("Starting values:\n")
#    print(start)
    
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
        print(cbind(start, agrad, ngrad))
        print(as.numeric(sum(fo)))
    }

    maxl <- maxLik(f, start = start,...)
    nb.iter <- maxl$iterations
    convergence.OK <- maxl$code <= 2
    coefficients <- maxl$estimate
    if (corr){
        nbeqs <- 1 + (! is.null(X1)) + (! is.null(X3))
        KR <- ifelse(nbeqs == 2, 1, 3)
    }
    else KR <- 0
#    if (robust){
    if (FALSE){
        N4 <- ifelse(is.null(X4), 1, ncol(X4))
        if (! is.null(corr)){
            firstCorr <- sum(c(ncol(X1), ncol(X2), ncol(X3), N4)) + 1
            posCorr <- firstCorr:(firstCorr + KR - 1)
            coefficients[posCorr] <- atan(coefficients[posCorr]) * 2 / pi
        }
        
        if (dist %in% c("bc2", "ln2")){
            if (dist == "bc2") posMu <- sum(c(ncol(X1), ncol(X2), ncol(X3), N4)) + KR + 2
            else  posMu <- sum(c(ncol(X1), ncol(X2), ncol(X3), N4)) + KR + 1
            coefficients[posMu] <- exp(coefficients[posMu])
        }
        f <- function(param) mhurdle.lnl(param, X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y,
                                         gradient = TRUE, fitted = FALSE,
                                         dist = dist, corr = corr,
                                         robust = FALSE)
        maxl <- maxLik(f, start = coefficients, iterlim = 0, ...)
   }



    
    if (fitted){
        fitted.values <- attr(mhurdle.lnl(coefficients, X1 = X1, X2 = X2, X3 = X3, X4 = X4, y = y,
                                          gradient = FALSE, fitted = TRUE, robust = FALSE,
                                          dist = dist, corr = corr), "fitted")
    }
    else fitted.values <- NULL

    # La ligne ci-dessous renvoie la contribution de chaque obs a la
    # vraisemblance et au gradient (en attribut)
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
    ocoef <- coefficients
    gtheta <- rep(1, length(coefficients))
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
    if (ncol(X2) > 1) class(result) <- c("mhurdle")
    result
}

