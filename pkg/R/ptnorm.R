punorm <- function(z){
    z <- as.double(z)
    lg <- as.integer(length(z))
    prob <- double(lg)
    z <- sanitize(z, string = c("z", "punorm"), replace = TRUE, verbal = FALSE)
    ans <- .Fortran("MYUNORM", prob, z, lg)[[1]]
    ans
}
    
pbnorm <- function(z1, z2, rho){
    z1 <- as.double(z1)
    z2 <- as.double(z2)
    al <- as.integer(length(z1))
    rho <- as.double(rho)
    prob <- double(al)
    z1 <- sanitize(z1, string = c("z1", "pbnorm"), replace = TRUE, verbal = FALSE)
    z2 <- sanitize(z2, string = c("z2", "pbnorm"), replace = TRUE, verbal = FALSE)
    ans <- .Fortran("MYBNORM", prob, z1, z2, rho, al)[[1]]
    ans
}

ptnorm <- function(z1, z2, z3, rho){
    # rho23 must be the largest corr coefficient in absolute value
    z <- cbind(z1, z2, z3)
    mxr <- which.max(abs(rho))
    if (mxr == 1){
        posh <- 1:3
        posr <- 1:3
    }
    if (mxr == 2){
        posh <- c(2, 1, 3)
        posr <- c(1, 3, 2)
    }
    if (mxr == 3){
        posh <- c(3, 1, 2)
        posr <- c(2, 3, 1)
    }
    z <- z[, posh, drop = FALSE]
    rho <- rho[posr]
    z1 <- as.double(z[, 1])
    z2 <- as.double(z[, 2])
    z3 <- as.double(z[, 3])
    z1 <- sanitize(z1, string = c("z1", "ptnorm"), verbal = FALSE, replace = TRUE)
    z2 <- sanitize(z2, string = c("z2", "ptnorm"), verbal = FALSE, replace = TRUE)
    z3 <- sanitize(z3, string = c("z3", "ptnorm"), verbal = FALSE, replace = TRUE)
    al <- as.integer(nrow(z))
    rho <- as.double(rho)
    prob <- double(al)
    ans <- .Fortran("MYTNORM", prob, z1, z2, z3, rho, al)[[1]]
    ans
    return(ans)
}

PHI2 <- function(z1, z2, rho){
    usebiv <- TRUE
    if (length(z1) == 1 && is.infinite(z1)){
        if (sign(z1) == 1){
            f <- punorm(z2)
            d1 <- 0 ; d2 <- dnorm(z2)
            dr <- 0
        }
        else{
            f <- 0
            d1 <- 0 ; d2 <- 0
            dr <- 0
        }
        usebiv <- FALSE
    }
    if (length(z2) == 1 && is.infinite(z2)){
        if (sign(z2) == 1){
            f <- pnorm(z1)
            d1 <- dnorm(z1) ; d2 <- 0
            dr <- 0
        }
        else{
            f <- 0
            d1 <- 0 ; d2 <- 0
            dr <- 0
        }
        usebiv <- FALSE
    }
    if (usebiv){
        if (rho == 0) f = pnorm(z1) * pnorm(z2) else f <- pbnorm(z1, z2, rho)
        d1 <- dnorm(z1) * pnorm( (z2 - rho * z1) / sqrt(1 - rho ^ 2) )
        d2 <- dnorm(z2) * pnorm( (z1 - rho * z2) / sqrt(1 - rho ^ 2) )
        dr <- dnorm(z1) * dnorm( (z2 - rho * z1) / sqrt(1 - rho ^ 2) )  / sqrt(1 - rho ^ 2)
    }
    list(f = f, d1 = d1, d2 = d2, dr = dr)
}

dbnorm <- function(z1, z2, rho)
    1 / (2 * pi * sqrt(1 - rho ^ 2)) *
        exp(- 0.5 * (z1 ^ 2 + z2 ^ 2 - 2 * rho * z1 * z2) / (1 - rho ^ 2))

PHI3 <- function(z1, z2, z3, rho){
    h1 <- ! (length(z1) == 1 && is.infinite(z1))
    h2 <- ! (length(z2) == 1 && is.infinite(z2))
    h3 <- ! (length(z3) == 1 && is.infinite(z3))

    s1 <- s2 <- s3 <- 0
    if (! h1) s1 <- sign(z1)
    if (! h2) s2 <- sign(z2)
    if (! h3) s3 <- sign(z3)
    
    # The complete case, all the arguments are finite
    if (all(h1, h2, h3)){
        if (all(rho == 0)) PT <- pnorm(z1) * pnorm(z2) * pnorm(z3)
        else PT <- ptnorm(z1, z2, z3, rho)
        d1 <- dnorm(z1) * pbnorm((z2 - rho[1] * z1) / sqrt( (1 - rho[1] ^ 2) ),
                                 (z3 - rho[2] * z1) / sqrt( (1 - rho[2] ^ 2) ),
                                 (rho[3] - rho[1] * rho[2]) / sqrt(1 - rho[1] ^ 2) /
                                     sqrt(1 - rho[2] ^ 2))
        d2 <- dnorm(z2) * pbnorm((z1 - rho[1] * z2) / sqrt( (1 - rho[1] ^ 2) ),
                                 (z3 - rho[3] * z2) / sqrt( (1 - rho[3] ^ 2) ),
                                 (rho[2] - rho[1] * rho[3]) / sqrt(1 - rho[1] ^ 2) /
                                     sqrt(1 - rho[3] ^ 2))
        d3 <- dnorm(z3) * pbnorm((z1 - rho[2] * z3) / sqrt( (1 - rho[2] ^ 2) ),
                                 (z2 - rho[3] * z3) / sqrt( (1 - rho[3] ^ 2) ),
                                 (rho[1] - rho[2] * rho[3]) / sqrt(1 - rho[2] ^ 2) /
                                     sqrt(1 - rho[3] ^ 2))
        S <- 1 - rho[1] ^ 2 - rho[2] ^ 2  - rho[3] ^ 2 + 2 * rho[1] * rho[2] * rho[3]
        if (S < 0) stop(paste("the coefficients of correlations are out of range",
                              paste(round(rho, 4), collapse = ","), sep = " "))
        dr1 <- dbnorm(z1, z2, rho[1]) *
            punorm( (
                (1 - rho[1] ^ 2) * z3 -
                    (rho[2] - rho[1] * rho[3]) * z1 -
                        (rho[3] - rho[1] * rho[2]) * z2) / sqrt( (1 - rho[1] ^ 2) * S))
        dr2 <- dbnorm(z1, z3, rho[2]) *
            punorm( (
                (1 - rho[2] ^ 2) * z2 -
                    (rho[1] - rho[2] * rho[3]) * z1 -
                        (rho[3] - rho[1] * rho[2]) * z3) / sqrt( (1 - rho[2] ^ 2) * S))
        dr3 <- dbnorm(z2, z3, rho[3]) *
            punorm( (
                (1 - rho[3] ^ 2) * z1 -
                    (rho[1] - rho[2] * rho[3]) * z2 -
                        (rho[2] - rho[1] * rho[3]) * z3) / sqrt( (1 - rho[3] ^ 2) * S))
        dr <- cbind(dr1, dr2, dr3)
    }
    else{
        # The trivial case, one of the bounds is - Infty, so that the result is 0
        if (any(c(s1, s2, s3) == -1)) PT <- d1 <- d2 <- d3 <- dr <- 0
        else{
            # The resulting cases, at least one of z1, z2, and z3 are + Inf
            if (h1){
                if (! h2){
                    if (h3){
                        # h1 and h3 -> bivariate
                        arho <- rho[2]
                        if (arho == 0) PT <- punorm(z1) * punorm(z3) else PT <- pbnorm(z1, z3, arho)
                        d1 <- dnorm(z1) * punorm( (z3 - arho * z1) / sqrt(1 - arho ^ 2) )
                        d3 <- dnorm(z3) * punorm( (z1 - arho * z3) / sqrt(1 - arho ^ 2) )
                        dr <- dnorm(z1) *  dnorm( (z3 - arho * z1) / sqrt(1 - arho ^ 2) ) / sqrt(1 - arho ^ 2)
                    }
                    else{
                        # h1 -> univariate
                        PT <- punorm(z1)
                        d1 <- dnorm(z1)
                        dr <- 0
                    }
                }
                else{
                    # h1 and h2 is present (so that h3 is not)
                    arho <- rho[1]
                    if (arho == 0) PT <- punorm(z1) * punorm(z2) else PT <- pbnorm(z1, z2, arho)
                    d1 <- dnorm(z1) * punorm( (z2 - arho * z1) / sqrt(1 - arho ^ 2) )
                    d2 <- dnorm(z2) * punorm( (z1 - arho * z2) / sqrt(1 - arho ^ 2) )
                    dr <- dnorm(z1) *  dnorm( (z2 - arho * z1) / sqrt(1 - arho ^ 2) ) / sqrt(1 - arho ^ 2)
                }
            }
            else{
                if (! h2){
                    # neither h1 and h2 -> univariate z3 distribution or 1
                    if (h3){
                        PT <- punorm(z3)
                        d3 <- dnorm(z3)
                        dr <- 0
                    }
                    else PT <- 1
                }
                else{
                    if (h3){
                        # h2 and h3 -> bivariate
                        arho <- rho[3]
                        if (arho == 0) PT <- punorm(z2) * punorm(z3) else PT <- pbnorm(z2, z3, arho)
                        d2 <- dnorm(z2) * punorm( (z3 - arho * z2) / sqrt(1 - arho ^ 2) )
                        d3 <- dnorm(z3) * punorm( (z2 - arho * z3) / sqrt(1 - arho ^ 2) )
                        dr <- dnorm(z2) *  dnorm( (z3 - arho * z2) / sqrt(1 - arho ^ 2) ) / sqrt(1 - arho ^ 2)
                    }                    
                    else{
                        # h2 > univariate
                        PT <- punorm(z2)
                        d2 <- dnorm(z2)
                        dr <- 0
                    }
                }
            }
        }
    }           
    if (! h1) d1 <- 0
    if (! h2) d2 <- 0
    if (! h3) d3 <- 0
    list(f = PT, d1 = d1, d2 = d2, d3 = d3, dr = dr)
}

