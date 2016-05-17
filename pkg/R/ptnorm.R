## system("R CMD SHLIB ~/Dropbox/Forge/mhurdle/pkg/src/ptnorm.f")
## dyn.load("~/Dropbox/Forge/mhurdle/pkg/src/ptnorm.so")

REPLACE = TRUE
verbal <- FALSE

punorm <- function(z){
    z <- as.double(z)
    lg <- as.integer(length(z))
    prob <- double(lg)
    if (REPLACE){
        z[is.na(z)] <- 0     ; if (verbal) cat("some NA values of z in punorm replaced by 0\n")
        z[z > 100] <- 100    ; if (verbal) cat("some values of z > 100 in punorm replaced by 100\n")
        z[z < - 100] <- -100 ; if (verbal) cat("some values of z < -100 in punorm replaced by -100\n")
    }
    ans <- .Fortran("MYUNORM", prob, z, lg)[[1]]
    return(ans)
}
    
pbnorm <- function(z1, z2, rho){
    # ROB
    z1 <- pmax(-10, z1)
    z2 <- pmax(-10, z2)
    z1 <- as.double(z1)
    z2 <- as.double(z2)
    al <- as.integer(length(z1))
    rho <- as.double(rho)
    prob <- double(al)
    if (REPLACE){
        z1[is.na(z1)] <- 0 ; if (verbal) cat("some NA values of z1 in pbnorm replaced by 0\n")
        z2[is.na(z2)] <- 0 ; if (verbal) cat("some NA values of z2 in pbnorm replaced by 0\n")
    }
    ans <- .Fortran("MYBNORM", prob, z1, z2, rho, al)[[1]]
    ans
    return(ans)
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
    z <- z[, posh]
    rho <- rho[posr]
    z1 <- as.double(z[, 1])
    z2 <- as.double(z[, 2])
    z3 <- as.double(z[, 3])
    if (REPLACE){
        if (any(is.na(z1))) z1[is.na(z1)] <- 0               ; if (verbal) cat("some NA values of z1 in pbnorm replaced by 0\n")
        if (any(is.na(z2))) z1[is.na(z2)] <- 0               ; if (verbal) cat("some NA values of z2 in pbnorm replaced by 0\n")
        if (any(is.na(z3))) z1[is.na(z3)] <- 0               ; if (verbal) cat("some NA values of z3 in pbnorm replaced by 0\n")
        if (any(is.infinite(z1))) z1[is.infinite(z1)] <- 0   ; if (verbal) cat("some Inf values of z1 in pbnorm replaced by 0\n")
        if (any(is.infinite(z2))) z2[is.infinite(z2)] <- 0   ; if (verbal) cat("some Inf values of z2 in pbnorm replaced by 0\n")
        if (any(is.infinite(z3))) z3[is.infinite(z3)] <- 0   ; if (verbal) cat("some Inf values of z3 in pbnorm replaced by 0\n")     
    }
    al <- as.integer(nrow(z))
    rho <- as.double(rho)
    prob <- double(al)
    ans <- .Fortran("MYTNORM", prob, z1, z2, z3, rho, al)[[1]]
    ans
    return(ans)
}

PHI2 <- function(z1, z2, rho){
    h1 <- ! is.null(z1)
    h2 <- ! is.null(z2)
    if (h1 & h2){
        if (rho == 0) f = pnorm(z1) * pnorm(z2) else f <- pbnorm(z1, z2, rho)
        d1 <- dnorm(z1) * pnorm( (z2 - rho * z1) / sqrt(1 - rho ^ 2) )
        d2 <- dnorm(z2) * pnorm( (z1 - rho * z2) / sqrt(1 - rho ^ 2) )
        dr <- dnorm(z1) * dnorm( (z2 - rho * z1) / sqrt(1 - rho ^ 2) )  / sqrt(1 - rho ^ 2)
    }
    if (! h1 & h2){
        f <- punorm(z2)
        d1 <- NULL
        d2 <- dnorm(z2)
        dr <- NULL
    }
    if (h1 & ! h2){
        f <- punorm(z1)
        d1 <- dnorm(z1)
        d2 <- NULL
        dr <- NULL
    }
    list(f = f, d1 = d1, d2 = d2, dr = dr)
}

dbnorm <- function(z1, z2, rho)
    1 / (2 * pi * sqrt(1 - rho ^ 2)) * exp(- 0.5 * (z1 ^ 2 + z2 ^ 2 - 2 * rho * z1 * z2) / (1 - rho ^ 2))

PHI3 <- function(z1, z2, z3, rho){
    h1 <- ! is.null(z1)
    h3 <- ! is.null(z3)
    if (h1 & h3){
        if (all(rho == 0)) PT <- pnorm(z1) * pnorm(z2) * pnorm(z3) else PT <- ptnorm(z1, z2, z3, rho)
        
        d1 <- dnorm(z1) * pbnorm((z2 - rho[1] * z1) / sqrt( (1 - rho[1] ^ 2) ),
                                 (z3 - rho[2] * z1) / sqrt( (1 - rho[2] ^ 2) ),
                                 (rho[3] - rho[1] * rho[2]) / sqrt(1 - rho[1] ^ 2) / sqrt(1 - rho[2] ^ 2))
        
        d2 <- dnorm(z2) * pbnorm((z1 - rho[1] * z2) / sqrt( (1 - rho[1] ^ 2) ),
                                 (z3 - rho[3] * z2) / sqrt( (1 - rho[3] ^ 2) ),
                                 (rho[2] - rho[1] * rho[3]) / sqrt(1 - rho[1] ^ 2) / sqrt(1 - rho[3] ^ 2))
        
        d3 <- dnorm(z3) * pbnorm((z1 - rho[2] * z3) / sqrt( (1 - rho[2] ^ 2) ),
                                 (z2 - rho[3] * z3) / sqrt( (1 - rho[3] ^ 2) ),
                                 (rho[1] - rho[2] * rho[3]) / sqrt(1 - rho[2] ^ 2) / sqrt(1 - rho[3] ^ 2))
        
        S <- 1 - rho[1] ^ 2 - rho[2] ^ 2  - rho[3] ^ 2 + 2 * rho[1] * rho[2] * rho[3]
        if (S < 0) stop(paste("the coefficients of correlations are out of range", paste(round(rho, 4), collapse = ","), sep = " "))

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
    if (h1 & ! h3){
        if (rho == 0) PT <- pnorm(z1) * pnorm(z2) else PT <- pbnorm(z1, z2, rho)
        d1 <- dnorm(z1) * punorm( (z2 - rho * z1) / sqrt(1 - rho ^ 2) )
        d2 <- dnorm(z2) * punorm( (z1 - rho * z2) / sqrt(1 - rho ^ 2) )
        d3 <- NULL
        dr <- dnorm(z1) *  dnorm( (z2 - rho * z1) / sqrt(1 - rho ^ 2) ) / sqrt(1 - rho ^ 2)
    }
    if (! h1 & h3){
        if (rho == 0) PT <- pnorm(z2) * pnorm(z3) else PT <- pbnorm(z2, z3, rho)
        d1 <- NULL
        d2 <- dnorm(z2) * punorm( (z3 - rho * z2) / sqrt(1 - rho ^ 2) )
        d3 <- dnorm(z3) * punorm( (z2 - rho * z3) / sqrt(1 - rho ^ 2) )
        dr <- dnorm(z3) *  dnorm( (z2 - rho * z3) / sqrt(1 - rho ^ 2) ) / sqrt(1 - rho ^ 2)
    }
    if (! h1 & ! h3){
        PT <- pnorm(z2)
        d1 <- NULL
        d2 <- dnorm(z2)
        d3 <- NULL
        dr <- NULL
    }
    list(f = PT, d1 = d1, d2 = d2, d3 = d3, dr = dr)
}

## 06 92 94 55 53
## 06 93 30 21 08
