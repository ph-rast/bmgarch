##' Simulates time series data from specified BEKK model.
##'
##' Simulates timeseries data from specified BEKK model.
##' Number of time series computed from the number of columns in C.
##' All matrices must be of the same dimension.
##' If ARMA parameters (phi, theta) unspecified (NULL), then assumes a constant mean of zero.
##' @title Simulate BEKK data.
##' @param N Integer. Length of time series.
##' @param C Numeric square matrix. Constant covariance matrix (C). Must be symmetric.
##' @param A Numeric square matrix. Moving average GARCH matrix (A).
##' @param B Numeric square matrix. Autoregressive ARCH matrix (B).
##' @param phi Numeric square matrix (Optional). Autoregressive coefficients (Phi).
##' @param theta Numeric square matrix (Optional). Moving average coefficients (Theta).
##' @return Matrix of observations.
##' @author Stephen R. Martin
##' @keywords internal
##' @importFrom stats rnorm
.sim.bekk <- function(N,C,A,B, phi = NULL, theta = NULL) {
    if(ncol(C) != nrow(C)){
        stop("C must be symmetric, square, PD.")
    }
    if(ncol(A) != nrow(A)){
        stop("A must be square.")
    }
    if(ncol(B) != nrow(B)){
        stop("B must be square.")
    }
    nt <- ncol(C)

    y <- array(0, dim = c(N, nt))
    y[1,] <- rnorm(nt, 0, sqrt(diag(C)))

    H <- array(0, dim = c(nt, nt, N))
    H[,,1] <- C

    for(i in 2:N) {
        H[,,i] <- C + t(A) %*% (t(y[i - 1,, drop = FALSE]) %*% y[i - 1,,drop = FALSE]) %*% A + t(B) %*% H[,,i-1] %*% B
        y[i,] <- MASS::mvrnorm(1, rep(0, nt), H[,,i])
    }

    if (!is.null(phi) & !is.null(theta)) {
        ## Assume phi0 (intercept) is zero.
        if (ncol(phi) != nrow(phi)) {
            stop("phi must be square [nt, nt].")
        }
        if (ncol(theta) != nrow(theta)) {
            stop("theta must be square [nt, nt].")
        }
        if (ncol(phi) != nt) {
            stop("phi must be square [nt, nt].")
        }
        if (ncol(theta) != nt) {
            stop("theta must be square [nt, nt].")
        }
        mu <- array(0, dim = c(N, nt))
        mu[1,] <- 0
        for(i in 2:N) {
            mu[i,] <- 10 + y[i - 1, , drop = FALSE] %*% phi + (y[i - 1, ,drop = FALSE] - mu[i - 1,,drop = FALSE])%*%theta
            y[i,] <- y[i,,drop = FALSE] + mu[i,,drop = FALSE]
        }
        ## y <- mu + y
    }

    return(y)
}
