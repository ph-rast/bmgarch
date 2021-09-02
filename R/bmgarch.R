#' Standardize input data to facilitate computation
#' 
#' @param data Time-series data
#' @param xC Numeric vector or matrix.
#' @param P Numeric.
#' @param Q Numeric.
#' @param standardize_data Logical.
#' @param distribution Character.
#' @param meanstructure Character.
#' @return bmgarch stan data list. 
#' @keywords internal
standat <- function(data, xC, P, Q, standardize_data, distribution, meanstructure){

    if(dim(data)[1] < dim(data)[2]) {
        data = t(data)
        warning("data is wider than it is long. Transposing...")
    }
    if ( is.null( colnames( data ) ) ) colnames( data ) = paste0('t', 1:ncol( data ) )

    ## Model for meanstructure
    if( meanstructure == "constant" | meanstructure == 0 ) {
        meanstructure <- 0
    } else if ( meanstructure == "arma" | meanstructure == 1 ){
        meanstructure <- 1
    } else {
        stop("meanstructure must be either 'constant' or 'arma'.")
    }

    ## Test that all data vectors have variance > 0
    ## Stop if a vector as zero variance
    dvar = apply(data, 2, var)
    if( sum(ifelse(dvar == 0, 1, 0)) > 0 ) stop(
                  paste0("Datavector ", names(which(dvar == 0 ) ), " has zero variance.") )
                                                       
    ## Tests on predictor
    ## Pass in a 0 matrix, so that stan does not complain
    if ( is.null(xC) ) {
        xC = matrix(0, nrow = nrow(data), ncol = ncol(data))
    }
    ## Match dimension of predictor to TS. If only one vector is given, it's assumed that it is the same for all TS's
    if ( is.null(ncol(xC)) ) {
        warning("xC is assumed constant across TS's")
        xC <- matrix(xC, nrow = nrow(data), ncol = ncol(data)) ## Risky, better to throw an error
    } else if ( dim( xC )[2] != dim( data )[2] ) { ## xC is not a vector  - check if it is of right dimension
        warning("xC is not of right dimension - adapt xC dimension to match number of TS")
    }

    if( standardize_data ) {
    ## Standardize time-series
        stdx <- scale(data)
        centered_data <- attr(stdx, "scaled:center")
        scaled_data <- attr(stdx, "scaled:scale")
        return_standat <- list(T = nrow(stdx),
                            rts = stdx,
                            xC = xC,
                            nt = ncol(stdx),
                            centered_data = centered_data,
                            scaled_data = scaled_data,
                            distribution = distribution,
                            P = P,
                            Q = Q,
                            meanstructure = meanstructure)
        } else {
        ## Unstandardized
        return_standat <- list(T = nrow(data),
                                rts = data,
                                xC = xC,
                                nt = ncol(data),
                                distribution = distribution,
                                P = P,
                                Q = Q,
                                meanstructure = meanstructure)
    }
    return(return_standat)
}

##' Draw samples from a specified multivariate GARCH model using 'Stan', given multivariate time-series. Currently supports CCC, DCC, BEKK, and pdBEKK model parameterizations.
##'
##' Four types of paramerizations are implemented. The constant conditional correlation (CCC) and the dynamic conditional correlation \insertCite{@DCC; Engle2002,Engle2001a}{bmgarch}, as well as  BEKK \insertCite{Engle1995}{bmgarch} and a BEKK model with positivity constraints on the diagonals of the ARCH and GARCH parameters "pdBEKK" \insertCite{Rast2020}{bmgarch}.
##'
##' The fitted models are 'rstan' objects and all posterior parameter estimates can be obtained and can be examined with either the 'rstan' toolbox, plotted and printed using generic functions  or passed to 'bmgarch' functions to 'forecast' or compute 'model_weights' or compute fit statistics based on leave-future-out cross-validation. 
##' 
##' @title Estimate Bayesian Multivariate GARCH
##' @param data Time-series or matrix object. A time-series or matrix object containing observations at the same interval.
##' @param xC Numeric vector or matrix. Covariates(s) for the constant variance terms in C, or c, used in a log-linear model on the constant variance terms \insertCite{Rast2020}{bmgarch}. If vector, then it acts as a covariate for all constant variance terms. If matrix, must have columns equal to number of time series, and each column acts as a covariate for the respective time series (e.g., column 1 predicts constant variance for time series 1).
##' @param parameterization Character (Default: "CCC"). The type of of parameterization. Must be one of "CCC", "DCC", "BEKK", or "pdBEKK".
##' @param P Integer. Dimension of GARCH component in MGARCH(P,Q).
##' @param Q Integer. Dimension of ARCH component in MGARCH(P,Q).
##' @param iterations Integer (Default: 2000). Number of iterations for each chain (including warmup).
##' @param chains Integer (Default: 4). The number of Markov chains.
##' @param standardize_data Logical (Default: FALSE). Whether data should be standardized to easy computations. 
##' @param distribution Character (Default: "Student_t"). Distribution of innovation: "Student_t"  or "Gaussian"
##' @param meanstructure Character (Default: "constant"). Defines model for means. Either 'constant'  or 'arma'. Currently arma(1,1) only.
##' @param ... Additional arguments can be ‘chain_id’, ‘init_r’, ‘test_grad’, ‘append_samples’, ‘refresh’, ‘enable_random_init’ etc. See the documentation in \code{\link[rstan]{stan}}.
##' @return \code{bmgarch} object.
##' @importFrom Rdpack reprompt
##' @author Philippe Rast, Stephen R. Martin
##' @references
##'    \insertAllCited()
##' @export
##' @examples
##' \dontrun{
##' data(panas)
##' # Fit BEKK(1,1) mgarch model with a ARMA(1,1) meanstructure,
##' # and student-t residual distribution
##' fit <- bmgarch(panas, parameterization = "BEKK",
##'                P = 1, Q = 1,
##'                meanstructure = "arma",
##'                distribution = "Student_t")
##'
##' # Summarize the parameters
##' summary(fit)
##'
##' # Forecast 5 ahead
##' fit.fc <- forecast(fit, ahead = 5)
##' print(fit.fc)
##'
##' # Plot mean forecasts
##' plot(fit.fc, type = "mean")
##' 
##' # Plot variance forecasts
##' plot(fit.fc, type = "var")
##' 
##' # Plot correlation forecasts
##' plot(fit.fc, type = "cor")
##'
##' # Plot modeled data ("backcasted values").
##' plot(fit, type = "mean")
##' 
##' # Save "backcasted" values
##' fit.bc <- fitted(fit)
##'
##' # Save estimated and forecasted data as a data.frame
##' df.fc <- as.data.frame(fit.fc)
##'
##' # Access rstan's model fit object
##' mf <- fit$model_fit
##'
##' # Return diagnostics and a plot of the first 10 parameters
##' rstan::check_hmc_diagnostics(mf)
##' rstan::plot(mf)
##' }
bmgarch <- function(data,
                   xC = NULL,
                   parameterization = "CCC",
                   P = 1,
                   Q = 1,
                   iterations = 2000,
                   chains = 4,
                   standardize_data = FALSE,
                   distribution = "Student_t",
                   meanstructure = "constant", ...) {
    if ( tolower(distribution) == "gaussian" ) {
        num_dist <- 0
    } else if ( tolower(distribution) == "student_t" ) {
        num_dist <- 1
    } else {
        stop( "\n\n Specify distribution: Gaussian or Student_t \n\n")
    }

    return_standat <- standat(data, xC, P, Q,  standardize_data, distribution = num_dist, meanstructure )
    stan_data <- return_standat[ c("T", "xC", "rts", "nt", "distribution", "P", "Q", "meanstructure")]

    stanmodel <- switch(parameterization,
                        CCC = stanmodels$CCCMGARCH,
                        DCC = stanmodels$DCCMGARCH,
                        BEKK = stanmodels$BEKKMGARCH,
                        pdBEKK = stanmodels$pdBEKKMGARCH,
                        NULL)
    if(is.null(stanmodel)) {
        stop("Not a valid model specification. ",
             parameterization,
             "must be one of: ",
             paste0(supported_models, collapse = ", "),
             ".")
    }

    model_fit <- rstan::sampling(stanmodel,
                                 data = stan_data,
                                 verbose = TRUE,
                                 iter = iterations,
                                 control = list(adapt_delta = .99),
                                 chains = chains,
                                 init_r = .05)

    ## Model fit is based on standardized values.
    mns <- return_standat$centered_data
    sds <- return_standat$scaled_data
    ## Values could be converted to original scale using something like this on the estimates
    ## orig_sd = stan_data$rts %*% diag(sds)
    ## orig_scale = orig_sd + array(rep(mns, each = aussi[[1]]$T), dim = c(aussi[[1]]$T, aussi[[1]]$nt) )
    return_fit <- list(model_fit = model_fit,
                       param = parameterization,
                       distribution = distribution,
                       num_dist = num_dist,
                       iter = iterations,
                       chains = chains,
                       elapsed_time = rstan::get_elapsed_time(model_fit),
                       date = date(),
                       nt = stan_data$nt,
                       TS_length = stan_data$T,
                       TS_names = colnames(stan_data$rts),
                       RTS_last = stan_data$rts[stan_data$T,],
                       RTS_full = stan_data$rts,
                       mgarchQ = stan_data$Q,
                       mgarchP = stan_data$P,
                       xC = stan_data$xC,
                       meanstructure = stan_data$meanstructure,
                       std_data = standardize_data)
    class(return_fit) <- "bmgarch"
    return(return_fit)
}


#' Models supported by bmgarch
#'
#' To be used when checking whether a parameterization or object type is a supported type.
#' May facilitate more parameterizations, as we only have to update these, and the switch statements.
#' @keywords internal
#' @author Stephen R. Martin
supported_models <- c("DCC", "CCC", "BEKK", "pdBEKK")
