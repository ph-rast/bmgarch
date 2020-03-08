##' Computes posterior summaries for all parameters of interest for bmgarch objects.
##'
##' @title Summary method for bmgarch objects.
##' @param object bmgarch object.
##' @param CrI Numeric vector (Default: \code{c(.025, .975)}). Lower and upper bound of predictive credible interval.
##' @param digits Integer (Default: 2, optional). Number of digits to round to when printing.
##' @return summary.bmgarch object. A named list containing "meta" and "model_summary". \code{model_summary} contains summary table for all model parameters.
##' @author Stephen R. Martin, Philippe Rast
##' @export
summary.bmgarch <- function(object, CrI = c(.025, .975), digits = 2) {
    
    # Parameters for each model
    common_params <- c("lp__", "nu")
    arma_params <- c("phi0", "phi", "theta")

    ccc_params <- c("c_h", "a_h", "b_h", "R", "beta", "c_h_var")
    dcc_params <- c("c_h", "a_h", "b_h", "beta", "c_h_var", "a_q", "b_q", "S")
    bekk_params <- c("C_R", "C_var", "A", "B", "beta1")
    
    # Meta-data needed for printing
    # TODO: Revisit this; some can be removed. Kitchen sink for now.
    metaNames <- c("param", "distribution", "num_dist", "iter", "chains", "elapsed_time", "date", "nt", "TS_names", "mgarchQ", "mgarchP", "meanstructure")
    meta <- with(object, mget(metaNames))
    meta$xC <- !all(object$xC == 0)
    out <- list()
    out$meta <- meta
    out$meta$digits <- digits

    # Get the model summaries. print.summary.bmgarch will process this + meta.
    params <- switch(object$param,
                     CCC = c(ccc_params, arma_params, common_params),
                     DCC = c(dcc_params, arma_params, common_params),
                     BEKK = c(bekk_params, arma_params, common_params),
                     pdBEKK = c(bekk_params, arma_params, common_params),
                     NULL
                     )
    if(is.null(params)) {
        stop("bmgarch object 'param' does not match a supported model. ",
             object$param, "is not one in ", paste0(supported_models, collapse = ", "), ".")
    }

    out$model_summary <- .get_stan_summary(object, params, CrI)

    class(out) <- "summary.bmgarch"
    return(out)
}

##' @title Get stan summaries.
##' @param object bmgarch object.
##' @param params Character vector. Names of params to pull from stan summary.
##' @param CrI Numeric vector (length 2).
##' @return Stan summary for parameters. Columns: mean, sd, mdn, and CrIs.
##' @author Stephen R. Martin, Philippe Rast
##' @keywords internal
.get_stan_summary <- function(object, params, CrI) {
    CrI <- c(.5, CrI)
    cols <- c("mean","sd",paste0(CrI*100, "%"), "n_eff", "Rhat")
    model_summary <- rstan::summary(object$model_fit, pars = params, probs = CrI)$summary[,cols]
    colnames(model_summary)[colnames(model_summary) == "50%"] <- "mdn"
    return(model_summary)
}

##' @title Print method for bmgarch.summary objects.
##' @param x summary.bmgarch object.
##' @param ... Not used.
##' @return x (invisible).
##' @author Stephen R. Martin
##' @export
print.summary.bmgarch <- function(x, ...) {
    if(x$meta$param == "CCC") {
        .print.summary.ccc(x)
    } else if(x$meta$param == "DCC") {
        .print.summary.dcc(x)
    } else if(x$meta$param %in% c("BEKK", "pdBEKK")) {
        .print.summary.bekk(x)
    }
    .newline(2)

    .print.summary.arma(x)
    .newline(2)

    if(x$meta$xC) {
        .print.summary.beta(x)
        .newline(2)
    }

    if(x$meta$num_dist == 1) {
        .print.summary.nu(x)
        .newline(2)
    }

    .print.summary.lp(x)
    .newline()

    return(invisible(x))
}
##' @title Print helper for CCC.
##' @param bmsum summary.bmgarch object.
##' @param x bmgarch.summary object.
##' @param ... Not used.
##' @return Void.
##' @author Stephen R. Martin, Philippe Rast
##' @keywords internal
.print.summary.ccc <- function(bmsum) {
    # Meta-data
    cat("Model:", paste0(bmsum$meta$param, "-MGARCH\n"))
    cat("Basic Specification: ")
    cat("H_t = D_t R D_t")
    .newline()
    .tab()
    cat("diag(D_t) = sqrt(h_[ii,t]) = c_h + a_h*y^2_[t-1] + b_h*h_[ii, t-1")
    .newline()

    # Sampling configuration
    .print.config(bmsum)

    # Get indices for params
    ms <- bmsum$model_summary
    ms <- ms[!grepl("c_h$", rownames(ms)),] # Remove c_h; will print c_h_var
    garch_h_index <- grep("_h", rownames(ms))
    cond_corr_index <- grep("R", rownames(ms))

    # Settings
    nt <- bmsum$meta$nt
    P <- bmsum$meta$mgarchP
    Q <- bmsum$meta$mgarchQ
    digits <- bmsum$meta$digits

    # Shortened TS names, if needed.
    short_names <- abbreviate(bmsum$meta$TS_names, minlength = 2)
    # Off-diagonals
    cormat_index <- matrix(1:(nt*nt), nrow = nt)
    corr_only <- cormat_index[lower.tri(cormat_index)]
    diag_only <- diag(cormat_index)
    ## obtain all combinations of TS varnames for A and B in BEKK
    full_varnames <- expand.grid( short_names, short_names)
    ## obtain off-diagonal TS varnames
    od_varnames <- full_varnames[corr_only, ]

    cat(paste0("GARCH(", P, ",", Q, ")"), " estimates for conditional variance:")
    .newline(2)


    #########
    # GARCH #
    #########

    # Change out GARCH param labels
    rn <- rownames(ms[garch_h_index,])
    for(i in 1:nt) {
        replace <- grep(paste0(as.character(i), "\\]"), rn)
        rn[replace] <- gsub(paste0(as.character(i), "\\]"), paste0(short_names[i]), rn)[replace]
        rn <- gsub("\\[", "_", rn)
    }
    garch_out <- ms[garch_h_index,]
    rownames(garch_out) <- rn

    ## Print GARCH params
    print(round(garch_out, digits = digits))
    .newline(2)

    #####
    # R #
    #####
    cat("Constant correlation (R) coefficients:")
    .newline(2)

    corr_out <- ms[cond_corr_index[corr_only],]
    if( nt == 2 ) {
        tmp <- matrix(corr_out, nrow = 1)
        rownames(tmp) <- paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
        colnames(tmp) <- names(corr_out)
        corr_out <- tmp
    } else {
        rownames(corr_out) <- paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
    }
    print(round(corr_out, digits = digits))
}

##' @title Print helper for DCC.
##' @param bmsum summary.bmgarch object.
##' @return Void.
##' @author Stephen R. Martin, Philippe Rast
##' @keywords internal
.print.summary.dcc <- function(bmsum) {
    # Meta-data
    cat("Model:", paste0(bmsum$meta$param, "-MGARCH\n"))
    cat("Basic Specification: ")
    cat("H_t = D_t R D_t")
    .newline()
    .tab()
    cat("diag(D_t) = sqrt(h_ii,t) = c_h + a_h*y^2_[t-1] + b_h*h_[ii,t-1]")
    .newline()
    .tab()
    cat("R_t = Q^[-1]_t Q_t Q^[-1]_t = ( 1 - a_q - b_q)S + a_q(u_[t-1]u'_[t-1]) + b_q(Q_[t-1])")
    .newline()

    # Sampling configuration
    .print.config(bmsum)

    # Get indices for params
    ms <- bmsum$model_summary
    ms <- ms[!grepl("c_h$", rownames(ms)),] # Remove c_h; will print c_h_var
    garch_h_index <- grep("_h", rownames(ms))
    garch_q_index  <- grep("_q", rownames(ms) )
    cond_corr_index <- grep("R", rownames(ms))
    S_index = grep("S", rownames(ms))

    # Settings
    nt <- bmsum$meta$nt
    P <- bmsum$meta$mgarchP
    Q <- bmsum$meta$mgarchQ
    digits <- bmsum$meta$digits

    # Shortened TS names, if needed.
    short_names <- abbreviate(bmsum$meta$TS_names, minlength = 2)
    # Off-diagonals
    cormat_index <- matrix(1:(nt*nt), nrow = nt)
    corr_only <- cormat_index[lower.tri(cormat_index)]
    diag_only <- diag(cormat_index)
    ## obtain all combinations of TS varnames for A and B in BEKK
    full_varnames <- expand.grid( short_names, short_names)
    ## obtain off-diagonal TS varnames
    od_varnames <- full_varnames[corr_only, ]

    #########
    # GARCH #
    #########
    cat(paste0("GARCH(", P, ",", Q, ")"), " estimates for conditional variance on D:")
    .newline(2)

    rn <- rownames(ms[garch_h_index,])
    for ( i in 1:nt ) {
        replace <- grep(paste0( as.character(i), "\\]"), rn) 
        rn[replace] <- gsub(paste0( as.character(i), "\\]" ), paste0(short_names[i]), rn)[replace]
    }
    rn <- gsub("\\[", "_", rn)

    ## Save into new object to change rownames
    garch_h_out <- ms[garch_h_index,]
    ## Assign new rownames
    rownames(garch_h_out) <- rn
    ## print garch parameters
    print(round( garch_h_out, digits = digits) )
    .newline(2)

    #####
    # Q #
    #####
    cat("GARCH(1,1) estimates for conditional variance on Q:")
    .newline(2)
    rn = rownames(ms[garch_q_index,])
    for ( i in 1:nt ) {
        replace = grep(paste("\\[", "\\]", sep=as.character(i)), rn)
        rn[replace] = gsub(paste("\\[", "\\]", sep=as.character(i)), paste0("_", short_names[i]), rn)[replace]
    }

    ## Save into new object to change rownames
    garch_q_out = ms[garch_q_index,]
    ## Assign new rownames
    rownames(garch_q_out) = rn
    ## print garch parameters
    print(round( garch_q_out, digits = digits) )
    .newline(2)

    cat("Unconditional correlation 'S' in Q:")
    .newline(2)
    S_out <- ms[S_index[corr_only],]
    if (nt == 2) {
        tmp <- matrix( S_out, nrow = 1 )
        rownames(tmp) <- paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
        colnames(tmp) <- names(S_out)
        S_out <- tmp 
    } else {
        rownames(S_out) <- paste( paste0("S_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
    }

    print(round(S_out, digits = digits))
}

##' @title Print helper for BEKK/pdBEKK.
##' @param bmsum summary.bmgarch object.
##' @return Void.
##' @author Stephen R. Martin, Philippe Rast
##' @keywords internal
.print.summary.bekk <- function(bmsum) {
    # Meta-data
    cat("Model:", paste0(bmsum$meta$param, "-MGARCH\n"))
    cat("Basic Specification: ")
    cat("H_t = D_t R D_t")
    .newline()
    cat("H_t = C + A'[y_(t-1)*y'_(t-1)]A + B'H_(t-1)B")
    .newline()

    # Sampling configuration
    .print.config(bmsum)

    # Get indices for params
    ms <- bmsum$model_summary
    ms <- ms[!grepl("c_h$", rownames(ms)),] # Remove c_h; will print c_h_var
    garch_h_index <- grep("_h", rownames(ms))
    garch_q_index  <- grep("_q", rownames(ms) )
    cond_corr_index <- grep("R", rownames(ms))
    S_index = grep("S", rownames(ms))

    # Settings
    nt <- bmsum$meta$nt
    P <- bmsum$meta$mgarchP
    Q <- bmsum$meta$mgarchQ
    digits <- bmsum$meta$digits

    # Shortened TS names, if needed.
    short_names <- abbreviate(bmsum$meta$TS_names, minlength = 2)
    # Off-diagonals
    cormat_index <- matrix(1:(nt*nt), nrow = nt)
    corr_only <- cormat_index[lower.tri(cormat_index)]
    diag_only <- diag(cormat_index)
    ## obtain all combinations of TS varnames for A and B in BEKK
    full_varnames <- expand.grid( short_names, short_names)
    ## obtain off-diagonal TS varnames
    od_varnames <- full_varnames[corr_only, ]

    garch_C_index <- grep("beta0", rownames(ms))
    garch_Cv_index <- grep("C_var", rownames(ms))
    garch_R_index <- grep("C_R", rownames(ms))
    garch_A_index <- grep("A", rownames(ms))
    garch_B_index <- grep("B", rownames(ms)) 

    ########
    # BEKK #
    ########
    .sep()
    cat("Constant correlation, R (diag[C]*R*diag[C]):")
    .newline(2)
    R_out <- ms[garch_R_index[corr_only],]
    if(nt == 2) {
        tmp <- matrix( R_out, nrow = 1 )
        rownames(tmp) <- paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
        colnames(tmp) <- names(R_out)
        R_out = tmp 
    } else {
        rownames(R_out) <- paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
    }

    print(round(R_out, digits = digits))
    .newline(2)

    cat("Constant variances (diag[C]):")
    .newline(2)
    C_out <- ms[garch_Cv_index,]
    if ( nt == 2 ) {
        tmp <- matrix( C_out, nrow = 2 )
        rownames(tmp) <- paste0("var_", short_names )
        colnames(tmp) <- colnames(C_out)
        C_out <- tmp
    } else {
        rownames(C_out) = paste0("var_", short_names )
    }
    print(round( C_out, digits = digits) )
    .newline(2)

    #####
    # A #
    #####
    cat(paste0(paste0("MGARCH(", P, ",", Q, ')')), "estimates for A:")
    .newline(2)

    a <- list()
    A_out <- ms[garch_A_index,]
    if (Q > 1) {
        for( q in 1:Q ){
            a[[q]] = paste( paste0( paste0("A_", q, "_"), full_varnames[ ,1] ), full_varnames[ ,2] , sep = '-')
        }
        rownames(A_out) = unlist( a )
    } else  {
        rownames(A_out) = paste( paste0("A_", full_varnames[ ,1] ), full_varnames[ ,2] , sep = '-')   
    }
    print(round( A_out, digits = digits) )
    .newline(2)

    #####
    # B #
    #####
    cat(paste0(paste0("MGARCH(", P, ",", Q, ')')), "estimates for B:")
    .newline(2)

    b <- list()
    B_out <- ms[garch_B_index,]
    if (P > 1) {
        for( p in 1:P ){
            b[[p]] <- paste( paste0( paste0("B_", p, "_"), full_varnames[ ,1] ), full_varnames[ ,2] , sep = '-')
        }
        rownames( B_out ) <- unlist( b )
    } else {
        rownames( B_out ) <- paste( paste0("B_", full_varnames[ ,1] ), full_varnames[ ,2] , sep = '-')   
    }
    print(round( B_out, digits = digits) )
}

##' @title Print helper for ARMA component.
##' @param bmsum summary.bmgarch object.
##' @return Void.
##' @author Stephen R. Martin, Philippe Rast
##' @keywords internal
.print.summary.arma <- function(bmsum) {
    ms <- bmsum$model_summary
    if(bmsum$meta$meanstructure == 0) {
        arma_index <- grep("phi0", rownames(ms))
    } else {
        arma_index <- grep("^phi|^theta", rownames(ms))
    }
    cat("ARMA(1,1) estimates on the location:")
    .newline(2)
    print(round(bmsum$model_summary[arma_index,], digits = bmsum$meta$digits))
}

##' @title Print helper for beta component.
##' @param bmsum summary.bmgarch object.
##' @return Void.
##' @author Stephen R. Martin, Philippe Rast
##' @keywords internal
.print.summary.beta <- function(bmsum) {

    ms <- bmsum$model_summary
    # Settings
    nt <- bmsum$meta$nt
    P <- bmsum$meta$mgarchP
    Q <- bmsum$meta$mgarchQ
    digits <- bmsum$meta$digits
    # Shortened TS names, if needed.
    short_names <- abbreviate(bmsum$meta$TS_names, minlength = 2)
    # Off-diagonals
    cormat_index <- matrix(1:(nt*nt), nrow = nt)
    corr_only <- cormat_index[lower.tri(cormat_index)]
    diag_only <- diag(cormat_index)
    ## obtain all combinations of TS varnames for A and B in BEKK
    full_varnames <- expand.grid( short_names, short_names)
    ## obtain off-diagonal TS varnames
    od_varnames <- full_varnames[corr_only, ]
    

    if(bmsum$meta$param %in% c("BEKK", "pdBEKK")) {
        cat("Exogenous predictor (beta1 on log scale: C = sRs with s = exp( x*beta ):")
        .newline(2)
        beta_index <- grep("beta1", rownames(ms))
        beta <- ms[beta_index,]
        rownames(beta) <- paste0("beta1_", short_names)
        print(round(beta, digits = digits))
    } else {
        cat("Exogenous predictor (beta1 on log scale: c = exp( beta_0 + beta_1*x ):")
        .newline(2)

        beta0_index <- grep("c_h\\[", rownames(ms))
        beta0 <- ms[beta0_index,]
        rownames(beta0) <- paste0("beta0_", short_names)

        beta_index <- grep("beta", rownames(ms))
        beta <- ms[beta_index,]
        rownames(beta) <- paste0("beta_", short_names)
        betas <- rbind(beta0, beta)
        print(round(betas, digits = digits))
    }
}

##' @title Print helper for nu component.
##' @param bmsum summary.bmgarch object.
##' @return Void.
##' @author Stephen R. Martin, Philippe Rast
##' @keywords internal
.print.summary.nu <- function(bmsum) {
    nu <- bmsum$model_summary["nu",]
    cat("Df constant student_t (nu):")
    .newline(2)
    print(round(nu, digits = bmsum$meta$digits))
}

##' @title Print helper for LP component.
##' @param bmsum summary.bmgarch object.
##' @return Void.
##' @author Stephen R. Martin, Philippe Rast
##' @keywords internal
.print.summary.lp <- function(bmsum) {
    cat("Log density posterior estimate:")
    .newline(2)
    print(round(bmsum$model_summary["lp__",], digits = bmsum$meta$digits))
}

##' @title Print helper - Return new line(s).
##' @param n Integer (Default: 1). Number of new lines.
##' @return Prints new lines.
##' @author Stephen R. Martin
##' @keywords internal
.newline <- function(n = 1) {
    for(i in 1:n) {
        cat("\n")
    }
}

##' @title Print helper - tab
##' @return Prints tab.
##' @author Stephen R. Martin
##' @keywords internal
.tab <- function() {
    cat("\t")
}

##' @title Print helper - Separator, new line
##' @return Prints "---" and a new line.
##' @author Stephen R. Martin
##' @keywords internal
.sep <- function() {
    cat("---")
    .newline()
}

##' @title Print helper for Sampling Config.
##' @param bmsum summary.bmgarch object.
##' @return Void.
##' @author Stephen R. Martin, Philippe Rast
##' @keywords internal
.print.config <- function(bmsum) {
    .newline()
    cat("Distribution: ", bmsum$meta$distribution)
    .newline()
    .sep()
    cat("Iterations: ", bmsum$meta$iter)
    .newline()
    cat("Chains: ", bmsum$meta$chains)
    .newline()
    cat("Date: ", bmsum$meta$date)
    .newline()
    cat("Elapsed time (min): ", round((max(bmsum$meta$elapsed_time[,1]) + max(bmsum$meta$elapsed_time[,2]))/60, 2))
    .newline(2)
    
}

##' @title Summary method for BMGARCH object
##' @param object bmgarch object.
##' @param CrI Numeric vector (Default: \code{c(.025, .975)}). Lower and upper bound of predictive credible interval.
##' Possible values are .025, .05, .10, .50, .90, .95, and .975
##' @return Summary object
##' @author Philippe Rast
##' @export
## summary.bmgarch <- function(object, CrI = c(0.025, 0.975), digits = 2 ) {
##     cat( "Model:", paste0(object$param, "-MGARCH\n"))
##     if( object$param == 'CCC') {
##         cat("Basic specification: H_t = D_t R D_t", "\n")
##             cat("               diag(D_t) = sqrt(h_ii,t) = c_h + a_h*y^2_[t-1] + b_h*h_[ii,t-1]", "\n")
##     } else {
##         if( object$param == 'DCC') {
##             cat("Basic specification: H_t = D_t R_t D_t", "\n")
##             cat("               diag(D_t) = sqrt(h_ii,t) = c_h + a_h*y^2_[t-1] + b_h*h_[ii,t-1]", "\n")
##             cat("               R_t = Q^[-1]_t Q_t Q^[-1]_t = ( 1 - a_q - b_q)S + a_q(u_[t-1]u'_[t-1]) + b_q(Q_[t-1])", "\n")
##         } else {
##             if( object$param == 'BEKK' | object$param == 'pdBEKK') {
##                 cat("Basic specification: H_t = C + A'[y_(t-1)*y'_(t-1)]A + B'H_(t-1)B", "\n")
##             }
##         }
##     }
##     cat("\n")
##     cat("Distribution: ", object$distribution, "\n")
##     cat("---\n")
##     cat("Iterations: ", object$iter, "\n")
##     cat("Chains: ", object$chains, "\n")
##     cat("Date: ", object$date, "\n")
##     cat("Elapsed time (min):",
##     round((max(object$elapsed_time[,1]) + max(object$elapsed_time[,2]))/60, 2) , "\n\n")

##     ## Shared objects
##     nt <- object$nt
##     short_names <- abbreviate(object$TS_names, minlength = 2)
##     ## extract and printonly off-diagonal elements
##     cormat_index <- matrix(1:(nt*nt), nrow = nt)
##     corr_only <- cormat_index[lower.tri(cormat_index)]
##     diag_only <- diag(cormat_index)
##     ## obtain all combinations of TS varnames for A and B in BEKK
##     full_varnames <- expand.grid( short_names, short_names)
##     ## obtain off-diagonal TS varnames
##     od_varnames <- full_varnames[corr_only, ]
##     P <- object$mgarchP
##     Q <- object$mgarchQ

##     ## depending on parameteriztion, different parameters will be returned:
##     if(object$param == 'CCC') {
##         model_summary <- rstan::summary(object$model_fit,
##                                        pars = c('c_h', 'a_h', 'b_h', 'R', 'phi0', 'phi', 'theta', 'lp__'),
##                                        probs = CrI)$summary[, -2 ] # drop se_mean with [,-2]
        
##         garch_h_index  = grep("_h", rownames(model_summary) )
##         cond_corr_index = grep("R", rownames(model_summary) )
##         if( object$meanstructure == 0 ) {
##             arma_index = grep("phi0", rownames(model_summary))
##         } else {
##             arma_index = grep("^phi|^theta", rownames(model_summary) )
##         }
          
##         ##############################
##         ## GARCH parameters of CCC  ##
##         ##############################
##         cat(paste0(paste0("GARCH(", P, ",", Q, ')')), "estimates for conditional variance:", "\n\n") 
        
##         rn = rownames(model_summary[garch_h_index,])
##         for ( i in 1:nt ) {
##            # replace = grep(paste("\\[", "\\]", sep=as.character(i)), rn)
##             replace = grep(paste0( as.character(i), "\\]"), rn) 
##         rn[replace] = gsub(paste0( as.character(i), "\\]" ), paste0(short_names[i]), rn)[replace]
##         }
##         rn = gsub("\\[", "_", rn)

##         ## Save into new object to change rownames
##         garch_out = model_summary[garch_h_index,]
##         ## Assign new rownames
##         rownames(garch_out) = rn
##         ## print garch parameters
##         print(round( garch_out, digits = digits) )
##         cat("\n\n")

##         ###############################
##         ## Constant Correlation R ##
##         ###############################
##         cat("Constant correlation (R) coefficient(s):", "\n\n")
        
##         corr_out = model_summary[cond_corr_index[corr_only],]
##         if ( nt == 2 ) {
##             tmp = matrix( corr_out, nrow = 1 )
##             rownames(tmp) = paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
##             colnames(tmp) = names(corr_out)
##             corr_out = tmp } else {
##                                rownames(corr_out) =
##                                    paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
##                            }
##         print(round( corr_out, digits = digits) )
##         cat("\n\n")
        
##         ###############################
##         ## Location parameters       ##
##         ###############################
##         cat("ARMA(1,1) estimates on the location:", "\n\n")
##         print(round( model_summary[arma_index,], digits = digits) )
##         cat("\n\n") } else {
                        
##     #########
##     ## DCC ##
##     #########
##     if(object$param == 'DCC') {
##         model_summary <- rstan::summary(object$model_fit,
##                                        pars = c('a_q', 'b_q', 'c_h_var', 'a_h', 'b_h', 'S', 'phi0', 'phi', 'theta','beta', 'lp__'),
##                                        probs = CrI)$summary[, -2]
            
##         garch_h_index <- grep("_h", rownames(model_summary) )
##         garch_q_index  = grep("_q", rownames(model_summary) )
##         S_index = grep("S", rownames(model_summary) )
##         if( object$meanstructure == 0 ) {
##             arma_index = grep("phi0", rownames(model_summary))
##         } else {
##             arma_index = grep("^phi|^theta", rownames(model_summary) )
##         }

##         ## ###########################
##         ## GARCH parameters of DCC  ##
##         ## ###########################
##         cat(paste0(paste0("GARCH(", P, ",", Q, ')')), "estimates for conditional variance on D:", "\n\n") 

##         rn <- rownames(model_summary[garch_h_index,])
##         for ( i in 1:nt ) {
##             replace <- grep(paste0( as.character(i), "\\]"), rn) 
##         rn[replace] <- gsub(paste0( as.character(i), "\\]" ), paste0(short_names[i]), rn)[replace]
##         }
##         rn <- gsub("\\[", "_", rn)

##         ## Save into new object to change rownames
##         garch_h_out <- model_summary[garch_h_index,]
##         ## Assign new rownames
##         rownames(garch_h_out) <- rn
##         ## print garch parameters
##         print(round( garch_h_out, digits = digits) )
##         cat("\n\n")

##         ##########################
##         ## GARCH estimates on Q ##
##         ##########################
##         cat("GARCH(1,1) estimates for conditional variance on Q:", "\n\n")
##         rn = rownames(model_summary[garch_q_index,])
##         for ( i in 1:nt ) {
##             replace = grep(paste("\\[", "\\]", sep=as.character(i)), rn)
##             rn[replace] = gsub(paste("\\[", "\\]", sep=as.character(i)), paste0("_", short_names[i]), rn)[replace]
##         }

##         ## Save into new object to change rownames
##         garch_q_out = model_summary[garch_q_index,]
##         ## Assign new rownames
##         rownames(garch_q_out) = rn
##         ## print garch parameters
##         print(round( garch_q_out, digits = digits) )
##         cat("\n\n")

##         cat("Unconditional correlation 'S' in Q:", "\n\n")
##         S_out <- model_summary[S_index[corr_only],]
##         if ( nt == 2 ) {
##             tmp <- matrix( S_out, nrow = 1 )
##             rownames(tmp) <- paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
##             colnames(tmp) <- names(S_out)
##             S_out <- tmp } else {
##                                rownames(S_out) <- 
##                                    paste( paste0("S_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
##                            }
##         print(round( S_out, digits = digits) )
##         cat("\n\n")

##         ###############################
##         ## Location parameters       ##
##         ###############################       
##         cat("ARMA(1,1) estimates on the location:", "\n\n")
##         print(round( model_summary[arma_index,], digits = digits) )
##         cat("\n\n") } else {

##     #####################              
##     ## BEKK and pdBEKK ##
##     #####################
##     if(object$param == 'BEKK' | object$param == 'pdBEKK') {
##         model_summary <- rstan::summary(object$model_fit,
##                                         pars = c('beta0', 'C_var', 'A', 'B', 'C_R', 'phi0', 'phi', 'theta', 'lp__'),
##                                         probs = CrI)$summary[, -2 ] 

##         garch_C_index <- grep("beta0", rownames(model_summary) )
##         garch_Cv_index <- grep("C_var", rownames(model_summary) )
##         garch_R_index <- grep("C_R", rownames(model_summary) )                            
##         garch_A_index <- grep("A", rownames(model_summary) )
##         garch_B_index <- grep("B", rownames(model_summary) )                        
##         if( object$meanstructure == 0 ) {
##             arma_index <- grep("phi0", rownames(model_summary))
##         } else {
##             arma_index <- grep("^phi|^theta", rownames(model_summary) )
##         }

##         #######################
##         ## Garch parameters  ##
##         #######################
        
##         cat("---\n\n")

##         cat("Constant correlation, R (diag[C]*R*diag[C]):", "\n\n")
##         R_out <- model_summary[garch_R_index[corr_only],]
##         if ( nt == 2 ) {
##             tmp <- matrix( R_out, nrow = 1 )
##             rownames(tmp) <- paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ),
##                                   substr(od_varnames[ ,2], 1, 2) , sep = '-')
##             colnames(tmp) <- names(R_out)
##             R_out = tmp } else {
##                                rownames(R_out) <- 
##                                    paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ),
##                                          substr(od_varnames[ ,2], 1, 2) , sep = '-')
##                            }
##         print(round( R_out, digits = digits) )
##         cat("\n\n")

##         cat("Constant variances (diag[C]):", "\n\n")
##         C_out = model_summary[garch_Cv_index,]
##         if ( nt == 2 ) {
##             tmp = matrix( C_out, nrow = 2 )
##             rownames(tmp) = paste0("var_", short_names )
##             colnames(tmp) = colnames(C_out)
##             C_out = tmp } else {
##                             rownames(C_out) = paste0("var_", short_names )
##                            }
##         print(round( C_out, digits = digits) )
##         cat("\n\n")

##         #############
##         ## A and B ##
##         #############

##         #######
##         ## A ##
##         #######
##         cat(paste0(paste0("MGARCH(", P, ",", Q, ')')), "estimates for A:", "\n\n")

##         a = list()
##         A_out = model_summary[garch_A_index,]
##             if (Q > 1) {
##                 for( q in 1:Q ){
##                     a[[q]] = paste( paste0( paste0("A_", q, "_"), full_varnames[ ,1] ), full_varnames[ ,2] , sep = '-')
##                 }
##                 rownames(A_out) = unlist( a )
##             } else rownames(A_out) = paste( paste0("A_", full_varnames[ ,1] ), full_varnames[ ,2] , sep = '-')
##         print(round( A_out, digits = digits) )
##         cat("\n\n")
        
##         #######
##         ## B ##
##         #######
##         cat(paste0(paste0("MGARCH(", P, ",", Q, ')')), "estimates for B:", "\n\n")

##         b <- list()
##         B_out <- model_summary[garch_B_index,]
##         ##if ( nt == 2 ) {
##         if (P > 1) {
##             for( p in 1:P ){
##                 b[[p]] <- paste( paste0( paste0("B_", p, "_"), full_varnames[ ,1] ), full_varnames[ ,2] , sep = '-')
##             }
##             rownames( B_out ) <- unlist( b )
##         } else rownames( B_out ) <- paste( paste0("B_", full_varnames[ ,1] ), full_varnames[ ,2] , sep = '-')
##         print(round( B_out, digits = digits) )
##         cat("\n\n")
        
##         #####################
##         ## ARMA parameters ##
##         #####################
##         cat(paste0(paste(paste0("ARMA(", 1), 1, sep = ','), ')'), "estimates on the location:", "\n\n")
##         print(round( model_summary[arma_index,], digits = digits) )
##         cat("\n\n")
##         }
##     }
##    }

##     ## ######################
##     ## Beta: Predictor on H
##     ## #####################
##     if(all(object$xC == 0)) {
##         if(object$param == 'BEKK' | object$param == 'pdBEKK') {
##         cat("Exogenous predictor (beta1 on log scale: C = sRs with s = exp( x*beta ):", "\n\n")
##         beta <- rstan::summary(object$model_fit, pars = c('beta1'), probs = CrI)$summary[, -2]
##         rownames(beta) <- paste0( "beta1_", short_names )
##         print(round(beta, digits = digits ) )
##         cat("\n\n" )
##         } else {
##             cat("Exogenous predictor (beta1 on log scale: c = exp( beta_0 + beta_1*x ):", "\n\n")
##             beta0 <- rstan::summary(object$model_fit, pars = c('c_h'), probs = CrI)$summary[, -2]
##             rownames(beta0) <- paste0( "beta0_", short_names )
##             beta <- rstan::summary(object$model_fit, pars = c('beta'), probs = CrI)$summary[, -2]
##             rownames(beta) <- paste0( "beta_", short_names )
##             print(round(rbind(beta0,beta), digits = digits ) )
##             cat("\n\n" )
##         }
##     }
        
##     nu <- rstan::summary(object$model_fit, pars = c('nu'))$summary[,'mean']
##     Lnu <- round( rstan::summary(object$model_fit, pars = c('nu'), probs = CrI)$summary[,4], 2)
##     Unu <- round( rstan::summary(object$model_fit, pars = c('nu'), probs = CrI)$summary[,5], 2)

    
##     if( object$num_dist == 1) {
##         cat("Df constant student_t: nu =", round( nu, digits = 2), "\n")
##         cat("                nu 95%CrI [", Lnu, ";", Unu, "]","\n\n")
##     }

##     cat("Log density posterior estimate:", "\n\n")
##     print(round( model_summary[grep("lp__", rownames(model_summary) ),], digits = digits) )
## }
