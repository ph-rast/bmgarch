##' @keywords internal
##' @importFrom stats sd
.colSDs <- function(x) {
    lapply(x, function(x) {
        dims <- dim(x)
        apply(x, 2:length(dims), sd)
    })
}

##' Obtain quantiles over columns in lists
##' @title Quantiles within lists
##' @param x 
##' @param probs Quantile(s). Inherits from \code{forecast} which defaults to \code{c(.025, .975)}.
##' @return Quantiles at the column level within lists
##' @author philippe
##' @keywords internal
.colQTs <- function(x, probs ) {
    lapply(x, function(x) {
        dims <- dim(x)
        apply(x, 2:length(dims), quantile, probs)
    })
}

##' Computes posterior summaries for all parameters of interest for bmgarch objects.
##'
##' @title Summary method for bmgarch objects.
##' @param object bmgarch object.
##' @param CrI Numeric vector (Default: \code{c(.025, .975)}). Lower and upper bound of predictive credible interval.
##' @param digits Integer (Default: 2, optional). Number of digits to round to when printing.
##' @param ... Not used.
##' @return summary.bmgarch object. A named list containing "meta" and "model_summary". \code{model_summary} contains summary table for all model parameters.
##' @author Stephen R. Martin, Philippe Rast
##' @export
summary.bmgarch <- function(object, CrI = c(.025, .975), digits = 2, ...) {
    
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

    out$model_summary <- .get_stan_summary(object$model_fit, params, CrI)

    class(out) <- "summary.bmgarch"
    return(out)
}

##' @title Get stan summaries.
##' @param model_fit stanfit object or list of stanfit objects.
##' @param params Character vector. Names of params to pull from stan summary.
##' @param CrI Numeric vector (length 2).
##' @param weights Numeric vector. Weights for each model in model_fit, if list.
##' @return Stan summary for parameters. Columns: mean, sd, mdn, and CrIs.
##' @author Stephen R. Martin, Philippe Rast
##' @keywords internal
.get_stan_summary <- function(model_fit, params, CrI, weights = NULL) {
    if(class(model_fit) == "stanfit" | (class(model_fit) == "list" & length(model_fit) == 1)) { # One model fit
        if(class(model_fit) == "list") {
            model_fit <- model_fit[[1]]
        }
        CrI <- c(.5, CrI)
        cols <- c("mean","sd",paste0(CrI*100, "%"), "n_eff", "Rhat")
        model_summary <- rstan::summary(model_fit, pars = params, probs = CrI)$summary[,cols]
        colnames(model_summary)[colnames(model_summary) == "50%"] <- "mdn"
        return(model_summary)
    } else { # List of stanfits
        if(length(model_fit) != length(weights)) {
            stop("Weights must be same length as model_fit list.")
        }

        samps_comb <- .weighted_samples(model_fit, params, weights)

        # Model summaries
        m <- Map(colMeans, samps_comb)
        SD <- .colSDs(samps_comb)
        qts <- .colQTs(samps_comb, c(.5, CrI))

        # Permutate
        ## Stan's summary() goes in right-fastest ([1, 1, 1] -> [1, 1, 2])
        ## Stan's extract/matrix/array goes in left-fastest ([1, 1, 1] -> [2, 1, 1])
        ## We want the wt'd method to return same order as stan's summary.
        ## R makes vectors in left-fastest; this reverses indexing, then fills a summary matrix.
        dims <- lapply(samps_comb, function(p) { # Each parameter
            dim(p)[2:length(dim(p))]
        })
        num_dims <- lapply(dims, length)

        out <- list()
        for(p in seq_len(length(num_dims))) { # Each parameter
            m[[p]] <- aperm(m[[p]], rev(1:num_dims[[p]])) # Right-most varies fastest
            SD[[p]] <- aperm(SD[[p]], rev(1:num_dims[[p]]))
            qts[[p]] <- aperm(qts[[p]], rev(1:(num_dims[[p]] + 1)))
            out[[p]] <- cbind(as.numeric(m[[p]]),
                              as.numeric(SD[[p]]),
                              matrix(as.numeric(qts[[p]]), ncol = 3),
                              n_eff = NA,
                              Rhat = NA)
            colnames(out[[p]]) <- c("mean", "sd", "mdn", paste0(CrI * 100, "%"), "n_eff", "Rhat")
        }
        names(out) <- names(m)

        # Recreate row names
        for(p in seq_len(length(num_dims))) {
            ind_mat <- do.call(expand.grid, rev(lapply(dims[[p]], seq_len)))
            ind_mat <- ind_mat[,rev(seq_len(ncol(ind_mat)))]
            prefix <- names(out)[p]
            inner <- apply(ind_mat, 1, paste0, collapse = ",")
            rownames(out[[p]]) <- paste0(prefix, "[", inner, "]")
        }
        # Combine all to single data.frame
        out <- as.matrix(do.call(rbind, out))
        return(out)
    }
}

.weighted_samples <- function(model_fit, params, weights) {
    ##################
    # Extract method #
    ##################
    samps <- lapply(model_fit, rstan::extract, pars = params)
    # Apply weights
    for(i in seq_len(length(samps))) { # Each model
        samps[[i]] <- lapply(samps[[i]], function(p) { # Each parameter
            p * weights[i] # Weight them
        })
    }
    # Reduce
    samps_comb <- lapply(params, function(p) { # For each parameter
        Reduce("+", lapply(samps, function(m) {m[[p]]})) # Sum samples together
    })
    names(samps_comb) <- params
    return(samps_comb)
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
.print.summary.ccc <- function(bmsum, ...) {
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
    ms <- ms[!grepl("c_h\\[", rownames(ms)),] # Remove c_h; will print c_h_var
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
    ms <- ms[!grepl("c_h\\[", rownames(ms)),] # Remove c_h; will print c_h_var
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
    ms <- ms[!grepl("c_h\\[", rownames(ms)),] # Remove c_h; will print c_h_var
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
    nt <- bmsum$meta$nt
    TS_names <- bmsum$meta$TS_names

    intercept_index <- grep("phi0", rownames(ms))
    ar_index <- grep("phi\\[", rownames(ms))
    ma_index <- grep("theta\\[", rownames(ms))

    ar_indices <- gsub("phi\\[([[:digit:]]+),([[:digit:]]+)]", "\\1,\\2", rownames(ms)[ar_index])
    ma_indices <- gsub("theta\\[([[:digit:]]+),([[:digit:]]+)]", "\\1,\\2", rownames(ms)[ma_index])
    ar_indices <- strsplit(ar_indices, ",")
    ma_indices <- strsplit(ma_indices, ",")
    ar_indices <- apply(do.call(rbind, ar_indices), 1:2, as.numeric)
    ma_indices <- apply(do.call(rbind, ma_indices), 1:2, as.numeric)

    rownames(ms)[intercept_index] <- paste0("(Intercept)_", bmsum$meta$TS_names)
    rownames(ms)[ar_index] <- paste0("Phi_", TS_names[ar_indices[,1]], "-", TS_names[ar_indices[,2]])
    rownames(ms)[ma_index] <- paste0("Theta_", TS_names[ma_indices[,1]], "-", TS_names[ma_indices[,2]])

    if(bmsum$meta$meanstructure == 0) {
        arma_index <- intercept_index
        msg <- "Intercept estimates on the location:"
    } else {
        arma_index <- c(intercept_index, ar_index, ma_index)
        msg <- "ARMA(1,1) estimates on the location:"
    }
    cat(msg)
    .newline(2)
    print(round(ms[arma_index,], digits = bmsum$meta$digits))
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
