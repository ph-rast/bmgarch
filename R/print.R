##' @title Summary for BMGARCH object
##' @param object 
##' @param CrI 
##' @return Summary object
##' @author philippe
##' @export
summary.bmgarch = function(object, CrI = c(0.05, 0.95), digits = 2 ) {
    cat( "Model:", paste0(object$param, "-MGARCH\n"))
    if( object$param == 'CCC') {
        cat("Basic specification: H_t = D_t R D_t", "\n")
            cat("               diag(D_t) = sqrt(h_ii,t) = c_h + a_h*y^2_[t-1] + b_h*h_[ii,t-1]", "\n")
    } else {
        if( object$param == 'DCC') {
            cat("Basic specification: H_t = D_t R_t D_t", "\n")
            cat("               diag(D_t) = sqrt(h_ii,t) = c_h + a_h*y^2_[t-1] + b_h*h_[ii,t-1]", "\n")
            cat("               R_t = Q^[-1]_t Q_t Q^[-1]_t = ( 1 - a_q - b_q)S + a_q(u_[t-1]u'_[t-1]) + b_q(Q_[t-1])", "\n")
        } else {
            if( object$param == 'BEKK' | object$param == 'pdBEKK') {
                cat("Basic specification: H_t = C + A'[y_(t-1)*y'_(t-1)]A + B'H_(t-1)B", "\n")
            }
        }
    }
    cat("\n")
    cat("Distribution: ", object$distribution, "\n")
    cat("---\n")
    cat("Iterations: ", object$iter, "\n")
    cat("Chains: ", object$chains, "\n")
    cat("Date: ", object$date, "\n")
    cat("Elapsed time (min):",
    round((max(object$elapsed_time[,1]) + max(object$elapsed_time[,2]))/60, 2) , "\n\n")

    ## Shared objects
    nt <- object$nt
    short_names <- abbreviate(object$TS_names, minlength = 2)
    ## extract and printonly off-diagonal elements
    cormat_index <- matrix(1:(nt*nt), nrow = nt)
    corr_only <- cormat_index[lower.tri(cormat_index)]
    diag_only <- diag(cormat_index)
    ## obtain all combinations of TS varnames for A and B in BEKK
    full_varnames <- expand.grid( short_names, short_names)
    ## obtain off-diagonal TS varnames
    od_varnames <- full_varnames[corr_only, ]
    P <- object$mgarchP
    Q <- object$mgarchQ

    ## depending on parameteriztion, different parameters will be returned:
    if(object$param == 'CCC') {
        model_summary <- rstan::summary(object$model_fit,
                                       pars = c('c_h', 'a_h', 'b_h', 'R', 'phi0', 'phi', 'theta', 'lp__'),
                                       probs = CrI)$summary[, -2 ] # drop se_mean with [,-2]
        
        garch_h_index  = grep("_h", rownames(model_summary) )
        cond_corr_index = grep("R", rownames(model_summary) )
        if( object$meanstructure == 0 ) {
            arma_index = grep("phi0", rownames(model_summary))
        } else {
            arma_index = grep("^phi|^theta", rownames(model_summary) )
        }
          
        ##############################
        ## GARCH parameters of CCC  ##
        ##############################
        cat(paste0(paste0("GARCH(", P, ",", Q, ')')), "estimates for conditional variance:", "\n\n") 
        
        rn = rownames(model_summary[garch_h_index,])
        for ( i in 1:nt ) {
           # replace = grep(paste("\\[", "\\]", sep=as.character(i)), rn)
            replace = grep(paste0( as.character(i), "\\]"), rn) 
        rn[replace] = gsub(paste0( as.character(i), "\\]" ), paste0(short_names[i]), rn)[replace]
        }
        rn = gsub("\\[", "_", rn)

        ## Save into new object to change rownames
        garch_out = model_summary[garch_h_index,]
        ## Assign new rownames
        rownames(garch_out) = rn
        ## print garch parameters
        print(round( garch_out, digits = digits) )
        cat("\n\n")

        ###############################
        ## Constant Correlation R ##
        ###############################
        cat("Constant correlation (R) coefficient(s):", "\n\n")
        
        corr_out = model_summary[cond_corr_index[corr_only],]
        if ( nt == 2 ) {
            tmp = matrix( corr_out, nrow = 1 )
            rownames(tmp) = paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
            colnames(tmp) = names(corr_out)
            corr_out = tmp } else {
                               rownames(corr_out) =
                                   paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
                           }
        print(round( corr_out, digits = digits) )
        cat("\n\n")
        
        ###############################
        ## Location parameters       ##
        ###############################
        cat("ARMA(1,1) estimates on the location:", "\n\n")
        print(round( model_summary[arma_index,], digits = digits) )
        cat("\n\n") } else {
                        
    #########
    ## DCC ##
    #########
    if(object$param == 'DCC') {
        model_summary <- rstan::summary(object$model_fit,
                                       pars = c('a_q', 'b_q', 'c_h_var', 'a_h', 'b_h', 'S', 'phi0', 'phi', 'theta','beta', 'lp__'),
                                       probs = CrI)$summary[, -2]
            
        garch_h_index <- grep("_h", rownames(model_summary) )
        garch_q_index  = grep("_q", rownames(model_summary) )
        S_index = grep("S", rownames(model_summary) )
        if( object$meanstructure == 0 ) {
            arma_index = grep("phi0", rownames(model_summary))
        } else {
            arma_index = grep("^phi|^theta", rownames(model_summary) )
        }

        ## ###########################
        ## GARCH parameters of DCC  ##
        ## ###########################
        cat(paste0(paste0("GARCH(", P, ",", Q, ')')), "estimates for conditional variance on D:", "\n\n") 

        rn <- rownames(model_summary[garch_h_index,])
        for ( i in 1:nt ) {
            replace <- grep(paste0( as.character(i), "\\]"), rn) 
        rn[replace] <- gsub(paste0( as.character(i), "\\]" ), paste0(short_names[i]), rn)[replace]
        }
        rn <- gsub("\\[", "_", rn)

        ## Save into new object to change rownames
        garch_h_out <- model_summary[garch_h_index,]
        ## Assign new rownames
        rownames(garch_h_out) <- rn
        ## print garch parameters
        print(round( garch_h_out, digits = digits) )
        cat("\n\n")

        ##########################
        ## GARCH estimates on Q ##
        ##########################
        cat("GARCH(1,1) estimates for conditional variance on Q:", "\n\n")
        rn = rownames(model_summary[garch_q_index,])
        for ( i in 1:nt ) {
            replace = grep(paste("\\[", "\\]", sep=as.character(i)), rn)
            rn[replace] = gsub(paste("\\[", "\\]", sep=as.character(i)), paste0("_", short_names[i]), rn)[replace]
        }

        ## Save into new object to change rownames
        garch_q_out = model_summary[garch_q_index,]
        ## Assign new rownames
        rownames(garch_q_out) = rn
        ## print garch parameters
        print(round( garch_q_out, digits = digits) )
        cat("\n\n")

        cat("Unconditional correlation 'S' in Q:", "\n\n")
        S_out <- model_summary[S_index[corr_only],]
        if ( nt == 2 ) {
            tmp <- matrix( S_out, nrow = 1 )
            rownames(tmp) <- paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
            colnames(tmp) <- names(S_out)
            S_out <- tmp } else {
                               rownames(S_out) <- 
                                   paste( paste0("S_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
                           }
        print(round( S_out, digits = digits) )
        cat("\n\n")

        ###############################
        ## Location parameters       ##
        ###############################       
        cat("ARMA(1,1) estimates on the location:", "\n\n")
        print(round( model_summary[arma_index,], digits = digits) )
        cat("\n\n") } else {

    #####################              
    ## BEKK and pdBEKK ##
    #####################
    if(object$param == 'BEKK' | object$param == 'pdBEKK') {
        model_summary <- rstan::summary(object$model_fit,
                                        pars = c('beta0', 'C_var', 'A', 'B', 'C_R', 'phi0', 'phi', 'theta', 'lp__'),
                                        probs = CrI)$summary[, -2 ] 

        garch_C_index <- grep("beta0", rownames(model_summary) )
        garch_Cv_index <- grep("C_var", rownames(model_summary) )
        garch_R_index <- grep("C_R", rownames(model_summary) )                            
        garch_A_index <- grep("A", rownames(model_summary) )
        garch_B_index <- grep("B", rownames(model_summary) )                        
        if( object$meanstructure == 0 ) {
            arma_index <- grep("phi0", rownames(model_summary))
        } else {
            arma_index <- grep("^phi|^theta", rownames(model_summary) )
        }

        #######################
        ## Garch parameters  ##
        #######################
        
        cat("---\n\n")

        cat("Constant correlation, R (diag[C]*R*diag[C]):", "\n\n")
        R_out <- model_summary[garch_R_index[corr_only],]
        if ( nt == 2 ) {
            tmp <- matrix( R_out, nrow = 1 )
            rownames(tmp) <- paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ),
                                  substr(od_varnames[ ,2], 1, 2) , sep = '-')
            colnames(tmp) <- names(R_out)
            R_out = tmp } else {
                               rownames(R_out) <- 
                                   paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ),
                                         substr(od_varnames[ ,2], 1, 2) , sep = '-')
                           }
        print(round( R_out, digits = digits) )
        cat("\n\n")

        cat("Constant variances (diag[C]):", "\n\n")
        C_out = model_summary[garch_Cv_index,]
        if ( nt == 2 ) {
            tmp = matrix( C_out, nrow = 2 )
            rownames(tmp) = paste0("var_", short_names )
            colnames(tmp) = colnames(C_out)
            C_out = tmp } else {
                            rownames(C_out) = paste0("var_", short_names )
                           }
        print(round( C_out, digits = digits) )
        cat("\n\n")

        #############
        ## A and B ##
        #############

        #######
        ## A ##
        #######
        cat(paste0(paste0("MGARCH(", P, ",", Q, ')')), "estimates for A:", "\n\n")

        a = list()
        A_out = model_summary[garch_A_index,]
            if (Q > 1) {
                for( q in 1:Q ){
                    a[[q]] = paste( paste0( paste0("A_", q, "_"), full_varnames[ ,1] ), full_varnames[ ,2] , sep = '-')
                }
                rownames(A_out) = unlist( a )
            } else rownames(A_out) = paste( paste0("A_", full_varnames[ ,1] ), full_varnames[ ,2] , sep = '-')
        print(round( A_out, digits = digits) )
        cat("\n\n")
        
        #######
        ## B ##
        #######
        cat(paste0(paste0("MGARCH(", P, ",", Q, ')')), "estimates for B:", "\n\n")

        b <- list()
        B_out <- model_summary[garch_B_index,]
        ##if ( nt == 2 ) {
        if (P > 1) {
            for( p in 1:P ){
                b[[p]] <- paste( paste0( paste0("B_", p, "_"), full_varnames[ ,1] ), full_varnames[ ,2] , sep = '-')
            }
            rownames( B_out ) <- unlist( b )
        } else rownames( B_out ) <- paste( paste0("B_", full_varnames[ ,1] ), full_varnames[ ,2] , sep = '-')
        print(round( B_out, digits = digits) )
        cat("\n\n")
        
        #####################
        ## ARMA parameters ##
        #####################
        cat(paste0(paste(paste0("ARMA(", 1), 1, sep = ','), ')'), "estimates on the location:", "\n\n")
        print(round( model_summary[arma_index,], digits = digits) )
        cat("\n\n")
        }
    }
   }

    ## ######################
    ## Beta: Predictor on H
    ## #####################
    if( sum(object$xC) != 0) {
        if(object$param == 'BEKK' | object$param == 'pdBEKK') {
        cat("Exogenous predictor (beta1 on log scale: C = sRs with s = exp( x*beta ):", "\n\n")
        beta <- rstan::summary(object$model_fit, pars = c('beta1'), probs = CrI)$summary[, -2]
        rownames(beta) <- paste0( "beta1_", short_names )
        print(round(beta, digits = digits ) )
        cat("\n\n" )
        } else {
            cat("Exogenous predictor (beta1 on log scale: c = exp( beta_0 + beta_1*x ):", "\n\n")
            beta0 <- rstan::summary(object$model_fit, pars = c('c_h'), probs = CrI)$summary[, -2]
            rownames(beta0) <- paste0( "beta0_", short_names )
            beta <- rstan::summary(object$model_fit, pars = c('beta'), probs = CrI)$summary[, -2]
            rownames(beta) <- paste0( "beta_", short_names )
            print(round(rbind(beta0,beta), digits = digits ) )
            cat("\n\n" )
        }
    }
        
    nu <- rstan::summary(object$model_fit, pars = c('nu'))$summary[,'mean']
    Lnu <- round( rstan::summary(object$model_fit, pars = c('nu'), probs = CrI)$summary[,4], 2)
    Unu <- round( rstan::summary(object$model_fit, pars = c('nu'), probs = CrI)$summary[,5], 2)

    
    if( object$num_dist == 1) {
        cat("Df constant student_t: nu =", round( nu, digits = 2), "\n")
        cat("                nu 95%CrI [", Lnu, ";", Unu, "]","\n\n")
    }

    cat("Log density posterior estimate:", "\n\n")
    print(round( model_summary[grep("lp__", rownames(model_summary) ),], digits = digits) )
}
