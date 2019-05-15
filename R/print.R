##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Summary for BMGARCH object
##' @param object 
##' @param CrI 
##' @return Summary object
##' @author philippe
summary.bmgarch = function(object, CrI = c(0.05, 0.95) ){
    cat( "Model:", paste0(object$param, "-MGARCH\n"))
    cat("---\n")
    cat("Iterations: ", object$iter, "\n")
    cat("Chains: ", object$chains, "\n")
    cat("Date: ", object$date, "\n")
    cat("Elapsed time (min):",
    round((max(object$elapsed_time[,1]) + max(object$elapsed_time[,2]))/60, 2) , "\n\n")

    ## Shared objects among models
    nt = object$nt
    short_names = substr(object$TS_names, 1, 2)
     ## extract and print only off-diagonal elements
    cormat_index = matrix(1:(nt*nt), nrow = nt)
    corr_only = cormat_index[lower.tri(cormat_index)]
    ## obtain off-diagonal TS varnames
    od_varnames = expand.grid( short_names, short_names)[corr_only, ]

    
    ## depending on parameteriztion, different parameters will be returned:
    if(object$param == 'CCC') {
        model_summary = rstan::summary(object$model_fit,
                                       pars = c('c_h', 'a_h', 'b_h', 'R', 'b0', 'b1', 'b2', 'b3', 'b4', 'lp__'),
                                       probs = CrI)$summary
        
        garch_h_index  = grep("_h", rownames(model_summary) )
        cond_corr_index = grep("R", rownames(model_summary) )
        arma_index = grep("b[0-9]", rownames(model_summary) )

        ##############################
        ## GARCH parameters of CCC  ##
        ##############################
        cat("GARCH(1,1) estimates for conditional variance:", "\n\n")
        
        rn = rownames(model_summary[garch_h_index,])
        for ( i in 1:nt ) {
            replace = grep(paste("\\[", "\\]", sep=as.character(i)), rn)
            rn[replace] = gsub(paste("\\[", "\\]", sep=as.character(i)), paste0("_", short_names[i]), rn)[replace]
        }

        ## Save into new object to change rownames
        garch_out = model_summary[garch_h_index,]
        ## Assign new rownames
        rownames(garch_out) = rn
        ## print garch parameters
        print(round( garch_out, 2) )
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
        print(round( corr_out, 2) )
        cat("\n\n")
        
        ###############################
        ## Location parameters       ##
        ###############################
        cat("ARMA(1,1) estimates on the location:", "\n\n")
        print(round( model_summary[arma_index,], 2) )
        cat("\n\n") } else {
                        
    #########
    ## DCC ##
    #########
    if(object$param == 'DCC') {
        model_summary = rstan::summary(object$model_fit,
                                       pars = c('a_q', 'b_q', 'c_h', 'a_h', 'b_h', 'S', 'b0', 'b1', 'b2', 'b3', 'b4', 'lp__'),
                                       probs = CrI)$summary
            
        garch_h_index  = grep("_h", rownames(model_summary) )
        garch_q_index  = grep("_q", rownames(model_summary) )
        S_index = grep("S", rownames(model_summary) )
        arma_index = grep("b[0-9]", rownames(model_summary) )

        ##############################
        ## GARCH parameters of DCC  ##
        ##############################
        cat("GARCH(1,1) estimates for conditional variance on D:", "\n\n")

        rn = rownames(model_summary[garch_h_index,])
        for ( i in 1:nt ) {
            replace = grep(paste("\\[", "\\]", sep=as.character(i)), rn)
            rn[replace] = gsub(paste("\\[", "\\]", sep=as.character(i)), paste0("_", short_names[i]), rn)[replace]
        }

        ## Save into new object to change rownames
        garch_h_out = model_summary[garch_h_index,]
        ## Assign new rownames
        rownames(garch_h_out) = rn
        ## print garch parameters
        print(round( garch_h_out, 2) )
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
        print(round( garch_q_out, 2) )
        cat("\n\n")

        cat("Unconditional correlation 'S' in Q:", "\n\n")
        S_out = model_summary[S_index[corr_only],]
        if ( nt == 2 ) {
            tmp = matrix( S_out, nrow = 1 )
            rownames(tmp) = paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
            colnames(tmp) = names(S_out)
            S_out = tmp } else {
                               rownames(S_out) =
                                   paste( paste0("S_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
                           }
        print(round( S_out, 2) )
        cat("\n\n")

        ###############################
        ## Location parameters       ##
        ###############################       
        cat("ARMA(1,1) estimates on the location:", "\n\n")
        print(round( model_summary[arma_index,], 2) )
        cat("\n\n") } else {

    ##########              
    ## BEKK ##
    ##########
    if(object$param == 'BEKK') {
        model_bekk = rstan::summary(object$model_fit,
                                    pars = c('Cnst', 'A', 'B', 'corC', 'b0', 'b1', 'b2', 'b3', 'b4', 'lp__'),
                                    probs = CrI)$summary }
                            
        garch_h_index  = grep("_h", rownames(model_summary) )
        cond_corr_index = grep("R", rownames(model_summary) )
        arma_index = grep("b[0-9]", rownames(model_summary) )
                            
        cat("GARCH(1,1) estimates for conditional variance:", "\n\n")
        print(round( model_summary[garch_h_index,], 2) )
        cat("\n\n")
                            
        cat("Conditional correlation(s):", "\n\n")
        print(round( model_summary[cond_corr_index,], 2) )
        cat("\n\n")
                           
        cat("ARMA(1,1) estimates on the location:", "\n\n")
        print(round( model_summary[arma_index,], 2) )
        cat("\n\n")
    }}
    
    cat("Log density posterior estimate:", "\n\n")
    print(round( model_summary[grep("lp__", rownames(model_summary) ),], 2) )
}
