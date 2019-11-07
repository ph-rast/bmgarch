##' @title Forecasting mean and variance 
##' @param object Fitted bmgarch object
##' @param ahead Periods to be forecasted ahead
##' @param type Plot conditional means (default) or (co)-variances. Takes arguments "means", "variances", "covariances"
##' @param CrI Lower and upper boundary of credibility interval. Default is "c( 0.5, .95)"
##' @return Forecasted mean and variance   
##' @author philippe
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 geom_errorbar
##' @importFrom ggplot2 aes_string
##' @export
forecast <- function(object, ahead =  1, type =  "means", CrI =  c(.05,  .95), seed = NA) {

    ## model_cppname <- "forecastGQ.stan"
    ## stanfit <- rstan::stanc("./src/stan_files/forecastGQ.stan",
    ##                         allow_undefined = TRUE, 
    ##                         obfuscate_model_name = FALSE)
    
    ## stanfit$model_cpp <- list(model_cppname = stanfit$model_name, 
    ##                           model_cppcode = stanfit$cppcode)

    ## return(do.call(methods::new, args = c(stanfit[-(1:3)], Class = "stanmodel", 
    ##              mk_cppmodule = function(x) get(paste0("model_", model_cppname)))))
    
    ## obtain data to be passed to gqs
    standat <-  list(T = object$TS_length,
                     nt = object$nt,
                     rts =  cbind(object$RTS_full),
                     xH = object$xH,
                     Q =  object$mgarchQ,
                     P =  object$mgarchP,
                     ahead =  ahead, 
                     meanstructure =  object$meanstructure,
                     distribution =  object$num_dist)

    ## Call forecasting functions
    if( object$param == 'DCC') {
        ## ################
        ## DCC Forecast ##
        ## ################
        forecasted <- rstan::gqs(stanmodels$forecastDCC,
                                 draws =  as.matrix(object$model_fit),
                                 data = standat,
                                 seed =  seed)
        ## ##########
        ## END DCC ##
        ## ##########
        
    } else {
        if ( object$param == 'CCC') {
            ## ###############
            ## CCC Forecast ##
            ## ###############
             forecasted <- rstan::gqs(stanmodels$forecastCCC,
                                 draws =  as.matrix(object$model_fit),
                                 data = standat,
                                 seed =  seed)
        
        ## ##########
        ## END CCC ##
        ## ###########
        
    } else {
        if( object$param == 'BEKK' ) {
            ## ###################
            ## BEKK Forecast ##
            ## ###################
             forecasted <- rstan::gqs(stanmodels$forecastBEKK,
                                 draws =  as.matrix(object$model_fit),
                                 data = standat,
                                 seed =  seed)
            ## ###########
            ## END BEKK ##
            ## ###########
        }
    }
    }
    
    ## We only need rts_p and H_p
    ## not that rts_p and H_p include past observations (given by P or Q)
    ## and forecasted elements given (given by ahead)

    ## colect rts and H elements
    rts_pred_mn <-  rstan::summary(forecasted, pars =  "rts_p",
                                   probs =  CrI )$summary[, 1]
    rts_pred_Lcri <-  rstan::summary(forecasted, pars =  "rts_p",
                                probs =  CrI )$summary[, 4]
    rts_pred_Ucri <-  rstan::summary(forecasted, pars =  "rts_p",
                                probs =  CrI )$summary[, 5]   
    
    H_pred_mn <-  rstan::summary(forecasted, pars =  "H_p",
                                probs =  CrI )$summary[, 1]
    H_pred_Lcri <-  rstan::summary(forecasted, pars =  "H_p",
                                probs =  CrI )$summary[, 4]
    H_pred_Ucri <-  rstan::summary(forecasted, pars =  "H_p",
                                probs =  CrI )$summary[, 5]
    
    ## Separate estimate (backcast) from forecast, given dimension of P, Q, and ahead
    backcast <-  max(object$mgarchP, object$mgarchQ)
    
    nt <- object$nt

    ## Dummy code periods: observed vs forecast
    rts_period <- rep( seq_len( length( rts_pred_mn )/nt ), each = nt )
    rts_frcst <- ifelse( rts_period < max(rts_period) - (ahead -1 ), 0, 1)
    
    rts_mn <- cbind(rts_pred_mn, rts_period, rts_frcst)
    rts_Lcri <- cbind(rts_pred_Lcri, rts_period, rts_frcst)
    rts_Ucri <- cbind(rts_pred_Ucri, rts_period, rts_frcst)

    ## Dummy code periods: observed vs forecast
    H_period <- rep( seq_len( length( H_pred_mn )/nt^2 ), each = nt^2 )
    H_frcst <- ifelse( H_period < max(H_period) - (ahead -1 ), 0, 1)
    
    H_mn <- cbind(H_pred_mn, H_period, H_frcst)
    H_Lcri <- cbind(H_pred_Lcri, H_period, H_frcst)
    H_Ucri <- cbind(H_pred_Ucri, H_period, H_frcst)

    ## ###########################
    ## Plot and print forecasts ##
    ## ###########################

    ## ##########
    ## Means:  ##
    ## ##########
    
    ## Get observed data
    rts_data <- object$RTS_full

    ## Get predicted data for means
    rts_ahead_means <- tail(matrix(rts_mn[, "rts_pred_mn"], ncol = nt, byrow = TRUE),
                            n = ahead)
    colnames(rts_ahead_means) <- colnames(rts_data )

    ## Get predicted CrI's
    rts_ahead_Lcri <- tail(matrix(rts_Lcri[, "rts_pred_Lcri"], ncol = nt, byrow = TRUE),
                            n = ahead)
    colnames(rts_ahead_Lcri) <- paste0( paste0("CrI_", CrI[1]), colnames(rts_data ))

    rts_ahead_Ucri <- tail(matrix(rts_Ucri[, "rts_pred_Ucri"], ncol = nt, byrow = TRUE),
                            n = ahead)
    colnames(rts_ahead_Ucri) <- paste0( paste0("CrI_", CrI[2]), colnames(rts_data ))
    
    ## combine observed with predicted
    ## create index for plotting
    index <- c(rep("observed", nrow(rts_data)), rep("forecast",  ahead))
    
    forecasted_data <- cbind( rts_ahead_means,  rts_ahead_Lcri, rts_ahead_Ucri )

    ## augment observed data matrix
    placeholder <- matrix(0, ncol = ncol(forecasted_data )-nt, nrow = nrow(rts_data) )

    observed_data <- cbind(rts_data,  placeholder )
    colnames(observed_data) <- colnames(forecasted_data )

    ## Check if xts format - if so, covert to matrix
    if(is.xts(observed_data )) observed_data <- data.frame(observed_data )
    
    rts_combined <- data.frame( rbind(observed_data,
                                      forecasted_data), index)

    ## Provide variable names
    names(rts_combined) <- c( colnames(forecasted_data ), "index" )
    rts_combined$period <- seq_len(nrow(rts_combined))

    if( type == "means" ) {
    ## Print:
    print( forecasted_data )

    ## Plots  
    
    for( i in seq_len(nt) ) {
        base <- ggplot(rts_combined,
                       aes_string(x = "period",
                                  y = colnames(rts_combined )[i] ,
                                  color =  "index") ) +
            geom_line() +
            geom_errorbar(aes_string( ymin = colnames(rts_ahead_Lcri )[i],
                                     ymax = colnames(rts_ahead_Ucri )[i] ), width = .1, alpha = .5) +
            geom_point(data = tail(rts_combined,  n =  ahead) )
        
        plot(base)
        if ( i == 1) {
            devAskNewPage( ask = TRUE )
        }
    }
    } else if(type == "covariances" ) {
        ## Use covariance covariance plot of bmgarch::plot

    }
}
