##' @title Forecasting mean and variance 
##' @param object Fitted bmgarch object
##' @param ahead Periods to be forecasted ahead
##' @return Forecasted mean and variance   
##' @author philippe
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 geom_errorbar
##' @export
forecast <- function(object, ahead =  1, seed = NA) {

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
             frcst_DCC <- rstan::gqs(stanmodels$forecastDCC,
                                 draws =  as.matrix(object$model_fit),
                                 data = standat,
                                 seed =  seed)
            
            return(frcst_DCC)

        ## ##########
        ## END DCC ##
        ## ##########
        
    } else {
        if ( object$param == 'CCC') {
            ## ###############
            ## CCC Forecast ##
            ## ###############
             frcst_CCC <- rstan::gqs(stanmodels$forecastCCC,
                                 draws =  as.matrix(object$model_fit),
                                 data = standat,
                                 seed =  seed)
            
            return(frcst_CCC)            
            ## CCC specific parameters
        
        ## ##########
        ## END CCC ##
        ## ###########
        
    } else {
        if( object$param == 'BEKK' ) {
            ## ###################
            ## BEKK Forecast ##
            ## ###################
             frcst_BEKK <- rstan::gqs(stanmodels$forecastBEKK,
                                 draws =  as.matrix(object$model_fit),
                                 data = standat,
                                 seed =  seed)
            
            return(frcst_BEKK)

            ## BEKK Specific parameters
            ## ###########
            ## END BEKK ##
            ## ###########
        }
    }
    }
    
    ## ###########################
    ## Plot and print forecasts ##
    ## ###########################
            
    ## ## #################################################
    ## ## compute means and CrI's from predicted values ##
    ## pred_means = array( NA, dim = c( ahead, object$nt ) )
    ## pred_ci = array( NA, dim = c( ahead, 2, object$nt ) )
    
    ## for ( p in 1:ahead){
    ##     pred_means[p, ] = colMeans( rts_p[, p, ] )
    ##     pred_ci[p, ,] =  apply( rts_p[, p, ], 2, FUN = bmgarch:::.qtile )
    ## }

    ## print( list( pred_means, pred_ci ) )
    

    ## ## ####################
    ## ## Plotting options ##
    ## ## ####################

    ## ## use plots from bgarch::plot() function and add prediction

    ##                                     #retro_plot = plot(object, type = 'means')

    ## rts_data = object$RTS_full

    ## ## TODO go into plot and define df as array, so that it saves the df for each TS?
    ##                                     #df = retro_plot$past_data

    
    ## cis = array(NA, dim = c( nrow(rts_data)+ahead, 2*object$nt ) )
    ## cis_names =  paste0(rep( object$TS_names, each = 2 ),
    ##                     rep (c('_lower', '_upper'), object$nt ))

    ## colnames(cis) = cis_names
    
    ## ## write CI limits into cis
    
    ## index = seq(from = nrow(rts_data)+1, nrow(rts_data)+ahead, by = 1)
    ## cis[ index,] = c(pred_ci) 

    ## ## append predictions
    ## df = data.frame ( rbind( as.matrix(rts_data), pred_means) )
    ## df$prediction = c( rep( 0, object$TS_length), rep( 1, ahead ) )
    ## df$prediction = factor(df$prediction, levels = c(0, 1), labels = c('Past', 'Forecast') )

    ## df$time = 1:nrow(df)

    ## ## add CI's to df
    ## df_p = data.frame(df, cis)
    
    
    ## colNames <- rlang::syms(names(df_p))
    
    
    ## for( i in 1:object$nt){
    ##     ## find relevant positions for current variable in df_p
    ##     var_pos = grep(object$TS_names[i], colNames )
        
    ##     base = ggplot( data = df_p, aes( x = time , y = !!colNames[[ var_pos[1] ]], color = prediction )) + geom_line()
    ##     base = base + geom_errorbar(aes( ymin = !!colNames[[var_pos[2]]], ymax = !!colNames[[var_pos[3]]]), width = .1) +
    ##         geom_point(data = df_p[ ( object$TS_length+1 ):( object$TS_length + ahead ), ], color = 'blue')
    ##     plot(base)
    ##     if ( i == 1){
    ##         devAskNewPage( ask = TRUE )
    ##     }
    ## }

    
######################
    ## ... Forecast     ##
######################
}
