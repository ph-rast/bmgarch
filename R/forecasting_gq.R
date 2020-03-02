##' %@title Forecasting mean and variance 
##' @param object Fitted bmgarch object
##' @param ahead Periods to be forecasted ahead
##' @param xC_p Predictor for H
##' @param type Plot conditional means, variances (default) or correlations. Takes arguments "mean", "var", "cor"
##' @param CrI Lower and upper bound of credible interval. Default is "c( 0.025, .975)".
##' Possible values are .025, .05, .10, .50, .90, .95, and .975
##' @param plot Should forecast be plotted. Logical argument, defaults to TRUE
##' @param last_t For plotting: Only show last t observations in plot
##' @return Forecasted conditional means, variances and correlations. Forecasts are shown in plots (unless plot = FALSE).
##' Conditional variance and correlation forecasts are appended to the corresponding plots with the estimates over the given period.
##' Mean forecasts are appended to the observed values if no meanstructure is estimated (default), otherwise they will be appended
##' to the esimated means. 
##' @author philippe
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 geom_errorbar
##' @importFrom ggplot2 aes_string
##' @export
forecast <- function(object,
                     ahead =  1,
                     type =  "var",
                     xC_p =  NULL,
                     CrI =  c(.025, .975),
                     seed = NA,
                     plot =  TRUE,
                     last_t =  100) {

    ## Define a 0 array as the stan models need some non-null matrix
    if( is.null( xC_p ) ) xC_p <- array( rep(0, object$nt),  dim = c(ahead,  object$nt ) )

    ## obtain data to be passed to gqs
    standat <-  list(T = object$TS_length,
                     nt = object$nt,
                     rts = cbind(object$RTS_full),
                     xC = object$xC,
                     Q =  object$mgarchQ,
                     P =  object$mgarchP,
                     ahead =  ahead, 
                     meanstructure =  object$meanstructure,
                     distribution =  object$num_dist,
                     xC_p =  xC_p)

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
            if( object$param == 'BEKK' | object$param == 'pdBEKK' ) {
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
    
    ## We only need rts_p and H_p and R_p
    ## note that rts_p and H_p include past observations (given by P or Q)
    ## and forecasted elements given (given by ahead)

    ## Global variables
    ## Separate estimate (backcast) from forecast, given dimension of P, Q, and ahead
    backcast <-  max(object$mgarchP, object$mgarchQ)
    nt <- object$nt

    ## ###########################
    ## Plot and print forecasts ##
    ## ###########################

    ## init plt
    plt <- list()
    
    if( type == "mean" ) {
        ## Two versions: When no meanstructure present (default), predicted means
        ## will be plotted with the observed values.
        ## If meanstructure is present ( != 0 ), forecast is added to plot with
        ## estimated means and CrI's

        ## collect rts elements
        rts_pred_mn <-  rstan::summary(forecasted, pars =  "rts_p",
                                       probs =  .5 )$summary[, 4]        
        
        rts_pred_Lcri <-  rstan::summary(forecasted, pars =  "rts_p",
                                         probs =  CrI[1] )$summary[, 4]
        rts_pred_Ucri <-  rstan::summary(forecasted, pars =  "rts_p",
                                         probs =  CrI[2] )$summary[, 4]
        
        ## Dummy code periods: observed vs forecast
        rts_period <- rep( seq_len( length( rts_pred_mn )/nt ), each = nt )
        rts_frcst <- ifelse( rts_period < max(rts_period) - (ahead -1 ), 0, 1)
        
        rts_mn <- cbind(rts_pred_mn, rts_period, rts_frcst)
        rts_Lcri <- cbind(rts_pred_Lcri, rts_period, rts_frcst)
        rts_Ucri <- cbind(rts_pred_Ucri, rts_period, rts_frcst)

        ## ##########
        ## Means:  ##
        ## ##########

        ## select only forecasted
        rts_mn_forecasted <- rts_mn[rts_mn[, 3] != 0, ]
        rts_Lcri_forecasted <- rts_Lcri[ rts_Lcri[, 3] != 0, ]
        rts_Ucri_forecasted <- rts_Ucri[ rts_Ucri[, 3] != 0, ]

        ## Write into new object
        forecasted_data_rts <- cbind(rts_mn_forecasted[, 1],
                                     rts_Lcri_forecasted[, 1],
                                     rts_Ucri_forecasted[, 1])

        colnames(forecasted_data_rts) <- c('rts_forecasted',
                                           paste0('CrI_', CrI[1]),
                                           paste0('CrI_', CrI[2]))

        forecasted_rts <- list()

        for(i in 1:nt ) {
            forecasted_rts[[i]] <- forecasted_data_rts[seq(1, ahead * nt , by = nt)+(i-1), ]            
            names(forecasted_rts[[i]])[1] <- paste(colnames(standat$rts)[i],
                                                      names(forecasted_rts[[i]])[1], sep = "_")
        }

        
        ## Plots:
        if( plot ) {
            if ( object$meanstructure != 0 ) {                
                ## Use conditional means plot of bmgarch::plot
                retro <- bmgarch::plot.bmgarch(object, type = 'means', askNewPage = FALSE, CrI = CrI )
                
                for( i in seq_len(nt) ) {
                    df <- array( NA, dim = c(standat$T, 3) )
                    df <- data.frame( rbind(df, forecasted_rts[[i]]) )
                    df$period <- seq_len( dim(df )[1] )
                    names(df )[1] <- object$TS_names[i]
                    
                    plt[[i]] <- retro$retro_plot[[i]] +
                        geom_line(data = df, aes_string( x = "period", y = names(df)[1] ) ) +                
                        geom_errorbar(data = df,  aes_string(x = "period", ymin = names(df)[2], ymax = names(df)[3]),
                                      inherit.aes = FALSE , alpha =  .3)+
                        geom_point(data = df, mapping = aes_string(x = "period", y = names(df)[1] ),
                                   inherit.aes = FALSE)+
                        coord_cartesian( xlim = c(dim(df )[1]-last_t,  dim(df )[1]) )
                    
                    plot(plt[[i]])
                    if ( i == 1) {
                        devAskNewPage( ask = TRUE )
                    }
                }
            } else if ( object$meanstructure == 0 ) {
                ## Plot on observed data
                rts_data <- object$RTS_full
                
                for( i in seq_len(nt) ) {
                    df <- array( NA, dim = c(standat$T, 3) )
                    df[, 1] <- rts_data[, i]
                    df <- data.frame( rbind(df, forecasted_rts[[i]]) )
                    df$Type <- factor( c( rep("Observed", standat$T),
                                         rep("Forecasted",  ahead)) )
                    
                    df$period <- seq_len( dim(df )[1] )
                    names(df )[1] <- object$TS_names[i]
                    
                    plt[[i]] <- ggplot(data = df, aes_string(x = "period",
                                                             y = names(df)[1],
                                                             color =  "Type") ) +
                        geom_line() +
                        geom_errorbar( aes_string(ymin = names(df)[2], ymax = names(df)[3]),
                                      alpha =  .3)+
                        geom_point(data = df[df$Type == "Forecasted", ], mapping = aes_string(x = "period", y = names(df)[1] ))+
                        coord_cartesian( xlim = c(dim(df )[1]-last_t,  dim(df )[1]) )
                    
                    plot(plt[[i]])
                    if ( i == 1) {
                        devAskNewPage( ask = TRUE )
                    }
                }  
            }
            
        }
        
        return( forecasted_rts )
        
    } else if(type == "var" ) {
        ## ############
        ## Variances ##
        ## ############
        ## Print variances
        ## obtain predicted variances:

        ## collect H elements
        H_pred_mn <-  rstan::summary(forecasted, pars =  "H_p",
                                     probs =  c(.5) )$summary[, 4]
        H_pred_Lcri <-  rstan::summary(forecasted, pars =  "H_p",
                                       probs =  CrI[1] )$summary[, 4]
        H_pred_Ucri <-  rstan::summary(forecasted, pars =  "H_p",
                                       probs =  CrI[2] )$summary[, 4]
        ## Dummy code periods: observed vs forecast
        H_period <- rep( seq_len( length( H_pred_mn )/nt^2 ), each = nt^2 )
        H_frcst <- ifelse( H_period < max(H_period) - (ahead -1 ), 0, 1)
        
        H_mn <- cbind(H_pred_mn, H_period, H_frcst)
        H_Lcri <- cbind(H_pred_Lcri, H_period, H_frcst)
        H_Ucri <- cbind(H_pred_Ucri, H_period, H_frcst)
        
        ## select only forecasted
        H_mn_forecasted <- H_mn[H_mn[, 3] != 0, ]
        H_Lcri_forecasted <- H_Lcri[ H_Lcri[, 3] != 0, ]
        H_Ucri_forecasted <- H_Ucri[ H_Ucri[, 3] != 0, ]
        
        ## find variances
        
        splitted <- strsplit( rownames(H_mn_forecasted), split = ",|]")
        row <- as.numeric( t( sapply( splitted, "[", 2:3) )[,1])
        col <- as.numeric( t( sapply( splitted, "[", 2:3) )[,2])
        
        ## obtain position of variances:
        var_pos <- which( row == col )

        ## Write into new object
        forecasted_data_H <- cbind(H_mn_forecasted[var_pos, 1],
                                   H_Lcri_forecasted[var_pos, 1],
                                   H_Ucri_forecasted[var_pos, 1])

        colnames(forecasted_data_H) <- c('H_forecasted',
                                         paste0('CrI_', CrI[1]),
                                         paste0('CrI_', CrI[2]))

        forecasted_H <- list()
        
        for(i in 1:nt ) {
            forecasted_H[[i]] <- forecasted_data_H[seq(1, ahead * nt , by = nt)+(i-1), ]            
            names(forecasted_H[[i]])[1] <- paste(colnames(standat$rts)[i],
                                                    names(forecasted_H[[i]])[1], sep = "_")
        }
        
        ## Plots
        if( plot ) {
            ## Use conditional variance plot of bmgarch::plot
            retro <- bmgarch::plot.bmgarch(object, type = 'cvar', askNewPage = FALSE, CrI = CrI)

            ## Loop through nt's
            for(i in seq_len(nt) ) {
                df <- array( NA, dim = c(standat$T, 3) )
                df <- data.frame( rbind(df, forecasted_H[[i]]) )
                df$period <- seq_len( dim(df )[1] )
                df$Type <- factor( c( rep("Estimated", standat$T),
                                     rep("Forecasted",  ahead)) )

                
                plt[[i]] <- retro$retro_plot[[i]] +
                    geom_line(data = df, aes_string( x = "period", y = names(df)[1],  color =  "Type" ) ) +
                    geom_errorbar(data = df,  aes_string(x = "period", ymin = names(df)[2], ymax = names(df)[3]),
                                  inherit.aes = FALSE , alpha =  .3)+
                    geom_point(data = df, mapping = aes_string(x = "period",
                                                               y = names(df)[1],
                                                               color =  "Type" ), inherit.aes = FALSE)+
                    coord_cartesian( xlim = c(dim(df )[1]-last_t,  dim(df )[1]) )
                
                plot(plt[[i]])
                if ( i == 1 ) {
                    devAskNewPage( ask = TRUE)
                }
            }
        }
        ## print
        return( forecasted_H )
    } else if(type == "cor" ) {
        ## ############
        ## Correlations  ##
        ## ############

        ## collect R elements
        R_pred_mn <-  rstan::summary(forecasted, pars =  "R_p",
                                     probs =  c(.5) )$summary[, 4]
        R_pred_Lcri <-  rstan::summary(forecasted, pars =  "R_p",
                                       probs =  CrI[1] )$summary[, 4]
        R_pred_Ucri <-  rstan::summary(forecasted, pars =  "R_p",
                                       probs =  CrI[2] )$summary[, 4]

        ## Dummy code periods: observed vs forecast
        R_period <- rep( seq_len( length( R_pred_mn )/nt^2 ), each = nt^2 )
        R_frcst <- ifelse( R_period < max(R_period) - (ahead -1 ), 0, 1)
        
        R_mn <- cbind(R_pred_mn, R_period, R_frcst)
        R_Lcri <- cbind(R_pred_Lcri, R_period, R_frcst)
        R_Ucri <- cbind(R_pred_Ucri, R_period, R_frcst)

        ## select only forecasted
        R_mn_forecasted <- R_mn[R_mn[, 3] != 0, ]
        R_Lcri_forecasted <- R_Lcri[ R_Lcri[, 3] != 0, ]
        R_Ucri_forecasted <- R_Ucri[ R_Ucri[, 3] != 0, ]
        
        ## find correlations
        splitted <- strsplit( rownames(R_mn_forecasted), split = ",|]")
        row <- as.numeric( t( sapply( splitted, "[", 2:3) )[,1])
        col <- as.numeric( t( sapply( splitted, "[", 2:3) )[,2])

        ## obtain position of correlations:
        cor_pos <- which( row != col  & row < col )
        
        row_col <- cbind(row, col )[cor_pos, ]

        ## Write int new object
        forecasted_data_R <- cbind(R_mn_forecasted[cor_pos, 1],
                                   R_Lcri_forecasted[cor_pos, 1],
                                   R_Ucri_forecasted[cor_pos, 1])

        colnames(forecasted_data_R) <- c('R_frcsted',
                                         paste0('CrI_', CrI[1]),
                                         paste0('CrI_', CrI[2]))

        forecasted_R <- list()
        corrs <- (nt^2 - nt)/2
        for(i in 1:corrs  ) {
            forecasted_R[[i]] <- forecasted_data_R[seq(1, ahead * corrs , by = corrs)+(i-1), ]
            if(corrs == 1 ) label_pos <-  1 else {
            label_pos <- row_col[seq(1, ahead * corrs , by = corrs)+(i-1), ][1, ]}
            names(forecasted_R[[i]])[1] <- paste( paste(abbreviate( colnames(standat$rts)[label_pos[1]]),
                                                        abbreviate( colnames(standat$rts)[label_pos[2]]),
                                                           sep =  "_"),
                                                    names(forecasted_R[[i]])[1], sep = "_")
        }
        
        ## Plots
        if( plot ) {
            ## Use conditional variance plot of bmgarch::plot
            retro <- plot(object, type = 'ccor', askNewPage = FALSE, CrI = CrI)

            ## FOR LOOP HERE with number of correlations (corrs)
            for(i in 1:corrs ) {
                df <- array( NA, dim = c(standat$T, 3) )
                df <- data.frame( rbind(df, forecasted_R[[i]]) )
                df$period <- seq_len( dim(df )[1] ) 
                df$Type <- factor( c( rep("Estimated", standat$T),
                                     rep("Forecasted",  ahead)) )
                
                plt[[i]] <- retro$retro_plot[[i]] +
                    geom_line(data = df, aes_string(x = "period",
                                                    y = names(df)[1],
                                                    color =  "Type") ) +                
                    geom_errorbar(data = df,  aes_string(x = "period",
                                                         ymin = names(df)[2],
                                                         ymax = names(df)[3]),
                                  inherit.aes = FALSE , alpha =  .3)+
                    geom_point(data = df[df$Type == "Forecasted", ],
                               mapping = aes_string(x = "period",
                                                    y = names(df)[1],
                                                    color =  "Type"), inherit.aes = FALSE) +
                    coord_cartesian( xlim = c(dim(df )[1]-last_t,  dim(df )[1]) )
                
                
                plot(plt[[i]])
                if ( i == 1 ) {
                    devAskNewPage( ask = TRUE)
                }
            }
        }
        ## print
        return( forecasted_R )
    }
}
