##' @title Forecasting mean and variance 
##' @param object Fitted bmgarch object
##' @param ahead Periods to be forecasted ahead
##' @return Forecasted mean and variance   
##' @author philippe
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 geom_errorbar
##' @export
forecast_old <- function(object, ahead){

    ## Global objects used in forecasts: Location object (mu) is the same for all parameterizations
    
    phi0 = rstan::extract(object$model_fit)[['phi0']]
    phi = rstan::extract(object$model_fit)[['phi']]
    theta = rstan::extract(object$model_fit)[['theta']]

    Q = object$mgarchQ
    P = object$mgarchP

    ## obtain max Q or P for lag in mu and rts
    max_lag =  max( Q,  P )
    mu = rstan::extract(object$model_fit)[['mu']][, (object$TS_length-1):object$TS_length, ]
    rts = t( tail(object$RTS_full, n =  1) )
 ##    array( object$RTS_last , dim = c(object$nt, 1) ) 

    ## Init mu
    mu_p = array(NA, dim = c( nrow( mu ), ahead, object$nt) )

    ## #####################
    ## Meanstructure
    ## #####################
    
    ## Initialize Moving Average array
    MA = array( NA, dim = c(nrow(theta), object$nt) )

    mu_p[, 1, ] =  phi0 + 
        if( is.null( phi )) 0 else {
                                  t( apply(phi, 1, FUN = function(x){ x %*% rts } ) )
                              } +  
        if(is.null(theta) ) 0 else {
                                  t( sapply(1:nrow(theta), function(i) bmgarch:::.f_MA( MA, theta, mu, rts, i ) ) )
                              }

            ## init rts_p 
            rts_p =  array(NA, dim = c( nrow( rstan::extract(object$model_fit)[['mu']] ), ahead, object$nt) )
            
            ## Call forecasting functions
            if( object$param == 'DCC') {
                ## ################
                ## DCC Forecast ##
                ## ################

                ## DCC specific parameters
                rr = t( apply( sweep( mu, 2, rts ), 1, bmgarch:::.square) )

                c_h = rstan::extract(object$model_fit)[['c_h']]
                a_h = rstan::extract(object$model_fit)[['a_h']]
                b_h = rstan::extract(object$model_fit)[['b_h']]

                ## Only last estimate
                D = rstan::extract(object$model_fit)[['D']][, object$TS_length, ]

                a_q = rstan::extract(object$model_fit)[['a_q']]
                b_q = rstan::extract(object$model_fit)[['b_q']]

                u = rstan::extract(object$model_fit)[['u']][, object$TS_length, ]
                S = rstan::extract(object$model_fit)[['S']]
                Q = rstan::extract(object$model_fit)[['Q']][, object$TS_length, ,]

                ## stan: for(d in 1:nt){
                ## stan: rr_p[1, d] = square( rts[T, d] - mu[T, d] )
                ## stan: D_p[1, d] = sqrt( c_h[d] + a_h[d]*rr_p[1, d] +  b_h[d]*D[T,d] );


                D_p = array(NA, dim = c( nrow(c_h), ahead, object$nt) )
                D_p[,1,] = sqrt( c_h + a_h * rr + b_h * D )

                ## stan: }
                ## stan: Q_p[1,] = (1 - a_q - b_q) * S + a_q * (u[T,] * u[T,]') + b_q * Q[T,];
                Q_p = array(NA, dim = c( nrow(c_h), ahead, object$nt, object$nt) )
                Q_sdi_p = Q_p
                R_p = Q_p
                H_p = Q_p
###rts_p = array(NA, dim = c( nrow(c_h), ahead,  object$nt) )
                for( i in 1:nrow(c_h) ){
                    Q_p[i, 1,,]  =  ( 1 - a_q - b_q )[i] * S[i,,] + a_q[i] * diag(u[i,]^2) +  b_q[i] * Q[i,,]
                    Q_sdi_p[i, 1,,] =    diag( 1/sqrt(diag(Q_p[i,1,,])) )
                    R_p[i, 1 ,,] = Q_sdi_p[i, 1,,] %*% Q_p[i, 1,,] %*% Q_sdi_p[i, 1,,]
                    H_p[i, 1 ,,] = diag( D_p[i, 1, ] ) %*% R_p[i, 1, ,] %*% diag( D_p[i, 1, ] )
                    rts_p[i,1,]  = mvtnorm::rmvnorm( 1, mean = mu_p[i, 1, ], sigma = H_p[i, 1 ,,] )
                }

                ## For more than 1 ahead    
                if(ahead >= 2){
                    for( p in 2:ahead ){
                        mu_p[, p, ] =
                            phi0  +
                            t( sapply(1:nrow( theta ), function(i) bmgarch:::.f_array_x_mat(mat_out = MA, array_obj = phi, mat_obj = rts_p[,p-1,], i) )  ) +
                            t( sapply(1:nrow( theta ), function(i) bmgarch:::.f_array_x_mat( MA, theta, ( rts_p[,p-1,] - mu_p[,p-1,]) , i ) ) )
                        
                        rr = (rts_p[,p-1,] - mu_p[, p-1, ])^2
                        D_p[ , p , ] = sqrt( c_h + a_h * rr + b_h * D_p[,p-1,] )
                        
                        ## Loop through iteartions

                        ## init u and overwrite for each p
                        u_p = array(NA, dim = c( nrow(c_h), object$nt) )

                        for( i in 1:nrow(c_h) ){

                            u_p[i,] = ( 1/D_p[i, p-1, ] ) * rts_p[i, p-1,] - mu_p[i, p-1,]
                            
                            Q_p[i, p,,]  = ( 1 - a_q - b_q )[i] * S[i,,] + a_q[i] * diag(u_p[i,]^2) +  b_q[i] * Q[i,,]
                            
                            Q_sdi_p[i, p,,] = diag( 1/sqrt(diag(Q_p[i,p,,])) )
                            
                            R_p[i, p ,,] = Q_sdi_p[i, p,,] %*% Q_p[i, p,,] %*% Q_sdi_p[i, p,,]
                            
                            H_p[i, p ,,] = diag( D_p[i,p, ] ) %*% R_p[i, p, ,] %*% diag( D_p[i,p, ] )
                            
                            ## todo -> this might be a temporary fix: Check forecast for sanity
                            if( is.infinite( H_p[i, p ,,] )[1,1] ) H_p[i, p ,,] = diag(object$nt)
                            
                            rts_p[i,p,]  = mvtnorm::rmvnorm( 1, mean = mu_p[i, p, ], sigma = H_p[i, p ,,] )
                        }
                    }
                }

                ## ##########
                ## END DCC ##
                ## ##########
                
            } else { if ( object$param == 'CCC') {
                         ## ###############
                         ## CCC Forecast ##
                         ## ###############

                         ## CCC specific parameters
                         rr = t( apply( sweep( mu, 2, rts ), 1, bmgarch:::.square) )

                         c_h = rstan::extract(object$model_fit)[['c_h']]
                         a_h = rstan::extract(object$model_fit)[['a_h']]
                         b_h = rstan::extract(object$model_fit)[['b_h']]

                         R = rstan::extract(object$model_fit)[['R']]
                         
                         ## Only last estimate
                         D = rstan::extract(object$model_fit)[['D']][, object$TS_length, ]

                         D_p = array(NA, dim = c( nrow(c_h), ahead, object$nt) )

                         D_p[,1,] = sqrt( c_h + a_h * rr + b_h * D )
                         
                         ## Init H_p
                         H_p = array(NA, dim = c( nrow(c_h), ahead, object$nt, object$nt) )
                         
                         for( i in 1:nrow(c_h) ){
                             H_p[i, 1 ,,] = diag( D_p[i, 1, ] ) %*% R[i, ,] %*% diag( D_p[i, 1, ] )
                             rts_p[i,1,]  = mvtnorm::rmvnorm( 1, mean = mu_p[i, 1, ], sigma = H_p[i, 1 ,,] )
                         }

                         ## For more than 1 ahead    
                         if(ahead >= 2){
                             for( p in 2:ahead ){
                                 mu_p[, p, ] =
                                     phi0  +
                                     t( sapply(1:nrow( theta ),
                                               function(i) bmgarch:::.f_array_x_mat(mat_out = MA, array_obj = phi, mat_obj = rts_p[,p-1,], i) )  ) +
                                     t( sapply(1:nrow( theta ),
                                               function(i) bmgarch:::.f_array_x_mat( MA, theta, ( rts_p[,p-1,] - mu_p[,p-1,]) , i ) ) )
                                 rr = (rts_p[,p-1,] - mu_p[, p-1, ])^2
                                 D_p[ , p , ] = sqrt( c_h + a_h * rr + b_h * D_p[,p-1,] )
                                 
                                 ## Loop through iterations
                                 for( i in 1:nrow(c_h) ){                   
                                     H_p[i, p ,,] = diag( D_p[i,p, ] ) %*% R[i, ,] %*% diag( D_p[i,p, ] )
                                     
                                     ## todo -> this might be a temporary fix: Check forecast for sanity
                                     if( is.infinite( H_p[i, p ,,] )[1,1] ) H_p[i, p ,,] = diag(object$nt)
                                     
                                     rts_p[i,p,]  = mvtnorm::rmvnorm( 1, mean = mu_p[i, p, ], sigma = H_p[i, p ,,] )
                                 }
                             }
                         }

                         ## ##########
                         ## END CCC ##
#############
                         
                     } else { if( object$param == 'BEKK' ) {
                                  ## ###################
                                  ## BEKK Forecast ##
                                  ## ###################

                                  ## BEKK Specific parameters
                                  
                                  Cnst = rstan::extract(object$model_fit)[['Cnst']]
                                  A = rstan::extract(object$model_fit)[['A']]
                                  B = rstan::extract(object$model_fit)[['B']]

                                  ## Obtain H for last time-point
                                  H = rstan::extract(object$model_fit)[['H']][ , object$TS_length , , ]

                                  ## Ahead == 1
                                  for( i in 1:nrow(A) ){
                                      rr = bmgarch:::.cp( sweep( mu, 2, rts )[i, ] )
                                      H_p = Cnst[i, ,] + t( A[i, ,] ) %*%  rr %*% A[i, ,] + t( B[i, ,] ) %*% H[i, ,] %*% B[i, ,]
                                      rts_p[i, 1, ] =  mvtnorm::rmvnorm( 1, mean = mu_p[i, 1, ], sigma = H_p )                             
                                  }

                                  ## Ahead > 1
                                  if(ahead >= 2){
                                      for( p in 2:ahead ){
                                          ## location 
                                          mu_p[, p, ] =
                                              phi0  +
                                              t( sapply(1:nrow( theta ),
                                                        function(i) bmgarch:::.f_array_x_mat(mat_out = MA, array_obj = phi,
                                                                                             mat_obj = rts_p[,p-1,], i) )  ) +
                                              t( sapply(1:nrow( theta ),
                                                        function(i) bmgarch:::.f_array_x_mat( MA, theta, ( rts_p[,p-1,] - mu_p[,p-1,]) , i ) ) )
                                          ## Scale
                                          for( i in 1:nrow(A) ){
                                              rr = bmgarch:::.cp( sweep( mu, 2, rts )[i, ] )
                                              H_p = Cnst[i, ,] + t( A[i, ,] ) %*%  rr %*% A[i, ,] + t( B[i, ,] ) %*% H[i, ,] %*% B[i, ,]
                                              rts_p[i, p, ] =  mvtnorm::rmvnorm( 1, mean = mu_p[i, p, ], sigma = H_p )                             
                                          }   
                                      }
                                  }
                              }
                         
                         ## ###########
                         ## END BEKK ##
                         ## ###########
                     }
            }
            
            ## ###########################
            ## Plot and print forecasts ##
            ## ###########################
            
###################################################
            ## compute means and CrI's from predicted values ##
            pred_means = array( NA, dim = c( ahead, object$nt ) )
            pred_ci = array( NA, dim = c( ahead, 2, object$nt ) )
            
            for ( p in 1:ahead){
                pred_means[p, ] = colMeans( rts_p[, p, ] )
                pred_ci[p, ,] =  apply( rts_p[, p, ], 2, FUN = bmgarch:::.qtile )
            }

            print( list( pred_means, pred_ci ) )
            

######################
            ## Plotting options ##
######################

            ## use plots from bgarch::plot() function and add prediction

                                        #retro_plot = plot(object, type = 'means')

            rts_data = object$RTS_full

            ## TODO go into plot and define df as array, so that it saves the df for each TS?
                                        #df = retro_plot$past_data

            
            cis = array(NA, dim = c( nrow(rts_data)+ahead, 2*object$nt ) )
            cis_names =  paste0(rep( object$TS_names, each = 2 ),
                                rep (c('_lower', '_upper'), object$nt ))

            colnames(cis) = cis_names
            
            ## write CI limits into cis
            
            index = seq(from = nrow(rts_data)+1, nrow(rts_data)+ahead, by = 1)
            cis[ index,] = c(pred_ci) 

            ## append predictions
            df = data.frame ( rbind( as.matrix(rts_data), pred_means) )
            df$prediction = c( rep( 0, object$TS_length), rep( 1, ahead ) )
            df$prediction = factor(df$prediction, levels = c(0, 1), labels = c('Past', 'Forecast') )

            df$time = 1:nrow(df)

            ## add CI's to df
            df_p = data.frame(df, cis)
            
            
            colNames <- rlang::syms(names(df_p))
            
            
            for( i in 1:object$nt){
                ## find relevant positions for current variable in df_p
                var_pos = grep(object$TS_names[i], colNames )
                
                base = ggplot( data = df_p, aes( x = time , y = !!colNames[[ var_pos[1] ]], color = prediction )) + geom_line()
                base = base + geom_errorbar(aes( ymin = !!colNames[[var_pos[2]]], ymax = !!colNames[[var_pos[3]]]), width = .1) +
                    geom_point(data = df_p[ ( object$TS_length+1 ):( object$TS_length + ahead ), ], color = 'blue')
                plot(base)
                if ( i == 1){
                    devAskNewPage( ask = TRUE )
                }
            }

            
######################
            ## ... Forecast     ##
######################
        }
