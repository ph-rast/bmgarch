##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Internal function to be used in sweep()
##' @param x Value to be squared
##' @return Squared value
##' @author philippe
.square = function(x){
    x^2
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Forecasting mean and variance 
##' @param object Fitted bmgarch object
##' @param ahead Periods to be forecasted ahead
##' @return Forecasted mean and variance   
##' @author philippe
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 geom_errorbar
forecast = function(object, ahead){
    b0 = rstan::extract(object$model_fit)[['b0']]
    b1 = rstan::extract(object$model_fit)[['b1']]
    b2 = rstan::extract(object$model_fit)[['b2']]
    b3 = rstan::extract(object$model_fit)[['b3']]
    b4 = rstan::extract(object$model_fit)[['b4']]

    ## for a GARCH 1,1 we only need last mu and last rts
    mu = rstan::extract(object$model_fit)[['mu']][, object$TS_length, ]
    rts = array( object$RTS_last ) ## TS object needs to be turned in to array for reversing
    rev = rev(rts)


    ##################
    ## DCC Forecast ##
    ##################
    
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

    mu_p = array(NA, dim = c( nrow(c_h), ahead, object$nt) )
    ## stan:  mu_p[1,] =  b0 + b1 .* (rts[T, ]-mu[T,]) + b2 .* (rev[T, ] - (sum(mu[T,]) - mu[T,]) ) + b3 .* rts[T,] +  b4 .* rev[T, ];    
    mu_p[, 1, ] = b0 + b1 * ( -sweep( mu, 2, rts ) ) + b2 * ( -sweep( mu[ ,rev(1:ncol(mu))], 2, rev )[1,] ) +
        sweep( b3, 2, rts, '*') + sweep( b4, 2, rev, '*')

    ## stan: for(d in 1:nt){
    ## stan: rr_p[1, d] = square( rts[T, d] - mu[T, d] )
    ## stan: D_p[1, d] = sqrt( c_h[d] + a_h[d]*rr_p[1, d] +  b_h[d]*D[T,d] );
    rr = t( apply( sweep( mu, 2, rts ), 1, .square) )

    D_p = array(NA, dim = c( nrow(c_h), ahead, object$nt) )
    D_p[,1,] = c_h + a_h * rr + b_h * D

    ## stan: }
    ## stan: Q_p[1,] = (1 - a_q - b_q) * S + a_q * (u[T,] * u[T,]') + b_q * Q[T,];
    Q_p = array(NA, dim = c( nrow(c_h), ahead, object$nt, object$nt) )
    Q_sdi_p = Q_p
    R_p = Q_p
    H_p = Q_p
    rts_p = array(NA, dim = c( nrow(c_h), ahead,  object$nt) )
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
            rev_p = rts_p[,p-1, rev(1:ncol( mu_p[,p-1,] ))]
            mu_p[, p, ] = b0 + b1 * ( rts_p[,p-1,] - mu_p[, p-1, ] ) + b2 * ( rev_p - mu_p[ , p-1, rev(1:ncol(mu_p[,p-1,]))] ) +
                b3 * rts_p[, p-1, ]  + b4*rev_p
            rr = (rts_p[,p-1,] - mu_p[, p-1, ])^2
            D_p[ , p , ] = c_h + a_h * rr + b_h * D_p[,p-1,]
            
            ## Loop through iteartions

            ## init u and overwrite for each p
            u_p = array(NA, dim = c( nrow(c_h), object$nt) )
            
            for( i in 1:nrow(c_h) ){

                u_p[i,] = ( 1/D_p[i, p-1, ] ) * rts_p[i, p-1,] - mu_p[i, p-1,]
                               
                Q_p[i, p,,]  = ( 1 - a_q - b_q )[i] * S[i,,] + a_q[i] * diag(u_p[i,]^2) +  b_q[i] * Q[i,,]
                
                Q_sdi_p[i, p,,] = diag( 1/sqrt(diag(Q_p[i,p,,])) )
                
                R_p[i, p ,,] = Q_sdi_p[i, p,,] %*% Q_p[i, p,,] %*% Q_sdi_p[i, p,,]
                
                H_p[i, p ,,] = diag( D_p[i,p, ] ) %*% R_p[i, p, ,] %*% diag( D_p[i,p, ] )
                
                rts_p[i,p,]  = mvtnorm::rmvnorm( 1, mean = mu_p[i, p, ], sigma = H_p[i, p ,,] )
            }
        }
    }

    ## compute means and CrI's from predicted values
    pred_means = array( NA, dim = c( ahead, object$nt ) )
    pred_ci = array( NA, dim = c( ahead, 2, object$nt ) )
    
    for ( p in 1:ahead){
        pred_means[p, ] = colMeans( rts_p[, p, ] )
        pred_ci[p, ,] =  apply( rts_p[, p, ], 2, FUN = .qtile )
    }

    print( list( pred_means, pred_ci ) )
    

    ######################
    ## Plotting options ##
    ######################

    ## use plots from bgarch::plot() function and add prediction

    #retro_plot = plot(object, type = 'means')

    rts_data = object$RTS_full

    ## TODO go into plot and define df as array, so that is saves the df for each TS?
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
    
    base = ggplot( data = df_p, aes( x = time , y = FB.Open, color = prediction )) + geom_line()
    base
    base + geom_errorbar(aes( ymin = FB.Open_lower, ymax = FB.Open_upper), width = .1) +
        geom_point(data = df_p[ ( object$TS_length+1 ):( object$TS_length + ahead ), ],
                   color = 'blue')


    
    ######################
    ## ... Forecast ##
    ####################
}
