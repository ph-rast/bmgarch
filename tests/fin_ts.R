## #####################
##  stocks
## #####################
library(quantmod)
library(plyr)
library(bmgarch)
library(bayesplot)

##########
# Stocks #
##########

getSymbols("TWTR")
TWTR$TWTR.Open[1:10]

TWTR$wday <- .indexwday( TWTR )
TWTR$mday <- ifelse( TWTR$wday == 1, 1, 0)

getSymbols("FB")
FB$FB.Open[1:10]
FB$wday <- .indexwday( FB )
FB$mday <- ifelse( FB$wday == 1, 1, 0)


getSymbols("GOOG")
GOOG$GOOG.Open[1:10]
GOOG$wday <- .indexwday( GOOG )
GOOG$mday <- ifelse( GOOG$wday == 1, 1, 0)


upper <- max(dim(cbind(FB$FB.Open,  TWTR$TWTR.Open,  GOOG$GOOG.Open)))
lower <- upper-100
leaveout <- 0
r2 <- cbind(FB$FB.Open, TWTR$TWTR.Open, GOOG$GOOG.Open)[lower:(upper-leaveout),] ## remove the last leaveout days

################
## BEKK       ##
################
## Forecasting

# If not using arma.
rlag <- scale( diff(r2, lag = 1,  log = FALSE )[-1, ]  )
colnames(rlag ) <- colnames(r2 )

# If using arma
sr2 <- scale(r2 )+5
sr2
r2

fit <- bmgarch(sr2,
               iterations = 1000,
               P = 1, Q = 1,
               meanstructure = "arma",
               standardize_data = FALSE,
               parameterization = 'DCC',
               xH = NULL,
               adapt_delta=0.85)
system("notify-send 'Done sampling' " )
summary(fit )

plot(fit, type = 'mean' )

forecast(fit, ahead = 3)
sr2[98:100,]

fit.constant <- bmgarch(rlag[,1:2],
                        iterations = 800,
                        P = 1, Q = 1,
                        meanstructure = "constant",
                        standardize_data = FALSE,
                        parameterization = "BEKK",
                        xH = NULL,
                        adapt_delta = .80)
system("notify-send 'Done sampling' " )
summary(fit.constant)
mcmc_parcoord(as.array(fit.constant$model_fit, pars = c("A","B","Cnst")), np = nuts_params(fit.constant$model_fit))

############
# Sim data #
############
sim.bekk <- function(N,C,A,B, phi = NULL, theta = NULL) {
    if(ncol(C) != nrow(C)){
        stop("C must be symmetric, square, PD.")
    }
    if(ncol(A) != nrow(A)){
        stop("A must be square.")
    }
    if(ncol(B) != nrow(B)){
        stop("B must be square.")
    }
    nt <- ncol(C)

    y <- array(0, dim = c(N, nt))
    y[1,] <- rnorm(nt, 0, sqrt(diag(C)))

    H <- array(0, dim = c(nt, nt, N))
    H[,,1] <- C

    for(i in 2:N) {
        H[,,i] <- C + t(A) %*% (t(y[i - 1,, drop = FALSE]) %*% y[i - 1,,drop = FALSE]) %*% A + t(B) %*% H[,,i-1] %*% B
        y[i,] <- MASS::mvrnorm(1, rep(0, nt), H[,,i])
    }

    if (!is.null(phi) & !is.null(theta)) {
        ## Assume phi0 (intercept) is zero.
        if (ncol(phi) != nrow(phi)) {
            stop("phi must be square [nt, nt].")
        }
        if (ncol(theta) != nrow(theta)) {
            stop("theta must be square [nt, nt].")
        }
        if (ncol(phi) != nt) {
            stop("phi must be square [nt, nt].")
        }
        if (ncol(theta) != nt) {
            stop("theta must be square [nt, nt].")
        }
        mu <- array(0, dim = c(N, nt))
        mu[1,] <- 0
        for(i in 2:N) {
            mu[i,] <- 10 + y[i - 1, , drop = FALSE] %*% phi + (y[i - 1, ,drop = FALSE] - mu[i - 1,,drop = FALSE])%*%theta
            y[i,] <- y[i,,drop = FALSE] + mu[i,,drop = FALSE]
        }
        ## y <- mu + y
    }

    return(y)
}

set.seed(13)

# nt = 2
N <-  200
C <-  matrix( c(1,  .3,  .3,  1 ) ,  ncol = 2)

A <-  matrix( c(.43,  -.07,  0.03,  .53 ) ,  ncol = 2, byrow = TRUE)
A <-  matrix( c(.34,  -.43,  -.35,  .16 ) ,  ncol = 2, byrow = TRUE)

B <-  matrix( c(.85,  -.11,  0.09,  .57 ) ,  ncol = 2, byrow = TRUE)
B <-  matrix(c(.81, -.23, .63, .31), ncol =  2, byrow = TRUE)

# nt = 3
set.seed(13)
N <- 100
nt <- 3
C_sd <- diag(rep(2, 3))
C <- C_sd %*% rethinking::rlkjcorr(1,3, 5) %*% C_sd
A <- matrix(runif(nt^2, -.5, .5), ncol=nt)
B <- matrix(runif(nt^2, -.5, .5), ncol=nt)

# ARMA(1,1)
phi <- matrix(runif(nt^2, -.5, .5), ncol = nt)
theta <- matrix(runif(nt^2, -.5, .5), ncol = nt)
phi <- matrix(0, ncol = nt, nrow = nt)
theta <- matrix(0, ncol = nt, nrow = nt)
diag(phi) <- rep(.8, nt)
diag(theta) <- rep(.5, nt)

y <- sim.bekk(N, C, A, B, phi = NULL, theta =  NULL)

fit <- bmgarch(y,
                iterations = 1000,
                P = 1, Q = 1,
                meanstructure = "constant",
                standardize_data = FALSE,
                parameterization = "pdBEKK",
                distribution = "Gaussian",
                xH = NULL,
                adapt_delta = .95)
system("notify-send 'Done sampling' " )
summary(fit)

mcmc_trace(as.array(fit$model_fit, pars = c("A","B","Cnst")))
mcmc_dens_overlay(as.array(fit$model_fit, pars = c("A","B","Cnst")))
mcmc_parcoord(as.array(fit$model_fit, pars = c("A","B","Cnst","beta0","beta1","phi","theta")), np = nuts_params(fit$model_fit))


library(bmgarch )

fit <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )],
               parameterization = "CCC", standardize_data = TRUE,
               iterations = 500)

fit1 <- bmgarch(data = stocks[1:100, c("toyota",  "nissan" )], 
                parameterization = "DCC", standardize_data = TRUE, #meanstructure = 'arma',
                iterations = 500)

fc <- forecast(fit, ahead = 1, newdata = stocks[101, c("toyota",  "nissan" )])
fc$forecast$log_lik

lfo <- loo(fit, mode = 'backward',  L = 50 )
lfo

bmgarch_objects <- bmgarch_list(fit, fit1 )
mw <-  model_weights(bmgarch_objects = bmgarch_objects, L =  80)
mw

bmgarch_objects[1]

forecast( bmgarch_objects[[1]], ahead = 1)
forecast( bmgarch_objects[[2]], ahead = 1)

class(mw )
           
fc <- forecast( bmgarch_objects, ahead = 1, weights = mw )
fc



## Figure out whther object is a list of models or just one model
## Check for nesting structure; If depth == 1, one object else list of models
.depth <- function(this,thisdepth=0){
  if(!is.list(this)){
    return(thisdepth)
  }else{
    return(max(unlist(lapply(this,.depth,thisdepth=thisdepth+1))))    
  }
}


object <- bmgarch_list(fit, fit1 )
    
.depth( fit )
n_mods <- .depth( bmgarch_objects )
n_mods

object <- bmgarch_objects
class( object )

object <- fit
object <- fit1 
##

## Define Weight for given model
wgt <- mw$wts
wgt <- c( 0.7, 0.3 )
wgt

seq_len(2)


for(i in 1:2 )

    i <- 2
object <- bmgarch_objects[[i]]
object$param
standat <- list(T = object$TS_length,
                nt = object$nt,
                rts = cbind(object$RTS_full),
                xC = object$xC,
                Q =  object$mgarchQ,
                P =  object$mgarchP,
                ahead =  ahead, 
                meanstructure =  object$meanstructure,
                distribution =  object$num_dist,
                xC_p =  xC,
                future_rts = newdata,
                compute_log_lik =  compute_log_lik)
standat

forecasted <- rstan::gqs( object = bmgarch:::stanmodels$forecastDCC,
                         draws = as.matrix(wgt_object$model_fit),data = standat)
##
forecasted


#rtsp <- rstan::extract(forecasted,  par = "rts_p" )
#Hp <- rstan::extract(forecasted,  par = "H_p" )
forecasted


## Try to recreate the structure of the list for the DCC model
## in order to be able to add the constant correlation
constant_r <- rstan::extract( fit$model_fit,  pars = "R" )

constant_r$R[,,2][,1] %*% array(1, dim = c(1, 3 ) )

lapply(constant_r, FUN = function(x ) {
    dims <- dim(x )
    apply(x, 2:length(dims), FUN = function(x) { mean(x )
#        x %*% array(1, dim = c(1, 3 ) )
    })})

lapply(constant_r, FUN = function(x ) {
    dims <- dim(x )
})

constant_r$R[,1,2]


bmgarch:::.get_stan_summary(forecasted,  "R_p",  c(.1, .9 ) )

parms <- rstan::extract(forecasted, par = c( "rts_p", "H_p") )
parms

corr_parms <- rstan::extract(forecasted, par = c( "R_p") )
corr_parms




if(i == 1 ) {
    ## write samples
    weighted_samp <-  Map("*", parms, wgt[i])    
} else if( i > 1 ) {
    ## sum over individual list elements given the weights
    weighted_tmp <- Map("*", parms, wgt[i])
    weighted_samp <-
        Map("+", weighted_samp, weighted_tmp)
}
#####
## Done combining samples
#####

#####
## Compute summary stats on weighted sample
####
weighted_samp

##
## Means 
weighted_means <- Map(colMeans,  weighted_samp)
weighted_means

c( t( weighted_means$rts_p ) )

bmgarch:::.get_stan_summary(forecasted, "rts_p",  c(.025, .975) )

## SD's 
.colSDs <- function(x) {
    lapply(x, function(x) {
        dims <- dim(x)
        apply(x, 2:length(dims), sd)
    })
}

weighted_sds <- .colSDs(weighted_samp )
weighted_sds

c( t( weighted_sds$rts_p ) )


## Quantiles
.colQTs <- function(x, probs = c(.025, .975 )) {
#    probs <- sort( c( probs, .5))
    lapply(weighted_samp, function(x) {
        dims <- dim(x)
        apply(x, 2:length(dims), quantile, probs)
    })
}

weighted_mdn <- .colQTs(weighted_samp,  probs = c(0.5) )
weighted_mdn

c( t( weighted_mdn$rts_p ) )

weighted_lower <- .colQTs(weighted_samp,  probs = c(0.025) )
c( t( weighted_lower$rts_p ) )

weighted_upper <- .colQTs(weighted_samp,  probs = c(0.095) )
c( t( weighted_upper$rts_p ) )


f.mean <- bmgarch:::.get_stan_summary(forecasted,  "rts_p",  c(0.025, .975 ) )
f.mean

f.mean[,"mean"] <- c( t( weighted_means$rts_p ) )
f.mean[,"sd"] <- c( t( weighted_sds$rts_p ) )
f.mean[,"mdn"] <- c( t( weighted_mdn$rts_p ) )
f.mean[,4] <- c( t( weighted_lower$rts_p ) )
f.mean[,5] <- c( t( weighted_upper$rts_p ) )
f.mean[,'n_eff'] <- NA
f.mean[,'Rhat'] <- NA
f.mean

colnames(f.mean )[4:5] <-
    c( paste0(min(0.025)*100, "%" ), paste0(max(0.095)*100, "%" )


f.var = bmgarch:::.get_stan_summary(forecasted,  "H_p",  c(0.025,  .975 ) )


.sort = function(x ) {
    c(apply( x, 1, FUN = function(x) {
        c(x )
    }))
}


f.var[,"mean"] <- .sort(weighted_means$H_p)
f.var[,"sd"] <- .sort(weighted_sds$H_p)
f.var[,"mdn"] <- .sort(weighted_mdn$H_p)
f.var[,4] <- .sort(weighted_lower$H_p)
f.var[,5] <- .sort(weighted_upper$H_p)
f.var[,'n_eff'] <- NA
f.var[,'Rhat'] <- NA
f.var





rstan::summary(weighted_forecasted, pars = "ar_d")$summary[, "mean"]
rstan::summary(forecasted, pars = "ar_d", probs = c(.1, .9))$summary

rstan::get_posterior_mean(forecasted)
rstan::get_posterior_mean(weighted_forecasted)
rstan::extract(forecasted, pars = "ar_d" )[[1]]



.get_stan_summary <- function(model_fit, params, CrI) {
    CrI <- c(.5, CrI)
    cols <- c("mean","sd",paste0(CrI*100, "%"), "n_eff", "Rhat")
    model_summary <- rstan::summary(model_fit, pars = params, probs = CrI)$summary[,cols]
    colnames(model_summary)[colnames(model_summary) == "50%"] <- "mdn"
    return(model_summary)
}
.get_stan_summary(forecasted, "ar_d",  c(.1, .92 ) )


.append_weight <- function(x ) {
    x$weight <- 
        }

bmgarch_objects[[1]]$weight <- 1

.weighted_forecast <- function(x, ahead = 1, xC = NULL, weight = 1) {
    fc <- forecast( x, ahead = 1, xC = xC , weight = weight)
    return(fc )
}

round( mw$wts[], 2)


## pass weight from model_weights to the forecasting function
## Map passes the corresponding elemnets of both lists to the
## forecasting function





rstan::summary(weighted_forecasted, pars = "rts_p")$summary[, "mean"]

rstan::summary(forecasted, pars = "rts_p")$summary[, "mean"]


exct <- lfocv( fit, mode = 'exact', L =  50)
awd <- lfocv( fit, mode = 'forward',  L = 50 )
bkwd <- lfocv( fit, mode = 'backward',  L =  50 )
sum(exct$out, na.rm = T)
sum(awd$out, na.rm = T)
sum(bkwd$out, na.rm = T)




log_lik_list <- list(  ll_dcc,  ll_bekk )

wts <- loo::loo_model_weights( log_lik_list, method = "stacking")
wts


