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
lower <- upper-300
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
sr2 <- scale(r2 )
sr2

fit <- bmgarch(sr2[,1:2],
               iterations = 800,
               P = 1, Q = 1,
               meanstructure = "arma",
               standardize_data = FALSE,
               parameterization = 'BEKK',
               xH = NULL,
               adapt_delta=0.80)
system("notify-send 'Done sampling' " )
summary(fit )

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

set.seed(13)
N <-  100
C <-  matrix( c(2,  0.5,  0.5,  2 ) ,  ncol = 2)
A <-  matrix( c(.4,  0.1,  -0.3,  .2 ) ,  ncol = 2)
B <-  matrix( c(.2,  0.1,  0.3,  .3 ) ,  ncol = 2)
y <- array(c(.5, -.5 ),  dim = c(N,  2) )
H <- array(C, dim = c(2, 2, N ) )
for( i in 2:N ) {
    H[,,i] <- C + t(A) %*% (y[i-1,] %*% t(y[i-1,])) %*% A + t(B) %*% H[,,i-1] %*% B
    y[i, ] <- MASS::mvrnorm(1,  c(0, 0 ),  H[,,i] )
}
y

fit.constant <- bmgarch(y,
                        iterations = 1000,
                        P = 1, Q = 1,
                        meanstructure = "constant",
                        standardize_data = FALSE,
                        parameterization = "BEKK",
                        distribution = "Gaussian",
                        xH = NULL,
                        adapt_delta = .95)
system("notify-send 'Done sampling' " )
summary(fit.constant)
A
B
C
## summary(fit.constant$model_fit, pars = c("A","B","Cnst"))$summary[,c("mean","2.5%","97.5%","Rhat")]

mcmc_trace(as.array(fit.constant$model_fit, pars = c("A","B","Cnst")))
mcmc_dens_overlay(as.array(fit.constant$model_fit, pars = c("A","B","Cnst")))
