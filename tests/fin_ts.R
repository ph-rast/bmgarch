## #####################
##  stocks
## #####################
library(quantmod)
library(plyr)
#symbols <- c("GOOG","MMM")
#symbols

#getQuote("AAPL")
#getQuote("FB",what=yahooQF(c("Bid","Ask")))
#standardQuote()


getSymbols("TWTR")
TWTR$TWTR.Open[1:10]
## tail( TWTR )

## class( TWTR )

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
upper
lower <- upper-116
leaveout <- 0
r2 <- cbind(FB$FB.Open, TWTR$TWTR.Open, GOOG$GOOG.Open)[lower:(upper-leaveout),] ## remove the last leaveout days
dim(r2)
head(r2 )

## Monday
mday <- cbind(tail( FB$mday, n = nrow(r2) ),
  tail( TWTR$mday, n = nrow(r2)),
  tail( GOOG$mday, n = nrow(r2)))


head(mday)

## Tweets from
## http://trumptwitterarchive.com/archive
goog <- read.csv( file = "./goog.csv" )
twtr <- read.csv( file = "./twtr.csv" )
fb <- read.csv( file = "./fb.csv" )

fb_tweet <- as.POSIXlt( as.character(twtr$created_at), format = "%m-%d-%Y %H:%M:%S")
twtr_tweet <- as.POSIXlt( as.character(twtr$created_at), format = "%m-%d-%Y %H:%M:%S")
goog_tweet <- as.POSIXlt( as.character(goog$created_at), format = "%m-%d-%Y %H:%M:%S")

ttweet <- data.frame( tweet, 1)
goog_tweet
twtr_tweet
fb_tweet
as.Date(goog_tweet)

mday[, 1:3] <- 0
head(mday)

mday[as.Date(fb_tweet), 'mday'] <- 1
mday[as.Date(twtr_tweet), 'mday.1'] <- 1
mday[as.Date(goog_tweet), 'mday.2'] <- 1

colnames(mday) <- c("fb",  "twtr",  "goog" )

head(mday )
head(r2 )
## ## write out for stata
## r2$t <- 1:nrow(r2)
## foreign::write.dta( as.data.frame(  r2  ), file = '~/Downloads/test.dta')

dim(mday )

nrow(mday )
## remove last day of mday and first row of rts, so that tweet dummy is at the
## day when stock opens. 
tweet <-  mday[-nrow(mday ), ]
rts <- r2[-1,]

dim(tweet ) == dim(rts )
################
## BEKK       ##
################
## Forecasting

#fit <- bmgarch(rts, iterations = 500, P = 1, Q = 1, meanstructure = "arma",
#               parameterization = 'DCC', xH = tweet)

#summary(fit)


#forecast(fit, ahead =  10)


## Lagged
## differences

rlag <- diff(r2, lag = 1,  log = TRUE )[-1, ]

rlag[1:5, ]
colnames(rlag ) <- colnames(r2 )

fit <- bmgarch(rlag, iterations = 500, P = 1, Q = 1, meanstructure = "constant",
               parameterization = 'DCC', xH = NULL)
summary(fit)

forecast( fit, ahead =  5, type = "cor", CrI =  c(.025, .975))

aussi <- plot(fit, type =  'cvar' )



## ###########
## Plotting
## ##########

dim(tweet )
df <- data.frame(tweet[,1:3], period =  seq_len( dim(tweet)[1]))
head(df)
df <- df[-1,]
df <- rbind(df, c(0, 0, 0 ) )
df$fb[df$fb == 0] <- NA
df$goog[df$goog == 0] <- NA
df$twtr[df$twtr == 0] <- NA

df$fb_inter <- df$fb*df$period

aussi$retro_plot[[1]] + geom_point(data = df,  aes(x = period, y = fb, size =  2), color =  'blue', shape =  18) + coord_cartesian( ylim = c(0,  4) )
aussi$retro_plot[[2]] + geom_point(data = df,  aes(x = period, y = goog, size =  2), color =  'blue', shape =  18) + coord_cartesian( ylim = c(0,  4) )
aussi$retro_plot[[3]] + geom_point(data = df,  aes(x = period, y = twtr, size =  2), color =  'blue', shape =  18) + coord_cartesian( ylim = c(0,  4) )

dfraw <- data.frame( rlag[,1:3], df)
head(dfraw )

ggplot(dfraw, aes(x = period, y = FB.Open ) ) + geom_line() + geom_point( aes(x = period, y = fb, size =  2), color =  'blue', shape =  18) + coord_cartesian( ylim = c(-10,  10) )

dfraw <- data.frame( rts[,1:3], df)

## Place diamond on actual TS
df <- data.frame(tweet[,1:3], period =  seq_len( dim(tweet)[1]))
realdt <- (rts[,1:3] * df)[,1:3]
realdt[realdt == 0] <- NA
head(realdt )

dfraw1 <- data.frame( rts[,1:3], realdt,  period =  seq_len( dim(tweet)[1]))
dfraw1

ggplot(dfraw1, aes(x = period, y = FB.Open ) ) + geom_line() + geom_point( aes(x = period, y = fb, size =  2), color =  'red', shape =  18)# + coord_cartesian( ylim = c(-10,  10) )



## ### without tweets
fit0 <- bmgarch(rlag, iterations = 500, P = 1, Q = 1, meanstructure = "constant",
                parameterization = 'DCC', xH = NULL)

summary(fit0 )

aussi0 <- plot(fit0, type =  'cvar' )

aussi0$retro_plot[[1]]+ geom_point(data = df,  aes(x = period, y = fb, size =  2), color =  'blue', shape =  18) + coord_cartesian( ylim = c(0,  4) )

library(loo )
log_lik_fit0 <- extract_log_lik(fit0$model_fit )
log_lik_fit <- extract_log_lik(fit$model_fit )

loo_fit0 <- loo(log_lik_fit0 )
loo_fit <- loo(log_lik_fit )
waic_fit0 <- waic(log_lik_fit0 )
waic_fit <- waic(log_lik_fit )

compare(loo_fit0, loo_fit)
compare(waic_fit0, waic_fit)



i <- 3
rt_rep <- rstan::extract(fit$model_fit )[['rts_out']][,i,]
bayesplot::color_scheme_set("red")
bayesplot::ppc_dens_overlay(y = array(scale(rlag[,i])), ## Depending on whether bmgarch used scaling or not - needs to be done here too
                            yrep = rt_rep[sample(1:1000, size = 150),]) +
    coord_cartesian( xlim = c(-10,  10)) +
    ggtitle(colnames(rlag)[i])


## Drop some objects
## fit and fit0 are almost identical
## try some models with different lags

fit22 <- bmgarch(rlag, iterations = 500, P = 2, Q = 2, 
               parameterization = 'DCC', xH = tweet)
fit33 <- bmgarch(rlag, iterations = 500, P = 3, Q = 3, 
               parameterization = 'DCC', xH = tweet)
fit44 <- bmgarch(rlag, iterations = 500, P = 4, Q = 4, 
               parameterization = 'DCC', xH = tweet)

####
## Compare forecasting performance
###
rowSums(tweet)
twtday <-  which(rowSums(tweet ) != 0)
twtday

sel <- which(twtday > 175 )
twtday
length(sel)

## Spquared predicton error
spe_0 <- array(NA, dim = c(length(sel), 9 ) )
spe_1 <- array(NA, dim = c(length(sel), 9 ) )


for(i in 1:length(sel) ) {

    upper <- twtday[sel[i]]-1
    lower <- upper-175

    fit0 <- bmgarch(rlag[ lower:upper, ], iterations = 500, P = 1, Q = 1, meanstructure = "constant",
                    parameterization = 'DCC', xH = NULL, standardize_data = TRUE)
    #summary(fit0 )

    fit1 <- bmgarch(rlag[ lower:upper, ], iterations = 500, P = 1, Q = 1, meanstructure = "constant",
                    parameterization = 'DCC', xH = tweet[ lower:upper, ] , standardize_data = TRUE)
    # summary(fit1)                                                                                                                                                                                                

    fc0 <- forecast(fit0, ahead =  1, plot = FALSE, xH_p =  NULL)

    fc1 <- forecast(fit1, ahead =  1, xH_p = tweet[upper+1,], plot = FALSE)

    ## Squared predicton error
    spe_0[i, 1:3] <- (fc0[,1:3] - rlag[ upper+1, ])^2
    spe_0[i, 4:6] <- (fc0[,4:6] - rlag[ upper+1, ])
    spe_0[i, 7:9] <- (fc0[,7:9] - rlag[ upper+1, ])
    
    spe_1[i, 1:3] <- (fc1[,1:3] - rlag[ upper+1, ])^2
    spe_1[i, 4:6] <- (fc1[,4:6] - rlag[ upper+1, ])
    spe_1[i, 7:9] <- (fc1[,7:9] - rlag[ upper+1, ])
}

## note, CrI's are differences to observed value
## of lower always negative and upper always postive,
## prediction was acutally in the CrI range
colnames(spe_0) <- colnames(fc0)
colnames(spe_1) <- colnames(fc1)
spe_0
spe_1
## Mean squared prediction error:                                                                                                                                                                         
mspe_0 <- colMeans(sqrt(spe_0[,1:3]))
mspe_1 <- colMeans(sqrt(spe_1[,1:3]))

mspe_0
mspe_1

## tweets lead to lower performance
