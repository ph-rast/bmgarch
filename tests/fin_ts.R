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
lower <- upper-300
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

rlag <- scale( diff(r2, lag = 1,  log = TRUE )[-1, ]  )

rlag[1:5, ]
colnames(rlag ) <- colnames(r2 )

library(bmgarch)

tweet
dim(rlag)

range(rlag[,1:2])

head(r2 )
tail(r2)

sr2 <- scale(r2 )
sr2

y <- r2[, 1:2]

ys1 <-  ( y[,1] - as.numeric(y[1,1]) ) / sd(y[, 1] )
ys2 <-  ( y[,2] - as.numeric(y[1,2]) ) / sd(y[, 2] )

ys <- cbind(ys1, ys2 )
ys


fit <- bmgarch(sr2[,1:2],
               iterations = 800,
               P = 1, Q = 1,
               meanstructure = "arma",
               standardize_data = FALSE,
               parameterization = 'DCC',
               xH = NULL,
               adapt_delta=0.80)
system("notify-send 'Done sampling' " )
summary(fit )

fit5 <- bmgarch(rlag[,1:2]*5,
               iterations = 6000,
               P = 1, Q = 1,
               meanstructure = "constant",
	       standardize_data = FALSE,
               parameterization = 'BEKK',
	       xH = NULL,
               adapt_delta=0.99)
system("notify-send 'Done sampling' " )
summary(fit5 )

rstan::traceplot(fit$model_fit,  pars =  c('phi0' ) ,  inc_warmup =  T)

forecast( fit, ahead =  25, type = "mean", CrI =  c(.025, .975), last_t = 50)
fit$RTS_full
fit$TS_names

plot(fit, type =  'ccor' )

library(astsa )
astsa::sarima(r2[,1] , p = 1,  q = 1, d = 1)

stata <- read.csv(file =  '~/Downloads/statadata.csv' )



stata_fit <- bmgarch(stata[1700:2000,1:2]*10,  parameterization = 'DCC',  standardize_data = FALSE )
summary(stata_fit )
plot(stata_fit,  type =  'ccor' )
forecast(stata_fit, type = "var", ahead =  25)
forecast(stata_fit, type = "cor", ahead =  25)
forecast(stata_fit, type = "mean", ahead =  25)
