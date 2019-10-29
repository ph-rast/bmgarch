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
lower <- upper-150
leaveout <- 0
r2 <- cbind(FB$FB.Open, TWTR$TWTR.Open, GOOG$GOOG.Open)[lower:(upper-leaveout),] ## remove the last leaveout days
dim(r2)

## Monday
mday <- cbind(tail( FB$mday, n = nrow(r2) )),
  tail( TWTR$mday, n = nrow(r2)),
  tail( GOOG$mday, n = nrow(r2)))

mday

## ## write out for stata
## r2$t <- 1:nrow(r2)
## foreign::write.dta( as.data.frame(  r2  ), file = '~/Downloads/test.dta')

################
## BEKK       ##
################
## Forecasting
library( bmgarch )

fit <- bmgarch(r2, iterations = 500, meanstructure = "constant", P = 1, Q = 1)
summary(fit)
plot(fit, type =  'means' )

forecast(fit, ahead =  3)


