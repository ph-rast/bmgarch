## #####################
##  stocks
## #####################
library(quantmod)
library(plyr)


## https://www.investing.com/search/?q=industrial
today <- format(Sys.time(), "%Y-%m-%d")
today <- "2019-11-05"

getSymbols(c("^DJI",    ## Dow Jones Industrial Average
             "^SSEI",   ## Shanghai industrial index SSEI
             "^STOXX"   ## STOXX Europe 600
             ), from = "2017-07-25", to = today)


tail(DJI$DJI.Open)
tail(SSEI$SSEI.Open)
tail(STOXX$STOXX.Open)


r2 <- cbind(DJI$DJI.Open, SSEI$SSEI.Open, STOXX$STOXX.Open)
dim(r2)
head(r2 )
r2

## Tweets on: China || billion || products || dollars || tariffs || trade
tweets <- read.csv( file = "./tariffs.csv" )
head(tweets )

tweet_time <- as.POSIXlt( as.character(tweets$created_at), format = "%m-%d-%Y %H:%M:%S")
tweet_time

## Count how many tweets in a given day - Tweets happen on WE as well. Push them over to Monday
tweet_day <- as.Date(tweet_time)
tweet_day
day <- unique(tweet_day)

## In order to push the WE tweets to Monday, we need a full time-index over the whole period
start <- range(index(r2 ))[1]
stop <- range(index(r2 ))[2]

## vector of all days
full_dates <- seq(from = start, to = stop, by = "days")

## match num tweets per day
tot_tweets <- data.frame(day, num_tweets =  0)
tot_tweets

for(i in seq_len(length(day) ) ) {    
    tot_tweets$num_tweets[i] <- length(tweet_day[tweet_day == day[i]])
}


## match all tweets to days
fd <- xts(rep(0,  length(full_dates ) ), order.by = full_dates)
ttw <- xts(tot_tweets$num_tweets, order.by = day)

fd[tot_tweets$day] <- ttw
fd
colnames(fd) <- 'tweets'


## find WE
fd$we <- 0
fd[ which(.indexwday( fd ) == 0 | .indexwday( fd ) == 6), 'we' ] <- 1
fd


fd[21:23]

for(i in seq_len(dim(fd )[1] ) ) {
    if(fd$we[i] == 1 & fd$tweets[i] != 0 ) fd$tweets[i+1] <- sum(fd$tweets[i+1], fd$tweets[i])
}

## drop weekends
fd_ww <- fd[fd$we == 0]

full <- na.omit( merge(fd_ww, r2) )
full

## predictor
pred <- full[,1]
pred_m <- cbind(pred,pred,pred)
dim(pred_m )

outcome_orig <- full[, 3:5]
outcome_orig
dim(outcome_orig )

tail(outcome_orig ,  n =  7)

outcome_log_lag <- diff(outcome_orig, lag = 1, log = TRUE)[-1,]
tail(outcome_log_lag )

dim(outcome_log_lag )
outcome_log_lag[1:3,]


## remove last day of mday and first row of rts, so that tweet dummy is at the
## day when stock opens. 
twt_pred <-  pred_m[-nrow(pred_m ), ]


dim(twt_pred )
dim(outcome_log_lag )

tail(outcome_log_lag )


## Rolling window: period of 100
## Predict next day volatility with tweet - or not

twtday <-  which(twt_pred[,1] != 0)
twtday

sel <- which(twtday > 100 )

sel[235]

## Spquared predicton error
spe_0 <- array(NA, dim = c(15, 9 ) )
spe_1 <- array(NA, dim = c(15, 9 ) )

for(i in 1:15 ) {
    
    upper <- twtday[sel[217 + i]]-1
    lower <- upper-100
    
    fit0 <- bmgarch(outcome_log_lag[ lower:upper, ], iterations = 500, P = 1, Q = 1, meanstructure = "constant",
                    parameterization = 'DCC', xH = NULL)
    #summary(fit0 )
    
    fit1 <- bmgarch(outcome_log_lag[ lower:upper, ], iterations = 500, P = 1, Q = 1, meanstructure = "constant",
                    parameterization = 'DCC', xH = twt_pred[ lower:upper, ] )
    #summary(fit1)

    fc0 <- forecast(fit0, ahead =  1, plot = FALSE, xH_p =  NULL)
    
    fc1 <- forecast(fit1, ahead =  1, xH_p = twt_pred[upper+1,], plot = FALSE)

    ## Squared predicton error
    spe_0[i, 1:3] <- (fc0[,1:3] - outcome_log_lag[ upper+1, ])^2
    spe_0[i, 4:6] <- (fc0[,4:6] - outcome_log_lag[ upper+1, ])^2
    spe_0[i, 7:9] <- (fc0[,7:9] - outcome_log_lag[ upper+1, ])^2

    spe_1[i, 1:3] <- (fc1[,1:3] - outcome_log_lag[ upper+1, ])^2
    spe_1[i, 4:6] <- (fc1[,4:6] - outcome_log_lag[ upper+1, ])^2
    spe_1[i, 7:9] <- (fc1[,7:9] - outcome_log_lag[ upper+1, ])^2
}

##
spe_0
colnames(spe_0) <- colnames(fc0)
colnames(spe_1) <- colnames(fc1)

## Mean squared prediction error:
mspe_0 <- colMeans(sqrt(spe_0))
mspe_1 <- colMeans(sqrt(spe_1))

mspe_0
mspe_1

plot(fit1, type =  'ccor' )

## in case X11 breaks due to ssh:
Sys.setenv("DISPLAY"=":11.0")

plot(spe_0[,1], type = 'l',  col = 'red',  ylim = c(0, max(spe_0) ) )
lines(spe_0[,2],  col = 'red',  lty =  2 )
lines(spe_0[,3],  col = 'red',  lty =  3 )

lines(spe_1[,1],  col = 'blue',  lty =  1 )
lines(spe_1[,2],  col = 'blue',  lty =  2 )
lines(spe_1[,3],  col = 'blue',  lty =  3 )


## 
rlag <- cbind(rl1, rl2, rl3 )[-1, ]
colnames(rlag ) <- colnames(r2 )
dim(rlag )
dim(tweet )

fit <- bmgarch(rlag, iterations = 500, P = 1, Q = 1, meanstructure = "constant",
               parameterization = 'DCC', xH = tweet)
summary(fit)

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
aussi$retro_plot[[1]] + geom_vline(data = df,  aes(xintercept = fb_inter), size =  .3, color =  'blue', alpha =  .5) + coord_cartesian( ylim = c(0,  4) , xlim = c(500,  700 ))
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
