## #####################
## ## actual stocks
## #####################
library(quantmod)
library(plyr)
symbols <- c("GOOG","MMM")
symbols

getQuote("AAPL")
getQuote("FB",what=yahooQF(c("Bid","Ask")))
standardQuote()


getSymbols("TWTR")
TWTR$TWTR.Open[1:10]

getSymbols("FB")
FB$FB.Open[1:10]

getSymbols("GOOG")
GOOG$GOOG.Open[1:10]

upper <- max(dim(cbind(FB$FB.Open,  TWTR$TWTR.Open,  GOOG$GOOG.Open)))
upper
lower <- upper-300
leaveout = 0
r2 <- cbind(FB$FB.Open, TWTR$TWTR.Open, GOOG$GOOG.Open)[lower:(upper-leaveout),] ## remove the last leaveout days

r2

actual <- cbind(FB$FB.Open, TWTR$TWTR.Open)[(upper-(leaveout-1)):upper,] 
actual

r2 <- na.omit(r2)
r2

names(r2)

################
## BEKK       ##
################
## Forecasting



fit1 = bmgarch(data = r2[,1:2], iterations = 1500, parameterization = 'DCC')

summary(fit1)
plot(fit1, type = 'cvar') 
plot(fit1, type = 'means' ) 
plot(fit1, type = 'ccor' ) 

forecast(fit1, ahead = 20)

## Cf. https://groups.google.com/forum/#!topic/stan-users/tWQdtndbSnA for failed initalization: init_r < 2
bekk_fit_0 <- sampling(bekk_mod, data = standat, verbose = TRUE, iter = 2000, control = list(adapt_delta = .99), init_r = 1, chains = 4)
dcc_fit <- sampling(dcc_mod, data = standat, verbose = TRUE, iter = 1900, control = list(adapt_delta = .99), init_r = 1, chains = 4)

## inits for ccc
initf2 <- function(chain_id = 1, L) {
  list(H = sigma1,
       R = cov2cor(sigma1))
}
n_chains <- 4
inits <- lapply(1:n_chains, function(id) initf2(chain_id=1:N, L))
str(inits)

inits

ccc_fit <- sampling(ccc_mod, data = standat, verbose = TRUE, iter = 1900, control = list(adapt_delta = .99), init_r = .1, chains = 4)

print(bekk_fit_0, pars = c('Cnst', 'A', 'B', 'corC', 'b0', 'b1', 'b2', 'b3', 'b4', 'lp__'), probs = c(0.05, 0.95))
print(dcc_fit, pars = c('a_q', 'b_q', 'c_h', 'a_h', 'b_h', 'S', 'b0', 'b1', 'b2', 'b3', 'b4', 'lp__'), probs = c(0.05, 0.95))
print(ccc_fit, pars = c('c_h', 'a_h', 'b_h', 'R', 'b0', 'b1', 'b2', 'b3', 'b4', 'lp__'), probs = c(0.05, 0.95))


bekk_fit = ccc_fit

pairs(bekk_fit, pars = c('A'))
traceplot(bekk_fit_0, pars = c('A', 'B'), inc_warmup = T)



######################33
## Check with rmgarch
library('rmgarch')

Dat <- r2
xspec = ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(garchOrder = c(1,1), model = 'sGARCH'), distribution.model = 'norm')
uspec = multispec(replicate(2, xspec))
spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
spec1a = dccspec(uspec = uspec, dccOrder = c(1, 1), model='aDCC', distribution = 'mvnorm')

cl = makePSOCKcluster(4)
multf = multifit(uspec, Dat, cluster = cl)

system.time({fit1 = dccfit(spec1, data = Dat, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)})
fit_adcc = dccfit(spec1, data = Dat, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
print(fit1)           
print(fit_adcc)

stopCluster(cl)


plot(dccforecast(fit = fit1, n.ahead = 5))



## forecasted params:
ahead <- standat$ahead
ahead

size <- dim(extract(bekk_fit)[['rts_p']])
size
## series 1
extract(bekk_fit)[['rts_p']][1:10,1,1]

## extract 
p1 <- array(extract(bekk_fit)[['rts_p']][,1:ahead,1], dim = c(size[1], size[2]))
p1

## series 2
p2 <- array(extract(bekk_fit)[['rts_p']][,1:ahead, 2], dim = c(size[1], size[2]))

p3 <- array(extract(bekk_fit)[['rts_p']][,1:ahead, 3], dim = c(size[1], size[2]))


muhat = colMeans(extract(bekk_fit)[['mu']][,,1])

qtile = function(x){
  cis = quantile(x, c(.025, .975) )
  return(cis)
}

qtile(extract(bekk_fit)[['mu']][1:10,1,1])

CIs = apply(extract(bekk_fit)[['mu']][1:10,,1], 2, qtile )

plot(as.numeric(r2[,1]), type = 'l')
lines(muhat, lty = 2)
lines(CIs[1,], lty = 3, col = rgb(1,0,0, .9))
lines(CIs[2,], lty = 3, col = rgb(1,0,0, .9))

##pdf(file = '../FIGURES/CorrInH.pdf', width = 5, height = 3)
plot(r2[,1])

plot(as.numeric(r2[,1]), type = 'l', ylim = c(-2.5, 5), xlim = c(0, (nrow(r2) + ahead + 15)))
CI.Lpct <- round(nrow(p1)*0.025, 0)  ## 95% CI
CI.Upct <- round(nrow(p1)*0.9755, 0)
for( i in 1:ahead){
    points((nrow(r2) + i), mean(p1[,i]), col = 'gray40')
    CI.L <- sort(p1[,i])[CI.Lpct]
    CI.U <- sort(p1[,i])[CI.Upct]
    lines(rep((nrow(r2) + i),2), c(CI.L, CI.U), col =  'gray70')
}

#
lines(as.numeric(r2[,2]), type = 'l', ylim = c(-3, 3), xlim = c(0, (nrow(r2) + ahead)), col = '#00b0e7')
for( i in 1:ahead ){
    points((nrow(r2) + i + .3), mean(p2[,i]), col = '#004bbc')
    CI.L <- sort(p2[,i])[CI.Lpct]
    CI.U <- sort(p2[,i])[CI.Upct]
    lines(rep((nrow(r2) + i + .5),2), c(CI.L, CI.U), col =  '#add8e6')
}

#
##op
lines(as.numeric(r2[,3]), type = 'l', ylim = c(-3, 3), xlim = c(0, (nrow(r2) + ahead)), col = '#0000e7')
for( i in 1:ahead ){
  points((nrow(r2) + i + .6), mean(p3[,i]), col = '#0000e790')
  CI.L <- sort(p3[,i])[CI.Lpct]
  CI.U <- sort(p3[,i])[CI.Upct]
  lines(rep((nrow(r2) + i + .6),2), c(CI.L, CI.U), col =  '#0000e750')
}


obs <- length(r2[,1])
obs

xact <- (obs+1):(obs+ahead)
xact

## rescale to r2 values using
## attr(r2, "scaled:center")
## attr(r2, "scaled:scale")

attr(r2, "scaled:center")
attr(r2, "scaled:scale")

actual

## obtain actual values and standardize them with same location and scale
actual.std = cbind( (actual[,1] - attr(r2, "scaled:center")[1])/attr(r2, "scaled:scale")[1],
                    (actual[,2] - attr(r2, "scaled:center")[2])/attr(r2, "scaled:scale")[2])
                    #(actual[,3] - attr(r2, "scaled:center")[3])/attr(r2, "scaled:scale")[3] )

actual.std
xact

points(xact[1]   , actual.std[1:ahead-1 ,1], pch = "*")
points(xact+.3, actual.std[1:ahead-1, 2], pch = "*", col = 'blue')
points(xact+.6, actual.std[1:ahead, 3], pch = "*", col = 'red')


prederror <- array(NA, dim = c(15,2))
prederror[1,1] <- mean(p1[,1]) - actual.std[1,1] ## if ahead is not == 1 then specify column in p1()
prederror[1,2] <- mean(p2[,1]) - actual.std[1,2]
prederror[1,3] <- mean(p3[,1]) - actual.std[1,3]


prederror

actual.std

dim(actual.std)[1]

for(i in 2:(dim(actual.std)[1])){
    ## attach actual day datapoint, to r2 and rerun
    ## append actual.std to the r2 vectors
    r2 <- cbind(rbind(r2[,1], actual.std[i-1,1]), rbind(r2[,2], actual.std[i-1,2]) )#, rbind(r2[,3], actual.std[i-1,3]) ) 
    standat <- list(T = nrow(r2), rts = r2, sigma1 = sigma1, nt = ncol(r2), ahead = 1)
    bekk_fit <- sampling(bekk_mod, data = standat, verbose = TRUE, iter = 1000, control = list(adapt_delta = .99), init_r = 1, chains = 4)
    p1 <- array(extract(bekk_fit)[['rts_p']][,1:ahead,1], dim = c(size[1], size[2]))
    p2 <- array(extract(bekk_fit)[['rts_p']][,1:ahead,2], dim = c(size[1], size[2]))
#    p3 <- array(extract(bekk_fit)[['rts_p']][,1:ahead,3], dim = c(size[1], size[2]))
  #
    for( k in 1:ahead){
        points((nrow(r2) + k), mean(p1[,k]), col = 'gray40')
        CI.L <- sort(p1[,k])[CI.Lpct]
        CI.U <- sort(p1[,k])[CI.Upct]
        lines(rep((nrow(r2) + k),2), c(CI.L, CI.U), col =  'gray70')
    }
    for( q in 1:ahead ){
        points((nrow(r2) + q + .3), mean(p2[,q]), col = '#004bbc')
        CI.L <- sort(p2[,q])[CI.Lpct]
        CI.U <- sort(p2[,q])[CI.Upct]
        lines(rep((nrow(r2) + q + .3),2), c(CI.L, CI.U), col =  '#add8e6')
#        points((nrow(r2) + q + .6), mean(p3[,q]), col = '#0000e790')
#        CI.L <- sort(p3[,q])[CI.Lpct]
#        CI.U <- sort(p3[,q])[CI.Upct]
#        lines(rep((nrow(r2) + q + .6),2), c(CI.L, CI.U), col =  '#0000e750')
    }
    obs <- length(r2[,1])
    xact <- (obs+1):(obs+ahead)
    points(xact, as.numeric(actual.std[i,1]), pch = "*")
    points(xact+.3, as.numeric(actual.std[i,2]), pch = "*", col = 'blue')
#    points(xact+.6, as.numeric(actual.std[i,3]), pch = "*", col = 'red')
    ## record difference
    prederror[i,1] <- mean(p1[,1]) - actual.std[i,1]
    prederror[i,2] <- mean(p2[,1]) - actual.std[i,2]
 #   prederror[i,3] <- mean(p3[,1]) - actual.std[i,3]
}

prederror

colSums(prederror^2)/ncol(prederror)
##[1] 0.5493340 0.7854347
sd1 = attr(r2, "scaled:scale")[1]
sd2 = attr(r2, "scaled:scale")[2]


prederror_orig <- cbind(prederror[,1]*sd1, prederror[,2]*sd2)
prederror_orig

colSums(prederror_orig)
#lines(187:200+.6,(actual[,2]-mn2)/sd2)


for ( k in 1:15){
    sel <- 187:201
    lines(c(sel[k],sel[k]), c( (actual[k,1]-mn1)/sd1, (actual[k,1]-mn1)/sd1+prederror[k,1]))
    lines(c(sel[k],sel[k]), c( (actual[k,2]-mn2)/sd2, (actual[k,2]-mn2)/sd2+prederror[k,2]), col = 'blue')
}

yrep <- extract(bekk_fit)[['rts_out']]

dim(yrep)

yrep[1,1,]
yrep[1,2,]
yrep[1, ,52]

n_ts <- ncol(r2)

op <- par(mfcol = c(1,n_ts))
for( k in 1:n_ts){
    plot(yrep[1,k,], type = 'l', ylim = c(-10, 10), col = "#81939920")
    for(i in sample(nrow(yrep), size = 500)){
        lines(yrep[i,k,], col =  '#81939920')
    }
    lines(
        r2[,k]
      , col = 'red', lwd = 4)   
}
op

library(bayesplot)

color_scheme_set("red")
ppc_dens_overlay(y=r2[,2],
                 yrep = yrep[,2,], trim = T)


op <- par(mfcol = c(1,n_ts))
for( k in 1:n_ts){
    plot(density(r2[,k]))
    for(i in sample(nrow(yrep))){
        lines(density(
            yrep[i,k,]
        ), col =  '#81939920')
    }
    lines(density(r2[,k]), col = 'red', lwd = 2)
}
op
