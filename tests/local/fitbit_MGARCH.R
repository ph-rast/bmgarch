## Project title: 
##    Created at: 
##        Author: Philippe Rast
##          Data: 
##       Summary: Reaction to Reviewer 2 with GARCH
## ---------------------------------------------------------------------- ##
rm(list=ls())


#####################
## Siwei's fitbit  ##
####################
fitbit <- read.csv(file = '~/UZH/Projekte/M_GARCH/WORK/-DATASETS/SOURCE/fit_bit_dat.csv')
fitbit <- read.csv(file = '~/UZH/Projekte/M_GARCH/WORK/-DATASETS/DERIVED/the_data.csv')
names(fitbit)
head(fitbit)

stem(fitbit$Mean_TemperatureF)
stem(fitbit$Max_TemperatureF)

variables = names(fitbit)[9:29]

length(variables)
#op <- par(mfcol = c(3,7))
#for(i in 1:length(variables)) {
#    plot(fitbit[fitbit$record_id ==1, variables[i]], type = 'l', main = variables[i])
#}
                                        #op


## Bivariate mgarch for one individual. Ts are steps and some variable
fitsub <- fitbit[,c( 'record_id', 'day', 'steps', 'couple_id', 'Max_TemperatureF', 'fitbitDate')]
fitsub


fitsub$NAf <-  scale(
    rowSums( fitbit[,c( "disinterested", "upset", "guilty", "scared", "hostile", "irritable",
                      "ashamed", "nervous", "jittery", "afraid", "stressed")] , na.rm = TRUE)
    )

fitsub$PAf <-scale(
    rowSums( fitbit[,c( "interested", "excited", "strong", "enthusiastic", "proud", "alert", "inspired",
                        "determined", "attentive", "active")], na.rm = TRUE)
    )

fitsub$PAf
fitsub$NAf

fitsub$stepsstd <- scale(fitbit$steps)

subs <- unique(fitsub$record_id)
subs

## Get Days of Week for WE effect in stressed
fitsub$weekday <- weekdays(as.Date(fitsub$fitbitDate, '%Y-%m-%d'))
fitsub$we_dummy = ifelse( fitsub$weekday == "Saturday" | fitsub$weekday == "Sunday", 1, 0)


## obtain max day count per person:
ids = unique(fitsub$record_id)
ids
for( i in 1:length(ids) ){
    fitsub[fitsub$record_id == ids[i], 'totday'] = length(fitsub[fitsub$record_id == ids[i], 'day'])
}

## individual 37 has a weirdly negative reading on PAf in row 48. Set that to zero - now it works

fitsub[fitsub$record_id==37,'PAf'][48] <- 0

head(fitsub)


########################################
## Loop through participants
########################################

########################################
## BEKK
########################################


# library(tcltk)

sel = unique(fitsub[fitsub$totday >= 90, "record_id"])
sel = sel[!sel %in% redo]

length(sel)

results <- array(NA, dim = c(length(sel),22))
resultsCI <- array(NA, dim = c(length(sel),44))
ahead = 1

names(fitsub)

## loop through sel participants
total <- length(sel)
# create progress bar
pb <- txtProgressBar(min = 0, max = total, style = 3)
for(i in 1:length(sel)){

#for(i in redo){

    ## Extract Time series per invididul   
    #i = 4#13
    #i = 10
    r <- fitsub[fitsub$record_id==sel[i],c('PAf', 'NAf')]
    xH <- fitsub[fitsub$record_id==sel[i],c('we_dummy')] 
    dim( r )
    r

    #plot(1:nrow(r), r[,1], type = 'l'); lines(1:nrow(r), r[,2], col = 'red')
    #rug(x = xH*fitsub[fitsub$record_id==sel[i], 'day'])

    xH = cbind( xH, xH )
    ##xH
    
    # rl = r
    # cbind(rl, r)
    # rl[,1]= quantmod::Lag(r[,1], 1)
    # rl[,2]= quantmod::Lag(r[,2], 1)
    # rl2 = r - rl
    # plot(1:nrow(rl2), rl2[,1], type = 'l')
    # lines(1:nrow(rl2), rl2[,2], col = 'red')
     
    fit1 = bmgarch(data = r, xH = xH, parameterization = "BEKK", iterations = 500,
                   distribution = 'student_t', Q = 1, P = 1, standardize_data = TRUE)
    #summary(fit1)
    
    #round(rstan::summary(fit1$model_fit, pars = c( 'beta', 'lp__'))$summary[,c('mean', '2.5%', '97.5%', 'Rhat')], 2)

    #class(fit1)
    #plot( fit1, type = 'means' )
    
    # summary(fit1, CrI = c(0.025, .975))
    # bmgarch::plot.bmgarch(fit1, type = 'cvar')
    # bmgarch::plot.bmgarch(fit1, type = 'means')
    # bmgarch::plot.bmgarch(fit1, type = 'ccor')
    # forecast(fit1, ahead = 3)
    
    pst <- rstan::summary(fit1$model_fit, pars = c('Cnst', 'A', 'B', 'corC',  'phi0', 'beta[1,1]', 'beta[1,2]', 'beta[2,2]', 'lp__'), probs = c(0.05, 0.95))$summary[,c('mean', '5%', '95%', 'Rhat')]
    dim(pst)
    ## Overwrite parameter if Rhat < 1.1
    results[i, ] <- pst[,1]*ifelse(pst[,'Rhat']>1.1, NA, 1)
    colnames(results) <- rownames(pst)
    ## Record lower and upper CI
    resultsCI[i,] <-  c(t(pst[,2:3]))
    low <- paste0(rownames(pst), '5%')
    hi <- paste0(rownames(pst), '95%')
    colnames(resultsCI)  <-   c(t(cbind(low, hi)))    
    ## Progress Bar
    setTxtProgressBar(pb, i)
}

close(pb)

## focus on b2 for now?

results
resultsCI[35,]

redo <- which(is.na(rowSums(results)))
length(redo)
redo
total

save(results, file = 'results_fitbit_bekk.RDat', compress = 'xz')
save(resultsCI, file = 'resultsCI_fitbit_bekk.RDat', compress = 'xz')

load('../results_fitbit_bekk.RDat')

#load(file = 'results_dyad.RDat')

## Create filter for results with posterior mass outside of zero
nonzero <- ifelse(sign(resultsCI[, seq(1, 44, by = 2)]) + sign(resultsCI[, seq(2, 44, by = 2)]) == 0, NA, 1)
nonzero

results[ !is.na( rowSums(results) ), ]

## Return nonzero results:
round(results * nonzero, 2) [ !is.na( rowSums(results) ), ]



## Predict some stuff with it.
sel

subpop = fitsub[fitsub$record_id %in% sel, ]

names(subpop)

subpop


########################################
## DCC
########################################


results_dcc <- array(NA, dim = c(length(sel),20))
resultsCI_dcc <- array(NA, dim = c(length(sel),40))
ahead = 1

# create progress bar
pb <- txtProgressBar(min = 0, max = total, style = 3)
for(i in 1:length(sel)){
  
for(i in redo){
    ## Extract Time series per invididul 
    r <- fitsub[fitsub$record_id==sel[i],c('PAf', 'NAf')]  
    r <- fitbit[fitbit$record_id==sel[i],c("disinterested", "upset", "guilty", "scared", "hostile")]
    r = na.omit(r)
    fit1 = bmgarch(data = r, parameterization = "DCC", iterations = 400, chains = 4)

    pst <-
        rstan::summary(fit1$model_fit, c('c_h', 'a_h', 'b_h', 'a_q', 'b_q', 'S[1,2]', 'phi0', 'phi', 'theta', 'lp__'))$summary[,c('mean', '2.5%', '97.5%', 'Rhat')]
    ## Overwrite parameter if Rhat < 1.1
    results_dcc[i, ] <- pst[,1]*ifelse(pst[,'Rhat']>1.1, NA, 1)
    colnames(results_dcc) <- rownames(pst)
    ## Record lower and upper CI
    resultsCI_dcc[i,] <-  c(t(pst[,2:3]))
    low <- paste0(rownames(pst), '2.5%')
    hi <- paste0(rownames(pst), '97.5%')
    colnames(resultsCI_dcc)  <-   c(t(cbind(low, hi)))    
    ## Progress Bar
    setTxtProgressBar(pb, i )
}

close(pb)

## focus on b2 for now?

results_dcc
resultsCI_dcc

redo <- which(is.na(rowSums(results_dcc)))
length(redo)
redo


cbind(1,rowSums(results_dcc))

save(results_dcc, file = 'results_fitbit_dcc.RDat', compress = 'xz')
save(resultsCI_dcc, file = 'resultsCI_fitbit_dcc.RDat', compress = 'xz')



## Create filter for results with posterior mass outside of zero
nonzero_dcc <- ifelse(sign(resultsCI_dcc[, seq(1, 40, by = 2)]) + sign(resultsCI_dcc[, seq(2, 40, by = 2)]) == 0, NA, 1)

nonzero

## Return nonzero results:
round(results_dcc * nonzero_dcc, 2)


## Look at some interesting cases. Eg. 4 and 3 with S[1,2] .80 and -.50
sel[4]

r <- fitsub[fitsub$record_id==sel[4],c('PAf', 'NAf')]
pdf(file = "no4_30.pdf")
op <- par(mfcol = c(2,1))
plot(r[,'PAf'], type = 'l', ylim = c(-3, 3.5) , main = "Individual with resid. corr: .80", ylab = 'outcome on z-scale', col = 'red')
lines(r[, 'NAf'], col = 'blue')
text(x = 80, y = -2.5, labels = 'Positive Affect', col = 'red')
text(x = 80, y = 3, labels = 'Negative Affect', col = 'blue')
r <- fitsub[fitsub$record_id==sel[3],c('PAf', 'NAf')] 
plot(r[,'PAf'], type = 'l', ylim = c(-3, 3.5) , main = "Individual with resid. corr.: -.50", ylab = 'outcome on z-scale', col = 'red')
lines(r[, 'NAf'], col = 'blue')
text(x = 80, y = -2.5, labels = 'Positive Affect', col = 'red')
text(x = 80, y = 3, labels = 'Negative Affect', col = 'blue')
op
dev.off()


## check out row 18, with respect to garch params (crazy data...)
r <- fitsub[fitsub$record_id==sel[18],c('PAf', 'NAf')]
pdf(file = "no18.pdf") 
plot(r[,'PAf'], type = 'l', ylim = c(-3, 3.5) , main = "a_h[1]:2.28, a_h[2]: 0.44; b_h's: 0.5",
     ylab = 'outcome on z-scale', col = 'red')
lines(r[, 'NAf'], col = 'blue')
dev.off()



results_dcc


######################
## Relate results to steps
#####################
## Create person median for steps:

med_steps = aggregate(steps ~ record_id, data = fitsub, FUN=median)

## subset with sel
md_steps = med_steps[sel,]

cor( cbind(results, md_steps[,2]), use = "pairwise.complete.obs")

cor( cbind(results_dcc, md_steps[,2]), use = "pairwise.complete.obs")


###################################
####### Pre March 25, 2019 ########
###################################

################
## BEKK       ##
################
## Forecasting

standat <- list(T = nrow(r2), rts = r2, sigma1 = sigma1, nt = ncol(r2), ahead = 1)#, X=x)
standat

## Cf. https://groups.google.com/forum/#!topic/stan-users/tWQdtndbSnA for failed initalization: init_r < 2
bekk_fit <- sampling(bekk_mod, data = standat, verbose = TRUE, iter = 1000, control = list(adapt_delta = .9), init_r = .1, chains = 4)
dcc_fit <-  sampling(dcc_mod,  data = standat, verbose = TRUE, iter = 1000, control = list(adapt_delta = .9), init_r = 1, chains = 4)

print(bekk_fit, pars = c('Const', 'A', 'B', 'b0', 'b1', 'b2', 'lp__'),
      probs = c(0.025, 0.975))

print(dcc_fit, pars = c('c_h', 'a_h', 'b_h', 'Rwt', 'S', 'b0', 'b1', 'b2', 'lp__'),
      probs = c(0.025, 0.975))


rstan::get_elapsed_time(bekk_fit)
rstan::get_elapsed_time(dcc_fit)

hist(extract(bekk_fit)$Const[,1,1])
hdi(extract(bekk_fit)$Const[,1,1])
hist(extract(bekk_fit)$Const[,2,1])
hdi(extract(bekk_fit)$Const[,2,1])
hist(extract(bekk_fit)$Const[,2,2])
hdi(extract(bekk_fit)$Const[,2,2])


library(loo)
log_lik_bekk <- extract_log_lik(bekk_fit)
log_lik_dcc <- extract_log_lik(dcc_fit)

loo_bekk <- loo(log_lik_bekk)
loo_dcc <- loo(log_lik_dcc)
waic_bekk <- waic(log_lik_bekk)
waic_dcc <- waic(log_lik_dcc)

print(loo_bekk)
print(loo_dcc)
print(waic_bekk)
print(waic_dcc)

compare(loo_bekk, loo_dcc)
compare(waic_bekk, waic_dcc)


ar1 <- c(mean(r2[,1]), r2[,1][-1]) - mean(r2[,1])

plot(r2[,1], type = 'l')
abline(h=mean(r2[,1]))

summary(lm(r2[,1] ~ 1 + ar1 + x))

sqrt(mean(colMeans(extract(bekk_fit)[['H']][,,1,1])))

library(xtable)


rownames(summary(bekk_fit, pars = c('Const', 'A', 'B', 'corC[1,2]', 'b0', 'b1', 'b2', 'lp__'), probs = c(0.05, 0.95))$summary)



cor(colMeans(extract(bekk_fit)[['corH']][,,2,1]), x)


dim(extract(bekk_fit)[['corH']])
covmatsize <- dim(extract(bekk_fit)[['corH']])[4]
covmatsize

#pdf(file =  paste0(paste0('../FIGURES/mcorr', subs[subj]),'.pdf'))
op <- par(mfcol = c(3,1))
for(s in 1:(covmatsize-1)){
    kl = s + 1
    for(m in kl:covmatsize){
        #if(s == 2 & m == 2) next
        corH <- extract(bekk_fit)[['corH']][,,m,s]
        corH
        dim(corH)
        plot(corH[1,], type = 'l', ylim = c(-1, 1), col = "#81939920", ylab = "Correlation in H", xlab = "Day")
        for(i in sample(nrow(corH), size = 500)){
            lines(corH[i, ], col = '#81939920')
        }
        lines(colMeans(corH[,]), col = 'red', lwd = 2)
        for(k in 1:ncol(corH)){
            points(k, sort(corH[,k])[0.025*nrow(corH)], col = 'blue', pch = '*', lwd = 4)
            points(k, sort(corH[,k])[0.975*nrow(corH)], col = 'blue', pch = '*', lwd = 4)
        }
        lower <- array(NA, ncol(corH))
        upper <- array(NA, ncol(corH))
        for(k in 1:ncol(corH)){
            lower[k] <- sort(corH[,k])[0.025*nrow(corH)]
            upper[k] <- sort(corH[,k])[0.975*nrow(corH)]
        }
        lines(lower, col = 'blue', lwd = 1)
        lines(upper, col = 'blue', lwd = 1)
        abline(h=0)
        lines(x)
    }
}


#dev.off()

}


##dev.off()
lines((as.numeric(r2[,2])+3.5)/2, type = 'l', ylim = c(-3, 3), xlim = c(0, (nrow(r2) + ahead)), col = '#00b0e7')
lines((as.numeric(r2[,1])+3.5)/2, type = 'l', ylim = c(-3, 3), xlim = c(0, (nrow(r2) + ahead)))


lines(scale(na[,1])*0.1-0.4)
lines(scale(na[,2])*0.1-0.4, col = 'red')


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

r

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



########### cp uso
## forecasted params:
ahead <- standat$ahead

dim(extract(bekk_fit)[['rts_p']])
## ts 1
r2[,1]
p1 <- extract(bekk_fit)[['rts_p']][,2:(ahead+1),1]
p1
## ts 2
p2 <- extract(bekk_fit)[['rts_p']][,2:(ahead+1),2]
p3 <- extract(bekk_fit)[['rts_p']][,2:(ahead+1),3]


## orig
ts1 <- r2[,1]*attributes(r2[,1])$'scaled:scale' #+ attributes(r2[,1])$'scaled:center'
ts2 <- r2[,2]*attributes(r2[,2])$'scaled:scale' #+ attributes(r2[,2])$'scaled:center'
ts3 <- r2[,3]*attributes(r2[,3])$'scaled:scale' #+ attributes(r2[,3])$'scaled:center'
pr1 <- p1*attributes(r2[,1])$'scaled:scale' #+ attributes(r2[,1])$'scaled:center'
pr2 <- p2*attributes(r2[,2])$'scaled:scale' #+ attributes(r2[,2])$'scaled:center'
pr3 <- p3*attributes(r2[,3])$'scaled:scale' #+ attributes(r2[,3])$'scaled:center'

op <- par(mfcol = c(3, 1))
plot(ts1, type = 'l', xlim = c(0, (nrow(r2) + ahead)))
CI.Lpct <- round(nrow(pr1)*0.025, 0)  ## 95% CI position in vector
CI.Upct <- round(nrow(pr1)*0.9755, 0)
for( i in 1:ahead){
    points((nrow(r2) + i), mean(pr1[,i]), col = 'gray40')
    CI.L <- sort(pr1[,i])[CI.Lpct]
    CI.U <- sort(pr1[,i])[CI.Upct]
    lines(rep((nrow(r2) + i),2), c(CI.L, CI.U), col =  'gray70')
}
#
plot(ts2, type = 'l',  xlim = c(0, (nrow(r2) + ahead)), col = '#00b0e7')
for( i in 1:ahead ){
    points((nrow(r2) + i + .5), mean(pr2[,i]), col = '#004bbc')
    CI.L <- sort(pr2[,i])[CI.Lpct]
    CI.U <- sort(pr2[,i])[CI.Upct]
    lines(rep((nrow(r2) + i + .5),2), c(CI.L, CI.U), col =  '#add8e6')
}
#
##op
plot(ts3, type = 'l', xlim = c(0, (nrow(r2) + ahead)), col = '#cc0000')
for( i in 1:ahead ){
    points((nrow(r2) + i + .5), mean(p3[,i]), col = '#cc0000')
    CI.L <- sort(pr3[,i])[CI.Lpct]
    CI.U <- sort(pr3[,i])[CI.Upct]
    lines(rep((nrow(r2) + i + .5),2), c(CI.L, CI.U), col =  '#cc000040')
}
op


##pdf(file = '../FIGURES/CorrInH.pdf', width = 5, height = 3)

dim(extract(bekk_fit)[['corH']])


nt = ncol(r2)





######################33
## Check with rmgarch
## library('mgarchBEKK')
## simu <- simulateBEKK(20, 200, order = c(1,1))
## simu
## Dat <- matrix(unlist(simu$eps), ncol = 20)

library('rmgarch')

dim(r2)
Dat <- rbind(r2, r2)  # add two entries to make go up to 100
Dat

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
