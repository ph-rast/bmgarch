## Project title: 
##    Created at: 
##        Author: Philippe Rast
##          Data: 
##       Summary: Reaction to Reviewer 2 with GARCH
## ---------------------------------------------------------------------- ##
rm(list=ls())


#############################
## stan
#############################
library('rstan')
library('shinystan')
## enable multicore computing
rstan_options(auto_write = TRUE) ## cf. https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
options(mc.cores = parallel::detectCores())


## compare divergence of H matrix with KL
## Define KL divergence
K_L = function(Theta,hatTheta, Mu, hatMu){
    p = ncol(Theta)   
    invTheta = solve(Theta,diag(1,p))
    kl  = 0.5 * ( sum(diag(invTheta%*%hatTheta)) - log(det(invTheta%*%hatTheta)) - p +
                  t((hatMu - Mu)) %*% solve(hatTheta) %*% (hatMu - Mu) )
    return(kl)
}



####################
## Simulate data ##
###################
library(MASS)
library(Matrix)

tslength = 750
n_ts = 2

## Intercept
b0 = array(c(0, 0, 0), dim = c(n_ts,1))
## AR(1) process
b1 = array(0.3, dim = c(n_ts,1))


##########
## BEKK ##
##########

simuMGARCH <- function(process = 'VECH', tslength, n_ts, b0, b1){
    if(process == 'BEKK') {
        ## GARCH parameters
        ## Constant Variance
        Cnst <- matrix(c(.3, 0.10, 0.001,
                         0.10, .3,   0.001,
                         0.001, 0.001, 0.01), ncol = 3)
        Cnst <- Cnst[1:2, 1:2]
        ## A matrix
        A <- matrix(c(0.1, 0.0,   0.0,
                      0.3, 0.1,  0.0,
                      0.0, 0.0, 0.01), ncol = 3)
        A <- A[1:2, 1:2]
        ## B matrix
        B <- matrix(c(0.3, -.1,   0.,
                      0, 0.5, .0,
                      0   , 0.0, 0.95), ncol = 3)
        B <- B[1:2, 1:2]

        ## Define parameters 
        ## Mean vector
        BEKK_mu <- array(NA, dim = c(n_ts, tslength))
        ## Conditional Covariance
        BEKK_H <- list(NA)
        ## Returns, or outcome
        r <- array(NA, dim = c(n_ts, tslength))

        ## Init mean model at t = 1
        BEKK_H[[1]] <- Cnst + diag(c(A[1,1], A[2,2])) + diag(c(B[1,1], B[2,2]))
        BEKK_mu[,1] <- b0 + b1*0
        r[,1] <- mvrnorm(mu = BEKK_mu[,1], Sigma = BEKK_H[[1]])

        ## Generat TS for t >= 2
        for(t in 2:tslength){
            BEKK_mu[,t] = b0 + b1 * (r[,t-1] - BEKK_mu[,t-1])
            rr = ( r[,t-1] - BEKK_mu[,t-1] )%*%t( r[,t-1] - BEKK_mu[,t-1] )
            BEKK_H[[t]] =  Cnst + t(A)%*%rr%*%A + t(B)%*%BEKK_H[[t-1]]%*%B 
            r[,t] <- mvrnorm(mu = BEKK_mu[,t], Sigma = BEKK_H[[t]])
        }

        BEKK_y = r
        return(list(ts_y=BEKK_y, ts_H=BEKK_H, ts_mu=BEKK_mu))
    } else if( process == 'VECH' ) {
##########
        ## VECH ##
##########


        ## GARCH parameters
        ## Constant Variance, define as Cnst and then VEC-H to Const
        Cnst <- matrix(c(.3,   NA, 
                         0.10, .3), ncol = 2, byrow = TRUE)
        Const = Cnst[lower.tri(Cnst, diag = TRUE)]

        ## A matrix
        A <- matrix(c(0.1,  0.0,  0.0,
                      0.3,  0.01,  0.0,
                      0.0,  0.0,  0.1), ncol = 3)
        A 
        ## B matrix
        B <- matrix(c(0.1, -.1,   0.,
                      0,   -.01,  .0,
                      0,   0.0, 0.2), ncol = 3)
        B 


        ## Define parameters 
        ## Mean vector
        VECH_mu <- array(NA, dim = c(n_ts, tslength))
        ## Conditional Covariance
        VECH_H <- list(NA)
        ## Returns, or outcome
        r <- array(NA, dim = c(n_ts, tslength))

        ## Init mean model at t = 1
        VECH_H[[1]][upper.tri(Cnst, diag = TRUE)]  <- Const
        VECH_H[[1]] <- forceSymmetric(matrix(VECH_H[[1]], ncol = n_ts))


        VECH_mu[,1] <-  b0 + b1*0
        r[,1] <-  mvrnorm(mu = VECH_mu[,1], Sigma = VECH_H[[1]])
        r

        ## Generate TS for t >= 2
        for(t in 2:tslength){
            ## Empty matrix for temp result
            H_tmp <- array(NA, dim = c(n_ts, n_ts))
            VECH_mu[,t] = b0 + b1 * (r[,t-1] - VECH_mu[,t-1])
            rr = ( r[,t-1] - VECH_mu[,t-1] )%*%t( r[,t-1] - VECH_mu[,t-1] )
            Vrr = rr[lower.tri(rr, diag = TRUE)]
            VECH_H_tmp = Const +  A%*%Vrr + B%*%VECH_H[[t-1]][lower.tri(VECH_H[[t-1]], diag = TRUE) ]
            H_tmp[upper.tri(H_tmp, diag = TRUE)] <- VECH_H_tmp
            VECH_H[[t]] = forceSymmetric(H_tmp)
            r[,t] <- mvrnorm(mu = VECH_mu[,t], Sigma = VECH_H[[t]])
        }

        VECH_y = r
        return(list(ts_y=VECH_y, ts_H=VECH_H, ts_mu=VECH_mu))
    } else {
##############
        ## DCC simu ##
##############

        ## variance
        c_h = c(0.2, 0.2)
        a_h = c(.2, .3)
        b_h = c(.2, .2)

        a_q = 0.3
        b_q = 0.4

        ## mean
        DCC_mu <- array(NA, dim = c(n_ts, tslength))
                                        #b0 = array(5, dim = c(n_ts,1))
                                        #b1 = array(0.3, dim = c(n_ts,1))

        y = array(b0, dim = c(n_ts, tslength))

        DCC_mu[,1] <- b0 + b1*0

        h <- array(c(0.1, 0.2), c(n_ts,tslength))
        h

        Q <- list(diag(n_ts))
        Qs <- list(diag(n_ts))
        R <- list(diag(n_ts))
        DCC_H <- list(diag(n_ts))


        S_sd <- diag(c(0.5, 0.8))
        S_cor  <- matrix(c(1, 0.2, 0.2, 1), ncol = 2)
        S <- S_sd %*% S_cor %*% S_sd
        S
        
        
        u <- array(c(0.1, 0.2), c(n_ts,tslength))


        for(t in 2:tslength){
            DCC_mu[,t] = b0 + b1 * (y[,t-1] - DCC_mu[,t-1])
            for(i in 1:n_ts){
                h[i,t] = sqrt(c_h[i] + a_h[i]*(y[i, t-1] - DCC_mu[i, t-1])^2 + b_h[i]*h[i,t-1])
            }
            u[,t] <- solve(diag(h[,t])) %*% (y[,t] - DCC_mu[,t-1])
            Q[[t]] <- (1 - a_q - b_q)*S+ a_q * (u[,t-1] %*% t(u[,t-1])) + b_q * Q[[t-1]]
            Qs[[t]] <- solve(diag(sqrt(diag(Q[[t]]))))
            R[[t]]= Qs[[t]] %*% Q[[t]] %*% Qs[[t]]
            DCC_H[[t]] = diag(h[,t])%*%R[[t]] %*%diag(h[,t])
            ##
            y[,t] <- mvrnorm(mu = DCC_mu[,t], Sigma = DCC_H[[t]])
        }

        DCC_y=y
        return(list(ts_y=DCC_y, ts_H=DCC_H, ts_mu=DCC_mu))
    }
}


data <- simuMGARCH(process = 'VECH', tslength = 100, n_ts = 2, b0 = c(3, 0), b1 = c(0.3, 0.3))
data

#####
r = data$ts_y

op <- par(mfcol = c(2,1))
plot(r[1,], type = 'l', ylim = c(-3, 6))
lines(r[2,], type = 'l', col = 'red')
abline(h = 0, col = "#81939980")
abline(h = 0, col = "#cc000050")

                                        #
series1 <- array(NA, dim = c(n_ts, 1))
series2 <- series1
Coseries <- series1
for(i in 1:tslength){
    series1[i] <- data$ts_H[[i]][1,1]
    series2[i] <- data$ts_H[[i]][2,2]
    Coseries[i] <- cov2cor(data$ts_H[[i]])[1,2]
}
plot(series1, lty = 2, col = "#819399", type ='l', ylim=c(0, max(series1)))
lines(series2, lty = 2, col = "#cc0000")
lines(Coseries, lty = 2, col = "#33cc33")
op


r2 <- t(r)
mean(r2[,1])

sigma1 <- cov(r2)
r2
##dev.off()



###############################
## Load Compiled Stan Models ##
###############################


###########
## VECH  ##
###########

load(file = 'VECH.rds')
                                        #
#vech_mod <- stan_model(file = 'VECH-MGARCH.stan', verbose = TRUE)
#save(vech_mod, file = 'VECH.rds')


#########
## DCC ##
#########
dcc_mod <- readRDS(file = 'DCC-BivariateMGARCH.rds')

dcc_mod <- stan_model(file = './DCC-BivariateMGARCH.stan', verbose = TRUE)
saveRDS(dcc_mod, file = 'DCC-BivariateMGARCH.rds')

##########
## CDCC ##
##########
## have not defined a specifc bivariate model with crossed AR(1) in b2
cdcc_mod <- readRDS(file = 'CDCC-MGARCH.rds')

cdcc_mod <- stan_model(file = './CDCC-MGARCH.stan', verbose = TRUE)
saveRDS(cdcc_mod, file = 'CDCC-MGARCH.rds')

###########
## BEKK  ##
##########
#bekk_mod <- readRDS(file = 'BEKK-BivariateGARCH.rds')

#bekk_mod <- stan_model(file = 'BEKK-BivariateGARCH.stan', verbose = TRUE)
#saveRDS(bekk_mod, file = 'BEKK_bivariate.rds')


bekk_mod <- stan_model(file = 'BEKK-MGARCH.stan', verbose = TRUE)
saveRDS(bekk_mod, file = 'BEKK_multivariate.rds')



runSim <- function(conditions){
    ## define big list KLiter
    KLiter <- list()
    for(rep in 1:conditions[2]){
        data <- simuMGARCH(process = 'BEKK', tslength = 200, n_ts = 2, b0 = c(0, 0), b1 = c(0.3, 0.3))    
        ## Data
        r2 = t(data$ts_y)
        standat <- list(T = nrow(r2), rts = r2, sigma1 = var(r2), nt = ncol(r2), ahead = 1)
        ## Cf. https://groups.google.com/forum/#!topic/stan-users/tWQdtndbSnA for failed initalization: init_r < 2
        dcc_fit  <- sampling(dcc_mod,  data = standat, verbose = TRUE, iter = 1000, control = list(adapt_delta = .90), init_r = 2, chains = 4)
        cdcc_fit <- sampling(cdcc_mod, data = standat, verbose = TRUE, iter = 1000, control = list(adapt_delta = .90), init_r = 2, chains = 4)
        bekk_fit <- sampling(bekk_mod, data = standat, verbose = TRUE, iter = 1000, control = list(adapt_delta = .90), init_r = 1, chains = 4)
        print(bekk_fit, pars = c('A', 'B', 'Cnst', 'b0', 'b1', 'b2'))
        pairs( bekk_fit, pars = c('A[1,1]', 'A[1,2]') )
        stan_dens( bekk_fit, pars =  'A', separate_chains = T)
        mu <- data$ts_mu
        H <- data$ts_H
        ##
        ## Define List for both model specs
        KLout <- list()
        ## compute KL's for all t's for dcc, bekk and cdcc and store them into list with 
        for(k in 1:3){
            mgarch_fit <- list(dcc_fit, cdcc_fit, bekk_fit)[[k]]
            ## construct H matrix (Hest) , estimated from DCC
            KL <- array(NA, dim = c(standat$T))
            ## extract medians (later could use full posterior)
            lh11 <- apply( extract(mgarch_fit)[['L_H']][,1:standat$T,1,1], 2 , median)
            lh21 <- apply( extract(mgarch_fit)[['L_H']][,1:standat$T,2,1], 2 , median)
            lh22 <- apply( extract(mgarch_fit)[['L_H']][,1:standat$T,2,2], 2 , median)
            mu1  <- apply( extract(mgarch_fit)[['mu']][, 1:standat$T, 1 ], 2, median)
            mu2  <- apply( extract(mgarch_fit)[['mu']][, 1:standat$T, 2 ], 2, median)
            for(i in 1:standat$T){
                L_H <- diag(2)
                L_H[1,1] <- lh11[i]
                L_H[2,1] <- lh21[i]
                L_H[2,2] <- lh22[i]
                Hest <- L_H%*%t(L_H)
                ## obtain KL with known H and mu from simulation and Hest and estMu from model
                ## Correction: H is from a MVN with mu =0, as it is for the residuals. Ie., mu's need to stay 0!
                KL[i] <- K_L(H[[i]], Hest, c(0, 0), c(mu1[i], mu2[i]))
            }
            KLout[[k]] = KL
        }
        names(KLout) <- c('DCC', 'CDCC', 'BEKK')
        KLiter[[rep]] = KLout        
    }
    print(KLiter)
}

## Run simulation
lengths <- 50#seq(50, 130, by = 20)
niter <- 20
conditions = expand.grid(lengths, niter)
conditions

system.time({
KLresults <- apply(conditions, 1, runSim)
})

KLresults

saveRDS(KLresults, file = 'KLresults_bekk_small_ns.Rds')
KLresults <- readRDS(file = 'KLresults_bekk_small_ns.Rds')


op <- par ( mfcol = c(2, 2) )
plot(density(KLresults[[1]][[1]]$DCC), ylim = c(0,10), xlim = c(0, 1), main = 'DCC')
for(i in 1:niter){
    lines(density(KLresults[[1]][[i]]$DCC))
}


plot(density(KLresults[[1]][[1]]$CDCC), ylim = c(0,10), xlim = c(0, 1), main = 'CDCC')
for(i in 1:niter){
  lines(density(KLresults[[1]][[i]]$CDCC))
}


plot(density(KLresults[[1]][[1]]$BEKK), ylim = c(0,10), xlim = c(0, 1), main = 'BEKK')
for(i in 1:niter){
    lines(density(KLresults[[1]][[i]]$BEKK))
}

op

tmp <- list()
tmp2 <- matrix(NA, nrow = niter, ncol = 3)

## get median without first element to account for
## guessing the first element 
medianc <- function(x){
    return(median(x[-1]))
}

for(i in 1:nrow(conditions)){
    for(j in 1:niter){
        tmp2[j,] <- unlist(lapply(KLresults[[i]][[j]], medianc))
        colnames(tmp2) <- c('DCC', 'CDCC', 'BEKK')
    }
    tmp[[i]] <- tmp2
}

tmp

results <- matrix(NA, ncol = 4, nrow = nrow(conditions))
for(i in 1:nrow(conditions)){
    results[i, 1:3] <- apply(tmp[[i]], 2, median)
    results[i, 4] <- lengths[i]
}

colnames(results) <- c('DCC', 'CDCC','BEKK', 'SampSize')
results

## DCC vs BEKK
plot(results[,4], results[,1], type = 'l')
lines(results[,4], results[,2], type = 'l', col = 'red')
lines(results[,4], results[,3], type = 'l', col = 'blue')


## Plot with uncertainty 
KLresults

KLresults[[1]]$DCC

## KLresults[["Sample Size"]][["nRep"]]$"Model Type"


boxplot()

mods = colnames(results)[-4]
x = NA##KLresults[[1]][[1]]$mods[1]
out = list()

for(j in 1:3){
  for(i in 1:niter){
    x = append(x, KLresults[[3]][[i]][j])
  }
  out[j] = list(unlist(x)[-1] )
  x = NA##KLresults[[1]][[1]]$mods[1]
}



## Create bnoxplots for different n conditions and Model types
n50dcc = data.frame(out[1][[1]], 'DCC')
names(n50dcc) <- c('KLD', 'Model')
n50cdcc = data.frame(out[2][[1]], 'CDCC')
names(n50cdcc) <- c('KLD', 'Model')
n50bekk = data.frame(out[3][[1]], 'BEKK')
names(n50bekk) <- c('KLD', 'Model')

n50 = rbind(n50dcc,
            n50cdcc,
            n50bekk)

ggplot(data = n50, aes(y = KLD, x = Model)) + geom_boxplot()



n50dcc = data.frame(out[1][[1]], 'DCC')
names(n50dcc) <- c('KLD', 'Model')
n50cdcc = data.frame(out[2][[1]], 'CDCC')
names(n50cdcc) <- c('KLD', 'Model')
n50bekk = data.frame(out[3][[1]], 'BEKK')
names(n50bekk) <- c('KLD', 'Model')

n50 = rbind(n50dcc,
            n50cdcc,
            n50bekk)



## https://stats.stackexchange.com/questions/14673/measures-of-similarity-or-distance-between-two-covariance-matrices


KL[-1]
median(KL[2:standat$T]) # 
plot(KL[-1], type = 'l')
plot(density(KL[-1]))


## compare KL divergence among models: https://en.wikipedia.org/wiki/Fisher_information_metric

## head(sort(summary(mgarch_fit)$summary[,'Rhat'],decreasing=TRUE,50))
## print(mgarch_fit, pars = c('c_h', 'a_h', 'b_h', 'Rwt', 'S', 'b0', 'b1','lp__'),  probs = c(0.05, 0.95))
## rstan::traceplot(mgarch_fit, pars = c('c_h', 'a_h', 'b_h', 'Rwt'), inc_warmup = T)


#################################################################################################################################

#################################################################################################################################
library('rmgarch')

Dat <- r2
xspec = ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(garchOrder = c(1,1), model = 'sGARCH'), distribution.model = 'norm')
uspec = multispec(replicate(4, xspec))
spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
spec1a = dccspec(uspec = uspec, dccOrder = c(1, 1), model='aDCC', distribution = 'mvnorm')

cl = makePSOCKcluster(4)
multf = multifit(uspec, Dat, cluster = cl)

system.time({fit1 = dccfit(spec1, data = Dat, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)})
print(fit1)           
stopCluster(cl)


#####################
## Siwei's fitbit  ##
####################
fitbit <- read.csv(file = '../../../WORK/-DATASETS/SOURCE/fit_bit_dat.csv')
names(fitbit)
head(fitbit)

fitbit[, 1:4]

## Bivariate mgarch for one individual. Ts are steps and some variable
fitsub <- fitbit[,c( 'record_id', 'day', 'steps', 'active', 'interested', 'excited', 'proud', 'couple_id', 'guilty', 'veryActiveMinutes', 'afraid', 'scared')]
fitsub

subs <- unique(fitsub$record_id)
subs

library(lme4)

fit  <- lmer(interested ~ day + (1 + day | record_id), data = fitbit, na.action = na.exclude)
summary(fit)

library(brms)
brm(bf(interested ~ day + (1 + day | record_id),
       sigma ~ 1 + steps  + (1 + steps | record_id)),
    data = fitbit, cores = 4, iter = 400)



fitbit$resid <- abs(residuals(fit))^2
head(fitbit)
str(fitbit)
cors <- round( cor(na.omit(fitbit[,c(8:28,48:60)])), 2)
cors

dim(cors)
cors[34,]

X <- fitbit[,c(8:28,48:60)]
names(X)

summary( lm(resid ~ steps + veryActiveMinutes + veryActiveDistance + totalDistance + sedentaryMinutes + moderatelyActiveDistance+loggedActivitiesDistance+lightlyActiveMinutes+lightlyActiveDistance+caloriesOut, data = fitbit) )

summary( lmer(resid ~ steps + sedentaryMinutes + loggedActivitiesDistance + lightlyActiveMinutes + (1 | record_id), data = fitbit) )



out <- fitsub[fitsub$record_id==subs[9], c( 'steps' ,'interested')]
## Set 0 steps to NA
out$steps[out$steps==0] <- NA
r2 <- na.omit(out)

dim(r2)


if(FALSE){ ## skip all code in here
cpl_id <- sort(unique(fitbit$couple_id))
## 15 couples
## select only couples
fitfull <- fitsub[fitsub$couple_id%in%cpl_id, ]
fitfull
#
r2 <- array(NA, dim = c(100,2))
r2
r2[,1] <-  fitfull[fitfull$record_id==232, 'interested']
r2[,2] <-  fitfull[fitfull$record_id==233, 'interested']
#
max(subs)
length(subs)
#
## select only those with XX days
maxday <- array(0, dim = c(length(subs), 1))
for(i in 1:length(subs) ){
    sel <- subs[i]
    maxday[i] <- max(fitsub[fitsub$record_id==sel, 'day'])
}
#
range(maxday)
## Keep those with 100
fullday_id <- subs[maxday==100]
#
## Subset ot only those with 100 days
fitfull <- fitsub[fitsub$record_id%in%fullday_id, ]
#
## get first two individuals
ids <- subs[1:2]
#
r2 <- array(NA, dim = c(100,2))
r2
r2[,1] <-  fitfull[fitfull$record_id==80, 'interested']
r2[,2] <-  fitfull[fitfull$record_id==81, 'interested']
#
r2 <- fitfull[, c('interested')]
r2[,1] <-  scale(r2[,1], scale = T)
r2[,2] <-  scale(r2[,2], scale = T)
r2[,3] <-  scale(r2[,3], scale = T)
r2
}


r2 <- na.omit(r2)
r2
## composite:
r2[,1] <- rowSums(r2[,1:3])/3

r2[,1] <-  scale(r2[,1], scale = T)
r2[,2] <-  scale(r2[,2], scale = T)
r2[,3] <-  scale(r2[,3], scale = T)
r2[,4] <-  scale(r2[,4], scale = T)

r2 <- scale(r2)


plot(r2[,1]/100, type = 'l')
lines(r2[,2], col = 'green')
lines(r2[,3], col = 'blue')
lines(r2[,4], col = 'gray')

r2 <- r2[,c(1,4)]

plot(rowSums(r2[,1:3])/3, type = 'l')

sigma1 <- matrix(c(var(r2[,1]), cov(r2[,1],r2[,4]),
                   cov(r2[,1],r2[,4]),var(r2[,4])), ncol = 2, byrow = T)
sigma1

sigma1 <- cov(r2)



#sigma1 <- diag(3)*c(var(r2[,1]),var(r2[,2]), var(r2[,3])) 
sigma1
cov2cor(sigma1)

chol(sigma1)

ncol(r2)
nrow(r2)

r2
ncol(r2)


## vech_mod <- stan_model(model_code = vech, verbose = TRUE)
##


##################
## VECH         ##
##################
##
standat <- list(T = nrow(r2), rts = r2, sigma1 = sigma1, nt = ncol(r2), d=ncol(r2)*(ncol(r2)+1)/2 )
vech_fit <- sampling(vech_mod, data = standat, verbose = TRUE, iter = 800, chains = 4)
## vech_fit
##

print(vech_fit, pars = c('Cnst', 'A', 'B', 'lp__'),  probs = c(0.05, 0.95))



################
## BEKK       ##
################
## Forecasting


length(subs)
N=194

## Initialize list
resout <- vector("list", N)

pb <- txtProgressBar(min = 0, max = N, style = 3)

for(i in 161:N){
    out <- fitsub[fitsub$record_id==subs[i], c( 'steps' ,'interested')]
    ## Set 0 steps to NA
    out$steps[out$steps==0] <- NA
    ## check that there are at least 50 observations
    if(dim(na.omit(out))[1] >= 50) {
    r2 <- scale(na.omit(out))
    sigma1 <- cov(r2)
    #
    standat <- list(T = nrow(r2), rts = r2, sigma1 = sigma1, nt = ncol(r2), ahead = 1)
    ## Cf. https://groups.google.com/forum/#!topic/stan-users/tWQdtndbSnA for failed initalization: init_r < 2
    bekk_fit <- sampling(bekk_mod, data = standat, verbose = TRUE, warmup = 250, iter = 500, control = list(adapt_delta = .90), init_r = 1, chains = 4)
    fit <- summary(bekk_fit, pars = c('Const', 'A', 'B', 'corC', 'b0', 'b1', 'b2', 'lp__'))$summary[,c('mean', '2.5%', '97.5%', 'Rhat')]
    ## is 0 excluded?
    sig <- ifelse(sign(fit[,2]) * sign(fit[,3]) == 1, 1, 0)
    sigest <- sig * fit[,1]
    ## params of interest
    resout[[i]][1] <- list( fit[,1] )
    resout[[i]][2] <- list(sig)
    resout[[i]][3] <- list( fit[, 4])
    resout[[i]][4] <- list( sigest )
    } else next
    # update progress bar
    setTxtProgressBar(pb, i)
}

resout
saveRDS(object = resout, file = 'StepsInterested.Rds')

## combine the list into a matrix with all 23 elements of each list
temp <- matrix( rapply(resout, function(x) head(x, 23)), ncol = 23, byrow = TRUE) 
temp

##  add id identifier
temp2 <- data.frame(id = rep(1:136, each = 4), temp)
names(temp2)  <- c('id', names(resout[[1]][1][[1]]))

temp2

## select only columns of interest c(2,6,710,11)
temp3 <- temp2[, c('id', "A[1,2]", "A[2,1]", "B[1,2]", "B[2,1]")]

## return only spill over effects - and if they were significant (4th line of each individual)
temp3[seq(4, dim(temp3)[1], by = 4), ]



print(bekk_fit, pars = c('corH'),  probs = c(0.05, 0.95))


dcc_mod <- readRDS(file = 'DCC-MGARCH.rds')

dcc_fit <-  sampling(dcc_mod, data = standat, verbose = TRUE, iter = 1000, control = list(adapt_delta = .90), init_r = 1, chains = 4)
cdcc_fit <-sampling(cdcc_mod, data = standat, verbose = TRUE, iter = 1000, control = list(adapt_delta = .90), init_r = 1, chains = 4)

print(bekk_fit, pars = c('Const', 'A', 'B', 'corC', 'b0', 'b1', 'b2', 'lp__'), probs = c(0.05, 0.95))
print(cdcc_fit, pars = c('c_h', 'a_h', 'b_h', 'a_q', 'b_q','S', 'b0', 'b1', 'lp__'), probs = c(0.05, 0.95))
print(dcc_fit,  pars = c('c_h', 'a_h', 'b_h', 'a_q', 'b_q', 'S', 'b0', 'b1', 'lp__'), probs = c(0.05, 0.95))

## traceplot(vech_fit, pars = c('lp__'), inc_warmup = T)

traceplot(bekk_fit, pars = c('Cnst','B[2,1]','lp__'), inc_warmup = T)
traceplot(dcc_fit, pars = c('a_q','lp__'), inc_warmup = T)
 
## forecasted params:
ahead <- standat$ahead
dim(extract(bekk_fit)[['H_p']])
plot(extract(bekk_fit)[['H_p']][,1,2,1]/sqrt(extract(bekk_fit)[['H_p']][,1,1,1]*extract(bekk_fit)[['H_p']][,1,2,2]))

dim(extract(bekk_fit)[['rts_p']])
## person 1
p1 <- extract(bekk_fit)[['rts_p']][,2:(ahead+1),1]
p1
## person 2
p2 <- extract(bekk_fit)[['rts_p']][,2:(ahead+1),2]
p3 <- extract(bekk_fit)[['rts_p']][,2:(ahead+1),3]


#op <- par(mfcol = c(2, 1))
as.numeric(r2[,1])
as.numeric(r2[,2])

plot(as.numeric(r2[,1]), type = 'l', ylim = c(-3, 8), xlim = c(0, (nrow(r2) + ahead)))
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
    points((nrow(r2) + i + .5), mean(p2[,i]), col = '#004bbc')
    CI.L <- sort(p2[,i])[CI.Lpct]
    CI.U <- sort(p2[,i])[CI.Upct]
    lines(rep((nrow(r2) + i + .5),2), c(CI.L, CI.U), col =  '#add8e6')
}
#
##op

lines(as.numeric(r2[,3]), type = 'l', ylim = c(-3, 8), xlim = c(0, (nrow(r2) + ahead)), col = '#01a0e3')
for( i in 1:ahead ){
    points((nrow(r2) + i + .5), mean(p3[,i]), col = '#004bbc')
    CI.L <- sort(p3[,i])[CI.Lpct]
    CI.U <- sort(p3[,i])[CI.Upct]
    lines(rep((nrow(r2) + i + .5),2), c(CI.L, CI.U), col =  '#add8e6')
}


obs <- length(r2[,1])
xact <- (obs+1):(obs+ahead)

lines(xact, as.numeric((actual[,1]-mn1)/sd1), type = 'l', lty = 2)
lines(xact, as.numeric((actual[,2]-mn2)/sd2), type = 'l', lty = 2, col = 'blue')

##pdf(file = '../FIGURES/CorrInH.pdf', width = 5, height = 3)

dim(extract(bekk_fit)[['corH']])




nt = ncol(r2)



#pdf(file = '../FIGURES/mcorr.pdf')
op <- par(mfcol = c(3,1))
for(s in 1:2){
    for(m in 2:3){
        if(s == 2 & m == 2) next
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
    }
}

dev.off()


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

#########################################
#########################################
H <- matrix(c(25, 5, 5, 20), ncol = 2)
cov2cor(H)

sH <- solve(H)
sH%*%H
sH



install.packages('ks')
library(ks)
x <- matrix(1:9, nrow=3, ncol=3)
y <- (x + t(x))/2
x
y
vech(x)


ks::vech


install.packages('expm')
library(expm)

V = matrix(c(1, 2, 3, 1, -1, 1, 2, 2, 1), ncol = 3)
V
D = diag(c(81, 9, 7))

## Create Matrix to be sqrt'ed
A <-V%*%D%*%solve(V)
A

A


## Solve by diagonalization
eigen(A)

V <- eigen(A)$vectors
V
D <- diag(eigen(A)$values)
D
S <- sqrt(D)
dim(S)
S

## Square root, R:
R <- V%*%S%*%solve(V)
R
R%*%R
A

## with package
sqrtm(A)



sqrtm

sqrtm(A)%*%sqrtm(A)

## Manually take square root of 2x2 matrix
tau <- sum(diag(A))
delta <- det(A)
s <- sqrt(delta)
s
t <-  sqrt(tau + 2*s)

1/t*(A + s*diag(2))


A = sigma1
1/sqrt(sum(diag(A)) + 2*sqrt(det(A)))*(A + sqrt(det(A))*diag(2))
## for stan, replace sum(diag(A)) with trace(A)


rev(r2)

r2


a  <-  array(c(4,7), dim = c(2,1))
a

sum(a) - a

r2[1,]

sum(r2[1,]) - r2[1,]
