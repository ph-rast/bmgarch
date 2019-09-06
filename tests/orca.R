#!/global/software/R-3.0.2/bin/Rscript
## Project title: CP with ORCA
##    Created at: 27-02-2015
##        Author: Philippe Rast
##          Data: orca
##       Summary:
## ---------------------------------------------------------------------- ##

rm(list=ls())

## load 'orca' object
load('~/UZH/Projekte/Hold/CPO_ChangePointInORCATECH/WORK/-DATASETS/DERIVED/orca.RDat')

head(orca)
names(orca)

###############################
## Data preparation for stan
###############################

## Step 1: Only converters
names(orca)
## SELECT ALL with <=
converters <- orca[orca$Donset<=1,c('oadc', 'day', 't2c.day', 'meanws',
                       'agebaseline.gc','onsetage.c', 'avgsess.cpu.time')]
## Remove NA from dependent
convert <-  converters[!is.na(converters$meanws),]
convert[1:900,]

id <- unique(convert$oadc)
length(id)
for(i in 1:length(id)){  ## enumarte from 1 thru N
    convert[convert$oadc==id[i], 'subject'] <- i
}

subject <- convert$subject
subject


names(convert)

## Dependent vairable: meanws
y.univ=convert$meanws
nobs=length(y.univ)

convert$t2c.day

## Rescale t2c.day to weekly scale
convert$t2c.week <- convert$t2c.day/7
convert$t2c.week <- convert$day/7

## Rescale t2c.day to monthly scale
convert$t2c.month <- convert$day/(365/12)

################################
## Aggregate data to two weeks.
## months are now increasing in 0.5 steps
################################
convert$t2c.month[1:20]
names(convert)

convert[1:100, c('oadc', 'meanws', 't2c.week')]
agr <- 1 ## number of weeks to aggregate
agr.vars <- c('meanws', 't2c.month','subject','agebaseline.gc', 'avgsess.cpu.time')

agr.orca <- matrix(NA, ncol = length(agr.vars))
agr.orca <- data.frame(agr.orca)
names(agr.orca) <- agr.vars
agr.orca

index=1

for(i in id){
    init <- convert[convert$oadc==i,'t2c.week'][1]
    end <- max(convert[convert$oadc==i,'t2c.week'])
    while(init<end){
        agr.orca[index,] <- colMeans(convert[convert$oadc==i &
                                                 convert$t2c.week>=init &
                                                     convert$t2c.week<init+agr,
                                             agr.vars], na.rm=T)
        ## keep month as scale
        init <- init+agr
        index <- index+1
    }
    index <- length(agr.orca[,1])+1
}

## omit NA's in dependent
agr.orca <- agr.orca[!is.na(agr.orca$meanws),]
agr.orca
names(agr.orca)

## IMPORTANT: Previously created objects are overwritten to match agr.orca
subject <- agr.orca$subject
subject
y.univ <- agr.orca$meanws
nobs=length(y.univ)

## write vectors that can be easily exchanged with different objects if necessary
t2c.month <- agr.orca$t2c.month
################################


## Construct BP baseline age vector
baselineage <-array(NA, dim = length(id))
baselineage

for(i in 1:length(id)){
    baselineage[i] <-  unique(convert[convert$oadc==id[i], c('agebaseline.gc')])
}
baselineage

## Construct BP baseline onset age vector
baselineonsetage <-array(NA, dim = length(id))
baselineonsetage

for(i in 1:length(id)){
    baselineonsetage[i] <- unique(convert[convert$oadc==id[i], 'onsetage.c'])}
baselineonsetage


## Note: orca is in long format
library(nlme)
mod.lme0 <- lme(meanws ~ 1 + agebaseline.gc+t2c.month,
                random=~1+t2c.month|subject,  data=agr.orca,
                na.action=na.omit)
summary(mod.lme0)
VarCorr(mod.lme0)


################################################################################
##        MODELING                                                            ##
################################################################################

test.dat = na.omit(agr.orca)

test.dat[test.dat$subject == 3, ]

plot(test.dat[test.dat$subject == 3, 't2c.month'], test.dat[test.dat$subject == 3, 'meanws'], type = 'l', ylim = c(0, 80) )
lines(test.dat[test.dat$subject == 3, 't2c.month'], test.dat[test.dat$subject == 3, 'avgsess.cpu.time'], type = 'l')


rt = test.dat[test.dat$subject == 3, c('meanws', 'avgsess.cpu.time')]
out = bmgarch(data = rt, parameterization = 'DCC')

rstan::traceplot(out$model_fit, pars = 'a_q', inc_warmup = T)

summary(out)

plot.bmgarch(object = out, type = 'ccor')

forecast( out, ahead = 3 )
