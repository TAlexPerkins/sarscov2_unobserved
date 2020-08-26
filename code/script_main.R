#=============================================================================#
# Authors: Alex Perkins, Sean Cavany, Sean Moore, Rachel Oidtman, Anita Lerch, Marya Poterek
# project: Estimating unobserved SARS-CoV-2 infections in the United States
# Year: 2020
# 
# Code to generate all figures and results from main text
#
#=============================================================================#
# set up workspace
#=============================================================================#

# load libraries
library(extraDistr)
library(doParallel)
library(mc2d)
library(MASS)
library(boot)

# load function to simulate autochthonous transmission
source('simOutbreak.R')

# set random number seed
set.seed(1234)

#=============================================================================#
# load in and process data
#=============================================================================#

# read in line list data for US
# updated 20200312
# data from https://github.com/midas-network/COVID-19/tree/master/data/cases/global/line_listings_nihfogarty
linelist = read.csv('../data/2020_03_12_1800EST_linelist_NIHFogarty.csv')
yesUS = subset(linelist, country=='USA')
# remove Diamond Princess repatriated cases
yesUS = yesUS[grep("Diamond",yesUS$summary,invert=T),]

# fit gamma parameters for symptom to report delay
data.delay = as.Date(yesUS$reporting.date) - as.Date(yesUS$symptom_onset)
data.delay = as.numeric(data.delay[which(!is.na(data.delay))])
delay.shape.baseline = MASS::fitdistr(data.delay,dgamma,start=list(shape=0.5,rate=0.5))$estimate[1]
delay.rate = MASS::fitdistr(data.delay,dgamma,start=list(shape=0.5,rate=0.5))$estimate[2]

# number of travelers that were cases or died
num.CF = c(
  nrow(subset(yesUS,international_traveler>0))-sum(subset(yesUS,international_traveler>0)$death>0,
                                                   na.rm=T),
  sum(subset(yesUS,international_traveler>0)$death>0,na.rm=T))

# read in case data internationally
# updated 20200307
# data from https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series
ts = read.csv('../data/time_series_19-covid-Confirmed.csv')
ts.natl = matrix(0,length(unique(ts$Country.Region)),ncol(ts)-4)
for(ii in 1:ncol(ts.natl)){
  ts.natl[,ii] = aggregate(ts[,4+ii],by=list(ts$Country.Region),FUN=sum)[,2]
}
row.names(ts.natl) = aggregate(ts[,4+ii],by=list(ts$Country.Region),FUN=sum)[,1]
for(ii in 1:nrow(ts.natl)){
  ts.natl[ii,-1] = pmax(0,diff(ts.natl[ii,]))
}

# correct for travel ban from China (for non-US citizens starting 2/2 at 5pm) - so 0 out starting 2/3
colnames(ts.natl) = 22:(ncol(ts.natl)+21)
which(colnames(ts.natl)==34)
ts.natl['China',which(colnames(ts.natl)==34):ncol(ts.natl)] = 0

# read in death data internationally
# updated 20200307
# data from https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series
tsd = read.csv('../data/time_series_19-covid-Deaths.csv')
tsd.natl = matrix(0,length(unique(tsd$Country.Region)),ncol(tsd)-4)
for(ii in 1:ncol(ts.natl)){
  tsd.natl[,ii] = aggregate(tsd[,4+ii],by=list(tsd$Country.Region),FUN=sum)[,2]
}
row.names(tsd.natl) = aggregate(tsd[,4+ii],by=list(tsd$Country.Region),FUN=sum)[,1]
for(ii in 1:nrow(tsd.natl)){
  tsd.natl[ii,-1] = pmax(0,diff(tsd.natl[ii,]))
}
colnames(tsd.natl) = 22:(ncol(tsd.natl)+21)

# count up local cases by day in the US
cases.US.total = c(rep(0,21),ts.natl['US',])
cases.US.imported =
  table(
    as.Date(as.character(subset(yesUS,international_traveler>0)$reporting.date)) -
    as.Date('2019-12-31'))
tmp = rep(0,length(cases.US.total))
tmp[as.numeric(names(cases.US.imported))] = cases.US.imported
cases.US.imported = tmp
rm(tmp)
cases.US.local = pmax(0, cases.US.total - cases.US.imported)

# count up local deaths by day in the US
deaths.US.total = c(rep(0,21),tsd.natl['US',])
deaths.US.imported =
  table(
    as.Date(as.character(subset(yesUS,international_traveler>0&death>0)$reporting.date)) -
      as.Date('2019-12-31'))
tmp = rep(0,length(deaths.US.total))
tmp[as.numeric(names(deaths.US.imported))] = deaths.US.imported
deaths.US.imported = tmp
rm(tmp)
deaths.US.local = pmax(0, deaths.US.total - deaths.US.imported)

#=============================================================================#
# simulate imported infections
#=============================================================================#
parameters = read.csv("../data/parameters.csv",sep=',',header=T,row.names=1)
# sample replicates of how many infections have been imported into the US
maxUS = 2e4
rangeUS = sum(yesUS$international_traveler>0,na.rm=T):maxUS
# estimate for asymptomatic proportion based on
# https://www.medrxiv.org/content/10.1101/2020.02.20.20025866v2
PrAsymptomatic = exp(optim(par=c(0,0),fn=function(par){
    sum((
        c(exp(par[1])/(exp(par[1])+exp(par[2])),
          qbeta(c(0.025,0.975),exp(par[1]),exp(par[2]))) -
        parameters["asymptomatic",1:3]) ^ 2)})$par)

# distribution parameters for death
PrDeathSymptom = c(parameters["cfr",1],
                   diff(as.numeric(parameters["cfr",2:3]))/1.96/2)

# set values of unknown parameters
# note that these values seem to maximize the probability of the cumulative
# deaths in the US as of March 12, 2020 predicted by the model
replicates = 1000
load("../results/sensitivity/param_estimates_posterior_1.rda", verbose=T)
indices = sample(1:length(PrCaseSymptom.trav_posterior), replicates, replace=TRUE)
PrCaseSymptom.trav = PrCaseSymptom.trav_posterior[indices]
asympRFraction = asympRFraction_posterior[indices]

# sample from uncertainty about proportions of infection outcomes
propns.ASCF = cbind(
  rbeta(replicates,PrAsymptomatic[1],PrAsymptomatic[2]),
  rnorm(replicates,PrDeathSymptom[1],PrDeathSymptom[2]))
propns.ASCF[propns.ASCF<0]=0
propns.ASCF = cbind(
  propns.ASCF[,1],
  (1-propns.ASCF[,1]) * (1-PrCaseSymptom.trav) * (1-propns.ASCF[,2]),
  (1-propns.ASCF[,1]) * PrCaseSymptom.trav * (1-propns.ASCF[,2]),
  (1-propns.ASCF[,1]) * propns.ASCF[,2])

# sample from uncertainty of parameters that we have uncertainty ranges of
R.reps = rnorm(replicates,parameters["R",1],
               diff(as.numeric(parameters["R",2:3]))/3.92)
R.reps[R.reps<0]=0

if (!any(is.na(parameters["incubation_shape",2:3]))) {
    ip.shape.reps = rnorm(replicates, parameters["incubation_shape",1],
                          diff(as.numeric(parameters["incubation_shape",2:3]))/3.92)
    ip.shape.reps[ip.shape.reps<0]=1e-10
} else {
    ip.shape.reps = rep(parameters["incubation_shape",1], replicates)
}
if (!any(is.na(parameters["incubation_scale",2:3]))) {
    ip.scale.reps = rnorm(replicates, parameters["incubation_scale",1],
                          diff(as.numeric(parameters["incubation_scale",2:3]))/3.92)
    ip.scale.reps[ip.shape.reps<0]=1e-10
} else {
    ip.scale.reps = rep(parameters["incubation_scale",1], replicates)
}

k.mean = parameters["k",1]
k.lower = parameters["k",2]
k.upper = parameters["k",3]
if (!(is.na(k.lower) | is.na(k.upper))) {
    k.meanlogs = seq(-10,log(k.mean),0.01)
    k.sdlogs = sqrt(2*(log(k.mean) - k.meanlogs))
    ind = which.min(sqrt((qlnorm(0.025,k.meanlogs,k.sdlogs)-k.lower)^2
                         + (qlnorm(0.975,k.meanlogs,k.sdlogs)-k.upper)^2))
    k.reps = rlnorm(replicates,k.meanlogs[ind],k.sdlogs[ind])
} else {
    k.reps = rep(k.mean,replicates)
}

if (!any(is.na(parameters["serial_mean",2:3]))) {
    serial.mean.reps = rnorm(replicates, parameters["serial_mean",1],
                          diff(as.numeric(parameters["serial_mean",2:3]))/3.92)
    serial.mean.reps[ip.shape.reps<0]=1e-10
} else {
    serial.mean.reps = rep(parameters["serial_mean",1], replicates)
}

if (!any(is.na(parameters["serial_sd",2:3]))) {
    serial.sd.reps = rnorm(replicates, parameters["serial_sd",1],
                          diff(as.numeric(parameters["serial_sd",2:3]))/3.92)
    serial.sd.reps[ip.shape.reps<0]=1e-10
} else {
    serial.sd.reps = rep(parameters["serial_sd",1], replicates)
}

if (!any(is.na(parameters["onset_death_mean",2:3]))) {
    onset_death.mean.reps = rnorm(replicates, parameters["onset_death_mean",1],
                          diff(as.numeric(parameters["onset_death_mean",2:3]))/3.92)
    onset_death.mean.reps[ip.shape.reps<0]=1e-10
} else {
    onset_death.mean.reps = rep(parameters["onset_death_mean",1], replicates)
}

if (!any(is.na(parameters["onset_death_sd",2:3]))) {
    onset_death.sd.reps = rnorm(replicates, parameters["onset_death_sd",1],
                          diff(as.numeric(parameters["onset_death_sd",2:3]))/3.92)
    onset_death.sd.reps[ip.shape.reps<0]=1e-10
} else {
    onset_death.sd.reps = rep(parameters["onset_death_sd",1], replicates)
}

# draw samples of the number of imported infections
imports = numeric(length=replicates)
for(ii in 1:replicates){
  PrImportedInfections =
    dmultinomial(
      x = cbind(
        0:(maxUS-sum(num.CF)),
        num.CF[1],num.CF[2]),
      prob = c(sum(propns.ASCF[ii,1:2]),propns.ASCF[ii,3:4]))
  imports[ii] =
    sample(
      sum(num.CF):maxUS,
      1,
      prob=PrImportedInfections,
      replace=T)
}

# draw samples of the number of local infections with simple extrapolation
num.CF.local = c(sum(cases.US.local),sum(deaths.US.local))
local.simple = rep(NA,replicates)
for(ii in 1:replicates){
  PrLocalSimpleInfections =
    dmultinomial(
      x = cbind(
        0:(maxUS-sum(num.CF.local)),
        num.CF.local[1],num.CF.local[2]),
      prob = c(sum(propns.ASCF[ii,1:2]),propns.ASCF[ii,3:4]))
  if(sum(PrLocalSimpleInfections) > 0){
    local.simple[ii] =
      sample(
        sum(num.CF.local):maxUS,
        1,
        prob=PrLocalSimpleInfections,
        replace=T)
  }
}

# draw samples of the day on which imported infections arrived
case.days=vector()
for(i in 1:length(cases.US.imported)){
  if(cases.US.imported[i]>0){
    if(length(case.days)==0){
      case.days=rep(i,cases.US.imported[i])
    }else{
      case.days=c(case.days,rep(i,cases.US.imported[i]))      
    }
  }
}
import.case.density = density(
  case.days,
  from = 1,
  to = length(cases.US.imported),
  n = length(cases.US.imported))$y

# estimate the day of the year on which imports occur
import.doy = list()
for(ii in 1:replicates){
  import.doy[[ii]] = sample(
    1:length(cases.US.imported),
    imports[ii],
    prob=import.case.density,
    replace=T)
}

#=============================================================================#
# simulate local transmission
#=============================================================================#

# simulate local transmission for each draw of imported infections
local = foreach(ii = 1:replicates) %do% {
  simOutbreak(
    timeImport = import.doy[[ii]], # timing of each imported infection
    R = R.reps[ii], # reproduction number
    k = k.reps[ii], # dispersion parameter
    si_mean = serial.mean.reps[ii], # mean of serial interval distribution
    si_sd = serial.sd.reps[ii], # standard deviation of serial interval distribution
    inc_shape = ip.shape.reps[ii], # shape parameter of incubation period distribution
    inc_scale = ip.scale.reps[ii], # scale parameter of incubation period distribution
    symp_to_death_mean = onset_death.mean.reps[ii], # mean of time between symptom onset and death
    symp_to_death_sd = onset_death.sd.reps[ii], # std. dev. of time between symptom onset and death
    report_delay_shape = delay.shape.baseline, # shape parameter for delay between symptom and reporting
    report_delay_rate = delay.rate, # rate parameter for delay between symptom and reporting
    stopSimulationDay = length(cases.US.imported), # day of year since Jan 1 when simulation stops
    asympProp = propns.ASCF[ii,1], # proportion of infections that are asymptomatic
    asympRFraction = asympRFraction[ii], # relative infectiousness of asymptomatics
    lnormFlag = T # toggles whether serial interval distribution is lognormal
  )
}

#simulate deaths out but turn transmission off on 12 March
local.predict = foreach(ii = 1:replicates) %do% {
  simOutbreakR0Change(
    timeImport = import.doy[[ii]], # timing of each imported infection
    R = R.reps[ii], # reproduction number
    k = k.reps[ii], # dispersion parameter
    si_mean = serial.mean.reps[ii], # mean of serial interval distribution
    si_sd = serial.sd.reps[ii], # standard deviation of serial interval distribution
    inc_shape = ip.shape.reps[ii], # shape parameter of incubation period distribution
    inc_scale = ip.scale.reps[ii], # scale parameter of incubation period distribution
    symp_to_death_mean = onset_death.mean.reps[ii], # mean of time between symptom onset and death
    symp_to_death_sd = onset_death.sd.reps[ii], # std. dev. of time between symptom onset and death
    report_delay_shape = delay.shape.baseline, # shape parameter for delay between symptom and reporting
    report_delay_rate = delay.rate, # rate parameter for delay between symptom and reporting
    stopSimulationDay = 180, # day of year since Jan 1 when simulation stops
    asympProp = propns.ASCF[ii,1], # proportion of infections that are asymptomatic
    asympRFraction = asympRFraction[ii], # relative infectiousness of asymptomatics
    lnormFlag = T, # toggles whether serial interval distribution is lognormal
    RChangeDay = length(cases.US.imported), # determines when R changes
    RChange = 0 # determines what R drops to at R0ChangeDay
  )
}

# load the following to generate the objects used to generate the figures in the paper
load("../results/baseline_projections.rda",verbose=T)

#=============================================================================#
# produce plots and results for all main text figures 
#=============================================================================#

# Infections and imports results - processing
local.mat = t(matrix(
  unlist(lapply(local, function(x) x$daily)),
  length(local[[1]]$daily),
  replicates))

# Infections and imports results - quantities
quantile(PrCaseSymptom.trav_posterior,c(0.025, 0.5, 0.975))
quantile(unlist(lapply(local,function(ll)ll$cum)),c(0.025,0.5,0.975))
quantile(local.mat[,ncol(local.mat)],c(0.025,0.5,0.975))
sum(sum(cases.US.local)<unlist(lapply(local,function(ll)ll$cum)))/nrow(local.mat)

# Figure 1
# plot distribution of cumulative infections
## plot all locally acquired infections over time
wdth = 9.5
hght = 5.5
pdf('../plots/figure_1_cumulative_infections_and_infections_daily.pdf',
    width=wdth/2.54,height=hght/2.54, pointsize=8)
par(mfrow=c(1,2),
    mar=c(4.6,4.1,2.1,0.6))
plot.temp = log(unlist(lapply(local,function(ll)ll$cum)),base=10)
hist(plot.temp,breaks=seq(2,10),plot=F)$density
plot.temp = plot.temp[plot.temp <=8]
hist(plot.temp,
     breaks=seq(2,8),col='gray',
     xlab=expression('Cumulative infections (log'[10]*')'),
     ylab='Proportion of simulations',main='',las=1,xaxt="n",freq=F)
axis(1,at=seq(2,8,2),labels=c(expression(10^2),expression(10^4),
                              expression(10^6),expression(10^8)))
lines(log(rep(sum(cases.US.local),2),base=10),
      c(0,replicates),
      col='red',lwd=2)
mtext("A",side=3,line=0, 
       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
       cex=1.2)
plot(
  as.Date('2019-12-31') + 31:ncol(local.mat),
  apply(local.mat,2,function(ii)median(ii,na.rm=T))[31:ncol(local.mat)],
  ylim=c(1,quantile(local.mat[,ncol(local.mat)],0.975)),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
  xlim=as.Date('2019-12-31') + c(31,ncol(local.mat)),
  xlab='Date',ylab='',main='', log="y",xaxt="n")
mtext("Infections",side=2,line=3.5)
axis(1,at=c(as.Date("2020-02-01"), as.Date("2020-02-15"),as.Date("2020-03-01")),
     labels=c("Feb 1","Feb 15","Mar 1"))
## for (ii in 1:nrow(local.mat)) {
##     lines(
##         as.Date('2019-12-31') + 31:ncol(local.mat),
##         local.mat[ii,31:ncol(local.mat)],
##         col=rgb(0,0,0,0.1),lwd=0.25,type='l')
## }
polygon(
  c(as.Date('2019-12-31') + 1:ncol(local.mat),
    rev(as.Date('2019-12-31') + 1:ncol(local.mat))),
  c(apply(local.mat,2,function(ii)quantile(ii,0.25,na.rm=T)),
    rev(apply(local.mat,2,function(ii)quantile(ii,0.75,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.25))
polygon(
  c(as.Date('2019-12-31') + 1:ncol(local.mat),
    rev(as.Date('2019-12-31') + 1:ncol(local.mat))),
  c(apply(local.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
    rev(apply(local.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.15))
mtext("B",side=3,line=0, 
       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
       cex=1.2)
dev.off()

# Cases results - processing
updateDaily = FALSE # turn on Bayesian daily updating
smoothSpline = TRUE # turn on smoothing spline
cases.mat = t(matrix(
  unlist(lapply(local, function(x) x$cases)),
  length(local[[1]]$cases),
  replicates))
p.mat = matrix(NA,nrow(cases.mat),ncol(cases.mat))
cases.mat.obs = rbinom(length(cases.mat), as.vector(cases.mat), rowSums(propns.ASCF[,2:4]))
cases.mat.obs = matrix(cases.mat.obs, replicates, ncol(cases.mat))
for(ii in 1:nrow(cases.mat)){
  alpha.old=1
  beta.old=1
  for(jj in 1:ncol(cases.mat)){
    if(cases.mat[ii,jj]){
      actual.cases = rbinom(1,cases.mat[ii,jj], sum(propns.ASCF[ii,2:4]))
      alpha.new = alpha.old+cases.US.local[jj]
      beta.new = beta.old+actual.cases-cases.US.local[jj]
      p.mat[ii,jj] =
          rbeta(1,alpha.new,max(1,beta.new))
      if (updateDaily) {
        alpha.old=alpha.new
        beta.old=beta.new
      }
    }
  }
  if (smoothSpline) {
    non.NA.indices = which(!is.na(p.mat[ii,]))
    if(length(non.NA.indices) > ncol(p.mat) / 3){
      temp.sp = smooth.spline((1:ncol(p.mat))[non.NA.indices],
                              logit(p.mat[ii,non.NA.indices]),
                              nknots=floor((ncol(p.mat) - non.NA.indices[1])/7 + 0.5))
      p.mat[ii,non.NA.indices[1]:ncol(p.mat)] = inv.logit(predict(temp.sp, non.NA.indices[1]:ncol(p.mat))$y)
    } else {
      p.mat[ii,non.NA.indices[1]:ncol(p.mat)] = NA
    }
  }
}
testing = read.csv("../data/testing_ts.csv")
incomplete.days = testing$Day[!testing$Complete]
low.p = which.min(apply(p.mat,2,function(ii)median(ii,na.rm=T)))

# Cases results - quantities
low.p+as.Date('2019-12-31')
apply(p.mat,2,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=T))[,low.p]
apply(p.mat,2,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=T))[,length(cases.US.local)]
quantile(rowSums(ifelse(p.mat < 0.2,1,0),na.rm=T), c(0.025,0.5,0.975))
quantile(rowSums(ifelse(p.mat < 0.1,1,0),na.rm=T), c(0.025,0.5,0.975))
# obtain correlations of p.mat and testing
totalTestsPerDay = c(rep(0, min(testing$Day)-1),
                     testing$Total[1:which(testing$Day==ncol(p.mat))])
start.day.cor = low.p
end.day.cor = max(testing$Day[testing$Complete])
correlation.reduced = cor(t(p.mat[,start.day.cor:end.day.cor]),
                          totalTestsPerDay[start.day.cor:end.day.cor], use="pairwise")
quantile(correlation.reduced,c(0.025,0.5,0.975),na.rm=T)

# Figure 2
# plot locally acquired symptomatic infections over time, alongside number of cases reported in the US
# plot proportion of locally acquired symptomatic infections reported over time
# alongside numbers of tests administered in the US
wdth = 11.5
hght = 10.5
pdf('../plots/figure_2_symptomatic_daily_and_symptomatic_detected.pdf',
    ##width=4.5,height=13.8, pointsize=14)
    width=wdth/2.54,height=hght/2.54, pointsize=8)
par(mfrow=c(2,2))
par(mar = c(4, 4, 2, 4) + 0.1,
    oma = c(0, 0, 0, 0))
plot(
  as.Date('2019-12-31') + 1:ncol(cases.mat),
  apply(cases.mat.obs,2,function(ii)median(ii,na.rm=T)),
  ylim=c(0,quantile(cases.mat.obs[,ncol(cases.mat)],0.75)),
  col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
  xlim=as.Date('2019-12-31') + c(31,ncol(cases.mat)),
  xlab='Date',ylab='Symptomatic infections',main='',xaxt="n")
axis(1,at=c(as.Date("2020-02-01"), as.Date("2020-02-15"),as.Date("2020-03-01")),
     labels=c("Feb 1","Feb 15","Mar 1"))
polygon(
  c(as.Date('2019-12-31') + 1:ncol(cases.mat),
    rev(as.Date('2019-12-31') + 1:ncol(cases.mat))),
  c(apply(cases.mat.obs,2,function(ii)quantile(ii,0.25,na.rm=T)),
    rev(apply(cases.mat.obs,2,function(ii)quantile(ii,0.75,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.25))
polygon(
  c(as.Date('2019-12-31') + 1:ncol(cases.mat),
    rev(as.Date('2019-12-31') + 1:ncol(cases.mat))),
  c(apply(cases.mat.obs,2,function(ii)quantile(ii,0.025,na.rm=T)),
    rev(apply(cases.mat.obs,2,function(ii)quantile(ii,0.975,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.15))
mtext("A",side=3,line=0, 
      at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.2)
lines(as.Date('2019-12-31') + 1:ncol(cases.mat),
      cases.US.local, col="red",lwd=2)
## par(new = TRUE)
## plot(as.Date('2019-12-31') + 1:ncol(cases.mat),
##      cases.US.local, type="l", col="red",
##      axes=F, bty = "n", xlab = "", ylab = "",
##      xlim=as.Date('2019-12-31') + c(31,ncol(cases.mat)), lwd=2,
##      xaxs='i',yaxs='i')
## axis(side=4, at = pretty(range(cases.US.local)), col="red", col.axis="red",las=1)
## mtext("Reported cases", side=4, line=3, lwd=2, col="red")        
legend("topleft", col=c("red", "black"), lty="solid",
       legend=rev(c("Model", "Data")),
       bty="n", lwd=c(2,2))
cases.mat.new = t(matrix(
  unlist(lapply(local.predict, function(x) x$cases[1:length(local[[1]]$cases)])),
  length(local[[1]]$cases),
  replicates))
det.cases.mat.obs = rbinom(length(cases.mat.new), as.vector(cases.mat.new), rowSums(propns.ASCF[,2:4])*p.mat)
det.cases.mat = matrix(det.cases.mat.obs, replicates, ncol(cases.mat.new))
plot(
  as.Date('2019-12-31') + 1:ncol(det.cases.mat),
  apply(det.cases.mat,2,function(ii)median(ii,na.rm=T)),
  ylim=c(0,apply(det.cases.mat,2,function(ii)quantile(ii,0.75,na.rm=T))[ncol(cases.mat)]),
  col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
  xlim=as.Date('2019-12-31') + c(31,ncol(det.cases.mat)),
  xlab='Date',ylab='Reported cases',main='',xaxt="n")
axis(1,at=c(as.Date("2020-02-01"), as.Date("2020-02-15"),as.Date("2020-03-01")),
     labels=c("Feb 1","Feb 15","Mar 1"))
polygon(
  c(as.Date('2019-12-31') + 1:ncol(det.cases.mat),
    rev(as.Date('2019-12-31') + 1:ncol(det.cases.mat))),
  c(apply(det.cases.mat,2,function(ii)quantile(ii,0.25,na.rm=T)),
    rev(apply(det.cases.mat,2,function(ii)quantile(ii,0.75,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.25))
polygon(
  c(as.Date('2019-12-31') + 1:ncol(det.cases.mat),
    rev(as.Date('2019-12-31') + 1:ncol(det.cases.mat))),
  c(apply(det.cases.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
    rev(apply(det.cases.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.15))
## for (ii in 1:nrow(det.cases.mat)) {
##     lines(
##         as.Date('2019-12-31') + 1:ncol(det.cases.mat),
##         det.cases.mat[ii,],
##         col=rgb(0,0,0,0.1),lwd=0.25,type='l')
## }
legend("topleft",lty=rep("solid",2),lwd=2,
       legend=c("Data", "Model"),col=c("red","black"),
       bty='n') 
lines(as.Date('2019-12-31') + 1:ncol(det.cases.mat),
      cases.US.local,
      type='l',lwd=2,col=2,xaxs='i',yaxs='i',
      xaxt='n',yaxt='n',xlab='',ylab='')
mtext("B",side=3,line=0, 
       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
       cex=1.2)
plot(
  as.Date('2019-12-31') + 1:ncol(p.mat),
  apply(p.mat,2,function(ii)median(ii,na.rm=T)),
  ylim=c(0,1),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
  xlim=as.Date('2019-12-31') + c(31,ncol(p.mat)),
  xlab='Date',ylab=expression('Symptomatics reporting ('*rho[local]*')'),
  main='',xaxt="n")
axis(1,at=c(as.Date("2020-02-01"), as.Date("2020-02-15"),as.Date("2020-03-01")),
     labels=c("Feb 1","Feb 15","Mar 1"))
polygon(
  c(as.Date('2019-12-31') + 1:ncol(p.mat),
    rev(as.Date('2019-12-31') + 1:ncol(p.mat))),
  c(apply(p.mat,2,function(ii)quantile(ii,0.25,na.rm=T)),
    rev(apply(p.mat,2,function(ii)quantile(ii,0.75,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.25))
polygon(
  c(as.Date('2019-12-31') + 1:ncol(p.mat),
    rev(as.Date('2019-12-31') + 1:ncol(p.mat))),
  c(apply(p.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
    rev(apply(p.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.15))
## for (ii in 1:nrow(p.mat)) {
##     lines(
##         as.Date('2019-12-31') + 1:ncol(p.mat),
##         p.mat[ii,],
##         col=rgb(0,0,0,0.1),lwd=0.25,type='l')
## }
par(new = TRUE)
plot(as.Date('2019-12-31') + testing$Day[testing$Complete],
     testing$Total[testing$Complete], type="l", col="red",
     axes=F, bty = "n", xlab = "", ylab = "",
     xlim=as.Date('2019-12-31') + c(31,ncol(cases.mat)), lwd=2,
     xaxs='i',yaxs='i')
axis(side=4, at = pretty(range(testing$Total)), col="red", col.axis="red",las=1)
mtext("Tests administered", side=4, line=3, lwd=2, col="red")        
legend("topleft", col=c("red", "black"), lty="solid",
       legend=c("Data", "Model"),
       bty="n", lwd=c(2,2))
mtext("C",side=3,line=0, 
       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.2)
hist(10*p.mat[,ncol(p.mat)],col='gray',
     xlab=expression('Symptomatics reporting ('*rho[local]*') on March 12'),
     ylab='Proportion of simulations',main='',las=1,freq=F,xaxt="n")
axis(1,at=seq(0,10,2),labels=seq(0,1,0.2))
mtext("D",side=3,line=0, 
      at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.2)
dev.off()

# Deaths results - processing
local.predict.death = t(matrix(unlist(lapply(local.predict, function(x) x$death)), length(local.predict[[1]]$death), replicates))
local.predict.death = rbinom(length(local.predict.death), as.vector(local.predict.death), propns.ASCF[,4])
local.predict.death = matrix(local.predict.death, replicates, length(local.predict[[1]]$death))
death.mat = t(matrix(
  unlist(lapply(local, function(x) x$death)),
  length(local[[1]]$death),
  replicates))
death.mat.obs = rbinom(length(death.mat), as.vector(death.mat), propns.ASCF[,4])
death.mat = matrix(death.mat.obs, replicates, ncol(death.mat))
all.death.mat = t(matrix(
  unlist(lapply(local, function(x) x$daily)),
  length(local[[1]]$cases),
  replicates))
all.death.mat.obs = rbinom(length(all.death.mat), as.vector(all.death.mat), propns.ASCF[,4])
all.death.mat = matrix(all.death.mat.obs, replicates, ncol(all.death.mat))
future.death.mat = all.death.mat-death.mat
death.day.min = which(deaths.US.local!=0)[1]
death.day.max = length(deaths.US.local)

quantile(rowSums(death.mat),c(0.025,0.5,0.975))
sum(deaths.US.local)
sum(death.mat[,death.day.min:death.day.max])/sum(death.mat)
quantile(rowSums(future.death.mat),c(0.025,0.5,0.975))
quantile(rowSums(future.death.mat)/sum(deaths.US.local),c(0.025,0.5,0.975),na.rm=T)

# Figure 3
# plot future deaths, assuming transmission stops on March 12
wdth = 9.5
hght = 5.5
pdf('../plots/figure_3_deaths_forecast.pdf',width=wdth/2.54,height=hght/2.54,pointsize=8)
par(mar=c(2.1,4.1,1.1,1.1))
time.max = 150#dim(local.predict.death)[2]
death.predict.median = apply(local.predict.death[,], 2, median)
death.predict.025 = apply(local.predict.death[,], 2, function(x)quantile(x,0.025))
death.predict.975 = apply(local.predict.death[,], 2, function(x)quantile(x,0.975))
death.predict.25 = apply(local.predict.death[,], 2, function(x)quantile(x,0.25))
death.predict.75 = apply(local.predict.death[,], 2, function(x)quantile(x,0.75))
times = seq(from=as.Date("2020-01-01"), by="1 day", length.out=time.max)
plot(times, c(deaths.US.local, rep(NA, time.max-length(deaths.US.local))),
     xlim = c(as.Date("2020-02-01"),as.Date("2020-05-15")), ylim = c(1,max(death.predict.975)),
     col="red", type="l", xlab="Month", ylab="Deaths", lwd=2, main="",xaxs='i',yaxs='i',las=1, xaxt = 'n',
     log='y')
month_starts = c(1, 30, 61, 91)
axis(side = 1, at = times[which(times %in% as.Date("2020-02-01"): as.Date("2020-05-15"))][month_starts], labels = F)
axis(side = 1, at = times[which(times %in% as.Date("2020-02-01"): as.Date("2020-05-15"))][month_starts + 15], tick = F,
     labels = c('Feb', 'Mar', 'Apr', 'May'))
lines(times, death.predict.median[1:time.max], col="black", lwd=2)
polygon(
  c(times, rev(times)),
  c(pmax(0.1,death.predict.75[1:time.max]),
    pmax(0.1,rev(death.predict.25[1:time.max]))),
  border=NA,col=rgb(0,0,0,alpha=0.25))
polygon(
  c(times, rev(times)),
  c(pmax(0.1,death.predict.975[1:time.max]),
    pmax(0.1,rev(death.predict.025[1:time.max]))),
  border=NA,col=rgb(0,0,0,alpha=0.15))
## for (ii in 1:nrow(local.predict.death)) {
##     lines(
##         times,
##         local.predict.death[ii,1:length(times)],
##         col=rgb(0,0,0,0.1),lwd=0.25,type='l')
## }
abline(v=as.Date("2020-03-12"), lty="dashed")
legend("topleft",lty=rep("solid",2),lwd=2,
       legend=c("Data", "Model"),col=c("red","black"),
       bty='n') 
dev.off()
