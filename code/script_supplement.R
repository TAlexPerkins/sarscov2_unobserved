#=============================================================================#
# Authors: Alex Perkins, Sean Cavany, Sean Moore, Rachel Oidtman, Anita Lerch, Marya Poterek
# project: Estimating unobserved SARS-CoV-2 infections in the United States
# Year: 2020
# 
# Code to generate all figures and results from supplementary text
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
# deaths in the US as of March 8, 2020 predicted by the model
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
    ip.shape.reps = rep(parameters["incubation_scale",1], replicates)
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
# produce plots and results for all supplementary text figures 
#=============================================================================#

# Figure S1 plot reporting delay distribution
pdf("../plots/gamma_reporting_delay.pdf", width=4, height=4, pointsize=10)
h = hist(data.delay, plot=F)
h$counts = h$counts / sum(h$counts)
plot(h, freq=TRUE, ylab="Relative Frequency",
     xlab="Delay (days)",
     xaxs='i',yaxs='i',las=1,
     xlim = c(0,16),main='')
lines(seq(0.1,16,0.1),dgamma(seq(0.1,16,0.1),
                             delay.shape.baseline, delay.rate),
      lwd=2)
dev.off()

# Figure S2 = smooth vs non-smoothed methods of pLocal
updateDaily = FALSE # turn on Bayesian daily updating
cases.mat = t(matrix(
  unlist(lapply(local, function(x) x$cases)),
  length(local[[1]]$cases),
  replicates))
pdf('../plots/smooth_vs_nonsmooth_pLocal.pdf',
    width=9,height=4.8, pointsize=14)
par(mfrow=c(1,2))
for (smoothSpline in c(FALSE, TRUE)) {
    p.mat = matrix(NA,nrow(cases.mat),ncol(cases.mat))
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
    plot(
        as.Date('2019-12-31') + 1:ncol(p.mat),
        apply(p.mat,2,function(ii)median(ii,na.rm=T)),
        ylim=c(0,1),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
        xlim=as.Date('2019-12-31') + c(31,ncol(p.mat)),
        xlab='Date',ylab='Symptomatics reporting',
        main='')
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
    mtext(ifelse(smoothSpline,"B","A"),side=3,line=0, 
       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
       cex=1.2)
}
dev.off()

# Figure S3 - plot alternative importation scenario
##Plots of import densities
alt.import.density=c(rep(0,21),
                     colSums(ts.natl[-which(row.names(ts.natl)=='US'),])/sum(ts.natl[-which(row.names(ts.natl)=='US'),]))
pdf('../plots/importation_patterns.pdf',
    width=7,height=7, pointsize=14)
par(mar = c(5, 4, 4, 4) + 0.3)
barplot(cases.US.imported,ylab='Imported cases',xlab='Date',ylim=c(0,35),
        xaxs='i',yaxs='i',las=1,cex.axis=1.2,cex.lab=1.3,width=0.85)
axis(1,at=seq(1,length(cases.US.imported),by=15),
     labels=seq(as.Date('2019-12-31'),as.Date('2019-12-31')+length(cases.US.imported)+1,by=15),
     cex.axis=1.2,cex.lab=1.3)
box()
par(new = TRUE)
plot(as.Date('2019-12-31') + 1:length(import.case.density),
     import.case.density, type="l", col="red",
     axes=F, bty = "n", xlab = "", ylab = "",
     lwd=2,
     xaxs='i',yaxs='i',cex.axis=1.2,cex.lab=1.3)
lines(as.Date('2019-12-31') + 1:length(import.case.density),
      alt.import.density, type="l", col="blue",lwd=2)
axis(side=4, at = pretty(range(import.case.density)), col.axis="black",las=1,cex.axis=1.2,cex.lab=1.3)
mtext("Importation timing", side=4, line=3.2, lwd=2,cex=1.3)        
legend("topleft", col=c("black","red", "blue"), lty=c(0,1,1),#fill=c("grey",NA,NA),
       legend=(c("Imported cases","Baseline importation", "Alternative importation")),cex=1.2,
       pch=c(22,NA,NA),pt.bg=c("grey",NA,NA),pt.cex=2,
       bty="n", lwd=2)
dev.off()

## Figure S4 p.mat histogram
pdf('../plots/rho_local_histogram.pdf',
    ##width=4.5,height=13.8, pointsize=14)
    width=4.8,height=4.8, pointsize=14)
hist(10*p.mat[,ncol(p.mat)],col='gray',
     xlab=expression('Symptomatics reporting ('*rho[local]*')'),
     ylab='Proportion of simulations',main='',las=1,freq=F,xaxt="n")
axis(1,at=seq(0,10,2),labels=seq(0,1,0.2))
dev.off()

## Compare predicted cases with actual
cases.mat.new = t(matrix(
  unlist(lapply(local.predict, function(x) x$cases[1:length(local[[1]]$cases)])),
  length(local[[1]]$cases),
  replicates))
det.cases.mat.obs = rbinom(length(cases.mat.new), as.vector(cases.mat.new), rowSums(propns.ASCF[,2:4])*p.mat)
det.cases.mat = matrix(det.cases.mat.obs, replicates, ncol(cases.mat.new))
quantile(apply(det.cases.mat,1,function(ii)sum(ii,na.rm=T)),c(0.025,0.5, 0.975))
