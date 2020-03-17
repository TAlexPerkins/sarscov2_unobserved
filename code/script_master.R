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

# read in line list data for US
# updated 20200307
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

# sample replicates of how many infections have been imported into the US
maxUS = 2e4
rangeUS = sum(yesUS$international_traveler>0,na.rm=T):maxUS
# estimate for asymptomatic proportion based on
# https://www.medrxiv.org/content/10.1101/2020.02.20.20025866v2
PrAsymptomatic = exp(optim(par=c(0,0),fn=function(par){
  sum((
    qbeta(c(0.5,0.025,0.975),exp(par[1]),exp(par[2])) -
      c(0.179,0.155,0.202)) ^ 2)})$par)
# estimate for proportion of symptomatic infections resulting in death based on
# http://weekly.chinacdc.cn/en/article/id/e53946e2-c6c4-41e9-9a9b-fea8db1a8f51
PrDeathSymptom = c(1+1023,1+44672-1023)

# set values of unknown parameters
# note that these values seem to maximize the probability of the cumulative
# deaths in the US as of March 8, 2020 predicted by the model
replicates = 1000
load("../results/param_estimates/param_estimates_posterior_1.rda", verbose=T)
indices = sample(1:length(PrCaseSymptom.trav_posterior), replicates, replace=TRUE)
PrCaseSymptom.trav = PrCaseSymptom.trav_posterior[indices]
asympRFraction = asympRFraction_posterior[indices]

# sample from uncertainty about proportions of infection outcomes
propns.ASCF = cbind(
  rbeta(replicates,PrAsymptomatic[1],PrAsymptomatic[2]),
  rbeta(replicates,PrDeathSymptom[1],PrDeathSymptom[2]))
propns.ASCF = cbind(
  propns.ASCF[,1],
  (1-propns.ASCF[,1]) * (1-PrCaseSymptom.trav) * (1-propns.ASCF[,2]),
  (1-propns.ASCF[,1]) * PrCaseSymptom.trav * (1-propns.ASCF[,2]),
  (1-propns.ASCF[,1]) * propns.ASCF[,2])

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
# mtext("A",side=3,line=0, 
#       at=par("usr")[1]+0.02*diff(par("usr")[1:2]),
#       cex=1.2)
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

import.doy = list()
for(ii in 1:replicates){
  import.doy[[ii]] = sample(
    1:length(cases.US.imported),
    imports[ii],
    prob=import.case.density,
    replace=T)
}

# simulate local transmission for each draw of imported infections
local = foreach(ii = 1:replicates) %do% {
  simOutbreak(
    timeImport = import.doy[[ii]], # timing of each imported infection
    R = 1.97, # reproduction number
    k = 1e3, # dispersion parameter
    si_mean = 4.56, # mean of serial interval distribution
    si_sd = 0.95, # standard deviation of serial interval distribution
    inc_shape = 1.88, # shape parameter of incubation period distribution
    inc_scale = 7.97, # scale parameter of incubation period distribution
    symp_to_death_mean = 14, # mean of time between symptom onset and death
    symp_to_death_sd = 5.37, # std. dev. of time between symptom onset and death
    report_delay_shape = delay.shape.baseline, # shape parameter for delay between symptom and reporting
    report_delay_rate = delay.rate, # rate parameter for delay between symptom and reporting
    stopSimulationDay = length(cases.US.imported), # day of year since Jan 1 when simulation stops
    asympProp = propns.ASCF[ii,1], # proportion of infections that are asymptomatic
    asympRFraction = asympRFraction[ii], # relative infectiousness of asymptomatics
    lnormFlag = F # toggles whether serial interval distribution is lognormal
  )
}

local.predict = foreach(ii = 1:replicates) %do% {
  simOutbreakR0Change(
    timeImport = import.doy[[ii]], # timing of each imported infection
    R = 1.97, # reproduction number
    k = 1e3, # dispersion parameter
    si_mean = 4.56, # mean of serial interval distribution
    si_sd = 0.95, # standard deviation of serial interval distribution
    inc_shape = 1.88, # shape parameter of incubation period distribution
    inc_scale = 7.97, # scale parameter of incubation period distribution
    symp_to_death_mean = 14, # mean of time between symptom onset and death
    symp_to_death_sd = 5.37, # std. dev. of time between symptom onset and death
    report_delay_shape = delay.shape.baseline, # shape parameter for delay between symptom and reporting
    report_delay_rate = delay.rate, # rate parameter for delay between symptom and reporting
    stopSimulationDay = 180, # day of year since Jan 1 when simulation stops
    asympProp = propns.ASCF[ii,1], # proportion of infections that are asymptomatic
    asympRFraction = asympRFraction[ii], # relative infectiousness of asymptomatics
    lnormFlag = F, # toggles whether serial interval distribution is lognormal
    RChangeDay = length(cases.US.imported), # determines when R changes
    RChange = 0 # determines what R drops to at R0ChangeDay
  )
}

save(local, local.predict, imports, import.doy, file="../data/transmission_model_outputs.RData")

#load("../data/transmission_model_outputs.RData", verbose = T)
load("../results/param_sensitivity/sensitivity_sims_1_upd_3-12.rda",verbose=T)

# set figure margins

par(mar=c(4,5,1,1))
# Figure 1
# plot distribution of cumulative infections
# plot all locally acquired infections over time
pdf('../plots/figure_1_cumulative_infections_and_infections_daily.pdf',
    width=9,height=5, pointsize=14)
par(mfrow=c(1,2))
hist(
  unlist(lapply(local,function(ll)ll$cum)),
  col='gray',xlab='Cumulative infections',
  ylab='Number of simulations',main='',las=1)
mtext("A",side=3,line=0, 
       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
       cex=1.2)
quantile(unlist(lapply(local,function(ll)ll$cum)),c(0.025,0.5,0.975))
local.mat = t(matrix(
  unlist(lapply(local, function(x) x$daily)),
  length(local[[1]]$daily),
  replicates))
quantile(local.mat[,ncol(local.mat)],c(0.025,0.5,0.975))
plot(
  as.Date('2019-12-31') + 1:ncol(local.mat),
  apply(local.mat,2,function(ii)median(ii,na.rm=T)),
  ylim=c(0,quantile(local.mat[,ncol(local.mat)],0.975)),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
  xlim=as.Date('2019-12-31') + c(31,ncol(local.mat)),
  xlab='Date',ylab='Infections',main='')
polygon(
  c(as.Date('2019-12-31') + 1:ncol(local.mat),
    rev(as.Date('2019-12-31') + 1:ncol(local.mat))),
  c(apply(local.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
    rev(apply(local.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.25))
mtext("B",side=3,line=0, 
       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
       cex=1.2)
dev.off()

# Figure 2
# plot locally acquired symptomatic infections over time, alongside number of cases reported in the US
# plot proportion of locally acquired symptomatic infections reported over time
# alongside numbers of tests administered in the US
updateDaily = FALSE # turn on Bayesian daily updating
smoothSpline = TRUE # turn on smoothing spline
pdf('../plots/figure_2_symptomatic_daily_and_symptomatic_detected.pdf',
    width=9,height=4.8, pointsize=14)
par(mfrow=c(1,2))
# figures with two vertical axes                                        
par(mar = c(5, 4, 4, 4) + 0.3)
cases.mat = t(matrix(
  unlist(lapply(local, function(x) x$cases)),
  length(local[[1]]$cases),
  replicates))
# cum.cases.vec = apply(cases.mat, 1, sum)
# quantile(cum.cases.vec,c(0.025, 0.5,0.975))
cases.mat.obs = rbinom(length(cases.mat), as.vector(cases.mat), rowSums(propns.ASCF[,2:3]))
cases.mat = matrix(cases.mat.obs, replicates, ncol(cases.mat))
plot(
  as.Date('2019-12-31') + 1:ncol(cases.mat),
  apply(cases.mat,2,function(ii)median(ii,na.rm=T)),
  ylim=c(0,quantile(cases.mat[,ncol(cases.mat)],0.975)),
  col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
  xlim=as.Date('2019-12-31') + c(31,ncol(cases.mat)),
  xlab='Date',ylab='Symptomatic infections',main='')
polygon(
  c(as.Date('2019-12-31') + 1:ncol(cases.mat),
    rev(as.Date('2019-12-31') + 1:ncol(cases.mat))),
  c(apply(cases.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
    rev(apply(cases.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.25))
mtext("A",side=3,line=0, 
       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
       cex=1.2)
par(new = TRUE)
plot(as.Date('2019-12-31') + 1:ncol(cases.mat),
     cases.US.local, type="l", col="red",
     axes=F, bty = "n", xlab = "", ylab = "",
     xlim=as.Date('2019-12-31') + c(31,ncol(cases.mat)), lwd=2,
     xaxs='i',yaxs='i')
axis(side=4, at = pretty(range(cases.US.local)), col="red", col.axis="red",las=1)
mtext("Reported cases", side=4, line=3, lwd=2, col="red")        
legend("topleft", col=c("red", "black"), lty="solid",
       legend=rev(c("Model", "Data")),
       bty="n", lwd=2)
cases.mat = t(matrix(
  unlist(lapply(local, function(x) x$cases)),
  length(local[[1]]$cases),
  replicates))
p.mat = matrix(NA,nrow(cases.mat),ncol(cases.mat))
for(ii in 1:nrow(cases.mat)){
  alpha.old=1
  beta.old=1
  for(jj in 1:ncol(cases.mat)){
    if(cases.mat[ii,jj]){
      actual.cases = rbinom(1,cases.mat[ii,jj], sum(propns.ASCF[ii,2:3]))
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
  c(apply(p.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
    rev(apply(p.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.25))
# read csv of test numbers
testing = read.csv("../data/testing_ts.csv")
incomplete.days = testing$Day[!testing$Complete]
## rect(as.Date('2019-12-31')+incomplete.days[1], 0,
##      as.Date('2019-12-31')+incomplete.days[length(incomplete.days)], 1,
##      col = rgb(1,0.5,0.5,0.25),border=NA)
# obtain low and high values of p and dates
low.p = which.min(apply(p.mat,2,function(ii)median(ii,na.rm=T)))
low.p+as.Date('2019-12-31')
apply(p.mat,2,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=T))[,low.p]
apply(p.mat,2,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=T))[,72]
par(new = TRUE)
plot(as.Date('2019-12-31') + testing$Day[testing$Complete],
     testing$Total[testing$Complete], type="l", col="red",
     axes=F, bty = "n", xlab = "", ylab = "",
     xlim=as.Date('2019-12-31') + c(31,ncol(cases.mat)), lwd=2,
     xaxs='i',yaxs='i')
axis(side=4, at = pretty(range(testing$Total)), col="red", col.axis="red",las=1)
mtext("Tests administered", side=4, line=3, lwd=2, col="red")        
legend("top", col=c("red", "black"), lty="solid",
       legend=c("Data", "Model"),
       bty="n", lwd=2)
mtext("B",side=3,line=0, 
       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
       cex=1.2)
dev.off()

# Figure 3
# plot future deaths, assuming transmission stops on March 12
pdf('../plots/figure_3_deaths_forecast.pdf',width=9,height=6,pointsize=14)
local.predict.death = t(matrix(unlist(lapply(local.predict, function(x) x$death)), length(local.predict[[1]]$death), replicates))
# estimate deaths from model. last col of propns.ACSF applies to infections only - asymps already removed at line 44
local.predict.death = rbinom(length(local.predict.death), as.vector(local.predict.death), propns.ASCF[,4])
local.predict.death = matrix(local.predict.death, replicates, length(local.predict[[1]]$death))
range = 1:replicates
time.max = 150#dim(local.predict.death)[2]
death.predict.median = apply(local.predict.death[range,], 2, median)
death.predict.025 = apply(local.predict.death[range,], 2, function(x)quantile(x,0.025))
death.predict.975 = apply(local.predict.death[range,], 2, function(x)quantile(x,0.975))
death.predict.975 = apply(local.predict.death[range,], 2, function(x) quantile(x,0.975))
death.predict.025 = apply(local.predict.death[range,], 2, function(x) quantile(x,0.025))
times = seq(from=as.Date("2020-01-01"), by="1 day", length.out=time.max)
plot(times, c(deaths.US.local, rep(NA, time.max-length(deaths.US.local))),
     xlim = c(as.Date("2020-02-01"),as.Date("2020-05-15")), ylim = c(0,max(death.predict.975)),
     col="red", type="l", xlab="Month", ylab="Deaths", lwd=2, main="",xaxs='i',yaxs='i',las=1, xaxt = 'n')
month_starts = c(1, 30, 61, 91)
axis(side = 1, at = times[which(times %in% as.Date("2020-02-01"): as.Date("2020-05-15"))][month_starts], labels = F)
axis(side = 1, at = times[which(times %in% as.Date("2020-02-01"): as.Date("2020-05-15"))][month_starts + 15], tick = F,
     labels = c('Feb', 'Mar', 'Apr', 'May'))
lines(times, death.predict.median[1:time.max], col="black", lwd=2)
polygon(
  c(times, rev(times)),
  c(death.predict.975[1:time.max],
    rev(death.predict.025[1:time.max])),
  border=NA,col=rgb(0,0,0,alpha=0.25))
abline(v=as.Date("2020-03-12"), lty="dashed")
legend("topleft",lty=rep("solid",2),lwd=2,
       legend=c("Data", "Model"),col=c("red","black"),
       bty='n') 
dev.off()

# obtain correlations of p.mat and testing
# correlation testing p.mat
totalTestsPerDay = c(rep(0, min(testing$Day)-1),
                     testing$Total[1:which(testing$Day==ncol(p.mat))])
start.day.plot = low.p
end.day.plot = max(testing$Day[testing$Complete])
correlation.reduced = cor(t(p.mat[,start.day.plot:end.day.plot]),
                          totalTestsPerDay[start.day.plot:end.day.plot], use="pairwise")
quantile(correlation.reduced,c(0.025,0.5,0.975))

# death numbers in manuscript
# process objects
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
# deaths so far
quantile(rowSums(death.mat),c(0.025,0.5,0.975))
sum(deaths.US.local)
# proportion of deaths in days when autochthonous deaths occurred
death.day.min = which(deaths.US.local!=0)[1]
death.day.max = length(deaths.US.local)
sum(death.mat[,death.day.min:death.day.max])/sum(death.mat)
# deaths which will happen vs deaths that have happened.
quantile(rowSums(future.death.mat),c(0.025,0.5,0.975))
quantile(rowSums(future.death.mat)/rowSums(death.mat),c(0.025,0.5,0.975))

# Other figures
# plot distribution of deaths expected to eventually result
# from infections through the end of the simulation
pdf('../plots/eventual_deaths.pdf',width=5,height=4)
hist(
  unlist(lapply(local,function(ll)ll$cum)) * propns.ASCF[,4],
  col='gray',xlab='Infections so far that will result in death',
  ylab='Number of simulations',main='',las=1)
dev.off()


# obtain correlations of p.mat and testing
# correlation testing p.mat
totalTestsPerDay = c(rep(0, min(testing$Day)-1),
                     testing$Total[1:which(testing$Day==ncol(p.mat))])
# whole timeseries
correlation.all = cor(t(p.mat), totalTestsPerDay, use="pairwise")
rank.correlation.all = cor(t(p.mat), totalTestsPerDay,
                           use="pairwise", method="spearman")
# reduced timeseries
start.day.plot = low.p
end.day.plot = max(testing$Day[testing$Complete])
correlation.reduced = cor(t(p.mat[,start.day.plot:end.day.plot]),
                          totalTestsPerDay[start.day.plot:end.day.plot], use="pairwise")
rank.correlation.reduced = cor(t(p.mat[,start.day.plot:end.day.plot]),
                               totalTestsPerDay[start.day.plot:end.day.plot],
                               use="pairwise", method="spearman")

quantile(correlation.all,c(0.025,0.5,0.975))
quantile(correlation.reduced,c(0.025,0.5,0.975))

pdf("../plots/scatter_test_pDetect.pdf",height=4, width=5)
plot(totalTestsPerDay[start.day.plot:end.day.plot],
     apply(p.mat[,start.day.plot:end.day.plot],2, median),
     xlim=c(0,max(totalTestsPerDay[start.day.plot:end.day.plot])),
     ylim=c(0,0.2),pch=20, bg="black",
     col=1,type='p',xaxs='i',yaxs='i',las=1,
     xlab='Tests administered',ylab='Symptomatics reporting',
     main='')
dev.off()

# plot deaths resulting from locally acquired infections over time
# note that this does not show deaths that are assumed to happen in
# the future as a result of recent infection 
pdf('../plots/deaths_daily.pdf',width=5,height=4)
plot(
  as.Date('2019-12-31') + 1:ncol(death.mat),
  apply(death.mat,2,function(ii)median(ii,na.rm=T)),
  ylim=c(0,max(c(deaths.US.local,death.mat))),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
  xlim=as.Date('2019-12-31') + c(31,ncol(death.mat)),
  xlab='Date',ylab='Deaths',main='')
polygon(
  c(as.Date('2019-12-31') + 1:ncol(death.mat),
    rev(as.Date('2019-12-31') + 1:ncol(death.mat))),
  c(apply(death.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
    rev(apply(death.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.25))
legend("topleft",lty=rep("solid",2),lwd=2,
       legend=c("Data", "Model"),col=c("red","black"),
       bty='n') 
par(new=T)
plot(deaths.US.local[31:length(deaths.US.local)],
     type='l',lwd=2,col=2,ylim=c(0,6),xaxs='i',yaxs='i',
     xaxt='n',yaxt='n',xlab='',ylab='')
dev.off()


# plot predicted cases and actual cases
pdf('../plots/cases_daily_resample.pdf',width=5,height=4)
cases.mat = t(matrix(
  unlist(lapply(local, function(x) x$cases)),
  length(local[[1]]$cases),
  replicates))
p.mat.resample = apply(p.mat,2,function(x) sample(x,length(x),replace=T))
det.cases.mat.obs = rbinom(length(cases.mat), as.vector(cases.mat), rowSums(propns.ASCF[,2:3])*p.mat.resample)
det.cases.mat = matrix(det.cases.mat.obs, replicates, ncol(cases.mat))
plot(
  as.Date('2019-12-31') + 1:ncol(det.cases.mat),
  apply(det.cases.mat,2,function(ii)median(ii,na.rm=T)),
  ylim=c(0,800),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
  xlim=as.Date('2019-12-31') + c(31,ncol(det.cases.mat)),
  xlab='Date',ylab='Reported cases',main='')
polygon(
  c(as.Date('2019-12-31') + 1:ncol(det.cases.mat),
    rev(as.Date('2019-12-31') + 1:ncol(det.cases.mat))),
  c(apply(det.cases.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
    rev(apply(det.cases.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
  border=NA,col=rgb(0,0,0,0.25))
legend("topleft",lty=rep("solid",2),lwd=2,
       legend=c("Data", "Model"),col=c("red","black"),
       bty='n') 
lines(as.Date('2019-12-31') + 1:ncol(det.cases.mat),
      cases.US.local,
      type='l',lwd=2,col=2,xaxs='i',yaxs='i',
      xaxt='n',yaxt='n',xlab='',ylab='')
quantile(rowSums(det.cases.mat, na.rm=T),c(0.025,0.5,0.975))
sum(cases.US.local)
dev.off()

# Old vs new methods of pLocal
updateDaily = FALSE # turn on Bayesian daily updating
pdf('../plots/smooth_vs_nonsmooth_pLocal.pdf',
    width=9,height=4.8, pointsize=14)
par(mfrow=c(1,2))
for (smoothSpline in c(FALSE, TRUE)) {
    cases.mat = t(matrix(
        unlist(lapply(local, function(x) x$cases)),
        length(local[[1]]$cases),
        replicates))
    p.mat = matrix(NA,nrow(cases.mat),ncol(cases.mat))
    for(ii in 1:nrow(cases.mat)){
        alpha.old=1
        beta.old=1
        for(jj in 1:ncol(cases.mat)){
            if(cases.mat[ii,jj]){
                actual.cases = rbinom(1,cases.mat[ii,jj], sum(propns.ASCF[ii,2:3]))
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
        c(apply(p.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
          rev(apply(p.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
        border=NA,col=rgb(0,0,0,0.25))
    mtext(ifelse(smoothSpline,"B","A"),side=3,line=0, 
       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
       cex=1.2)
}
dev.off()

# Comparing to Shaman et al.
# weighted p.mat
p.mat.weighted = apply(p.mat*local.mat,1,function(x) sum(x, na.rm=T))/rowSums(local.mat,na.rm=T)
quantile(p.mat.weighted,c(0.025,0.5,0.975))
undocumented.prop.ave = propns.ASCF[,1] + (1-p.mat.weighted)*rowSums(propns.ASCF[,2:3])
undocumented.alpha.ave = (asympRFraction*propns.ASCF[,1] + (1-p.mat.weighted)*rowSums(propns.ASCF[,2:3]))/undocumented.prop.ave
quantile(undocumented.prop.ave,c(0.025,0.5,0.975))
quantile(undocumented.alpha.ave,c(0.025,0.5,0.975))

#daily p.mat
undocumented.prop = propns.ASCF[,1] + (1-p.mat)*rowSums(propns.ASCF[,2:3])
undocumented.alpha = (asympRFraction*propns.ASCF[,1] + (1-p.mat)*rowSums(propns.ASCF[,2:3]))/undocumented.prop
quantile(undocumented.prop[,low.p],c(0.025,0.5,0.975))
quantile(undocumented.prop[,ncol(p.mat)],c(0.025,0.5,0.975))
quantile(undocumented.alpha[,low.p],c(0.025,0.5,0.975))
quantile(undocumented.alpha[,ncol(p.mat)],c(0.025,0.5,0.975))
