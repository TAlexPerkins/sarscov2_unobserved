# load libraries
library(extraDistr)
library(doParallel)
library(mc2d)
#setwd("..") #Set to main code directory
# load function to simulate autochthonous transmission
source('simOutbreak.R')

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least the replicate argument must be supplied", call.=FALSE)
}

set.seed(1234)

# set up parallel cluster
cl = makeCluster(20)
registerDoParallel(cl)

##Parameter sensitivity scenario to run
run_num=as.numeric(args[1])

replicates=1000
##Parameter file
param_grid=read.csv(file="../data/sensitivity/covid_params_estimate.csv",header=T,as.is=T)
param_run=param_grid[run_num,]

##Load PrCaseSymptom.trav and rAsymp posteriors
load(paste('../results/sensitivity/param_estimates_posterior_',run_num,'_gam.rda',sep=""))
param_reps=sample(1:length(PrCaseSymptom.trav_posterior),
                  size=replicates,replace=T)
PrCaseSymptom.trav=PrCaseSymptom.trav_posterior[param_reps]
asympRFraction=asympRFraction_posterior[param_reps]


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

print(warnings())
# correct for travel ban from China (for non-US citizens starting 2/2 at 5pm) - so 0 out starting 2/3
colnames(ts.natl) = 22:(ncol(ts.natl)+21)
#which(colnames(ts.natl)==34)
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
rangeUS = sum(yesUS$international_traveler==1,na.rm=T):maxUS

# sample from uncertainty about proportions of infection outcomes
propns.ASCF = cbind(
  rbeta(replicates,param_run$PrAsymptomatic_alpha,param_run$PrAsymptomatic_beta),
  rnorm(replicates,param_run$PrDeathSymptom_mean,param_run$PrDeathSymptom_sd))
propns.ASCF[propns.ASCF<0]=0

propns.ASCF = cbind(
    propns.ASCF[,1],
    (1-propns.ASCF[,1]) * (1- PrCaseSymptom.trav) * (1-propns.ASCF[,2]),
    (1-propns.ASCF[,1]) * PrCaseSymptom.trav * (1-propns.ASCF[,2]),
    (1-propns.ASCF[,1]) * propns.ASCF[,2])

R.reps = rnorm(replicates,param_run$R_mean,
               param_run$R_sd)
R.reps[R.reps<0]=0
if (!is.na(param_run$inc_shape_sd)) {
  ip.shape.reps = rnorm(replicates, param_run$inc_shape_mean,
                        param_run$inc_shape_sd)
  ip.shape.reps[ip.shape.reps<0]=1e-10
} else {
  ip.shape.reps = rep(param_run$inc_shape_mean, replicates)
}
if (!is.na(param_run$inc_scale_sd)) {
  ip.scale.reps = rnorm(replicates, param_run$inc_scale_mean,
                        param_run$inc_scale_sd)
  ip.scale.reps[ip.scale.reps<0]=1e-10
} else {
  ip.scale.reps = rep(param_run$inc_scale_mean, replicates)
}
k.mean = param_run$k_mean
k.lower = param_run$k_lower
k.upper = param_run$k_upper
if (!(is.na(k.lower) | is.na(k.upper))) {
  k.meanlogs = seq(-10,log(k.mean),0.01)
  k.sdlogs = sqrt(2*(log(k.mean) - k.meanlogs))
  ind = which.min(sqrt((qlnorm(0.025,k.meanlogs,k.sdlogs)-k.lower)^2
                       + (qlnorm(0.975,k.meanlogs,k.sdlogs)-k.upper)^2))
  k.reps = rlnorm(replicates,k.meanlogs[ind],k.sdlogs[ind])
} else {
  k.reps = rep(k.mean,replicates)
}

# sample imported infections
imports = numeric(length=replicates)
import.doy = list()

#Scenario where we still use international incidence to estimate imports
if(param_run$Variable=="Import"&param_run$Scenario=="Incid"){
  import.doy=foreach(icnt = 1:replicates,.packages='mc2d') %dopar%{
    PrImportedInfections = 
      dmultinomial(
        x = cbind(
          0:(maxUS-sum(num.CF)),
          num.CF[1],num.CF[2]),
        prob = c(sum(propns.ASCF[icnt,1:2]),propns.ASCF[icnt,3:4]))
    imports =
      sample(
        sum(num.CF):maxUS,
        1,
        prob=PrImportedInfections,
        replace=T)
    import.when=sample(
      1:ncol(ts.natl),
      imports,
      prob=colSums(ts.natl[-which(row.names(ts.natl)=='US'),]),
      replace=T)
    
    return(as.numeric(colnames(ts.natl)[1]) - 1 + import.when)
  }
}else{
  import.doy=foreach(icnt = 1:replicates,.packages='mc2d') %dopar%{
    PrImportedInfections = 
      dmultinomial(
        x = cbind( 
         0:(maxUS-sum(num.CF)),
          num.CF[1],num.CF[2]),
        prob = c(sum(propns.ASCF[icnt,1:2]),propns.ASCF[icnt,3:4]))
    imports =
      sample(
        sum(num.CF):maxUS,
        1,
        prob=PrImportedInfections,
        replace=T)
    # if(param_run$Scenario=="DP"){
    #   #Need to exclude DP cases from shape of curve. Have already been included in number of imports above
    #   cases.US.imported =
    #     table(
    #       as.Date(as.character(subset(yesUS_noDP,traveler>0)$reporting.date)) -
    #         as.Date('2019-12-31'))
    #   tmp = rep(0,length(cases.US.total))
    #   tmp[as.numeric(names(cases.US.imported))] = cases.US.imported
    #   cases.US.imported = tmp
    # }
    ## Have imported infection curve drawn from distribution of reported imported cases rather than international incidence
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
    import.case.density=density(case.days,from = 1,to=length(cases.US.imported),n=length(cases.US.imported))$y
    
    return(sample(
      1:length(cases.US.imported),
      imports,
      prob=import.case.density,
      replace=T))
  }  
}

# simulate local transmission for each imported case draw
local = foreach(icnt = 1:replicates) %dopar% {
  simOutbreak(import.doy[[icnt]],
              R=R.reps[icnt],
              k=k.reps[icnt],
              si_mean=param_run$serial_int_mean,
              si_sd = param_run$serial_int_sd,
              report_delay_shape = param_run$rep_delay_shape,
              report_delay_rate = param_run$rep_delay_rate,
              inc_shape = ip.shape.reps[icnt],
              inc_scale = ip.scale.reps[icnt],
              symp_to_death_mean = param_run$time_death_mean,
              symp_to_death_sd = param_run$time_death_sd,
              asympProp=propns.ASCF[icnt,1],
              asympRFraction=asympRFraction[icnt],
              stopSimulationDay = length(cases.US.imported),
              lnormFlag = param_run$serial_int_lnorm)
  
}

stopCluster(cl)

save(local,import.doy,propns.ASCF,file=paste("../results/sensitivity/sensitivity_sims_",run_num,"_gam.rda",sep=""))
