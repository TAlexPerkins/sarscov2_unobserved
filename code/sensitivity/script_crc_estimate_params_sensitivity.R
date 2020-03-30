# load libraries
library(extraDistr)
library(doParallel)
library(mc2d)
setwd("..") #Set to main code directory
# load function to simulate autochthonous transmission
source('simOutbreak.R')

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least the replicate argument must be supplied", call.=FALSE)
}

# set up parallel cluster
cl = makeCluster(48)
registerDoParallel(cl)

##Parameter sensitivity scenario to run
run_num=as.numeric(args[1])

replicates=200
##Parameter file
param_grid=read.csv(file="../data/sensitivity/covid_params_estimate.csv",header=T,as.is=T)
param_run=param_grid[run_num,]

# negative log likleihood
ll.localDeaths = function(local_in, propns.ASCF_in){
  
  log(max(1e-300,mean(dbinom(x = cum.deaths.US.total,
                  size = unlist(lapply(local_in,sum)),
                  prob = propns.ASCF_in, log = F))))
}

# ll.localDeaths.pois = function(local_in, propns.ASCF_in){
#   
#   log(mean(dpois(x = deaths.US.total,
#                  lambda = min(unlist(lapply(local_in,sum))*propns.ASCF_in, log = F)))
# }

# read in case data for US
# updated 20200312
# data from https://github.com/midas-network/COVID-19/tree/master/data/cases/global/line_listings_nihfogarty
linelist = read.csv('../data/2020_03_12_1800EST_linelist_NIHFogarty.csv')
yesUS = subset(linelist, country=='USA')
#Only keep Diamond Princess if Import with DP included scenario
yesUS_noDP = yesUS[grep("Diamond",yesUS$summary,invert=T),]
if(param_run$Scenario!="DP"){
  yesUS = yesUS_noDP
}

# number of travelers that were cases or died
num.CF = c(
  nrow(subset(yesUS,international_traveler>0))-sum(subset(yesUS,international_traveler>0)$death>0,na.rm=T),
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

# Exclude China incidence starting 2/3 to exclude from imports
colnames(ts.natl) = 22:(ncol(ts.natl)+21)
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
cum.deaths.US.total = sum(deaths.US.local)

# sample replicates of how many infections have been imported into the US
maxUS = 2e4
rangeUS = sum(yesUS$international_traveler==1,na.rm=T):maxUS

##Param range
PrCaseSymptom.trav = c(0.01,seq(0.05, 1, by = 0.05))
asympRFraction = seq(0, 1, by = 0.05)

params.gridded = expand.grid(PrCaseSymptom.trav, asympRFraction)
names(params.gridded) = c('PrCaseSymptom.trav', 'asympRFraction')

local = list()
length(local) = nrow(params.gridded)
ll.out = rep(NA, nrow = nrow(params.gridded))
propns.ASCF.list=list()
for(ii in 1:nrow(params.gridded)){
  # sample from uncertainty about proportions of infection outcomes
  propns.ASCF = cbind(
    rbeta(replicates,param_run$PrAsymptomatic_alpha,param_run$PrAsymptomatic_beta),
    rbeta(replicates,param_run$PrDeathSymptom_alpha,param_run$PrDeathSymptom_beta))
  
  propns.ASCF = cbind(
    propns.ASCF[,1],
    (1-propns.ASCF[,1]) * (1- params.gridded$PrCaseSymptom.trav[ii]) * (1-propns.ASCF[,2]),
    (1-propns.ASCF[,1]) * params.gridded$PrCaseSymptom.trav[ii] * (1-propns.ASCF[,2]),
    (1-propns.ASCF[,1]) * propns.ASCF[,2])
  propns.ASCF.list[[ii]]=propns.ASCF
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
  local[[ii]] = foreach(icnt = 1:replicates) %dopar% {
    simOutbreak(import.doy[[icnt]],
                R=param_run$r_val,
                k=param_run$k_val,
                si_mean=param_run$serial_int_mean,
                si_sd = param_run$serial_int_sd,
                report_delay_shape = param_run$rep_delay_shape,
                report_delay_rate = param_run$rep_delay_rate,
                inc_shape = param_run$inc_shape,
                inc_scale = param_run$inc_scale,
                symp_to_death_mean = param_run$time_death_mean,
                symp_to_death_sd = param_run$time_death_sd,
                asympProp=propns.ASCF[icnt,1],
                asympRFraction=params.gridded$asympRFraction[ii],
                stopSimulationDay = length(cases.US.imported),
                lnormFlag = param_run$serial_int_lnorm)

  }
  ll.out[ii] = ll.localDeaths(local_in = lapply(local[[ii]],function(x) x$death), 
                                  propns.ASCF_in = propns.ASCF[,4])
}

stopCluster(cl)

save(local,params.gridded,ll.out,propns.ASCF.list,file=paste("../results/sensitivity/mle_local_",run_num,".rda",sep=""))
