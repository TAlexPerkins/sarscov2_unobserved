
#=============================================================================#
# set up workspace
#=============================================================================#

# load libraries
library(extraDistr)
library(doParallel)
library(mc2d)
library(BayesianTools) #,lib="~/myRlibs-3.6.2")
library(akima)
setwd("..")
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least the replicate argument must be supplied", call.=FALSE)
}

run_num=as.numeric(args[1])

param_grid=read.csv(file="../data/sensitivity/covid_params_estimate.csv",header=T,as.is=T)
param_run=param_grid[run_num,]
#=============================================================================#
# necessary functions
#=============================================================================#

logit = function(x){
  return(log(x / (1-x)))
}

inv_logit = function(x){
  return(exp(x) / (1 + exp(x)))
}

##Param range
PrCaseSymptom.trav = c(0.01,seq(0.05, 1, by = 0.05))
asympRFraction = seq(0, 1, by = 0.05)

params.gridded = expand.grid(PrCaseSymptom.trav, asympRFraction)
names(params.gridded) = c('PrCaseSymptom.trav', 'asympRFraction')

load(paste('../results/sensitivity/mle_local_',run_num,'.rda',sep=""))

params.gridded_new = cbind(params.gridded, ll.out)
colnames(params.gridded_new)[3] = 'll'
ll.out.mat=matrix((ll.out),nrow=length(PrCaseSymptom.trav),ncol=length(asympRFraction),byrow = F)

mod=bicubic.grid(x=PrCaseSymptom.trav,
             y=asympRFraction,
             z=ll.out.mat,
             xlim=c(0,1),ylim=c(0,1), 
             dx=0.001,dy=0.001)

mod_out=as.vector(t(mod$z))
mod_x=unlist(foreach(i=1:length(mod$x))%do%{rep(mod$x[i],length(mod$y))})
mod_y=rep(mod$y,length(mod$x))

sample_out=sample(1:length(mod_out),
       size=10000,
       prob=exp(mod_out),replace = T)

PrCaseSymptom.trav_posterior = mod_x[sample_out]
asympRFraction_posterior = mod_y[sample_out]

quantile(PrCaseSymptom.trav_posterior,c(0.025,.5,.975))
quantile(asympRFraction_posterior,c(0.025,.5,.975))

save(PrCaseSymptom.trav_posterior, asympRFraction_posterior,
     file = paste('../results/sensitivity/param_estimates_posterior_',run_num,'.rda',sep=""))

