##
## Combine sim results
##
#library(doParallel)
library(boot)
setwd("../results/sensitivity")
rfiles=list.files(pattern="sensitivity_sims_")

params_dat=read.csv(file="../../data/sensitivity/covid_params_estimate.csv",header=T,as.is = T)
dsim=1 #Default sim

sim1=unlist(lapply(strsplit(rfiles,"_"), function(x) x[3]))
sim_nums=as.numeric(unlist(lapply(strsplit(sim1,".rda"),function(x) x[1])))

paste(params_dat$Variable, params_dat$Scenario)[sim_nums]
sim.names = c("Baseline",
              "Alternative importation",
              "Asymptomatic probability - Low", "Asymptomatic probability - High",
              expression('R'[0]*" - Low"), expression('R'[0]*" - High"),
              "Overdispersion - Low", "Overdispersion - Moderate",
              "Reporting delay - Low", "Reporting delay - High",
              "CFR - Low", "CFR - High",
              "Serial interval - Low", "Serial interval - High",
              "Incubation period - Low","Incubation period - High",
              "Time to death - Low", "Time to death - High")

params_dat=params_dat[sim_nums,]
# get cases.US.local
# get local cases
linelist = read.csv('../../data/2020_03_12_1800EST_linelist_NIHFogarty.csv')
yesUS = subset(linelist, country=='USA')
# remove Diamond Princess repatriated cases
yesUS = yesUS[grep("Diamond",yesUS$summary,invert=T),]
ts = read.csv('../../data/time_series_19-covid-Confirmed.csv')
ts.natl = matrix(0,length(unique(ts$Country.Region)),ncol(ts)-4)
for(ii in 1:ncol(ts.natl)){
  ts.natl[,ii] = aggregate(ts[,4+ii],by=list(ts$Country.Region),FUN=sum)[,2]
}
row.names(ts.natl) = aggregate(ts[,4+ii],by=list(ts$Country.Region),FUN=sum)[,1]
for(ii in 1:nrow(ts.natl)){
  ts.natl[ii,-1] = pmax(0,diff(ts.natl[ii,]))
}
colnames(ts.natl) = 22:(ncol(ts.natl)+21)
which(colnames(ts.natl)==34)
ts.natl['China',which(colnames(ts.natl)==34):ncol(ts.natl)] = 0
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

# read csv of test numbers
testing = read.csv("../../data/testing_ts.csv")
incomplete.days = testing$Day[!testing$Complete]


first.batch.end=12
updateDaily = FALSE
smoothSpline = TRUE
pdf("../../plots/sensitivity/combined_pDetect_plots1.pdf",width = 4.5, height = 6,pointsize=8)
par(mfrow=c(4,3))
par(oma=c(3, 3, 1, 1), mar=c(4.1,3.1,3.1,1.1))
for(ii in 1:first.batch.end){
  load(rfiles[ii])
  cases.mat = matrix(unlist(lapply(local, function(x) t(x$cases))),
                     ncol=ncol(local[[1]]$cases),byrow=T)
  p.mat = matrix(NA,nrow(cases.mat),ncol(cases.mat))
  for (jj in 1:nrow(cases.mat)) {
    alpha.old=1
    beta.old=1
    for (kk in 1:ncol(cases.mat)){
      if(cases.mat[jj,kk]){
        actual.cases = rbinom(1,cases.mat[jj,kk], sum(propns.ASCF[ii,2:3]))
        alpha.new = alpha.old+cases.US.local[kk]
        beta.new = beta.old+actual.cases-cases.US.local[kk]
        p.mat[jj,kk] = rbeta(1,alpha.new,max(1,beta.new))
        if (updateDaily) {
          alpha.old=alpha.new
          beta.old=beta.new
        }
      }
    }
    if (smoothSpline) {
      non.NA.indices = which(!is.na(p.mat[jj,]))
      if(length(non.NA.indices) > ncol(p.mat) / 3){
        temp.sp = smooth.spline((1:ncol(p.mat))[non.NA.indices],
                                logit(p.mat[jj,non.NA.indices]),
                                nknots=floor((ncol(p.mat) - non.NA.indices[1])/7 + 0.5))
        p.mat[jj,non.NA.indices[1]:ncol(p.mat)] = inv.logit(predict(temp.sp, non.NA.indices[1]:ncol(p.mat))$y)
      }else{
        # print('ELSE')
        p.mat[jj,non.NA.indices[1]:ncol(p.mat)] = NA
      }
    }
  }
  plot(
    as.Date('2019-12-31') + 1:ncol(p.mat),
    apply(p.mat,2,function(ii)median(ii,na.rm=T)),
    ylim=c(0,1),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
    xlim=as.Date('2019-12-31') + c(31,ncol(p.mat)),
    xlab="", ylab="",
    main=sim.names[ii])
  polygon(
    c(as.Date('2019-12-31') + 1:ncol(p.mat),
      rev(as.Date('2019-12-31') + 1:ncol(p.mat))),
    c(apply(p.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
      rev(apply(p.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
    border=NA,col=rgb(0,0,0,0.25))
}
mtext("Date", 1, 1, outer=TRUE)
mtext("Symptomatics detected", 2, 1, outer=TRUE, las=0)
dev.off()


pdf("../../plots/sensitivity/combined_pDetect_plots2.pdf",width = 4.5, height = 3.19,pointsize=8)
par(mfrow=c(2,3))
par(oma=c(3, 3, 1, 1), mar=c(4.1,3.1,3.1,1.1))
for(ii in ((first.batch.end+1):length(rfiles))){
  load(rfiles[ii])
  cases.mat = matrix(unlist(lapply(local, function(x) t(x$cases))),
                     ncol=ncol(local[[1]]$cases),byrow=T)
  p.mat = matrix(NA,nrow(cases.mat),ncol(cases.mat))
  for (jj in 1:nrow(cases.mat)) {
    alpha.old=1
    beta.old=1
    for (kk in 1:ncol(cases.mat)){
      if(cases.mat[jj,kk]){
        actual.cases = rbinom(1,cases.mat[jj,kk], sum(propns.ASCF[ii,2:3]))
        alpha.new = alpha.old+cases.US.local[kk]
        beta.new = beta.old+actual.cases-cases.US.local[kk]
        p.mat[jj,kk] = rbeta(1,alpha.new,max(1,beta.new))
        if (updateDaily) {
          alpha.old=alpha.new
          beta.old=beta.new
        }
      }
    }
    if (smoothSpline) {
      non.NA.indices = which(!is.na(p.mat[jj,]))
      if(length(non.NA.indices) > ncol(p.mat) / 3){
        temp.sp = smooth.spline((1:ncol(p.mat))[non.NA.indices],
                                logit(p.mat[jj,non.NA.indices]),
                                nknots=floor((ncol(p.mat) - non.NA.indices[1])/7 + 0.5))
        p.mat[jj,non.NA.indices[1]:ncol(p.mat)] = inv.logit(predict(temp.sp, non.NA.indices[1]:ncol(p.mat))$y)
      }else{
        # print('ELSE')
        p.mat[jj,non.NA.indices[1]:ncol(p.mat)] = NA
      }
    }
  }
  plot(
    as.Date('2019-12-31') + 1:ncol(p.mat),
    apply(p.mat,2,function(ii)median(ii,na.rm=T)),
    ylim=c(0,1),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
    xlim=as.Date('2019-12-31') + c(31,ncol(p.mat)),
    xlab="", ylab="",
    main=sim.names[ii])
  polygon(
    c(as.Date('2019-12-31') + 1:ncol(p.mat),
      rev(as.Date('2019-12-31') + 1:ncol(p.mat))),
    c(apply(p.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
      rev(apply(p.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
    border=NA,col=rgb(0,0,0,0.25))
}
mtext("Date", 1, 1, outer=TRUE)
mtext("Symptomatics detected", 2, 1, outer=TRUE, las=0)
dev.off()
