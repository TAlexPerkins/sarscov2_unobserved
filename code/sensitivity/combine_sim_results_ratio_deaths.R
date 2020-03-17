##
## Combine sim results
##
setwd("../../results/sensitivity")
rfiles=list.files(pattern="sensitivity_sims")

params_dat=read.csv(file="../../data/sensitivity/covid_params_estimate.csv",header=T,as.is = T)
dsim=1 #Default sim

sim1=unlist(lapply(strsplit(rfiles,"_"), function(x) x[3]))
sim_nums=as.numeric(unlist(lapply(strsplit(sim1,".rda"),function(x) x[1])))

params_dat=params_dat[sim_nums,]

potential.deaths=list()
potential.deaths.all=list()
for(ii in 1:length(rfiles)){
  load(rfiles[ii])
  cum.deaths.all = apply(matrix(unlist(lapply(local, function(x) t(x$daily))),
                                ncol=ncol(local[[1]]$death),byrow=T),
                         1,sum)
  cum.deaths = apply(matrix(unlist(lapply(local, function(x) t(x$death))),
                            ncol=ncol(local[[1]]$death),byrow=T),
                     1,sum)
  cum.deaths.future = cum.deaths.all - cum.deaths
  potential.deaths.all[[ii]] = rbinom(length(cum.deaths.future), as.vector(cum.deaths.future), propns.ASCF[,4]) 
  potential.deaths[[ii]] = rbinom(length(cum.deaths), as.vector(cum.deaths), propns.ASCF[,4]) 
}

output = matrix(unlist(potential.deaths), ncol = 1000, byrow = TRUE)
output.all = matrix(unlist(potential.deaths.all), ncol = 1000, byrow = TRUE)
names(output)=paste("sim",1:1000,sep="")
names(output.all)=paste("sim",1:1000,sep="")
output_dat=cbind(params_dat,output)
output_dat.all=cbind(params_dat,output.all)
names(output_dat)[19:1018]=paste("sim",1:1000,sep="")
names(output_dat.all)[19:1018]=paste("sim",1:1000,sep="")


output_dat[,19:1018] = output_dat.all[,19:1018]/output_dat[,19:1018]

##Table of results
d_dat=data.frame(sim=sim_nums,Parameter=params_dat$Variable,Scenario=params_dat$Level,
                 infs_m=round(apply(output_dat[,19:1018],1,function(x) quantile(x,0.5)),1),
                 infs_l=round(apply(output_dat[,19:1018],1,function(x) quantile(x,0.025)),1),
                 infs_u=round(apply(output_dat[,19:1018],1,function(x) quantile(x,0.975)),1))
d_dat=d_dat[order(d_dat$sim),]
write.table(d_dat,file="parameter_sensitivity_death_ratios.csv",sep=",",row.names = F,col.names = T)

library(tidyr)
output_r0=output_dat[which(output_dat$Variable=="R0"|is.na(output_dat$Variable)),]
output_r0_long=gather(output_r0,scenario,value,sim1:sim1000,factor_key = T)

output_k=output_dat[which(output_dat$Variable=="k"|is.na(output_dat$Variable)),]
output_k_long=gather(output_k,scenario,value,sim1:sim1000,factor_key = T)
##Default overdispersion was low
output_k_long$Level[output_k_long$Scenario=="Default"]="Low"

output_pr_asymp=output_dat[which(output_dat$Variable=="PrAsymp"|is.na(output_dat$Variable)),]
output_pr_asymp_long=gather(output_pr_asymp,scenario,value,sim1:sim1000,factor_key = T)

output_pr_death=output_dat[which(output_dat$Variable=="PrDeath"|is.na(output_dat$Variable)),]
output_pr_death_long=gather(output_pr_death,scenario,value,sim1:sim1000,factor_key = T)
# 
output_delay=output_dat[which(output_dat$Variable=="Rep_delay"|is.na(output_dat$Variable)),]
output_delay_long=gather(output_delay,scenario,value,sim1:sim1000,factor_key = T)

output_si=output_dat[which(output_dat$Variable=="Serial_Int"|is.na(output_dat$Variable)),]
output_si_long=gather(output_si,scenario,value,sim1:sim1000,factor_key = T)

output_inc_per=output_dat[which(output_dat$Variable=="Inc_Period"|is.na(output_dat$Variable)),]
output_inc_per_long=gather(output_inc_per,scenario,value,sim1:sim1000,factor_key = T)

output_dd=output_dat[which(output_dat$Variable=="Time_Death"|is.na(output_dat$Variable)),]
output_dd_long=gather(output_dd,scenario,value,sim1:sim1000,factor_key = T)

output_import=output_dat[which(output_dat$Variable=="Import"|is.na(output_dat$Variable)),]
output_import_long=gather(output_import,scenario,value,sim1:sim1000,factor_key = T)

library(ggplot2)

##Combo figure
output_k_long$Variable="Overdispersion"
output_delay_long$Variable="Reporting delay"
output_pr_death_long$Variable="Case fatality risk"
output_pr_asymp_long$Variable="Asymptomatic prob."
output_si_long$Variable="Serial interval"
output_r0_long$Variable="R0"
output_inc_per_long$Variable="Incubation period"
output_dd_long$Variable="Time to death"
output_import_long$Variable="Importation timing"

output_long=rbind(output_r0_long,
                  output_k_long,
                  output_delay_long,
                  output_pr_death_long,
                  output_si_long,
                  output_pr_asymp_long,
                  output_inc_per_long,
                  output_dd_long,
                  output_import_long)
output_long$Level[output_long$Level=="Middle"]="Mid"
output_long$Level=factor(output_long$Level,levels=c("Low","Mid","High"))
output_long$Scenario[output_long$Scenario!="Default"]="Sensitivity"


pdf(file="../../plots/sensitivity/sensitivity_death_ratio_facet_plot.pdf",height=12,width=12,pointsize=14,useDingbats = F)
ggplot(output_long[which(!is.na(output_long$Variable)),],aes(x=as.factor(Level),y=value,fill=Scenario))+
  geom_violin()+scale_fill_manual(values=c("#2166ac","#d1e5f0"))+
  facet_wrap(~Variable,nrow=3)+
  
  theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
        axis.title=element_text(size=16),legend.text = element_text(size=14),strip.text.x = element_text(size=14),
        panel.background = element_rect(fill="white",color="grey50"),legend.position = "",#c(0.08,0.94),
        legend.title=element_blank())+
    scale_y_log10()+xlab("")+ylab("Ratio of expected deaths to deaths through March 12")

dev.off()
