
#=============================================================================#
# set up workspace
#=============================================================================#

# load libraries
library(akima)
library(hexbin)
library(RColorBrewer)
library(ggplot2)
library(data.table)
setwd("..")
param_grid=read.csv(file="../data/sensitivity/covid_params_scenario_names.csv",header=T,as.is=T)

bin_list=list()
dat_list=list()
pdat=list()
for(ii in 1:18){
  run_num=ii
  load(paste('../results/sensitivity/param_estimates_posterior_',run_num,'.rda',sep=""))
  
  dat_list[[ii]]=data.frame(run=ii,pr_det=PrCaseSymptom.trav_posterior,r_asymp=asympRFraction_posterior,
                            name=param_grid$Name[ii])
  # joint posterior distribution on transformed parameters
  bin_list[[ii]] = hexbin(PrCaseSymptom.trav_posterior,
               asympRFraction_posterior,
               xbins=40)

  pp=ggplot(dat_list[[ii]],aes(x=pr_det,y=r_asymp))+geom_hex(bins=30)+
    theme(axis.text.x=element_text(size=14),axis.text.y=element_text(size=14),
          axis.title=element_text(size=16),legend.text = element_text(size=14),strip.text.x = element_text(size=14),
          panel.background = element_rect(fill="white",color="grey50"),
          legend.title=element_text(size=16),plot.title = element_text(size=16,hjust=0.5))+
    scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0,1))+
    xlab("")+ylab("")+ggtitle(param_grid$Name[ii])

  ggsave(filename=paste('../plots/sensitivity/param_estimates_posterior_hexbin',run_num,'_ggplot.pdf',sep=""),
         plot=pp,height = 6, width = 6)
  pdat[[ii]]=data.frame(scenario=param_grid$Name[ii],
             det_m=quantile(PrCaseSymptom.trav_posterior,.5),
             det_l=quantile(PrCaseSymptom.trav_posterior,.025),
             det_u=quantile(PrCaseSymptom.trav_posterior,.975),
  r_asymp_m=quantile(asympRFraction_posterior,.5),
  r_asymp_l=quantile(asympRFraction_posterior,.025),
  r_asymp_u=quantile(asympRFraction_posterior,.975))
}

#post_dat=rbindlist(dat_list)
post_s=rbindlist(pdat)

write.table(post_s,file="../results/sensitivity/parameter_estimates_ppi.csv",sep=",",row.names = F,col.names = T)

