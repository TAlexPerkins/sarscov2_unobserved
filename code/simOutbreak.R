## function to simulate autochthonous transmission

simOutbreak = function(
  timeImport,
  R = 1.97,
  k = 1e3,
  si_mean = 4.56,
  si_sd = 0.95,
  inc_shape = 1.88,
  inc_scale = 7.97,
  symp_to_death_mean = 14,
  symp_to_death_sd = 5.37,
  report_delay_shape = 3.434,
  report_delay_rate = 0.5724,
  stopSimulationDay = 68,
  asympProp = 0.35,
  asympRFraction = 0.5,
  lnormFlag = T)
{
  # decide what to keep track of across replicate simulations
  daily.mat = matrix(0,1,stopSimulationDay)
    
  case.mat = matrix(0,1,stopSimulationDay)
  death.mat = matrix(0,1,stopSimulationDay)
  
  cum.vec = rep(0, 1)
  
  
  # calculate meanlog and sdlogs for lnorm distributions
  if(lnormFlag){
    si_sdlog = sqrt(log((si_sd/si_mean)^2 + 1))
    si_meanlog = log(si_mean) - 0.5*si_sdlog^2
  }
  
  symp_to_death_sdlog = sqrt(log((symp_to_death_sd/symp_to_death_mean)^2 + 1))
  symp_to_death_meanlog = log(symp_to_death_mean) - 0.5*symp_to_death_sdlog^2
  
  # initiate vector to store daily incidence of autochthonous infections
  dailyIncidence = rep(0,stopSimulationDay)

  # initiate vector to store daily reporting of autochthonous infections
  dailyCases = rep(0,stopSimulationDay)
  
  # initiate vector to store daily mortality of autochthonous infections
  dailyMortality = rep(0,stopSimulationDay)
      
  # initiate vector to store timing of parent infections
  time.exp = timeImport
    
  # loop through generations of transmission until extinct or time is up
  while(length(time.exp) > 0){
    
    # draw a number of offspring for each parent  
    R_vec = ifelse(runif(length(time.exp),0,1) < asympProp, R*asympRFraction, R)
    number.offspring = rnbinom(n=length(time.exp), mu=R_vec, size=k)
      
    if(sum(number.offspring > 0)){
        
      # remove parents with zero offspring
      time.exp = time.exp[number.offspring>0]
      number.offspring = number.offspring[number.offspring>0]

      # determine the time of exposure of each offspring
      if(!lnormFlag) {
        time.exp =
          rep(time.exp,times=number.offspring) +
          rnorm(sum(number.offspring), mean=si_mean, sd=si_sd)
      } else {
        time.exp =
          rep(time.exp,times=number.offspring) +
          rlnorm(sum(number.offspring), meanlog=si_meanlog, sdlog=si_sdlog)
      }
        
      # retain only those infections that occur before time is up
      time.exp =
        time.exp[floor(time.exp) <= stopSimulationDay]

      # determine incubation period for each offspring
      incubation_periods = rweibull(length(time.exp), shape=inc_shape, scale=inc_scale)
          
      # determine the time of infection detection of each offspring
      time.det =
        time.exp +
        incubation_periods + 
        rgamma(length(time.exp),
               shape=report_delay_shape, rate=report_delay_rate)
          
      # retain only those detected infections that occur before time is up
      time.det =
        time.det[floor(time.det) <= stopSimulationDay]
  
      # update vector of daily incidence of detected infections      
      dailyIncidence[as.numeric(names(table(floor(time.exp))))] =
        dailyIncidence[as.numeric(names(table(floor(time.exp))))] +
        table(floor(time.exp))
         
      # update vector of daily case reporting
      dailyCases[as.numeric(names(table(floor(time.det))))] =
        dailyCases[as.numeric(names(table(floor(time.det))))] +
        table(floor(time.det))

      # determine the time of death of each offspring
      time.death =
        time.exp +
        incubation_periods +
        rlnorm(length(time.exp),meanlog=symp_to_death_meanlog,
               sdlog=symp_to_death_sdlog)
            
      # retain only deaths that occur before time is up
      time.death =
        time.death[floor(time.death) <= stopSimulationDay]    
          
      # update vector of daily mortality
      dailyMortality[as.numeric(names(table(floor(time.death))))] =
        dailyMortality[as.numeric(names(table(floor(time.death))))] +
          table(floor(time.death))
    } 
  }
  # keep track of things that need to be kept track of
  daily.mat[1,] = dailyIncidence
  case.mat[1,] = dailyCases
  death.mat[1,] = dailyMortality    
  cum.vec[1] = sum(dailyIncidence)
  
  # return things that need to be returned
  return(list(daily=daily.mat,death=death.mat,cases=case.mat,cum=cum.vec))
}

## function to simulate autochthonous transmission with a step change in R0
## function to simulate autochthonous transmission
simOutbreakR0Change = function(
  timeImport,
  R = 1.97,
  k = 1e3,
  si_mean = 4.56,
  si_sd = 0.95,
  inc_shape = 1.88,
  inc_scale = 7.97,
  symp_to_death_mean = 14,
  symp_to_death_sd = 5.37,
  report_delay_shape = 3.434,
  report_delay_rate = 0.5724,
  stopSimulationDay = 68,
  asympProp = 0.35,
  asympRFraction = 0.5,
  RChangeDay = Inf,
  RChange = R,
  lnormFlag = T)
{
  # decide what to keep track of across replicate simulations
  daily.mat = matrix(0,1,stopSimulationDay)
    
  case.mat = matrix(0,1,stopSimulationDay)
  death.mat = matrix(0,1,stopSimulationDay)
  
  cum.vec = rep(0, 1)
  
  
  # calculate meanlog and sdlogs for lnorm distributions
  if(lnormFlag){
    si_sdlog = sqrt(log((si_sd/si_mean)^2 + 1))
    si_meanlog = log(si_mean) - 0.5*si_sdlog^2
  }
  
  symp_to_death_sdlog = sqrt(log((symp_to_death_sd/symp_to_death_mean)^2 + 1))
  symp_to_death_meanlog = log(symp_to_death_mean) - 0.5*symp_to_death_sdlog^2
  
  # initiate vector to store daily incidence of autochthonous infections
  dailyIncidence = rep(0,stopSimulationDay)

  # initiate vector to store daily reporting of autochthonous infections
  dailyCases = rep(0,stopSimulationDay)
  
  # initiate vector to store daily mortality of autochthonous infections
  dailyMortality = rep(0,stopSimulationDay)
      
  # initiate vector to store timing of parent infections
  time.exp = timeImport
    
  # loop through generations of transmission until extinct or time is up
  while(length(time.exp) > 0){
    
    # draw a number of offspring for each parent  
    R_vec = ifelse(time.exp < RChangeDay, R, RChange)
    R_vec = ifelse(runif(length(time.exp),0,1) < asympProp, asympRFraction, 1)*R_vec
    number.offspring = rnbinom(n=length(time.exp), mu=R_vec, size=k)
      
        
    # remove parents with zero offspring
    time.exp = time.exp[number.offspring>0]
    number.offspring = number.offspring[number.offspring>0]

    # remove parents exposed after date when transmission stops
    number.offspring = number.offspring[floor(time.exp) < RChangeDay]
    time.exp = time.exp[floor(time.exp) < RChangeDay]

    if(sum(number.offspring > 0)){
      
      # determine the time of exposure of each offspring
      if(!lnormFlag) {
        time.exp =
          rep(time.exp,times=number.offspring) +
          rnorm(sum(number.offspring), mean=si_mean, sd=si_sd)
      } else {
        time.exp =
          rep(time.exp,times=number.offspring) +
          rlnorm(sum(number.offspring), meanlog=si_meanlog, sdlog=si_sdlog)
      }
        
      # retain only those infections that occur before time is up
      time.exp =
        time.exp[floor(time.exp) <= stopSimulationDay]
      # retain only those infections that are exposed before r0 changes
      number.offspring = number.offspring[time.exp <= RChangeDay]
      time.exp = time.exp[floor(time.exp) <= RChangeDay]


      # determine incubation period for each offspring
      incubation_periods = rweibull(length(time.exp), shape=inc_shape, scale=inc_scale)
          
      # determine the time of infection detection of each offspring
      time.det =
        time.exp +
        incubation_periods + 
        rgamma(length(time.exp),
               shape=report_delay_shape, rate=report_delay_rate)
          
      # retain only those detected infections that occur before time is up
      time.det =
        time.det[floor(time.det) <= stopSimulationDay]
  
      # update vector of daily incidence of detected infections      
      dailyIncidence[as.numeric(names(table(floor(time.exp))))] =
        dailyIncidence[as.numeric(names(table(floor(time.exp))))] +
        table(floor(time.exp))
         
      # update vector of daily case reporting
      dailyCases[as.numeric(names(table(floor(time.det))))] =
        dailyCases[as.numeric(names(table(floor(time.det))))] +
        table(floor(time.det))

      # determine the time of death of each offspring
      time.death =
        time.exp +
        incubation_periods +
        rlnorm(length(time.exp),meanlog=symp_to_death_meanlog,
               sdlog=symp_to_death_sdlog)
            
      # retain only deaths that occur before time is up
      time.death =
        time.death[floor(time.death) <= stopSimulationDay]    
          
      # update vector of daily mortality
      dailyMortality[as.numeric(names(table(floor(time.death))))] =
        dailyMortality[as.numeric(names(table(floor(time.death))))] +
          table(floor(time.death))
    } 
  }
  # keep track of things that need to be kept track of
  daily.mat[1,] = dailyIncidence
  case.mat[1,] = dailyCases
  death.mat[1,] = dailyMortality    
  cum.vec[1] = sum(dailyIncidence)
  
  # return things that need to be returned
  return(list(daily=daily.mat,death=death.mat,cases=case.mat,cum=cum.vec))
}
