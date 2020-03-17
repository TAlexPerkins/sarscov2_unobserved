# sarscov2_unobserved[a]


This repository contains code to reproduce results from the following pre-print


**Estimating unobserved SARS-CoV-2 infections in the United States**
TA Perkins<sup>&#42;</sup>, SM Cavany<sup>&#42;</sup>, SM Moore<sup>&#42;</sup>, RJ Oidtman, A Lerch, M Poterek
http://perkinslab.weebly.com/uploads/2/5/6/2/25629832/perkins_etal_sarscov2.pdf


### License


This code is being released under the MIT License.
### Overview
Contents of this repository are organized according to code, data, results, and plots. All plots and results shown in the main text and supplementary material are generated in script_main.R and script_supplement.R respectively, with the exception of the sensitivity analyses. Folders for code, data, results, and plots related to the sensitivity analysis described in the supplemental text are located in a sensitivity subfolder within each of those folders.
All code for the analysis featured in the main text was written in the R language (version 3.5.2) and was executed on a Lenovo Thinkpad running Linux (PureOS 9.0). The sensitivity analyses were written in R (version 3.6.2) and run on a Linux cluster. At the time the research was done, all R packages used in this analysis were available on CRAN and were straightforward to install and load.
### Sensitivity analysis
The results and plots for the sensitivity analysis can be recreated by running the following scripts in the code/sensitivity subfolder:
1. script_crc_estimate_parameter_sensitivity.R: This script runs a sweep of the AsympRfraction and parameters for each of 18 different parameter scenarios and calculate log likelihoods.
   1. Input: scenario run number
   2. Output: mle_local_run-number.rda object containing likelihood values
2. parameter_profiles_sensitivity_analysis.R: This script uses the log likelihoods from each parameter sweep, uses bic.grid to create a modified log likelihood surface, and then samples from that surface to get posteriors for the two estimated parameters.
   1. Input: mle_local_run-number.rda
   2. Output: parameter_estimates_posterior_run-number.rda
3. bic_profile_plots.R: Uses output files from step 2 to create hexbin plots of joint estimated parameter posteriors
   1. Input: parameter_estimates_posterior_run-number.rda
4. Script_crc_params_sensitivity_fit. R: Script using the parameter posteriors to simulate infections, deaths, etc for each of 18 scenarios
   1. Input: parameter_estimates_posterior_run-number.rda
   2. Output: sensitivity_sims_run-number.rda
5. combine_sim_results_fit_import.R: Creates summary statistics for cumulative infections and plots for all 18 parameter sensitivity scenarios
   1. Input: sensitivity_sims_run-number.rda
6. combine_sim_results_ratio_deaths.R: Creates summary statistics and plots of deaths compared expected after versus before March 13 for all 18 parameter sensitivity scenarios
   1. Input: sensitivity_sims_run-number.rda
7. multiple_pDetect_ts.R: Creates plots of the probability of detecting a local symptomatic infection over time for all 18 parameter sensitivity scenarios
   1. Input: sensitivity_sims_run-number.rda[b]






[a]Something about like this is what I would like to see


https://github.com/TAlexPerkins/TimeSeriesSpatialScale
[b]You should see what will be required to make it appear how you want it in markdown +smoore15@nd.edu