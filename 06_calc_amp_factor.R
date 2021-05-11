####################################################################################################################################################################
# 
# Author: AUTHOR
# Purpose: Calculate Amp Factor for STGPR based upon amplitude and data variance
#
####################################################################################################################################################################

rm(list=ls())

library(crosswalk, lib.loc = "FILEPATH")
library(openxlsx)
source("FILEPATH/get_crosswalk_version.R")

# set model options -----------------------------------------------------------------------------------------------------------------------------------------------

root <- 'FILEPATH'

run_id <- 179774
xv_id <- 34760
config_path <- paste0('FILEPATH')

# first, load and transform crosswalk version data to calculate nsd ----------------------------------------------------------------------------------------------

# load in data
dt <- as.data.table(read.xlsx(paste0('FILEPATH')))

#fill in standard error before logit transform
dt[is.na(standard_error), standard_error := sqrt(variance)]

#get logit calcs using the delta transform package
logit_means <- as.data.table(delta_transform(mean=dt$val, sd=dt$standard_error, transformation='linear_to_logit'))
setnames(logit_means, c('mean_logit', 'sd_logit'), c('logit_mean', 'logit_se'))
dt <- cbind(dt, logit_means)

#calculate the median of the non-sampling deviation * 1.96
nsd <- 1.96*median(dt[is_outlier!=1, logit_se], na.rm=T)

# function to calculate new amp factor, based on data variance --------------------------------------------------------------------------------------------------------------------------------

source("FILEPATH/utility.r")

amp_factor_calculator <- function(run_id=NULL) {
  
  amp <- model_load(run_id = run_id, "amp_nsv")
  parameters <- model_load(run_id = run_id, 'parameters')
  amp_fac <- parameters[, gpr_amp_factor]
  me_name <- parameters[, me_name]
  print(paste0('me_name = ', me_name))
  
  # the mean amplitude when amp_factor==1
  amp_mean <- mean(amp[,st_amp])/amp_fac
  
  print(paste0('estimated amplitude = ', amp_mean))
  
  # calculate the amp_factor
  print(nsd)
  amp_factor <- nsd/amp_mean
  amp.list <- list("amp_factor"=amp_factor, "nsd"=nsd)
  
  print(paste0('new_amp_factor = ', amp_factor))
  
  return(amp.list)
}

# run the function ------------------------------------------------------------------------------------------------------------------------------------------------------

amp_factor <- amp_factor_calculator(run_id)

# read in config and change amp_factor value accordingly -----------------------------------------------------------------------------------------------------------------

config <- fread(paste0(config_path, '.csv'))
config[, gpr_amp_factor:=amp_factor$amp_factor]
config[, crosswalk_version_id:=xv_id]
config[, description:='SDG STGPR run with calc amp factor & gpr_scale=25, 3.1.21']
write.csv(config, paste0(config_path, '_calculated_ampfactor.csv'), row.names=F)


