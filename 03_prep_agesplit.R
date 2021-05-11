####################################################################################################################################################################
# 
# Author: AUTHOR
# Purpose: Prepare SDG IPV data for agesplit: 
#          (1) Split into train and test data 
#          (2) Launch ST-GPR model on test data
#
####################################################################################################################################################################

rm(list=ls())

library(data.table)
library(dplyr)
library(openxlsx)

# source central functions
source("FILEPATH/save_crosswalk_version.R")
source("FILEPATH/get_bundle_version.R")
source("FILEPATH/get_age_metadata.R")
source("FILEPATH/get_location_metadata.R")
source("FILEPATH/get_population.R")
library(crosswalk, lib.loc = "FILEPATH")
library(logitnorm, lib.loc='FILEPATH')

bv <- 37034

## Function to calculate standard error
get_se <- function(raw_dt){
  dt <- copy(raw_dt)
  dt[is.na(standard_error) & !is.na(lower) & !is.na(upper), standard_error := (upper-lower)/3.92]
  z <- qnorm(0.975)
  dt[is.na(standard_error) & measure %in% c("prevalence","proportion"), standard_error := sqrt(mean*(1-mean)/sample_size + z^2/(4*sample_size^2))]
  dt[is.na(standard_error) & measure == "incidence" & cases < 5, standard_error := ((5-mean*sample_size)/sample_size+mean*sample_size*sqrt(5/sample_size^2))/5]
  dt[is.na(standard_error) & measure == "incidence" & cases >= 5, standard_error := sqrt(mean/sample_size)]
  return(dt)
}

## Function to fill out mean/cases/sample sizes
get_cases_sample_size <- function(raw_dt){
  dt <- copy(raw_dt)
  dt[is.na(mean), mean := cases/sample_size]
  dt[is.na(cases) & !is.na(sample_size), cases := mean * sample_size]
  dt[is.na(sample_size) & !is.na(cases), sample_size := cases / mean]
  return(dt)
}

## Function to get cases if they are missing
calculate_cases_fromse <- function(raw_dt){
  dt <- copy(raw_dt)
  dt[is.na(cases) & is.na(sample_size) & measure %in% c("prevalence","proportion"), sample_size := (mean*(1-mean)/standard_error^2)]
  dt[is.na(cases) & is.na(sample_size) & measure == "incidence", sample_size := mean/standard_error^2]
  dt[is.na(cases), cases := mean * sample_size]
  return(dt)
}
# read in input data ----------------------------------------------------------------------------------------------------------------------------------------------

data <- as.data.table(read.xlsx(paste0('FILEPATH')))

ages <- get_age_metadata(19)
ages[, age_group_years_end:=age_group_years_end-1]
ages[age_group_id==235, age_group_years_end:=99]
setnames(ages, c('age_group_years_start', 'age_group_years_end'), c('age_start', 'age_end'))
locs <- get_location_metadata(22)

# divide data into test and train to be used for ST-GPR --------------------------------------------------------------------------------------------------------------------

#formatting
setnames(data, 'mean', 'val')
data[, measure_id:=18]
data[, sex_id:=2]
data[, variance:=standard_error^2]

#split into train and test datasets
train <- data[age_group_id!=999]
test <- fsetdiff(data, train)

#subset to data added in 2020 only for training set
train <- train[is.na(cv_sdg_ipv_no_physical)]

#format
col.stgpr <- c("nid","location_id","year_id","age_group_id","sex_id","val","variance","sample_size","measure","is_outlier", "measure_id")
train.in.stgpr <- train[,col.stgpr,with = F]
train.in.stgpr[age_group_id==30, age_group_id:=21]

#write datasets
write.csv(train.in.stgpr, 'FILEPATH')
write.csv(test, 'FILEPATH')

# run ST-GPR model on the train data -----------------------------------------------------------------------------------------------------------------------------

central_root <- 'FILEPATH'
setwd(central_root)
source('FILEPATH/register.R')
source('FILEPATH/sendoff.R')

run_id <- register_stgpr_model('FILEPATH')
stgpr_sendoff(run_id, 'proj_team', log_path = 'FILEPATH')

