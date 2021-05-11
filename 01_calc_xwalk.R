####################################################################################################################################################################
# 
# Author: USERNAME
# Purpose: Crosswalking Network Meta-analysis for IPV alternate definitions - Sub definition approach
#
####################################################################################################################################################################

rm(list=ls())

# set-up environment -----------------------------------------------------------------------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(crosswalk, lib.loc = "FILEPATH")

source("FILEPATH/get_bundle_data.R")
source("FILEPATH/get_bundle_version.R")
source("FILEPATH/save_bundle_version.R")
source("FILEPATH/save_crosswalk_version.R")
source("FILEPATH/get_age_metadata.R")
source("FILEPATH/upload_bundle_data.R")

## Function to fill out mean/cases/sample sizes
get_cases_sample_size <- function(raw_dt){
  dt <- copy(raw_dt)
  dt[is.na(mean), mean := cases/sample_size]
  dt[is.na(cases) & !is.na(sample_size), cases := mean * sample_size]
  dt[is.na(sample_size) & !is.na(cases), sample_size := cases / mean]
  return(dt)
}

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

# Read in data and get into correct matched format for within-study pairs ------------------------------------------------------------------------------------------

#Read in collapsed microdata
data <- get_bundle_data(8921, 'iterative', 7)

#filter to 2020 data only
old <- read.xlsx('FILEPATH')
data <- data[!nid %in% unique(old$nid)]

#sex info
data[, sex_id:=2]

#age info
data[age_start==999, age_start:=age_start_orig]
data[age_end==999, age_end:=age_end_orig]
data[, age_midpt:=(as.numeric(age_start)+as.numeric(age_end))/2]

# Set zero mean adjustment value depending on version: 
mean_adj_val <- quantile(data[val!=0]$val, probs=0.025)

#renaming
setnames(data, c('standard_error', 'val'), c('se', 'mean'), skip_absent = T)

#remove unnecessary columns
data <- data[, c('age_start', 'age_end', 'age_midpt', 'ihme_loc_id', 'nid', 'year_id', 'survey_name', 'var', 'mean', 'se')]

data_orig <- copy(data) #save an unmatched version of the data to apply adjustments to 

# Drop subnat locations
data <- data[!ihme_loc_id %like% '_']

#var subsetting
data <- data[var %in% c('any_sexpv_pastyr', 'any_physpv_pastyr', 'anysex_anyphys_ipv_pastyr')]

#get subsets of gs and alternate defs
ref_subset <- copy(data) #for network, refvar can be any of the alternative definitions
alts_subset <- data[!var=='anysex_anyphys_ipv_pastyr'] #altvar is any of the alternative defs only

#set names accordingly
setnames(ref_subset, c('mean', 'se', 'var'), c('ref_mean', 'ref_se', 'refvar'))
setnames(alts_subset, c('mean', 'se', 'var'), c('alt_mean', 'alt_se', 'altvar'))

#merge back onto each other
matched <- merge(ref_subset, alts_subset, by=c('age_start', 'age_end', 'ihme_loc_id', 'nid', 'year_id', 'survey_name', 'age_midpt'), allow.cartesian = T)

#get rid of rows that match the same def
onepoint_data <- matched[altvar==refvar]
matched <- matched[!altvar==refvar]

#remove duplicate indirect comparisons (B:C == C:B)
data <- copy(matched)

alt_defs <- unique(alts_subset$altvar)

for (i in 1:length(alt_defs)){
  for (j in 1:length(alt_defs)){
    data[refvar==alt_defs[i] & altvar==alt_defs[j], comparison_pair:=paste0(alt_defs[i], ' to ', alt_defs[j])]
    data[refvar==alt_defs[j] & altvar==alt_defs[i], comparison_pair:=paste0(alt_defs[i], ' to ', alt_defs[j])]
  }
}

data[!is.na(comparison_pair), duplicate_pair:=duplicated(comparison_pair), by=c('nid', 'age_start', 'age_end')]
data <- data[is.na(duplicate_pair) | duplicate_pair==FALSE]

data[, comparison_pair:=NULL]
data[, duplicate_pair:=NULL]

matched <- copy(data)

# Change zeros to adj_value set above
matched <- matched[ref_mean<=mean_adj_val, ref_mean:=mean_adj_val]
matched <- matched[alt_mean<=mean_adj_val, alt_mean:=mean_adj_val]

#get logit calcs using the delta transform package
logit_alt_means <- as.data.table(delta_transform(mean=matched$alt_mean, sd=matched$alt_se, transformation='linear_to_logit'))
setnames(logit_alt_means, c('mean_logit', 'sd_logit'), c('logit_alt_mean', 'logit_alt_se'))
logit_ref_means <- as.data.table(delta_transform(mean=matched$ref_mean, sd=matched$ref_se, transformation='linear_to_logit'))
setnames(logit_ref_means, c('mean_logit', 'sd_logit'), c('logit_ref_mean', 'logit_ref_se'))

#bind back onto main data table
matched <- cbind(matched, logit_alt_means)
matched <- cbind(matched, logit_ref_means)

matched[, c("logit_diff", "logit_diff_se")] <- calculate_diff(
  df = matched, 
  alt_mean = "logit_alt_mean", alt_sd = "logit_alt_se",
  ref_mean = "logit_ref_mean", ref_sd = "logit_alt_se" )

# Specify the data for model fitting -------------------------------------------------------------------------------------------------------------------------------------------------------

df <- CWData(
  df = matched,             # dataset for metaregression
  obs = "logit_diff",       # column name for the observation mean
  obs_se = "logit_diff_se", # column name for the observation standard error
  alt_dorms = "altvar",     # column name of the variable indicating the alternative method
  ref_dorms = "refvar",     # column name of the variable indicating the reference method
  covs = list(),
  study_id = 'nid',          # name of the column indicating group membership, usually the matching groups
  add_intercept = TRUE
)

# Fit the model ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

fit_ipv_network <- CWModel(
  cwdata = df,            # object returned by `CWData()`
  obs_type = "diff_logit", # "diff_log" or "diff_logit" depending on whether bounds are [0, Inf) or [0, 1]
  cov_models = list(       # specifying predictors in the model; see help(CovModel)
    CovModel("intercept")),
  gold_dorm = "anysex_anyphys_ipv_pastyr",  # the level of `ref_dorms` that indicates it's the gold standard
  order_prior = list(c('any_sexpv_pastyr', 'anysex_anyphys_ipv_pastyr'), 
                     c('any_physpv_pastyr', 'anysex_anyphys_ipv_pastyr')),
  max_iter=100L, #default is 100
  inlier_pct=0.9 #set trimming to 10%
)

# Save results ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

version <- 'final_noagecov'
py_save_object(object = fit_ipv_network, filename = 'FILEPATH', pickle = "dill")


