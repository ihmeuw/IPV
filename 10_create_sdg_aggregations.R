###################################################################################################################################################
#
# Author: AUTHOR
# Purpose: Aggregate and age-standardize  (15-49) IPV SDG indicator estimates
#
###################################################################################################################################################

rm(list=ls())

###################################################################################################################################################
library(data.table)
library(dplyr)
library(parallel)
library(ggplot2)

source("FILEPATH/get_model_results.R")
source("FILEPATH/get_location_metadata.R")
source("FILEPATH/get_age_metadata.R")
source("FILEPATH/get_outputs.R")
source("FILEPATH/get_population.R")
source("FILEPATH/get_covariate_estimates.R")
source("FILEPATH/get_draws.R")
source("FILEPATH/get_bundle_data.R")
source("FILEPATH/get_crosswalk_version.R")
source("FILEPATH/get_bundle_version.R")
source("FILEPATH/interpolate.R")
source("FILEPATH/utility.r")

out_dir <- 'FILEPATH'

age_metadata <- get_age_metadata(gbd_round_id=6, age_group_set_id = 12)
locs <- get_location_metadata(35, gbd_round_id = 7)
location_set <- c(locs[level==3,location_id], 354)
years <- seq(1990, 2019, 1)

# function to make draws long
draws_to_long <- function(dt, value_name = "value"){
  out <- melt.data.table(dt,
                         measure.vars = patterns("draw_"),
                         variable.name = "draw",
                         value.name = value_name)
  
  out[, draw := tstrsplit(draw, "_", keep = 2)]
  out[, draw := as.integer(draw)]
}


get_stgpr_draws <- function(run_id, location_ids=0){
  #' @description My own get draws function so I don't have to upload to EPI upload to get draws for a specific age-sex-loc-yr
  #' @param run_id int. The STGPR run_id that you want to get draws from
  #' @param location_ids numeric. A numeric vector
  #' @return a data.table of all N number of draws from the STGPR output
  
  path <- paste0("/share/covariates/ubcov/model/output/",run_id,"/draws_temp_0/")
  if(location_ids==0){files <- list.files(path)}else{files <- paste0(location_ids, ".csv")}
  if(F){
    # for testing
    files <- paste0("/share/covariates/ubcov/model/output/",run_id,"/draws_temp_0/18.csv")
  }
  
  vars <- paste0("draw_",seq(0,99,1))
  read_draw <- function(file){
    message(paste("Reading", file))
    data <- fread(file)
    message(paste("Done Reading", file))
    data
  }
  
  #' Internal parallelization to read draws from each csv (csvs divided up by location)
  mclapply(paste0(path, files), read_draw, mc.cores = 40) %>% rbindlist(use.names = T)
}


# 1. Read in and melt draws for ipv sdg -------------------------------------------------------------------------------------------------------------------------------------------------------------

#get final prevalence draws from epi
sdg <- get_draws(gbd_id_type='modelable_entity_id',
                 gbd_id=10820,
                 source='epi',
                 measure_id=18,
                 metric_id=3,
                 gbd_round_id = 7,
                 decomp_step='iterative',
                 status='best',
                 year_id=years,
                 age_group_id=seq(8,14,1),
                 location_id = location_set)

dt <- draws_to_long(sdg)
dt <- dt[, .(age_group_id, year_id, sex_id, draw, value, location_id)]

# 2. Get population weight for IPV ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

pop <- get_population(age_group_id = unique(dt$age_group_id),
                      location_id = unique(dt$location_id),
                      sex_id = unique(dt$sex_id),
                      year_id = unique(dt$year_id),
                      gbd_round_id = 6,
                      decomp_step = "step5")

setnames(pop,"population", "weight")
pop[, run_id := NULL]
dt <- merge(dt, pop, all.x=T)

# 3. Create age-std and all age estimates for 15-49 -------------------------------------------------------------------------------------------------------

dt <- merge(dt,age_metadata[,.(age_group_id, age_group_weight_value)], by="age_group_id", all.x=T)

all_age_dt <- dt[,.(value = weighted.mean(value,weight), age_standardized = 0), 
                 by = setdiff(names(dt), c("age_group_id", "value", "weight", "age_group_weight_value", "age_group_name", "age_group_years_start", "age_group_years_end"))]

age_std_dt <- dt[,.(value = weighted.mean(value,age_group_weight_value), age_standardized = 1), 
                 by = setdiff(names(dt), c("age_group_id", "value", "weight", "age_group_weight_value", "age_group_name", "age_group_years_end", "age_group_years_start"))]

dt <- rbind(all_age_dt,age_std_dt)

# 4. Collapse draws to mean/upper/lower -------------------------------------------------------------------------------------------------------------------------------

final <- dt[, .(value = mean(value),
                upper = quantile(value, .975),
                lower = quantile(value, .025)), by=.(location_id, sex_id, year_id, age_standardized)]

final <- merge(final, locs[, c('location_id', 'location_name')], by='location_id')
final[, sex:='Female']
final[, metric_name:='Rate']
final[, measure_name:='Prevalence']
final[, age_group:='15 to 49']
final[, model_name:='SDG Indicator: prevalence of intimate partner physical or sexual violence within the last 12 months (female only)']

write.csv(final, paste0(out_dir, 'sdg_ipv_15to49_aggregated_2020results.csv'))

# 5 year age bins ---------------------------------------------------------------------------------------------------------------------

drawvars <- paste0("draw_",0:999)

sdg_long <- melt.data.table(sdg, id.vars = names(sdg)[!grepl("draw", names(sdg))], measure.vars = patterns("draw"),
                            variable.name = "draw.id", value.name = "draw.value")  

fiveyr <- sdg_long[, .(val = mean(draw.value),
                      upper = quantile(draw.value, .975),
                      lower = quantile(draw.value, .025)), by=.(age_group_id, location_id, sex_id, year_id)]

#create some helpful fields (age group name, location name, measure name, metric name, sex)
fiveyr <- merge(fiveyr, age_metadata[, c('age_group_id', 'age_group_name')])
fiveyr <- merge(fiveyr, locs[, c('location_id', 'location_name')], by='location_id')
fiveyr[, sex:='Female']
fiveyr[, metric_name:='Rate']
fiveyr[, measure_name:='Prevalence']
fiveyr[, model_name:='SDG Indicator: prevalence of intimate partner physical or sexual violence within the last 12 months (female only)']

#Save
write.csv(fiveyr, paste0(out_dir, 'sdg_ipv_5yrbin_2020results.csv'))

#format into one file --------------------------------------------------------------------------------------------------------------------------------------------------

agg <- copy(final)
yr <- copy(fiveyr)
names(agg)
names(yr)

agg[, c('age_group_id', 'sex_id'):=NULL]
setnames(yr, 'age_group_name', 'age_group')
yr[, c('age_group_id', 'sex_id'):=NULL]
setnames(yr, 'val', 'value')
yr[, age_standardized:=0]

#rbind
all <- rbind(agg, yr)

#sort
setorder(all, year_id, location_name, location_id, age_group, age_standardized)

#save
write.csv(all, paste0(out_dir, 'sdg_ipv_2020results.csv'))

#format dataset for GHDx record
ghdx <- fread(paste0(out_dir, 'sdg_ipv_2020results.csv'))[, V1:=NULL]
ghdx[, `:=` (estimate_type='past', indicator_id=1047, indicator_short='Int Partner Viol', 
             ihme_indicator_description='Indicator 5.2.1: Age-standardised prevalence of ever-partnered women aged 15 years and older who experienced physical or sexual violence by a current or former intimate partner in the last 12 months (%)',
             indicator_outline='5.2.1', indicator_unit='0-100%', 
             target_description='Eliminate all forms of violence against all women and girls in the public and private spheres, including trafficking and sexual and other types of exploitation.', 
             goal_description='Achieve gender equality and empower all women and girls.', 
             indicator_scale_type='negative', 
             indicator_target='Equal to or less than 0.5% by 2030')]
values <- c('value', 'upper', 'lower')
ghdx[, (values):=lapply(.SD, function(x) x*100), .SDcols=values]
ghdx <- ghdx[age_standardized==1 & age_group=='15 to 49']
ghdx[, c('age_standardized', 'age_group'):=NULL]

write.csv(ghdx, 'FILEPATH', row.names=F)
