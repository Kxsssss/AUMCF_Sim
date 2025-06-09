#setwd("~/Documents/AUMCF_Sim/jacc")
setwd('~/Documents/GitHub/AUMCF_Sim')
rm(list = ls())

# Packages.
library(optparse)
library(MCC)
library(parallel)
library(dplyr)

#source("../JACC_methods.R")
source("JACC_methods.R")

my_reps <- 2000
params <- list(
  n = 200,
  time = 4,
  censor = 0.2,
  frailtyVar = 6,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.4,
  BaseEvent1 = 1.0,
  BetaDeath = 0,
  BetaEvent = 0,
  reps = my_reps,
  adjusted = 0,
  tvr = 1.2,
  tvd = 0.2,
  experiment = 2,
  out = "Test/"
)
# Option parsing.
t0 <- proc.time()
# Output stem.
out_suffix <- paste0(
  "N", params$n, 
  "_T", params$time,
  "_l0", params$BaseEvent0,
  "_l1", params$BaseEvent1,
  "_adj", params$adjusted,
  ".rds"
)

source('test_sim_comparison.R')

summarize_results <- function(sim_augmented){
  cbind(sim_augmented %>%
          group_by(type) %>%
          summarise(Pt_est = mean(value, na.rm = TRUE)),
        sim_augmented %>%
          group_by(type) %>%
          summarise(Prob_reject_H0 = mean(p_value < 0.05, na.rm = TRUE)), 
        sim_augmented %>%
          group_by(type) %>%
          summarise(ASE = mean(se, na.rm = TRUE)),
        sim_augmented %>%
          group_by(type) %>%
          summarise(ESE = sd(value, na.rm = TRUE)),
        #sim_augmented %>%
         # group_by(type) %>%
          # summarise(CP = mean( (lower < 1/1.4) * (upper > 1/1.4), na.rm = TRUE)),
        sim_augmented %>%
          group_by(type) %>%
          summarize(n_NA = sum(is.na(value) | is.na(se) | is.na(p_value)))
  )[, -c(3,5,7,9)]
}


print(dim(sim_augmented)) 
print(summarize_results(sim_augmented))





# setting for different death rates
# params <- list(
#   n = 100,
#   time = 4,
#   censor = 0.2,
#   frailtyVar = 6,
#   BaseDeath0 = 0.5,
#   BaseDeath1 = 0.1,
#   BaseEvent0 = 1.0,
#   BaseEvent1 = 1.0,
#   BetaDeath = 0,
#   BetaEvent = 0,
#   reps = 500,
#   adjusted = 0,
#   tvr = 1.2,
#   tvd = 0.2,
#   experiment = 3,
#   out = "Test/"
# )
