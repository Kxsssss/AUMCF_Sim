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

params <- list(
  n = 100,
  time = 5,
  censor = 0.5,
  frailtyVar = 1,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.6,
  BaseEvent0 = 1.0,
  BaseEvent1 = 1.0,
  BetaDeath = 0,
  BetaEvent = 0,
  reps = 500,
  adjusted = 0,
  tvr = 1.2,
  tvd = 0.2,
  experiment = 1,
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


print(dim(sim_augmented)) # should be 2000 * 9 = 18000


print(cbind(sim_augmented %>%
              group_by(type) %>%
              summarise(Pt_est = mean(value)),
  sim_augmented %>%
              group_by(type) %>%
              summarise(Prob_reject_H0 = mean(p_value < 0.05)), 
  sim_augmented %>%
        group_by(type) %>%
        summarise(ASE = mean(se)),
        sim_augmented %>%
          group_by(type) %>%
          summarise(ESE = sd(value))
        )[, -c(3,5,7)])
