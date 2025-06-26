setwd('~/Documents/GitHub/AUMCF_Sim')
rm(list = ls())

# Packages.
library(optparse)
library(MCC)
library(parallel)
library(dplyr)

# Functions for comparison methods.
source("comparison_methods.R")

# Functions for data generating with time-varying treatment effect.
# Note: Not currently in MCC package.
source("data_generation.R")

# Test for Case 1a. 
my_reps <- 100

params <- list(
  n = 200, 
  censor = 0.2,
  time = 4,
  frailtyVar = 0,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.0,
  BaseEvent1 = 1.0,
  BetaDeath = 0,
  BetaEvent = 0, 
  adjusted = 0,
  TV_effect = 0,
  experiment = 1,
  reps = my_reps,
  out = "Test/"
)
source('run_simulation.R')

# Test for Case 3a. 
my_reps <- 500

params <- list(
  n = 200, 
  censor = 0.2,
  time = 4,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  TV_effect = 0,
  adjusted = 0,
  experiment = 3,
  reps = my_reps,
  out = "Test/"
)
source('run_simulation.R')

# Test for Case 3b. 
my_reps <- 500

params <- list(
  n = 200, 
  censor = 0.2,
  time = 4,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  TV_effect = log(0.6),
  adjusted = 0,
  experiment = 3,
  reps = my_reps,
  out = "Test/"
)
source('run_simulation.R')

# Results for this run. 
# type        bias prob_reject_H0 cov_p        ase        ese
# 1   aucmf_diff 0.041782368          0.692 0.954 0.49328382 0.48764647
# 2  aucmf_ratio 0.011900395          0.694 0.958 0.07125965 0.07042072
# 3        cox_p 0.040281840          0.122 0.926 0.10973810 0.10574539
# 4      frailty 0.009709693          0.836 0.948 0.07769361 0.06260279
# 5         lwyy 0.009503208          0.932 0.956 0.07593875 0.05838590
# 6           nb 0.010281907          0.920 0.960 0.07742137 0.05835420
# 7       wr_FWR 0.234321079          0.268 0.732 0.18535262 0.19364306
# 8       wr_LWR 0.236309103          0.272 0.728 0.18438422 0.19146985
# 9       wr_NWR 0.331599309          0.344 0.656 0.23137935 0.23931787
# 10      wr_STD 0.036698754          0.052 0.948 0.16231245 0.16923089
# p_value true_value
# 1  0.07837816 -1.2436419
# 2  0.07810764  0.7982747
# 3  0.40054231  0.8846132
# 4  0.03748982  0.7894341
# 5  0.01739252  0.7609708
# 6  0.01880018  0.7613297
# 7  0.27272200  1.0000000
# 8  0.26870924  1.0000000
# 9  0.22302242  1.0000000
# 10 0.50110607  1.0000000