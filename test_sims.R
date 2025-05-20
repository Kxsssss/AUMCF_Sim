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

# A: NULL SCENARIO
# A1: High censoring, tau = 4 
params <- list(
  n = 100,
  time = 4,
  censor = 0.4,
  frailtyVar = 0,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.0,
  BaseEvent1 = 1.0,
  BetaDeath = 0,
  BetaEvent = 0,
  reps = 5000,
  adjusted = 0,
  tvr = 1.2,
  tvd = 0.2,
  experiment = 3,
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


# A2: Low censoring, tau = 4 
params <- list(
  n = 100,
  time = 4,
  censor = 0.2,
  frailtyVar = 0,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.0,
  BaseEvent1 = 1.0,
  BetaDeath = 0,
  BetaEvent = 0,
  reps = 5000,
  adjusted = 0,
  tvr = 1.2,
  tvd = 0.2,
  experiment = 3,
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


# A3: High censoring, tau = 1 
params <- list(
  n = 100,
  time = 1,
  censor = 0.4,
  frailtyVar = 0,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.0,
  BaseEvent1 = 1.0,
  BetaDeath = 0,
  BetaEvent = 0,
  reps = 5000,
  adjusted = 0,
  tvr = 1.2,
  tvd = 0.2,
  experiment = 3,
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


# A4: Low censoring, tau = 1 
params <- list(
  n = 100,
  time = 1,
  censor = 0.2,
  frailtyVar = 0,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.0,
  BaseEvent1 = 1.0,
  BetaDeath = 0,
  BetaEvent = 0,
  reps = 5000,
  adjusted = 0,
  tvr = 1.2,
  tvd = 0.2,
  experiment = 3,
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


# B: NON-NULL SCENARIO
# B1: High censoring, tau = 4 
params <- list(
  n = 100,
  time = 4,
  censor = 0.4,
  frailtyVar = 0,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.4,
  BaseEvent1 = 1.0,
  BetaDeath = 0,
  BetaEvent = 0,
  reps = 5000,
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


# B2: Low censoring, tau = 4 
params <- list(
  n = 100,
  time = 4,
  censor = 0.2,
  frailtyVar = 0,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.4,
  BaseEvent1 = 1.0,
  BetaDeath = 0,
  BetaEvent = 0,
  reps = 5000,
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


# B3: High censoring, tau = 1 
params <- list(
  n = 100,
  time = 1,
  censor = 0.4,
  frailtyVar = 0,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.4,
  BaseEvent1 = 1.0,
  BetaDeath = 0,
  BetaEvent = 0,
  reps = 5000,
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


# B4: Low censoring, tau = 1 
params <- list(
  n = 100,
  time = 4,
  censor = 0.2,
  frailtyVar = 0,
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.0,
  BaseEvent1 = 1.0,
  BetaDeath = 0,
  BetaEvent = 0,
  reps = 500,
  adjusted = 0,
  tvr = 1.2,
  tvd = 0.2,
  experiment = 3,
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



