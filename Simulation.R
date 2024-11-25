#setwd("~/Documents/AUMCF_Sim")
rm(list = ls())

# Packages.
library(optparse)
library(MCC)
library(parallel)
library(dplyr)

# -----------------------------------------------------------------------------
# Command line arguments.
# -----------------------------------------------------------------------------

# Command line options.
opt_list <- list()

# Sample size.
opt <- make_option(c("--n"), type = "integer", help = "Patients", default = 50)
opt_list <- c(opt_list, opt)

# Truncation time (tau).
opt <- make_option(c("--time"), type = "numeric", help = "Patients", default = 2)
opt_list <- c(opt_list, opt)

# Censoring rate.
opt <- make_option(c("--censor"), type = "numeric", help = "Censoring", default = 0.2)
opt_list <- c(opt_list, opt)

# Frailty variance
opt <- make_option(c("--frailtyVar"), type = "numeric", 
                   help = "Frailty Variance", default = 0)
opt_list <- c(opt_list, opt)


# Base death rate for the reference(0) and treatment(1) arm
opt <- make_option(c("--BaseDeath0"), type = "numeric", 
                   help = "Base death rate for the reference arm", default = 0.2)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--BaseDeath1"), type = "numeric",
                   help = "Base death rate for the treatment arm", default = 0.2)
opt_list <- c(opt_list, opt)

# Beta death rate 
opt <- make_option(c("--BetaDeath"), type = "numeric", 
                   help = "Beta death rate", default = 0)
opt_list <- c(opt_list, opt)


# Base event rate for the reference(0) and treatment(1) arm
opt <- make_option(c("--BaseEvent0"), type = "numeric", 
                   help = "Base event rate for the reference arm", default = 1)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--BaseEvent1"), type = "numeric", 
                   help = "Base event rate for the treatment arm", default = 1)
opt_list <- c(opt_list, opt)

# Beta event rate 
opt <- make_option(c("--BetaEvent"), type = "numeric", 
                   help = "Beta death rate", default = 0)
opt_list <- c(opt_list, opt)


# Simulation replicates.
opt <- make_option(c("--reps"), type = "integer", help = "MC replicates", default = 1000)
opt_list <- c(opt_list, opt)

# Adjustment
opt <- make_option(c("--adjusted"), type = "integer", help = "Indicator of adjustment", default = 0)
opt_list <- c(opt_list, opt)

# True value 
opt <- make_option(c("--tv"), type = "numeric", help = "True value", default = 0)
opt_list <- c(opt_list, opt)

# Experience index
opt <- make_option(c("--ei"), type = "integer", help = "Index of experience", default = 1)
opt_list <- c(opt_list, opt)

# Output directory.
opt <- make_option(c("--out"), type = "character", help = "Output stem", default = "Test/")
opt_list <- c(opt_list, opt)

# Option parsing.
t0 <- proc.time()
parsed_opts <- OptionParser(option_list = opt_list)
params <- parse_args(object = parsed_opts)

# Output stem.
out_suffix <- paste0(
  "N", params$n, 
  "_T", params$time,
  "_f", params$frailtyVar,
  "_bb",params$BetaDeath, params$BetaEvent,
  "_Adj", params$adjusted,
  ".rds"
)

# -----------------------------------------------------------------------------
# Simulation.
# -----------------------------------------------------------------------------
# Data generation for unadjusted cases
Gen_data <- function(params){
  # E1, E4
  if(params$ei == 1 | params$ei == 4){
    beta_d <- params$BetaDeath
    beta_e <- params$BetaEvent
    covariate <- data.frame(arm = c(rep(0, params$n), rep(1, params$n)))
  }
  # E3
  if(params$ei == 3){
    # if BetaEvent is set to be 0, then run the null case
    if(params$BetaEvent == 0){
      beta_e <- c(log(1), log(1))
    }else{
      beta_e <- c(log(0.5), log(2))
    }
    beta_d <- c(0,0)
    covariate <- data.frame(arm = c(rep(0, params$n), rep(1, params$n)),
                            covar = stats::rnorm(params$n*2))
  }
  
  if(params$ei == 2){
    beta_d <- params$BetaDeath
    beta_e <- params$BetaEvent
    covariate0 <- data.frame(arm = rep(0, params$n))
    covariate1 <- data.frame(arm = rep(1, params$n))
    
    # Generate data for each arm based on different based death rate  
    # and base event rate
    data0 <- MCC::GenData(
      n = params$n,
      censoring_rate = params$censor,
      base_death_rate = params$BaseDeath0,
      base_event_rate = params$BaseEvent0,
      beta_death = beta_d,
      beta_event = beta_e,
      covariates = covariate0,
      frailty_variance = params$frailtyVar,
      tau = params$time
    )
    
    data1 <- MCC::GenData(
      n = params$n,
      censoring_rate = params$censor,
      base_death_rate = params$BaseDeath1,
      base_event_rate = params$BaseEvent1,
      beta_death = beta_d,
      beta_event = beta_e,
      covariates = covariate1,
      frailty_variance = params$frailtyVar,
      tau = params$time
    )
    
    # Combine the data for control arm and treatment arm
    data1$idx <- data1$idx + params$n
    data <- rbind(data0, data1)
    return(data)
    
  }else{
    data <- MCC::GenData(
      n = params$n * 2,
      censoring_rate = params$censor,
      base_death_rate = params$BaseDeath0,
      base_event_rate = params$BaseEvent0,
      beta_death = beta_d,
      beta_event = beta_e,
      covariates = covariate,
      frailty_variance = params$frailtyVar,
      tau = params$time
    )
    return(data)
  }
}

Loop <- function(i) {
  data <- Gen_data(params)
  # if adjustment is needed
  if(params$adjusted == 1){
    boot <- try(
      MCC::CompareAUCs(data, tau = params$time, 
                       covars = data %>% dplyr::select(covar))
    )
  }else{ # others
    boot <- try(
      MCC::CompareAUCs(data, tau = params$time)
    )
  }
  
  if (class(boot) != "try-error") {
    sim_results <- c(
      "n" = params$n,
      "tau" = params$time,
      "observed" = boot@CIs$observed[1],
      "asymptotic_se" = boot@CIs$se[1],
      "cover_ind" = boot@CIs$lower[1] <= params$tv & params$tv <= boot@CIs$upper[1],
      "p_value" = boot@Pvals$p[1]
    )
    return(sim_results)
  }
}

numCores <- detectCores()
sim <- parallel::mclapply(seq_len(params$reps), Loop, mc.cores = numCores)
sim <- do.call(rbind, sim)

# -----------------------------------------------------------------------------
# Summarize.
# -----------------------------------------------------------------------------
# Summarized simulation results.
out <- data.frame(
  "bias" = mean(sim[, 3]) - params$tv,
  "empirical_se" = sd(sim[, 3]),
  "asymptotic_se" = mean(sim[, 4]),
  "CP" = mean(sim[,5]),     # coverage probability
  "p_value"= mean(sim[,6]),
  "MSE" = mean(sum((sim[, 3] - params$tv)^2))
)

# Store simulation settings.
out$n <- params$n
out$time <- params$time

# save the summary table
out_stem <- paste0(params$out)
if (!dir.exists(out_stem)) {dir.create(out_stem, recursive = TRUE)}
out_file <- paste0(out_stem, out_suffix)
saveRDS(object = out, file = out_file)
# save the points estimates
sim_file <- paste0(out_stem, "sim",
                   "_N", params$n,
                   "_T", params$time,
                   "_Adj", params$adjusted,
                   "_bb", params$BetaDeath, params$BetaEvent,
                   "_f", params$frailtyVar,
                   ".rds")
saveRDS(object = sim, file = sim_file)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------
t1 <- proc.time()
cat(t1-t0, "\n")
