#setwd("~/Documents/AUMCF_Sim")
rm(list = ls())

# Packages.
library(optparse)
library(MCC)
library(parallel)

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
opt <- make_option(c("--censor"), type = "numeric", help = "Censoring", default = 1)
opt_list <- c(opt_list, opt)

# Frailty variance
opt <- make_option(c("--frailtyVar"), type = "numeric", 
                   help = "Frailty Variance", default = 0.2)
opt_list <- c(opt_list, opt)


# Base death rate for the reference(0) and treatment(1) arm
opt <- make_option(c("--BaseDeath0"), type = "numeric", 
                   help = "Base death rate for the reference arm", default = 0.3)
opt_list <- c(opt_list, opt)

opt <- make_option(c("--BaseDeath1"), type = "numeric",
                   help = "Base death rate for the treatment arm", default = 1)
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
                   help = "Base event rate for the treatment arm", default = 3)
opt_list <- c(opt_list, opt)

# Beta event rate 
opt <- make_option(c("--BetaEvent"), type = "numeric", 
                   help = "Beta death rate", default = 0)
opt_list <- c(opt_list, opt)


# Simulation replicates.
opt <- make_option(c("--reps"), type = "integer", help = "MC replicates", default = 1000)
opt_list <- c(opt_list, opt)

# Permutation replicates.
opt <- make_option(c("--boot"), type = "integer", help = "Bootstrap replicates", default = 500)
opt_list <- c(opt_list, opt)

# Adjustment
opt <- make_option(c("--adjusted"), type = "integer", help = "Indicator of adjustment", default = 0)
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
  "_B", params$boot,
  "_R", params$reps,
  "_Adj", params$adjusted,
  ".rds"
)

# -----------------------------------------------------------------------------
# Simulation.
# -----------------------------------------------------------------------------
Gen_data <- function(params){
  # Check if need to be adjusted
  # NO
  beta_d <- params$BetaDeath
  beta_e <- params$BetaEvent
  covariate0 <- data.frame(arm = rep(0, params$n))
  covariate1 <- data.frame(arm = rep(1, params$n))
  # YES
  if(params$adjusted == 1){
    beta_d <- rep(params$BetaDeath, 2)
    beta_e <- rep(params$BetaEvent, 2)
    covariate0 <- data.frame(arm = rep(0, params$n),
                              covar = stats::rnorm(params$n))
    covariate1 <- data.frame(arm = rep(1, params$n),
                              covar = stats::rnorm(params$n))
  }
  
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
}

Loop <- function(i) {
  data <- Gen_data(params)
  
  boot <- try(
    MCC::CompareAUCs(
      data,
      tau = params$time,
      boot = TRUE,
      reps = params$boot
    )
  )
  
  if (class(boot) != "try-error") {
    out <- c(
      "observed" = boot@CIs$observed[1],
      "asymptotic_se" = boot@CIs$se[1],
      "boot_se" = boot@CIs$se[2]
    )
    return(out)
  }
}

numCores <- detectCores()
#sim <- parallel::mclapply(1:5, Loop, mc.cores = numCores)
sim <- parallel::mclapply(seq_len(params$reps), Loop, mc.cores = numCores)
sim <- do.call(rbind, sim)

# -----------------------------------------------------------------------------
# Summarize.
# -----------------------------------------------------------------------------

# Summarized simulation results.
out <- data.frame(
  "empirical_var" = var(sim[, 1]),
  "asymptotic_var" = mean(sim[, 2]^2),
  "bootstrap_var" = mean(sim[, 3]^2)
)

# Store simulation settings.
out$n <- params$n
out$time <- params$time
out$censor <- params$censor
out$BaseDeath0 <- params$BaseDeath0
out$BaseEvent0 <- params$BaseEvent0
out$BetaDeath <- params$BetaDeath
out$BetaEvent <- params$BetaEvent
out$reps <- params$reps
out$boot <- params$boot


out_stem <- paste0(params$out)
if (!dir.exists(out_stem)) {dir.create(out_stem, recursive = TRUE)}
out_file <- paste0(out_stem, out_suffix)
saveRDS(object = out, file = out_file)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------
t1 <- proc.time()
cat(t1-t0, "\n")
