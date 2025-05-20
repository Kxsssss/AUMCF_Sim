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
# -----------------------------------------------------------------------------
# Command line arguments.
# -----------------------------------------------------------------------------

# Command line options.
opt_list <- list()

# Sample size.
opt <- make_option(c("--n"), type = "integer", help = "Patients", default = 100)
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
                   help = "Beta death rate", default = 1)
opt_list <- c(opt_list, opt)


# Simulation replicates.
opt <- make_option(c("--reps"), type = "integer", help = "MC replicates", default = 10)
opt_list <- c(opt_list, opt)

# Adjustment
opt <- make_option(c("--adjusted"), type = "integer", help = "Indicator of adjustment", default = 1)
opt_list <- c(opt_list, opt)

# True value 
opt <- make_option(c("--tvr"), type = "numeric", help = "True value (ratio)", default = 1)
opt_list <- c(opt_list, opt)
opt <- make_option(c("--tvd"), type = "numeric", help = "True value (difference)", default = 0)
opt_list <- c(opt_list, opt)

# Experience index
opt <- make_option(c("--experiment"), type = "integer", help = "Index of experience", default = 3)
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
  "_l0", params$BaseEvent0,
  "_l1", params$BaseEvent1,
  "_adj", params$adjusted,
  ".rds"
)

# -----------------------------------------------------------------------------
# Simulation.
# -----------------------------------------------------------------------------
# Data generation for unadjusted cases
Gen_data <- function(params){
  # E1, E4
  if(params$experiment == 1 | params$experiment == 4){
    beta_d <- params$BetaDeath
    beta_e <- params$BetaEvent
    covariate <- data.frame(arm = c(rep(0, params$n), rep(1, params$n)))
  }
  # E3
  if(params$experiment == 3){
    # if BetaEvent is set to be 0, then run the null case
    if(params$BetaEvent == 0){
      beta_e <- c(log(1), log(1))
      beta_d <- c(log(1), log(1))
    }else{
      beta_e <- c(0, log(2))
      beta_d <- c(0, log(0.5))
    }
    covariate <- data.frame(arm = c(rep(0, params$n), rep(1, params$n)),
                            covar = stats::rnorm(params$n*2))
  }
  
  if(params$experiment == 2){
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
  set.seed(i)
  data <- Gen_data(params)
  boot <- try(
      MCC::CompareAUCs(data, tau = params$time)
    )
  if (class(boot) != "try-error") {
    # aucmf
    aucmf <- data.frame(
      value = boot@CIs$observed[2],
      se = boot@CIs$se[2],
      lower = boot@CIs$lower[2],
      upper = boot@CIs$upper[2],
      p_value = boot@Pvals$p[2],
      #z_value = (boot@CIs$observed[2]-1)/boot@CIs$se[2],
      type = "aucmf_ratio"
    )
    
    coxp <- cox_prop(data)
    lwyy <- lwyy(data)
    nb <- nb(data)
    frailty <- frailty_(data)
    wr <- wr(data)
    
    aucmf_diff <- data.frame(
      value = boot@CIs$observed[1],
      se = boot@CIs$se[1],
      lower = boot@CIs$lower[1],
      upper = boot@CIs$upper[1],
      p_value = boot@Pvals$p[1],
      type = "aucmf_diff"
    )
    
    results <- rbind(aucmf, coxp, lwyy, nb, frailty, wr, aucmf_diff)
    
    # if need to compare to the adjusted case 
    if(params$adjusted == 1){
      boot_adj <- try(
        MCC::CompareAUCs(data, tau = params$time, 
                         covars = data %>% dplyr::select(covar))
      )
      if (class(boot_adj) != "try-error") {
        aucmf_diff_adj <- data.frame(
          value = boot_adj@CIs$observed[1],
          se = boot_adj@CIs$se[1],
          lower = boot_adj@CIs$lower[1],
          upper = boot_adj@CIs$upper[1],
          p_value = boot_adj@Pvals$p[1],
          type = "aucmf_diff_adj"
        )
        results <- rbind(results, aucmf_diff_adj)
        }
      }
    rownames(results) <- 1:nrow(results)
    
    return(results)
  }
}

numCores <- detectCores()
#sim_l <- parallel::mclapply(seq_len(params$reps), Loop, mc.cores = numCores)
sim_l <- parallel::mclapply(seq_len(5), Loop, mc.cores = numCores)
sim_ul <- do.call(rbind, sim_l)
#sim$rep <- rep(seq_along(sim_l), each = nrow(sim_l[[1]]))

char_cols <- c("type")
sim <- sim_ul %>%
  mutate(across(
    .cols = -all_of(char_cols),
    .fns = ~ suppressWarnings(as.numeric(.))
  )) %>%
  filter(complete.cases(.))

extra_run <- 0
baseline <- if_else(params$adjusted == 0, 7, 8)
extra_i <- 0
while (nrow(sim) < baseline*params$reps) {
  extra_run <- (baseline*params$reps - nrow(sim))/baseline
  sim_l1 <- parallel::mclapply(seq_len(extra_run) + params$reps + extra_i,
                               Loop, mc.cores = numCores)
  sim_ul1 <- do.call(rbind, sim_l1)
  
  char_cols <- c("type")
  sim1 <- sim_ul1 %>%
    mutate(across(
      .cols = -all_of(char_cols),
      .fns = ~ suppressWarnings(as.numeric(.))
    )) %>%
    filter(complete.cases(.))
  sim <- rbind(sim, sim1)
  extra_i = extra_i + extra_run
}

# -----------------------------------------------------------------------------
# Summarize.
# -----------------------------------------------------------------------------
# Summarized simulation results.
tv_lookup <- tibble(
  type = unique(sim$type),
  true_value = if_else(type == "aucmf_diff" | type == "aucmf_diff_adj", 
                       params$tvd, params$tvr)
)

sim_augmented <- sim %>%
  left_join(tv_lookup, by = "type")

summary_table <- sim_augmented %>%
  group_by(type) %>%
  summarise(
    bias = mean(value) - first(true_value),
    lower = mean(lower),
    upper = mean(upper),
    cov_p = mean(lower <= first(true_value) & first(true_value) <= upper),
    ase = mean(se),
    ese = sd(value),
    p_value = mean(p_value),
    .groups = "drop"
  )

summary_table_tv <- summary_table %>%
  left_join(tv_lookup, by = "type")

out <- data.frame(summary_table_tv)
# Store simulation settings.
out$n <- params$n
out$time <- params$time
out$rep <- params$reps

# save the summary table
out_stem <- paste0(params$out)
if (!dir.exists(out_stem)) {dir.create(out_stem, recursive = TRUE)}
out_file <- paste0(out_stem, out_suffix)
saveRDS(object = out, file = out_file)

# save the points estimates
sim_file <- paste0(out_stem, "sim",
                   "N", params$n, 
                   "_T", params$time,
                   "_l0", params$BaseEvent0,
                   "_l1", params$BaseEvent1,
                   ".rds")
saveRDS(object = sim_augmented, file = sim_file)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------
t1 <- proc.time()
cat(t1-t0, "\n")
