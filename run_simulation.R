# -----------------------------------------------------------------------------
# Main script to run the simulation settings.
# -----------------------------------------------------------------------------

# Convenience function for data generation.

# Generates data under the scenarios in the paper, including:
  ## Case 1: Constant treatment effect (Under the Null).
    ## Case 1a: No covariate effect.
    ## Case 1b: Covariate effect.
  ## Case 2: Constant treatment effect (Under the Alternative).
  ## Case 3: Time-varying treatment effect.
    ## Case 3a: Under the Null.
    ## Case 3b: Under the Alternative.

SimulateData <- function(params, calc_truth = FALSE){

  
  if (!calc_truth){
    
    # Keep censoring and n rate if not calculating the true parameter value.
    censoring_rate <- params$censor
    n <- params$n
    
  }else{
    
    # Set censoring rate to 0 and increase n if calculating the true parameter value.
    censoring_rate <- 0
    n <- 10000
    
  }
  
  
  if(params$experiment == 1){
    
    # Generate for Case 1: Constant treatment effect (Under the Null).
    
    if(params$BetaEvent == 0){
      
      # Case 1a: No covariate effect. 
      beta_e <- c(log(1), log(1))
      beta_d <- c(log(1), log(1))
      
    }else{
      
      # Case 1b: Covariate effect.
      beta_e <- c(0, log(2))
      beta_d <- c(0, log(0.5))
      
    }
    
    covariate <- data.frame(arm = c(rep(0, n), rep(1, n)),
                            covar = stats::rnorm(n*2))
    
    data <- MCC::GenData(
      n = n * 2,
      censoring_rate = censoring_rate,
      base_death_rate = params$BaseDeath0,
      base_event_rate = params$BaseEvent0,
      beta_death = beta_d,
      beta_event = beta_e,
      covariates = covariate,
      frailty_variance = params$frailtyVar,
      tau = params$time
    )
    
    return(data)
    
    }else if(params$experiment == 2){
      
    # Generate for Case 2: Constant treatment effect (Under the Alternative).
      
    beta_d <- params$BetaDeath
    beta_e <- params$BetaEvent
    covariate0 <- data.frame(arm = rep(0, n))
    covariate1 <- data.frame(arm = rep(1, n))
    
    # Control arm.
    data0 <- MCC::GenData(
      n = n,
      censoring_rate = censoring_rate,
      base_death_rate = params$BaseDeath0,
      base_event_rate = params$BaseEvent0,
      beta_death = beta_d,
      beta_event = beta_e,
      covariates = covariate0,
      frailty_variance = params$frailtyVar,
      tau = params$time
    )
    
    # Treatment arm. 
    data1 <- MCC::GenData(
      n = n,
      censoring_rate = censoring_rate,
      base_death_rate = params$BaseDeath1,
      base_event_rate = params$BaseEvent1,
      beta_death = beta_d,
      beta_event = beta_e,
      covariates = covariate1,
      frailty_variance = params$frailtyVar,
      tau = params$time
    )
    
    # Combine the data for control arm and treatment arm.
    data1$idx <- data1$idx + n
    data <- rbind(data0, data1)
    
    return(data)
    
    }else if(params$experiment == 3){
    
    # Generate for Case 3: Time-varying treatment effect.
    data <- SimData(n = n,
                    lambda_cens = params$censor,
                    lambda_death = params$BaseDeath0,
                    tau = params$time, 
                    beta = params$TV_effect)
    
    return(data)
    }

}
  
  
# Convenience function to loop across simulation replicates.
t0 <- proc.time()
SimulationLoop <- function(i) {
  
  set.seed(i)
  data <- SimulateData(params)
  
  boot <- try(
    MCC::CompareAUCs(data, tau = params$time)
  )
  
  if (class(boot) != "try-error") {
    
    # Comparison Methods. 
    coxp <- coxPHmodel(data)
    lwyy <- LWYYmodel(data)
    nb <- NegBmodel(data)
    frailty <- FRYmodel(data)
    wr <- WRstat(data)
    
    # AUMCF (Difference).
    aumcf_diff <- data.frame(
      value = boot@CIs$observed[1],
      se = boot@CIs$se[1],
      lower = boot@CIs$lower[1],
      upper = boot@CIs$upper[1],
      p_value = boot@Pvals$p[1],
      type = "aucmf_diff"
    )
    
    # AUMCF (Ratio).
    aumcf_ratio <- data.frame(
      value = boot@CIs$observed[2],
      se = boot@CIs$se[2],
      lower = boot@CIs$lower[2],
      upper = boot@CIs$upper[2],
      p_value = boot@Pvals$p[2],
      type = "aucmf_ratio"
    )
    
    # Return results.
    results <- rbind(aumcf_diff, aumcf_ratio, coxp, lwyy, nb, frailty, wr)
    
    # If needed, run the covariate adjusted analysis. 
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

# Run the simulation and format the output. 
output <- lapply(1:params$reps, SimulationLoop) 
sim <- do.call(rbind, output)

# -----------------------------------------------------------------------------
# Compute the true parameter values.  
# -----------------------------------------------------------------------------

method_names <-unique(sim$type)

if(params$experiment == 1 | (params$experiment == 3 & params$TV_effect == 0) ){
  
  # Null Case
  truth_values <- c(0, rep(1, length(method_names) - 1))
  truth <- setNames(truth_values, method_names)

}else{
  
  # Non-Null Case
  big_data <- SimulateData(params, calc_truth = TRUE)
  
  boot <- try(
    MCC::CompareAUCs(big_data, tau = params$time)
  )
  
  if (class(boot) != "try-error") {
    
    # Comparison Methods. 
    coxp <- coxPHmodel(big_data)
    lwyy <- LWYYmodel(big_data)
    nb <- NegBmodel(big_data)
    frailty <- FRYmodel(big_data)
    # wr <- WRstat(big_data) NOTE issue here
    
    # AUMCF (Difference).
    aumcf_diff <- data.frame(
      value = boot@CIs$observed[1],
      se = boot@CIs$se[1],
      lower = boot@CIs$lower[1],
      upper = boot@CIs$upper[1],
      p_value = boot@Pvals$p[1],
      type = "aucmf_diff"
    )
    
    # AUMCF (Ratio).
    aumcf_ratio <- data.frame(
      value = boot@CIs$observed[2],
      se = boot@CIs$se[2],
      lower = boot@CIs$lower[2],
      upper = boot@CIs$upper[2],
      p_value = boot@Pvals$p[2],
      type = "aucmf_ratio"
    )
    
    # Return results.
    truth_results <- rbind(aumcf_diff, aumcf_ratio, coxp, lwyy, nb, frailty)
    # Note: Filler for WR methods.
    truth_values <- c(truth_results[, "value"], rep(1, 4))
    truth <- setNames(truth_values, method_names)
  }

}

# -----------------------------------------------------------------------------
# Summarize the results. 
# -----------------------------------------------------------------------------
# Add true value to the simulation.
sim_augmented <- sim %>%
  mutate(true_value = truth[type])

# Calculate summary stats.
summary_table <- sim_augmented %>%
  group_by(type) %>%
  summarise(
    bias = mean(value) - first(true_value),
    prob_reject_H0 = mean(p_value < 0.05, na.rm = TRUE),
    cov_p = mean(lower <= first(true_value) & first(true_value) <= upper),
    ase = mean(se),
    ese = sd(value),
    p_value = mean(p_value),
    .groups = "drop"
  )

# Save summary stats + true value.
summary_table_tv <- summary_table %>%
  mutate(true_value = truth[type])
out <- data.frame(summary_table_tv)
print(out)

# Store simulation settings.
out$n <- params$n
out$time <- params$time
out$rep <- params$reps

# Output directory. 
out_stem <- paste0(params$out)


# Save results. 
if(params$experiment != 3 ){
  
sim_file <- paste0(out_stem, "sim",
                   "N", params$n, 
                   "_T", params$time,
                   "_l0", params$BaseEvent0,
                   "_l1", params$BaseEvent1,
                   "_F", params$frailtyVar,
                   "_c", params$censor,
                   ".rds")
}else{
  
  sim_file <- paste0(out_stem, "sim",
                     "N", params$n, 
                     "_T", params$time,
                     "_c", params$censor,
                     "_TV", params$TV_effect,
                     ".rds")
  
}

saveRDS(object = sim_augmented, file = sim_file)

# -----------------------------------------------------------------------------
# Runtime.
# -----------------------------------------------------------------------------
t1 <- proc.time()
cat(t1-t0, "\n")



