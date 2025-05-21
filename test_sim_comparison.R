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
    # If BetaEvent is set to be 0, then run the null case
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
    wr_rec <- wr_rec(data)
    wr <- wr(data)
    
    aucmf_diff <- data.frame(
      value = boot@CIs$observed[1],
      se = boot@CIs$se[1],
      lower = boot@CIs$lower[1],
      upper = boot@CIs$upper[1],
      p_value = boot@Pvals$p[1],
      type = "aucmf_diff"
    )
    
    results <- rbind(aucmf, coxp, lwyy, nb, frailty, wr, wr_rec, aucmf_diff)
    #print(results)
    
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

output <- lapply(1 : params$reps, Loop) 
sim <- do.call(rbind, output)
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



# Sanity check
# print(dim(sim_augmented)) # should be 2000 * 9 = 18000
# print(sim_augmented %>% 
#   group_by(type) %>%
#   summarise(mean_pvalue = mean(p_value < 0.05)))
# 
# print(sim_augmented %>% 
#         group_by(type) %>%
#         summarise(mean_val = mean(value)))
# 
# print(head(sim_augmented, 9))