#setwd("~/Documents/AUMCF_Sim")
#rm(list = ls())

##########################################################################
##### All functions will return a data frame which contains:
#   value    se   lower    upper   p_value   type
#   (ratio)       (95% ci       )            (name of the model)
##########################################################################
setwd('~/Documents/GitHub/AUMCF_Sim')
# Cox Proportional Hazards Model (First Event Only)
library(dplyr)
library(survival)
cox_prop <- function(data){
  first_event_data <- data %>%
    group_by(idx) %>%
    arrange(time) %>%
    slice(1) %>%
    ungroup()
  
  fit_coxph <- coxph(Surv(time, status != 0) ~ arm, 
                     data = first_event_data,
                     control = coxph.control(timefix = FALSE)
                     )
  
  s_coxp <- summary(fit_coxph)
  coxp_hr <- s_coxp$coefficients[1,2]
  coxp_se <- s_coxp$coefficients[1,3]
  coxp_ci_l <- s_coxp$conf.int[1,3]
  coxp_ci_u <- s_coxp$conf.int[1,4]
  coxp_p <- s_coxp$coefficients[1,5]
  #coxp_z <- s_coxp$coefficients[1,4]
  
  result <- data.frame(value = coxp_hr,
                       se = coxp_se,
                       lower = coxp_ci_l,
                       upper = coxp_ci_u,
                       p_value = coxp_p,
                       #z_value = coxp_z,
                       type  = "cox_p")
  return(result)
}


# LWYY
#install.packages("reReg")
library(reReg)
lwyy <- function(data){
  fit_lwyy <- reReg(
    Recur(time, idx, status == 1) ~ arm,
    data = data,
    model = "cox.LWYY"
  )
  
  s_lwyy <- summary(fit_lwyy)
  lwyy_coef <- s_lwyy$coefficients.rec[1,1]
  lwyy_se <- s_lwyy$coefficients.rec[1,2]
  lwyy_hr <- exp(lwyy_coef)
  lwyy_ci_l <- exp(lwyy_coef - 1.96 * lwyy_se)
  lwyy_ci_u <- exp(lwyy_coef + 1.96 * lwyy_se)
  lwyy_p <- s_lwyy$coefficients[1,4]
  #lwyy_z <- s_lwyy$coefficients[1,3]
  
  result <- data.frame(value = lwyy_hr,
                       se = lwyy_se,
                       lower = lwyy_ci_l,
                       upper = lwyy_ci_u,
                       p_value = lwyy_p,
                       #z_value = lwyy_z,
                       type  = "lwyy")
  return(result)
}


# Negative Binomial Rate Ratio
library(dplyr)
library(MASS)

nb <- function(data){
  nb_data <- data %>%
    group_by(idx, arm) %>%
    summarise(
      count = sum(status == 1),       # number of recurrent events
      followup = max(time),           # total follow-up time
      .groups = "drop"
    )
  
  fit_nb <- glm.nb(count ~ arm + offset(log(followup)), data = nb_data)
  s_nb <- summary(fit_nb)
  nb_coef <- s_nb$coefficients[2,1]
  nb_se <- s_nb$coefficients[2,2]
  nb_rr <- exp(nb_coef)
  nb_ci_l <- exp(nb_coef - 1.96 * nb_se)
  nb_ci_u <- exp(nb_coef + 1.96 * nb_se)
  nb_p <- s_nb$coefficients[2,4]
  #nb_z <- s_nb$coefficients[2,3]
  
  result <- data.frame(value = nb_rr,
                       se = nb_se,
                       lower = nb_ci_l,
                       upper = nb_ci_u,
                       p_value = nb_p,
                       #z_value = nb_z,
                       type  = "nb")
  return(result)
}


# Joint Frailty Model
#install.packages("frailtypack")
library(frailtypack)

frailty_ <- function(data){
  df_frailty <- data %>%
    filter(status %in% c(0, 1)) %>%
    arrange(idx, time) %>%
    group_by(idx) %>%
    mutate(start = lag(time, default = 0)) %>%
    ungroup()
  
  fit_frailty <- frailtyPenal(
    formula = Surv(start, time, status) ~ cluster(idx) + arm,
    data = df_frailty,
    recurrentAG = TRUE, 
    n.knots = 10,
    kappa = 1e-2
  )
  
  frailty_coef <- fit_frailty$coef["arm"]
  frailty_se <- sqrt(fit_frailty$varH)
  frailty_hr <- exp(frailty_coef)
  frailty_ci_l <- exp(frailty_coef - 1.96 * frailty_se)
  frailty_ci_u <- exp(frailty_coef + 1.96 * frailty_se)
  frailty_p <- fit_frailty$beta_p.value["arm"]
  #frailty_z <- frailty_coef/frailty_se
  
  result <- data.frame(value = frailty_hr,
                       se = frailty_se,
                       lower = frailty_ci_l,
                       upper = frailty_ci_u,
                       p_value = frailty_p,
                       #z_value = frailty_z,
                       type  = "frailty")
  return(result)
}



# win ratio
# Load packages
library(dplyr)
library(survival)
library(WR)

wr <- function(data){
  # Step 1: Collapse long format to one row per subject
  df_wr <- data %>%
    group_by(idx, arm) %>%
    summarise(
      time = max(time),                        # Time to last observed event (death/censoring)
      event = as.integer(any(status == 2)),    # 1 = death, 0 = censored or recurrent only
      recurrent_count = sum(status == 1),      # Number of recurrent events
      .groups = "drop"
    )
  
  # Step 2: Generate all possible treatment vs. control pairs
  treated <- df_wr %>% filter(arm == 1)
  control <- df_wr %>% filter(arm == 0)
  
  # Initialize counters
  wins <- 0
  losses <- 0
  ties <- 0
  
  # Step 3: Pairwise comparisons
  for (i in 1:nrow(treated)) {
    for (j in 1:nrow(control)) {
      
      ti <- treated[i, ]
      cj <- control[j, ]
      
      # Step 1: Compare death outcome
      if (ti$event == 1 && cj$event == 1) {
        # both died → compare death times
        if (ti$time > cj$time) {
          wins <- wins + 1
        } else if (ti$time < cj$time) {
          losses <- losses + 1
        } else {
          # tie in death time → go to recurrence
          if (ti$recurrent_count < cj$recurrent_count) {
            wins <- wins + 1
          } else if (ti$recurrent_count > cj$recurrent_count) {
            losses <- losses + 1
          } else {
            ties <- ties + 1
          }
        }
      } else if (ti$event == 1 && cj$event == 0) {
        # treated died, control didn't
        losses <- losses + 1
      } else if (ti$event == 0 && cj$event == 1) {
        # treated alive, control died
        wins <- wins + 1
      } else {
        # neither died → compare recurrence
        if (ti$recurrent_count < cj$recurrent_count) {
          wins <- wins + 1
        } else if (ti$recurrent_count > cj$recurrent_count) {
          losses <- losses + 1
        } else {
          ties <- ties + 1
        }
      }
    }
  }
  
  # Step 4: Calculate Win Ratio
  win_ratio <- wins / losses
  
  # log scale
  log_wr <- log(win_ratio)
  se_log_wr <- sqrt(1/wins + 1/losses)
  
  # CI
  wr_ci_l<- exp(log_wr - 1.96 * se_log_wr)
  wr_ci_u <- exp(log_wr + 1.96 * se_log_wr)
  
  # z and p (null hypothesis WR = 1)
  wr_z <- log_wr / se_log_wr
  wr_p <- 2 * (1 - pnorm(abs(wr_z)))
  
  # try it with the packages.
  #wr_result <- WRrec(ID, time, status, trt, strata = NULL, naive = FALSE)
  
  
  result <- data.frame(value = win_ratio,
                       se = win_ratio * se_log_wr,
                       lower = wr_ci_l,
                       upper = wr_ci_u,
                       p_value = wr_p,
                       #z_value = wr_z,
                       type  = "wr")
  return(result)
}


wr_rec <- function(data){

  wr_rec_all <- WRrec(data[, "idx"],
          data[, "time"],
          data[, "status"],
          data[, "arm"],
          strata = NULL, 
          naive = TRUE)
  
  result_LWR <- data.frame(value = exp(wr_rec_all$log.WR),
                           se = wr_rec_all$se * exp(wr_rec_all$log.WR),
                           lower = exp(wr_rec_all$log.WR - 1.96 * wr_rec_all$se),
                           upper = exp(wr_rec_all$log.WR + 1.96 * wr_rec_all$se),
                           p_value = wr_rec_all$pval,
                           type  = "wr_LWR")
  
  wr_z_FWR <- wr_rec_all$log.WR.FI/wr_rec_all$se.FI
  p_val_FWR <- 2 * (1 - pnorm(abs(wr_z_FWR)))
  result_FWR <- data.frame(value = exp(wr_rec_all$log.WR.FI),
                           se = wr_rec_all$se.FI * exp(wr_rec_all$log.WR.FI),
                           lower = exp(wr_rec_all$log.WR.FI - 1.96 * wr_rec_all$se.FI),
                           upper = exp(wr_rec_all$log.WR.FI + 1.96 * wr_rec_all$se.FI),
                           p_value =  p_val_FWR,
                           type  = "wr_FWR")
  
  wr_z_NWR <- wr_rec_all$log.WR.naive/wr_rec_all$se.naive
  p_val_NWR <- 2 * (1 - pnorm(abs(wr_z_NWR)))
  result_NWR <- data.frame(value = exp(wr_rec_all$log.WR.naive),
                           se = wr_rec_all$se.naive * exp(wr_rec_all$log.WR.naive),
                           lower = exp(wr_rec_all$log.WR.naive - 1.96 * wr_rec_all$se.naive),
                           upper = exp(wr_rec_all$log.WR.naive + 1.96 * wr_rec_all$se.naive),
                           p_value =  p_val_NWR,
                           type  = "wr_NWR")
  
  result <- rbind(result_LWR, result_FWR, result_NWR)
  
  return(result)

}