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
coxPH_model <- function(data){
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
  
  fit_lwyy <- tryCatch(reReg(
    Recur(time, idx, status == 1) ~ arm,
    data = data,
    model = "cox.LWYY"
  ), error=function(e) NULL)
  
  
  if(!is.null(fit_lwyy)){
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
  }else{
    
    result <- data.frame(value = NA,
                         se = NA,
                         lower = NA,
                         upper = NA,
                         p_value = NA,
                         #z_value = lwyy_z,
                         type  = "lwyy")
    
  }     
  
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
    n.knots = 7,
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
# win ratio
# Load packages
library(dplyr)
library(survival)
library(WinRatio)

wr <- function(data){
  
  # Creation of dataset 'df' with 2 outcomes:
  #   Outcome 1: death (survival event)
  #   Outcome 2: recurrence (repeated survival event)
  data1 <- data  %>%  
    group_by(idx) %>% 
    mutate(death = if_else(max(status) == 2, 1, 0),
           t2death = max(time)) %>% 
    ungroup() %>%
    dplyr::select(idx, arm, death, t2death) %>% 
    unique()
  
  data2 <- data %>% 
    filter(status == 1) %>% 
    dplyr::select(idx, t2recurr = time) %>%
    mutate(recurr = 1) %>% 
    arrange(idx, t2recurr) %>% 
    group_by(idx) %>% 
    mutate(nrecurr = row_number()) %>% 
    ungroup() %>% 
    full_join((data1 %>% dplyr::select(idx)), by = "idx") %>% 
    complete(idx, nrecurr) %>%
    filter(!is.na(nrecurr)) %>% 
    full_join((data1 %>% dplyr::select(idx, t2death)), by = "idx") %>% 
    mutate(
      recurr = replace(recurr, is.na(recurr), 0),
      t2recurr = if_else(is.na(t2recurr), t2death, t2recurr)
    ) %>% 
    dplyr::select(idx, nrecurr, recurr, t2recurr)
  
  
  data3 <- data2 %>% 
    pivot_wider(id_cols = "idx", names_from =  nrecurr, 
                values_from = recurr, names_prefix  = "recurr") 
  
  data4 <- data2 %>% 
    pivot_wider(id_cols = "idx", names_from =  nrecurr, 
                values_from = t2recurr, names_prefix  = "t2recurr") 
  
  data5 <- full_join(data3, data4, by = "idx")
  
  df <- full_join(data1, data5, by = "idx")
  
  max_k <- max(data2$nrecurr)
  recurr_vars <- paste0("recurr", 1:max_k)
  t2recurr_vars <- paste0("t2recurr", 1:max_k)
  
  # Calculate the win ratio
  wr <- winratio(id = "idx", trt = "arm", active = 1, 
                 outcomes = list(outc1 = c("death", "s", "t2death"),
                                 outc2 =  list(recurr_vars, "r", t2recurr_vars)),
                 fu = "t2death", data = df)
  
  wr_se <- (log(wr$wr.upper) - log(wr$wr.lower)) / (2 * 1.96)
  
  result <- data.frame(value = wr$wr,
                       se = wr_se,
                       lower = wr$wr.lower,
                       upper = wr$wr.upper,
                       p_value = wr$p.value,
                       type  = "wr")
  return(result)
}

library(WR)
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
  
  # Calculate the standard win ratio.
  o <- order(data$idx, data$time)
  dat_o <- data[o,]
  dat_oc <-dat_o[!duplicated(dat_o[c("idx","status")]),]
  
  wr_rec_std <- WRrec(dat_oc[, "idx"],
                      dat_oc[, "time"],
                      dat_oc[, "status"],
                      dat_oc[, "arm"],
                      strata = NULL, 
                      naive = FALSE)
  result_STD <- data.frame(value = exp(wr_rec_std$log.WR),
                           se = wr_rec_std$se * exp(wr_rec_std$log.WR),
                           lower = exp(wr_rec_std$log.WR - 1.96 * wr_rec_std$se),
                           upper = exp(wr_rec_std$log.WR + 1.96 * wr_rec_std$se),
                           p_value = wr_rec_std$pval,
                           type  = "wr_STD")
  
  result <- rbind(result_LWR, result_FWR, result_NWR, result_STD)
  
  return(result)
  
}


# Ghosh + Lin - gives error.
#install.packages("reReg")
library(reReg)
gl <- function(data){
  data$event <- (data$status == 1) * 1
  data$status_gl <- (data$status == 2) * 1
  fm <- Recur(time, idx, event, status_gl) ~ arm
  fit_gl <- reReg(
    fm,
    data = data,
    model = "cox.GL"
  )
  
  s_gl <- summary(fit_gl)
  gl_coef <- s_gl$coefficients.rec[1,1]
  gl_se <- s_gl$coefficients.rec[1,2]
  gl_hr <- exp(gl_coef)
  gl_ci_l <- exp(gl_coef - 1.96 * gl_se)
  gl_ci_u <- exp(gl_coef + 1.96 * gl_se)
  gl_p <- s_gl$coefficients[1,4]
  #gl_z <- s_gl$coefficients[1,3]
  
  result <- data.frame(value = gl_hr,
                       se = gl_se,
                       lower = gl_ci_l,
                       upper = gl_ci_u,
                       p_value = gl_p,
                       #z_value = gl_z,
                       type  = "gl")
  return(result)
}





library(reReg)
hh <- function(data){
  data$event <- (data$status == 1) * 1
  data$status_gl <- (data$status == 2) * 1
  fm <- Recur(time, idx, event, status_gl) ~ arm
  fit_gl <- reReg(
    fm,
    data = data,
    model = "cox.HH"
  )
  
  s_gl <- summary(fit_gl)
  gl_coef <- s_gl$coefficients.rec[1,1]
  gl_se <- s_gl$coefficients.rec[1,2]
  gl_hr <- exp(gl_coef)
  gl_ci_l <- exp(gl_coef - 1.96 * gl_se)
  gl_ci_u <- exp(gl_coef + 1.96 * gl_se)
  gl_p <- s_gl$coefficients[1,4]
  #gl_z <- s_gl$coefficients[1,3]
  
  result <- data.frame(value = gl_hr,
                       se = gl_se,
                       lower = gl_ci_l,
                       upper = gl_ci_u,
                       p_value = gl_p,
                       #z_value = gl_z,
                       type  = "gl")
  return(result)
}


