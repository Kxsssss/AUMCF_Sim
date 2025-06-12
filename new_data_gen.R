# Purpose: Data generation.
# Updated: 2022-05-15

#' Simulate Recurrent Events Data
#' 
#' Simulate recurrent events data using exponential censoring, gap, and death 
#' times. Status is coded as 0 for censoring, 1 for event, 2 for death.
#' \itemize{
#'   \item The subject-specific death rate is calculated as frailty x base_death_rate x 
#'   exp(beta_death x covariates).
#'   \item The subject-specific event rate is calculated as frailty x base_event_rate x
#'   exp(beta_event x covariates).
#' }
#'
#' @param base_death_rate Baseline arrival rate for the terminal event.
#' @param base_event_rate Baseline arrival rate for recurrent events.
#' @param beta_death Numeric vector of log rate ratios for the death rate.
#' @param beta_event Numeric vector of log rate ratios for the event rate.
#' @param censoring_rate Arrival rate for the censoring time.
#' @param covariates Numeric design matrix.
#' @param frailty_variance Variance of the gamma frailty.
#' @param min_death_rate Minimum subject-specific event rate. Must be positive.
#' @param min_event_rate Minimum subject-specific event rate. Must be positive.
#' @param n Number of subjects. Overwritten by `nrow(covariates)` if covariates are provided.
#' @param tau Truncation time.
#' @return Data.frame, containing:
#' \itemize{
#'  \item The subject identifier `idx` (index).
#'  \item `time` and `status` of the event: 0 for censoring, 1 for an 
#'     event, 2 for death.
#'  \item The `true_death_rate`, `true_event_rate`, and `frailty`, which are subject-specific.
#' }
#' @export
GenData <- function(
    base_death_rate = 0.25,
    base_event_rate = 1.0,
    beta_death = NULL,
    beta_event = NULL,
    censoring_rate = 0.25,
    covariates = NULL,
    frailty_variance = 0.0,
    min_death_rate = 0.05,
    min_event_rate = 0.05,
    n = 100,
    tau = 4,
    time_threshold = 2,
    post_threshold_beta_event = log(0.5)
) {
  
  if (is.null(covariates)) {
    covariates <- data.matrix(rep(1, n))
  } else {
    covariates <- data.matrix(covariates)
    n <- nrow(covariates)
  }
  df <- data.frame(idx = seq_len(n), covariates)
  
  if (is.null(beta_death)) beta_death <- rep(0, ncol(covariates))
  if (is.null(beta_event)) beta_event <- rep(0, ncol(covariates))
  
  df$cens_rate <- censoring_rate
  df$death_rate <- base_death_rate * exp(covariates %*% beta_death)
  df$death_rate <- pmax(df$death_rate, min_death_rate)
  
  if (frailty_variance > 0) {
    theta <- 1 / frailty_variance
    df$frailty <- stats::rgamma(n = n, shape = theta, rate = theta)
  } else {
    df$frailty <- 1
  }
  
  df$death_rate <- df$death_rate * df$frailty
  df$cens_time <- rexp(n, rate = df$cens_rate)
  df$death_time <- rexp(n, rate = df$death_rate)
  
  sim_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    id <- df$idx[i]
    t <- 0
    times <- c()
    statuses <- c()
    
    x_i <- covariates[i, , drop = FALSE]
    
    while (TRUE) {
      # Determine current beta_event based on time
      current_beta_event <- if (t < time_threshold) beta_event else rep(post_threshold_beta_event, length(beta_event))
      
      # Compute current event rate
      lambda <- base_event_rate * exp(x_i %*% current_beta_event) * df$frailty[i]
      lambda <- pmax(lambda, min_event_rate)
      
      gap_time <- rexp(1, rate = lambda)
      t_next <- t + gap_time
      
      # Check stopping conditions
      t_cens <- df$cens_time[i]
      t_death <- df$death_time[i]
      
      if (t_next >= t_cens & t_cens <= t_death & t_cens <= tau) {
        times <- c(times, t_cens)
        statuses <- c(statuses, 0)
        break
      } else if (t_next >= t_death & t_death <= t_cens & t_death <= tau) {
        times <- c(times, t_death)
        statuses <- c(statuses, 2)
        break
      } else if (t_next >= tau) {
        times <- c(times, tau)
        statuses <- c(statuses, 0)
        break
      } else {
        times <- c(times, t_next)
        statuses <- c(statuses, 1)
        t <- t_next
      }
    }
    
    sim_list[[i]] <- data.frame(
      idx = id,
      time = times,
      status = statuses
    )
  }
  
  data <- do.call(rbind, sim_list)
  out <- merge(data, df, by = "idx")
  return(out)
}