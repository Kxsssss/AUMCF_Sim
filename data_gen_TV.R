# Purpose: Simulate data with a time-varying treatment effect.
# Updated: 2025-06-12

#' Simulate Observation Time and Status
#'
#' @param lambda_cens Rate for the time-to censoring.
#' @param lambda_death Rate for the time-to death.
#' @param tau Truncation time.
#' @return Data.frame.
.ObsTimeStatus <- function(
    lambda_cens,
    lambda_death,
    tau
) {
  
  # Time to censoring.
  if (lambda_cens > 0) {
    time_cens <- stats::rexp(n = 1, rate = lambda_cens)
    time_cens <- min(time_cens, tau)
  } else {
    time_cens <- tau
  }
  
  # Time to death.
  if (lambda_death > 0) {
    time_death <- stats::rexp(n = 1, rate = lambda_death)
  } else {
    time_death <- Inf
  }
  
  # Final status.
  obs_time <- ifelse(time_death <= time_cens, time_death, time_cens)
  status <- ifelse(obs_time == time_death, 2, 0)
  out <- data.frame(
    time = obs_time,
    status = status
  )
  return(out)
}


#' Find Max Lambda
#'
#' @param lambda_fn A *function* of time that returns the rate for the 
#'   recurrent event process.
#' @param tau Truncation time.
#' @param n_grid Number of grid points.
#' @return Scalar maximum event rate.
.LambdaMax <- function(lambda_fn, tau, n_grid = 100) {
  grid <- seq(from = 0, to = tau, length.out = n_grid)
  rates <- lambda_fn(grid)
  max_lambda <- max(rates)
  return(max_lambda)
}


#' Simulate for a Single Subject
#' 
#' Simulates data from a non-homogeneous Poisson process by thinning
#' a homogeneous Poisson process with rate = lambda_max.
#'
#' @param lambda_fn A *function* of time that returns the rate for the 
#'   recurrent event process.
#' @param lambda_max Maximum value of event rate.
#' @param lambda_cens Rate for the time-to censoring.
#' @param lambda_death Rate for the time-to death.
#' @param tau The truncation time.
#' @return Data.frame.
SimSubj <- function(
    lambda_fn,
    lambda_max,
    lambda_cens = 0.25,
    lambda_death = 0.25,
    tau = 4
) {
  
  # Time and status.
  df <- .ObsTimeStatus(
    lambda_cens = lambda_cens, 
    lambda_death = lambda_death,
    tau = tau
  )
  
  # Poisson process.
  event_times <- numeric(0)
  t <- 0
  while (t < df$time) {
    
    # Simulate gap time.
    gap <- stats::rexp(n = 1, rate = lambda_max)
    
    # Increment homogeneous process.
    t <- t + gap
    if (t > df$time) {break}
    
    # Decide whether to accept proposal.
    pi <- min(lambda_fn(t) / lambda_max, 1)
    accept <- stats::rbinom(n = 1, size = 1, prob = pi)
    if (accept == 1) {
      event_times <- c(event_times, t)
    }
    
  }
  
  # Output.
  if (length(event_times) > 0) {
    out <- data.frame(
      time = event_times,
      status = 1
    )
    out <- rbind(out, df)
  } else {
    out <- df
  }
  return(out)
}


#' Simulate arm
#' 
#' @param lambda_fn A *function* of time that returns the rate for the 
#'   recurrent event process.
#' @param n Sample size.
#' @param lambda_cens Rate for the time-to censoring.
#' @param lambda_death Rate for the time-to death.
#' @param lambda_max Maximum value of event rate.
#' @param tau The truncation time.
#' @return Data.frame.
SimArm <- function(
    lambda_fn,
    n,
    lambda_max = NULL,
    lambda_cens = 0.25,
    lambda_death = 0.25,
    tau = 4
) {
  
  lambda_fn <- Vectorize(lambda_fn)
  
  # Find lambda max.
  if (is.null(lambda_max)) {
    lambda_max <- .LambdaMax(lambda_fn, tau)
  }
  
  # Generate data.
  data <- lapply(seq_len(n), function(idx) {
    
    df <- SimSubj(
      lambda_fn = lambda_fn,
      lambda_max = lambda_max,
      lambda_cens = lambda_cens,
      lambda_death = lambda_death,
      tau = tau
    )
    df$idx <- idx
    return(df)
    
  })
  data <- do.call(rbind, data)
  return(data)
}


#' Simulate data
#' 
#' Set the rate functions for the treatment and control arms within this function.
#' 
#' @param n Sample size.
#' @param lambda_cens Rate for the time-to censoring.
#' @param lambda_death Rate for the time-to death.
#' @param tau The truncation time.
#' @return Data.frame.
SimData <- function(
    n,
    lambda_cens = 0.25,
    lambda_death = 0.25,
    tau = 4.0
) {
  
  # Base event rate.
  lambda_base <- 1.0
  
  # Treatment effect.
  beta <- log(0.5)
  
  # Break point, after which the rate for the treatment arm switches
  # from lambda_base to lambda_base * exp(beta).
  bp <- 1.0
  
  # Rate for treatment arm.
  rate_1 <- function(time) {
    if (time <= bp) {
      out <- lambda_base
    } else {
      out <- lambda_base * exp(beta)
    }
    return(out)
  }
  
  # Rate for control arm.
  rate_0 <- function(time) {
    out <- lambda_base
    return(out)
  }
  
  # Simulate data for treatment arm.
  data_1 <- SimArm(
    lambda_fn = rate_1,
    n = n,
    lambda_cens = lambda_cens,
    lambda_death = lambda_death,
    tau = tau
  )
  data_1$arm <- 1
  
  # Simulate data for control arm.
  data_0 <- SimArm(
    lambda_fn = rate_0,
    n = n,
    lambda_cens = lambda_cens,
    lambda_death = lambda_death,
    tau = tau
  )
  data_0$arm <- 0
  data_0$idx <- data_0$idx + n
  
  # Final data set.
  data <- rbind(data_1, data_0)
  data <- data[, c("idx", "arm", "time", "status")]
  return(data)
  
}

# Simulate data.
data <- SimData(n = 100)


