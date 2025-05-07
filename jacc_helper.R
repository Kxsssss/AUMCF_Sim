library(dplyr)
library(stringr)
library(purrr)

# Step 1: Generate all expected filenames
n_vals <- c(50, 100, 200, 400)
t_vals <- 1:4
t_vals <- 2
expected_files <- expand.grid(N = n_vals, T = t_vals) %>%
  mutate(
    filename = paste0("N", N, "_T", T, "_l01_l11_adj1.rds"),
    full_path = file.path("~/Documents/AUMCF_Sim/jacc/e3/blog/", filename)
  )

# Step 2: Define processing function
process_file <- function(file_path, filename) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  df <- readRDS(file_path)
  return(df)
}

# Step 3: Apply across all files and combine
summary_table <- map2_dfr(expected_files$full_path, expected_files$filename, process_file)

updated_table <- summary_table %>%
  mutate(
    bias = if_else(type == "wr", bias + 0.5 - 2, bias),
    true_value = if_else(type == "wr", true_value + 1.5, true_value)
  )

# Step 4: Save the table
setwd("~/Documents/AUMCF_Sim/jacc")
out_suffix <- paste0("informative_Cov",".rds")
out_file <- paste0(out_suffix)
saveRDS(object = summary_table, file = out_file)

###########################################sim####################################################
expected_files <- expand.grid(N = n_vals, T = t_vals) %>%
  mutate(
    filename = paste0("simN", N, "_T", T, "_l01_l11.rds"),
    full_path = file.path("~/Documents/AUMCF_Sim/jacc/e3/blog/", filename)
  )

# Step 2: Define processing function
process_file <- function(file_path, filename) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  df <- readRDS(file_path)
  
  # Extract N and T from filename
  n_val <- as.numeric(str_extract(filename, "(?<=N)\\d+"))
  t_val <- as.numeric(str_extract(filename, "(?<=T)\\d+"))
  
  # Add time and sample size
  df %>%
    mutate(
      time = t_val,
      n = n_val
    ) %>%
    dplyr::select(type, time, n, value, se, lower, upper, p_value, true_value)
}

# Step 3: Apply across all files and combine
summary_table <- map2_dfr(expected_files$full_path, expected_files$filename, process_file)

updated_table <- summary_table %>%
  mutate(
    true_value = if_else(type == "wr", true_value + 0.5, true_value)
  )

# Step 4: Save the table
setwd("~/Documents/AUMCF_Sim/jacc")
out_suffix <- paste0("informative_Cov_sim",".rds")
out_file <- paste0(out_suffix)
saveRDS(object = summary_table, file = out_file)

###########################################tv####################################################
# Step 1: Generate all filenames
t_vals <- 2
rep_vals <- 1:5

file_names <- expand.grid(T = t_vals, R = rep_vals) %>%
  mutate(
    filename = paste0("N10000_T", T, "_l01_l11_be1_", R, ".rds")
  ) %>%
  pull(filename)

# Step 2: Read and combine with correct full path
all_data <- map_dfr(file_names, function(f) {
  full_path <- file.path("~/Documents/AUMCF_Sim/jacc/e3/blog/tv/", f)
  
  if (!file.exists(full_path)) {
    warning(paste("Missing:", full_path))
    return(NULL)
  }
  
  df <- readRDS(full_path)
  t_val <- as.numeric(str_extract(f, "(?<=T)\\d+"))
  df$time <- t_val
  df
})

# Step 3: Summarise averages
summary_table <- all_data %>%
  group_by(time) %>%
  summarise(
    mean_tv_d = mean(tv_d, na.rm = TRUE),
    mean_tv_r = mean(tv_r, na.rm = TRUE),
    .groups = "drop"
  )

# Step 4: View or export
print(summary_table)

setwd("~/Documents/AUMCF_Sim/jacc")
out_suffix <- paste0("power_log_tv",".rds")
out_file <- paste0(out_suffix)
saveRDS(object = summary_table, file = out_file)



###############################################################################################
t <- informative_Cov_sim
summary_table <- t %>%
  group_by(type, n, time) %>%
  summarise(
    bias = mean(value) - first(true_value),
    cov_p = mean(lower <= true_value & true_value <= upper, na.rm = TRUE),
    lower = mean(lower, na.rm = TRUE),
    upper = mean(upper, na.rm = TRUE),
    ase = mean(se, na.rm = TRUE),
    ese = sd(value, na.rm = TRUE),
    p_value = mean(p_value, na.rm = TRUE),
    true_value = mean(true_value, na.rm = TRUE),  
    rep = n(),  
    .groups = "drop"
  ) %>%
  dplyr::select(type, bias, lower, upper, cov_p, ase, ese, p_value, true_value, n, time, rep)

# Step 4: Save the table
setwd("~/Documents/AUMCF_Sim/jacc")
out_suffix <- paste0("informative_Cov",".rds")
out_file <- paste0(out_suffix)
saveRDS(object = summary_table, file = out_file)
