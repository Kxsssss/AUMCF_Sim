# Purpose: Plot recurrent events paper simulation results.
# Updated: 2024-12-02

library(dplyr)
library(ggplot2)

# -------------------------------------------------------------------------

# Theme options.
gg_theme <- theme_bw() + 
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank(),
    axis.line = element_line(color = "black")
  )

# -------------------------------------------------------------------------
# Bias plots.
# -------------------------------------------------------------------------

data <- readRDS(file = "data/E1_10k_sim.rds")
data$bias <- data$observed - data$true_value

data$n_factor <- factor(
  x = data$n,
  levels = c(50, 100, 200, 400),
  ordered = TRUE
)

data$tau_factor <- factor(
  x = data$tau,
  levels = c(1, 2, 3, 4),
  ordered = TRUE
)

# Bias plot.
q <- ggplot(data = data) + 
  gg_theme + 
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray"
  ) + 
  geom_boxplot(
    aes(x = n_factor, y = bias, fill = n_factor),
    show.legend = FALSE
  ) +
  facet_grid(
    . ~ tau_factor,
    labeller = label_bquote(cols = tau == .(tau_factor))
  ) +
  scale_x_discrete(
    name = "Sample Size"
  ) + 
  scale_y_continuous(
    name = "Bias"
  ) + 
  scale_fill_brewer(
    name = "N",
    palette = "Blues"
  )
show(q)
q_box <- q

if (FALSE) {
  ggsave(
    plot = q,
    file = "figures/null_bias_plots.pdf",
    width = 9.0,
    height = 4.5
  )
}


# -------------------------------------------------------------------------
# QQ plots.
# -------------------------------------------------------------------------

#' Prepare QQ plotting frame
#' @param pvals P-values.
#' @param p_max Maximum p-value to include in plot.
#' @return Data.frame.
PrepQQDF <- function(pvals, p_max = 0.05) {
  
  n_orig <- length(pvals)
  pvals <- pvals[pvals <= p_max]
  n <- length(pvals)
  obs <- -log10(sort(pvals))
  exp <- -log10(seq_len(n) / (n_orig + 1))
  df = data.frame(obs = obs, exp = exp)
  
  return(df)
}

#' Prepare Confidence Bands frame
#' @param n_sim Total number of simulations at each setting.
#' @param p_max Maximum p-value to include in plot.
#' @return Data.frame.
PrepBandsDF <- function(n_sim = 1e4, p_max = 0.05) {
  df <- data.frame(coord = seq_len(n_sim))
  df$exp <- -log10(seq_len(n_sim) / (n_sim + 1))
  df$lower <- -log10(
    stats::qbeta(0.975, df$coord, n_sim - df$coord + 1)
  )
  df$upper <- -log10(
    stats::qbeta(0.025, df$coord, n_sim - df$coord + 1)
  )
  df <- df %>%
    dplyr::filter(exp >= -log10(0.05))
  return(df)
}


# Reducing to a single time point because there are not enough
# dimensions to present multiple.
df_plot <- data %>%
  dplyr::group_by(n_factor, tau_factor) %>%
  dplyr::reframe(
    PrepQQDF(pvals = p_value)
  )
df_plot

# Confidence bands frame.
df_bands <- PrepBandsDF()

q <- ggplot(data = df_plot) +
  gg_theme + 
  facet_grid(
    . ~ tau_factor,
    labeller = label_bquote(cols = tau == .(tau_factor))
  ) +
  geom_abline(
    intercept = 0, 
    slope = 1, 
    lty = "dashed",
  ) +
  geom_ribbon(
    data = df_bands,
    aes(x = exp, ymin = lower, ymax = upper),
    alpha = 0.2,
    fill = "gray"
  ) + 
  geom_point(
    aes(y = obs, x = exp, color = n_factor)
  ) + 
  scale_x_continuous(
    name = expression(Expected~-log[10](p))
  ) + 
  scale_y_continuous(
    name = expression(Observed~-log[10](p))
  ) + 
  scale_color_brewer(
    name = "N",
    palette = "Blues"
  )
show(q)
q_qq <- q

if (FALSE) {
  ggsave(
    plot = q,
    file = "figures/null_qq_plots.pdf",
    width = 9.0,
    height = 4.5
  )
}

if (FALSE) {
  q <- cowplot::plot_grid(
    labels = c("A", "B"),
    plotlist = list(q_box, q_qq),
    ncol = 1
  )
  ggsave(
    plot = q,
    file = "figures/t1e_overall.pdf",
    width = 9.0,
    height = 9.0
  )
}


# -------------------------------------------------------------------------
# Coverage: type I error setting.
# -------------------------------------------------------------------------

z <- stats::qnorm(0.975)
df <- data %>%
  dplyr::group_by(n_factor, tau_factor) %>%
  dplyr::summarise(
    p = mean(cover_ind),
    p_se = sqrt(p * (1 - p) / dplyr::n()),
    p_l = p - z * p_se,
    p_u = p + z * p_se,
    ase = sqrt(mean(asymptotic_se^2)),
    ese = sd(observed)
  ) %>%
  dplyr::mutate(
    p = sprintf("%.1f", 100 * p),
    p_l = sprintf("%.1f", 100 * p_l),
    p_u = sprintf("%.1f", 100 * p_u),
    p_ci = glue::glue("{p} ({p_l}-{p_u})"),
    ase = sprintf("%.3f", ase),
    ese = sprintf("%.3f", ese)
  )

df <- df %>% 
  dplyr::select(n_factor, tau_factor, p_ci, ase, ese)
show(df)

t1e_tab_df <- df
tab <- xtable::xtable(df)
xtable::print.xtable(tab, include.rownames = FALSE)


# -------------------------------------------------------------------------
# Experiment 2: Presence of a true difference
# -------------------------------------------------------------------------

# data <- readRDS(file = "data/E2_sim.rds")
data <- readRDS(file = "data/update/E2_sim.rds")
data$bias <- data$observed - data$true_value

data$n_factor <- factor(
  x = data$n,
  levels = c(50, 100, 200, 400),
  ordered = TRUE
)

data$tau_factor <- factor(
  x = data$tau,
  levels = c(1, 2, 3, 4),
  ordered = TRUE
)

# Bias.
q <- ggplot(data = data) + 
  gg_theme + 
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray"
  ) + 
  geom_boxplot(
    aes(x = n_factor, y = bias, fill = n_factor),
    show.legend = FALSE
  ) +
  facet_grid(
    . ~ tau_factor,
    labeller = label_bquote(cols = tau == .(tau_factor))
  ) +
  scale_x_discrete(
    name = "Sample Size"
  ) + 
  scale_y_continuous(
    name = "Bias"
  ) + 
  scale_fill_brewer(
    name = "N",
    palette = "Blues"
  )
show(q)
q_box <- q


# Average chi2 statistic.
z <- stats::qnorm(0.975)
df <- data %>%
  dplyr::group_by(n_factor, tau_factor) %>%
  dplyr::mutate(x2 = stats::qchisq(p_value, df = 1, lower.tail = FALSE)) %>%
  dplyr::summarise(
    
    mu = mean(x2),
    se = sqrt(var(x2) / dplyr::n()),
    lower = mu - z * se,
    upper = mu + z * se
  )

q <- ggplot(data = df) + 
  gg_theme + 
  geom_col(
    aes(x = n_factor, y = mu, fill = n_factor)
  ) +
  geom_errorbar(
    aes(x = n_factor, ymin = lower, ymax = upper, group = tau_factor),
    width = 0.5
  ) + facet_grid(
    . ~ tau_factor,
    labeller = label_bquote(cols = tau == .(tau_factor))
  ) +
  scale_x_discrete(
    name = "Sample Size"
  ) + 
  scale_y_continuous(
    name = expression(E(chi^2))
  ) + 
  scale_fill_brewer(
    name = "N",
    palette = "Blues"
  )
show(q)
q_pwr <- q

if (FALSE) {
  q <- cowplot::plot_grid(
    labels = c("A", "B"),
    plotlist = list(q_box, q_pwr),
    ncol = 1
  )
  ggsave(
    plot = q,
    file = "figures/power_overall.pdf",
    width = 9.0,
    height = 9.0
  )
}

# Table.
tab <- df %>%
  dplyr::mutate(
    mu = sprintf("%.1f", mu),
    mu_l = sprintf("%.1f", lower),
    mu_u = sprintf("%.1f", upper),
    mu_ci = glue::glue("{mu} ({mu_l}-{mu_u})"),
  ) %>%
  dplyr::select(n_factor, tau_factor, mu_ci)


# -------------------------------------------------------------------------
# Coverage in the power setting. 
# -------------------------------------------------------------------------

z <- stats::qnorm(0.975)
df <- data %>%
  dplyr::group_by(n_factor, tau_factor) %>%
  dplyr::summarise(
    p = mean(cover_ind),
    p_se = sqrt(p * (1 - p) / dplyr::n()),
    p_l = p - z * p_se,
    p_u = p + z * p_se,
    ase = sqrt(mean(asymptotic_se^2)),
    ese = sd(observed)
  ) %>%
  dplyr::mutate(
    p = sprintf("%.1f", 100 * p),
    p_l = sprintf("%.1f", 100 * p_l),
    p_u = sprintf("%.1f", 100 * p_u),
    p_ci = glue::glue("{p} ({p_l}-{p_u})"),
    ase = sprintf("%.3f", ase),
    ese = sprintf("%.3f", ese)
  )

df <- df %>% 
  dplyr::select(n_factor, tau_factor, p_ci, ase, ese)
show(df)

power_tab_df <- df
tab <- t1e_tab_df %>% dplyr::inner_join(power_tab_df, by = c("n_factor", "tau_factor"))
tab <- xtable::xtable(tab)
xtable::print.xtable(tab, include.rownames = FALSE)

# -------------------------------------------------------------------------
# Experiment 3: Augmentation.
# -------------------------------------------------------------------------

df0 <- readRDS(file = "data/E3_b0b0_sim.rds")
df0$covar_effect <- 0
df1 <- readRDS(file = "data/update/E3_log0.5log2_sim.rds")
df1$covar_effect <- 1
data <- rbind(df0, df1)

data$n_factor <- factor(
  x = data$n,
  levels = c(50, 100, 200, 400),
  ordered = TRUE
)

df <- data %>%
  dplyr::group_by(n_factor, covar_effect, adjusted) %>%
  dplyr::summarise(
    ese = sd(observed),
    ase = sqrt(mean(asymptotic_se^2))
  )

df$adj_factor <- factor(
  x = df$adjusted,
  levels = c(0, 1),
  labels = c("Baseline", "Augmented")
)

df$covar_factor <- factor(
  x = df$covar_effect,
  levels = c(0, 1),
  labels = c("Uninformative Covariate", "Informative Covariate")
)

# df <- df %>%
#   tidyr::pivot_longer(
#     ese:ase,
#     names_to = "method",
#     values_to = "se"
#   )
# df$method <- factor(
#   x = df$method,
#   levels = c("ase", "ese"),
#   labels = c("Asymptotic", "Empirical")
# )

q <- ggplot(data = df) + 
  gg_theme + 
  geom_col(
    aes(x = n_factor, y = ese, fill = adj_factor),
    position = position_dodge()
  ) +
  facet_grid(
    . ~ covar_factor
  ) +
  scale_x_discrete(
    name = "Sample Size"
  ) + 
  scale_y_continuous(
    name = "Empirical SE",
    limits = c(0, 0.55),
    breaks = seq(from = 0.0, to = 0.5, by = 0.1)
  ) + 
  ggsci::scale_fill_nejm(
    name = "Estimator"
  )
show(q)
q0 <- q

if (FALSE) {
  ggsave(
    plot = q,
    file = "figures/augmentation.pdf",
    width = 9.0,
    height = 4.5
  )
}

tab <- df %>%
  dplyr::ungroup() %>%
  dplyr::select(n_factor, adj_factor, covar_factor, ese, ase)
tab <- xtable::xtable(tab, digits = 3)
xtable::print.xtable(tab, include.rownames = FALSE)

#' Jackknife variance
#' @param df Data.frame for a single configuration.
Jackknife <- function(df) {
  y0 <- df %>% 
    dplyr::filter(adjusted == 0) %>%
    dplyr::pull(observed)
  y1 <- df %>%
    dplyr::filter(adjusted == 1) %>%
    dplyr::pull(observed)
  n <- length(df0)
  jk <- sapply(seq_len(n), function(i) {
    v0 <- var(y0[-i])
    v1 <- var(y1[-i])
    rel_eff <- v0 / v1
    return(rel_eff)
  })
  jk_var <- (n - 1)^2 / n * var(jk)
  out <- data.frame(
    v0 = var(y0),
    v1 = var(y1),
    rel_eff = var(y0) / var(y1),
    se_rel_eff = sqrt(jk_var)
  )
  return(out)
}

#' Relative efficiency.
#' @param data Data.frame.
CalcRE <- function(data) {
  config <- data %>%
    dplyr::group_by(n_factor, covar_effect) %>%
    dplyr::summarise(n = dplyr::n())
  n_config <- nrow(config)
  
  out <- lapply(seq_len(n_config), function(i) {
    which_n <- config$n_factor[i]
    which_covar <- config$covar_effect[i]
    df <- data %>% dplyr::filter(n_factor == which_n, covar_effect == which_covar)
    out <- Jackknife(df)
    out$n_factor <- which_n
    out$covar_effect <- which_covar
    return(out)
  })
  out <- do.call(rbind, out)
  
  z <- stats::qnorm(0.975)
  out <- out %>%
    dplyr::mutate(
      lower = rel_eff - z * se_rel_eff,
      upper = rel_eff + z * se_rel_eff
    )
  
  return(out)
}

df_rel <- CalcRE(data)
df_rel$n_factor <- factor(
  x = df_rel$n,
  levels = c(50, 100, 200, 400),
  ordered = TRUE
)
df_rel$covar_effect <- factor(
  x = df_rel$covar_effect,
  levels = c(0, 1),
  labels = c("Uninformative Covariate", "Informative Covariate")
)

q <- ggplot(data = df_rel) + 
  gg_theme + 
  geom_col(
    aes(x = n_factor, y = rel_eff),
    fill = "#0072B5FF",
    alpha = 0.8
  ) +
  geom_errorbar(
    aes(x = n_factor, ymin = lower, ymax = upper),
    width = 0.5
  ) + 
  geom_hline(
    yintercept = 1.0,
    linetype = "dashed",
    color = "gray",
    linewidth = 1.3
  ) + 
  facet_grid(
    . ~ covar_effect
  ) +
  scale_x_discrete(
    name = "Sample Size"
  ) + 
  scale_y_continuous(
    name = "Relative Efficiency"
  ) + 
  ggsci::scale_fill_nejm(
    name = "Estimator"
  )
show(q)
q1 <- q

if (FALSE) {
  q <- cowplot::plot_grid(
    labels = c("A", "B"),
    plotlist = list(q0, q1),
    ncol = 1
  )
  ggsave(
    plot = q,
    file = "figures/augmentation_re.pdf",
    width = 9.0,
    height = 9.0
  )
}

tab_rel <- df_rel %>%
  dplyr::mutate(
    re = sprintf("%.2f", rel_eff),
    lo = sprintf("%.2f", lower),
    hi = sprintf("%.2f", upper),
    re_ci = glue::glue("{re} ({lo}-{hi})"),
  ) %>%
  dplyr::select(n_factor, covar_effect, re_ci)
tab_rel <- xtable::xtable(tab_rel)
xtable::print.xtable(tab_rel, include.rownames = FALSE)

# -------------------------------------------------------------------------
# Experiment 4: Frailty.
# -------------------------------------------------------------------------

data <- readRDS(file = "data/E4_sim.rds")

data$n_factor <- factor(
  x = data$n,
  levels = c(50, 100, 200, 400),
  ordered = TRUE
)

data$frailty_factor <- factor(
  x = data$frailty,
  levels = c(0, 2, 4, 8),
  ordered = TRUE
)

df <- data %>%
  dplyr::group_by(n_factor, frailty_factor) %>%
  dplyr::summarise(
    ese = sd(observed),
    ase = sqrt(mean(asymptotic_se^2))
  )
show(df)

