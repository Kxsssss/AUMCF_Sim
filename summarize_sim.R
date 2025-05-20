my_data <- readRDS("~/Desktop/power_sim.rds")

print(my_data %>% 
  group_by(type, n) %>%
  summarise(mean_pvalue = mean(p_value < 0.05)), n = 50)

my_data <- readRDS("~/Desktop/validity_sim.rds")

print(my_data %>% 
        group_by(type, n) %>%
        summarise(mean_pvalue = mean(p_value < 0.05)), n = 50)


dim(my_data %>% filter(type == 'wr'))

my_data %>% filter(type == 'wr') %>%
  filter(p_value <= 0.05 & upper >= 1 & lower <= 1)