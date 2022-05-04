## variant proportions

### load packages
library(tidyverse)
library(ggplot2)
library(stats)
library(zoo)
library(imputeTS)

library(devtools)
install_github("covid-19-Re/estimateR")

library(estimateR)

# import data
df_ww <- read_csv(file = "data/variant_proportions_wastewater.csv") %>% 
  select(- ...1) %>% 
  rename(date = index)

# reshape data from wide to long for better plotting
df_ww_long <- df_ww %>% 
  gather(variant, proportion, -c(date, city))

# plot data
df_ww_long %>% 
  ggplot() +
  geom_line(aes(x = date, y = proportion, color = variant)) +
  facet_wrap(~ city, ncol = 3)

# compute max. proportion per variant x city to filter out irrelavant variants
max_prop <- df_ww_long %>% 
  group_by(city, variant) %>% 
  arrange(date) %>% 
  summarise(max_prop = max(proportion, na.rm = TRUE))

# create a filter vector
filter_vec <- max_prop %>% 
  filter(max_prop >= 0.75) %>% 
  mutate(filter_vec = paste0(city,"_",variant)) %>% 
  ungroup() %>% 
  select(filter_vec) %>% 
  pull()

# filter data for plotting
df_plot <- df_ww_long %>% 
  filter(paste0(city,"_",variant) %in% filter_vec)

df_plot %>% 
  ggplot() +
  geom_line(aes(x = date, y = proportion, color = variant)) +
  facet_wrap(~ city, ncol = 3)


  