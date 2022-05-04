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

all_files <- list.files(path = "data", pattern = ".csv$")
files <- all_files[2:(length(all_files)-1)]
file_path <- paste0("data/", files)

data <- file_path %>%
  purrr::set_names(nm = (basename(.) %>% tools::file_path_sans_ext())) %>%
  purrr::map_df(function(x){read.csv(x, sep = ";")}, .id = "FileName") %>% 
  mutate(rna_v = str_extract(FileName, "_v\\d$")) %>% 
  mutate(location = str_extract(FileName,"(?<=processed_normed_data_)\\w*(?=_?)")) %>% 
  mutate(location = str_replace(FileName, "_v\\d$", "")) %>% 
  mutate(location = str_replace(location, "processed_normed_data_", "")) %>% 
  mutate(rna_v = str_replace(rna_v, "_", "")) %>% 
  filter(location != "lausanne") %>% 
  mutate(date = as.Date(X)) %>% 
  select(-c(FileName, X)) 


# clean variable names
names(data) <- str_extract(names(data), "^\\w*")


# plot data
coeff <- max(data$median_7d_sars_cov2_rna, na.rm = TRUE)/max(data$median_7d_new_cases, na.rm = TRUE)


# reproducing plots with missing values

ggplot(data, aes(x = date)) +
  geom_line(aes(y = median_7d_sars_cov2_rna/coeff, colour = rna_v)) +
  geom_line(aes(y = median_7d_new_cases, colour = "New cases")) +
  scale_y_continuous(
    
    # Features of the first axis
    name = "7 day Median of new infections",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="RNA concentration in Wasterwater")
  ) +
  facet_wrap(~location, nrow = 3)
  

# now, impute missing values for new_cases with 7-day moving average

data <- data %>% 
  arrange(rna_v, location, date) %>% 
  group_by(location, rna_v) %>% 
  mutate(new_cases_imp = na_ma(new_cases, k = 1)) %>% 
  mutate(cases_7d_med = rollmedian(new_cases_imp,7, fill = NA, align = "right"))


ggplot(data, aes(x = date)) +
  geom_line(aes(y = cases_7d_med, colour = "Imputed")) +
  geom_line(aes(y = median_7d_new_cases, colour = "New cases")) +
  facet_wrap(~location, nrow = 3)


data_imp2 <- data %>% 
  arrange(rna_v, location, date) %>% 
  group_by(location, rna_v) %>% 
  mutate(cases_7d_med = rollmedian(new_cases_imp,7, fill = NA, align = "right")) %>% 
  mutate(cases_7d_med = na_ma(cases_7d_med, k = 1)) 




test <- data %>% 
  select(date, location, rna_v, new_cases, new_cases_imp, median_7d_new_cases, cases_7d_med)



ggplot(data, aes(x = date)) +
  geom_line(aes(y = median_7d_sars_cov2_rna/coeff, colour = rna_v)) +
  geom_line(aes(y = median_7d_new_cases, colour = "New cases")) +
  scale_y_continuous(
    
    # Features of the first axis
    name = "7 day Median of new infections",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="RNA concentration in Wasterwater")
  ) +
  facet_wrap(~location, nrow = 3)



# counting missing values and range

data_v1 <- data %>% 
  filter(rna_v == "v1")


data_v2 <- data %>% 
  filter(rna_v == "v2")

data_v1 %>% 
  group_by(location) %>% 
  summarize(nobs_cases = sum(!is.na(new_cases)),
            nobs_rna = sum(!is.na(sars_cov2_rna)),
            start_date = min(date, na.rm = TRUE),
            end_date = max(date, na.rm = TRUE))

data_v2 %>% 
  group_by(location) %>% 
  summarize(nobs_cases = sum(!is.na(new_cases)),
            nobs_rna = sum(!is.na(sars_cov2_rna)),
            start_date = min(date, na.rm = TRUE),
            end_date = max(date, na.rm = TRUE))



data_v1 %>% 
  filter(!is.na(sars_cov2_rna)) %>% 
  group_by(location) %>% 
  summarize(nobs_cases = sum(!is.na(new_cases)),
            nobs_rna = sum(!is.na(sars_cov2_rna)),
            start_date = min(date, na.rm = TRUE),
            end_date = max(date, na.rm = TRUE))

data_v2 %>% 
  filter(!is.na(sars_cov2_rna)) %>% 
  group_by(location) %>% 
  summarize(nobs_cases = sum(!is.na(new_cases)),
            nobs_rna = sum(!is.na(sars_cov2_rna)),
            start_date = min(date, na.rm = TRUE),
            end_date = max(date, na.rm = TRUE))



# Plot lagged 7 day median of new cases and new cases

ggplot(data_v1, aes(x = date)) +
  geom_segment( aes(x = date, xend = date, y=0, yend=new_cases,colour = "New cases"), size = 0.1) +
  #geom_point(aes(y = new_cases, colour = "New cases"), size = 0.1) +
  geom_line(aes(y = median_7d_new_cases, colour = "New cases 7 day median (lagged)")) +
  facet_wrap(~location, nrow = 3)

smooth_data <- data_v1 %>% 
  mutate(center_med = rollmedian(new_cases, k = 7, align = "center", fill = NA))

ggplot(smooth_data, aes(x = date)) +
  geom_segment( aes(x = date, xend = date, y=0, yend=new_cases,colour = "New cases"), size = 0.1) +
  #geom_point(aes(y = new_cases, colour = "New cases"), size = 0.1) +
  geom_line(aes(y = median_7d_new_cases, colour = "New cases 7 day median (lagged)"), size = 0.4) +
  geom_line(aes(y = center_med, colour = "New cases 7 day median (centered)"), size = 0.4, linetype = "twodash") +
  facet_wrap(~location, nrow = 3)


# Zurich only
smooth_data %>% 
  filter(location == "zurich") %>% 
  arrange(date) %>% 
  mutate(supersmooth = supsmu(x = date, y = new_cases, span = "cv")) %>% 
  ggplot(aes(x = date)) +
  geom_segment( aes(x = date, xend = date, y=0, yend=new_cases,colour = "New cases"), size = 0.3) +
  #geom_point(aes(y = new_cases, colour = "New cases"), size = 0.1) +
  geom_line(aes(y = median_7d_new_cases, colour = "New cases 7 day median (lagged)"), size = 0.6) +
  geom_line(aes(y = center_med, colour = "New cases 7 day median (centered)"), size = 0.4, linetype = "twodash") +
  geom_line(aes(y = supersmooth, colour = "Super Smoother"), size = 0.4, linetype = "twodash") +
  facet_wrap(~location) +
  scale_color_manual(values=c("#999999", "red", "green"))


