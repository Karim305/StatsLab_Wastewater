# Imputing missing values

# load packages
library(imputeTS)
library(tidyverse)
library(xts)

# load data
all_files <- list.files(path = "data", pattern = ".csv$")
files <- all_files[2:length(all_files)]
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
  select(-c(FileName, X)) %>% 
  filter(rna_v == "v1")

# clean variable names
names(data) <- str_extract(names(data), "^\\w*")

# obtain time series for Zurich
data_zh <- data %>% 
  filter(location == "zurich") %>% 
  select(date, new_cases)

data_zh %>% 
  filter(is.na(new_cases))

data_zh %>% 
  filter(!is.na(new_cases)) %>% 
  arrange(date) %>% 
  mutate(gap = date - lag(date,1)) %>% 
  arrange(desc(gap))

data_zh$new_cases[1:length(data_zh$new_cases) -1] %>% 
  is.na() %>% 
  any()

set.seed(1)
NA_rows <- sample(10:500, 20, replace=FALSE) %>% 
  unique() %>% 
  sort()

data_zh <- data_zh %>% 
  filter(!is.na(new_cases)) %>% 
  mutate(new_cases_gaps = new_cases)

data_zh$new_cases_gaps[NA_rows] <- NA

zh_gaps <- ts(data_zh$new_cases_gaps, frequency = 7) 

zh_ts <- data_zh$new_cases[!is.na(data_zh$new_cases)]


# This is the Case data with some observations removed
plot(data_zh[,c("date","new_cases_gaps")], main="New cases ZH with missing values")
statsNA(data_zh$new_cases_gaps)

## non-TS methods

par(mfrow=c(2,2))
# Mean Imputation
plot(na_mean(zh_gaps, option = "mean") - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "Mean")
b1 <- mean((na_mean(zh_gaps, option = "mean") - zh_ts)^2)

# Median Imputation
plot(na_mean(zh_gaps, option = "median") - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "Median")
b2 <- mean((na_mean(zh_gaps, option = "median") - zh_ts)^2)

# Mode Imputation
plot(na_mean(zh_gaps, option = "mode") - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "Mode")
b3 <- mean((na_mean(zh_gaps, option = "mode") - zh_ts)^2)    

# Random Imputation
plot(na_random(zh_gaps) - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "Random")
b4 <- mean((na_random(zh_gaps) - zh_ts)^2)

NoTS <- data.frame(methods=c('Mean Imputation', 'Median Imputation', 'Mode Imputation', 'Random Imputation'), MSE=c(b1, b2, b3, b4))
### Time-Series specific method      
#- Last observation carried forward (LOCF)     
#- Next observation carried backward (NOCB)   
#- Linear interpolation     
#- Spline interpolation
#### These methods rely on the assumption that adjacent observations are similar to one another. These methods do not work well when this assumption is not valid, especially when the presence of strong seasonality.      


par(mfrow=c(2,2))
# Last Observartion Carried Forward
plot(na_locf(zh_gaps, option = "locf") - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "LOCF")
m1 <- mean((na_locf(zh_gaps, option = "locf") - zh_ts)^2)

# Next Observartion Carried Backward
plot(na_locf(zh_gaps, option = "nocb") - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "NOCB")
m2 <- mean((na_locf(zh_gaps, option = "nocb") - zh_ts)^2)

# Linear Interpolation
plot(na_interpolation(zh_gaps, option = "linear") - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "Linear")
m3 <- mean((na_interpolation(zh_gaps, option = "linear") - zh_ts)^2)

# Spline Interpolation
plot(na_interpolation(zh_gaps, option = "spline") - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "Spline")
m4 <- mean((na_interpolation(zh_gaps, option = "spline") - zh_ts)^2)

TS1 <- data.frame(methods=c('LOCF', 'NOCB', 'Linear', 'Spline'), MSE=c(m1, m2, m3, m4))



### The Combination of Seasonal Adjustment and other methods      
#### This approach is very effective when seasonality is present.      
# We de-seasonlize the data first, and then do interpolation on the data.    
# Once the mising values are inputed, we need to re-seasonalize the data.      

par(mfrow=c(3,2))
# Seasonal Adjustment then Random
plot(na_seadec(zh_gaps, algorithm = "random") - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "Seas-Adj -> Random")
ma1 <- mean((na_seadec(zh_gaps, algorithm = "random") - zh_ts)^2)

# Seasonal Adjustment then Mean
plot(na_seadec(zh_gaps, algorithm = "mean") - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "Seas-Adj -> Mean")
ma2 <- mean((na_seadec(zh_gaps, algorithm = "mean") - zh_ts)^2)

# Seasonal Adjustment then LOCF
plot(na_seadec(zh_gaps, algorithm = "locf") - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "Seas-Adj -> LOCF")
ma3 <- mean((na_seadec(zh_gaps, algorithm = "locf") - zh_ts)^2)

# Seasonal Adjustment then Linear Interpolation
plot(na_seadec(zh_gaps, algorithm = "interpolation") - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "Seas-Adj -> Linear")
ma4 <- mean((na_seadec(zh_gaps, algorithm = "interpolation") - zh_ts)^2)

# Seasonal Adjustment then Kalman smoothing (State space model)
plot(na_seadec(zh_gaps, algorithm = "kalman") - zh_ts, ylim = c(-mean(zh_ts), mean(zh_ts)), ylab = "Difference", main = "Seas-Adj -> Kalman smoothing")
ma5 <- mean((na_seadec(zh_gaps, algorithm = "kalman") - zh_ts)^2)

TS2 <- data.frame(methods=c("Seas-Adj+Random", "Seas-Adj+Mean", "Seas-Adj+LOCF","Seas-Adj+Linear", "Seas-Adj+Kalman"),
                  MSE=c(ma1, ma2, ma3, ma4, ma5))



