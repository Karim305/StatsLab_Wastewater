# Simulate performance of seasonal imputation algorithms

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


zh_ts <- data_zh$new_cases[!is.na(data_zh$new_cases)]




## Define simulations imputations that should be simulated using replicate() as function

sim_RMSE <-function(num_NA = 20, data = zh_ts){
  
  NA_rows <- sample(10:500, num_NA, replace=FALSE) %>% 
    unique() %>% 
    sort()
  
  gaps_vec <- data 
  gaps_vec[NA_rows] <- NA
  
  zh_gaps <- ts(gaps_vec, frequency = 7) 
  

  ### Time-Series specific method      
  #- Last observation carried forward (LOCF)     
  #- Next observation carried backward (NOCB)   
  #- Linear interpolation     
  #- Spline interpolation
  #### These methods rely on the assumption that adjacent observations are similar to one another. These methods do not work well when this assumption is not valid, especially when the presence of strong seasonality.
  
  # Last Observartion Carried Forward
  m1 <- mean((na_locf(zh_gaps, option = "locf") - data)^2)
  
  # Next Observartion Carried Backward
  m2 <- mean((na_locf(zh_gaps, option = "nocb") - data)^2)
  
  # Linear Interpolation
  m3 <- mean((na_interpolation(zh_gaps, option = "linear") - data)^2)
  
  # Spline Interpolation
  m4 <- mean((na_interpolation(zh_gaps, option = "spline") - data)^2)
  
  TS1 <- t(data.frame(row.names =c('LOCF', 'NOCB', 'Linear', 'Spline'), MSE=c(m1, m2, m3, m4)))
  
  
  
  ### The Combination of Seasonal Adjustment and other methods      
  #### This approach is very effective when seasonality is present.      
  # We de-seasonlize the data first, and then do interpolation on the data.    
  # Once the mising values are inputed, we need to re-seasonalize the data. 
  
  # Seasonal Adjustment then Random
  ma1 <- mean((na_seadec(zh_gaps, algorithm = "random") - data)^2)
  
  # Seasonal Adjustment then Mean
  ma2 <- mean((na_seadec(zh_gaps, algorithm = "mean") - data)^2)
  
  # Seasonal Adjustment then LOCF
  ma3 <- mean((na_seadec(zh_gaps, algorithm = "locf") - data)^2)
  
  # Seasonal Adjustment then Linear Interpolation
  ma4 <- mean((na_seadec(zh_gaps, algorithm = "interpolation") - data)^2)
  
  # Seasonal Adjustment then Kalman smoothing (State space model)
  ma5 <- mean((na_seadec(zh_gaps, algorithm = "kalman") - data)^2)
  
  TS2 <- t(data.frame(row.names = c("Seas-Adj+Random", "Seas-Adj+Mean", "Seas-Adj+LOCF","Seas-Adj+Linear", "Seas-Adj+Kalman"),
                    MSE=c(ma1, ma2, ma3, ma4, ma5)))
  
  cbind(TS1,TS2)
  
}

test <- replicate(100, sim_RMSE(num_NA = 50, data = zh_ts), simplify = "matrix")

df <- t(data.frame(test,
                 row.names = c('LOCF', 'NOCB', 'Linear', 'Spline', "Seas-Adj+Random", "Seas-Adj+Mean", "Seas-Adj+LOCF","Seas-Adj+Linear", "Seas-Adj+Kalman")))

# get data frame in long format to plot histograms with ggplot
for (i in seq_along(colnames(df))){
  var_name <- colnames(df)[i]
  if(i == 1){
    
    df_long <- NULL
    
  }
  
  df_i <- data.frame(MSE = df[,i],
                       Method = var_name,
                       row.names = NULL)
  
  df_long <- rbind(df_long, df_i)
  
}

(summary <- df_long %>% 
  group_by(Method) %>% 
  summarize(mean = mean(MSE)))

mean_plot <- summary %>% 
  filter(Method != "Seas-Adj+Random") %>%
  filter(Method != "Seas-Adj+Mean") 

df_long %>% 
  filter(Method != "Seas-Adj+Random") %>%
  filter(Method != "Seas-Adj+Mean") %>%
  ggplot() + 
  geom_histogram(aes(MSE, fill = Method), bins = 100) +
  facet_wrap(~Method) +
  geom_vline(data = mean_plot, aes(xintercept = mean))


df_long %>% 
  filter(Method != "Seas-Adj+Random") %>%
  filter(Method != "Seas-Adj+Mean") %>%
  ggplot() + 
  geom_density(aes(MSE, fill = Method)) +
  facet_wrap(~Method) +
  geom_vline(data = mean_plot, aes(xintercept = mean))

