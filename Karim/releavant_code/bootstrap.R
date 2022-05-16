get_bootstrap_replicates <- function(
  data, 
  block_size = 10, 
  days_incl = 21,
  integer_results = FALSE,
  days_in_past = 30
  ) { # to sample one time series based on LOESS and block bootstrap
  
  wastewater_data <- data %>% 
    filter(datatype == "wastewater")
  
  used_dates <- wastewater_data %>% 
    pull(date)
  
  first_date_inc <- used_dates[1] - days_in_past
  
  if(days_in_past == 0){
    
    dates_in_past <- c()
    
  } else {
    
    dates_in_past <- seq(first_date_inc, used_dates[1] - 1, by = "days")
    
  }
  
  
  
  all_dates <- c(dates_in_past,used_dates) %>% 
    unique()
  
  incidence_data <- data %>% 
    filter(datatype == "confirmed cases",
           date %in% all_dates)
  
  # incidence_data_all <- data %>% 
  #   filter(datatype == "confirmed cases",
  #          date %in% all_dates)
  
  # incidence_data_past <- data %>% 
  #   filter(datatype == "confirmed cases",
  #          date %in% dates_in_past)
  
  tmp1 <- incidence_data
  tmp1$log_value <- log(tmp1$value + 1)
  
  tmp2 <- wastewater_data
  tmp2$log_value <- log(tmp2$value + 1)
  
  # Smoothing
  smoothed_incidence_data_all <- tmp1 %>%
    complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(log_value = 0)) %>%
    mutate(log_loess = getLOESSCases(dates = date, count_data = log_value, days_incl), # getLOESSCases returns a vector
           log_diff = log_value - log_loess)
  
  smoothed_incidence_data <- smoothed_incidence_data_all %>% 
    filter(date %in% used_dates)
  
  smoothed_wastewater_data <- tmp2 %>%
    complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(log_value = 0)) %>%
    mutate(log_loess = getLOESSCases(dates = date, count_data = log_value, days_incl), # getLOESSCases returns a vector
           log_diff = log_value - log_loess)
  
  ts_case <- smoothed_incidence_data$log_diff
  ts_ww <- smoothed_wastewater_data$log_diff
  
  
  log_smoothed_cases <- smoothed_incidence_data$log_loess
  log_smoothed_ww <- smoothed_wastewater_data$log_loess
  
  
  ## Past:
  smoothed_incidence_data_past <- smoothed_incidence_data_all %>% 
    filter(date %in% dates_in_past)
  
  ts_case_past <- smoothed_incidence_data_past$log_diff
  log_smoothed_cases_past <- smoothed_incidence_data_past$log_loess
  ###
    
    
    # get the weekdays for each position of ts
    weekdays_index <- (1:length(ts_case)) %% 7
    weekdays_index[which(weekdays_index==0)] <- 7
    
    
    
    
    ts_case_boot <-c()
    ts_ww_boot <-c()
    last_day_index <- 7
    
    ###### get the ts_boot: make sure glue wrt the correct days
    while(length(ts_case_boot) < length(ts_case)){
      start_index <- sample(1:(length(ts_case)-block_size+1), 1)
      sampled_index <- start_index:(start_index+block_size-1)
      sampled_weekdays <- weekdays_index[sampled_index] 
      
      # make sure the day related to the first sample is after the previous ts_boot
      first_day_index <- which(sampled_weekdays==last_day_index)[1] + 1
      ts_boot_index <- sampled_index[first_day_index:block_size]
      
      last_day_index <- tail(weekdays_index[ts_boot_index],1)
      
      ts_case_boot <- c(ts_case_boot, ts_case[ts_boot_index])
      ts_ww_boot <- c(ts_ww_boot, ts_ww[ts_boot_index])
    }
    
    # take the same length as previous ts
    
    ts_case_boot <- ts_case_boot[1:length(ts_case)]
    ts_ww_boot <- ts_ww_boot[1:length(ts_ww)]
    
    
    
    
    # return(list("cases" = ts_case_boot, "wastewater" = ts_ww_boot))  
    
    # Sampling 
    # log_diff_boot <- block_boot_overlap_func(smoothed_incidence_data$log_diff, block_size)
    
    
    ts_case_boot_out <- exp(ts_case_boot + log_smoothed_cases) -1
    ts_case_boot_out[ts_case_boot_out<0] <- 0
    
    
    ts_ww_boot_out <- exp(ts_ww_boot + log_smoothed_ww) -1
    ts_ww_boot_out[ts_ww_boot_out<0] <- 0
    

    
    if (days_in_past != 0){
      ##### Now, past
      weekdays_index_past <- ((1:length(ts_case_past)) - length(ts_case_past) %% 7) %% 7
      weekdays_index_past[which(weekdays_index_past==0)] <- 7
      
      ts_case_boot_past <-c()
      
      last_day_index_past <- ifelse(weekdays_index_past[1]-1 == 0,
                                    7,
                                    weekdays_index_past[1]-1)
      
      
      
      
      while(length(ts_case_boot_past) < length(ts_case_past)){
        start_index <- sample(1:(length(ts_case_past)-block_size+1), 1)
        sampled_index <- start_index:(start_index+block_size-1)
        sampled_weekdays <- weekdays_index_past[sampled_index] 
        
        # make sure the day related to the first sample is after the previous ts_boot
        first_day_index <- which(sampled_weekdays==last_day_index_past)[1] + 1
        ts_boot_index <- sampled_index[first_day_index:block_size]
        
        last_day_index_past <- tail(weekdays_index_past[ts_boot_index],1)
        
        ts_case_boot_past <- c(ts_case_boot_past, ts_case_past[ts_boot_index])
        
      }
      
      
      # take the same length as previous ts
      ts_case_boot_past <- ts_case_boot_past[1:length(ts_case_past)]
      
    
      
      ts_case_boot_past_out <- exp(ts_case_boot_past + log_smoothed_cases_past) -1
      ts_case_boot_past_out[ts_case_boot_past_out<0] <- 0
      
      ts_case_boot_out <- c(ts_case_boot_past_out,ts_case_boot_out) 
      
    }
    
    
    if (integer_results){
      
      ts_case_boot_out <- round(ts_case_boot_out)  
      ts_ww_boot_out <- round(ts_ww_boot_out)  
      
    }
    
    replicate_cases <- incidence_data %>%
      complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0)) %>%
      dplyr::mutate(value = ts_case_boot_out) %>%
      arrange(date)
    
    
    replicate_wastewater <- wastewater_data %>%
      complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0)) %>%
      dplyr::mutate(value = ts_ww_boot_out) %>%
      arrange(date)
    
    results <- bind_rows(replicate_cases, replicate_wastewater)
    
    if((replicate_wastewater %>% 
          pull(date) 
              == 
        replicate_cases[(1+days_in_past):length(replicate_cases$date),"date"] %>% 
          pull(date)) %>% 
        all()) {
      
      return(results) 
      
    }
  
}

"""
## To Do: 
- reply to Martin
- potentially make the blockbootstrapping modular function again...
- what data do we want to use?
- absolute numbers and then log??
- 


to do:

- multiply incidence by 100,000
- deliver full wastewater data


"""


data.frame(x = 1:50) %>% 
  slice((10 -5 + 1):n())






















######


# get_bootstrap_replicates <- function(
#   incidence_data, wastewater_data, 
#   block_size = 10, 
#   days_incl = 21,
#   integer_results = FALSE) { # to sample one time series based on LOESS and block bootstrap
#   
#   tmp1 <- incidence_data
#   tmp1$log_value <- log(tmp1$value + 1)
#   
#   tmp2 <- wastewater_data
#   tmp2$log_value <- log(tmp2$value + 1)
#   
#   # Smoothing
#   smoothed_incidence_data <- tmp1 %>%
#     complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(log_value = 0)) %>%
#     mutate(log_loess = getLOESSCases(dates = date, count_data = log_value, days_incl), # getLOESSCases returns a vector
#            log_diff = log_value - log_loess)
#   
#   smoothed_wastewater_data <- tmp2 %>%
#     complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(log_value = 0)) %>%
#     mutate(log_loess = getLOESSCases(dates = date, count_data = log_value, days_incl), # getLOESSCases returns a vector
#            log_diff = log_value - log_loess)
#   
#   ts_case <- smoothed_incidence_data$log_diff
#   ts_ww <- smoothed_wastewater_data$log_diff
#   
#   log_smoothed_cases <- smoothed_incidence_data$log_loess
#   log_smoothed_ww <- smoothed_wastewater_data$log_loess
#   
#   
#   
#   
#   if(length(ts_case) == length(ts_ww)) {
#   
#     
#     # get the weekdays for each position of ts
#   weekdays_index <- (1:length(ts_case)) %% 7
#   weekdays_index[which(weekdays_index==0)] <- 7
#     
#     ts_case_boot <-c()
#     ts_ww_boot <-c()
#     last_day_index <- 7
#     
#     ###### get the ts_boot: make sure glue wrt the correct days
#     while(length(ts_case_boot) < length(ts_case)){
#       start_index <- sample(1:(length(ts_case)-block_size+1), 1)
#       sampled_index <- start_index:(start_index+block_size-1)
#       sampled_weekdays <- weekdays_index[sampled_index] 
#       
#       # make sure the day related to the first sample is after the previous ts_boot
#       first_day_index <- which(sampled_weekdays==last_day_index)[1] + 1
#       ts_boot_index <- sampled_index[first_day_index:block_size]
#       
#       last_day_index <- tail(weekdays_index[ts_boot_index],1)
#       
#       ts_case_boot <- c(ts_case_boot, ts_case[ts_boot_index])
#       ts_ww_boot <- c(ts_ww_boot, ts_ww[ts_boot_index])
#     }
#     
#     # take the same length as previous ts
#     
#     ts_case_boot <- ts_case_boot[1:length(ts_case)]
#     ts_ww_boot <- ts_ww_boot[1:length(ts_ww)]
#     
#     # return(list("cases" = ts_case_boot, "wastewater" = ts_ww_boot))  
#     
#   # Sampling 
#   # log_diff_boot <- block_boot_overlap_func(smoothed_incidence_data$log_diff, block_size)
#   
#   
#   ts_case_boot_out <- exp(ts_case_boot + log_smoothed_cases) -1
#   ts_case_boot_out[ts_case_boot_out<0] <- 0
#   
#   
#   ts_ww_boot_out <- exp(ts_ww_boot + log_smoothed_ww) -1
#   ts_ww_boot_out[ts_ww_boot_out<0] <- 0
#   
#   
#   if (integer_results){
#     
#     ts_case_boot_out <- round(ts_case_boot_out)  
#     ts_ww_boot_out <- round(ts_ww_boot_out)  
#     
#   }
#   
#   replicate_cases <- incidence_data %>%
#     complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0)) %>%
#     dplyr::mutate(value = ts_case_boot_out) %>%
#     arrange(date)
#   
#   
#   replicate_wastewater <- wastewater_data %>%
#     complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0)) %>%
#     dplyr::mutate(value = ts_case_boot_out) %>%
#     arrange(date)
#   
#   if((replicate_wastewater %>% pull(date) == replicate_cases %>% pull(date)) %>% 
#      all()) {
#     
#     return(list("cases" = replicate_cases, "wastewater" = replicate_wastewater)) 
#     
#   }
#   
#   
#   
#   
#   } else {
#     
#     stop("Input time series are not of same length!")
#     
#   }
#   
# }
# 
# 
# 
# 
# 
# 
