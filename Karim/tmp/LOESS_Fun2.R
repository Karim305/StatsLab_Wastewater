get_bootstrap_replicate <- function(original_time_series, block_size = 10, days_incl = 21) {
  tmp <- original_time_series
  
  # Change introduced after meeting on 19.1
  #tmp$log_value <- ifelse(tmp$value != 0, log(tmp$value), 0)
  tmp$log_value <- log(tmp$value + 1)
  
  smoothed_incidence_data <- tmp %>%
    complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(log_value = 0)) %>%
    mutate(log_loess = getLOESSCases(dates = date, count_data = log_value, days_incl),
           log_diff = log_value - log_loess)
  
  # The bootstrap is performed on log_diff = log_value - log_loess. The assumption that the deviations from the smoother are iid is quite reasonable
  log_diff_boot <- block_boot_overlap_func(smoothed_incidence_data$log_diff, block_size)
  log_smoothed_data <- smoothed_incidence_data$log_loess
  
  ts_boot <- exp(log_diff_boot + log_smoothed_data) -1
  ts_boot[ts_boot<0] <- 0
  ts_boot <- round(ts_boot)
  
  replicate <- original_time_series %>%
    complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0)) %>%
    dplyr::mutate(value = ts_boot) %>%
    arrange(date)
  
  return(replicate)
}

### 
block_boot_overlap_func <- function(ts, block_size = 10){
  
  # get the weekdays for each position of ts
  weekdays_index <- (1:length(ts)) %% 7
  weekdays_index[which(weekdays_index==0)] <- 7
  
  ts_boot <-c()
  last_day_index <- 7
  
  ###### get the ts_boot: make sure glue wrt the correct days
  while(length(ts_boot) < length(ts)){
    start_index <- sample(1:(length(ts)-block_size+1), 1) # sample start index somewhere in the data
    sampled_index <- start_index:(start_index+block_size-1) 
    sampled_weekdays <- weekdays_index[sampled_index]
    
    # make sure the day related to the first sample is after the previous ts_boot
    first_day_index <- which(sampled_weekdays==last_day_index)[1] + 1
    ts_boot_index <- sampled_index[first_day_index:block_size]
    
    last_day_index <- tail(weekdays_index[ts_boot_index],1)
    
    ts_boot <- c(ts_boot, ts[ts_boot_index])
  }
  
  # take the same length as previous ts
  ts_boot <- ts_boot[1:length(ts)]
  
  return(ts_boot)
}
