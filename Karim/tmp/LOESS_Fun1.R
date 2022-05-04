## FUN 2

get_infection_incidence_by_deconvolution <- function(
  data_subset,
  constant_delay_distribution,
  constant_delay_distribution_incubation = c(),
  is_onset_data = F,
  is_local_cases = T,
  smooth_incidence = T,
  days_incl = 21,
  empirical_delays  = tibble(),
  n_bootstrap = 5,
  days_further_in_the_past = 30,
  days_further_in_the_past_incubation = 5,
  max_iterations = 100,
  verbose = FALSE) {
  
  #TODO make the days_further_in_the_past type specific
  
  if(nrow(data_subset) == 0) {
    return(tibble())
  }
  
  data_type_subset <- unique(data_subset$data_type)[1]
  
  # exclude leading zeroes
  data_subset <- data_subset %>%
    arrange(date) %>%
    filter(cumsum(value) > 0)
  
  if(nrow(data_subset) == 0) {
    return(tibble())
  }
  
  minimal_date <- min(data_subset$date) - days_further_in_the_past
  maximal_date <- max(data_subset$date)
  all_dates <- seq(minimal_date, maximal_date, by = "days")
  
  is_empirical = (nrow(empirical_delays) > 0)
  
  if(verbose && is_empirical) {
    cat("\tEmpirical delay distribution available\n")
  }
  
  if( is_onset_data ) {
    delay_distribution_matrix_incubation <- get_matrix_constant_waiting_time_distr(
      constant_delay_distribution_incubation,
      all_dates)
    
    initial_delta_incubation <- min(which(cumsum(constant_delay_distribution_incubation) > 0.5)) - 1 # take median value (-1 because index 1 corresponds to zero days)
    
    
    if(unique(data_subset$region)[1] != "ESP") { # hack to workaround weirdness of Spanish data
      
      # account for additional right-truncation of onset data (needs to be reported first)
      if(is_empirical) {
        delay_distribution_matrix_onset_to_report <- get_matrix_empirical_waiting_time_distr(
          empirical_delays,
          seq.Date(min(data_subset$date), max(data_subset$date), by = "days"))
      } else {
        delay_distribution_matrix_onset_to_report <- get_matrix_constant_waiting_time_distr(
          constant_delay_distribution,
          seq.Date(min(data_subset$date), max(data_subset$date), by = "days"))
      }
      
      data_subset <- data_subset %>%
        complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0))
      
      Q_vector_onset_to_report <- apply(delay_distribution_matrix_onset_to_report, MARGIN = 2, sum)
      
      #TODO remove
      # if(unique(data_subset$region)[1] == "ESP") { # hack to work around spanish data between symptom onset dates only
      #   right_truncation <- 3
      #   # need to offset the Q vector by how many days were truncated off originally
      #   Q_vector_onset_to_report <- c(rep(1, right_truncation), Q_vector_onset_to_report[1:(length(Q_vector_onset_to_report) - right_truncation)] )
      # }
      
      data_subset <- data_subset %>%
        mutate(value = value / Q_vector_onset_to_report) %>% 
        mutate(value = if_else(value == Inf, 0, value))
      
    }
    
  } else {
    if(is_empirical) {
      delay_distribution_matrix_onset_to_report <- get_matrix_empirical_waiting_time_distr(
        empirical_delays,
        all_dates[(days_further_in_the_past_incubation + 1):length(all_dates)])
      
      delay_distribution_matrix_incubation <- get_matrix_constant_waiting_time_distr(
        constant_delay_distribution_incubation,
        all_dates)
      
      initial_delta_incubation <- min(which(cumsum(constant_delay_distribution_incubation) > 0.5)) - 1 # take median value (-1 because index 1 corresponds to zero days)
      initial_delta_report <-  median(empirical_delays$delay, na.rm = T)
    } else {
      delay_distribution_matrix <- get_matrix_constant_waiting_time_distr(
        constant_delay_distribution,
        all_dates)
      
      initial_delta <- min(which(cumsum(constant_delay_distribution) > 0.5)) - 1 # take median value (-1 because index 1 corresponds to zero days)
    }
  }
  
  
  
  results <- list(tibble())
  
  for (bootstrap_replicate_i in 0:n_bootstrap) {
    
    if (verbose == T) {
      cat("    Bootstrap replicate: ", bootstrap_replicate_i, "\n")
    }
    
    if (bootstrap_replicate_i == 0) {
      time_series <- data_subset
    } else {
      time_series <- get_bootstrap_replicate(data_subset)
    }
    
    if (smooth_incidence == T) {
      smoothed_incidence_data <- time_series %>%
        complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0)) %>% 
        mutate(value = getLOESSCases(dates = date, count_data = value, days_incl))
      
      raw_total_incidence <- sum(time_series$value, na.rm = TRUE)
      smoothed_total_incidence <- sum(smoothed_incidence_data$value, na.rm = T)
      
      if (smoothed_total_incidence > 0) {
        smoothed_incidence_data <- smoothed_incidence_data %>%
          mutate(value = value * raw_total_incidence / smoothed_total_incidence)
      }
      
    } else {
      smoothed_incidence_data <- time_series  %>%
        complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0))
    }
    
    
    if (is_onset_data) {
      deconvolved_infections <-  do_deconvolution(smoothed_incidence_data,
                                                  delay_distribution_matrix = delay_distribution_matrix_incubation,
                                                  days_further_in_the_past = days_further_in_the_past,
                                                  initial_delta = initial_delta_incubation,
                                                  max_iterations = max_iterations,
                                                  verbose = verbose)
    } else {
      if(is_empirical) {
        # perform the deconvolution in two steps
        deconvolved_symptom_onsets <- do_deconvolution(smoothed_incidence_data,
                                                       delay_distribution_matrix = delay_distribution_matrix_onset_to_report,
                                                       days_further_in_the_past = days_further_in_the_past - days_further_in_the_past_incubation,
                                                       initial_delta = initial_delta_report,
                                                       max_iterations = max_iterations,
                                                       verbose = verbose)
        
        deconvolved_infections <- do_deconvolution(deconvolved_symptom_onsets,
                                                   delay_distribution_matrix = delay_distribution_matrix_incubation,
                                                   days_further_in_the_past = days_further_in_the_past_incubation,
                                                   initial_delta = initial_delta_incubation,
                                                   max_iterations = max_iterations,
                                                   verbose = verbose)
      } else {
        deconvolved_infections <-  do_deconvolution(smoothed_incidence_data,
                                                    delay_distribution_matrix = delay_distribution_matrix,
                                                    days_further_in_the_past = days_further_in_the_past,
                                                    initial_delta = initial_delta,
                                                    max_iterations = max_iterations,
                                                    verbose = verbose)
      }
    }
    
    
    deconvolved_infections <- deconvolved_infections %>% slice((days_further_in_the_past -5 + 1):n())
    
    data_type_name <- paste0("infection_", data_type_subset)
    
    ## dataframe containing results
    deconvolved_infections <- tibble(
      date = deconvolved_infections$date,
      region = unique(time_series$region)[1],
      country = unique(time_series$country)[1],
      source = unique(time_series$source)[1],
      local_infection = is_local_cases,
      data_type = data_type_name,
      replicate = bootstrap_replicate_i,
      value = deconvolved_infections$value
    )
    
    results <- c(results, list(deconvolved_infections))
  }
  
  return(bind_rows(results))
}
