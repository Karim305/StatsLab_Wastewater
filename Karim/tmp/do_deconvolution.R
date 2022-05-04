# do_deconvolution.R function

function(
  incidence_data,
  days_further_in_the_past = 30,
  verbose = FALSE,
  delay_distribution_matrix,
  initial_delta,
  max_iterations = 100
) {
  
  # use mode of 'constant_delay_distribution'. -1 because indices are offset by one as the delay can be 0.
  
  first_guess_delay <- ceiling(initial_delta)
  
  if (verbose) {
    cat("\tDelay on first guess: ", first_guess_delay, "\n")
  }
  
  first_recorded_incidence <-  with(filter(incidence_data, cumsum(value) > 0), value[which.min(date)])
  last_recorded_incidence <- with(incidence_data, value[which.max(date)])
  
  minimal_date <- min(incidence_data$date) - days_further_in_the_past
  maximal_date <- max(incidence_data$date)
  
  first_guess <- incidence_data %>%
    mutate(date = date - first_guess_delay) %>%
    complete(date = seq.Date(minimal_date, min(date), by = "days"),
             fill = list(value = first_recorded_incidence)) %>% # left-pad with first recorded value
    complete(date = seq.Date(max(date), maximal_date, by = "days"),
             fill = list(value = last_recorded_incidence)) %>% # right-pad with last recorded value
    arrange(date) %>% 
    filter(date >=  minimal_date)
  
  original_incidence <- incidence_data %>% 
    complete(date = seq.Date(minimal_date, maximal_date, by = "days"),
             fill = list(value = 0)) %>% 
    pull(value)
  
  final_estimate <- iterate_RL(
    first_guess$value,
    original_incidence,
    delay_distribution_matrix = delay_distribution_matrix,
    max_delay = days_further_in_the_past,
    max_iterations = max_iterations,
    verbose = verbose)
  
  deconvolved_dates <- first_guess %>% pull(date)
  
  result <- tibble(date = deconvolved_dates, value = final_estimate)
  
  result <- result %>%
    filter(date <= maximal_date - first_guess_delay)
  
  return(result)
}
