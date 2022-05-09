
library(tidyverse)
library(lubridate)
library(patchwork)
library(viridis)

###########################################################################33
#########################################################################33
# Nested functions 1, don't need editing

### nested functions 1, in get_infection_incidence_by_deconvolution
get_matrix_constant_waiting_time_distr <- function(waiting_time_distr,all_dates) { 
  # waiting_time_distr is a 'pdf', all_dates
  
  #### initialization
  N <- length(all_dates) # length of time series
  if(length(all_dates) >= length(waiting_time_distr)) { # 我们的waiting_time_distr有200天
    waiting_time_distr <- c(waiting_time_distr, rep(0, times = N - length(waiting_time_distr))) # 用0补上, assune probabiliyt = 0
  }
  
  #### draw matrix
  delay_distribution_matrix <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    delay_distribution_matrix[, i ] <-  c(rep(0, times = i - 1 ), waiting_time_distr[1:(N - i + 1)])
  }
  
  return(delay_distribution_matrix)
}

###################################
# we dont even need this anymore, find gamma parameters from mean/sd of distribution
getGammaParams <- function(meanParam, sdParam){
  shapeParam <- meanParam^2 / (sdParam^2)
  scaleParam <- (sdParam^2) / meanParam
  return(list(shape = shapeParam, scale = scaleParam))
}

#########################################################################
#########################################################################
### Nested functions 2, bootstrapping and smoothing related functions
# get_bootstrap_replicate - get_infection_incidence_by_deconvolution
# days_incl is a parameter of getLOESSCases
# block_size is a parameter of block_boot_overlap_funcc
get_bootstrap_replicate <- function(
  original_time_series, 
  block_size = 10, 
  days_incl = 21) { # to sample one time series based on LOESS and block bootstrap
  
  tmp <- original_time_series
  tmp$log_value <- log(tmp$value + 1)
  
  # Smoothing
  smoothed_incidence_data <- tmp %>%
    complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(log_value = 0)) %>%
    mutate(log_loess = getLOESSCases(dates = date, count_data = log_value, days_incl), # getLOESSCases returns a vector
           log_diff = log_value - log_loess)
  
  # Sampling 
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

##################################################################
# smooth time series with LOESS method, get_bootstrap_replicate -  get_infection_incidence_by_deconvolution
getLOESSCases <- function(dates, count_data, days_incl = 21, degree = 1, truncation = 0) {
  
  if (truncation != 0) {
    dates <- dates[1:(length(dates) - truncation)]
    count_data <- count_data[1:(length(count_data) - truncation)]
  }
  
  n_points <- length(unique(dates))
  sel_span <- days_incl / n_points
  
  n_pad <- round(length(count_data) * sel_span * 0.5)
  
  c_data <- data.frame(value = c(rep(0, n_pad), count_data),
                       date_num = c(seq(as.numeric(dates[1]) - n_pad, as.numeric(dates[1]) - 1),
                                    as.numeric(dates)))
  c_data.lo <- loess(value ~ date_num, data = c_data, span = sel_span, degree = degree)
  smoothed <- predict(c_data.lo)
  smoothed[smoothed < 0] <- 0
  raw_smoothed_counts <- smoothed[(n_pad + 1):length(smoothed)]
  normalized_smoothed_counts <-
    raw_smoothed_counts * sum(count_data, na.rm = T) / sum(raw_smoothed_counts, na.rm = T)
  
  if (truncation != 0) {
    normalized_smoothed_counts <- append(normalized_smoothed_counts, rep(NA, truncation))
  }
  return(normalized_smoothed_counts)
}

###########################################################################33
#########################################################################33
### Nested functions 3, convolution related functions
# do_deconvolution - get_infection_incidence_by_deconvolution
do_deconvolution <- function(
  incidence_data,# smoothed incidence data
  days_further_in_the_past, # this is the dimension of the deconvolution function, i.e. how many days after the infection
  #verbose = FALSE, # we dont need this
  delay_distribution_matrix,# calculated matrix of 'probability'
  initial_delta, # this is the mode of the delay distribution, 
  max_iterations = 100, # cap on the number of iterations of RL algorithm
  threshold_chi_squared = threshold_chi_squared
) {
  ######### Initialization
  first_guess_delay <- ceiling(initial_delta) # fixed delay for the initial value
  
  first_recorded_incidence <-  with(filter(incidence_data, cumsum(value) > 0), value[which.min(date)]) 
  last_recorded_incidence <- with(incidence_data, value[which.max(date)])
  
  minimal_date <- min(incidence_data$date) - days_further_in_the_past
  maximal_date <- max(incidence_data$date)
  
  first_guess <- incidence_data %>% # assume fixed delay, and directly shift data, not using deconvolution, described in Deconvolution paper
    mutate(date = date - first_guess_delay) %>% # shift date column
    complete(date = seq.Date(minimal_date, min(date), by = "days"), 
             fill = list(value = first_recorded_incidence)) %>% # left-pad with first recorded value: fill the missing values before the srart of original data set with a fixed number first_recorded_incidence
    complete(date = seq.Date(max(date), maximal_date, by = "days"),
             fill = list(value = last_recorded_incidence)) %>% # right-pad with last recorded value: fill the missing values after the end of original data set
    arrange(date) %>% 
    filter(date >=  minimal_date)
  
  original_incidence <- incidence_data %>% # impute missing values with 0
    complete(date = seq.Date(minimal_date, maximal_date, by = "days"),
             fill = list(value = 0)) %>% 
    pull(value)
  
  ############################ run the RL algorithm
  # nested nested function, iterate_RL
  final_estimate <- iterate_RL( 
    first_guess$value, # shifted smoothed data
    original_incidence,
    delay_distribution_matrix = delay_distribution_matrix,# matrix of pdf
    max_delay = days_further_in_the_past,# 30
    max_iterations = max_iterations,
    threshold_chi_squared = threshold_chi_squared)# 100
  #verbose = verbose)
  
  ###### output the deconvolved data
  deconvolved_dates <- first_guess %>% pull(date)
  result <- tibble(date = deconvolved_dates, value = final_estimate)
  result <- result %>%
    filter(date <= maximal_date - first_guess_delay) # truncate
  return(result)
}


##################################################################################3
#### nested function - do_deconvolution - get_infection_incidence_by_deconvolution
iterate_RL <- function( # Richardson-Lucy EM Algorithm, described in paper. it returns deconvolved incidence data
  initial_estimate,
  original_incidence,
  delay_distribution_matrix,
  threshold_chi_squared,
  max_iterations,
  max_delay) { # hyper parameter, number of days in the delay distribution = days_further_in_the_past
  #verbose = FALSE # we dont need this parameter
  
  current_estimate <- initial_estimate
  N <- length(current_estimate)
  N0 <- N - max_delay
  chi_squared <- Inf
  count <- 1
  
  delay_distribution_matrix <- delay_distribution_matrix[1:length(current_estimate), 1:length(current_estimate)] # truncate delay distribution according to the dimension of data
  truncated_delay_distribution_matrix <- delay_distribution_matrix[(1 + max_delay):NROW(delay_distribution_matrix),, drop = F] # trancate again, row starting from a different date
  
  Q_vector <- apply(truncated_delay_distribution_matrix, MARGIN = 2, sum) # summing each columns and return a vector
  
  while(chi_squared > threshold_chi_squared & count <= max_iterations) {
    
    #if (verbose) {
    #  cat("\t\tStep: ", count, " - Chi squared: ", chi_squared, "\n")
    #}
    
    E <- as.vector(delay_distribution_matrix %*% current_estimate) # matrix multiplication, n*n * n*1 = n*1
    B <- replace_na(original_incidence/E, 0) # ? wny missing values arise
    
    current_estimate <- current_estimate / Q_vector *  as.vector(crossprod(B, delay_distribution_matrix)) # ? why
    current_estimate <- replace_na(current_estimate, 0) #? why missing values arise
    
    chi_squared <- 1/N0 * sum((E[(max_delay + 1): length(E)] - original_incidence[(max_delay + 1) : length(original_incidence)])^2/E[(max_delay + 1): length(E)], na.rm = T) # what's this metric?
    count <- count + 1
  }
  
  return(current_estimate)
}

#########################################################################
#########################################################################
##### Nested functions 4 to construct delay distribution !
###  build usable delay distribution (parameter name constant_delay_distributions in deconvolve Incidence)
### get_vector_constant_waiting_time_distr. It build an empirical pdf from mix
get_vector_constant_waiting_time_distr <- function(shape_incubation,
                                                   scale_incubation,
                                                   shape_onset_to_report,
                                                   scale_onset_to_report,
                                                   length_out,
                                                   hypothesis,
                                                   days_further_in_the_past){ 
  F_h <- make_ecdf_from_gammas(shape = c(shape_incubation, shape_onset_to_report), 
                               scale = c(scale_incubation, scale_onset_to_report), 
                               hypothesis = hypothesis,
                               days_further_in_the_past = days_further_in_the_past)
  
  f <- Vectorize(function(x){
    if(x < 0) {
      return(0)
    } else if(x < 0.5) {
      return(F_h(0.5))
    } else {
      return(F_h(round(x + 1E-8) + 0.5) - F_h(round(x + 1E-8) - 0.5))
    }
  })
  
  x <- 0:(length_out - 1)
  
  return(f(x))
}


########################################## 
# Build empirical CDF from draws summing samples from two gamma distributions
make_ecdf_from_gammas <- function(shape, scale, 
                                  numberOfSamples = 1E6, 
                                  hypothesis,
                                  days_further_in_the_past) {
  if ( hypothesis == 'gamma') {
    draws <- rgamma(numberOfSamples, shape = shape[1], scale = scale[1]) + 
      rgamma(numberOfSamples, shape = shape[2], scale = scale[2])
    return(Vectorize(ecdf(draws)))
  } 
  else if ( hypothesis == 'beta') {
    draws <- rgamma(numberOfSamples, shape = shape[1], scale = scale[1]) + 
      rbeta(numberOfSamples, shape[2], scale[2]) * days_further_in_the_past
    return(Vectorize(ecdf(draws)))
  } 
  else if ( hypothesis == 'weibull') {
    draws <-
      rgamma(numberOfSamples, shape = shape[1], scale = scale[1]) + 
      rweibull(numberOfSamples, shape = shape[2], scale = scale[2])
    return(Vectorize(ecdf(draws)))
  } 
  #else {
  #  statement4}
}


###########################################################################33
#########################################################################33
### Main function to do deconvolution, get_infection_incidence_by_deconvolution.
###  Version, with boostrap
get_infection_incidence_by_deconvolution <- function(
  data_subset,# unsmoothed data with two columns, value and date
  constant_delay_distribution, # pdf for incubation + onset
  #constant_delay_distribution_incubation = c(), # pdf for incubation only
  #is_onset_data = F,# 不需要 and default
  #is_local_cases = T, # 不需要 and default
  #smooth_incidence = TRUE,# depends on the version of get_bootstrap_replicate
  days_incl,# hyperparameter
  #empirical_delays  = tibble(),# default
  n_bootstrap,# hyperparameter and 50
  days_further_in_the_past, # hyperparameter and default. this one is used in the data preprocessing as well as iterateRL
  #days_further_in_the_past_incubation = 5, # hyperparaeter and default. what's this for: its not in our case
  threshold_chi_squared,
  is_sampling, # TODO
  max_iterations = 100) {# default verbose = FALSE
  
  ####### Initialization
  data_subset <- data_subset %>%   # exclude leading zeroes
    arrange(date) %>%
    filter(cumsum(value) > 0)
  
  #data_type_subset <- unique(data_subset$data_type)[1] # 我们是n1或者n2, 只用在命名col里
  #data_type_name <- paste0("infection_", data_type_subset) # "infection_n1" or "infection_n2", 只用在命名col里
  
  minimal_date <- min(data_subset$date) - days_further_in_the_past
  maximal_date <- max(data_subset$date)
  all_dates <- seq(minimal_date, maximal_date, by = "days")
  
  ##### nested function that doesn't need to change, get_matrix_constant_waiting_time_distr
  delay_distribution_matrix <- get_matrix_constant_waiting_time_distr(constant_delay_distribution,all_dates)
  initial_delta <- min(which(cumsum(constant_delay_distribution) > 0.5)) - 1 # take median value (-1 because index 1 corresponds to zero days)
  
  ######### Output only one deconvolution       #no bootstrap for value? might need to add it. 
  ###### nested functions that we dont need to change, getLOESSCases
  #added bootstrap for value
  results <- list(tibble())
  for (bootstrap_replicate_i in 0:n_bootstrap) {
    
    if (bootstrap_replicate_i == 0) {
      time_series <- data_subset
    } else {
      time_series <- get_bootstrap_replicate(data_subset)
    }
    smoothed_incidence_data <- time_series %>%
      complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0)) %>% 
      mutate(value = getLOESSCases(dates = date, count_data = value, days_incl)) 
    
    ##### Scale the smoothed data, do we need this?
    raw_total_incidence <- sum(time_series$value, na.rm = TRUE)
    smoothed_total_incidence <- sum(smoothed_incidence_data$value, na.rm = T)
    
    if (smoothed_total_incidence > 0) {
      smoothed_incidence_data <- smoothed_incidence_data %>%
        mutate(value = value * raw_total_incidence / smoothed_total_incidence)
    }
    ###########################################
    # Deconvolution on smoothed data
    deconvolved_infections <-  do_deconvolution(smoothed_incidence_data,
                                                delay_distribution_matrix = delay_distribution_matrix,
                                                days_further_in_the_past = days_further_in_the_past,# 30
                                                initial_delta = initial_delta,# median value of gamma mixture
                                                max_iterations = max_iterations,
                                                threshold_chi_squared = threshold_chi_squared)  
    deconvolved_infections <- deconvolved_infections %>% slice((days_further_in_the_past -5 + 1):n())
    
    ## dataframe containing results
    deconvolved_infections <- tibble(
      date = deconvolved_infections$date,
      value = deconvolved_infections$value,
      replicate = bootstrap_replicate_i
    )
    results <- c(results, list(deconvolved_infections))
    
 } # end of function. tibble with two columns date and value.
  return(bind_rows(results))
}    





#try it with zurich alpha 


zurich_alpha<-data[which(data$city=='Zürich (ZH)'&data$variant=='B.1.1.7 (Alpha)'),]      

zurich_alpha_case<-zurich_alpha[c(2,8)]
colnames(zurich_alpha_case) <- c('date','value')

oneboot_data <- deconvolveIncidence(zurich_alpha_case, 
                                           #incidence_var = inc_var,
                                           getGammaParams(5.3, 3.2), # incubation
                                           getGammaParams(2.83, 2.96), # when its zero, it means consider only one delay
                                           n_bootstrap = 1,# by changing this from 0-50 we can get different 
                                           days_incl = 21, # smoothing parameter, get_infection_incidence_by_deconvolution - get_bootstrap_replicate - getLOESSCases
                                           days_further_in_the_past = 30, # deconvolution parameter specifying the maximum delay possible, get_infection_incidence_by_deconvolution 
                                           # also used to specify the range of beta distribution
                                           threshold_chi_squared = 1,
                                           hypothesis = "gamma" #
                                           #num_delays = 1, # number of delay distributions considered in the model
                                           #smooth_param = TRUE, 
)


#produce one original and one replicate. the dataset output is three column, one date, one value, one replicate numbner. 

# ggplot()+
#   geom_line(data=fiftyboot_data%>%filter(replicate<=10),aes(x=date,y=value,color=factor(replicate)))

zurich_alpha_ww<-zurich_alpha[c(2,4)]
colnames(zurich_alpha_ww) <- c('date','value')

#combine dataset
oneboot_data$datatype<-'confirmed cases'
zurich_alpha_ww$datatype<-'wastewater'
combine<-bind_rows(zurich_alpha_ww, oneboot_data)



set.seed(1234) #should I put it here? 

gammawholeloss<-function(df, par){
  
  waste<-df[which(df$datatype=='wastewater'),]      
  caseinfect<-df[which(df$datatype=='confirmed cases'),]
  
  if(par[1]<0|par[2]<0){   #the parameter in optim(lower=c(0,0) works the same)
    return(Inf) }
  ww_infect = deconvolveIncidence(waste, 
                                  #incidence_var = inc_var,
                                  IncubationParams=list(shape = 0, scale = 0),
                                  #IncubationParams = getGammaParams(5.3, 3.2),
                                  OnsetToCountParams = getGammaParams(par[1], par[2]), # incubation
                                  n_bootstrap = 1,# try with 1 bootstrap
                                  days_incl = 21, # smoothing parameter, get_infection_incidence_by_deconvolution - get_bootstrap_replicate - getLOESSCases
                                  days_further_in_the_past = 30, # deconvolution parameter specifying the maximum delay possible, get_infection_incidence_by_deconvolution 
                                  threshold_chi_squared = 1,
                                  hypothesis="gamma"
                                  #smooth_param = TRUE, 
  )
  result = compareTracesRMSE(ww_infect, caseinfect)   #call our rmse lose 
  return(result)
}

###### nested function in the loss function
compareTracesRMSE <- function(infect_i, infect_j){
  compare_df = infect_i %>%
    left_join(infect_j, by = 'date', suffix = c(".i", ".j")) %>%
    mutate(se = (value.i - value.j)^2)
  
  se = compare_df %>% pull(se)
  rmse = sqrt(sum(se, na.rm = T)/length(infect_i$date))
  
  return(rmse)
}

###########################################################################33
#########################################################################33
### gamma Optimization----get mean and sd from optimiztion 
wholeloss(combine,par)

par<-as.matrix(c(5.686912,4.067996))
#optim(par=par, fn=gammawholeloss, df=combine,method="SANN")    #run a long time. 





##########   beta. directly get shape and scale from optimization

betawholeloss<-function(df, par){
  
  waste<-df[which(df$datatype=='wastewater'),]      
  caseinfect<-df[which(df$datatype=='confirmed cases'),]
  
  if(par[1]<0|par[2]<0){   #the parameter in optim(lower=c(0,0) works the same)
    return(Inf) }
  ww_infect = deconvolveIncidence(waste, 
                                  #incidence_var = inc_var,
                                  IncubationParams=list(shape = 0, scale = 0),
                                  #IncubationParams = getGammaParams(5.3, 3.2),
                                  OnsetToCountParams = list(shape = par[1], scale = par[2]), # incubation
                                  n_bootstrap = 1,# try it with one bootstrap
                                  days_incl = 21, # smoothing parameter, get_infection_incidence_by_deconvolution - get_bootstrap_replicate - getLOESSCases
                                  days_further_in_the_past = 30, # deconvolution parameter specifying the maximum delay possible, get_infection_incidence_by_deconvolution 
                                  threshold_chi_squared = 1,
                                  hypothesis="beta",
                                  is_sampling = FALSE
                                  #smooth_param = TRUE, 
  )
  result = compareTracesRMSE(ww_infect, caseinfect)   #call our rmse lose 
  return(result)
}

par<-as.matrix(c(2,3))
optim(par=par, fn=betawholeloss, df=combine) #method="L-BFGS-B", lower=c(0, 0))







#grid search for one boot~ with mean and sd. what if we directly try with shape and scales? for now this one works. 

      meanOpts = seq(0.5, 10, 0.5)
      sdOpts = seq(0.5, 5, 0.5)

      deconv_results = cbind(expand_grid(meanOpts, sdOpts),
                             'rmse_cc' = NA)

      for (row_id in 1:nrow(deconv_results)){
        #set.seed=1234   should we set seed here? 
        deconv_config = try(deconvolveIncidence(zurich_alpha_ww, 
                                                IncubationParams=list(shape = 0, scale = 0),
                                                getGammaParams(deconv_results[row_id, 'meanOpts'],
                                                               deconv_results[row_id, 'sdOpts']), #can change for shape and scale. 
                                                hypothesis="gamma", n_boot = 1 ))

        if('try-error' %in% class(deconv_config)){
          deconv_results[row_id, c('rmse_cc')] = c(Inf)
          next
        }

        deconv_results[row_id, c('rmse_cc')] = compareTracesRMSE(deconv_config, oneboot_data#we get from decovole case data
)
      }

      write_csv(deconv_results, '/Users/meiyilong/Downloads/wastewaterRe-main/scan/deconv_try.csv')




#write an r markdown file that explains everything 





#f2 can directly get to the function 
#f1 can open help. 

