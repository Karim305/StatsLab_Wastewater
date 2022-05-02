
#only use rmse for the lose

compareTracesRSE <- function(Re_i, Re_j){
  compare_df = Re_i %>%
    left_join(Re_j, by = 'date', suffix = c(".i", ".j")) %>%
    mutate(se = (median_R_mean.i - median_R_mean.j)^2)
  
  se = compare_df %>% pull(se)
  rmse = sqrt(sum(se, na.rm = T)/length(Re_i$date))
  
  return(rmse)
}


#combine the dataset. first part is teh wastewater RNA, second part is the confirmedcase Re we get from previous coding. because they have different column names I just created a bunch of NAs..


wwtry<-ww_data[c('date','flow','region','norm_n1','norm_n2')]
wwtry$country<-NA
wwtry$datatype<-NA
wwtry$estimate_type<-NA
wwtry$median_R_mean<-NA
wwtry$median_R_highHPD<-NA
wwtry$median_R_lowHPD<-NA
wwtry$countrylso3<-NA
wwtry$data_type<-'wastewater'

caseRetry<-plotData #plotData is case confirmed Re
caseRetry$flow<-NA
caseRetry$norm_n1<-NA
caseRetry$norm_n2<-NA

combine<-bind_rows(wwtry, caseRetry)

#this combine is the dataset we will use in optim()


#this is our lose function 

wholeloss<-function(df, par){
  
  waste<-df[which(df$data_type=='wastewater'),]      
  caseRe<-df[which(df$data_type=='Confirmed cases'),]
  
  #deconvolution for wastewater
  config_df = expand.grid("region" = c('ZH'),  
                          'incidence_var' = c('norm_n1'),
                          'FirstGamma' = 'incubation'    
  )
  
  
  deconv_ww_data <- data.frame()
  Re_ww <- data.frame()
  
  for(row_i in 1:nrow(config_df)){
    new_deconv_data = deconvolveIncidence(waste, #I deleted zurich filter here. no need
                                          incidence_var = config_df[row_i, 'incidence_var'],
                                          getCountParams(as.character(config_df[row_i, 'FirstGamma'])), 
                                          getGammaParams(mean=par[1], sd=par[2]), 
                                          smooth_param = TRUE, n_boot = 50) 
    
    new_deconv_data <- new_deconv_data %>%
      mutate(incidence_var = config_df[row_i, 'incidence_var'])
    
    ##### Get Re #####
    new_Re_ww = getReBootstrap(new_deconv_data)
    new_Re_ww <- new_Re_ww %>%
      mutate(variable = config_df[row_i, 'incidence_var'],
             region = config_df[row_i, 'region'])
    
    deconv_ww_data <- bind_rows(deconv_ww_data, new_deconv_data)
    Re_ww = bind_rows(Re_ww, new_Re_ww)
  }
  

  #combine wastewater Re and case Re dataset
  
  all_Re <- as_tibble(Re_ww) %>%
    select(region, data_type, date, median_R_mean, 
           median_R_highHPD, median_R_lowHPD, source) %>%
    mutate(data_type = recode(data_type,
                              'infection_norm_n1' = 'N1')) %>%
    bind_rows(caseRe) #deleted zurich filter
  
  ## Compare Rww and Rcc in all combinations ####
  data_types = unique(all_Re$data_type)
  rmse_matrix = matrix(data = NA, nrow = length(data_types), ncol = length(data_types),
                       dimnames = list(data_types, data_types))
  coverage_matrix = matrix(data = NA, nrow = length(data_types), ncol = length(data_types),
                           dimnames = list(data_types, data_types))
  mape_matrix = matrix(data = NA, nrow = length(data_types), ncol = length(data_types),
                       dimnames = list(data_types, data_types))
  
  for(ind_i in seq_along(data_types)){
    for (ind_j in seq_along(data_types)){
      if (ind_i == ind_j){next}
      
      subRe_i = all_Re %>% 
        dplyr::filter(data_type == data_types[ind_i])
      subRe_j = all_Re %>% 
        dplyr::filter(data_type == data_types[ind_j])
      
      result = compareTracesRSE(subRe_i, subRe_j)   #call our rmse lose 
      
      rmse_matrix[ind_i, ind_j] = result
    }
  }
  
  return(result)
  
}


install.packages("optimParallel")
install.packages("parallel")
library("optimParallel")
library(dbplyr)

# initialize par

par<-as.matrix(c(7.5,0.5))

#do the optimization 
optim(par=c(1,2), fn=wholeloss, df=combine, method = "L-BFGS-B")

# Error in ecdf(draws) : 'x' must have 1 or more non-missing values
# In addition: Warning message:
#   In rgamma(numberOfSamples, shape = shape[2], scale = scale[2]) :
#   
#   Error in ecdf(draws) : 'x' must have 1 or more non-missing values 
#even if it is 0,0 for shape2 and scale2, draws is not 0 because it plus the incubation....how chould ecdf draws don't have value?


#deal with iteration. I think the bug happened when updating parameters. the amount of iteration needed for the bug changes each time. it must be the parameter update problem. what kind of parameter will make this bug?



#piss me off! how could the value be 0.09263511? maybe if we change to infection, and change the otpimization algorism, it can be better. if we do it since incubation, if might also be better




