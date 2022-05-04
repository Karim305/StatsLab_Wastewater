## StatsLab Checking Jana's Code

library(knitr)
library(tidyverse)
library(lubridate)
library(patchwork)
library(viridis)
library(EpiEstim)
library(magrittr)

dir <-"/Users/Karim/Library/CloudStorage/OneDrive-Persönlich/ETH/T2/Stats_Lab/JanaRepo"

app_location = paste0(dir,'/covid-19-re-shiny-app')
source(paste0(app_location,'/app/otherScripts/2_utils_getInfectionIncidence.R'))
source(paste0(app_location,'/app/otherScripts/3_utils_doReEstimates.R'))
source(paste0(dir,'/wastewaterRe/code/wastewater_functions.R'))

plot_dir = paste0(dir,'/figures2')

theme_set(theme_minimal() +
            theme(
              strip.text = element_text(size=20),
              axis.text= element_text(size=17),
              axis.title =  element_text(size=20),
              legend.text= element_text(size=17),
              legend.title= element_text(size=20)
            ))


ZH_flow_url = "http://parsivel-eawag.ch/sarscov2/__data__/ARA%20Werdhoelzli_flow_cases.csv"
ZH_genes_url = "http://parsivel-eawag.ch/sarscov2/__data__/ARA%20Werdhoelzli_genes.csv"

raw_flow_data_ZH <- read_delim(ZH_flow_url, delim = ';',
                               col_names = c('date', 'cases', 'cases_smooth', 
                                             'flow', 'n1_smooth', 'n2_smooth'),
                               col_types = cols(date = col_date(format = '')),
                               skip = 1) 

raw_gene_data_ZH <- read_delim(ZH_genes_url, delim = ';',
                               col_names = c('date', 'n1', 'n2'),
                               col_types = cols(date = col_date(format = '')),
                               skip = 1) 


"""
# Summary of code used in WastewaterRe Paper by Huisman et. al


# - wastewater_Zurich.R

## - Initial data manipulation (outside of complicated, nested functions)

### - Wastewater ("Flow") data

- flow data is loaded with col_names = c('date', 'cases', 'cases_smooth','flow', 'n1_smooth', 'n2_smooth')

- gene data is loaded with col_names = c('date', 'n1', 'n2'),

- missing values in flow data are imputed by linear interpolation:
"""  

raw_data_ZH <- raw_flow_data_ZH %>%
  left_join(raw_gene_data_ZH, c('date')) %>%
  filter(!is.na(n1),
         date >= as_date("2020-09-01"),
         date <= as_date("2021-01-20"),
         date != as_date("2020-10-29")) %>%
  mutate(orig_data = TRUE) %>%
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) %>%
  mutate(region = 'ZH')

"""
### - Case data

- Catchment level clinical case data: missing values are set to 0

- Intuition: a positive case is most likely reported for the next day if the actual day was missing 
"""

orig_cases <- read_csv(paste0(dir,'/wastewaterRe/data/ZH_case_incidence_data.csv')) %>%
  filter(date >= as_date("2020-08-15"),
         date <= as_date("2021-01-20")) %>%
  mutate(across(c(-date), ~ ifelse(is.na(.), 0, .)))  %>%
  select(date, confirmed)



plot_orig_cases <- orig_cases %>% 
  pivot_longer(cols = c(-date))


case_data_plot <- ggplot(plot_orig_cases) +
  geom_bar(aes(x=date, y= value, fill = name), alpha = 0.5,
           position = 'identity', stat = 'identity', show.legend = F) +
  labs(x = 'Date' , y='Cases per day') +
  scale_x_date(limits = c(as_date('2020-08-15'), as_date('2021-01-20')) ) +
  scale_fill_manual(values = c(viridis(4)[1:2])) + 
  labs(colour = 'Variable')

case_data_plot

"""
## - Deconvolution of Case data 

- Let $C$ denote the observed cases and $D$ denote the delay distribution ("transfer function") from true underlying incidence $I_{cc}$ to $C$. 
- Then, $C = I_{cc} \ast D$ is the convolution of true incidence $I_{cc}$ and delay function $D$.
- Thus, in order to obtain an estimate for the underlying true incidence, which we want to use to obtain our Shedding Load Distribution (SLD), we need to do a deconvolution of observed cases and delay distribution.
"""

deconv_cases <- data.frame()
Re_cases <- data.frame()
for(inc_var in c('confirmed')){
  new_deconv_data = deconvolveIncidence(orig_cases, 
                                        incidence_var = inc_var,
                                        getCountParams('incubation'), 
                                        getCountParams(paste0(inc_var, '_zh')),
                                        smooth_param = TRUE, n_boot = 50)
  new_Re = getReBootstrap(new_deconv_data) 
  
  deconv_cases <- bind_rows(deconv_cases, new_deconv_data)
  Re_cases = bind_rows(Re_cases, new_Re)
}

"""

# - deconvolveIncidence()

- Let's have a closer look at the called function wrapper:
"""
inc_var <- "confirmed"
IncubationParams <- getCountParams('incubation') 
OnsetToCountParams <- getCountParams(paste0(inc_var, '_zh'))

constant_delay_distributions <- list("Simulated" = get_vector_constant_waiting_time_distr(
    IncubationParams$shape, IncubationParams$scale,
    OnsetToCountParams$shape, OnsetToCountParams$scale),
    "Symptoms" = get_vector_constant_waiting_time_distr(
      IncubationParams$shape, IncubationParams$scale,
      0, 0))
  
constant_delay_distributions

sum(constant_delay_distributions[[1]])
sum(constant_delay_distributions[[2]])


par(mfrow=c(2,1))
plot(1:30,constant_delay_distributions[[1]][1:30], type = "l")
plot(1:30,constant_delay_distributions[[2]][1:30], type = "l")


ggplot() +
  geom_line(x = 1:30, y= constant_delay_distributions[[1]][1:30], color = 'Incubation')

delay_distribution_inc <- tibble(t = 1:30,
                              y = constant_delay_distributions[[1]][1:30],
                              type = 'Incubation')


delay_distribution_sym <- tibble(t = 1:30,
                                 y = constant_delay_distributions[[2]][1:30],
                                 type = 'Symptom Onset')

c(list(delay_distribution_inc),list(delay_distribution_sym)) %>% 
  dplyr::bind_rows()

c(delay_distribution_inc,delay_distribution_sym) %>% 
  rowbind()