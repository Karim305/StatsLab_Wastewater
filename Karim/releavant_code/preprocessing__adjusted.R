### Data preprocessing
# previously, we kept only days where the selected variants are dominant for cases as well.
# change that
### load packages
library(tidyverse)
library(ggplot2)
library(stats)
library(zoo)
library(imputeTS)

library(devtools)
install_github("covid-19-Re/estimateR")

dir <-"/Users/Karim/Library/CloudStorage/OneDrive-Persönlich/ETH/T2/Stats_Lab/JanaRepo"

app_location = paste0(dir,'/covid-19-re-shiny-app')
source(paste0(app_location,'/app/otherScripts/2_utils_getInfectionIncidence.R'))
source(paste0(app_location,'/app/otherScripts/3_utils_doReEstimates.R'))

source(paste0(dir,'/wastewaterRe/code/wastewater_functions.R'))

library(estimateR)

## import data ##
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

# cities with enough data
cities2 <- c("Altenrhein (SG)", "Chur (GR)", "Laupen (BE)", "Lugano (TI)", "Zürich (ZH)")
variants <- c("B.1.1.7 (Alpha)", "B.1.617.2 (Delta)")

data %>% 
  distinct(location)


## Preprocessing for new_cases
df_cases <- data %>% 
  mutate(city = case_when(
    
    location == "altenrhein" ~ "Altenrhein (SG)",
    location == "chur" ~ "Chur (GR)",
    location == "laupen" ~ "Laupen (BE)",
    location == "lugano" ~ "Lugano (TI)",
    location == "zurich" ~ "Zürich (ZH)"
    
  )) %>% 
  filter(!is.na(city)) %>% 
  distinct(city, date, new_cases) %>% 
  group_by(city) %>% 
  arrange(date) %>% 
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>% 
  mutate(new_cases = ifelse(is.na(new_cases),0,new_cases))

# save case data
write_csv(df_cases, "data/processed/cases.csv")


## preprocessing for wastewater data  
df_flow <- data %>% 
  mutate(city = case_when(
    
    location == "altenrhein" ~ "Altenrhein (SG)",
    location == "chur" ~ "Chur (GR)",
    location == "laupen" ~ "Laupen (BE)",
    location == "lugano" ~ "Lugano (TI)",
    location == "zurich" ~ "Zürich (ZH)"
    
  )) %>% 
  select(city, date, rna = sars_cov2_rna, flow ,rna_v) %>% 
  filter(!is.na(city), !is.na(rna)) 

df_flow %>% 
  group_by(city, date) %>% 
  mutate(nobs = n_distinct(rna)) %>% 
  arrange(desc(nobs))

df_v1 <- df_flow %>% 
  filter(rna_v == "v1") %>% 
  group_by(city) %>% 
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>% 
  mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) %>% 
  mutate(rna_v = "v1")

df_v2 <- df_flow %>% 
  filter(rna_v == "v2") %>% 
  group_by(city) %>% 
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>% 
  mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) %>% 
  mutate(rna_v = "v2")

# variant proportions
# import data
df_ww <- read_csv(file = "data/variant_proportions_wastewater.csv") %>% 
  select(- ...1) %>% 
  rename(date = index)

# reshape data from wide to long for better plotting
df_ww_long <- df_ww %>% 
  gather(variant, proportion, -c(date, city)) %>% 
  mutate(variant = ifelse(variant == "B.1.1.7", "B.1.1.7 (Alpha)",
                          ifelse(variant == "B.1.617.2", "B.1.617.2 (Delta)", variant)))

# plot data
df_ww_long %>% 
  ggplot() +
  geom_line(aes(x = date, y = proportion, color = variant)) +
  facet_wrap(~ city, ncol = 3)

city_vec <- df_ww_long %>% 
  distinct(city) %>% 
  pull()

variant_vec <- df_ww_long %>% 
  distinct(variant) %>% 
  pull()



cities2 %in% city_vec
variants %in% variant_vec



data_prop <- subset(df_ww_long, (city %in% cities2) & (variant %in% variants))



# threshold 
thresh <- 0.90
min_days <- 20

df_props <-data_prop %>% 
  filter(proportion >= thresh) %>% 
  group_by(city, variant) %>% 
  arrange(date)

df_select <- df_props %>% 
  group_by(city, variant) %>% 
  summarize(n_days = n_distinct(date, na.rm = TRUE),
            first_day = min(date, na.rm = TRUE),
            last_day = max(date, na.rm = TRUE)) 

df_int <- df_ww_long %>% 
  inner_join(df_select, by = c("city", "variant")) %>% 
  filter(date >= first_day, date <= last_day) %>% 
  filter(n_days >= min_days) %>% 
  group_by(city, variant) %>% 
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>% 
  #mutate(rna = ifelse(is.na(rna),0,rna)) %>% 
  mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) %>%
  select(-n_days) #%>% 
#mutate(cv = paste0(city,'_',variant))


df_int %>% 
  ggplot() +
  geom_line(aes(x = date, y = proportion, color = variant)) +
  facet_wrap(~ city, ncol = 3)

df_out <- df_v1 %>% 
  inner_join(df_int, by = c("city", "date")) %>% 
  group_by(city, variant) %>% 
  mutate(rna_norm = rna / min(rna, na.rm = T)) %>% 
  ungroup() 

df_out %>% 
  ggplot() +
  geom_line(aes(x=date,y=rna_norm,colour = variant)) +
  facet_wrap(~city, ncol = 2,scales = "free_y")


df_out %>% 
  filter(rna_v == "v1") %>% 
  select(-rna_v) %>% 
  write_csv("data/processed/ww.csv")



