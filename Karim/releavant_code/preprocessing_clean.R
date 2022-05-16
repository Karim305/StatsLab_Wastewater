### Data Preprocessing


# setup -------------------------------------------------------------------

### load packages
library(tidyverse)
library(ggplot2)

### figure directory
fig_dir <- "/Users/Karim/Library/CloudStorage/OneDrive-Persönlich/ETH/T2/Stats_Lab/StatsLab_Wastewater/Karim/figures/"

### load data
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
cities2 <- c("Altenrhein (SG)", "Chur (GR)", "Laupen (BE)", "Lugano (TI)", "Zürich (ZH)", "Genève (GE)")
variants <- c("B.1.1.7 (Alpha)", "B.1.617.2 (Delta)", "BA.1 / BA.2 (Omicron)")


# case data ---------------------------------------------------------------

df_cases <- data %>% 
  mutate(city = case_when(
    
    location == "altenrhein" ~ "Altenrhein (SG)",
    location == "chur" ~ "Chur (GR)",
    location == "laupen" ~ "Laupen (BE)",
    location == "lugano" ~ "Lugano (TI)",
    location == "zurich" ~ "Zürich (ZH)",
    location == "geneve" ~ "Genève (GE)"
    
  )) %>% 
  # filter(!is.na(city)) %>% 
  distinct(city, date, new_cases) %>% 
  mutate(new_cases = round(new_cases * 100000)) %>% 
  group_by(city) %>% 
  arrange(date) %>% 
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>% 
  mutate(new_cases = ifelse(is.na(new_cases),0,new_cases)) 

df_cases %>% 
  write_csv("data/processed/cases.csv")  



# wastewater data ---------------------------------------------------------

## preprocessing for wastewater data  
df_flow <- data %>% 
  mutate(city = case_when(
    
    location == "altenrhein" ~ "Altenrhein (SG)",
    location == "chur" ~ "Chur (GR)",
    location == "laupen" ~ "Laupen (BE)",
    location == "lugano" ~ "Lugano (TI)",
    location == "zurich" ~ "Zürich (ZH)",
    location == "geneve" ~ "Genève (GE)"
    
  )) %>% 
  select(city, date, rna = sars_cov2_rna, flow ,rna_v) %>% 
  # filter(!is.na(city), !is.na(rna)) 
  filter(!is.na(rna)) 

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

## variant proportions

# import data
df_ww <- read_csv(file = "data/variant_proportions_wastewater.csv") %>% 
  select(- ...1) %>% 
  rename(date = index)


# reshape data from wide to long for better plotting

df_ww_long <- df_ww %>% 
  gather(variant, proportion, -c(date, city)) %>% 
  mutate(variant = case_when(variant == "B.1.1.7" ~ "B.1.1.7 (Alpha)",
                             variant == "B.1.617.2" ~ "B.1.617.2 (Delta)",
# create new omikron category that includes both subvariants
                             (variant == "BA.1" | variant == "BA.2") ~ "BA.1 / BA.2 (Omicron)",
                             TRUE ~ variant)) %>% 
  group_by(city, date, variant) %>% 
# add up proportions for omicron subvariants  
  summarize(proportion = sum(proportion), .groups = "drop")


# select relevant periods -------------------------------------------------

data_prop <- subset(df_ww_long, (city %in% cities2) & (variant %in% variants))

# threshold 
thresh <- 0.90
min_days <- 20

df_select <- data_prop %>% 
  filter(proportion >= thresh) %>% 
  group_by(city, variant) %>% 
  summarize(n_days = n_distinct(date, na.rm = TRUE),
            first_day = min(date, na.rm = TRUE),
            last_day = max(date, na.rm = TRUE),
            .groups = "drop") 

# dataframe with relevant periods only
df_int <- df_ww_long %>% 
  inner_join(df_select, by = c("city", "variant")) %>% 
  filter(date >= first_day & date <= last_day & n_days >= min_days) %>% 
  group_by(city, variant) %>% 
  complete(date = seq.Date(min(date), max(date), by = 'days')) %>% 
  mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = F) )) 


  
dominant_variants <- subset(df_ww_long, (city %in% cities2)) %>% 
  group_by(city,date) %>% 
  arrange(desc(proportion)) %>%
  filter(row_number() == 1) 

# problem: we consider a variant dominant from the first date to the last date 
# it has a proportion above 0.9.
# nonetheless, we want to keep proportions for the other periods as well, so we
# set proportion and variant to the proportion and the variant that dominates 
# over the used periods even if undetermined is dominant at some date withing that
# period. For the periods we don't include in our analysis, proportion and variant
# correspond to the proportion and name of the variant with the highest proportion
# daywise.
df_variants_out <- dominant_variants %>% 
  left_join(df_int, by = c("city", "date"),
            suffix = c("_d", "_i")) %>% 
  mutate(variant = variant_i,
         proportion = proportion_i,
         include = !is.na(variant_i)) %>% 
  mutate(variant = ifelse(is.na(variant),variant_d,variant),
         proportion = ifelse(is.na(proportion),proportion_d, proportion)) %>% 
  ungroup() %>% 
  select(city, date, variant, proportion, include, first_day, last_day)



# plotting proportions ----------------------------------------------------

# selected periods
df_int %>% 
  ggplot() +
  geom_line(aes(x = date, y = proportion, color = variant)) +
  facet_wrap(~ city, ncol = 2) + 
  ggtitle("Variant proportions for selected variants and periods")

ggsave(paste0(fig_dir,"variant_proportions_selected.png"))

# all periods and variants
subset(df_ww_long, (city %in% cities2)) %>% 
  group_by(city, variant) %>% 
  ggplot() +
  # highlight
  geom_rect(aes(xmin=first_day, xmax = last_day, fill = variant), ymin = -0.5, ymax = 1,
            alpha = 0.5, data = df_variants_out %>% filter(include)) +
  
  geom_area(aes(x = date, y = proportion, fill = variant), alpha = 0.8) +

  geom_text(
    aes(x = first_day, y = -0.1, label = "used"), color = "white",
    data = df_variants_out %>% filter(include),
    size = 3, vjust = 0, hjust = 0, nudge_x = 2, check_overlap = F
  ) +

    facet_wrap(~city) +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle("Variant Proportions (boxes denote periods used in our analysis)")

ggsave(paste0(fig_dir,"variant_proportions_all.png"))




# combining proportions and wastewater ------------------------------------

## we treat v1 and v2 measures as distinct data 

df_int_v1 <- df_v1 %>%
  left_join(df_variants_out, by = c("city", "date")) %>% 
  mutate(include = ifelse(is.na(include), FALSE, include))

df_int_v2 <- df_v2 %>% 
  left_join(df_variants_out, by = c("city", "date")) %>% 
  mutate(include = ifelse(is.na(include), FALSE, include))

df_int_v1 %>% 
  select(city, date, rna, flow, rna_v, variant, proportion, include) %>% 
  write_csv("data/processed/ww_v1.csv")

df_int_v2 %>% 
  select(city, date, rna, flow, rna_v, variant, proportion, include) %>% 
  write_csv("data/processed/ww_v2.csv")  


# plotting wastewater -----------------------------------------------------

df_version_plot <- df_int_v1 %>% 
  bind_rows(df_int_v2) %>% 
  mutate(type = ifelse(
    include,
    paste0("Version: ", rna_v,", Variant: ", variant), 
    "excluded") 
    )

df_version_plot %>% 
  ggplot() + 
  geom_line(aes(x = date, y = rna, color = rna_v)) +
  facet_wrap(~city, ncol = 2) + 
  scale_y_log10() +
  ggtitle("All wastewater data by datatype (log-scale)") 

ggsave(paste0(fig_dir,"wastewater_by_datatype_log.png"))

df_version_plot %>% 
  ggplot() + 
  geom_line(aes(x = date, y = rna, color = rna_v)) +
  facet_wrap(~city, ncol = 2) + 
  # scale_y_log10() +
  ggtitle("All wastewater data by datatype") 

ggsave(paste0(fig_dir,"wastewater_by_datatype.png"))



df_version_plot %>%   
  filter(include) %>% 
  ggplot() + 
  geom_line(aes(x = date, y = rna, color = type)) +
  facet_wrap(~city, ncol = 2) + 
  scale_y_log10() +
  ggtitle("Used wastewater data (log-scale)") 

ggsave(paste0(fig_dir,"wastewater_used_log.png"))

df_version_plot %>%   
  filter(include) %>% 
  ggplot() + 
  geom_line(aes(x = date, y = rna, color = type)) +
  facet_wrap(~city, ncol = 2) + 
  ggtitle("Used wastewater data (log-scale)") 

ggsave(paste0(fig_dir,"wastewater_used.png"))



