#Import Libraries
library(conflicted)
library(tidyverse)
library(ggalluvial)
pacman::p_load(rio, janitor)
library(lubridate)
library(tidyplots)
library(sf)
library(cartography)

setwd("~/Documents/Taf Stats/KlebGen_Presentation_Kampala")

#Import data
merged_data <- import("./data/Klebsiella pneumoniae sequencing tracking Africa CDC_Cameroon.xlsx", which = "LNSP-MergeData")
sociodemo_data <- import("./data/Klebsiella pneumoniae sequencing tracking Africa CDC_Cameroon.xlsx", which = "LNSP-sociodemo") %>% clean_names()
lab_geoloc <- import("./data/Klebsiella pneumoniae sequencing tracking Africa CDC_Cameroon.xlsx", which = 'Lab-GeoLoc')

#View data
View(merged_data)
view(sociodemo_data)
View(lab_geoloc)

#Data summary
skimr::skim(merged_data)
skimr::skim(sociodemo_data)

#Filter data with low completude 
#Clean columns names
complete_merged_data <- merged_data %>%
  dplyr::filter(completude>0.9) %>% clean_names()

skimr::skim(complete_merged_data)

#Change age in year
complete_merged_data <- complete_merged_data %>%
  mutate(age_in_year = case_match(unite_age_annee_a_mois_m_jour_j,
                                  "J" ~ age/365,
                                  "M" ~ age/12,
                                  .default = age))

#Select important columns 
complete_merged_data <- complete_merged_data %>%
  select(sequencing_sample_id, sample_collection_date_yyyy_mm_dd, strain, mlst, dna_extract_concentration_ng_u_l, starting_dna_amount_for_library_prep_ng, final_library_concentration_ng_u_l,
         age_in_year, gender, patient_in_out, geo_loc_name_admin1_e_g_region, geo_loc_name_admin1_e_g_location, specimen_type, origine_lab_organism_specif_ic)

#Rename columns 
names(complete_merged_data)
complete_merged_data <- complete_merged_data %>%
  rename(sample_id = sequencing_sample_id,
         collection_date = sample_collection_date_yyyy_mm_dd,
         dna_extract_concentration = dna_extract_concentration_ng_u_l,
         starting_dna = starting_dna_amount_for_library_prep_ng, 
         final_concentration = final_library_concentration_ng_u_l,
         region = geo_loc_name_admin1_e_g_region,
         laboratory_name  = geo_loc_name_admin1_e_g_location,
         original_strain = origine_lab_organism_specif_ic)

complete_merged_data <- complete_merged_data %>%
  rename(patient_type = patient_in_out)

complete_merged_data <- complete_merged_data %>%
  mutate(gender = case_match(gender,
                             "M" ~ "Male",
                             "F" ~ "Female",
                             .default = gender),
         patient_type = str_to_lower(patient_type),
         region = str_to_lower(region),
         laboratory_name = case_match(laboratory_name,
                                      "SW-06-RHL" ~ "RHL",
                                      "LT-HGOPED" ~ "HGOPED",
                                      "HRBAF" ~ "HR BAFOUSSAM",
                                      "HRBUEA" ~ "HR BUEA",
                                      "HRGAROUA"~ "HR GAROUA",
                                      .default = laboratory_name))
    
table(complete_merged_data$region)

complete_merged_data <- complete_merged_data %>%
  mutate(region = case_match(region,
                             "extreme-nord" ~ "extreme nord",
                             "south-west"~ "sud ouest",
                             .default = region),
         specimen_type = str_to_lower(specimen_type),
         original_strain = str_to_title(original_strain))

names(sociodemo_data)
sociodemo_data <- sociodemo_data %>%
  rename(sample_id = origine_sample_id,
         nphl_arrival_date = date_arrivee_lnsp, 
         collection_date = sample_collection_date_yyyy_mm_dd,
         region = geo_loc_name_admin1_e_g_region, 
         laboratory = geo_loc_name_admin1_e_g_location,
         strain_abrv = origine_lab_organism,
         strain = origine_lab_organism_specif_ic,
         patient_type = patient_in_out) %>%
  mutate(age_in_year = case_match(unite_age_annee_a_mois_m_jour_j,
                                  "J" ~ age/365,
                                  "M" ~ age/12,
                                  .default = age)) %>%
  select(sample_id, nphl_arrival_date, collection_date, region, laboratory, age_in_year, gender,patient_type,specimen_type,
         strain_abrv, strain) %>%
  mutate(gender = case_match(gender,
                             "M" ~ "Male",
                             "F" ~ "Female"),
         patient_type = str_to_lower(patient_type),
         region = str_to_lower(region),
         laboratory = case_match(laboratory,
                                      "SW-06-RHL" ~ "RHL",
                                      "LT-HGOPED" ~ "HGOPED",
                                      "HRBAF" ~ "HR BAFOUSSAM",
                                      "HRBUEA" ~ "HR BUEA",
                                      "HRGAROUA"~ "HR GAROUA",
                                      .default = laboratory),
         region = case_match(region,
                             "extreme-nord" ~ "extreme nord",
                             "south-west"~ "sud ouest",
                             .default = region),
         specimen_type = str_to_lower(specimen_type),
         strain = str_to_title(strain),
         strain_abrv = str_to_lower(strain_abrv)
         )



#Distribution of collected samples per month of collection and region
sociodemo_data %>%
  select(collection_date, region) %>%
  drop_na()%>%
  mutate(collection_date = as.Date(collection_date),
    collection_month = paste(year(collection_date), month(collection_date), sep = "-"),
    collection_month = case_match(collection_month,
                                  "2026-3" ~ "2025-3",
                                  "2823-8"~"2023-8",
                                  .default = collection_month)) %>%
  ggplot() +
  geom_bar(aes(x=collection_month, fill = region), color="black")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "collection month", y = "Total samples collected")

#Distribution of collected samples per month of arrival at NPHL
sociodemo_data %>%
  select(nphl_arrival_date, region) %>%
  drop_na()%>%
  mutate(nphl_arrival_date = as.Date(nphl_arrival_date),
         collection_month = paste(year(nphl_arrival_date), month(nphl_arrival_date), sep = "-")) %>%
  ggplot() +
  geom_bar(aes(x=collection_month), color="black", fill="grey")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "NPHL arrival month", y = "Total samples collected")

#Distribution of collected samples per laboratories 
sociodemo_data %>%
  ggplot() +
  geom_bar(aes(x=laboratory), color="black", fill="grey")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "Laboratory name", y = "Total samples collected")

#Distribution of collected samples per laboratories and strain
sociodemo_data %>%
  tidyplot(x=laboratory, color=strain_abrv) %>%
  add_barstack_relative() %>%
  adjust_x_axis(title = "Laboratory name", rotate_labels = 45) %>%
  adjust_y_axis_title(title = "Total samples collected") %>%
  adjust_legend_title("Strain") %>%
  adjust_size(width = 180, height = 180)

#Distribution of collected samples per region and strain
sociodemo_data %>%
  tidyplot(x=region, color=strain_abrv) %>%
  add_barstack_relative() %>%
  adjust_x_axis(title = "Region", rotate_labels = 45) %>%
  adjust_y_axis_title(title = "Total samples collected") %>%
  adjust_legend_title("Strain") %>%
  adjust_size(width = 180, height = 180)

#Distribution of collected samples per region and strain
sociodemo_data %>%
  select(nphl_arrival_date, region) %>%
  drop_na()%>%
  mutate(nphl_arrival_date = as.Date(nphl_arrival_date),
         collection_month = paste(year(nphl_arrival_date), month(nphl_arrival_date), sep = "-")) %>%
  tidyplot(x=collection_month, color=region) %>%
  add_barstack_relative() %>%
  adjust_x_axis(title = "NPHL arrival month", rotate_labels = 45) %>%
  adjust_y_axis_title(title = "Total samples collected") %>%
  adjust_legend_title("Region") %>%
  adjust_size(width = 180, height = 180)


sociodemo_data %>%
  select(strain_abrv) %>%
  drop_na() %>%
  tidyplot(color=strain_abrv) %>%
  add_pie() %>%
  adjust_legend_title("Region") %>%
  adjust_size(width = 150, height = 150)


#Mapping of region 
shp_path <- "./data/shape_2022/region_cam_2022.shp"
mtq <- st_read(shp_path)
plot(st_geometry(mtq))

df_region <- sociodemo_data %>%
  group_by(region) %>%
  summarise(total_sample = n()) %>%
  mutate(region = str_to_title(region))

code_rs <- c("CEN", "EXT", "LIT", "NOR", "OUE", "SUW")
df_region$code_rs <- code_rs

mtq <- mtq %>% 
  left_join(df_region, by = c("Code_RS" = "code_rs"))


choroLayer(x = mtq, var = "total_sample",
           method = "quantile", nclass = 4,
           border = "grey40",
           legend.pos = "topright", legend.values.rnd = 1,
           legend.title.txt = "Total samples collected")
labelLayer(x = mtq, txt = "Code_RS",
           halo = TRUE, overlap = FALSE)


#Sequencing result 
complete_merged_data %>% 
  ggplot() +
  geom_bar(aes(y = strain)) +
  labs(x = "Total sequences", y = "Strain")


complete_merged_data %>% 
  ggplot() +
  geom_bar(aes(y = strain, fill = gender)) +
  labs(x = "Total sequences", y = "Strain", fill="Gender")

complete_merged_data %>% 
  ggplot() +
  geom_boxplot(aes(y = age_in_year, x=strain)) +
  labs(y = "Age (year)", x = "Strain")

complete_merged_data %>% 
  dplyr::filter(strain %in% c("Klebsiella Pneumoniae", "Klebsiella quasipneumoniae", "Klebsiella Variicola")) %>%
  ggplot() +
  geom_bar(aes(x=mlst))+
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  labs(x = "MLST", y = "Total sequences", title = "KPSC") 
  
table(complete_merged_data$strain, complete_merged_data$region)

