#########################################
#                                       #
#       GMA Project - Data Process      #
#                                       #
#########################################

# Sys.setlocale("LC_ALL", "Portuguese_Brazil.1252")

### Choose which data frame to save and work on the bottom line

#packages
library(compositions)
library(stats)
library(sp)
library(gstat)
library(tidyverse)

## Data loading
setwd("~/Desktop/Doutorado UFV/FIT 799 - Pesquisa/Scripts/GMA Scripts")
df_irr <- read.csv("Analise_de_Geracao-TomateAnao-REV00-DENIZ-PC - Uniformidade.csv",
                   header = T)
df_raw <- read.csv("Analise_de_Geracao-TomateAnao-REV48-DENIZ-PC - ColetaDados.csv",
                   header = T)

##  Raw data filtering
df_irr <- data.frame(df_irr[1:248,1:3])
df_raw <- as_tibble(df_raw[1:600,1:198])
df_clean_vars <- df_raw %>% 
  select_if(~ !all(is.na(.))) %>% 
  arrange(Ind)

##  Setting index data.frame
df_names <- cbind(df_clean_vars[,6:10], rep(seq(1,31, length.out = 75),8))
names(df_names) <- c("ind", "generation", "gen_color", "line", "plant_num", "plant_pos")

df_names <- df_names %>%
  as_tibble() %>% 
  mutate(across(c(ind:plant_num), as.factor)) %>% 
  mutate(generation = factor(df_names$gen_color,
                             # levels = c("VERDE", "ROSA", "ROXO", "AMARELO", "VERMELHO", "AZUL"),
                             # labels = c("F2", "BC1", "BC2", "P1", "F1", "P2")))
                             levels = c("AMARELO", "AZUL", "VERMELHO", "VERDE", "ROSA", "ROXO"),
                             labels = c("P1", "P2", "F1", "F2", "BC1", "BC2")))

df_names <- df_names %>%
  mutate(dam = factor(df_names$generation,
                      levels = c("F2", "BC1", "BC2", "P1", "F1", "P2"),
                      labels = c("F1F0", "F1F0", "P2F9", "P1F9", "P2F9", "P2F9"))) %>%
  mutate(sire = factor(df_names$generation,
                       levels = c("F2", "BC1", "BC2", "P1", "F1", "P2"),
                       labels = c("F1F0", "P1F9", "F1F0", "P1F9", "P1F9", "P2F9"))) %>%
  mutate(geno = factor(if_else(df_names$generation %in% c("F1", "P1", "P2"),
                               paste0(df_names$generation),
                               paste0(df_names$generation, "-", df_names$ind)))) %>% 
  relocate(c(generation, geno, dam, sire), .after = gen_color)

df_names <- df_names %>% 
  mutate(rep = ave(seq_along(df_names$geno), df_names$geno, FUN = seq_along)) %>% 
  relocate(c(rep), .after = sire)

##  Irrigation data interpolation
# Setting grid
grid <- df_names %>%
  select(ind, line, plant_pos) %>% 
  mutate(line = as.numeric(line))
coordinates(grid) <- ~ line + plant_pos
gridded(grid) = TRUE

# Interpolating
irrig <- df_irr %>% 
  mutate(Valor = (0.12*(as.numeric(Valor))+1))
coordinates(irrig) <- ~ Linha + Medida
irrig.idw <- idw(Valor ~ 1, irrig, grid, idp=1, nmax=2)

spplot(irrig.idw, zcol = "var1.pred")

# Binding irrigation covariate to data.frame
df_names <- cbind(df_names,
                  irr_flw = irrig.idw$var1.pred)

##  Phenotypic data processing
# Main phenotypic variables df
df_phenoMain <- df_clean_vars[,c(16,25,30,35,36,85,89,91:107)]
df_phenoMain <- df_phenoMain %>%
  as_tibble() %>%
  mutate(across(where(is.character), ~ as.numeric(gsub("N/A", "", .)))) %>%
  mutate(across(c(1:length(df_phenoMain)), as.numeric))
names(df_phenoMain) <- c("num_sideshoot",
                         "total_leaf_sideshoot",
                         "total_cluster_sideshoot",
                         "sideshoot_length",
                         "hypocotyl_diam",
                         "num_clusters_to_det",
                         "num_leaves_to_firstcluster",
                         "num_leaves_from_firstcluster_to_det",
                         "num_leaflets_per_leaf",
                         "petiole_length_cm",
                         "leaf_length_cm",
                         "main_leaflet_length_cm",
                         "main_leaflet_width_cm",
                         "side_leaflet_length_cm",
                         "side_leaflet_width_cm",
                         "petiolule_length_cm",
                         "rachilla_length_cm",
                         "secondary_leaflet_length",
                         "secondary_leaflet_blade_length",
                         "secondary_leaflet_width",
                         "height_to_det",
                         "three_internode_length",
                         "leaf_angle",
                         "leaflet_angle")
df_phenoMain <- cbind(df_names, df_phenoMain) %>% 
  arrange(line, plant_num)

#All phenotypic variables df
df_pheno <- df_clean_vars[,c(16:36,84:107)]
df_pheno <- df_pheno %>%
  as_tibble() %>%
  mutate(Acamamento = case_when(
    Acamamento == "S" ~ 1,
    Acamamento == "N" ~ 0)) %>% 
  mutate(Indeterminada = case_when(
    Indeterminada == "S" ~ 1,
    Indeterminada == "N" ~ 0)) %>% 
  mutate(across(where(is.character), ~ as.numeric(na_if(., "N/A")))) %>%
  mutate(across(c(1:length(df_pheno)), as.numeric))
names(df_pheno) <- c("num_sideshoot",
                     "desuckering_29_02_2024",
                     "desuckering_13_03_2024",
                     "desuckering_22_03_2024",
                     "desuckering_03_04_2024",
                     "desuckering_13_04_2024",
                     "desuckering_21_04_2024",
                     "desuckering_03_05_2024",
                     "desuckering_25_06_2024",
                     "total_leaf_sideshoot",
                     "leaf_desuckering_25_06_2024",
                     "leaf_desuckering_1",
                     "leaf_desuckering_2",
                     "leaf_desuckering_3",
                     "total_cluster_sideshoot",
                     "cluster_desuckering_25_06_2024",
                     "cluster_desuckering_1",
                     "cluster_desuckering_2",
                     "cluster_desuckering_3",
                     "sideshoot_length",
                     "hypocotyl_diam",
                     "lodging",
                     "num_clusters_to_det",
                     "leaf_color_L",
                     "leaf_color_a",
                     "leaf_color_b",
                     "num_leaves_to_firstcluster",
                     "growth_habit",                       #sim = indeterminate
                     "num_leaves_from_firstcluster_to_det",
                     "num_leaflets",
                     "petiole_length_cm",
                     "leaf_length_cm",
                     "main_leaflet_length_cm",
                     "main_leaflet_width_cm",
                     "side_leaflet_length_cm",
                     "side_leaflet_width_cm",
                     "petiolule_length_cm",
                     "rachilla_length_cm",
                     "secondary_leaflet_length",
                     "secondary_leaflet_blade_length",
                     "secondary_leaflet_width",
                     "height_to_det",
                     "three_internode_length",
                     "leaf_angle",
                     "leaflet_angle")

df_pheno <- cbind(df_names, df_pheno) %>% 
  arrange(line, plant_num) %>%
  mutate(across(total_leaf_sideshoot, ~ replace_na(., 0))) %>% 
  mutate(across(total_cluster_sideshoot, ~ replace_na(., 0))) %>% 
  mutate(internode_length = three_internode_length/3) %>% 
  mutate(height_to_det = height_to_det/100) %>% 
  mutate(leaf_area = (0.347*(leaf_length_cm*side_leaflet_length_cm*2) - 10.7)/100) %>% 
  #Blanco, F. F., & Folegatti, M. V. (2003).
  # A new method for estimating the leaf 
  # area index of cucumberand tomato plants.
  # Horticultura Brasileira, 21, 666-669.
  mutate(lai = leaf_area*num_leaves_from_firstcluster_to_det/(4*10))


#Yield related variables df
df_yield <- df_clean_vars[,c(108:124,127)]
names(df_yield) <- c("ripe_fruits_harvested",
                     "harvest_08_04_2024",
                     "harvest_11_04_2024",
                     "harvest_17_04_2024",
                     "harvest_20_04_2024",
                     "harvest_24_04_2024",
                     "harvest_26_04_2024",
                     "harvest_01_05_2024",
                     "harvest_05_05_2024",
                     "harvest_13_05_2024",
                     "harvest_20_05_2024",
                     "harvest_27_05_2024",
                     "harvest_03_06_2024",
                     "harvest_10_06_2024",
                     "harvest_18_06_2024",
                     "harvest_24_07_2024",
                     "harvest_02_07_2024",
                     "weight_six_fruits")
df_yield <- cbind(df_names, df_yield)

#Post harvested fruit qualities variables df
df_fruitqual <- df_clean_vars[,127:152]
names(df_fruitqual) <- c("weight_6_fruits",
                         "fruit_color_L_1",
                         "fruit_color_a_1",
                         "fruit_color_b_1",
                         "fruit_color_L_2",
                         "fruit_color_a_2",
                         "fruit_color_b_2",
                         "fruit_color_L_3",
                         "fruit_color_a_3",
                         "fruit_color_b_3",
                         "fruit_polar_diameter_1",
                         "fruit_equatorial_diameter_1",
                         "fruit_polar_diameter_2",
                         "fruit_equatorial_diameter_2",
                         "fruit_polar_diameter_3",
                         "fruit_equatorial_diameter_3",
                         "penetrometer_fruit_1",
                         "penetrometer_fruit_2",
                         "penetrometer_fruit_3",
                         "num_locules",
                         "brix",
                         "ph",
                         "sample_weight",
                         "ml_naoh",
                         "naoh_fc",
                         "viscosity")
df_fruitqual <- cbind(df_names, df_fruitqual)

#DNA extraction QC checks and DNA Index variables df
df_dna <- df_clean_vars[,c(159:162,164:166,175:189)]
names(df_dna) <- c("plate",
                   "well_line",
                   "well_column",
                   "well",
                   "eppendorf_num",
                   "dna_extraction_date",
                   "qubit_ng_per_ul",
                   "absorbance_230nm",
                   "absorbance_260nm",
                   "absorbance_280nm",
                   "absorbance_320nm",
                   "absorbance_minuscontrol_230nm",
                   "absorbance_minuscontrol_260nm",
                   "absorbance_minuscontrol_280nm",
                   "absorbance_minuscontrol_320nm",
                   "absorbance_900nm",
                   "absorbance_975nm",
                   "absorbance_sub320nm",
                   "dna_yield_ng_per_ul",
                   "260-280_ratio",
                   "260-230_ratio",
                   "electrophoresis_gel_quality")
df_dna <- df_dna %>%
  mutate(electrophoresis_gel_quality = factor(electrophoresis_gel_quality,
                                              levels = c("Ótimo", "Bom", "Fraco", "", "Ruim"),
                                              labels = c("Great", "Good", "Weak", NA, "Bad")))
df_dna <- cbind(df_names, df_dna)

#Stress variables df
df_stress <- cbind(df_names, df_clean_vars[,37:83])


trait <- c(
  "height_to_det",
  "num_sideshoot",
  "total_leaf_sideshoot",
  "total_cluster_sideshoot",
  "sideshoot_length",
  "hypocotyl_diam",
  "num_clusters_to_det",
  "num_leaves_to_firstcluster",
  "num_leaves_from_firstcluster_to_det",
  "num_leaflets",
  "internode_length",
  "leaf_area",
  "lai",
  "leaf_angle",
  "leaflet_angle"
)

name <- c(
  "Height to last leaf ($m$)",
  "N° of side shoots",
  "N° of leaf side shoots",
  "N° of cluster side shoots",
  "Side shoot length ($cm$)",
  "Hypocotyl diameter ($mm$)",
  "N° of clusters to last leaf",
  "N° of leaves up to first cluster",
  "N° of leaves after first cluster",
  "N° of Leaflets",
  "Internode length ($cm$)",
  "Leaf area ($dm^{2}$)",
  "Leaf Area Index (LAI)",
  "Leaf angle ($^{\\circ}$)",
  "Leaflet angle ($^{\\circ}$)"
)

traits <- data.frame(cbind(trait, name))

x <- acomp(df_pheno[traits[,"trait"]])
pcx <- princomp(x)
sum(pcx$sdev[1:2]^2)/sum(pcx$sdev^2)

pheno_fa <- factanal(na.omit(df_pheno[traits[,"trait"]]+0.00001), factors = 7)
pheno_fa

save(df_names, df_pheno, traits, file = "./PhenotypicData.RData")
save(pcx, irrig.idw, file = "./PlotData.RData")
