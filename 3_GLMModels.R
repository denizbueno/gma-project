#########################################
#                                       #
#          GMA Project - Models         #
#                                       #
#########################################

setwd("~/Desktop/Doutorado UFV/FIT 799 - Pesquisa/Scripts/GMA Scripts")
load("PhenotypicData.RData")
load("GenotypicData.RData")
load("GLMMs.RData")
source("./0_Formulas.R")

#packages
library(asreml)
library(tidyverse)

##Model Fitting

traits <- c(
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
  "height_to_det",
  "internode_length",
  "leaf_area",
  "lai",
  "leaf_angle",
  "leaflet_angle"
  )

df <- df_pheno %>% select(names(df_names), all_of(traits))

# lrts <- list()
# asremls <- list()

asreml.options(gammaPar = TRUE,
               ai.sing = FALSE)

for (trait in traits){
  df <- df_pheno
  
  tom1.lm <- glm(as.formula(paste0(trait," ~ geno + line + plant_num +
                                       line:irr_flw")),
                 data = df)
  
  #MODEL 1
  tom1.asr <- asreml(as.formula(paste0(trait," ~ geno + line + plant_num")),
                     residual = ~ units,
                     maxit = 20,
                     workspace = 2.048e9,
                     data = df)

  #MODEL 2
  tom2.asr <- asreml(as.formula(paste0("", trait, " ~ geno + line + plant_num")),
                     random = ~ geno:(irr_flw-mean(irr_flw)),
                     residual = ~ar1v(line):ar1(plant_num),
                     maxit = 20,
                     na.action = na.method(y = "include"),
                     workspace = 2.048e9,
                     data = df)
  
  asremls[[paste0("tom1_",trait,".lm")]] <- tom1.lm
  asremls[[paste0("tom1_",trait,".asr")]] <- tom1.asr
  asremls[[paste0("tom2_",trait,".asr")]] <- tom2.asr
  
  lrt <- lrt(asremls[[paste0("tom1_",trait,".asr")]],
             asremls[[paste0("tom2_",trait,".asr")]])
  
  lrts[[paste0("lrt_",trait)]] <- data.frame("trait" = trait,
                                             "df" = lrt$df,
                                             "LR-statistic" = lrt$"LR-statistic",
                                             "Pr(Chisq)" = lrt$"Pr(Chisq)")
  
  # # #REML log-likelihood, random components and wald statistics from the fit
  # tom2.asr$loglik
  # summary(tom2.asr)$varcomp
  # summary(tom2.asr)$bic
  # wald(tom2.asr)
}


save(asremls, lrts, file = paste0("./GLMMs.RData"))

add_object_to_rda(trait, "./PlotData.RData", overwrite = TRUE)
add_object_to_rda(asremls, "./PlotData.RData", overwrite = TRUE)