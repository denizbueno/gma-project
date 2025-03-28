##########################################
#                                        #
#        GMA Project - Means Tables      #
#                                        #
##########################################

setwd("~/Desktop/Doutorado UFV/FIT 799 - Pesquisa/Scripts/GMA Scripts")
load("PhenotypicData.RData")
load("GLMMs.RData")
source("./0_Formulas.R")


#packages
library(tidyverse)

# traits <- c("height_to_det")

meanstable <- list()

for (trait in traits$trait){
  asr <- asremls[[paste0("tom2_", trait, ".asr")]]
  
  pred <- summary(asr, coef = TRUE)$coef.fixed +
    summary(asr, coef = TRUE)$coef.fixed["(Intercept)","solution"]
  
  pred <- pred[startsWith(rownames(pred), "geno") &
                 !endsWith(rownames(pred), ":irr_flw"), ]
  
  pred <- data.frame(pred) %>%
    rownames_to_column("id") %>%  # Convert rownames to a column named "id"
    separate(
      col = id,                       # Column to split
      into = c("generation", "ind"),  # New column names
      sep = "-",                      # Split on ":"
      remove = TRUE,
      fill = "right") %>%              # Keep original "id" column (set to TRUE to drop it)
    mutate(generation = str_remove(generation, "^geno_")) %>% 
    mutate(solution = ifelse(is.na(std.error), NA, solution))
  
  predseg <- pred %>% filter(!is.na(ind))
  predpar <- pred %>% filter(is.na(ind)) %>% select(-ind)
  
  # res <- data.frame(res = residuals(asr))
  # sd <- sqrt(summary(asr)$
  #               varcomp["line:plant_num!R","component"])
  # set.seed(1234)     #Reset seed for each trait
  # res <- data.frame(res = rnorm(600, sd = sd, mean = 0))
  
  # Perform two left joins and coalesce the results
  df <- df_pheno %>%
    select(c(names(df_names), all_of(trait))) %>% 
    left_join(predseg, by = c("generation", "ind")) %>%
    left_join(predpar, by = "generation", suffix = c("", "_gen")) %>% 
    mutate(solution = coalesce(solution, solution_gen)) %>%
    # bind_cols(res) %>%
    # mutate(solution = if_else(
    #   generation %in% c("P1", "P2", "F1"),
    #   solution + res,
    #   solution)) %>%
    # select(c(names(df_names), all_of(trait), solution, res))
    select(c(names(df_names), all_of(trait), solution))
  
  # Define names for the list elements
  rawvar <- meanstable[[paste0("rawvar_", trait)]]
  predvar <- meanstable[[paste0("predvar_", trait)]]
  
  rawvar <- df[!is.na(df[[trait]]),] %>% 
    group_by(generation) %>% 
    na.omit() %>% 
    summarise(mean = mean(.data[[trait]]), var = var(.data[[trait]]), n = n()) %>% 
    arrange(factor(generation, levels = c("P1", "P2", "F1", "F2", "BC1", "BC2"))) %>% 
    mutate(meanvar = var/n)
  
  predvar <- df[!is.na(df$solution ),] %>% 
    group_by(generation) %>% 
    na.omit() %>% 
    summarise(mean = mean(solution ), var = var(solution ), n = n()) %>% 
    arrange(factor(generation, levels = c("P1", "P2", "F1", "F2", "BC1", "BC2"))) %>% 
    mutate(meanvar = var/n)
  
  for (g in c("P1", "P2", "F1")) {
    predvar[predvar$generation == g,]$n <- rawvar[rawvar$generation == g,]$n
    predvar[predvar$generation == g,]$var <- pred[pred$generation == g,]$std.error^2 /
      predvar[predvar$generation == g,]$n
    predvar[predvar$generation == g,]$meanvar <- predvar[predvar$generation == g,]$var /
      predvar[predvar$generation == g,]$n
  }
  
  # Define names for the list elements
  meanstable[[paste0("rawvar_", trait)]] <- rawvar
  meanstable[[paste0("predvar_", trait)]] <- predvar
}

# Step 1: Empty list to store rows
variance_list <- list()

# Step 2: Loop through each data frame

for (trait in traits$trait){
  for (v in c("rawvar", "predvar")) {
    dftab <- meanstable[[paste0(v,"_",trait)]]
    
    # Convert values into named vectors
    sig2 <- setNames(dftab$var, dftab$generation)
    mu <- setNames(dftab$mean, dftab$generation)
    n <- setNames(dftab$n, dftab$generation)
    
    # Estimate the environmental variance for the F2 population
    sig2env <- (1/4) * (sig2[["P1"]] + sig2[["P2"]] + 2 * sig2[["F1"]])
    sig2gF2 <- sig2[["F2"]] - sig2env
    
    # Additive and Dominance Variance on F2
    sig2a <- 2 * sig2[["F2"]] - (sig2[["BC1"]] + sig2[["BC2"]])
    sig2d <- sig2gF2 - sig2a
    
    # Heritability Estimates
    h2b <- (100 * sig2gF2) / sig2[["F2"]]
    h2n <- (100 * sig2a) / sig2[["F2"]]
    
    # Step 3: Store results as a data frame
    new_row <- data.frame(
      Trait = as.factor(trait),  # Indicate whether from rawvar or predvar
      Dataset = as.factor(v),  # Indicate whether from rawvar or predvar
      sig2env = sig2env,
      sig2gF2 = sig2gF2,
      sig2a = sig2a,
      sig2d = sig2d,
      h2b = h2b,
      h2n = h2n
    )
    
    # Append the new row to the list
    variance_list[[paste0(trait,"_",v)]] <- new_row
  }
}

# Step 4: Combine list into a data frame
variance <- bind_rows(variance_list)
variance[3:length(variance)] <- round(variance[3:length(variance)], 4)
variance[7:length(variance)] <- round(variance[7:length(variance)], 1)

res_list <- list()

for (trait in traits$trait){
   
  varcomp <- summary(asremls[[paste0("tom2_", trait, ".asr")]])$varcomp
  varcomp <- varcomp %>% rownames_to_column()
  
  res_list[[trait]] <- cbind(
   trait,
   varcomp
  )
}

residuals <- bind_rows(res_list)

lrtest <- do.call(rbind.data.frame, lrts)

save(meanstable, variance, traits, file = "./MeansTable.RData")
add_object_to_rda(trait, "./PlotData.RData", overwrite = TRUE)