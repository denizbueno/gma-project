##########################################
#                                        #
#        GMA Project - Means Tables      #
#                                        #
##########################################

setwd("~/Desktop/Doutorado UFV/FIT 799 - Pesquisa/Scripts/GMA Scripts/")
load("PhenotypicData.RData")
load("modelsGLMMs.RData")
source("./0_Formulas.R")


#packages
library(asreml)
library(tidyverse)


#AIC models comparison
aics <- enframe(map(asremls, ~summary(.x$model)$aic)) |>
  unnest_wider(value, names_sep = "_") |>
  separate(name, sep = "_",
           into=c("abrv", "fam", "try", "sufix")) |>
  select(!c(sufix, fam)) |>
  pivot_wider(names_from = try, values_from = value_1) |>
  select(abrv, sprintf("%02d", 1:length(mods))) |>
  select(!c("01","07"))
n <- 2:ncol(aics)
aics$min_val <- apply(aics[,n], 1, function(x) min(x, na.rm = TRUE))
aics$min_col <- colnames(aics[,n])[apply(aics[,n], 1, which.min)]
aics$fam <- traits$fam


#BIC models comparison
bics <- enframe(map(asremls, ~.x$bics)) |>
  unnest_wider(value, names_sep = "_") |>
  separate(name, sep = "_",
           into=c("abrv", "fam", "try", "sufix")) |>
  select(!c(sufix, fam)) |>
  pivot_wider(names_from = try, values_from = value_1) |>
  select(abrv, sprintf("%02d", 1:length(mods))) |>
  select(!c("01","07"))
n <- 2:ncol(bics)
bics$min_val <- apply(bics[,n], 1, function(x) min(x, na.rm = TRUE))
bics$min_col <- colnames(bics[,n])[apply(bics[,n], 1, which.min)]
bics$fam <- traits$fam

# lrts <- list()
# for (i in 1:nrow(traits)) {
#   abrv <- traits$abrv[i]
#   t <- traits$trait[i]
#   f <- traits$fam[i]
#   
#   n <- try - 1
#   
#   for (j in 1:(n-1)) {
#     for (k in (j+1):n) {
#       
#       mod1 <- paste(abrv, f, sprintf("%02d", j), "tom.asr", sep = "_")
#       mod2 <- paste(abrv, f, sprintf("%02d", k), "tom.asr", sep = "_")
#       
#       # Handle LRT with error catching
#       llrt <- tryCatch({
#         lrt(asremls[[mod1]]$model, asremls[[mod2]]$model)
#       }, error = function(e) {
#         return(NULL)
#       })
#       
#       # Create entry with real results or NA if error
#       lrts[[paste0(abrv, "_lrt", j, "x", k)]] <- if (is.null(llrt)) {
#         data.frame(
#           trait = t,
#           lrt = paste0(j, "x", k),
#           df = NA,
#           LR.statistic = NA,
#           Pr.Chisq = NA
#         )
#       } else {
#         data.frame(
#           trait = t,
#           lrt = paste0(j, "x", k),
#           df = llrt$df,
#           LR.statistic = llrt$"LR-statistic",
#           Pr.Chisq = llrt$"Pr(Chisq)"
#         )
#       }
#     }
#   }
# }
# 
# 
# lrtss <- enframe(lrts) |> 
#   unnest_wider(value, names_sep = "_") |> 
#   rename_with(~ sub("^value_", "", .x)) |> 
#   na.omit()
# lrtss$name[lrtss$Pr.Chisq < 0.05]

traits <- read.csv("traitsGLMMres.csv", header = TRUE)

#Residuals of best models
res_list <- list()
for (i in 1:nrow(traits)){
  abrv <- traits$abrv[i]
  t <- traits$trait[i]
  m <- traits$mod[i]
  f <- traits$fam[i]
  # # Debug  
  #   print(paste0("tom", m, "_", t,".asr"))
  # }
  asr <- asremls[[paste(abrv,f, sprintf("%02d", m),"tom.asr", sep = "_")]]$model
  
  varcomp <- summary(asr)$varcomp
  varcomp <- varcomp |> rownames_to_column()
  
  res_list[[abrv]] <- cbind(
    t,
    varcomp
  )
}
residuals <- bind_rows(res_list)

# res_trait <- list()
# for (i in 1:4){
#   
#   t <- traits$trait[10]
#   m <- paste0("tom",traits$mod[i],"_", t, ".asr")
#   
#   varcomp <- summary(asremls[[m]])$varcomp
#   varcomp <- varcomp |> rownames_to_column()
#   
#   res_trait[[i]] <- cbind(
#     t,
#     varcomp
#   )
# }
# residuals <- bind_rows(res_trait)

meanstable <- list()

for (i in 1:nrow(traits)){
  abrv <- traits$abrv[i]
  t <- traits$trait[i]
  # m <- traits$mod[i]
  m <- 12
  f <- traits$fam[i]
# # Debug  
#   print(paste(abrv,f, sprintf("%02d", m),"tom.asr", sep = "_")))
# }
  
  asr <- asremls[[paste(abrv,f, sprintf("%02d", m),"tom.asr", sep = "_")]]$model
  
  pred <- summary(asr, coef = TRUE)$coef.fixed
  # if (f %in% c(2,4)){
  #   pred = exp(pred)
  #   } else {
  #     pred = pred
  #   }
  pred[,"solution"] <- pred[,"solution"] + pred["(Intercept)","solution"]
  
  pred <- pred[startsWith(rownames(pred), "geno") &
                 !endsWith(rownames(pred), ":irr_flw"), ]
  
  pred <- data.frame(pred) |>
    rownames_to_column("id") |>  # Convert rownames to a column named "id"
    separate(
      col = id,                       # Column to split
      into = c("generation", "ind"),  # New column names
      sep = "-",                      # Split on ":"
      remove = TRUE,
      fill = "right") |>              # Keep original "id" column (set to TRUE to drop it)
    mutate(generation = str_remove(generation, "^geno_")) |> 
    mutate(solution = ifelse(is.na(std.error), NA, solution))
  
  predseg <- pred |> filter(!is.na(ind))
  predpar <- pred |> filter(is.na(ind)) |> select(-ind)
  
  df <- df_pheno |>
    select(c(names(df_names), any_of(t))) |> 
    left_join(predseg, by = c("generation", "ind")) |>
    left_join(predpar, by = "generation", suffix = c("", "_gen")) |> 
    mutate(solution = coalesce(solution, solution_gen)) |> 
    # bind_cols(res) |>
    # mutate(solution = if_else(
    #   generation %in% c("P1", "P2", "F1"),
    #   solution + res,
    #   solution)) |>
    # select(c(names(df_names), t, solution, res))
    select(c(names(df_names), all_of(t), solution))
  
  # if (f %in% c(2,4)){
  #   df[t] = log(df[t])
  #   } else {
  #     df[t] = df[t]
  #   }
  
  # Define names for the list elements
  rawvar <- meanstable[[paste0("rawvar_", t)]]
  predvar <- meanstable[[paste0("predvar_", t)]]
  
  rawvar <- df[!is.na(df[[t]]),] |> 
    group_by(generation) |> 
    na.omit() |> 
    summarise(mean = mean(.data[[t]]), var = var(.data[[t]]), n = n()) |> 
    arrange(factor(generation, levels = c("P1", "P2", "F1", "F2", "BC1", "BC2"))) |> 
    mutate(meanvar = var/n)
  
  predvar <- df[!is.na(df$solution),] |> 
    group_by(generation) |> 
    na.omit() |> 
    summarise(mean = mean(solution), var = var(solution), n = n()) |> 
    arrange(factor(generation, levels = c("P1", "P2", "F1", "F2", "BC1", "BC2"))) |> 
    mutate(meanvar = var/n)
  
  sd <- setNames(pred$std.error, pred$generation)
  n <- setNames(predvar$n, predvar$generation)
  
  for (g in c("P1", "P2", "F1")) {
    # predvar[predvar$generation == g, ]$var <- sd[[g]]^2 / n[[g]]
    predvar[predvar$generation == g, ]$meanvar <- sd[[g]]^2 / n[[g]]
  }
  # predvar <- predvar |> mutate(meanvar = var / n)
  predvar <- predvar |> mutate(var = meanvar * n)
  
  # Define names for the list elements
  meanstable[[paste0("rawvar_", t)]] <- rawvar
  meanstable[[paste0("predvar_", t)]] <- predvar
}

# Step 1: Empty list to store rows
variance_list <- list()

# Step 2: Loop through each data frame

for (i in 1:nrow(traits)){
  
  t <- traits$trait[i]
  
  for (v in c("rawvar", "predvar")) {
    dftab <- meanstable[[paste0(v,"_",t)]]
    
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
    #   #Weber and Moorthy (1952)
    # h2b <- 100*(sig2[["F2"]] - (sig2[["P1"]]*sig2[["P2"]]*sig2[["F1"]])^(1/3))/sig2[["F2"]]
    #   #Burton (1951)
    # h2b <- 100*(sig2[["F2"]] - sig2[["F1"]]) / sig2[["F2"]]
    #   #Mahmud and Kramer (1951)
    # h2b <- 100*(sig2[["F2"]] - (sig2[["P1"]]*sig2[["P2"]])) / sig2[["F2"]]
    #   #Warner (1952)
    # h2n <- 100*(2*sig2[["F2"]]-(sig2[["BC1"]]+sig2[["BC2"]])) / sig2[["F2"]]
    
    
    # Step 3: Store results as a data frame
    new_row <- data.frame(
      Trait = as.factor(t),  # Indicate whether from rawvar or predvar
      Dataset = as.factor(v),  # Indicate whether from rawvar or predvar
      sig2env = sig2env,
      sig2gF2 = sig2gF2,
      sig2a = sig2a,
      sig2d = sig2d,
      h2b = h2b,
      h2n = h2n
    )
    
    # Append the new row to the list
    variance_list[[paste0(t,"_",v)]] <- new_row
  }
}


# Step 4: Combine list into a data frame
variance <- bind_rows(variance_list)
variance[3:length(variance)] <- round(variance[3:length(variance)], 4)
variance[7:length(variance)] <- round(variance[7:length(variance)], 1)
options(scipen = 99999)

save(meanstable, variance, traits, file = "./MeansTable.RData")