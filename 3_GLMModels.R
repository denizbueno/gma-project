#########################################
#                                       #
#          GMA Project - Models         #
#                                       #
#########################################

setwd("~/Desktop/Doutorado UFV/FIT 799 - Pesquisa/Scripts/GMA Scripts")
load("PhenotypicData.RData")
load("GenotypicData.RData")
load("modelsGLMMs.RData")
source("./0_Formulas.R")

#packages
library(asreml)
library(tidyverse)


# lrts <- list()
# asremls <- list()
# bics <- list()

traits <- read.csv("traitsGLMM.csv", header = TRUE)

asreml.options(gammaPar = TRUE,
               ai.sing = FALSE)

start_time <- Sys.time()
# for (i in 1:nrow(traits)){
for (i in 15){
  df <- df_pheno
  
  abrv <- traits$abrv[i]
  t <- traits$trait[i]
  f <- traits$fam[i]
  
  print(paste(abrv, t, i, Sys.time(), sep = " | "))

  #LinkFunctions
  asrfam <-switch(traits$fam[i],
               "asr_gaussian(link = 'identity', dispersion = NA)",
               "asr_gaussian(link = 'log', dispersion = NA)",
               "asr_poisson(link = 'identity', dispersion = NA)",
               "asr_poisson(link = 'log', dispersion = NA)",
               "asr_negative.binomial(link = 'log', dispersion = NA, phi = 1)")
  
  #MODEL
  tom.asr <- asreml(fixed = as.formula(paste0(t, " ~ geno")),
                    random = ~ line + plant_num + geno:irr_flw + irr_flw,
                    residual = ~ar1v(line):ar1(plant_num),
                    na.action = na.method(y = "include"),
                    family = eval(parse(text = asrfam)),
                    data = df, maxit = 500, workspace = 2.048e8)

  # for (k in 1:length(mods)){
  #   tom.asr <- eval(parse(text = mods[k]))
  #   # try <- length(mods) + k
  #   try <- k

    asremls[[paste(abrv,f, sprintf("%02d", try),"tom.asr", sep = "_")]] <- list(
      model = tom.asr,
      rundate = Sys.time(),
      bics= summary(tom.asr)$bic
      )
  # }
}

try = try+1

# #REML log-likelihood, random components and wald statistics from the fit
# plot(tom.asr)
# wald(tom.asr)
# plot(varioGram(tom.asr))

asremls_bin <- c(asremls, asremls_bin)
asremls_bin <- asremls_bin[union(names(asreml),names(asremls_bin))]
mods <- unique(map(asremls_bin[sort(names(asremls_bin))], ~.x$model$call))

# bin <- enframe(map(asremls, ~.x$model$converge)) |> filter(value == FALSE)
# bin <- data.frame(models = names(asremls_bin)) |>
#   separate(models,
#            into = c("abrv", "fam", "mods", "sufix"),
#            sep = "_") |>
#   # arrange(abrv, mods) |>
#   # filter(!(mods %in% sprintf("%02d", c(3,9)))) |>
#   unite(c(abrv:sufix), col = "name", sep = "_")

asremls <- asremls[bin$name]

  
save(asremls, asremls_bin, traits, try, mods, file = paste0("./modelsGLMMs.RData"))
add_object_to_rda(asremls, "./PlotData.RData", overwrite = TRUE)

end_time <- Sys.time()
time_diff <- difftime(end_time, start_time, units = "min")
print(time_diff)
