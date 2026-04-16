##########################################
#                                        #
#        GMA Project - Means Tables      #
#                                        #
##########################################

setwd("/home/deniz/workspace/GMA_Scripts/")
load("PhenotypicData.RData")
load("modelsGLMMs.RData")
source("./0_Formulas.R")


#packages
library(asreml)
library(tidyverse)

traits <- read.csv("traitsGLMMres.csv", header = TRUE)

# Model filtering
model_sel <- list(
  wospcor = traits$mod - 1,
  woirr = traits$woirr,
  full = traits$mod
)

ics <- map(names(model_sel), function(model_name) {
  s <- model_sel[[model_name]]
  model_files <- paste(
    traits$abrv,
    traits$fam,
    sprintf("%02d", s),
    "tom.asr",
    sep = "_"
  )

  models <- keep_at(asremls_bin, model_files)
  names(models) <- traits$abrv
  return(models)
}) |>
  setNames(names(model_sel)) |>
  modify_depth(3, summary)


tests <- c("aic", "bic")
ics_df <- map(tests, function(test) {
  table <- ics |>
    map_depth(2, \(x) x$model[[test]][1]) |>
    enframe() |>
    unnest_longer(value) |>
    pivot_wider(names_from = name, values_from = value) |>
    mutate(
      difirr = full - woirr,
      perirr = 100 * (woirr / abs(woirr)) * difirr / woirr,
      difspcor = full - wospcor,
      perspcor = 100 * (wospcor / abs(wospcor)) * difspcor / wospcor
    )
  return(table)
}) |>
  setNames(tests)

ics_df


#LRT statistic
lrts <- list()
for (i in 1:nrow(traits)) {
  abrv <- traits$abrv[i]
  t <- traits$trait[i]
  f <- traits$fam[i]

  n <- try - 1

  for (j in 1:(n - 1)) {
    for (k in (j + 1):n) {
      mod1 <- paste(abrv, f, sprintf("%02d", j), "tom.asr", sep = "_")
      mod2 <- paste(abrv, f, sprintf("%02d", k), "tom.asr", sep = "_")

      # Handle LRT with error catching
      llrt <- tryCatch(
        {
          lrt(asremls[[mod1]]$model, asremls[[mod2]]$model)
        },
        error = function(e) {
          return(NULL)
        }
      )

      # Create entry with real results or NA if error
      lrts[[paste0(abrv, "_lrt", j, "x", k)]] <- if (is.null(llrt)) {
        data.frame(
          trait = t,
          lrt = paste0(j, "x", k),
          df = NA,
          LR.statistic = NA,
          Pr.Chisq = NA
        )
      } else {
        data.frame(
          trait = t,
          lrt = paste0(j, "x", k),
          df = llrt$df,
          LR.statistic = llrt$"LR-statistic",
          Pr.Chisq = llrt$"Pr(Chisq)"
        )
      }
    }
  }
}


lrtss <- enframe(lrts) |>
  unnest_wider(value, names_sep = "_") |>
  rename_with(~ sub("^value_", "", .x)) |>
  na.omit()
lrtss$name[lrtss$Pr.Chisq < 0.05]

#Residuals of best models
res_list <- list()
for (i in 1:nrow(traits)) {
  abrv <- traits$abrv[i]
  t <- traits$trait[i]
  m <- traits$mod[i]
  f <- traits$fam[i]
  # # Debug
  #   print(paste0("tom", m, "_", t,".asr"))
  # }
  asr <- asremls[[
    paste(abrv, f, sprintf("%02d", m), "tom.asr", sep = "_")
  ]]$model

  varcomp <- summary(asr)$varcomp
  varcomp <- varcomp |> rownames_to_column()

  res_list[[abrv]] <- cbind(
    abrv,
    varcomp
  )
}
residuals <- bind_rows(res_list, .id = "abrv") |>
  arrange(rowname)

meanstable <- list()
df <- list()

for (i in 1:nrow(traits)) {
  abrv <- traits$abrv[i]
  t <- traits$trait[i]
  m <- traits$mod[i]
  f <- traits$fam[i]
  # # Debug
  # 	print(paste(abrv,f, sprintf("%02d", m),"tom.asr", sep = "_")))
  # }

  asrfam <- switch(
    f,
    "asr_gaussian(link = 'identity', dispersion = NA)",
    "asr_gaussian(link = 'log', dispersion = NA)",
    "asr_poisson(link = 'identity', dispersion = NA)",
    "asr_poisson(link = 'log', dispersion = NA)",
    "asr_poisson(link = 'sqrt', dispersion = NA)"
  )

  asr <- asremls[[
    paste(abrv, f, sprintf("%02d", m), "tom.asr", sep = "_")
  ]]$model

  pred <- summary(asr, coef = TRUE)$coef.fixed
  # if (f %in% c(2, 4)) {
  #   pred = exp(pred)
  # } else {
  #   pred = pred
  # }
  pred[, "solution"] <- pred[, "solution"] + pred["(Intercept)", "solution"]

  pred <- pred[
    startsWith(rownames(pred), "geno") & !endsWith(rownames(pred), ":irr_flw"),
  ]

  pred <- data.frame(pred) |>
    rownames_to_column("id") |> # Convert rownames to a column named "id"
    separate(
      col = id, # Column to split
      into = c("generation", "ind"), # New column names
      sep = "-", # Split on ":"
      remove = TRUE,
      # Keep original "id" column (set to TRUE to drop it)
      fill = "right"
    ) |>
    mutate(generation = str_remove(generation, "^geno_")) |>
    mutate(solution = ifelse(is.na(std.error), NA, solution))

  predseg <- pred |> filter(!is.na(ind))
  predpar <- pred |> filter(is.na(ind)) |> select(-ind)

  res <- data.frame(res = residuals(asr))

  df[[abrv]] <- df_pheno |>
    select(c(names(df_names), any_of(t))) |>
    left_join(predseg, by = c("generation", "ind")) |>
    left_join(predpar, by = "generation", suffix = c("", "_gen")) |>
    mutate(solution = coalesce(solution, solution_gen)) |>
    rename(value = t) |>
    bind_cols(res) |>
    # mutate(solution = if_else(
    # 			  generation %in% c("P1", "P2", "F1"),
    # 			  solution + res,
    # 			  solution)) |>
    select(c(names(df_names), value, solution, res))

  if (f %in% c(2, 4)) {
    df[[abrv]]$value = log(df[[abrv]]$value)
  }

  # Define names for the list elements
  rawvar <- meanstable[[paste0("rawvar_", t)]]
  predvar <- meanstable[[paste0("predvar_", t)]]

  rawvar <- df[[abrv]] |>
    na.omit() |>
    group_by(generation) |>
    # filter(between(value, quantile(value, 0.01), quantile(value, 0.99))) |>
    summarise(mean = mean(value), var = var(value), n = n()) |>
    arrange(factor(
      generation,
      levels = c("P1", "P2", "F1", "F2", "BC1", "BC2")
    )) |>
    mutate(meanvar = var / n)

  predvar <- df[[abrv]] |>
    group_by(generation) |>
    na.omit() |>
    summarise(mean = mean(solution), var = var(solution), n = n()) |>
    arrange(factor(
      generation,
      levels = c("P1", "P2", "F1", "F2", "BC1", "BC2")
    )) |>
    mutate(meanvar = var / n)

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

df_tot <- bind_rows(df, .id = "abrv") |>
  pivot_wider(
    names_from = abrv,
    values_from = c(value, solution, res),
    names_vary = "slowest"
  )

meanstb <- bind_rows(meanstable, .id = "ds_trait") |>
  separate_wider_delim(
    ds_trait,
    delim = "_",
    names = c("Dataset", "trait"),
    too_many = "merge"
  )

# Step 1: Empty list to store rows
variance_list <- list()

# Step 2: Loop through each data frame

for (i in 1:nrow(traits)) {
  t <- traits$trait[i]

  for (v in c("rawvar", "predvar")) {
    dftab <- meanstable[[paste0(v, "_", t)]]

    # Convert values into named vectors
    sig2 <- setNames(dftab$var, dftab$generation)
    mu <- setNames(dftab$mean, dftab$generation)
    n <- setNames(dftab$n, dftab$generation)

    # Estimate the environmental variance for the F2 population
    sig2env <- (1 / 4) * (sig2[["P1"]] + sig2[["P2"]] + 2 * sig2[["F1"]])
    sig2gF2 <- sig2[["F2"]] - sig2env
    sig2gF2 <- ifelse(sig2gF2 < 0, 0, sig2gF2)

    # Additive and Dominance Variance on F2
    sig2a <- 2 * sig2[["F2"]] - (sig2[["BC1"]] + sig2[["BC2"]])
    sig2a <- ifelse(sig2a < 0, 0, sig2a)

    sig2d <- sig2gF2 - sig2a
    sig2d <- ifelse(sig2d < 0, 0, sig2d)

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
      trait = as.factor(t), # Indicate whether from rawvar or predvar
      Dataset = as.factor(v), # Indicate whether from rawvar or predvar
      sig2env = sig2env,
      # sig2P1 = sig2[["P1"]],
      # sig2P2 = sig2[["P2"]],
      # sig2F1 = sig2[["F1"]],
      # sig2F2 = sig2[["F2"]],
      sig2gF2 = sig2gF2,
      sig2a = sig2a,
      sig2d = sig2d,
      h2b = h2b,
      h2n = h2n
    )

    # Append the new row to the list
    variance_list[[paste0(t, "_", v)]] <- new_row
  }
}


# Step 4: Combine list into a data frame
variance <- bind_rows(variance_list)
variance[3:length(variance)] <- round(variance[3:length(variance)], 4)
variance[7:length(variance)] <- round(variance[7:length(variance)], 1)
options(scipen = 99)

save(
  df_names,
  df_tot,
  meanstable,
  variance,
  traits,
  ics_df,
  file = "./MeansTable.RData"
)

add_object_to_rda(asremls, "./PlotData.RData", overwrite = TRUE)
