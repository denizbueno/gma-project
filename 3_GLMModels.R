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
# bics <- list()
# asremls <- list()
# asremls_bin <- list()

# try = 1

# traits <- data.frame(
#   abrv = c("nfh"),
#   trait = c("ripe_fruits_harvested"),
#   name = c("Number of fruits harvested"),
#   fam = c(4),
#   mod = c(12)
# )
# traits <- data.frame(
#   abrv = c("colL", "colH", "colC"),
#   trait = c("leaf_color_L", "leaf_color_H", "leaf_color_C"),
#   name = c("Leaf color (L)", "Leaf color (Hue)", "Leaf color (Chroma)"),
#   fam = c(1, 1, 1),
#   mod = c(12, 12, 12)
# )
traits <- read.csv("traitsGLMMres.csv", header = TRUE)

asreml.options(gammaPar = TRUE, ai.sing = FALSE)

# Define model options with paired fixed and random effects
model_options <- list(
  # list(
  #   fixed = ~ line + plant_pos + irr_flw + geno,
  #   random = NULL
  # ),
  # list(
  #   fixed = ~ line + plant_pos*irr_flw + geno,
  #   random = NULL
  # ),
  # list(
  #   fixed = ~ line + plant_pos + geno,
  #   random = ~ plant_pos:irr_flw + irr_flw
  # ),
  # list(
  #   fixed = ~ geno,
  #   random = ~ line + plant_pos + irr_flw
  # ),
  # list(
  #   fixed = ~ geno,
  #   random = ~ line + plant_pos*irr_flw
  # ),
  list(
    fixed = ~ line + plant_pos + geno,
    random = NULL
  ),
  list(
    fixed = ~geno,
    random = ~ line + plant_pos
  )
)

# Define residual options
residual_options <- list(
  # ~ units,
  ~ ar1v(line):ar1(plant_pos)
)

# Create the model list
mods <- list()
for (i in seq_along(model_options)) {
  for (j in seq_along(residual_options)) {
    fixed_expr <- model_options[[i]]$fixed
    random_expr <- model_options[[i]]$random
    residual_expr <- residual_options[[j]]

    # Build the base model call without random effects
    base_call <- bquote(
      asreml(
        fixed = as.formula(paste0(t, " ", .(deparse(fixed_expr)))),
        residual = .(residual_expr),
        na.action = na.method(y = "include"),
        family = eval(parse(text = asrfam)),
        data = df,
        maxit = 100,
        workspace = 2048000000
      )
    )

    # Add random effects if specified
    if (!is.null(random_expr)) {
      # Convert to list, add random, then convert back to call
      call_list <- as.list(base_call)
      call_list$random <- random_expr

      # Reorder the elements to place random right after fixed
      # Get the names of all elements
      arg_names <- names(call_list)

      # Find the position of 'fixed'
      fixed_idx <- which(arg_names == "fixed")

      # Create a new order: all elements up to and including 'fixed',
      # then 'random', then the rest
      new_order <- c(
        1:fixed_idx,
        which(arg_names == "random"),
        setdiff((fixed_idx + 1):length(call_list), which(arg_names == "random"))
      )

      # Reorder the call list
      call_list <- call_list[new_order]

      model_call <- as.call(call_list)
    } else {
      model_call <- base_call
    }

    mods <- c(mods, list(model_call))
  }
}

start_time <- Sys.time()
for (i in 1:nrow(traits)) {
  # for (f in 1:4){
  df <- df_pheno

  abrv <- traits$abrv[i]
  t <- traits$trait[i]
  f <- traits$fam[i]

  df$plant_pos <- as.factor(df$plant_pos)

  #LinkFunctions
  asrfam <- switch(
    f,
    "asr_gaussian(link = 'identity', dispersion = NA)",
    "asr_gaussian(link = 'log', dispersion = NA)",
    "asr_poisson(link = 'identity', dispersion = NA)",
    "asr_poisson(link = 'log', dispersion = NA)",
    "asr_negative.binomial(link = 'log', dispersion = NA, phi = 1)"
  )

  #model
  tom.asr <- asreml(
    fixed = as.formula(
      paste0(t, " ~ line:irr_flw + plant_num:irr_flw + irr_flw + geno")
    ),
    random = ~ line + plant_num,
    residual = ~ ar1v(line):ar1(plant_pos),
    na.action = na.method(y = "include"),
    family = eval(parse(text = asrfam)),
    data = df,
    maxit = 100,
    workspace = 2.048e9
  )

  for (k in 1:length(mods)) {
    print(
      paste(
        abrv,
        t,
        i,
        paste0("fam: ", f),
        paste0("mod: ", try + k - 1),
        Sys.time(),
        sep = " | "
      )
    )

    tryCatch(
      {
        tom.asr <- eval(parse(text = mods[k]))
        asr <- paste(
          abrv,
          f,
          sprintf("%02d", try + k - 1),
          "tom.asr",
          sep = "_"
        )

        #store only sucessful
        if (tom.asr$converge == TRUE) {
          asremls[[asr]] <- list(
            model = tom.asr,
            rundate = Sys.time(),
            bics = summary(tom.asr)$bic
          )
        }
      },
      error = function(e) {
        message("Error at trait ", t, ", model ", k, ": ", conditionMessage(e))
      }
    )
    # try <- k
  }
  # }
}

try <- try + k

# #REML log-likelihood, random components and wald statistics from the fit
# plot(tom.asr)
# wald(tom.asr)
# plot(varioGram(tom.asr))

# asremls <- discard_at(
#   asremls,
#   paste("llag", 2, sprintf("%02d", c(1:15)), "tom.asr", sep = "_")
# )

asremls_bin <- c(asremls, asremls_bin)
asremls_bin <- asremls_bin[union(names(asreml), names(asremls_bin))]
# mods <- unique(map(asremls_bin[sort(names(asremls_bin))], ~.x$model$call))

save(
  asremls,
  asremls_bin,
  traits,
  try,
  mods,
  file = paste0("./modelsGLMMs.RData")
)
add_object_to_rda(asremls, "./PlotData.RData", overwrite = TRUE)
add_object_to_rda(mods, "./PlotData.RData", overwrite = TRUE)

end_time <- Sys.time()
time_diff <- difftime(end_time, start_time, units = "min")
print(time_diff)
