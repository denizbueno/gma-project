##########################################
#                                        #
# GMA Project - Generation Mean Analysis #
#                                        #
##########################################

setwd("~/Desktop/Doutorado UFV/FIT 799 - Pesquisa/Scripts/GMA Scripts")
load("MeansTable.RData")
source("./0_Formulas.R")

#packages
library(tidyverse)

# dftab <- data.frame(generation = c("P1", "P2", "F1", "F2", "BC1", "BC2"),
#                     mean = c(44.0029, 32.2236, 34.3404, 35.8047, 42.9800, 32.5364),
#                     var = c(34.2978, 11.1927, 11.8896, 57.5939, 54.5265, 24.6219),
#                     n = c(15, 15, 16, 144, 50, 47))

gma_list <- list()
r2AD_list <- list()

for (i in 1:nrow(traits)){
  
  trait <- traits$trait[i]
  abrv <- traits$abrv[i]
  gma <- list()
  
  for (v in c("rawvar", "predvar")) {
    dftab <- meanstable[[paste0(v,"_",trait)]]
        
    # Convert values into named vectors
    sig2 <- setNames(dftab$var, dftab$generation)
    mu <- setNames(dftab$mean, dftab$generation)
    n <- setNames(dftab$n, dftab$generation)
    
    ##Model Y = m + a + d
    #Incidence Matrix
    im <- matrix(c(1, 1, 0, 1, 0, 0,
                   1, -1, 0, 1, 0, 0,
                   1, 0, 1, 0, 0, 1,
                   1, 0, 1/2, 0, 0, 1/4,
                   1, 1/2, 1/2, 1/4, 1/4, 1/4,
                   1, -1/2, 1/2, 1/4, -1/4, 1/4), 
                 nrow = 6, ncol = 6, byrow = TRUE)
    #Parameter Incidence Matrix
    parhat <- matrix(c(1/2, 1/2, 0, 4, -2, -2, 
                       1/2, -1/2, 0, 0, 0, 0,
                       -3/2, -3/2, -1, -8, 6, 6,
                       0, 0, 0, -4, 2, 2,
                       -1, 1, 0, 0, 2, -2,
                       1, 1, 2, 4, -4, -4), 
                     nrow = 6, ncol = 6, byrow = TRUE)
    
    dinv <- solve(diag(sig2/n))
    
    imAD <- im[,1:3]
    
    dinvAD <- solve(diag(sig2/n))
    
    xdx <- t(imAD)%*%dinv%*%imAD
    xdy <- t(imAD)%*%dinv%*%mu
    parhat <- solve(xdx,xdy)
    xdxinv <- solve(xdx)
    
    gmaAD <- data.frame(parameters = as.factor(c("mhat", "ahat", "dhat")),
                        estimates = parhat,
                        varcomp = diag(xdxinv),
                        sdcomp = sqrt(diag(xdxinv)),
                        tstat = parhat / sqrt(diag(xdxinv)),
                        df = sum(n-1),
                        pval = round(2 * (1 - pt(abs(parhat / sqrt(diag(xdxinv))),
                                                 df = sum(n-1))),8))
    
    muhat <- imAD %*% parhat
    r2AD <- cor(mu,muhat)^2
    te <- solve(xdx)
    w <- solve(xdx,xdy)
    
    r2AD_list[[paste0(v,"_",trait)]] <- r2AD
    
    sq <- w^2/diag(te)
    r2 <- ((w^2/diag(te))*100)/sum(w^2/diag(te))
    
    gmaAD <- data.frame(cbind(model = as.factor(round(r2AD*100,3)), 
                              gmaAD,
                              SQ = sq,
                              R2 = r2))
    
    #Model Y = m + a + d + aa + ad + dd
    
    scamat <- matrix(c(-1, 0, -1, 0, 2, 0,
                       0, -1, -1, 0, 0, 2,
                       -1, -1, -2, 4, 0, 0,
                       0, 0, 0, 2, -1, -1),
                     nrow = 4, ncol = 6, byrow = TRUE)
    
    scaling <- data.frame(parameters = c("A", "B", "C", "D"),
                          scau = scamat %*% mu,
                          scavar = scamat^2 %*% (sig2/n))
    
    parhat <- matrix(c(1/2, 1/2, 0, 4, -2, -2, 
                       1/2, -1/2, 0, 0, 0, 0,
                       -3/2, -3/2, -1, -8, 6, 6,
                       0, 0, 0, -4, 2, 2,
                       -1, 1, 0, 0, 2, -2,
                       1, 1, 2, 4, -4, -4), 
                     nrow = 6, ncol = 6, byrow = TRUE)
    
    degfree <- matrix(c(1, 1, 0, 1, 1, 1, 
                        1, 1, 0, 0, 0, 0,
                        1, 1, 1, 1, 1, 1,
                        0, 0, 0, 1, 1, 1,
                        1, 1, 0, 0, 1, 1,
                        1, 1, 1, 1, 1, 1), 
                      nrow = 6, ncol = 6, byrow = TRUE)
    
    gmaE <- data.frame(parameters = as.factor(c("mhat", "ahat", "dhat", 
                                    "aahat", "adhat", "ddhat")),
                      estimates = parhat %*% mu,
                      varcomp = parhat^2 %*% (sig2/n),
                      sdcomp = sqrt(parhat^2 %*% (sig2/n)),
                      tstat = (parhat %*% mu) / sqrt(parhat^2 %*% (sig2/n)),
                      df = degfree%*%(n-1),
                      pval = round(2 * (1 - pt(abs((parhat %*% mu) / sqrt(parhat^2 %*% (sig2/n))),
                                               df = degfree%*%(n-1))),8))
    
    # for (i in row.names(scaling)) {
    #   cat(i, "[", scaling[i,"scau"] - sqrt(scaling[i,"scavar"]),
    #       "to", scaling[i,"scau"] + sqrt(scaling[i,"scavar"]), "] ","\n" , sep = " ")
    # }
    
    im <- matrix(c(1, 1, 0, 1, 0, 0,
                   1, -1, 0, 1, 0, 0,
                   1, 0, 1, 0, 0, 1,
                   1, 0, 1/2, 0, 0, 1/4,
                   1, 1/2, 1/2, 1/4, 1/4, 1/4,
                   1, -1/2, 1/2, 1/4, -1/4, 1/4), 
                 nrow = 6, ncol = 6, byrow = TRUE)
    dinv <- solve(diag(sig2/n))
    
    xdx <- t(im)%*%dinv%*%im
    xdy <- t(im)%*%dinv%*%mu
    te <- solve(xdx)
    w <- solve(xdx,xdy)
    
    sq <- w^2/diag(te)
    r2 <- ((w^2/diag(te))*100)/sum(w^2/diag(te))
    
    gmaE <- data.frame(cbind(model = as.factor("epi"), gmaE, SQ = sq, R2 = r2))
    
    # Rounding digits for columns 3 to 10:
    round_digits <- c(3, 3, 3, 3, 0, 4, 3, 1)
    
    # Loop over columns 2 to 9
    for(j in seq_along(round_digits)) {
      col_index <- j + 2  # because our columns start at 3
      gmaAD[, col_index] <- round(gmaAD[, col_index], round_digits[j])
      gmaE[, col_index] <- round(gmaE[, col_index], round_digits[j])
    }
    
    gma[[v]] <- data.frame(
      Trait = as.factor(abrv),
      Dataset = as.factor(v),
      rbind(
        gmaAD,
        # "a+d" = c(rep(NA,9), round(r2AD*100, 1)),
        # Epi = rep(NA,9),
        gmaE))
    
    # gma_list[[paste0(trait,"_",v)]] <- gma
  }
  gma_list[[trait]] <- bind_rows(gma)
}

raw <- do.call(rbind.data.frame, r2AD_list[startsWith(names(r2AD_list), "raw")])
pred <- do.call(rbind.data.frame, r2AD_list[startsWith(names(r2AD_list), "pred")])
r2AD <- do.call(rbind.data.frame, r2AD_list)

delta <- pred |> 
  bind_cols(raw, delta = pred - raw) |> 
  rownames_to_column() |> 
  setNames(c("trait", "pred","raw","delta")) |> 
  mutate(trait = gsub("predvar_", "", trait))

r2AD <- do.call(rbind.data.frame, r2AD_list) |>
  rownames_to_column() |> 
  separate(
    col = rowname,                   # Column to split
    into = c("Dataset", "Trait"),    # New column names
    sep = "_",                       # Split on ":"
    extra = "merge",                 # Merge after
    remove = TRUE,
    fill = "right") |>              # Keep original "id" col
  right_join(variance, by = c("Trait", "Dataset")) |> 
  rename(r2 = "V1") |> 
  relocate(Trait) |> 
  relocate(r2, .after = h2n)

delta_list <- list()

for (trait in traits$trait){
  dgma <- gma_list[[trait]]
  pred <- dgma[dgma$Dataset == "predvar",5:12]
  raw <- dgma[dgma$Dataset == "rawvar",5:12]
  
  d <- pred - raw
  dgma <- dgma |> 
    filter(Dataset == "predvar") |> 
    select("Trait", "model", "parameters") |> 
    bind_cols(d)
  
  delta_list[[trait]] <- dgma
}

add_object_to_rda(gma_list, "./MeansTable.RData", overwrite = TRUE)
add_object_to_rda(r2AD_list, "./MeansTable.RData", overwrite = TRUE)