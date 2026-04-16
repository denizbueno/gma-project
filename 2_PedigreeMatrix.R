#########################################
#                                       #
#   GMA Project - Relationship Matrix   #
#                                       #
#########################################
setwd("~/Desktop/Doutorado UFV/FIT 799 - Pesquisa/Scripts/GMA Scripts")
load("./PhenotypicData.RData")

#packages
library(AGHmatrix)
library(tidyverse)
# source("./0_Formulas.R")

df <- df_names

## Generations as matrix
# Identity-by-Desent

PedTable <- df |> distinct(geno, dam, sire)

Endogamy <- data.frame(
  geno = factor(c(
    "P1F0",
    "P2F0",
    "P1F1",
    "P2F1",
    "P1F2",
    "P2F2", #Pre-breeding
    "P1F3",
    "P2F3",
    "P1F4",
    "P2F4",
    "P1F5",
    "P2F5", #Pre-breeding
    "P1F6",
    "P2F6",
    "P1F7",
    "P2F7",
    "P1F8",
    "P2F8", #Pre-breeding
    "P1F9",
    "P2F9",
    "F1F0" #Seed source
  )),
  dam = factor(c(
    0,
    0,
    "P1F0",
    "P2F0",
    "P1F1",
    "P2F1", #Pre-breeding
    "P1F2",
    "P2F2",
    "P1F3",
    "P2F3",
    "P1F4",
    "P2F4", #Pre-breeding
    "P1F5",
    "P2F5",
    "P1F6",
    "P2F6",
    "P1F7",
    "P2F7", #Pre-breeding
    "P1F8",
    "P2F8",
    "P2F8" #Seed source
  )),
  sire = factor(c(
    0,
    0,
    "P1F0",
    "P2F0",
    "P1F1",
    "P2F1", #Pre-breeding
    "P1F2",
    "P2F2",
    "P1F3",
    "P2F3",
    "P1F4",
    "P2F4", #Pre-breeding
    "P1F5",
    "P2F5",
    "P1F6",
    "P2F6",
    "P1F7",
    "P2F7", #Pre-breeding
    "P1F8",
    "P2F8",
    "P1F8" #Seed source
  ))
)

GenAMatrix <- rbind(Endogamy, PedTable)
# GenAMatrix <- ainverse(GenAMatrix)
GenAMatrix <- data.frame(Amatrix(GenAMatrix, ploidy = 2))
GenAMatrix <- GenAMatrix[-(1:21), -(1:21)]

save(GenAMatrix, file = "./GenotypicData.RData")
