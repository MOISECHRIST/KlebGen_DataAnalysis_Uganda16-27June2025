#Import Libraries
library(conflicted)
library(tidyverse)
library(ggalluvial)
pacman::p_load(rio)

setwd("~/Documents/Taf Stats/KlebGen_Presentation_Kampala")

#Import data
merged_data <- import("./data/Klebsiella pneumoniae sequencing tracking Africa CDC_Cameroon.xlsx", which = "LNSP-MergeData")

