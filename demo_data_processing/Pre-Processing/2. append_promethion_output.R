################################################################################
# Promethion Data Processing
#
# This script processes Promethion metabolic data and merges it with the metadata
# from a separate file. It calculates the mean values of various metabolic
# parameters (kcal/hr, RQ, VO2, VCO2) for each timepoint and adds them to the
# metadata.
#
# Dependencies:
#   - tidyverse
#   - readxl
#   - janitor
#
# Author: Nelly Amenyogbe, Adrien Eynaud & Danny ___
# Date: Wednesday 6th March 2024
################################################################################

# Load required packages
library(tidyverse)
library(readxl)
library(janitor)

# Source helper functions (if applicable)
source("Pre-processing/helper_functions.R")

# Set working directory (adjust as needed)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Define the parent folder containing the Promethion data
FOLDER <- "./Promethion_Outputs/"

# Set the experiment date (used for matching Promethion IDs)
EXP_DATE <- "300322"

# Data import and naming
# Format is a nested list with dataframes at the bottom
data_list <- lapply(list.dirs(FOLDER, recursive = FALSE), data_import)
names(data_list) <- str_sub(list.dirs(FOLDER, recursive = FALSE), -5, -1)  # Assumes box numbers are always length 5

# Omit any rows containing NA and generate means
means_list <- list()
kcal_means <- list()
rq_means <- list()
vo2_means <- list()
vco2_means <- list()

for (i in seq_along(data_list)) {
  data_list[[i]] <- lapply(data_list[[i]], na.omit)  # Omit NA rows
  means_list[[i]] <- lapply(data_list[[i]], function(x) summarise_all(x, mean))
  kcal_means[[i]] <- lapply(means_list[[i]], function(x) x[, grepl("kcal", colnames(x))])
  rq_means[[i]] <- lapply(means_list[[i]], function(x) x[, grepl("rq", colnames(x))])
  vo2_means[[i]] <- lapply(means_list[[i]], function(x) x[, grepl("vo2", colnames(x))])
  vco2_means[[i]] <- lapply(means_list[[i]], function(x) x[, grepl("vco2", colnames(x))])
}

# Convert means to named numeric vectors
kcal_means <- lapply(kcal_means, unlist)
kcal_means <- unlist(kcal_means)
rq_means <- lapply(rq_means, unlist)
rq_means <- unlist(rq_means)
vo2_means <- lapply(vo2_means, unlist)
vo2_means <- unlist(vo2_means)
vco2_means <- lapply(vco2_means, unlist)
vco2_means <- unlist(vco2_means)

# Reformat names to match metadata
names(kcal_means) <- gsub(".kcal_hr", "", names(kcal_means))
names(rq_means) <- gsub(".rq", "", names(rq_means))
names(vo2_means) <- gsub(".vo2", "", names(vo2_means))
names(vco2_means) <- gsub(".vco2", "", names(vco2_means))

# Import metadata and clean
metadata <- read.csv("Result/intermediate_metadata.csv")
metadata <- remove_empty(data.frame(metadata), which = "rows")

# Add unique identifier column
metadata <- metadata %>%
  mutate(cage_expdate_visit_pup = paste(cage, EXP_DATE, visit, pup_no, sep = "_"))

# Add metabolic data columns to metadata
metadata$kcal_hr <- kcal_means[metadata$cage_expdate_visit_pup]
metadata$rq <- rq_means[metadata$cage_expdate_visit_pup]
metadata$vo2 <- vo2_means[metadata$cage_expdate_visit_pup]
metadata$vco2 <- vco2_means[metadata$cage_expdate_visit_pup]

# Create kcal_hr_g and avg_score columns
metadata <- metadata %>%
  mutate(
    kcal_hr_g = as.numeric(kcal_hr) / as.numeric(weight),
    avg_score = ifelse(is.na(score.1), NA_real_, ceiling((score.1 + score.2) / 2))
  )

# Write result file
write.csv(metadata[, 2:ncol(metadata)], "Result/result_export.csv", row.names = FALSE)
