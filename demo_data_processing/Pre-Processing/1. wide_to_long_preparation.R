################################################################################
# Mouse Metadata Preparation: Wide to Long Data Transformation
#
# This script reads in a wide-format mouse metadata set and transforms it into a
# long format suitable for survival curves or other analyses. It calculates the hours post-challenge
# for experimental animals and assigns survival outcomes based on sacrifice types.
#
# The input data is expected to be in an Excel file named 'data_entry_template.xlsx'
# located in the same directory as the script.
#
# Dependencies:
#   - plyr
#   - dplyr
#   - lubridate
#   - tidyr
#   - readxl
#
# Author: Nelly Amenyogbe, Adrien Eynaud & Danny ___
# Date: Wednesday 6th of April 2024
################################################################################

# Load required packages
library(plyr)
library(dplyr)
library(lubridate)
library(tidyr)
library(readxl)

# Set working directory (adjust as needed)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')

# Load data
data_file <- "data_entry_template.xlsx"
data <- read_excel(data_file, skip = 2)

# Prepare long data format
metadata_vars <- c("pup_ID", "pup_no", "cage", "sex", "DOB_dam.dd.mm.yyyy", "DOB_pup.dd.mm.yyyy",
                   "Exclude", "Exclude_Notes", "AEC", "Experiment_type", "Experiment_name",
                   "challenge_stat", "challenge_dose", "sacrifice_type", "sacrifice_time.hh.mm",
                   "sacrifice_date.dd.mm.yyyy", "group", "treatment_visit", "treatment_dose",
                   "challenge_visit", "challenge_dose")

measure_vars <- setdiff(colnames(data), metadata_vars)

visit_cols <- data.frame(measure_vars) %>%
  separate(measure_vars, c("visit", "variable"), sep = "_")

visits <- unique(visit_cols$visit)
variables <- unique(visit_cols$variable)

data_long <- ldply(seq_along(visits), function(i) {
  vis_columns <- paste0(visits[i], "_", variables)
  dat_vis <- data[, c(metadata_vars, vis_columns)]
  colnames(dat_vis) <- c(metadata_vars, variables)
  dat_vis$visit <- visits[i]
  dat_vis
})

# Add visit type column
data_long$visit_type <- NA
for (i in seq_along(data_long$visit)) {
  chal_vis <- data_long$challenge_visit[i]
  treat_vis <- data_long$treatment_visit[i]
  
  if (is.na(chal_vis)) {
    chal_vis <- "NA"
  }
  
  if (is.na(treat_vis)) {
    treat_vis <- "NA"
  }
  
  if (data_long$visit[i] == chal_vis) {
    data_long$visit_type[i] <- "challenge"
  } else if (data_long$visit[i] == treat_vis) {
    data_long$visit_type[i] <- "treatment"
  } else if (chal_vis == "NA") {
    data_long$visit_type[i] <- NA
  } else {
    data_long$visit_type[i] <- "monitor"
  }
}

# Generate hours post challenge
data_long$date.dd.mm.yyyy <- gsub(".21", ".2021", data_long$date.dd.mm.yyyy, fixed = TRUE)
data_long$date.dd.mm.yyyy <- gsub("\\b.\\b", "/", as.character(data_long$date.dd.mm.yyyy))
data_long$time <- format(round(data_long$time.hh.mm, digits = 2), nsmall = 2)
data_long$time <- gsub("\\.", ":", data_long$time)
data_long$visit_date <- paste(data_long$date.dd.mm.yyyy, data_long$time)
data_long$visit_date <- dmy_hm(data_long$visit_date) # visit_date will fail to parse for some visits, as expected.

# NOTE: For the code below, when visit_date contains NA, lubridate::interval() will emit
#   “In with_tz.default(starts, tzone): Unrecognized time zone 'hours'”.
# This warning reflects the absence of a defined timestamp (e.g. an unmonitored visit)
# and gives an expected NA result; it does not indicate an underlying error
# and may therefore be disregarded.
data_long <- ddply(data_long, .(pup_ID, cage),
                   transform,
                   hours_post_challenge = interval(visit_date[visit_type == "challenge"],
                                                   visit_date, "hours") %>%
                     time_length(unit = "hour") %>%
                     round(digits = 2))

# Prepare data for survival analysis
data_long$sacrifice_date.dd.mm.yyyy <- gsub("\\b.\\b", "/", as.character(data_long$sacrifice_date.dd.mm.yyyy))
data_long$sacrifice_time.hh.mm <- sprintf("%.2f", data_long$sacrifice_time)
data_long$sacrifice_time.hh.mm <- as.character(data_long$sacrifice_time.hh.mm)
data_long$sacrifice_time.hh.mm <- gsub("\\.", ":", data_long$sacrifice_time)
data_long$sac_dt <- paste(data_long$sacrifice_date.dd.mm.yyyy, data_long$sacrifice_time.hh.mm)
data_long$sac_dt <- dmy_hm(data_long$sac_dt) # again, visit_date will fail to parse for some visits, as expected.

# As above: NA visit_date triggers expected lubridate tz warning
data_long <- ddply(data_long, .(pup_ID, cage),
                   transform,
                   sac_hpc = interval(visit_date[visit_type == "challenge"],
                                      sac_dt, "hours") %>%
                     time_length(unit = "hour") %>%
                     round(digits = 2))

data_long$outcome <- ifelse(data_long$sacrifice_type %in% c("HE", "FD"), 1,
                            ifelse(data_long$sacrifice_type == "EE", 0, data_long$sacrifice_type))

unique(data_long$outcome)  # Should be a vector containing 1, 0, or NA for dams.

# Assign "monitor" to visit_type for healthy pups
for (i in seq_along(data_long$visit)) {
  group_na <- data_long$group[i]
  if (is.na(group_na)) {
    next
  } else if (group_na == "HL") {
    data_long$visit_type[i] <- "monitor"
  }
}

# Export data for analysis
meta_vars <- c("pup_ID", "pup_no", "cage", "sex", "DOB_dam.dd.mm.yyyy", "DOB_pup.dd.mm.yyyy",
               "Exclude", "Exclude_Notes", "AEC", "Experiment_type", "Experiment_name",
               "challenge_stat", "challenge_dose", "challenge_visit", "treatment_dose",
               "treatment_visit", "sacrifice_type", "sacrifice_time.hh.mm", "sacrifice_date.dd.mm.yyyy",
               "group", "outcome", "sac_hpc")

measure_vars <- c("visit", "visit_type", variables, "visit_date", "hours_post_challenge")

data_export <- data_long[, c(meta_vars, measure_vars)]
data_export$expdate <- "300322"  # Adjust this value as needed

write.csv(data_export, "Result/intermediate_metadata.csv", row.names = FALSE)







