# 2020-August-07
# Created by Danny Harbeson
# Adapted 2021-Mar by Adrien Eynaud
### Helper functions for Metabolic chamber data analysis


# Imports all csvs in a given folder as a list of dataframes then assigns names according to filename. Folder must be formated as "./FOLDER" (also acceptable: "./FOLDER/"). Default behavior excludes 'datetime' column.
data_import <- function(FOLDER, INCLUDE_DATETIME = F){
  require(tidyverse)
  # Generate vector of all .csv files in folder, with paths
  CSVS = list.files(FOLDER, full.names = TRUE)[grepl("csv",list.files(FOLDER, full.names = TRUE))] 
  # Call names of files and then remove the '.csv' string, will be used to name dfs in list
  FILENAMES = list.files(FOLDER)[grepl("csv",list.files(FOLDER))] 
  for(i in 1:length(FILENAMES)){
    FILENAMES[i] <- substr(FILENAMES[i], 0, gregexpr(".csv", FILENAMES)[[i]]-1)
  }
  # Read csvs, name as described
  out = lapply(CSVS, function(x) read.csv(x, stringsAsFactors = F))
  names(out) = FILENAMES
  if(INCLUDE_DATETIME == T){
    return(out)
  } else {
    out <- lapply(out, function(x) x <- x[,-1])
    return(out)
  }
}

# Simple function to run substr backwards
right = function(text, num_char) {
  substr(text, nchar(text) - (num_char-1), nchar(text))
}

# Merge list of identical dataframes together into singular dataframe. Will only work if all columns are identical
list_to_dataframe <- function(lst){
  len = length(lst)
  out = rbind(lst[[1]], lst[[2]])
  
  i = 3
  while(i <= len){
    out <- rbind(out, lst[[i]])
    i <- i+1
  }
  return(out)
}
