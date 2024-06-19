#Bioconductor packages:
library(dada2)
library(decontam)
library(ggplot2)
library(microbiome)
library(phyloseq)

#Base R packages:
library(stringr)
library(readr)


## COMMON TOOLS
drop_columns <- function (df, columns) {
  return(df[, !(colnames(df) %in% columns)])
}

rename_column <- function(df, old_name, new_name) {
  colnames(df)[colnames(df) == old_name] <- new_name
  return(df)
}

replace_na <- function (df, to_replace) {
  df[is.na(df)] <- to_replace
  return(df)
}

rename_columns <- function(df, old_names, new_names) {
  if (!is.character(old_names)) {
    stop("old_names must be a character atomic vector")
  }
  if (!is.character(new_names)) {
    stop("new_names must be a character atomic vector")
  }
  if (length(old_names) != length(new_names)) {
    stop("old_names and new_names must have same length")
  }
  if (length(old_names) != 0) {
    for (i in 1:length(old_names)) {
      old_name <- old_names[i]
      new_name <- new_names[i]
      colnames(df)[colnames(df) == old_name] <- new_name
    }
  }
  return(df)
}

rename_rows <- function(df, old_names, new_names) {
  if (!is.character(old_names)) {
    stop("old_names must be a character atomic vector")
  }
  if (!is.character(new_names)) {
    stop("new_names must be a character atomic vector")
  }
  if (length(old_names) != length(new_names)) {
    stop("old_names and new_names must have same length")
  }
  if (length(old_names) != 0) {
    for (i in 1:length(old_names)) {
      old_name <- old_names[i]
      new_name <- new_names[i]
      rownames(df)[rownames(df) == old_name] <- new_name
    }
  }
  return(df)
}

replace_values <- function(char_vector, old_names, new_names) {
  if (!is.atomic(old_names)) {
    stop("old_names must be a atomic vector")
  }
  if (!is.atomic(new_names)) {
    stop("new_names must be a atomic vector")
  }
  if (!is.atomic(char_vector)) {
    stop("char_vector must be a atomic vector")
  }
  if (!is.atomic(new_names)) {
    stop("new_names must be a atomic vector")
  }
  if (length(old_names) != length(new_names)) {
    stop("old_names and new_names must have same length")
  }
  if (length(old_names) != 0) {
    for (i in 1:length(old_names)) {
      old_name <- old_names[i]
      new_name <- new_names[i]
      char_vector[char_vector == old_name] <- new_name
    }
  }
  return(char_vector)
}


are_vals_in_df <- function(df = ? data.frame, vals){
  bool_df <- FALSE
  for(val in vals){
    bool_df <- bool_df | (df == val)
  }
  return(bool_df)
}


check_df_consistency <- function(df1 = ? data.frame, df2 = ? data.frame){
  return(
    all(row.names(df1) == row.names(df2)) & all(colnames(df1) == colnames(df2))
  )
}





