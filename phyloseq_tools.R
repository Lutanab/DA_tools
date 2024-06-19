#Bioconductor packages:
library(dada2)
library(decontam)
library(ggplot2)
library(microbiome)
library(phyloseq)

#Base R packages:
library(stringr)
library(readr)

library(networkD3)
library(htmlwidgets)
library(htmltools)

project_dir <- "/home/gukin_eg/projects/DA_tools"
source(paste0(project_dir, "/common_tools.R"))

prepare_count_table <- function (ps_obj, taxa_name_func) {
  otu_matrix <- as(otu_table(ps_obj), 'matrix')
  taxa_matrix <- as(tax_table(ps_obj), 'matrix')
  taxa_matrix <- cbind(ASV=rownames(taxa_matrix), taxa_matrix)  
  taxa_matrix_good_names <- apply(taxa_matrix, MARGIN=1, taxa_name_func)
  # ps_meta <- as(sample_data(ps_obj), 'matrix')
  colnames(otu_matrix) <- taxa_matrix_good_names
  return(otu_matrix)
}

add_metadata_to_ps <- function (ps_obj, meta_df) {
  additional_meta_fields <- colnames(meta_df)
  for (sample_name in rownames(ps_obj@sam_data)) {
    if (sample_name %in% rownames(meta_df)) {
      ps_obj@sam_data[sample_name, additional_meta_fields] <- meta_df[sample_name, additional_meta_fields]
    }
  }
  return(ps_obj)
}

remove_metadata_columns <- function (ps_obj, columns) {
  ps_obj@sam_data <- drop_columns(ps_obj@sam_data, columns)
  return(ps_obj)
}

rename_metadata_columns <- function (ps_obj, old_names, new_names) {
  ps_obj@sam_data <- rename_columns(ps_obj@sam_data, old_names, new_names)
  return(ps_obj)
}


get_decotam_table <- function(ps, dividing_column_name, concentration_column_name, use_labels=NULL) {
  if (!class(ps) == "phyloseq"){
    stop("ps must be an phyloseq object")
  }
  if (!dividing_column_name %in% colnames(ps@sam_data)){
    stop(paste0("Column '", dividing_column_name, "' must be in sample data table"))
  }
  if (!concentration_column_name %in% colnames(ps@sam_data)){
    stop(paste0("Column '", concentration_column_name, "' must be in sample data table"))
  }
  if (! is.numeric(ps@sam_data[[concentration_column_name]])){
    warning(paste0("Column '", concentration_column_name, "' (as concentration column) must be numeric"))
    ps@sam_data[[concentration_column_name]] <- as.numeric(ps@sam_data[[concentration_column_name]])
  }
  if (! is.character(ps@sam_data[[dividing_column_name]])){
    stop(paste0("Column '", dividing_column_name, "' (as column of dividing features) must be character"))
  }
  
  label_values <- NULL
  possible_label_values <- unique(ps@sam_data[[dividing_column_name]])
  if (is.null(use_labels)){
    label_values <- possible_label_values
  }
  else if(!is.character(use_labels)){
    stop(paste0("'use_labels' must be a character atomic vector"))
  }
  else {
    for (label in use_labels){
      if (label %in% possible_label_values){
        label_values <- c(label_values, label)
      }
      else{
        warning(paste0("Value '", label, "' not fount among column '", dividing_column_name, "'s values"))
      }
    }
  }
  
  # print(label_values)
  decotam_table <- data.frame(row.names=taxa_names(ps))
  for (label in label_values){
    tmp_ps <- ps
    tmp_ps@sam_data$subsetting_condition <- (tmp_ps@sam_data[[dividing_column_name]] == label)
    tmp_ps@sam_data$subsetting_condition <- (tmp_ps@sam_data$subsetting_condition) & !(is.na(tmp_ps@sam_data[[concentration_column_name]]))
    tmp_ps@sam_data$subsetting_condition <- (tmp_ps@sam_data$subsetting_condition) & !(is.null(tmp_ps@sam_data[[concentration_column_name]]))
    if (!any(tmp_ps@sam_data$subsetting_condition)){
      next
    }
    tmp_ps <- subset_samples(tmp_ps, subsetting_condition)
    tmp_ps@sam_data[[concentration_column_name]] <- as.numeric(tmp_ps@sam_data[[concentration_column_name]])
    contamdf.freq <-isContaminant(tmp_ps, method="freq", conc=concentration_column_name, normalize=T)
    decotam_table[[paste0(label, "_is_contam")]] <- contamdf.freq$contaminant
    print(label)
  }
  return(decotam_table)
}

## SANKEY_DIAGRAM
plot_sankey_diagram <- function(links,
                                fontSize = 15
) {
  nodes <- data.frame(
    name=unique(c(
      as.character(links$source), 
      as.character(links$target)
    ))
  )
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  p <- sankeyNetwork(Links = links, Nodes = nodes,
                     Source = "IDsource", Target = "IDtarget",
                     Value = "value", NodeID = "name", 
                     sinksRight=FALSE, fontSize=fontSize)
  return(p)
}