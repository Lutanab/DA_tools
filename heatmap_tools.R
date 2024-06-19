library(phyloseq)
library(microbiome)
library(pheatmap)
library(ggplot2)

#TODO Нормальный импорт зависимостей
source("/home/gukin_eg/projects/TEETH_2024/scripts/common_tools.R")
source("/home/gukin_eg/projects/TEETH_2024/scripts/common_meta_transform_class.R")
source("/home/gukin_eg/projects/TEETH_2024/scripts/diff_abundance_class.R")
source("/home/gukin_eg/projects/TEETH_2024/scripts/diff_abundance_tools.R")


## HEATMAP HELPERS  ##

get_taxa_names_aligner <- function(local_ps = ? phyloseq){
  names_aligner <- data.frame(row.names = row.names(local_ps@tax_table))
  tmp_tax_df <- data.frame(local_ps@tax_table)
  names_aligner$name <- paste0(row.names(tmp_tax_df), '__', tmp_tax_df$Family,'.g__' ,gsub('/', '', gsub("-", "", tmp_tax_df$Genus)),'__', tmp_tax_df$Species)
  return(names_aligner)
} ? data.frame


get_significant_taxa <- function(table = ? data.frame, unsignificance_label = ""){
  relevant_taxa <- FALSE
  for(feature in colnames(table)){
    relevant_taxa <- relevant_taxa | (table[[feature]] != unsignificance_label)
  }
  return(relevant_taxa)
}


filter_result_by_values <- function(diff_abundance_data = ? DiffAbundanceComputationResult, columns_to_filter = NA_character_ ? character,
                                    filter_lfc_less_then = +Inf ? numeric, filter_lfc_more_then = -Inf ? numeric){
  if(is.na(columns_to_filter)){
    columns_to_filter <- diff_abundance_data$columns
  }
  if(any(!columns_to_filter %in% diff_abundance_data$columns)){
    unmatched_columns <- columns_to_filter[!columns_to_filter %in% diff_abundance_data$columns]
    stop(
      paste0(paste0(unmatched_columns, collapse=", "), ": these columns are not present among idff abundance columns")
      )
  }
  appropriate_taxa <- FALSE
  for(column in columns_to_filter){
    current_appropriate_taxa <- (diff_abundance_data$lfc_table[[column]] <= filter_lfc_less_then) | (diff_abundance_data$lfc_table[[column]] >= filter_lfc_more_then)
    current_appropriate_taxa[is.na(current_appropriate_taxa)] <- FALSE
    appropriate_taxa <- appropriate_taxa | current_appropriate_taxa
  }
  return(appropriate_taxa)
}

preprocess_heatmap_data <- function(
    local_ps = ? phyloseq, diff_abundance_data = ? DiffAbundanceComputationResult,
    unsignificance_label = "" ? character, REMOVE_IRRELEVANT_TAXA = TRUE ? logical, FULL_NAMES = TRUE ? logical,
    columns_to_filter = NA_character_ ? character, filter_lfc_less_then = +Inf ? numeric, filter_lfc_more_then = -Inf ? numeric
    ){
  if(FULL_NAMES){
    names_aligner <- get_taxa_names_aligner(local_ps)
    diff_abundance_data$rename_taxa(
      names_aligner[row.names(diff_abundance_data$lfc_table), 'name']
    )
  }
  if(REMOVE_IRRELEVANT_TAXA){
    diff_abundance_data$prune_taxa(
      get_significant_taxa(diff_abundance_data$significance_labels, unsignificance_label=unsignificance_label)
    )
  }
  if((filter_lfc_more_then != -Inf) | (filter_lfc_less_then != +Inf)){
    diff_abundance_data$prune_taxa(
      filter_result_by_values(diff_abundance_data, columns_to_filter = columns_to_filter, 
                              filter_lfc_less_then = filter_lfc_less_then, filter_lfc_more_then = filter_lfc_more_then)
    )
  }
  return(diff_abundance_data)
}


get_heatmap_breaks <- function(heatmap_values = ? data.frame, color_palette = ? character){
  palette_length <- length(color_palette)
  relevant_values <- !is.na(heatmap_values)
  
  custom_breaks <- c(seq(min(heatmap_values[relevant_values]), 0, length.out=ceiling(palette_length/2) + 1),
                     seq(max(heatmap_values[relevant_values])/palette_length, max(heatmap_values[relevant_values]), length.out=floor(palette_length/2)))
  return(custom_breaks)
}


## For common taxa abundance

prepare_taxa_abundance_heatmap_table <- function(
    local_ps = ? phyloseq, dividing_feature = ? character, features_of_interest = ? character,
    FULL_FEATURE_NAMES = TRUE ? logical, CLUSSTERING_WITHIN_DIVIDING_FEATURE = FALSE ? logical
)
  {
  names_aligner <- get_taxa_names_aligner(local_ps)
  
  count_mtrx <- data.frame(local_ps@otu_table)
  sample_meta <- as(sample_data(local_ps), 'data.frame')
  features_of_interest <- c(dividing_feature, features_of_interest)
  sample_meta <- sample_meta[, features_of_interest]
  
  feature_vals <- unique(sample_meta[[dividing_feature]])
  feature_vals <- feature_vals[order(feature_vals)]
  
  if (CLUSSTERING_WITHIN_DIVIDING_FEATURE){
    sub_matrixes <- list()
    i = 1
    breaks <- c()
    for (val in feature_vals) {
      sub_counts <- count_mtrx[row.names(sample_meta[sample_meta[,dividing_feature] == val,]),]
      sub_meta <- sample_meta[row.names(sub_counts),]
      #h clustering
      hc2  <- stats::hclust(stats::dist(sub_counts, method='euclidean'), method = 'ward.D2')
      sub_counts <- sub_counts[hc2$order,]
      sub_matrixes[[i]] <- sub_counts
      if (is.null(breaks[[i-1]])) { prev <- 0 }
      else { prev <- breaks[[i-1]] }
      breaks <- append(breaks, prev + dim(sub_counts)[1])
      i = i + 1
    }
    ordered_counts <- do.call(rbind, sub_matrixes)
  } else {ordered_counts <- count_mtrx}
  
  # ==============================================================================
  sample_meta <- sample_meta[row.names(ordered_counts),]
  #taxa_meta <- decontam_table[colnames(ordered_counts),]
  
  if (FULL_FEATURE_NAMES){
    colnames(ordered_counts) <- names_aligner[colnames(ordered_counts), 'name']
    #row.names(taxa_meta) <- names_aligner[row.names(taxa_meta), 'name']
  }
  
  return(list(
    "heat_matrix" = t(ordered_counts),
    "sample_meta" = sample_meta
  ))
}


## PLOT FUNCTIONS ##
plot_one_tool_on_multiple_features_heatmap <- function(
    local_ps = ? phyloseq, diff_abundance_data = ? DiffAbundanceComputationResult,
    filepath = ? character, relevant_colnames = ? character, relevant_colnames_rename_to = ? character,
    color_palette = colorRampPalette((c("red", "white", "green")))(100) ? character, cluster_taxa = FALSE ? logical, cutree_taxa = 1 ? integer,
    columns_to_lfc_filter = NA_character_, filter_lfc_less_then = +Inf ? numeric, filter_lfc_more_then = -Inf ? numeric,
    REMOVE_IRRELEVANT_TAXA = TRUE ? logical, FULL_NAMES = TRUE ? logical
    ){
  diff_abundance_data <- diff_abundance_data$clone()
  diff_abundance_data$lfc_table <- diff_abundance_data$lfc_table[,relevant_colnames]
  diff_abundance_data$pval_table <- NULL
  if(is.null(diff_abundance_data$significance_labels)){
    stop("Please set pval stars table")
  }
  diff_abundance_data$prune_columns(relevant_colnames)
  diff_abundance_data$rename_columns(relevant_colnames, relevant_colnames_rename_to)
  preprocess_heatmap_data(local_ps, diff_abundance_data, 
                          REMOVE_IRRELEVANT_TAXA = REMOVE_IRRELEVANT_TAXA, FULL_NAMES = FULL_NAMES,
                          columns_to_filter=columns_to_lfc_filter, filter_lfc_less_then=filter_lfc_less_then, filter_lfc_more_then=filter_lfc_more_then)
  
  
  if(length(diff_abundance_data$taxa) == 0){
    print("0 taxa were found after all processings (if was)")
    return(NULL)
  }
  breaks <- get_heatmap_breaks(diff_abundance_data$lfc_table, color_palette)
  p1 <- pheatmap(diff_abundance_data$lfc_table,
                 color = color_palette, breaks = breaks,
                 cluster_rows = cluster_taxa, cutree_rows = cutree_taxa,
                 cluster_cols = F,
                 show_rownames = TRUE,
                 show_colnames = TRUE,
                 filename=filepath,
                 cellwidth = 30,
                 cellheight = 30,
                 border_color=NA,
                 display_numbers = diff_abundance_data$significance_labels,
                 fontsize_number = 10
  )
  return(p1)
}

plot_multiple_tools_on_one_feature_heatmap <- function(local_ps = ? phyloseq, common_results = ? ExpandedDiffAbundanceComputationResult,
                                                 target_feature = ? character, heatmap_title = ? character,
                                                 filename = ? character,  REMOVE_IRRELEVANT_TAXA = TRUE ? logical, FULL_NAMES = TRUE ? logical,
                                                 color_palette = colorRampPalette((c("red", "white", "green")))(100) ? character,
                                                 cluster_taxa = FALSE ? logical, cutree_taxa = 1 ? integer,
                                                 columns_to_lfc_filter = NA_character_, filter_lfc_less_then = +Inf ? numeric, filter_lfc_more_then = -Inf ? numeric
){
  common_lfc_table <- data.frame()
  common_pval_stars_table <- data.frame()
  
  for(tool in common_results$tools){
    tmp_lfc_table <- as.data.frame(common_results[[tool]]$lfc_table)
    tmp_pval_stars_table <- as.data.frame(common_results[[tool]]$significance_labels)
    
    avalible_taxa <- row.names(tmp_lfc_table)
    common_lfc_table[avalible_taxa, tool] <- tmp_lfc_table[avalible_taxa, target_feature]
    common_pval_stars_table[avalible_taxa, tool] <- tmp_pval_stars_table[avalible_taxa, target_feature]
  }
  common_pval_stars_table <- replace_na(common_pval_stars_table, "")
  result <- DiffAbundanceComputationResult$from_significance_labels(common_lfc_table, common_pval_stars_table)
  
  result <- preprocess_heatmap_data(local_ps, result, REMOVE_IRRELEVANT_TAXA = REMOVE_IRRELEVANT_TAXA, FULL_NAMES = FULL_NAMES,
                          columns_to_filter=columns_to_lfc_filter, filter_lfc_less_then=filter_lfc_less_then, filter_lfc_more_then=filter_lfc_more_then)
  
  if(length(result$taxa) == 0){
    print("0 taxa were found after all processings (if was)")
    return(NULL)
  }
  breaks <- get_heatmap_breaks(result$lfc_table, color_palette)
  p1 <- pheatmap(result$lfc_table,
                 color = color_palette, breaks = breaks,
                 cluster_rows = cluster_taxa, cutree_rows = cutree_taxa,
                 cluster_cols = F,
                 show_rownames = TRUE,
                 show_colnames = TRUE,
                 filename=filename,
                 main=heatmap_title,
                 cellwidth = 30,
                 cellheight = 30,
                 border_color=NA,
                 display_numbers = result$significance_labels,
                 fontsize_number = 10
  )
  return(p1)
}


plot_multiple_tools_on_multiple_features_heatmap <- function(
    local_ps = ? phyloseq, common_results = ? ExpandedDiffAbundanceComputationResult, filename = ? character, relevant_meta_features = ? character, 
    relevant_meta_features_rename_to = ? character, n_stars_detection = 1 ? integer, REMOVE_IRRELEVANT_TAXA = TRUE ? logical, FULL_NAMES = TRUE ? logical,
    EMPHASIZE_INCONSISTENT_CASES = TRUE ? logical, heatmap_title = "Number of tools that give significant LFC" ? character,
    color_palette = colorRampPalette((c("red", "white", "green")))(100) ? character, cluster_taxa = FALSE ? logical, cutree_taxa = 1 ? integer,
    filter_lfc_less_then = +Inf ? numeric, filter_lfc_more_then = -Inf ? numeric
){
  if((length(n_stars_detection) != 1) | !(n_stars_detection %in% c(1, 2, 3, 4))){
    stop("'n_stars_detection' must be 1, 2, 3 or 4")
  }
  
  tools_count_table <- data.frame()
  tools_count_table[common_results$taxa, common_results$columns] <- 0
  avg_lfc_table <- data.frame()
  avg_lfc_table[common_results$taxa, common_results$columns] <- 0
  if(EMPHASIZE_INCONSISTENT_CASES){
    lfc_all_plus_sign <- data.frame()
    lfc_all_plus_sign[common_results$taxa, common_results$columns] <- TRUE
    lfc_all_minus_sign <- data.frame()
    lfc_all_minus_sign[common_results$taxa, common_results$columns] <- TRUE
  }
  
  for(tool in common_results$tools){
    tmp_pval_stars_table <- as.data.frame(common_results[[tool]]$significance_labels)
    tmp_lfc_table <- as.data.frame(common_results[[tool]]$lfc_table)
    if(is.null(tmp_pval_stars_table)){
      stop(paste0("Please define significance labels. Tool: ", tool))
    }
    
    if(n_stars_detection == 1){
      is_enough_stars <- (tmp_pval_stars_table != "")
    } else if(n_stars_detection == 2){
      is_enough_stars <- !(are_vals_in_df(tmp_pval_stars_table, c("", "*")))
    } else if(n_stars_detection == 3){
      is_enough_stars <- (are_vals_in_df(tmp_pval_stars_table, c("**\n*", "**\n**")))
    } else {
      is_enough_stars <- (tmp_pval_stars_table == "**\n**")
    }
    
    in_lfc_filter_range <- (tmp_lfc_table <= filter_lfc_less_then) | (tmp_lfc_table >= filter_lfc_more_then)
    
    appropriate_lfcs <- is_enough_stars & in_lfc_filter_range
    
    tools_count_table <- tools_count_table + as.integer(appropriate_lfcs)
    avg_lfc_table[appropriate_lfcs] <- avg_lfc_table[appropriate_lfcs] + tmp_lfc_table[appropriate_lfcs]
    if(EMPHASIZE_INCONSISTENT_CASES){
      lfc_all_plus_sign[appropriate_lfcs] <- (lfc_all_plus_sign[appropriate_lfcs]) & (tmp_lfc_table >= 0)[appropriate_lfcs]
      lfc_all_minus_sign[appropriate_lfcs] <- (lfc_all_minus_sign[appropriate_lfcs]) & (tmp_lfc_table < 0)[appropriate_lfcs]
    }
  }
  
  significant_cases <- tools_count_table != 0
  avg_lfc_table[significant_cases] <- avg_lfc_table[significant_cases] / tools_count_table[significant_cases]
  
  if(EMPHASIZE_INCONSISTENT_CASES){
    inconsistent_cases <- !((lfc_all_minus_sign) | (lfc_all_plus_sign))
    avg_lfc_table[significant_cases & inconsistent_cases] <- NA
  }
  
  for(column in colnames(tools_count_table)){
    tools_count_table[[column]] <- as.character(tools_count_table[[column]])
  }
  
  result <- DiffAbundanceComputationResult$from_significance_labels(avg_lfc_table, tools_count_table)
  preprocess_heatmap_data(local_ps, result, unsignificance_label = "0", REMOVE_IRRELEVANT_TAXA = REMOVE_IRRELEVANT_TAXA, FULL_NAMES = FULL_NAMES)
  result$prune_columns(relevant_meta_features)
  result$rename_columns(relevant_meta_features, relevant_meta_features_rename_to)
  
  if(length(result$taxa) == 0){
    print("0 taxa were found after all processings (if was)")
    return(NULL)
  }
  breaks <- get_heatmap_breaks(result$lfc_table, color_palette)
  p1 <- pheatmap::pheatmap(result$lfc_table,
                           color = color_palette, breaks = breaks,
                           cluster_rows = cluster_taxa, cutree_rows = cutree_taxa,
                           cluster_cols = F,
                           show_rownames = TRUE,
                           show_colnames = TRUE,
                           filename=filename,
                           main = heatmap_title,
                           cellwidth = 30,
                           cellheight = 30,
                           border_color=NA,
                           display_numbers = result$significance_labels,
                           fontsize_number = 10
  )
  return(p1)
}


