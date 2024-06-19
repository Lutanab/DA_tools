library(phyloseq)
library(microbiome)
library(dplyr)
library(tidyr)
library(randomcoloR)
library(types)
library(ggplot2)


project_dir <- "/home/gukin_eg/projects/DA_tools"
source(paste0(project_dir, "/common_tools.R"))
source(paste0(project_dir, "/common_meta_transform_class.R"))
source(paste0(project_dir, "/diff_abundance_class.R"))
source(paste0(project_dir, "/diff_abundance_tools.R"))
source(paste0(project_dir, "/heatmap_tools.R"))


get_time_table <- function(
    local_ps = ? phyloseq, time_point_column = ? character, subject_uid_column = ? character
){
  time_table_local <- data.frame()
  tmp_df <- data.frame(local_ps@sam_data)
  NA_SUBJECTS <- FALSE
  NA_TIMEPOINTS <- FALSE
  
  for(current_sample in row.names(tmp_df)){
    subj <- tmp_df[current_sample, subject_uid_column]
    if(is.na(subj) & !NA_SUBJECTS){
      warning("Subject_uid_column contain NA values. Dont consider such cases")
      NA_SUBJECTS <- TRUE
      next
    }
    
    timepoint <- tmp_df[current_sample, time_point_column]
    if(is.na(timepoint) & !NA_TIMEPOINTS){
      warning("Subject_uid_column contain NA values. Dont consider such cases")
      NA_TIMEPOINTS <- TRUE
      next
    }
    
    if(!is.null(time_table_local[subj, timepoint])){
      if(!is.na(time_table_local[subj, timepoint])){
        stop(paste0("Phyloseq contain samples with the same timepoint and subject_uid.\nSample: ", current_sample, " Timepoint: ", timepoint, " Subject: ", subj))
      }
    }
    time_table_local[subj, timepoint] <- current_sample
  }
  time_table_local <- time_table_local[, order(colnames(time_table_local))]
  return(time_table_local)
}


plot_sample_count_per_timepoint <- function(time_table = ? data.frame, theme = ? theme){
  tmp_df <- data.frame(colSums(!is.na(time_table)))
  colnames(tmp_df) <- c("N_subjects")
  tmp_df$timepoint <- row.names(tmp_df)
  
  g <- ggplot2::ggplot(tmp_df, aes(x=timepoint, y=N_subjects)) +
    ggplot2::geom_bar(stat = "identity")
  return(g)
}


plot_timepoint_arrower_pca <- function(
    local_ps =? phyloseq, time_table = ? data.frame, time_point_column = ? character, subject_uid_column = ? character,
    alpha = 0.35 ? numeric
){
  # Validators
  bad_timepoints <- setdiff(colnames(time_table), timepoint_order)
  if(length(bad_timepoints) != 0){
    stop(paste0(
      "There are timepoints, present in 'time_table', but absent in 'timepoint_order': ",
      paste0(bad_timepoints, collapse = ", ")
    ))
  }
  
  bad_timepoints <- setdiff(timepoint_order, colnames(time_table))
  if(length(bad_timepoints) != 0){
    stop(paste0(
      "There are timepoints, present in 'timepoint_order', but absent in 'time_table': ",
      paste0(bad_timepoints, collapse = ", ")
    ))
  }
  
  # Base PCA
  dist = phyloseq::distance(
    microbiome::transform(current_ps, "compositional"),
    method="bray")
  ordination = phyloseq::ordinate(current_ps, method="PCoA", distance=dist)
  tmp_df <- (phyloseq::plot_ordination(current_ps, ordination, color=time_point_column))$data
  tmp_df <- tmp_df[, c("Axis.1", "Axis.2", subject_uid_column, time_point_column)]
  tmp_df$time_point <- factor(tmp_df[[time_point_column]], levels = timepoint_order)
  tmp_df <- tmp_df[order(tmp_df[[subject_uid_column]], tmp_df[[time_point_column]]),]
  
  ## assembling arrow table
  for(i in 1:(nrow(tmp_df)-1)){
    if(tmp_df[i,subject_uid_column] == tmp_df[i+1,subject_uid_column]){
      tmp_df[i, 'x_end'] <- tmp_df[i+1, 'Axis.1']
      tmp_df[i, 'y_end'] <- tmp_df[i+1, 'Axis.2']
    }
  }
  
  # Plot
  g <- ggplot2::ggplot() + ggplot2::geom_segment(data=tmp_df, aes(x=Axis.1, y=Axis.2, xend=x_end, yend=y_end), alpha = alpha, arrow = arrow(length=unit(0.30,"cm"), angle=10, type = "closed"))
  g <- g + ggplot2::geom_point(data = tmp_df, size=2, aes(x=Axis.1, y=Axis.2, colour=tmp_df[[time_point_column]])) +
    ggplot2::guides(color = guide_legend(title = time_point_column))
  
  return(g)
}



get_glom_count_table <- function(local_ps = ? phyloseq, tax_rank = ? character)
  {
  if((!tax_rank %in% colnames(local_ps@tax_table)) & (tax_rank != "ASV")){
    stop(paste0(
      "Tax_rank '", tax_rank, "' is not present in phyloseq object"
    ))
  }
  
  if(tax_rank != "ASV"){local_ps <- tax_glom(local_ps, taxrank=tax_rank)}
  glom_count_table <- data.frame(local_ps@otu_table)
  if(!local_ps@otu_table@taxa_are_rows){glom_count_table <- data.frame(t(glom_count_table))}
  
  # Naming
  tmp_tax_table <- data.frame(local_ps@tax_table)
  n_tax_rank <- which(colnames(tmp_tax_table)==tax_rank)
  names_mapper <- data.frame(row.names=taxa_names(local_ps))
  if(n_tax_rank == 1){
    names_mapper[row.names(tmp_tax_table) ,"name"] <- paste0(colnames(tmp_tax_table)[1], ": ", tmp_tax_table[[1]])
  }else{
    names_mapper[row.names(tmp_tax_table) ,"name"] <- paste0(colnames(tmp_tax_table)[n_tax_rank-1], ": ", tmp_tax_table[[n_tax_rank-1]],
                                                             "; ", colnames(tmp_tax_table)[n_tax_rank], ": ", tmp_tax_table[[n_tax_rank]])
  }
  row.names(glom_count_table) <- names_mapper[row.names(glom_count_table), 'name']
  return(glom_count_table)
}


get_taxa_order <- function(glom_count_table = ? data.frame){
  taxa_order <- rowSums(glom_count_table, na.rm=TRUE)
  taxa_order <- taxa_order[order(taxa_order, decreasing=TRUE)]
  taxa_order <- c("Other", names(taxa_order))
  return(taxa_order)
}

generate_taxa_colors <- function(taxa_order = ? character){
  taxa_colors <- c('grey', randomcoloR::distinctColorPalette(length(taxa_order) - 1))
  names(taxa_colors) <- taxa_order
  return(taxa_colors)
}



plot_taxa_barchart_over_time <- function(
    glom_count_table = ? data.frame, time_table = ? data.frame, subjects = ? character,
    taxa_order = NULL ? character, taxa_colors = NULL ? character, max_taxa = 8 ? integer
){
  # fill taxa order and values if null
  if(is.null(taxa_order)){
    taxa_order <- get_taxa_order(glom_count_table)
  }
  if(is.null(taxa_colors)){
    taxa_colors <- generate_taxa_colors(taxa_order)
  }
  
  # Validators
  bad_timepoints <- setdiff(colnames(time_table), timepoint_order)
  if(length(bad_timepoints) != 0){
    stop(paste0(
      "There are timepoints, present in 'time_table', but absent in 'timepoint_order': ",
      paste0(bad_timepoints, collapse = ", ")
    ))
  }
  bad_timepoints <- setdiff(timepoint_order, colnames(time_table))
  if(length(bad_timepoints) != 0){
    stop(paste0(
      "There are timepoints, present in 'timepoint_order', but absent in 'time_table': ",
      paste0(bad_timepoints, collapse = ", ")
    ))
  }
  bad_subjects <- setdiff(subjects, row.names(time_table))
  if(length(bad_subjects) != 0){
    stop(paste0(
      "There are subjects, present in 'subjects' variable, but absent in 'time_table': ",
      paste0(bad_subjects, collapse = ", ")
    ))
  }
  
  if(
    length(setdiff(union(taxa_order, row.names(glom_count_table)), union(taxa_order, row.names(glom_count_table)))) != 0
    ){
    stop("'taxa_order' and 'glom_count_table' must have tha same taxa")
  }
  if(taxa_order[1] != "Other"){stop("1st element of taxa_order must be 'Other'")}
  if(
    length(setdiff(union(taxa_order, names(taxa_colors)), union(taxa_order, names(taxa_colors)))) != 0
  ){
    stop("'taxa_colors' must be a named charactre vector with names same as 'taxa_order'")
  }
  
  # Splitting by time points
  res_tmp_df <- data.frame(row.names=row.names(glom_count_table))
  timepoint_order <- colnames(time_table)
  for(timepoint in timepoint_order){
    tmp_samples <- time_table[subjects, timepoint]
    tmp_samples <- tmp_samples[!is.na(tmp_samples)]
    if(length(tmp_samples) == 0){next}
    tmp_samples <- rowSums((glom_count_table %>% dplyr::select(all_of(tmp_samples))), na.rm=TRUE)
    tmp_read_sum <- sum(tmp_samples)
    if(tmp_read_sum==0){next}
    res_tmp_df[names(tmp_samples), timepoint] <- as.numeric(tmp_samples)/tmp_read_sum
  }
  timepoint_order <- colnames(res_tmp_df)
  if(length(timepoint_order) <= 1){
    stop("Given subjects have eq/less 1 timepoints, which have summary !=0 abundance. Barchart over time is irrelevant")
  }
  
  #First top n taxa:
  tmp_samples <- time_table[subjects, timepoint_order]
  tmp_samples <- unlist(tmp_samples)
  tmp_samples <- tmp_samples[!is.na(tmp_samples)]
  
  top_n_taxa <-  rowSums(glom_count_table %>% dplyr::select(all_of(tmp_samples)))
  top_n_taxa <- top_n_taxa[order(top_n_taxa, decreasing=TRUE)]
  top_n_taxa <- names(top_n_taxa)
  top_n_taxa <- top_n_taxa[1:min(length(top_n_taxa), max_taxa)]
  
  bad_taxa <- setdiff(row.names(res_tmp_df), top_n_taxa)
  if(length(bad_taxa)!=0){
    bad_taxa_sums <- colSums(res_tmp_df[bad_taxa,])
    res_tmp_df["Other", names(bad_taxa_sums)] <- as.numeric(bad_taxa_sums)
    res_tmp_df <- res_tmp_df[c(top_n_taxa, "Other"),]
  }
  
  # df transformations for ggplot
  res_tmp_df_compressed <- data.frame()
  for(timepoint in colnames(res_tmp_df)){
    nums <- seq(nrow(res_tmp_df_compressed)+1, nrow(res_tmp_df_compressed)+nrow(res_tmp_df))
    res_tmp_df_compressed[nums, "taxon"] <- row.names(res_tmp_df)
    res_tmp_df_compressed[nums, "abundance"] <- res_tmp_df[[timepoint]]
    res_tmp_df_compressed[nums, "time_point"] <- as.numeric(which(timepoint_order == timepoint))
  }
  res_tmp_df_compressed$taxon <- factor(res_tmp_df_compressed$taxon, levels=taxa_order)
  
  #plot
  p <- ggplot2::ggplot(res_tmp_df_compressed, aes(x=time_point, y=abundance, fill = taxon)) +
    ggplot2::geom_area(position = 'stack') + 
    ggplot2::scale_fill_manual(values=taxa_colors) + 
    ggplot2::scale_x_continuous(name='time_point', breaks=(1:length(timepoint_order)), labels=(timepoint_order))
  return(p)
}

  
  
  
  
  
  
  
  
  
















































##############
# plot_taxa_barchart_over_time <- function(
#     local_ps = ? phyloseq, time_table = ? data.frame, subjects = ? character,
#     tax_rank = ? character, max_taxa = 8 ? integer
# ){
#   # Validators
#   bad_timepoints <- setdiff(colnames(time_table), timepoint_order)
#   if(length(bad_timepoints) != 0){
#     stop(paste0(
#       "There are timepoints, present in 'time_table', but absent in 'timepoint_order': ",
#       paste0(bad_timepoints, collapse = ", ")
#     ))
#   }
#   bad_timepoints <- setdiff(timepoint_order, colnames(time_table))
#   if(length(bad_timepoints) != 0){
#     stop(paste0(
#       "There are timepoints, present in 'timepoint_order', but absent in 'time_table': ",
#       paste0(bad_timepoints, collapse = ", ")
#     ))
#   }
#   if((!tax_rank %in% colnames(local_ps@tax_table)) & (tax_rank != "ASV")){
#     stop(paste0(
#       "Tax_rank '", tax_rank, "' is not present in phyloseq object"
#     ))
#   }
#   bad_subjects <- setdiff(subjects, row.names(time_table))
#   if(length(bad_subjects) != 0){
#     stop(paste0(
#       "There are subjects, present in 'subjects' variable, but absent in 'time_table': ",
#       paste0(bad_subjects, collapse = ", ")
#     ))
#   }
#   
#   if(tax_rank != "ASV"){local_ps <- tax_glom(local_ps, taxrank=tax_rank)}
#   tmp_df <- data.frame(local_ps@otu_table)
#   if(!local_ps@otu_table@taxa_are_rows){tmp_df <- data.frame(t(tmp_df))}
#   
#   # Naming
#   tmp_tax_table <- data.frame(local_ps@tax_table)
#   n_tax_rank <- which(colnames(tmp_tax_table)==tax_rank)
#   names_mapper <- data.frame(row.names=taxa_names(local_ps))
#   if(n_tax_rank == 1){
#     names_mapper[row.names(tmp_tax_table) ,"name"] <- paste0(colnames(tmp_tax_table)[1], ": ", tmp_tax_table[[1]])
#   }else{
#     names_mapper[row.names(tmp_tax_table) ,"name"] <- paste0(colnames(tmp_tax_table)[n_tax_rank-1], ": ", tmp_tax_table[[n_tax_rank-1]],
#                                                              "; ", colnames(tmp_tax_table)[n_tax_rank], ": ", tmp_tax_table[[n_tax_rank]])
#   }
#   row.names(tmp_df) <- names_mapper[row.names(tmp_df), 'name']
#   
#   
#   # Splitting by time points
#   res_tmp_df <- data.frame(row.names=row.names(tmp_df))
#   timepoint_order <- colnames(time_table)
#   for(timepoint in timepoint_order){
#     tmp_samples <- time_table[subjects, timepoint]
#     tmp_samples <- tmp_samples[!is.na(tmp_samples)]
#     if(length(tmp_samples) == 0){next}
#     tmp_samples <- rowSums((tmp_df %>% dplyr::select(all_of(tmp_samples))), na.rm=TRUE)
#     tmp_read_sum <- sum(tmp_samples)
#     if(tmp_read_sum==0){next}
#     res_tmp_df[names(tmp_samples), timepoint] <- as.numeric(tmp_samples)/tmp_read_sum
#   }
#   timepoint_order <- colnames(res_tmp_df)
#   if(length(timepoint_order) <= 1){
#     stop("Given subjects have eq/less 1 timepoints, which have summary !=0 abundance. Barchart over time is irrelevant")
#   }
#   
#   #First top n_taxa:
#   top_taxa <- rowSums(tmp_df, na.rm=TRUE)
#   top_taxa <- top_taxa[order(top_taxa, decreasing=TRUE)]
#   top_taxa <- top_taxa[1:min(length(top_taxa), max_taxa)]
#   top_taxa <- names(top_taxa)
#   
#   bad_taxa <- setdiff(row.names(res_tmp_df), top_taxa)
#   if(length(bad_taxa)!=0){
#     bad_taxa_sums <- colSums(res_tmp_df[bad_taxa,])
#     res_tmp_df["Other", names(bad_taxa_sums)] <- as.numeric(bad_taxa_sums)
#     res_tmp_df <- res_tmp_df[c(top_taxa, "Other"),]
#   }
#   
#   # df transformations for ggplot
#   res_tmp_df_compressed <- data.frame()
#   for(timepoint in colnames(res_tmp_df)){
#     nums <- seq(nrow(res_tmp_df_compressed)+1, nrow(res_tmp_df_compressed)+nrow(res_tmp_df))
#     res_tmp_df_compressed[nums, "taxon"] <- row.names(res_tmp_df)
#     res_tmp_df_compressed[nums, "abundance"] <- res_tmp_df[[timepoint]]
#     res_tmp_df_compressed[nums, "time_point"] <- as.numeric(which(timepoint_order == timepoint))
#   }
#   # res_tmp_df_compressed[['time_point']] <- as.numeric(res_tmp_df_compressed[['time_point']])
#   
#   res_tmp_df_compressed$taxon <- factor(res_tmp_df_compressed$taxon)
#   if(length(bad_taxa)!=0){
#     res_tmp_df_compressed$taxon <- factor(res_tmp_df_compressed$taxon)
#     res_tmp_df_compressed$taxon <- relevel(res_tmp_df_compressed$taxon, "Other")
#   }else{
#     res_tmp_df_compressed$taxon <- factor(res_tmp_df_compressed$taxon)
#   }
#   
#   if(length(bad_taxa)!=0){
#     barchart_colors <- c('grey', randomcoloR::distinctColorPalette(length(top_taxa) + 1))
#   }else{
#     barchart_colors <- randomcoloR::distinctColorPalette(length(top_taxa))
#   }
#   View(res_tmp_df_compressed)
#   View(res_tmp_df)
#   View(tmp_df)
#   
#   #plot
#   p <- ggplot2::ggplot(res_tmp_df_compressed, aes(x=time_point, y=abundance, fill = taxon)) +
#     ggplot2::geom_area(position = 'stack') + 
#     ggplot2::scale_fill_manual(values=barchart_colors) + 
#     ggplot2::scale_x_continuous(name='time_point', breaks=(1:length(timepoint_order)), labels=(timepoint_order))
#   return(p)
# }
