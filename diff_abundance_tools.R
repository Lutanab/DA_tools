library(phyloseq)
library(microbiome)
library(DESeq2)
library(ALDEx2)
library(ANCOMBC)
library(apeglm)
library(ashr)
library(mia)

#TODO Нормальный импорт зависимостей
project_dir <- "/home/gukin_eg/projects/DA_tools"
source(paste0(project_dir, "/common_meta_transform_class.R"))
source(paste0(project_dir, "/diff_abundance_class.R"))
source(paste0(project_dir, "/common_tools.R"))

##
# Подярок TARGET VS REF. Таргет - второй, референс - третий
# Разница это Target - Ref
##

## DESEQ 2

ps_to_deseq2 <- function(
    local_ps = ? phyloseq, design_formula = ? formula, ref_values = list() ? list)
  {
  count_data <- t(data.frame(local_ps@otu_table))
  meta_data <- data.frame(local_ps@sam_data)
  for(feature in names(ref_values)){
    meta_data[[feature]] <- set_reference_value(meta_data[[feature]], ref_values[[feature]])
  }
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = meta_data,
    design = design_formula)
  dds <- DESeq(dds)
  return(dds)
} ? DESeqDataSet



prepare_deseq2_feature_table <- function(
    dds = ? DESeqDataSet, USE_ADJUSTED_PVAL = TRUE ? logical, contrasts_table = NULL ? data.frame,
    lfc_shrink_type = "apeglm" ? character)
{
  if (USE_ADJUSTED_PVAL){
    pval_col <- "padj"
  }
  else {
    pval_col <- "pvalue"
  }
  lfc_table <- data.frame()
  pval_table <- data.frame()
  if(is.null(contrasts_table)){
    features <- resultsNames(dds)
    features <- features[2: length(features)]
    for(feature in features){
      tmp_res <- data.frame(lfcShrink(dds, coef = feature, type = lfc_shrink_type))
      lfc_table[row.names(tmp_res), feature] <- tmp_res[row.names(tmp_res), "log2FoldChange"]
      pval_table[row.names(tmp_res), feature] <- tmp_res[row.names(tmp_res), pval_col]
    }
  }
  else{
    for(column_name in c("condition", "ref", "target")){
      if(!column_name %in% colnames(contrasts_table)){
        stop(paste0("Custom contrasts table must contain '", column_name, "' column"))
      }
    }
    for(n_contrast in row.names(contrasts_table)){
      contrast_arg <- c(contrasts_table[n_contrast, "condition"],
                         contrasts_table[n_contrast, "target"],
                         contrasts_table[n_contrast, "ref"])
      feature <- paste0(contrast_arg[1], "_", contrast_arg[2], "_vs_", contrast_arg[3])
      tmp_res <- data.frame(lfcShrink(dds, contrast = contrast_arg, type = lfc_shrink_type))
      lfc_table[row.names(tmp_res), feature] <- tmp_res[row.names(tmp_res), "log2FoldChange"]
      pval_table[row.names(tmp_res), feature] <- tmp_res[row.names(tmp_res), pval_col]
    }
  }
  result <- DiffAbundanceComputationResult$new(
    lfc_table = as.data.frame(lfc_table),
    pval_table = as.data.frame(pval_table)
  )
  return(result)
} ? DiffAbundanceComputationResult

## ALDEX2

prepare_aldex2_feature_table <- function(
    local_ps = ? phyloseq, design_formula = ? formula,
    ref_values = ? list, USE_ADJUSTED_PVAL = TRUE ? logical)
{
  count_table <- t(data.frame(local_ps@otu_table))
  meta_data <- data.frame(local_ps@sam_data)
  
  mm <- get_custom_categorial_design_matrix(meta_data, ref_values, design_formula)
  feature_names <- colnames(mm)
  feature_names <- feature_names[2: length(feature_names)]
  
  result <- aldex.clr(count_table, mm, denom="all")
  result <- aldex.glm(result)
  
  if (USE_ADJUSTED_PVAL){
    pval_suffix <- ":pval.holm"
  }
  else {
    pval_suffix <- ":pval"
  }
  lfc_table <- data.frame()
  pval_table <- data.frame()
  
  lfc_table[row.names(result), feature_names] <- result[row.names(result), paste0(feature_names, ":Est")]
  pval_table[row.names(result), feature_names] <- result[row.names(result), paste0(feature_names, pval_suffix)]
  
  result <- DiffAbundanceComputationResult$new(
    lfc_table = as.data.frame(lfc_table),
    pval_table = as.data.frame(pval_table)
  )
  return(result)
} ? list



## ANCOMBC2

prepare_ancombc2_feature_table <- function(
    local_ps = ? phyloseq, design_formula_char = ? character,
    ref_values = ? list, p_adj_method = "holm" ? character)
  {
  tmp_sam_data <- local_ps@sam_data[,names(ref_values)]
  local_ps@sam_data <- tmp_sam_data
  
  ref_values_data <- set_reference_values(df=local_ps@sam_data, ref_values=ref_values)
  local_ps@sam_data <- ref_values_data[["df"]]
  old_column_names <- ref_values_data[["old_column_names"]]
  new_column_names <- ref_values_data[["new_column_names"]]
  
  tse <- makeTreeSummarizedExperimentFromPhyloseq(local_ps)
  result <- ancombc2(data = tse, assay_name = "counts",
                     fix_formula = design_formula_char,
                     p_adj_method = p_adj_method)
  result <- result$res
  
  lfc_table <- data.frame()
  pval_table <- data.frame()
  
  lfc_table[result[[1]], new_column_names] <- result[row.names(result), paste0("lfc_", old_column_names)]
  pval_table[result[[1]], new_column_names] <- result[row.names(result), paste0("p_", old_column_names)]
  
  result <- DiffAbundanceComputationResult$new(
    lfc_table = as.data.frame(lfc_table),
    pval_table = as.data.frame(pval_table)
  )
  return(result)
} ? list


## RAW GLM ###

prepare_raw_glm_feature_table <- function(
    local_ps = ? phyloseq, desing_formula = ? formula, ref_values = ? list,
    count_transform = "clr" ? character, USE_ADJUSTED_PVAL = TRUE ? logical)
{
  local_ps <- microbiome::transform(local_ps, count_transform)
  
  tmp_count <- data.frame(local_ps@otu_table)
  tmp_meta <- data.frame(local_ps@sam_data)
  tmp_meta <- set_reference_values(tmp_meta, ref_values)[["df"]]
  
  # Тот цирк, что ниже, нужен для того, чтобы норм работать с df c одной колонкой
  tmp_meta_x2 <- data.frame()
  tmp_meta_x2[row.names(tmp_meta), names(ref_values)] <- tmp_meta[row.names(tmp_meta), names(ref_values)]
  tmp_meta <- tmp_meta_x2
  
  # tmp_meta <- tmp_meta[,names(ref_values)]
  tmp_meta_x2 <- get_custom_categorial_design_matrix(tmp_meta, ref_values, desing_formula)
  tmp_meta <- data.frame()
  tmp_meta[row.names(tmp_meta_x2), colnames(tmp_meta_x2)[2:ncol(tmp_meta_x2)]] <- tmp_meta_x2[row.names(tmp_meta_x2), 2:ncol(tmp_meta_x2)]
  remove(tmp_meta_x2)
  
  feature_names <- colnames(tmp_meta)
  tmp_formula <- formula(paste0("target ~ ", paste(feature_names, collapse = " + ")))
  tmp_meta <- as.data.frame(tmp_meta)
  
  lfc_table <- data.frame()
  pval_table <- data.frame()
  
  for(current_taxon in colnames(tmp_count)){
    tmp_meta[['target']] <- tmp_count[[current_taxon]]
    tmp_model <- glm(tmp_formula, data = tmp_meta)
    tmp_model <- coef(summary(tmp_model))
    
    lfc_table[feature_names, current_taxon] <- tmp_model[feature_names, "Estimate"]
    pval_table[feature_names, current_taxon] <- tmp_model[feature_names, "Pr(>|t|)"]
    if(USE_ADJUSTED_PVAL){
      pval_table[feature_names, current_taxon] <- p.adjust(pval_table[feature_names, current_taxon])
    }
  }
  
  result <- DiffAbundanceComputationResult$new(
    lfc_table = as.data.frame(t(lfc_table)),
    pval_table = as.data.frame(t(pval_table))
  )
  return(result)
} ? list

## COMMON ##

set_reference_value <- function(
    categorial_vector, ref_value = ? character)
  {
  categorial_vector <- factor(categorial_vector)
  categorial_vector <- relevel(categorial_vector, ref_value)
  return(categorial_vector)
}



set_reference_values <- function(
    df = ? data.frame, ref_values = ? list)
  {
  old_column_names <- c()
  new_column_names <- c()
  
  for(feature in names(ref_values)){
    ref_value <- ref_values[[feature]]
    df[[feature]] <- set_reference_value(df[[feature]], ref_value)
    target_values <- levels(df[[feature]])
    if(length(target_values) >= 2){
      target_values <- target_values[2:length(target_values)]
    }
    else{
      target_values <- c()
    }
    
    current_old_column_names <- paste0(feature, target_values)
    old_column_names <- c(old_column_names,current_old_column_names)
    current_new_column_names <- paste(feature, target_values, "vs", ref_value, sep = "_")
    new_column_names <- c(new_column_names, current_new_column_names)
  }
  return(list(
    "df" = df,
    "old_column_names" = old_column_names,
    "new_column_names" = new_column_names
  ))
}



get_custom_categorial_design_matrix <- function(
    df = ? data.frame, ref_values = ? list, desing_formula = ? formula)
  {
  ref_values_data <- set_reference_values(df=df, ref_values=ref_values)
  df <- ref_values_data[["df"]]
  old_column_names <- ref_values_data[["old_column_names"]]
  new_column_names <- ref_values_data[["new_column_names"]]
  
  mm <- model.matrix(design_formula, df)
  mm <- rename_columns(mm, old_column_names, new_column_names)
  return(mm)
}


prepare_feature_table <- function(
    tool = ? character, local_ps = ? phyloseq, design_formula = ? formula,
    ref_values = ? list, USE_ADJUSTED_PVAL = TRUE ? logical){
  if(tool == "DESeq2"){
    dds <- ps_to_deseq2(local_ps, design_formula = design_formula, ref_values = ref_values)
    results <- prepare_deseq2_feature_table(dds, USE_ADJUSTED_PVAL = USE_ADJUSTED_PVAL)
  }
  else if(tool == "ALDEx2"){
    results <- prepare_aldex2_feature_table(local_ps, design_formula,
                                            ref_values, USE_ADJUSTED_PVAL=USE_ADJUSTED_PVAL)
  }
  else if(tool == "ANCOMBC"){
    if(USE_ADJUSTED_PVAL){p_adj_method <- "holm"}
    else{p_adj_method <- "none"}
    results <- prepare_ancombc2_feature_table(local_ps, as.character(design_formula)[2],
                                              ref_values, p_adj_method=p_adj_method)
  }
  else if(tool == "glm_raw"){
    results <- prepare_raw_glm_feature_table(local_ps, design_formula,
                                             ref_values, USE_ADJUSTED_PVAL=USE_ADJUSTED_PVAL)
  }
  else{
    stop("'tool' must be 'DESeq2', 'ALDEx2', 'ANCOMBC' or 'glm_raw'")
  }
  return(results)
} ? DiffAbundanceComputationResult



get_multiple_tools_feature_table <- function(
    local_ps = ? phyloseq, design_formula = ? formula,
    ref_values = ? list, tools = c("DESeq2", "ALDEx2", "ANCOMBC", "glm_raw")? character, USE_ADJUSTED_PVAL = TRUE ? logical)
  {
  common_results <- list()
  for(tool in tools){
    common_results[[tool]] <- prepare_feature_table(tool, local_ps, design_formula, ref_values)
  }
  common_results <- ExpandedDiffAbundanceComputationResult$new(common_results)
  common_results$conform_results()
  return(common_results)
} ? ExpandedDiffAbundanceComputationResult





