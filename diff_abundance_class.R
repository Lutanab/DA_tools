library(types)
library(R6)



DiffAbundanceComputationResult <- R6Class(
  "DiffAbundanceComputationResult",
  public = list(
    lfc_table = NULL ? data.frame,
    pval_table = NULL ? data.frame,
    significance_labels = NULL ? data.frame,
    initialize = function(lfc_table = ? data.frame, pval_table = ? data.frame){
      self$lfc_table <- lfc_table
      self$pval_table <- pval_table
      if(!(is.null(lfc_table) & is.null(pval_table))){
        if(!check_df_consistency(lfc_table, pval_table)){
          stop("lfc and pval tables have different rows/columns (or in differernt order)")
        }
      }
    },
    pval_to_stars = function(one_star_treshold = 0.05 ? numeric, two_stars_treshold = 0.01 ? numeric,
                              three_stars_treshold = 0.005 ? numeric, four_stars_treshold = 0.001 ? numeric)
      {
      stars_table <- data.frame(self$pval_table)
      stars_table[,] <- ""
      stars_table[self$pval_table <= one_star_treshold] <- "*"
      stars_table[self$pval_table <= two_stars_treshold] <- "**"
      stars_table[self$pval_table <= three_stars_treshold] <- "**\n*"
      stars_table[self$pval_table <= four_stars_treshold] <- "**\n**"
      self$significance_labels <- stars_table
    },
    prune_taxa = function(condition = ? logical){
      self$lfc_table <- self$lfc_table[condition,]
      if(!is.null(self$pval_table)){
        self$pval_table <- self$pval_table[condition,]
      }
      if(!is.null(self$significance_labels)){
        self$significance_labels <- self$significance_labels[condition,]
      }
    },
    prune_columns = function(relevant_columns = ? character){
      self$lfc_table <- self$lfc_table[,relevant_columns]
      if(!is.null(self$pval_table)){
        self$pval_table <- self$pval_table[,relevant_columns]
      }
      if(!is.null(self$significance_labels)){
        self$significance_labels <- self$significance_labels[,relevant_columns]
      }
    },
    rename_taxa = function(to_rename = ? character){
      row.names(self$lfc_table) <- to_rename
      if(!is.null(self$pval_table)){
        row.names(self$pval_table) <- to_rename
      }
      if(!is.null(self$significance_labels)){
        row.names(self$significance_labels) <- to_rename
      }
    },
    rename_columns = function(old_columns_names = ? character, new_columns_names = ? character){
      self$lfc_table <- rename_columns(
        self$lfc_table, old_columns_names, new_columns_names
        )
      if(!is.null(self$pval_table)){
        self$pval_table <- rename_columns(
          self$pval_table, old_columns_names, new_columns_names
        )
      }
      if(!is.null(self$significance_labels)){
        self$significance_labels <- rename_columns(
          self$significance_labels, old_columns_names, new_columns_names
        )
      }
    }
    # get_columns = function(){return(colnames(self$lfc_table))},
    # get_taxa = function(){return(row.names(self$lfc_table))}
  ),
  active = list(
    columns = function(){return(colnames(self$lfc_table))},
    taxa = function(){return(row.names(self$lfc_table))}
  )
)

DiffAbundanceComputationResult$from_significance_labels <- function(
    lfc_table = ? data.frame, significance_labels = ? data.frame
){
  obj <- DiffAbundanceComputationResult$new(NULL, NULL)
  obj$lfc_table <- lfc_table
  obj$significance_labels <- significance_labels
  if(!check_df_consistency(lfc_table, significance_labels)){
    stop("lfc and significance labels tables have different rows/columns (or in differernt order)")
  }
  return(obj)
}



ExpandedDiffAbundanceComputationResult <- R6Class(
  "ExpandedDiffAbundanceComputationResult",
  lock_objects = FALSE,
  public = list(
    tools = NA_character_ ? character,
    conform_results = function(){
      all_taxa <- c()
      all_columns <- c()
      
      for(tool in self$tools){
        all_taxa <- union(all_taxa, self[[tool]]$taxa)
        all_columns <- union(all_columns, self[[tool]]$columns)
      }
      
      for(tool in self$tools){
        absent_taxa <- setdiff(all_taxa, self[[tool]]$taxa)
        absent_columns <- setdiff(all_columns, self[[tool]]$columns)
        
        self[[tool]]$lfc_table[absent_taxa, ] <- NA_real_
        self[[tool]]$lfc_table[, absent_columns] <- NA_real_
        if(!is.null(self[[tool]]$pval_table)){
          self[[tool]]$pval_table[absent_taxa, ] <- NA_real_
          self[[tool]]$pval_table[, absent_columns] <- NA_real_
        }
        if(!is.null(self[[tool]]$significance_labels)){
          self[[tool]]$significance_labels[absent_taxa, ] <- NA_real_
          self[[tool]]$significance_labels[, absent_columns] <- NA_real_
        }
      }
    },
    initialize = function(result_list = ? list){
      if(length(result_list) == 0){
        stop("Input must be not empty list")
      }
      self$tools <- names(result_list)
      
      for(key in names(result_list)){
        if(!is(result_list[[key]], "DiffAbundanceComputationResult")){
          stop("Input must be a list like 'tool name' = DiffAbundanceComputationResult object")
        }
        self[[key]] <- result_list[[key]]
        # self$conform_results()
        # TODO Понять чому они по отдельности работают хорошо, а когда запускаю из под инита пздц
      }
    }
  ),
  active = list(
    columns = function(){return(self[[self$tools[1]]]$columns)},
    taxa = function(){return(self[[self$tools[1]]]$taxa)}
  )
)
