library(phyloseq)
library(microbiome)
library(types)
library(R6)

CommonMetaTransform <- R6Class('CommonMetaTransform',
                               public = list(
                                 min_reads=0,
                                 detection=NA_integer_,
                                 prevalence=0.0,
                                 tax_glom_level=NA_character_,
                                 count_transformation=NA_character_,
                                 initialize = function(min_reads = 0 ? integer, detection = NA_integer_ ? integer,
                                                       prevalence = 0.0 ? numeric, tax_glom_level = NA_character_ ? character,
                                                       count_transformation = NA_character_ ? character){
                                   self$min_reads <- as.integer(min_reads)
                                   self$detection <- as.integer(detection)
                                   self$prevalence <- prevalence
                                   self$tax_glom_level <- tax_glom_level
                                   self$count_transformation <- count_transformation
                                 },
                                 execute = function(local_ps = ? phyloseq) {
                                   if(!is.na(self$detection)){
                                     taxas <- core_members(local_ps, detection = self$detection, prevalence = self$prevalence)
                                     local_ps <- prune_taxa(taxas, local_ps)
                                   }
                                   local_ps <- prune_samples(sample_sums(local_ps) > self$min_reads, local_ps)
                                   if(!is.na(self$tax_glom_level)){
                                     local_ps <- tax_glom(local_ps, taxrank=self$tax_glom_level)
                                   }
                                   if(!is.na(self$count_transformation)){
                                     local_ps <- microbiome::transform(local_ps, self$count_transformation)
                                   }
                                   return(local_ps)
                                 } ? phyloseq
                               )
)


