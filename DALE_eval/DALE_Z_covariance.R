#!/usr/bin/env Rscript
library(argparse)
library(dplyr)

parser <- ArgumentParser(add_help = F)
parser$add_argument(
  "--benchmarking_objs", '-o',
  type = "character",
  nargs = '+', 
  required = FALSE,
  help = 'Benchmarking obj name (accepts multiple values)'
)

args <- parser$parse_args()
benchmarking_objs = args$benchmarking_objs

if(is.null(benchmarking_objs)){
  benchmarking_objs = list.files('../Benchmarking_obj/')
}

reference_mapping = read.delim('../Benchmarking_obj/reference_mapping.csv',sep = ',')

source('../scripts/evalu.R')
source('../scripts/helpers.R')

# by default, will only consider these methods
methods = c('bMIND',
            'cbsx',
            'ENIGMAtrace',
            'ENIGMAtrace(unnormalized)',
            'ENIGMAL2',
            'ENIGMAL2(unnormalized)',
            'EPICunmix',
            'InstaPrism',
            'TCA',
            'Unico') 


for(benchmarking_obj in benchmarking_objs){
  
  print(paste('############### Start covariance summary for',benchmarking_obj,'###############'))
  
  obj_dir = paste0('../Benchmarking_obj/',benchmarking_obj)
  export_dir = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvSummary/')
  
  indep_ref = reference_mapping$indep_ref[reference_mapping$datasets == benchmarking_obj]
  limma_top_genes = read.delim(paste0('../Indep_scReference/',indep_ref,'/limma_top_genes.csv'),sep = ',')
  indep_ref_mapping = read.delim(paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/indep_ref/indep_ref_mapping.csv'),sep = ',')
  limma_top_genes = indep_colnames_rename(limma_top_genes,indep_ref_mapping)

  Z_truth = readRDS(paste0('../Benchmarking_obj/',benchmarking_obj,'/init/Z_truth.RDS'))
  
  # will only considered covariance structure for matched cell types!
  offDiag_list = list()
  offDiag_list[['Z_truth']] = get_offDiagCor(Z_truth[,,colnames(limma_top_genes)])
  
  for(method in methods){
    Z_indep_path = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/indep_ref/',method,'_Z_inferred.RDS')
    
    if(file.exists(Z_indep_path)){
      Z_inferred = readRDS(Z_indep_path)
      Z_inferred = Z_inferred[,,indep_ref_mapping$maxCorName]
      dimnames(Z_inferred)[[3]] = indep_ref_mapping$cell_type
      stopifnot(setequal(colnames(limma_top_genes),dimnames(Z_inferred)[[3]]))
      
      offDiag_list[[method]] = get_offDiagCor(Z_inferred)
    }
  }
  
  all_genes <- unique(unlist(lapply(offDiag_list, names)))
  
  offDiag_cov_summary <- data.frame(gene = all_genes)
  for (i in seq_along(offDiag_list)) {
    offDiag_cov_summary[[names(offDiag_list)[i]]] <- offDiag_list[[i]][all_genes]
  }
  
  top_1000_marker_list = get_marker_list(limma_top_genes,select_method = 'top_n',n = 1000)
  
  offDiag_cov_summary$gene_group = ifelse(offDiag_cov_summary$gene %in% do.call(c,top_1000_marker_list),'marker','others')
  offDiag_cov_summary$dataset = benchmarking_obj
  saveRDS(offDiag_cov_summary,file = paste0(export_dir,'offDiag_cov_summary.RDS'))
}

print('done')



