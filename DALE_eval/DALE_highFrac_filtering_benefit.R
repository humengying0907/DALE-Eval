#!/usr/bin/env Rscript
library(argparse)
library(dplyr)

# must run DALE_Z_cor.R and DALE_Z_cor_highFrac_samples.R first before running this script!
# note: by default will only consider n_vec = c(10,30,100,300,1000,3000,10000) and spearman correlation!

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

source('../scripts/evalu.R')
source('../scripts/helpers.R')

if(is.null(benchmarking_objs)){
  benchmarking_objs = list.files('../Benchmarking_obj/')
}

reference_mapping = read.delim('../Benchmarking_obj/reference_mapping.csv',sep = ',')

get_benefit_summary = function(dataset,
                               refType = 'indep', 
                               n_vec = c(10,30,100,300,1000,3000,10000)){

  obj_dir = paste0('../Benchmarking_obj/',dataset)
  cor_res_list = readRDS(paste0(obj_dir,'/deconvSummary/cor_res_',refType,'_ref.RDS'))
  highFrac_samples_cor_res_list = readRDS(paste0(obj_dir,'/deconvSummary/cor_res_highFrac_samples_',refType,'_ref.RDS'))
  
  if(refType == 'self'){
    Z_limma_statistics =  read.delim(paste0(obj_dir,'/reference/self/limma_top_genes.csv'),sep = ',')
    frac_inferred = read.delim(paste0(obj_dir,'/deconvRes/self_ref/InstaPrismFrac.csv'),sep = ',')
    
  }else if(refType == 'indep'){
    indep_d = reference_mapping$indep_ref[reference_mapping$datasets == basename(normalizePath(obj_dir))]
    Z_limma_statistics = read.delim(paste0('../Indep_scReference/',indep_d,'/limma_top_genes.csv'),sep = ',')
    indep_mapping = read.delim(paste0(obj_dir,'/deconvRes/indep_ref/indep_ref_mapping.csv'),sep = ',')
    Z_limma_statistics =  indep_colnames_rename(Z_limma_statistics,indep_mapping)
    Z_limma_statistics = Z_limma_statistics[rownames(Z_limma_statistics) %in% rownames(cor_res_list$bulk),]
    
    frac_inferred = read.delim(paste0(obj_dir,'/deconvRes/indep_ref/InstaPrismFrac.csv'),sep = ',')
    frac_inferred = indep_colnames_rename(frac_inferred,indep_mapping)
  }
  
  highFrac_nSamples = colSums(frac_inferred>0.1) # please make sure this is consistent with min_frac used in DALE_highFrac_filtering_benefit.R
  inferred_ct_abundance = colMeans(frac_inferred)
  
  a = get_by_specificity_summary(cor_res_list,Z_limma_statistics,'top_n', n_vec, NULL, T)
  b = get_by_specificity_summary(highFrac_samples_cor_res_list,Z_limma_statistics,'top_n', n_vec, NULL, T)
  
  summary_df = a[,c('method','cell_type','group','n_genes','avg_cor')]
  colnames(summary_df)[5] = 'avg_cor_all_samples'
  
  summary_df$avg_cor_highFrac_samples = NA
  
  key1 = paste0(a$method,a$cell_type,a$group)
  key2 = paste0(b$method,b$cell_type,b$group)
  
  summary_df$avg_cor_highFrac_samples = b$avg_cor[match(key1,key2)]
  
  summary_df$sample_size = nrow(frac_inferred)
  summary_df$highFrac_sample_size = highFrac_nSamples[match(summary_df$cell_type,names(highFrac_nSamples))]
  summary_df$inferred_ct_abundance = inferred_ct_abundance[match(summary_df$cell_type,names(inferred_ct_abundance))]
  
  summary_df$refType = refType
  summary_df$dataset = dataset
  
  return(summary_df)
}

for(benchmarking_obj in benchmarking_objs){
  
  print(benchmarking_obj)
  export_dir = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvSummary/')
  
  self_summary = get_benefit_summary(benchmarking_obj,'self')
  indep_summary = get_benefit_summary(benchmarking_obj,'indep')
  
  out_df = rbind(self_summary,indep_summary)
  saveRDS(out_df,file = paste0(export_dir,'highFrac_filtering_benefit.RDS'))
}
print('done')
