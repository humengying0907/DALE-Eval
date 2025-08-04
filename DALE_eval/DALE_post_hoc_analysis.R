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

# prepare a vector indicating the scales of Z_inferred
Z_scales = list()
Z_scales[['linear']] = c('BayesPrism','cbsx','ENIGMAL2','ENIGMAtrace',
                         'ENIGMAL2(unnormalized)','ENIGMAtrace(unnormalized)',
                         'Unico','InstaPrism','InstaPrismUpdated')
Z_scales[['log2']] = c('bMIND','EPICunmix','TCA')
Z_scales_vec = unlist(Z_scales)
names(Z_scales_vec) <- Z_scales_vec
Z_scales_vec[] <- rep(names(Z_scales), lengths(Z_scales)) 


for(benchmarking_obj in benchmarking_objs){
  
  print(paste('############### Start post hoc analysis for',benchmarking_obj,'###############'))
  
  obj_dir = paste0('../Benchmarking_obj/',benchmarking_obj)
  export_dir = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvSummary/post_hoc_statisitics/')
  
  bulk_expr = read.delim(paste0('../Benchmarking_obj/',benchmarking_obj,'/init/bulk_expr.csv'),sep = ',')
  indep_mapping = read.delim(paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/indep_ref/indep_ref_mapping.csv'),sep = ',')
  
  if(!dir.exists(export_dir)){
    dir.create(export_dir,recursive = T)
  }
  
  subdirs <- c("Z_limma_statistics", "Z_mean_fc", "Z_partial_r2")
  for (ref in  c("self_ref", "indep_ref")) {
    for (sub in subdirs) {
      dir.create(file.path(export_dir, ref, sub), recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  for(method in methods){
    print(paste('**********',method,'**********'))
    
    for(refType in c('self','indep')){
      Z_path = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/',refType,'_ref/',method,'_Z_inferred.RDS')
      
      if(file.exists(Z_path)){
        
        Z_inferred = readRDS(Z_path)
        
        Z_limma_statistics = get_Z_limma_statistics(Z_inferred,Z_scales_vec[method]) %>% round(3)
        Z_means = compute_Z_mean(Z_inferred)
        Z_means_fc = compute_Z_mean_fc(Z_means,Z_scales_vec[method]) %>% round(3)
        Z_partial_r2 = compute_partial_r2(Z_inferred,bulk_expr)
        
        if(refType == 'indep'){
          
          Z_limma_statistics = indep_colnames_rename(Z_limma_statistics,indep_mapping)
          Z_means_fc = indep_colnames_rename(Z_means_fc,indep_mapping)
          Z_partial_r2 = indep_colnames_rename(Z_partial_r2,indep_mapping)

        }
        
        write.table(Z_limma_statistics,file = paste0(export_dir,refType,'_ref/Z_limma_statistics/',method,'_Z_limma_statistics.csv'),sep = ',')
        write.table(Z_means_fc,file = paste0(export_dir,refType,'_ref/Z_mean_fc/',method,'_Z_mean_fc.csv'),sep = ',')
        write.table(Z_partial_r2,file = paste0(export_dir,refType,'_ref/Z_partial_r2/',method,'_Z_partial_r2.csv'),sep = ',')

      }
    }
  }  
}

print('done')




