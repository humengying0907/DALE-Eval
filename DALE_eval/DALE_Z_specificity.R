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
parser$add_argument("--add_self_ref", type = "integer", required = FALSE, default = 0,choices = c(0, 1),
                    help = "Use 1 to enable specificity evaluation on deconvRes from self ref or 0 to skip (FALSE). Default 0")


args <- parser$parse_args()
benchmarking_objs = args$benchmarking_objs
add_self_ref = as.logical(args$add_self_ref)

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

top_n = 100 # by default, will consider top 100 marker genes per cell type

Z_scales = list()
Z_scales[['linear']] = c('BayesPrism','cbsx','ENIGMAL2','ENIGMAtrace','ENIGMAL2(unnormalized)','ENIGMAtrace(unnormalized)',
                         'Unico','InstaPrism','InstaPrismUpdated')
Z_scales[['log2']] = c('bMIND','EPICunmix','TCA','bMIND(EPICunmix)','EPICunmix(EPICunmix)')
Z_scales_vec = unlist(Z_scales)
names(Z_scales_vec) <- Z_scales_vec
Z_scales_vec[] <- rep(names(Z_scales), lengths(Z_scales)) 


for(benchmarking_obj in benchmarking_objs){
  
  print(paste('############### Start specificity analysis for',benchmarking_obj,'###############'))
  
  obj_dir = paste0('../Benchmarking_obj/',benchmarking_obj)
  export_dir = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvSummary/')
  
  if(!dir.exists(export_dir)){
    dir.create(export_dir,recursive = T)
  }
  
  Z_truth_limma_statistics = read.delim(paste0('../Benchmarking_obj/',benchmarking_obj,'/init/Z_limma_statistics.csv'),sep = ',')
  indep_ref_mapping = read.delim(paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/indep_ref/indep_ref_mapping.csv'),sep = ',')
  Z_truth = readRDS(paste0('../Benchmarking_obj/',benchmarking_obj,'/init/Z_truth.RDS'))
  
  Z_truth_var_speicificity = get_variation_specificity(Z_truth,Z_truth_limma_statistics,NULL,top_n)
  Z_truth_expr_speicificity = get_expr_specificity(Z_truth,Z_truth_limma_statistics,NULL,top_n)
  
  var_specificity_self = list(Z_truth = Z_truth_var_speicificity)
  var_specificity_indep = list(Z_truth = Z_truth_var_speicificity)
  
  expr_specificity_self = list(Z_truth = Z_truth_expr_speicificity)
  expr_specificity_indep = list(Z_truth = Z_truth_expr_speicificity)
  
  for(method in methods){
    print(method)
    
    if(add_self_ref){
      Z_self_path = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/self_ref/',method,'_Z_inferred.RDS')
      if(file.exists(Z_self_path)){
        Z = readRDS(Z_self_path)
        var_specificity_self[[method]] = get_variation_specificity(Z,Z_truth_limma_statistics,NULL,top_n) 
        expr_specificity_self[[method]] = get_expr_specificity(Z,Z_truth_limma_statistics,NULL,top_n,Z_scales_vec[method])
      }
      saveRDS(var_specificity_self,file = paste0(export_dir,'/var_specificity_self_ref.RDS'))
      saveRDS(expr_specificity_self,file = paste0(export_dir,'/expr_specificity_self_ref.RDS'))
    }
    
    Z_indep_path = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/indep_ref/',method,'_Z_inferred.RDS')
    if(file.exists(Z_indep_path)){
      Z = readRDS(Z_indep_path)
      var_specificity_indep[[method]] = get_variation_specificity(Z,Z_truth_limma_statistics,indep_ref_mapping,top_n) 
      expr_specificity_indep[[method]] = get_expr_specificity(Z,Z_truth_limma_statistics,indep_ref_mapping,top_n,Z_scales_vec[method])
    }
  }  
  
  saveRDS(var_specificity_indep,file = paste0(export_dir,'/var_specificity_indep_ref.RDS'))
  saveRDS(expr_specificity_indep,file = paste0(export_dir,'/expr_specificity_indep_ref.RDS'))
}

print('done')




