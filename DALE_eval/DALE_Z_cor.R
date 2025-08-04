#!/usr/bin/env Rscript
library(argparse)
library(dplyr)

# important note:
# must have indep_ref_mapping.csv available in deconvRes/indep_ref first before running Z-cor module!

parser <- ArgumentParser(add_help = F)
parser$add_argument(
  "--benchmarking_objs", '-o',
  type = "character",
  nargs = '+', 
  required = FALSE,
  help = 'Benchmarking obj name (accepts multiple values)'
)
parser$add_argument("--cor_method",'-r',type="character",required = F, default = 'spearman',choices = c("spearman", "pearson"))

args <- parser$parse_args()
benchmarking_objs = args$benchmarking_objs
cor_method = args$cor_method

if(cor_method == 'spearman'){
  prefix = 'cor_res_'
}else if(cor_method == 'pearson'){
  prefix = 'pearson_cor_res_'
}

if(is.null(benchmarking_objs)){
  benchmarking_objs = list.files('../Benchmarking_obj/')
}

source('../scripts/evalu.R')
source('../scripts/helpers.R')

get_Z_cor = function(obj_dir = './',
                     refType = 'self',
                     cor_method = 'spearman'){
  
  # please note that Pearson correlation is depreciated for Z_inferred in log scale
  # since transforming either Z_truth to log scale or Z_inferred by exponentiation will still result in incomparable values
  # in such case, we will directly compare Z_inferred with Z_truth without further transformation
  
  require(dplyr)
  Z_truth = readRDS(paste0(obj_dir,'/init/Z_truth.RDS'))
  bulk_expr = read.delim(paste0(obj_dir,'/init/bulk_expr.csv'),sep = ',') %>% as.matrix()
  truthFrac = read.delim(paste0(obj_dir,'/init/truthFrac.csv'),sep = ',') %>% as.matrix()

  
  res_dir = paste0(obj_dir,'/deconvRes/',paste0(refType,'_ref'))
  InstaPrismFrac = read.delim(paste0(res_dir,'/InstaPrismFrac.csv'),sep = ',') %>% as.matrix()
  if(refType == 'indep'){
    indep_mapping = read.delim(paste0(res_dir,'/indep_ref_mapping.csv'),sep = ',')
  }
  
  bulk_regressed = linear_regress_out(bulk_expr, model.matrix(~0+InstaPrismFrac))
  
  res_files = list.files(res_dir)
  res_files <- res_files[grepl('Z_inferred.RDS', res_files)]
  
  bulk_Z = get_3d_bulk(bulk_expr,dimnames(Z_truth)[[3]])
  bulk_regressed_Z = get_3d_bulk(bulk_regressed,dimnames(Z_truth)[[3]])
  
  print('start baseline correlation calculation ...')
  cor_res_list = list()
  cor_res_list[['bulk']] = get_per_gene_cor(bulk_Z,Z_truth,truthFrac,margin = 'row',cor_method = cor_method)
  cor_res_list[['bulk_regressed']] = get_per_gene_cor(bulk_regressed_Z,Z_truth,truthFrac,margin = 'row',cor_method = cor_method)
  
  for(res_file in res_files){
    Z_inferred = readRDS(paste0(res_dir,'/',res_file))
    method = sub('_Z_inferred.RDS','',res_file)
    print(paste('Start correlation calculation for',method))
    
    # map cell_type names in Z_inferred
    if(refType == 'indep'){
      Z_inferred = Z_inferred[,,indep_mapping$maxCorName]
      dimnames(Z_inferred)[[3]] = indep_mapping$cell_type
    }
    
    cor_res_list[[method]] = get_per_gene_cor(Z_inferred,Z_truth,truthFrac,margin = 'row',cor_method = cor_method)
  }
  
  return(cor_res_list)
}

for(obj in benchmarking_objs){
  
  print(paste('################### Start Z-cor calculation for',obj,'###################'))
  
  export_dir = paste0('../Benchmarking_obj/',obj,'/deconvSummary/')
  
  if(!dir.exists(export_dir)){
    dir.create(export_dir)
  }
  
  print('*** self ref ***')
  cor_res_list_self = get_Z_cor(obj_dir = paste0('../Benchmarking_obj/',obj),refType = 'self',cor_method = cor_method)
  saveRDS(cor_res_list_self,file = paste0(export_dir,prefix,'self_ref.RDS'))
  
  print('*** indep ref ***')
  cor_res_list_indep = get_Z_cor(obj_dir = paste0('../Benchmarking_obj/',obj),refType = 'indep',cor_method = cor_method)
  saveRDS(cor_res_list_indep,file = paste0(export_dir,prefix,'indep_ref.RDS'))
  
}

print('done')

