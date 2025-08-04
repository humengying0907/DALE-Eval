#!/usr/bin/env Rscript
library(argparse)
library(InstaPrism)

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

for(benchmarking_obj in benchmarking_objs){
  
  print(paste('################### start processing',benchmarking_obj,'###################'))
  
  truthFrac = read.delim(paste0('../Benchmarking_obj/',benchmarking_obj,'/init/truthFrac.csv'),sep = ',')
  indep_frac = read.delim(paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/indep_ref/InstaPrismFrac.csv'),sep = ',')
  
  l = InstaPrism::frac_evalu(truthFrac,indep_frac)
  cor_summ = l$summ %>% rownames_to_column('cell_type')
  
  poor_matching = cor_summ$cell_type[cor_summ$cor < 0.6]
  if(length(poor_matching)>0){
    print(paste(
      'cell types with poor frac correlations:',
      paste(poor_matching, collapse=", ")
    ))
  }
  
  unused_indep_cell_types = colnames(indep_frac)[!colnames(indep_frac) %in% cor_summ$maxCorName]
  if(length(unused_indep_cell_types)>0){
    print(paste(
      'cell types in indep reference without a matched cell type:',
      paste(unused_indep_cell_types, collapse=", ")
    ))
  }
  
  used_indep_cell_types = gsub("\\.[0-9]+", "", cor_summ$maxCorName) # reverse the make.names() effect
  if(any(duplicated(used_indep_cell_types))){
    print(paste(
      'cell types in indep reference matched to multiple ground truth cell types:',
      paste(used_indep_cell_types[duplicated(used_indep_cell_types)], collapse=", ")
    ))
  }
  
  if ((length(poor_matching) > 0) | (length(unused_indep_cell_types) > 0) | any(duplicated(used_indep_cell_types))) {
    print('a raw matching table is exported in deconvRes/indep_ref dir, please manually adjust it to include matched cell-types only, and save it as indep_ref_mapping.csv')
    write.table(cor_summ,file = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/indep_ref/indep_ref_mapping_raw.csv'),sep = ',')
  }else{
    write.table(cor_summ,file = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/indep_ref/indep_ref_mapping.csv'),sep = ',')
  }
}

