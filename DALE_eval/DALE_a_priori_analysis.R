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


for(benchmarking_obj in benchmarking_objs){
  
  print(paste('############### Start a priori analysis for',benchmarking_obj,'###############'))
  
  obj_dir = paste0('../Benchmarking_obj/',benchmarking_obj)
  export_dir = paste0('../Benchmarking_obj/',benchmarking_obj,'/init/')
  
  bulk_expr = read.delim(paste0('../Benchmarking_obj/',benchmarking_obj,'/init/bulk_expr.csv'),sep = ',')
  Z_truth = readRDS(paste0('../Benchmarking_obj/',benchmarking_obj,'/init/Z_truth.RDS'))

  Z_limma_statistics = get_Z_limma_statistics(Z_truth) %>% round(3)
  write.table(Z_limma_statistics,file = paste0(export_dir,'Z_limma_statistics.csv'),sep = ',')
}

print('done')




