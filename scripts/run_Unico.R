#!/usr/bin/env Rscript
library(argparse)
library(dplyr)

parser <- ArgumentParser(add_help = F)
parser$add_argument("--refType",'-r',type="character")
parser$add_argument("--obj_dir",'-o',type="character",required = F,default = './',help = 'path to benchmarking objs directory')
parser$add_argument("--n_core", type = "integer", required = FALSE, default = 15)
parser$add_argument("--indep_ref_dir",type="character",required = F) 

args <- parser$parse_args()
refType = args$refType
obj_dir = args$obj_dir
n_core = args$n_core
indep_ref_dir = args$indep_ref_dir

contents = list.files(obj_dir)

if(!'init' %in% contents | !'reference' %in% contents){
  stop(paste('invalid obj_dir provided, unable to find `init` or `reference` folder from the obj_dir provided:',obj_dir))
}
  
export_dir = paste0(obj_dir,'/deconvRes/',paste0(refType,'_ref'))
if(!dir.exists(export_dir)){
  stop('deconvRes folder with InstaPrism inferred frac is required to run Unico')
}

bulk_expr = read.delim(paste0(obj_dir,'/init/bulk_expr.csv'),sep = ',')
bulk_expr = as.matrix(bulk_expr)

if(refType == 'self'){
  
  frac = read.delim(paste0(export_dir,'/InstaPrismFrac.csv'),sep = ',') %>% as.matrix()
  limma_top_genes = read.delim(paste0(obj_dir,'reference/self/limma_top_genes.csv'),sep = ',')
  
}else if(refType == 'indep'){

  frac = read.delim(paste0(export_dir,'/InstaPrismFrac.csv'),sep = ',',check.names = F) %>% as.matrix()
  limma_top_genes = read.delim(paste0(indep_ref_dir,'/limma_top_genes.csv'),sep = ',',check.names = F)
  
}else{
  stop('invalid refType argument provided')
}

dataset = basename(normalizePath(obj_dir))

# X must not include features with standard deviation less than 1e-04
gene_var = matrixStats::rowSds(bulk_expr)
keep_id = which(gene_var > 10^-4)
bulk_expr = bulk_expr[keep_id,]

library(Unico)

print('start running Unico ...')
start.time = Sys.time()

params.hat = Unico(bulk_expr,frac,C1 = NULL, C2 = NULL, num_cores = n_core)
Z.hat = Unico::tensor(bulk_expr,frac,C1=NULL,C2 =NULL,params.hat, num_cores = n_core)
Z.hat = aperm(Z.hat,c(2,3,1))

end.time = Sys.time()

saveRDS(Z.hat,file = paste0(export_dir,'/Unico_Z_inferred.RDS'))

# run time for Unico
execution_time = end.time - start.time
units(execution_time) = "mins"
execution_time = paste0(round(as.numeric(execution_time),2),'mins')

# write logs
current_log = data.frame(
  dataset = dataset,
  method = 'Unico',
  refType = refType,
  run_time = execution_time
)

if(!dir.exists(paste0(obj_dir,'/logs'))){
  dir.create(paste0(obj_dir,'/logs'))
}

log_file_path <- paste0(obj_dir,'/logs/runtime.txt')

if(file.exists(log_file_path)){
  existing_data <- read.delim(log_file_path,sep = '\t')
  new_data <- rbind(existing_data, current_log)
  write.table(new_data, log_file_path,sep = '\t',row.names = F,quote = F)
} else {
  write.table(current_log, log_file_path,sep = '\t',row.names = F,quote = F)
}

print('done')
