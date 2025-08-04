#!/usr/bin/env Rscript
library(argparse)
library(dplyr)
library(ENIGMA)

parser <- ArgumentParser(add_help = F)
parser$add_argument("--refType",'-r',type="character")
parser$add_argument("--obj_dir",'-o',type="character",required = F,default = './',help = 'path to benchmarking objs directory')
parser$add_argument("--n_core", type = "integer", required = FALSE, default = 15)
parser$add_argument("--indep_ref_dir",type="character",required = F) 
parser$add_argument("--ENIGMAmode",'-m',type="character",required = T,default = 'L2')

args <- parser$parse_args()
refType = args$refType
obj_dir = args$obj_dir
n_core = args$n_core
indep_ref_dir = args$indep_ref_dir
ENIGMAmode = args$ENIGMAmode

contents = list.files(obj_dir)

if(!'init' %in% contents | !'reference' %in% contents){
  stop(paste('invalid obj_dir provided, unable to find `init` or `reference` folder from the obj_dir provided:',obj_dir))
}
  
stopifnot(ENIGMAmode %in% c('L2','trace'))

export_dir = paste0(obj_dir,'/deconvRes/',paste0(refType,'_ref'))
if(!dir.exists(export_dir)){
  stop('deconvRes folder with InstaPrism inferred frac is required to run ENIGMA')
}

bulk_expr = read.delim(paste0(obj_dir,'/init/bulk_expr.csv'),sep = ',')
bulk_expr = as.matrix(bulk_expr)

if(refType == 'self'){
  
  frac = read.delim(paste0(export_dir,'/InstaPrismFrac.csv'),sep = ',') %>% as.matrix()
  ENIGMA_ref = read.delim(paste0(obj_dir,'reference/self/rowMeans_sig.csv'),sep = ',')  %>% as.matrix()
  limma_top_genes = read.delim(paste0(obj_dir,'reference/self/limma_top_genes.csv'),sep = ',')
  
}else if(refType == 'indep'){

  frac = read.delim(paste0(export_dir,'/InstaPrismFrac.csv'),sep = ',',check.names = F) %>% as.matrix()
  ENIGMA_ref = read.delim(paste0(indep_ref_dir,'/rowMeans_sig.csv'),sep = ',')  %>% as.matrix()
  limma_top_genes = read.delim(paste0(indep_ref_dir,'/limma_top_genes.csv'),sep = ',',check.names = F)
  
}else{
  stop('invalid refType argument provided')
}

dataset = basename(normalizePath(obj_dir))

cms = intersect(rownames(bulk_expr),rownames(ENIGMA_ref))

bulk_expr = bulk_expr[cms,]
ENIGMA_ref = ENIGMA_ref[cms,]

# build ENIGMA object
egm = create_ENIGMA(bulk = bulk_expr, ref = ENIGMA_ref,ref_type = "aggre")
egm@result_cell_proportion = frac

print('start running ENIGMA ...')
start.time = Sys.time()

if(ENIGMAmode == 'L2'){
  method = 'ENIGMAL2'
  egm = ENIGMA_L2_max_norm(egm, alpha = 0.1,model_tracker = TRUE,model_name = "log", preprocess = "log")
  cse_normalized = sce2array(egm,norm_output = T, model_name = "log") 
  cse_unnormalized = sce2array(egm,norm_output = F, model_name = "log") 
  
}else if(ENIGMAmode == 'trace'){
  method = 'ENIGMAtrace'
  egm = ENIGMA_trace_norm(egm, alpha = 0.1,model_tracker = TRUE,model_name = "log", preprocess = "log")
  cse_normalized = sce2array(egm,norm_output = T, model_name = "log") 
  cse_unnormalized = sce2array(egm,norm_output = F, model_name = "log") 
  
}

end.time = Sys.time()

saveRDS(cse_normalized,file = paste0(export_dir,'/',method,'_Z_inferred.RDS'))
saveRDS(cse_unnormalized,file = paste0(export_dir,'/',method,'(unnormalized)_Z_inferred.RDS'))


# run time for ENIGMA
execution_time = end.time - start.time
units(execution_time) = "mins"
execution_time = paste0(round(as.numeric(execution_time),2),'mins')

# write logs
current_log = data.frame(
  dataset = dataset,
  method = method,
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
