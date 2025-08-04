#!/usr/bin/env Rscript
library(argparse)

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
  stop('deconvRes folder with InstaPrism inferred frac is required to run TCA')
}

bulk_expr = read.delim(paste0(obj_dir,'/init/bulk_expr.csv'),sep = ',')
bulk_expr = as.matrix(bulk_expr)

if(refType == 'self'){
  
  frac = read.delim(paste0(export_dir,'/InstaPrismFrac.csv'),sep = ',')
  limma_top_genes = read.delim(paste0(obj_dir,'reference/self/limma_top_genes.csv'),sep = ',')
  
}else if(refType == 'indep'){

  frac = read.delim(paste0(export_dir,'/InstaPrismFrac.csv'),sep = ',',check.names = F)
  limma_top_genes = read.delim(paste0(indep_ref_dir,'/limma_top_genes.csv'),sep = ',',check.names = F)
  
}else{
  stop('invalid refType argument provided')
}

dataset = basename(normalizePath(obj_dir))

# exclude genes with var less then 1e-08
X = log2(bulk_expr+1)
gene_var = matrixStats::rowVars(X)
keep_id = which(gene_var > 10^-8)
X = X[keep_id,]

library(TCA)

print('start running TCA ...')
start.time = Sys.time()
TCA_res <- tca(X = X, W = frac, parallel = T, num_cores = n_core)
Z_hat <- TCA::tensor(X = X, tca.mdl = TCA_res)
stopifnot(all.equal(colnames(TCA_res[["mus_hat"]]),colnames(frac)))
Z_array = array(NA,dim=c(nrow(Z_hat[[1]]),ncol(Z_hat[[1]]),length(Z_hat)))
dimnames(Z_array)[[1]] = rownames(Z_hat[[1]])
dimnames(Z_array)[[2]] = colnames(Z_hat[[1]])
dimnames(Z_array)[[3]] = colnames(frac)

for(i in 1:length(Z_hat)){
  Z_array[,,i] = Z_hat[[i]]
}

file.remove('TCA.log')
end.time = Sys.time()

saveRDS(Z_array,file = paste0(export_dir,'/TCA_Z_inferred.RDS'))

# run time for TCA
execution_time = end.time - start.time
units(execution_time) = "mins"
execution_time = paste0(round(as.numeric(execution_time),2),'mins')

# write logs
current_log = data.frame(
  dataset = dataset,
  method = 'TCA',
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
