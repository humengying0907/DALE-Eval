#!/usr/bin/env Rscript
library(argparse)
library(dplyr)
library(tibble)

parser <- ArgumentParser(add_help = F)
parser$add_argument("--refType",'-r',type="character")
parser$add_argument("--obj_dir",'-o',type="character",required = F,default = './',help = 'path to benchmarking objs directory')
parser$add_argument("--indep_ref_dir",type="character",required = F)
parser$add_argument("--n_core", type = "integer", required = FALSE, default = 1) # set this to 1 to avoid Segmentation fault (core dumped)
parser$add_argument("--username",type="character",required = T)
parser$add_argument("--token",type="character",required = T) 
parser$add_argument("--singularity_container_path",type="character",required = T)

args <- parser$parse_args()

refType = args$refType
obj_dir = args$obj_dir
indep_ref_dir = args$indep_ref_dir
n_core = args$n_core
username = args$username
token = args$token
singularity_container_path = args$singularity_container_path
  
  
export_dir = paste0(obj_dir,'/deconvRes/',paste0(refType,'_ref'))
if(!dir.exists(export_dir)){
  stop('deconvRes folder with InstaPrism inferred frac is required to run bMIND')
}


bulk = read.delim(paste0(obj_dir,'/init/bulk_expr.csv'),sep = ',')

if(refType == 'self'){
  frac = read.delim(paste0(export_dir,'/InstaPrismFrac.csv'),sep = ',') 
  cbsx_sig= read.delim(paste0(obj_dir,'/reference/self/cbsx_sig.txt'),sep = '\t',check.names = F)
  
}else if(refType == 'indep'){
  frac = read.delim(paste0(export_dir,'/InstaPrismFrac.csv'),sep = ',',check.names = F) 
  cbsx_sig = read.delim( paste0(indep_ref_dir,'/cbsx_sig.txt'),sep = '\t',check.names = F)
  
}else{
  stop('invalid refType argument provided')
}


# transform to cbsx required format and store in tempdir
bulk = cbind(data.frame(GeneSymbol = rownames(bulk)),bulk)
frac = frac %>% rownames_to_column('Mixture')
cbsx_sig = cbsx_sig[cbsx_sig$GeneSymbol %in% rownames(bulk),]

# export bulk expr to a temp cbsx_path
temp_dir <- tempdir()
write.table(bulk,file = paste0(temp_dir,'/bulk.txt'),sep = '\t',row.names = F,quote = F)
write.table(frac,file = paste0(temp_dir,'/frac.txt'),sep = '\t',row.names = F,quote = F)
write.table(cbsx_sig,file = paste0(temp_dir,'/cbsx_sig.txt'),sep = '\t',row.names = F,quote = F)

print(list.files(temp_dir))

library(omnideconv)

# run cibersortx_hires
cbsx_hires_command = function(cbsx_path,username,token,singularity_container_path){
  command = paste0('singularity exec -c -B ',cbsx_path, 
                   '/:/src/data -B ',cbsx_path,
                   '/:/src/outdir ',
                   singularity_container_path, 
                   ' /src/CIBERSORTxHiRes --mixture bulk.txt --sigmatrix cbsx_sig.txt --cibresults frac.txt --label HiRes --threads ',n_core,
                   ' --heatmap FALSE --username ',username,' --token ',token)
  return(command)
}

command_to_run = cbsx_hires_command(temp_dir,username,token,singularity_container_path)
verbose = T
start.time = Sys.time()
code <- system(command_to_run, ignore.stdout = !verbose, ignore.stderr = !verbose) 


hires_res_files = list.files(temp_dir)
cell_types = colnames(cbsx_sig)[-1]

if(length(hires_res_files) < 5){
  stop('cbsx singularity job failed')
}

g = nrow(bulk)
n = ncol(bulk)-1

cbsx_inferred = array(NA, dim = c(g,n,length(cell_types)))
dimnames(cbsx_inferred)[[3]] = cell_types

for(i in 1:length(cell_types)){
  ct = cell_types[i]
  print(ct)
  cts_out = read.delim(paste0(temp_dir,'/',hires_res_files[grepl(paste0(ct,'_Window'),hires_res_files)])) %>% 
    column_to_rownames('GeneSymbol') %>% as.matrix()
  
  print(dim(cts_out)) # even though some genes are empty, cbsx still export lines for all genes
  stopifnot(all.equal(dim(cts_out),dim(bulk[,-1])))
  cbsx_inferred[,,i] = cts_out
}

dimnames(cbsx_inferred)[[1]] = rownames(bulk)
dimnames(cbsx_inferred)[[2]] = colnames(bulk)[-1]

end.time = Sys.time()

saveRDS(cbsx_inferred,file = paste0(export_dir,'/cbsx_Z_inferred.RDS'))



# run time for cibersortx
execution_time = end.time - start.time
units(execution_time) = "mins"
execution_time = paste0(round(as.numeric(execution_time),2),'mins')

# write logs
current_log = data.frame(
  dataset = basename(normalizePath(obj_dir)),
  method = 'CIBERSORTx',
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
