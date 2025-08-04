#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(add_help = F)
parser$add_argument("--refType",'-r',type="character")
parser$add_argument("--obj_dir",'-o',type="character",required = F,default = './',help = 'path to benchmarking objs directory')
parser$add_argument("--indep_ref_dir",type="character",required = F)
parser$add_argument("--n_core", type = "integer", required = FALSE, default = 15)

args <- parser$parse_args()
refType = args$refType
obj_dir = args$obj_dir
indep_ref_dir = args$indep_ref_dir
n_core = args$n_core

contents = list.files(obj_dir)

if(!'init' %in% contents | !'reference' %in% contents){
  stop(paste('invalid obj_dir provided, unable to find `init` or `reference` folder from the obj_dir provided:',obj_dir))
}

dataset = basename(normalizePath(obj_dir))

export_dir = paste0(obj_dir,'/deconvRes/',paste0(refType,'_ref'))
if(!dir.exists(export_dir)){
  stop('deconvRes folder with InstaPrism inferred frac is required to run bMIND')
}

bulk_expr = read.delim(paste0(obj_dir,'/init/bulk_expr.csv'),sep = ',')
bulk_expr = as.matrix(bulk_expr)

if(refType == 'self'){
  profile = read.delim(paste0(obj_dir,'/reference/self/bMIND_profile.csv'),sep = ',')
  frac = read.delim(paste0(export_dir,'/InstaPrismFrac.csv'),sep = ',')
  stopifnot(setequal(colnames(profile),colnames(frac)))
  frac = frac[,colnames(profile)]
  
  if(file.exists(paste0(obj_dir,'/reference/self/bMIND_covariance.RDS'))){
    bMIND_covariance = readRDS(paste0(obj_dir,'/reference/self/bMIND_covariance.RDS'))
    stopifnot(all.equal(colnames(profile),dimnames(bMIND_covariance)[[1]]))
  }else{
    bMIND_covariance = NULL
  }
  
  limma_top_genes = read.delim(paste0(obj_dir,'reference/self/limma_top_genes.csv'),sep = ',')
  
}else if(refType == 'indep'){
  if(is.null(indep_ref_dir)){
    stop('please provide indep_ref_dir')
  }
  
  profile = read.delim(paste0(indep_ref_dir,'/bMIND_profile.csv'),sep = ',',check.names = F)
  frac = read.delim(paste0(export_dir,'/InstaPrismFrac.csv'),sep = ',',check.names = F)
  stopifnot(setequal(colnames(profile),colnames(frac)))
  frac = frac[,colnames(profile)]
  
  if(file.exists(paste0(indep_ref_dir,'/bMIND_covariance.RDS'))){
    bMIND_covariance = readRDS(paste0(indep_ref_dir,'/bMIND_covariance.RDS'))
    stopifnot(all.equal(rownames(profile),dimnames(bMIND_covariance)[[1]]))
    stopifnot(all.equal(colnames(profile),dimnames(bMIND_covariance)[[2]]))
    
  }else{
    bMIND_covariance = NULL
  }
  
  limma_top_genes = read.delim(paste0(indep_ref_dir,'/limma_top_genes.csv'),sep = ',',check.names = F)
  
}else{
  stop('invalid refType argument provided')
}

cms = intersect(rownames(bulk_expr),rownames(profile))
print(paste('total genes for deconvolution:',length(cms)))

bulk_expr = bulk_expr[cms,]
profile = profile[cms,]

if(!is.null(bMIND_covariance)){
  bMIND_covariance = bMIND_covariance[cms,,]
}

library(EPICunmix)
stopifnot(all.equal(rownames(bulk_expr),rownames(profile)))

print('start running bMIND ...')
start.time = Sys.time()
posterior = MIND::bMIND(log2(bulk_expr + 1), frac, profile = profile, covariance = bMIND_covariance, ncore = n_core)
Z = aperm(posterior$A, c(1, 3, 2))
dimnames(Z)[[2]] = colnames(bulk_expr)
dimnames(Z)[[3]] = colnames(frac)
end.time = Sys.time()

saveRDS(Z,file = paste0(export_dir,'/bMIND_Z_inferred.RDS'))
# saveRDS(posterior,file = paste0(export_dir,'/bMIND_posterior.RDS')) # export posterior as well in case EPICunmix failed

# run time for bMIND
execution_time1 = end.time - start.time
units(execution_time1) = "mins"
execution_time1 = paste0(round(as.numeric(execution_time1),2),'mins')

# EPICunmix
print('start running EPICunmix using posterior from bMIND ...')
start.time = Sys.time()
epic_unmix = EPICunmix::run_epic_unmix(log2(bulk_expr + 1), frac, posterior, outf = F, ncore = n_core)
Z2 = aperm(epic_unmix$A, c(1, 3, 2))
dimnames(Z2)[[3]] = dimnames(Z)[[3]] # reverse the effect of name changes
end.time = Sys.time()
saveRDS(Z2,file = paste0(export_dir,'/EPICunmix_Z_inferred.RDS'))

# document run time for EPICunmix
execution_time2 = end.time - start.time
units(execution_time2) = "mins"
execution_time2 = as.numeric(execution_time2) + as.numeric(gsub('mins','',execution_time1))
execution_time2 = paste0(round(as.numeric(execution_time2),2),'mins')

# write logs
current_log = data.frame(
  dataset = rep(dataset,2),
  method = c('bMIND','EPICunmix'),
  refType = rep(refType,2),
  run_time = c(execution_time1,execution_time2)
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
