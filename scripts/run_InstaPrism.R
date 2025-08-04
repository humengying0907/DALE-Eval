#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(add_help = F)
parser$add_argument("--refType",'-r',type="character")
parser$add_argument("--obj_dir",'-o',type="character",required = F,default = './',help = 'path to benchmarking objs directory')
parser$add_argument("--indep_ref_dir",type="character",required = F)
parser$add_argument("--n_core", type = "integer", required = FALSE, default = 15)
parser$add_argument("--return_raw", type = "integer", required = FALSE, default = 0,choices = c(0, 1),
                    help = "Use 1 to export BayesPrism_Z_inferred.RDS (TRUE) or 0 to skip (FALSE). Default 0")
parser$add_argument("--InstaPrism_update", type = "integer", required = FALSE, default = 0,choices = c(0, 1),
                    help = "Use 1 to enable InstaPrism update and export InstaPrismUpdated_Z_inferred.RDS (TRUE) or 0 to skip (FALSE). Default 0")

args <- parser$parse_args()
refType = args$refType
obj_dir = args$obj_dir
indep_ref_dir = args$indep_ref_dir
n_core = args$n_core
return_raw = as.logical(args$return_raw)
InstaPrism_update = as.logical(args$InstaPrism_update)

################### functions for InstaPrism normalization and frac matching #################
normalize_lapCorrected = function(slice,libSize,Nconstant,Dconstant){
  division_value = (colSums(slice) + Nconstant)/(libSize  + Dconstant)
  scaled_slice <- sweep(slice, 2, division_value, '/')
  return(scaled_slice)
}

InstaPrism_normalize = function(data, normalize_method = 'lapCorrected'){
  stopifnot(normalize_method %in% c('lapCorrected','fracDiv'))

  scaled_array_3d <- array(dim = dim(data))
  dimnames(scaled_array_3d) <- dimnames(data)
  libSize = apply(data, MARGIN = 2, sum) # bulk libsize
  
  if((range(libSize)[2] - range(libSize)[1])> 0.01 * mean(libSize)){
    warning('inconsistent libsize between reconstructed bulk profiles')
  }
  
  k = dim(data)[[3]]
  Nconstant = 0.1*libSize/k 
  Dconstant = 0.1*libSize/k
  
  for(cell_type in dimnames(data)[[3]]) {
    slice <- data[,,cell_type]
    
    if(normalize_method == 'lapCorrected'){
      scaled_slice = normalize_lapCorrected(slice,libSize,Nconstant,Dconstant)
    }else if(normalize_method == 'fracDiv'){
      scaled_slice = normalize_lapCorrected(slice,libSize,Nconstant = 0,Dconstant = 0)
    }
    scaled_array_3d[,,cell_type] <- scaled_slice
  }
  return(scaled_array_3d)
}


##################
library(InstaPrism)
key_table = data.frame(dataset = c('BRCA_Bassez2021','BRCA_Wu2021',
                                   'CRC_Lee2020','CRC_Pelka2021',
                                   'LUAD_Kim2020','NSCLC_Wu2021','LUAD_Laughney2020'),
                       key = c('Cancer_cell','Cancer_Epithelial',
                               'Malignant','Malignant',
                               'Malignant','Malignant','Malignant'))

stopifnot(refType %in% c('self','indep'))

dataset = basename(normalizePath(obj_dir))

if(refType == 'self'){
  contents = list.files(obj_dir)
  if(!'init' %in% contents | !'reference' %in% contents){
    stop(paste('invalid obj_dir provided, unable to find `init` or `reference` folder from the obj_dir provided:',obj_dir))
  }
  refPhi = readRDS(paste0(obj_dir,'/reference/self/refPhi.RDS'))
  
  if(grepl('PBMC',dataset) | grepl('ROSMAP',dataset)){
    key = NA
  }else{
    key = key_table$key[key_table$dataset == dataset]
  }
  
}else if(refType == 'indep'){
  if(grepl("/$", indep_ref_dir)){
    indep_ref_dir <- sub("/$", "", indep_ref_dir)
  }
  indep_refName = basename(indep_ref_dir)
  
  if(grepl('PBMC',dataset) | grepl('ROSMAP',dataset)){
    key = NA
  }else{
    key = key_table$key[key_table$dataset == indep_refName]
  }
  
  refPhi = readRDS(paste0(indep_ref_dir,'/refPhi.RDS'))
  
}else{
  stop('invalid refType argument provided')
}

export_dir = paste0(obj_dir,'/deconvRes/',paste0(refType,'_ref'))
if(!dir.exists(export_dir)){
  dir.create(export_dir,recursive = TRUE)
}

bulk_expr = read.delim(paste0(obj_dir,'/init/bulk_expr.csv'),sep = ',')
bulk_expr = as.matrix(bulk_expr)

start.time = Sys.time()
InstaPrism_res = InstaPrism(bulk_Expr = bulk_expr,refPhi_cs = refPhi, n.iter = 500, filter = F,n.core = n_core) 
inferred_frac= t(InstaPrism_res@Post.ini.ct@theta)
Z = get_Z_array(InstaPrism_res,resolution = 'ct',n.core = n_core)
Z = aperm(Z,c(2,1,3))
# apply normalization
Z_normalized = InstaPrism_normalize(Z,normalize_method = 'lapCorrected')
end.time = Sys.time()

write.table(inferred_frac,file = paste0(export_dir,'/InstaPrismFrac.csv'),sep = ',')
saveRDS(Z_normalized,file = paste0(export_dir,'/InstaPrism_Z_inferred.RDS'))

if(return_raw){
  saveRDS(Z,file = paste0(export_dir,'/BayesPrism_Z_inferred.RDS'))
}

# document run time
execution_time1 = end.time - start.time
units(execution_time1) = "mins"
execution_time1 = paste0(round(as.numeric(execution_time1),2),'mins')

# get InstaPrism updated results
if(InstaPrism_update){
  start.time = Sys.time()
  InstaPrism_res_updated = InstaPrism_update(InstaPrism_res,bulk_expr,key = key, n.iter = 500, n.core = n_core)
  Z_updated = get_Z_array(InstaPrism_res_updated,resolution = 'ct',n.core = n_core)
  Z_updated = aperm(Z_updated,c(2,1,3))
  updated_frac = t(InstaPrism_res_updated@theta)
  Z_updated_normalized = InstaPrism_normalize(Z_updated,normalize_method = 'lapCorrected')
  end.time = Sys.time()
  
  saveRDS(Z_updated_normalized,file = paste0(export_dir,'/InstaPrismUpdated_Z_inferred.RDS'))
  
  # document run time for InstaPrismUpdate
  execution_time2 = end.time - start.time
  units(execution_time2) = "mins"
  execution_time2 = as.numeric(execution_time2) + as.numeric(gsub('mins','',execution_time1))
  execution_time2 = paste0(round(as.numeric(execution_time2),2),'mins')
  # write logs
  current_log = data.frame(
    dataset = rep(dataset,2),
    method = c('InstaPrism','InstaPrismUpdated'),
    refType = rep(refType,2),
    run_time = c(execution_time1,execution_time2)
  )
}else{
  # write logs
  current_log = data.frame(
    dataset = dataset,
    method = 'InstaPrism',
    refType = refType,
    run_time = execution_time1
  )
}

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
