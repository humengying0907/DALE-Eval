#!/usr/bin/env Rscript
library(argparse)
library(dplyr)

# only consider hallmark genes for sample correlation
# will exclude samples with frac < 0.001 for sample correlation
read.gmt<-function(gmt.file){
  g<-GSA::GSA.read.gmt(gmt.file)
  genelist<-g$genesets
  names(genelist)<-g$geneset.names
  rm(g)
  return(genelist)
}

hallmark_genesets = read.gmt('../source_data/h.all.v7.5.1.symbols.gmt')
hallmark_genes = do.call(c,hallmark_genesets)
hallmark_genes = hallmark_genes[!duplicated(hallmark_genes)]

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
  prefix = 'sample_cor_res_'
}else if(cor_method == 'pearson'){
  prefix = 'sample_pearson_cor_res_'
}

if(is.null(benchmarking_objs)){
  benchmarking_objs = list.files('../Benchmarking_obj/')
}

source('../scripts/evalu.R')
source('../scripts/helpers.R')

reference_mapping = read.delim('../Benchmarking_obj/reference_mapping.csv',sep = ',')

get_Z_sample_cor = function(obj_dir = './',
                            refType = 'self',
                            included_genes = NULL,
                            cor_method = 'spearman',
                            methods = c('bMIND',
                                        'cbsx',
                                        'ENIGMAtrace',
                                        'ENIGMAtrace(unnormalized)',
                                        'ENIGMAL2',
                                        'ENIGMAL2(unnormalized)',
                                        'EPICunmix',
                                        'InstaPrism',
                                        'TCA',
                                        'Unico')){
  
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
  
  
  bulk_Z = get_3d_bulk(bulk_expr,dimnames(Z_truth)[[3]])
  bulk_regressed_Z = get_3d_bulk(bulk_regressed,dimnames(Z_truth)[[3]])
  
  
  # only consider inter-sample correlation for provided gene list
  if(!is.null(included_genes)){
    cms = intersect(dimnames(Z_truth)[[1]],included_genes)
    Z_truth = Z_truth[cms,,]
    
    print(paste('only consider inter-sample gene correlation for',length(cms),'genes'))
  }
  
  print('start baseline correlation calculation ...')
  cor_res_list = list()
  cor_res_list[['bulk']] = get_per_gene_cor(bulk_Z,Z_truth,truthFrac,margin = 'column',cor_method = cor_method)
  cor_res_list[['bulk_regressed']] = get_per_gene_cor(bulk_regressed_Z,Z_truth,truthFrac,margin = 'column',cor_method = cor_method)
  
  for(method in methods){
    print(paste('Start correlation calculation for',method))
    
    Z_path = paste0(obj_dir,'/deconvRes/',refType,'_ref/',method,'_Z_inferred.RDS')
    if (!file.exists(Z_path)) {
      next  
    }
    Z_inferred <- readRDS(Z_path)
    
    # map cell_type names in Z_inferred
    if(refType == 'indep'){
      Z_inferred = Z_inferred[,,indep_mapping$maxCorName]
      dimnames(Z_inferred)[[3]] = indep_mapping$cell_type
    }
    
    cor_res_list[[method]] = get_per_gene_cor(Z_inferred,Z_truth,truthFrac,margin = 'column',cor_method = cor_method)
  }
  
  return(cor_res_list)
}

for(obj in benchmarking_objs){
  
  print(paste('################### Start Z sample cor calculation for',obj,'###################'))
  export_dir = paste0('../Benchmarking_obj/',obj,'/deconvSummary/')
  
  if(!dir.exists(export_dir)){
    dir.create(export_dir)
  }
  
  cor_res_list_indep = get_Z_sample_cor(obj_dir = paste0('../Benchmarking_obj/',obj),refType = 'indep',included_genes = hallmark_genes,cor_method = cor_method)
  saveRDS(cor_res_list_indep,file = paste0(export_dir,prefix,'indep_ref.RDS'))
}

print('done')

