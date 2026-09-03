#!/usr/bin/env Rscript
# Expression specificity scored against Z-truth, as defined in the Methods: per gene,
# Lin's CCC between the inferred and the ground-truth log-scale mean cell type profiles.
#
# This complements DALE_Z_specificity.R, which reports the logFC contrast of each profile
# on its own. That contrast never looks at Z-truth, so a gene predicted to be sharply
# specific to the WRONG cell type still scores well on it; the CCC here penalises that.
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
parser$add_argument("--add_self_ref", type = "integer", required = FALSE, default = 0, choices = c(0, 1),
                    help = "Use 1 to enable CCC evaluation on deconvRes from self ref or 0 to skip (FALSE). Default 0")
parser$add_argument("--top_n", type = "integer", required = FALSE, default = 100,
                    help = "Restrict to the top n marker genes per cell type. Use 0 for all genes. Default 100")

args <- parser$parse_args()
benchmarking_objs = args$benchmarking_objs
add_self_ref = as.logical(args$add_self_ref)
top_n = args$top_n

if(is.null(benchmarking_objs)){
  benchmarking_objs = list.files('../Benchmarking_obj/')
}

source('../scripts/evalu.R')
source('../scripts/helpers.R')

# by default, will only consider these methods
methods = c('bMIND',
            'cbsx',
            'ENIGMAtrace',
            'ENIGMAtrace(unnormalized)',
            'ENIGMAL2',
            'ENIGMAL2(unnormalized)',
            'EPICunmix',
            'InstaPrism',
            'TCA',
            'Unico')

Z_scales = list()
Z_scales[['linear']] = c('BayesPrism','cbsx','ENIGMAL2','ENIGMAtrace','ENIGMAL2(unnormalized)','ENIGMAtrace(unnormalized)',
                         'Unico','InstaPrism','InstaPrismUpdated')
Z_scales[['log2']] = c('bMIND','EPICunmix','TCA','bMIND(EPICunmix)','EPICunmix(EPICunmix)')
Z_scales_vec = unlist(Z_scales)
names(Z_scales_vec) <- Z_scales_vec
Z_scales_vec[] <- rep(names(Z_scales), lengths(Z_scales))


for(benchmarking_obj in benchmarking_objs){

  print(paste('############### Start expression CCC analysis for',benchmarking_obj,'###############'))

  obj_dir = paste0('../Benchmarking_obj/',benchmarking_obj)
  export_dir = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvSummary/')

  if(!dir.exists(export_dir)){
    dir.create(export_dir,recursive = T)
  }

  Z_truth_limma_statistics = read.delim(paste0('../Benchmarking_obj/',benchmarking_obj,'/init/Z_limma_statistics.csv'),sep = ',')
  indep_ref_mapping = read.delim(paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/indep_ref/indep_ref_mapping.csv'),sep = ',')
  Z_truth = readRDS(paste0('../Benchmarking_obj/',benchmarking_obj,'/init/Z_truth.RDS'))

  if(top_n > 0){
    marker_genes = unique(unlist(get_marker_list(Z_truth_limma_statistics,select_method = 'top_n',n = top_n)))
  }else{
    marker_genes = NULL
  }

  expr_ccc_self = list()
  expr_ccc_indep = list()

  for(method in methods){
    print(method)

    if(add_self_ref){
      Z_self_path = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/self_ref/',method,'_Z_inferred.RDS')
      if(file.exists(Z_self_path)){
        Z = readRDS(Z_self_path)
        expr_ccc_self[[method]] = get_expr_ccc(Z,Z_truth,NULL,Z_scales_vec[method],'linear',marker_genes)
      }
      saveRDS(expr_ccc_self,file = paste0(export_dir,'/expr_ccc_self_ref.RDS'))
    }

    Z_indep_path = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvRes/indep_ref/',method,'_Z_inferred.RDS')
    if(file.exists(Z_indep_path)){
      Z = readRDS(Z_indep_path)
      expr_ccc_indep[[method]] = get_expr_ccc(Z,Z_truth,indep_ref_mapping,Z_scales_vec[method],'linear',marker_genes)
    }
  }

  saveRDS(expr_ccc_indep,file = paste0(export_dir,'/expr_ccc_indep_ref.RDS'))

  ccc_summary = do.call(rbind,lapply(names(expr_ccc_indep),function(m){
    data.frame(method = m,
               n_genes = nrow(expr_ccc_indep[[m]]),
               mean_ccc = mean(expr_ccc_indep[[m]]$ccc,na.rm = T),
               mean_pearson = mean(expr_ccc_indep[[m]]$pearson,na.rm = T),
               mean_Cb = mean(expr_ccc_indep[[m]]$Cb,na.rm = T))
  }))
  print(ccc_summary)
}

print('done')
