#!/usr/bin/env Rscript
library(argparse)
library(dplyr)

# must run DALE_Z_cor.R first before running this script!
# note: by default will only consider n_vec = c(10,30,100,300,1000,3000,10000) and cutoff_vec = c(0,1,2,3,4,5)

parser <- ArgumentParser(add_help = F)
parser$add_argument(
  "--benchmarking_objs", '-o',
  type = "character",
  nargs = '+', 
  required = FALSE,
  help = 'Benchmarking obj name (accepts multiple values)'
)
parser$add_argument("--cor_method",'-r',type="character",required = F, default = 'spearman',choices = c("spearman", "pearson"))
parser$add_argument("--add_post_hoc_derived_summary", type = "integer", required = FALSE, default = 0,choices = c(0, 1),
                    help = "Use 1 to enable gene prioritization with post hoc derived statistics or 0 to skip (FALSE). Default 0.
                    If enabled, need to have DALE_post_hoc_analysis.R run first!")
reference_mapping = read.delim('../Benchmarking_obj/reference_mapping.csv',sep = ',')

args <- parser$parse_args()
benchmarking_objs = args$benchmarking_objs
cor_method = args$cor_method
add_post_hoc = as.logical(args$add_post_hoc_derived_summary)

if(cor_method == 'spearman'){
  export_name = 'by_limma_cor_summary.RDS'
}else if(cor_method == 'pearson'){
  export_name = 'by_limma_pearson_cor_summary.RDS'
}
if(is.null(benchmarking_objs)){
  benchmarking_objs = list.files('../Benchmarking_obj/')
}

source('../scripts/evalu.R')
source('../scripts/helpers.R')

get_by_limma_specificity_summary = function(obj_dir = './',
                                            indep_ref_dir = './', 
                                            n_vec = c(10,30,100,300,1000,3000,10000),
                                            cutoff_vec = c(0,1,2,3,4,5),
                                            cor_method = 'spearman',
                                            add_post_hoc = T){
  
  abbre1 = sub('_.*','',basename(normalizePath(obj_dir)))
  abbre2 = sub('_.*','',basename(normalizePath(indep_ref_dir)))
  
  stopifnot(abbre1 == abbre2 | (abbre1 == "LUAD" & abbre2 == "NSCLC"))
  
  if(cor_method == 'spearman'){
    cor_res_list_self = readRDS(paste0(obj_dir,'/deconvSummary/cor_res_self_ref.RDS'))
    cor_res_list_indep = readRDS(paste0(obj_dir,'/deconvSummary/cor_res_indep_ref.RDS'))
  }else if(cor_method == 'pearson'){
    cor_res_list_self = readRDS(paste0(obj_dir,'/deconvSummary/pearson_cor_res_self_ref.RDS'))
    cor_res_list_indep = readRDS(paste0(obj_dir,'/deconvSummary/pearson_cor_res_indep_ref.RDS'))
  }else{
    stop('invalid cor_method provided')
  }
  
  # 01. reference derived specificity
  # 01a. self ref
  self_limma_top_genes = read.delim(paste0(obj_dir,'/reference/self/limma_top_genes.csv'),sep = ',')
  summary1a = bind_rows(
    get_by_specificity_summary(cor_res_list_self, self_limma_top_genes, 'top_n', n_vec, cutoff_vec),
    get_by_specificity_summary(cor_res_list_self, self_limma_top_genes, 'cutoff', n_vec, cutoff_vec)
  ) %>%
    mutate(refType = 'self', specificity = 'reference')
  
  # 01b. indep ref
  indep_limma_top_genes = read.delim(paste0(indep_ref_dir,'/limma_top_genes.csv'),sep = ',')
  indep_mapping = read.delim(paste0(obj_dir,'/deconvRes/indep_ref/indep_ref_mapping.csv'),sep = ',')
  indep_limma_top_genes = indep_colnames_rename(indep_limma_top_genes,indep_mapping)
  # note: find overlapping genes first, then select top_n genes!
  indep_limma_top_genes = indep_limma_top_genes[rownames(indep_limma_top_genes) %in% rownames(cor_res_list_indep$bulk),]
  
  summary1b = bind_rows(
    get_by_specificity_summary(cor_res_list_indep,indep_limma_top_genes,'top_n', n_vec, cutoff_vec),
    get_by_specificity_summary(cor_res_list_indep,indep_limma_top_genes,'cutoff', n_vec, cutoff_vec)
  ) %>% 
    mutate(refType = 'indep', specificity = 'reference')
  
  # 02. post-hoc derived specificity (Z-mean logFC over second best)
  if(add_post_hoc){
    self_res = list()
    indep_res = list()
    for(method in methods){
      post_hoc_limma_self_path = paste0(obj_dir,'/deconvSummary/post_hoc_statisitics/self_ref/Z_mean_fc/',method,'_Z_mean_fc.csv')
      if(file.exists(post_hoc_limma_self_path)){
        Z_mean_fc = read.delim(post_hoc_limma_self_path,sep = ',')
        self_res[[method]] = bind_rows(
          get_by_specificity_summary(cor_res_list_self[method],Z_mean_fc,'top_n', n_vec, cutoff_vec),
          get_by_specificity_summary(cor_res_list_self[method],Z_mean_fc,'cutoff', n_vec, cutoff_vec)
        ) 
      }
      post_hoc_limma_indep_path = paste0(obj_dir,'/deconvSummary/post_hoc_statisitics/indep_ref/Z_mean_fc/',method,'_Z_mean_fc.csv')
      if(file.exists(post_hoc_limma_indep_path)){
        Z_mean_fc = read.delim(post_hoc_limma_indep_path,sep = ',')
        indep_res[[method]] = bind_rows(
          get_by_specificity_summary(cor_res_list_indep[method],Z_mean_fc,'top_n', n_vec, cutoff_vec) ,
          get_by_specificity_summary(cor_res_list_indep[method],Z_mean_fc,'cutoff', n_vec, cutoff_vec) 
        )
      }
    }
    
    summary2a = do.call(rbind,self_res) %>% mutate(refType = 'self', specificity = 'post_hoc_fc',rownames = NULL)
    summary2b = do.call(rbind,indep_res) %>% mutate(refType = 'indep', specificity = 'post_hoc_fc',rownames = NULL)
    
    by_limma_summary = do.call(rbind,list(summary1a,summary1b,
                                          summary2a,summary2b))    
  }else{
    by_limma_summary = do.call(rbind,list(summary1a,summary1b))  
  }

  rownames(by_limma_summary) = NULL
  
  truthFrac = read.delim(paste0(obj_dir,'/init/truthFrac.csv'),sep = ',')  %>% as.matrix()
  avg_frac = colMeans(truthFrac)
  avg_frac = avg_frac[order(avg_frac,decreasing = T)]
  by_limma_summary$cell_type = factor(by_limma_summary$cell_type,levels = names(avg_frac))
  
  return(by_limma_summary)
}


for(benchmarking_obj in benchmarking_objs){
  
  print(paste('################### start processing',benchmarking_obj,'###################'))
  
  indep_ref = reference_mapping$indep_ref[reference_mapping$datasets == benchmarking_obj]
  export_dir = paste0('../Benchmarking_obj/',benchmarking_obj,'/deconvSummary/')
  
  by_limma_top_genes_summary = get_by_limma_specificity_summary(obj_dir = paste0('../Benchmarking_obj/',benchmarking_obj),
                                                                indep_ref_dir = paste0('../Indep_scReference/',indep_ref),
                                                                n_vec = c(10,30,100,300,1000,3000,10000),
                                                                cutoff_vec = c(0,1,2,3,4,5),
                                                                cor_method = cor_method,
                                                                add_post_hoc = add_post_hoc)
  by_limma_top_genes_summary$dataset = benchmarking_obj
  
  saveRDS(by_limma_top_genes_summary,file = paste0(export_dir,export_name))
}

print('done')

