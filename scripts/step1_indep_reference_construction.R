# Canonical independent-reference construction workflow.
#
# Paired with scripts/step1_indep_reference_construction.ipynb.
# This file is intentionally linear, like scripts/step2_benchmarking_input.ipynb:
# each dataset has a divider and then the full R construction code. The older
# per-reference files under Indep_scReference/<ref>/ are kept unchanged for now.
#
# Usage from scripts/:
#   Rscript step1_indep_reference_construction.R
#
# All outputs are written explicitly to ../Indep_scReference/<ref>/.
# CIBERSORTx sections require credentials in the environment:
#   export CIBERSORTX_USERNAME='...'
#   export CIBERSORTX_TOKEN='...'

library(dplyr)

source("../DALE_Eval/modules/reference_prep_helpers.R")

get_cibersortx_username <- function() {
  value <- Sys.getenv("CIBERSORTX_USERNAME", unset = "")
  if (!nzchar(value)) {
    stop("CIBERSORTx username is required. Set CIBERSORTX_USERNAME before running CIBERSORTx reference sections.")
  }
  value
}

get_cibersortx_token <- function() {
  value <- Sys.getenv("CIBERSORTX_TOKEN", unset = "")
  if (!nzchar(value)) {
    stop("CIBERSORTx token is required. Set CIBERSORTX_TOKEN before running CIBERSORTx reference sections.")
  }
  value
}

configure_cibersortx_credentials <- function() {
  if (!requireNamespace("omnideconv", quietly = TRUE)) {
    stop("CIBERSORTx reference construction requires the omnideconv package.")
  }
  suppressPackageStartupMessages(library(omnideconv))
  set_cibersortx_credentials(get_cibersortx_username(), get_cibersortx_token())
}

validate_reference_outputs <- function(ref_name, required = c("refPhi.RDS", "cbsx_sig.txt", "bMIND_profile.csv", "rowMeans_sig.csv", "limma_top_genes.csv")) {
  ref_dir <- paste0("../Indep_scReference/", ref_name)
  status <- data.frame(
    indep_ref = ref_name,
    file = required,
    exists = file.exists(paste0(ref_dir, "/", required)),
    stringsAsFactors = FALSE
  )
  print(status)
  invisible(status)
}


############ BRCA_Wu2021 ############

message('=== BRCA_Wu2021 ===')

library(dplyr)
# source handled by scripts/step1_indep_reference_construction.R

seurat_obj = readRDS('../../../InstaPrismExtension/curated_dataset/BRCA/BRCA_SP1039/seurat_object_UMAP_withMeta.rds')

scExpr =  seurat_obj@assays[["RNA"]]@counts
scMeta = seurat_obj@meta.data
rm(seurat_obj)

# exclude ribosomal and mitochondrial genes
category.matrix <- BayesPrism:::assign.category(input.genes = rownames(scExpr), species= 'hs')
gene.group = c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY")
category.matrix <- category.matrix[, gene.group, drop=F]
print(colSums(category.matrix), na.rm=TRUE)

exclude.idx <- rowSums(category.matrix)>0	
scExpr = scExpr[!exclude.idx,]

dim(scExpr) # 28600 * 41514

cell.type.labels = scMeta$celltype_major
cell.type.labels[cell.type.labels=='Plasmablasts'] = 'B-cells'
cell.type.labels = gsub(' ','_',cell.type.labels)
cell.type.labels = gsub('-','_',cell.type.labels)



################# InstaPrism reference ###############
require(InstaPrism)
# build cell_type_labels
cell.type.labels = scMeta$celltype_major

# for malignant cells, find subclusters
scExpr_mal = scExpr[,cell.type.labels =='Cancer Epithelial']
scMeta_mal = scMeta[scMeta$celltype_major == 'Cancer Epithelial',]

mal_subcluster = deconvBenchmarking:::get_subcluster(scExpr_mal,scMeta_mal$celltype_minor,min.subcluster.size = 50)
scMeta_mal$mal_subcluster = mal_subcluster

# overwrite celltype_subset for malignant
for(i in 1:nrow(scMeta_mal)){
  matched_id = which(rownames(scMeta) == rownames(scMeta_mal)[i])
  scMeta$celltype_subset[matched_id] = scMeta_mal$mal_subcluster[i]
}

# merge plasmablast with B at cell type level
cell.type.labels[cell.type.labels=='Plasmablasts'] = 'B-cells'

brca_refPhi = refPrepare(scExpr,
                         cell.type.labels = cell.type.labels, 
                         cell.state.labels = scMeta$celltype_subset)

# avoid any blank in cell_type names
names(brca_refPhi@map) = gsub(' ','_',names(brca_refPhi@map))

# replace '-' with '_' in cell_type names
names(brca_refPhi@map) = gsub('-','_',names(brca_refPhi@map))


# ensure that each cell state only map to one cell type
length(do.call(c,brca_refPhi@map)) == ncol(brca_refPhi@phi.cs) # TRUE

saveRDS(brca_refPhi,file = '../Indep_scReference/BRCA_Wu2021/refPhi.RDS')

################ cell type marker genes from pseudobulk-DE analysis ##################
cpm_train = cpm_normalization(scExpr,target_sum = 1e4)
all.equal(colnames(cpm_train),rownames(scMeta)) # TRUE

cell.type.labels = scMeta$celltype_major
cell.type.labels[cell.type.labels=='Plasmablasts'] = 'B-cells'
cell.type.labels = gsub(' ','_',cell.type.labels)
cell.type.labels = gsub('-','_',cell.type.labels)

cell_types = unique(cell.type.labels)

ct_labels = c()
pseuodbulk_ct_all = data.frame(matrix(NA,ncol = 0,nrow = nrow(cpm_train)))

for(ct in cell_types){
  print(ct)
  ct_indices = which(cell.type.labels == ct)
  scExpr_sub = cpm_train[,ct_indices]
  pseudobulk_ct = get_mu_by_group(scExpr_sub,scMeta[ct_indices,]$orig.ident)
  
  pseudobulk_ct = as.data.frame(pseudobulk_ct)
  colnames(pseudobulk_ct) = paste0(as.character(ct),'_',colnames(pseudobulk_ct))
  ct_labels = c(ct_labels,rep(as.character(ct),ncol(pseudobulk_ct)))
  pseuodbulk_ct_all = cbind(pseuodbulk_ct_all,pseudobulk_ct)
}

#limma_statistics = deconvBenchmarking:::get_limma_statistics(pseuodbulk_ct_all, ct_labels, hv_genes)
#limma_top_genes = limma_statistics[[1]][["coefficients"]]

limma_top_genes = compute_limma_statistics(log2(pseuodbulk_ct_all+1),ct_labels)
write.table(limma_top_genes,file = '../Indep_scReference/BRCA_Wu2021/limma_top_genes.csv',sep = ',',row.names = T)

############### cbsx reference ###########
# will only consider top 10000 for cbsx reference construction
hv_genes = get_top_variable_genes(pseuodbulk_ct_all,top_n = 10000)

limma_top_genes = read.delim('../Indep_scReference/BRCA_Wu2021/limma_top_genes.csv',sep = ',')
hv_genes = rownames(limma_top_genes)

cell_type_annotations = cell.type.labels

# save(cpm_train,cell_type_annotations,file = '../Indep_scReference/BRCA_Wu2021/cbsx_required_input.RData')

# library(omnideconv) # must restart R session first!
# load('../Indep_scReference/BRCA_Wu2021/cbsx_required_input.RData')
configure_cibersortx_credentials() 
cbsx_sig = build_model_cibersortx(as.matrix(cpm_train[hv_genes,]),  # require a matrix format
                                  cell_type_annotations,
                                  container  = "singularity",
                                  container_path = '../../../deconvolution/modules/cibersortx/',
                                  verbose = T) 
# export cbsx_sig into required format
cbsx_sig_export=cbind(data.frame(GeneSymbol=rownames(cbsx_sig)),cbsx_sig)
write.table(cbsx_sig_export,file = '../Indep_scReference/BRCA_Wu2021/cbsx_sig.txt',sep = '\t',row.names = F,quote = F)


######### bNIND referernce ##########
EPICunmix_prior = EPICunmix::get_prior(sc = as.matrix(scExpr), # takes single-cell count matrix
                                       sample = scMeta$orig.ident,
                                       cell_type = cell.type.labels,
                                       filter_pd = T)

# 10073 genes available, while all other genes are filtered
# since this value is large enough, will export mean and covariance at the same time
bMIND_profile = EPICunmix_prior$profile 
bMIND_covariance = EPICunmix_prior$covariance

write.table(bMIND_profile,file = '../Indep_scReference/BRCA_Wu2021/bMIND_profile(depreicated).csv',sep = ',',row.names = T) # only 10073 genes available
saveRDS(bMIND_covariance,file = '../Indep_scReference/BRCA_Wu2021/bMIND_covariance(depreicated).RDS') # only 10073 genes available, deconvRes related is in bMIND(with_cov)_Z_inferred.RDS

# 0320 update: disable filter_pd to ensure enough genes in the bMIND reference
profile = log2(cpm_normalization(get_mu_by_group(scExpr,
                                                 cell.type.labels))+1) %>% as.matrix()
write.table(profile,file = '../Indep_scReference/BRCA_Wu2021/bMIND_profile.csv',sep = ',',row.names = T)

######### ENIGMA reference ###########
cpm_train = cpm_normalization(scExpr,target_sum = 1e6)

mu = get_mu_by_group(cpm_train,cell.type.labels)
write.table(mu,file = '../Indep_scReference/BRCA_Wu2021/rowMeans_sig.csv',sep = ',')

validate_reference_outputs('BRCA_Wu2021')


############ CRC_Lee2020 ############

message('=== CRC_Lee2020 ===')

library(dplyr)
library(tibble)
# source handled by scripts/step1_indep_reference_construction.R

dataset_path = '../../../InstaPrismExtension/curated_dataset/CRC/Data_Lee2020_Colorectal/'

scExpr = Seurat::ReadMtx(mtx = paste0(dataset_path,'Exp_data_UMIcounts.mtx'),
                         cells = paste0(dataset_path,'Cells.csv'),
                         features = paste0(dataset_path,'Genes.txt'),
                         feature.column = 1,skip.cell = 1,cell.sep = ',') 
scMeta = read.delim(paste0(dataset_path,'Cells.csv'),sep =',')
scMeta = scMeta %>% column_to_rownames('cell_name')
scMeta$nCount = Matrix::colSums(scExpr)
scMeta$nGenes = Matrix::colSums(scExpr>0)

print(dim(scExpr)) # 22276 * 21657

# remove cells with no annotations/ also remove 'Mast' since there's only 1 mast cell
keep_id = which(scMeta$cell_type!='' & scMeta$cell_type!='Mast' & scMeta$cell_type!= 'Epithelial'
                & scMeta$nCount >= 1000 & scMeta$nGenes >= 500) 
print(paste(length(keep_id)/nrow(scMeta),'of cells are kept')) # "0.946760862538671 of cells are kept"

scExpr = scExpr[,keep_id]
scMeta = scMeta[keep_id,]

# quality control of genes
v = Matrix::rowSums(scExpr > 0)
scExpr = scExpr[v >= 100,]

print(dim(scExpr)) # 14675 * 20504




cpm_train = cpm_normalization(scExpr,target_sum = 1e6)
save(cpm_train,scMeta,file = '../Indep_scReference/CRC_Lee2020/scRNA_training.RData')


################# InstaPrism reference ###############
require(InstaPrism)

# for malignant cells, use sampleID as cell_states
# for all other cells, use InstaPrism::get_subcluster() to find subclusters

scMeta$cell_state = NA
mal_id = which(scMeta$cell_type == 'Malignant')
scMeta$cell_state[mal_id] = paste0(scMeta$cell_type[mal_id],'_',scMeta$sample[mal_id])

recluster_id = which(is.na(scMeta$cell_state))

SCISSORS_clusters = InstaPrism::get_subcluster(scExpr = cpm_train[,recluster_id],
                                               cell_type_labels = scMeta$cell_type[recluster_id],
                                               subcluster_method = 'SCISSORS',
                                               CalculateSilhouette = F)
scMeta$cell_state[recluster_id] = SCISSORS_clusters
refPhi = InstaPrism::refPrepare(cpm_train,scMeta$cell_type,scMeta$cell_state)


save(scMeta,file = '../Indep_scReference/CRC_Lee2020/scMeta_annotated.RData')
saveRDS(refPhi,file = '../Indep_scReference/CRC_Lee2020/refPhi.RDS')

################ cell type marker genes from pseudobulk-DE analysis ##################
load('../Indep_scReference/CRC_Lee2020/scRNA_training.RData')
all.equal(colnames(cpm_train),rownames(scMeta)) # TRUE

cell.type.labels = scMeta$cell_type
cell_types = unique(cell.type.labels)

ct_labels = c()
pseuodbulk_ct_all = data.frame(matrix(NA,ncol = 0,nrow = nrow(cpm_train)))

for(ct in cell_types){
  print(ct)
  ct_indices = which(cell.type.labels == ct)
  scExpr_sub = cpm_train[,ct_indices]
  pseudobulk_ct = get_mu_by_group(scExpr_sub,scMeta[ct_indices,]$sample)
  
  pseudobulk_ct = as.data.frame(pseudobulk_ct)
  colnames(pseudobulk_ct) = paste0(as.character(ct),'_',colnames(pseudobulk_ct))
  ct_labels = c(ct_labels,rep(as.character(ct),ncol(pseudobulk_ct)))
  pseuodbulk_ct_all = cbind(pseuodbulk_ct_all,pseudobulk_ct)
}


#limma_statistics = deconvBenchmarking:::get_limma_statistics(pseuodbulk_ct_all, ct_labels, hv_genes)
#limma_top_genes = limma_statistics[[1]][["coefficients"]]
limma_top_genes = compute_limma_statistics(log2(pseuodbulk_ct_all+1),ct_labels)

write.table(limma_top_genes,file = '../Indep_scReference/CRC_Lee2020/limma_top_genes.csv',sep = ',',row.names = T)


############### cbsx reference ###########
# will only consider top 10000 for cbsx reference construction
#limma_top_genes = read.delim('../Indep_scReference/CRC_Lee2020/limma_top_genes.csv',sep = ',')
#hv_genes = rownames(limma_top_genes)
hv_genes = get_top_variable_genes(pseuodbulk_ct_all,top_n = 10000)

# library(omnideconv) # must restart R session first!
configure_cibersortx_credentials() 
cbsx_sig = build_model_cibersortx(as.matrix(cpm_train[hv_genes,]),  # require a matrix format
                                  scMeta$cell_type,
                                  container  = "singularity",
                                  container_path = '../../../deconvolution/modules/cibersortx/',
                                  verbose = T) 
# export cbsx_sig into required format
cbsx_sig_export=cbind(data.frame(GeneSymbol=rownames(cbsx_sig)),cbsx_sig)
write.table(cbsx_sig_export,file = '../Indep_scReference/CRC_Lee2020/cbsx_sig.txt',sep = '\t',row.names = F,quote = F)

######### bNIND referernce ##########
# 1431 genes are filtered out because cell-type covariance matrix is not positive-definite (PD)
EPICunmix_prior = EPICunmix::get_prior(sc = as.matrix(scExpr), # takes single-cell count matrix
                                       sample = scMeta$sample,
                                       cell_type = scMeta$cell_type,
                                       filter_pd = T)

# 13244 genes available; since this value is large enough, will export mean and covariance at the same time
bMIND_profile = EPICunmix_prior$profile 
bMIND_covariance = EPICunmix_prior$covariance

write.table(bMIND_profile,file = '../Indep_scReference/CRC_Lee2020/bMIND_profile.csv',sep = ',',row.names = T)
saveRDS(bMIND_covariance,file = '../Indep_scReference/CRC_Lee2020/bMIND_covariance.RDS')

######### ENIGMA reference ###########
mu = get_mu_by_group(cpm_train,scMeta$cell_type)
write.table(mu,file = '../Indep_scReference/CRC_Lee2020/rowMeans_sig.csv',sep = ',')

validate_reference_outputs('CRC_Lee2020')


############ LUAD_Laughney2020 ############

message('=== LUAD_Laughney2020 ===')

library(dplyr)
library(tibble)
# source handled by scripts/step1_indep_reference_construction.R

# note: Laughney_Massague_2020 from Salcher2022 only contains 26403 cells! but from 3ca resource we have 40505 cells
# therefore we will only consider scRNA data from 3ca

ad = anndata::read_h5ad('../scRNA_datasets/LUAD_Laughney2020/LUAD_Laughney2020_processed.h5ad')
patient_info = read.delim('../../../InstaPrismExtension/curated_dataset/LUAD/Data_Laughney2020/Meta-data.csv',sep = ',')

scExpr = Matrix::t(ad$X) # 19222* 40505
scMeta = ad$obs

scMeta$site = patient_info$site[match(scMeta$sample,patient_info$sample)]

# exclude Neutrophil
keep_id = which(scMeta$cell_type!='Neutrophil'
                & scMeta$total_counts >= 1000 & scMeta$n_genes_by_counts >= 500)

scExpr = scExpr[,keep_id] 
scMeta = scMeta[keep_id,]

# quality control of genes
v = Matrix::rowSums(scExpr > 0)
scExpr = scExpr[v >= 100,]

print(dim(scExpr)) # 15508 * 39614

cpm_train = cpm_normalization(scExpr,target_sum = 1e6)


save(cpm_train,scMeta,file = '../Indep_scReference/LUAD_Laughney2020/scRNA_training.RData')

################# InstaPrism reference ###############
require(InstaPrism)

scMeta$site[scMeta$site == 'Primary '] = 'Primary'
scMeta$site[grepl('M\\(',scMeta$site)] = 'Metastasis'

scMeta$cell_state = paste0(scMeta$cell_subtype,'_',scMeta$site)

mal_id = which(scMeta$cell_type == 'Malignant')
scMeta$cell_state[mal_id] = paste0(scMeta$cell_type[mal_id],'_',scMeta$sample[mal_id])

refPhi = InstaPrism::refPrepare(cpm_train,scMeta$cell_type,scMeta$cell_state)

all.equal(length(do.call(c,refPhi@map)),ncol(refPhi@phi.cs)) # TRUE
saveRDS(refPhi,file = '../Indep_scReference/LUAD_Laughney2020/refPhi.RDS')

################# bMIND reference #################
EPICunmix_prior = EPICunmix::get_prior(sc = as.matrix(scExpr), # takes single-cell count matrix
                                       sample = scMeta$sample,
                                       cell_type = scMeta$cell_type,
                                       filter_pd = T) # unable to calculate covariance: Error in if (!is.symmetric.matrix(x)) stop("argument x is not a symmetric matrix")

# will manually calculate profile only (using code from EPICunmix::get_prior)
sc = as.matrix(scExpr)
K = length(sort(unique(scMeta$cell_type)))
cell_type = sort(unique(scMeta$cell_type))
profile = matrix(NA, nrow(sc), K)
rownames(profile) = rownames(sc)
colnames(profile) = cell_type
for (i in cell_type) {
  profile[, i] = log2(edgeR::cpm(rowMeans(sc[, scMeta$cell_type == 
                                               i])) + 1)
}

write.table(profile,file = '../Indep_scReference/LUAD_Laughney2020/bMIND_profile.csv',sep = ',',row.names = T)

# or equivalently
profile = log2(cpm_normalization(get_mu_by_group(scExpr,scMeta$cell_type))+1) %>% as.matrix()


######### ENIGMA reference ###########
mu = get_mu_by_group(cpm_train,scMeta$cell_type)
write.table(mu,file = '../Indep_scReference/LUAD_Laughney2020/rowMeans_sig.csv',sep = ',')

################ cell type marker genes from pseudobulk-DE analysis ##################
all.equal(colnames(cpm_train),rownames(scMeta)) # TRUE

cell.type.labels = scMeta$cell_type
cell_types = unique(cell.type.labels)

ct_labels = c()
pseuodbulk_ct_all = data.frame(matrix(NA,ncol = 0,nrow = nrow(cpm_train)))

for(ct in cell_types){
  print(ct)
  ct_indices = which(cell.type.labels == ct)
  scExpr_sub = cpm_train[,ct_indices]
  pseudobulk_ct = get_mu_by_group(scExpr_sub,scMeta[ct_indices,]$sample)
  
  pseudobulk_ct = as.data.frame(pseudobulk_ct)
  colnames(pseudobulk_ct) = paste0(as.character(ct),'_',colnames(pseudobulk_ct))
  ct_labels = c(ct_labels,rep(as.character(ct),ncol(pseudobulk_ct)))
  pseuodbulk_ct_all = cbind(pseuodbulk_ct_all,pseudobulk_ct)
}


limma_top_genes = compute_limma_statistics(log2(pseuodbulk_ct_all+1),ct_labels)
write.table(limma_top_genes,file = '../Indep_scReference/LUAD_Laughney2020/limma_top_genes.csv',sep = ',',row.names = T)


#### cbsx signatures ####
# will only consider top 10000 for cbsx reference construction
#limma_top_genes = read.delim('../Indep_scReference/LUAD_Laughney2020/limma_top_genes.csv',sep = ',')
#hv_genes = rownames(limma_top_genes)

hv_genes = get_top_variable_genes(pseuodbulk_ct_all,top_n = 10000)

# library(omnideconv) # must restart R session first!
load('../Indep_scReference/LUAD_Laughney2020/scRNA_training.RData',verbose = T)
configure_cibersortx_credentials() 
# cbsx takes CPM as input
cbsx_sig = build_model_cibersortx(as.matrix(cpm_train[hv_genes,]), 
                                  as.character(scMeta$cell_type),
                                  container  = "singularity",
                                  container_path = '../../../deconvolution/modules/cibersortx/',
                                  verbose = T) 

# export cbsx_sig into required format
cbsx_sig_export=cbind(data.frame(GeneSymbol=rownames(cbsx_sig)),cbsx_sig)
write.table(cbsx_sig_export,file = '../Indep_scReference/LUAD_Laughney2020/cbsx_sig.txt',sep = '\t',row.names = F,quote = F)

validate_reference_outputs('LUAD_Laughney2020')


############ PBMC_AIDA2024 / PBMC_refined_AIDA2024 ############

message('=== PBMC_AIDA2024 / PBMC_refined_AIDA2024 ===')

library(dplyr)

# source handled by scripts/step1_indep_reference_construction.R

# AnnData object with n_obs × n_vars = 1058909 × 36266
ad = anndata::read_h5ad('../../../InstaPrismExtension/scRNA_resource/Cellxgene/PBMC_AIDA/85ca63ad-10c2-4c9e-9f76-a52568b294f1.h5ad')
meta = ad$obs
var = ad$var
obsm = ad$obsm

# subset to donor_id from South_Korea only for reference construction!

ad <- ad[ad$obs$Country == "South_Korea", ]  # 386150 × 36266
var = ad$var
meta = ad$obs

max(ad$X) # 8.988204

# meta = readRDS('../../../InstaPrismExtension/scRNA_resource/Cellxgene/PBMC_AIDA/PBMC_meta.RDS')
# var = readRDS('../../../InstaPrismExtension/scRNA_resource/Cellxgene/PBMC_AIDA/var.RDS')

# gene filtering (keep genes expressed in at least 100 cells)
nCells = Matrix::colSums(ad$X > 0)

keep_gene_id = which(nCells >= 100)
sum(duplicated(var$feature_name[keep_gene_id])) # 0

table(meta$author_cell_type)

#B          CD4+_T       CD4+_T_cm      CD4+_T_cyt       CD4+_T_em 
#2862            7002           38777            5299            8383 
#CD4+_T_naive          CD8+_T    CD8+_T_GZMB+    CD8+_T_GZMK+    CD8+_T_naive 
#51057             728           22164           13889           26733 
#CD14+_Monocyte  CD16+_Monocyte        CD16+_NK        CD56+_NK              DC 
#57783           14072           53512            2163              63 
#IGHMhi_memory_B IGHMlo_memory_B             ILC            MAIT        Monocyte 
#5058            6399              63            5699            7760 
#NK        Plasma_B        Platelet             RBC               T 
#5147             574            8372             104           14576 
#Treg      atypical_B             cDC            cDC1            cDC2 
#3279            1639             186             224            5490 
#dnT             gdT         naive_B             pDC 
#237            7802            7278            1776 


to_remove_ct = c('ILC','Platelet','RBC','pDC')
all_ct = unique(meta$author_cell_type)
to_keep_ct = setdiff(all_ct,to_remove_ct) 


to_keep_ct
#[1] "CD14+_Monocyte"  "T"               "CD16+_NK"        "Monocyte"       
#[5] "Treg"            "CD4+_T_cm"       "CD8+_T_naive"    "CD4+_T_naive"   
#[9] "B"               "IGHMhi_memory_B" "CD16+_Monocyte"  "cDC1"           
#[13] "CD8+_T_GZMK+"    "gdT"             "IGHMlo_memory_B" "CD8+_T_GZMB+"   
#[17] "CD4+_T_em"       "naive_B"         "NK"              "MAIT"           
#[21] "cDC2"            "CD4+_T_cyt"      "CD4+_T"          "atypical_B"     
#[25] "Plasma_B"        "dnT"             "CD56+_NK"        "CD8+_T"         
#[29] "cDC"             "DC"             


keep_id = which(meta$nCount_RNA >= 1000 & meta$nFeature_RNA >=500 & meta$author_cell_type %in% to_keep_ct) # 375621 cells kept
ad_raw = ad
ad = ad[keep_id,keep_gene_id] # 375621 × 20738

X_original <- exp(ad$X) - 1
count = Matrix::t(X_original)
rownames(count) = ad$var$feature_name

count[1:5,1:2]
#5 x 2 Matrix of class "dgeMatrix"
#CTGATAGAGTACGACG-KR_B2_L1 CACATTTCAGCTGGCT-KR_B2_L1
#RP11-34P13.7                          0                         0
#RP11-34P13.13                         0                         0
#LINC01409                             0                         0
#FAM87B                                0                         0
#LINC01128                             0                         0

scMeta = ad$obs
all.equal(rownames(scMeta) %>% as.character(),colnames(count)%>% as.character()) # TRUE

scMeta$author_cell_type = as.vector(scMeta$author_cell_type)
scMeta$donor_id = as.vector(scMeta$donor_id)

length(unique(scMeta$donor_id)) # 170



################## BLUE/scTAPE processed h5ad export ##############
# The PBMC BLUE/scTAPE processed h5ad export is intentionally maintained in
# scripts/step1_indep_reference_construction.ipynb as a Python-only cell.
#
# Reason: the PBMC source h5ad is loaded through Python/anndata, and converting
# large sparse/log matrices through R/reticulate can produce PyCapsule or object
# dtype conversion errors. The notebook cell builds the AnnData object directly
# in Python and writes:
#   ../scRNA_datasets/PBMC_AIDA2024/PBMC_AIDA2024_processed.h5ad
#   ../Indep_scReference/PBMC_AIDA2024/celltype_mapping.yaml
#   ../Indep_scReference/PBMC_refined_AIDA2024/celltype_mapping.yaml
#
# Run that notebook cell before BLUE/scTAPE PBMC runs. The remaining R sections
# below still generate the non-BLUE/scTAPE independent-reference files.

ct_recode <- function(x, mapping_list) {
  lookup_table <- unlist(mapping_list)
  names(lookup_table) <- rep(names(mapping_list), lengths(mapping_list))
  recoded <- names(lookup_table)[match(x, lookup_table)]
  return(recoded)
}

################## PBMC_AIDA2024 ##############
if(!dir.exists('../Indep_scReference/PBMC_AIDA2024/')){
  dir.create('../Indep_scReference/PBMC_AIDA2024/')
}

ct_map_coarse = list(
  T_cell = c('CD4+_T', 'CD4+_T_cm','CD4+_T_cyt','CD4+_T_em','CD4+_T_naive',
             'CD8+_T','CD8+_T_GZMB+','CD8+_T_GZMK+','CD8+_T_naive','MAIT',
             'dnT','gdT','T','Treg'),
  B_cell = c('atypical_B','B','IGHMhi_memory_B','IGHMlo_memory_B','naive_B','Plasma_B'),
  myeloid = c('CD14+_Monocyte','CD16+_Monocyte','Monocyte',
              'cDC','cDC1','cDC2','DC'),
  NK = c('CD16+_NK','CD56+_NK','NK')
)

length(do.call(c,ct_map_coarse))
length(unique(scMeta$author_cell_type))

scMeta$coarse_cell_type = ct_recode(scMeta$author_cell_type,ct_map_coarse)

### 01. bMIND reference ###
profile = log2(cpm_normalization(get_mu_by_group(count,scMeta$coarse_cell_type))+1) %>% as.matrix()
write.table(profile,file = '../Indep_scReference/PBMC_AIDA2024/bMIND_profile.csv',sep = ',',row.names = T)

### 02. ENIGMA reference ###
cpm = cpm_normalization(count)
mu = get_mu_by_group(cpm,scMeta$coarse_cell_type)
write.table(mu,file = '../Indep_scReference/PBMC_AIDA2024/rowMeans_sig.csv',sep = ',')

### 03. cell type marker genes from pseudobulk-DE analysis ###
cpm = cpm_normalization(count)
all.equal(colnames(cpm) %>% as.character(),rownames(scMeta)%>% as.character()) # TRUE

cell.type.labels = scMeta$coarse_cell_type
cell_types = unique(cell.type.labels)

ct_labels = c()
pseuodbulk_ct_all = data.frame(matrix(NA,ncol = 0,nrow = nrow(cpm)))


for(ct in cell_types){
  print(ct)
  ct_indices = which(cell.type.labels == ct)
  scExpr_sub = cpm[,ct_indices]
  pseudobulk_ct = get_mu_by_group(scExpr_sub,scMeta[ct_indices,]$donor_id)
  
  pseudobulk_ct = as.data.frame(pseudobulk_ct)
  colnames(pseudobulk_ct) = paste0(as.character(ct),'_',colnames(pseudobulk_ct))
  ct_labels = c(ct_labels,rep(as.character(ct),ncol(pseudobulk_ct)))
  pseuodbulk_ct_all = cbind(pseuodbulk_ct_all,pseudobulk_ct)
}


limma_top_genes = compute_limma_statistics(log2(pseuodbulk_ct_all+1),ct_labels)
# rownames(limma_top_genes) = var$feature_name[match(rownames(limma_top_genes),rownames(var))]

write.table(limma_top_genes,file = '../Indep_scReference/PBMC_AIDA2024/limma_top_genes.csv',sep = ',',row.names = T)

### 04. InstaPrism reference ###
refPhi = InstaPrism::refPrepare(cpm,scMeta$coarse_cell_type,scMeta$author_cell_type)
all.equal(length(do.call(c,refPhi@map)),ncol(refPhi@phi.cs)) 
saveRDS(refPhi,file = '../Indep_scReference/PBMC_AIDA2024/refPhi.RDS')

### 05. cbsx reference ###
# random subset to around 50,000 cells for reference construction
set.seed(123)

hv_genes = get_top_variable_genes(pseuodbulk_ct_all,top_n = 10000)
cbsx_id = sample(nrow(scMeta),50000,replace = F)
cbsx_cpm = cpm[hv_genes,cbsx_id]
cell_type_anno = scMeta$coarse_cell_type[cbsx_id]
save(cbsx_cpm,cell_type_anno,file = '../Indep_scReference/PBMC_AIDA2024/cbsx_required_input.RData')


library(omnideconv) # must restart R session first!
load('../Indep_scReference/PBMC_AIDA2024/cbsx_required_input.RData',verbose = T)

random_id = sample(1:50000, 45000,replace = F)

configure_cibersortx_credentials() 
# cbsx takes CPM as input
cbsx_sig = build_model_cibersortx(as.matrix(cbsx_cpm[,random_id]), 
                                  cell_type_anno[random_id],
                                  container  = "singularity",
                                  container_path = '../../../deconvolution/modules/cibersortx/',
                                  verbose = T) 
cbsx_sig_export=cbind(data.frame(GeneSymbol=rownames(cbsx_sig)),cbsx_sig)

write.table(cbsx_sig_export,file = '../Indep_scReference/PBMC_AIDA2024/cbsx_sig.txt',sep = '\t',row.names = F,quote = F)

################## PBMC_refined_AIDA2024 ##############
# a refined PBMC reference 

ct_map_refined = list(
  CD4_T = c('CD4+_T', 'CD4+_T_cm','CD4+_T_cyt','CD4+_T_em','CD4+_T_naive','Treg'),
  CD8_T = c('CD8+_T','CD8+_T_GZMB+','CD8+_T_GZMK+','CD8+_T_naive','MAIT'),
  gd_T = 'gdT',
  B_cell = c('atypical_B','B','IGHMhi_memory_B','IGHMlo_memory_B','naive_B'),
  Plasma_B = 'Plasma_B',
  CD14_Mono = 'CD14+_Monocyte',
  CD16_Mono = 'CD16+_Monocyte',
  DC = c('cDC','cDC1','cDC2','DC'),
  NK = c('CD16+_NK','CD56+_NK','NK')
)


refined_id = which(scMeta$author_cell_type %in% do.call(c,ct_map_refined)) 
scMeta_refined = scMeta[refined_id,]
count_refined = count[,refined_id]
cpm_refined = cpm[,refined_id]


scMeta_refined$refined_cell_type = ct_recode(scMeta_refined$author_cell_type,ct_map_refined)


### 01. bMIND reference ###
profile = log2(cpm_normalization(get_mu_by_group(count_refined,scMeta_refined$refined_cell_type))+1) %>% as.matrix()
write.table(profile,file = '../Indep_scReference/PBMC_refined_AIDA2024/bMIND_profile.csv',sep = ',',row.names = T)

### 02. ENIGMA reference ###
mu = get_mu_by_group(cpm_refined,scMeta_refined$refined_cell_type)
write.table(mu,file = '../Indep_scReference/PBMC_refined_AIDA2024/rowMeans_sig.csv',sep = ',')

### 03. cell type marker genes from pseudobulk-DE analysis ###
all.equal(colnames(cpm_refined) %>% as.character(),rownames(scMeta_refined) %>% as.character()) 

cell.type.labels = scMeta_refined$refined_cell_type
cell_types = unique(cell.type.labels)

ct_labels = c()
pseuodbulk_ct_all = data.frame(matrix(NA,ncol = 0,nrow = nrow(cpm_refined)))

for(ct in cell_types){
  print(ct)
  ct_indices = which(cell.type.labels == ct)
  scExpr_sub = cpm_refined[,ct_indices]
  pseudobulk_ct = get_mu_by_group(scExpr_sub,scMeta_refined[ct_indices,]$donor_id)
  
  pseudobulk_ct = as.data.frame(pseudobulk_ct)
  colnames(pseudobulk_ct) = paste0(as.character(ct),'_',colnames(pseudobulk_ct))
  ct_labels = c(ct_labels,rep(as.character(ct),ncol(pseudobulk_ct)))
  pseuodbulk_ct_all = cbind(pseuodbulk_ct_all,pseudobulk_ct)
}


limma_top_genes = compute_limma_statistics(log2(pseuodbulk_ct_all+1),ct_labels)
write.table(limma_top_genes,file = '../Indep_scReference/PBMC_refined_AIDA2024/limma_top_genes.csv',sep = ',',row.names = T)


### 04. InstaPrism reference ###
refPhi = InstaPrism::refPrepare(cpm_refined,scMeta_refined$refined_cell_type,scMeta_refined$author_cell_type)
all.equal(length(do.call(c,refPhi@map)),ncol(refPhi@phi.cs)) 
saveRDS(refPhi,file = '../Indep_scReference/PBMC_refined_AIDA2024/refPhi.RDS')

### 05. cbsx reference ###

cell_types <- scMeta_refined$refined_cell_type
cell_barcodes <- rownames(scMeta_refined)
cell_type_counts <- table(cell_types)

#cell_type_counts
#cell_types
#B_cell CD14_Mono CD16_Mono     CD4_T     CD8_T        DC      gd_T        NK
#23236     57643     14063    113797     69213      5963      7802     60822
#Plasma_B
#574
set.seed(123) 

cells_to_keep <- unlist(lapply(names(cell_type_counts), function(ct) {
  cells_of_type <- cell_barcodes[cell_types == ct]
  n_cells <- length(cells_of_type)
  
  if (n_cells < 1000) {
    return(cells_of_type)
  } else {
    n_downsample <- round(n_cells * 0.12)
    return(sample(cells_of_type, n_downsample))
  }
})) # 42880 cells kept

cpm_cbsx = cpm_refined[hv_genes,cells_to_keep]

cell_type_anno = scMeta_refined[cells_to_keep,'refined_cell_type']
save(cpm_cbsx,cell_type_anno,file = '../Indep_scReference/PBMC_refined_AIDA2024/cbsx_required_input.RData')

load('../Indep_scReference/PBMC_refined_AIDA2024/cbsx_required_input.RData',verbose = T)
# library(omnideconv) # must restart R session first!

configure_cibersortx_credentials() 
# cbsx takes CPM as input
cbsx_sig = build_model_cibersortx(as.matrix(cpm_cbsx), 
                                  cell_type_anno,
                                  container  = "singularity",
                                  container_path = '../../../deconvolution/modules/cibersortx/',
                                  verbose = T) 
cbsx_sig_export=cbind(data.frame(GeneSymbol=rownames(cbsx_sig)),cbsx_sig)

write.table(cbsx_sig_export,file = '../Indep_scReference/PBMC_refined_AIDA2024/cbsx_sig.txt',sep = '\t',row.names = F,quote = F)

validate_reference_outputs('PBMC_AIDA2024')
validate_reference_outputs('PBMC_refined_AIDA2024')


############ ROSMAP_MultiregionAD_Mathys2024 ############

message('=== ROSMAP_MultiregionAD_Mathys2024 ===')

library(Seurat)
# source handled by scripts/step1_indep_reference_construction.R

# downloaded from https://www.synapse.org/#!Synapse:syn52293442 --> https://www.synapse.org/Synapse:syn52383412
# check umap obj from scRNA_datasets/ROSMAP_MultiregionAD_Mathys2024

seurat_obj = readRDS('../../../InstaPrismExtension/scRNA_resource/ROSMAP/MultiregionAD_Mathys2024/Prefrontal_cortex.rds') # 33538*254721

table(seurat_obj$major_cell_type)
#Ast    CAMs    CPEC     End     Epd     Exc     Fib     Inh     Mic     Oli 
#17004     495       1     865       1  112143     635   40290    9263   63192 
#OPC     Per     SMC T cells 
#9781     614     140     297 

length(unique(seurat_obj$projid)) # 48

# keep 7 major cell types only for reference construction
to_keep_ct = c('Ast',
               'OPC',
               'Exc',
               'Inh',
               'Mic','T cells', 'CAMs',
               'Oli',
               'Per','SMC','End','Fib')

seurat_obj_main <- subset(seurat_obj, subset = major_cell_type %in% to_keep_ct) 
seurat_obj_main <- subset(seurat_obj_main, features = rownames(seurat_obj_main)[Matrix::rowSums(seurat_obj_main@assays$RNA@counts > 0) >= 1000]) # 20997 * 254719


scMeta = seurat_obj_main@meta.data

ct_map = list(
  Excitatory_neuron = 'Exc',
  Inhibitory_neuron = 'Inh',
  Immune = c('Mic','T cells','CAMs'),
  Astrocyte = 'Ast',
  Oligodendrocyte_precursor_cell = 'OPC',
  Oligodendrocyte = 'Oli',
  Vascular = c('Per','SMC','End','Fib')
)

x = unique(scMeta$major_cell_type)
x[!x %in% do.call(c,ct_map)]

ct_recode <- function(x, mapping_list) {
  lookup_table <- unlist(mapping_list)
  names(lookup_table) <- rep(names(mapping_list), lengths(mapping_list))
  recoded <- names(lookup_table)[match(x, lookup_table)]
  return(recoded)
}

scMeta$cell_type = ct_recode(scMeta$major_cell_type,ct_map)
scMeta$sample = paste0('sample_', scMeta$projid)

table(scMeta$cell_type)

################## BLUE/scTAPE processed h5ad export ##############
# BLUE and scTAPE need a processed h5ad with count-like X, gene symbols in
# var_names, obs['cell_type'] as the fine labels used by the mapping YAML, and
# obs['sample'] as the library/sample label.
if (!dir.exists('../scRNA_datasets/ROSMAP_MultiregionAD_Mathys2024/')) {
  dir.create('../scRNA_datasets/ROSMAP_MultiregionAD_Mathys2024/', recursive = TRUE)
}

scExpr_for_h5ad = seurat_obj_main@assays$RNA@counts
obs_for_h5ad = scMeta
var_for_h5ad = data.frame(gene_symbol = rownames(scExpr_for_h5ad))
rownames(var_for_h5ad) = rownames(scExpr_for_h5ad)

ad_rosmap = anndata::AnnData(
  X = Matrix::t(scExpr_for_h5ad),
  obs = obs_for_h5ad,
  var = var_for_h5ad
)
ad_rosmap$write_h5ad('../scRNA_datasets/ROSMAP_MultiregionAD_Mathys2024/ROSMAP_MultiregionAD_Mathys2024_processed.h5ad')
rm(ad_rosmap, scExpr_for_h5ad, obs_for_h5ad, var_for_h5ad)

################## BLUE/scTAPE celltype mapping ##############
# Keys are BLUE/scTAPE output classes. Values are labels in obs['cell_type'].
writeLines(c(
  'Excitatory_neuron:',
  '  - Excitatory_neuron',
  'Inhibitory_neuron:',
  '  - Inhibitory_neuron',
  'Immune:',
  '  - Immune',
  'Astrocyte:',
  '  - Astrocyte',
  'Oligodendrocyte_precursor_cell:',
  '  - Oligodendrocyte_precursor_cell',
  'Oligodendrocyte:',
  '  - Oligodendrocyte',
  'Vascular:',
  '  - Vascular'
), con = '../Indep_scReference/ROSMAP_MultiregionAD_Mathys2024/celltype_mapping.yaml')

#Astrocyte              Excitatory_neuron 
#17004                         112143 
#Immune              Inhibitory_neuron 
#10055                          40290 
#Oligodendrocyte Oligodendrocyte_precursor_cell 
#63192                           9781 
#Vascular 
#2254 

################# bMIND reference #################
scExpr = seurat_obj_main@assays$RNA@counts
# equivalent to profile[, i] = log2(edgeR::cpm(rowMeans(sc[,meta_sc$cell_type == i])) + 1) from EPICunmix::get_prior()
profile = log2(cpm_normalization(get_mu_by_group(scExpr,scMeta$cell_type))+1) %>% as.matrix()
write.table(profile,file = '../Indep_scReference/ROSMAP_MultiregionAD_Mathys2024/bMIND_profile.csv',sep = ',',row.names = T)

################ ENIGMA reference ##############
cpm = cpm_normalization(scExpr)
mu = get_mu_by_group(cpm,scMeta$cell_type)
write.table(mu,file = '../Indep_scReference/ROSMAP_MultiregionAD_Mathys2024/rowMeans_sig.csv',sep = ',')

################ cell type marker genes from pseudobulk-DE analysis ##################
scExpr = seurat_obj_main@assays$RNA@counts
cpm = cpm_normalization(scExpr)
all.equal(colnames(cpm),rownames(scMeta)) # TRUE

cell.type.labels = scMeta$cell_type
cell_types = unique(cell.type.labels)

ct_labels = c()
pseuodbulk_ct_all = data.frame(matrix(NA,ncol = 0,nrow = nrow(cpm)))

for(ct in cell_types){
  print(ct)
  ct_indices = which(cell.type.labels == ct)
  scExpr_sub = cpm[,ct_indices]
  pseudobulk_ct = get_mu_by_group(scExpr_sub,scMeta[ct_indices,]$sample)
  
  pseudobulk_ct = as.data.frame(pseudobulk_ct)
  colnames(pseudobulk_ct) = paste0(as.character(ct),'_',colnames(pseudobulk_ct))
  ct_labels = c(ct_labels,rep(as.character(ct),ncol(pseudobulk_ct)))
  pseuodbulk_ct_all = cbind(pseuodbulk_ct_all,pseudobulk_ct)
}

limma_top_genes = compute_limma_statistics(log2(pseuodbulk_ct_all+1),ct_labels)
write.table(limma_top_genes,file = '../Indep_scReference/ROSMAP_MultiregionAD_Mathys2024/limma_top_genes.csv',sep = ',',row.names = T)

#### cbsx signatures ####
hv_genes = get_top_variable_genes(pseuodbulk_ct_all,top_n = 10000)

# subset to around 1/10 for cbsx reference construction
cell_types <- scMeta$cell_type
cell_barcodes <- colnames(seurat_obj_main)
cell_type_counts <- table(cell_types)
set.seed(123) 
cells_to_keep <- unlist(lapply(names(cell_type_counts), function(ct) {
  cells_of_type <- cell_barcodes[cell_types == ct]
  n_cells <- length(cells_of_type)
  
  if (n_cells < 1000) {
    return(cells_of_type)
  } else {
    n_downsample <- round(n_cells * 0.18)
    return(sample(cells_of_type, n_downsample))
  }
}))

length(cells_to_keep) # 45851


cpm_cbsx = cpm[hv_genes,cells_to_keep] # 10000 * 45851
cell_type_labels = scMeta[cells_to_keep,'cell_type']
table(cell_type_labels)
#Astrocyte              Excitatory_neuron 
#3061                          20186 
#Immune              Inhibitory_neuron 
#1810                           7252 
#Oligodendrocyte Oligodendrocyte_precursor_cell 
#311375                           1761 
#Vascular 
#406

save(cpm_cbsx,cell_type_labels,file = '../Indep_scReference/ROSMAP_MultiregionAD_Mathys2024/cbsx_required_input.RData')

# library(omnideconv) # must restart R session first!
load('../Indep_scReference/ROSMAP_MultiregionAD_Mathys2024/cbsx_required_input.RData',verbose = T)

configure_cibersortx_credentials() 
# cbsx takes CPM as input
cbsx_sig = build_model_cibersortx(as.matrix(cpm_cbsx), 
                                  cell_type_labels,
                                  container  = "singularity",
                                  container_path = '../../../deconvolution/modules/cibersortx/',
                                  verbose = T) 
cbsx_sig_export=cbind(data.frame(GeneSymbol=rownames(cbsx_sig)),cbsx_sig)
write.table(cbsx_sig_export,file = '../Indep_scReference/ROSMAP_MultiregionAD_Mathys2024/cbsx_sig.txt',sep = '\t',row.names = F,quote = F)


############# InstaPrism reference ##############
scMeta$cell_state = paste0(scMeta$cell_type,'_',scMeta$cell_type_high_resolution)

cs_counts = table(scMeta$cell_state)
to_keep_cs = names(cs_counts)[cs_counts>5]

keep_id = which(scMeta$cell_state %in% to_keep_cs)

cpm_InstaPrism = cpm[,keep_id]
scMeta_InstaPrism = scMeta[keep_id,]

all.equal(rownames(scMeta_InstaPrism),colnames(cpm_InstaPrism)) # TRUE

refPhi = InstaPrism::refPrepare(cpm_InstaPrism,scMeta_InstaPrism$cell_type,scMeta_InstaPrism$cell_state)
all.equal(length(do.call(c,refPhi@map)),ncol(refPhi@phi.cs)) # TRUE
saveRDS(refPhi,file = '../Indep_scReference/ROSMAP_MultiregionAD_Mathys2024/refPhi.RDS')

validate_reference_outputs('ROSMAP_MultiregionAD_Mathys2024')
