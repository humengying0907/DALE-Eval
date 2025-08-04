# DALE-Eval

DALE-Eval (Deconvolution Assessment across muLtiple Environments) is a comprehensive benchmarking framework for **cell type–specific expression (CTSE)** deconvolution methods.

## 🚀Introduction
While cell type **fraction** deconvolution has been extensively developed and benchmarked, the next frontier—cell type–specific **expression** (CTSE) deconvolution—remains largely underexplored. **DALE-Eval** standardizes the evaluation of CTSE deconvolution accuracy, featuring:

- **Reproducible benchmark data construction:** end-to-end pipeline to build benchmark objects, including bulk input, CTSE ground truth, reference profiles, and gene-level specificity annotations.
- **Multi-method evaluation:** supports methods such as BayesPrism/InstaPrism, bMIND, EPICunmix, Unico, TCA, CIBERSORTx, and ENIGMA.  
- **Multi-aspect assessment:** quantifies performance along multiple dimensions, including gene expression prediction accuracy, cell type-specificity and more

## Directory Layout (suggested)


### 1. Generating a Benchmarking Object

Benchmark objects bundle:
- Bulk expression
- Reference profiles (e.g., from scRNA-seq)
- Ground-truth fractions
- Ground-truth CTSE profiles

Example:

```bash
cat Benchmarking_obj/BRCA_example/data_initialization.R
```

### 2. Running Deconvolution



```bash
cd Benchmarking_obj/BRCA_example
cat deconv_commands.txt
```



Expected outputs per method:
- Estimated CTSE expression profiles (gene × sample × cell type)
- Estimated fractions (if applicable)
- Log (runtime)

### 3. Performance evaluation with DALE-eval

Compare inferred profiles to truth and baselines:

```bash
cd DALE_eval
Rscript DALE_indep_ref_matching.R -o BRCA_example
Rscript DALE_Z_cor.R -o BRCA_example 
Rscript DALE_post_hoc_analysis.R -o BRCA_example
Rscript DALE_by_limma_cor_summary.R -o BRCA_example
Rscript DALE_a_priori_analysis.R -o BRCA_example
Rscript DALE_Z_specificity.R -o BRCA_example
```

Core metrics (configurable):
- Gene-level correlation 
- Covariance/Expression specificity 
- Post-hoc prioritized gene specificity
- Performance summary table

The resulting summary files are available in the deconvSummary/ folder of each benchmarking object directory.

## Citation
If you use this framework or its data, please cite:

> DALE-Eval: A comprehensive cell type-specific expression deconvolution benchmark for transcriptomics data
Mengying Hu, Maria Chikina, Martin Jinye Zhang
bioRxiv 2025.07.31.667984; doi: https://doi.org/10.1101/2025.07.31.667984



