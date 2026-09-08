# DALE-Eval

**DALE-Eval** is a comprehensive benchmark for evaluating **cell type-specific expression (CTSE) deconvolution** methods from bulk transcriptomic data.

## Introduction

While cell type **fraction** deconvolution has been extensively developed and benchmarked, the next frontier—**cell type-specific expression (CTSE)** deconvolution—remains largely underexplored. Here, we introduce DALE-Eval, a large-scale, biologically grounded benchmark of 10 state-of-the-art CTSE methods across diverse tissues, disease contexts, evaluation settings, and downstream tasks.

This GitHub repository is intended to provide the reproducible components of the benchmark, including:

- source code and shared utility functions;
- benchmark and method configuration files;
- method runners;
- evaluation and benchmark-summary scripts; and
- scripts used to construct benchmark inputs and reproduce the analysis workflow.

The **complete DALE-Eval benchmark release**, including the corresponding source data, benchmark objects, deconvolution outputs, intermediate objects, evaluation results, and other generated materials, is available separately at:

**https://doi.org/10.5281/zenodo.16649010**

The directory layout below describes the **complete benchmark workspace**. Because of repository size limitations, benchmarking obj, large data files, intermediate objects, and generated results shown in this layout are not stored in this GitHub repository.


---


## Comprehensive DALE-Eval benchmark directory layout

**Note:** The directory structure below describes the **complete DALE-Eval benchmark workspace** and is provided as a guide to the organization of the full source data, intermediate objects, and benchmark results. This GitHub repository contains only the relevant **scripts, configuration files, and documentation** shown within this layout; large data files and generated benchmark outputs are provided separately as part of the complete benchmark release.

### Main components

| Component | Location |
|---|---|
| Dataset-specific inputs, truth, CTSE estimates, and metrics | `Benchmarking_obj/<dataset>/` |
| Independent-reference inputs and cell-type mappings | `Indep_scReference/` |
| Benchmark configurations, shared helpers, runners, and evaluators | `DALE_Eval/` |
| Reference construction, input preparation, and workflow commands | `scripts/` |
| Auxiliary gene sets and transcript-to-gene mapping | `other_source_data/` |
| Benchmark-wide summaries and reusable summary objects | `benchmark_summary/` |
| Figure code and supplied figure products | `visualization/` |

### Per-dataset benchmark objects

```text
Benchmarking_obj/
└── <dataset>/
    ├── metadata/
    │   ├── sample_meta.txt
    │   └── ctse_meta.txt
    ├── bulk_input/
    │   ├── wcpm.txt
    │   └── sumcount.txt
    ├── ctse_truth/
    │   ├── meancpm/
    │   │   └── <cell_type>.txt.gz
    │   ├── sumcount/
    │   │   └── <cell_type>.txt.gz
    │   └── sumcount_cpm/
    │       └── <cell_type>.txt.gz
    ├── frac_input/
    │   ├── truth_cellfrac.txt
    │   └── truth_transcriptfrac.txt
    ├── self_reference/
    │   ├── sample_split.txt
    │   ├── limma_top_genes.csv
    │   ├── refPhi.RDS
    │   ├── cbsx_sig.txt
    │   ├── bMIND_profile.csv
    │   └── rowMeans_sig.csv
    ├── deconv_res/
    │   └── <config_slug>/
    │       └── <method>/
    │           ├── <cell_type>.txt.gz
    │           ├── InstaPrismfrac.txt       # InstaPrism only
    │           ├── BLUEfrac.txt             # BLUE only
    │           ├── scTAPEfrac.txt           # scTAPE only
    │           └── bMIND_posterior.RDS      # bMIND only, when posterior export is enabled
    ├── regressed_bulk/
    │   └── <baseline_id>.txt
    ├── evaluation_metadata/
    │   ├── cor_na_acceptable_mask/
    │   │   └── cor_na_acceptable_mask_truth-<truth_type>.txt
    │   └── sample_mask/
    │       └── InstaPrismfrac_<config_id>_gt0.1.txt
    ├── deconv_performance/
    │   ├── <config_slug>__truth-<truth_type>/
    │   │   ├── spearman_cor/
    │   │   │   └── <method>.txt
    │   │   ├── spearman_cor_highFrac_sample/
    │   │   │   └── <method>.txt
    │   │   ├── sample_cor_all_genes/
    │   │   │   └── <method>.txt
    │   │   ├── sample_cor_hallmark_genes/
    │   │   │   └── <method>.txt
    │   │   ├── sample_logmean/
    │   │   │   └── <method>.txt
    │   │   ├── celltype_pairwise_pearson_cor/
    │   │   │   └── <method>.txt
    │   │   ├── sample_mean/
    │   │   │   └── <method>[_libnorm].txt
    │   │   ├── normalized_variance/
    │   │   │   └── <method>[_libnorm].txt
    │   │   └── sample_cor_metadata.txt
    │   ├── truth_<truth_type>/
    │   │   ├── celltype_pairwise_pearson_cor/
    │   │   │   └── Z_truth.txt
    │   │   ├── sample_mean/
    │   │   │   └── Z_truth.txt
    │   │   └── normalized_variance/
    │   │       └── Z_truth.txt
    │   └── bulk_baseline/
    │       ├── sample_cor_all_genes/
    │       │   ├── truth-meancpm__bulk-wcpm.txt
    │       │   └── truth-meancpm__bulk-wcpm__regress-InstaPrismfrac_config01.txt
    │       ├── sample_cor_hallmark_genes/
    │       │   ├── truth-meancpm__bulk-wcpm.txt
    │       │   └── truth-meancpm__bulk-wcpm__regress-InstaPrismfrac_config01.txt
    │       └── spearman_cor/
    │           └── truth-<truth_type>__<baseline_id>.txt
    ├── gene_prioritization/
    │   └── <config_slug>/
    │       ├── post_hoc_logmean_contrast/
    │       │   └── <method>.txt.gz
    │       ├── post_hoc_logmean_margin/
    │       │   └── <method>.txt.gz
    │       ├── post_hoc_meanlog_margin/
    │       │   └── <method>[_libnorm].txt.gz
    │       ├── post_hoc_meanlog_mean_contrast/
    │       │   └── <method>[_libnorm].txt.gz
    │       ├── post_hoc_partial_r2/
    │       │   └── <method>[_libnorm].txt.gz
    │       └── genesigtest_adapted_fdr/
    │           └── GeneSigTest_adapted.txt.gz
    ├── scITD_res/                          # available for PBMC_Perez2022 and ROSMAP_AD430_Mathys2023
    │   ├── truth_sumcount/
    │   │   ├── sample_scores.txt
    │   │   └── gene_celltype_loading_Factor<K>.txt
    │   └── <config_slug>/
    │       └── <method>[_libnorm]/
    │           ├── sample_scores.txt
    │           ├── gene_celltype_loading_Factor<K>.txt
    │           └── negative_gene_shifts.txt
    └── DE_res/                             # available for PBMC_Perez2022 and ROSMAP_AD430_Mathys2023
        ├── truth_sumcount/
        │   └── <cell_type>.txt.gz
        ├── <config_slug>/
        │   └── <method>/
        │       └── <cell_type>.txt.gz
        ├── baseline/
        │   ├── bulk_limma/
        │   │   └── bulk.txt.gz
        │   └── bulk_truthfrac_adjusted_limma/
        │       └── bulk.txt.gz
        ├── direct_ctsDEG_config02_bulk-sumcount/
        │   ├── bMIND/
        │   │   └── <cell_type>.txt.gz
        │   ├── ENIGMAL2/
        │   │   └── <cell_type>.txt.gz
        │   ├── ENIGMAtrace/
        │   │   └── <cell_type>.txt.gz
        │   └── runtime.txt
        └── run_summaries/
            └── <run_id>.txt
```

Each dataset directory collects the benchmark inputs, ground-truth CTSE, fraction inputs, reference objects, deconvolution outputs, evaluation products, and downstream-analysis results for that dataset. Each `<cell_type>.txt.gz` CTSE file is a compressed, tab-delimited gene-by-sample matrix for one cell type. Fraction files are sample-by-cell-type matrices. Performance files store method-level metrics for a specific dataset, configuration, and truth definition.

#### Evaluation results in `deconv_performance/`



| Folder | Main content and purpose |
|---|---|
| `spearman_cor/` | Gene-by-cell-type Spearman correlations between inferred and truth CTSE across eligible samples; the main CTSE-accuracy measure. |
| `spearman_cor_highFrac_sample/` | The same gene-level comparison after additionally retaining samples with the relevant InstaPrism fraction greater than 0.1. |
| `sample_cor_all_genes/` | Sample-by-cell-type Spearman correlations between inferred and truth expression across all overlapping genes. |
| `sample_cor_hallmark_genes/` | The same within-sample comparison restricted to overlapping Hallmark genes. |
| `celltype_pairwise_pearson_cor/` | For each gene, Pearson correlations across samples between pairs of cell types; measures cross-cell-type expression dependence within a method or within truth. |
| `sample_mean/` | Gene-by-cell-type averages across samples after method-specific preprocessing; used to compare average expression patterns. |
| `sample_logmean/` | Log-mean expression summaries for gene prioritization, with method-specific handling of already-transformed estimates. |
| `normalized_variance/` | Gene variability relative to the fitted mean–variance trend, calculated separately for each cell type using scITD. |
| `truth_<truth_type>/` | Truth-only counterparts of cell-type-pairwise correlation, sample mean, and normalized variance, stored as `Z_truth.txt` in the corresponding metric subfolders. |
| `bulk_baseline/` | Raw and fraction-regressed bulk comparisons against CTSE truth; `spearman_cor/` contains gene-level correlations, while the two `sample_cor_*` folders contain within-sample correlations. |

The last two folders are directly under `deconv_performance/`. `sample_cor_metadata.txt` records method-level gene and sample counts from the first matched cell type. Calculation details are in benchmark_summary/README.md.


### Independent-reference inputs

```text
Indep_scReference/
├── cell_type_mapping.txt
├── refPhi_map.xlsx
├── refPhi_summary.R
└── <indep_ref>/
    ├── celltype_mapping.yaml
    ├── limma_top_genes.csv
    ├── refPhi.RDS
    ├── cbsx_sig.txt
    ├── bMIND_profile.csv
    └── rowMeans_sig.csv
```

Independent-reference directories contain the processed reference inputs required by the evaluated methods, together with cell-type mappings and marker-gene summaries where applicable. The `<indep_ref>` options are `BRCA_Wu2021`, `CRC_Lee2020`, `LUAD_Laughney2020`, `PBMC_AIDA2024`, `PBMC_refined_AIDA2024`, and `ROSMAP_MultiregionAD_Mathys2024`. 

### Evaluation framework and runners

```text
DALE_Eval/
├── configs/
│   ├── benchmark_dataset_info.txt
│   ├── benchmark_ref_assignment.txt
│   ├── deconv_configs.txt
│   ├── eval_configs.txt
│   ├── method_configs.txt
│   ├── method_default_extra_configs.txt
│   ├── scITD_configs.txt
│   ├── scITD_method_preprocessing.txt
│   ├── ctse_de_configs.txt
│   ├── BLUE_hyperparameters.yaml
│   └── scTAPE_hyperparameters.yaml
├── modules/
│   ├── config_helpers.R
│   ├── runner_helpers.R
│   ├── runner_helpers.py
│   ├── marker_helpers.R
│   ├── mapping_helpers.R
│   ├── performance_summary_helpers.R
│   ├── reference_prep_helpers.R
│   ├── scITD_helpers.R
│   ├── ctse_de_helpers.R
│   ├── direct_ctse_de_helpers.R
│   ├── pseudobulk.R
│   ├── pseudobulk.py
│   ├── evalu.R
│   ├── celltype_pairwise_cor_helpers.R
│   ├── sample_mean_nv_helpers.R
│   ├── gene_prioritization_helpers.R
│   ├── utils.py
│   └── visual_helpers.py
├── eval/
│   ├── evaluate_ctse_performance.R
│   ├── evaluate_ctse_sample_cor.R
│   ├── evaluate_method_celltype_pairwise_cor.R
│   ├── evaluate_truth_celltype_pairwise_cor.R
│   ├── evaluate_method_sample_mean_nv.R
│   ├── evaluate_truth_sample_mean_nv.R
│   └── build_gene_prioritization_scores.R
└── run/
    ├── run_InstaPrism.R
    ├── run_BLUE.py
    ├── run_scTAPE.py
    ├── run_CIBERSORTx.R
    ├── run_ENIGMA.R
    ├── run_bMIND.R
    ├── run_EPICunmix_with_posterior.R
    ├── run_TCA.R
    └── run_Unico.R
```

`DALE_Eval/` contains the core benchmarking framework. `configs/` defines dataset assignments, deconvolution settings, evaluation settings, downstream-analysis settings, and method-specific parameters. `modules/` contains shared implementation code, `run/` provides method-specific deconvolution entry points, and `eval/` provides evaluation entry points.


#### Note

Compared with archived DALE-Eval, the EPICunmix runner now processes genes in configurable chunks (1,000 genes by default) to limit memory use during fitting. This enables successful deconvolution of large datasets such as PBMC_1k1k, for which previous runs failed.

#### Deconvolution configuration options

Deconvolution outputs are organized by configuration. The base result-directory slug is:

```text
<config_id>_bulk-<bulk_input>__frac-<frac_input>__ref-<refType>
```

For non-default normalization, `__norm-<bulk_normalization>` is appended.

| Config ID | Bulk input | Scale | Normalization | Fraction input | Reference |
|---|---|---|---|---|---|
| `config01` | `wcpm` | CPM | CPM | `InstaPrismfrac` | independent |
| `config02` | `sumcount` | counts | CPM | `InstaPrismfrac` | independent |
| `config03` (exploratory) | `wcpm_centered_noise` | CPM | CPM | `InstaPrismfrac` | independent |
| `config04` (exploratory) | `wcpm_mushift_noise` | CPM | CPM | `InstaPrismfrac` | independent |
| `config05` (exploratory) | `sumcount_centered_noise` | counts | CPM | `InstaPrismfrac` | independent |
| `config06` (exploratory) | `sumcount_mushift_noise` | counts | CPM | `InstaPrismfrac` | independent |
| `config07` | `wcpm` | CPM | CPM | `InstaPrismfrac` | self |
| `config08` | `sumcount` | counts | CPM | `InstaPrismfrac` | self |
| `config09` | `wcpm` | CPM | CPM | `truth_cellfrac` | self |
| `config02_tmm`  (exploratory) | `sumcount` | counts | TMM | `InstaPrismfrac` | independent |
| `config02_uq`  (exploratory) | `sumcount` | counts | upper quartile | `InstaPrismfrac` | independent |
| `real01`  (exploratory)| `dissociated_polyA` | counts | CPM | `InstaPrismfrac` | self |

`config03`–`config06` are exploratory noise settings outside the paper analyses. Their noise inputs, QC files, per-dataset deconvolution and evaluation outputs, and baseline matrices are excluded from this release; the configuration definitions remain available for optional reruns.

`real01` is an additional real-bulk recipe without a standard CTSE evaluation entry or supplied deconvolution results.

### Project workflow scripts

```text
scripts/
├── downstream_evaluation_commands.sh
├── step1_indep_reference_construction.R
├── step1_indep_reference_construction.ipynb
├── step2_benchmarking_input.ipynb
├── step3_deconvolution_commands.sh
├── step4_run_add_bulk_noise.R # optional
├── step5_scITD_benchmark.R
└── step6_ctse_DE_benchmark.R
```

These scripts document the project-level workflow from reference construction and benchmark-input preparation through deconvolution and downstream evaluation. `step3_deconvolution_commands.sh` records deconvolution commands and optional exploratory noise runs; their noise inputs can be regenerated with `step4_run_add_bulk_noise.R`. `downstream_evaluation_commands.sh` records evaluation preparation, benchmark evaluation, scITD analysis, and gene-prioritization commands.

### Auxiliary source data

```text
other_source_data/
├── h.all.v7.5.1.symbols.gmt
└── t2g.RData
```

This directory contains auxiliary resources used by downstream project workflows. `h.all.v7.5.1.symbols.gmt` contains hallmark gene-set definitions, and `t2g.RData` contains the transcript-to-gene mapping.

### Benchmark-wide summaries

```text
benchmark_summary/
├── README.md
├── preparation/
│   ├── build_benchmark_dataset_info.R
│   ├── build_bulk_baseline_cor.R
│   ├── build_bulk_baseline_sample_cor.R
│   ├── build_cor_na_acceptable_mask.R
│   └── build_highFrac_sample_masks.R
├── ctse_spearman_cor_summary/
│   ├── summarize_spearman_by_top_genes.R
│   ├── summarize_interpretable_gene_coverage.R
│   ├── spearman_by_top_genes_summary_list.RDS
│   └── interpretable_gene_coverage_summary_list.RDS
├── ref_type_comparison_summary/
│   ├── ref_type_comparison_self_top_genes.R
│   ├── config01_vs_config07_self_top_genes_summary_list.RDS
│   └── config02_vs_config08_self_top_genes_summary_list.RDS
├── sample_filter_effect_summary/
│   ├── summarize_all_vs_highFrac_samples.R
│   ├── summarize_highFrac_sample_metadata.R
│   ├── all_vs_highFrac_samples_summary_list.RDS
│   └── highFrac_sample_metadata_list.RDS
├── covariance_specificity_summary/
│   ├── summarize_gene_celltype_covariance_specificity.R
│   ├── config01__gene_celltype_specificity.RDS
│   └── config01__gene_celltype_specificity_method_summary.RDS
├── gene_prioritization_summary/
│   ├── summarize_posthoc_vs_indep_ref.R
│   └── posthoc_vs_indep_ref_summary_list.RDS
├── sample_mean_comparision_summary/
│   ├── sample_mean_comparision_summary.R
│   ├── ccc_summary_list.RDS
│   ├── truth_weighted_ordering_score_summary_list.RDS
│   └── intra_cell_type_cor_summary.txt
├── scITD_factor_matching_summary/
│   ├── summarize_scITD_factor_matching.R
│   └── scITD_factor_matching_summary.RDS
└── DE_summary/
    ├── DEG_call_agreement_summary.R
    ├── DEG_call_metric_spec.tsv
    ├── DEG_call_agreement_summary.tsv
    ├── DEG_call_agreement_summary_both.tsv
    ├── DEG_call_agreement_summary_up.tsv
    └── DEG_call_agreement_summary_down.tsv
```

This directory contains benchmark-wide aggregation code and the summary objects used by downstream analyses and figure-generation scripts.

### Figure code and retained figure products

```text
visualization/
├── FigureS1_runtime_summary/
│   ├── FigureS1_runtime_summary.R
│   └── FigureS1_runtime_summary_config01.tsv
├── FigureS2_gene_specificity_R2/
│   ├── FigureS2_gene_specificity_R2.R
│   └── config01_gene_specificity_R2.RDS
├── FigureS3_top_n_gene_summary/
│   └── FigureS3_top_n_gene_summary.R
├── FigureS4_frac_accuracy_scatter/
│   └── FigureS4_frac_accuracy_scatter.R
├── FigureS5_celltype_pairwise_cor_boxplot/
│   └── FigureS5_celltype_pairwise_cor_boxplot.R
├── FigureS6_expression_specificity_by_top_n/
│   └── FigureS6_expression_specificity_by_top_n.R
├── FigureS7_csDEG_truth_overlap/
│   └── FigureS7_csDEG_truth_overlap.R
├── FigureS8_sample_filtering_benefit_lineplot/
│   ├── FigureS8_sample_filtering_benefit_all_datasets.R
│   └── FigureS8_sample_filtering_benefit_selected_cell_types.R
├── fig2a_UpSetR_like/
│   ├── fig2a.R
│   └── fig2a_updated.R
├── fig2b_cor_by_ct_abundance/
│   └── fig2b.R
├── fig2c_cor_by_gene_specificity/
│   └── fig2c.R
├── fig2d_cor_by_group_boxplot/
│   ├── fig2d.R
│   └── fig2d_per_group_comparison.R
├── fig3b_data_generation_comparison_scatter/
│   └── fig3b.R
├── fig3c_data_generation_comparison/
│   ├── fig3c_boxplot.R
│   ├── fig3c_cell_type_boxplot.R
│   └── fig3c_lineplot.R
├── fig3d_data_generation_method_ranking_heatmap/
│   ├── fig3d.R
│   └── fig3d_tissue_abundance_rank_heatmap.R
├── fig4a_ref_type_comparison_scatter/
│   ├── fig4a_by_cell_type.R
│   └── fig4a_by_gene_group.R
├── fig4b_frac_input_comparison_scatter/
│   └── fig4b_scatter.R
├── fig4c_frac_input_comparison_boxplot/
│   └── fig4c_boxplot.R
├── fig5b_overall_celltype_pairwise_cor_boxplot/
│   └── fig5b_overall_celltype_pairwise_cor_boxplot.R
├── fig5c_celltype_pairwise_cor_boxplot/
│   ├── fig5c_celltype_pairwise_cor_boxplot.R
│   └── fig5c_other_datasets.R
├── fig5d_pairwise_cor_heatmap/
│   └── fig5d_pairwise_cor_heatmap.R
├── fig5f_sample_mean_concordance/
│   ├── fig5f_sample_mean_concordance.R
│   └── fig5f_sample_mean_concordance_blue_overlap.R
├── fig5g_sample_mean_heatmap/
│   └── fig5g_sample_mean_heatmap.R
├── fig6b_DE_funky_heatmap/
│   ├── fig6b_DE_funky_heatmap.R
│   ├── fig6b_pre_funky_table.R
│   ├── funky_data.RDS
│   └── pre_funky_table_PBMC_Perez2022.tsv
├── fig7b_scITD_gene_loading_truth/
│   ├── fig7b_scITD_gene_loading_truth.R
│   └── fig7b_scITD_selected_gsea.tsv
├── fig7c_scITD_loading_boxplot/
│   ├── fig7c_scITD_SLE_score_boxplot.R
│   └── fig7c_scITD_truth_SLE_score_boxplot.R
├── fig7e_scITD_loading_scatter/
│   └── fig7e_scITD_loading_score_scatter.R
├── fig7f_scITD_score_summary/
│   └── fig7f_scITD_score_summary.R
├── fig8c_highFrac_filtering_benefit/
│   └── fig8c.R
├── fig8d_highFrac_method_performance/
│   └── fig8d.R
├── fig8e_posthoc_gene_prioritization/
│   ├── fig8e.R
│   └── fig8e_gene_ranking_lineplot.R
└── fig9_overall_method_performance_funky_heatmap/
    ├── fig9_overall_method_performance_funky_heatmap.R
    ├── fig9_metric_dataset_coverage.tsv
    ├── fig9_overall_method_performance_summary.tsv
    └── tmp_fig9_plot_from_summary.R
```

## Citation

If you use DALE-Eval, its framework, or its benchmark data, please cite:

> **DALE-Eval: A comprehensive cell type-specific expression deconvolution benchmark for transcriptomics data**  
> Mengying Hu, Maria Chikina, Martin Jinye Zhang  
> bioRxiv 2025.07.31.667984; https://doi.org/10.1101/2025.07.31.667984

