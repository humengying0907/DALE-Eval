#!/usr/bin/env bash
# CTSE benchmark command catalog: downstream evaluation and derived analyses
#
# Run selected commands from scripts/. This file is a reproducibility record,
# not a single end-to-end job: commands can refresh existing benchmark outputs.
# Complete the applicable step 3 deconvolution commands before these commands.

# =============================================================================
# EVALUATION PREPARATION — ALL DATASETS
# Build truth-evaluability masks before the dataset annotation, then prepare
# the high-fraction sample masks and bulk-expression baseline correlations.
# =============================================================================
Rscript ../benchmark_summary/preparation/build_cor_na_acceptable_mask.R
Rscript ../benchmark_summary/preparation/build_benchmark_dataset_info.R
Rscript ../benchmark_summary/preparation/build_highFrac_sample_masks.R
Rscript ../benchmark_summary/preparation/build_bulk_baseline_cor.R

# =============================================================================
# STANDARD CTSE SPEARMAN EVALUATION — ALL RETAINED DATASET/CONFIG RESULTS
# Evaluate every dataset/config combination currently represented in deconv_res/.
# =============================================================================
# BRCA_Bassez2021 — standard CTSE Spearman evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config02 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config02_tmm --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config02_uq --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config03 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config04 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config05 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config06 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config07 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config08 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config09 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
# CRC_Pelka2021 — standard CTSE Spearman evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset CRC_Pelka2021 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset CRC_Pelka2021 --config_id config02 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset CRC_Pelka2021 --config_id config07 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset CRC_Pelka2021 --config_id config08 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset CRC_Pelka2021 --config_id config09 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
# LUAD_Kim2020 — standard CTSE Spearman evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset LUAD_Kim2020 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset LUAD_Kim2020 --config_id config02 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset LUAD_Kim2020 --config_id config07 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset LUAD_Kim2020 --config_id config08 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset LUAD_Kim2020 --config_id config09 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
# PBMC_1k1k — standard CTSE Spearman evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_1k1k --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_1k1k --config_id config07 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
# PBMC_Perez2022 — standard CTSE Spearman evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_Perez2022 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_Perez2022 --config_id config02 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_Perez2022 --config_id config03 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_Perez2022 --config_id config04 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_Perez2022 --config_id config05 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_Perez2022 --config_id config06 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_Perez2022 --config_id config07 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_Perez2022 --config_id config08 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_Perez2022 --config_id config09 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
# PBMC_refined_1k1k — standard CTSE Spearman evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_refined_1k1k --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_refined_1k1k --config_id config07 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
# PBMC_refined_Perez2022 — standard CTSE Spearman evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_refined_Perez2022 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_refined_Perez2022 --config_id config07 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_refined_Perez2022 --config_id config09 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
# ROSMAP_AD430_Mathys2023 — standard CTSE Spearman evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD430_Mathys2023 --config_id config07 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD430_Mathys2023 --config_id config08 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
# ROSMAP_AD92_Xiong2023 — standard CTSE Spearman evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD92_Xiong2023 --config_id config03 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD92_Xiong2023 --config_id config04 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD92_Xiong2023 --config_id config05 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD92_Xiong2023 --config_id config06 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD92_Xiong2023 --config_id config07 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD92_Xiong2023 --config_id config08 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --min_n_sample 10 --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD92_Xiong2023 --config_id config09 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --min_n_sample 10 --digits 6

# =============================================================================
# STANDARD TOP-GENE SPEARMAN SUMMARY — ALL DATASETS
# Run after the preparation and standard method-evaluation commands above.
# Output: benchmark_summary/ctse_spearman_cor_summary/spearman_by_top_genes_summary_list.RDS
# =============================================================================
Rscript ../benchmark_summary/ctse_spearman_cor_summary/summarize_spearman_by_top_genes.R

# =============================================================================
# HIGH-FRACTION-SAMPLE CTSE SPEARMAN EVALUATION
# Evaluate the retained config01/config02 result sets using InstaPrism >0.1 masks.
# =============================================================================
# BRCA_Bassez2021 — high-fraction sample evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --sample_mask Benchmarking_obj/BRCA_Bassez2021/evaluation_metadata/sample_mask/InstaPrismfrac_config01_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset BRCA_Bassez2021 --config_id config02 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --sample_mask Benchmarking_obj/BRCA_Bassez2021/evaluation_metadata/sample_mask/InstaPrismfrac_config02_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
# CRC_Pelka2021 — high-fraction sample evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset CRC_Pelka2021 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --sample_mask Benchmarking_obj/CRC_Pelka2021/evaluation_metadata/sample_mask/InstaPrismfrac_config01_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset CRC_Pelka2021 --config_id config02 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --sample_mask Benchmarking_obj/CRC_Pelka2021/evaluation_metadata/sample_mask/InstaPrismfrac_config02_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
# LUAD_Kim2020 — high-fraction sample evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset LUAD_Kim2020 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --sample_mask Benchmarking_obj/LUAD_Kim2020/evaluation_metadata/sample_mask/InstaPrismfrac_config01_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset LUAD_Kim2020 --config_id config02 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --sample_mask Benchmarking_obj/LUAD_Kim2020/evaluation_metadata/sample_mask/InstaPrismfrac_config02_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
# PBMC_1k1k — high-fraction sample evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_1k1k --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --sample_mask Benchmarking_obj/PBMC_1k1k/evaluation_metadata/sample_mask/InstaPrismfrac_config01_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
# PBMC_Perez2022 — high-fraction sample evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_Perez2022 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --sample_mask Benchmarking_obj/PBMC_Perez2022/evaluation_metadata/sample_mask/InstaPrismfrac_config01_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_Perez2022 --config_id config02 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --sample_mask Benchmarking_obj/PBMC_Perez2022/evaluation_metadata/sample_mask/InstaPrismfrac_config02_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
# PBMC_refined_1k1k — high-fraction sample evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_refined_1k1k --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --sample_mask Benchmarking_obj/PBMC_refined_1k1k/evaluation_metadata/sample_mask/InstaPrismfrac_config01_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
# PBMC_refined_Perez2022 — high-fraction sample evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset PBMC_refined_Perez2022 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --sample_mask Benchmarking_obj/PBMC_refined_Perez2022/evaluation_metadata/sample_mask/InstaPrismfrac_config01_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
# ROSMAP_AD430_Mathys2023 — high-fraction sample evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --sample_mask Benchmarking_obj/ROSMAP_AD430_Mathys2023/evaluation_metadata/sample_mask/InstaPrismfrac_config01_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --sample_mask Benchmarking_obj/ROSMAP_AD430_Mathys2023/evaluation_metadata/sample_mask/InstaPrismfrac_config02_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
# ROSMAP_AD92_Xiong2023 — high-fraction sample evaluation
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --truth_type meancpm --methods all --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001 --sample_mask Benchmarking_obj/ROSMAP_AD92_Xiong2023/evaluation_metadata/sample_mask/InstaPrismfrac_config01_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --truth_type sumcount_cpm --methods all --metrics spearman_cor --filter_frac truth_transcriptfrac --min_frac 0.001 --sample_mask Benchmarking_obj/ROSMAP_AD92_Xiong2023/evaluation_metadata/sample_mask/InstaPrismfrac_config02_gt0.1.txt --min_n_sample 10 --output_tag highFrac_sample --digits 6



# =============================================================================
# CELL-TYPE PAIRWISE CORRELATION — CONFIG01 AND MEAN-CPM TRUTH
# Compare cell-type profiles within each method and within the truth matrices.
# =============================================================================
# BRCA_Bassez2021 — method and truth cell-type pairwise correlations
Rscript ../DALE_Eval/eval/evaluate_method_celltype_pairwise_cor.R --dataset BRCA_Bassez2021 --config_id config01 --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_celltype_pairwise_cor.R --dataset BRCA_Bassez2021 --truth_type meancpm --digits 6
# CRC_Pelka2021 — method and truth cell-type pairwise correlations
Rscript ../DALE_Eval/eval/evaluate_method_celltype_pairwise_cor.R --dataset CRC_Pelka2021 --config_id config01 --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_celltype_pairwise_cor.R --dataset CRC_Pelka2021 --truth_type meancpm --digits 6
# LUAD_Kim2020 — method and truth cell-type pairwise correlations
Rscript ../DALE_Eval/eval/evaluate_method_celltype_pairwise_cor.R --dataset LUAD_Kim2020 --config_id config01 --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_celltype_pairwise_cor.R --dataset LUAD_Kim2020 --truth_type meancpm --digits 6
# PBMC_1k1k — method and truth cell-type pairwise correlations
Rscript ../DALE_Eval/eval/evaluate_method_celltype_pairwise_cor.R --dataset PBMC_1k1k --config_id config01 --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_celltype_pairwise_cor.R --dataset PBMC_1k1k --truth_type meancpm --digits 6
# PBMC_Perez2022 — method and truth cell-type pairwise correlations
Rscript ../DALE_Eval/eval/evaluate_method_celltype_pairwise_cor.R --dataset PBMC_Perez2022 --config_id config01 --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_celltype_pairwise_cor.R --dataset PBMC_Perez2022 --truth_type meancpm --digits 6
# ROSMAP_AD430_Mathys2023 — method and truth cell-type pairwise correlations
Rscript ../DALE_Eval/eval/evaluate_method_celltype_pairwise_cor.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_celltype_pairwise_cor.R --dataset ROSMAP_AD430_Mathys2023 --truth_type meancpm --digits 6
# ROSMAP_AD92_Xiong2023 — method and truth cell-type pairwise correlations
Rscript ../DALE_Eval/eval/evaluate_method_celltype_pairwise_cor.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_celltype_pairwise_cor.R --dataset ROSMAP_AD92_Xiong2023 --truth_type meancpm --digits 6



# =============================================================================
# scITD BENCHMARK — SELECTED PBMC AND ROSMAP RESULTS
# Build scITD outputs for five methods across config01 and config02.
# =============================================================================
# PBMC_Perez2022 — scITD benchmark for config01/config02
Rscript step5_scITD_benchmark.R PBMC_Perez2022 config01 InstaPrism false
Rscript step5_scITD_benchmark.R PBMC_Perez2022 config01 CIBERSORTx false
Rscript step5_scITD_benchmark.R PBMC_Perez2022 config01 ENIGMAL2 false
Rscript step5_scITD_benchmark.R PBMC_Perez2022 config01 ENIGMAtrace false
Rscript step5_scITD_benchmark.R PBMC_Perez2022 config01 Unico false
Rscript step5_scITD_benchmark.R PBMC_Perez2022 config02 InstaPrism false
Rscript step5_scITD_benchmark.R PBMC_Perez2022 config02 CIBERSORTx false
Rscript step5_scITD_benchmark.R PBMC_Perez2022 config02 ENIGMAL2 false
Rscript step5_scITD_benchmark.R PBMC_Perez2022 config02 ENIGMAtrace false
Rscript step5_scITD_benchmark.R PBMC_Perez2022 config02 Unico false
# ROSMAP_AD430_Mathys2023 — scITD benchmark for config01/config02
Rscript step5_scITD_benchmark.R ROSMAP_AD430_Mathys2023 config01 InstaPrism false
Rscript step5_scITD_benchmark.R ROSMAP_AD430_Mathys2023 config01 CIBERSORTx false
Rscript step5_scITD_benchmark.R ROSMAP_AD430_Mathys2023 config01 ENIGMAL2 false
Rscript step5_scITD_benchmark.R ROSMAP_AD430_Mathys2023 config01 ENIGMAtrace false
Rscript step5_scITD_benchmark.R ROSMAP_AD430_Mathys2023 config01 Unico false
Rscript step5_scITD_benchmark.R ROSMAP_AD430_Mathys2023 config02 InstaPrism false
Rscript step5_scITD_benchmark.R ROSMAP_AD430_Mathys2023 config02 CIBERSORTx false
Rscript step5_scITD_benchmark.R ROSMAP_AD430_Mathys2023 config02 ENIGMAL2 false
Rscript step5_scITD_benchmark.R ROSMAP_AD430_Mathys2023 config02 ENIGMAtrace false
Rscript step5_scITD_benchmark.R ROSMAP_AD430_Mathys2023 config02 Unico false




# =============================================================================
# SAMPLE-MEAN AND NOISE-VARIANCE EVALUATION
# Build method summaries alongside the corresponding truth summaries.
# =============================================================================
# BRCA_Bassez2021 — method and truth sample-mean/noise-variance summaries
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset BRCA_Bassez2021 --config_id config01 --methods all --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset BRCA_Bassez2021 --config_id config02 --methods all --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset BRCA_Bassez2021 --truth_type sumcount --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset BRCA_Bassez2021 --truth_type meancpm --digits 6
# CRC_Pelka2021 — method and truth sample-mean/noise-variance summaries
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset CRC_Pelka2021 --config_id config01 --methods all --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset CRC_Pelka2021 --config_id config02 --methods all --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset CRC_Pelka2021 --truth_type sumcount --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset CRC_Pelka2021 --truth_type meancpm --digits 6
# LUAD_Kim2020 — method and truth sample-mean/noise-variance summaries
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset LUAD_Kim2020 --config_id config01 --methods all --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset LUAD_Kim2020 --config_id config02 --methods all --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset LUAD_Kim2020 --truth_type sumcount --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset LUAD_Kim2020 --truth_type meancpm --digits 6
# PBMC_1k1k — method and truth sample-mean/noise-variance summaries
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset PBMC_1k1k --config_id config01 --methods all --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset PBMC_1k1k --truth_type meancpm --digits 6
# PBMC_Perez2022 — method and truth sample-mean/noise-variance summaries
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset PBMC_Perez2022 --config_id config01 --methods all --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset PBMC_Perez2022 --config_id config02 --methods BLUE,CIBERSORTx,ENIGMAL2,ENIGMAtrace,EPICunmix,InstaPrism,TCA,Unico,bMIND,scTAPE --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset PBMC_Perez2022 --truth_type sumcount --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset PBMC_Perez2022 --truth_type meancpm --digits 6
# ROSMAP_AD430_Mathys2023 — method and truth sample-mean/noise-variance summaries
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --methods all --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --methods all --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset ROSMAP_AD430_Mathys2023 --truth_type sumcount --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset ROSMAP_AD430_Mathys2023 --truth_type meancpm --digits 6
# ROSMAP_AD92_Xiong2023 — method and truth sample-mean/noise-variance summaries
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --methods all --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_method_sample_mean_nv.R --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --methods all --lib_norm false --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset ROSMAP_AD92_Xiong2023 --truth_type sumcount --digits 6
Rscript ../DALE_Eval/eval/evaluate_truth_sample_mean_nv.R --dataset ROSMAP_AD92_Xiong2023 --truth_type meancpm --digits 6

# =============================================================================
# WITHIN-SAMPLE CROSS-GENE SPEARMAN CORRELATION — CONFIG01
# =============================================================================
# Bulk and InstaPrism-regressed bulk: prepare regressed matrices with the
# build_bulk_baseline_cor.R command above before running this generator.
# This writes all-gene and Hallmark baselines from current inputs, restricted
# to the declared test samples; see benchmark_summary/README.md.
Rscript ../benchmark_summary/preparation/build_bulk_baseline_sample_cor.R

# BRCA_Bassez2021 — within-sample cross-gene correlation
Rscript ../DALE_Eval/eval/evaluate_ctse_sample_cor.R --dataset BRCA_Bassez2021 --config_id config01 --truth_type meancpm --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_sample_cor.R --dataset CRC_Pelka2021 --config_id config01 --truth_type meancpm --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_sample_cor.R --dataset LUAD_Kim2020 --config_id config01 --truth_type meancpm --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_sample_cor.R --dataset PBMC_1k1k --config_id config01 --truth_type meancpm --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_sample_cor.R --dataset PBMC_Perez2022 --config_id config01 --truth_type meancpm --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_sample_cor.R --dataset PBMC_refined_1k1k --config_id config01 --truth_type meancpm --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_sample_cor.R --dataset PBMC_refined_Perez2022 --config_id config01 --truth_type meancpm --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_sample_cor.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --truth_type meancpm --methods all --digits 6
Rscript ../DALE_Eval/eval/evaluate_ctse_sample_cor.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --truth_type meancpm --methods all --digits 6

# =============================================================================
# TRUTH-INDEPENDENT GENE PRIORITIZATION
# Build scores for available config01/config02 sample-mean results.
# =============================================================================
# BRCA_Bassez2021 — gene-prioritization scores
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset BRCA_Bassez2021 --config_id config01 --scores all --methods all --lib_norm false --n_core 15 --digits 6
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset BRCA_Bassez2021 --config_id config02 --scores all --methods all --lib_norm false --n_core 15 --digits 6
# CRC_Pelka2021 — gene-prioritization scores
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset CRC_Pelka2021 --config_id config01 --scores all --methods all --lib_norm false --n_core 15 --digits 6
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset CRC_Pelka2021 --config_id config02 --scores all --methods all --lib_norm false --n_core 15 --digits 6
# LUAD_Kim2020 — gene-prioritization scores
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset LUAD_Kim2020 --config_id config01 --scores all --methods all --lib_norm false --n_core 15 --digits 6
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset LUAD_Kim2020 --config_id config02 --scores all --methods all --lib_norm false --n_core 15 --digits 6
# PBMC_1k1k — gene-prioritization scores
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset PBMC_1k1k --config_id config01 --scores all --methods all --lib_norm false --n_core 15 --digits 6
# PBMC_Perez2022 — gene-prioritization scores
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset PBMC_Perez2022 --config_id config01 --scores all --methods all --lib_norm false --n_core 15 --digits 6
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset PBMC_Perez2022 --config_id config02 --scores all --methods all --lib_norm false --n_core 15 --digits 6
# ROSMAP_AD430_Mathys2023 — gene-prioritization scores
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset ROSMAP_AD430_Mathys2023 --config_id config01 --scores all --methods all --lib_norm false --n_core 15 --digits 6
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset ROSMAP_AD430_Mathys2023 --config_id config02 --scores all --methods all --lib_norm false --n_core 15 --digits 6
# ROSMAP_AD92_Xiong2023 — gene-prioritization scores
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset ROSMAP_AD92_Xiong2023 --config_id config01 --scores all --methods all --lib_norm false --n_core 15 --digits 6
Rscript ../DALE_Eval/eval/build_gene_prioritization_scores.R --dataset ROSMAP_AD92_Xiong2023 --config_id config02 --scores all --methods all --lib_norm false --n_core 15 --digits 6
