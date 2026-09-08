# Benchmark Summary Definitions

This folder contains the main benchmark-summary workflows. The summary scripts
read completed evaluation results under Benchmarking_obj/ and save reusable
RDS summary lists beside the corresponding scripts. They do not generate CTSE
method estimates.

## Core Definitions

| Term | Definition |
|---|---|
| group | Limma-ranked marker-gene group: top_10, top_30, top_100, top_300, top_1000, top_3000, top_10000, or all_genes. Top N is selected before later overlap and eligibility filtering. |
| Acceptable-NA mask | A truth-side mask identifying genes for which a correlation NA is expected because fewer than 10 eligible test samples are available or truth CTSE is constant over the eligible test samples. These genes are excluded from the standard summary denominator. |
| n_genes | Number of selected marker genes remaining after the required gene overlap and truth-evaluability filtering. It can be smaller than the nominal top N. |
| n_constant_genes | Among n_genes, the number with Spearman_cor = NA. Because truth-side non-evaluable genes were already removed, this generally represents constant method estimates under the benchmark contract that methods export the complete gene-by-sample result. More generally, it records any remaining non-finite method-side correlation. |
| n_eval_samples | Number of test samples actually eligible for the correlation for the relevant truth type and cell type. This is not the total dataset sample size. |
| avg_cor | Mean of finite gene-level Spearman correlations among n_genes. |
| avg_cor_with_NA_penalty | Mean over the same n_genes after replacing remaining non-finite correlations with zero. This penalizes method-side constant or otherwise undefined correlations. |

## Evaluation Sample Definitions

For meancpm truth, the standard evaluation samples are:

~~~text
test samples
intersect samples present in meancpm CTSE truth
intersect truth_cellfrac > 0.001
~~~

For sumcount_cpm truth, the standard evaluation samples are:

~~~text
test samples
intersect samples present in sumcount_cpm CTSE truth
intersect truth_transcriptfrac > 0.001
~~~

The canonical counts are stored as n_eval_samples_meancpm and
n_eval_samples_sumcount_cpm in benchmark_dataset_info.txt. A summary row maps
its truth type to the corresponding count. Under the benchmark contract that
each method exports all test samples, the count is shared across methods and
top-gene groups for the same dataset, truth type, and cell type.

The acceptable-NA masks are built using these same test-sample and
truth-fraction rules. A cell type with fewer than 10 eligible samples has no
truth-evaluable genes for correlation.

The dataset annotation also stores gene-universe counts for cross-setting
comparisons. `n_truth_genes` is the common CTSE truth gene count;
`n_genes_indep_ref` is the gene count in the assigned reference's
`rowMeans_sig.csv`; and `n_genes_truth_indep_ref_overlap` is the overlap of
those complete gene universes. These dataset-level values are repeated across
cell types. `n_eval_genes_meancpm` and `n_eval_genes_sumcount_cpm` count
zeros in the corresponding acceptable-NA mask for each cell type, meaning genes
that remain evaluable for correlation. An unavailable truth type is recorded
as `NA`.

## Standard Top-Gene Spearman Summary

Script:

~~~text
ctse_spearman_cor_summary/summarize_spearman_by_top_genes.R
~~~

Output:

~~~text
ctse_spearman_cor_summary/spearman_by_top_genes_summary_list.RDS
~~~

For the standard summary, marker rankings follow the reference type of each
configuration. Independent-reference configurations use the
independent-reference limma statistics, while self-reference configurations
use the self-reference limma statistics. Consequently, two configurations
with different reference types may summarize different top-N genes.

Within each method, cell type, and group, the script selects the original top N
genes, intersects them with the correlation result, removes genes marked
truth-side non-evaluable by the acceptable-NA mask, and reports the resulting
n_genes, correlation averages, remaining n_constant_genes, and n_eval_samples.

The standard summary also reports three baseline columns:

~~~text
avg_cor_bulk
avg_cor_bulk_truthFrac_regressed
avg_cor_bulk_InstaPrismFrac_regressed
~~~

For one dataset, config, cell type, and group, the config-specific top-gene
ranking and acceptable-NA mask are applied consistently to all three baseline
matrices. These matrices have the same gene and cell-type universe under the
current baseline contract, so the three baseline averages in that row are
calculated over the same genes. The resulting baseline values are summarized
once and repeated across all method rows with the same config, cell type, and
group.

Baseline genes are selected independently from the method correlation file.
Consequently, the genes used for the three repeated baseline averages can
differ from the genes represented by the method row's n_genes, avg_cor, and
avg_cor_with_NA_penalty. The baseline columns should therefore be interpreted
as repeated config/cell-type/group annotations, not as method-gene-matched
averages.

Config01 and config07 read the same raw-wcpm and truth-cell-fraction-regressed
Spearman matrices; config02 and config08 similarly read the same raw
CPM-normalized-sumcount and truth-transcript-fraction-regressed matrices. Each
config still summarizes those shared matrices using its own reference-specific
top-gene ranking. The InstaPrism-fraction-regressed baseline is config-specific:
config01, config02, config07, and config08 read matrices regressed using their
matching InstaPrism fraction estimates.

## Interpretable Gene Coverage

Script:

~~~text
ctse_spearman_cor_summary/summarize_interpretable_gene_coverage.R
~~~

Output:

~~~text
ctse_spearman_cor_summary/interpretable_gene_coverage_summary_list.RDS
~~~

This reference-aware summary has one row per dataset, independent-reference
configuration, expected truth type, cell type, and CTSE method. It measures the
fraction of biologically evaluable genes for which a method delivers a finite
gene-level Spearman correlation.

For each dataset, truth type, and cell type, the denominator is:

~~~text
truth/reference/variable gene universe =
  truth genes
  intersect genes in the assigned independent reference
  intersect genes marked truth-evaluable by the acceptable-NA mask
~~~

The acceptable-NA mask removes genes with a truth-side reason for an undefined
correlation, including insufficient eligible samples or constant truth CTSE.
Within this denominator, the summary reports:

| Column | Definition |
|---|---|
| `n_truth_indep_ref_variable_genes` | Number of genes in the truth/reference/variable denominator. |
| `n_method_genes_in_universe` | Denominator genes present in the method correlation result, whether their correlations are finite or not. |
| `n_missing_method_genes` | Denominator genes absent from the method result. |
| `n_nonfinite_method_genes` | Present denominator genes with a non-finite method correlation. Because truth-side non-evaluable genes were removed first, these generally reflect constant or otherwise uninterpretable method estimates. |
| `n_interpretable_genes` | Present denominator genes with a finite method correlation. |
| `interpretable_gene_coverage` | `n_interpretable_genes / n_truth_indep_ref_variable_genes`. |

Missing method genes and non-finite method correlations both reduce
`interpretable_gene_coverage`. The value is therefore in `[0, 1]`, where `1`
means that the method delivers a finite correlation for every gene in the
truth/reference/variable universe.

The standard summary `all_genes` counts are not restricted to the assigned
independent-reference gene universe. They must not be divided by
`n_truth_indep_ref_variable_genes`; numerator and denominator would refer to
different gene universes and the resulting ratio could exceed one.

## Reference-Type Comparisons

Script:

~~~text
ref_type_comparison_summary/ref_type_comparison_self_top_genes.R
~~~

Outputs:

~~~text
ref_type_comparison_summary/config01_vs_config07_self_top_genes_summary_list.RDS
ref_type_comparison_summary/config02_vs_config08_self_top_genes_summary_list.RDS
~~~

These paired comparisons ask how config01/config07 and config02/config08
perform over the same genes. Top-N genes are always defined by the matching
self-reference limma ranking. For each CTSE method, the selected genes are
intersected with the two method correlation matrices and the truth-evaluable
gene set. The `bulk_InstaPrismFrac_regressed` comparison is represented as a
standalone method: its selected genes are intersected only with the two
matching config-specific bulk-regressed matrices and the truth-evaluable gene
set. It is not intersected with any CTSE method matrix.

The exact resulting gene set is used for both configurations within a row.
Therefore, differences between paired values are not caused by different
marker genes. Each output contains only the paired method statistics:

~~~text
method
cell_type
group
n_genes
n_eval_samples
n_constant_genes_<indep_config>
n_constant_genes_<self_config>
avg_cor_<indep_config>
avg_cor_<self_config>
avg_cor_with_NA_penalty_<indep_config>
avg_cor_with_NA_penalty_<self_config>
~~~

## All-Sample Versus High-Fraction Sample Comparison

Script:

~~~text
sample_filter_effect_summary/summarize_all_vs_highFrac_samples.R
~~~

Output:

~~~text
sample_filter_effect_summary/all_vs_highFrac_samples_summary_list.RDS
~~~

This comparison is performed separately for config01 and config02. Top-N genes
are always defined using the independent-reference limma statistics.

The filtered evaluation samples are recalculated as:

~~~text
config01:
test samples
intersect samples present in meancpm CTSE truth
intersect truth_cellfrac > 0.001
intersect config01 InstaPrismfrac > 0.1 mask

config02:
test samples
intersect samples present in sumcount_cpm CTSE truth
intersect truth_transcriptfrac > 0.001
intersect config02 InstaPrismfrac > 0.1 mask
~~~

Filtered-sample truth eligibility is redetermined rather than reusing the
all-sample acceptable-NA mask. A gene is filtered-sample eligible only when at
least 10 filtered samples are available and truth CTSE is nonconstant over
those samples.

For each row, the script selects the original independent-reference top N and
then intersects it with:

~~~text
filtered-sample truth-evaluable genes
intersect genes in the all-sample method correlation
intersect genes in the filtered-sample method correlation
intersect genes in both all-sample bulk baselines
intersect genes in both filtered-sample bulk baselines
~~~

This one exact gene set defines n_genes and is applied to all method and
baseline columns in the row. The output records both n_eval_samples_all and
n_eval_samples_filtered. If fewer than 10 filtered samples are available, the
row is retained with n_genes = 0 and correlation summaries equal to NA.

## High-Fraction Sample Metadata

Script:

~~~text
sample_filter_effect_summary/summarize_highFrac_sample_metadata.R
~~~

Output:

~~~text
sample_filter_effect_summary/highFrac_sample_metadata_list.RDS
~~~

This method-independent metadata has one row per dataset, independent-reference
configuration, expected truth type, and mapped cell type. It is shared across
CTSE methods and top-gene groups. High-fraction samples are test samples with a
config-specific InstaPrism fraction greater than `instaprism_min_frac`, which
defaults to 0.1.

The sample-retention and abundance columns use two distinct sample universes:

| Column | Definition |
|---|---|
| `n_test_samples_with_inferred_fraction` | Test samples present in the config-specific InstaPrism fraction result and sample mask. |
| `n_highFrac_samples` | Those test samples passing the high-fraction mask. |
| `sample_fraction_left` / `sample_percent_left` | `n_highFrac_samples / n_test_samples_with_inferred_fraction`, stored as a fraction and percentage. |
| `n_eval_samples_all` | Standard truth-evaluation samples after intersecting test samples, truth CTSE samples, and the config-specific truth-fraction threshold. |
| `n_eval_samples_filtered` | Standard truth-evaluation samples also passing the high-fraction mask. |
| `eval_sample_fraction_left` / `eval_sample_percent_left` | `n_eval_samples_filtered / n_eval_samples_all`, stored as a fraction and percentage. |
| `inferred_ct_abundance_all_test_samples` | Mean config-specific InstaPrism fraction over test samples with an inferred fraction; this is the abundance used by the adaptive top-N rule. |
| `inferred_ct_abundance_after_filtering` | Mean config-specific InstaPrism fraction over retained high-fraction samples. |
| `passes_min_n_sample` | Whether `n_eval_samples_filtered` meets the configuration's minimum evaluation-sample count. |

The output RDS is a named list with one data frame per dataset. The script also
leaves the combined `highFrac_sample_metadata` table inspectable in the R
session.

## Within-Sample CTSE Correlations

`DALE_Eval/eval/evaluate_ctse_sample_cor.R` measures how well a method
reconstructs the expression profile across genes within each sample and cell
type. It calculates Spearman correlation between the truth and inferred gene
vectors, producing a sample-by-cell-type matrix for each method.

Inputs are the per-cell-type `.txt.gz` matrices under
`Benchmarking_obj/<dataset>/ctse_truth/<truth_type>/` and
`Benchmarking_obj/<dataset>/deconv_res/<config_slug>/<method>/`, plus the
Hallmark GMT file. Cell types are matched by filename. Within each matched
cell type, the evaluator intersects truth and method genes and sample IDs.

The two gene scopes are all overlapping genes and their intersection with the
union of Hallmark genes. The default GMT is
`other_source_data/h.all.v7.5.1.symbols.gmt`; `--hallmark_gmt` selects another
file. Input `NA` values become zero. The optional `--truth_transform` and
`--estimate_transform` settings independently accept `none` (default) or
`log2p1`; the latter applies `log2(x + 1)` and rejects negative inputs. There
is no automatic method-specific preprocessing.

All truth/method-overlapping samples are evaluated. The script applies no
fraction filter, sample mask, or explicit `sample_split.txt` restriction;
testing-only evaluation relies on the samples stored in the method output.
There is no minimum-sample-count rule because each correlation is calculated
across genes. Fewer than two overlapping genes or a constant truth or inferred
gene vector gives `NA`. Output rows combine the sample IDs across cell types,
with `NA` where a sample is unavailable for a cell type.

For example, run from `benchmark_summary/` after preparing truth and method
CTSE outputs:

~~~bash
Rscript ../DALE_Eval/eval/evaluate_ctse_sample_cor.R \
  --dataset LUAD_Kim2020 --config_id config01 --truth_type meancpm \
  --methods all --digits 6
~~~

`--methods` accepts `all` or a comma-separated list. The output directory is
`Benchmarking_obj/<dataset>/deconv_performance/<config_slug>__truth-<truth_type>/`:

~~~text
sample_cor_all_genes/<method>.txt
sample_cor_hallmark_genes/<method>.txt
sample_cor_metadata.txt
~~~

Correlations are rounded to six decimal places by default. Rerunning replaces
the selected method tables and writes metadata for the methods processed in
that run. The metadata has one row per method and gene scope, recording input
gene counts, overlap counts, genes used, and common samples from the first
alphabetically matched cell type; it is not a per-cell-type count table.

The bulk-baseline sample correlations described below use the same correlation
helper, but additionally apply the declared test split and truth-fraction
threshold. A comparison with those baselines therefore needs the same sample
and cell-type eligibility applied to the method tables during analysis.

## Bulk Baselines

The standard top-gene summary uses three baseline types:

| Baseline | Definition |
|---|---|
| Bulk | Prepared bulk expression used as the same gene-expression estimate for each cell type. |
| Bulk truth-fraction regressed | Prepared bulk expression after regressing out the matching truth fraction. |
| Bulk InstaPrism-fraction regressed | Prepared bulk expression after regressing out the InstaPrism fractions matching the evaluated config. Config07 and config08 use their self-reference fraction estimates. |

The all-sample versus high-fraction workflow uses only the first two baseline
types. In both workflows, baseline correlations use samples marked `test` in
`Benchmarking_obj/<dataset>/self_reference/sample_split.txt`, intersected with
samples present in the truth and baseline estimate and passing the matching
truth-fraction threshold (`> 0.001`). The high-fraction workflow additionally
applies its cell-type-specific sample mask. Here, "all samples" means all
eligible test samples.

Regression fitting uses all samples shared by the prepared bulk matrix and its
regression-fraction matrix. Test-sample selection is applied to correlation
evaluation after fitting. The builder also writes a local baseline manifest
recording `sample_policy=test_samples` and, for regressed estimates,
`regression_sample_policy=bulk_fraction_overlap`. Bulk-baseline manifests are
omitted from the public bundle and are not required by the main Spearman
summary.

To refresh these baseline correlations and the main Spearman summary using the
supplied benchmark inputs, method evaluations, and eligibility masks, run from
`benchmark_summary/`:

~~~bash
Rscript preparation/build_bulk_baseline_cor.R
Rscript ctse_spearman_cor_summary/summarize_spearman_by_top_genes.R
~~~

The baseline builder defaults to all available datasets and bulk recipes and
also refreshes high-fraction baseline correlations where the required masks
are supplied. It rewrites the corresponding regressed expression matrices and
baseline manifests using the existing regression-fit policy. The summary
command then replaces
`ctse_spearman_cor_summary/spearman_by_top_genes_summary_list.RDS`.

To also refresh the all-sample versus high-fraction summary, which contains the
same corrected bulk baselines, run afterward:

~~~bash
Rscript sample_filter_effect_summary/summarize_all_vs_highFrac_samples.R
~~~

For high-fraction evaluation, the regression is not refitted on the filtered
samples. The existing bulk and truth-fraction-regressed matrices are reused,
and only the samples used to calculate gene correlations are filtered. This
keeps the estimator fixed and isolates the effect of sample filtering.

The all-sample versus filtered summary reports four baseline columns:

~~~text
avg_cor_bulk_all_samples
avg_cor_bulk_filtered_samples
avg_cor_bulk_truthFrac_regressed_all_samples
avg_cor_bulk_truthFrac_regressed_filtered_samples
~~~

All four are calculated over the same genes used for the paired method
comparison.

### Within-sample bulk correlations

`preparation/build_bulk_baseline_sample_cor.R` calculates Spearman correlation
across genes for each sample and cell type, using prepared bulk expression or
bulk expression regressed on the configuration's InstaPrism fractions. These
sample-by-cell-type outputs differ from the gene-by-cell-type correlations
used by the top-gene summaries above.

Run from `scripts/` in this order:

```sh
Rscript ../benchmark_summary/preparation/build_bulk_baseline_cor.R
Rscript ../benchmark_summary/preparation/build_bulk_baseline_sample_cor.R
```

The first script creates the regressed bulk matrices. The second reads the
selected configuration from `DALE_Eval/configs/deconv_configs.txt` and its
truth, fraction-filter, threshold, and baseline identifiers from
`DALE_Eval/configs/eval_configs.txt`. Its editable defaults select all nine
datasets, `config01`, and both all-gene and Hallmark correlations. For another
configuration, select datasets with the corresponding truth and bulk inputs
and prepare its regressed matrix first. Configurations without defined bulk
baselines are rejected.

Both baselines are restricted to samples marked `test` in
`Benchmarking_obj/<dataset>/self_reference/sample_split.txt` using the shared
`restrict_to_test_samples()` helper. This excludes training samples and samples
absent from the split. For each cell type, the calculation then intersects
truth, baseline, and fraction sample IDs and truth/baseline gene IDs.
Hallmark evaluation further intersects the
gene IDs with `other_source_data/h.all.v7.5.1.symbols.gmt`. Truth and expression
input NA values become zero, following the shared CTSE reader. Correlations
remain NA for undefined profiles and samples whose truth fraction is missing
or at most the configured threshold (0.001 for config01). There is no
minimum-sample-count requirement for a correlation across genes. The regressed
bulk is reused without refitting on the fraction-filtered samples.

The supplied correlation tables are under
`Benchmarking_obj/<dataset>/deconv_performance/bulk_baseline/`:

```text
sample_cor_all_genes/truth-<truth_type>__<baseline_id>.txt
sample_cor_hallmark_genes/truth-<truth_type>__<baseline_id>.txt
```

For config01, the two baseline IDs are `bulk-wcpm` and
`bulk-wcpm__regress-InstaPrismfrac_config01`. Running the generator replaces the
selected tables. Generated correlations are rounded to the editable `digits`
setting (default 6). The generator also writes `sample_cor_manifest.tsv` as
local provenance; this file is omitted from the public bundle and is not
required to read the supplied correlation tables. It records the scope,
fraction filter, threshold, source, dimensions, output filename, and sample
policy. Selected entries are
updated to `result_source=current_inputs` and `sample_policy=test_samples`;
entries for other tables remain. Its `metric_file` paths are relative to the
bulk-baseline directory and `source_file` paths are repository-relative.

The supplied all-gene and Hallmark tables were calculated from current inputs,
restricted to the declared test samples, and rounded to six decimal places.
Both baselines use the same test-sample set as the config01 method results in
all nine datasets. BRCA has 68 samples; ROSMAP AD430 has 80. Undefined
correlations and fraction filtering can still give different numbers of
finite values across cell types and methods.

Baseline fraction masks are stored as NA values in the tables. To change the
threshold, update the config01 entry in `DALE_Eval/configs/eval_configs.txt`
and regenerate the baseline tables.

## Sample Mean, Sample Log-Mean, and Normalized Variance Inputs

The gene-by-cell-type inferred sample-mean, sample-log-mean, and
normalized-variance matrices, together with truth sample-mean and
normalized-variance matrices, are generated by:

~~~text
DALE_Eval/eval/evaluate_method_sample_mean_nv.R
DALE_Eval/eval/evaluate_truth_sample_mean_nv.R
~~~

For each cell type `c`, preprocessing produces a sample-by-gene matrix `Y_c`.
For gene `g`, the stored sample mean is calculated across the samples in that
matrix:

~~~text
sample_mean[g, c] = mean_s(Y_c[s, g])
~~~

For inferred CTSE, the standard `lib_norm=false` evaluation also stores a
separate `sample_logmean` matrix for gene prioritization. It is calculated
directly from the aligned original CTSE values rather than from `sample_mean`.
Let

~~~text
raw_sample_mean[g, c] = mean_s(X_c[g, s])
~~~

For methods in the `log1p` preprocessing group, define one shift per gene over
the benchmark-matched cell-type means and then use the archived floor:

~~~text
raw_mean_shift[g] = max(0, -min_c(raw_sample_mean[g, c]))

sample_logmean[g, c]
  = log2(max(raw_sample_mean[g, c] + raw_mean_shift[g], 1e-10))
~~~

For methods marked `already_transformed`, `sample_logmean` is the arithmetic
mean of the original transformed CTSE values, without another logarithm or
the sample-level negative-value shift used for `sample_mean`. Consequently,
`sample_logmean` and `sample_mean` are identical for an already-transformed
gene/cell-type profile only when that profile did not require the latter
shift. The retained `DALE_Eval/eval/evaluate_method_sample_mean_nv.R` workflow
writes `sample_logmean` during standard `lib_norm=false` evaluations.

Normalized variance is estimated from the same `Y_c` by
`scITD::get_normalized_variance()`, independently for each cell type. scITD
calculates gene means and variances across samples, fits the cell-type-specific
mean-variance trend across genes using a GAM of log variance against log mean,
and measures each gene's departure from its expected variance. The stored
value is scITD's `norm_variances` value, not ordinary sample variance and not
the later `NV^2` tensor weight.

For inferred CTSE, the calculation uses the benchmark-matched cell types, all
genes shared across those cell types, and every aligned sample stored by the
method. All stored samples must belong to the testing split; the method is not
intersected with truth samples or truth genes. If a gene/cell-type profile has
negative values, one constant is added to all its samples so its minimum is
zero. The configured method preprocessing is then applied:

| Method preprocessing | Methods | Definition of `Y` |
|---|---|---|
| `log1p` | CIBERSORTx, ENIGMAL2, ENIGMAtrace, InstaPrism, Unico | Apply `log1p` directly after any negative-value shift. No library normalization is used in the standard `lib_norm=false` evaluation. |
| `already_transformed` | BLUE, bMIND, EPICunmix, scTAPE, TCA | Retain the shifted method scale without an additional log or library normalization. |

The optional `lib_norm=true` mode applies library normalization followed by
`log1p` only to methods in the `log1p` group and writes sample-mean and
normalized-variance files with a `_libnorm` suffix. It does not write
`sample_logmean`; the log-mean branch is currently standard-only.

Truth calculations use all truth cell types, all genes shared across them, and
every testing sample. Truth preprocessing is:

| Truth type | Definition of `Y` |
|---|---|
| `sumcount` | Within each cell type, calculate TMM factors from positive-library profiles, normalize to the configured scale factor of 10,000, and apply `log1p`. |
| `meancpm` | Apply `log1p` directly, without TMM or scale-factor rescaling. |
| `sumcount_cpm` | Apply `log1p` directly, without TMM or scale-factor rescaling. |

An all-zero truth sample/cell-type profile remains zero in `Y`. Neither method
nor truth generation applies fraction filtering, a fraction mask, all-zero
profile removal, variable-gene selection, or normalized-variance QC. All genes
in the common input universe remain in the output; non-finite normalized
variances are written as `NA`.

These matrices have genes as rows and cell types as columns and are written to
six decimal places:

~~~text
Benchmarking_obj/<dataset>/deconv_performance/<config_slug>__truth-<expected_truth_type>/sample_mean/<method>[_libnorm].txt
Benchmarking_obj/<dataset>/deconv_performance/<config_slug>__truth-<expected_truth_type>/sample_logmean/<method>.txt
Benchmarking_obj/<dataset>/deconv_performance/<config_slug>__truth-<expected_truth_type>/normalized_variance/<method>[_libnorm].txt
Benchmarking_obj/<dataset>/deconv_performance/<truth_output_folder>/sample_mean/Z_truth.txt
Benchmarking_obj/<dataset>/deconv_performance/<truth_output_folder>/normalized_variance/Z_truth.txt
~~~

For method outputs, `expected_truth_type` comes from `eval_configs.txt` for the
requested config and affects only the output path. Truth CTSE is not used in
the method calculation.

Truth output folders are `truth_sumcount`, `truth_meancpm`, and
`truth_sumcount_cpm` for the corresponding truth types.

## `scITD_factor_matching_summary`

The notebook-style script
`scITD_factor_matching_summary/summarize_scITD_factor_matching.R` matches
independently fitted inferred scITD factors to truth factors. By default it
summarizes config01 and config02 results for PBMC Perez and ROSMAP. It includes
available standard results and `_libnorm` results for methods configured with
`preprocess_mode=log1p`: CIBERSORTx, ENIGMAL2, ENIGMAtrace, InstaPrism, and
Unico.

For each dataset, config, and inferred result, the script intersects truth and
method genes and cell types and gives both loading matrices the same row and
column order. Each factor gene-by-cell-type loading matrix is flattened in R
column-major order: all genes for the first cell type, then all genes for the
second cell type, and so on. Pearson correlation between every truth and
method loading vector produces a truth-factor-by-method-factor correlation
matrix.

The retained factor match is the one-to-one assignment that maximizes the sum
of absolute loading correlations across factors. With the current five
factors, the script evaluates all 5! possible assignments. Every truth factor
is paired with one inferred factor, and every inferred factor is used once.

This global match can differ from independently selecting the largest absolute
correlation in each truth-factor row. Independent row-wise selection can reuse
the same inferred factor for multiple truth factors and leave other inferred
factors unused. It therefore does not define a coherent factor permutation and
could reuse one inferred sample-score vector in several truth-factor
comparisons. For this reason, row-wise-best matches are not retained. The stored
`sample_score_pearson` uses the one-to-one match and multiplies the matched
method score by `sign_flip`; factors are not rematched using sample scores.

The script writes one RDS:

~~~text
scITD_factor_matching_summary/scITD_factor_matching_summary.RDS
~~~

The RDS has two elements:

| Element | Contents |
|---|---|
| `loading_correlation_matrix_list` | Complete truth-factor-by-method-factor Pearson correlation matrix for each available result. Entries are named `<dataset>__<config_id>__<method>`, with `_libnorm` retained in the method name where applicable. |
| `factor_matching` | Matched-factor data frames grouped by `<dataset>__<config_id>`. Unavailable method results are absent. |

Each `factor_matching` data frame contains:

| Column | Meaning |
|---|---|
| `method` | Inferred CTSE result whose scITD factors were fitted; library-normalized variants retain the `_libnorm` suffix. |
| `truth_factor` | Truth factor being matched. |
| `matched_method_factor` | Inferred factor selected by the global one-to-one assignment. |
| `loading_pearson` | Signed Pearson correlation between the flattened truth and matched inferred loading matrices. |
| `abs_loading_pearson` | Absolute value of `loading_pearson`, used as the matched-pair strength. |
| `sign_flip` | `-1` for a negative matched correlation and `1` otherwise. Multiply inferred loadings and sample scores by this value before aligned comparisons. |
| `sample_score_pearson` | Pearson correlation across shared samples between the truth-factor score and the globally matched inferred-factor score after applying `sign_flip` to the inferred score. |

### scITD loading FDR

The publication-style fig7b/fig7d loading workflow reconstructs scITD's loading
significance test separately for each truth or inferred scITD result. For every
gene, cell type, and sample factor, it tests the association across samples
between the preprocessed cell-type-specific expression values and the sample
factor score. The p-values from all genes, cell types, and factors in that one
scITD result are adjusted together with `stats::p.adjust(method = "fdr")`, which
is the Benjamini-Hochberg procedure in R. Thus, the adjustment is global within
one scITD result but is not pooled across truth and inferred methods.

The editable cutoff in
`visualization/fig7b_scITD_gene_loading_truth/fig7b_scITD_gene_loading_truth.R`
is:

~~~text
loading_fdr_threshold = 0.02
~~~

An entry is treated as significant when its adjusted FDR is strictly below the
cutoff. At the default 0.02 cutoff, the procedure is intended to limit the
expected false-discovery proportion among declared gene/cell-type/factor
associations to approximately 2%, under the assumptions of the adjustment. An
FDR value is not the probability that one specific loading is false, and it
does not measure the loading's direction or magnitude.

For the loading heatmap, truth genes are retained when the selected truth
factor has FDR below the cutoff in at least one cell type. For truth and each
matched inferred factor, available gene/cell-type entries with FDR greater than
or equal to the cutoff are displayed as zero. Truth-panel genes absent from an
inferred method remain `NA` and are displayed in gray. The gene rows are
ordered once from the filtered truth matrix and that exact order is reused in
every inferred panel.

Increasing the cutoff retains more truth genes and more nonzero inferred
entries; decreasing it produces a stricter, sparser heatmap. Because the truth
gene panel is defined using this cutoff, changing it can alter both the plotted
gene set and its truth-derived row order.

## `sample_mean_comparision_summary`

This section defines the complete config01 sample-mean comparison summary and
the formulas used by its three outputs. The notebook-style script
`sample_mean_comparision_summary/sample_mean_comparision_summary.R` compares
config01 inferred and truth sample-mean matrices. It creates two gene-level RDS
lists and one cell-type-level tab-separated table:

~~~text
sample_mean_comparision_summary/truth_weighted_ordering_score_summary_list.RDS
sample_mean_comparision_summary/ccc_summary_list.RDS
sample_mean_comparision_summary/intra_cell_type_cor_summary.txt
~~~

| Output | Statistical unit | Comparison axis |
|---|---|---|
| Truth-weighted, tie-aware Kendall-style ordering score | One value per dataset, gene, and method | Matched cell types within a gene |
| Concordance correlation coefficient (CCC) | One value per dataset, gene, and method | Matched cell types within a gene |
| Intra-cell-type Spearman correlation | One value per dataset, method, and cell type | Shared genes within a cell type |

Each RDS is a named list with one element per dataset. Each element is a data
frame with genes in rows and methods in columns. A method gene that is not
reported is `NA`. The text table has the columns `method`, `cell_type`,
`n_genes`, `intra_cell_type_cor`, and `dataset`.

The gene universe for a dataset is the intersection of genes in the truth
sample-mean matrix and its assigned independent-reference signature. Cell
types are the mapped target cell types present in both the truth and the
assigned signature. Gene-level scores require a method to provide every
matched cell type and finite truth/inferred values for that gene; otherwise
that method-gene score is `NA`.

The first two scores answer a gene-level question: across cell types, did the
inferred mean profile preserve the truth profile? The third score changes the
axis of comparison and asks, within one cell type, whether the inferred values
rank genes similarly to truth.

### Truth-weighted, tie-aware Kendall-style ordering score

This is the full name of the score stored in
`truth_weighted_ordering_score_summary_list.RDS`. It is Kendall-style because
it evaluates pairwise ordering across cell types, but it is not ordinary
Kendall's tau: pair contributions are weighted by the corresponding absolute
truth difference and inferred ties receive zero credit.

For one dataset, method, and gene `g`, let `T[g,c]` be the truth sample mean
and `I[g,c]` the inferred sample mean in matched cell type `c`. For every
matched cell-type pair `(a, b)`, define:

~~~text
d_truth[g, a, b] = truth[g, a] - truth[g, b]
d_inferred[g, a, b] = inferred[g, a] - inferred[g, b]
weight[g, a, b] = abs(d_truth[g, a, b])

truth-weighted ordering score[g] =
  sum_{a < b} weight[g, a, b] *
              sign(d_truth[g, a, b]) *
              sign(d_inferred[g, a, b])
  / sum_{a < b} weight[g, a, b]
~~~

Equivalently, in mathematical notation:

For gene *g* and every cell-type pair *a*, *b*:

<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">
  <msub>
    <mi>S</mi>
    <mi>g</mi>
  </msub>
  <mo>=</mo>
  <mfrac>
    <mrow>
      <munder>
        <mo>∑</mo>
        <mrow>
          <mi>a</mi>
          <mo>&lt;</mo>
          <mi>b</mi>
        </mrow>
      </munder>
      <mo>|</mo>
      <msub>
        <mi>T</mi>
        <mrow><mi>g</mi><mi>a</mi></mrow>
      </msub>
      <mo>−</mo>
      <msub>
        <mi>T</mi>
        <mrow><mi>g</mi><mi>b</mi></mrow>
      </msub>
      <mo>|</mo>
      <mspace width="0.4em"/>
      <mi mathvariant="normal">sign</mi>
      <mo>(</mo>
      <msub>
        <mi>T</mi>
        <mrow><mi>g</mi><mi>a</mi></mrow>
      </msub>
      <mo>−</mo>
      <msub>
        <mi>T</mi>
        <mrow><mi>g</mi><mi>b</mi></mrow>
      </msub>
      <mo>)</mo>
      <mspace width="0.4em"/>
      <mi mathvariant="normal">sign</mi>
      <mo>(</mo>
      <msub>
        <mi>I</mi>
        <mrow><mi>g</mi><mi>a</mi></mrow>
      </msub>
      <mo>−</mo>
      <msub>
        <mi>I</mi>
        <mrow><mi>g</mi><mi>b</mi></mrow>
      </msub>
      <mo>)</mo>
    </mrow>
    <mrow>
      <munder>
        <mo>∑</mo>
        <mrow>
          <mi>a</mi>
          <mo>&lt;</mo>
          <mi>b</mi>
        </mrow>
      </munder>
      <mo>|</mo>
      <msub>
        <mi>T</mi>
        <mrow><mi>g</mi><mi>a</mi></mrow>
      </msub>
      <mo>−</mo>
      <msub>
        <mi>T</mi>
        <mrow><mi>g</mi><mi>b</mi></mrow>
      </msub>
      <mo>|</mo>
    </mrow>
  </mfrac>
</math>

- **S<sub>g</sub> = 1**: all high/low cell-type relationships are preserved.
- **S<sub>g</sub> = 0**: no trend is recovered, including a constant inferred profile.
- **S<sub>g</sub> = −1**: the ordering is completely reversed.
- Large truth differences receive more weight than negligible differences.
- Shifts and positive rescaling of a method do not affect the score.



The summary contains one `S_g` for each available dataset-gene-method
combination.

A correctly ordered pair contributes its positive truth weight, a reversed
pair contributes the negative weight, and an inferred tie contributes zero.
A truth tie has zero weight, so it neither helps nor hurts the score. The
score ranges from `-1` to `1`: `1` means that every truth ordering is
preserved, `0` means no net ordering information, and `-1` means that every
non-tied truth ordering is reversed. If the truth is constant across all
matched cell types, the denominator is zero and the score is `NA`.

For example, suppose one gene has truth means `Ast = 10`, `Exc = 6`, and
`Mic = 2`. The truth weights for `(Ast, Exc)`, `(Ast, Mic)`, and `(Exc, Mic)`
are `4`, `8`, and `4`, giving a denominator of `16`.

| Inferred values: Ast, Exc, Mic | Weighted numerator | Score | Interpretation |
|---|---:|---:|---|
| `8, 5, 1` | `4 + 8 + 4` | `1` | All three relationships are preserved. |
| `8, 8, 1` | `0 + 8 + 4` | `0.75` | Ast and Exc are tied; the other two relationships are preserved. |
| `5, 5, 5` | `0 + 0 + 0` | `0` | No cell-type ordering is inferred. |
| `1, 5, 9` | `-4 - 8 - 4` | `-1` | All three relationships are reversed. |
| `7, 9, 1` | `-4 + 8 + 4` | `0.5` | Ast versus Exc is reversed, while the two larger truth separations are preserved. |

This score deliberately evaluates ordering rather than absolute calibration.
It is also truth weighted: reversing a cell-type pair with a large truth
difference is penalized more than reversing a nearly tied pair.

### Concordance correlation coefficient (CCC)

For each gene, Lin's CCC compares its truth and inferred profiles across the
matched cell types. Using population moments across those cell types:

~~~text
CCC[g] = 2 * covariance(truth[g, ], inferred[g, ]) /
         (variance(truth[g, ]) + variance(inferred[g, ]) +
          (mean(truth[g, ]) - mean(inferred[g, ]))^2)
~~~

CCC ranges from `-1` to `1`. Unlike the truth-weighted ordering score, it
penalizes disagreement in both scale and mean level as well as association.
An exact match has CCC `1`; a completely reversed profile with the same mean
and variance has CCC `-1`. A variable truth profile compared with a constant
inferred profile has CCC `0`. If the denominator is zero, including when both
profiles are the same constant, CCC is `NA`.

As a worked example, use truth values `(10, 6, 2)` and inferred values
`(12, 8, 4)`. Their means are `6` and `8`; both population variances and their
population covariance are `32/3`. Therefore:

~~~text
CCC = (2 * 32/3) / (32/3 + 32/3 + (6 - 8)^2)
    = 16/19
    = 0.842
~~~

The cell-type trend is perfect, but the inferred profile is shifted upward by
two, so CCC is below one. For the same truth profile, inferred `(10, 6, 2)`
has CCC `1`, inferred `(2, 6, 10)` has CCC `-1`, and inferred `(5, 5, 5)` has
CCC `0`.

### Intra-cell-type Spearman correlation

For each dataset, method, and matched cell type, `intra_cell_type_cor` is the
Spearman correlation between truth and inferred sample means across their
finite shared genes:

~~~text
intra_cell_type_cor[c] =
  cor_g(truth[g, c], inferred[g, c], method = "spearman")
~~~

`n_genes` is the number of finite shared gene pairs used. The correlation is
`NA` when fewer than two genes are available or either vector is constant.
For a toy cell type with truth values `(1, 5, 9)` for genes `(G1, G2, G3)`,
inferred values `(10, 20, 30)` give identical ranks and Spearman correlation
`1`; inferred values `(30, 20, 10)` give `-1`; inferred values `(4, 4, 4)`
give `NA`.

## Covariance Specificity

`covariance_specificity_summary/summarize_gene_celltype_covariance_specificity.R`
summarizes config01 gene-level Pearson correlations across samples between
pairs of cell types. For each gene and focal cell type, it compares the median
absolute correlation with partner cell types in the method and truth results.
Genes come from the truth/independent-reference overlap and must also be
present in the method result. Both medians use the same partner pairs shared
by truth and the method; the default `require_all_partner_pairs = TRUE`
requires finite values for every one of those shared pairs.

The score is defined as:

~~~text
truth_raw_specificity  = 1 - median(abs(truth partner correlations))
method_raw_specificity = 1 - median(abs(method partner correlations))

relative_covariance_specificity
  = min(1, max(0, method_raw_specificity / truth_raw_specificity))
~~~

When `truth_raw_specificity` is zero, the relative score is set to 1. Otherwise,
a score below 1 indicates greater median absolute correlation between cell
types than in truth. A score of 1 includes both truth-matching and lower
absolute correlations, so it does not imply exact recovery of the truth
correlation structure.

The two outputs are:

~~~text
covariance_specificity_summary/config01__gene_celltype_specificity.RDS
covariance_specificity_summary/config01__gene_celltype_specificity_method_summary.RDS
~~~

The first is a data frame with one row per scored dataset, method, gene, and
focal cell type. The second summarizes scores by dataset, method, and focal
cell type, including their mean and median, `n_genes_scored`, and
`gene_score_coverage`. Coverage divides the scored gene count by the full
truth/independent-reference gene count; missing or unevaluable genes reduce
coverage and are omitted from score averages. Figure 9 uses the detailed file,
selects its marker genes, and averages scores within cell type and then within
dataset.

## Truth-Independent Gene Prioritization Inputs

The config-level gene-prioritization matrices are generated by:

~~~text
DALE_Eval/eval/build_gene_prioritization_scores.R
~~~

All outputs have genes as rows and benchmark cell types as columns. The five
method-specific scores use only inferred CTSE from the requested deconvolution
config; they do not read truth CTSE. Method cell types are first restricted and
mapped to the benchmark cell types.

Let `X_c[g, s]` be the inferred value for gene `g`, sample `s`, and cell type
`c` in the method's CTSE output. Before sample-level values are used, each
gene/cell-type profile containing negative values is shifted by one constant
across all samples:

~~~text
X_shifted_c[g, s] = X_c[g, s] + max(0, -min_s(X_c[g, s]))
~~~

The sample-by-gene matrix `Y_c` then follows the method preprocessing table
above: `log1p(X_shifted_c)` for methods in the `log1p` group, or unchanged
`X_shifted_c` for methods marked `already_transformed`. With
`lib_norm=true`, supported methods instead receive per-sample library
normalization to the configured scale factor followed by `log1p`.

### Mean-log contrast rankings

`post_hoc_meanlog_margin` and `post_hoc_meanlog_mean_contrast` directly read
the existing `sample_mean/<method>[_libnorm].txt` matrix. Their input is
therefore the stored six-decimal value

~~~text
M[g, c] = mean_s(Y_c[s, g])
~~~

rather than the original sample-level CTSE. The ranking step applies no
additional log transformation, normalization, standardization, or truth
comparison. In particular, `meanlog` does not imply another log step for
methods whose configured preprocessing is `already_transformed`.

~~~text
post_hoc_meanlog_margin[g, c]
  = M[g, c] - max_{d != c}(M[g, d])

post_hoc_meanlog_mean_contrast[g, c]
  = M[g, c] - mean_{d != c}(M[g, d])
~~~

Larger positive values indicate that the gene's average inferred expression is
more specific to the target cell type than to the other benchmark cell types.

### Log-mean contrast rankings

`post_hoc_logmean_margin` and `post_hoc_logmean_contrast` read the standard
`sample_logmean/<method>.txt` matrix described above. Let its value be
`L[g, c]`. The scores are

~~~text
post_hoc_logmean_margin[g, c]
  = L[g, c] - max_{d != c}(L[g, d])

post_hoc_logmean_contrast[g, c]
  = L[g, c] - mean_{d != c}(L[g, d])
~~~

For a nonnegative method in the `log1p` preprocessing group, the margin is the
archived log2 ratio between the target cell-type arithmetic mean and the
largest competing cell-type arithmetic mean. The contrast compares the target
mean with the geometric mean of the competing cell-type means. For methods
marked `already_transformed`, both scores are differences on the original
transformed scale. These rankings are available only for `lib_norm=false`.

### Partial-R2 ranking

`post_hoc_partial_r2` does not use the stored sample-mean matrix. It rereads
the sample-level CTSE output and rebuilds each `Y_c` using the shift and method
preprocessing above. It uses the samples stored by the method, verifies that
they are testing samples, and aligns them to the bulk samples. Genes must occur
in the bulk matrix and every selected CTSE cell type; a gene is retained only
when its config-prepared bulk value is positive in at least two aligned
samples.

Let `B_config[g, s]` be bulk expression after the requested config's bulk
preparation. `bulk_normalization=none` keeps the input values; `cpm` converts
count input to CPM or keeps CPM input; and `tmm`/`uq` convert count input to
edgeR-normalized CPM. Partial R2 applies natural-log `log1p` to `B_config` and
fits an intercept-containing joint model:

~~~text
log1p(B_config[g, ]) ~ Y_1[, g] + ... + Y_C[, g]

post_hoc_partial_r2[g, c]
  = t_c^2 / (t_c^2 + residual_df)
~~~

Here `t_c` is the coefficient t-statistic for cell type `c`. The score is the
additional variance associated with that predictor after accounting for all
other cell types. It ranges from 0 to 1 and discards coefficient direction;
larger is better. Fractions are not used for this score.

### Adapted GeneSigTest ranking

`genesigtest_adapted_fdr` is config-level and method-independent: it does not
use inferred CTSE. It uses `B_config` on its linear prepared scale, without
the extra `log1p` used for partial R2. On testing samples shared with the
configured fraction matrix `F`, it fits every bulk gene with

~~~text
B_config[g, ] ~ F[, 1] + ... + F[, C] - 1
~~~

For each cell-type coefficient, the script calculates a one-sided p-value for
a positive association and applies BH adjustment across genes separately
within each cell type. It does not apply the partial-R2 positive-in-two-samples
filter. Independent-reference fraction columns are mapped to benchmark target
cell types for output. The stored value is raw adjusted FDR rounded to
`--digits` significant digits, so smaller is better. This is the deterministic,
non-refined approximation of archived `GeneSigTest_adapted()`; it does not call
ENIGMA's optional refinement or bootstrap procedure.

Outputs are stored without a truth suffix:

~~~text
Benchmarking_obj/<dataset>/gene_prioritization/<config_slug>/post_hoc_meanlog_margin/<method>[_libnorm].txt.gz
Benchmarking_obj/<dataset>/gene_prioritization/<config_slug>/post_hoc_meanlog_mean_contrast/<method>[_libnorm].txt.gz
Benchmarking_obj/<dataset>/gene_prioritization/<config_slug>/post_hoc_logmean_margin/<method>.txt.gz
Benchmarking_obj/<dataset>/gene_prioritization/<config_slug>/post_hoc_logmean_contrast/<method>.txt.gz
Benchmarking_obj/<dataset>/gene_prioritization/<config_slug>/post_hoc_partial_r2/<method>[_libnorm].txt.gz
Benchmarking_obj/<dataset>/gene_prioritization/<config_slug>/genesigtest_adapted_fdr/GeneSigTest_adapted.txt.gz
~~~

## Post-hoc Versus Independent-Reference Gene-Prioritization Summary

Script:

~~~text
gene_prioritization_summary/summarize_posthoc_vs_indep_ref.R
~~~

Output:

~~~text
gene_prioritization_summary/posthoc_vs_indep_ref_summary_list.RDS
~~~

This summary evaluates whether genes prioritized from inferred CTSE recover
genes with strong truth-CTSE correlation as effectively as genes prioritized
from the assigned independent reference. The default run uses `config01`, all
datasets with gene-prioritization outputs, all methods with a matching
Spearman-correlation matrix, and top-N groups of 10, 30, 100, 300, 1,000, and
3,000 genes.

The default alternative rankings are:

~~~text
meanlog_margin
meanlog_mean_contrast
logmean_margin
logmean_contrast
partial_r2
genesigtest_adapted_fdr
~~~

The first five are method-specific post-hoc rankings. Adapted GeneSigTest is
config-level and method-independent, but its selected genes are evaluated
separately with each CTSE method's gene-level Spearman correlations. Larger
values rank first for the post-hoc and independent-reference limma scores;
smaller adjusted FDR values rank first for GeneSigTest.

The comparison is pairwise. For each alternative ranking, cell type, and
method, the alternative and mapped independent-reference limma ranking use the
same truth-independent candidate universe:

~~~text
finite alternative-score genes
intersect finite mapped independent-reference limma genes
intersect genes in the configured bulk input
~~~

This pairwise design prevents the more restrictive partial-R2 gene universe
from changing the candidate genes used for the mean-log, log-mean, GeneSigTest,
or other comparisons. Independent-reference rows are therefore repeated once
for each `comparison_score`; each repeated row is the matched reference
baseline for that alternative ranking.

Top N is selected from the matched candidate universe before intersecting with
the method correlation result and before removing truth-side acceptable NAs.
The remaining non-finite correlations are omitted from `avg_cor` and replaced
with zero in `avg_cor_with_NA_penalty`.

The main ranking columns are:

| Column | Definition |
|---|---|
| `ranking_source` | Broad source family: `indep_ref`, `posthoc`, or `genesigtest`. |
| `ranking_score` | Exact ranking used for the row: `indep_ref_limma`, `meanlog_margin`, `meanlog_mean_contrast`, `logmean_margin`, `logmean_contrast`, `partial_r2`, or `genesigtest_adapted_fdr`. |
| `comparison_score` | Alternative score defining the pairwise comparison. For an `indep_ref_limma` row, this identifies the alternative ranking whose candidate universe is used. |
| `ranking_label` | Plot-ready label for the ranking. |
| `group` / `top_n` | Requested top-gene group and its numeric N. |

Gene-availability and evaluation columns are:

| Column | Definition |
|---|---|
| `n_candidate_genes` | Genes in the pair-specific truth-independent candidate universe. |
| `n_selected_genes` | Candidate genes selected for the requested top N; it can be smaller than N when the candidate universe is smaller. |
| `n_cor_available_genes` | Selected genes present in the CTSE method's correlation matrix. |
| `n_evaluated_genes` | Correlation-available genes remaining after truth-side acceptable-NA removal. |
| `n_constant_genes` | Evaluated genes with a remaining non-finite method correlation. |
| `avg_cor` | Mean of the finite correlations among evaluated genes. |
| `avg_cor_with_NA_penalty` | Mean over evaluated genes after replacing remaining non-finite method correlations with zero. |

The output RDS is a named list with one data frame per dataset. The script also
leaves `summary_df`, the combined data frame, and `availability_summary`, the
dataset/config availability table, inspectable in the R session. Only
`summary_list` is written to the RDS.

This summary intentionally reads the ordinary unsuffixed post-hoc files, such
as `<method>.txt.gz`. It has no `lib_norm` setting and does not use
`<method>_libnorm.txt.gz` outputs.

## DE Effect-Size and Reported-q Call Agreement

`DE_summary/DEG_call_agreement_summary.R` compares cell-type-specific DE
results with the truth-CTSE DE reference. The upregulated and downregulated
summaries evaluate calls separately by direction. A positive call requires
finite `effect_size`, `p_value`, and reported `q_value`, with the default rules:

~~~text
up:   q_value < 0.05 and effect_size > 0
down: q_value < 0.05 and effect_size < 0
~~~

The workflow uses each result's reported q-values without recalculating or
harmonizing them. It selects the direction after applying the reported
significance threshold; it does not recalculate one-sided tests. Methods may
have adjusted over different tested gene families, and these summaries
preserve their reported decisions. In each directional summary, TP already
means that truth and result are significant in the same selected direction.

The workflow also calculates effect-size Spearman correlation over
truth-q-evaluable genes in the assigned independent reference that are tested
by the result. This comparison includes all paired genes regardless of DEG
significance or effect direction.

### Current exported scope

The reusable directional files are:

- [`DEG_call_agreement_summary_up.tsv`](DE_summary/DEG_call_agreement_summary_up.tsv)
- [`DEG_call_agreement_summary_down.tsv`](DE_summary/DEG_call_agreement_summary_down.tsv)

Each currently contains 253 rows and 33 columns:

- Datasets: `PBMC_Perez2022` (92 rows) and `ROSMAP_AD430_Mathys2023` (161 rows).
- PBMC cell types: B cell, T cell, NK, and Myeloid.
- ROSMAP cell types: Ast, Exc, Inh, Mic, Oli, Opc, and Vas.
- Two-stage family: 10 CTSE methods followed by DE testing, with two
  truth-universe rows per method and cell type.
- Direct family: bMIND, ENIGMAL2, and ENIGMAtrace, with one truth-universe row
  per method and cell type.

Rows are identified by `dataset`, `cell_type`, `compared_family`,
`compared_method`, `compared_result`, and `truth_universe_type`. The filename
identifies the call direction. Figure 6b and Figure S7 currently read the
upregulated file. Selecting one two-stage universe and the single direct
universe gives 52 PBMC method-family-cell-type comparisons, avoiding duplicate
weighting of two-stage methods.

The directional tables remain available as `deg_call_agreement_up` and
`deg_call_agreement_down` in the R session. `result_index` lists existing
method/cell-type result files, and `dataset_result_list` contains the
per-dataset results. The current setting is `write_outputs = TRUE`, which
writes the agreement TSVs and `DEG_call_metric_spec.tsv` beside the script.
Additional combined-direction exports are described at the end of this section.

### Truth-universe row design

A gene is tested when `effect_size` and `p_value` are finite, and q-evaluable
when its reported `q_value` is also finite. Let `Q` be the truth-q-evaluable
genes, `I` the assigned independent-reference genes, and `C` the genes exported
by a two-stage CTSE method.

| Family | `truth_universe_type` | Gene universe |
|---|---|---|
| Direct | `truth_q_evaluable_indep_ref` | `Q` intersect `I`. |
| Two-stage | `truth_q_evaluable_indep_ref` | `Q` intersect `I`. |
| Two-stage | `truth_q_evaluable_indep_ref_ctse` | `Q` intersect `I` intersect `C`. |

The CTSE DE workflow retains every CTSE input gene as a DE-result row,
including untested genes with NA statistics. Its result-table gene identifiers
therefore provide `C`. The broader independent-reference universe evaluates
the full two-stage workflow, including truth positives absent from CTSE output.
The CTSE-overlap universe evaluates DE agreement conditional on genes exported
by CTSE. Genes unavailable in the independent reference are excluded from both
universes.

### Coverage column reference

| Statistic | Meaning |
|---|---|
| `n_truth_input_genes` | Gene rows present in the truth DE table, including genes with NA DE statistics. |
| `n_result_input_genes` | Gene rows present in the compared result table. |
| `n_common_input_genes` | Gene identifiers present in both input tables. |
| `n_truth_tested` | Truth genes with finite `effect_size` and `p_value`. |
| `n_result_tested` | Result genes with finite `effect_size` and `p_value`. |
| `n_truth_indep_ref_q_evaluable` | Truth-q-evaluable genes in the assigned independent reference: `Q` intersect `I`. |
| `n_result_indep_ref_input` | Result rows whose genes are in `I`. |
| `n_result_indep_ref_tested` | Result rows in `I` with finite `effect_size` and `p_value`. |
| `n_truth_indep_ref_tested_by_result` | Genes in `Q` intersect `I` tested by the result; these are the paired genes for effect-size Spearman. |
| `result_indep_ref_tested_rate` | `n_result_indep_ref_tested / n_result_indep_ref_input`. |
| `truth_indep_ref_tested_coverage` | `n_truth_indep_ref_tested_by_result / n_truth_indep_ref_q_evaluable`. |
| `coverage_warning` | Flags a reference-aware coverage rate below 0.50 without removing the comparison. |

These coverage statistics are independent of call direction and are repeated
across the two universe rows for a two-stage result. A coverage rate is `NA`
when its denominator is zero. `result_indep_ref_tested_rate` describes
testability among reference-eligible rows present in the result;
`truth_indep_ref_tested_coverage` additionally accounts for missing result
genes, using all of `Q` intersect `I` as its denominator.

### Why truth DE tables can contain NA genes

Truth DE tables retain the complete input-gene list even when some genes
cannot be tested. Genes excluded by expression filtering or with insufficient
information for the DE model remain as rows with unavailable statistics.
They contribute to `n_truth_input_genes` but not to the q-evaluable truth
universe. An NA row therefore indicates an unevaluable gene rather than a
failure of the whole cell-type analysis.

### DEG-call count column reference

The following columns apply to both the up and down files. Positive and
negative refer to calls in the direction selected by the file.

| Statistic | Meaning |
|---|---|
| `truth_universe_type` | Definition of the truth universe used by this row. |
| `n_truth_universe` | Number of genes in that universe: `TP + FP + FN + TN`. |
| `n_truth_sig` | Truth positives within the selected universe and direction: `TP + FN`. |
| `n_result_sig_total` | All result positives in the selected direction, including those outside the truth universe. |
| `n_result_sig_in_truth_universe` | Result positives inside the selected universe: `TP + FP`. |
| `n_result_sig_outside_truth_universe` | `n_result_sig_total - n_result_sig_in_truth_universe`; recorded but excluded from FP and the agreement metrics. |
| `n_shared_deg` | TP: truth and result are both positive in the selected direction within the universe. |
| `n_result_only_deg` | FP: result positive but truth not positive in the selected direction within the universe. |
| `n_missed_truth_deg` | FN: truth positive but result not positive in the selected direction within the universe. |

TN is calculated internally but has no exported count column. It can be
recovered as:

~~~text
TN = n_truth_universe - n_shared_deg - n_result_only_deg - n_missed_truth_deg
~~~

For two-stage results, the difference in `n_missed_truth_deg` between the two
universe rows counts reference-eligible truth positives absent from CTSE output.

### Reference universe and no-call rule

Within the selected universe, a missing result gene, unavailable DE statistics,
a nonsignificant q-value, or an effect outside the selected direction is a
negative/no-call. Such a gene is an FN when truth is positive and a TN when
truth is negative. A truth-negative gene can include a significant truth DEG
in the opposite direction.

For example, a significant truth-up/result-down gene is an FN in the up
summary and an FP in the down summary, provided it belongs to the selected
universe. The CTSE-overlap universe excludes genes absent from `C` before
counting calls; the broader independent-reference universe retains them as
possible FNs or TNs.

### Heatmap-ready performance column reference

| Statistic | Interpretation |
|---|---|
| `effect_size_spearman` | Per-cell-type Spearman correlation over genes in `Q` intersect `I` tested by the result, regardless of DEG significance or direction. |
| `precision_vs_truth` | Fraction of result positives inside the universe also positive in truth: `TP / (TP + FP)`. |
| `recall_vs_truth` | Fraction of truth positives recovered: `TP / (TP + FN)`. |
| `f1_vs_truth` | Balance of precision and recall: `2*TP / (2*TP + FP + FN)`. |
| `call_jaccard` | Intersection divided by union of positive calls: `TP / (TP + FP + FN)`. |
| `specificity_vs_truth` | Fraction of truth negatives receiving no positive result call: `TN / (TN + FP)`. |
| `false_positive_rate_vs_truth` | Fraction of truth negatives called positive: `FP / (TN + FP)`. |

Ratios are `NA` when their denominator is zero. In particular, precision is
`NA` when there are no result positives, recall is `NA` when there are no truth
positives, and F1/Jaccard are `NA` when neither has positives. If one side has
positives and the other has none, F1/Jaccard are zero. Figure 6b's scatter view
uses zero to display undefined precision for no-call results; the exported
TSVs retain `NA`.

For an upregulated-call example with ten universe genes, suppose four are
truth-up DEGs and the result calls three genes up, of which two are truth-up
DEGs. Then `TP = 2`, `FP = 1`, `FN = 2`, and `TN = 5`:

| Metric | Calculation | Value |
|---|---:|---:|
| Precision | `2 / 3` | `0.667` |
| Recall | `2 / 4` | `0.500` |
| F1 | `4 / 7` | `0.571` |
| Jaccard | `2 / 5` | `0.400` |
| Specificity | `5 / 6` | `0.833` |
| False-positive rate | `1 / 6` | `0.167` |

The same definitions apply to downregulated calls with the direction reversed.

### Effect-size rank agreement

For each dataset, cell type, and method, the paired gene set `G` consists of
genes in `Q` intersect `I` tested by the result:

~~~text
effect_size_spearman = Spearman correlation(
  truth effect_size[G], result effect_size[G]
)
~~~

`n_truth_indep_ref_tested_by_result` records the paired gene count. This
correlation is the same in the up/down files and the two two-stage universe
rows. It is `NA` with fewer than two paired genes or a constant effect-size
vector. It measures agreement in effect ranking across genes without requiring
comparable effect magnitudes or scales.

### Coverage, convergence, and interpretation

Precision, recall, and F1 describe agreement with the empirical truth-CTSE DE
reference in the selected direction and universe. Their interpretation depends
on gene coverage: a method can have high precision among a small tested subset
while missing many reference-eligible truth positives. The two coverage rates
and the broader truth-universe view expose this difference.

The result index contains existing files only. Missing method/cell-type
combinations are absent from the summaries, so their completion status is
separate from the reported agreement metrics. A missing file alone does not
establish whether a method was unrun or failed to converge.

Specificity and false-positive rate are supporting diagnostics because
missing or untestable result genes count as no-calls and can inflate apparent
true negatives. Jaccard and F1 contain redundant confusion-count information.
Two-stage and direct DE families remain labeled separately because their
models and testing targets may differ. Agreement precision is relative to the
empirical truth reference, rather than an absolute biological true discovery
rate.

### Additional combined-direction exports

The generator also writes `DEG_call_agreement_summary_both.tsv` and its
identical unsuffixed copy, `DEG_call_agreement_summary.tsv`. These count
significant calls regardless of effect sign and retain additional
direction-agreement columns absent from the 33-column up/down files. The in-session
`deg_call_agreement` alias and `deg_call_agreement_long` also refer to this
combined summary. `DEG_call_metric_spec.tsv` includes combined-summary metrics
and is therefore broader than the up/down column references above.
