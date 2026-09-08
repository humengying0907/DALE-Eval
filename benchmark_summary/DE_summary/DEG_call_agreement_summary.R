# Effect-size and reported-q DEG-call agreement between inferred/direct results
# and truth CTSE.
#
# Notebook-style usage:
# 1. Run from the repository root, benchmark_summary/, or
#    benchmark_summary/DE_summary/.
# 2. Edit the settings below.
# 3. Inspect dataset_result_list, result_index, deg_call_agreement_both,
#    deg_call_agreement_up, deg_call_agreement_down, deg_call_agreement_long,
#    deg_call_metric_spec, and paper_metric_notes.
# 4. Set write_outputs = TRUE only when reusable TSV exports are wanted.
#
# Important: this script never recalculates BH q-values. Each result is called
# significant from the q_value reported in its own DE table.


## settings ----

dataset_include = c(
  "PBMC_Perez2022",
  "ROSMAP_AD430_Mathys2023"
)
de_q_cutoff = 0.05

truth_folder = "truth_sumcount"
two_stage_folder =
  "config02_bulk-sumcount__frac-InstaPrismfrac__ref-indep"
direct_folder = "direct_ctsDEG_config02_bulk-sumcount"

families_include = c("Two-stage", "Direct")
call_directions = c("both", "up", "down")

# These thresholds create warnings/labels only. They do not remove results.
result_indep_ref_tested_rate_warning = 0.50
truth_indep_ref_tested_coverage_warning = 0.50

write_outputs = T


## paths ----

working_dir = normalizePath(getwd(), mustWork = TRUE)

if (
  dir.exists(file.path(working_dir, "Benchmarking_obj")) &&
    dir.exists(file.path(working_dir, "benchmark_summary"))
) {
  repo_root = working_dir
} else if (basename(working_dir) == "benchmark_summary") {
  repo_root = normalizePath(file.path(working_dir, ".."), mustWork = TRUE)
} else if (
  basename(working_dir) == "DE_summary" &&
    basename(dirname(working_dir)) == "benchmark_summary"
) {
  repo_root = normalizePath(file.path(working_dir, "../.."), mustWork = TRUE)
} else {
  stop(
    "Run this script from the repository root, benchmark_summary/, or ",
    "benchmark_summary/DE_summary/."
  )
}

output_dir = file.path(repo_root, "benchmark_summary", "DE_summary")

ref_assignment = read.delim(
  file.path(
    repo_root,
    "DALE_Eval",
    "configs",
    "benchmark_ref_assignment.txt"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)


## helpers ----

read_de_result = function(path, dataset, family, method, cell_type) {
  x = read.delim(
    gzfile(path),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_columns = c("gene", "effect_size", "p_value", "q_value")
  missing_columns = setdiff(required_columns, names(x))

  if (length(missing_columns) > 0L) {
    stop(
      "Missing columns in ", path, ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if (anyDuplicated(x$gene)) {
    stop("Duplicated gene names in ", path)
  }

  x = x[, required_columns]
  x$dataset = dataset
  x$family = family
  x$method = method
  x$result = if (family == "Truth") {
    "Truth"
  } else {
    paste(family, method, sep = ": ")
  }
  x$cell_type = cell_type
  x$tested = is.finite(x$effect_size) & is.finite(x$p_value)
  x$q_available = x$tested & is.finite(x$q_value)
  x$reported_significant = x$q_available & x$q_value < de_q_cutoff

  x
}

safe_ratio = function(numerator, denominator) {
  if (denominator > 0) numerator / denominator else NA_real_
}

direction_positive = function(effect_size, direction) {
  switch(
    direction,
    both = rep(TRUE, length(effect_size)),
    up = is.finite(effect_size) & effect_size > 0,
    down = is.finite(effect_size) & effect_size < 0
  )
}

compare_reported_deg_calls = function(
    truth,
    result,
    indep_ref_genes,
    direction) {
  truth_q_reference = truth[truth$q_available, , drop = FALSE]
  truth_indep_reference = truth_q_reference[
    truth_q_reference$gene %in% indep_ref_genes,
    ,
    drop = FALSE
  ]
  result_indep_reference = result[
    result$gene %in% indep_ref_genes,
    ,
    drop = FALSE
  ]

  result_tested_genes = result$gene[result$tested]
  common_tested_indep_ref_genes = truth_indep_reference$gene[
    truth_indep_reference$gene %in% result_tested_genes
  ]
  truth_effect = truth_indep_reference$effect_size[
    match(common_tested_indep_ref_genes, truth_indep_reference$gene)
  ]
  result_effect_common = result$effect_size[
    match(common_tested_indep_ref_genes, result$gene)
  ]
  effect_size_spearman = if (
    length(common_tested_indep_ref_genes) > 1L
  ) {
    suppressWarnings(cor(
      truth_effect,
      result_effect_common,
      method = "spearman"
    ))
  } else {
    NA_real_
  }

  truth_universe_list = if (result$family[1] == "Two-stage") {
    list(
      truth_q_evaluable_indep_ref = truth_indep_reference,
      truth_q_evaluable_indep_ref_ctse = truth_indep_reference[
        truth_indep_reference$gene %in% result$gene,
        ,
        drop = FALSE
      ]
    )
  } else {
    list(truth_q_evaluable_indep_ref = truth_indep_reference)
  }

  n_common_input_genes = length(intersect(truth$gene, result$gene))
  n_result_indep_ref_tested = sum(result_indep_reference$tested)
  n_truth_indep_ref_tested_by_result =
    length(common_tested_indep_ref_genes)
  n_result_sig_total = sum(
    result$reported_significant &
      direction_positive(result$effect_size, direction)
  )

  do.call(rbind, lapply(names(truth_universe_list), function(
      truth_universe_type) {
    truth_universe = truth_universe_list[[truth_universe_type]]
    result_match = match(truth_universe$gene, result$gene)
    result_present = !is.na(result_match)

    result_positive = rep(FALSE, nrow(truth_universe))
    result_effect = rep(NA_real_, nrow(truth_universe))
    result_positive[result_present] =
      result$reported_significant[result_match[result_present]] &
      direction_positive(
        result$effect_size[result_match[result_present]],
        direction
      )
    result_effect[result_present] =
      result$effect_size[result_match[result_present]]

    truth_positive = truth_universe$reported_significant &
      direction_positive(truth_universe$effect_size, direction)
    shared_deg = truth_positive & result_positive
    result_only_deg = !truth_positive & result_positive
    missed_truth_deg = truth_positive & !result_positive
    shared_non_deg = !truth_positive & !result_positive

    direction_match = shared_deg &
      is.finite(result_effect) &
      sign(truth_universe$effect_size) == sign(result_effect)

    n_truth_positive = sum(truth_positive)
    n_result_positive = sum(result_positive)
    n_shared_deg = sum(shared_deg)
    n_result_only_deg = sum(result_only_deg)
    n_missed_truth_deg = sum(missed_truth_deg)
    n_shared_non_deg = sum(shared_non_deg)
    n_signed_shared_deg = sum(direction_match)

    precision = safe_ratio(n_shared_deg, n_result_positive)
    recall = safe_ratio(n_shared_deg, n_truth_positive)
    f1 = safe_ratio(
      2 * n_shared_deg,
      2 * n_shared_deg + n_result_only_deg + n_missed_truth_deg
    )

    signed_precision = safe_ratio(n_signed_shared_deg, n_result_positive)
    signed_recall = safe_ratio(n_signed_shared_deg, n_truth_positive)
    signed_f1 = safe_ratio(
      2 * n_signed_shared_deg,
      n_result_positive + n_truth_positive
    )

    data.frame(
      dataset = truth$dataset[1],
      cell_type = truth$cell_type[1],
      compared_family = result$family[1],
      compared_method = result$method[1],
      compared_result = result$result[1],
      truth_universe_type = truth_universe_type,
      n_truth_universe = nrow(truth_universe),
      n_truth_input_genes = nrow(truth),
      n_result_input_genes = nrow(result),
      n_common_input_genes = n_common_input_genes,
      n_truth_tested = sum(truth$tested),
      n_result_tested = sum(result$tested),
      n_truth_indep_ref_q_evaluable = nrow(truth_indep_reference),
      n_result_indep_ref_input = nrow(result_indep_reference),
      n_result_indep_ref_tested = n_result_indep_ref_tested,
      n_truth_indep_ref_tested_by_result =
        n_truth_indep_ref_tested_by_result,
      result_indep_ref_tested_rate = safe_ratio(
        n_result_indep_ref_tested,
        nrow(result_indep_reference)
      ),
      truth_indep_ref_tested_coverage = safe_ratio(
        n_truth_indep_ref_tested_by_result,
        nrow(truth_indep_reference)
      ),
      effect_size_spearman = effect_size_spearman,
      n_truth_sig = n_truth_positive,
      n_result_sig_total = n_result_sig_total,
      n_result_sig_in_truth_universe = n_result_positive,
      n_result_sig_outside_truth_universe =
        n_result_sig_total - n_result_positive,
      n_shared_deg = n_shared_deg,
      n_result_only_deg = n_result_only_deg,
      n_missed_truth_deg = n_missed_truth_deg,
      n_shared_deg_direction_match = n_signed_shared_deg,
      precision_vs_truth = precision,
      recall_vs_truth = recall,
      f1_vs_truth = f1,
      signed_precision_vs_truth = signed_precision,
      signed_recall_vs_truth = signed_recall,
      signed_f1_vs_truth = signed_f1,
      shared_deg_direction_agreement = safe_ratio(
        n_signed_shared_deg,
        n_shared_deg
      ),
      call_jaccard = safe_ratio(
        n_shared_deg,
        n_shared_deg + n_result_only_deg + n_missed_truth_deg
      ),
      specificity_vs_truth = safe_ratio(
        n_shared_non_deg,
        n_shared_non_deg + n_result_only_deg
      ),
      false_positive_rate_vs_truth = safe_ratio(
        n_result_only_deg,
        n_shared_non_deg + n_result_only_deg
      ),
      stringsAsFactors = FALSE
    )
  }))
}


## dataset summaries ----

build_dataset_summary = function(dataset_name) {
  indep_ref = ref_assignment$indep_ref[
    match(dataset_name, ref_assignment$dataset)
  ]
  indep_ref_signature = read.csv(
    file.path(
      repo_root,
      "Indep_scReference",
      indep_ref,
      "rowMeans_sig.csv"
    ),
    row.names = 1,
    check.names = FALSE
  )
  indep_ref_genes = rownames(indep_ref_signature)

  de_root = file.path(
    repo_root,
    "Benchmarking_obj",
    dataset_name,
    "DE_res"
  )
  truth_de_dir = file.path(de_root, truth_folder)
  two_stage_de_dir = file.path(de_root, two_stage_folder)
  direct_de_dir = file.path(de_root, direct_folder)

  stopifnot(
    dir.exists(truth_de_dir),
    dir.exists(two_stage_de_dir),
    dir.exists(direct_de_dir)
  )

  cell_types = sort(sub(
    "\\.txt\\.gz$",
    "",
    list.files(truth_de_dir, pattern = "\\.txt\\.gz$")
  ))

  result_index = data.frame(
    dataset = character(),
    family = character(),
    method = character(),
    cell_type = character(),
    path = character(),
    stringsAsFactors = FALSE
  )

  family_dirs = list(
    `Two-stage` = two_stage_de_dir,
    Direct = direct_de_dir
  )

  for (family in intersect(families_include, names(family_dirs))) {
    family_dir = family_dirs[[family]]
    methods = sort(list.dirs(
      family_dir,
      recursive = FALSE,
      full.names = FALSE
    ))

    for (method in methods) {
      for (cell_type in cell_types) {
        path = file.path(
          family_dir,
          method,
          paste0(cell_type, ".txt.gz")
        )

        if (file.exists(path)) {
          result_index = rbind(
            result_index,
            data.frame(
              dataset = dataset_name,
              family = family,
              method = method,
              cell_type = cell_type,
              path = path,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
  }

  truth_results = setNames(lapply(cell_types, function(cell_type) {
    read_de_result(
      file.path(truth_de_dir, paste0(cell_type, ".txt.gz")),
      dataset = dataset_name,
      family = "Truth",
      method = "Truth",
      cell_type = cell_type
    )
  }), cell_types)

  compared_results = lapply(seq_len(nrow(result_index)), function(i) {
    read_de_result(
      result_index$path[i],
      dataset = dataset_name,
      family = result_index$family[i],
      method = result_index$method[i],
      cell_type = result_index$cell_type[i]
    )
  })

  deg_call_agreement = setNames(lapply(call_directions, function(direction) {
    do.call(rbind, lapply(seq_along(compared_results), function(i) {
      result = compared_results[[i]]
      truth = truth_results[[result$cell_type[1]]]
      compare_reported_deg_calls(
        truth,
        result,
        indep_ref_genes,
        direction
      )
    }))
  }), call_directions)

  list(
    result_index = result_index,
    deg_call_agreement = deg_call_agreement
  )
}


## combine dataset summaries ----

dataset_result_list = setNames(
  lapply(dataset_include, build_dataset_summary),
  dataset_include
)

result_index = do.call(rbind, lapply(
  dataset_result_list,
  `[[`,
  "result_index"
))
rownames(result_index) = NULL

deg_call_agreement_list = setNames(lapply(call_directions, function(direction) {
  x = do.call(rbind, lapply(dataset_result_list, function(dataset_result) {
    dataset_result$deg_call_agreement[[direction]]
  }))
  rownames(x) = NULL

  x$coverage_warning =
    x$result_indep_ref_tested_rate <
      result_indep_ref_tested_rate_warning |
    x$truth_indep_ref_tested_coverage <
      truth_indep_ref_tested_coverage_warning

  x[order(
    x$dataset,
    x$cell_type,
    x$compared_family,
    x$compared_method,
    x$truth_universe_type
  ), ]
}), call_directions)

deg_call_agreement_both = deg_call_agreement_list$both

direction_redundant_columns = c(
  "n_shared_deg_direction_match",
  "signed_precision_vs_truth",
  "signed_recall_vs_truth",
  "signed_f1_vs_truth",
  "shared_deg_direction_agreement"
)

deg_call_agreement_up = deg_call_agreement_list$up[, setdiff(
  names(deg_call_agreement_list$up),
  direction_redundant_columns
)]
deg_call_agreement_down = deg_call_agreement_list$down[, setdiff(
  names(deg_call_agreement_list$down),
  direction_redundant_columns
)]

# Backward-compatible notebook object and unsuffixed export.
deg_call_agreement = deg_call_agreement_both


## metric specification for later heatmap work ----

deg_call_metric_spec = data.frame(
  id = c(
    "result_indep_ref_tested_rate",
    "truth_indep_ref_tested_coverage",
    "effect_size_spearman",
    "precision_vs_truth",
    "recall_vs_truth",
    "f1_vs_truth",
    "signed_precision_vs_truth",
    "signed_recall_vs_truth",
    "signed_f1_vs_truth",
    "shared_deg_direction_agreement",
    "specificity_vs_truth",
    "false_positive_rate_vs_truth",
    "call_jaccard"
  ),
  group = c(
    rep("Coverage", 2),
    "Effect-size agreement",
    rep("Reported-q call agreement", 7),
    rep("Paper diagnostics", 2),
    "Redundant diagnostic"
  ),
  label = c(
    "Result tested rate within independent reference",
    "Truth independent-reference tested coverage",
    "Effect-size Spearman",
    "Precision / TDR",
    "Recall / sensitivity",
    "F1",
    "Signed precision",
    "Signed recall",
    "Signed F1",
    "Shared-DEG direction",
    "Specificity",
    "False-positive rate",
    "DEG Jaccard"
  ),
  direction = c(
    rep("higher", 11),
    "lower",
    "higher"
  ),
  heatmap_role = c(
    "essential QC",
    "essential QC",
    "primary rank agreement",
    "primary",
    "primary",
    "primary unsigned summary",
    "supporting",
    "supporting",
    "candidate ranking metric",
    "supporting",
    "diagnostic only",
    "diagnostic only",
    "omit when F1 is shown"
  ),
  stringsAsFactors = FALSE
)

metric_ids = deg_call_metric_spec$id

deg_call_agreement_long = do.call(rbind, lapply(metric_ids, function(metric) {
  data.frame(
    dataset = deg_call_agreement$dataset,
    cell_type = deg_call_agreement$cell_type,
    compared_family = deg_call_agreement$compared_family,
    compared_method = deg_call_agreement$compared_method,
    compared_result = deg_call_agreement$compared_result,
    truth_universe_type = deg_call_agreement$truth_universe_type,
    n_truth_universe = deg_call_agreement$n_truth_universe,
    coverage_warning = deg_call_agreement$coverage_warning,
    metric = metric,
    value = deg_call_agreement[[metric]],
    stringsAsFactors = FALSE
  )
}))

deg_call_agreement_long = merge(
  deg_call_agreement_long,
  deg_call_metric_spec,
  by.x = "metric",
  by.y = "id",
  all.x = TRUE,
  sort = FALSE
)


## metrics from Meng et al. (Briefings in Bioinformatics, 2023) ----

# Paper: "A comprehensive assessment of cell type-specific differential
# expression methods in bulk data", doi:10.1093/bib/bbac516.
#
# Its simulation truth is known and all methods are evaluated over controlled
# gene universes. Here, truth CTSE DE is an empirical reference and method gene
# testability differs. The mappings below therefore distinguish directly useful
# call-agreement metrics from metrics needing a fixed-universe design.

paper_metric_notes = data.frame(
  paper_metric = c(
    "True discovery rate (TDR)",
    "Sensitivity",
    "Specificity",
    "False-positive rate",
    "ROC curve / AUROC",
    "TDR at fixed top-N genes",
    "Expression-stratified sensitivity",
    "Runtime"
  ),
  current_mapping = c(
    "precision_vs_truth",
    "recall_vs_truth",
    "specificity_vs_truth",
    "false_positive_rate_vs_truth",
    "not calculated",
    "not calculated",
    "future stratified summary",
    "outside DEG-call agreement"
  ),
  recommendation = c(
    "Primary; use reported-q calls in the selected truth-universe type",
    "Primary; untested method genes remain no-calls and can be missed truth DEGs",
    "Diagnostic only; untested method genes inflate apparent true negatives",
    "Diagnostic only; inherits the specificity/testability limitation",
    "Defer until a fixed gene universe and missing-q score policy are prespecified",
    "Do not use with variable available gene counts; it can reward low coverage",
    "Relevant future analysis because the paper finds strong expression-level effects",
    "Useful operational metric, but not evidence of biological DEG agreement"
  ),
  stringsAsFactors = FALSE
)


## inspect ----

agreement_preview_columns = c(
  "dataset",
  "cell_type",
  "compared_family",
  "compared_method",
  "truth_universe_type",
  "n_truth_universe",
  "n_result_tested",
  "n_truth_indep_ref_q_evaluable",
  "n_result_indep_ref_input",
  "n_result_indep_ref_tested",
  "n_truth_indep_ref_tested_by_result",
  "result_indep_ref_tested_rate",
  "truth_indep_ref_tested_coverage",
  "effect_size_spearman",
  "n_truth_sig",
  "n_result_sig_in_truth_universe",
  "n_shared_deg",
  "n_missed_truth_deg",
  "precision_vs_truth",
  "recall_vs_truth",
  "f1_vs_truth",
  "signed_f1_vs_truth",
  "coverage_warning"
)


print(
  deg_call_agreement[, agreement_preview_columns],
  row.names = FALSE,
  digits = 3
)


## optional exports ----

if (write_outputs) {
  write.table(
    deg_call_agreement_both,
    file.path(output_dir, "DEG_call_agreement_summary.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    deg_call_agreement_both,
    file.path(output_dir, "DEG_call_agreement_summary_both.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    deg_call_agreement_up,
    file.path(output_dir, "DEG_call_agreement_summary_up.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    deg_call_agreement_down,
    file.path(output_dir, "DEG_call_agreement_summary_down.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    deg_call_metric_spec,
    file.path(output_dir, "DEG_call_metric_spec.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}
