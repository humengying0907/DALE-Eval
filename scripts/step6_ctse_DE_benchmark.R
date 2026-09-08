# Notebook-style runner for config02 cell type-specific DE benchmarking.
# Run from scripts/ with: Rscript step6_ctse_DE_benchmark.R

# ============================================================
# 1. Settings
# ============================================================

datasets <- c(
  "PBMC_Perez2022",
  "ROSMAP_AD430_Mathys2023"
)

methods <- c(
  "InstaPrism",
  "BLUE",
  "scTAPE",
  "CIBERSORTx",
  "ENIGMAL2",
  "ENIGMAtrace",
  "bMIND",
  "TCA",
  "Unico",
  "EPICunmix"
)

run_truth <- T
run_inferred <- T
run_bMIND_direct <- F
run_ENIGMAL2_direct <- F
run_ENIGMAtrace_direct <- F
run_bulk_limma_baseline <- TRUE
resume_bMIND_chunks <- F

overwrite_outputs <- T
n_core <- 12

de_covariate_cols <- list(
  PBMC_Perez2022 = c("Sex", "pop_cov"),
  ROSMAP_AD430_Mathys2023 = c("msex")
)

# ============================================================
# 2. Source helpers
# ============================================================

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
step6_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
} else {
  normalizePath(getwd())
}

source(file.path(step6_dir, "..", "DALE_Eval", "modules", "ctse_de_helpers.R"))
repo_root <- ctse_de_find_repo_root(step6_dir)
source(file.path(repo_root, "DALE_Eval", "modules", "config_helpers.R"))
source(file.path(repo_root, "DALE_Eval", "modules", "mapping_helpers.R"))
source(file.path(repo_root, "DALE_Eval", "modules", "direct_ctse_de_helpers.R"))

step6_attach_de_covariates <- function(context) {
  columns <- de_covariate_cols[[context$dataset]]
  if (is.null(columns)) {
    stop("No step6 DE covariates configured for dataset: ", context$dataset)
  }
  context$covariates <- ctse_de_read_sample_covariates(context, columns)
  context
}

truth_summary_rows <- stats::setNames(vector("list", length(datasets)), datasets)
inferred_summary_rows <- stats::setNames(vector("list", length(datasets)), datasets)
direct_summary_rows <- stats::setNames(vector("list", length(datasets)), datasets)

# ============================================================
# 3. Truth CTSE DE: TMM -> voom -> limma, one cell type at a time
# ============================================================

if (run_truth) {
  for (dataset in datasets) {
    message("\n================ TRUTH: ", dataset, " ================")
    config <- read_ctse_de_config(dataset, repo_root)
    context <- ctse_de_read_sample_context(dataset, config, repo_root)
    context <- step6_attach_de_covariates(context)
    ctse_meta <- ctse_de_read_truth_meta(context)
    truth_files <- ctse_de_truth_files(dataset, config, repo_root)

    for (cell_type in names(truth_files)) {
      message("[truth] ", dataset, " / ", cell_type)
      output_path <- file.path(
        context$paths$truth_output_dir,
        ctse_de_safe_cell_type_filename(cell_type)
      )
      job <- tryCatch(
        {
          fitted <- ctse_de_run_truth_cell_type(
            dataset = dataset,
            cell_type = cell_type,
            config = config,
            context = context,
            ctse_meta = ctse_meta
          )
          ctse_de_write_limma_result(
            fitted$result,
            output_path,
            overwrite = overwrite_outputs
          )
          fitted$summary$output_file <- output_path
          fitted$summary
        },
        error = function(error) {
          message("[skip truth] ", conditionMessage(error))
          ctse_de_summary_row(
            dataset = dataset,
            source = "truth",
            method = "truth",
            cell_type = cell_type,
            status = "skipped",
            preprocessing = "TMM + voom",
            effect_scale = "voom log2-expression difference",
            problem = conditionMessage(error),
            output_file = output_path
          )
        }
      )
      truth_summary_rows[[dataset]][[cell_type]] <- job
    }
  }
}

# ============================================================
# 4. Inferred CTSE DE: explicit method and cell-type loops
# ============================================================

if (run_inferred) {
  for (dataset in datasets) {
    message("\n================ INFERRED: ", dataset, " ================")
    config <- read_ctse_de_config(dataset, repo_root)
    context <- ctse_de_read_sample_context(dataset, config, repo_root)
    context <- step6_attach_de_covariates(context)

    for (method in methods) {
      message("\n[inferred method] ", dataset, " / ", method)
      method_dir <- file.path(context$paths$deconv_dir, method)
      specification_problem <- ""
      specification <- tryCatch(
        ctse_de_detect_method_preprocessing(method, method_dir),
        error = function(error) {
          specification_problem <<- conditionMessage(error)
          message("[skip inferred method] ", specification_problem)
          NULL
        }
      )
      if (is.null(specification)) {
        key <- paste(method, "method", sep = "\r")
        inferred_summary_rows[[dataset]][[key]] <- ctse_de_summary_row(
          dataset = dataset,
          source = "inferred",
          method = method,
          status = "skipped",
          problem = specification_problem
        )
        next
      }

      message("  contains_negative: ", specification$contains_negative)
      message("  max_ctse_value: ", signif(specification$max_ctse_value, 6))
      message("  preprocessing: ", specification$preprocessing)

      for (cell_type in names(specification$files)) {
        message("[inferred cell type] ", method, " / ", cell_type)
        output_path <- file.path(
          context$paths$inferred_output_dir,
          method,
          ctse_de_safe_cell_type_filename(cell_type)
        )
        job <- tryCatch(
          {
            fitted <- ctse_de_run_inferred_cell_type(
              dataset = dataset,
              method = method,
              cell_type = cell_type,
              specification = specification,
              context = context
            )
            ctse_de_write_limma_result(
              fitted$result,
              output_path,
              overwrite = overwrite_outputs
            )
            fitted$summary$output_file <- output_path
            fitted$summary
          },
          error = function(error) {
            message("[skip inferred cell type] ", conditionMessage(error))
            ctse_de_summary_row(
              dataset = dataset,
              source = "inferred",
              method = method,
              cell_type = cell_type,
              status = "skipped",
              contains_negative = specification$contains_negative,
              max_ctse_value = specification$max_ctse_value,
              preprocessing = specification$preprocessing,
              effect_scale = specification$effect_scale,
              problem = conditionMessage(error),
              output_file = output_path
            )
          }
        )
        key <- paste(method, cell_type, sep = "\r")
        inferred_summary_rows[[dataset]][[key]] <- job
      }
    }
  }
}

# ============================================================
# 5. Native bMIND DE, in resumable gene chunks
# Output columns: gene, effect_size, p_value, q_value (no t column)
# ============================================================

if (run_bMIND_direct) {
  for (dataset in datasets) {
    message("\n================ DIRECT bMIND: ", dataset, " ================")
    config <- read_ctse_de_config(dataset, repo_root)
    context <- ctse_de_read_sample_context(dataset, config, repo_root)
    context <- step6_attach_de_covariates(context)

    input_problem <- ""
    prepared <- tryCatch(
      direct_ctse_de_prepare_inputs(dataset, config, context),
      error = function(error) {
        input_problem <<- conditionMessage(error)
        message("[skip direct bMIND input] ", input_problem)
        NULL
      }
    )
    if (is.null(prepared)) {
      direct_summary_rows[[dataset]][["bMIND\rmethod"]] <- ctse_de_summary_row(
        dataset = dataset,
        source = "direct",
        method = "bMIND",
        status = "skipped",
        problem = input_problem
      )
      next
    }

    plan_problem <- ""
    bmind_plan <- tryCatch(
      direct_ctse_de_prepare_bmind_plan(
        prepared,
        chunk_size = config$bmind_de_chunk_size
      ),
      error = function(error) {
        plan_problem <<- conditionMessage(error)
        message("[skip direct bMIND plan] ", plan_problem)
        NULL
      }
    )
    if (is.null(bmind_plan)) {
      direct_summary_rows[[dataset]][["bMIND\rmethod"]] <- ctse_de_summary_row(
        dataset = dataset,
        source = "direct",
        method = "bMIND",
        status = "skipped",
        n_input_samples = length(prepared$samples),
        n_input_genes = length(prepared$genes),
        problem = plan_problem
      )
      next
    }

    runtime_entry <- direct_ctse_de_start_runtime(
      paths = context$paths,
      dataset = dataset,
      method = "bMIND",
      config_id = config$config_id,
      n_core = n_core,
      n_input_samples = length(prepared$samples),
      n_input_genes = length(prepared$genes),
      chunk_size = bmind_plan$chunk_size,
      n_chunks = bmind_plan$n_chunks
    )

    native_problem <- ""
    native_error <- NULL
    native_job <- tryCatch(
      run_bmind_ctse_de(
        prepared = prepared,
        plan = bmind_plan,
        config = config,
        paths = context$paths,
        n_core = n_core,
        resume = resume_bMIND_chunks,
        run_id = runtime_entry$run_id
      ),
      interrupt = function(condition) {
        direct_ctse_de_finish_runtime(
          runtime_entry,
          status = "failed",
          problem = conditionMessage(condition)
        )
        stop(condition)
      },
      error = function(error) {
        native_error <<- error
        native_problem <<- conditionMessage(error)
        message("[skip direct bMIND] ", native_problem)
        NULL
      }
    )
    if (is.null(native_job)) {
      n_tested <- if (is.null(native_error$n_tested_genes)) {
        NA_integer_
      } else {
        native_error$n_tested_genes
      }
      n_completed <- if (is.null(native_error$n_completed_chunks)) {
        0L
      } else {
        native_error$n_completed_chunks
      }
      n_failed <- if (is.null(native_error$n_failed_chunks)) {
        0L
      } else {
        native_error$n_failed_chunks
      }
      direct_ctse_de_finish_runtime(
        runtime_entry,
        status = "failed",
        n_tested_genes = n_tested,
        n_completed_chunks = n_completed,
        n_failed_chunks = n_failed,
        problem = native_problem
      )
      direct_summary_rows[[dataset]][["bMIND\rmethod"]] <- ctse_de_summary_row(
        dataset = dataset,
        source = "direct",
        method = "bMIND",
        status = "skipped",
        n_input_samples = length(prepared$samples),
        n_control = prepared$n_control,
        n_case = prepared$n_case,
        n_input_genes = length(prepared$genes),
        n_tested_genes = n_tested,
        problem = native_problem
      )
      next
    }

    if (native_job$resumed) {
      message("[direct bMIND] resumed compatible temporary chunks")
    }
    output_problems <- character()
    for (cell_type in names(native_job$results)) {
      message("[direct bMIND cell type] ", cell_type)
      output_path <- file.path(
        context$paths$direct_output_dir,
        "bMIND",
        ctse_de_safe_cell_type_filename(cell_type)
      )
      summary_row <- native_job$summaries[[cell_type]]
      write_error <- tryCatch(
        {
          ctse_de_write_native_result(
            native_job$results[[cell_type]],
            output_path,
            overwrite = overwrite_outputs
          )
          NULL
        },
        error = function(error) error
      )
      if (is.null(write_error)) {
        summary_row$output_file <- output_path
      } else {
        problem <- conditionMessage(write_error)
        message("[skip direct bMIND output] ", problem)
        output_problems <- c(output_problems, problem)
        summary_row$status <- "skipped"
        summary_row$problem <- problem
        summary_row$output_file <- output_path
      }
      key <- paste("bMIND", cell_type, sep = "\r")
      direct_summary_rows[[dataset]][[key]] <- summary_row
    }

    if (length(output_problems) > 0L) {
      direct_ctse_de_finish_runtime(
        runtime_entry,
        status = "failed",
        n_tested_genes = native_job$n_tested_genes,
        n_completed_chunks = native_job$n_completed_chunks,
        n_failed_chunks = native_job$n_failed_chunks,
        problem = paste(unique(output_problems), collapse = " | ")
      )
      message("[direct bMIND] temporary chunks retained after output failure")
      next
    }

    cleanup_error <- tryCatch(
      {
        direct_ctse_de_cleanup_bmind_chunks(native_job$chunk_dir)
        NULL
      },
      error = function(error) error
    )
    if (!is.null(cleanup_error)) {
      direct_ctse_de_finish_runtime(
        runtime_entry,
        status = "failed",
        n_tested_genes = native_job$n_tested_genes,
        n_completed_chunks = native_job$n_completed_chunks,
        n_failed_chunks = native_job$n_failed_chunks,
        problem = conditionMessage(cleanup_error)
      )
      message("[direct bMIND cleanup failed] ", conditionMessage(cleanup_error))
      next
    }

    direct_ctse_de_finish_runtime(
      runtime_entry,
      status = "completed",
      n_tested_genes = native_job$n_tested_genes,
      n_completed_chunks = native_job$n_completed_chunks,
      n_failed_chunks = native_job$n_failed_chunks
    )
  }
}

# ============================================================
# 6. Native ENIGMA L2 DE
# Output columns: gene, effect_size, p_value, q_value (no t column)
# ============================================================

if (run_ENIGMAL2_direct) {
  for (dataset in datasets) {
    message("\n================ DIRECT ENIGMA L2: ", dataset, " ================")
    config <- read_ctse_de_config(dataset, repo_root)
    context <- ctse_de_read_sample_context(dataset, config, repo_root)
    context <- step6_attach_de_covariates(context)

    input_problem <- ""
    prepared_common <- tryCatch(
      direct_ctse_de_prepare_inputs(dataset, config, context),
      error = function(error) {
        input_problem <<- conditionMessage(error)
        NULL
      }
    )
    prepared <- if (is.null(prepared_common)) {
      NULL
    } else {
      tryCatch(
        direct_ctse_de_prepare_enigma_inputs(
          prepared_common,
          config,
          repo_root = repo_root
        ),
        error = function(error) {
          input_problem <<- conditionMessage(error)
          NULL
        }
      )
    }
    if (is.null(prepared)) {
      message("[skip direct ENIGMAL2 input] ", input_problem)
      direct_summary_rows[[dataset]][["ENIGMAL2\rmethod"]] <-
        ctse_de_summary_row(
          dataset = dataset,
          source = "direct",
          method = "ENIGMAL2",
          status = "skipped",
          problem = input_problem
        )
      next
    }

    runtime_entry <- direct_ctse_de_start_runtime(
      paths = context$paths,
      dataset = dataset,
      method = "ENIGMAL2",
      config_id = config$config_id,
      n_core = n_core,
      n_input_samples = length(prepared_common$samples),
      n_input_genes = length(prepared_common$genes)
    )
    native_problem <- ""
    native_job <- tryCatch(
      run_enigma_ctse_de(prepared = prepared, mode = "L2"),
      interrupt = function(condition) {
        direct_ctse_de_finish_runtime(
          runtime_entry,
          status = "failed",
          problem = conditionMessage(condition)
        )
        stop(condition)
      },
      error = function(error) {
        native_problem <<- conditionMessage(error)
        message("[skip direct ENIGMAL2] ", native_problem)
        NULL
      }
    )
    if (is.null(native_job)) {
      direct_ctse_de_finish_runtime(
        runtime_entry,
        status = "failed",
        problem = native_problem
      )
      direct_summary_rows[[dataset]][["ENIGMAL2\rmethod"]] <-
        ctse_de_summary_row(
          dataset = dataset,
          source = "direct",
          method = "ENIGMAL2",
          status = "skipped",
          n_input_samples = length(prepared_common$samples),
          n_control = prepared_common$n_control,
          n_case = prepared_common$n_case,
          n_input_genes = length(prepared_common$genes),
          problem = native_problem
        )
      next
    }

    output_problems <- character()
    for (cell_type in names(native_job$results)) {
      message("[direct ENIGMAL2 cell type] ", cell_type)
      output_path <- file.path(
        context$paths$direct_output_dir,
        "ENIGMAL2",
        ctse_de_safe_cell_type_filename(cell_type)
      )
      summary_row <- native_job$summaries[[cell_type]]
      write_error <- tryCatch(
        {
          ctse_de_write_native_result(
            native_job$results[[cell_type]],
            output_path,
            overwrite = overwrite_outputs
          )
          NULL
        },
        error = function(error) error
      )
      if (is.null(write_error)) {
        summary_row$output_file <- output_path
      } else {
        problem <- conditionMessage(write_error)
        message("[skip direct ENIGMAL2 output] ", problem)
        output_problems <- c(output_problems, problem)
        summary_row$status <- "skipped"
        summary_row$problem <- problem
        summary_row$output_file <- output_path
      }
      key <- paste("ENIGMAL2", cell_type, sep = "\r")
      direct_summary_rows[[dataset]][[key]] <- summary_row
    }

    runtime_status <- if (length(output_problems) == 0L) "completed" else "failed"
    runtime_problem <- if (length(output_problems) == 0L) {
      ""
    } else {
      paste(unique(output_problems), collapse = " | ")
    }
    direct_ctse_de_finish_runtime(
      runtime_entry,
      status = runtime_status,
      n_tested_genes = native_job$n_tested_genes,
      n_completed_chunks = NA_integer_,
      n_failed_chunks = NA_integer_,
      problem = runtime_problem
    )
  }
}

# ============================================================
# 7. Native ENIGMA trace DE
# Output columns: gene, effect_size, p_value, q_value (no t column)
# ============================================================

if (run_ENIGMAtrace_direct) {
  for (dataset in datasets) {
    message("\n================ DIRECT ENIGMA trace: ", dataset, " ================")
    config <- read_ctse_de_config(dataset, repo_root)
    context <- ctse_de_read_sample_context(dataset, config, repo_root)
    context <- step6_attach_de_covariates(context)

    input_problem <- ""
    prepared_common <- tryCatch(
      direct_ctse_de_prepare_inputs(dataset, config, context),
      error = function(error) {
        input_problem <<- conditionMessage(error)
        NULL
      }
    )
    prepared <- if (is.null(prepared_common)) {
      NULL
    } else {
      tryCatch(
        direct_ctse_de_prepare_enigma_inputs(
          prepared_common,
          config,
          repo_root = repo_root
        ),
        error = function(error) {
          input_problem <<- conditionMessage(error)
          NULL
        }
      )
    }
    if (is.null(prepared)) {
      message("[skip direct ENIGMAtrace input] ", input_problem)
      direct_summary_rows[[dataset]][["ENIGMAtrace\rmethod"]] <-
        ctse_de_summary_row(
          dataset = dataset,
          source = "direct",
          method = "ENIGMAtrace",
          status = "skipped",
          problem = input_problem
        )
      next
    }

    runtime_entry <- direct_ctse_de_start_runtime(
      paths = context$paths,
      dataset = dataset,
      method = "ENIGMAtrace",
      config_id = config$config_id,
      n_core = n_core,
      n_input_samples = length(prepared_common$samples),
      n_input_genes = length(prepared_common$genes)
    )
    native_problem <- ""
    native_job <- tryCatch(
      run_enigma_ctse_de(prepared = prepared, mode = "trace"),
      interrupt = function(condition) {
        direct_ctse_de_finish_runtime(
          runtime_entry,
          status = "failed",
          problem = conditionMessage(condition)
        )
        stop(condition)
      },
      error = function(error) {
        native_problem <<- conditionMessage(error)
        message("[skip direct ENIGMAtrace] ", native_problem)
        NULL
      }
    )
    if (is.null(native_job)) {
      direct_ctse_de_finish_runtime(
        runtime_entry,
        status = "failed",
        problem = native_problem
      )
      direct_summary_rows[[dataset]][["ENIGMAtrace\rmethod"]] <-
        ctse_de_summary_row(
          dataset = dataset,
          source = "direct",
          method = "ENIGMAtrace",
          status = "skipped",
          n_input_samples = length(prepared_common$samples),
          n_control = prepared_common$n_control,
          n_case = prepared_common$n_case,
          n_input_genes = length(prepared_common$genes),
          problem = native_problem
        )
      next
    }

    output_problems <- character()
    for (cell_type in names(native_job$results)) {
      message("[direct ENIGMAtrace cell type] ", cell_type)
      output_path <- file.path(
        context$paths$direct_output_dir,
        "ENIGMAtrace",
        ctse_de_safe_cell_type_filename(cell_type)
      )
      summary_row <- native_job$summaries[[cell_type]]
      write_error <- tryCatch(
        {
          ctse_de_write_native_result(
            native_job$results[[cell_type]],
            output_path,
            overwrite = overwrite_outputs
          )
          NULL
        },
        error = function(error) error
      )
      if (is.null(write_error)) {
        summary_row$output_file <- output_path
      } else {
        problem <- conditionMessage(write_error)
        message("[skip direct ENIGMAtrace output] ", problem)
        output_problems <- c(output_problems, problem)
        summary_row$status <- "skipped"
        summary_row$problem <- problem
        summary_row$output_file <- output_path
      }
      key <- paste("ENIGMAtrace", cell_type, sep = "\r")
      direct_summary_rows[[dataset]][[key]] <- summary_row
    }

    runtime_status <- if (length(output_problems) == 0L) "completed" else "failed"
    runtime_problem <- if (length(output_problems) == 0L) {
      ""
    } else {
      paste(unique(output_problems), collapse = " | ")
    }
    direct_ctse_de_finish_runtime(
      runtime_entry,
      status = runtime_status,
      n_tested_genes = native_job$n_tested_genes,
      n_completed_chunks = NA_integer_,
      n_failed_chunks = NA_integer_,
      problem = runtime_problem
    )
  }
}

# ============================================================
# 8. Non-cell-type-specific bulk limma baselines
#
# Baseline 1:
# raw sum counts -> filterByExpr -> TMM -> voom ->
# configured DE covariates + disease
#
# Baseline 2:
# raw sum counts -> filterByExpr -> TMM -> voom ->
# K-1 truth transcript fractions + configured DE covariates + disease
#
# Neither baseline receives a direct CTS-DE runtime row.
# ============================================================

if (run_bulk_limma_baseline) {
  for (dataset in datasets) {
    config <- read_ctse_de_config(dataset, repo_root)
    context <- ctse_de_read_sample_context(dataset, config, repo_root)
    context <- step6_attach_de_covariates(context)
    
    # --------------------------------------------------------
    # 8a. Unadjusted bulk limma
    # --------------------------------------------------------
    
    message(
      "\n================ BULK LIMMA BASELINE: ",
      dataset,
      " ================"
    )
    
    bulk_output_path <- file.path(
      context$paths$baseline_output_dir,
      "bulk_limma",
      "bulk.txt.gz"
    )
    
    bulk_job <- tryCatch(
      {
        prepared_bulk <- direct_bulk_de_prepare_inputs(
          dataset,
          context
        )
        
        fitted_bulk <- run_bulk_limma_de(prepared_bulk)
        
        ctse_de_write_limma_result(
          fitted_bulk$result,
          bulk_output_path,
          overwrite = overwrite_outputs
        )
        
        fitted_bulk
      },
      error = function(error) {
        message(
          "[skip bulk limma baseline] ",
          conditionMessage(error)
        )
        NULL
      }
    )
    
    if (!is.null(bulk_job)) {
      message(
        "[bulk limma baseline] saved ",
        bulk_output_path,
        " (",
        bulk_job$n_tested_genes,
        " tested genes; control=",
        bulk_job$n_control,
        ", case=",
        bulk_job$n_case,
        ")"
      )
    }
    
    # --------------------------------------------------------
    # 8b. Truth-fraction-adjusted bulk limma
    # --------------------------------------------------------
    
    message(
      "\n================ BULK TRUTH-FRACTION-ADJUSTED LIMMA: ",
      dataset,
      " ================"
    )
    
    adjusted_output_path <- file.path(
      context$paths$baseline_output_dir,
      "bulk_truthfrac_adjusted_limma",
      "bulk.txt.gz"
    )
    
    adjusted_job <- tryCatch(
      {
        prepared_adjusted <- direct_bulk_truthfrac_de_prepare_inputs(
          dataset = dataset,
          config = config,
          context = context
        )
        
        fitted_adjusted <- run_bulk_truthfrac_adjusted_limma(
          prepared_adjusted
        )
        
        ctse_de_write_limma_result(
          fitted_adjusted$result,
          adjusted_output_path,
          overwrite = overwrite_outputs
        )
        
        fitted_adjusted
      },
      error = function(error) {
        message(
          "[skip bulk truth-fraction-adjusted limma] ",
          conditionMessage(error)
        )
        NULL
      }
    )
    
    if (!is.null(adjusted_job)) {
      message(
        "[bulk truth-fraction-adjusted limma] saved ",
        adjusted_output_path,
        " (",
        adjusted_job$n_tested_genes,
        " tested genes; control=",
        adjusted_job$n_control,
        ", case=",
        adjusted_job$n_case,
        "; reference fraction=",
        adjusted_job$reference_fraction,
        ")"
      )
    }
  }
}

# ============================================================
# 9. Persist per-run analysis summaries
# ============================================================

summary_run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

for (dataset in datasets) {
  summary_rows <- c(
    unname(truth_summary_rows[[dataset]]),
    unname(inferred_summary_rows[[dataset]]),
    unname(direct_summary_rows[[dataset]])
  )
  
  if (length(summary_rows) == 0L) {
    next
  }
  
  summary_table <- ctse_de_bind_summaries(summary_rows)
  
  summary_path <- file.path(
    repo_root,
    "Benchmarking_obj",
    dataset,
    "DE_res",
    "run_summaries",
    paste0("step6_", summary_run_id, ".txt")
  )
  
  ctse_de_write_table(
    summary_table,
    summary_path,
    overwrite = FALSE
  )
  
  message(
    "[analysis summary] saved ",
    summary_path,
    " (",
    nrow(summary_table),
    " rows)"
  )
}

message("\nstep6 complete")
