# scITD downstream benchmark for truth and inferred CTSE.
#
# Without command-line arguments, edit the values in section 2 and run the
# driver interactively. The positional command-line interface runs one method
# job at a time. Reusable logic lives in DALE_Eval/modules/scITD_helpers.R.

library(scITD)
library(Matrix)

# ============================================================
# 1. Locate the repository and load modules
# ============================================================


source('../DALE_Eval/modules/config_helpers.R')
repo_root <- find_repo_root()

source(file.path(
  repo_root,
  "DALE_Eval",
  "modules",
  "scITD_helpers.R"
))


# ============================================================
# 2. Editable benchmark selection
# ============================================================

cli_args <- commandArgs(trailingOnly = TRUE)
if (
  length(cli_args) > 0L &&
    !length(cli_args) %in% c(3L, 4L)
) {
  stop(
    "Usage: Rscript step5_scITD_benchmark.R ",
    "<dataset> <config_id> <method> [lib_norm:true|false]"
  )
}

parse_cli_boolean <- function(value) {
  normalized <- tolower(trimws(as.character(value)))
  if (!normalized %in% c("true", "false")) {
    stop("lib_norm must be true or false")
  }
  identical(normalized, "true")
}

datasets <- c(
  "PBMC_Perez2022",
  "ROSMAP_AD430_Mathys2023"
)
selected_config_ids <- NULL
selected_methods <- NULL
lib_norm_methods <- character()

scitd_config_id <- "tutorial_v1"

# TRUE creates the canonical truth decomposition before method runs.
# Set FALSE to reuse truth outputs already present under
# scITD_res/truth_sumcount/.
run_truth <- TRUE

# Existing managed scITD output files are protected unless this is TRUE.
overwrite_outputs <- FALSE

if (length(cli_args) > 0L) {
  datasets <- cli_args[[1]]
  selected_config_ids <- cli_args[[2]]
  selected_methods <- cli_args[[3]]
  if (length(cli_args) == 4L) {
    lib_norm <- parse_cli_boolean(cli_args[[4]])
    lib_norm_methods <- if (lib_norm) selected_methods else character()
  }
  run_truth <- FALSE
}


# ============================================================
# 3. Load the documented scITD benchmark config
# ============================================================

scitd_config <- read_scitd_config(
  scitd_config_id = scitd_config_id,
  repo_root = repo_root
)
method_preprocessing <- read_scitd_method_preprocessing(
  repo_root = repo_root
)

if (is.null(selected_config_ids)) {
  selected_config_ids <- scitd_config$deconv_config_ids
}
unsupported_config_ids <- setdiff(
  selected_config_ids,
  scitd_config$deconv_config_ids
)
if (length(unsupported_config_ids) > 0L) {
  stop("Unsupported scITD config: ", paste(unsupported_config_ids, collapse = ", "))
}
if (!is.null(selected_methods)) {
  unknown_methods <- setdiff(selected_methods, method_preprocessing$method)
  if (length(unknown_methods) > 0L) {
    stop("Unknown scITD method: ", paste(unknown_methods, collapse = ", "))
  }
}

message(
  "scITD config: ", scitd_config$scitd_config_id,
  "; deconvolution configs: ",
  paste(scitd_config$deconv_config_ids, collapse = ", "),
  "; method preprocessing mappings: ",
  nrow(method_preprocessing)
)


# ============================================================
# 4. Generate truth and inferred scITD results
# ============================================================

truth_datasets_completed <- 0L
method_jobs_completed <- 0L
method_jobs_skipped <- 0L

for (dataset in datasets) {
  message("\n################ ", dataset, " ################")

  truth_output_dir <- scitd_truth_result_dir(
    dataset,
    repo_root = repo_root
  )
  protected_truth_outputs <- scitd_managed_output_files(truth_output_dir)
  protected_truth_outputs <- protected_truth_outputs[
    file.exists(protected_truth_outputs)
  ]

  truth_spec <- tryCatch(
    {
      if (
        run_truth &&
          length(protected_truth_outputs) > 0 &&
          !overwrite_outputs
      ) {
        message(
          "Reusing protected truth scITD outputs: ",
          truth_output_dir
        )
        read_scitd_truth_spec(truth_output_dir)
      } else if (run_truth) {
        truth_result <- run_scitd_truth(
          dataset = dataset,
          scitd_config = scitd_config,
          repo_root = repo_root
        )
        write_scitd_outputs(
          result = truth_result,
          output_dir = truth_output_dir,
          overwrite = overwrite_outputs
        )
        specification <- scitd_truth_spec_from_result(truth_result)
        rm(truth_result)
        gc(verbose = FALSE)
        specification
      } else {
        message("Reusing saved truth scITD outputs: ", truth_output_dir)
        read_scitd_truth_spec(truth_output_dir)
      }
    },
    error = function(error) {
      message(
        "[SKIP] ", dataset, " truth: ",
        conditionMessage(error)
      )
      NULL
    }
  )

  if (is.null(truth_spec)) {
    if (length(cli_args) > 0L) {
      stop("Truth scITD result is unavailable for: ", dataset)
    }
    message("[SKIP] No inferred runs attempted because truth is unavailable")
    next
  }
  truth_datasets_completed <- truth_datasets_completed + 1L

  jobs <- discover_scitd_method_jobs(
    dataset = dataset,
    deconv_config_ids = selected_config_ids,
    methods = selected_methods,
    lib_norm_methods = lib_norm_methods,
    repo_root = repo_root
  )
  if (nrow(jobs) == 0) {
    if (length(cli_args) > 0L) {
      stop("Requested scITD method job is unavailable")
    }
    message("[SKIP] No available method jobs for: ", dataset)
    next
  }

  for (job_index in seq_len(nrow(jobs))) {
    job <- jobs[job_index, , drop = FALSE]
    protected_method_outputs <- scitd_managed_output_files(
      job$output_dir[[1]]
    )
    protected_method_outputs <- protected_method_outputs[
      file.exists(protected_method_outputs)
    ]
    if (
      length(protected_method_outputs) > 0 &&
        !overwrite_outputs
    ) {
      message(
        "[SKIP] ",
        job$dataset[[1]], " / ",
        job$config_id[[1]], " / ",
        job$method[[1]],
        ": protected scITD outputs already exist"
      )
      method_jobs_skipped <- method_jobs_skipped + 1L
      next
    }

    completed <- tryCatch(
      {
        inferred_result <- run_scitd_inferred(
          dataset = job$dataset[[1]],
          config_id = job$config_id[[1]],
          method = job$method[[1]],
          truth_spec = truth_spec,
          scitd_config = scitd_config,
          method_preprocessing = method_preprocessing,
          repo_root = repo_root,
          lib_norm = job$lib_norm[[1]],
          input_dir = job$input_dir[[1]]
        )
        write_scitd_outputs(
          result = inferred_result,
          output_dir = job$output_dir[[1]],
          overwrite = overwrite_outputs
        )
        rm(inferred_result)
        gc(verbose = FALSE)
        TRUE
      },
      error = function(error) {
        message(
          "[SKIP] ",
          job$dataset[[1]], " / ",
          job$config_id[[1]], " / ",
          job$method[[1]], ": ",
          conditionMessage(error)
        )
        FALSE
      }
    )

    if (completed) {
      method_jobs_completed <- method_jobs_completed + 1L
    } else {
      if (length(cli_args) > 0L) {
        stop("Requested scITD method job failed")
      }
      method_jobs_skipped <- method_jobs_skipped + 1L
    }
  }
}


# ============================================================
# 5. Completion summary
# ============================================================

message(
  "\nscITD workflow complete. Truth datasets available: ",
  truth_datasets_completed,
  "; method jobs completed: ",
  method_jobs_completed,
  "; method jobs skipped: ",
  method_jobs_skipped
)

# PBMC_Perez2022: config02 TCA does not retain usable B-cell patient-to-patient variation among genes passing variance QC


