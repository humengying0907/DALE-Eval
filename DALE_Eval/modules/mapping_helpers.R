read_indep_ref_cell_type_mapping <- function(paths, repo_root) {
  mapping_path <- file.path(repo_root, "Indep_scReference", "cell_type_mapping.txt")
  if (!file.exists(mapping_path)) {
    stop("Independent-reference cell type mapping file not found: ", mapping_path)
  }

  mapping <- read.delim(mapping_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("dataset", "indep_ref", "target_cell_type", "indep_ref_cell_type")
  missing_cols <- setdiff(required_cols, colnames(mapping))
  if (length(missing_cols) > 0) {
    stop("cell_type_mapping.txt missing columns: ", paste(missing_cols, collapse = ", "))
  }

  ref_name <- basename(paths$ref_dir)
  mapping[mapping$dataset == paths$dataset & mapping$indep_ref == ref_name, , drop = FALSE]
}

apply_indep_ref_cell_type_mapping <- function(z_array, paths, repo_root = find_repo_root()) {
  if (paths$refType != "indep") {
    return(z_array)
  }

  cell_types <- dimnames(z_array)[[3]]
  if (is.null(cell_types)) {
    stop("Cannot map CTSE array without cell-type dimnames on dimension 3")
  }

  mapping_path <- file.path(repo_root, "Indep_scReference", "cell_type_mapping.txt")
  if (!file.exists(mapping_path)) {
    message("Independent-reference cell type mapping file not found; exporting original reference cell type names: ", mapping_path)
    return(z_array)
  }

  mapping <- read_indep_ref_cell_type_mapping(paths, repo_root)
  if (nrow(mapping) == 0) {
    message("No cell type mapping found for dataset=", paths$dataset, " indep_ref=", basename(paths$ref_dir), "; exporting original reference cell type names")
    return(z_array)
  }

  duplicated_ref <- unique(mapping$indep_ref_cell_type[duplicated(mapping$indep_ref_cell_type)])
  if (length(duplicated_ref) > 0) {
    stop("Duplicate indep_ref_cell_type mappings found: ", paste(duplicated_ref, collapse = ", "))
  }

  map_to_target <- setNames(mapping$target_cell_type, mapping$indep_ref_cell_type)
  mapped_cell_types <- ifelse(cell_types %in% names(map_to_target), map_to_target[cell_types], cell_types)

  duplicated_out <- unique(mapped_cell_types[duplicated(mapped_cell_types)])
  if (length(duplicated_out) > 0) {
    stop("Cell type mapping creates duplicated output names: ", paste(duplicated_out, collapse = ", "))
  }

  dimnames(z_array)[[3]] <- unname(mapped_cell_types)
  z_array
}

map_truth_fraction_to_indep_ref <- function(frac, reference_cell_types, paths, repo_root) {
  mapping <- read_indep_ref_cell_type_mapping(paths, repo_root)
  if (nrow(mapping) == 0) {
    stop("No cell type mapping found for dataset=", paths$dataset, " indep_ref=", basename(paths$ref_dir))
  }

  missing_targets <- setdiff(colnames(frac), mapping$target_cell_type)
  if (length(missing_targets) > 0) {
    stop("Truth fraction cell types missing independent-reference mapping: ", paste(missing_targets, collapse = ", "))
  }

  mapping <- mapping[match(colnames(frac), mapping$target_cell_type), , drop = FALSE]
  duplicated_ref <- unique(mapping$indep_ref_cell_type[duplicated(mapping$indep_ref_cell_type)])
  if (length(duplicated_ref) > 0) {
    stop("Two truth fraction cell types map to the same independent-reference cell type: ", paste(duplicated_ref, collapse = ", "))
  }

  missing_reference <- setdiff(mapping$indep_ref_cell_type, reference_cell_types)
  if (length(missing_reference) > 0) {
    stop("Mapped independent-reference cell types missing from reference: ", paste(missing_reference, collapse = ", "))
  }

  colnames(frac) <- mapping$indep_ref_cell_type
  list(
    frac = frac[, mapping$indep_ref_cell_type, drop = FALSE],
    reference_cell_types = mapping$indep_ref_cell_type
  )
}

align_fraction_and_reference_cell_types <- function(frac, reference_cell_types, paths, repo_root) {
  frac_input <- paths$config$frac_input[[1]]
  if (setequal(reference_cell_types, colnames(frac))) {
    ordered_cell_types <- reference_cell_types
    return(list(
      frac = frac[, ordered_cell_types, drop = FALSE],
      reference_cell_types = ordered_cell_types
    ))
  }

  if (paths$refType == "indep" && frac_input %in% c("truth_cellfrac", "truth_transcriptfrac")) {
    return(map_truth_fraction_to_indep_ref(frac, reference_cell_types, paths, repo_root))
  }

  missing_in_frac <- setdiff(reference_cell_types, colnames(frac))
  missing_in_reference <- setdiff(colnames(frac), reference_cell_types)
  stop(
    "Reference and fraction cell types do not match. ",
    "Missing in fraction: ", paste(missing_in_frac, collapse = ", "),
    "; missing in reference: ", paste(missing_in_reference, collapse = ", ")
  )
}
