import os
import sys
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

import glob
import random
import numpy as np
import pandas as pd
from natsort import natsorted
from math import floor


def aggregate_by_group(Expr, 
                       identifier_labels,  # A list or 1D NumPy array containing labels used for grouping. 
                       gene_names=None, 
                       min_nCells_per_pseudobulk=10):
    
    if Expr.shape[1] != len(identifier_labels):
        raise ValueError("Number of columns in Expr must equal the length of identifier_labels")
        
    if isinstance(identifier_labels, list):
        identifier_labels = np.array(identifier_labels)  
    elif isinstance(identifier_labels, np.ndarray) and identifier_labels.ndim != 1:
        raise ValueError("identifier_labels must be a 1D array if it's a NumPy array.")

    # Get unique labels and filter by min_nCells_per_pseudobulk
    labels, counts = np.unique(identifier_labels, return_counts=True)
    valid_labels = [label for label, count in zip(labels, counts) if count >= min_nCells_per_pseudobulk]

    if len(valid_labels) == 0:
        print(f"No identifiers with more than {min_nCells_per_pseudobulk} cells to aggregate. Returning an empty DataFrame.")
        return pd.DataFrame(), []

    # Group indices by valid labels
    groups = {}
    nCells_per_pseudobulk = []
    
    groups = {label: np.where(identifier_labels == label)[0] for label in valid_labels}
    nCells_per_pseudobulk.extend([len(groups[label]) for label in valid_labels])

    # Aggregate expression values for each group
    C = np.column_stack([np.mean(Expr[:, indices], axis=1) for indices in groups.values()])

    # Create DataFrame for aggregated results
    df_C = pd.DataFrame(C, columns=valid_labels)
    if gene_names is not None:
        df_C.index = gene_names

    return df_C, nCells_per_pseudobulk


def adata_qc(adata, min_nCells_per_pseudobulk=10, min_cells=100):    
    error_messages = []
    
    # Check required columns in adata.obs
    required_columns = ["cell_type", "sample"]
    missing_columns = [col for col in required_columns if col not in adata.obs.columns]
    if missing_columns:
        error_messages.append(f"QC failed: Missing columns in adata.obs: {', '.join(missing_columns)}.")
    
    # Check for unlabeled cells (NaN in 'cell_type')
    if "cell_type" in adata.obs.columns:
        n_unlabeled = adata.obs['cell_type'].isna().sum()
        if n_unlabeled > 0:
            error_messages.append(f"QC failed: {n_unlabeled} cells are unlabeled (contain NaN in 'cell_type' column).")
        
        # Check for blank cell_type names
        if adata.obs["cell_type"].str.strip().eq("").any():
            error_messages.append("QC failed: Some cell_type names are blank.")
    
    # Check if `adata.X` is log-transformed
    if np.max(adata.X) < 50:
        error_messages.append("QC failed: adata.X appears to be log-transformed (maximum value < 50).")
    
    # Check for genes expressed in fewer than `min_cells` cells
    genes_expressed_per_cell = (adata.X > 0).sum(axis=0).A1  # Nonzero values per gene
    if np.sum(genes_expressed_per_cell < min_cells) > 0:
        print(f"Warning: Some genes are expressed in fewer than {min_cells} cells.")

    colname_of_cellType = "cell_type"
    colname_of_sample = "sample"

    for cell_type, filtered_obs in adata.obs.groupby(colname_of_cellType):
        group_sizes = filtered_obs.groupby(colname_of_sample).size()

        # At least 30% of samples have >= min_nCells_per_pseudobulk
        if (group_sizes >= min_nCells_per_pseudobulk).sum() / len(group_sizes) < 0.3:
            error_messages.append(f"QC failed: Cell type '{cell_type}' does not meet the sample threshold criteria (>= 30% of samples with >= {min_nCells_per_pseudobulk} cells).")

    if error_messages:
        for message in error_messages:
            print(message)
    else:
        print(
            f"QC passed: All checks satisfied\n"
            f"Final adata shape: {adata.shape}\n"
            f"Final samples present in adata: {adata.obs['sample'].nunique()}\n"
        )


def pseudobulk_nCells_by_ct(obs,
                            colname_of_cellType = 'cell_type',
                            colname_of_cellType_heterogeneity = 'sample',
                            min_nCells_per_pseudobulk = 10):
    
    results = []
    for cell_type, filtered_obs in obs.groupby(colname_of_cellType):
        group_sizes = filtered_obs.groupby(colname_of_cellType_heterogeneity).size()
        min_ncells = group_sizes.min()
        max_ncells = group_sizes.max()
        results.append({'cell_type': cell_type, 
                        'min_nCells': min_ncells, 
                        'max_nCells': max_ncells,
                        'nSamples' : (group_sizes >= min_nCells_per_pseudobulk).sum()})

    result_df = pd.DataFrame(results)

    return result_df


def train_test_split(obs,
                     training_ratio=0.2,
                     colname_of_cellType='cell_type',
                     colname_of_sample='sample',
                     min_nCells_per_ct_reference=100,
                     min_nSamples_per_cell_type=10,
                     max_iterations=20):

    unique_sampleIDs = obs[colname_of_sample].unique()
    iterations = 0  

    # Initialize training_obs
    while iterations < max_iterations:
        iterations += 1

        train_sampleIDs = np.random.choice(unique_sampleIDs, 
                                           size=floor(len(unique_sampleIDs) * training_ratio), 
                                           replace=False)
        
        training_obs = obs[obs[colname_of_sample].isin(train_sampleIDs)]
        test_obs = obs[~obs[colname_of_sample].isin(train_sampleIDs)]

        # Count cells for each cell type in the training set
        n = training_obs[colname_of_cellType].value_counts()

        # Check the number of samples for each cell type in the test set
        nSamples = test_obs.groupby(colname_of_cellType)[colname_of_sample].nunique()

        if all(n >= min_nCells_per_ct_reference) and all(nSamples >= min_nSamples_per_cell_type):
            test_sampleIDs = [sample for sample in unique_sampleIDs if sample not in train_sampleIDs]
            return train_sampleIDs.tolist(), test_sampleIDs

    # If the maximum number of iterations is reached and no valid split found
    raise ValueError(
        "Unable to satisfy conditions for training and test splits. "
        "Please check if the cell types and sample sizes meet the requirements."
    )

