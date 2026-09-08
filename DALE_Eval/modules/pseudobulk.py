from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc


PSEUDOBULK_TYPES = {"count", "cpm"}


def _matrix_values(X, max_values=100000):
    values = X.data if sp.issparse(X) else np.asarray(X).ravel()
    if values.size > max_values:
        step = max(1, values.size // max_values)
        values = values[::step][:max_values]
    return values


def _check_matrix_basic(X):
    values = _matrix_values(X)
    if values.size == 0:
        raise ValueError("adata.X is empty or all zero")
    if not np.isfinite(values).all():
        raise ValueError("adata.X contains non-finite values")
    if np.nanmin(values) < 0:
        raise ValueError("adata.X contains negative values")

    lib_sizes = np.asarray(X.sum(axis=1)).ravel()
    if lib_sizes.size == 0 or not np.isfinite(lib_sizes).any() or np.all(lib_sizes == 0):
        raise ValueError("adata.X has empty/non-finite/all-zero cell libraries")
    return values, lib_sizes


def _detect_expression_state(X):
    values, lib_sizes = _check_matrix_basic(X)
    max_value = float(np.nanmax(values))
    if max_value < 50:
        return "log"

    integer_like = np.allclose(values, np.round(values), rtol=0, atol=1e-6)
    lib_mean = float(np.nanmean(lib_sizes))
    lib_cv = float(np.nanstd(lib_sizes) / (lib_mean + 1e-12))
    cpm_like = lib_cv < 0.05 and lib_mean > 1e5
    if integer_like and not cpm_like:
        return "count"
    return "cpm"


def _require_columns(obs, columns):
    missing = [col for col in columns if col not in obs.columns]
    if missing:
        raise KeyError(f"Missing required metadata columns: {missing}")
    for col in columns:
        if obs[col].isna().any():
            raise ValueError(f"{col!r} must not contain NA values")


def compute_cellfrac(scMeta, sample_col="sample", ct_col="cell_type"):
    """
    Compute sample-by-cell-type cell count fractions.
    """
    _require_columns(scMeta, [sample_col, ct_col])
    df = (
        scMeta
        .groupby([sample_col, ct_col], observed=True)
        .size()
        .reset_index(name="n")
    )
    df["frac"] = df["n"] / df.groupby(sample_col, observed=True)["n"].transform("sum")
    frac_matrix = df.pivot(index=sample_col, columns=ct_col, values="frac").fillna(0)
    frac_matrix.index.name = None
    frac_matrix.columns.name = None
    return frac_matrix


def compute_transcriptfrac(adata, sample_col="sample", ct_col="cell_type"):
    """
    Compute sample-by-cell-type raw UMI/transcript fractions from adata.X.
    """
    X = adata.X.tocsr() if sp.issparse(adata.X) else np.asarray(adata.X)
    if _detect_expression_state(X) != "count":
        raise ValueError("truth_transcriptfrac requires raw count-like adata.X")

    _require_columns(adata.obs, [sample_col, ct_col])
    cell_totals = np.asarray(X.sum(axis=1)).ravel()
    df = pd.DataFrame({
        sample_col: adata.obs[sample_col].astype(str).to_numpy(),
        ct_col: adata.obs[ct_col].astype(str).to_numpy(),
        "transcripts": cell_totals,
    })
    transcript_sum = (
        df.groupby([sample_col, ct_col], observed=True)["transcripts"]
        .sum()
        .unstack(fill_value=0)
    )
    frac_matrix = transcript_sum.div(transcript_sum.sum(axis=1), axis=0).fillna(0)
    frac_matrix.index.name = None
    frac_matrix.columns.name = None
    return frac_matrix


def build_ct_pseudobulk_profiles(
    adata,
    ct_group_keys,
    pseudobulk_type="count",
    cpm_target_sum=1e6,
    verbose=True,
):
    """
    Build one pseudobulk assumption from adata.X.

    Parameters
    ----------
    adata : AnnData
        Single-cell object with expression in adata.X.
    ct_group_keys : str or list[str]
        One key returns a gene x group DataFrame. Two keys return a dictionary
        keyed by the first grouping variable, with gene x second-key DataFrames.
        For CTSE benchmark truth, use [ct_col, sample_col].
    pseudobulk_type : {'count', 'cpm'}
        'count' returns raw count sums and requires count-like adata.X.
        'cpm' returns mean per-cell CPM; count-like adata.X is normalized first,
        while CPM-like adata.X is averaged directly.
    cpm_target_sum : float
        Target sum for per-cell CPM when pseudobulk_type='cpm' and adata.X is counts.
    verbose : bool
        If True, print input-state and aggregation messages.
    """
    if pseudobulk_type not in PSEUDOBULK_TYPES:
        raise ValueError(f"pseudobulk_type must be one of {sorted(PSEUDOBULK_TYPES)}")

    if isinstance(ct_group_keys, str):
        ct_group_keys = [ct_group_keys]
    ct_group_keys = list(ct_group_keys)
    if len(ct_group_keys) not in (1, 2):
        raise ValueError("ct_group_keys must contain one or two columns")
    _require_columns(adata.obs, ct_group_keys)

    X = adata.X.tocsr() if sp.issparse(adata.X) else np.asarray(adata.X)
    expression_state = _detect_expression_state(X)

    if verbose:
        print(f"adata.X detected as: {expression_state}")

    if expression_state == "log":
        raise ValueError("adata.X appears log transformed; cannot build pseudobulk")

    if pseudobulk_type == "count" and expression_state != "count":
        raise ValueError(
            "pseudobulk_type='count' requires raw count-like adata.X; "
            f"detected {expression_state!r}"
        )

    adata_use = adata.copy()
    adata_use.X = X

    if pseudobulk_type == "count":
        aggre_func = "sum"
        if verbose:
            print("Using aggregation: sum")
            print("Output meaning: raw count sum per pseudobulk group")
    else:
        aggre_func = "mean"
        if expression_state == "count":
            if verbose:
                print(f"Normalizing each cell to CPM with target_sum={cpm_target_sum:g}")
                print("Using aggregation: mean")
                print("Output meaning: mean CPM per pseudobulk group")
            adata_use = sc.pp.normalize_total(
                adata_use,
                target_sum=cpm_target_sum,
                copy=True,
            )
        else:
            if verbose:
                print("adata.X appears CPM-like; using it directly")
                print("Using aggregation: mean")
                print("Output meaning: mean CPM per pseudobulk group")

    combined_key = ct_group_keys[0]
    if len(ct_group_keys) == 2:
        outer_key, inner_key = ct_group_keys
        combined_key = "__pseudobulk_group__"

        observed_groups = (
            adata_use.obs[[outer_key, inner_key]]
            .astype(str)
            .drop_duplicates()
            .reset_index(drop=True)
        )
        observed_groups[combined_key] = [
            f"pseudobulk_group_{i}" for i in range(observed_groups.shape[0])
        ]

        adata_use.obs[combined_key] = (
            adata_use.obs[[outer_key, inner_key]]
            .astype(str)
            .merge(observed_groups, on=[outer_key, inner_key], how="left")[combined_key]
            .to_numpy()
        )

        if verbose:
            print(
                f"Two-key mode: returning dict keyed by {outer_key!r}, "
                f"with columns from {inner_key!r}"
            )
    else:
        if verbose:
            print(f"One-key mode: returning gene x {ct_group_keys[0]} DataFrame")

    ref_ct = sc.get.aggregate(adata_use, by=combined_key, func=aggre_func)
    mat = ref_ct.layers[aggre_func]
    if sp.issparse(mat):
        mat = mat.toarray()
    mat = np.asarray(mat)

    groups_order = list(ref_ct.obs_names)
    pseudobulk_df = pd.DataFrame(mat, index=groups_order, columns=ref_ct.var_names).T
    pseudobulk_df.index.name = None
    pseudobulk_df.columns.name = None

    if len(ct_group_keys) == 1:
        counts = (
            adata_use.obs[ct_group_keys[0]]
            .value_counts()
            .reindex(groups_order)
            .fillna(0)
            .astype(int)
        )
        pseudobulk_meta = pd.DataFrame({
            "name": groups_order,
            "anno_level": ct_group_keys[0],
            "cluster_name": groups_order,
            "nCells": counts.to_numpy(),
        })

        if verbose:
            print(f"Output pseudobulk shape: {pseudobulk_df.shape[0]} genes x {pseudobulk_df.shape[1]} groups")

        return {
            "pseudobulk": pseudobulk_df,
            "pseudobulk_meta": pseudobulk_meta,
        }

    outer_order = list(pd.unique(adata_use.obs[outer_key].astype(str)))
    inner_order = list(pd.unique(adata_use.obs[inner_key].astype(str)))

    group_lookup = {
        (row[outer_key], row[inner_key]): row[combined_key]
        for _, row in observed_groups.iterrows()
    }

    pseudobulk = {}
    mapping_rows = []
    counts = adata_use.obs.groupby([outer_key, inner_key], observed=True).size()

    for outer in outer_order:
        out = pd.DataFrame(0.0, index=pseudobulk_df.index, columns=inner_order)

        for inner in inner_order:
            n_cells = int(counts.get((outer, inner), 0))
            mapping_rows.append({
                outer_key: outer,
                inner_key: inner,
                "nCells": n_cells,
            })

            group_id = group_lookup.get((outer, inner))
            if group_id in pseudobulk_df.columns:
                out[inner] = pseudobulk_df[group_id]

        out.index.name = None
        out.columns.name = None
        pseudobulk[outer] = out

    group_counts = pd.DataFrame(mapping_rows)

    if verbose:
        print(f"Output contains {len(pseudobulk)} {outer_key} groups")
        print(f"Each matrix shape: {pseudobulk_df.shape[0]} genes x {len(inner_order)} {inner_key} groups")

    return {
        "pseudobulk": pseudobulk,
        "pseudobulk_meta": group_counts,
    }


def compute_weighted_bulk_from_ct_profiles(meancpm_by_ct, cellfrac, log2_transform=False):
    """
    Compute weighted bulk CPM from CTSE meancpm and cell fractions.

    The log2_transform option is retained for isolated wlogcpm experiments under
    scripts/wlogcpm_test; active DALE_Eval runners do not use wlogcpm inputs.
    """
    if not meancpm_by_ct:
        raise ValueError("meancpm_by_ct is empty")
    first = next(iter(meancpm_by_ct.values()))
    genes = first.index
    samples = first.columns
    out = pd.DataFrame(0.0, index=genes, columns=samples)

    for ct, mat in meancpm_by_ct.items():
        if not mat.index.equals(genes) or not mat.columns.equals(samples):
            raise ValueError(f"Matrix for cell type {ct!r} has inconsistent genes/samples")
        weights = cellfrac.reindex(index=samples).get(ct)
        if weights is None:
            weights = pd.Series(0.0, index=samples)
        values = np.log2(mat + 1) if log2_transform else mat
        out = out + values.mul(weights.reindex(samples).fillna(0).to_numpy(), axis=1)

    out.index.name = None
    out.columns.name = None
    return out


def normalize_counts_to_cpm(counts, target_sum=1e6):
    """
    Column-wise CPM normalization for a genes x samples count matrix.

    Columns with zero library size are returned as all zero. Input may be a
    pandas DataFrame or an array-like object; DataFrame index/columns are
    preserved.
    """
    is_df = isinstance(counts, pd.DataFrame)
    values = counts.to_numpy(dtype=float, copy=True) if is_df else np.asarray(counts, dtype=float).copy()
    values[~np.isfinite(values)] = 0
    values[values < 0] = 0

    lib_sizes = values.sum(axis=0)
    out = np.zeros_like(values, dtype=float)
    keep = np.isfinite(lib_sizes) & (lib_sizes > 0)
    out[:, keep] = values[:, keep] / lib_sizes[keep] * float(target_sum)

    if is_df:
        return pd.DataFrame(out, index=counts.index.copy(), columns=counts.columns.copy())
    return out


def normalize_ctse_sumcount_to_cpm(ctse_sumcount_by_ct, target_sum=1e6):
    """
    CPM-normalize each CTSE sumcount matrix independently.

    Parameters
    ----------
    ctse_sumcount_by_ct : dict[str, pd.DataFrame]
        Cell-type keyed genes x samples count matrices.
    target_sum : float
        CPM target sum, usually 1e6.
    """
    return {
        ct: normalize_counts_to_cpm(mat, target_sum=target_sum)
        for ct, mat in ctse_sumcount_by_ct.items()
    }


def _rowwise_pearson_fast(a, b, min_samples=3):
    """
    Fast row-wise Pearson correlation for two same-shaped arrays.

    Non-finite paired positions are ignored. Rows with fewer than
    ``min_samples`` paired finite values or zero variance return NaN.
    """
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    if x.shape != y.shape:
        raise ValueError("a and b must have the same shape")

    keep = np.isfinite(x) & np.isfinite(y)
    xx = np.where(keep, x, 0.0)
    yy = np.where(keep, y, 0.0)

    n = keep.sum(axis=1).astype(float)
    sx = xx.sum(axis=1)
    sy = yy.sum(axis=1)
    sxx = (xx * xx).sum(axis=1)
    syy = (yy * yy).sum(axis=1)
    sxy = (xx * yy).sum(axis=1)

    vx = n * sxx - sx * sx
    vy = n * syy - sy * sy
    out = np.full(x.shape[0], np.nan, dtype=float)

    good = (
        (n >= float(min_samples))
        & np.isfinite(vx)
        & np.isfinite(vy)
        & (vx > 0)
        & (vy > 0)
    )
    out[good] = (n[good] * sxy[good] - sx[good] * sy[good]) / np.sqrt(vx[good] * vy[good])
    out[~np.isfinite(out)] = np.nan
    return out, n.astype(int)


def _rowwise_spearman_fast(a_df, b_df, min_samples=3):
    """
    Fast row-wise Spearman correlation as Pearson correlation of row ranks.
    """
    keep = np.isfinite(a_df.to_numpy(dtype=float)) & np.isfinite(b_df.to_numpy(dtype=float))
    a_rank = a_df.where(keep).rank(axis=1, method="average", na_option="keep")
    b_rank = b_df.where(keep).rank(axis=1, method="average", na_option="keep")
    return _rowwise_pearson_fast(a_rank.to_numpy(), b_rank.to_numpy(), min_samples=min_samples)[0]


def compute_ctse_gene_correlations(
    ctse_a_by_ct,
    ctse_b_by_ct,
    a_label="sumcount_cpm",
    b_label="meancpm",
    min_samples=3,
):
    """
    Compute gene-wise Pearson and Spearman correlations between CTSE matrices.

    Each cell type is compared on overlapping genes and samples. Correlations
    are computed across samples for each gene using vectorized row-wise matrix
    operations, matching the fast Pearson/Spearman approach in evalu.R.
    """
    frames = []
    common_cts = sorted(set(ctse_a_by_ct).intersection(ctse_b_by_ct))
    if not common_cts:
        raise ValueError("No overlapping cell types between CTSE inputs")

    for ct in common_cts:
        a = ctse_a_by_ct[ct]
        b = ctse_b_by_ct[ct]
        common_genes = a.index.intersection(b.index)
        common_samples = a.columns.intersection(b.columns)
        if len(common_genes) == 0 or len(common_samples) == 0:
            continue

        a_use = a.loc[common_genes, common_samples].astype(float)
        b_use = b.loc[common_genes, common_samples].astype(float)

        pearson, n_samples = _rowwise_pearson_fast(
            a_use.to_numpy(),
            b_use.to_numpy(),
            min_samples=min_samples,
        )
        spearman = _rowwise_spearman_fast(
            a_use,
            b_use,
            min_samples=min_samples,
        )

        df = pd.DataFrame({
            "cell_type": ct,
            "gene": common_genes.to_numpy(),
            "pearson": pearson,
            "spearman": spearman,
            "n_samples": n_samples,
            "a_label": a_label,
            "b_label": b_label,
        })
        df = df[np.isfinite(df["pearson"]) | np.isfinite(df["spearman"])]
        if not df.empty:
            frames.append(df)

    if not frames:
        return pd.DataFrame(
            columns=["cell_type", "gene", "pearson", "spearman", "n_samples", "a_label", "b_label"]
        )
    return pd.concat(frames, axis=0, ignore_index=True)

def read_ctse_truth_dir(ctse_dir):
    ctse_dir = Path(ctse_dir)
    out = {}

    for path in sorted(ctse_dir.glob("*.txt.gz")):
        ct = path.name.replace(".txt.gz", "")
        out[ct] = pd.read_csv(path, sep="\t", index_col=0)

    if not out:
        raise FileNotFoundError(f"No .txt.gz files found in {ctse_dir}")

    return out

