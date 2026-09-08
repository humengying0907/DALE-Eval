import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_frac_comparison(
    truth_frac,
    estimated_frac,
    title=None,
    nrow=1,
    xlabel="true fraction",
    ylabel="estimated fraction",
    point_size=10,
    return_data=False,
):
    """
    Compare two sample-by-cell-type fraction matrices.

    Parameters
    ----------
    truth_frac, estimated_frac : pd.DataFrame
        Samples as index and cell types as columns.
    return_data : bool
        If True, return (fig, summary, plot_data). Otherwise return fig.
    """
    common_samples = sorted(set(truth_frac.index).intersection(estimated_frac.index))
    if len(common_samples) == 0:
        raise ValueError("No common samples found between truth_frac and estimated_frac.")

    common_celltypes = sorted(set(truth_frac.columns).intersection(estimated_frac.columns))
    if len(common_celltypes) == 0:
        raise ValueError("No common cell types found between truth_frac and estimated_frac.")

    truth = truth_frac.loc[common_samples, common_celltypes]
    estimate = estimated_frac.loc[common_samples, common_celltypes]

    keep_ct = (truth.sum(axis=0, skipna=True) > 0) & (estimate.sum(axis=0, skipna=True) > 0)
    truth = truth.loc[:, keep_ct]
    estimate = estimate.loc[:, keep_ct]

    if truth.shape[1] == 0:
        raise ValueError("No common non-zero cell types remain after filtering.")

    plot_rows = []
    summary_rows = []

    for ct in truth.columns:
        x = truth[ct].astype(float)
        y = estimate[ct].astype(float)
        valid = x.notna() & y.notna()

        if valid.sum() >= 2 and x[valid].nunique() > 1 and y[valid].nunique() > 1:
            cor = float(np.corrcoef(x[valid], y[valid])[0, 1])
        else:
            cor = np.nan

        rmse = float(np.sqrt(np.mean((x[valid] - y[valid]) ** 2))) if valid.any() else np.nan
        lim = float(np.nanmax([x.max(skipna=True), y.max(skipna=True), 0]))

        summary_rows.append({
            "cell_type": ct,
            "cor": cor,
            "RMSE": rmse,
            "lim": lim,
        })

        plot_rows.append(pd.DataFrame({
            "sample": common_samples,
            "cell_type": ct,
            "true_frac": x.to_numpy(),
            "estimate": y.to_numpy(),
        }))

    plot_data = pd.concat(plot_rows, axis=0, ignore_index=True)
    summary = pd.DataFrame(summary_rows)

    n_ct = len(summary)
    nrow = max(1, int(nrow))
    ncol = int(np.ceil(n_ct / nrow))

    fig, axes = plt.subplots(
        nrow,
        ncol,
        figsize=(3.2 * ncol, 3.2 * nrow),
        squeeze=False,
    )

    for ax, (_, row) in zip(axes.ravel(), summary.iterrows()):
        ct = row["cell_type"]
        df = plot_data[plot_data["cell_type"] == ct]
        lim = row["lim"]

        if not np.isfinite(lim) or lim <= 0:
            lim = 1.0

        ax.scatter(df["true_frac"], df["estimate"], s=point_size)
        ax.plot([0, lim], [0, lim], color="red", linestyle=":", linewidth=1)

        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(ct)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        ax.text(
            0.02 * lim,
            0.98 * lim,
            f"cor = {row['cor']:.2f}\nRMSE = {row['RMSE']:.2f}",
            ha="left",
            va="top",
            fontsize=9,
        )

    for ax in axes.ravel()[n_ct:]:
        ax.set_visible(False)

    if title is not None:
        fig.suptitle(title)

    fig.tight_layout()

    if return_data:
        return fig, summary, plot_data
    return fig


def plot_ctse_gene_correlation_histograms(
    gene_correlations,
    metrics=("pearson", "spearman"),
    bins=40,
    ncol=3,
    title="CTSE gene-wise correlation",
    alpha=0.75,
    stacked=True,
    xlim=None,
    reference_zero=False,
):
    """
    Plot histograms of gene-wise CTSE correlations by cell type.

    Parameters
    ----------
    gene_correlations : pd.DataFrame
        Output from ``compute_ctse_gene_correlations`` with columns
        cell_type, pearson, and spearman.
    metrics : tuple[str]
        Correlation columns to show.
    bins : int
        Histogram bin count. When xlim is None, bins are chosen separately for
        each cell type using that cell type's observed correlation range.
    ncol : int
        Number of facet columns.
    title : str or None
        Optional figure title.
    alpha : float
        Histogram transparency.
    stacked : bool
        If True, stack Pearson and Spearman histograms. If False, overlay them.
    xlim : tuple[float, float] or None
        Optional shared x-axis limits. Use None for free x-axis scaling.
    reference_zero : bool
        If True, draw a vertical line at zero. Defaults to False so plots are
        not visually anchored around zero when all correlations are high.
    """
    if gene_correlations is None or gene_correlations.empty:
        raise ValueError("gene_correlations is empty")

    metrics = tuple(metrics)
    missing = [m for m in metrics if m not in gene_correlations.columns]
    if missing:
        raise KeyError(f"Missing metric columns: {missing}")
    if "cell_type" not in gene_correlations.columns:
        raise KeyError("gene_correlations must contain a 'cell_type' column")

    cell_types = list(pd.unique(gene_correlations["cell_type"]))
    n_ct = len(cell_types)
    ncol = max(1, int(ncol))
    nrow = int(np.ceil(n_ct / ncol))

    fig, axes = plt.subplots(
        nrow,
        ncol,
        figsize=(4.0 * ncol, 3.0 * nrow),
        squeeze=False,
        sharex=xlim is not None,
        sharey=False,
    )

    colors = dict(pearson="#4C78A8", spearman="#F58518")

    for ax, ct in zip(axes.ravel(), cell_types):
        df = gene_correlations[gene_correlations["cell_type"] == ct]
        metric_values = []
        metric_labels = []
        metric_colors = []

        for metric in metrics:
            values = pd.to_numeric(df[metric], errors="coerce")
            values = values[np.isfinite(values)]
            if values.empty:
                continue
            metric_values.append(values.to_numpy())
            metric_labels.append(metric)
            metric_colors.append(colors.get(metric))

        if metric_values:
            combined = np.concatenate(metric_values)
            if xlim is None:
                lo = float(np.nanmin(combined))
                hi = float(np.nanmax(combined))
                if not np.isfinite(lo) or not np.isfinite(hi):
                    lo, hi = -1.0, 1.0
                elif lo == hi:
                    pad = max(0.01, abs(lo) * 0.02)
                    lo, hi = lo - pad, hi + pad
                else:
                    pad = 0.04 * (hi - lo)
                    lo, hi = max(-1.0, lo - pad), min(1.0, hi + pad)
                edges = np.linspace(lo, hi, int(bins) + 1)
            else:
                lo, hi = xlim
                edges = np.linspace(lo, hi, int(bins) + 1)

            if stacked:
                ax.hist(
                    metric_values,
                    bins=edges,
                    stacked=True,
                    label=metric_labels,
                    color=metric_colors,
                    alpha=alpha,
                    edgecolor="white",
                    linewidth=0.3,
                )
            else:
                for values, metric, color in zip(metric_values, metric_labels, metric_colors):
                    ax.hist(
                        values,
                        bins=edges,
                        alpha=alpha,
                        label=metric,
                        color=color,
                        edgecolor="white",
                        linewidth=0.3,
                    )
            ax.set_xlim(lo, hi)
        elif xlim is not None:
            ax.set_xlim(*xlim)

        medians = []
        for metric in metrics:
            values = pd.to_numeric(df[metric], errors="coerce")
            values = values[np.isfinite(values)]
            if not values.empty:
                medians.append(f"{metric} med={values.median():.2f}")
        subtitle = "\n".join(medians)
        ax.set_title(f"{ct}\n{subtitle}" if subtitle else ct, fontsize=10)
        if reference_zero:
            ax.axvline(0, color="black", linestyle=":", linewidth=0.8)
        ax.set_xlabel("Gene-wise correlation across samples")
        ax.set_ylabel("Genes")

    for ax in axes.ravel()[n_ct:]:
        ax.set_visible(False)

    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right", frameon=False)
    if title is not None:
        fig.suptitle(title)
    fig.tight_layout(rect=(0, 0, 0.96, 0.96) if title is not None else None)
    return fig

