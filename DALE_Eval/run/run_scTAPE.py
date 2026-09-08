#!/usr/bin/env python3
from __future__ import annotations

import argparse
import atexit
import copy
import csv
import gzip
import json
import os
import random
import shlex
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

RUN_DIR = Path(__file__).resolve().parent
MODULE_DIR = RUN_DIR.parent / "modules"
sys.path.insert(0, str(MODULE_DIR))

from runner_helpers import (  # noqa: E402
    apply_method_default_extra_args,
    extra_arg,
    find_repo_root,
    finish_runtime_log,
    format_deconv_dim,
    method_message_log_path,
    parse_extra_args,
    parse_simple_blue_hyperparameters,
    prepare_bulk_dataframe_for_deconv,
    read_bulk_dataframe,
    read_bulk_header,
    read_indep_ref_cell_type_mapping,
    read_test_samples,
    resolve_deconv_paths,
    start_runtime_log,
    validate_extra_args,
)

METHOD = "scTAPE"
INDEP_REF_SCRNA_SOURCE = {"PBMC_refined_AIDA2024": "PBMC_AIDA2024"}


class TeeStream:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data: str) -> int:
        for stream in self.streams:
            stream.write(data)
            stream.flush()
        return len(data)

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()

    def isatty(self) -> bool:
        return any(getattr(stream, "isatty", lambda: False)() for stream in self.streams)

    def __getattr__(self, name: str):
        return getattr(self.streams[0], name)


def start_message_log(log_path: Path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_file = log_path.open("a")
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    sys.stdout = TeeStream(old_stdout, log_file)
    sys.stderr = TeeStream(old_stderr, log_file)
    print(f"message_log: {log_path}")
    print("command: " + shlex.join(["python", *sys.argv]))

    def cleanup() -> None:
        sys.stdout = old_stdout
        sys.stderr = old_stderr
        log_file.close()

    return cleanup


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run scTAPE through the benchmark config interface.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--config_id", required=True)
    parser.add_argument("--n_core", type=int, default=15, help="Accepted for benchmark interface consistency; scTAPE does not use it directly.")
    parser.add_argument("--extra_args", default="")
    return parser.parse_args()


def bool_extra(value: Any, key: str) -> bool:
    if isinstance(value, bool):
        return value
    raise ValueError(f"{key} must be true or false")


def positive_int_extra(value: Any, key: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise ValueError(f"{key} must be a positive integer")
    return value


def reference_scRNA_source(reference_name: str) -> str:
    return INDEP_REF_SCRNA_SOURCE.get(reference_name, reference_name)


def find_reference_h5ad(repo_root: Path, reference_name: str) -> Path:
    source_name = reference_scRNA_source(reference_name)
    ref_sc_dir = repo_root / "scRNA_datasets" / source_name
    if not ref_sc_dir.is_dir():
        raise FileNotFoundError(f"scRNA dataset directory not found for scTAPE independent reference {reference_name}: {ref_sc_dir}")

    candidate = ref_sc_dir / f"{source_name}_processed.h5ad"
    if candidate.exists():
        return candidate

    globbed = sorted(ref_sc_dir.glob("*_processed.h5ad"))
    if len(globbed) == 1:
        return globbed[0]
    if not globbed:
        raise FileNotFoundError(f"No *_processed.h5ad found for scTAPE independent reference {reference_name} in {ref_sc_dir}")
    raise ValueError(f"Multiple candidate scTAPE independent-reference h5ads found in {ref_sc_dir}: {globbed}")


def load_sctape_hyperparameters(repo_root: Path) -> dict[str, dict[str, Any]]:
    hyper_path = repo_root / "DALE_Eval" / "configs" / "scTAPE_hyperparameters.yaml"
    return parse_simple_blue_hyperparameters(hyper_path)


def resolve_repo_path(repo_root: Path, value: Any, label: str) -> Path:
    p = Path(str(value)).expanduser()
    if not p.is_absolute():
        p = repo_root / p
    if not p.exists():
        raise FileNotFoundError(f"{label} not found: {p}")
    return p


def parse_celltype_mapping(path: Path) -> dict[str, str]:
    fine_to_coarse: dict[str, str] = {}
    current: str | None = None
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        if not line.startswith(" "):
            if not line.endswith(":"):
                raise ValueError(f"Invalid mapping line in {path}: {raw_line}")
            current = line[:-1].strip().strip("'\"")
            if not current:
                raise ValueError(f"Empty coarse cell type in {path}")
            continue
        if current is None:
            continue
        stripped = line.strip()
        if not stripped.startswith("-"):
            continue
        fine = stripped[1:].strip().strip("'\"")
        if not fine:
            continue
        if fine in fine_to_coarse:
            raise ValueError(f"Fine cell type {fine!r} appears multiple times in {path}")
        fine_to_coarse[fine] = current
    if not fine_to_coarse:
        raise ValueError(f"No fine cell types found in mapping YAML: {path}")
    return fine_to_coarse


def resolve_celltype_mapping(
    repo_root: Path,
    reference_mapping_dir: Path,
    hyper: dict[str, Any],
    extra_mapping_yaml: Any,
) -> tuple[Path, str, dict[str, str]]:
    config_mapping = hyper.get("cell_types", {}).get("mapping_yaml")
    mapping_value = extra_mapping_yaml if extra_mapping_yaml is not None else config_mapping
    if mapping_value is None:
        mapping_path = reference_mapping_dir / "celltype_mapping.yaml"
        if not mapping_path.exists():
            raise FileNotFoundError(f"scTAPE cell type mapping YAML not found: {mapping_path}")
        return mapping_path, "default", parse_celltype_mapping(mapping_path)
    mapping_path = resolve_repo_path(repo_root, mapping_value, "scTAPE cell type mapping YAML")
    return mapping_path, "extra_args" if extra_mapping_yaml is not None else "config", parse_celltype_mapping(mapping_path)



def maybe_subset_reference_h5ad(
    reference_h5ad: Path,
    workdir: Path,
    cell_type_col: str,
    sample_col: str,
    celltype_mapping: dict[str, str],
    ref_sample_subset: bool,
    ref_sample_subset_n: int | None,
) -> tuple[Path, dict[str, Any] | None]:
    if not ref_sample_subset:
        return reference_h5ad, None
    if ref_sample_subset_n is None:
        raise ValueError("ref_sample_subset=true requires ref_sample_subset_n=<positive integer>")

    import anndata as ad

    subset_path = workdir / "reference_subset.h5ad"
    report_path = workdir / "reference_subset_report.json"
    adata = ad.read_h5ad(reference_h5ad, backed="r")
    try:
        if cell_type_col not in adata.obs.columns:
            raise ValueError(f"input h5ad missing obs[{cell_type_col!r}]: {reference_h5ad}")
        if sample_col not in adata.obs.columns:
            raise ValueError(f"ref_sample_subset=true requires obs[{sample_col!r}] in h5ad: {reference_h5ad}")

        sample_values = adata.obs[sample_col].astype(str)
        unique_samples = list(dict.fromkeys(sample_values.tolist()))
        if ref_sample_subset_n > len(unique_samples):
            raise ValueError(
                f"ref_sample_subset_n={ref_sample_subset_n} but only {len(unique_samples)} unique samples exist in obs[{sample_col!r}]"
            )
        selected_samples = unique_samples[:ref_sample_subset_n]
        selected_set = set(selected_samples)
        mask = sample_values.isin(selected_set).to_numpy()
        if int(mask.sum()) == 0:
            raise ValueError("Reference sample subset selected zero cells")

        selected_cell_types = sorted(set(adata.obs.loc[mask, cell_type_col].astype(str).tolist()))
        selected_fine = set(selected_cell_types)
        if celltype_mapping:
            coarse_to_fine: dict[str, list[str]] = {}
            for fine, coarse in celltype_mapping.items():
                coarse_to_fine.setdefault(coarse, []).append(fine)
            missing_classes = {coarse: values for coarse, values in coarse_to_fine.items() if not (selected_fine & set(values))}
            if missing_classes:
                formatted = "; ".join(f"{k}: {','.join(v)}" for k, v in missing_classes.items())
                raise ValueError("Reference sample subset is missing mapped output classes: " + formatted)

        subset = adata[mask, :].to_memory()
        subset_path.parent.mkdir(parents=True, exist_ok=True)
        subset.write_h5ad(subset_path)
    finally:
        if hasattr(adata, "file") and adata.file is not None:
            adata.file.close()

    report = {
        "original_h5ad": str(reference_h5ad),
        "subset_h5ad": str(subset_path),
        "sample_col": sample_col,
        "cell_type_col": cell_type_col,
        "requested_samples": ref_sample_subset_n,
        "available_samples": len(unique_samples),
        "selected_samples": selected_samples,
        "selected_sample_count": len(selected_samples),
        "selected_cell_count": int(mask.sum()),
        "selected_cell_types": selected_cell_types,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True))
    return subset_path, report


def top_variable_genes_from_reference(adata, bulk_genes: set[str], max_genes: Any) -> list[str] | None:
    if max_genes is None:
        return None
    if not isinstance(max_genes, int) or max_genes <= 0:
        raise ValueError("scTAPE genes.max_genes must be null or a positive integer")

    import numpy as np
    import scipy.sparse as sp

    overlap = [g for g in adata.var_names.astype(str) if g in bulk_genes]
    if len(overlap) <= max_genes:
        return overlap
    sub = adata[:, overlap]
    X = sub.X
    if sp.issparse(X):
        mean = np.asarray(X.mean(axis=0)).ravel()
        mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
        var = mean_sq - mean * mean
    else:
        var = np.asarray(X).var(axis=0)
    order = np.argsort(var)[::-1][:max_genes]
    return [overlap[i] for i in order]


def prepare_reference_adata(
    reference_h5ad: Path,
    celltype_mapping: dict[str, str],
    bulk_genes: set[str],
    max_genes: Any,
    cell_type_col: str,
):
    import anndata as ad
    import scipy.sparse as sp

    adata = ad.read_h5ad(reference_h5ad)
    if cell_type_col not in adata.obs.columns:
        raise ValueError(f"input h5ad missing obs[{cell_type_col!r}]: {reference_h5ad}")

    source = adata.obs[cell_type_col].astype(str)
    if celltype_mapping:
        missing = sorted(set(source.unique()) - set(celltype_mapping))
        if missing:
            print("WARNING: scTAPE dropping fine labels not in mapping: " + ", ".join(missing))
        keep = source.isin(celltype_mapping).to_numpy()
        if int(keep.sum()) == 0:
            raise ValueError("scTAPE cell type mapping kept zero cells")
        adata = adata[keep].copy()
        source = source.loc[adata.obs_names]
        adata.obs["CellType"] = source.map(celltype_mapping).astype(str).values
    else:
        adata.obs["CellType"] = source.values

    selected_genes = top_variable_genes_from_reference(adata, bulk_genes, max_genes)
    if selected_genes is not None:
        adata = adata[:, selected_genes].copy()

    if sp.issparse(adata.X):
        adata.X = adata.X.astype("float32")
    else:
        adata.X = adata.X.astype("float32", copy=False)
    return adata, selected_genes


def sanitize_cell_type_filename(cell_type: str) -> str:
    return cell_type.replace("/", "_") + ".txt.gz"


def export_sctape_outputs(cell_type_sigm: dict[str, Any], pred_frac, output_dir: Path, mapping: dict[str, str]) -> None:
    method_dir = output_dir / METHOD
    if method_dir.exists() and any(method_dir.iterdir()):
        print(f"  overwrite note: existing {METHOD} outputs found for this config; files in {method_dir} may be overwritten")
    method_dir.mkdir(parents=True, exist_ok=True)

    cell_types = list(cell_type_sigm.keys())
    mapped = [mapping.get(ct, ct) for ct in cell_types]
    dupes = sorted({ct for ct in mapped if mapped.count(ct) > 1})
    if dupes:
        raise ValueError("Cell type mapping creates duplicated scTAPE output names: " + ", ".join(dupes))

    for source_ct, out_ct in zip(cell_types, mapped):
        df = cell_type_sigm[source_ct].copy()  # samples x genes
        out_path = method_dir / sanitize_cell_type_filename(out_ct)
        with gzip.open(out_path, "wt", newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(["", *[str(x) for x in df.index]])
            for gene in df.columns:
                values = [round(float(v), 2) if float(v) != 0 else 0 for v in df[gene].values]
                writer.writerow([gene, *values])

    frac = pred_frac.loc[:, cell_types].copy()
    frac.columns = mapped
    frac.round(6).to_csv(method_dir / "scTAPEfrac.txt", sep="\t", index_label="")
    print(f"exported {len(cell_types)} CTSE files to {method_dir}")
    print(f"exported fraction file to {method_dir / 'scTAPEfrac.txt'}")



def save_sctape_loss_history(loss: list[Any], recon_loss: list[Any], workdir: Path) -> Path:
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    loss_path = workdir / "loss_history.csv"
    plot_path = workdir / "loss_history.png"
    pred_loss = [float(x) for x in loss]
    reconstruction_loss = [float(x) for x in recon_loss]
    df = pd.DataFrame(
        {
            "iteration": list(range(1, len(pred_loss) + 1)),
            "prediction_loss": pred_loss,
            "reconstruction_loss": reconstruction_loss,
        }
    )
    df.to_csv(loss_path, index=False)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(df["iteration"], df["prediction_loss"], linewidth=1.2, label="prediction_loss")
    ax.plot(df["iteration"], df["reconstruction_loss"], linewidth=1.2, label="reconstruction_loss")
    ax.set_xlabel("training iteration")
    ax.set_ylabel("loss")
    ax.set_title("scTAPE training convergence")
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=160)
    plt.close(fig)
    print(f"  loss_history: {loss_path}")
    print(f"  loss_plot: {plot_path}")
    return plot_path


def run_log_loss_plot_path(message_log_path: Path) -> Path:
    return message_log_path.with_name(f"{message_log_path.stem}_loss.png")


def copy_loss_plot_to_run_logs(loss_plot_path: Path | None, run_log_plot_path: Path) -> None:
    if loss_plot_path is None or not loss_plot_path.exists():
        print("  run_log_loss_plot: not available")
        return
    run_log_plot_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(loss_plot_path, run_log_plot_path)
    print(f"  run_log_loss_plot: {run_log_plot_path}")


def cleanup_temp_workdir(workdir: Path, keep_temp: bool) -> None:
    if keep_temp:
        print(f"  temp_workdir: kept at {workdir}")
        return
    shutil.rmtree(workdir, ignore_errors=True)
    print(f"  temp_workdir: removed {workdir}")


def run_sctape(
    reference_h5ad: Path,
    bulk_df,
    hyper: dict[str, dict[str, Any]],
    celltype_mapping: dict[str, str],
    workdir: Path,
    cell_type_col: str,
):
    import numpy as np
    import torch
    from torch.optim import Adam

    tape_root = RUN_DIR.parent / "external_modules" / "TAPE"
    sys.path.insert(0, str(tape_root))

    from torch.utils.data import DataLoader

    from TAPE import utils as tape_utils  # noqa: E402
    from TAPE.model import AutoEncoder, device, simdatset  # noqa: E402
    from TAPE.simulation import generate_simulated_data  # noqa: E402
    from TAPE.train import adaptive_stage, reproducibility, training_stage  # noqa: E402

    tape_utils.plt.show = lambda *args, **kwargs: None

    sim_cfg = hyper.get("simulation", {})
    prep_cfg = hyper.get("preprocess", {})
    train_cfg = hyper.get("training", {})
    pred_cfg = hyper.get("prediction", {})
    gene_cfg = hyper.get("genes", {})

    seed = int(train_cfg.get("seed", 0))
    reproducibility(seed)
    random.seed(seed)
    np.random.seed(seed)

    adata, selected_genes = prepare_reference_adata(
        reference_h5ad,
        celltype_mapping,
        set(map(str, bulk_df.columns)),
        gene_cfg.get("max_genes"),
        cell_type_col,
    )
    if selected_genes is not None:
        bulk_df = bulk_df.loc[:, [g for g in selected_genes if g in bulk_df.columns]]

    print("Generating scTAPE simulated pseudobulk training data")
    print(f"  sc_reference: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"  scTAPE cell types ({adata.obs['CellType'].nunique()}): {sorted(adata.obs['CellType'].astype(str).unique().tolist())}")
    adapted_ref_path = workdir / "scTAPE_reference.h5ad"
    adata.write_h5ad(adapted_ref_path)
    print(f"  scTAPE_reference: {adapted_ref_path}")
    simudata = generate_simulated_data(
        sc_data=str(adapted_ref_path),
        samplenum=int(sim_cfg.get("samplenum", 5000)),
        d_prior=None,
        n=int(sim_cfg.get("n_cells", 500)),
        random_state=sim_cfg.get("random_state"),
        sparse=bool(sim_cfg.get("sparse", True)),
        sparse_prob=float(sim_cfg.get("sparse_prob", 0.5)),
        rare=bool(sim_cfg.get("rare", False)),
        rare_percentage=float(sim_cfg.get("rare_percentage", 0.4)),
    )
    sim_path = workdir / "simulated_training.h5ad"
    simudata.write_h5ad(sim_path)
    print(f"  simulated_training: {simudata.n_obs} samples x {simudata.n_vars} genes")

    print("Preparing scTAPE model matrices")
    train_x, train_y, test_x, genename, celltypes, samplename = tape_utils.ProcessInputData(
        train_data=simudata,
        test_data=bulk_df,
        sep="\t",
        datatype=str(prep_cfg.get("datatype", "counts")),
        genelenfile=None,
        variance_threshold=float(prep_cfg.get("variance_threshold", 0.98)),
        scaler=str(prep_cfg.get("scaler", "mms")),
    )
    print("training data shape is ", train_x.shape, "\ntest data shape is ", test_x.shape)

    batch_size = int(train_cfg.get("batch_size", 128))
    epochs = int(train_cfg.get("epochs", 128))
    print("Start training")
    print(f"  batch_size: {batch_size}")
    print(f"  epochs: {epochs}")
    train_loader = DataLoader(simdatset(train_x, train_y), batch_size=batch_size, shuffle=True)
    model = AutoEncoder(train_x.shape[1], train_y.shape[1]).to(device)
    optimizer = Adam(model.parameters(), lr=1e-4)
    model, loss, recon_loss = training_stage(model, train_loader, optimizer, epochs=epochs)
    print("Training is done")
    loss_plot_path = save_sctape_loss_history(loss, recon_loss, workdir)

    adaptive = bool(pred_cfg.get("adaptive", True))
    mode = str(pred_cfg.get("mode", "high-resolution"))
    if not adaptive or mode != "high-resolution":
        raise ValueError("The benchmark scTAPE runner currently requires prediction.adaptive=true and mode=high-resolution for CTSE export")

    step = int(pred_cfg.get("adaptive_step", 300))
    max_iter = int(pred_cfg.get("adaptive_max_iter", 3))
    print("Starting scTAPE high-resolution adaptive prediction")
    print(f"  adaptive_step: {step}")
    print(f"  adaptive_max_iter: {max_iter}")
    test_sigm_list = np.zeros((test_x.shape[0], len(celltypes), len(genename)), dtype=np.float32)
    test_pred = np.zeros((test_x.shape[0], len(celltypes)), dtype=np.float32)
    for i in range(len(test_x)):
        print(f"  sample {i + 1}/{len(test_x)}: {samplename[i]}")
        sample_model = copy.deepcopy(model)
        decoder_parameters = [{"params": [p for n, p in sample_model.named_parameters() if "decoder" in n]}]
        encoder_parameters = [{"params": [p for n, p in sample_model.named_parameters() if "encoder" in n]}]
        optimizer_d = Adam(decoder_parameters, lr=1e-4)
        optimizer_e = Adam(encoder_parameters, lr=1e-4)
        sigm, _, pred = adaptive_stage(sample_model, test_x[i, :].reshape(1, -1), optimizer_d, optimizer_e, step=step, max_iter=max_iter)
        test_sigm_list[i, :, :] = sigm.astype(np.float32)
        test_pred[i, :] = pred.astype(np.float32)

    import pandas as pd

    pred_frac = pd.DataFrame(test_pred, columns=celltypes, index=samplename)
    cell_type_sigm = {}
    for ci, cell_type in enumerate(celltypes):
        cell_type_sigm[str(cell_type)] = pd.DataFrame(test_sigm_list[:, ci, :], columns=genename, index=samplename)
    return cell_type_sigm, pred_frac, len(genename), loss_plot_path


def main() -> None:
    args = parse_args()
    repo_root = find_repo_root(RUN_DIR)
    extra = apply_method_default_extra_args(parse_extra_args(args.extra_args), "scTAPE", repo_root)
    validate_extra_args(extra, {"map_cell_types", "use_test_samples", "cell_type_col", "sample_col", "celltype_mapping_yaml", "keep_temp", "ref_sample_subset", "ref_sample_subset_n"}, METHOD)

    dataset = args.dataset
    config_id = args.config_id
    n_core = args.n_core
    map_cell_types = bool_extra(extra_arg(extra, "map_cell_types", True), "map_cell_types")
    use_test_samples = bool_extra(extra_arg(extra, "use_test_samples", True), "use_test_samples")
    keep_temp = bool_extra(extra_arg(extra, "keep_temp", False), "keep_temp")
    ref_sample_subset = bool_extra(extra_arg(extra, "ref_sample_subset", False), "ref_sample_subset")
    ref_sample_subset_n_value = extra_arg(extra, "ref_sample_subset_n", None)
    if ref_sample_subset_n_value is not None and not ref_sample_subset:
        raise ValueError("ref_sample_subset_n was provided but ref_sample_subset is not true")
    ref_sample_subset_n = positive_int_extra(ref_sample_subset_n_value, "ref_sample_subset_n") if ref_sample_subset_n_value is not None else None
    cell_type_col = str(extra_arg(extra, "cell_type_col", "cell_type"))
    sample_col = str(extra_arg(extra, "sample_col", "sample"))
    if not cell_type_col:
        raise ValueError("cell_type_col must be a non-empty obs column name")
    if not sample_col:
        raise ValueError("sample_col must be a non-empty obs column name")
    run_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    paths = resolve_deconv_paths(dataset, config_id, repo_root)
    if paths.ref_type == "self":
        raise ValueError(
            "scTAPE does not support refType=self in this benchmark runner. "
            "scTAPE trains from scRNA-derived pseudobulks, and self-reference splits may leave too few cells for stable training. "
            "Use an independent reference config instead."
        )
    ref_name = paths.ref_dir.name
    reference_h5ad = find_reference_h5ad(repo_root, ref_name)
    original_reference_h5ad = reference_h5ad

    workdir = paths.obj_dir / "logs" / f"scTAPE_temp_{config_id}_{run_stamp}"
    workdir.mkdir(parents=True, exist_ok=True)
    message_log_path = method_message_log_path(paths, METHOD, stamp=run_stamp)
    cleanup_message_log = start_message_log(message_log_path)
    atexit.register(cleanup_message_log)

    first_col, all_samples = read_bulk_header(paths.bulk_path)
    del first_col
    selected_samples: list[str] | None = None
    if use_test_samples:
        test_samples = read_test_samples(paths)
        selected_samples = [sample for sample in all_samples if sample in set(test_samples)]
        if not selected_samples:
            raise ValueError("use_test_samples=true but no test samples overlap with bulk samples")

    bulk_df_gene_sample = read_bulk_dataframe(paths.bulk_path)
    bulk_prep = prepare_bulk_dataframe_for_deconv(bulk_df_gene_sample, paths, METHOD)
    bulk_df_gene_sample = bulk_prep["bulk_expr"]
    if selected_samples is not None:
        bulk_df_gene_sample = bulk_df_gene_sample.loc[:, selected_samples]
    bulk_df = bulk_df_gene_sample.T
    bulk_input_dim = format_deconv_dim(bulk_df.shape[1], bulk_df.shape[0])

    hyper = load_sctape_hyperparameters(repo_root)
    mapping_path, mapping_source, celltype_mapping = resolve_celltype_mapping(
        repo_root,
        paths.ref_dir,
        hyper,
        extra_arg(extra, "celltype_mapping_yaml", None),
    )
    reference_h5ad, subset_report = maybe_subset_reference_h5ad(
        reference_h5ad=reference_h5ad,
        workdir=workdir,
        cell_type_col=cell_type_col,
        sample_col=sample_col,
        celltype_mapping=celltype_mapping,
        ref_sample_subset=ref_sample_subset,
        ref_sample_subset_n=ref_sample_subset_n,
    )

    print("Running scTAPE")
    print(f"  dataset: {dataset}")
    print(f"  config_id: {config_id}")
    print(f"  bulk_input: {paths.config['bulk_input']} ({bulk_input_dim})")
    print(f"  bulk_scale: {paths.config['bulk_scale']}")
    print(f"  bulk_normalization: {paths.config['bulk_normalization']}")
    print(f"  bulk_preparation: {bulk_prep['action']}")
    print(f"  config_frac_input: {paths.config['frac_input']}")
    print(f"  refType: {paths.ref_type}")
    print(f"  reference_h5ad_original: {original_reference_h5ad}")
    print(f"  reference_h5ad_effective: {reference_h5ad}")
    print(f"  h5ad_cell_type_col: {cell_type_col}")
    print(f"  h5ad_sample_col: {sample_col}")
    print(f"  ref_sample_subset: {ref_sample_subset}")
    if subset_report is not None:
        print(f"  ref_sample_subset_n: {subset_report['requested_samples']}")
        print(f"  ref_sample_subset_selected: {subset_report['selected_sample_count']} / {subset_report['available_samples']} samples")
        print(f"  ref_sample_subset_cells: {subset_report['selected_cell_count']}")
        print(f"  ref_sample_subset_cell_types: {', '.join(subset_report['selected_cell_types'])}")
    print(f"  scTAPE_mapping: {mapping_path}")
    print(f"  scTAPE_mapping_source: {mapping_source}")
    print(f"  workdir: {workdir}")
    print(f"  map_cell_types: {map_cell_types}")
    print(f"  use_test_samples: {use_test_samples}")
    print(f"  keep_temp: {keep_temp}")
    print("  note: scTAPE does not consume frac_input; it estimates fractions from bulk + scRNA reference")
    print("Starting deconvolution")
    print(f"  deconv_input: {bulk_input_dim}")

    runtime_entry = start_runtime_log(paths, METHOD, n_core=n_core, deconv_input=bulk_input_dim, log_path=str(message_log_path))
    try:
        cell_type_sigm, pred_frac, n_model_genes, loss_plot_path = run_sctape(reference_h5ad, bulk_df, hyper, celltype_mapping, workdir, cell_type_col)
        export_mapping = read_indep_ref_cell_type_mapping(paths, repo_root) if map_cell_types else {}
        export_sctape_outputs(cell_type_sigm, pred_frac, paths.output_dir, export_mapping)
    except KeyboardInterrupt:
        finish_runtime_log(paths, runtime_entry, "killed")
        raise
    except BaseException:
        finish_runtime_log(paths, runtime_entry, "failed")
        raise
    else:
        finish_runtime_log(paths, runtime_entry, "completed")
        copy_loss_plot_to_run_logs(loss_plot_path, run_log_loss_plot_path(message_log_path))
        cleanup_temp_workdir(workdir, keep_temp)

    print("Output saved")
    print(f"  model_genes: {n_model_genes}")
    print(f"  scTAPE ctse: {paths.output_dir / METHOD / '<cell_type>.txt.gz'}")
    print(f"  scTAPE fractions: {paths.output_dir / METHOD / 'scTAPEfrac.txt'}")
    if keep_temp:
        print(f"  scTAPE workdir: {workdir}")


if __name__ == "__main__":
    main()
