#!/usr/bin/env python3
from __future__ import annotations

import argparse
import atexit
import json
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

RUN_DIR = Path(__file__).resolve().parent
MODULE_DIR = RUN_DIR.parent / "modules"
sys.path.insert(0, str(MODULE_DIR))

from runner_helpers import (  # noqa: E402
    apply_method_default_extra_args,
    DeconvPaths,
    extra_arg,
    find_repo_root,
    finish_runtime_log,
    format_deconv_dim,
    parse_extra_args,
    method_message_log_path,
    parse_simple_blue_hyperparameters,
    read_bulk_header,
    read_indep_ref_cell_type_mapping,
    read_test_samples,
    resolve_deconv_paths,
    start_runtime_log,
    validate_extra_args,
    write_prepared_bulk_tsv_for_deconv,
    write_simple_yaml,
)


METHOD = "BLUE"
CTKEY = "cell.type.for.deconv"
INDEP_REF_SCRNA_SOURCE = {"PBMC_refined_AIDA2024": "PBMC_AIDA2024"}


def resolve_uv_command(extra_path: str | None = None) -> list[str]:
    if extra_path:
        if extra_path == "python -m uv":
            return [sys.executable, "-m", "uv"]
        if extra_path.startswith("python:"):
            return [extra_path.split(":", 1)[1], "-m", "uv"]
        p = Path(extra_path).expanduser()
        if p.exists() and os.access(p, os.X_OK):
            return [str(p)]

    found = shutil.which("uv")
    if found:
        return [found]

    for candidate in [Path.home() / ".local" / "bin" / "uv", Path.home() / ".cargo" / "bin" / "uv"]:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return [str(candidate)]

    try:
        subprocess.run(
            [sys.executable, "-m", "uv", "--version"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return [sys.executable, "-m", "uv"]
    except Exception:
        pass

    raise FileNotFoundError(
        "Could not find uv. Run `command -v uv` or `python -m uv --version`; "
        "then either add uv to PATH, pass `uv_path=/absolute/path/to/uv`, "
        "or pass `uv_path=python:/absolute/path/to/python` for a Python where uv is installed."
    )


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
    parser = argparse.ArgumentParser(description="Run BLUE through the benchmark config interface.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--config_id", required=True)
    parser.add_argument("--n_core", type=int, default=15, help="Accepted for benchmark interface consistency; BLUE does not use it directly.")
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
        raise FileNotFoundError(f"scRNA dataset directory not found for BLUE independent reference {reference_name}: {ref_sc_dir}")

    candidate = ref_sc_dir / f"{source_name}_processed.h5ad"
    if candidate.exists():
        return candidate

    globbed = sorted(ref_sc_dir.glob("*_processed.h5ad"))
    if len(globbed) == 1:
        return globbed[0]
    if not globbed:
        raise FileNotFoundError(f"No *_processed.h5ad found for BLUE independent reference {reference_name} in {ref_sc_dir}")
    raise ValueError(f"Multiple candidate BLUE independent-reference h5ads found in {ref_sc_dir}: {globbed}")


def validate_reference_h5ad_columns(
    blue_dir: Path,
    uv_cmd: list[str],
    h5ad_path: Path,
    workdir: Path,
    cell_type_col: str,
    sample_col: str,
) -> None:
    code = r"""
from pathlib import Path
import sys
import anndata as ad

h5ad = Path(sys.argv[1])
cell_type_col = sys.argv[2]
sample_col = sys.argv[3]
a = ad.read_h5ad(h5ad, backed="r")
if cell_type_col not in a.obs.columns:
    raise SystemExit(f"input h5ad missing obs[{cell_type_col!r}]: {h5ad}")
if sample_col not in a.obs.columns:
    print(f"warning: input h5ad missing obs[{sample_col!r}]; BLUE will fall back to the h5ad filename as the library id")
if hasattr(a, "file") and a.file is not None:
    a.file.close()
"""
    run_uv_python(
        blue_dir,
        uv_cmd,
        code,
        [str(h5ad_path), cell_type_col, sample_col],
        log_path=workdir / "logs" / "validate_reference_h5ad.log",
    )


def maybe_subset_reference_h5ad(
    blue_dir: Path,
    uv_cmd: list[str],
    reference_h5ad: Path,
    workdir: Path,
    cell_type_col: str,
    sample_col: str,
    mapping_path: Path,
    ref_sample_subset: bool,
    ref_sample_subset_n: int | None,
) -> tuple[Path, dict[str, Any] | None]:
    if not ref_sample_subset:
        return reference_h5ad, None
    if ref_sample_subset_n is None:
        raise ValueError("ref_sample_subset=true requires ref_sample_subset_n=<positive integer>")

    subset_path = workdir / "reference_subset.h5ad"
    report_path = workdir / "reference_subset_report.json"
    code = r"""
from pathlib import Path
import json
import sys

import anndata as ad

h5ad = Path(sys.argv[1])
subset_path = Path(sys.argv[2])
cell_type_col = sys.argv[3]
sample_col = sys.argv[4]
n_keep = int(sys.argv[5])
mapping_path = Path(sys.argv[6])
report_path = Path(sys.argv[7])

def parse_mapping_yaml(path):
    out = {}
    current = None
    for raw in path.read_text().splitlines():
        line = raw.split('#', 1)[0].rstrip()
        if not line.strip():
            continue
        if not line.startswith(' '):
            if not line.endswith(':'):
                raise SystemExit(f'Invalid mapping line in {path}: {raw}')
            current = line[:-1].strip().strip('\"\'')
            out[current] = []
            continue
        if current is None:
            continue
        stripped = line.strip()
        if stripped.startswith('-'):
            value = stripped[1:].strip().strip('\"\'')
            if value:
                out[current].append(value)
    if not out:
        raise SystemExit(f'No cell-type mapping entries found: {path}')
    return out

adata = ad.read_h5ad(h5ad, backed='r')
try:
    if cell_type_col not in adata.obs.columns:
        raise SystemExit(f'input h5ad missing obs[{cell_type_col!r}]: {h5ad}')
    if sample_col not in adata.obs.columns:
        raise SystemExit(f'ref_sample_subset=true requires obs[{sample_col!r}] in h5ad: {h5ad}')

    sample_values = adata.obs[sample_col].astype(str)
    unique_samples = list(dict.fromkeys(sample_values.tolist()))
    if n_keep > len(unique_samples):
        raise SystemExit(f'ref_sample_subset_n={n_keep} but only {len(unique_samples)} unique samples exist in obs[{sample_col!r}]')
    selected_samples = unique_samples[:n_keep]
    selected_set = set(selected_samples)
    mask = sample_values.isin(selected_set).to_numpy()
    if int(mask.sum()) == 0:
        raise SystemExit('Reference sample subset selected zero cells')

    selected_cell_types = sorted(set(adata.obs.loc[mask, cell_type_col].astype(str).tolist()))
    selected_fine = set(selected_cell_types)
    mapping = parse_mapping_yaml(mapping_path)
    missing_classes = {coarse: values for coarse, values in mapping.items() if not (selected_fine & set(values))}
    if missing_classes:
        formatted = '; '.join(f'{k}: {",".join(v)}' for k, v in missing_classes.items())
        raise SystemExit('Reference sample subset is missing mapped output classes: ' + formatted)

    subset = adata[mask, :].to_memory()
    subset_path.parent.mkdir(parents=True, exist_ok=True)
    subset.write_h5ad(subset_path)
finally:
    if hasattr(adata, 'file') and adata.file is not None:
        adata.file.close()

report = {
    'original_h5ad': str(h5ad),
    'subset_h5ad': str(subset_path),
    'sample_col': sample_col,
    'cell_type_col': cell_type_col,
    'requested_samples': n_keep,
    'available_samples': len(unique_samples),
    'selected_samples': selected_samples,
    'selected_sample_count': len(selected_samples),
    'selected_cell_count': int(mask.sum()),
    'selected_cell_types': selected_cell_types,
}
report_path.write_text(json.dumps(report, indent=2, sort_keys=True))
print(json.dumps(report, sort_keys=True))
"""
    result = run_uv_python(
        blue_dir,
        uv_cmd,
        code,
        [str(reference_h5ad), str(subset_path), cell_type_col, sample_col, str(ref_sample_subset_n), str(mapping_path), str(report_path)],
    )
    if result is not None and result.stdout:
        print(result.stdout.strip())
    report = json.loads(report_path.read_text())
    return subset_path, report


def run_logged(cmd: list[str], cwd: Path, log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print("  command:", " ".join(cmd))
    with log_path.open("w") as log:
        log.write("$ " + " ".join(cmd) + "\n\n")
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="")
            log.write(line)
        rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, cmd)


def run_uv_python(blue_dir: Path, uv_cmd: list[str], code: str, args: list[str], log_path: Path | None = None) -> subprocess.CompletedProcess[str] | None:
    cmd = [*uv_cmd, "run", "python", "-c", code, *args]
    if log_path is not None:
        run_logged(cmd, blue_dir, log_path)
        return None
    return subprocess.run(cmd, cwd=str(blue_dir), check=True, text=True, capture_output=True)


def load_blue_hyperparameters(repo_root: Path, hyperparameters_yaml: Any = None) -> tuple[dict[str, dict[str, Any]], Path]:
    if hyperparameters_yaml is None:
        hyper_path = repo_root / "DALE_Eval" / "configs" / "BLUE_hyperparameters.yaml"
    else:
        hyper_path = Path(str(hyperparameters_yaml)).expanduser()
        if not hyper_path.is_absolute():
            hyper_path = repo_root / hyper_path
    if not hyper_path.exists():
        raise FileNotFoundError(f"BLUE hyperparameters YAML not found: {hyper_path}")
    return parse_simple_blue_hyperparameters(hyper_path), hyper_path


def resolve_blue_celltype_mapping(
    repo_root: Path,
    reference_h5ad: Path,
    reference_mapping_dir: Path,
    blue_dir: Path,
    uv_cmd: list[str],
    workdir: Path,
    hyper: dict[str, Any],
    extra_mapping_yaml: Any = None,
) -> tuple[Path, str]:
    config_mapping = hyper.get("cell_types", {}).get("mapping_yaml")
    mapping_value = extra_mapping_yaml if extra_mapping_yaml is not None else config_mapping
    if mapping_value is None:
        mapping_path = reference_mapping_dir / "celltype_mapping.yaml"
        if not mapping_path.exists():
            raise FileNotFoundError(f"BLUE cell type mapping YAML not found: {mapping_path}")
        return mapping_path, "default"

    mapping_path = Path(str(mapping_value)).expanduser()
    if not mapping_path.is_absolute():
        mapping_path = repo_root / mapping_path
    if not mapping_path.exists():
        raise FileNotFoundError(f"BLUE cell type mapping YAML not found: {mapping_path}")
    return mapping_path, "extra_args" if extra_mapping_yaml is not None else "config"


def _blue_gene_path(gene_lists_dir: Path, gene_set: str) -> Path:
    if gene_set == "deg":
        return gene_lists_dir / "all_DEG_lst.txt"
    if gene_set == "common":
        return gene_lists_dir / "all_common_gene_lst.txt"
    raise ValueError(f"Unsupported BLUE gene set: {gene_set}")


def _normalize_gene_set(value: Any, name: str, allowed: set[str], allow_null: bool = False) -> str | None:
    if value is None:
        if allow_null:
            return None
        raise ValueError(f"BLUE genes.{name} cannot be null")
    value = str(value).lower()
    if value not in allowed:
        raise ValueError(f"BLUE genes.{name} must be one of {', '.join(sorted(allowed))}; got {value!r}")
    return value


def resolve_blue_gene_paths(hyper: dict[str, Any], gene_lists_dir: Path) -> dict[str, str | None]:
    genes = hyper.setdefault("genes", {})
    predict_gene_set = _normalize_gene_set(
        genes.get("predict_gene_set", "deg"),
        "predict_gene_set",
        {"deg", "common"},
    )
    input_gene_set = _normalize_gene_set(
        genes.get("input_gene_set", None),
        "input_gene_set",
        {"deg", "common"},
        allow_null=True,
    )
    pseudobulk_gene_set = _normalize_gene_set(
        genes.get("pseudobulk_gene_set", "auto"),
        "pseudobulk_gene_set",
        {"auto", "deg", "common"},
    )

    effective_input_gene_set = input_gene_set or "deg"
    if pseudobulk_gene_set == "auto":
        pseudobulk_gene_set = "common" if "common" in {effective_input_gene_set, predict_gene_set} else "deg"

    if pseudobulk_gene_set == "deg" and "common" in {effective_input_gene_set, predict_gene_set}:
        raise ValueError(
            "BLUE genes.pseudobulk_gene_set=deg cannot cover common input/output genes; "
            "use pseudobulk_gene_set=auto or common."
        )

    genes["predict_gene_set"] = predict_gene_set
    genes["input_gene_set"] = input_gene_set
    genes["pseudobulk_gene_set"] = pseudobulk_gene_set

    return {
        "input_gene_list": str(_blue_gene_path(gene_lists_dir, input_gene_set)) if input_gene_set else None,
        "output_gene_list": str(_blue_gene_path(gene_lists_dir, predict_gene_set)),
        "pseudobulk_gene_list": str(_blue_gene_path(gene_lists_dir, pseudobulk_gene_set)),
    }

def write_generated_blue_config(
    paths: DeconvPaths,
    reference_h5ad: Path,
    mapping_path: Path,
    workdir: Path,
    bulk_tsv: Path,
    hyper: dict[str, Any],
    cell_type_col: str,
    sample_col: str,
) -> Path:
    gene_lists_dir = workdir / "deconv_training_data"
    gene_paths = resolve_blue_gene_paths(hyper, gene_lists_dir)

    config = {
        "paths": {
            "sc_h5ad_dir": str(reference_h5ad.parent),
            "sc_h5ad_glob": reference_h5ad.name,
            "bulk_tsv": str(bulk_tsv),
            "bulk_metadata_tsv": None,
            "combined_sc_h5ad": str(workdir / "processed_data" / "combined_labeled.h5ad"),
            "celltype_mapping_yaml": str(mapping_path),
            "gene_lists_dir": str(gene_lists_dir),
            "pseudobulk_dir": str(gene_lists_dir / "pseudobulks"),
            "ckpt_root": str(workdir / "deconv_ckpt"),
            "predictions_dir": str(workdir / "deconv_predictions"),
            **gene_paths,
        },
        "sc_schema": {
            "fine_celltype_col": cell_type_col,
            "library_id_col": sample_col,
            "raw_counts_layer": None,
            "ensembl_var_names": False,
        },
        **hyper,
    }

    config_path = workdir / "pipeline_config_BLUE.yaml"
    write_simple_yaml(config_path, config)
    return config_path



def generate_blue_loss_plot(blue_dir: Path, uv_cmd: list[str], workdir: Path, run_stamp: str) -> Path | None:
    loss_csv = workdir / "deconv_ckpt" / run_stamp / "loss_history.csv"
    loss_png = loss_csv.with_suffix(".png")
    if not loss_csv.exists():
        print(f"  loss_history: not found at {loss_csv}")
        return None

    code = r"""
from pathlib import Path
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

loss_csv = Path(sys.argv[1])
loss_png = Path(sys.argv[2])
df = pd.read_csv(loss_csv)
x = df['epoch'] if 'epoch' in df.columns else range(1, len(df) + 1)
fig, ax = plt.subplots(figsize=(7, 4.5))
for col in ['train', 'val_total', 'val_prop', 'val_gep']:
    if col in df.columns:
        ax.plot(x, df[col], linewidth=1.5, marker='o', markersize=3, label=col)
ax.set_xlabel('epoch')
ax.set_ylabel('loss')
ax.set_title('BLUE training convergence')
ax.grid(alpha=0.25, linewidth=0.6)
if ax.lines:
    ax.legend(frameon=False)
fig.tight_layout()
loss_png.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(loss_png, dpi=160)
plt.close(fig)
"""
    print(f"  loss_history: {loss_csv}")
    run_uv_python(
        blue_dir,
        uv_cmd,
        code,
        [str(loss_csv), str(loss_png)],
        workdir / "logs" / "plot_loss_history.log",
    )
    print(f"  loss_plot: {loss_png}")
    return loss_png


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


def export_blue_outputs(
    blue_dir: Path,
    uv_cmd: list[str],
    predictions_dir: Path,
    output_dir: Path,
    mapping: dict[str, str],
    workdir: Path,
) -> None:
    runs = sorted([p for p in predictions_dir.iterdir() if p.is_dir()], key=lambda p: p.stat().st_mtime)
    if not runs:
        raise FileNotFoundError(f"No BLUE prediction run directories found under {predictions_dir}")
    pred_run = runs[-1]
    h5ads = sorted(
        pred_run.glob("predicted_ctGEP_ep*.h5ad"),
        key=lambda p: int(p.stem.replace("predicted_ctGEP_ep", "")),
    )
    props = sorted(
        pred_run.glob("predicted_proportions_ep*.csv"),
        key=lambda p: int(p.stem.replace("predicted_proportions_ep", "")),
    )
    if not h5ads or not props:
        raise FileNotFoundError(f"Missing BLUE prediction h5ad/csv under {pred_run}")

    method_dir = output_dir / METHOD
    if method_dir.exists() and any(method_dir.iterdir()):
        print(f"  overwrite note: existing {METHOD} outputs found for this config; files in {method_dir} may be overwritten")
    method_dir.mkdir(parents=True, exist_ok=True)

    mapping_json = workdir / "export_cell_type_mapping.json"
    mapping_json.write_text(json.dumps(mapping, indent=2, sort_keys=True))
    code = r"""
from pathlib import Path
import json
import sys

import anndata as ad
import numpy as np
import pandas as pd

h5ad_path = Path(sys.argv[1])
prop_path = Path(sys.argv[2])
out_dir = Path(sys.argv[3])
mapping = json.loads(Path(sys.argv[4]).read_text())

adata = ad.read_h5ad(h5ad_path)
cell_types = list(adata.layers.keys())
mapped = [mapping.get(ct, ct) for ct in cell_types]
dupes = sorted({ct for ct in mapped if mapped.count(ct) > 1})
if dupes:
    raise SystemExit("Cell type mapping creates duplicated BLUE output names: " + ", ".join(dupes))

genes = [str(x) for x in adata.var_names]
samples = [str(x) for x in adata.obs_names]
for source_ct, out_ct in zip(cell_types, mapped):
    safe_name = out_ct.replace("/", "_") + ".txt.gz"
    matrix = np.asarray(adata.layers[source_ct]).T
    df = pd.DataFrame(matrix, index=genes, columns=samples)
    df.round(2).to_csv(out_dir / safe_name, sep="\t", compression="gzip", index_label="")

frac = pd.read_csv(prop_path, index_col=0)
missing = [ct for ct in cell_types if ct not in frac.columns]
if missing:
    raise SystemExit("Predicted fraction file missing columns: " + ", ".join(missing))
frac = frac.loc[:, cell_types].copy()
frac.columns = mapped
frac.round(6).to_csv(out_dir / "BLUEfrac.txt", sep="\t", index_label="")

print(f"exported {len(cell_types)} CTSE files to {out_dir}")
print(f"exported fraction file to {out_dir / 'BLUEfrac.txt'}")
"""
    run_uv_python(
        blue_dir,
        uv_cmd,
        code,
        [str(h5ads[-1]), str(props[-1]), str(method_dir), str(mapping_json)],
        log_path=workdir / "logs" / "export_outputs.log",
    )


def main() -> None:
    args = parse_args()
    repo_root = find_repo_root(RUN_DIR)
    extra = apply_method_default_extra_args(parse_extra_args(args.extra_args), "BLUE", repo_root)
    validate_extra_args(
        extra,
        {
            "map_cell_types",
            "use_test_samples",
            "cell_type_col",
            "sample_col",
            "uv_path",
            "celltype_mapping_yaml",
            "hyperparameters_yaml",
            "samplenum_per_ct",
            "val_samplenum_per_patient",
            "keep_temp",
            "ref_sample_subset",
            "ref_sample_subset_n",
        },
        METHOD,
    )

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

    blue_dir = repo_root / "DALE_Eval" / "external_modules" / "BLUE"
    paths = resolve_deconv_paths(dataset, config_id, repo_root)
    if paths.ref_type == "self":
        raise ValueError(
            "BLUE does not support refType=self in this benchmark runner. "
            "BLUE trains from scRNA-derived pseudobulks, and self-reference splits may leave too few cells for stable training. "
            "Use an independent reference config instead."
        )
    uv_cmd = resolve_uv_command(extra_arg(extra, "uv_path", None))
    ref_name = paths.ref_dir.name
    reference_h5ad = find_reference_h5ad(repo_root, ref_name)
    original_reference_h5ad = reference_h5ad

    logs_dir = paths.obj_dir / "logs"
    workdir = logs_dir / f"BLUE_temp_{config_id}_{run_stamp}"
    workdir.mkdir(parents=True, exist_ok=True)
    validate_reference_h5ad_columns(blue_dir, uv_cmd, reference_h5ad, workdir, cell_type_col, sample_col)
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

    bulk_tsv = workdir / "bulk_input_BLUE.tsv"
    n_genes, bulk_samples, bulk_prep = write_prepared_bulk_tsv_for_deconv(
        paths.bulk_path,
        bulk_tsv,
        paths,
        selected_samples=selected_samples,
        method=METHOD,
    )
    deconv_input = format_deconv_dim(n_genes, len(bulk_samples))

    hyper, hyper_path = load_blue_hyperparameters(repo_root, extra_arg(extra, "hyperparameters_yaml", None))
    pseudobulk_hyper = hyper.setdefault("pseudobulk", {})
    for key in ("samplenum_per_ct", "val_samplenum_per_patient"):
        value = extra_arg(extra, key, None)
        if value is not None:
            pseudobulk_hyper[key] = positive_int_extra(value, key)
    mapping_path, mapping_source = resolve_blue_celltype_mapping(
        repo_root=repo_root,
        reference_h5ad=reference_h5ad,
        reference_mapping_dir=paths.ref_dir,
        blue_dir=blue_dir,
        uv_cmd=uv_cmd,
        workdir=workdir,
        hyper=hyper,
        extra_mapping_yaml=extra_arg(extra, "celltype_mapping_yaml", None),
    )
    reference_h5ad, subset_report = maybe_subset_reference_h5ad(
        blue_dir=blue_dir,
        uv_cmd=uv_cmd,
        reference_h5ad=reference_h5ad,
        workdir=workdir,
        cell_type_col=cell_type_col,
        sample_col=sample_col,
        mapping_path=mapping_path,
        ref_sample_subset=ref_sample_subset,
        ref_sample_subset_n=ref_sample_subset_n,
    )

    blue_config = write_generated_blue_config(
        paths=paths,
        reference_h5ad=reference_h5ad,
        mapping_path=mapping_path,
        workdir=workdir,
        bulk_tsv=bulk_tsv,
        hyper=hyper,
        cell_type_col=cell_type_col,
        sample_col=sample_col,
    )

    bulk_input = paths.config["bulk_input"]
    frac_input = paths.config["frac_input"]
    print("Running BLUE")
    print(f"  dataset: {dataset}")
    print(f"  config_id: {config_id}")
    print(f"  bulk_input: {bulk_input} ({deconv_input})")
    print(f"  bulk_scale: {paths.config['bulk_scale']}")
    print(f"  bulk_normalization: {paths.config['bulk_normalization']}")
    print(f"  bulk_preparation: {bulk_prep['action']}")
    print(f"  config_frac_input: {frac_input}")
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
    print(f"  BLUE_mapping: {mapping_path}")
    print(f"  BLUE_mapping_source: {mapping_source}")
    print(f"  BLUE_hyperparameters: {hyper_path}")
    print(f"  BLUE_samplenum_per_ct: {pseudobulk_hyper.get('samplenum_per_ct')}")
    print(f"  BLUE_val_samplenum_per_patient: {pseudobulk_hyper.get('val_samplenum_per_patient')}")
    print(f"  workdir: {workdir}")
    print(f"  map_cell_types: {map_cell_types}")
    print(f"  use_test_samples: {use_test_samples}")
    print(f"  keep_temp: {keep_temp}")
    print("  note: BLUE does not consume frac_input; it estimates fractions from bulk + scRNA reference")
    print("Starting deconvolution")
    print(f"  deconv_input: {deconv_input}")

    commands = [
        ("00_combine_sc_h5ads", [*uv_cmd, "run", "python", "scripts/deconv/00_combine_sc_h5ads.py", "--config", str(blue_config)]),
        ("01_build_gene_lists", [*uv_cmd, "run", "python", "scripts/deconv/01_build_gene_lists.py", "--config", str(blue_config)]),
        ("02_simulate_pseudobulk", [*uv_cmd, "run", "python", "scripts/deconv/02_simulate_pseudobulk.py", "--config", str(blue_config)]),
        ("03_preprocess_bulk", [*uv_cmd, "run", "python", "scripts/deconv/03_preprocess_bulk.py", "--config", str(blue_config)]),
    ]
    train_cmd = [*uv_cmd, "run", "python", "scripts/deconv/04_train_unet.py", "--config", str(blue_config), "--run-stamp", run_stamp]
    predict_cmd = [*uv_cmd, "run", "python", "scripts/deconv/05_predict_bulk.py", "--config", str(blue_config)]
    commands.extend([
        ("04_train_unet", train_cmd),
        ("05_predict_bulk", predict_cmd),
    ])

    runtime_entry = start_runtime_log(paths, METHOD, n_core=n_core, deconv_input=deconv_input, log_path=str(message_log_path))
    loss_plot_path: Path | None = None
    try:
        for step_name, cmd in commands:
            print(f"\nBLUE step: {step_name}")
            run_logged(cmd, blue_dir, workdir / "logs" / f"{step_name}.log")
            if step_name == "04_train_unet":
                loss_plot_path = generate_blue_loss_plot(blue_dir, uv_cmd, workdir, run_stamp)

        export_mapping = read_indep_ref_cell_type_mapping(paths, repo_root) if map_cell_types else {}
        export_blue_outputs(
            blue_dir=blue_dir,
            uv_cmd=uv_cmd,
            predictions_dir=workdir / "deconv_predictions",
            output_dir=paths.output_dir,
            mapping=export_mapping,
            workdir=workdir,
        )
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

    print(f"  BLUE ctse: {paths.output_dir / METHOD / '<cell_type>.txt.gz'}")
    print(f"  BLUE fractions: {paths.output_dir / METHOD / 'BLUEfrac.txt'}")
    if keep_temp:
        print(f"  BLUE workdir: {workdir}")


if __name__ == "__main__":
    main()
