from __future__ import annotations

import csv
import gzip
import math
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any


def find_repo_root(start: str | Path | None = None) -> Path:
    cur = Path(start or Path.cwd()).resolve()
    if cur.is_file():
        cur = cur.parent
    while True:
        if (
            (cur / "DALE_Eval" / "configs" / "deconv_configs.txt").is_file()
            and (cur / "Benchmarking_obj").is_dir()
            and (cur / "DALE_Eval").is_dir()
        ):
            return cur
        if cur.parent == cur:
            raise FileNotFoundError(f"Could not find repo root from: {start or Path.cwd()}")
        cur = cur.parent


def read_tsv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def config_bulk_normalization(config: dict[str, str]) -> str:
    return config["bulk_normalization"].strip().lower()


def config_slug(config: dict[str, str]) -> str:
    norm = config_bulk_normalization(config)
    norm_suffix = "" if norm == "cpm" else f"__norm-{norm}"
    return (
        f"{config['config_id']}_bulk-{config['bulk_input']}"
        f"__frac-{config['frac_input']}__ref-{config['refType']}"
        f"{norm_suffix}"
    )


def read_deconv_config(config_id: str, repo_root: Path) -> dict[str, str]:
    config_path = repo_root / "DALE_Eval" / "configs" / "deconv_configs.txt"
    rows = read_tsv_dicts(config_path)
    hits = [row for row in rows if row.get("config_id") == config_id]
    if len(hits) != 1:
        raise ValueError(f"Expected exactly one deconv config for config_id={config_id}, found {len(hits)}")
    required = {"config_id", "bulk_input", "bulk_scale", "bulk_normalization", "frac_input", "refType"}
    missing = required.difference(hits[0])
    if missing:
        raise ValueError(f"deconv config missing columns: {', '.join(sorted(missing))}")
    return hits[0]


def resolve_indep_ref_name(dataset: str, repo_root: Path) -> str:
    assignment_path = repo_root / "DALE_Eval" / "configs" / "benchmark_ref_assignment.txt"
    rows = read_tsv_dicts(assignment_path)
    hits = [row for row in rows if row.get("dataset") == dataset]
    if len(hits) != 1:
        raise ValueError(f"Expected exactly one independent reference assignment for dataset={dataset}")
    return hits[0]["indep_ref"]


@dataclass
class DeconvPaths:
    dataset: str
    config_id: str
    config: dict[str, str]
    ref_type: str
    obj_dir: Path
    ref_dir: Path
    bulk_path: Path
    frac_path: Path
    output_dir: Path
    runtime_path: Path


def resolve_deconv_paths(dataset: str, config_id: str, repo_root: Path) -> DeconvPaths:
    config = read_deconv_config(config_id, repo_root)
    ref_type = config["refType"]
    if ref_type not in {"indep", "self"}:
        raise ValueError(f"Invalid refType in deconv config {config_id}: {ref_type}")

    obj_dir = repo_root / "Benchmarking_obj" / dataset
    if ref_type == "self":
        ref_dir = obj_dir / "self_reference"
    else:
        ref_dir = repo_root / "Indep_scReference" / resolve_indep_ref_name(dataset, repo_root)

    return DeconvPaths(
        dataset=dataset,
        config_id=config_id,
        config=config,
        ref_type=ref_type,
        obj_dir=obj_dir,
        ref_dir=ref_dir,
        bulk_path=obj_dir / "bulk_input" / f"{config['bulk_input']}.txt",
        frac_path=obj_dir / "frac_input" / f"{config['frac_input']}.txt",
        output_dir=obj_dir / "deconv_res" / config_slug(config),
        runtime_path=obj_dir / "logs" / "deconv_runs.txt",
    )


def parse_extra_args(extra_args: str | None) -> dict[str, Any]:
    if not extra_args:
        return {}
    out: dict[str, Any] = {}
    for piece in extra_args.split(";"):
        piece = piece.strip()
        if not piece:
            continue
        if "=" not in piece:
            raise ValueError(f"Invalid extra_args piece: {piece}. Expected key=value")
        key, value = [x.strip() for x in piece.split("=", 1)]
        lower = value.lower()
        if lower in {"true", "false"}:
            parsed: Any = lower == "true"
        else:
            try:
                parsed = int(value)
            except ValueError:
                try:
                    parsed = float(value)
                except ValueError:
                    parsed = value
        out[key] = parsed
    return out



def _parse_method_default_value(value: str, value_type: str) -> Any:
    value_type = value_type.lower()
    if value_type == "null":
        return None
    if value_type == "logical":
        lower = value.lower()
        if lower not in {"true", "false"}:
            raise ValueError(f"Invalid logical default value: {value}")
        return lower == "true"
    if value_type == "integer":
        return int(value)
    if value_type == "numeric":
        return float(value)
    if value_type == "character":
        return value
    raise ValueError(f"Unsupported method default value_type: {value_type}")


def read_method_default_extra_configs(repo_root: Path) -> list[dict[str, str]]:
    path = repo_root / "DALE_Eval" / "configs" / "method_default_extra_configs.txt"
    if not path.exists():
        raise FileNotFoundError(f"Method default extra config file not found: {path}")
    rows = read_tsv_dicts(path)
    required = {"method", "extra_arg", "default_value", "value_type"}
    for row in rows:
        missing = required.difference(row)
        if missing:
            raise ValueError(f"method_default_extra_configs.txt missing columns: {', '.join(sorted(missing))}")
    return rows


def apply_method_default_extra_args(extra: dict[str, Any], method: str, repo_root: Path) -> dict[str, Any]:
    merged: dict[str, Any] = {}
    for row in read_method_default_extra_configs(repo_root):
        if row["method"] != method:
            continue
        merged[row["extra_arg"]] = _parse_method_default_value(row["default_value"], row["value_type"])
    merged.update(extra)
    return merged


def validate_extra_args(extra: dict[str, Any], allowed: set[str], method: str) -> None:
    unknown = sorted(set(extra).difference(allowed))
    if unknown:
        allowed_msg = ", ".join(sorted(allowed)) if allowed else "none"
        raise ValueError(
            f"{method} received unsupported extra_args: {', '.join(unknown)}. "
            f"Allowed extra_args: {allowed_msg}"
        )


def extra_arg(extra: dict[str, Any], key: str, default: Any = None) -> Any:
    return extra[key] if key in extra else default


def read_test_samples(paths: DeconvPaths) -> list[str]:
    sample_split_path = paths.obj_dir / "self_reference" / "sample_split.txt"
    if not sample_split_path.exists():
        raise FileNotFoundError(f"use_test_samples=true but sample split file not found: {sample_split_path}")
    rows = read_tsv_dicts(sample_split_path)
    return [
        row["sampleIDs"]
        for row in rows
        if row.get("group", "").lower() == "test" and row.get("sampleIDs")
    ]


def read_bulk_header(path: Path) -> tuple[str, list[str]]:
    with path.open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
    if len(header) < 2:
        raise ValueError(f"Bulk file has fewer than 2 columns: {path}")
    return header[0], header[1:]


def read_bulk_dataframe(path: Path, selected_samples: list[str] | None = None):
    import pandas as pd

    df = pd.read_csv(path, sep="	", index_col=0)
    if selected_samples is not None:
        missing = [sample for sample in selected_samples if sample not in df.columns]
        if missing:
            raise ValueError("Selected bulk samples missing from bulk file: " + ", ".join(missing))
        df = df.loc[:, selected_samples]
    return df


def _counts_to_cpm_dataframe(df):
    col_sums = df.sum(axis=0)
    if (col_sums < 0).any():
        raise ValueError("Count bulk input contains negative library sizes")
    out = df.copy()
    positive = col_sums > 0
    out.loc[:, ~positive] = 0
    out.loc[:, positive] = out.loc[:, positive].div(col_sums[positive], axis=1) * 1_000_000.0
    return out


def prepare_bulk_dataframe_for_deconv(df, paths: DeconvPaths, method: str):
    bulk_scale, bulk_normalization, action = _python_bulk_prep_config(paths, method)

    out = df
    if bulk_scale == "counts" and bulk_normalization == "cpm":
        out = _counts_to_cpm_dataframe(df)

    return {
        "bulk_expr": out,
        "bulk_scale": bulk_scale,
        "bulk_normalization": bulk_normalization,
        "action": action,
        "note": (
            f"{method} bulk preparation; config_id={paths.config_id} "
            f"bulk_input={paths.config['bulk_input']} bulk_scale={bulk_scale} "
            f"bulk_normalization={bulk_normalization}; {action}"
        ),
    }


def _python_bulk_prep_config(paths: DeconvPaths, method: str) -> tuple[str, str, str]:
    bulk_scale = paths.config["bulk_scale"].strip().lower()
    bulk_normalization = paths.config["bulk_normalization"].strip().lower()
    if bulk_scale not in {"counts", "cpm"}:
        raise ValueError(f"Unsupported bulk_scale={bulk_scale} for config_id={paths.config_id}")
    if bulk_normalization not in {"cpm", "tmm", "uq", "none"}:
        raise ValueError(f"Unsupported bulk_normalization={bulk_normalization} for config_id={paths.config_id}")
    if bulk_normalization == "none":
        return bulk_scale, bulk_normalization, "used input without runner normalization"
    if bulk_normalization == "cpm":
        if bulk_scale == "counts":
            return bulk_scale, bulk_normalization, "converted counts to CPM"
        return bulk_scale, bulk_normalization, "used CPM input without runner normalization"
    raise ValueError(
        f"{method} does not yet support bulk_normalization={bulk_normalization}. "
        "edgeR-based normalization is currently supported only by the R runners."
    )


def write_prepared_bulk_tsv_for_deconv(
    source_bulk: Path,
    out_path: Path,
    paths: DeconvPaths,
    selected_samples: list[str] | None,
    method: str,
) -> tuple[int, list[str], dict[str, str]]:
    bulk_scale, bulk_normalization, action = _python_bulk_prep_config(paths, method)
    first_col, samples = read_bulk_header(source_bulk)
    del first_col
    if selected_samples is None:
        selected_samples = samples
    missing = [sample for sample in selected_samples if sample not in samples]
    if missing:
        raise ValueError(f"Selected bulk samples missing from bulk file: {', '.join(missing)}")
    keep_idx = [samples.index(sample) for sample in selected_samples]

    lib_sizes = [0.0 for _ in selected_samples]
    if bulk_scale == "counts" and bulk_normalization == "cpm":
        with source_bulk.open(newline="") as src:
            reader = csv.reader(src, delimiter="	")
            next(reader)
            for row in reader:
                if not row:
                    continue
                for out_i, src_i in enumerate(keep_idx):
                    lib_sizes[out_i] += float(row[src_i + 1] or 0)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_genes = 0
    with source_bulk.open(newline="") as src, out_path.open("w", newline="") as dst:
        reader = csv.reader(src, delimiter="	")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        next(reader)
        writer.writerow(["gene_id", *selected_samples])
        for row in reader:
            if not row:
                continue
            gene = row[0]
            values: list[str] = []
            for out_i, src_i in enumerate(keep_idx):
                value = float(row[src_i + 1] or 0)
                if bulk_scale == "counts" and bulk_normalization == "cpm":
                    lib_size = lib_sizes[out_i]
                    value = 0.0 if lib_size <= 0 else value / lib_size * 1_000_000.0
                values.append(f"{value:.6g}")
            writer.writerow([gene, *values])
            n_genes += 1

    prep = {
        "bulk_scale": bulk_scale,
        "bulk_normalization": bulk_normalization,
        "action": action,
        "note": (
            f"{method} bulk preparation; config_id={paths.config_id} "
            f"bulk_input={paths.config['bulk_input']} bulk_scale={bulk_scale} "
            f"bulk_normalization={bulk_normalization}; {action}"
        ),
    }
    return n_genes, selected_samples, prep


def format_deconv_dim(n_genes: int, n_samples: int) -> str:
    return f"{n_genes} genes x {n_samples} samples"


RUNTIME_COLUMNS = [
    "run_id",
    "dataset",
    "method",
    "deconv_config",
    "bulk_input",
    "frac_input",
    "refType",
    "start_time",
    "n_core",
    "deconv_input",
    "status",
    "run_time",
]


def method_message_log_path(paths: DeconvPaths, method: str, stamp: str | None = None) -> Path:
    stamp = stamp or datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_method = "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in method)
    safe_config = "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in paths.config_id)
    return paths.obj_dir / "logs" / "run_logs" / f"{safe_method}_{safe_config}_{stamp}.log"


def make_run_id(method: str, config_id: str, stamp: str | None = None) -> str:
    stamp = stamp or datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_method = "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in method)
    safe_config = "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in config_id)
    return f"{safe_method}_{safe_config}_{stamp}"


def run_id_from_log_path(log_path: str, method: str, config_id: str) -> str:
    if log_path:
        stem = Path(log_path).stem
        safe_method = "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in method)
        safe_config = "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in config_id)
        if stem.startswith(f"{safe_method}_{safe_config}_"):
            return stem
    return make_run_id(method, config_id)


METHODS_WITHOUT_FRAC_INPUT = {"BLUE", "scTAPE", "InstaPrism", "InstaPrismUpdated"}
METHODS_WITHOUT_N_CORE = {"BLUE", "scTAPE", "CIBERSORTx", "ENIGMAL2", "ENIGMAtrace"}


def runtime_frac_input(paths: DeconvPaths, method: str) -> str:
    if method in METHODS_WITHOUT_FRAC_INPUT:
        return "NA"
    return paths.config["frac_input"]


def runtime_n_core(method: str, n_core: int | None) -> str:
    if method in METHODS_WITHOUT_N_CORE or n_core is None:
        return "NA"
    return str(n_core)

def _runtime_timestamp() -> str:
    return datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S.%f")[:-3] + " " + time.tzname[0]


def _elapsed_mins(start_clock: float) -> str:
    return f"{round((time.time() - start_clock) / 60, 2)}mins"


def read_runtime_log(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = []
        for row in reader:
            for col in RUNTIME_COLUMNS:
                row.setdefault(col, "")
            rows.append({col: row[col] for col in RUNTIME_COLUMNS})
        return rows


def write_runtime_log(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=RUNTIME_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


@dataclass
class RuntimeEntry:
    run_id: str
    method: str
    start_time: str
    start_clock: float
    n_core: str
    deconv_input: str


def start_runtime_log(paths: DeconvPaths, method: str, n_core: int | None, deconv_input: str, log_path: str = "") -> RuntimeEntry:
    entry = RuntimeEntry(
        run_id=run_id_from_log_path(log_path, method, paths.config_id),
        method=method,
        start_time=_runtime_timestamp(),
        start_clock=time.time(),
        n_core=runtime_n_core(method, n_core),
        deconv_input=deconv_input,
    )
    rows = read_runtime_log(paths.runtime_path)
    rows.append(
        {
            "run_id": entry.run_id,
            "dataset": paths.dataset,
            "method": method,
            "deconv_config": paths.config_id,
            "bulk_input": paths.config["bulk_input"],
            "frac_input": runtime_frac_input(paths, method),
            "refType": paths.ref_type,
            "start_time": entry.start_time,
            "n_core": entry.n_core,
            "deconv_input": deconv_input,
            "status": "running",
            "run_time": "",
        }
    )
    write_runtime_log(paths.runtime_path, rows)
    return entry


def finish_runtime_log(paths: DeconvPaths, entry: RuntimeEntry, status: str) -> None:
    rows = read_runtime_log(paths.runtime_path)
    run_time = _elapsed_mins(entry.start_clock)
    for row in reversed(rows):
        if row["run_id"] == entry.run_id:
            row["status"] = status
            row["run_time"] = run_time
            break
    else:
        rows.append(
            {
                "run_id": entry.run_id,
                "dataset": paths.dataset,
                "method": entry.method,
                "deconv_config": paths.config_id,
                "bulk_input": paths.config["bulk_input"],
                "frac_input": runtime_frac_input(paths, entry.method),
                "refType": paths.ref_type,
                "start_time": entry.start_time,
                "n_core": entry.n_core,
                "deconv_input": entry.deconv_input,
                "status": status,
                "run_time": run_time,
            }
        )
    write_runtime_log(paths.runtime_path, rows)


def read_indep_ref_cell_type_mapping(paths: DeconvPaths, repo_root: Path) -> dict[str, str]:
    mapping_path = repo_root / "Indep_scReference" / "cell_type_mapping.txt"
    if not mapping_path.exists():
        return {}
    rows = read_tsv_dicts(mapping_path)
    ref_name = paths.ref_dir.name
    mapping: dict[str, str] = {}
    for row in rows:
        if row.get("dataset") == paths.dataset and row.get("indep_ref") == ref_name:
            mapping[row["indep_ref_cell_type"]] = row["target_cell_type"]
    return mapping


def _safe_float(value: str) -> float:
    try:
        out = float(value)
    except ValueError:
        return 0.0
    return out if math.isfinite(out) else 0.0


def gzip_write_matrix(path: Path, row_names: list[str], col_names: list[str], values: list[list[float]], digits: int = 2) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["", *col_names])
        for row_name, row_values in zip(row_names, values):
            writer.writerow([row_name, *[round(v, digits) if v != 0 else 0 for v in row_values]])


def simple_yaml_value(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, list):
        return "[" + ", ".join(simple_yaml_value(v) for v in value) + "]"
    text = str(value)
    if not text or any(ch in text for ch in ":#[]{}&,*!|>'\"%@`"):
        return repr(text)
    return text


def write_simple_yaml(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for key, value in data.items():
        if isinstance(value, dict):
            lines.append(f"{key}:")
            for subkey, subvalue in value.items():
                lines.append(f"  {subkey}: {simple_yaml_value(subvalue)}")
        else:
            lines.append(f"{key}: {simple_yaml_value(value)}")
    path.write_text("\n".join(lines) + "\n")


def parse_simple_blue_hyperparameters(path: Path) -> dict[str, dict[str, Any]]:
    """Parse the benchmark-owned BLUE defaults YAML.

    This intentionally supports only the simple two-level structure used by
    DALE_Eval/configs/BLUE_hyperparameters.yaml so run_BLUE.py can be launched
    with the system Python even when PyYAML is only installed in BLUE's uv env.
    """
    out: dict[str, dict[str, Any]] = {}
    current: str | None = None
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        if not line.startswith(" "):
            key = line.rstrip(":")
            out[key] = {}
            current = key
            continue
        if current is None or ":" not in line:
            continue
        key, value = [x.strip() for x in line.split(":", 1)]
        out[current][key] = parse_simple_scalar(value)
    return out


def parse_simple_scalar(value: str) -> Any:
    if value == "null":
        return None
    if value in {"true", "false"}:
        return value == "true"
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1].strip()
        if not inner:
            return []
        return [parse_simple_scalar(piece.strip()) for piece in inner.split(",")]
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value
