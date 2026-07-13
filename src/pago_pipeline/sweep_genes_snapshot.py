from __future__ import annotations

import importlib
import shutil
import sys
import tempfile
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, TypeAlias

import numpy as np
import pandas as pd
import psutil
from Bio import SeqIO

from src.pago_pipeline.ncbi_fasta_snapshot import (
    load_fasta_snapshot_by_directory,
    load_latest_fasta_snapshot,
)
from src.pago_pipeline.ncbi_snapshot import (
    SnapshotMode,
    _coerce_snapshot_mode,
    _replace_latest_directory,
    build_snapshot_directory_name,
    get_most_recent_snapshot_directory,
    list_saved_snapshot_directories,
)
from src.pago_pipeline.storage import (
    read_json_file,
    sha256_of_file,
    write_csv_rows_atomic,
    write_json_atomic,
)

PathLike: TypeAlias = str | Path

DEFAULT_SWEEP_GENES_MASKS: tuple[tuple[int, ...], ...] = (
    (1, 1, 0, 0, 1),
    (1, 1, 1, 0, 0),
    (1, 0, 1, 0, 1),
    (1, 1, 0, 1, 0),
)
DEFAULT_PROJECTED_DIMENSIONS_PER_MASK = 700
DEFAULT_SEQUENCE_METADATA_FILE_NAME = "sequence_metadata.csv"
DEFAULT_PROFILING_LOG_FILE_NAME = "profiling_log.csv"


@dataclass(frozen=True)
class SweepGenesSnapshotResult:
    snapshot_directory: Path
    manifest_file_path: Path
    embeddings_file_path: Path
    sequence_metadata_file_path: Path
    profiling_log_file_path: Path
    source_fasta_snapshot_directory: Path


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _current_utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _build_embeddings_file_name(total_embedding_dimensions: int) -> str:
    if total_embedding_dimensions <= 0:
        raise ValueError("total_embedding_dimensions must be a positive integer.")
    return f"sweep_genes_embeddings_{total_embedding_dimensions}D.npy"


def _get_process_cpu_seconds(process: psutil.Process) -> float:
    cpu_times = process.cpu_times()
    return float(cpu_times.user + cpu_times.system)


def _get_process_rss_memory_mb(process: psutil.Process) -> float:
    return float(process.memory_info().rss / (1024**2))


def _finalize_profiling_row(
    *,
    row_payload: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        key: (
            round(value, 6)
            if isinstance(value, float)
            else value
        )
        for key, value in row_payload.items()
    }


def _normalize_masks(
    masks: Sequence[Sequence[int]],
) -> list[list[int]]:
    if not masks:
        raise ValueError("masks must contain at least one binary mask.")

    normalized_masks = [list(mask) for mask in masks]
    mask_lengths = {len(mask) for mask in normalized_masks}
    if len(mask_lengths) != 1:
        raise ValueError("All masks must have the same length.")

    for mask in normalized_masks:
        if not mask:
            raise ValueError("Masks must not be empty.")
        if any(mask_value not in {0, 1} for mask_value in mask):
            raise ValueError("Masks must contain only 0 and 1 values.")
        if not any(mask):
            raise ValueError("Each mask must keep at least one position.")

    return normalized_masks


def _load_vendored_sweep_package(
    *,
    scripts_sweep_root_directory: PathLike,
) -> Any:
    resolved_scripts_sweep_root_directory = _as_path(scripts_sweep_root_directory)
    if not resolved_scripts_sweep_root_directory.exists():
        raise FileNotFoundError(
            "Vendored sweep package root directory was not found: "
            f"{resolved_scripts_sweep_root_directory}."
        )
    if not (resolved_scripts_sweep_root_directory / "sweep").exists():
        raise FileNotFoundError(
            "Vendored sweep package directory was not found: "
            f"{resolved_scripts_sweep_root_directory / 'sweep'}."
        )

    resolved_package_root = str(resolved_scripts_sweep_root_directory.resolve())
    if resolved_package_root not in sys.path:
        sys.path.insert(0, resolved_package_root)

    return importlib.import_module("sweep")


def _write_sequence_metadata_csv(
    *,
    sequence_records: Sequence[Any],
    output_file_path: Path,
) -> Path:
    output_file_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_dataframe = pd.DataFrame(
        {
            "sequence_index": np.arange(len(sequence_records), dtype=np.int64),
            "record_id": [record.id for record in sequence_records],
            "description": [record.description for record in sequence_records],
            "sequence_length": [len(record.seq) for record in sequence_records],
        }
    )

    with tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        dir=output_file_path.parent,
        encoding="utf-8",
        newline="",
        suffix=".csv",
    ) as temporary_file:
        resolved_temporary_file_path = Path(temporary_file.name)

    metadata_dataframe.to_csv(resolved_temporary_file_path, index=False)
    resolved_temporary_file_path.replace(output_file_path)
    return output_file_path


def _write_profiling_log_csv(
    *,
    profiling_rows: Sequence[Mapping[str, Any]],
    output_file_path: Path,
) -> Path:
    output_file_path.parent.mkdir(parents=True, exist_ok=True)
    return write_csv_rows_atomic(
        rows=[dict(row) for row in profiling_rows],
        output_file_path=output_file_path,
    )


def _build_sweep_genes_snapshot_manifest(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    embeddings_file_path: Path,
    sequence_metadata_file_path: Path,
    profiling_log_file_path: Path,
    sequence_count: int,
    masks: Sequence[Sequence[int]],
    projected_dimensions_per_mask: int,
    composition: str,
    projection: bool,
    random_seed: int,
    chunk_size: int,
    n_jobs: int,
    scripts_sweep_root_directory: Path,
    source_fasta_snapshot_payload: Mapping[str, object],
    profiling_rows: Sequence[Mapping[str, Any]],
) -> dict[str, object]:
    source_fasta_manifest = source_fasta_snapshot_payload.get("manifest")
    if not isinstance(source_fasta_manifest, Mapping):
        raise RuntimeError(
            "Resolved FASTA snapshot payload is missing a valid source manifest."
        )

    total_embedding_dimensions = len(masks) * projected_dimensions_per_mask
    total_profiling_rows = len(profiling_rows)
    snapshot_total_row = next(
        (
            row
            for row in profiling_rows
            if row.get("work_unit_kind") == "snapshot_total"
        ),
        None,
    )
    return {
        "snapshot_format_version": "1.1",
        "artifact_type": "sweep_genes_embedding_snapshot",
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "manifest_file_name": "manifest.json",
        "embeddings_file_name": embeddings_file_path.name,
        "sequence_metadata_file_name": sequence_metadata_file_path.name,
        "profiling_log_file_name": profiling_log_file_path.name,
        "sequence_count": sequence_count,
        "embedding_dimension": total_embedding_dimensions,
        "projected_dimensions_per_mask": projected_dimensions_per_mask,
        "mask_count": len(masks),
        "mask_length": len(masks[0]),
        "masks": [list(mask) for mask in masks],
        "fasta_type": "AA",
        "composition": composition,
        "projection": projection,
        "random_seed": random_seed,
        "mask_projection_seeds": [
            random_seed + index for index in range(len(masks))
        ],
        "chunk_size": chunk_size,
        "n_jobs": n_jobs,
        "dtype": "float32",
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "source_fasta_snapshot_relative_path": source_fasta_manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_fasta_snapshot_directory_name": source_fasta_manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_fasta_file_name": source_fasta_manifest.get("fasta_file_name"),
        "source_fasta_file_sha256": source_fasta_manifest.get("fasta_file_sha256"),
        "source_fasta_manifest_sha256": sha256_of_file(
            input_file_path=Path(source_fasta_snapshot_payload["manifest_file_path"])
        ),
        "source_metadata_snapshot_relative_path": source_fasta_manifest.get(
            "source_metadata_snapshot_relative_path"
        ),
        "source_xml_snapshot_relative_path": source_fasta_manifest.get(
            "source_xml_snapshot_relative_path"
        ),
        "search_query": source_fasta_manifest.get("search_query"),
        "translated_query": source_fasta_manifest.get("translated_query"),
        "sweep_package_root_directory": str(scripts_sweep_root_directory),
        "sweep_orthbase_source": str(
            scripts_sweep_root_directory / "sweep" / "orthbase.py"
        ),
        "embeddings_file_sha256": sha256_of_file(input_file_path=embeddings_file_path),
        "sequence_metadata_file_sha256": sha256_of_file(
            input_file_path=sequence_metadata_file_path
        ),
        "profiling_log_file_sha256": sha256_of_file(input_file_path=profiling_log_file_path),
        "profiling_row_count": total_profiling_rows,
        "profiling_summary": {
            "snapshot_total_elapsed_seconds": (
                snapshot_total_row.get("elapsed_seconds")
                if snapshot_total_row is not None
                else None
            ),
            "snapshot_total_cpu_seconds": (
                snapshot_total_row.get("cpu_seconds")
                if snapshot_total_row is not None
                else None
            ),
            "snapshot_total_rss_memory_delta_mb": (
                snapshot_total_row.get("rss_memory_delta_mb")
                if snapshot_total_row is not None
                else None
            ),
            "mask_work_unit_count": sum(
                1 for row in profiling_rows if row.get("work_unit_kind") == "mask_total"
            ),
            "chunk_work_unit_count": sum(
                1 for row in profiling_rows if row.get("work_unit_kind") == "mask_chunk"
            ),
        },
    }


def save_sweep_genes_snapshot(
    *,
    snapshot_root_directory: PathLike,
    source_fasta_snapshot_directory: PathLike,
    scripts_sweep_root_directory: PathLike,
    masks: Sequence[Sequence[int]] = DEFAULT_SWEEP_GENES_MASKS,
    projected_dimensions_per_mask: int = DEFAULT_PROJECTED_DIMENSIONS_PER_MASK,
    composition: str = "binary",
    projection: bool = True,
    random_seed: int = 42,
    chunk_size: int = 256,
    n_jobs: int = 1,
    update_latest_directory: bool = True,
) -> SweepGenesSnapshotResult:
    if not projection:
        raise ValueError(
            "SWeeP Genes snapshot creation requires projection=True to materialize the compact embedding matrix."
        )

    normalized_masks = _normalize_masks(masks)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_fasta_snapshot_directory = _as_path(source_fasta_snapshot_directory)
    resolved_scripts_sweep_root_directory = _as_path(scripts_sweep_root_directory)
    total_embedding_dimensions = len(normalized_masks) * projected_dimensions_per_mask
    embeddings_file_name = _build_embeddings_file_name(total_embedding_dimensions)

    source_fasta_snapshot_payload = load_fasta_snapshot_by_directory(
        snapshot_directory=resolved_source_fasta_snapshot_directory,
    )
    fasta_file_path = Path(source_fasta_snapshot_payload["fasta_file_path"])
    source_fasta_manifest = source_fasta_snapshot_payload.get("manifest")
    if not isinstance(source_fasta_manifest, Mapping):
        raise RuntimeError(
            "Resolved FASTA snapshot payload is missing a valid source manifest."
        )

    with fasta_file_path.open("r", encoding="utf-8") as fasta_file_handle:
        sequence_records = list(SeqIO.parse(fasta_file_handle, "fasta"))
    if not sequence_records:
        raise RuntimeError("The FASTA snapshot does not contain any sequences.")

    sweep = _load_vendored_sweep_package(
        scripts_sweep_root_directory=resolved_scripts_sweep_root_directory,
    )
    process = psutil.Process()

    snapshot_created_at_utc = _current_utc_timestamp()
    snapshot_started_at_utc = _current_utc_timestamp()
    snapshot_elapsed_started_at = time.perf_counter()
    snapshot_cpu_seconds_before = _get_process_cpu_seconds(process)
    snapshot_rss_memory_before_mb = _get_process_rss_memory_mb(process)
    search_query = str(
        source_fasta_manifest.get("search_query")
        or source_fasta_manifest.get("translated_query")
        or fasta_file_path.name
    )
    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=snapshot_created_at_utc,
        search_query=search_query,
    )
    immutable_snapshot_directory = (
        resolved_snapshot_root_directory / "snapshots" / snapshot_directory_name
    )
    immutable_snapshot_directory.mkdir(parents=True, exist_ok=False)
    immutable_snapshot_complete = False

    try:
        embeddings_file_path = immutable_snapshot_directory / embeddings_file_name
        sequence_metadata_file_path = (
            immutable_snapshot_directory / DEFAULT_SEQUENCE_METADATA_FILE_NAME
        )
        profiling_log_file_path = (
            immutable_snapshot_directory / DEFAULT_PROFILING_LOG_FILE_NAME
        )
        manifest_file_path = immutable_snapshot_directory / "manifest.json"
        profiling_rows: list[dict[str, Any]] = []

        _write_sequence_metadata_csv(
            sequence_records=sequence_records,
            output_file_path=sequence_metadata_file_path,
        )

        embeddings_memmap = np.lib.format.open_memmap(
            embeddings_file_path,
            mode="w+",
            dtype=np.float32,
            shape=(len(sequence_records), total_embedding_dimensions),
        )

        for mask_index, mask in enumerate(normalized_masks):
            mask_started_at_utc = _current_utc_timestamp()
            mask_elapsed_started_at = time.perf_counter()
            mask_cpu_seconds_before = _get_process_cpu_seconds(process)
            mask_rss_memory_before_mb = _get_process_rss_memory_mb(process)
            orth_mat = sweep.orthbase(
                sweep.calc_proj_mat_size(mask, "AA"),
                projected_dimensions_per_mask,
                seed=random_seed + mask_index,
            ).astype(np.float32, copy=False)
            column_start = mask_index * projected_dimensions_per_mask
            column_end = column_start + projected_dimensions_per_mask
            mask_label = "".join(str(mask_value) for mask_value in mask)

            def _record_mask_chunk(chunk_metrics: Mapping[str, Any]) -> None:
                profiling_rows.append(
                    _finalize_profiling_row(
                        row_payload={
                            "work_unit_kind": "mask_chunk",
                            "mask_index": mask_index,
                            "mask_label": mask_label,
                            "mask_seed": random_seed + mask_index,
                            "mask_values": " ".join(str(mask_value) for mask_value in mask),
                            "column_start_index": column_start,
                            "column_end_index_exclusive": column_end,
                            "embedding_dimension": projected_dimensions_per_mask,
                            **dict(chunk_metrics),
                        }
                    )
                )

            projected_embeddings = sweep.fas2sweep(
                sequence_records,
                orth_mat=orth_mat,
                mask=mask,
                composition=composition,
                chunk_size=chunk_size,
                projection=True,
                fasta_type="AA",
                n_jobs=n_jobs,
                dtype=np.float32,
                work_unit_callback=_record_mask_chunk,
            )
            embeddings_memmap[:, column_start:column_end] = projected_embeddings
            mask_completed_at_utc = _current_utc_timestamp()
            mask_cpu_seconds_after = _get_process_cpu_seconds(process)
            mask_rss_memory_after_mb = _get_process_rss_memory_mb(process)
            profiling_rows.append(
                _finalize_profiling_row(
                    row_payload={
                        "work_unit_kind": "mask_total",
                        "mask_index": mask_index,
                        "mask_label": mask_label,
                        "mask_seed": random_seed + mask_index,
                        "mask_values": " ".join(str(mask_value) for mask_value in mask),
                        "column_start_index": column_start,
                        "column_end_index_exclusive": column_end,
                        "embedding_dimension": projected_dimensions_per_mask,
                        "row_count": len(sequence_records),
                        "started_at_utc": mask_started_at_utc,
                        "completed_at_utc": mask_completed_at_utc,
                        "elapsed_seconds": float(
                            time.perf_counter() - mask_elapsed_started_at
                        ),
                        "cpu_seconds": float(
                            mask_cpu_seconds_after - mask_cpu_seconds_before
                        ),
                        "rss_memory_before_mb": mask_rss_memory_before_mb,
                        "rss_memory_after_mb": mask_rss_memory_after_mb,
                        "rss_memory_delta_mb": float(
                            mask_rss_memory_after_mb - mask_rss_memory_before_mb
                        ),
                    }
                )
            )

        embeddings_memmap.flush()
        del embeddings_memmap
        snapshot_completed_at_utc = _current_utc_timestamp()
        snapshot_cpu_seconds_after = _get_process_cpu_seconds(process)
        snapshot_rss_memory_after_mb = _get_process_rss_memory_mb(process)
        profiling_rows.append(
            _finalize_profiling_row(
                row_payload={
                    "work_unit_kind": "snapshot_total",
                    "mask_index": None,
                    "mask_label": None,
                    "mask_seed": None,
                    "mask_values": None,
                    "column_start_index": 0,
                    "column_end_index_exclusive": total_embedding_dimensions,
                    "embedding_dimension": total_embedding_dimensions,
                    "row_count": len(sequence_records),
                    "started_at_utc": snapshot_started_at_utc,
                    "completed_at_utc": snapshot_completed_at_utc,
                    "elapsed_seconds": float(
                        time.perf_counter() - snapshot_elapsed_started_at
                    ),
                    "cpu_seconds": float(
                        snapshot_cpu_seconds_after - snapshot_cpu_seconds_before
                    ),
                    "rss_memory_before_mb": snapshot_rss_memory_before_mb,
                    "rss_memory_after_mb": snapshot_rss_memory_after_mb,
                    "rss_memory_delta_mb": float(
                        snapshot_rss_memory_after_mb - snapshot_rss_memory_before_mb
                    ),
                }
            )
        )
        _write_profiling_log_csv(
            profiling_rows=profiling_rows,
            output_file_path=profiling_log_file_path,
        )

        immutable_snapshot_relative_path = str(
            Path("snapshots") / snapshot_directory_name
        )
        manifest_payload = _build_sweep_genes_snapshot_manifest(
            snapshot_created_at_utc=snapshot_created_at_utc,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
            embeddings_file_path=embeddings_file_path,
            sequence_metadata_file_path=sequence_metadata_file_path,
            profiling_log_file_path=profiling_log_file_path,
            sequence_count=len(sequence_records),
            masks=normalized_masks,
            projected_dimensions_per_mask=projected_dimensions_per_mask,
            composition=composition,
            projection=projection,
            random_seed=random_seed,
            chunk_size=chunk_size,
            n_jobs=n_jobs,
            scripts_sweep_root_directory=resolved_scripts_sweep_root_directory,
            source_fasta_snapshot_payload=source_fasta_snapshot_payload,
            profiling_rows=profiling_rows,
        )
        write_json_atomic(
            payload=manifest_payload,
            output_file_path=manifest_file_path,
        )
        immutable_snapshot_complete = True

        if update_latest_directory:
            _replace_latest_directory(
                latest_directory=resolved_snapshot_root_directory / "latest",
                files_to_copy=[
                    (embeddings_file_path, embeddings_file_path.name),
                    (sequence_metadata_file_path, sequence_metadata_file_path.name),
                    (profiling_log_file_path, profiling_log_file_path.name),
                    (manifest_file_path, manifest_file_path.name),
                ],
            )

    except Exception:
        if not immutable_snapshot_complete and immutable_snapshot_directory.exists():
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

    return SweepGenesSnapshotResult(
        snapshot_directory=immutable_snapshot_directory,
        manifest_file_path=manifest_file_path,
        embeddings_file_path=embeddings_file_path,
        sequence_metadata_file_path=sequence_metadata_file_path,
        profiling_log_file_path=profiling_log_file_path,
        source_fasta_snapshot_directory=resolved_source_fasta_snapshot_directory,
    )


def _validate_loaded_sweep_genes_snapshot_payload(
    *,
    snapshot_directory: PathLike,
    manifest_payload: Mapping[str, object],
    require_profiling_log: bool = False,
    require_artifact_hashes: bool = False,
) -> tuple[Path, Path, Optional[Path]]:
    resolved_snapshot_directory = _as_path(snapshot_directory)

    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != "sweep_genes_embedding_snapshot":
        raise RuntimeError(
            "Saved SWeeP Genes snapshot manifest artifact_type mismatch. "
            f"Expected 'sweep_genes_embedding_snapshot', got {artifact_type!r}."
        )

    embeddings_file_name = manifest_payload.get("embeddings_file_name")
    if not isinstance(embeddings_file_name, str) or not embeddings_file_name.strip():
        raise RuntimeError(
            "Saved SWeeP Genes snapshot manifest must define a non-empty embeddings_file_name."
        )

    sequence_metadata_file_name = manifest_payload.get("sequence_metadata_file_name")
    if (
        not isinstance(sequence_metadata_file_name, str)
        or not sequence_metadata_file_name.strip()
    ):
        raise RuntimeError(
            "Saved SWeeP Genes snapshot manifest must define a non-empty sequence_metadata_file_name."
        )

    resolved_embeddings_file_path = resolved_snapshot_directory / embeddings_file_name
    resolved_sequence_metadata_file_path = (
        resolved_snapshot_directory / sequence_metadata_file_name
    )
    profiling_log_file_name = manifest_payload.get("profiling_log_file_name")
    resolved_profiling_log_file_path: Optional[Path] = None
    if isinstance(profiling_log_file_name, str) and profiling_log_file_name.strip():
        resolved_profiling_log_file_path = (
            resolved_snapshot_directory / profiling_log_file_name
        )

    if not resolved_embeddings_file_path.exists():
        raise FileNotFoundError(
            f"Saved SWeeP Genes embeddings file not found: {resolved_embeddings_file_path}."
        )
    if not resolved_sequence_metadata_file_path.exists():
        raise FileNotFoundError(
            "Saved SWeeP Genes sequence metadata file not found: "
            f"{resolved_sequence_metadata_file_path}."
        )
    if require_profiling_log and resolved_profiling_log_file_path is None:
        raise RuntimeError(
            "Saved SWeeP Genes snapshot manifest must define a profiling_log_file_name."
        )
    if (
        resolved_profiling_log_file_path is not None
        and not resolved_profiling_log_file_path.exists()
    ):
        raise FileNotFoundError(
            "Saved SWeeP Genes profiling log file not found: "
            f"{resolved_profiling_log_file_path}."
        )

    expected_embeddings_file_sha256 = manifest_payload.get("embeddings_file_sha256")
    if require_artifact_hashes and expected_embeddings_file_sha256 is None:
        raise RuntimeError(
            "Saved SWeeP Genes snapshot manifest must define embeddings_file_sha256."
        )
    if expected_embeddings_file_sha256 is not None:
        actual_embeddings_file_sha256 = sha256_of_file(
            input_file_path=resolved_embeddings_file_path
        )
        if actual_embeddings_file_sha256 != expected_embeddings_file_sha256:
            raise RuntimeError(
                "Saved SWeeP Genes embeddings file hash mismatch. "
                f"Expected {expected_embeddings_file_sha256}, got {actual_embeddings_file_sha256}."
            )

    expected_sequence_metadata_file_sha256 = manifest_payload.get(
        "sequence_metadata_file_sha256"
    )
    if require_artifact_hashes and expected_sequence_metadata_file_sha256 is None:
        raise RuntimeError(
            "Saved SWeeP Genes snapshot manifest must define "
            "sequence_metadata_file_sha256."
        )
    if expected_sequence_metadata_file_sha256 is not None:
        actual_sequence_metadata_file_sha256 = sha256_of_file(
            input_file_path=resolved_sequence_metadata_file_path
        )
        if actual_sequence_metadata_file_sha256 != expected_sequence_metadata_file_sha256:
            raise RuntimeError(
                "Saved SWeeP Genes sequence metadata file hash mismatch. "
                "Expected "
                f"{expected_sequence_metadata_file_sha256}, got "
                f"{actual_sequence_metadata_file_sha256}."
            )

    expected_profiling_log_file_sha256 = manifest_payload.get(
        "profiling_log_file_sha256"
    )
    if (
        require_artifact_hashes
        and resolved_profiling_log_file_path is not None
        and expected_profiling_log_file_sha256 is None
    ):
        raise RuntimeError(
            "Saved SWeeP Genes snapshot manifest must define profiling_log_file_sha256."
        )
    if (
        expected_profiling_log_file_sha256 is not None
        and resolved_profiling_log_file_path is not None
    ):
        actual_profiling_log_file_sha256 = sha256_of_file(
            input_file_path=resolved_profiling_log_file_path
        )
        if actual_profiling_log_file_sha256 != expected_profiling_log_file_sha256:
            raise RuntimeError(
                "Saved SWeeP Genes profiling log file hash mismatch. "
                f"Expected {expected_profiling_log_file_sha256}, got "
                f"{actual_profiling_log_file_sha256}."
            )

    return (
        resolved_embeddings_file_path,
        resolved_sequence_metadata_file_path,
        resolved_profiling_log_file_path,
    )


def load_sweep_genes_snapshot_manifest_by_directory(
    *,
    snapshot_directory: PathLike,
    require_profiling_log: bool = False,
    require_artifact_hashes: bool = False,
) -> dict[str, object]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / "manifest.json"
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved SWeeP Genes snapshot manifest must deserialize into a JSON object."
        )

    embeddings_file_path, sequence_metadata_file_path, profiling_log_file_path = (
        _validate_loaded_sweep_genes_snapshot_payload(
            snapshot_directory=resolved_snapshot_directory,
            manifest_payload=manifest_payload,
            require_profiling_log=require_profiling_log,
            require_artifact_hashes=require_artifact_hashes,
        )
    )

    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "embeddings_file_path": embeddings_file_path,
        "sequence_metadata_file_path": sequence_metadata_file_path,
        "profiling_log_file_path": profiling_log_file_path,
    }


def load_sweep_genes_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> dict[str, object]:
    manifest_only_payload = load_sweep_genes_snapshot_manifest_by_directory(
        snapshot_directory=snapshot_directory,
    )
    resolved_snapshot_directory = Path(manifest_only_payload["snapshot_directory"])
    manifest_file_path = Path(manifest_only_payload["manifest_file_path"])
    manifest_payload = manifest_only_payload["manifest"]
    embeddings_file_path = Path(manifest_only_payload["embeddings_file_path"])
    sequence_metadata_file_path = Path(
        manifest_only_payload["sequence_metadata_file_path"]
    )
    profiling_log_payload = manifest_only_payload["profiling_log_file_path"]
    profiling_log_file_path = (
        Path(profiling_log_payload)
        if profiling_log_payload is not None
        else None
    )

    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": manifest_payload,
        "embeddings_file_path": embeddings_file_path,
        "sequence_metadata_file_path": sequence_metadata_file_path,
        "profiling_log_file_path": profiling_log_file_path,
        "embeddings": np.load(embeddings_file_path, mmap_mode="r"),
        "sequence_metadata": pd.read_csv(sequence_metadata_file_path),
        "profiling_log": (
            pd.read_csv(profiling_log_file_path)
            if profiling_log_file_path is not None
            else None
        ),
    }


def load_latest_sweep_genes_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return load_sweep_genes_snapshot_by_directory(
        snapshot_directory=resolved_snapshot_root_directory / "latest",
    )


def latest_sweep_genes_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
) -> bool:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"
    manifest_file_path = latest_directory / "manifest.json"
    if not latest_directory.exists() or not manifest_file_path.exists():
        return False

    try:
        load_sweep_genes_snapshot_manifest_by_directory(
            snapshot_directory=latest_directory,
            require_profiling_log=True,
            require_artifact_hashes=True,
        )
    except (FileNotFoundError, RuntimeError, OSError, ValueError):
        return False

    return True


def _normalize_snapshot_path_identity(path_value: object) -> Optional[str]:
    if not isinstance(path_value, (str, Path)) or not str(path_value).strip():
        return None
    return str(Path(path_value).expanduser().resolve(strict=False))


def _sweep_genes_snapshot_manifest_matches_request(
    *,
    manifest_payload: Mapping[str, object],
    source_fasta_snapshot_payload: Mapping[str, object],
    scripts_sweep_root_directory: PathLike,
    masks: Sequence[Sequence[int]],
    projected_dimensions_per_mask: int,
    composition: str,
    projection: bool,
    random_seed: int,
    chunk_size: int,
    n_jobs: int,
) -> bool:
    source_fasta_manifest = source_fasta_snapshot_payload.get("manifest")
    if not isinstance(source_fasta_manifest, Mapping):
        return False

    source_fasta_manifest_file_path = source_fasta_snapshot_payload.get(
        "manifest_file_path"
    )
    if not isinstance(source_fasta_manifest_file_path, Path):
        source_fasta_manifest_file_path = Path(str(source_fasta_manifest_file_path))

    expected_masks = _normalize_masks(masks)
    expected_source_values = {
        "source_fasta_snapshot_relative_path": source_fasta_manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_fasta_snapshot_directory_name": source_fasta_manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_fasta_file_name": source_fasta_manifest.get("fasta_file_name"),
        "source_fasta_file_sha256": source_fasta_manifest.get("fasta_file_sha256"),
        "source_fasta_manifest_sha256": sha256_of_file(
            input_file_path=source_fasta_manifest_file_path
        ),
    }
    for key, expected_value in expected_source_values.items():
        if manifest_payload.get(key) != expected_value:
            return False

    expected_sweep_root_identity = _normalize_snapshot_path_identity(
        scripts_sweep_root_directory
    )
    actual_sweep_root_identity = _normalize_snapshot_path_identity(
        manifest_payload.get("sweep_package_root_directory")
    )
    if actual_sweep_root_identity != expected_sweep_root_identity:
        return False

    expected_total_embedding_dimensions = (
        len(expected_masks) * projected_dimensions_per_mask
    )
    return (
        manifest_payload.get("artifact_type") == "sweep_genes_embedding_snapshot"
        and manifest_payload.get("masks") == expected_masks
        and manifest_payload.get("projected_dimensions_per_mask")
        == projected_dimensions_per_mask
        and manifest_payload.get("embedding_dimension")
        == expected_total_embedding_dimensions
        and manifest_payload.get("composition") == composition
        and manifest_payload.get("projection") == projection
        and manifest_payload.get("random_seed") == random_seed
        and manifest_payload.get("chunk_size") == chunk_size
        and manifest_payload.get("n_jobs") == n_jobs
    )


def _find_matching_immutable_sweep_genes_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
    source_fasta_snapshot_payload: Mapping[str, object],
    scripts_sweep_root_directory: PathLike,
    masks: Sequence[Sequence[int]],
    projected_dimensions_per_mask: int,
    composition: str,
    projection: bool,
    random_seed: int,
    chunk_size: int,
    n_jobs: int,
) -> Optional[Path]:
    for snapshot_directory in reversed(
        list_saved_sweep_genes_snapshot_directories(
            snapshot_root_directory=snapshot_root_directory
        )
    ):
        try:
            manifest_only_payload = load_sweep_genes_snapshot_manifest_by_directory(
                snapshot_directory=snapshot_directory,
                require_profiling_log=True,
                require_artifact_hashes=True,
            )
        except (FileNotFoundError, RuntimeError, OSError, ValueError):
            continue

        manifest_payload = manifest_only_payload.get("manifest")
        if not isinstance(manifest_payload, Mapping):
            continue
        if _sweep_genes_snapshot_manifest_matches_request(
            manifest_payload=manifest_payload,
            source_fasta_snapshot_payload=source_fasta_snapshot_payload,
            scripts_sweep_root_directory=scripts_sweep_root_directory,
            masks=masks,
            projected_dimensions_per_mask=projected_dimensions_per_mask,
            composition=composition,
            projection=projection,
            random_seed=random_seed,
            chunk_size=chunk_size,
            n_jobs=n_jobs,
        ):
            return snapshot_directory

    return None


def _publish_sweep_genes_snapshot_as_latest(
    *,
    snapshot_root_directory: PathLike,
    snapshot_directory: PathLike,
) -> None:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    manifest_only_payload = load_sweep_genes_snapshot_manifest_by_directory(
        snapshot_directory=snapshot_directory,
        require_profiling_log=True,
        require_artifact_hashes=True,
    )
    _replace_latest_directory(
        latest_directory=resolved_snapshot_root_directory / "latest",
        files_to_copy=[
            (
                Path(manifest_only_payload["embeddings_file_path"]),
                Path(manifest_only_payload["embeddings_file_path"]).name,
            ),
            (
                Path(manifest_only_payload["sequence_metadata_file_path"]),
                Path(manifest_only_payload["sequence_metadata_file_path"]).name,
            ),
            (
                Path(manifest_only_payload["profiling_log_file_path"]),
                Path(manifest_only_payload["profiling_log_file_path"]).name,
            ),
            (
                Path(manifest_only_payload["manifest_file_path"]),
                Path(manifest_only_payload["manifest_file_path"]).name,
            ),
        ],
    )


def resolve_sweep_genes_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    source_fasta_snapshot_root_directory: PathLike,
    scripts_sweep_root_directory: PathLike,
    masks: Sequence[Sequence[int]] = DEFAULT_SWEEP_GENES_MASKS,
    projected_dimensions_per_mask: int = DEFAULT_PROJECTED_DIMENSIONS_PER_MASK,
    composition: str = "binary",
    projection: bool = True,
    random_seed: int = 42,
    chunk_size: int = 256,
    n_jobs: int = 1,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_fasta_snapshot_root_directory = _as_path(
        source_fasta_snapshot_root_directory
    )

    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        latest_is_available = latest_sweep_genes_snapshot_is_available(
            snapshot_root_directory=resolved_snapshot_root_directory
        )
        if not latest_is_available:
            raise FileNotFoundError(
                "No latest SWeeP Genes snapshot directory was found. Run the workflow once "
                "with snapshot_mode='create_new' before using 'reuse_latest'."
            )

        latest_snapshot_payload = load_latest_sweep_genes_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory
        )
        return latest_snapshot_payload

    if resolved_snapshot_mode not in {
        SnapshotMode.create_new,
        SnapshotMode.reuse_latest_or_create,
    }:
        raise ValueError(
            "Invalid snapshot_mode. Expected one of: "
            "'create_new', 'reuse_latest', 'reuse_latest_or_create'."
        )

    try:
        source_fasta_snapshot_payload = load_latest_fasta_snapshot(
            snapshot_root_directory=resolved_source_fasta_snapshot_root_directory
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "No reusable source FASTA snapshot was found for the SWeeP Genes workflow. "
            "Create the upstream FASTA snapshot before running the SWeeP Genes step."
        ) from exc

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_sweep_genes_snapshot_is_available(
            snapshot_root_directory=resolved_snapshot_root_directory
        )
    ):
        latest_manifest_only_payload = load_sweep_genes_snapshot_manifest_by_directory(
            snapshot_directory=resolved_snapshot_root_directory / "latest",
            require_profiling_log=True,
            require_artifact_hashes=True,
        )
        latest_manifest_payload = latest_manifest_only_payload.get("manifest")
        if isinstance(
            latest_manifest_payload, Mapping
        ) and _sweep_genes_snapshot_manifest_matches_request(
            manifest_payload=latest_manifest_payload,
            source_fasta_snapshot_payload=source_fasta_snapshot_payload,
            scripts_sweep_root_directory=scripts_sweep_root_directory,
            masks=masks,
            projected_dimensions_per_mask=projected_dimensions_per_mask,
            composition=composition,
            projection=projection,
            random_seed=random_seed,
            chunk_size=chunk_size,
            n_jobs=n_jobs,
        ):
            print("Latest SWeeP Genes snapshot is compatible. Reusing frozen snapshot.")
            latest_snapshot_payload = load_latest_sweep_genes_snapshot(
                snapshot_root_directory=resolved_snapshot_root_directory
            )
            return latest_snapshot_payload

    if resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create:
        matching_snapshot_directory = (
            _find_matching_immutable_sweep_genes_snapshot_directory(
                snapshot_root_directory=resolved_snapshot_root_directory,
                source_fasta_snapshot_payload=source_fasta_snapshot_payload,
                scripts_sweep_root_directory=scripts_sweep_root_directory,
                masks=masks,
                projected_dimensions_per_mask=projected_dimensions_per_mask,
                composition=composition,
                projection=projection,
                random_seed=random_seed,
                chunk_size=chunk_size,
                n_jobs=n_jobs,
            )
        )
        if matching_snapshot_directory is not None:
            print(
                "Found a compatible immutable SWeeP Genes snapshot. Reusing frozen snapshot."
            )
            if update_latest_directory:
                _publish_sweep_genes_snapshot_as_latest(
                    snapshot_root_directory=resolved_snapshot_root_directory,
                    snapshot_directory=matching_snapshot_directory,
                )
            return load_sweep_genes_snapshot_by_directory(
                snapshot_directory=matching_snapshot_directory
            )

    saved_snapshot = save_sweep_genes_snapshot(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_fasta_snapshot_directory=Path(
            source_fasta_snapshot_payload["snapshot_directory"]
        ),
        scripts_sweep_root_directory=scripts_sweep_root_directory,
        masks=masks,
        projected_dimensions_per_mask=projected_dimensions_per_mask,
        composition=composition,
        projection=projection,
        random_seed=random_seed,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        update_latest_directory=update_latest_directory,
    )

    return load_sweep_genes_snapshot_by_directory(
        snapshot_directory=saved_snapshot.snapshot_directory
    )


def list_saved_sweep_genes_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    return list_saved_snapshot_directories(snapshot_root_directory=snapshot_root_directory)


def get_most_recent_sweep_genes_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory
    )
