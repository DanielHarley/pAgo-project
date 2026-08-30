from __future__ import annotations

import shutil
import tempfile
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, TypeAlias

from src.pago_pipeline.ncbi_esearch_preflight import (
    DEFAULT_MAX_UID_COUNT,
    DEFAULT_SAMPLE_SIZE,
    EsearchPreflightResult,
    run_ncbi_esearch_preflight,
)
from src.pago_pipeline.ncbi_snapshot import (
    SnapshotMode,
    _coerce_snapshot_mode,
    _replace_latest_directory,
    build_snapshot_directory_name,
    get_most_recent_snapshot_directory,
    list_saved_snapshot_directories,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path

ARTIFACT_TYPE = "ncbi_esearch_preflight"
SNAPSHOT_FORMAT_VERSION = "1.0"
DEFAULT_MANIFEST_FILE_NAME = "manifest.json"
DEFAULT_PREFLIGHT_REPORT_FILE_NAME = "preflight_report.json"
DEFAULT_SAMPLE_UID_FILE_NAME = "sample_protein_uids.txt"


@dataclass(frozen=True)
class NCBIEsearchPreflightSnapshotResult:
    snapshot_directory: Path
    snapshot_root_directory: Path
    manifest_file_path: Path
    preflight_report_file_path: Path
    sample_uid_file_path: Path
    result_count: int
    exceeds_max_uid_count: bool


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _current_utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _write_text_lines_atomic(*, text_lines: list[str], output_file_path: Path) -> None:
    output_file_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        dir=output_file_path.parent,
        encoding="utf-8",
        newline="\n",
        suffix=".txt",
    ) as temporary_file:
        resolved_temporary_file_path = Path(temporary_file.name)
        for text_line in text_lines:
            temporary_file.write(f"{text_line}\n")
    resolved_temporary_file_path.replace(output_file_path)


def _build_manifest_payload(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    preflight_result: EsearchPreflightResult,
    output_file_path_by_key: Mapping[str, Path],
) -> dict[str, object]:
    output_files = {
        key: {
            "file_name": file_path.name,
            "path": str(file_path),
            "sha256": sha256_of_file(input_file_path=file_path),
        }
        for key, file_path in output_file_path_by_key.items()
    }
    return {
        "artifact_type": ARTIFACT_TYPE,
        "snapshot_format_version": SNAPSHOT_FORMAT_VERSION,
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "manifest_file_name": DEFAULT_MANIFEST_FILE_NAME,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "search_query": preflight_result.search_query,
        "translated_query": preflight_result.translated_query,
        "result_count": int(preflight_result.result_count),
        "max_uid_count": int(preflight_result.max_uid_count),
        "exceeds_max_uid_count": bool(preflight_result.exceeds_max_uid_count),
        "history_web_env": preflight_result.history_web_env,
        "history_query_key": preflight_result.history_query_key,
        "retrieved_at_utc": preflight_result.retrieved_at_utc,
        "sample_requested_count": int(preflight_result.sample_requested_count),
        "sample_record_count": int(preflight_result.sample_record_count),
        "sample_records_with_sequence": int(
            preflight_result.sample_records_with_sequence
        ),
        "sample_records_missing_sequence": int(
            preflight_result.sample_records_missing_sequence
        ),
        "sample_records_with_extractable_uid": int(
            preflight_result.sample_records_with_extractable_uid
        ),
        "sample_fetch_error": preflight_result.sample_fetch_error,
        "python_version": preflight_result.python_version,
        "biopython_version": preflight_result.biopython_version,
        "output_files": output_files,
    }


def save_ncbi_esearch_preflight_snapshot(
    *,
    snapshot_root_directory: PathLike,
    search_query: str,
    ncbi_email: str,
    ncbi_api_key: Optional[str] = None,
    max_uid_count: int = DEFAULT_MAX_UID_COUNT,
    sample_size: int = DEFAULT_SAMPLE_SIZE,
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
    update_latest_directory: bool = True,
) -> NCBIEsearchPreflightSnapshotResult:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    preflight_result = run_ncbi_esearch_preflight(
        search_query=search_query,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        max_uid_count=max_uid_count,
        sample_size=sample_size,
        ssl_ca_file=ssl_ca_file,
        ssl_ca_directory=ssl_ca_directory,
    )

    snapshot_created_at_utc = _current_utc_timestamp()
    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=snapshot_created_at_utc,
        search_query=search_query,
    )
    immutable_snapshot_directory = (
        resolved_snapshot_root_directory / "snapshots" / snapshot_directory_name
    )
    immutable_snapshot_directory.mkdir(parents=True, exist_ok=False)
    immutable_snapshot_relative_path = str(
        Path("snapshots") / snapshot_directory_name
    )
    immutable_snapshot_complete = False

    try:
        output_file_path_by_key = {
            "preflight_report_file": (
                immutable_snapshot_directory / DEFAULT_PREFLIGHT_REPORT_FILE_NAME
            ),
            "sample_uid_file": (
                immutable_snapshot_directory / DEFAULT_SAMPLE_UID_FILE_NAME
            ),
        }
        write_json_atomic(
            payload=preflight_result.as_report_payload(),
            output_file_path=output_file_path_by_key["preflight_report_file"],
        )
        _write_text_lines_atomic(
            text_lines=list(preflight_result.sample_uid_list),
            output_file_path=output_file_path_by_key["sample_uid_file"],
        )

        manifest_file_path = (
            immutable_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
        )
        manifest_payload = _build_manifest_payload(
            snapshot_created_at_utc=snapshot_created_at_utc,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
            preflight_result=preflight_result,
            output_file_path_by_key=output_file_path_by_key,
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
                    (
                        output_file_path_by_key["preflight_report_file"],
                        DEFAULT_PREFLIGHT_REPORT_FILE_NAME,
                    ),
                    (
                        output_file_path_by_key["sample_uid_file"],
                        DEFAULT_SAMPLE_UID_FILE_NAME,
                    ),
                    (manifest_file_path, DEFAULT_MANIFEST_FILE_NAME),
                ],
            )
    except Exception:
        if (
            not immutable_snapshot_complete
            and immutable_snapshot_directory.exists()
        ):
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

    return NCBIEsearchPreflightSnapshotResult(
        snapshot_directory=immutable_snapshot_directory,
        snapshot_root_directory=resolved_snapshot_root_directory,
        manifest_file_path=manifest_file_path,
        preflight_report_file_path=output_file_path_by_key["preflight_report_file"],
        sample_uid_file_path=output_file_path_by_key["sample_uid_file"],
        result_count=int(preflight_result.result_count),
        exceeds_max_uid_count=bool(preflight_result.exceeds_max_uid_count),
    )


def _validate_loaded_ncbi_esearch_preflight_payload(
    *,
    snapshot_directory: Path,
    manifest_payload: Mapping[str, object],
) -> dict[str, Path]:
    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != ARTIFACT_TYPE:
        raise RuntimeError(
            "Saved ESearch preflight snapshot artifact_type mismatch. "
            f"Expected {ARTIFACT_TYPE!r}, got {artifact_type!r}."
        )
    snapshot_format_version = manifest_payload.get("snapshot_format_version")
    if snapshot_format_version != SNAPSHOT_FORMAT_VERSION:
        raise RuntimeError(
            "Saved ESearch preflight snapshot snapshot_format_version mismatch. "
            f"Expected {SNAPSHOT_FORMAT_VERSION!r}, got {snapshot_format_version!r}."
        )
    output_files = manifest_payload.get("output_files")
    if not isinstance(output_files, Mapping):
        raise RuntimeError(
            "Saved ESearch preflight snapshot manifest must define output_files."
        )
    resolved: dict[str, Path] = {}
    for key, entry in output_files.items():
        if not isinstance(entry, Mapping):
            raise RuntimeError(f"output_files entry {key!r} must be an object.")
        file_name = entry.get("file_name")
        expected_sha256 = entry.get("sha256")
        if not isinstance(file_name, str) or not file_name.strip():
            raise RuntimeError(f"output_files entry {key!r} must define file_name.")
        if not isinstance(expected_sha256, str) or not expected_sha256.strip():
            raise RuntimeError(f"output_files entry {key!r} must define sha256.")
        file_path = snapshot_directory / file_name
        if not file_path.exists():
            raise FileNotFoundError(
                f"Saved ESearch preflight snapshot file not found: {file_path}."
            )
        if sha256_of_file(input_file_path=file_path) != expected_sha256:
            raise RuntimeError(
                f"Saved ESearch preflight snapshot hash mismatch for {key!r}."
            )
        resolved[key] = file_path
    return resolved


def load_ncbi_esearch_preflight_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved ESearch preflight snapshot manifest must deserialize into a "
            "JSON object."
        )
    resolved = _validate_loaded_ncbi_esearch_preflight_payload(
        snapshot_directory=resolved_snapshot_directory,
        manifest_payload=manifest_payload,
    )
    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "preflight_report_file_path": resolved["preflight_report_file"],
        "sample_uid_file_path": resolved["sample_uid_file"],
        "preflight_report": read_json_file(
            input_file_path=resolved["preflight_report_file"]
        ),
    }


def load_latest_ncbi_esearch_preflight_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return load_ncbi_esearch_preflight_snapshot_by_directory(
        snapshot_directory=resolved_snapshot_root_directory / "latest"
    )


def latest_ncbi_esearch_preflight_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
    search_query: Optional[str] = None,
) -> bool:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"
    latest_manifest_file_path = latest_directory / DEFAULT_MANIFEST_FILE_NAME
    if not latest_directory.exists() or not latest_manifest_file_path.exists():
        return False
    try:
        manifest_payload = read_json_file(input_file_path=latest_manifest_file_path)
        if not isinstance(manifest_payload, Mapping):
            return False
        _validate_loaded_ncbi_esearch_preflight_payload(
            snapshot_directory=latest_directory,
            manifest_payload=manifest_payload,
        )
        if search_query is not None and manifest_payload.get("search_query") != (
            search_query
        ):
            return False
    except (FileNotFoundError, RuntimeError, OSError, ValueError):
        return False
    return True


def _raise_if_exceeds_max_uid_count(
    *,
    manifest_payload: Mapping[str, object],
    allow_exceeds_max_uid_count: bool,
) -> None:
    if allow_exceeds_max_uid_count:
        return
    if bool(manifest_payload.get("exceeds_max_uid_count")):
        result_count = manifest_payload.get("result_count")
        max_uid_count = manifest_payload.get("max_uid_count")
        raise RuntimeError(
            "NCBI ESearch preflight reported "
            f"{result_count} protein UIDs, above the configured max_uid_count "
            f"of {max_uid_count}. The preflight report is materialized for "
            "audit, but the full retrieval is blocked. Raise max_uid_count "
            "explicitly (allow_exceeds_max_uid_count=True) only after "
            "reviewing translated_query and result_count."
        )


def resolve_ncbi_esearch_preflight_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    search_query: Optional[str] = None,
    ncbi_email: Optional[str] = None,
    ncbi_api_key: Optional[str] = None,
    max_uid_count: int = DEFAULT_MAX_UID_COUNT,
    sample_size: int = DEFAULT_SAMPLE_SIZE,
    allow_exceeds_max_uid_count: bool = False,
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    latest_is_available = latest_ncbi_esearch_preflight_snapshot_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory,
        search_query=search_query,
    )

    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            raise FileNotFoundError(
                "No latest ESearch preflight snapshot was found for this query. "
                "Run the workflow once with snapshot_mode='create_new' before "
                "using 'reuse_latest'."
            )
        payload = load_latest_ncbi_esearch_preflight_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory
        )
        _raise_if_exceeds_max_uid_count(
            manifest_payload=payload["manifest"],
            allow_exceeds_max_uid_count=allow_exceeds_max_uid_count,
        )
        return payload

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        print(
            "Latest ESearch preflight snapshot is available for this query. "
            "Reusing frozen snapshot."
        )
        payload = load_latest_ncbi_esearch_preflight_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory
        )
        _raise_if_exceeds_max_uid_count(
            manifest_payload=payload["manifest"],
            allow_exceeds_max_uid_count=allow_exceeds_max_uid_count,
        )
        return payload

    if resolved_snapshot_mode not in {
        SnapshotMode.create_new,
        SnapshotMode.reuse_latest_or_create,
    }:
        raise ValueError(
            "Invalid snapshot_mode. Expected one of: "
            "'create_new', 'reuse_latest', 'reuse_latest_or_create'."
        )

    if not ncbi_email:
        raise ValueError(
            "ncbi_email is required when an ESearch preflight must be run."
        )
    if not search_query or not search_query.strip():
        raise ValueError(
            "search_query is required when an ESearch preflight must be run."
        )

    saved_result = save_ncbi_esearch_preflight_snapshot(
        snapshot_root_directory=resolved_snapshot_root_directory,
        search_query=search_query,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        max_uid_count=max_uid_count,
        sample_size=sample_size,
        ssl_ca_file=ssl_ca_file,
        ssl_ca_directory=ssl_ca_directory,
        update_latest_directory=update_latest_directory,
    )
    payload = load_ncbi_esearch_preflight_snapshot_by_directory(
        snapshot_directory=saved_result.snapshot_directory
    )
    _raise_if_exceeds_max_uid_count(
        manifest_payload=payload["manifest"],
        allow_exceeds_max_uid_count=allow_exceeds_max_uid_count,
    )
    return payload


def list_saved_ncbi_esearch_preflight_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    return list_saved_snapshot_directories(
        snapshot_root_directory=snapshot_root_directory
    )


def get_most_recent_ncbi_esearch_preflight_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory
    )
