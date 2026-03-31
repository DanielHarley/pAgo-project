from __future__ import annotations

import hashlib
import shutil
from dataclasses import asdict
from enum import StrEnum
from pathlib import Path
from typing import Any, Dict, Optional, Union

from src.pago_pipeline.ncbi_api import (
    NCBIProteinUidFetchResult,
    fetch_ncbi_protein_uid_snapshot,
)
from src.pago_pipeline.storage import (
    read_json_file,
    read_text_lines_from_file,
    write_json_atomic,
)

PathLike = Union[str, Path]


class SnapshotMode(StrEnum):
    create_new = "create_new"
    reuse_latest = "reuse_latest"
    reuse_latest_or_create = "reuse_latest_or_create"


def _as_path(path_like: PathLike) -> Path:
    """
    Convert PathLike input into pathlib.Path explicitly.

    Keeping this conversion in one helper improves readability and helps
    static type checkers understand that downstream variables are true Path
    objects rather than Union[str, Path].
    """
    return Path(path_like)


def _coerce_snapshot_mode(
    snapshot_mode: SnapshotMode | str,
) -> SnapshotMode:
    """
    Normalize a runtime snapshot mode value into SnapshotMode.
    """
    if isinstance(snapshot_mode, SnapshotMode):
        return snapshot_mode

    try:
        return SnapshotMode(snapshot_mode)
    except ValueError as error:
        raise ValueError(
            "Invalid snapshot_mode. Expected one of: "
            "'create_new', 'reuse_latest', 'reuse_latest_or_create'."
        ) from error


def _sanitize_utc_timestamp_for_path(utc_timestamp: str) -> str:
    """
    Convert an ISO-like UTC timestamp into a filesystem-friendly string.

    Example:
        2026-03-21T18:42:03Z -> 2026-03-21T18-42-03Z
    """
    return utc_timestamp.replace(":", "-")


def _build_query_hash(search_query: str, hash_length: int = 12) -> str:
    """
    Build a short deterministic hash from the search query to reduce
    snapshot-directory name collisions and preserve provenance.
    """
    if hash_length <= 0:
        raise ValueError("hash_length must be a positive integer.")

    full_hash = hashlib.sha256(search_query.encode("utf-8")).hexdigest()
    return full_hash[:hash_length]


def build_snapshot_directory_name(
    *,
    retrieved_at_utc: str,
    search_query: str,
) -> str:
    """
    Build a deterministic, filesystem-friendly snapshot directory name.
    """
    safe_timestamp = _sanitize_utc_timestamp_for_path(retrieved_at_utc)
    query_hash = _build_query_hash(search_query=search_query)

    return f"{safe_timestamp}__q_{query_hash}"


def _write_text_lines(
    *,
    text_lines: list[str],
    output_file_path: PathLike,
) -> None:
    """
    Write one text line per entry with a trailing newline.
    """
    resolved_output_file_path = _as_path(output_file_path)
    resolved_output_file_path.parent.mkdir(parents=True, exist_ok=True)

    with resolved_output_file_path.open(
        "w",
        encoding="utf-8",
        newline="\n",
    ) as file_handle:
        for text_line in text_lines:
            file_handle.write(f"{text_line}\n")


def _build_snapshot_manifest(
    *,
    fetch_result: NCBIProteinUidFetchResult,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
) -> Dict[str, Any]:
    """
    Build the manifest payload persisted alongside protein_uids.txt.
    """
    manifest_payload = asdict(fetch_result)
    manifest_payload.pop("protein_uids", None)

    manifest_payload["snapshot_format_version"] = "1.0"
    manifest_payload["snapshot_file_name"] = "protein_uids.txt"
    manifest_payload["manifest_file_name"] = "manifest.json"
    manifest_payload["immutable_snapshot_directory_name"] = (
        immutable_snapshot_directory_name
    )
    manifest_payload["immutable_snapshot_relative_path"] = (
        immutable_snapshot_relative_path
    )

    return manifest_payload


def save_ncbi_protein_uid_snapshot(
    *,
    fetch_result: NCBIProteinUidFetchResult,
    snapshot_root_directory: PathLike,
    update_latest_directory: bool = True,
) -> Path:
    """
    Persist one immutable local snapshot of an NCBI protein-UID retrieval event.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=fetch_result.retrieved_at_utc,
        search_query=fetch_result.search_query,
    )

    immutable_snapshot_directory = (
        resolved_snapshot_root_directory / "snapshots" / snapshot_directory_name
    )
    immutable_snapshot_directory.mkdir(parents=True, exist_ok=False)

    immutable_snapshot_relative_path = str(
        Path("snapshots") / snapshot_directory_name
    )

    protein_uids_file_path = immutable_snapshot_directory / "protein_uids.txt"
    manifest_file_path = immutable_snapshot_directory / "manifest.json"

    _write_text_lines(
        text_lines=fetch_result.protein_uids,
        output_file_path=protein_uids_file_path,
    )

    manifest_payload = _build_snapshot_manifest(
        fetch_result=fetch_result,
        immutable_snapshot_directory_name=snapshot_directory_name,
        immutable_snapshot_relative_path=immutable_snapshot_relative_path,
    )
    write_json_atomic(
        payload=manifest_payload,
        output_file_path=manifest_file_path,
    )

    if update_latest_directory:
        latest_directory = resolved_snapshot_root_directory / "latest"
        latest_directory.mkdir(parents=True, exist_ok=True)

        shutil.copy2(protein_uids_file_path, latest_directory / "protein_uids.txt")
        shutil.copy2(manifest_file_path, latest_directory / "manifest.json")

    return immutable_snapshot_directory


def load_snapshot_manifest(
    *,
    manifest_file_path: PathLike,
) -> Dict[str, Any]:
    """
    Load a snapshot manifest from an explicit manifest path.
    """
    return read_json_file(input_file_path=manifest_file_path)


def load_snapshot_protein_uids(
    *,
    protein_uids_file_path: PathLike,
) -> list[str]:
    """
    Load protein UIDs from an explicit snapshot protein_uids.txt path.
    """
    return read_text_lines_from_file(input_file_path=protein_uids_file_path)


def load_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> Dict[str, Any]:
    """
    Load both manifest metadata and protein UIDs from a snapshot directory.
    """
    resolved_snapshot_directory = _as_path(snapshot_directory)

    manifest_file_path = resolved_snapshot_directory / "manifest.json"
    protein_uids_file_path = resolved_snapshot_directory / "protein_uids.txt"

    manifest_payload = load_snapshot_manifest(
        manifest_file_path=manifest_file_path,
    )
    protein_uids = load_snapshot_protein_uids(
        protein_uids_file_path=protein_uids_file_path,
    )

    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "protein_uids_file_path": protein_uids_file_path,
        "manifest": manifest_payload,
        "protein_uids": protein_uids,
    }


def load_latest_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> Dict[str, Any]:
    """
    Load the convenience latest snapshot copy.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"

    return load_snapshot_by_directory(snapshot_directory=latest_directory)


def get_latest_snapshot_manifest_path(
    *,
    snapshot_root_directory: PathLike,
) -> Path:
    """
    Return the manifest path for the convenience latest snapshot copy.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return resolved_snapshot_root_directory / "latest" / "manifest.json"


def get_latest_snapshot_protein_uids_path(
    *,
    snapshot_root_directory: PathLike,
) -> Path:
    """
    Return the protein_uids.txt path for the convenience latest snapshot copy.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return resolved_snapshot_root_directory / "latest" / "protein_uids.txt"


def get_snapshot_manifest_path(
    *,
    snapshot_directory: PathLike,
) -> Path:
    """
    Return the manifest path inside a specific immutable snapshot directory.
    """
    resolved_snapshot_directory = _as_path(snapshot_directory)
    return resolved_snapshot_directory / "manifest.json"


def get_snapshot_protein_uids_path(
    *,
    snapshot_directory: PathLike,
) -> Path:
    """
    Return the protein_uids.txt path inside a specific immutable snapshot directory.
    """
    resolved_snapshot_directory = _as_path(snapshot_directory)
    return resolved_snapshot_directory / "protein_uids.txt"


def list_saved_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    """
    List saved immutable snapshot directories in lexical order.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    snapshots_directory = resolved_snapshot_root_directory / "snapshots"

    if not snapshots_directory.exists():
        return []

    return sorted(
        path
        for path in snapshots_directory.iterdir()
        if path.is_dir()
    )


def get_most_recent_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    """
    Return the most recent immutable snapshot directory, or None if none exist.
    """
    snapshot_directories = list_saved_snapshot_directories(
        snapshot_root_directory=snapshot_root_directory,
    )

    if not snapshot_directories:
        return None

    return snapshot_directories[-1]


def latest_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
) -> bool:
    """
    Return True only if the convenience latest snapshot copy is complete.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    latest_directory = resolved_snapshot_root_directory / "latest"
    latest_manifest_file_path = latest_directory / "manifest.json"
    latest_protein_uids_file_path = latest_directory / "protein_uids.txt"

    return (
        latest_directory.exists()
        and latest_manifest_file_path.exists()
        and latest_protein_uids_file_path.exists()
    )


def resolve_ncbi_protein_uid_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    search_query: str,
    deduplicate_uids: bool = True,
    sort_uids: bool = True,
    page_size: int = 1000,
    max_retry_attempts: int = 5,
    request_delay_seconds: Optional[float] = None,
    ncbi_email: Optional[str] = None,
    ncbi_api_key: Optional[str] = None,
    update_latest_directory: bool = True,
) -> Dict[str, Any]:
    """
    Resolve the active snapshot payload according to the requested mode.
    """
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    latest_is_available = latest_snapshot_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory,
    )

    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            latest_directory = resolved_snapshot_root_directory / "latest"
            latest_manifest_file_path = latest_directory / "manifest.json"
            latest_protein_uids_file_path = latest_directory / "protein_uids.txt"

            if not latest_directory.exists():
                raise FileNotFoundError(
                    "No latest snapshot directory was found. Run the workflow "
                    "once with snapshot_mode='create_new' before using "
                    "'reuse_latest'."
                )

            if not latest_manifest_file_path.exists():
                raise FileNotFoundError(
                    f"Latest snapshot manifest not found: "
                    f"{latest_manifest_file_path}. Run the workflow once with "
                    f"snapshot_mode='create_new' to create it."
                )

            raise FileNotFoundError(
                f"Latest snapshot protein UID file not found: "
                f"{latest_protein_uids_file_path}. Run the workflow once with "
                f"snapshot_mode='create_new' to create it."
            )

        return load_latest_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory,
        )

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        print("Latest snapshot is available. Reusing frozen snapshot.")
        return load_latest_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory,
        )

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
            "ncbi_email is required when snapshot creation from NCBI is needed."
        )

    fetch_result = fetch_ncbi_protein_uid_snapshot(
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        query=search_query,
        deduplicate_uids=deduplicate_uids,
        sort_uids=sort_uids,
        page_size=page_size,
        max_retry_attempts=max_retry_attempts,
        request_delay_seconds=request_delay_seconds,
    )

    saved_snapshot_directory = save_ncbi_protein_uid_snapshot(
        fetch_result=fetch_result,
        snapshot_root_directory=resolved_snapshot_root_directory,
        update_latest_directory=update_latest_directory,
    )

    return load_snapshot_by_directory(
        snapshot_directory=saved_snapshot_directory,
    )