from __future__ import annotations
import hashlib
import shutil
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Literal, Optional, Union
from src.pago_pipeline.ncbi_api import NCBIProteinIdFetchResult, fetch_ncbi_protein_id_snapshot
from src.pago_pipeline.storage import (
    read_json_file,
    read_text_lines_from_file,
    save_ncbi_protein_ids_as_txt,
    write_json_atomic,
)


PathLike = Union[str, Path]


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


def _build_snapshot_manifest(
    *,
    fetch_result: NCBIProteinIdFetchResult,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
) -> Dict[str, Any]:
    """
    Build the manifest payload persisted alongside protein_ids.txt.

    The manifest intentionally duplicates key provenance metadata so the
    snapshot can be audited and reused without depending on live NCBI state.

    Important:
    - immutable_snapshot_directory_name identifies the immutable snapshot folder
      inside snapshots/
    - immutable_snapshot_relative_path preserves the original physical origin
      even when this manifest is later copied into latest/
    """
    manifest_payload = asdict(fetch_result)
    manifest_payload.pop("protein_ids", None)

    manifest_payload["snapshot_format_version"] = "1.0"
    manifest_payload["snapshot_file_name"] = "protein_ids.txt"
    manifest_payload["manifest_file_name"] = "manifest.json"
    manifest_payload["immutable_snapshot_directory_name"] = (
        immutable_snapshot_directory_name
    )
    manifest_payload["immutable_snapshot_relative_path"] = (
        immutable_snapshot_relative_path
    )

    return manifest_payload


def save_ncbi_protein_id_snapshot(
    *,
    fetch_result: NCBIProteinIdFetchResult,
    snapshot_root_directory: PathLike,
    update_latest_directory: bool = True,
) -> Path:
    """
    Persist one immutable local snapshot of an NCBI protein-ID retrieval event.

    Directory structure:
        <snapshot_root_directory>/
            snapshots/
                <timestamp>__q_<hash>/
                    protein_ids.txt
                    manifest.json
            latest/
                protein_ids.txt
                manifest.json

    Important:
    - The directory inside `snapshots/` is immutable and must not already exist.
    - The `latest/` directory is only a convenience copy.
    """
    snapshot_root_directory = Path(snapshot_root_directory)

    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=fetch_result.retrieved_at_utc,
        search_query=fetch_result.search_query,
    )

    immutable_snapshot_directory = (
        snapshot_root_directory / "snapshots" / snapshot_directory_name
    )
    immutable_snapshot_directory.mkdir(parents=True, exist_ok=False)

    immutable_snapshot_relative_path = str(
        Path("snapshots") / snapshot_directory_name
    )

    protein_ids_file_path = immutable_snapshot_directory / "protein_ids.txt"
    manifest_file_path = immutable_snapshot_directory / "manifest.json"

    save_ncbi_protein_ids_as_txt(
        ncbi_protein_id_list=fetch_result.protein_ids,
        output_txt_file_path=protein_ids_file_path,
        deduplicate_ids=False,
        sort_ids=False,
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
        latest_directory = snapshot_root_directory / "latest"
        latest_directory.mkdir(parents=True, exist_ok=True)

        shutil.copy2(protein_ids_file_path, latest_directory / "protein_ids.txt")
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


def load_snapshot_protein_ids(
    *,
    protein_ids_file_path: PathLike,
) -> list[str]:
    """
    Load protein IDs from an explicit snapshot protein_ids.txt path.
    """
    return read_text_lines_from_file(input_file_path=protein_ids_file_path)


def load_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> Dict[str, Any]:
    """
    Load both manifest metadata and protein IDs from a snapshot directory.

    Returns a dictionary with:
    - manifest
    - protein_ids
    - snapshot_directory
    - manifest_file_path
    - protein_ids_file_path
    """
    snapshot_directory = Path(snapshot_directory)

    manifest_file_path = snapshot_directory / "manifest.json"
    protein_ids_file_path = snapshot_directory / "protein_ids.txt"

    manifest_payload = load_snapshot_manifest(
        manifest_file_path=manifest_file_path,
    )
    protein_ids = load_snapshot_protein_ids(
        protein_ids_file_path=protein_ids_file_path,
    )

    return {
        "snapshot_directory": snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "protein_ids_file_path": protein_ids_file_path,
        "manifest": manifest_payload,
        "protein_ids": protein_ids,
    }


def load_latest_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> Dict[str, Any]:
    """
    Load the convenience 'latest' snapshot copy.
    """
    snapshot_root_directory = Path(snapshot_root_directory)
    latest_directory = snapshot_root_directory / "latest"

    return load_snapshot_by_directory(snapshot_directory=latest_directory)


def get_latest_snapshot_manifest_path(
    *,
    snapshot_root_directory: PathLike,
) -> Path:
    """
    Return the manifest path for the convenience 'latest' snapshot copy.
    """
    snapshot_root_directory = Path(snapshot_root_directory)
    return snapshot_root_directory / "latest" / "manifest.json"


def get_latest_snapshot_protein_ids_path(
    *,
    snapshot_root_directory: PathLike,
) -> Path:
    """
    Return the protein_ids.txt path for the convenience 'latest' snapshot copy.
    """
    snapshot_root_directory = Path(snapshot_root_directory)
    return snapshot_root_directory / "latest" / "protein_ids.txt"


def get_snapshot_manifest_path(
    *,
    snapshot_directory: PathLike,
) -> Path:
    """
    Return the manifest path inside a specific immutable snapshot directory.
    """
    snapshot_directory = Path(snapshot_directory)
    return snapshot_directory / "manifest.json"


def get_snapshot_protein_ids_path(
    *,
    snapshot_directory: PathLike,
) -> Path:
    """
    Return the protein_ids.txt path inside a specific immutable snapshot directory.
    """
    snapshot_directory = Path(snapshot_directory)
    return snapshot_directory / "protein_ids.txt"


def list_saved_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    """
    List saved immutable snapshot directories in lexical order.

    Because the directory name begins with the UTC timestamp, lexical order
    is also chronological order.
    """
    snapshot_root_directory = Path(snapshot_root_directory)
    snapshots_directory = snapshot_root_directory / "snapshots"

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

from enum import StrEnum

class SnapshotMode(StrEnum):
    create_new = "create_new"
    reuse_latest = "reuse_latest"
    reuse_latest_or_create = "reuse_latest_or_create"


def latest_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
) -> bool:
    """
    Return True only if the convenience latest snapshot copy is complete.

    A complete latest snapshot requires:
    - latest/
    - latest/manifest.json
    - latest/protein_ids.txt
    """
    snapshot_root_directory = Path(snapshot_root_directory)

    latest_directory = snapshot_root_directory / "latest"
    latest_manifest_file_path = latest_directory / "manifest.json"
    latest_protein_ids_file_path = latest_directory / "protein_ids.txt"

    return (
        latest_directory.exists()
        and latest_manifest_file_path.exists()
        and latest_protein_ids_file_path.exists()
    )


def resolve_ncbi_protein_id_snapshot(
    *,
    snapshot_mode: SnapshotMode,
    snapshot_root_directory: PathLike,
    search_query: str,
    deduplicate_ids: bool = True,
    sort_ids: bool = True,
    page_size: int = 1000,
    max_retry_attempts: int = 5,
    request_delay_seconds: Optional[float] = None,
    ncbi_email: Optional[str] = None,
    ncbi_api_key: Optional[str] = None,
    update_latest_directory: bool = True,
) -> Dict[str, Any]:
    """
    Resolve the active snapshot payload according to the requested mode.

    Modes:
    - create_new:
        Always query NCBI, save a new immutable snapshot, and return it.
    - reuse_latest:
        Reuse the latest convenience copy only. Fail if it is incomplete.
    - reuse_latest_or_create:
        Reuse the latest convenience copy if complete; otherwise create a new
        snapshot from NCBI.

    Returns:
        Same payload structure as load_snapshot_by_directory(...) and
        load_latest_snapshot(...), with keys:
        - snapshot_directory
        - manifest_file_path
        - protein_ids_file_path
        - manifest
        - protein_ids
    """
    snapshot_root_directory = Path(snapshot_root_directory)

    latest_is_available = latest_snapshot_is_available(
        snapshot_root_directory=snapshot_root_directory,
    )

    if snapshot_mode == "reuse_latest":
        if not latest_is_available:
            latest_directory = snapshot_root_directory / "latest"
            latest_manifest_file_path = latest_directory / "manifest.json"
            latest_protein_ids_file_path = latest_directory / "protein_ids.txt"

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
                f"Latest snapshot protein ID file not found: "
                f"{latest_protein_ids_file_path}. Run the workflow once with "
                f"snapshot_mode='create_new' to create it."
            )

        return load_latest_snapshot(
            snapshot_root_directory=snapshot_root_directory,
        )

    if snapshot_mode == "reuse_latest_or_create" and latest_is_available:
        print("Latest snapshot is available. Reusing frozen snapshot.")
        return load_latest_snapshot(
            snapshot_root_directory=snapshot_root_directory,
        )

    if snapshot_mode not in {"create_new", "reuse_latest_or_create"}:
        raise ValueError(
            "Invalid snapshot_mode. Expected one of: "
            "'create_new', 'reuse_latest', 'reuse_latest_or_create'."
        )

    if not ncbi_email:
        raise ValueError(
            "ncbi_email is required when snapshot creation from NCBI is needed."
        )

    fetch_result = fetch_ncbi_protein_id_snapshot(
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        query=search_query,
        deduplicate_ids=deduplicate_ids,
        sort_ids=sort_ids,
        page_size=page_size,
        max_retry_attempts=max_retry_attempts,
        request_delay_seconds=request_delay_seconds,
    )

    saved_snapshot_directory = save_ncbi_protein_id_snapshot(
        fetch_result=fetch_result,
        snapshot_root_directory=snapshot_root_directory,
        update_latest_directory=update_latest_directory,
    )

    return load_snapshot_by_directory(
        snapshot_directory=saved_snapshot_directory,
    )