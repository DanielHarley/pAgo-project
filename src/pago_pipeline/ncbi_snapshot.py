from __future__ import annotations

import hashlib
import shutil
from dataclasses import asdict
from enum import StrEnum
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, Optional, Union

from src.pago_pipeline.ncbi_api import (
    NCBIProteinUidFetchResult,
    NCBIProteinXmlFetchResult,
    fetch_ncbi_protein_uid_snapshot,
    fetch_ncbi_protein_xml_batches,
)
from src.pago_pipeline.storage import (
    read_json_file,
    read_text_lines_from_file,
    save_ncbi_protein_uids_as_txt,
    sha256_of_file,
    write_bytes_atomic,
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


def _build_consolidated_xml_payload(
    *,
    fetch_result: NCBIProteinXmlFetchResult,
) -> tuple[bytes, int]:
    """
    Merge all XML batches into one valid XML document with a single root.

    Returns:
        tuple[bytes, int]:
            - consolidated XML payload as UTF-8 bytes with XML declaration
            - total number of child records appended under the root
    """
    if not fetch_result.xml_batches:
        raise ValueError("fetch_result.xml_batches must contain at least one batch.")

    consolidated_root: Optional[ET.Element] = None
    expected_root_tag: Optional[str] = None
    total_record_count = 0

    for xml_batch in fetch_result.xml_batches:
        try:
            batch_root = ET.fromstring(xml_batch.xml_payload_bytes)
        except ET.ParseError as error:
            raise RuntimeError(
                f"Failed to parse XML batch {xml_batch.batch_index}: {error}"
            ) from error

        if consolidated_root is None:
            consolidated_root = ET.Element(batch_root.tag, batch_root.attrib)
            expected_root_tag = batch_root.tag
        else:
            if batch_root.tag != expected_root_tag:
                raise RuntimeError(
                    "XML batch root tag mismatch during consolidation. "
                    f"Expected {expected_root_tag}, got {batch_root.tag} "
                    f"in batch {xml_batch.batch_index}."
                )

        batch_children = list(batch_root)
        for child in batch_children:
            consolidated_root.append(child)

        total_record_count += len(batch_children)

    if consolidated_root is None:
        raise RuntimeError("Failed to create consolidated XML root.")

    consolidated_xml_payload = ET.tostring(
        consolidated_root,
        encoding="utf-8",
        xml_declaration=True,
    )

    return consolidated_xml_payload, total_record_count


def _validate_saved_consolidated_xml_snapshot(
    *,
    xml_file_path: PathLike,
    expected_root_tag: str = "GBSet",
    expected_record_tag: str = "GBSeq",
    expected_record_count: Optional[int] = None,
) -> int:
    """
    Parse a saved consolidated XML snapshot and validate basic structural invariants.

    Returns:
        int:
            Number of direct child records found under the root element.
    """
    resolved_xml_file_path = _as_path(xml_file_path)

    try:
        parsed_tree = ET.parse(resolved_xml_file_path)
    except ET.ParseError as error:
        raise RuntimeError(
            f"Saved consolidated XML snapshot is not well-formed: {error}"
        ) from error

    root_element = parsed_tree.getroot()
    if root_element.tag != expected_root_tag:
        raise RuntimeError(
            "Saved consolidated XML snapshot root tag mismatch. "
            f"Expected {expected_root_tag}, got {root_element.tag}."
        )

    child_elements = list(root_element)
    record_count = len(child_elements)

    invalid_child_tags = sorted(
        {child.tag for child in child_elements if child.tag != expected_record_tag}
    )
    if invalid_child_tags:
        raise RuntimeError(
            "Saved consolidated XML snapshot contains unexpected child tags "
            f"under {expected_root_tag}: {invalid_child_tags}."
        )

    if (
        expected_record_count is not None
        and record_count != expected_record_count
    ):
        raise RuntimeError(
            "Saved consolidated XML snapshot record count mismatch. "
            f"Expected {expected_record_count}, got {record_count}."
        )

    return record_count


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


def _build_xml_snapshot_manifest(
    *,
    fetch_result: NCBIProteinXmlFetchResult,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    source_uid_snapshot_manifest: Dict[str, Any],
    source_uid_snapshot_manifest_sha256: str,
    consolidated_xml_file_name: str,
    consolidated_xml_file_sha256: str,
    consolidated_record_count: int,
    batch_metadata_records: list[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Build the manifest payload persisted alongside one consolidated XML file.
    """
    return {
        "snapshot_format_version": "1.0",
        "artifact_type": "ncbi_protein_xml_snapshot",
        "database_name": fetch_result.database_name,
        "identifier_type": fetch_result.identifier_type,
        "search_query": source_uid_snapshot_manifest["search_query"],
        "translated_query": source_uid_snapshot_manifest.get("translated_query"),
        "retrieved_at_utc": fetch_result.retrieved_at_utc,
        "retmode": fetch_result.retmode,
        "requested_protein_uid_count": fetch_result.requested_protein_uid_count,
        "normalized_protein_uid_count": fetch_result.normalized_protein_uid_count,
        "protein_uids_sha256": fetch_result.protein_uids_sha256,
        "batch_size": fetch_result.batch_size,
        "batch_count": fetch_result.batch_count,
        "request_delay_seconds": fetch_result.request_delay_seconds,
        "max_retry_attempts": fetch_result.max_retry_attempts,
        "python_version": fetch_result.python_version,
        "biopython_version": fetch_result.biopython_version,
        "manifest_file_name": "manifest.json",
        "protein_uids_file_name": "protein_uids.txt",
        "xml_file_name": consolidated_xml_file_name,
        "xml_file_sha256": consolidated_xml_file_sha256,
        "consolidated_record_count": consolidated_record_count,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "source_uid_snapshot_relative_path": source_uid_snapshot_manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_uid_snapshot_directory_name": source_uid_snapshot_manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_uid_snapshot_manifest_sha256": source_uid_snapshot_manifest_sha256,
        "source_uid_count": source_uid_snapshot_manifest.get(
            "normalized_protein_uid_count"
        ),
        "source_uid_sha256": source_uid_snapshot_manifest.get("protein_uids_sha256"),
        "batches": batch_metadata_records,
    }


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


def save_ncbi_protein_xml_snapshot(
    *,
    fetch_result: NCBIProteinXmlFetchResult,
    snapshot_root_directory: PathLike,
    source_uid_snapshot_manifest: Dict[str, Any],
    source_uid_snapshot_manifest_file_path: PathLike,
    protein_uids: list[str],
    update_latest_directory: bool = True,
) -> Path:
    """
    Persist one immutable local snapshot of consolidated NCBI protein XML.

    XML is fetched in batches upstream for robustness, then merged into one
    valid XML document with a single root and persisted as one file.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_source_uid_snapshot_manifest_file_path = _as_path(
        source_uid_snapshot_manifest_file_path
    )

    if not protein_uids:
        raise ValueError("protein_uids must contain at least one UID.")

    if fetch_result.normalized_protein_uid_count != len(protein_uids):
        raise ValueError(
            "The XML fetch result UID count does not match the provided "
            "protein_uids list length."
        )

    source_uid_snapshot_manifest_sha256 = sha256_of_file(
        input_file_path=resolved_source_uid_snapshot_manifest_file_path,
    )

    source_uid_sha256 = source_uid_snapshot_manifest.get("protein_uids_sha256")
    if (
        source_uid_sha256 is not None
        and source_uid_sha256 != fetch_result.protein_uids_sha256
    ):
        raise ValueError(
            "The XML fetch result does not match the source UID snapshot. "
            "protein_uids_sha256 values differ."
        )

    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=fetch_result.retrieved_at_utc,
        search_query=source_uid_snapshot_manifest["search_query"],
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
    consolidated_xml_file_name = "protein_records.xml"
    consolidated_xml_file_path = (
        immutable_snapshot_directory / consolidated_xml_file_name
    )

    save_ncbi_protein_uids_as_txt(
        ncbi_protein_uid_list=protein_uids,
        output_txt_file_path=protein_uids_file_path,
        deduplicate_uids=False,
        sort_uids=False,
    )

    batch_metadata_records: list[Dict[str, Any]] = []
    for xml_batch in fetch_result.xml_batches:
        batch_metadata_records.append(
            {
                "batch_index": xml_batch.batch_index,
                "batch_start_index": xml_batch.batch_start_index,
                "batch_end_index": xml_batch.batch_end_index,
                "protein_uid_count": xml_batch.protein_uid_count,
                "xml_payload_sha256": xml_batch.xml_payload_sha256,
            }
        )

    consolidated_xml_payload, consolidated_record_count = (
        _build_consolidated_xml_payload(fetch_result=fetch_result)
    )

    write_bytes_atomic(
        binary_payload=consolidated_xml_payload,
        output_file_path=consolidated_xml_file_path,
    )

    validated_record_count = _validate_saved_consolidated_xml_snapshot(
        xml_file_path=consolidated_xml_file_path,
        expected_record_count=consolidated_record_count,
    )

    consolidated_xml_file_sha256 = sha256_of_file(
        input_file_path=consolidated_xml_file_path,
    )

    manifest_payload = _build_xml_snapshot_manifest(
        fetch_result=fetch_result,
        immutable_snapshot_directory_name=snapshot_directory_name,
        immutable_snapshot_relative_path=immutable_snapshot_relative_path,
        source_uid_snapshot_manifest=source_uid_snapshot_manifest,
        source_uid_snapshot_manifest_sha256=source_uid_snapshot_manifest_sha256,
        consolidated_xml_file_name=consolidated_xml_file_name,
        consolidated_xml_file_sha256=consolidated_xml_file_sha256,
        consolidated_record_count=validated_record_count,
        batch_metadata_records=batch_metadata_records,
    )

    write_json_atomic(
        payload=manifest_payload,
        output_file_path=manifest_file_path,
    )

    if update_latest_directory:
        latest_directory = resolved_snapshot_root_directory / "latest"

        if latest_directory.exists():
            shutil.rmtree(latest_directory)

        latest_directory.mkdir(parents=True, exist_ok=False)

        shutil.copy2(protein_uids_file_path, latest_directory / "protein_uids.txt")
        shutil.copy2(manifest_file_path, latest_directory / "manifest.json")
        shutil.copy2(
            consolidated_xml_file_path,
            latest_directory / consolidated_xml_file_name,
        )

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


def get_latest_snapshot_xml_file_path(
    *,
    snapshot_root_directory: PathLike,
) -> Path:
    """
    Return the consolidated XML file path for the convenience latest snapshot copy.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return resolved_snapshot_root_directory / "latest" / "protein_records.xml"


def get_snapshot_xml_file_path(
    *,
    snapshot_directory: PathLike,
) -> Path:
    """
    Return the consolidated XML file path inside a specific immutable snapshot directory.
    """
    resolved_snapshot_directory = _as_path(snapshot_directory)
    return resolved_snapshot_directory / "protein_records.xml"


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


def latest_xml_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
) -> bool:
    """
    Return True only if the convenience latest XML snapshot copy is complete.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    latest_directory = resolved_snapshot_root_directory / "latest"
    latest_manifest_file_path = latest_directory / "manifest.json"
    latest_protein_uids_file_path = latest_directory / "protein_uids.txt"
    latest_xml_file_path = latest_directory / "protein_records.xml"

    return (
        latest_directory.exists()
        and latest_manifest_file_path.exists()
        and latest_protein_uids_file_path.exists()
        and latest_xml_file_path.exists()
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


def resolve_ncbi_protein_xml_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    source_uid_snapshot_root_directory: PathLike,
    search_query: str,
    source_uid_snapshot_mode: SnapshotMode | str = SnapshotMode.reuse_latest_or_create,
    deduplicate_uids: bool = True,
    sort_uids: bool = True,
    uid_page_size: int = 1000,
    xml_batch_size: int = 100,
    max_retry_attempts: int = 5,
    uid_request_delay_seconds: Optional[float] = None,
    xml_request_delay_seconds: Optional[float] = None,
    ncbi_email: Optional[str] = None,
    ncbi_api_key: Optional[str] = None,
    update_latest_directory: bool = True,
) -> Dict[str, Any]:
    """
    Resolve the active raw XML snapshot payload according to the requested mode.

    If snapshot creation is needed, the XML retrieval is performed from an upstream
    frozen UID snapshot.
    """
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    latest_is_available = latest_xml_snapshot_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory,
    )

    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            latest_directory = resolved_snapshot_root_directory / "latest"
            latest_manifest_file_path = latest_directory / "manifest.json"
            latest_protein_uids_file_path = latest_directory / "protein_uids.txt"
            latest_xml_file_path = latest_directory / "protein_records.xml"

            if not latest_directory.exists():
                raise FileNotFoundError(
                    "No latest XML snapshot directory was found. Run the workflow "
                    "once with snapshot_mode='create_new' before using "
                    "'reuse_latest'."
                )

            if not latest_manifest_file_path.exists():
                raise FileNotFoundError(
                    f"Latest XML snapshot manifest not found: "
                    f"{latest_manifest_file_path}. Run the workflow once with "
                    f"snapshot_mode='create_new' to create it."
                )

            if not latest_protein_uids_file_path.exists():
                raise FileNotFoundError(
                    f"Latest XML snapshot protein UID file not found: "
                    f"{latest_protein_uids_file_path}. Run the workflow once with "
                    f"snapshot_mode='create_new' to create it."
                )

            raise FileNotFoundError(
                f"Latest XML snapshot consolidated XML file not found: "
                f"{latest_xml_file_path}. Run the workflow once with "
                f"snapshot_mode='create_new' to create it."
            )

        return load_latest_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory,
        )

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        print("Latest XML snapshot is available. Reusing frozen snapshot.")
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

    source_uid_snapshot_payload = resolve_ncbi_protein_uid_snapshot(
        snapshot_mode=source_uid_snapshot_mode,
        snapshot_root_directory=source_uid_snapshot_root_directory,
        search_query=search_query,
        deduplicate_uids=deduplicate_uids,
        sort_uids=sort_uids,
        page_size=uid_page_size,
        max_retry_attempts=max_retry_attempts,
        request_delay_seconds=uid_request_delay_seconds,
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        update_latest_directory=update_latest_directory,
    )

    protein_uids = source_uid_snapshot_payload["protein_uids"]
    source_uid_snapshot_manifest = source_uid_snapshot_payload["manifest"]
    source_uid_snapshot_manifest_file_path = source_uid_snapshot_payload[
        "manifest_file_path"
    ]

    fetch_result = fetch_ncbi_protein_xml_batches(
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        protein_uids=protein_uids,
        batch_size=xml_batch_size,
        max_retry_attempts=max_retry_attempts,
        request_delay_seconds=xml_request_delay_seconds,
    )

    saved_snapshot_directory = save_ncbi_protein_xml_snapshot(
        fetch_result=fetch_result,
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_uid_snapshot_manifest=source_uid_snapshot_manifest,
        source_uid_snapshot_manifest_file_path=source_uid_snapshot_manifest_file_path,
        protein_uids=protein_uids,
        update_latest_directory=update_latest_directory,
    )

    return load_snapshot_by_directory(
        snapshot_directory=saved_snapshot_directory,
    )