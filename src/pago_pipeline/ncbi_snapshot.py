from __future__ import annotations

import hashlib
import platform
import re
import shutil
import tempfile
from dataclasses import asdict
from datetime import datetime, timezone
from enum import StrEnum
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, Optional, Union
from collections import Counter

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
    sha256_of_lines,
    write_bytes_atomic,
    write_json_atomic,
)

PathLike = Union[str, Path]
DEFAULT_LOCAL_EXCEL_IDENTIFIER_COLUMN = "Accession"


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


def _normalize_local_protein_identifiers(
    *,
    protein_identifiers: list[Any],
    deduplicate_identifiers: bool,
    sort_identifiers: bool,
) -> list[str]:
    """
    Normalize locally provided protein identifiers for snapshot persistence.
    """
    normalized_identifiers = [
        str(protein_identifier).strip()
        for protein_identifier in protein_identifiers
        if str(protein_identifier).strip()
    ]

    if deduplicate_identifiers:
        seen_identifiers = set()
        deduplicated_identifiers: list[str] = []

        for protein_identifier in normalized_identifiers:
            if protein_identifier not in seen_identifiers:
                seen_identifiers.add(protein_identifier)
                deduplicated_identifiers.append(protein_identifier)

        normalized_identifiers = deduplicated_identifiers

    if sort_identifiers:
        normalized_identifiers = sorted(normalized_identifiers)

    return normalized_identifiers


def _read_excel_identifier_values(
    *,
    excel_file_path: PathLike,
    sheet_name: str | int,
    identifier_column: str,
) -> list[Any]:
    """
    Read one identifier column from a local Excel workbook.
    """
    resolved_excel_file_path = _as_path(excel_file_path)

    if not resolved_excel_file_path.exists():
        raise FileNotFoundError(f"Excel source file not found: {resolved_excel_file_path}")

    if not resolved_excel_file_path.is_file():
        raise ValueError(f"Excel source path is not a file: {resolved_excel_file_path}")

    try:
        import pandas as pd
    except ImportError as error:
        raise ImportError(
            "pandas is required to read local Excel protein identifier sources."
        ) from error

    try:
        excel_dataframe = pd.read_excel(
            resolved_excel_file_path,
            sheet_name=sheet_name,
        )
    except ImportError as error:
        raise ImportError(
            "Reading .xls sources requires the optional dependency xlrd. "
            "Install project dependencies again with requirements.lock.txt."
        ) from error

    if identifier_column not in excel_dataframe.columns:
        available_columns = [str(column) for column in excel_dataframe.columns]
        raise ValueError(
            f"Identifier column {identifier_column!r} was not found in "
            f"{resolved_excel_file_path}. Available columns: {available_columns}."
        )

    identifier_series = excel_dataframe[identifier_column].dropna()
    return identifier_series.tolist()


def fetch_local_excel_protein_identifier_snapshot(
    *,
    excel_file_path: PathLike,
    sheet_name: str | int = 0,
    identifier_column: str = DEFAULT_LOCAL_EXCEL_IDENTIFIER_COLUMN,
    deduplicate_identifiers: bool = True,
    sort_identifiers: bool = True,
) -> NCBIProteinUidFetchResult:
    """
    Build a snapshot-ready protein identifier payload from a local Excel file.

    The downstream XML stage can submit accession.version values to NCBI EFetch
    in the same id list field historically used for numeric protein UIDs.
    """
    resolved_excel_file_path = _as_path(excel_file_path)
    raw_identifier_values = _read_excel_identifier_values(
        excel_file_path=resolved_excel_file_path,
        sheet_name=sheet_name,
        identifier_column=identifier_column,
    )

    normalized_identifiers = _normalize_local_protein_identifiers(
        protein_identifiers=raw_identifier_values,
        deduplicate_identifiers=deduplicate_identifiers,
        sort_identifiers=sort_identifiers,
    )

    if not normalized_identifiers:
        raise ValueError(
            "The local Excel source did not contain any non-empty protein "
            f"identifiers in column {identifier_column!r}."
        )

    protein_identifiers_sha256 = sha256_of_lines(
        text_lines=normalized_identifiers,
        deduplicate_lines_preserving_order=False,
        sort_lines=False,
    )
    excel_file_sha256 = sha256_of_file(input_file_path=resolved_excel_file_path)
    search_query = (
        f"local_excel:{resolved_excel_file_path.name}:"
        f"sheet={sheet_name}:column={identifier_column}"
    )

    return NCBIProteinUidFetchResult(
        database_name="protein",
        search_query=search_query,
        translated_query=None,
        identifier_type="accession.version",
        retrieved_at_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        ncbi_reported_result_count=len(raw_identifier_values),
        protein_uids=normalized_identifiers,
        raw_protein_uid_count=len(raw_identifier_values),
        normalized_protein_uid_count=len(normalized_identifiers),
        deduplicate_uids=deduplicate_identifiers,
        sort_uids=sort_identifiers,
        protein_uids_sha256=protein_identifiers_sha256,
        page_size=0,
        request_delay_seconds=0.0,
        max_retry_attempts=0,
        history_web_env=None,
        history_query_key=None,
        python_version=platform.python_version(),
        biopython_version="not_used_for_local_excel_source",
        source_type="local_excel",
        source_file_path=str(resolved_excel_file_path),
        source_file_sha256=excel_file_sha256,
        source_sheet_name=str(sheet_name),
        source_identifier_column=identifier_column,
    )


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


def _read_xml_batch_payload(
    *,
    xml_batch: NCBIProteinXmlBatchFetchResult,
) -> bytes:
    if xml_batch.xml_payload_bytes:
        return xml_batch.xml_payload_bytes

    if xml_batch.xml_payload_file_path:
        return _as_path(xml_batch.xml_payload_file_path).read_bytes()

    raise ValueError(
        f"XML batch {xml_batch.batch_index} has neither in-memory payload nor "
        "xml_payload_file_path."
    )


def _write_consolidated_xml_payload_streaming(
    *,
    fetch_result: NCBIProteinXmlFetchResult,
    output_file_path: PathLike,
) -> int:
    """
    Consolidate XML batches one at a time to keep memory bounded.
    """
    if not fetch_result.xml_batches:
        raise ValueError("fetch_result.xml_batches must contain at least one batch.")

    resolved_output_file_path = _as_path(output_file_path)
    resolved_output_file_path.parent.mkdir(parents=True, exist_ok=True)

    temporary_output_file = tempfile.NamedTemporaryFile(
        mode="wb",
        delete=False,
        dir=resolved_output_file_path.parent,
        suffix=".tmp",
    )
    temporary_output_file_path = Path(temporary_output_file.name)
    total_record_count = 0
    expected_root_tag: Optional[str] = None

    try:
        with temporary_output_file:
            temporary_output_file.write(b'<?xml version="1.0" encoding="utf-8"?>\n')
            temporary_output_file.write(b"<GBSet>\n")

            for xml_batch in fetch_result.xml_batches:
                xml_payload_bytes = _read_xml_batch_payload(xml_batch=xml_batch)

                try:
                    batch_root = ET.fromstring(xml_payload_bytes)
                except ET.ParseError as error:
                    raise RuntimeError(
                        f"Failed to parse XML batch {xml_batch.batch_index}: {error}"
                    ) from error

                if expected_root_tag is None:
                    expected_root_tag = batch_root.tag
                elif batch_root.tag != expected_root_tag:
                    raise RuntimeError(
                        "XML batch root tag mismatch during consolidation. "
                        f"Expected {expected_root_tag}, got {batch_root.tag} "
                        f"in batch {xml_batch.batch_index}."
                    )

                if batch_root.tag != "GBSet":
                    raise RuntimeError(
                        "XML batch root tag mismatch during consolidation. "
                        f"Expected GBSet, got {batch_root.tag}."
                    )

                batch_children = list(batch_root)
                for child in batch_children:
                    if child.tag != "GBSeq":
                        raise RuntimeError(
                            "XML batch contains unexpected child tag during "
                            f"consolidation: {child.tag}."
                        )

                    temporary_output_file.write(ET.tostring(child, encoding="utf-8"))
                    temporary_output_file.write(b"\n")

                total_record_count += len(batch_children)
                batch_root.clear()

            temporary_output_file.write(b"</GBSet>\n")

        temporary_output_file_path.replace(resolved_output_file_path)
    except Exception:
        if temporary_output_file_path.exists():
            temporary_output_file_path.unlink(missing_ok=True)
        raise

    return total_record_count


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

    record_count = 0
    invalid_child_tags: set[str] = set()
    root_seen = False
    depth = 0

    try:
        for event, element in ET.iterparse(
            resolved_xml_file_path,
            events=("start", "end"),
        ):
            if event == "start":
                if not root_seen:
                    root_seen = True
                    if element.tag != expected_root_tag:
                        raise RuntimeError(
                            "Saved consolidated XML snapshot root tag mismatch. "
                            f"Expected {expected_root_tag}, got {element.tag}."
                        )
                else:
                    depth += 1
                    if depth == 1 and element.tag != expected_record_tag:
                        invalid_child_tags.add(element.tag)
            elif event == "end":
                if depth == 1 and element.tag == expected_record_tag:
                    record_count += 1
                if depth <= 1:
                    element.clear()
                if depth > 0:
                    depth -= 1
    except ET.ParseError as error:
        raise RuntimeError(
            f"Saved consolidated XML snapshot is not well-formed: {error}"
        ) from error

    if not root_seen:
        raise RuntimeError("Saved consolidated XML snapshot is empty.")

    invalid_child_tags = sorted(invalid_child_tags)
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


def _extract_ncbi_protein_uid_from_gbseq_element(
    *,
    gbseq_element: ET.Element,
) -> str:
    """
    Extract one protein UID from a GBSeq record using its GBSeqid fields.
    """
    uid_candidates: set[str] = set()

    for gbseqid_element in gbseq_element.findall(".//GBSeqid"):
        gbseqid_text = (gbseqid_element.text or "").strip()
        if not gbseqid_text:
            continue

        uid_match = re.fullmatch(r"gi\|(\d+)\|?", gbseqid_text)
        if uid_match is not None:
            uid_candidates.add(uid_match.group(1))

    if not uid_candidates:
        raise RuntimeError(
            "Failed to extract a protein UID from GBSeqid fields in one XML record."
        )

    if len(uid_candidates) != 1:
        raise RuntimeError(
            "Found multiple conflicting protein UID candidates in one XML record: "
            f"{sorted(uid_candidates)}."
        )

    return next(iter(uid_candidates))


def _extract_ncbi_protein_identifier_candidates_from_gbseq_element(
    *,
    gbseq_element: ET.Element,
) -> set[str]:
    """
    Extract equivalent protein identifiers from one GBSeq record.
    """
    identifier_candidates: set[str] = set()
    ignored_pipe_tokens = {
        "bbs",
        "dbj",
        "emb",
        "gb",
        "gi",
        "gnl",
        "lcl",
        "pdb",
        "pir",
        "prf",
        "ref",
        "sp",
    }

    for gbseqid_element in gbseq_element.findall(".//GBSeqid"):
        gbseqid_text = (gbseqid_element.text or "").strip()
        if not gbseqid_text:
            continue

        identifier_candidates.add(gbseqid_text)

        uid_match = re.fullmatch(r"gi\|(\d+)\|?", gbseqid_text)
        if uid_match is not None:
            identifier_candidates.add(uid_match.group(1))

        for token in gbseqid_text.split("|"):
            normalized_token = token.strip()
            if (
                normalized_token
                and normalized_token.lower() not in ignored_pipe_tokens
            ):
                identifier_candidates.add(normalized_token)

    for xpath in (
        "./GBSeq_accession-version",
        "./GBSeq_primary-accession",
        ".//GBSeq_other-seqids/GBSeqid",
    ):
        for identifier_element in gbseq_element.findall(xpath):
            identifier_text = (identifier_element.text or "").strip()
            if identifier_text:
                identifier_candidates.add(identifier_text)

    return identifier_candidates


def _validate_xml_record_uids(
    *,
    xml_file_path: PathLike,
    expected_protein_uids: list[str],
    expected_root_tag: str = "GBSet",
    expected_record_tag: str = "GBSeq",
) -> list[str]:
    """
    Validate that the XML record UIDs match the requested protein UIDs.
    """
    resolved_xml_file_path = _as_path(xml_file_path)
    expected_uid_counter = Counter(expected_protein_uids)
    matched_protein_uids: list[str] = []
    unmatched_record_candidates: list[list[str]] = []
    invalid_child_tags: set[str] = set()
    root_seen = False
    depth = 0

    try:
        for event, element in ET.iterparse(
            resolved_xml_file_path,
            events=("start", "end"),
        ):
            if event == "start":
                if not root_seen:
                    root_seen = True
                    if element.tag != expected_root_tag:
                        raise RuntimeError(
                            "Saved consolidated XML snapshot root tag mismatch "
                            "during UID validation. "
                            f"Expected {expected_root_tag}, got {element.tag}."
                        )
                else:
                    depth += 1
                    if depth == 1 and element.tag != expected_record_tag:
                        invalid_child_tags.add(element.tag)
            elif event == "end":
                if depth == 1 and element.tag == expected_record_tag:
                    identifier_candidates = (
                        _extract_ncbi_protein_identifier_candidates_from_gbseq_element(
                            gbseq_element=element,
                        )
                    )
                    matching_expected_identifiers = [
                        expected_identifier
                        for expected_identifier, expected_count
                        in expected_uid_counter.items()
                        if (
                            expected_count > 0
                            and expected_identifier in identifier_candidates
                        )
                    ]

                    if len(matching_expected_identifiers) > 1:
                        raise RuntimeError(
                            "Saved consolidated XML snapshot record matched "
                            "multiple expected protein identifiers: "
                            f"{matching_expected_identifiers}."
                        )

                    if len(matching_expected_identifiers) == 1:
                        matched_identifier = matching_expected_identifiers[0]
                        expected_uid_counter[matched_identifier] -= 1
                        matched_protein_uids.append(matched_identifier)
                    else:
                        unmatched_record_candidates.append(
                            sorted(identifier_candidates)
                        )

                if depth <= 1:
                    element.clear()
                if depth > 0:
                    depth -= 1
    except ET.ParseError as error:
        raise RuntimeError(
            f"Saved consolidated XML snapshot is not well-formed: {error}"
        ) from error

    if not root_seen:
        raise RuntimeError("Saved consolidated XML snapshot is empty.")

    invalid_child_tags = sorted(invalid_child_tags)
    if invalid_child_tags:
        raise RuntimeError(
            "Saved consolidated XML snapshot contains unexpected child tags "
            f"during UID validation: {invalid_child_tags}."
        )

    if any(expected_uid_counter.values()) or unmatched_record_candidates:
        missing_uids = sorted(expected_uid_counter.elements())
        raise RuntimeError(
            "Saved consolidated XML snapshot record UIDs do not match the "
            "expected protein UIDs. "
            f"Missing: {missing_uids[:5]}; "
            f"Unexpected record identifier candidates: "
            f"{unmatched_record_candidates[:5]}."
        )

    return matched_protein_uids

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


def _replace_latest_directory(
    *,
    latest_directory: PathLike,
    files_to_copy: list[tuple[PathLike, str]],
) -> None:
    """
    Replace the convenience latest directory from a fully prepared temp copy.

    This avoids exposing a partially populated latest/ directory if copying one
    of the snapshot artifacts fails midway through the publish step.
    """
    resolved_latest_directory = _as_path(latest_directory)
    resolved_latest_directory.parent.mkdir(parents=True, exist_ok=True)

    temporary_latest_directory = Path(
        tempfile.mkdtemp(
            prefix=f"{resolved_latest_directory.name}_tmp_",
            dir=resolved_latest_directory.parent,
        )
    )

    try:
        for source_file_path, destination_file_name in files_to_copy:
            shutil.copy2(
                _as_path(source_file_path),
                temporary_latest_directory / destination_file_name,
            )

        if resolved_latest_directory.exists():
            shutil.rmtree(resolved_latest_directory)

        temporary_latest_directory.replace(resolved_latest_directory)
    except Exception:
        if temporary_latest_directory.exists():
            shutil.rmtree(temporary_latest_directory, ignore_errors=True)
        raise


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
        "batch_workspace_directory": fetch_result.batch_workspace_directory,
        "socket_idle_timeout_seconds": fetch_result.socket_idle_timeout_seconds,
        "batch_deadline_seconds": fetch_result.batch_deadline_seconds,
        "circuit_breaker_failure_threshold": (
            fetch_result.circuit_breaker_failure_threshold
        ),
        "circuit_breaker_cooldown_seconds": (
            fetch_result.circuit_breaker_cooldown_seconds
        ),
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
    immutable_snapshot_complete = False

    try:
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
        immutable_snapshot_complete = True

        if update_latest_directory:
            _replace_latest_directory(
                latest_directory=resolved_snapshot_root_directory / "latest",
                files_to_copy=[
                    (protein_uids_file_path, "protein_uids.txt"),
                    (manifest_file_path, "manifest.json"),
                ],
            )
    except Exception:
        if not immutable_snapshot_complete and immutable_snapshot_directory.exists():
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

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
    immutable_snapshot_complete = False

    try:
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
                    "xml_payload_file_path": xml_batch.xml_payload_file_path,
                    "response_byte_count": xml_batch.response_byte_count,
                    "open_latency_seconds": xml_batch.open_latency_seconds,
                    "read_latency_seconds": xml_batch.read_latency_seconds,
                    "total_latency_seconds": xml_batch.total_latency_seconds,
                    "average_throughput_bytes_per_second": (
                        xml_batch.average_throughput_bytes_per_second
                    ),
                    "attempt_count": xml_batch.attempt_count,
                }
            )

        expected_record_count = fetch_result.normalized_protein_uid_count
        has_persisted_xml_batches = any(
            xml_batch.xml_payload_file_path for xml_batch in fetch_result.xml_batches
        )

        if has_persisted_xml_batches:
            consolidated_record_count = _write_consolidated_xml_payload_streaming(
                fetch_result=fetch_result,
                output_file_path=consolidated_xml_file_path,
            )
        else:
            consolidated_xml_payload, consolidated_record_count = (
                _build_consolidated_xml_payload(fetch_result=fetch_result)
            )

        if consolidated_record_count != expected_record_count:
            raise RuntimeError(
                "Consolidated XML snapshot record count does not match the "
                "expected protein UID count. "
                f"Expected {expected_record_count}, got "
                f"{consolidated_record_count}."
            )

        if not has_persisted_xml_batches:
            write_bytes_atomic(
                binary_payload=consolidated_xml_payload,
                output_file_path=consolidated_xml_file_path,
            )

        validated_record_count = _validate_saved_consolidated_xml_snapshot(
            xml_file_path=consolidated_xml_file_path,
            expected_record_count=expected_record_count,
        )
        _validate_xml_record_uids(
            xml_file_path=consolidated_xml_file_path,
            expected_protein_uids=protein_uids,
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
        immutable_snapshot_complete = True

        if update_latest_directory:
            _replace_latest_directory(
                latest_directory=resolved_snapshot_root_directory / "latest",
                files_to_copy=[
                    (protein_uids_file_path, "protein_uids.txt"),
                    (manifest_file_path, "manifest.json"),
                    (consolidated_xml_file_path, consolidated_xml_file_name),
                ],
            )
    except Exception:
        if not immutable_snapshot_complete and immutable_snapshot_directory.exists():
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

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


def _validate_loaded_xml_snapshot_payload(
    *,
    snapshot_directory: PathLike,
    manifest_payload: Dict[str, Any],
    protein_uids: list[str],
) -> Path:
    """
    Validate a saved XML snapshot before it is reused as pipeline input.
    """
    resolved_snapshot_directory = _as_path(snapshot_directory)

    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != "ncbi_protein_xml_snapshot":
        raise RuntimeError(
            "Saved XML snapshot manifest artifact_type mismatch. "
            f"Expected 'ncbi_protein_xml_snapshot', got {artifact_type!r}."
        )

    xml_file_name = manifest_payload.get("xml_file_name")
    if not isinstance(xml_file_name, str) or not xml_file_name.strip():
        raise RuntimeError(
            "Saved XML snapshot manifest must define a non-empty xml_file_name."
        )

    resolved_xml_file_path = resolved_snapshot_directory / xml_file_name
    if not resolved_xml_file_path.exists():
        raise FileNotFoundError(
            f"Saved XML snapshot file not found: {resolved_xml_file_path}."
        )

    expected_uid_count = manifest_payload.get("normalized_protein_uid_count")
    if expected_uid_count is not None:
        if not isinstance(expected_uid_count, int):
            raise RuntimeError(
                "Saved XML snapshot manifest normalized_protein_uid_count "
                "must be an integer."
            )

        if len(protein_uids) != expected_uid_count:
            raise RuntimeError(
                "Saved XML snapshot protein UID count mismatch. "
                f"Expected {expected_uid_count}, got {len(protein_uids)}."
            )

    expected_uid_sha256 = manifest_payload.get("protein_uids_sha256")
    if expected_uid_sha256 is not None:
        actual_uid_sha256 = sha256_of_lines(
            text_lines=protein_uids,
            deduplicate_lines_preserving_order=False,
            sort_lines=False,
        )
        if actual_uid_sha256 != expected_uid_sha256:
            raise RuntimeError(
                "Saved XML snapshot protein UID hash mismatch. "
                f"Expected {expected_uid_sha256}, got {actual_uid_sha256}."
            )

    expected_record_count = manifest_payload.get("consolidated_record_count")
    if expected_record_count is not None and not isinstance(expected_record_count, int):
        raise RuntimeError(
            "Saved XML snapshot manifest consolidated_record_count must be an integer."
        )

    validated_record_count = _validate_saved_consolidated_xml_snapshot(
        xml_file_path=resolved_xml_file_path,
        expected_record_count=expected_record_count,
    )
    _validate_xml_record_uids(
        xml_file_path=resolved_xml_file_path,
        expected_protein_uids=protein_uids,
    )

    if (
        expected_uid_count is not None
        and validated_record_count != expected_uid_count
    ):
        raise RuntimeError(
            "Saved XML snapshot record count does not match the saved protein "
            "UID count. "
            f"Expected {expected_uid_count}, got {validated_record_count}."
        )

    expected_xml_file_sha256 = manifest_payload.get("xml_file_sha256")
    if expected_xml_file_sha256 is not None:
        actual_xml_file_sha256 = sha256_of_file(
            input_file_path=resolved_xml_file_path,
        )
        if actual_xml_file_sha256 != expected_xml_file_sha256:
            raise RuntimeError(
                "Saved XML snapshot file hash mismatch. "
                f"Expected {expected_xml_file_sha256}, got {actual_xml_file_sha256}."
            )

    return resolved_xml_file_path


def load_xml_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> Dict[str, Any]:
    """
    Load and validate an XML snapshot directory before reuse.
    """
    snapshot_payload = load_snapshot_by_directory(
        snapshot_directory=snapshot_directory,
    )
    xml_file_path = _validate_loaded_xml_snapshot_payload(
        snapshot_directory=snapshot_directory,
        manifest_payload=snapshot_payload["manifest"],
        protein_uids=snapshot_payload["protein_uids"],
    )
    snapshot_payload["xml_file_path"] = xml_file_path
    return snapshot_payload


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


def load_latest_xml_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> Dict[str, Any]:
    """
    Load and validate the convenience latest XML snapshot copy.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"

    return load_xml_snapshot_by_directory(snapshot_directory=latest_directory)


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


def _local_excel_source_matches_manifest(
    *,
    manifest_payload: Dict[str, Any],
    excel_file_path: PathLike,
    sheet_name: str | int,
    identifier_column: str,
) -> bool:
    """
    Return True when a saved UID snapshot matches the requested Excel source.
    """
    resolved_excel_file_path = _as_path(excel_file_path)

    if not resolved_excel_file_path.exists():
        return False

    return (
        manifest_payload.get("source_type") == "local_excel"
        and manifest_payload.get("source_file_path") == str(resolved_excel_file_path)
        and manifest_payload.get("source_file_sha256")
        == sha256_of_file(input_file_path=resolved_excel_file_path)
        and manifest_payload.get("source_sheet_name") == str(sheet_name)
        and manifest_payload.get("source_identifier_column") == identifier_column
    )


def _source_uid_snapshot_matches_xml_manifest(
    *,
    xml_manifest_payload: Dict[str, Any],
    source_uid_snapshot_payload: Dict[str, Any],
) -> bool:
    """
    Return True when an XML snapshot was built from the given UID snapshot.
    """
    source_uid_manifest = source_uid_snapshot_payload.get("manifest", {})
    source_uid_sha256 = source_uid_manifest.get("protein_uids_sha256")

    if source_uid_sha256 is None:
        return False

    return xml_manifest_payload.get("source_uid_sha256") == source_uid_sha256


def resolve_ncbi_protein_uid_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    search_query: Optional[str] = None,
    local_excel_file_path: Optional[PathLike] = None,
    local_excel_sheet_name: str | int = 0,
    local_excel_identifier_column: str = DEFAULT_LOCAL_EXCEL_IDENTIFIER_COLUMN,
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
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

    When snapshot_mode is 'reuse_latest', only snapshot_mode and
    snapshot_root_directory are required; search_query, ncbi_email,
    and the other creation-only parameters are ignored.
    """
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    use_local_excel_source = local_excel_file_path is not None

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

        latest_snapshot_payload = load_latest_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory,
        )

        if use_local_excel_source and not _local_excel_source_matches_manifest(
            manifest_payload=latest_snapshot_payload["manifest"],
            excel_file_path=local_excel_file_path,
            sheet_name=local_excel_sheet_name,
            identifier_column=local_excel_identifier_column,
        ):
            raise RuntimeError(
                "The latest source UID snapshot does not match the requested "
                "local Excel source."
            )

        return latest_snapshot_payload

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        latest_snapshot_payload = load_latest_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory,
        )

        if not use_local_excel_source or _local_excel_source_matches_manifest(
            manifest_payload=latest_snapshot_payload["manifest"],
            excel_file_path=local_excel_file_path,
            sheet_name=local_excel_sheet_name,
            identifier_column=local_excel_identifier_column,
        ):
            print("Latest snapshot is available. Reusing frozen snapshot.")
            return latest_snapshot_payload

        print(
            "Latest snapshot is available but does not match the requested "
            "local Excel source. Creating a new frozen snapshot."
        )

    if resolved_snapshot_mode not in {
        SnapshotMode.create_new,
        SnapshotMode.reuse_latest_or_create,
    }:
        raise ValueError(
            "Invalid snapshot_mode. Expected one of: "
            "'create_new', 'reuse_latest', 'reuse_latest_or_create'."
        )

    if use_local_excel_source:
        fetch_result = fetch_local_excel_protein_identifier_snapshot(
            excel_file_path=local_excel_file_path,
            sheet_name=local_excel_sheet_name,
            identifier_column=local_excel_identifier_column,
            deduplicate_identifiers=deduplicate_uids,
            sort_identifiers=sort_uids,
        )
    elif not ncbi_email:
        raise ValueError(
            "ncbi_email is required when snapshot creation from NCBI is needed."
        )
    elif not search_query or not search_query.strip():
        raise ValueError(
            "search_query is required when snapshot creation from NCBI is needed."
        )
    else:
        fetch_result = fetch_ncbi_protein_uid_snapshot(
            ncbi_email=ncbi_email,
            ncbi_api_key=ncbi_api_key,
            query=search_query,
            ssl_ca_file=ssl_ca_file,
            ssl_ca_directory=ssl_ca_directory,
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
    source_uid_snapshot_payload: Optional[Dict[str, Any]] = None,
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
    xml_batch_size: int = 50,
    max_retry_attempts: int = 5,
    xml_request_delay_seconds: Optional[float] = None,
    socket_idle_timeout_seconds: float = 90.0,
    batch_deadline_seconds: float = 300.0,
    circuit_breaker_failure_threshold: int = 3,
    circuit_breaker_cooldown_seconds: float = 180.0,
    batch_workspace_root_directory: Optional[PathLike] = None,
    ncbi_email: Optional[str] = None,
    ncbi_api_key: Optional[str] = None,
    update_latest_directory: bool = True,
) -> Dict[str, Any]:
    """
    Resolve the active raw XML snapshot payload according to the requested mode.

    If snapshot creation is needed, the XML retrieval is performed strictly from an
    already-frozen upstream UID snapshot. This function does not create the source
    UID snapshot; that upstream step must have been run previously.
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

        latest_xml_snapshot_payload = load_latest_xml_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory,
        )

        if (
            source_uid_snapshot_payload is not None
            and not _source_uid_snapshot_matches_xml_manifest(
                xml_manifest_payload=latest_xml_snapshot_payload["manifest"],
                source_uid_snapshot_payload=source_uid_snapshot_payload,
            )
        ):
            raise RuntimeError(
                "The latest XML snapshot does not match the requested source "
                "UID snapshot."
            )

        return latest_xml_snapshot_payload

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        latest_xml_snapshot_payload = load_latest_xml_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory,
        )

        if (
            source_uid_snapshot_payload is None
            or _source_uid_snapshot_matches_xml_manifest(
                xml_manifest_payload=latest_xml_snapshot_payload["manifest"],
                source_uid_snapshot_payload=source_uid_snapshot_payload,
            )
        ):
            print("Latest XML snapshot is available. Reusing frozen snapshot.")
            return latest_xml_snapshot_payload

        print(
            "Latest XML snapshot is available but does not match the requested "
            "source UID snapshot. Creating a new frozen XML snapshot."
        )

    if resolved_snapshot_mode not in {
        SnapshotMode.create_new,
        SnapshotMode.reuse_latest_or_create,
    }:
        raise ValueError(
            "Invalid snapshot_mode. Expected one of: "
            "'create_new', 'reuse_latest', 'reuse_latest_or_create'."
        )

    if source_uid_snapshot_payload is None:
        try:
            source_uid_snapshot_payload = resolve_ncbi_protein_uid_snapshot(
                snapshot_mode=SnapshotMode.reuse_latest,
                snapshot_root_directory=source_uid_snapshot_root_directory,
            )
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                "No reusable source UID snapshot was found for the XML workflow. "
                "Create the upstream protein UID snapshot before running the XML "
                "snapshot step."
            ) from exc

    if not ncbi_email:
        raise ValueError(
            "ncbi_email is required when snapshot creation from NCBI is needed."
        )

    protein_uids = source_uid_snapshot_payload["protein_uids"]
    source_uid_snapshot_manifest = source_uid_snapshot_payload["manifest"]
    source_uid_snapshot_manifest_file_path = source_uid_snapshot_payload[
        "manifest_file_path"
    ]
    source_uid_sha256 = source_uid_snapshot_manifest.get(
        "protein_uids_sha256",
        sha256_of_lines(
            text_lines=protein_uids,
            deduplicate_lines_preserving_order=False,
            sort_lines=False,
        ),
    )
    resolved_batch_workspace_root_directory = (
        _as_path(batch_workspace_root_directory)
        if batch_workspace_root_directory is not None
        else resolved_snapshot_root_directory / "workspaces"
    )
    batch_workspace_directory = (
        resolved_batch_workspace_root_directory
        / f"uids_{str(source_uid_sha256)[:16]}__batch_{xml_batch_size}__xml"
    )

    fetch_result = fetch_ncbi_protein_xml_batches(
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        protein_uids=protein_uids,
        ssl_ca_file=ssl_ca_file,
        ssl_ca_directory=ssl_ca_directory,
        batch_size=xml_batch_size,
        max_retry_attempts=max_retry_attempts,
        request_delay_seconds=xml_request_delay_seconds,
        batch_workspace_directory=batch_workspace_directory,
        socket_idle_timeout_seconds=socket_idle_timeout_seconds,
        batch_deadline_seconds=batch_deadline_seconds,
        circuit_breaker_failure_threshold=circuit_breaker_failure_threshold,
        circuit_breaker_cooldown_seconds=circuit_breaker_cooldown_seconds,
    )

    saved_snapshot_directory = save_ncbi_protein_xml_snapshot(
        fetch_result=fetch_result,
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_uid_snapshot_manifest=source_uid_snapshot_manifest,
        source_uid_snapshot_manifest_file_path=source_uid_snapshot_manifest_file_path,
        protein_uids=protein_uids,
        update_latest_directory=update_latest_directory,
    )

    return load_xml_snapshot_by_directory(
        snapshot_directory=saved_snapshot_directory,
    )
