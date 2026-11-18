from __future__ import annotations

import hashlib
import json
import re
import shutil
import tempfile
from dataclasses import asdict
from datetime import datetime, timezone
from enum import StrEnum
import xml.etree.ElementTree as ET
from pathlib import Path, PureWindowsPath
from typing import Any, Dict, Optional, Union
from collections import Counter

from src.pago_pipeline.ncbi_api import (
    NCBIProteinIdentifierType,
    NCBIProteinUidFetchResult,
    NCBIProteinXmlFetchResult,
    _normalize_protein_identifier_list,
    coerce_ncbi_protein_identifier_type,
    fetch_ncbi_protein_uid_snapshot,
    fetch_ncbi_protein_xml_batches,
    get_default_ncbi_xml_circuit_breaker,
    validate_ncbi_protein_identifiers,
)
from src.pago_pipeline.storage import (
    read_json_file,
    read_text_lines_from_file,
    sha256_of_file,
    sha256_of_lines,
    write_bytes_atomic,
    write_json_atomic,
)

PathLike = Union[str, Path]
NCBI_SNAPSHOT_FORMAT_VERSION = "1.1"
NCBI_PROTEIN_IDENTIFIERS_FILE_NAME = "protein_identifiers.txt"
NCBI_LEGACY_PROTEIN_UIDS_FILE_NAME = "protein_uids.txt"


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


def _manifest_field_value(
    manifest_payload: Dict[str, Any],
    *field_names: str,
    default: Any = None,
) -> Any:
    """Select the first present field, without bypassing invalid generic values."""
    for field_name in field_names:
        if field_name in manifest_payload:
            return manifest_payload[field_name]
    return default


def _snapshot_format_version(manifest_payload: Dict[str, Any]) -> str:
    version = manifest_payload.get("snapshot_format_version", "1.0")
    if version not in ("1.0", NCBI_SNAPSHOT_FORMAT_VERSION):
        raise ValueError(f"Unsupported snapshot_format_version: {version!r}.")
    return version


def _identifier_type_from_snapshot_manifest(
    manifest_payload: Dict[str, Any],
) -> NCBIProteinIdentifierType:
    """Read legacy metadata without changing it; never use this default for writing."""
    version = _snapshot_format_version(manifest_payload)
    if "identifier_type" not in manifest_payload:
        if version == "1.0":
            return NCBIProteinIdentifierType.uid
        raise ValueError(
            "Snapshot manifest identifier_type is required for "
            f"snapshot_format_version={version!r}."
        )
    if manifest_payload["identifier_type"] is None:
        raise ValueError("Snapshot manifest identifier_type must not be null.")
    return coerce_ncbi_protein_identifier_type(manifest_payload["identifier_type"])


def _identifier_type_for_new_snapshot(
    identifier_type: NCBIProteinIdentifierType | str,
) -> NCBIProteinIdentifierType:
    """Require an explicit valid domain when creating a new snapshot manifest."""
    if identifier_type is None:
        raise ValueError("identifier_type must be explicit when creating a new NCBI snapshot.")
    return coerce_ncbi_protein_identifier_type(identifier_type)


def _validate_xml_fetch_identifier_type(
    *,
    fetch_result: NCBIProteinXmlFetchResult,
    expected_identifier_type: NCBIProteinIdentifierType,
) -> NCBIProteinIdentifierType:
    """Require an explicitly typed fetch matching its frozen source domain."""
    fetch_identifier_type = _identifier_type_for_new_snapshot(fetch_result.identifier_type)
    if fetch_identifier_type != expected_identifier_type:
        raise ValueError(
            "The XML fetch result identifier_type does not match the source "
            "snapshot identifier_type. "
            f"Expected {expected_identifier_type.value}, got {fetch_identifier_type.value}."
        )
    return fetch_identifier_type


def _validate_source_identifier_list(
    *,
    protein_identifiers: list[str],
    source_manifest: Dict[str, Any],
) -> str:
    """Check the actual list against the selected source hash and count."""
    _validate_snapshot_integrity_fields(source_manifest)
    hash_field = (
        "protein_identifiers_sha256"
        if "protein_identifiers_sha256" in source_manifest
        else "protein_uids_sha256"
    )
    if hash_field not in source_manifest:
        raise ValueError(
            "Source identifier snapshot manifest must define "
            "protein_identifiers_sha256 or protein_uids_sha256."
        )
    expected_sha256 = source_manifest[hash_field]
    if (
        not isinstance(expected_sha256, str)
        or re.fullmatch(r"[0-9a-fA-F]{64}", expected_sha256) is None
    ):
        raise ValueError(
            f"Source identifier snapshot manifest {hash_field} "
            "must be a SHA-256 hash with 64 hexadecimal characters."
        )

    actual_sha256 = sha256_of_lines(
        text_lines=protein_identifiers,
        deduplicate_lines_preserving_order=False,
        sort_lines=False,
    )
    if actual_sha256 != expected_sha256:
        raise ValueError(
            "The protein identifier list does not match the source identifier "
            "snapshot. protein_identifiers_sha256 values differ."
        )

    expected_count = _manifest_field_value(
        source_manifest,
        "normalized_protein_identifier_count",
        "protein_identifier_count",
        "normalized_protein_uid_count",
    )
    if type(expected_count) is not int:
        raise ValueError(
            "Source identifier snapshot manifest must define an integer "
            "protein identifier count."
        )
    if len(protein_identifiers) != expected_count:
        raise ValueError(
            "Source identifier snapshot protein identifier count mismatch. "
            f"Expected {expected_count}, got {len(protein_identifiers)}."
        )

    return actual_sha256


def _validate_snapshot_integrity_fields(
    manifest_payload: Dict[str, Any],
    *,
    include_xml: bool = False,
) -> None:
    """Require canonical integrity metadata for 1.1; retain legacy fallbacks."""
    if _snapshot_format_version(manifest_payload) != NCBI_SNAPSHOT_FORMAT_VERSION:
        return

    hash_fields = ["protein_identifiers_sha256"]
    count_fields = ["normalized_protein_identifier_count"]
    if include_xml:
        hash_fields.append("xml_file_sha256")
        count_fields.append("consolidated_record_count")

    for field_name in (*hash_fields, *count_fields):
        if field_name not in manifest_payload:
            raise RuntimeError(
                f"Saved snapshot manifest {field_name} is required for "
                "snapshot_format_version=1.1."
            )

    for field_name in hash_fields:
        value = manifest_payload[field_name]
        if not isinstance(value, str) or not value.strip():
            raise RuntimeError(
                f"Saved snapshot manifest {field_name} must be a non-empty string."
            )
        if re.fullmatch(r"[0-9a-f]{64}", value) is None:
            raise RuntimeError(
                f"Saved snapshot manifest {field_name} must be a SHA-256 hash "
                "with 64 lowercase hexadecimal characters."
            )

    for field_name in count_fields:
        value = manifest_payload[field_name]
        if type(value) is not int:
            raise RuntimeError(f"Saved snapshot manifest {field_name} must be an integer.")
        if value < 0:
            raise RuntimeError(f"Saved snapshot manifest {field_name} must be non-negative.")

    if (
        include_xml
        and manifest_payload["consolidated_record_count"]
        != manifest_payload["normalized_protein_identifier_count"]
    ):
        raise RuntimeError(
            "Saved XML snapshot record count mismatch. "
            f"Expected {manifest_payload['normalized_protein_identifier_count']}, "
            f"got {manifest_payload['consolidated_record_count']}."
        )


def _validate_saved_identifier_list(
    *,
    manifest_payload: Dict[str, Any],
    protein_identifiers: list[str],
) -> Optional[int]:
    """Validate saved identifiers using canonical 1.1 or optional legacy metadata."""
    _validate_snapshot_integrity_fields(manifest_payload)
    count_field = next(
        (
            field_name for field_name in (
                "normalized_protein_identifier_count",
                "protein_identifier_count",
                "normalized_protein_uid_count",
            )
            if field_name in manifest_payload
        ),
        None,
    )
    expected_count = manifest_payload[count_field] if count_field else None
    if count_field is not None:
        if type(expected_count) is not int:
            raise RuntimeError(f"Saved snapshot manifest {count_field} must be an integer.")
        if len(protein_identifiers) != expected_count:
            raise RuntimeError(
                "Saved snapshot protein identifier count mismatch. "
                f"Expected {expected_count}, got {len(protein_identifiers)}."
            )

    hash_field = (
        "protein_identifiers_sha256"
        if "protein_identifiers_sha256" in manifest_payload
        else "protein_uids_sha256"
    )
    if hash_field in manifest_payload:
        expected_sha256 = manifest_payload[hash_field]
        if not isinstance(expected_sha256, str) or not expected_sha256.strip():
            raise RuntimeError(f"Saved snapshot manifest {hash_field} must be a non-empty string.")
        actual_sha256 = sha256_of_lines(
            text_lines=protein_identifiers,
            deduplicate_lines_preserving_order=False,
            sort_lines=False,
        )
        if actual_sha256 != expected_sha256:
            raise RuntimeError(
                "Saved snapshot protein identifier hash mismatch. "
                f"Expected {expected_sha256}, got {actual_sha256}."
            )
    if manifest_payload.get("artifact_type") == "ncbi_protein_xml_snapshot":
        _validate_snapshot_integrity_fields(manifest_payload, include_xml=True)
    return expected_count


def _snapshot_member_path(snapshot_directory: Path, relative_path: str) -> Path:
    """Resolve a portable relative path confined to its snapshot directory."""
    member = PureWindowsPath(relative_path)
    if (
        member.anchor or not member.parts or ".." in member.parts
        or ":" in relative_path or "\0" in relative_path
    ):
        raise RuntimeError(f"Snapshot path must remain inside the snapshot directory: {relative_path!r}.")
    path = snapshot_directory.joinpath(*member.parts)
    if not path.resolve().is_relative_to(snapshot_directory.resolve()):
        raise RuntimeError(f"Snapshot path must remain inside the snapshot directory: {relative_path!r}.")
    return path


def _protein_identifier_file_path_from_manifest(
    snapshot_directory: Path,
    manifest_payload: Dict[str, Any],
) -> Path:
    """Resolve the canonical identifier file, then the legacy names if absent."""
    default_file_name = (
        NCBI_LEGACY_PROTEIN_UIDS_FILE_NAME
        if _snapshot_format_version(manifest_payload) == "1.0"
        else NCBI_PROTEIN_IDENTIFIERS_FILE_NAME
    )
    file_name = _manifest_field_value(
        manifest_payload,
        "protein_identifiers_file_name", "protein_uids_file_name", "snapshot_file_name",
        default=default_file_name,
    )
    if not isinstance(file_name, str) or not file_name.strip():
        raise RuntimeError(
            "Snapshot manifest must define a non-empty protein identifier file name."
        )
    return _snapshot_member_path(snapshot_directory, file_name)


def _xml_file_path_from_manifest(
    snapshot_directory: Path,
    manifest_payload: Dict[str, Any],
) -> Path:
    xml_file_name = manifest_payload.get("xml_file_name")
    if not isinstance(xml_file_name, str) or not xml_file_name.strip():
        raise RuntimeError(
            "Saved XML snapshot manifest must define a non-empty xml_file_name."
        )
    return _snapshot_member_path(snapshot_directory, xml_file_name)


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
    timestamp_pattern = (
        r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}"
        r"(?:\.[0-9]{1,6})?(?:Z|[+-](?:[01][0-9]|2[0-3]):[0-5][0-9])"
    )
    if not isinstance(utc_timestamp, str) or re.fullmatch(timestamp_pattern, utc_timestamp) is None:
        raise ValueError("retrieved_at_utc must be an ISO 8601 timestamp with an explicit UTC offset.")
    try:
        timestamp = datetime.fromisoformat(utc_timestamp).astimezone(timezone.utc)
    except (ValueError, OverflowError) as error:
        raise ValueError(f"Invalid retrieved_at_utc timestamp: {utc_timestamp!r}.") from error
    return timestamp.isoformat().removesuffix("+00:00").replace(":", "-") + "Z"


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


def _extract_ncbi_protein_uid_from_gbseq_element(
    *,
    gbseq_element: ET.Element,
) -> str:
    """
    Extract one distinct UID from gi|<number> candidates in GBSeqid fields.

    Inspect direct and nested GBSeqid elements only. Match each complete value
    after trimming surrounding whitespace; retain the legacy optional trailing
    pipe. Bare numbers, other identifier domains and malformed values are not
    UID candidates. Missing candidates or distinct conflicting UIDs are errors.
    Repeating the same UID within a record still identifies that single record.
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


def _extract_ncbi_protein_accession_version_from_gbseq_element(
    *,
    gbseq_element: ET.Element,
) -> str:
    """
    Extract one distinct accession.version from the record's supported fields.

    A present GBSeq_accession-version must be valid and agree with any valid
    accession-bearing GBSeqid, direct or inside GBSeq_other-seqids. GI and local
    IDs are not accession candidates. Primary-accession and sequence-version
    never establish an exact revision.
    """
    primary_candidates = {
        (element.text or "").strip()
        for element in gbseq_element.findall("GBSeq_accession-version")
    }
    if primary_candidates:
        try:
            validate_ncbi_protein_identifiers(
                protein_identifiers=sorted(primary_candidates),
                identifier_type=NCBIProteinIdentifierType.accession_version,
            )
        except ValueError as error:
            raise RuntimeError(
                "Invalid accession.version in GBSeq_accession-version in one "
                f"XML record: {sorted(primary_candidates)}."
            ) from error

    seqid_candidates: set[str] = set()
    seqid_elements = (
        gbseq_element.findall("GBSeqid")
        + gbseq_element.findall("GBSeq_other-seqids/GBSeqid")
    )
    for element in seqid_elements:
        seqid_text = (element.text or "").strip()
        # Match complete accession-bearing Seq-ids, excluding gi/lcl/gnl.
        accession_match = re.fullmatch(
            r"(?:ref|gb|emb|dbj|tpg|tpe|tpd|sp|tr|pir|prf)\|([^|]+)\|?",
            seqid_text,
        )
        if accession_match is None:
            continue

        candidate = accession_match.group(1)
        try:
            validate_ncbi_protein_identifiers(
                protein_identifiers=[candidate],
                identifier_type=NCBIProteinIdentifierType.accession_version,
            )
        except ValueError:
            # An unversioned or malformed Seq-id cannot confirm a revision.
            continue
        seqid_candidates.add(candidate)

    accession_candidates = primary_candidates | seqid_candidates
    if not accession_candidates:
        raise RuntimeError(
            "Failed to extract an accession.version from GBSeqid fallback "
            "fields in one XML record."
        )

    if len(accession_candidates) != 1:
        raise RuntimeError(
            "Found multiple conflicting accession.version candidates in one "
            f"XML record. GBSeq_accession-version: {sorted(primary_candidates)}. "
            f"GBSeqid: {sorted(seqid_candidates)}."
        )

    return next(iter(accession_candidates))


def _validate_xml_record_identifiers(
    *,
    xml_file_path: PathLike,
    expected_protein_identifiers: list[str],
    identifier_type: NCBIProteinIdentifierType | str,
    expected_root_tag: str = "GBSet",
    expected_record_tag: str = "GBSeq",
) -> list[str]:
    """
    Validate XML structure and compare record identifiers with a requested multiset.

    The identifier domain is explicit; it is never inferred from XML contents.
    UID extraction retains the legacy behavior. Accession.version extraction
    requires supported versioned fields in each record to agree.
    Return identifiers in XML record order, preserving duplicate occurrences.
    """
    resolved_identifier_type = coerce_ncbi_protein_identifier_type(identifier_type)
    if resolved_identifier_type == NCBIProteinIdentifierType.uid:
        extract_identifier = _extract_ncbi_protein_uid_from_gbseq_element
        identifier_label = "UID"
    elif resolved_identifier_type == NCBIProteinIdentifierType.accession_version:
        extract_identifier = _extract_ncbi_protein_accession_version_from_gbseq_element
        identifier_label = "accession.version identifier"
    else:
        raise AssertionError(f"Unhandled identifier_type: {resolved_identifier_type}")

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
            "Saved consolidated XML snapshot root tag mismatch during "
            f"{identifier_label} "
            f"validation. Expected {expected_root_tag}, got {root_element.tag}."
        )

    child_elements = list(root_element)
    invalid_child_tags = sorted(
        {child.tag for child in child_elements if child.tag != expected_record_tag}
    )
    if invalid_child_tags:
        raise RuntimeError(
            "Saved consolidated XML snapshot contains unexpected child tags "
            f"during {identifier_label} validation: {invalid_child_tags}."
        )

    extracted_protein_identifiers = [
        extract_identifier(gbseq_element=child_element)
        for child_element in child_elements
    ]

    # Match record multiplicities, independently of order. Repeated candidates
    # within one record are resolved by the extractor; repeated records count.
    expected_identifier_counter = Counter(expected_protein_identifiers)
    extracted_identifier_counter = Counter(extracted_protein_identifiers)
    if extracted_identifier_counter != expected_identifier_counter:
        missing_identifiers = sorted(
            (
                expected_identifier_counter - extracted_identifier_counter
            ).elements()
        )
        unexpected_identifiers = sorted(
            (
                extracted_identifier_counter - expected_identifier_counter
            ).elements()
        )
        raise RuntimeError(
            f"Saved consolidated XML snapshot record {identifier_label}s do not match the "
            f"expected protein {identifier_label}s. "
            f"Missing: {missing_identifiers[:5]}; "
            f"Unexpected: {unexpected_identifiers[:5]}."
        )

    return extracted_protein_identifiers


def _validate_xml_record_uids(
    *,
    xml_file_path: PathLike,
    expected_protein_uids: list[str],
    expected_root_tag: str = "GBSet",
    expected_record_tag: str = "GBSeq",
) -> list[str]:
    """
    Legacy UID wrapper preserving the original parameters and return value.
    """
    return _validate_xml_record_identifiers(
        xml_file_path=xml_file_path,
        expected_protein_identifiers=expected_protein_uids,
        identifier_type="uid",
        expected_root_tag=expected_root_tag,
        expected_record_tag=expected_record_tag,
    )


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


def _build_identifier_snapshot_directory_name(
    *,
    retrieved_at_utc: str,
    source_label: str,
) -> str:
    """
    Build a deterministic directory name for a generic identifier snapshot.
    """
    return build_snapshot_directory_name(
        retrieved_at_utc=retrieved_at_utc,
        search_query=source_label,
    )


def _build_identifier_snapshot_source_label(
    *,
    identifier_type: str,
    source_type: str,
    protein_identifiers_sha256: str,
    search_query: Optional[str] = None,
    source_description: Optional[str] = None,
) -> str:
    """
    Build the stable label used to name generic identifier snapshot directories.
    """
    if search_query and search_query.strip():
        return search_query

    source_label_parts = [
        source_type,
        identifier_type,
        protein_identifiers_sha256,
    ]
    if source_description and source_description.strip():
        source_label_parts.append(source_description.strip())

    return (
        ":".join(source_label_parts)
    )


def _source_label_from_identifier_snapshot_manifest(
    manifest_payload: Dict[str, Any],
) -> str:
    """
    Resolve the label used when a downstream snapshot needs source provenance.
    """
    search_query = manifest_payload.get("search_query")
    if isinstance(search_query, str) and search_query.strip():
        return search_query

    identifier_type = _identifier_type_from_snapshot_manifest(manifest_payload).value
    source_type = manifest_payload.get("source_type")
    if not isinstance(source_type, str) or not source_type.strip():
        source_type = "unknown_identifier_source"

    protein_identifiers_sha256 = _manifest_field_value(
        manifest_payload, "protein_identifiers_sha256", "protein_uids_sha256",
    )
    if not isinstance(protein_identifiers_sha256, str) or not (
        protein_identifiers_sha256.strip()
    ):
        protein_identifiers_sha256 = "missing_identifier_hash"

    source_description = manifest_payload.get("source_description")
    if not isinstance(source_description, str):
        source_description = None

    return _build_identifier_snapshot_source_label(
        identifier_type=identifier_type,
        source_type=source_type,
        protein_identifiers_sha256=protein_identifiers_sha256,
        source_description=source_description,
    )


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


def _write_snapshot_identifier_files(
    *,
    snapshot_directory: Path,
    protein_identifiers: list[str],
) -> list[tuple[Path, str]]:
    """Write the canonical 1.1 file plus an identical legacy compatibility copy."""
    canonical_path = snapshot_directory / NCBI_PROTEIN_IDENTIFIERS_FILE_NAME
    legacy_path = snapshot_directory / NCBI_LEGACY_PROTEIN_UIDS_FILE_NAME
    _write_text_lines(text_lines=protein_identifiers, output_file_path=canonical_path)
    shutil.copy2(canonical_path, legacy_path)
    return [
        (canonical_path, NCBI_PROTEIN_IDENTIFIERS_FILE_NAME),
        (legacy_path, NCBI_LEGACY_PROTEIN_UIDS_FILE_NAME),
    ]


def _build_snapshot_manifest(
    *,
    fetch_result: NCBIProteinUidFetchResult,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
) -> Dict[str, Any]:
    """
    Build the 1.1 identifier manifest with generic fields and legacy aliases.
    """
    manifest_payload = asdict(fetch_result)
    manifest_payload.pop("protein_uids", None)

    manifest_payload["snapshot_format_version"] = NCBI_SNAPSHOT_FORMAT_VERSION
    manifest_payload["identifier_type"] = _identifier_type_for_new_snapshot(
        fetch_result.identifier_type,
    ).value
    manifest_payload["artifact_type"] = "ncbi_protein_identifier_snapshot"
    manifest_payload["source_type"] = "entrez_search"
    manifest_payload["snapshot_file_name"] = "protein_uids.txt"
    manifest_payload["protein_identifiers_file_name"] = NCBI_PROTEIN_IDENTIFIERS_FILE_NAME
    manifest_payload["protein_uids_file_name"] = "protein_uids.txt"
    manifest_payload["raw_protein_identifier_count"] = (
        fetch_result.raw_protein_identifier_count
    )
    manifest_payload["normalized_protein_identifier_count"] = (
        fetch_result.normalized_protein_identifier_count
    )
    manifest_payload["protein_identifier_count"] = (
        fetch_result.protein_identifier_count
    )
    manifest_payload["protein_identifiers_sha256"] = (
        fetch_result.protein_identifiers_sha256
    )
    manifest_payload["deduplicate_identifiers"] = fetch_result.deduplicate_uids
    manifest_payload["sort_identifiers"] = fetch_result.sort_uids
    manifest_payload["manifest_file_name"] = "manifest.json"
    manifest_payload["immutable_snapshot_directory_name"] = (
        immutable_snapshot_directory_name
    )
    manifest_payload["immutable_snapshot_relative_path"] = (
        immutable_snapshot_relative_path
    )

    return manifest_payload


def _build_provided_identifier_snapshot_manifest(
    *,
    protein_identifiers: list[str],
    identifier_type: str,
    source_type: str,
    retrieved_at_utc: str,
    raw_protein_identifier_count: int,
    deduplicate_identifiers: bool,
    sort_identifiers: bool,
    protein_identifiers_sha256: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    source_description: Optional[str] = None,
    source_name: Optional[str] = None,
    source_file_path: Optional[PathLike] = None,
) -> Dict[str, Any]:
    """
    Build the manifest for a generic identifier-list snapshot.
    """
    manifest_payload: Dict[str, Any] = {
        "snapshot_format_version": NCBI_SNAPSHOT_FORMAT_VERSION,
        "artifact_type": "ncbi_protein_identifier_snapshot",
        "source_type": source_type,
        "identifier_type": _identifier_type_for_new_snapshot(identifier_type).value,
        "retrieved_at_utc": retrieved_at_utc,
        "raw_protein_identifier_count": raw_protein_identifier_count,
        "normalized_protein_identifier_count": len(protein_identifiers),
        "protein_identifier_count": len(protein_identifiers),
        "deduplicate_identifiers": deduplicate_identifiers,
        "sort_identifiers": sort_identifiers,
        "protein_identifiers_sha256": protein_identifiers_sha256,
        "protein_identifiers_file_name": NCBI_PROTEIN_IDENTIFIERS_FILE_NAME,
        "manifest_file_name": "manifest.json",
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        # Compatibility fields consumed by existing UID snapshot readers.
        "snapshot_file_name": "protein_uids.txt",
        "protein_uids_file_name": "protein_uids.txt",
        "normalized_protein_uid_count": len(protein_identifiers),
        "protein_uids_sha256": protein_identifiers_sha256,
    }

    if source_description is not None:
        manifest_payload["source_description"] = source_description
    elif source_name is not None:
        manifest_payload["source_description"] = source_name

    if source_name is not None:
        manifest_payload["source_name"] = source_name

    if source_file_path is not None:
        manifest_payload["source_file_name"] = _as_path(source_file_path).name

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
    source_identifier_count = _manifest_field_value(
        source_uid_snapshot_manifest,
        "normalized_protein_identifier_count", "protein_identifier_count", "normalized_protein_uid_count",
    )
    source_identifier_sha256 = _manifest_field_value(
        source_uid_snapshot_manifest, "protein_identifiers_sha256", "protein_uids_sha256",
    )
    return {
        "snapshot_format_version": NCBI_SNAPSHOT_FORMAT_VERSION,
        "artifact_type": "ncbi_protein_xml_snapshot",
        "source_type": "identifier_snapshot",
        "database_name": fetch_result.database_name,
        "identifier_type": _identifier_type_for_new_snapshot(fetch_result.identifier_type).value,
        "search_query": source_uid_snapshot_manifest.get("search_query"),
        "translated_query": source_uid_snapshot_manifest.get("translated_query"),
        "retrieved_at_utc": fetch_result.retrieved_at_utc,
        "retmode": fetch_result.retmode,
        "requested_protein_identifier_count": fetch_result.requested_protein_identifier_count,
        "normalized_protein_identifier_count": fetch_result.normalized_protein_identifier_count,
        "protein_identifier_count": fetch_result.protein_identifier_count,
        "protein_identifiers_sha256": fetch_result.protein_identifiers_sha256,
        # Keep legacy aliases for existing consumers of XML snapshots.
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
        "protein_identifiers_file_name": NCBI_PROTEIN_IDENTIFIERS_FILE_NAME,
        "protein_uids_file_name": "protein_uids.txt",
        "xml_file_name": consolidated_xml_file_name,
        "xml_file_sha256": consolidated_xml_file_sha256,
        "consolidated_record_count": consolidated_record_count,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "source_identifier_snapshot_relative_path": source_uid_snapshot_manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_identifier_snapshot_directory_name": source_uid_snapshot_manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_identifier_snapshot_manifest_sha256": source_uid_snapshot_manifest_sha256,
        "source_uid_snapshot_relative_path": source_uid_snapshot_manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_uid_snapshot_directory_name": source_uid_snapshot_manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_uid_snapshot_manifest_sha256": source_uid_snapshot_manifest_sha256,
        "source_identifier_snapshot_artifact_type": source_uid_snapshot_manifest.get(
            "artifact_type"
        ),
        "source_identifier_snapshot_source_type": source_uid_snapshot_manifest.get(
            "source_type"
        ),
        "source_identifier_count": source_identifier_count,
        "source_identifier_sha256": source_identifier_sha256,
        "source_uid_count": source_identifier_count,
        "source_uid_sha256": source_identifier_sha256,
        "batches": batch_metadata_records,
    }


def save_ncbi_protein_identifier_snapshot(
    *,
    protein_identifiers: list[Any],
    identifier_type: str,
    snapshot_root_directory: PathLike,
    source_type: str = "provided_identifier_list",
    source_description: Optional[str] = None,
    source_name: Optional[str] = None,
    source_file_path: Optional[PathLike] = None,
    deduplicate_identifiers: bool = True,
    sort_identifiers: bool = True,
    retrieved_at_utc: Optional[str] = None,
    update_latest_directory: bool = True,
) -> Path:
    """
    Persist an immutable generic NCBI protein identifier snapshot.

    This entry point is for identifier lists supplied directly by a caller,
    workflow, or file. It does not call NCBI.
    """
    if not source_type or not source_type.strip():
        raise ValueError("source_type must be a non-empty string.")

    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_identifier_type = _identifier_type_for_new_snapshot(
        identifier_type,
    )
    normalized_protein_identifiers = _normalize_protein_identifier_list(
        protein_identifiers=protein_identifiers,
        identifier_type=resolved_identifier_type,
        deduplicate_identifiers=deduplicate_identifiers,
        sort_identifiers=sort_identifiers,
    )

    if not normalized_protein_identifiers:
        raise ValueError(
            "protein_identifiers must contain at least one non-empty "
            "identifier."
        )

    resolved_retrieved_at_utc = (
        retrieved_at_utc
        if retrieved_at_utc is not None
        else datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    )
    protein_identifiers_sha256 = sha256_of_lines(
        text_lines=normalized_protein_identifiers,
        deduplicate_lines_preserving_order=False,
        sort_lines=False,
    )
    source_label = _build_identifier_snapshot_source_label(
        identifier_type=resolved_identifier_type.value,
        source_type=source_type,
        protein_identifiers_sha256=protein_identifiers_sha256,
        source_description=(
            source_description if source_description is not None else source_name
        ),
    )
    snapshot_directory_name = _build_identifier_snapshot_directory_name(
        retrieved_at_utc=resolved_retrieved_at_utc,
        source_label=source_label,
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

        manifest_file_path = immutable_snapshot_directory / "manifest.json"

        identifier_files = _write_snapshot_identifier_files(
            snapshot_directory=immutable_snapshot_directory,
            protein_identifiers=normalized_protein_identifiers,
        )

        manifest_payload = _build_provided_identifier_snapshot_manifest(
            protein_identifiers=normalized_protein_identifiers,
            identifier_type=resolved_identifier_type.value,
            source_type=source_type,
            retrieved_at_utc=resolved_retrieved_at_utc,
            raw_protein_identifier_count=len(protein_identifiers),
            deduplicate_identifiers=deduplicate_identifiers,
            sort_identifiers=sort_identifiers,
            protein_identifiers_sha256=protein_identifiers_sha256,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
            source_description=source_description,
            source_name=source_name,
            source_file_path=source_file_path,
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
                    *identifier_files,
                    (manifest_file_path, "manifest.json"),
                ],
            )
    except Exception:
        if not immutable_snapshot_complete and immutable_snapshot_directory.exists():
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

    return immutable_snapshot_directory


def save_ncbi_protein_identifier_snapshot_from_file(
    *,
    identifier_file_path: PathLike,
    identifier_type: str,
    snapshot_root_directory: PathLike,
    source_type: str = "provided_identifier_list",
    source_description: Optional[str] = None,
    source_name: Optional[str] = None,
    deduplicate_identifiers: bool = True,
    sort_identifiers: bool = True,
    retrieved_at_utc: Optional[str] = None,
    update_latest_directory: bool = True,
) -> Path:
    """
    Persist a generic identifier snapshot from a local text file.
    """
    resolved_identifier_file_path = _as_path(identifier_file_path)
    protein_identifiers = read_text_lines_from_file(
        input_file_path=resolved_identifier_file_path,
    )

    return save_ncbi_protein_identifier_snapshot(
        protein_identifiers=protein_identifiers,
        identifier_type=identifier_type,
        snapshot_root_directory=snapshot_root_directory,
        source_type=source_type,
        source_description=source_description,
        source_name=source_name,
        source_file_path=resolved_identifier_file_path,
        deduplicate_identifiers=deduplicate_identifiers,
        sort_identifiers=sort_identifiers,
        retrieved_at_utc=retrieved_at_utc,
        update_latest_directory=update_latest_directory,
    )


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
    resolved_identifier_type = _identifier_type_for_new_snapshot(
        fetch_result.identifier_type,
    )
    validate_ncbi_protein_identifiers(
        protein_identifiers=fetch_result.protein_uids,
        identifier_type=resolved_identifier_type,
        parameter_name="fetch_result.protein_uids",
    )

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

        manifest_file_path = immutable_snapshot_directory / "manifest.json"

        identifier_files = _write_snapshot_identifier_files(
            snapshot_directory=immutable_snapshot_directory,
            protein_identifiers=fetch_result.protein_identifiers,
        )

        manifest_payload = _build_snapshot_manifest(
            fetch_result=fetch_result,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
        )
        for identifier_file_path, _ in identifier_files:
            _validate_saved_identifier_list(
                manifest_payload=manifest_payload,
                protein_identifiers=load_snapshot_protein_uids(
                    protein_uids_file_path=identifier_file_path,
                ),
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
                    *identifier_files,
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
    source_uid_snapshot_manifest_sha256: Optional[str] = None,
    source_uid_snapshot_manifest_bytes: Optional[bytes] = None,
) -> Path:
    """
    Persist one immutable local snapshot of consolidated NCBI protein XML.

    XML is fetched in batches upstream for robustness, then merged into one
    valid XML document with a single root and persisted as one file.
    The fetch must explicitly declare its identifier domain and match the source
    manifest. Only an unversioned or 1.0 source manifest without identifier_type
    defaults to UID. XML record identity is validated in that same domain.
    A valid source identifier SHA-256 is required. The generic hash takes
    precedence; its legacy alias is read only when the generic field is absent.
    Legacy parameter and file names are retained for both supported domains.
    The resolver supplies the manifest bytes captured when selecting the source.
    Direct callers can omit them to read and verify the supplied file.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    if not protein_uids:
        raise ValueError(
            "protein_uids must contain at least one protein identifier."
        )

    resolved_source_identifier_type = _identifier_type_from_snapshot_manifest(
        source_uid_snapshot_manifest,
    )
    resolved_fetch_identifier_type = _validate_xml_fetch_identifier_type(
        fetch_result=fetch_result,
        expected_identifier_type=resolved_source_identifier_type,
    )
    validate_ncbi_protein_identifiers(
        protein_identifiers=protein_uids,
        identifier_type=resolved_fetch_identifier_type,
        parameter_name="protein_uids",
    )

    if (
        type(fetch_result.normalized_protein_identifier_count) is not int
        or fetch_result.normalized_protein_identifier_count != len(protein_uids)
    ):
        raise ValueError(
            "The XML fetch result identifier count does not match the provided "
            "protein_uids list length."
        )

    source_identifier_sha256 = _validate_source_identifier_list(
        protein_identifiers=protein_uids,
        source_manifest=source_uid_snapshot_manifest,
    )
    if source_identifier_sha256 != fetch_result.protein_identifiers_sha256:
        raise ValueError(
            "The XML fetch result does not match the source identifier snapshot. "
            "protein_identifiers_sha256 values differ."
        )

    source_manifest_bytes = source_uid_snapshot_manifest_bytes
    if source_manifest_bytes is None:
        source_manifest_bytes = _as_path(source_uid_snapshot_manifest_file_path).read_bytes()
    if json.dumps(json.loads(source_manifest_bytes), sort_keys=True) != json.dumps(
        source_uid_snapshot_manifest, sort_keys=True,
    ):
        raise ValueError("The source manifest file differs from the selected source manifest.")
    source_manifest_sha256 = hashlib.sha256(source_manifest_bytes).hexdigest()
    if (
        source_uid_snapshot_manifest_sha256 is not None
        and source_uid_snapshot_manifest_sha256 != source_manifest_sha256
    ):
        raise ValueError("The source manifest SHA-256 does not match the source manifest bytes.")

    snapshot_directory_name = _build_identifier_snapshot_directory_name(
        retrieved_at_utc=fetch_result.retrieved_at_utc,
        source_label=_source_label_from_identifier_snapshot_manifest(
            source_uid_snapshot_manifest,
        ),
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

        manifest_file_path = immutable_snapshot_directory / "manifest.json"
        consolidated_xml_file_name = "protein_records.xml"
        consolidated_xml_file_path = (
            immutable_snapshot_directory / consolidated_xml_file_name
        )

        identifier_files = _write_snapshot_identifier_files(
            snapshot_directory=immutable_snapshot_directory,
            protein_identifiers=protein_uids,
        )
        for identifier_file_path, _ in identifier_files:
            saved_identifiers = load_snapshot_protein_uids(
                protein_uids_file_path=identifier_file_path,
            )
            saved_identifier_sha256 = _validate_source_identifier_list(
                protein_identifiers=saved_identifiers,
                source_manifest=source_uid_snapshot_manifest,
            )
            if saved_identifier_sha256 != fetch_result.protein_identifiers_sha256:
                raise ValueError(
                    "The saved protein identifier list does not match the XML "
                    "fetch result. protein_identifiers_sha256 values differ."
                )

        batch_metadata_records: list[Dict[str, Any]] = []
        for xml_batch in fetch_result.xml_batches:
            batch_metadata_records.append(
                {
                    "batch_index": xml_batch.batch_index,
                    "batch_start_index": xml_batch.batch_start_index,
                    "batch_end_index": xml_batch.batch_end_index,
                    "protein_identifier_count": xml_batch.protein_identifier_count,
                    "protein_uid_count": xml_batch.protein_uid_count,
                    "xml_payload_sha256": xml_batch.xml_payload_sha256,
                }
            )

        consolidated_xml_payload, consolidated_record_count = (
            _build_consolidated_xml_payload(fetch_result=fetch_result)
        )
        expected_record_count = fetch_result.normalized_protein_identifier_count

        if consolidated_record_count != expected_record_count:
            raise RuntimeError(
                "Consolidated XML snapshot record count does not match the "
                "expected protein UID count. "
                f"Expected {expected_record_count}, got "
                f"{consolidated_record_count}."
            )

        write_bytes_atomic(
            binary_payload=consolidated_xml_payload,
            output_file_path=consolidated_xml_file_path,
        )

        validated_record_count = _validate_saved_consolidated_xml_snapshot(
            xml_file_path=consolidated_xml_file_path,
            expected_record_count=expected_record_count,
        )
        _validate_xml_record_identifiers(
            xml_file_path=consolidated_xml_file_path,
            expected_protein_identifiers=protein_uids,
            identifier_type=resolved_fetch_identifier_type,
        )

        consolidated_xml_file_sha256 = sha256_of_file(
            input_file_path=consolidated_xml_file_path,
        )

        manifest_payload = _build_xml_snapshot_manifest(
            fetch_result=fetch_result,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
            source_uid_snapshot_manifest=source_uid_snapshot_manifest,
            source_uid_snapshot_manifest_sha256=source_manifest_sha256,
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
                    *identifier_files,
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
    Load identifiers and validate the versioned manifest's integrity metadata.
    """
    resolved_snapshot_directory = _as_path(snapshot_directory)

    manifest_file_path = resolved_snapshot_directory / "manifest.json"
    manifest_bytes = manifest_file_path.read_bytes()
    manifest_payload = json.loads(manifest_bytes)
    protein_uids_file_path = _protein_identifier_file_path_from_manifest(
        resolved_snapshot_directory, manifest_payload,
    )
    protein_uids = load_snapshot_protein_uids(
        protein_uids_file_path=protein_uids_file_path,
    )
    resolved_identifier_type = _identifier_type_from_snapshot_manifest(manifest_payload)
    validate_ncbi_protein_identifiers(
        protein_identifiers=protein_uids,
        identifier_type=resolved_identifier_type,
        parameter_name="protein_uids",
    )
    _validate_saved_identifier_list(
        manifest_payload=manifest_payload,
        protein_identifiers=protein_uids,
    )

    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest_file_bytes": manifest_bytes,
        "manifest_file_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
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

    Use the manifest's identifier domain, defaulting to UID only for legacy
    manifests without that field. Format 1.1 requires canonical identifier and
    XML hashes and counts. Older manifests retain their legacy field fallbacks.
    """
    resolved_snapshot_directory = _as_path(snapshot_directory)

    _validate_snapshot_integrity_fields(manifest_payload, include_xml=True)
    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != "ncbi_protein_xml_snapshot":
        raise RuntimeError(
            "Saved XML snapshot manifest artifact_type mismatch. "
            f"Expected 'ncbi_protein_xml_snapshot', got {artifact_type!r}."
        )

    resolved_identifier_type = _identifier_type_from_snapshot_manifest(manifest_payload)
    validate_ncbi_protein_identifiers(
        protein_identifiers=protein_uids,
        identifier_type=resolved_identifier_type,
        parameter_name="protein_uids",
    )

    resolved_xml_file_path = _xml_file_path_from_manifest(
        resolved_snapshot_directory, manifest_payload,
    )
    if not resolved_xml_file_path.exists():
        raise FileNotFoundError(
            f"Saved XML snapshot file not found: {resolved_xml_file_path}."
        )

    expected_identifier_count = _validate_saved_identifier_list(
        manifest_payload=manifest_payload,
        protein_identifiers=protein_uids,
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
    _validate_xml_record_identifiers(
        xml_file_path=resolved_xml_file_path,
        expected_protein_identifiers=protein_uids,
        identifier_type=resolved_identifier_type,
    )

    if (
        expected_identifier_count is not None
        and validated_record_count != expected_identifier_count
    ):
        raise RuntimeError(
            "Saved XML snapshot record count does not match the saved protein "
            "identifier count. "
            f"Expected {expected_identifier_count}, got {validated_record_count}."
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
    return _with_validated_xml_file(
        load_snapshot_by_directory(snapshot_directory=snapshot_directory),
    )


def _with_validated_xml_file(snapshot_payload: Dict[str, Any]) -> Dict[str, Any]:
    snapshot_payload["xml_file_path"] = _validate_loaded_xml_snapshot_payload(
        snapshot_directory=snapshot_payload["snapshot_directory"],
        manifest_payload=snapshot_payload["manifest"],
        protein_uids=snapshot_payload["protein_uids"],
    )
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
    Return the preferred identifier file path, retaining the legacy API name.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    return get_snapshot_protein_uids_path(
        snapshot_directory=resolved_snapshot_root_directory / "latest",
    )


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
    Return the manifest-selected identifier file, or the unversioned legacy path.
    """
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / "manifest.json"
    if not manifest_file_path.exists():
        return resolved_snapshot_directory / NCBI_LEGACY_PROTEIN_UIDS_FILE_NAME
    manifest_payload = load_snapshot_manifest(manifest_file_path=manifest_file_path)
    return _protein_identifier_file_path_from_manifest(resolved_snapshot_directory, manifest_payload)


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


def _load_available_snapshot_manifest(snapshot_directory: Path) -> Optional[Dict[str, Any]]:
    """Read the manifest of a possible reuse candidate."""
    manifest_file_path = snapshot_directory / "manifest.json"
    if not manifest_file_path.is_file():
        return None
    try:
        manifest_payload = load_snapshot_manifest(manifest_file_path=manifest_file_path)
    except (OSError, ValueError, RuntimeError, TypeError):
        # Incomplete/unreadable copies cannot be offered as reusable snapshots.
        return None
    return manifest_payload if isinstance(manifest_payload, dict) else None


def _snapshot_files_are_available(
    snapshot_directory: Path,
    *,
    include_xml: bool = False,
    manifest_payload: Optional[Dict[str, Any]] = None,
) -> bool:
    """Probe files selected by a manifest; full integrity is checked on load."""
    if manifest_payload is None:
        manifest_payload = _load_available_snapshot_manifest(snapshot_directory)
    if manifest_payload is None:
        return False
    identifier_file_path = _protein_identifier_file_path_from_manifest(snapshot_directory, manifest_payload)
    if not identifier_file_path.is_file():
        return False
    if include_xml:
        return _xml_file_path_from_manifest(snapshot_directory, manifest_payload).is_file()
    return True


def latest_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
) -> bool:
    """
    Return True only if the convenience latest snapshot copy is complete.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    return _snapshot_files_are_available(
        resolved_snapshot_root_directory / "latest",
    )


def latest_xml_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
) -> bool:
    """
    Return True only if the convenience latest XML snapshot copy is complete.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    return _snapshot_files_are_available(
        resolved_snapshot_root_directory / "latest", include_xml=True,
    )


def resolve_ncbi_protein_uid_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    search_query: Optional[str] = None,
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

    latest_is_available = (
        resolved_snapshot_mode != SnapshotMode.create_new
        and latest_snapshot_is_available(snapshot_root_directory=resolved_snapshot_root_directory)
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

            latest_protein_uids_file_path = get_snapshot_protein_uids_path(
                snapshot_directory=latest_directory,
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

    if not search_query or not search_query.strip():
        raise ValueError(
            "search_query is required when snapshot creation from NCBI is needed."
        )

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


def _load_xml_source_snapshot(source_root_directory: Path) -> Dict[str, Any]:
    """Select the immutable identifier snapshot referenced by the current latest copy."""
    latest = resolve_ncbi_protein_uid_snapshot(
        snapshot_mode=SnapshotMode.reuse_latest,
        snapshot_root_directory=source_root_directory,
    )
    manifest = latest["manifest"]
    relative_path = manifest.get("immutable_snapshot_relative_path")
    directory_name = manifest.get("immutable_snapshot_directory_name")
    if relative_path is None and directory_name is None:
        if _snapshot_format_version(manifest) == "1.0":
            return latest
        raise RuntimeError("The source manifest must identify its immutable snapshot directory.")

    if relative_path is None:
        relative_path = f"snapshots/{directory_name}"
    if not isinstance(relative_path, str):
        raise RuntimeError("The source manifest must identify a directory under snapshots/.")
    parts = PureWindowsPath(relative_path).parts
    if (
        len(parts) != 2 or parts[0] != "snapshots"
        or (directory_name is not None and parts[1] != directory_name)
    ):
        raise RuntimeError("The source manifest must identify a directory under snapshots/.")
    snapshot_directory = _snapshot_member_path(source_root_directory, relative_path)
    if snapshot_directory.resolve().parent != source_root_directory.resolve() / "snapshots":
        raise RuntimeError("The source manifest must identify an immutable directory under snapshots/.")

    source = load_snapshot_by_directory(snapshot_directory=snapshot_directory)
    if source["manifest_file_sha256"] != latest["manifest_file_sha256"]:
        raise RuntimeError("The latest source manifest differs from its immutable snapshot manifest.")
    return source


def resolve_ncbi_protein_xml_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    source_uid_snapshot_root_directory: Optional[PathLike] = None,
    source_identifier_snapshot_root_directory: Optional[PathLike] = None,
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
    xml_batch_size: int = 100,
    max_retry_attempts: int = 5,
    xml_request_delay_seconds: Optional[float] = None,
    fetch_timeout_seconds: float = 30.0,
    batch_deadline_seconds: float = 300.0,
    retry_backoff_initial_seconds: Optional[float] = None,
    retry_backoff_multiplier: float = 2.0,
    retry_backoff_max_seconds: float = 30.0,
    circuit_breaker_failure_threshold: int = 3,
    circuit_breaker_cooldown_seconds: float = 60.0,
    ncbi_email: Optional[str] = None,
    ncbi_api_key: Optional[str] = None,
    update_latest_directory: bool = True,
) -> Dict[str, Any]:
    """
    Resolve the active raw XML snapshot payload according to the requested mode.

    If snapshot creation is needed, the XML retrieval is performed strictly from an
    already-frozen upstream identifier snapshot, using its declared identifier_type.
    This function does not create the source snapshot.

    source_identifier_snapshot_root_directory is the generic source-root name;
    source_uid_snapshot_root_directory is retained as a compatible alias. Creation
    requires either argument. When both are supplied, they must resolve to the
    same directory. In reuse_latest mode the source arguments are ignored. In
    reuse_latest_or_create mode a supplied source must match the latest XML's
    provenance and identifiers before that XML is opened and validated for
    reuse. Reuse needs no NCBI credentials.

    Circuit-breaker state is retained by the NCBI API infrastructure across
    repeated snapshot workflow invocations for the same threshold and cooldown.
    Callers configure the policy here but do not need to own its mutable state.
    """
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"

    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_xml_snapshot_is_available(
            snapshot_root_directory=resolved_snapshot_root_directory,
        ):
            latest_manifest_file_path = latest_directory / "manifest.json"
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

            latest_protein_uids_file_path = get_snapshot_protein_uids_path(
                snapshot_directory=latest_directory,
            )
            latest_manifest = load_snapshot_manifest(manifest_file_path=latest_manifest_file_path)
            if not latest_protein_uids_file_path.exists():
                raise FileNotFoundError(
                    f"Latest XML snapshot protein UID file not found: "
                    f"{latest_protein_uids_file_path}. Run the workflow once with "
                    f"snapshot_mode='create_new' to create it."
                )

            latest_xml_file_path = _xml_file_path_from_manifest(latest_directory, latest_manifest)
            raise FileNotFoundError(
                f"Latest XML snapshot consolidated XML file not found: "
                f"{latest_xml_file_path}. Run the workflow once with "
                f"snapshot_mode='create_new' to create it."
            )

        return load_latest_xml_snapshot(
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

    if (
        source_uid_snapshot_root_directory is not None
        and source_identifier_snapshot_root_directory is not None
        and _as_path(source_uid_snapshot_root_directory).resolve()
        != _as_path(source_identifier_snapshot_root_directory).resolve()
    ):
        raise ValueError(
            "source_identifier_snapshot_root_directory and "
            "source_uid_snapshot_root_directory must refer to the same source directory."
        )
    source_root_directory = (
        source_identifier_snapshot_root_directory
        if source_identifier_snapshot_root_directory is not None
        else source_uid_snapshot_root_directory
    )
    if source_root_directory is None:
        if (
            resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
            and latest_xml_snapshot_is_available(snapshot_root_directory=resolved_snapshot_root_directory)
        ):
            latest_xml_snapshot = load_latest_xml_snapshot(
                snapshot_root_directory=resolved_snapshot_root_directory,
            )
            print("Latest XML snapshot is available. Reusing frozen snapshot.")
            return latest_xml_snapshot
        raise ValueError(
            "source_identifier_snapshot_root_directory (or legacy "
            "source_uid_snapshot_root_directory) is required when creating an XML snapshot."
        )

    try:
        source_identifier_snapshot_payload = _load_xml_source_snapshot(_as_path(source_root_directory))
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "No reusable source identifier snapshot was found for the XML workflow. "
            "Create the upstream protein identifier snapshot before running the XML "
            "snapshot step."
        ) from exc

    protein_uids = source_identifier_snapshot_payload["protein_uids"]
    source_identifier_snapshot_manifest = source_identifier_snapshot_payload["manifest"]
    source_identifier_snapshot_manifest_file_path = source_identifier_snapshot_payload[
        "manifest_file_path"
    ]
    source_identifier_type = _identifier_type_from_snapshot_manifest(source_identifier_snapshot_manifest)
    _validate_source_identifier_list(
        protein_identifiers=protein_uids,
        source_manifest=source_identifier_snapshot_manifest,
    )
    source_manifest_sha256 = source_identifier_snapshot_payload["manifest_file_sha256"]
    source_manifest_bytes = source_identifier_snapshot_payload["manifest_file_bytes"]
    latest_manifest = (
        _load_available_snapshot_manifest(latest_directory)
        if resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        else None
    )
    if latest_manifest is not None:
        latest_source_sha256 = _manifest_field_value(
            latest_manifest,
            "source_identifier_snapshot_manifest_sha256",
            "source_uid_snapshot_manifest_sha256",
        )
        source_hash_matches = latest_source_sha256 == source_manifest_sha256
        if not source_hash_matches and _snapshot_format_version(source_identifier_snapshot_manifest) == "1.0":
            # Git can translate line endings in committed legacy manifests.
            lf_bytes = source_manifest_bytes.replace(b"\r\n", b"\n")
            crlf_bytes = lf_bytes.replace(b"\n", b"\r\n")
            source_hash_matches = latest_source_sha256 in (
                hashlib.sha256(lf_bytes).hexdigest(),
                hashlib.sha256(crlf_bytes).hexdigest(),
            )
        if (
            source_hash_matches
            and _identifier_type_from_snapshot_manifest(latest_manifest) == source_identifier_type
            and _snapshot_files_are_available(
                latest_directory, include_xml=True, manifest_payload=latest_manifest,
            )
        ):
            latest_snapshot = load_snapshot_by_directory(snapshot_directory=latest_directory)
            # latest may have been replaced since the provenance check.
            if (
                latest_snapshot["manifest"] == latest_manifest
                and latest_snapshot["protein_uids"] == protein_uids
            ):
                latest_xml_snapshot = _with_validated_xml_file(latest_snapshot)
                print("Latest XML snapshot matches the selected source. Reusing frozen snapshot.")
                return latest_xml_snapshot

    if not ncbi_email:
        raise ValueError(
            "ncbi_email is required when snapshot creation from NCBI is needed."
        )

    fetch_result = fetch_ncbi_protein_xml_batches(
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        protein_uids=protein_uids,
        identifier_type=source_identifier_type,
        ssl_ca_file=ssl_ca_file,
        ssl_ca_directory=ssl_ca_directory,
        batch_size=xml_batch_size,
        max_retry_attempts=max_retry_attempts,
        request_delay_seconds=xml_request_delay_seconds,
        fetch_timeout_seconds=fetch_timeout_seconds,
        batch_deadline_seconds=batch_deadline_seconds,
        retry_backoff_initial_seconds=retry_backoff_initial_seconds,
        retry_backoff_multiplier=retry_backoff_multiplier,
        retry_backoff_max_seconds=retry_backoff_max_seconds,
        circuit_breaker_failure_threshold=circuit_breaker_failure_threshold,
        circuit_breaker_cooldown_seconds=circuit_breaker_cooldown_seconds,
        circuit_breaker=get_default_ncbi_xml_circuit_breaker(
            failure_threshold=circuit_breaker_failure_threshold,
            cooldown_seconds=circuit_breaker_cooldown_seconds,
        ),
    )

    _validate_xml_fetch_identifier_type(
        fetch_result=fetch_result,
        expected_identifier_type=source_identifier_type,
    )

    saved_snapshot_directory = save_ncbi_protein_xml_snapshot(
        fetch_result=fetch_result,
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_uid_snapshot_manifest=source_identifier_snapshot_manifest,
        source_uid_snapshot_manifest_file_path=source_identifier_snapshot_manifest_file_path,
        source_uid_snapshot_manifest_sha256=source_manifest_sha256,
        source_uid_snapshot_manifest_bytes=source_manifest_bytes,
        protein_uids=protein_uids,
        update_latest_directory=update_latest_directory,
    )

    return load_xml_snapshot_by_directory(
        snapshot_directory=saved_snapshot_directory,
    )
