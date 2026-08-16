from __future__ import annotations

import hashlib
import shutil
import tempfile
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from enum import StrEnum
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, Optional, Union

from src.pago_pipeline.ncbi_api import (
    NCBI_PROTEIN_GBSEQ_XML_RETTYPE,
    NCBIProteinUidFetchResult,
    NCBIProteinXmlFetchResult,
    fetch_ncbi_protein_uid_snapshot,
    fetch_ncbi_protein_xml_batches,
    get_default_ncbi_xml_circuit_breaker,
)
from src.pago_pipeline.ncbi_batch_workspace import purge_batch_workspace_directory
from src.pago_pipeline.ncbi_telemetry import NCBIRetrievalTelemetry
from src.pago_pipeline.ncbi_xml_stream import (
    XmlBatchPayloadSource,
    extract_protein_uid_from_gbseq_element,
    validate_consolidated_xml_file,
    write_consolidated_xml_document,
)
from src.pago_pipeline.storage import (
    read_json_file,
    read_text_lines_from_file,
    save_ncbi_protein_uids_as_txt,
    sha256_of_file,
    sha256_of_lines,
    write_json_atomic,
)

PathLike = Union[str, Path]

UID_SNAPSHOT_FORMAT_VERSION = "1.1"
XML_SNAPSHOT_FORMAT_VERSION = "1.1"
SUPPORTED_XML_SNAPSHOT_FORMAT_VERSIONS = frozenset({"1.0", "1.1"})
BATCH_WORKSPACE_DIRECTORY_NAME = ".batch_workspace"


class SnapshotMode(StrEnum):
    create_new = "create_new"
    reuse_latest = "reuse_latest"
    reuse_latest_or_create = "reuse_latest_or_create"


@contextmanager
def _measure_local_phase(
    *,
    telemetry: Optional[NCBIRetrievalTelemetry],
    phase_name: str,
):
    """
    Measure one local processing phase when a telemetry collector is present.
    """
    if telemetry is None:
        yield
        return

    with telemetry.measure_local_phase(phase_name=phase_name):
        yield


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


def _build_xml_batch_payload_sources(
    *,
    fetch_result: NCBIProteinXmlFetchResult,
) -> list[XmlBatchPayloadSource]:
    """
    Order fetched batches by plan index and expose them as consolidation inputs.

    Ordering is taken from the batch index rather than from arrival order,
    because concurrent retrieval completes batches out of order while the
    consolidated snapshot hash is defined by plan order alone.
    """
    if not fetch_result.xml_batches:
        raise ValueError("fetch_result.xml_batches must contain at least one batch.")

    return [
        XmlBatchPayloadSource(
            batch_index=xml_batch.batch_index,
            payload_bytes=xml_batch.xml_payload_bytes,
            payload_file_path=xml_batch.xml_payload_file_path,
        )
        for xml_batch in sorted(
            fetch_result.xml_batches,
            key=lambda xml_batch: xml_batch.batch_index,
        )
    ]


def _extract_ncbi_protein_uid_from_gbseq_element(
    *,
    gbseq_element: ET.Element,
) -> str:
    """
    Extract one protein UID from a GBSeq record using its GBSeqid fields.
    """
    return extract_protein_uid_from_gbseq_element(gbseq_element=gbseq_element)


def _validate_saved_consolidated_xml_snapshot(
    *,
    xml_file_path: PathLike,
    expected_root_tag: str = "GBSet",
    expected_record_tag: str = "GBSeq",
    expected_record_count: Optional[int] = None,
) -> int:
    """
    Validate the structure of a saved consolidated XML snapshot in one pass.

    Returns:
        int:
            Number of direct child records found under the root element.
    """
    return validate_consolidated_xml_file(
        xml_file_path=xml_file_path,
        expected_root_tag=expected_root_tag,
        expected_record_tag=expected_record_tag,
        expected_record_count=expected_record_count,
    ).record_count


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
    return validate_consolidated_xml_file(
        xml_file_path=xml_file_path,
        expected_root_tag=expected_root_tag,
        expected_record_tag=expected_record_tag,
        expected_protein_uids=expected_protein_uids,
        validation_context="during UID validation",
    ).protein_uids


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

    manifest_payload["snapshot_format_version"] = UID_SNAPSHOT_FORMAT_VERSION
    manifest_payload["snapshot_file_name"] = "protein_uids.txt"
    manifest_payload["manifest_file_name"] = "manifest.json"
    manifest_payload["immutable_snapshot_directory_name"] = (
        immutable_snapshot_directory_name
    )
    manifest_payload["immutable_snapshot_relative_path"] = (
        immutable_snapshot_relative_path
    )

    return manifest_payload


SUPERSEDED_LATEST_DIRECTORY_PREFIX_SUFFIX = "_superseded_"
STAGED_LATEST_DIRECTORY_PREFIX_SUFFIX = "_staged_"


def _rename_directory_with_retries(
    *,
    source_directory: Path,
    destination_directory: Path,
    attempt_count: int = 5,
    initial_retry_delay_seconds: float = 0.2,
) -> bool:
    """
    Rename a directory, retrying while the filesystem still holds it.

    On Windows a directory containing a file another process has open cannot be
    renamed or removed until that handle closes. A virus scanner or indexer
    reading a freshly written multi-hundred megabyte artifact routinely holds
    one for a short while, so a bounded retry converts a hard failure into a
    brief wait.
    """
    for attempt_index in range(attempt_count):
        try:
            source_directory.replace(destination_directory)
            return True
        except OSError:
            if attempt_index == attempt_count - 1:
                return False

            time.sleep(initial_retry_delay_seconds * (2 ** attempt_index))

    return False


def _reserve_sibling_directory_name(
    *,
    parent_directory: Path,
    name_prefix: str,
) -> Path:
    """
    Reserve an unused sibling directory name for a rename target.
    """
    reserved_directory = Path(
        tempfile.mkdtemp(prefix=name_prefix, dir=parent_directory)
    )
    reserved_directory.rmdir()
    return reserved_directory


def _sweep_superseded_latest_directories(
    *,
    parent_directory: Path,
    latest_directory_name: str,
) -> None:
    """
    Remove leftovers a previous interrupted publication could not delete.
    """
    for leftover_directory in parent_directory.glob(
        f"{latest_directory_name}"
        f"{SUPERSEDED_LATEST_DIRECTORY_PREFIX_SUFFIX}*"
    ):
        shutil.rmtree(leftover_directory, ignore_errors=True)


def _replace_latest_directory(
    *,
    latest_directory: PathLike,
    files_to_copy: list[tuple[PathLike, str]],
) -> None:
    """
    Replace the convenience latest directory from a fully prepared temp copy.

    The previous directory is moved aside rather than deleted in place. Deleting
    first means a failure part way through leaves no latest/ at all, which is
    exactly what happens when a large artifact is briefly locked: the manifest
    is removed, the payload is not, and the directory is left in a state no
    reader can validate. Moving aside fails before anything is destroyed, and
    the old directory is restored if the new one cannot be put in place.
    """
    resolved_latest_directory = _as_path(latest_directory)
    parent_directory = resolved_latest_directory.parent
    parent_directory.mkdir(parents=True, exist_ok=True)

    _sweep_superseded_latest_directories(
        parent_directory=parent_directory,
        latest_directory_name=resolved_latest_directory.name,
    )

    staged_latest_directory = Path(
        tempfile.mkdtemp(
            prefix=(
                f"{resolved_latest_directory.name}"
                f"{STAGED_LATEST_DIRECTORY_PREFIX_SUFFIX}"
            ),
            dir=parent_directory,
        )
    )
    superseded_latest_directory: Optional[Path] = None

    try:
        for source_file_path, destination_file_name in files_to_copy:
            shutil.copy2(
                _as_path(source_file_path),
                staged_latest_directory / destination_file_name,
            )

        if resolved_latest_directory.exists():
            superseded_latest_directory = _reserve_sibling_directory_name(
                parent_directory=parent_directory,
                name_prefix=(
                    f"{resolved_latest_directory.name}"
                    f"{SUPERSEDED_LATEST_DIRECTORY_PREFIX_SUFFIX}"
                ),
            )
            if not _rename_directory_with_retries(
                source_directory=resolved_latest_directory,
                destination_directory=superseded_latest_directory,
            ):
                superseded_latest_directory = None
                raise RuntimeError(
                    "Unable to move the existing latest snapshot directory "
                    f"aside: {resolved_latest_directory}. The published "
                    "immutable snapshot is unaffected and the existing "
                    "latest/ directory was left intact. Close any process "
                    "holding those files open and retry the publication."
                )

        if not _rename_directory_with_retries(
            source_directory=staged_latest_directory,
            destination_directory=resolved_latest_directory,
        ):
            if superseded_latest_directory is not None:
                _rename_directory_with_retries(
                    source_directory=superseded_latest_directory,
                    destination_directory=resolved_latest_directory,
                )
                superseded_latest_directory = None

            raise RuntimeError(
                "Unable to publish the latest snapshot directory: "
                f"{resolved_latest_directory}. The published immutable "
                "snapshot is unaffected."
            )
    except Exception:
        if staged_latest_directory.exists():
            shutil.rmtree(staged_latest_directory, ignore_errors=True)
        raise
    finally:
        if (
            superseded_latest_directory is not None
            and superseded_latest_directory.exists()
        ):
            # Best effort: latest/ is already correct, so a leftover that is
            # still locked is swept by the next publication instead of failing
            # one that has otherwise succeeded.
            shutil.rmtree(superseded_latest_directory, ignore_errors=True)


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
    retrieval_telemetry: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Build the manifest payload persisted alongside one consolidated XML file.
    """
    return {
        "snapshot_format_version": XML_SNAPSHOT_FORMAT_VERSION,
        "artifact_type": "ncbi_protein_xml_snapshot",
        "database_name": fetch_result.database_name,
        "identifier_type": fetch_result.identifier_type,
        "search_query": source_uid_snapshot_manifest["search_query"],
        "translated_query": source_uid_snapshot_manifest.get("translated_query"),
        "retrieved_at_utc": fetch_result.retrieved_at_utc,
        "retmode": fetch_result.retmode,
        "rettype": fetch_result.rettype,
        "requested_protein_uid_count": fetch_result.requested_protein_uid_count,
        "normalized_protein_uid_count": fetch_result.normalized_protein_uid_count,
        "protein_uids_sha256": fetch_result.protein_uids_sha256,
        "batch_size": fetch_result.batch_size,
        "batch_count": fetch_result.batch_count,
        "reused_batch_count": fetch_result.reused_batch_count,
        "fetched_batch_count": fetch_result.fetched_batch_count,
        "request_delay_seconds": fetch_result.request_delay_seconds,
        "max_retry_attempts": fetch_result.max_retry_attempts,
        "request_policy": {
            "fetch_timeout_seconds": fetch_result.fetch_timeout_seconds,
            "batch_deadline_seconds": fetch_result.batch_deadline_seconds,
            "retry_backoff_initial_seconds": (
                fetch_result.retry_backoff_initial_seconds
            ),
            "retry_backoff_multiplier": fetch_result.retry_backoff_multiplier,
            "retry_backoff_max_seconds": fetch_result.retry_backoff_max_seconds,
            "rate_limit_backoff_seconds": fetch_result.rate_limit_backoff_seconds,
            "circuit_breaker_failure_threshold": (
                fetch_result.circuit_breaker_failure_threshold
            ),
            "circuit_breaker_cooldown_seconds": (
                fetch_result.circuit_breaker_cooldown_seconds
            ),
            "max_concurrent_requests": fetch_result.max_concurrent_requests,
            "max_request_starts_per_second": (
                fetch_result.max_request_starts_per_second
            ),
            "reuse_http_connection": fetch_result.reuse_http_connection,
        },
        "retrieval_telemetry": retrieval_telemetry,
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


@dataclass(frozen=True)
class SavedXmlSnapshot:
    """
    One published XML snapshot together with the state validated while writing.

    The publication path already parsed, hashed, and UID-checked the artifact it
    just wrote. Returning that state lets callers skip re-reading a multi-hundred
    megabyte file only to reconstruct facts that are already known.
    """
    snapshot_directory: Path
    manifest_file_path: Path
    protein_uids_file_path: Path
    xml_file_path: Path
    manifest: Dict[str, Any]
    protein_uids: list[str]
    consolidated_record_count: int
    xml_file_sha256: str

    def as_snapshot_payload(self) -> Dict[str, Any]:
        return {
            "snapshot_directory": self.snapshot_directory,
            "manifest_file_path": self.manifest_file_path,
            "protein_uids_file_path": self.protein_uids_file_path,
            "manifest": self.manifest,
            "protein_uids": self.protein_uids,
            "xml_file_path": self.xml_file_path,
        }


def save_ncbi_protein_xml_snapshot(
    *,
    fetch_result: NCBIProteinXmlFetchResult,
    snapshot_root_directory: PathLike,
    source_uid_snapshot_manifest: Dict[str, Any],
    source_uid_snapshot_manifest_file_path: PathLike,
    protein_uids: list[str],
    update_latest_directory: bool = True,
    telemetry: Optional[NCBIRetrievalTelemetry] = None,
) -> SavedXmlSnapshot:
    """
    Persist one immutable local snapshot of consolidated NCBI protein XML.

    XML is fetched in batches upstream for robustness, then streamed into one
    valid XML document with a single root and persisted as one file. The
    document is never held in memory: records are re-serialized straight into
    the destination file, the SHA-256 is accumulated over the same bytes as they
    are written, and one streaming pass afterwards validates structure, record
    count, and the exact protein UID multiset.
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
        for xml_batch in sorted(
            fetch_result.xml_batches,
            key=lambda batch: batch.batch_index,
        ):
            batch_metadata_records.append(
                {
                    "batch_index": xml_batch.batch_index,
                    "batch_start_index": xml_batch.batch_start_index,
                    "batch_end_index": xml_batch.batch_end_index,
                    "protein_uid_count": xml_batch.protein_uid_count,
                    "xml_payload_sha256": xml_batch.xml_payload_sha256,
                    "reused_from_workspace": xml_batch.reused_from_workspace,
                }
            )

        expected_record_count = fetch_result.normalized_protein_uid_count

        with _measure_local_phase(
            telemetry=telemetry,
            phase_name="consolidated_xml_write",
        ):
            consolidated_write_result = write_consolidated_xml_document(
                batch_payload_sources=_build_xml_batch_payload_sources(
                    fetch_result=fetch_result,
                ),
                output_file_path=consolidated_xml_file_path,
            )

        if consolidated_write_result.record_count != expected_record_count:
            raise RuntimeError(
                "Consolidated XML snapshot record count does not match the "
                "expected protein UID count. "
                f"Expected {expected_record_count}, got "
                f"{consolidated_write_result.record_count}."
            )

        with _measure_local_phase(
            telemetry=telemetry,
            phase_name="consolidated_xml_validation",
        ):
            consolidated_validation_result = validate_consolidated_xml_file(
                xml_file_path=consolidated_xml_file_path,
                expected_record_count=expected_record_count,
                expected_protein_uids=protein_uids,
            )

        validated_record_count = consolidated_validation_result.record_count
        consolidated_xml_file_sha256 = consolidated_write_result.sha256

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
            retrieval_telemetry=(
                telemetry.build_summary()
                if telemetry is not None
                else fetch_result.retrieval_telemetry
            ),
        )

        write_json_atomic(
            payload=manifest_payload,
            output_file_path=manifest_file_path,
        )
        immutable_snapshot_complete = True

        if update_latest_directory:
            with _measure_local_phase(
                telemetry=telemetry,
                phase_name="latest_directory_copy",
            ):
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

    return SavedXmlSnapshot(
        snapshot_directory=immutable_snapshot_directory,
        manifest_file_path=manifest_file_path,
        protein_uids_file_path=protein_uids_file_path,
        xml_file_path=consolidated_xml_file_path,
        manifest=manifest_payload,
        protein_uids=list(protein_uids),
        consolidated_record_count=validated_record_count,
        xml_file_sha256=consolidated_xml_file_sha256,
    )


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

    snapshot_format_version = manifest_payload.get("snapshot_format_version")
    if snapshot_format_version not in SUPPORTED_XML_SNAPSHOT_FORMAT_VERSIONS:
        raise RuntimeError(
            "Saved XML snapshot manifest format version is not supported. "
            f"Expected one of "
            f"{sorted(SUPPORTED_XML_SNAPSHOT_FORMAT_VERSIONS)}, got "
            f"{snapshot_format_version!r}."
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

    validated_record_count = validate_consolidated_xml_file(
        xml_file_path=resolved_xml_file_path,
        expected_record_count=expected_record_count,
        expected_protein_uids=protein_uids,
    ).record_count

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


def resolve_ncbi_protein_uid_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    search_query: Optional[str] = None,
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
    deduplicate_uids: bool = True,
    sort_uids: bool = True,
    page_size: int = 10_000,
    max_retry_attempts: int = 5,
    request_delay_seconds: Optional[float] = None,
    fetch_timeout_seconds: float = 30.0,
    request_deadline_seconds: float = 300.0,
    retry_backoff_initial_seconds: Optional[float] = None,
    retry_backoff_multiplier: float = 2.0,
    retry_backoff_max_seconds: float = 30.0,
    rate_limit_backoff_seconds: float = 5.0,
    reuse_http_connection: bool = False,
    verbose_logging: bool = False,
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
        fetch_timeout_seconds=fetch_timeout_seconds,
        request_deadline_seconds=request_deadline_seconds,
        retry_backoff_initial_seconds=retry_backoff_initial_seconds,
        retry_backoff_multiplier=retry_backoff_multiplier,
        retry_backoff_max_seconds=retry_backoff_max_seconds,
        rate_limit_backoff_seconds=rate_limit_backoff_seconds,
        reuse_http_connection=reuse_http_connection,
        verbose_logging=verbose_logging,
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
    ssl_ca_file: Optional[str] = None,
    ssl_ca_directory: Optional[str] = None,
    xml_batch_size: int = 500,
    max_retry_attempts: int = 5,
    xml_request_delay_seconds: Optional[float] = None,
    fetch_timeout_seconds: float = 30.0,
    batch_deadline_seconds: float = 300.0,
    retry_backoff_initial_seconds: Optional[float] = None,
    retry_backoff_multiplier: float = 2.0,
    retry_backoff_max_seconds: float = 30.0,
    rate_limit_backoff_seconds: float = 5.0,
    circuit_breaker_failure_threshold: int = 3,
    circuit_breaker_cooldown_seconds: float = 60.0,
    xml_rettype: Optional[str] = NCBI_PROTEIN_GBSEQ_XML_RETTYPE,
    max_concurrent_requests: int = 1,
    max_request_starts_per_second: Optional[float] = None,
    reuse_http_connection: bool = False,
    batch_workspace_directory: Optional[PathLike] = None,
    enable_batch_resume: bool = True,
    purge_batch_workspace_on_success: bool = True,
    validate_batch_protein_uids: bool = True,
    verbose_batch_logging: bool = False,
    ncbi_email: Optional[str] = None,
    ncbi_api_key: Optional[str] = None,
    update_latest_directory: bool = True,
) -> Dict[str, Any]:
    """
    Resolve the active raw XML snapshot payload according to the requested mode.

    If snapshot creation is needed, the XML retrieval is performed strictly from an
    already-frozen upstream UID snapshot. This function does not create the source
    UID snapshot; that upstream step must have been run previously.

    Circuit-breaker state is retained by the NCBI API infrastructure across
    repeated snapshot workflow invocations for the same threshold and cooldown.
    Callers configure the policy here but do not need to own its mutable state.

    Fetched batches are written to a workspace directory under the snapshot root.
    A run that fails part way through therefore re-fetches only the batches it is
    missing on the next attempt, and the workspace is removed once a snapshot has
    been published successfully.
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

        return load_latest_xml_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory,
        )

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        print("Latest XML snapshot is available. Reusing frozen snapshot.")
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

    resolved_batch_workspace_directory = (
        _as_path(batch_workspace_directory)
        if batch_workspace_directory is not None
        else resolved_snapshot_root_directory / BATCH_WORKSPACE_DIRECTORY_NAME
    )
    retrieval_telemetry = NCBIRetrievalTelemetry()

    fetch_result = fetch_ncbi_protein_xml_batches(
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        protein_uids=protein_uids,
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
        rate_limit_backoff_seconds=rate_limit_backoff_seconds,
        circuit_breaker_failure_threshold=circuit_breaker_failure_threshold,
        circuit_breaker_cooldown_seconds=circuit_breaker_cooldown_seconds,
        circuit_breaker=get_default_ncbi_xml_circuit_breaker(
            failure_threshold=circuit_breaker_failure_threshold,
            cooldown_seconds=circuit_breaker_cooldown_seconds,
        ),
        rettype=xml_rettype,
        max_concurrent_requests=max_concurrent_requests,
        max_request_starts_per_second=max_request_starts_per_second,
        reuse_http_connection=reuse_http_connection,
        batch_workspace_directory=resolved_batch_workspace_directory,
        enable_batch_resume=enable_batch_resume,
        validate_batch_protein_uids=validate_batch_protein_uids,
        verbose_batch_logging=verbose_batch_logging,
        telemetry=retrieval_telemetry,
    )

    saved_xml_snapshot = save_ncbi_protein_xml_snapshot(
        fetch_result=fetch_result,
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_uid_snapshot_manifest=source_uid_snapshot_manifest,
        source_uid_snapshot_manifest_file_path=source_uid_snapshot_manifest_file_path,
        protein_uids=protein_uids,
        update_latest_directory=update_latest_directory,
        telemetry=retrieval_telemetry,
    )

    if purge_batch_workspace_on_success:
        purge_batch_workspace_directory(
            workspace_directory=fetch_result.batch_workspace_directory,
        )

    # The snapshot was parsed, hashed, and UID-checked while it was written, so
    # reloading it here would only repeat those passes over the published file.
    resolved_snapshot_payload = saved_xml_snapshot.as_snapshot_payload()

    # The manifest is written before the latest/ copy runs, so the manifest copy
    # of the telemetry cannot describe publication itself. The payload carries
    # the complete run, including the phases that follow the manifest write.
    resolved_snapshot_payload["retrieval_telemetry"] = (
        retrieval_telemetry.build_summary()
    )

    return resolved_snapshot_payload
