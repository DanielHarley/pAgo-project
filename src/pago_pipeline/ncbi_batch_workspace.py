from __future__ import annotations

import hashlib
import shutil
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

from src.pago_pipeline.ncbi_xml_stream import validate_xml_batch_payload
from src.pago_pipeline.storage import (
    read_json_file,
    sha256_of_lines,
    write_bytes_atomic,
    write_json_atomic,
)

PathLike = Union[str, Path]

BATCH_PLAN_FILE_NAME = "batch_plan.json"
BATCH_PLAN_FORMAT_VERSION = "1.0"
BATCH_PAYLOAD_DIRECTORY_NAME = "batches"

_PLAN_IDENTITY_FIELD_NAMES = (
    "protein_uids_sha256",
    "protein_uid_count",
    "batch_size",
    "database_name",
    "retmode",
    "rettype",
)


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _utc_now_string() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


@dataclass(frozen=True)
class PlannedXmlBatch:
    """
    One immutable unit of the XML retrieval plan.
    """
    batch_index: int
    batch_start_index: int
    batch_end_index: int
    protein_uids: List[str]
    protein_uids_sha256: str

    @property
    def protein_uid_count(self) -> int:
        return len(self.protein_uids)


@dataclass(frozen=True)
class ReusableXmlBatch:
    """
    One previously fetched batch that revalidated against the current plan.
    """
    payload_file_path: Path
    xml_payload_sha256: str
    record_count: int
    response_byte_count: int
    round_trip_latency_seconds: Optional[float]
    attempt_count: Optional[int]


def build_xml_batch_plan(
    *,
    protein_uids: Sequence[str],
    batch_size: int,
) -> List[PlannedXmlBatch]:
    """
    Split a frozen protein UID list into deterministic, index-addressed batches.
    """
    if batch_size <= 0:
        raise ValueError("batch_size must be a positive integer.")

    planned_batches: List[PlannedXmlBatch] = []

    for batch_index, batch_start_index in enumerate(
        range(0, len(protein_uids), batch_size),
        start=1,
    ):
        current_batch_protein_uids = list(
            protein_uids[batch_start_index: batch_start_index + batch_size]
        )
        planned_batches.append(
            PlannedXmlBatch(
                batch_index=batch_index,
                batch_start_index=batch_start_index,
                batch_end_index=(
                    batch_start_index + len(current_batch_protein_uids) - 1
                ),
                protein_uids=current_batch_protein_uids,
                protein_uids_sha256=sha256_of_lines(
                    text_lines=current_batch_protein_uids,
                    deduplicate_lines_preserving_order=False,
                    sort_lines=False,
                ),
            )
        )

    return planned_batches


class NCBIXmlBatchWorkspace:
    """
    Own the on-disk lifecycle of fetched XML batches for one retrieval run.

    The workspace serves two purposes that the plan keeps deliberately separate.
    It is the spill target that keeps fetched payloads out of process memory,
    and, when resume is enabled, it is the durable record that lets a failed run
    re-fetch only the batches it is actually missing. A batch is reused only
    when its recorded identity, its stored payload hash, and a full structural
    and UID revalidation of the payload all agree with the current plan.
    """

    def __init__(
        self,
        *,
        workspace_directory: Optional[PathLike],
        plan_identity: Dict[str, Any],
        request_policy: Optional[Dict[str, Any]] = None,
        enable_resume: bool = True,
        ephemeral_parent_directory: Optional[PathLike] = None,
    ) -> None:
        self._plan_identity = dict(plan_identity)
        self._request_policy = dict(request_policy or {})
        self._is_ephemeral = workspace_directory is None
        self._enable_resume = enable_resume and not self._is_ephemeral

        if workspace_directory is not None:
            self._workspace_directory = _as_path(workspace_directory).expanduser()
        else:
            resolved_parent_directory = (
                _as_path(ephemeral_parent_directory)
                if ephemeral_parent_directory is not None
                else None
            )
            if resolved_parent_directory is not None:
                resolved_parent_directory.mkdir(parents=True, exist_ok=True)

            self._workspace_directory = Path(
                tempfile.mkdtemp(
                    prefix="ncbi_xml_batches_",
                    dir=(
                        str(resolved_parent_directory)
                        if resolved_parent_directory is not None
                        else None
                    ),
                )
            )

        self._payload_directory = (
            self._workspace_directory / BATCH_PAYLOAD_DIRECTORY_NAME
        )
        self.reused_batch_count = 0
        self.purged_batch_count = 0
        self.plan_was_reset = False

    @property
    def workspace_directory(self) -> Path:
        return self._workspace_directory

    @property
    def resume_enabled(self) -> bool:
        return self._enable_resume

    @property
    def is_ephemeral(self) -> bool:
        return self._is_ephemeral

    @property
    def batch_plan_file_path(self) -> Path:
        return self._workspace_directory / BATCH_PLAN_FILE_NAME

    def batch_payload_file_paths(self, *, batch_index: int) -> tuple[Path, Path]:
        batch_file_stem = f"batch_{batch_index:06d}"
        return (
            self._payload_directory / f"{batch_file_stem}.xml",
            self._payload_directory / f"{batch_file_stem}.json",
        )

    def _current_plan_identity(self) -> Dict[str, Any]:
        return {
            field_name: self._plan_identity.get(field_name)
            for field_name in _PLAN_IDENTITY_FIELD_NAMES
        }

    def _saved_plan_is_compatible(self, *, saved_plan: Dict[str, Any]) -> bool:
        if saved_plan.get("batch_plan_format_version") != BATCH_PLAN_FORMAT_VERSION:
            return False

        saved_identity = saved_plan.get("plan_identity")
        if not isinstance(saved_identity, dict):
            return False

        return saved_identity == self._current_plan_identity()

    def open(self, *, planned_batches: Sequence[PlannedXmlBatch]) -> None:
        """
        Prepare the workspace, discarding any state incompatible with this plan.
        """
        self._workspace_directory.mkdir(parents=True, exist_ok=True)

        if self._enable_resume and self.batch_plan_file_path.exists():
            try:
                saved_plan = read_json_file(
                    input_file_path=self.batch_plan_file_path,
                )
            except Exception:
                saved_plan = None

            if not isinstance(saved_plan, dict) or not self._saved_plan_is_compatible(
                saved_plan=saved_plan,
            ):
                self._reset_payload_directory()
                self.plan_was_reset = True

        self._payload_directory.mkdir(parents=True, exist_ok=True)

        if not self._enable_resume:
            return

        write_json_atomic(
            payload={
                "batch_plan_format_version": BATCH_PLAN_FORMAT_VERSION,
                "created_at_utc": _utc_now_string(),
                "plan_identity": self._current_plan_identity(),
                "request_policy": self._request_policy,
                "batch_count": len(planned_batches),
                "batches": [
                    {
                        "batch_index": planned_batch.batch_index,
                        "batch_start_index": planned_batch.batch_start_index,
                        "batch_end_index": planned_batch.batch_end_index,
                        "protein_uid_count": planned_batch.protein_uid_count,
                        "protein_uids_sha256": planned_batch.protein_uids_sha256,
                    }
                    for planned_batch in planned_batches
                ],
            },
            output_file_path=self.batch_plan_file_path,
        )

    def _reset_payload_directory(self) -> None:
        if self._payload_directory.exists():
            shutil.rmtree(self._payload_directory, ignore_errors=True)

    def _discard_batch_files(self, *, batch_index: int) -> None:
        payload_file_path, metadata_file_path = self.batch_payload_file_paths(
            batch_index=batch_index,
        )
        removed_any_file = False
        for stale_file_path in (payload_file_path, metadata_file_path):
            if stale_file_path.exists():
                stale_file_path.unlink(missing_ok=True)
                removed_any_file = True

        if removed_any_file:
            self.purged_batch_count += 1

    def load_reusable_batch(
        self,
        *,
        planned_batch: PlannedXmlBatch,
        validate_protein_uids: bool = True,
    ) -> Optional[ReusableXmlBatch]:
        """
        Return a previously fetched batch only if it fully revalidates.
        """
        if not self._enable_resume:
            return None

        payload_file_path, metadata_file_path = self.batch_payload_file_paths(
            batch_index=planned_batch.batch_index,
        )
        if not payload_file_path.exists() or not metadata_file_path.exists():
            return None

        try:
            batch_metadata = read_json_file(input_file_path=metadata_file_path)
        except Exception:
            self._discard_batch_files(batch_index=planned_batch.batch_index)
            return None

        if not isinstance(batch_metadata, dict) or (
            batch_metadata.get("batch_index") != planned_batch.batch_index
            or batch_metadata.get("batch_start_index")
            != planned_batch.batch_start_index
            or batch_metadata.get("batch_end_index") != planned_batch.batch_end_index
            or batch_metadata.get("protein_uid_count")
            != planned_batch.protein_uid_count
            or batch_metadata.get("protein_uids_sha256")
            != planned_batch.protein_uids_sha256
        ):
            self._discard_batch_files(batch_index=planned_batch.batch_index)
            return None

        try:
            xml_payload_bytes = payload_file_path.read_bytes()
            xml_payload_sha256 = hashlib.sha256(xml_payload_bytes).hexdigest()
            if batch_metadata.get("xml_payload_sha256") != xml_payload_sha256:
                self._discard_batch_files(batch_index=planned_batch.batch_index)
                return None

            validation_result = validate_xml_batch_payload(
                xml_payload_bytes=xml_payload_bytes,
                expected_protein_uids=planned_batch.protein_uids,
                validate_protein_uids=validate_protein_uids,
            )
        except Exception:
            self._discard_batch_files(batch_index=planned_batch.batch_index)
            return None

        if validation_result.record_count != planned_batch.protein_uid_count:
            self._discard_batch_files(batch_index=planned_batch.batch_index)
            return None

        self.reused_batch_count += 1

        return ReusableXmlBatch(
            payload_file_path=payload_file_path,
            xml_payload_sha256=xml_payload_sha256,
            record_count=validation_result.record_count,
            response_byte_count=len(xml_payload_bytes),
            round_trip_latency_seconds=batch_metadata.get(
                "round_trip_latency_seconds"
            ),
            attempt_count=batch_metadata.get("attempt_count"),
        )

    def store_batch(
        self,
        *,
        planned_batch: PlannedXmlBatch,
        xml_payload_bytes: bytes,
        xml_payload_sha256: str,
        record_count: int,
        round_trip_latency_seconds: Optional[float] = None,
        attempt_count: Optional[int] = None,
        fetch_timeout_seconds: Optional[float] = None,
        batch_deadline_seconds: Optional[float] = None,
    ) -> Path:
        """
        Persist one validated batch payload atomically and record its identity.
        """
        payload_file_path, metadata_file_path = self.batch_payload_file_paths(
            batch_index=planned_batch.batch_index,
        )

        write_bytes_atomic(
            binary_payload=xml_payload_bytes,
            output_file_path=payload_file_path,
        )

        if self._enable_resume:
            write_json_atomic(
                payload={
                    "batch_index": planned_batch.batch_index,
                    "batch_start_index": planned_batch.batch_start_index,
                    "batch_end_index": planned_batch.batch_end_index,
                    "protein_uid_count": planned_batch.protein_uid_count,
                    "protein_uids_sha256": planned_batch.protein_uids_sha256,
                    "xml_payload_sha256": xml_payload_sha256,
                    "record_count": record_count,
                    "response_byte_count": len(xml_payload_bytes),
                    "round_trip_latency_seconds": round_trip_latency_seconds,
                    "attempt_count": attempt_count,
                    "fetch_timeout_seconds": fetch_timeout_seconds,
                    "batch_deadline_seconds": batch_deadline_seconds,
                    "recorded_at_utc": _utc_now_string(),
                },
                output_file_path=metadata_file_path,
            )

        return payload_file_path

    def purge(self) -> None:
        """
        Remove the whole workspace directory tree.
        """
        if self._workspace_directory.exists():
            shutil.rmtree(self._workspace_directory, ignore_errors=True)


def purge_batch_workspace_directory(
    *,
    workspace_directory: Optional[PathLike],
) -> None:
    """
    Remove a batch workspace directory if it still exists.
    """
    if workspace_directory is None:
        return

    resolved_workspace_directory = _as_path(workspace_directory)
    if resolved_workspace_directory.exists():
        shutil.rmtree(resolved_workspace_directory, ignore_errors=True)
