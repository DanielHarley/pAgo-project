"""Persist and reuse query-recall measurements as immutable snapshots (ADR-0002).

The computation lives in ``query_reference_recall``; this module only pins it.
A snapshot records what was measured, against which panel and which retrieval,
under which methodology, with a hash for every written file.

Five fields decide whether a stored snapshot may be reused:
``panel_lock_sha256``, ``matching_strategy_sha256``,
``snapshot_format_version``, ``source_metadata_manifest_sha256`` and
``source_metadata_csv_sha256``. The panel's lock covers the panel CSV and
curation notes transitively, because the panel loader verifies both against the
lock on every load; the other two panel hashes are recorded as provenance.

The source metadata snapshot is resolved **once** per operation and the
resulting immutable directory is used throughout, so the compatibility check
and the computation cannot see different sources.

``latest/`` is a convenience copy. Its publication is **not atomic**: the
shared helper populates a temporary directory, removes the existing
``latest/`` and only then renames. A failure in that window leaves ``latest/``
absent, never partially populated. The immutable snapshot under
``snapshots/<name>/`` is already complete, hashed and manifested before
``latest/`` is touched, so it remains the source of truth and the next run
simply recomputes the convenience copy.

Writing goes through ``storage.write_csv_rows_atomic`` with explicit field
names and pre-formatted cells: the stdlib ``csv`` writer's line terminator is
platform-independent, whereas ``DataFrame.to_csv`` follows the platform and
renders a missing value as an empty field, which would make a hashed artifact
differ between Linux and Windows.
"""

from __future__ import annotations

import csv
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, TypeAlias

from src.pago_pipeline.ncbi_metadata_snapshot import (
    load_latest_metadata_snapshot,
    load_metadata_snapshot_by_directory,
)
from src.pago_pipeline.ncbi_snapshot import (
    SnapshotMode,
    _coerce_snapshot_mode,
    _replace_latest_directory,
    build_snapshot_directory_name,
)
from src.pago_pipeline.query_recall_reference_panel import (
    DEFAULT_PANEL_DIRECTORY,
    LOCK_FILE_NAME,
    QueryRecallReferencePanelValidation,
    validate_query_recall_reference_panel,
)
from src.pago_pipeline.query_reference_recall import (
    DEFAULT_RETRIEVED_ACCESSION_COLUMN,
    DEFAULT_RETRIEVED_PROTEIN_UID_COLUMN,
    DEFAULT_RETRIEVED_SEQUENCE_COLUMN,
    DETAIL_FIELDNAMES,
    SUMMARY_FIELDNAMES,
    QueryReferenceRecallResult,
    build_matching_strategy_payload,
    build_matching_strategy_sha256,
    compute_query_reference_recall,
)
from src.pago_pipeline.storage import (
    read_json_file,
    sha256_of_file,
    write_csv_rows_atomic,
    write_json_atomic,
)

PathLike: TypeAlias = str | Path

ARTIFACT_TYPE = "query_reference_recall"
# 2.0 is the first format written against the versioned query-recall panel
# resource. Historical 1.x snapshots described a different reference schema and
# different strata and are never reused.
SNAPSHOT_FORMAT_VERSION = "2.0"
DEFAULT_MANIFEST_FILE_NAME = "manifest.json"
DEFAULT_SUMMARY_FILE_NAME = "reference_recall_summary.csv"
DEFAULT_DETAIL_FILE_NAME = "reference_recall_detail.csv"

_SNAPSHOT_DIRECTORY_QUERY_LITERAL = "query_reference_recall"

_REQUIRED_OUTPUT_FILE_KEYS: tuple[str, ...] = ("summary_file", "detail_file")


@dataclass(frozen=True)
class QueryReferenceRecallSnapshotResult:
    snapshot_directory: Path
    snapshot_root_directory: Path
    manifest_file_path: Path
    summary_file_path: Path
    detail_file_path: Path
    reference_count: int
    same_accession_recovered_count: int
    equivalent_recovered_count: int
    sequence_evidence_complete: bool


@dataclass(frozen=True)
class _ResolvedSourceMetadata:
    """One resolution of the upstream metadata snapshot, reused throughout."""

    manifest: Mapping[str, object]
    manifest_file_path: Path
    csv_file_path: Path
    manifest_sha256: str
    csv_sha256: str
    csv_sha256_declared_upstream: bool


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _current_utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _csv_cell(value: object) -> str:
    """Fix the serialized form so written bytes are determined.

    Missing becomes an empty field, booleans become ``true``/``false``, floats
    use ``repr`` (shortest round-tripping form), and everything else is its
    text. ``bool`` is checked before the numeric branch because it is an
    ``int`` subclass.
    """
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return repr(value)
    return str(value)


def _write_rows(
    *,
    rows: Sequence[Mapping[str, object]],
    fieldnames: Sequence[str],
    output_file_path: Path,
) -> Path:
    return write_csv_rows_atomic(
        rows=[
            {name: _csv_cell(row.get(name)) for name in fieldnames} for row in rows
        ],
        fieldnames=list(fieldnames),
        output_file_path=output_file_path,
    )


def _resolve_source_metadata(
    *, source_metadata_snapshot_root_directory: PathLike
) -> _ResolvedSourceMetadata:
    """Resolve the upstream metadata snapshot once.

    The upstream loader verifies the metadata CSV against the manifest's
    ``csv_file_sha256`` only when that field is declared, so an upstream
    manifest without it is accepted here and recorded as such; the independent
    comparison below is what protects reuse in that case.
    """
    try:
        payload = load_latest_metadata_snapshot(
            snapshot_root_directory=_as_path(source_metadata_snapshot_root_directory)
        )
    except FileNotFoundError as error:
        raise FileNotFoundError(
            "No reusable metadata snapshot was found for query-recall "
            "computation. Run the upstream metadata extraction stage first."
        ) from error

    manifest = payload.get("manifest")
    manifest_file_path = payload.get("manifest_file_path")
    csv_file_path = payload.get("csv_file_path")
    if (
        not isinstance(manifest, Mapping)
        or not isinstance(manifest_file_path, (str, Path))
        or not isinstance(csv_file_path, (str, Path))
    ):
        raise RuntimeError(
            "Resolved metadata snapshot payload is missing its manifest or CSV "
            "path."
        )

    # latest/ only *selects* the source. It is a convenience copy another
    # process may republish, so the operation is pinned to the immutable
    # snapshot it names — and that snapshot is loaded and validated in its own
    # right, by the same per-directory loader, rather than inheriting latest/'s
    # validation. Otherwise a corrupt immutable CSV would be computed over
    # while latest/ still validated, and provenance fields could come from one
    # manifest while the recorded hash came from another.
    immutable_relative_path = manifest.get("immutable_snapshot_relative_path")
    if not isinstance(immutable_relative_path, str) or not immutable_relative_path:
        raise RuntimeError(
            "The selected metadata snapshot does not declare "
            "immutable_snapshot_relative_path, so no immutable source can be "
            "pinned for this measurement. Computing over the mutable latest/ "
            "copy is not an acceptable substitute."
        )

    immutable_directory = (
        _as_path(source_metadata_snapshot_root_directory) / immutable_relative_path
    )
    try:
        immutable_payload = load_metadata_snapshot_by_directory(
            snapshot_directory=immutable_directory
        )
    except FileNotFoundError as error:
        raise FileNotFoundError(
            "The immutable metadata snapshot named by latest/ is missing: "
            f"{immutable_directory}. Query recall is measured against an "
            "immutable source, never against the latest/ copy alone."
        ) from error

    immutable_manifest = immutable_payload.get("manifest")
    immutable_manifest_file_path = immutable_payload.get("manifest_file_path")
    immutable_csv_file_path = immutable_payload.get("csv_file_path")
    if (
        not isinstance(immutable_manifest, Mapping)
        or not isinstance(immutable_manifest_file_path, (str, Path))
        or not isinstance(immutable_csv_file_path, (str, Path))
    ):
        raise RuntimeError(
            "The immutable metadata snapshot payload is missing its manifest or "
            f"CSV path: {immutable_directory}."
        )

    # latest/ is written as a copy of the immutable snapshot, so any difference
    # means the selected identity is not the identity of the files about to be
    # read. Rejecting is the only way to avoid recording a provenance that
    # mixes two sources.
    if dict(immutable_manifest) != dict(manifest):
        raise RuntimeError(
            "The latest/ metadata manifest and the immutable manifest it names "
            f"({immutable_directory}) differ, so the selected snapshot's "
            "identity cannot be established. Republish latest/ from the "
            "immutable snapshot before measuring query recall."
        )

    resolved_manifest_file_path = _as_path(immutable_manifest_file_path)
    resolved_csv_file_path = _as_path(immutable_csv_file_path)
    return _ResolvedSourceMetadata(
        manifest=immutable_manifest,
        manifest_file_path=resolved_manifest_file_path,
        csv_file_path=resolved_csv_file_path,
        manifest_sha256=sha256_of_file(input_file_path=resolved_manifest_file_path),
        csv_sha256=sha256_of_file(input_file_path=resolved_csv_file_path),
        csv_sha256_declared_upstream=(
            immutable_manifest.get("csv_file_sha256") is not None
        ),
    )


def _reserve_immutable_snapshot_directory(
    *, snapshots_directory: Path, preferred_name: str
) -> tuple[Path, str]:
    """Claim a fresh immutable snapshot directory, never reusing an existing one.

    ``build_snapshot_directory_name`` resolves to the second, so two snapshots
    created within the same second would collide. An existing immutable
    snapshot must never be written into, so the name is disambiguated with a
    suffix instead; ``mkdir(exist_ok=False)`` is the atomic reservation. The
    directory name is one of the manifest fields that legitimately varies
    between runs, so this does not affect reproducibility of the measured
    output.
    """
    snapshots_directory.mkdir(parents=True, exist_ok=True)
    attempt = 1
    while True:
        candidate_name = (
            preferred_name if attempt == 1 else f"{preferred_name}__{attempt}"
        )
        candidate_directory = snapshots_directory / candidate_name
        try:
            candidate_directory.mkdir(exist_ok=False)
        except FileExistsError:
            attempt += 1
            continue
        return candidate_directory, candidate_name


def _panel_lock_sha256(panel_validation: QueryRecallReferencePanelValidation) -> str:
    return sha256_of_file(
        input_file_path=panel_validation.panel_directory / LOCK_FILE_NAME
    )


def _build_manifest_payload(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    panel_validation: QueryRecallReferencePanelValidation,
    source_metadata: _ResolvedSourceMetadata,
    recall_result: QueryReferenceRecallResult,
    output_file_path_by_key: Mapping[str, Path],
) -> dict[str, object]:
    lock_payload = panel_validation.lock_payload
    return {
        "artifact_type": ARTIFACT_TYPE,
        "snapshot_format_version": SNAPSHOT_FORMAT_VERSION,
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "manifest_file_name": DEFAULT_MANIFEST_FILE_NAME,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "matching_strategy": build_matching_strategy_payload(),
        "matching_strategy_sha256": build_matching_strategy_sha256(),
        # Panel identity. The lock hash is the invalidating one: the panel
        # loader verifies the CSV and the curation notes against the lock on
        # every load, so those two are provenance.
        "panel_lock_sha256": _panel_lock_sha256(panel_validation),
        "panel_file_sha256": panel_validation.panel_sha256,
        "curation_notes_sha256": panel_validation.curation_notes_sha256,
        "intended_use": lock_payload.get("intended_use"),
        "detector_validation_eligible": lock_payload.get(
            "detector_validation_eligible"
        ),
        "source_metadata_snapshot_relative_path": source_metadata.manifest.get(
            "immutable_snapshot_relative_path"
        ),
        "source_metadata_snapshot_directory_name": source_metadata.manifest.get(
            "immutable_snapshot_directory_name"
        ),
        "source_metadata_manifest_sha256": source_metadata.manifest_sha256,
        "source_metadata_csv_sha256": source_metadata.csv_sha256,
        "source_metadata_csv_sha256_declared_upstream": (
            source_metadata.csv_sha256_declared_upstream
        ),
        "search_query": source_metadata.manifest.get("search_query"),
        "translated_query": source_metadata.manifest.get("translated_query"),
        "reference_count": recall_result.reference_count,
        "same_accession_recovered_count": (
            recall_result.same_accession_recovered_count
        ),
        "equivalent_recovered_count": recall_result.equivalent_recovered_count,
        # None (JSON null), never 0.0, for a stratum with no references. Only
        # these two maps carry null: the status map holds strings and the
        # lower-bound map holds booleans.
        "stratum_same_accession_recall": dict(
            recall_result.stratum_same_accession_recall
        ),
        "stratum_equivalent_recall": dict(recall_result.stratum_equivalent_recall),
        "stratum_recall_status": dict(recall_result.stratum_recall_status),
        "stratum_equivalent_is_lower_bound": dict(
            recall_result.stratum_equivalent_is_lower_bound
        ),
        "retrieved_row_count": recall_result.retrieved_row_count,
        "retrieved_rows_with_sequence_count": (
            recall_result.retrieved_rows_with_sequence_count
        ),
        "sequence_evidence_complete": recall_result.sequence_evidence_complete,
        "retrieved_inadmissible_accession_count": (
            recall_result.retrieved_inadmissible_accession_count
        ),
        "retrieved_unattributable_row_count": (
            recall_result.retrieved_unattributable_row_count
        ),
        "duplicate_reference_sequence_groups": [
            list(group) for group in recall_result.duplicate_reference_sequence_groups
        ],
        "output_files": {
            key: {
                "file_name": file_path.name,
                "path": str(file_path),
                "sha256": sha256_of_file(input_file_path=file_path),
            }
            for key, file_path in output_file_path_by_key.items()
        },
    }


def _read_retrieved_metadata(csv_file_path: Path):
    import pandas as pd

    return pd.read_csv(
        csv_file_path,
        low_memory=False,
        usecols=lambda column_name: column_name
        in {
            DEFAULT_RETRIEVED_PROTEIN_UID_COLUMN,
            DEFAULT_RETRIEVED_ACCESSION_COLUMN,
            DEFAULT_RETRIEVED_SEQUENCE_COLUMN,
        },
        dtype=str,
    )


def save_query_reference_recall_snapshot(
    *,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: PathLike,
    panel_directory: PathLike = DEFAULT_PANEL_DIRECTORY,
    update_latest_directory: bool = True,
    source_metadata: Optional[_ResolvedSourceMetadata] = None,
    panel_validation: Optional[QueryRecallReferencePanelValidation] = None,
) -> QueryReferenceRecallSnapshotResult:
    """Compute and write one immutable snapshot.

    ``source_metadata`` and ``panel_validation`` let a caller that already
    resolved them pass them in, so one operation never resolves the source
    twice.
    """
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    if panel_validation is None:
        panel_validation = validate_query_recall_reference_panel(
            panel_directory=panel_directory
        )
    if source_metadata is None:
        source_metadata = _resolve_source_metadata(
            source_metadata_snapshot_root_directory=(
                source_metadata_snapshot_root_directory
            )
        )

    snapshot_created_at_utc = _current_utc_timestamp()
    (
        immutable_snapshot_directory,
        snapshot_directory_name,
    ) = _reserve_immutable_snapshot_directory(
        snapshots_directory=resolved_snapshot_root_directory / "snapshots",
        preferred_name=build_snapshot_directory_name(
            retrieved_at_utc=snapshot_created_at_utc,
            search_query=_SNAPSHOT_DIRECTORY_QUERY_LITERAL,
        ),
    )
    immutable_snapshot_relative_path = str(Path("snapshots") / snapshot_directory_name)
    immutable_snapshot_complete = False

    try:
        recall_result = compute_query_reference_recall(
            panel_validation=panel_validation,
            retrieved_metadata_dataframe=_read_retrieved_metadata(
                source_metadata.csv_file_path
            ),
        )

        output_file_path_by_key = {
            "summary_file": immutable_snapshot_directory / DEFAULT_SUMMARY_FILE_NAME,
            "detail_file": immutable_snapshot_directory / DEFAULT_DETAIL_FILE_NAME,
        }
        _write_rows(
            rows=recall_result.summary_rows,
            fieldnames=SUMMARY_FIELDNAMES,
            output_file_path=output_file_path_by_key["summary_file"],
        )
        _write_rows(
            rows=recall_result.detail_rows,
            fieldnames=DETAIL_FIELDNAMES,
            output_file_path=output_file_path_by_key["detail_file"],
        )

        manifest_file_path = immutable_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
        write_json_atomic(
            payload=_build_manifest_payload(
                snapshot_created_at_utc=snapshot_created_at_utc,
                immutable_snapshot_directory_name=snapshot_directory_name,
                immutable_snapshot_relative_path=immutable_snapshot_relative_path,
                panel_validation=panel_validation,
                source_metadata=source_metadata,
                recall_result=recall_result,
                output_file_path_by_key=output_file_path_by_key,
            ),
            output_file_path=manifest_file_path,
        )
        immutable_snapshot_complete = True

        if update_latest_directory:
            _replace_latest_directory(
                latest_directory=resolved_snapshot_root_directory / "latest",
                files_to_copy=[
                    (
                        output_file_path_by_key["summary_file"],
                        DEFAULT_SUMMARY_FILE_NAME,
                    ),
                    (output_file_path_by_key["detail_file"], DEFAULT_DETAIL_FILE_NAME),
                    (manifest_file_path, DEFAULT_MANIFEST_FILE_NAME),
                ],
            )
    except Exception:
        if not immutable_snapshot_complete and immutable_snapshot_directory.exists():
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

    return QueryReferenceRecallSnapshotResult(
        snapshot_directory=immutable_snapshot_directory,
        snapshot_root_directory=resolved_snapshot_root_directory,
        manifest_file_path=manifest_file_path,
        summary_file_path=output_file_path_by_key["summary_file"],
        detail_file_path=output_file_path_by_key["detail_file"],
        reference_count=recall_result.reference_count,
        same_accession_recovered_count=(
            recall_result.same_accession_recovered_count
        ),
        equivalent_recovered_count=recall_result.equivalent_recovered_count,
        sequence_evidence_complete=recall_result.sequence_evidence_complete,
    )


def _validate_loaded_payload(
    *,
    snapshot_directory: Path,
    manifest_payload: Mapping[str, object],
) -> dict[str, Path]:
    """Integrity of a saved result: artifact kind, methodology, and file hashes."""
    artifact_type = manifest_payload.get("artifact_type")
    if artifact_type != ARTIFACT_TYPE:
        raise RuntimeError(
            "Saved query-recall snapshot artifact_type mismatch. Expected "
            f"{ARTIFACT_TYPE!r}, got {artifact_type!r}."
        )
    snapshot_format_version = manifest_payload.get("snapshot_format_version")
    if snapshot_format_version != SNAPSHOT_FORMAT_VERSION:
        raise RuntimeError(
            "Saved query-recall snapshot snapshot_format_version mismatch. "
            f"Expected {SNAPSHOT_FORMAT_VERSION!r}, got "
            f"{snapshot_format_version!r}. A snapshot written under an earlier "
            "output schema is not reused."
        )
    if manifest_payload.get("matching_strategy_sha256") != (
        build_matching_strategy_sha256()
    ):
        raise RuntimeError(
            "Saved query-recall snapshot matching_strategy_sha256 mismatch; the "
            "matching methodology changed."
        )

    output_files = manifest_payload.get("output_files")
    if not isinstance(output_files, Mapping):
        raise RuntimeError(
            "Saved query-recall snapshot manifest must define output_files."
        )
    # Both outputs are required before any of them is validated: a manifest
    # that declares neither, or only one, describes an incomplete result and
    # must fail validation rather than surface later as a missing key.
    missing_entries = [key for key in _REQUIRED_OUTPUT_FILE_KEYS if key not in output_files]
    if missing_entries:
        raise RuntimeError(
            "Saved query-recall snapshot manifest output_files is missing the "
            f"required entries {missing_entries}; a result without both "
            f"{list(_REQUIRED_OUTPUT_FILE_KEYS)} is incomplete."
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
                f"Saved query-recall snapshot file not found: {file_path}."
            )
        if sha256_of_file(input_file_path=file_path) != expected_sha256:
            raise RuntimeError(
                f"Saved query-recall snapshot hash mismatch for {key!r}: "
                f"{file_path}."
            )
        resolved[key] = file_path
    return resolved


def _read_csv_rows(csv_file_path: Path) -> list[dict[str, str]]:
    """Read a written CSV back as text rows.

    ``csv.DictReader`` keeps every cell a string, so identifiers and hashes
    cannot be coerced into numbers on the way back in.
    """
    with csv_file_path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def load_query_reference_recall_snapshot_by_directory(
    *, snapshot_directory: PathLike
) -> dict[str, object]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved query-recall snapshot manifest must deserialize into a JSON "
            "object."
        )
    resolved = _validate_loaded_payload(
        snapshot_directory=resolved_snapshot_directory,
        manifest_payload=manifest_payload,
    )
    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "summary_file_path": resolved["summary_file"],
        "detail_file_path": resolved["detail_file"],
        "summary_rows": _read_csv_rows(resolved["summary_file"]),
        "detail_rows": _read_csv_rows(resolved["detail_file"]),
    }


def load_latest_query_reference_recall_snapshot(
    *, snapshot_root_directory: PathLike
) -> dict[str, object]:
    return load_query_reference_recall_snapshot_by_directory(
        snapshot_directory=_as_path(snapshot_root_directory) / "latest"
    )


def _latest_is_compatible(
    *,
    snapshot_root_directory: Path,
    source_metadata: Optional[_ResolvedSourceMetadata],
    panel_validation: Optional[QueryRecallReferencePanelValidation],
) -> bool:
    latest_directory = snapshot_root_directory / "latest"
    latest_manifest_file_path = latest_directory / DEFAULT_MANIFEST_FILE_NAME
    if not latest_directory.exists() or not latest_manifest_file_path.exists():
        return False
    try:
        manifest_payload = read_json_file(input_file_path=latest_manifest_file_path)
        if not isinstance(manifest_payload, Mapping):
            return False
        _validate_loaded_payload(
            snapshot_directory=latest_directory,
            manifest_payload=manifest_payload,
        )
        if panel_validation is not None:
            if manifest_payload.get("panel_lock_sha256") != _panel_lock_sha256(
                panel_validation
            ):
                return False
        if source_metadata is not None:
            if (
                manifest_payload.get("source_metadata_manifest_sha256")
                != source_metadata.manifest_sha256
            ):
                return False
            # Compared independently of whether the upstream manifest declared
            # its own csv_file_sha256: when it did not, this is the only check
            # that notices a changed metadata CSV.
            if (
                manifest_payload.get("source_metadata_csv_sha256")
                != source_metadata.csv_sha256
            ):
                return False
    except (FileNotFoundError, RuntimeError, OSError, ValueError):
        return False
    return True


def latest_query_reference_recall_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: Optional[PathLike] = None,
    panel_directory: Optional[PathLike] = None,
) -> bool:
    """Whether the stored latest result is intact and compatible.

    A mismatch in any invalidating field makes the stored result unusable; the
    requested snapshot mode, not this function, decides what happens next.
    """
    source_metadata = None
    if source_metadata_snapshot_root_directory is not None:
        try:
            source_metadata = _resolve_source_metadata(
                source_metadata_snapshot_root_directory=(
                    source_metadata_snapshot_root_directory
                )
            )
        except (FileNotFoundError, RuntimeError, OSError, ValueError):
            return False
    panel_validation = None
    if panel_directory is not None:
        panel_validation = validate_query_recall_reference_panel(
            panel_directory=panel_directory
        )
    return _latest_is_compatible(
        snapshot_root_directory=_as_path(snapshot_root_directory),
        source_metadata=source_metadata,
        panel_validation=panel_validation,
    )


def resolve_query_reference_recall_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    source_metadata_snapshot_root_directory: PathLike,
    panel_directory: PathLike = DEFAULT_PANEL_DIRECTORY,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    """Reuse a compatible stored result, or create one, according to the mode.

    ``create_new`` always creates. ``reuse_latest`` requires a compatible
    stored result and raises when there is none — it never falls back to
    creating. ``reuse_latest_or_create`` reuses when compatible and creates
    otherwise.

    An invalid input raises before any of that: the panel preconditions and the
    upstream metadata snapshot are validated first, in every mode.
    """
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    if resolved_snapshot_mode not in {
        SnapshotMode.create_new,
        SnapshotMode.reuse_latest,
        SnapshotMode.reuse_latest_or_create,
    }:
        raise ValueError(
            "Invalid snapshot_mode. Expected one of: 'create_new', "
            "'reuse_latest', 'reuse_latest_or_create'."
        )
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)

    panel_validation = validate_query_recall_reference_panel(
        panel_directory=panel_directory
    )
    source_metadata = _resolve_source_metadata(
        source_metadata_snapshot_root_directory=(
            source_metadata_snapshot_root_directory
        )
    )

    if resolved_snapshot_mode in {
        SnapshotMode.reuse_latest,
        SnapshotMode.reuse_latest_or_create,
    }:
        if _latest_is_compatible(
            snapshot_root_directory=resolved_snapshot_root_directory,
            source_metadata=source_metadata,
            panel_validation=panel_validation,
        ):
            if resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create:
                print(
                    "Latest query-recall snapshot is compatible. Reusing frozen "
                    "snapshot."
                )
            return load_latest_query_reference_recall_snapshot(
                snapshot_root_directory=resolved_snapshot_root_directory
            )
        if resolved_snapshot_mode == SnapshotMode.reuse_latest:
            raise FileNotFoundError(
                "No compatible latest query-recall snapshot was found for the "
                "current panel and metadata snapshot. Run the workflow once "
                "with snapshot_mode='create_new' before using 'reuse_latest'."
            )

    saved_result = save_query_reference_recall_snapshot(
        snapshot_root_directory=resolved_snapshot_root_directory,
        source_metadata_snapshot_root_directory=(
            source_metadata_snapshot_root_directory
        ),
        panel_directory=panel_directory,
        update_latest_directory=update_latest_directory,
        source_metadata=source_metadata,
        panel_validation=panel_validation,
    )
    return load_query_reference_recall_snapshot_by_directory(
        snapshot_directory=saved_result.snapshot_directory
    )
