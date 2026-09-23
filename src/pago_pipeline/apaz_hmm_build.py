from __future__ import annotations

import csv
import hashlib
import io
import json
import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Protocol, TypeAlias

from src.pago_pipeline.apaz_split_groups import (
    validate_apaz_partition_invariants,
    validate_apaz_seed_consistency,
)
from src.pago_pipeline.storage import (
    read_json_file,
    sha256_of_file,
    write_bytes_atomic,
)

PathLike: TypeAlias = str | Path

APAZ_BUILD_INPUT_FORMAT_VERSION = "1.0"
APAZ_SEED_LOCK_ARTIFACT_TYPE = "apaz_seed_lock"
APAZ_SEED_LOCK_FORMAT_VERSION = "2.0"
# APAZ HMM construction is pinned to a single PyHMMER release so that the
# serialized profiles stay bit-for-bit reproducible across machines.
EXPECTED_PYHMMER_VERSION = "0.12.3"
APAZ_BUILD_PROTOCOL = "apaz_hmm_build/1.0 pyhmmer.plan7.Builder.build_msa"
APAZ_SPLIT_GROUPS_DIRECTORY_NAME = "split_groups"
APAZ_BUILD_PARTITION = "BUILD"
APAZ_ALLOWED_PARTITIONS = frozenset(
    {"BUILD", "CALIBRATION", "FINAL_HOLDOUT"}
)
APAZ_PARTITIONS_FILE_NAME = "apaz_partitions.csv"
APAZ_SEED_LOCK_FILE_NAME = "seeds_lock.json"
DEFAULT_APAZ_SEED_DIRECTORY = (
    Path(__file__).resolve().parent / "resources" / "apaz_seed"
)


@dataclass(frozen=True)
class ApazModelSpec:
    model_id: str
    lock_subgroup: str
    seed_file_name: str
    hmm_file_name: str
    hmm_name: str


APAZ_MODEL_SPECS = (
    ApazModelSpec(
        model_id="global",
        lock_subgroup="GLOBAL",
        seed_file_name="apaz_global.sto",
        hmm_file_name="apaz_global.hmm",
        hmm_name="APAZ_GLOBAL",
    ),
    ApazModelSpec(
        model_id="Ia",
        lock_subgroup="IA",
        seed_file_name="apaz_Ia.sto",
        hmm_file_name="apaz_Ia.hmm",
        hmm_name="APAZ_IA",
    ),
    ApazModelSpec(
        model_id="Ib",
        lock_subgroup="IB",
        seed_file_name="apaz_Ib.sto",
        hmm_file_name="apaz_Ib.hmm",
        hmm_name="APAZ_IB",
    ),
    ApazModelSpec(
        model_id="IIa",
        lock_subgroup="IIA",
        seed_file_name="apaz_IIa.sto",
        hmm_file_name="apaz_IIa.hmm",
        hmm_name="APAZ_IIA",
    ),
    ApazModelSpec(
        model_id="IIb",
        lock_subgroup="IIB",
        seed_file_name="apaz_IIb.sto",
        hmm_file_name="apaz_IIb.hmm",
        hmm_name="APAZ_IIB",
    ),
    ApazModelSpec(
        model_id="III",
        lock_subgroup="III",
        seed_file_name="apaz_III.sto",
        hmm_file_name="apaz_III.hmm",
        hmm_name="APAZ_III",
    ),
)

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_STOCKHOLM_COORDINATE_SUFFIX_PATTERN = re.compile(r"/\d+-\d+$")


class HMMBuildCallable(Protocol):
    def __call__(
        self,
        *,
        seed_alignment_file_path: Path,
        output_hmm_file_path: Path,
        model_name: str,
        random_seed: int,
    ) -> object: ...


@dataclass(frozen=True)
class StockholmAlignmentSummary:
    sequence_identifiers: tuple[str, ...]
    accessions: tuple[str, ...]
    sequence_count: int
    alignment_length: int
    aligned_sequence_by_accession: Mapping[str, str]


@dataclass(frozen=True)
class ApazSeedArtifact:
    model_spec: ApazModelSpec
    seed_file_path: Path
    seed_sha256: str
    alignment_summary: StockholmAlignmentSummary


@dataclass(frozen=True)
class ApazBuildInputValidationResult:
    reference_directory: Path
    seed_lock_file_path: Path
    seed_lock_sha256: str
    partitions_file_path: Path
    partitions_sha256: str
    partition_by_accession: Mapping[str, str]
    homology_cluster_by_accession: Mapping[str, str]
    seed_artifacts: tuple[ApazSeedArtifact, ...]


@dataclass(frozen=True)
class HmmStructuralReport:
    match_state_count: int
    model_seed_sequence_count: int | None
    hmm_checksum: int | None
    file_size_bytes: int


@dataclass(frozen=True)
class ApazHmmArtifact:
    model_spec: ApazModelSpec
    seed_file_path: Path
    seed_sha256: str
    hmm_file_path: Path
    hmm_sha256: str
    sequence_count: int
    alignment_length: int
    match_state_count: int
    model_seed_sequence_count: int | None
    hmm_checksum: int | None
    hmm_file_size_bytes: int


@dataclass(frozen=True)
class ApazHmmBuildResult:
    input_validation: ApazBuildInputValidationResult
    output_directory: Path
    builder_identity: str
    random_seed: int
    build_policy_sha256: str
    hmm_artifacts: tuple[ApazHmmArtifact, ...]


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _require_builder_identity(builder_identity: object) -> str:
    if (
        not isinstance(builder_identity, str)
        or not builder_identity.strip()
        or builder_identity != builder_identity.strip()
    ):
        raise ValueError("builder_identity must be a non-empty trimmed string.")
    return builder_identity


def _require_deterministic_random_seed(random_seed: object) -> int:
    if type(random_seed) is not int or not 1 <= random_seed <= 0xFFFFFFFF:
        raise ValueError(
            "random_seed must be an integer from 1 through 4294967295 so HMM "
            "construction remains deterministic."
        )
    return random_seed


def _require_sha256(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
        raise RuntimeError(
            f"{field_name} must be a lowercase 64-character SHA-256 value."
        )
    return value


def _require_canonical_relative_path(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or not value:
        raise RuntimeError(f"{field_name} must be a non-empty relative path.")
    pure_path = PurePosixPath(value)
    if (
        pure_path.is_absolute()
        or ".." in pure_path.parts
        or pure_path.as_posix() != value
    ):
        raise RuntimeError(
            f"{field_name} must be a canonical POSIX relative path, got {value!r}."
        )
    return value


def _canonical_stockholm_accession(sequence_identifier: str) -> str:
    return _STOCKHOLM_COORDINATE_SUFFIX_PATTERN.sub("", sequence_identifier)


def _ungapped_sequence(aligned_sequence: str) -> str:
    return re.sub(r"[-._~]", "", aligned_sequence).upper()


def read_stockholm_alignment_summary(
    *,
    alignment_file_path: PathLike,
) -> StockholmAlignmentSummary:
    resolved_alignment_file_path = _as_path(alignment_file_path)
    sequence_fragments: dict[str, list[str]] = {}
    first_content_line: str | None = None
    terminator_seen = False

    with resolved_alignment_file_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            stripped_line = raw_line.strip()
            if not stripped_line:
                continue
            if first_content_line is None:
                first_content_line = stripped_line
            if stripped_line == "//":
                if terminator_seen:
                    raise RuntimeError(
                        "Stockholm alignment contains more than one // terminator: "
                        f"{resolved_alignment_file_path}."
                    )
                terminator_seen = True
                continue
            if terminator_seen:
                raise RuntimeError(
                    "Stockholm alignment contains content after its // terminator "
                    f"at {resolved_alignment_file_path}:{line_number}."
                )
            if stripped_line.startswith("#"):
                continue
            fields = stripped_line.split()
            if len(fields) != 2:
                raise RuntimeError(
                    "Invalid Stockholm sequence row at "
                    f"{resolved_alignment_file_path}:{line_number}."
                )
            sequence_identifier, aligned_fragment = fields[0], fields[1]
            if not sequence_identifier or not aligned_fragment:
                raise RuntimeError(
                    "Stockholm sequence identifiers and fragments must be "
                    f"non-empty at {resolved_alignment_file_path}:{line_number}."
                )
            sequence_fragments.setdefault(sequence_identifier, []).append(
                aligned_fragment
            )

    if first_content_line != "# STOCKHOLM 1.0":
        raise RuntimeError(
            f"Alignment is not a Stockholm document: {resolved_alignment_file_path}."
        )
    if not terminator_seen:
        raise RuntimeError(
            f"Stockholm alignment has no // terminator: {resolved_alignment_file_path}."
        )
    if not sequence_fragments:
        raise RuntimeError(
            f"Stockholm alignment contains no sequences: {resolved_alignment_file_path}."
        )

    aligned_sequences = {
        identifier: "".join(fragments)
        for identifier, fragments in sequence_fragments.items()
    }
    alignment_lengths = {len(sequence) for sequence in aligned_sequences.values()}
    if len(alignment_lengths) != 1:
        raise RuntimeError(
            "Stockholm alignment contains unequal aligned sequence lengths in "
            f"{resolved_alignment_file_path}."
        )

    sequence_identifiers = tuple(aligned_sequences)
    accessions = tuple(
        _canonical_stockholm_accession(identifier)
        for identifier in sequence_identifiers
    )
    if len(set(accessions)) != len(accessions):
        raise RuntimeError(
            "Stockholm alignment contains duplicate canonical accessions after "
            f"coordinate suffix removal: {resolved_alignment_file_path}."
        )
    aligned_sequence_by_accession = {
        _canonical_stockholm_accession(identifier): sequence
        for identifier, sequence in aligned_sequences.items()
    }
    return StockholmAlignmentSummary(
        sequence_identifiers=sequence_identifiers,
        accessions=accessions,
        sequence_count=len(sequence_identifiers),
        alignment_length=next(iter(alignment_lengths)),
        aligned_sequence_by_accession=aligned_sequence_by_accession,
    )


def _load_partition_table(
    *,
    partitions_file_path: Path,
) -> tuple[dict[str, str], dict[str, str]]:
    """Return (partition_by_accession, split_group_id_by_accession).

    The v2 table carries an explicit ``split_group_id`` (the MMseqs2 90/80
    connected component). ``homology_cluster_id`` is kept only for backward
    compatibility and must equal ``split_group_id`` on every row. Neither a
    ``split_group_id`` nor a ``sequence_sha256`` may appear in more than one
    partition.
    """
    with partitions_file_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required_columns = {
            "accession", "homology_cluster_id", "split_group_id",
            "partition", "sequence_sha256",
        }
        if reader.fieldnames is None or not required_columns.issubset(
            set(reader.fieldnames)
        ):
            raise RuntimeError(
                "APAZ partition table must contain accession, homology_cluster_id, "
                "split_group_id, partition, and sequence_sha256 columns."
            )

        partition_by_accession: dict[str, str] = {}
        split_group_by_accession: dict[str, str] = {}
        partitions_by_split_group: dict[str, set[str]] = {}
        partitions_by_sequence_sha256: dict[str, set[str]] = {}
        for row_number, row in enumerate(reader, start=2):
            accession = str(row.get("accession", "")).strip()
            legacy_cluster_id = str(row.get("homology_cluster_id", "")).strip()
            split_group_id = str(row.get("split_group_id", "")).strip()
            partition = str(row.get("partition", "")).strip()
            sequence_sha256 = str(row.get("sequence_sha256", "")).strip()
            if not accession or not split_group_id or not partition or not sequence_sha256:
                raise RuntimeError(
                    "APAZ partition rows require non-empty accession, split_group_id, "
                    f"partition, and sequence_sha256 values at row {row_number}."
                )
            if legacy_cluster_id != split_group_id:
                raise RuntimeError(
                    "APAZ partition table requires homology_cluster_id == "
                    f"split_group_id (row {row_number}: {legacy_cluster_id!r} != "
                    f"{split_group_id!r})."
                )
            if partition not in APAZ_ALLOWED_PARTITIONS:
                raise RuntimeError(
                    f"Unsupported APAZ partition {partition!r} at row {row_number}."
                )
            if accession in partition_by_accession:
                raise RuntimeError(
                    f"Duplicate APAZ partition accession {accession!r}."
                )
            partition_by_accession[accession] = partition
            split_group_by_accession[accession] = split_group_id
            partitions_by_split_group.setdefault(split_group_id, set()).add(partition)
            partitions_by_sequence_sha256.setdefault(sequence_sha256, set()).add(partition)

    if not partition_by_accession:
        raise RuntimeError("APAZ partition table contains no records.")
    missing_partitions = sorted(
        APAZ_ALLOWED_PARTITIONS - set(partition_by_accession.values())
    )
    if missing_partitions:
        raise RuntimeError(
            "APAZ partition table must represent every required partition. "
            "Missing=" + ", ".join(missing_partitions)
        )
    leaked_groups = sorted(
        split_group_id
        for split_group_id, partitions in partitions_by_split_group.items()
        if len(partitions) != 1
    )
    if leaked_groups:
        raise RuntimeError(
            "APAZ split_group partition leakage detected for split groups: "
            + ", ".join(leaked_groups)
        )
    leaked_sequences = sorted(
        sequence_sha256
        for sequence_sha256, partitions in partitions_by_sequence_sha256.items()
        if len(partitions) != 1
    )
    if leaked_sequences:
        raise RuntimeError(
            "APAZ sequence_sha256 partition leakage detected (an identical protein "
            "sequence in more than one partition): " + ", ".join(leaked_sequences)
        )
    return partition_by_accession, split_group_by_accession


def _lock_file_entries_by_path(
    *,
    lock_payload: Mapping[str, object],
) -> dict[str, Mapping[str, object]]:
    raw_entries = lock_payload.get("files")
    if not isinstance(raw_entries, Sequence) or isinstance(
        raw_entries, (str, bytes)
    ):
        raise RuntimeError("APAZ seeds lock files must be a JSON array.")
    entries_by_path: dict[str, Mapping[str, object]] = {}
    for index, raw_entry in enumerate(raw_entries):
        if not isinstance(raw_entry, Mapping):
            raise RuntimeError(f"APAZ seeds lock files[{index}] must be an object.")
        relative_path = _require_canonical_relative_path(
            raw_entry.get("path"),
            field_name=f"files[{index}].path",
        )
        if relative_path in entries_by_path:
            raise RuntimeError(
                f"Duplicate APAZ seed lock path {relative_path!r}."
            )
        entries_by_path[relative_path] = raw_entry
    return entries_by_path


def validate_apaz_hmm_build_inputs(
    *,
    reference_directory: PathLike = DEFAULT_APAZ_SEED_DIRECTORY,
) -> ApazBuildInputValidationResult:
    resolved_reference_directory = _as_path(reference_directory)
    seed_lock_file_path = resolved_reference_directory / APAZ_SEED_LOCK_FILE_NAME
    partitions_file_path = resolved_reference_directory / APAZ_PARTITIONS_FILE_NAME
    lock_payload = read_json_file(input_file_path=seed_lock_file_path)
    if not isinstance(lock_payload, Mapping):
        raise RuntimeError("APAZ seeds lock must deserialize into a JSON object.")
    if lock_payload.get("artifact_type") != APAZ_SEED_LOCK_ARTIFACT_TYPE:
        raise RuntimeError(
            "APAZ seeds lock artifact_type mismatch. Expected "
            f"{APAZ_SEED_LOCK_ARTIFACT_TYPE!r}."
        )
    if lock_payload.get("lock_format_version") != APAZ_SEED_LOCK_FORMAT_VERSION:
        raise RuntimeError(
            "APAZ seeds lock lock_format_version mismatch. Expected "
            f"{APAZ_SEED_LOCK_FORMAT_VERSION!r}."
        )
    curation_protocol_version = lock_payload.get("curation_protocol_version")
    if (
        not isinstance(curation_protocol_version, str)
        or not curation_protocol_version.strip()
        or curation_protocol_version != curation_protocol_version.strip()
    ):
        raise RuntimeError(
            "APAZ seeds lock curation_protocol_version must be a non-empty "
            "trimmed string."
        )

    entries_by_path = _lock_file_entries_by_path(lock_payload=lock_payload)
    expected_seed_file_names = {spec.seed_file_name for spec in APAZ_MODEL_SPECS}
    if set(entries_by_path) != expected_seed_file_names:
        missing = sorted(expected_seed_file_names - set(entries_by_path))
        undeclared = sorted(set(entries_by_path) - expected_seed_file_names)
        raise RuntimeError(
            "APAZ seeds lock must declare exactly the six required seed "
            f"alignments. Missing={missing}; undeclared={undeclared}."
        )

    raw_partitions_entry = lock_payload.get("partitions_file")
    if not isinstance(raw_partitions_entry, Mapping):
        raise RuntimeError("APAZ seeds lock must define partitions_file.")
    locked_partitions_path = _require_canonical_relative_path(
        raw_partitions_entry.get("path"),
        field_name="partitions_file.path",
    )
    if locked_partitions_path != APAZ_PARTITIONS_FILE_NAME:
        raise RuntimeError(
            "APAZ partitions_file.path must be "
            f"{APAZ_PARTITIONS_FILE_NAME!r}."
        )
    expected_partitions_sha256 = _require_sha256(
        raw_partitions_entry.get("sha256"),
        field_name="partitions_file.sha256",
    )
    actual_partitions_sha256 = sha256_of_file(
        input_file_path=partitions_file_path
    )
    if actual_partitions_sha256 != expected_partitions_sha256:
        raise RuntimeError(
            "APAZ partitions file hash mismatch. Expected "
            f"{expected_partitions_sha256}, got {actual_partitions_sha256}."
        )

    physical_seed_inputs = {
        path.name
        for path in resolved_reference_directory.iterdir()
        if path.is_file() and path.suffix.lower() == ".sto"
    }
    undeclared_seed_inputs = sorted(physical_seed_inputs - expected_seed_file_names)
    if undeclared_seed_inputs:
        raise RuntimeError(
            "Undeclared APAZ seed alignment inputs are present: "
            + ", ".join(undeclared_seed_inputs)
        )

    partition_by_accession, cluster_by_accession = _load_partition_table(
        partitions_file_path=partitions_file_path
    )

    # When the frozen B2 split-group resources are present (the committed
    # reference directory), defer to the canonical B2 library as the single
    # source of truth for every partition and seed-consistency invariant.
    # Synthetic fixtures that ship no split_groups/ tree skip this and rely on
    # the structural checks below.
    split_groups_directory = (
        resolved_reference_directory / APAZ_SPLIT_GROUPS_DIRECTORY_NAME
    )
    if split_groups_directory.is_dir():
        validate_apaz_partition_invariants(
            partitions_csv_path=partitions_file_path,
            split_groups_directory=split_groups_directory,
        )
        validate_apaz_seed_consistency(
            resource_directory=resolved_reference_directory
        )

    seed_artifacts: list[ApazSeedArtifact] = []
    for model_spec in APAZ_MODEL_SPECS:
        seed_file_path = resolved_reference_directory / model_spec.seed_file_name
        lock_entry = entries_by_path[model_spec.seed_file_name]
        expected_seed_sha256 = _require_sha256(
            lock_entry.get("sha256"),
            field_name=f"files[{model_spec.seed_file_name}].sha256",
        )
        actual_seed_sha256 = sha256_of_file(input_file_path=seed_file_path)
        if actual_seed_sha256 != expected_seed_sha256:
            raise RuntimeError(
                f"APAZ seed hash mismatch for {model_spec.seed_file_name}. "
                f"Expected {expected_seed_sha256}, got {actual_seed_sha256}."
            )
        if lock_entry.get("alphabet") != "amino":
            raise RuntimeError(
                f"APAZ seed {model_spec.seed_file_name} must declare alphabet "
                "as 'amino'."
            )
        if lock_entry.get("subgroup") != model_spec.lock_subgroup:
            raise RuntimeError(
                f"APAZ seed {model_spec.seed_file_name} subgroup mismatch."
            )

        alignment_summary = read_stockholm_alignment_summary(
            alignment_file_path=seed_file_path
        )
        if lock_entry.get("sequence_count") != alignment_summary.sequence_count:
            raise RuntimeError(
                f"APAZ seed {model_spec.seed_file_name} sequence_count mismatch."
            )
        if lock_entry.get("alignment_length") != alignment_summary.alignment_length:
            raise RuntimeError(
                f"APAZ seed {model_spec.seed_file_name} alignment_length mismatch."
            )

        missing_partition_accessions = sorted(
            set(alignment_summary.accessions) - set(partition_by_accession)
        )
        if missing_partition_accessions:
            raise RuntimeError(
                f"APAZ seed {model_spec.seed_file_name} contains accessions absent "
                "from apaz_partitions.csv: "
                + ", ".join(missing_partition_accessions)
            )
        non_build_accessions = sorted(
            accession
            for accession in alignment_summary.accessions
            if partition_by_accession[accession] != APAZ_BUILD_PARTITION
        )
        if non_build_accessions:
            raise RuntimeError(
                f"APAZ seed {model_spec.seed_file_name} contains non-BUILD "
                "accessions: "
                + ", ".join(non_build_accessions)
            )
        seed_artifacts.append(
            ApazSeedArtifact(
                model_spec=model_spec,
                seed_file_path=seed_file_path,
                seed_sha256=actual_seed_sha256,
                alignment_summary=alignment_summary,
            )
        )

    artifact_by_model_id = {
        artifact.model_spec.model_id: artifact for artifact in seed_artifacts
    }
    subgroup_artifacts = [
        artifact
        for artifact in seed_artifacts
        if artifact.model_spec.model_id != "global"
    ]
    subgroup_by_accession: dict[str, list[str]] = {}
    for artifact in subgroup_artifacts:
        for accession in artifact.alignment_summary.accessions:
            subgroup_by_accession.setdefault(accession, []).append(
                artifact.model_spec.model_id
            )
    overlapping_subgroup_accessions = {
        accession: model_ids
        for accession, model_ids in subgroup_by_accession.items()
        if len(model_ids) != 1
    }
    if overlapping_subgroup_accessions:
        raise RuntimeError(
            "APAZ subgroup seed alignments must be accession-disjoint. "
            f"Overlaps={overlapping_subgroup_accessions}."
        )

    global_summary = artifact_by_model_id["global"].alignment_summary
    global_accessions = set(global_summary.accessions)
    subgroup_accessions = set(subgroup_by_accession)
    if global_accessions != subgroup_accessions:
        raise RuntimeError(
            "APAZ global seed accession set must equal the exact union of the "
            "five subgroup seeds. "
            f"Only global={sorted(global_accessions - subgroup_accessions)}; "
            f"only subgroups={sorted(subgroup_accessions - global_accessions)}."
        )
    for artifact in subgroup_artifacts:
        for accession, aligned_sequence in (
            artifact.alignment_summary.aligned_sequence_by_accession.items()
        ):
            global_sequence = global_summary.aligned_sequence_by_accession[accession]
            if _ungapped_sequence(global_sequence) != _ungapped_sequence(
                aligned_sequence
            ):
                raise RuntimeError(
                    "APAZ global and subgroup seeds contain different ungapped "
                    f"sequences for accession {accession!r}."
                )

    return ApazBuildInputValidationResult(
        reference_directory=resolved_reference_directory,
        seed_lock_file_path=seed_lock_file_path,
        seed_lock_sha256=sha256_of_file(input_file_path=seed_lock_file_path),
        partitions_file_path=partitions_file_path,
        partitions_sha256=actual_partitions_sha256,
        partition_by_accession=partition_by_accession,
        homology_cluster_by_accession=cluster_by_accession,
        seed_artifacts=tuple(seed_artifacts),
    )


def get_pyhmmer_builder_identity() -> str:
    try:
        import pyhmmer
    except ImportError as error:
        raise RuntimeError(
            "PyHMMER is required for production HMM construction. Install the "
            "project's pinned runtime dependencies before creating HMMs."
        ) from error
    if pyhmmer.__version__ != EXPECTED_PYHMMER_VERSION:
        raise RuntimeError(
            "APAZ HMM construction is pinned to pyhmmer "
            f"{EXPECTED_PYHMMER_VERSION}; the active interpreter provides "
            f"pyhmmer {pyhmmer.__version__}."
        )
    return f"pyhmmer/{pyhmmer.__version__}"


def describe_hmm_builder_configuration(*, random_seed: int) -> dict[str, object]:
    """Return the PyHMMER ``Builder`` configuration used for APAZ construction.

    Recorded verbatim in the build snapshot manifest so that the exact profile
    HMM construction parameters travel with the artifacts.
    """
    resolved_random_seed = _require_deterministic_random_seed(random_seed)
    try:
        import pyhmmer
    except ImportError as error:
        raise RuntimeError(
            "PyHMMER is required to describe the HMM builder configuration."
        ) from error
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet, seed=resolved_random_seed)
    effective_number = builder.effective_number
    return {
        "tool": "pyhmmer.plan7.Builder.build_msa",
        "pyhmmer_version": pyhmmer.__version__,
        "alphabet": "amino",
        "architecture": builder.architecture,
        "weighting": builder.weighting,
        "effective_number": (
            effective_number
            if isinstance(effective_number, str)
            else float(effective_number)
        ),
        "prior_scheme": builder.prior_scheme,
        "gap_open_probability": float(builder.popen),
        "gap_extend_probability": float(builder.pextend),
        "random_seed": resolved_random_seed,
        "seed_msa_name_forced": True,
        "creation_time_cleared": True,
    }


def _verify_hmm_artifact_on_disk(
    *,
    hmm_file_path: Path,
    expected_model_name: str,
    structural_tolerance: float = 1e-4,
) -> HmmStructuralReport:
    """Reopen a freshly written HMM from disk and prove it is a single, well
    formed, amino-alphabet model carrying the expected name.

    The in-memory HMM object the builder produced is never trusted: a profile
    is accepted only once it survives a full parse from its closed on-disk
    artifact.
    """
    if not hmm_file_path.is_file():
        raise RuntimeError(f"HMM artifact was not written: {hmm_file_path}.")
    file_size_bytes = hmm_file_path.stat().st_size
    if file_size_bytes == 0:
        raise RuntimeError(f"HMM artifact is empty: {hmm_file_path}.")

    try:
        import pyhmmer
    except ImportError as error:
        raise RuntimeError(
            "PyHMMER is required to verify a constructed HMM artifact."
        ) from error

    try:
        with pyhmmer.plan7.HMMFile(str(hmm_file_path)) as hmm_file:
            first_model = hmm_file.read()
            if first_model is None:
                raise RuntimeError(
                    f"HMM artifact does not parse into a model: {hmm_file_path}."
                )
            if hmm_file.read() is not None:
                raise RuntimeError(
                    "HMM artifact contains more than one model: "
                    f"{hmm_file_path}."
                )
    except (ValueError, EOFError, OSError) as error:
        raise RuntimeError(
            f"HMM artifact failed to parse from disk ({hmm_file_path}): {error}"
        ) from error

    if not first_model.alphabet.is_amino():
        raise RuntimeError(
            f"HMM artifact is not an amino-alphabet model: {hmm_file_path}."
        )
    raw_name = first_model.name
    actual_name = (
        raw_name.decode("ascii")
        if isinstance(raw_name, (bytes, bytearray))
        else str(raw_name)
    )
    if actual_name != expected_model_name:
        raise RuntimeError(
            f"HMM artifact name mismatch for {hmm_file_path}: expected "
            f"{expected_model_name!r}, found {actual_name!r}."
        )
    try:
        first_model.validate(tolerance=structural_tolerance)
    except Exception as error:
        raise RuntimeError(
            "HMM artifact violates the Plan7 structural constraints "
            f"({hmm_file_path}): {error}"
        ) from error
    if first_model.M < 1:
        raise RuntimeError(
            f"HMM artifact has no match states: {hmm_file_path}."
        )

    return HmmStructuralReport(
        match_state_count=int(first_model.M),
        model_seed_sequence_count=(
            int(first_model.nseq) if first_model.nseq is not None else None
        ),
        hmm_checksum=(
            int(first_model.checksum)
            if first_model.checksum is not None
            else None
        ),
        file_size_bytes=int(file_size_bytes),
    )


def build_hmm_with_pyhmmer(
    *,
    seed_alignment_file_path: Path,
    output_hmm_file_path: Path,
    model_name: str,
    random_seed: int,
) -> None:
    try:
        import pyhmmer
    except ImportError as error:
        raise RuntimeError(
            "PyHMMER is required for production HMM construction. Install the "
            "project's pinned runtime dependencies before creating HMMs."
        ) from error

    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.MSAFile(
        str(seed_alignment_file_path),
        digital=True,
        alphabet=alphabet,
    ) as alignment_file:
        alignment = alignment_file.read()
    if alignment is None:
        raise RuntimeError(
            f"PyHMMER could not read an MSA from {seed_alignment_file_path}."
        )
    alignment.name = model_name.encode("ascii")

    builder = pyhmmer.plan7.Builder(alphabet, seed=int(random_seed))
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(alignment, background)
    hmm.name = model_name.encode("ascii")
    # Strip every machine- and invocation-specific provenance field so the
    # serialized profile is bit-for-bit reproducible: creation_time yields a
    # DATE line, command_line (sys.argv) yields a COM line.
    for volatile_attribute in ("creation_time", "command_line"):
        try:
            setattr(hmm, volatile_attribute, None)
        except (AttributeError, TypeError):
            pass

    serialized_hmm = io.BytesIO()
    hmm.write(serialized_hmm)
    hmm_payload = serialized_hmm.getvalue()
    if not hmm_payload:
        raise RuntimeError(
            f"PyHMMER serialized an empty HMM for {seed_alignment_file_path}."
        )
    for volatile_line_prefix in (b"\nCOM ", b"\nDATE "):
        if volatile_line_prefix in hmm_payload:
            raise RuntimeError(
                "PyHMMER serialized a non-reproducible provenance line "
                f"({volatile_line_prefix.strip().decode('ascii')}) into the HMM "
                f"for {seed_alignment_file_path}."
            )
    write_bytes_atomic(
        binary_payload=hmm_payload,
        output_file_path=output_hmm_file_path,
    )


def build_apaz_hmm_policy_payload(
    *,
    builder_identity: str,
    random_seed: int,
) -> dict[str, object]:
    resolved_builder_identity = _require_builder_identity(builder_identity)
    resolved_random_seed = _require_deterministic_random_seed(random_seed)
    return {
        "build_input_format_version": APAZ_BUILD_INPUT_FORMAT_VERSION,
        "builder_identity": resolved_builder_identity,
        "random_seed": resolved_random_seed,
        "required_partition": APAZ_BUILD_PARTITION,
        "models": [
            {
                "model_id": spec.model_id,
                "seed_file_name": spec.seed_file_name,
                "hmm_file_name": spec.hmm_file_name,
                "hmm_name": spec.hmm_name,
            }
            for spec in APAZ_MODEL_SPECS
        ],
    }


def build_apaz_hmm_policy_sha256(
    *,
    builder_identity: str,
    random_seed: int,
) -> str:
    serialized_payload = json.dumps(
        build_apaz_hmm_policy_payload(
            builder_identity=builder_identity,
            random_seed=random_seed,
        ),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    )
    return hashlib.sha256(serialized_payload.encode("utf-8")).hexdigest()


def build_apaz_hmms(
    *,
    output_directory: PathLike,
    reference_directory: PathLike = DEFAULT_APAZ_SEED_DIRECTORY,
    hmm_builder: HMMBuildCallable | None = None,
    builder_identity: str | None = None,
    random_seed: int = 42,
) -> ApazHmmBuildResult:
    resolved_random_seed = _require_deterministic_random_seed(random_seed)
    if hmm_builder is None:
        resolved_hmm_builder: HMMBuildCallable = build_hmm_with_pyhmmer
        resolved_builder_identity = _require_builder_identity(
            builder_identity or get_pyhmmer_builder_identity()
        )
    else:
        resolved_hmm_builder = hmm_builder
        if builder_identity is None:
            raise ValueError(
                "builder_identity is required when an injected HMM builder is used."
            )
        resolved_builder_identity = _require_builder_identity(builder_identity)

    input_validation = validate_apaz_hmm_build_inputs(
        reference_directory=reference_directory
    )
    resolved_output_directory = _as_path(output_directory)
    resolved_output_directory.mkdir(parents=True, exist_ok=True)
    target_paths = {
        artifact.model_spec.model_id: (
            resolved_output_directory / artifact.model_spec.hmm_file_name
        )
        for artifact in input_validation.seed_artifacts
    }
    preexisting_targets = [path for path in target_paths.values() if path.exists()]
    if preexisting_targets:
        raise FileExistsError(
            "Refusing to overwrite existing APAZ HMM outputs: "
            + ", ".join(str(path) for path in preexisting_targets)
        )

    created_paths: list[Path] = []
    hmm_artifacts: list[ApazHmmArtifact] = []
    try:
        for seed_artifact in input_validation.seed_artifacts:
            output_hmm_file_path = target_paths[
                seed_artifact.model_spec.model_id
            ]
            resolved_hmm_builder(
                seed_alignment_file_path=seed_artifact.seed_file_path,
                output_hmm_file_path=output_hmm_file_path,
                model_name=seed_artifact.model_spec.hmm_name,
                random_seed=resolved_random_seed,
            )
            if (
                not output_hmm_file_path.exists()
                or output_hmm_file_path.stat().st_size == 0
            ):
                raise RuntimeError(
                    "HMM builder did not create a non-empty output for "
                    f"{seed_artifact.model_spec.model_id}."
                )
            created_paths.append(output_hmm_file_path)
            # Mandatory post-build validation: reopen the closed artifact from
            # disk and prove it is a single, structurally valid, amino-alphabet
            # model with the expected name before it is recorded. The SHA-256 is
            # taken only after this check, over the same closed file.
            structural_report = _verify_hmm_artifact_on_disk(
                hmm_file_path=output_hmm_file_path,
                expected_model_name=seed_artifact.model_spec.hmm_name,
            )
            hmm_artifacts.append(
                ApazHmmArtifact(
                    model_spec=seed_artifact.model_spec,
                    seed_file_path=seed_artifact.seed_file_path,
                    seed_sha256=seed_artifact.seed_sha256,
                    hmm_file_path=output_hmm_file_path,
                    hmm_sha256=sha256_of_file(
                        input_file_path=output_hmm_file_path
                    ),
                    sequence_count=(
                        seed_artifact.alignment_summary.sequence_count
                    ),
                    alignment_length=(
                        seed_artifact.alignment_summary.alignment_length
                    ),
                    match_state_count=structural_report.match_state_count,
                    model_seed_sequence_count=(
                        structural_report.model_seed_sequence_count
                    ),
                    hmm_checksum=structural_report.hmm_checksum,
                    hmm_file_size_bytes=structural_report.file_size_bytes,
                )
            )
    except Exception:
        for created_path in created_paths:
            created_path.unlink(missing_ok=True)
        for target_path in target_paths.values():
            if target_path.exists() and target_path not in created_paths:
                target_path.unlink(missing_ok=True)
        raise

    return ApazHmmBuildResult(
        input_validation=input_validation,
        output_directory=resolved_output_directory,
        builder_identity=resolved_builder_identity,
        random_seed=resolved_random_seed,
        build_policy_sha256=build_apaz_hmm_policy_sha256(
            builder_identity=resolved_builder_identity,
            random_seed=resolved_random_seed,
        ),
        hmm_artifacts=tuple(hmm_artifacts),
    )
