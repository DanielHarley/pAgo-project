from __future__ import annotations

import shutil
import re
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, TypeAlias

from src.pago_pipeline.apaz_hmm_build import (
    APAZ_BUILD_PROTOCOL,
    APAZ_MODEL_SPECS,
    DEFAULT_APAZ_SEED_DIRECTORY,
    HMMBuildCallable,
    ApazHmmBuildResult,
    build_apaz_hmm_policy_payload,
    build_apaz_hmm_policy_sha256,
    build_apaz_hmms,
    describe_hmm_builder_configuration,
    get_pyhmmer_builder_identity,
    validate_apaz_hmm_build_inputs,
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

ARTIFACT_TYPE = "apaz_hmm_build"
SNAPSHOT_FORMAT_VERSION = "1.0"
DEFAULT_MANIFEST_FILE_NAME = "manifest.json"
_SNAPSHOT_DIRECTORY_QUERY_LITERAL = "apaz_hmm_build"
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class ApazHmmBuildSnapshotResult:
    snapshot_directory: Path
    snapshot_root_directory: Path
    manifest_file_path: Path
    hmm_file_path_by_model_id: Mapping[str, Path]
    build_result: ApazHmmBuildResult


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _current_utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _resolve_builder_identity(
    *,
    hmm_builder: HMMBuildCallable | None,
    builder_identity: str | None,
) -> str:
    if hmm_builder is not None:
        if builder_identity is None or not builder_identity.strip():
            raise ValueError(
                "builder_identity is required when an injected HMM builder is used."
            )
        return builder_identity
    return builder_identity or get_pyhmmer_builder_identity()


def _manifest_payload_from_build_result(
    *,
    build_result: ApazHmmBuildResult,
    snapshot_created_at_utc: str,
    snapshot_directory_name: str,
) -> dict[str, object]:
    models = []
    output_files: dict[str, dict[str, str]] = {}
    for artifact in build_result.hmm_artifacts:
        model_id = artifact.model_spec.model_id
        output_key = f"{model_id}_hmm"
        output_files[output_key] = {
            "file_name": artifact.hmm_file_path.name,
            "sha256": artifact.hmm_sha256,
            "file_size_bytes": artifact.hmm_file_size_bytes,
        }
        models.append(
            {
                "model_id": model_id,
                "hmm_name": artifact.model_spec.hmm_name,
                "seed_file_name": artifact.seed_file_path.name,
                "seed_sha256": artifact.seed_sha256,
                "hmm_file_name": artifact.hmm_file_path.name,
                "hmm_sha256": artifact.hmm_sha256,
                "hmm_file_size_bytes": artifact.hmm_file_size_bytes,
                "sequence_count": artifact.sequence_count,
                "alignment_length": artifact.alignment_length,
                "match_state_count": artifact.match_state_count,
                "model_seed_sequence_count": artifact.model_seed_sequence_count,
                "hmm_checksum": artifact.hmm_checksum,
            }
        )
    validation = build_result.input_validation
    uses_pyhmmer_builder = build_result.builder_identity.startswith("pyhmmer/")
    return {
        "artifact_type": ARTIFACT_TYPE,
        "snapshot_format_version": SNAPSHOT_FORMAT_VERSION,
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "immutable_snapshot_directory_name": snapshot_directory_name,
        "immutable_snapshot_relative_path": str(
            Path("snapshots") / snapshot_directory_name
        ),
        "manifest_file_name": DEFAULT_MANIFEST_FILE_NAME,
        "build_protocol": APAZ_BUILD_PROTOCOL,
        "source_seed_directory_name": validation.reference_directory.name,
        "source_seed_lock_file_name": validation.seed_lock_file_path.name,
        "source_seed_lock_sha256": validation.seed_lock_sha256,
        "source_partitions_file_name": validation.partitions_file_path.name,
        "source_partitions_sha256": validation.partitions_sha256,
        "required_partition": "BUILD",
        "builder_identity": build_result.builder_identity,
        "pyhmmer_version": (
            build_result.builder_identity.split("/", 1)[-1]
            if uses_pyhmmer_builder
            else None
        ),
        "builder_configuration": (
            describe_hmm_builder_configuration(
                random_seed=build_result.random_seed
            )
            if uses_pyhmmer_builder
            else None
        ),
        "random_seed": build_result.random_seed,
        "build_policy": build_apaz_hmm_policy_payload(
            builder_identity=build_result.builder_identity,
            random_seed=build_result.random_seed,
        ),
        "build_policy_sha256": build_result.build_policy_sha256,
        "model_count": len(models),
        "models": models,
        "output_files": output_files,
    }


def save_apaz_hmm_build_snapshot(
    *,
    snapshot_root_directory: PathLike,
    reference_directory: PathLike = DEFAULT_APAZ_SEED_DIRECTORY,
    hmm_builder: HMMBuildCallable | None = None,
    builder_identity: str | None = None,
    random_seed: int = 42,
    update_latest_directory: bool = True,
) -> ApazHmmBuildSnapshotResult:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    resolved_builder_identity = _resolve_builder_identity(
        hmm_builder=hmm_builder, builder_identity=builder_identity
    )
    validate_apaz_hmm_build_inputs(reference_directory=reference_directory)
    snapshot_created_at_utc = _current_utc_timestamp()
    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=snapshot_created_at_utc,
        search_query=_SNAPSHOT_DIRECTORY_QUERY_LITERAL,
    )
    immutable_snapshot_directory = (
        resolved_snapshot_root_directory / "snapshots" / snapshot_directory_name
    )
    immutable_snapshot_directory.mkdir(parents=True, exist_ok=False)
    snapshot_complete = False
    try:
        build_result = build_apaz_hmms(
            output_directory=immutable_snapshot_directory,
            reference_directory=reference_directory,
            hmm_builder=hmm_builder,
            builder_identity=resolved_builder_identity,
            random_seed=random_seed,
        )
        manifest_file_path = (
            immutable_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
        )
        write_json_atomic(
            payload=_manifest_payload_from_build_result(
                build_result=build_result,
                snapshot_created_at_utc=snapshot_created_at_utc,
                snapshot_directory_name=snapshot_directory_name,
            ),
            output_file_path=manifest_file_path,
        )
        snapshot_complete = True
        if update_latest_directory:
            _replace_latest_directory(
                latest_directory=resolved_snapshot_root_directory / "latest",
                files_to_copy=[
                    *[
                        (artifact.hmm_file_path, artifact.hmm_file_path.name)
                        for artifact in build_result.hmm_artifacts
                    ],
                    (manifest_file_path, DEFAULT_MANIFEST_FILE_NAME),
                ],
            )
    except Exception:
        if not snapshot_complete and immutable_snapshot_directory.exists():
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

    return ApazHmmBuildSnapshotResult(
        snapshot_directory=immutable_snapshot_directory,
        snapshot_root_directory=resolved_snapshot_root_directory,
        manifest_file_path=manifest_file_path,
        hmm_file_path_by_model_id={
            artifact.model_spec.model_id: artifact.hmm_file_path
            for artifact in build_result.hmm_artifacts
        },
        build_result=build_result,
    )


def _validate_loaded_apaz_hmm_build_payload(
    *,
    snapshot_directory: Path,
    manifest_payload: Mapping[str, object],
) -> dict[str, Path]:
    if manifest_payload.get("artifact_type") != ARTIFACT_TYPE:
        raise RuntimeError("Saved APAZ HMM snapshot artifact_type mismatch.")
    if manifest_payload.get("snapshot_format_version") != SNAPSHOT_FORMAT_VERSION:
        raise RuntimeError(
            "Saved APAZ HMM snapshot snapshot_format_version mismatch."
        )
    if manifest_payload.get("manifest_file_name") != DEFAULT_MANIFEST_FILE_NAME:
        raise RuntimeError("Saved APAZ HMM snapshot manifest_file_name mismatch.")
    if manifest_payload.get("source_seed_lock_file_name") != "seeds_lock.json":
        raise RuntimeError(
            "Saved APAZ HMM snapshot source_seed_lock_file_name mismatch."
        )
    if manifest_payload.get("source_partitions_file_name") != (
        "apaz_partitions.csv"
    ):
        raise RuntimeError(
            "Saved APAZ HMM snapshot source_partitions_file_name mismatch."
        )
    for field_name in (
        "source_seed_lock_sha256",
        "source_partitions_sha256",
        "build_policy_sha256",
    ):
        value = manifest_payload.get(field_name)
        if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
            raise RuntimeError(
                f"Saved APAZ HMM snapshot {field_name} is invalid."
            )
    if manifest_payload.get("required_partition") != "BUILD":
        raise RuntimeError("Saved APAZ HMM snapshot required_partition mismatch.")
    build_protocol = manifest_payload.get("build_protocol")
    if not isinstance(build_protocol, str) or not build_protocol.strip():
        raise RuntimeError("Saved APAZ HMM snapshot build_protocol is missing.")
    builder_identity = manifest_payload.get("builder_identity")
    if not isinstance(builder_identity, str) or not builder_identity.strip():
        raise RuntimeError("Saved APAZ HMM snapshot builder_identity is missing.")
    random_seed = manifest_payload.get("random_seed")
    if type(random_seed) is not int:
        raise RuntimeError("Saved APAZ HMM snapshot random_seed must be an integer.")
    expected_build_policy = build_apaz_hmm_policy_payload(
        builder_identity=builder_identity,
        random_seed=random_seed,
    )
    if manifest_payload.get("build_policy") != expected_build_policy:
        raise RuntimeError(
            "Saved APAZ HMM snapshot build_policy is inconsistent with its "
            "builder metadata."
        )
    expected_build_policy_sha256 = build_apaz_hmm_policy_sha256(
        builder_identity=builder_identity,
        random_seed=random_seed,
    )
    if (
        manifest_payload.get("build_policy_sha256")
        != expected_build_policy_sha256
    ):
        raise RuntimeError(
            "Saved APAZ HMM snapshot build_policy_sha256 is inconsistent."
        )
    models = manifest_payload.get("models")
    if not isinstance(models, list) or len(models) != len(APAZ_MODEL_SPECS):
        raise RuntimeError("Saved APAZ HMM snapshot must define six models.")
    if manifest_payload.get("model_count") != len(APAZ_MODEL_SPECS):
        raise RuntimeError("Saved APAZ HMM snapshot model_count mismatch.")
    expected_model_id_order = [spec.model_id for spec in APAZ_MODEL_SPECS]
    if [entry.get("model_id") for entry in models if isinstance(entry, Mapping)] != (
        expected_model_id_order
    ):
        raise RuntimeError(
            "Saved APAZ HMM snapshot models are not in canonical order."
        )
    model_entry_by_id: dict[str, Mapping[str, object]] = {}
    for entry in models:
        if not isinstance(entry, Mapping):
            raise RuntimeError("Saved APAZ model entries must be objects.")
        model_id = entry.get("model_id")
        if not isinstance(model_id, str) or model_id in model_entry_by_id:
            raise RuntimeError("Saved APAZ model IDs must be unique strings.")
        model_entry_by_id[model_id] = entry
    expected_model_ids = {spec.model_id for spec in APAZ_MODEL_SPECS}
    if set(model_entry_by_id) != expected_model_ids:
        raise RuntimeError("Saved APAZ HMM snapshot model set mismatch.")

    output_files = manifest_payload.get("output_files")
    if not isinstance(output_files, Mapping):
        raise RuntimeError("Saved APAZ HMM snapshot must define output_files.")
    expected_output_keys = {f"{model_id}_hmm" for model_id in expected_model_ids}
    if set(output_files) != expected_output_keys:
        raise RuntimeError("Saved APAZ HMM output file set mismatch.")

    resolved_paths: dict[str, Path] = {}
    spec_by_id = {spec.model_id: spec for spec in APAZ_MODEL_SPECS}
    for model_id in expected_model_ids:
        output_entry = output_files[f"{model_id}_hmm"]
        if not isinstance(output_entry, Mapping):
            raise RuntimeError("Saved APAZ HMM output entries must be objects.")
        file_name = output_entry.get("file_name")
        expected_sha256 = output_entry.get("sha256")
        if file_name != spec_by_id[model_id].hmm_file_name:
            raise RuntimeError(f"Saved APAZ HMM file name mismatch for {model_id}.")
        if not isinstance(expected_sha256, str):
            raise RuntimeError(f"Saved APAZ HMM hash is missing for {model_id}.")
        model_entry = model_entry_by_id[model_id]
        spec = spec_by_id[model_id]
        if (
            model_entry.get("hmm_file_name") != file_name
            or model_entry.get("hmm_sha256") != expected_sha256
        ):
            raise RuntimeError(
                f"Saved APAZ model and output metadata disagree for {model_id}."
            )
        if (
            model_entry.get("hmm_name") != spec.hmm_name
            or model_entry.get("seed_file_name") != spec.seed_file_name
        ):
            raise RuntimeError(
                f"Saved APAZ model identity metadata mismatch for {model_id}."
            )
        seed_sha256 = model_entry.get("seed_sha256")
        if (
            not isinstance(seed_sha256, str)
            or _SHA256_PATTERN.fullmatch(seed_sha256) is None
        ):
            raise RuntimeError(
                f"Saved APAZ seed hash is invalid for {model_id}."
            )
        if (
            type(model_entry.get("sequence_count")) is not int
            or model_entry.get("sequence_count", 0) < 1
            or type(model_entry.get("alignment_length")) is not int
            or model_entry.get("alignment_length", 0) < 1
        ):
            raise RuntimeError(
                f"Saved APAZ alignment dimensions are invalid for {model_id}."
            )
        if (
            type(model_entry.get("match_state_count")) is not int
            or model_entry.get("match_state_count", 0) < 1
        ):
            raise RuntimeError(
                f"Saved APAZ match_state_count is invalid for {model_id}."
            )
        if (
            type(model_entry.get("hmm_file_size_bytes")) is not int
            or model_entry.get("hmm_file_size_bytes", 0) < 1
        ):
            raise RuntimeError(
                f"Saved APAZ hmm_file_size_bytes is invalid for {model_id}."
            )
        recorded_output_size = output_entry.get("file_size_bytes")
        if recorded_output_size is not None and (
            recorded_output_size != model_entry.get("hmm_file_size_bytes")
        ):
            raise RuntimeError(
                f"Saved APAZ output and model file size disagree for {model_id}."
            )
        if (
            not isinstance(expected_sha256, str)
            or _SHA256_PATTERN.fullmatch(expected_sha256) is None
        ):
            raise RuntimeError(f"Saved APAZ HMM hash is invalid for {model_id}.")
        file_path = snapshot_directory / file_name
        if not file_path.exists():
            raise FileNotFoundError(f"Saved APAZ HMM file not found: {file_path}.")
        actual_sha256 = sha256_of_file(input_file_path=file_path)
        if actual_sha256 != expected_sha256:
            raise RuntimeError(
                f"Saved APAZ HMM hash mismatch for {model_id}. Expected "
                f"{expected_sha256}, got {actual_sha256}."
            )
        if file_path.stat().st_size < 1:
            raise RuntimeError(f"Saved APAZ HMM file is empty for {model_id}.")
        resolved_paths[model_id] = file_path
    return resolved_paths


def load_apaz_hmm_build_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
) -> dict[str, object]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError("Saved APAZ HMM manifest must be a JSON object.")
    paths_by_model = _validate_loaded_apaz_hmm_build_payload(
        snapshot_directory=resolved_snapshot_directory,
        manifest_payload=manifest_payload,
    )
    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "hmm_file_path_by_model_id": paths_by_model,
    }


def load_latest_apaz_hmm_build_snapshot(
    *,
    snapshot_root_directory: PathLike,
) -> dict[str, object]:
    return load_apaz_hmm_build_snapshot_by_directory(
        snapshot_directory=_as_path(snapshot_root_directory) / "latest"
    )


def latest_apaz_hmm_build_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
    reference_directory: PathLike = DEFAULT_APAZ_SEED_DIRECTORY,
    builder_identity: str | None = None,
    random_seed: int = 42,
) -> bool:
    latest_directory = _as_path(snapshot_root_directory) / "latest"
    manifest_file_path = latest_directory / DEFAULT_MANIFEST_FILE_NAME
    if not latest_directory.exists() or not manifest_file_path.exists():
        return False
    try:
        payload = load_apaz_hmm_build_snapshot_by_directory(
            snapshot_directory=latest_directory
        )
        manifest = payload["manifest"]
        if not isinstance(manifest, Mapping):
            return False
        validation = validate_apaz_hmm_build_inputs(
            reference_directory=reference_directory
        )
        if manifest.get("source_seed_lock_sha256") != validation.seed_lock_sha256:
            return False
        if manifest.get("source_partitions_sha256") != validation.partitions_sha256:
            return False
        current_seed_hashes = {
            artifact.model_spec.model_id: artifact.seed_sha256
            for artifact in validation.seed_artifacts
        }
        saved_models = manifest.get("models")
        if not isinstance(saved_models, list):
            return False
        saved_seed_hashes = {
            str(entry.get("model_id")): entry.get("seed_sha256")
            for entry in saved_models
            if isinstance(entry, Mapping)
        }
        if saved_seed_hashes != current_seed_hashes:
            return False
        if builder_identity is not None:
            expected_policy_sha256 = build_apaz_hmm_policy_sha256(
                builder_identity=builder_identity,
                random_seed=random_seed,
            )
            if manifest.get("build_policy_sha256") != expected_policy_sha256:
                return False
    except (FileNotFoundError, RuntimeError, OSError, ValueError):
        return False
    return True


def resolve_apaz_hmm_build_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    reference_directory: PathLike = DEFAULT_APAZ_SEED_DIRECTORY,
    hmm_builder: HMMBuildCallable | None = None,
    builder_identity: str | None = None,
    random_seed: int = 42,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    resolved_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_builder_identity = _resolve_builder_identity(
        hmm_builder=hmm_builder, builder_identity=builder_identity
    )
    latest_available = latest_apaz_hmm_build_snapshot_is_available(
        snapshot_root_directory=snapshot_root_directory,
        reference_directory=reference_directory,
        builder_identity=resolved_builder_identity,
        random_seed=random_seed,
    )
    if resolved_mode == SnapshotMode.reuse_latest:
        if not latest_available:
            raise FileNotFoundError(
                "No compatible latest APAZ HMM build snapshot is available."
            )
        return load_latest_apaz_hmm_build_snapshot(
            snapshot_root_directory=snapshot_root_directory
        )
    if resolved_mode == SnapshotMode.reuse_latest_or_create and latest_available:
        return load_latest_apaz_hmm_build_snapshot(
            snapshot_root_directory=snapshot_root_directory
        )
    if resolved_mode not in {
        SnapshotMode.create_new,
        SnapshotMode.reuse_latest_or_create,
    }:
        raise ValueError("Invalid snapshot_mode for APAZ HMM build.")
    result = save_apaz_hmm_build_snapshot(
        snapshot_root_directory=snapshot_root_directory,
        reference_directory=reference_directory,
        hmm_builder=hmm_builder,
        builder_identity=resolved_builder_identity,
        random_seed=random_seed,
        update_latest_directory=update_latest_directory,
    )
    return load_apaz_hmm_build_snapshot_by_directory(
        snapshot_directory=result.snapshot_directory
    )


def list_saved_apaz_hmm_build_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    return list_saved_snapshot_directories(
        snapshot_root_directory=snapshot_root_directory
    )


def get_most_recent_apaz_hmm_build_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory
    )
