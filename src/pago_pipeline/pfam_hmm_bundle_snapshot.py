from __future__ import annotations

import hashlib
import importlib
import json
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, TypeAlias

from src.pago_pipeline.ncbi_snapshot import (
    SnapshotMode,
    _coerce_snapshot_mode,
    _replace_latest_directory,
    build_snapshot_directory_name,
    get_most_recent_snapshot_directory,
    list_saved_snapshot_directories,
)
from src.pago_pipeline.pfam_hmm_bundle import (
    CONCATENATION_POLICY_VERSION,
    DEFAULT_BUNDLE_FILE_NAME,
    DEFAULT_REGISTRY_FILE_PATH,
    PRESSED_FILE_SUFFIXES,
    PfamHmmBundleBuildResult,
    PfamHmmReferenceValidation,
    PfamHmmRegistryEntry,
    build_pfam_hmm_bundle,
    build_pfam_hmm_bundle_derivation_identity_sha256,
    load_pfam_hmm_bundle_registry,
    validate_pfam_hmm_bundle_file,
    validate_pfam_hmm_reference_data,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path

ARTIFACT_TYPE = "pfam_hmm_bundle_snapshot"
SNAPSHOT_FORMAT_VERSION = "1.0"
DEFAULT_MANIFEST_FILE_NAME = "manifest.json"
_SNAPSHOT_DIRECTORY_QUERY_LITERAL = "pfam_hmm_bundle"

_OUTPUT_KEY_BY_SUFFIX = {
    ".h3f": "pressed_h3f",
    ".h3i": "pressed_h3i",
    ".h3m": "pressed_h3m",
    ".h3p": "pressed_h3p",
}
_REQUIRED_OUTPUT_KEYS = frozenset({"bundle_hmm", *_OUTPUT_KEY_BY_SUFFIX.values()})
_EXPECTED_OUTPUT_FILE_NAME_BY_KEY = {
    "bundle_hmm": DEFAULT_BUNDLE_FILE_NAME,
    **{
        key: DEFAULT_BUNDLE_FILE_NAME + suffix
        for suffix, key in _OUTPUT_KEY_BY_SUFFIX.items()
    },
}
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_PFAM_ACCESSION_PATTERN = re.compile(r"^PF[0-9]{5}(?:\.[0-9]+)?$")


@dataclass(frozen=True)
class PfamHmmBundleSnapshotResult:
    snapshot_directory: Path
    snapshot_root_directory: Path
    manifest_file_path: Path
    bundle_file_path: Path
    pressed_file_paths: tuple[Path, ...]
    source_registry_sha256: str
    source_release: str
    model_count: int


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _current_utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _source_model_manifest_entries(
    *,
    validation: PfamHmmReferenceValidation,
) -> list[dict[str, object]]:
    return [
        {
            "canonical_order": model.registry_entry.canonical_order,
            "file_name": model.registry_entry.file_name,
            "accession": model.registry_entry.accession,
            "name": model.registry_entry.name,
            "alphabet": model.registry_entry.alphabet,
            "model_count": model.registry_entry.model_count,
            "sha256": model.registry_entry.sha256,
            "gathering_cutoffs_required": (
                model.registry_entry.gathering_cutoffs_required
            ),
        }
        for model in validation.validated_models
    ]


def _output_file_path_by_key(
    *,
    build_result: PfamHmmBundleBuildResult,
) -> dict[str, Path]:
    if len(build_result.pressed_file_paths) != len(PRESSED_FILE_SUFFIXES):
        raise RuntimeError(
            "A production Pfam HMM bundle snapshot requires all four hmmpress "
            "output files."
        )

    output_file_paths: dict[str, Path] = {
        "bundle_hmm": build_result.bundle_file_path,
    }
    for pressed_file_path in build_result.pressed_file_paths:
        matching_suffixes = [
            suffix
            for suffix in PRESSED_FILE_SUFFIXES
            if str(pressed_file_path).endswith(suffix)
        ]
        if len(matching_suffixes) != 1:
            raise RuntimeError(
                f"Unexpected hmmpress output file name: {pressed_file_path.name}."
            )
        suffix = matching_suffixes[0]
        output_file_paths[_OUTPUT_KEY_BY_SUFFIX[suffix]] = pressed_file_path

    if set(output_file_paths) != _REQUIRED_OUTPUT_KEYS:
        raise RuntimeError("The Pfam HMM bundle pressed output set is incomplete.")
    return output_file_paths


def _build_manifest_payload(
    *,
    snapshot_created_at_utc: str,
    immutable_snapshot_directory_name: str,
    immutable_snapshot_relative_path: str,
    build_result: PfamHmmBundleBuildResult,
    output_file_path_by_key: Mapping[str, Path],
) -> dict[str, object]:
    validation = build_result.validation
    source_models = _source_model_manifest_entries(validation=validation)
    output_files = {
        key: {
            "file_name": file_path.name,
            "path": str(file_path),
            "sha256": sha256_of_file(input_file_path=file_path),
            "size_bytes": file_path.stat().st_size,
        }
        for key, file_path in sorted(output_file_path_by_key.items())
    }
    return {
        "artifact_type": ARTIFACT_TYPE,
        "snapshot_format_version": SNAPSHOT_FORMAT_VERSION,
        "snapshot_created_at_utc": snapshot_created_at_utc,
        "manifest_file_name": DEFAULT_MANIFEST_FILE_NAME,
        "immutable_snapshot_directory_name": immutable_snapshot_directory_name,
        "immutable_snapshot_relative_path": immutable_snapshot_relative_path,
        "source_registry_file_name": validation.registry.registry_file_path.name,
        "source_registry_path": str(validation.registry.registry_file_path),
        "source_registry_sha256": validation.registry.registry_sha256,
        "source_database": validation.registry.source_database,
        "source_release": validation.registry.source_release,
        "source_model_count": validation.registry.model_count,
        "source_model_files": source_models,
        "canonical_model_file_names": [
            model.registry_entry.file_name for model in validation.validated_models
        ],
        "canonical_model_accessions": [
            model.registry_entry.accession for model in validation.validated_models
        ],
        "canonical_model_names": [
            model.registry_entry.name for model in validation.validated_models
        ],
        "canonical_model_alphabet": "amino",
        "bundle_model_count": validation.registry.model_count,
        "bundle_file_name": build_result.bundle_file_path.name,
        "bundle_file_sha256": build_result.bundle_sha256,
        "pressed_output_count": len(build_result.pressed_file_paths),
        "gathering_cutoffs_policy": "required_for_all_models",
        "concatenation_policy_version": CONCATENATION_POLICY_VERSION,
        "pyhmmer_version": validation.pyhmmer_version,
        "derivation_identity_sha256": (
            build_pfam_hmm_bundle_derivation_identity_sha256(
                validation=validation,
                pressed_outputs_required=True,
            )
        ),
        "output_files": output_files,
    }


def _manifest_models(
    *,
    manifest_payload: Mapping[str, object],
) -> tuple[PfamHmmRegistryEntry, ...]:
    source_model_files = manifest_payload.get("source_model_files")
    if not isinstance(source_model_files, Sequence) or isinstance(
        source_model_files, (str, bytes)
    ):
        raise RuntimeError(
            "Saved Pfam HMM bundle manifest must define source_model_files as an "
            "array."
        )

    models: list[PfamHmmRegistryEntry] = []
    for index, raw_entry in enumerate(source_model_files):
        if not isinstance(raw_entry, Mapping):
            raise RuntimeError(
                f"Saved Pfam HMM source_model_files[{index}] must be an object."
            )
        try:
            canonical_order = raw_entry["canonical_order"]
            file_name = raw_entry["file_name"]
            accession = raw_entry["accession"]
            name = raw_entry["name"]
            alphabet = raw_entry["alphabet"]
            model_count = raw_entry["model_count"]
            expected_sha256 = raw_entry["sha256"]
            gathering_cutoffs_required = raw_entry["gathering_cutoffs_required"]
        except KeyError as error:
            raise RuntimeError(
                f"Saved Pfam HMM source_model_files[{index}] is incomplete."
            ) from error

        if type(canonical_order) is not int or canonical_order != index + 1:
            raise RuntimeError(
                "Saved Pfam HMM source_model_files are not in contiguous canonical "
                "order."
            )
        if not all(
            isinstance(value, str) and bool(value)
            for value in (file_name, accession, name, alphabet, expected_sha256)
        ):
            raise RuntimeError(
                f"Saved Pfam HMM source_model_files[{index}] contains invalid text "
                "fields."
            )
        if Path(file_name).name != file_name or Path(file_name).suffix != ".hmm":
            raise RuntimeError(
                f"Saved Pfam HMM source_model_files[{index}] has an invalid "
                "file_name."
            )
        if not _PFAM_ACCESSION_PATTERN.fullmatch(accession):
            raise RuntimeError(
                f"Saved Pfam HMM source_model_files[{index}] has an invalid "
                "accession."
            )
        if alphabet != "amino":
            raise RuntimeError(
                f"Saved Pfam HMM source_model_files[{index}] has a non-amino "
                "alphabet."
            )
        if not _SHA256_PATTERN.fullmatch(expected_sha256):
            raise RuntimeError(
                f"Saved Pfam HMM source_model_files[{index}] has an invalid "
                "SHA-256."
            )
        if (
            type(model_count) is not int
            or model_count != 1
            or gathering_cutoffs_required is not True
        ):
            raise RuntimeError(
                f"Saved Pfam HMM source_model_files[{index}] violates the one-model "
                "or gathering-cutoff invariant."
            )
        models.append(
            PfamHmmRegistryEntry(
                canonical_order=canonical_order,
                file_name=file_name,
                accession=accession,
                name=name,
                alphabet=alphabet,
                model_count=model_count,
                sha256=expected_sha256,
                gathering_cutoffs_required=True,
            )
        )

    if not models:
        raise RuntimeError(
            "Saved Pfam HMM bundle source_model_files must contain at least one "
            "model."
        )

    expected_model_count = manifest_payload.get("source_model_count")
    if type(expected_model_count) is not int or expected_model_count != len(models):
        raise RuntimeError(
            "Saved Pfam HMM bundle source_model_count does not match "
            "source_model_files."
        )
    if manifest_payload.get("bundle_model_count") != len(models):
        raise RuntimeError(
            "Saved Pfam HMM bundle_model_count does not match source_model_files."
        )

    for field_name, values in {
        "file_name": [model.file_name for model in models],
        "accession": [model.accession for model in models],
        "name": [model.name for model in models],
    }.items():
        if len(values) != len(set(values)):
            raise RuntimeError(
                f"Saved Pfam HMM source model {field_name} values must be unique."
            )

    expected_accessions = [model.accession for model in models]
    expected_names = [model.name for model in models]
    expected_file_names = [model.file_name for model in models]
    for manifest_field, expected_value in (
        ("canonical_model_accessions", expected_accessions),
        ("canonical_model_names", expected_names),
        ("canonical_model_file_names", expected_file_names),
    ):
        if manifest_payload.get(manifest_field) != expected_value:
            raise RuntimeError(
                f"Saved Pfam HMM bundle {manifest_field} does not match the source "
                "model registry entries."
            )
    return tuple(models)


def _derivation_identity_from_manifest(
    *,
    manifest_payload: Mapping[str, object],
    expected_models: Sequence[PfamHmmRegistryEntry],
) -> str:
    payload = {
        "source_registry_sha256": manifest_payload.get(
            "source_registry_sha256"
        ),
        "source_database": manifest_payload.get("source_database"),
        "source_release": manifest_payload.get("source_release"),
        "canonical_models": [
            {
                "canonical_order": model.canonical_order,
                "file_name": model.file_name,
                "accession": model.accession,
                "name": model.name,
                "alphabet": model.alphabet,
                "model_count": model.model_count,
                "sha256": model.sha256,
                "gathering_cutoffs_required": (
                    model.gathering_cutoffs_required
                ),
            }
            for model in expected_models
        ],
        "concatenation_policy_version": manifest_payload.get(
            "concatenation_policy_version"
        ),
        "pyhmmer_version": manifest_payload.get("pyhmmer_version"),
        "pressed_outputs_required": True,
    }
    canonical_json = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical_json.encode("utf-8")).hexdigest()


def _validate_loaded_pfam_hmm_bundle_payload(
    *,
    snapshot_directory: Path,
    manifest_payload: Mapping[str, object],
    pyhmmer_module: object | None = None,
) -> dict[str, Path]:
    if manifest_payload.get("artifact_type") != ARTIFACT_TYPE:
        raise RuntimeError(
            "Saved Pfam HMM bundle snapshot manifest artifact_type mismatch."
        )
    if manifest_payload.get("snapshot_format_version") != SNAPSHOT_FORMAT_VERSION:
        raise RuntimeError(
            "Saved Pfam HMM bundle snapshot manifest snapshot_format_version "
            "mismatch."
        )
    if manifest_payload.get("concatenation_policy_version") != (
        CONCATENATION_POLICY_VERSION
    ):
        raise RuntimeError(
            "Saved Pfam HMM bundle concatenation policy version mismatch."
        )
    if manifest_payload.get("gathering_cutoffs_policy") != (
        "required_for_all_models"
    ):
        raise RuntimeError("Saved Pfam HMM bundle gathering-cutoff policy mismatch.")
    if manifest_payload.get("pressed_output_count") != len(PRESSED_FILE_SUFFIXES):
        raise RuntimeError("Saved Pfam HMM bundle pressed output count mismatch.")
    if manifest_payload.get("source_database") != "Pfam":
        raise RuntimeError("Saved Pfam HMM bundle source_database mismatch.")
    if manifest_payload.get("manifest_file_name") != DEFAULT_MANIFEST_FILE_NAME:
        raise RuntimeError("Saved Pfam HMM bundle manifest_file_name mismatch.")
    source_registry_file_name = manifest_payload.get("source_registry_file_name")
    if (
        not isinstance(source_registry_file_name, str)
        or Path(source_registry_file_name).name != source_registry_file_name
        or not source_registry_file_name.endswith(".json")
    ):
        raise RuntimeError(
            "Saved Pfam HMM bundle source_registry_file_name is invalid."
        )
    if manifest_payload.get("canonical_model_alphabet") != "amino":
        raise RuntimeError(
            "Saved Pfam HMM bundle canonical_model_alphabet mismatch."
        )
    source_release = manifest_payload.get("source_release")
    if not isinstance(source_release, str) or not source_release:
        raise RuntimeError("Saved Pfam HMM bundle source_release is missing.")
    source_registry_sha256 = manifest_payload.get("source_registry_sha256")
    if not isinstance(source_registry_sha256, str) or not _SHA256_PATTERN.fullmatch(
        source_registry_sha256
    ):
        raise RuntimeError(
            "Saved Pfam HMM bundle source_registry_sha256 is invalid."
        )
    derivation_identity_sha256 = manifest_payload.get("derivation_identity_sha256")
    if not isinstance(
        derivation_identity_sha256, str
    ) or not _SHA256_PATTERN.fullmatch(derivation_identity_sha256):
        raise RuntimeError(
            "Saved Pfam HMM bundle derivation_identity_sha256 is invalid."
        )
    pyhmmer_version = manifest_payload.get("pyhmmer_version")
    if not isinstance(pyhmmer_version, str) or not pyhmmer_version.strip():
        raise RuntimeError("Saved Pfam HMM bundle pyhmmer_version is missing.")

    expected_models = _manifest_models(manifest_payload=manifest_payload)
    expected_derivation_identity_sha256 = _derivation_identity_from_manifest(
        manifest_payload=manifest_payload,
        expected_models=expected_models,
    )
    if derivation_identity_sha256 != expected_derivation_identity_sha256:
        raise RuntimeError(
            "Saved Pfam HMM bundle derivation identity is inconsistent with "
            "its source and build-policy metadata."
        )
    output_files = manifest_payload.get("output_files")
    if not isinstance(output_files, Mapping):
        raise RuntimeError(
            "Saved Pfam HMM bundle manifest must define output_files."
        )
    if set(output_files) != _REQUIRED_OUTPUT_KEYS:
        raise RuntimeError(
            "Saved Pfam HMM bundle manifest output_files has an incomplete or "
            "unexpected key set."
        )

    resolved_output_file_path_by_key: dict[str, Path] = {}
    for key, raw_entry in output_files.items():
        if not isinstance(raw_entry, Mapping):
            raise RuntimeError(
                f"Saved Pfam HMM output_files entry {key!r} must be an object."
            )
        file_name = raw_entry.get("file_name")
        expected_sha256 = raw_entry.get("sha256")
        expected_size_bytes = raw_entry.get("size_bytes")
        if not isinstance(file_name, str) or Path(file_name).name != file_name:
            raise RuntimeError(
                f"Saved Pfam HMM output_files entry {key!r} has an invalid file_name."
            )
        if file_name != _EXPECTED_OUTPUT_FILE_NAME_BY_KEY[str(key)]:
            raise RuntimeError(
                f"Saved Pfam HMM output_files entry {key!r} has an unexpected "
                "file_name."
            )
        if not isinstance(expected_sha256, str) or not _SHA256_PATTERN.fullmatch(
            expected_sha256
        ):
            raise RuntimeError(
                f"Saved Pfam HMM output_files entry {key!r} has an invalid SHA-256."
            )
        if type(expected_size_bytes) is not int or expected_size_bytes < 1:
            raise RuntimeError(
                f"Saved Pfam HMM output_files entry {key!r} has an invalid size."
            )

        output_file_path = snapshot_directory / file_name
        if not output_file_path.is_file():
            raise FileNotFoundError(
                f"Saved Pfam HMM bundle file was not found: {output_file_path}."
            )
        if output_file_path.stat().st_size != expected_size_bytes:
            raise RuntimeError(
                f"Saved Pfam HMM bundle size mismatch for {key!r}."
            )
        actual_sha256 = sha256_of_file(input_file_path=output_file_path)
        if actual_sha256 != expected_sha256:
            raise RuntimeError(
                f"Saved Pfam HMM bundle hash mismatch for {key!r}. Expected "
                f"{expected_sha256}, got {actual_sha256}."
            )
        resolved_output_file_path_by_key[str(key)] = output_file_path

    bundle_file_path = resolved_output_file_path_by_key["bundle_hmm"]
    if manifest_payload.get("bundle_file_name") != bundle_file_path.name:
        raise RuntimeError("Saved Pfam HMM bundle_file_name mismatch.")
    if manifest_payload.get("bundle_file_sha256") != sha256_of_file(
        input_file_path=bundle_file_path
    ):
        raise RuntimeError("Saved Pfam HMM bundle_file_sha256 mismatch.")

    validate_pfam_hmm_bundle_file(
        bundle_file_path=bundle_file_path,
        expected_models=expected_models,
        pyhmmer_module=pyhmmer_module,
    )
    resolved_pyhmmer_module = (
        pyhmmer_module
        if pyhmmer_module is not None
        else importlib.import_module("pyhmmer")
    )
    current_pyhmmer_version = str(
        getattr(resolved_pyhmmer_module, "__version__", "unknown")
    )
    if current_pyhmmer_version != pyhmmer_version:
        raise RuntimeError(
            "Saved Pfam HMM bundle PyHMMER version mismatch. Expected "
            f"{pyhmmer_version!r}, got {current_pyhmmer_version!r}."
        )
    return resolved_output_file_path_by_key


def save_pfam_hmm_bundle_snapshot(
    *,
    snapshot_root_directory: PathLike,
    registry_file_path: PathLike = DEFAULT_REGISTRY_FILE_PATH,
    hmm_directory: PathLike | None = None,
    pyhmmer_module: object | None = None,
    update_latest_directory: bool = True,
) -> PfamHmmBundleSnapshotResult:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    registry = load_pfam_hmm_bundle_registry(
        registry_file_path=registry_file_path
    )
    snapshot_created_at_utc = _current_utc_timestamp()
    snapshot_directory_name = build_snapshot_directory_name(
        retrieved_at_utc=snapshot_created_at_utc,
        search_query=(
            f"{_SNAPSHOT_DIRECTORY_QUERY_LITERAL}:{registry.registry_sha256}"
        ),
    )
    immutable_snapshot_directory = (
        resolved_snapshot_root_directory / "snapshots" / snapshot_directory_name
    )
    immutable_snapshot_directory.mkdir(parents=True, exist_ok=False)
    immutable_snapshot_relative_path = str(Path("snapshots") / snapshot_directory_name)
    immutable_snapshot_complete = False

    try:
        build_result = build_pfam_hmm_bundle(
            output_directory=immutable_snapshot_directory,
            registry_file_path=registry_file_path,
            hmm_directory=hmm_directory,
            bundle_file_name=DEFAULT_BUNDLE_FILE_NAME,
            pyhmmer_module=pyhmmer_module,
            press_bundle=True,
        )
        output_file_path_by_key = _output_file_path_by_key(build_result=build_result)
        manifest_file_path = immutable_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
        manifest_payload = _build_manifest_payload(
            snapshot_created_at_utc=snapshot_created_at_utc,
            immutable_snapshot_directory_name=snapshot_directory_name,
            immutable_snapshot_relative_path=immutable_snapshot_relative_path,
            build_result=build_result,
            output_file_path_by_key=output_file_path_by_key,
        )
        write_json_atomic(
            payload=manifest_payload,
            output_file_path=manifest_file_path,
        )
        _validate_loaded_pfam_hmm_bundle_payload(
            snapshot_directory=immutable_snapshot_directory,
            manifest_payload=manifest_payload,
            pyhmmer_module=pyhmmer_module,
        )
        immutable_snapshot_complete = True

        if update_latest_directory:
            files_to_copy = [
                (file_path, file_path.name)
                for file_path in output_file_path_by_key.values()
            ]
            files_to_copy.append((manifest_file_path, DEFAULT_MANIFEST_FILE_NAME))
            _replace_latest_directory(
                latest_directory=resolved_snapshot_root_directory / "latest",
                files_to_copy=files_to_copy,
            )
    except Exception:
        if not immutable_snapshot_complete and immutable_snapshot_directory.exists():
            shutil.rmtree(immutable_snapshot_directory, ignore_errors=True)
        raise

    return PfamHmmBundleSnapshotResult(
        snapshot_directory=immutable_snapshot_directory,
        snapshot_root_directory=resolved_snapshot_root_directory,
        manifest_file_path=manifest_file_path,
        bundle_file_path=build_result.bundle_file_path,
        pressed_file_paths=build_result.pressed_file_paths,
        source_registry_sha256=build_result.validation.registry.registry_sha256,
        source_release=build_result.validation.registry.source_release,
        model_count=build_result.validation.registry.model_count,
    )


def load_pfam_hmm_bundle_snapshot_by_directory(
    *,
    snapshot_directory: PathLike,
    pyhmmer_module: object | None = None,
) -> dict[str, object]:
    resolved_snapshot_directory = _as_path(snapshot_directory)
    manifest_file_path = resolved_snapshot_directory / DEFAULT_MANIFEST_FILE_NAME
    manifest_payload = read_json_file(input_file_path=manifest_file_path)
    if not isinstance(manifest_payload, Mapping):
        raise RuntimeError(
            "Saved Pfam HMM bundle manifest must deserialize into a JSON object."
        )
    output_file_paths = _validate_loaded_pfam_hmm_bundle_payload(
        snapshot_directory=resolved_snapshot_directory,
        manifest_payload=manifest_payload,
        pyhmmer_module=pyhmmer_module,
    )
    pressed_file_paths = tuple(
        output_file_paths[_OUTPUT_KEY_BY_SUFFIX[suffix]]
        for suffix in PRESSED_FILE_SUFFIXES
    )
    return {
        "snapshot_directory": resolved_snapshot_directory,
        "manifest_file_path": manifest_file_path,
        "manifest": dict(manifest_payload),
        "bundle_file_path": output_file_paths["bundle_hmm"],
        "pressed_file_paths": pressed_file_paths,
        "pressed_file_path_by_suffix": {
            suffix: output_file_paths[_OUTPUT_KEY_BY_SUFFIX[suffix]]
            for suffix in PRESSED_FILE_SUFFIXES
        },
    }


def load_latest_pfam_hmm_bundle_snapshot(
    *,
    snapshot_root_directory: PathLike,
    pyhmmer_module: object | None = None,
) -> dict[str, object]:
    return load_pfam_hmm_bundle_snapshot_by_directory(
        snapshot_directory=_as_path(snapshot_root_directory) / "latest",
        pyhmmer_module=pyhmmer_module,
    )


def _source_reference_identity_matches(
    *,
    manifest_payload: Mapping[str, object],
    registry_file_path: PathLike,
    hmm_directory: PathLike | None,
    pyhmmer_module: object | None,
) -> bool:
    current_validation = validate_pfam_hmm_reference_data(
        registry_file_path=registry_file_path,
        hmm_directory=hmm_directory,
        pyhmmer_module=pyhmmer_module,
    )
    expected_derivation_identity = (
        build_pfam_hmm_bundle_derivation_identity_sha256(
            validation=current_validation,
            pressed_outputs_required=True,
        )
    )
    return (
        manifest_payload.get("source_registry_sha256")
        == current_validation.registry.registry_sha256
        and manifest_payload.get("source_release")
        == current_validation.registry.source_release
        and manifest_payload.get("derivation_identity_sha256")
        == expected_derivation_identity
        and manifest_payload.get("source_model_files")
        == _source_model_manifest_entries(validation=current_validation)
    )


def latest_pfam_hmm_bundle_snapshot_is_available(
    *,
    snapshot_root_directory: PathLike,
    registry_file_path: PathLike = DEFAULT_REGISTRY_FILE_PATH,
    hmm_directory: PathLike | None = None,
    pyhmmer_module: object | None = None,
) -> bool:
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_directory = resolved_snapshot_root_directory / "latest"
    latest_manifest_file_path = latest_directory / DEFAULT_MANIFEST_FILE_NAME
    if not latest_directory.is_dir() or not latest_manifest_file_path.is_file():
        return False

    try:
        payload = load_pfam_hmm_bundle_snapshot_by_directory(
            snapshot_directory=latest_directory,
            pyhmmer_module=pyhmmer_module,
        )
        manifest_payload = payload["manifest"]
        if not isinstance(manifest_payload, Mapping):
            return False
        if not _source_reference_identity_matches(
            manifest_payload=manifest_payload,
            registry_file_path=registry_file_path,
            hmm_directory=hmm_directory,
            pyhmmer_module=pyhmmer_module,
        ):
            return False
    except (FileNotFoundError, RuntimeError, OSError, TypeError, ValueError):
        return False
    return True


def resolve_pfam_hmm_bundle_snapshot(
    *,
    snapshot_mode: SnapshotMode | str,
    snapshot_root_directory: PathLike,
    registry_file_path: PathLike = DEFAULT_REGISTRY_FILE_PATH,
    hmm_directory: PathLike | None = None,
    pyhmmer_module: object | None = None,
    update_latest_directory: bool = True,
) -> dict[str, object]:
    resolved_snapshot_mode = _coerce_snapshot_mode(snapshot_mode)
    resolved_snapshot_root_directory = _as_path(snapshot_root_directory)
    latest_is_available = latest_pfam_hmm_bundle_snapshot_is_available(
        snapshot_root_directory=resolved_snapshot_root_directory,
        registry_file_path=registry_file_path,
        hmm_directory=hmm_directory,
        pyhmmer_module=pyhmmer_module,
    )

    if resolved_snapshot_mode == SnapshotMode.reuse_latest:
        if not latest_is_available:
            raise FileNotFoundError(
                "No latest Pfam HMM bundle snapshot matches the committed reference "
                "registry and local models."
            )
        return load_latest_pfam_hmm_bundle_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory,
            pyhmmer_module=pyhmmer_module,
        )

    if (
        resolved_snapshot_mode == SnapshotMode.reuse_latest_or_create
        and latest_is_available
    ):
        print("Latest Pfam HMM bundle snapshot is available. Reusing frozen snapshot.")
        return load_latest_pfam_hmm_bundle_snapshot(
            snapshot_root_directory=resolved_snapshot_root_directory,
            pyhmmer_module=pyhmmer_module,
        )

    if resolved_snapshot_mode not in {
        SnapshotMode.create_new,
        SnapshotMode.reuse_latest_or_create,
    }:
        raise ValueError(
            "Invalid snapshot_mode. Expected one of: 'create_new', "
            "'reuse_latest', 'reuse_latest_or_create'."
        )

    saved_result = save_pfam_hmm_bundle_snapshot(
        snapshot_root_directory=resolved_snapshot_root_directory,
        registry_file_path=registry_file_path,
        hmm_directory=hmm_directory,
        pyhmmer_module=pyhmmer_module,
        update_latest_directory=update_latest_directory,
    )
    return load_pfam_hmm_bundle_snapshot_by_directory(
        snapshot_directory=saved_result.snapshot_directory,
        pyhmmer_module=pyhmmer_module,
    )


def list_saved_pfam_hmm_bundle_snapshot_directories(
    *,
    snapshot_root_directory: PathLike,
) -> list[Path]:
    return list_saved_snapshot_directories(
        snapshot_root_directory=snapshot_root_directory
    )


def get_most_recent_pfam_hmm_bundle_snapshot_directory(
    *,
    snapshot_root_directory: PathLike,
) -> Optional[Path]:
    return get_most_recent_snapshot_directory(
        snapshot_root_directory=snapshot_root_directory
    )
