from __future__ import annotations

import hashlib
import importlib
import json
import re
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TypeAlias

from src.pago_pipeline.storage import sha256_of_file, write_bytes_atomic

PathLike: TypeAlias = str | Path

REGISTRY_ARTIFACT_TYPE = "pfam_hmm_bundle_lock"
REGISTRY_FORMAT_VERSION = "1.0"
CONCATENATION_POLICY_VERSION = "1.0"
DEFAULT_BUNDLE_FILE_NAME = "pago_domains.hmm"
PRESSED_FILE_SUFFIXES: tuple[str, ...] = (".h3f", ".h3i", ".h3m", ".h3p")

DEFAULT_REFERENCE_DIRECTORY = Path(__file__).resolve().parent / "resources" / "pfam_hmm"
DEFAULT_REGISTRY_FILE_PATH = (
    DEFAULT_REFERENCE_DIRECTORY / "pfam_hmm_bundle_lock.json"
)

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_PFAM_ACCESSION_PATTERN = re.compile(r"^PF[0-9]{5}\.[0-9]+$")


@dataclass(frozen=True)
class PfamHmmRegistryEntry:
    canonical_order: int
    file_name: str
    accession: str
    name: str
    alphabet: str
    model_count: int
    sha256: str
    gathering_cutoffs_required: bool


@dataclass(frozen=True)
class PfamHmmBundleRegistry:
    registry_file_path: Path
    registry_sha256: str
    source_database: str
    source_release: str
    model_count: int
    models: tuple[PfamHmmRegistryEntry, ...]


@dataclass(frozen=True)
class InspectedHmmModel:
    accession: str
    name: str
    alphabet: str
    has_gathering_cutoffs: bool


@dataclass(frozen=True)
class ValidatedPfamHmmModel:
    registry_entry: PfamHmmRegistryEntry
    source_file_path: Path
    inspected_model: InspectedHmmModel


@dataclass(frozen=True)
class PfamHmmReferenceValidation:
    registry: PfamHmmBundleRegistry
    hmm_directory: Path
    validated_models: tuple[ValidatedPfamHmmModel, ...]
    pyhmmer_version: str


@dataclass(frozen=True)
class PfamHmmBundleBuildResult:
    validation: PfamHmmReferenceValidation
    bundle_file_path: Path
    pressed_file_paths: tuple[Path, ...]
    bundle_sha256: str


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _require_nonempty_string(*, value: object, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise RuntimeError(f"Pfam HMM registry field {field_name!r} must be non-empty.")
    if value != value.strip():
        raise RuntimeError(
            f"Pfam HMM registry field {field_name!r} must not contain surrounding "
            "whitespace."
        )
    return value


def _require_exact_integer(*, value: object, field_name: str) -> int:
    if type(value) is not int:
        raise RuntimeError(f"Pfam HMM registry field {field_name!r} must be an integer.")
    return value


def _validate_registry_entry(
    *,
    entry_payload: object,
    entry_index: int,
) -> PfamHmmRegistryEntry:
    if not isinstance(entry_payload, Mapping):
        raise RuntimeError(
            f"Pfam HMM registry models[{entry_index}] must be a JSON object."
        )

    canonical_order = _require_exact_integer(
        value=entry_payload.get("canonical_order"),
        field_name=f"models[{entry_index}].canonical_order",
    )
    file_name = _require_nonempty_string(
        value=entry_payload.get("file_name"),
        field_name=f"models[{entry_index}].file_name",
    )
    accession = _require_nonempty_string(
        value=entry_payload.get("accession"),
        field_name=f"models[{entry_index}].accession",
    )
    name = _require_nonempty_string(
        value=entry_payload.get("name"),
        field_name=f"models[{entry_index}].name",
    )
    alphabet = _require_nonempty_string(
        value=entry_payload.get("alphabet"),
        field_name=f"models[{entry_index}].alphabet",
    )
    model_count = _require_exact_integer(
        value=entry_payload.get("model_count"),
        field_name=f"models[{entry_index}].model_count",
    )
    expected_sha256 = _require_nonempty_string(
        value=entry_payload.get("sha256"),
        field_name=f"models[{entry_index}].sha256",
    )
    gathering_cutoffs_required = entry_payload.get("gathering_cutoffs_required")

    if canonical_order < 1:
        raise RuntimeError("Pfam HMM canonical_order values must start at one.")
    if Path(file_name).name != file_name or Path(file_name).suffix != ".hmm":
        raise RuntimeError(
            "Every Pfam HMM registry file_name must be a local leaf path ending "
            "in '.hmm'."
        )
    if not _PFAM_ACCESSION_PATTERN.fullmatch(accession):
        raise RuntimeError(
            f"Pfam HMM accession {accession!r} must be a version-pinned Pfam "
            "accession in the PFxxxxx.version format (for example 'PF02171.24'); "
            "an unversioned accession is not accepted."
        )
    if alphabet != "amino":
        raise RuntimeError(
            f"Pfam HMM registry alphabet must be exactly 'amino', got {alphabet!r}."
        )
    if model_count != 1:
        raise RuntimeError(
            "Each committed Pfam HMM source file must contain exactly one model."
        )
    if not _SHA256_PATTERN.fullmatch(expected_sha256):
        raise RuntimeError(
            f"Pfam HMM registry SHA-256 for {file_name!r} must be 64 lowercase "
            "hexadecimal characters."
        )
    if gathering_cutoffs_required is not True:
        raise RuntimeError(
            f"Pfam HMM registry entry {accession!r} must require gathering cutoffs."
        )

    return PfamHmmRegistryEntry(
        canonical_order=canonical_order,
        file_name=file_name,
        accession=accession,
        name=name,
        alphabet=alphabet,
        model_count=model_count,
        sha256=expected_sha256,
        gathering_cutoffs_required=True,
    )


def load_pfam_hmm_bundle_registry(
    *,
    registry_file_path: PathLike = DEFAULT_REGISTRY_FILE_PATH,
) -> PfamHmmBundleRegistry:
    """Load and structurally validate the committed Pfam bundle lock file."""
    resolved_registry_file_path = _as_path(registry_file_path)
    try:
        with resolved_registry_file_path.open("r", encoding="utf-8") as file_handle:
            registry_payload = json.load(file_handle)
    except json.JSONDecodeError as error:
        raise RuntimeError(
            f"Pfam HMM registry is not valid JSON: {resolved_registry_file_path}."
        ) from error

    if not isinstance(registry_payload, Mapping):
        raise RuntimeError("Pfam HMM registry must deserialize into a JSON object.")
    if registry_payload.get("artifact_type") != REGISTRY_ARTIFACT_TYPE:
        raise RuntimeError(
            "Pfam HMM registry artifact_type mismatch. Expected "
            f"{REGISTRY_ARTIFACT_TYPE!r}."
        )
    if registry_payload.get("lock_format_version") != REGISTRY_FORMAT_VERSION:
        raise RuntimeError(
            "Pfam HMM registry lock_format_version mismatch. Expected "
            f"{REGISTRY_FORMAT_VERSION!r}."
        )

    source_database = _require_nonempty_string(
        value=registry_payload.get("source_database"),
        field_name="source_database",
    )
    source_release = _require_nonempty_string(
        value=registry_payload.get("source_release"),
        field_name="source_release",
    )
    model_count = _require_exact_integer(
        value=registry_payload.get("model_count"),
        field_name="model_count",
    )
    if source_database != "Pfam":
        raise RuntimeError(
            f"Pfam HMM registry source_database must be 'Pfam', got {source_database!r}."
        )
    if model_count < 1:
        raise RuntimeError("Pfam HMM registry must declare at least one model.")

    model_payloads = registry_payload.get("models")
    if not isinstance(model_payloads, Sequence) or isinstance(
        model_payloads, (str, bytes)
    ):
        raise RuntimeError("Pfam HMM registry models must be a JSON array.")
    models = tuple(
        _validate_registry_entry(entry_payload=payload, entry_index=index)
        for index, payload in enumerate(model_payloads)
    )
    if len(models) != model_count:
        raise RuntimeError(
            "Pfam HMM registry model_count does not equal the number of model "
            f"entries: declared {model_count}, found {len(models)}."
        )

    observed_orders = [entry.canonical_order for entry in models]
    expected_orders = list(range(1, model_count + 1))
    if observed_orders != expected_orders:
        raise RuntimeError(
            "Pfam HMM registry models must already be listed in contiguous canonical "
            f"order {expected_orders}; found {observed_orders}."
        )

    for field_name, values in {
        "file_name": [entry.file_name for entry in models],
        "accession": [entry.accession for entry in models],
        "name": [entry.name for entry in models],
    }.items():
        if len(values) != len(set(values)):
            raise RuntimeError(
                f"Pfam HMM registry model {field_name} values must be unique."
            )

    return PfamHmmBundleRegistry(
        registry_file_path=resolved_registry_file_path,
        registry_sha256=sha256_of_file(input_file_path=resolved_registry_file_path),
        source_database=source_database,
        source_release=source_release,
        model_count=model_count,
        models=models,
    )


def _load_pyhmmer_module(pyhmmer_module: object | None) -> object:
    if pyhmmer_module is not None:
        return pyhmmer_module
    try:
        return importlib.import_module("pyhmmer")
    except ModuleNotFoundError as error:
        raise RuntimeError(
            "PyHMMER is required to validate and press the committed Pfam HMM bundle."
        ) from error


def _decode_hmm_text(value: object, *, field_name: str) -> str:
    if isinstance(value, bytes):
        try:
            return value.decode("utf-8")
        except UnicodeDecodeError as error:
            raise RuntimeError(
                f"HMM {field_name} is not valid UTF-8 text."
            ) from error
    if isinstance(value, str):
        return value
    if value is None:
        return ""
    return str(value)


def _canonical_alphabet_name(alphabet: object) -> str:
    alphabet_name = getattr(alphabet, "name", None)
    if callable(alphabet_name):
        alphabet_name = alphabet_name()
    if alphabet_name is not None:
        decoded_name = _decode_hmm_text(alphabet_name, field_name="alphabet").lower()
        if decoded_name in {"amino", "protein", "amino acid", "amino acids"}:
            return "amino"
        if decoded_name in {"dna", "rna"}:
            return decoded_name

    for canonical_name, method_name in (
        ("amino", "is_amino"),
        ("dna", "is_dna"),
        ("rna", "is_rna"),
    ):
        predicate = getattr(alphabet, method_name, None)
        if callable(predicate) and bool(predicate()):
            return canonical_name

    return str(alphabet).strip().lower()


def _has_gathering_cutoffs(hmm_model: object) -> bool:
    cutoffs = getattr(hmm_model, "cutoffs", None)
    if cutoffs is None:
        return False
    gathering_available = getattr(cutoffs, "gathering_available", None)
    if callable(gathering_available):
        return bool(gathering_available())
    gathering = getattr(cutoffs, "gathering", None)
    if gathering is None:
        return False
    if isinstance(gathering, Sequence) and not isinstance(gathering, (str, bytes)):
        return len(gathering) > 0 and all(value is not None for value in gathering)
    return True


def _inspect_hmm_file(
    *,
    hmm_file_path: Path,
    pyhmmer_module: object,
) -> tuple[InspectedHmmModel, ...]:
    plan7 = getattr(pyhmmer_module, "plan7", None)
    hmm_file_factory = getattr(plan7, "HMMFile", None)
    if not callable(hmm_file_factory):
        raise RuntimeError("The provided PyHMMER module does not expose plan7.HMMFile.")

    try:
        with hmm_file_factory(str(hmm_file_path)) as hmm_file:
            hmm_models = list(hmm_file)
    except Exception as error:
        raise RuntimeError(f"Unable to parse HMM file {hmm_file_path}.") from error

    return tuple(
        InspectedHmmModel(
            accession=_decode_hmm_text(
                getattr(hmm_model, "accession", None),
                field_name="accession",
            ),
            name=_decode_hmm_text(
                getattr(hmm_model, "name", None),
                field_name="name",
            ),
            alphabet=_canonical_alphabet_name(getattr(hmm_model, "alphabet", None)),
            has_gathering_cutoffs=_has_gathering_cutoffs(hmm_model),
        )
        for hmm_model in hmm_models
    )


def _validate_exact_hmm_inventory(
    *,
    hmm_directory: Path,
    registry: PfamHmmBundleRegistry,
) -> None:
    expected_file_names = {entry.file_name for entry in registry.models}
    actual_file_names = {
        file_path.name
        for file_path in hmm_directory.iterdir()
        if file_path.is_file() and file_path.suffix == ".hmm"
    }
    missing_file_names = sorted(expected_file_names - actual_file_names)
    unexpected_file_names = sorted(actual_file_names - expected_file_names)
    if missing_file_names or unexpected_file_names:
        raise RuntimeError(
            "Pfam HMM directory inventory does not exactly match the registry. "
            f"Missing={missing_file_names}; unexpected={unexpected_file_names}."
        )


def validate_pfam_hmm_reference_data(
    *,
    registry_file_path: PathLike = DEFAULT_REGISTRY_FILE_PATH,
    hmm_directory: PathLike | None = None,
    pyhmmer_module: object | None = None,
) -> PfamHmmReferenceValidation:
    """Validate the committed lock and every local Pfam HMM without mutation."""
    registry = load_pfam_hmm_bundle_registry(registry_file_path=registry_file_path)
    resolved_hmm_directory = (
        _as_path(hmm_directory)
        if hmm_directory is not None
        else registry.registry_file_path.parent
    )
    if not resolved_hmm_directory.is_dir():
        raise FileNotFoundError(
            f"Pfam HMM reference directory was not found: {resolved_hmm_directory}."
        )

    _validate_exact_hmm_inventory(
        hmm_directory=resolved_hmm_directory,
        registry=registry,
    )
    resolved_pyhmmer_module = _load_pyhmmer_module(pyhmmer_module)
    validated_models: list[ValidatedPfamHmmModel] = []

    for registry_entry in registry.models:
        source_file_path = resolved_hmm_directory / registry_entry.file_name
        actual_sha256 = sha256_of_file(input_file_path=source_file_path)
        if actual_sha256 != registry_entry.sha256:
            raise RuntimeError(
                f"Pfam HMM SHA-256 mismatch for {registry_entry.file_name!r}. "
                f"Expected {registry_entry.sha256}, got {actual_sha256}."
            )

        inspected_models = _inspect_hmm_file(
            hmm_file_path=source_file_path,
            pyhmmer_module=resolved_pyhmmer_module,
        )
        if len(inspected_models) != registry_entry.model_count:
            raise RuntimeError(
                f"Pfam HMM model count mismatch for {registry_entry.file_name!r}. "
                f"Expected {registry_entry.model_count}, got {len(inspected_models)}."
            )
        inspected_model = inspected_models[0]
        for field_name, expected_value, actual_value in (
            ("accession", registry_entry.accession, inspected_model.accession),
            ("name", registry_entry.name, inspected_model.name),
            ("alphabet", registry_entry.alphabet, inspected_model.alphabet),
        ):
            if actual_value != expected_value:
                raise RuntimeError(
                    f"Pfam HMM {field_name} mismatch for "
                    f"{registry_entry.file_name!r}. Expected {expected_value!r}, "
                    f"got {actual_value!r}."
                )
        if (
            registry_entry.gathering_cutoffs_required
            and not inspected_model.has_gathering_cutoffs
        ):
            raise RuntimeError(
                f"Pfam HMM {registry_entry.accession!r} is missing required "
                "gathering cutoffs."
            )

        validated_models.append(
            ValidatedPfamHmmModel(
                registry_entry=registry_entry,
                source_file_path=source_file_path,
                inspected_model=inspected_model,
            )
        )

    return PfamHmmReferenceValidation(
        registry=registry,
        hmm_directory=resolved_hmm_directory,
        validated_models=tuple(validated_models),
        pyhmmer_version=str(getattr(resolved_pyhmmer_module, "__version__", "unknown")),
    )


def validate_pfam_hmm_bundle_file(
    *,
    bundle_file_path: PathLike,
    expected_models: Sequence[PfamHmmRegistryEntry],
    pyhmmer_module: object | None = None,
) -> tuple[InspectedHmmModel, ...]:
    """Validate the count, order, identity, alphabet, and GA data of a bundle."""
    resolved_pyhmmer_module = _load_pyhmmer_module(pyhmmer_module)
    inspected_models = _inspect_hmm_file(
        hmm_file_path=_as_path(bundle_file_path),
        pyhmmer_module=resolved_pyhmmer_module,
    )
    if len(inspected_models) != len(expected_models):
        raise RuntimeError(
            "Concatenated Pfam HMM bundle model count mismatch. Expected "
            f"{len(expected_models)}, got {len(inspected_models)}."
        )

    for model_index, (expected_model, inspected_model) in enumerate(
        zip(expected_models, inspected_models, strict=True),
        start=1,
    ):
        for field_name, expected_value, actual_value in (
            ("accession", expected_model.accession, inspected_model.accession),
            ("name", expected_model.name, inspected_model.name),
            ("alphabet", expected_model.alphabet, inspected_model.alphabet),
        ):
            if actual_value != expected_value:
                raise RuntimeError(
                    "Concatenated Pfam HMM bundle canonical-order mismatch at "
                    f"model {model_index} for {field_name}. Expected "
                    f"{expected_value!r}, got {actual_value!r}."
                )
        if (
            expected_model.gathering_cutoffs_required
            and not inspected_model.has_gathering_cutoffs
        ):
            raise RuntimeError(
                "Concatenated Pfam HMM bundle model "
                f"{expected_model.accession!r} is missing gathering cutoffs."
            )

    return inspected_models


def _build_concatenated_hmm_bytes(
    *,
    validated_models: Sequence[ValidatedPfamHmmModel],
) -> bytes:
    model_payloads: list[bytes] = []
    for validated_model in validated_models:
        source_payload = validated_model.source_file_path.read_bytes()
        if not source_payload:
            raise RuntimeError(
                f"Pfam HMM file is empty: {validated_model.source_file_path}."
            )
        model_payloads.append(source_payload.rstrip(b"\r\n") + b"\n")
    return b"".join(model_payloads)


def _iter_source_hmm_models(
    *,
    validation: PfamHmmReferenceValidation,
    pyhmmer_module: object,
):
    plan7 = getattr(pyhmmer_module, "plan7", None)
    hmm_file_factory = getattr(plan7, "HMMFile", None)
    if not callable(hmm_file_factory):
        raise RuntimeError("The provided PyHMMER module does not expose plan7.HMMFile.")
    for validated_model in validation.validated_models:
        with hmm_file_factory(str(validated_model.source_file_path)) as hmm_file:
            hmm_models = tuple(hmm_file)
        for hmm_model in hmm_models:
            yield hmm_model


def build_pfam_hmm_bundle(
    *,
    output_directory: PathLike,
    registry_file_path: PathLike = DEFAULT_REGISTRY_FILE_PATH,
    hmm_directory: PathLike | None = None,
    bundle_file_name: str = DEFAULT_BUNDLE_FILE_NAME,
    pyhmmer_module: object | None = None,
    press_bundle: bool = True,
) -> PfamHmmBundleBuildResult:
    """Validate, concatenate, and optionally press the canonical Pfam HMM bundle."""
    if Path(bundle_file_name).name != bundle_file_name or not bundle_file_name.endswith(
        ".hmm"
    ):
        raise ValueError("bundle_file_name must be a local leaf path ending in '.hmm'.")

    resolved_pyhmmer_module = _load_pyhmmer_module(pyhmmer_module)
    validation = validate_pfam_hmm_reference_data(
        registry_file_path=registry_file_path,
        hmm_directory=hmm_directory,
        pyhmmer_module=resolved_pyhmmer_module,
    )
    resolved_output_directory = _as_path(output_directory)
    resolved_output_directory.mkdir(parents=True, exist_ok=True)
    bundle_file_path = resolved_output_directory / bundle_file_name
    canonical_output_paths = (
        bundle_file_path,
        *(
            Path(str(bundle_file_path) + suffix)
            for suffix in PRESSED_FILE_SUFFIXES
        ),
    )
    preexisting_output_paths = [
        output_path for output_path in canonical_output_paths if output_path.exists()
    ]
    if preexisting_output_paths:
        raise FileExistsError(
            "Refusing to overwrite pre-existing Pfam HMM bundle outputs: "
            + ", ".join(str(path) for path in preexisting_output_paths)
        )

    staging_directory = Path(
        tempfile.mkdtemp(
            prefix=".pfam_hmm_bundle_staged_",
            dir=resolved_output_directory,
        )
    )
    staged_bundle_file_path = staging_directory / bundle_file_name
    staged_pressed_file_paths = tuple(
        Path(str(staged_bundle_file_path) + suffix)
        for suffix in PRESSED_FILE_SUFFIXES
    )
    pressed_file_paths: tuple[Path, ...] = tuple()

    try:
        write_bytes_atomic(
            binary_payload=_build_concatenated_hmm_bytes(
                validated_models=validation.validated_models
            ),
            output_file_path=staged_bundle_file_path,
        )

        validate_pfam_hmm_bundle_file(
            bundle_file_path=staged_bundle_file_path,
            expected_models=validation.registry.models,
            pyhmmer_module=resolved_pyhmmer_module,
        )

        if press_bundle:
            hmmer = getattr(resolved_pyhmmer_module, "hmmer", None)
            hmmpress = getattr(hmmer, "hmmpress", None)
            if not callable(hmmpress):
                raise RuntimeError(
                    "The provided PyHMMER module does not expose hmmer.hmmpress."
                )
            hmmpress(
                _iter_source_hmm_models(
                    validation=validation,
                    pyhmmer_module=resolved_pyhmmer_module,
                ),
                staged_bundle_file_path,
            )
            for pressed_file_path in staged_pressed_file_paths:
                if not pressed_file_path.is_file() or pressed_file_path.stat().st_size < 1:
                    raise RuntimeError(
                        "PyHMMER hmmpress did not produce a complete non-empty "
                        f"pressed bundle: {pressed_file_path}."
                    )

        published_file_paths: list[Path] = []
        try:
            staged_bundle_file_path.replace(bundle_file_path)
            published_file_paths.append(bundle_file_path)
            if press_bundle:
                final_pressed_file_paths: list[Path] = []
                for staged_pressed_file_path in staged_pressed_file_paths:
                    final_pressed_file_path = Path(
                        str(bundle_file_path)
                        + staged_pressed_file_path.name.removeprefix(bundle_file_name)
                    )
                    staged_pressed_file_path.replace(final_pressed_file_path)
                    published_file_paths.append(final_pressed_file_path)
                    final_pressed_file_paths.append(final_pressed_file_path)
                pressed_file_paths = tuple(final_pressed_file_paths)
        except Exception:
            for published_file_path in published_file_paths:
                published_file_path.unlink(missing_ok=True)
            raise
    finally:
        shutil.rmtree(staging_directory, ignore_errors=True)

    return PfamHmmBundleBuildResult(
        validation=validation,
        bundle_file_path=bundle_file_path,
        pressed_file_paths=pressed_file_paths,
        bundle_sha256=sha256_of_file(input_file_path=bundle_file_path),
    )


def build_pfam_hmm_bundle_derivation_identity_sha256(
    *,
    validation: PfamHmmReferenceValidation,
    pressed_outputs_required: bool = True,
) -> str:
    payload: dict[str, Any] = {
        "source_registry_sha256": validation.registry.registry_sha256,
        "source_database": validation.registry.source_database,
        "source_release": validation.registry.source_release,
        "canonical_models": [
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
        ],
        "concatenation_policy_version": CONCATENATION_POLICY_VERSION,
        "pyhmmer_version": validation.pyhmmer_version,
        "pressed_outputs_required": bool(pressed_outputs_required),
    }
    canonical_json = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical_json.encode("utf-8")).hexdigest()
