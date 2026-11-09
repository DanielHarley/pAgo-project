from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.pago_pipeline.ncbi_api import (
    validate_ncbi_protein_identifiers,
)
from src.pago_pipeline.ncbi_snapshot import (
    _identifier_type_from_snapshot_manifest,
    _protein_identifier_file_path_from_manifest,
    _snapshot_format_version,
    _validate_loaded_xml_snapshot_payload,
    _validate_saved_identifier_list,
    _validate_snapshot_integrity_fields,
    _xml_file_path_from_manifest,
)
from src.pago_pipeline.storage import sha256_of_file


LFS_POINTER_VERSION_LINE = "version https://git-lfs.github.com/spec/v1"
LFS_POINTER_MAX_BYTES = 1024


def _read_lfs_pointer_sha256(path: Path) -> str | None:
    if path.stat().st_size > LFS_POINTER_MAX_BYTES:
        return None

    try:
        pointer_text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return None

    pointer_lines = [
        line.strip() for line in pointer_text.splitlines() if line.strip()
    ]
    if not pointer_lines or pointer_lines[0] != LFS_POINTER_VERSION_LINE:
        return None

    oid_prefix = "oid sha256:"
    for line in pointer_lines[1:]:
        if not line.startswith(oid_prefix):
            continue

        candidate_sha256 = line[len(oid_prefix) :].strip().lower()
        if len(candidate_sha256) != 64:
            return None
        if any(character not in "0123456789abcdef" for character in candidate_sha256):
            return None
        return candidate_sha256

    return None


def _read_nonempty_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _has_identifier_data(manifest: dict[str, object]) -> bool:
    """Recognize generic or legacy identifier metadata, including partial data."""
    return manifest.get("snapshot_format_version") == "1.1" or manifest.get("artifact_type") in (
        "ncbi_protein_identifier_snapshot", "ncbi_protein_uid_snapshot", "ncbi_protein_xml_snapshot",
    ) or any(field in manifest for field in (
        "protein_identifiers_sha256", "protein_uids_sha256",
        "protein_identifiers_file_name", "protein_uids_file_name", "snapshot_file_name",
    ))


def _has_xml_data(manifest: dict[str, object]) -> bool:
    return (
        manifest.get("artifact_type") == "ncbi_protein_xml_snapshot"
        or "xml_file_sha256" in manifest
        or "xml_file_name" in manifest
    )


def _verify_uid_file(snapshot_directory: Path, manifest: dict[str, object]) -> list[str]:
    """Verify one identifier file in either domain; retain the legacy function name."""
    if not _has_identifier_data(manifest):
        return []
    hash_field = (
        "protein_identifiers_sha256"
        if "protein_identifiers_sha256" in manifest
        else "protein_uids_sha256"
    )
    if hash_field not in manifest:
        return [
            f"Missing protein identifier SHA-256 metadata: {snapshot_directory}: "
            "must define protein_identifiers_sha256 or protein_uids_sha256."
        ]
    expected_sha256 = manifest[hash_field]
    if not isinstance(expected_sha256, str) or not expected_sha256.strip():
        return [
            f"Invalid protein identifier SHA-256 metadata: {snapshot_directory}: "
            f"{hash_field} must be a non-empty string."
        ]

    try:
        identifier_file_path = _protein_identifier_file_path_from_manifest(snapshot_directory, manifest)
    except (OSError, RuntimeError, ValueError) as error:
        return [f"Invalid identifier file name: {snapshot_directory}: {error}"]

    if not identifier_file_path.is_file():
        return [f"Missing protein identifier file: {identifier_file_path}"]

    try:
        protein_identifiers = _read_nonempty_lines(identifier_file_path)
    except (OSError, UnicodeError) as error:
        return [f"Failed to read protein identifier file: {identifier_file_path}: {error}"]
    try:
        identifier_type = _identifier_type_from_snapshot_manifest(manifest)
        validate_ncbi_protein_identifiers(
            protein_identifiers=protein_identifiers,
            identifier_type=identifier_type,
            parameter_name="protein_identifiers",
        )
    except ValueError as error:
        return [
            "Identifier type validation failed: "
            f"{identifier_file_path}: {error}"
        ]

    try:
        _validate_saved_identifier_list(
            manifest_payload=manifest,
            protein_identifiers=protein_identifiers,
        )
    except RuntimeError as error:
        return [f"Invalid identifier snapshot: {snapshot_directory}: {error}"]

    return []


def _verify_xml_file(snapshot_directory: Path, manifest: dict[str, object]) -> list[str]:
    if not _has_xml_data(manifest):
        return []
    try:
        _validate_snapshot_integrity_fields(manifest, include_xml=True)
    except (RuntimeError, ValueError) as error:
        return [f"Invalid XML snapshot metadata: {snapshot_directory}: {error}"]
    expected_sha256 = manifest.get("xml_file_sha256")

    if not isinstance(expected_sha256, str) or not expected_sha256.strip():
        return [f"Invalid XML SHA-256 metadata: {snapshot_directory}: xml_file_sha256 must be a non-empty string."]
    try:
        xml_file_path = _xml_file_path_from_manifest(snapshot_directory, manifest)
    except (OSError, RuntimeError, ValueError) as error:
        return [f"Invalid XML file name: {snapshot_directory}: {error}"]
    if not xml_file_path.is_file():
        return [f"Missing XML file: {xml_file_path}"]

    lfs_pointer_sha256 = _read_lfs_pointer_sha256(xml_file_path)
    if lfs_pointer_sha256 is not None:
        if lfs_pointer_sha256 != expected_sha256:
            return [
                "XML LFS pointer SHA-256 mismatch: "
                f"{xml_file_path} expected {expected_sha256}, got {lfs_pointer_sha256}"
            ]
        return []

    if _snapshot_format_version(manifest) == "1.1":
        try:
            identifier_file_path = _protein_identifier_file_path_from_manifest(snapshot_directory, manifest)
            _validate_loaded_xml_snapshot_payload(
                snapshot_directory=snapshot_directory,
                manifest_payload=manifest,
                protein_uids=_read_nonempty_lines(identifier_file_path),
            )
        except (OSError, UnicodeError, ValueError, RuntimeError) as error:
            return [f"Invalid XML snapshot: {snapshot_directory}: {error}"]
        return []

    actual_sha256 = sha256_of_file(input_file_path=xml_file_path)
    if actual_sha256 != expected_sha256:
        return [
            "XML SHA-256 mismatch: "
            f"{xml_file_path} expected {expected_sha256}, got {actual_sha256}"
        ]

    return []


def verify_raw_data(raw_data_root: Path) -> int:
    manifest_file_paths = sorted(raw_data_root.rglob("manifest.json"))
    if not manifest_file_paths:
        print(f"No manifest.json files found under {raw_data_root}", file=sys.stderr)
        return 1

    errors: list[str] = []
    checked_file_count = 0

    for manifest_file_path in manifest_file_paths:
        with manifest_file_path.open("r", encoding="utf-8") as file_handle:
            manifest = json.load(file_handle)

        try:
            _snapshot_format_version(manifest)
        except ValueError as error:
            errors.append(f"Invalid snapshot manifest: {manifest_file_path}: {error}")
            continue

        snapshot_directory = manifest_file_path.parent
        snapshot_errors = [
            *_verify_uid_file(snapshot_directory, manifest),
            *_verify_xml_file(snapshot_directory, manifest),
        ]

        if snapshot_errors:
            errors.extend(snapshot_errors)
        else:
            # Both naming schemes describe one identity, not two separate files.
            # Each requested check above must succeed before it can be counted.
            checked_file_count += int(_has_identifier_data(manifest))
            checked_file_count += int(_has_xml_data(manifest))

    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1

    print(
        f"Verified {checked_file_count} raw data identities "
        f"from {len(manifest_file_paths)} manifests."
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify raw NCBI snapshot files against manifest SHA-256 hashes.",
    )
    parser.add_argument(
        "--raw-data-root",
        type=Path,
        default=Path("data") / "01-raw",
        help="Raw data root directory to verify.",
    )
    args = parser.parse_args()

    return verify_raw_data(args.raw_data_root)


if __name__ == "__main__":
    raise SystemExit(main())
