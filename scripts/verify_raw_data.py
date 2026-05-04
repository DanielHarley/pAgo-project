from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.pago_pipeline.storage import sha256_of_file, sha256_of_lines


def _read_nonempty_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _verify_uid_file(snapshot_directory: Path, manifest: dict[str, object]) -> list[str]:
    expected_sha256 = manifest.get("protein_uids_sha256")
    uid_file_name = manifest.get("protein_uids_file_name") or manifest.get("snapshot_file_name")

    if not isinstance(expected_sha256, str) or not isinstance(uid_file_name, str):
        return []

    uid_file_path = snapshot_directory / uid_file_name
    if not uid_file_path.exists():
        return [f"Missing UID file: {uid_file_path}"]

    actual_sha256 = sha256_of_lines(
        text_lines=_read_nonempty_lines(uid_file_path),
        deduplicate_lines_preserving_order=False,
        sort_lines=False,
    )
    if actual_sha256 != expected_sha256:
        return [
            "UID SHA-256 mismatch: "
            f"{uid_file_path} expected {expected_sha256}, got {actual_sha256}"
        ]

    return []


def _verify_xml_file(snapshot_directory: Path, manifest: dict[str, object]) -> list[str]:
    expected_sha256 = manifest.get("xml_file_sha256")
    xml_file_name = manifest.get("xml_file_name")

    if not isinstance(expected_sha256, str) or not isinstance(xml_file_name, str):
        return []

    xml_file_path = snapshot_directory / xml_file_name
    if not xml_file_path.exists():
        return [f"Missing XML file: {xml_file_path}"]

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

        snapshot_directory = manifest_file_path.parent
        snapshot_errors = [
            *_verify_uid_file(snapshot_directory, manifest),
            *_verify_xml_file(snapshot_directory, manifest),
        ]

        if snapshot_errors:
            errors.extend(snapshot_errors)
        else:
            checked_file_count += int("protein_uids_sha256" in manifest)
            checked_file_count += int("xml_file_sha256" in manifest)

    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1

    print(
        f"Verified {checked_file_count} raw data file hashes "
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
