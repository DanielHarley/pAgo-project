from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from scripts.verify_raw_data import _verify_xml_file


class VerifyRawDataTests(unittest.TestCase):
    def test_materialized_xml_is_verified_from_file_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_directory = Path(temporary_directory_name)
            xml_file_path = snapshot_directory / "protein_records.xml"
            xml_bytes = b"<GBSet><GBSeq /></GBSet>\n"
            xml_file_path.write_bytes(xml_bytes)

            expected_sha256 = hashlib.sha256(xml_bytes).hexdigest()
            manifest = {
                "xml_file_name": xml_file_path.name,
                "xml_file_sha256": expected_sha256,
            }

            self.assertEqual(_verify_xml_file(snapshot_directory, manifest), [])

    def test_matching_lfs_pointer_is_verified_without_materialized_xml(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_directory = Path(temporary_directory_name)
            xml_file_path = snapshot_directory / "protein_records.xml"
            expected_sha256 = (
                "ada0563932d20f54f68c614df3c40d4fa3520038ad438a41cc471a1eacb970bf"
            )
            xml_file_path.write_text(
                "version https://git-lfs.github.com/spec/v1\n"
                f"oid sha256:{expected_sha256}\n"
                "size 363034184\n",
                encoding="utf-8",
            )
            manifest = {
                "xml_file_name": xml_file_path.name,
                "xml_file_sha256": expected_sha256,
            }

            self.assertEqual(_verify_xml_file(snapshot_directory, manifest), [])

    def test_mismatching_lfs_pointer_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_directory = Path(temporary_directory_name)
            xml_file_path = snapshot_directory / "protein_records.xml"
            expected_sha256 = "a" * 64
            actual_sha256 = "b" * 64
            xml_file_path.write_text(
                "version https://git-lfs.github.com/spec/v1\n"
                f"oid sha256:{actual_sha256}\n"
                "size 363034184\n",
                encoding="utf-8",
            )
            manifest = {
                "xml_file_name": xml_file_path.name,
                "xml_file_sha256": expected_sha256,
            }

            errors = _verify_xml_file(snapshot_directory, manifest)

            self.assertEqual(len(errors), 1)
            self.assertIn("XML LFS pointer SHA-256 mismatch", errors[0])
            self.assertIn(expected_sha256, errors[0])
            self.assertIn(actual_sha256, errors[0])


if __name__ == "__main__":
    unittest.main()
