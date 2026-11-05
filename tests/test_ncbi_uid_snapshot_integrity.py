
import json
import tempfile
import unittest
from pathlib import Path

from src.pago_pipeline.ncbi_snapshot import (
    load_snapshot_by_directory,
)
from src.pago_pipeline.storage import sha256_of_lines


class TestNCBIUidSnapshotIntegrity(unittest.TestCase):

    def _create_test_snapshot(
        self,
        directory,
        uids,
    ):
        snapshot_directory = Path(directory)

        uid_hash = sha256_of_lines(
            text_lines=uids,
            deduplicate_lines_preserving_order=False,
            sort_lines=False,
        )

        manifest = {
            "normalized_protein_uid_count": len(uids),
            "protein_uids_sha256": uid_hash,
        }

        manifest_path = snapshot_directory / "manifest.json"
        uids_path = snapshot_directory / "protein_uids.txt"

        manifest_path.write_text(
            json.dumps(manifest),
            encoding="utf-8",
        )

        uids_path.write_text(
            "\n".join(uids) + "\n",
            encoding="utf-8",
        )

        return snapshot_directory, uids_path

    def test_accept_valid_uid_snapshot(self):
        with tempfile.TemporaryDirectory() as temp_dir:

            expected_uids = ["101", "102", "103", "104"]

            snapshot_directory, _ = self._create_test_snapshot(
                temp_dir,
                expected_uids,
            )

            result = load_snapshot_by_directory(
                snapshot_directory=snapshot_directory,
            )

            self.assertEqual(
                result["protein_uids"],
                expected_uids,
            )

            self.assertEqual(
                result["manifest"]["normalized_protein_uid_count"],
                4,
            )

    def test_reject_uid_snapshot_with_count_mismatch(self):
        with tempfile.TemporaryDirectory() as temp_dir:

            expected_uids = ["101", "102", "103", "104"]

            snapshot_directory, uids_path = self._create_test_snapshot(
                temp_dir,
                expected_uids,
            )

            # Simulate a corrupted snapshot:
            # one UID is removed after manifest creation.
            uids_path.write_text(
                "101\n102\n103\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                RuntimeError,
                "UID snapshot integrity validation failed",
            ):
                load_snapshot_by_directory(
                    snapshot_directory=snapshot_directory,
                )

    def test_reject_uid_snapshot_with_sha256_mismatch(self):
        with tempfile.TemporaryDirectory() as temp_dir:

            expected_uids = ["101", "102", "103", "104"]

            snapshot_directory, uids_path = self._create_test_snapshot(
                temp_dir,
                expected_uids,
            )

            # Change one UID while preserving the original count.
            uids_path.write_text(
                "101\n102\n999\n104\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                RuntimeError,
                "SHA-256 mismatch",
            ):
                load_snapshot_by_directory(
                    snapshot_directory=snapshot_directory,
                )

    def test_reject_missing_or_invalid_uid_count(self):

        invalid_counts = [None, "4", -1, True, 4.0]

        for invalid_count in invalid_counts:
            with self.subTest(invalid_count=invalid_count):

                with tempfile.TemporaryDirectory() as temp_dir:

                    snapshot_directory, _ = self._create_test_snapshot(
                        temp_dir,
                        ["101", "102", "103", "104"],
                    )

                    manifest_path = (
                        snapshot_directory / "manifest.json"
                    )

                    manifest = json.loads(
                        manifest_path.read_text(encoding="utf-8")
                    )

                    if invalid_count is None:
                        manifest.pop("normalized_protein_uid_count")
                    else:
                        manifest["normalized_protein_uid_count"] = (
                            invalid_count
                        )

                    manifest_path.write_text(
                        json.dumps(manifest),
                        encoding="utf-8",
                    )

                    with self.assertRaisesRegex(
                        RuntimeError,
                        "normalized_protein_uid_count",
                    ):
                        load_snapshot_by_directory(
                            snapshot_directory=snapshot_directory,
                        )

    def test_reject_missing_or_invalid_sha256(self):

        invalid_hashes = [
            None,
            "abc",
            "g" * 64,
            123,
            "A" * 64,
        ]

        for invalid_hash in invalid_hashes:
            with self.subTest(invalid_hash=invalid_hash):

                with tempfile.TemporaryDirectory() as temp_dir:

                    snapshot_directory, _ = self._create_test_snapshot(
                        temp_dir,
                        ["101", "102", "103", "104"],
                    )

                    manifest_path = (
                        snapshot_directory / "manifest.json"
                    )

                    manifest = json.loads(
                        manifest_path.read_text(encoding="utf-8")
                    )

                    if invalid_hash is None:
                        manifest.pop("protein_uids_sha256")
                    else:
                        manifest["protein_uids_sha256"] = invalid_hash

                    manifest_path.write_text(
                        json.dumps(manifest),
                        encoding="utf-8",
                    )

                    with self.assertRaisesRegex(
                        RuntimeError,
                        "protein_uids_sha256",
                    ):
                        load_snapshot_by_directory(
                            snapshot_directory=snapshot_directory,
                        )

    def test_reject_uid_snapshot_with_reordered_uids(self):
        with tempfile.TemporaryDirectory() as temp_dir:

            original_uids = ["101", "102", "103", "104"]

            snapshot_directory, uids_path = self._create_test_snapshot(
                temp_dir,
                original_uids,
            )

            # Keep all UIDs but change their original order.
            reordered_uids = ["102", "101", "103", "104"]

            uids_path.write_text(
                "\n".join(reordered_uids) + "\n",
                encoding="utf-8",
            )

            # The UID count is unchanged, but the canonical
            # SHA-256 must differ because the order changed.
            with self.assertRaisesRegex(
                RuntimeError,
                "SHA-256 mismatch",
            ):
                load_snapshot_by_directory(
                    snapshot_directory=snapshot_directory,
                )

    def test_accept_uid_snapshot_with_different_line_endings(self):
        expected_uids = ["101", "102", "103", "104"]

        # Test Unix/Linux (LF) and Windows (CRLF).
        for line_ending in ("\n", "\r\n"):
            with self.subTest(line_ending=repr(line_ending)):

                with tempfile.TemporaryDirectory() as temp_dir:

                    snapshot_directory, uids_path = (
                        self._create_test_snapshot(
                            temp_dir,
                            expected_uids,
                        )
                    )

                    # Write the same UIDs using the selected
                    # line ending, preserving the exact bytes.
                    file_content = (
                        line_ending.join(expected_uids)
                        + line_ending
                    ).encode("utf-8")

                    uids_path.write_bytes(file_content)

                    # Confirm the intended bytes were written.
                    self.assertEqual(
                        uids_path.read_bytes(),
                        file_content,
                    )

                    # Both representations must be accepted.
                    result = load_snapshot_by_directory(
                        snapshot_directory=snapshot_directory,
                    )

                    self.assertEqual(
                        result["protein_uids"],
                        expected_uids,
                    )

                    self.assertEqual(
                        result["manifest"][
                            "normalized_protein_uid_count"
                        ],
                        4,
                    )

                    expected_sha256 = sha256_of_lines(
                        text_lines=expected_uids,
                        deduplicate_lines_preserving_order=False,
                        sort_lines=False,
                    )

                    self.assertEqual(
                        result["manifest"]["protein_uids_sha256"],
                        expected_sha256,
                    )

if __name__ == "__main__":
    unittest.main()
