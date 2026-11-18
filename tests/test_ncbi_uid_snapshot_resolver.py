
import json
import tempfile
import unittest

from pathlib import Path
from unittest.mock import patch

from src.pago_pipeline.ncbi_snapshot import (
    SnapshotMode,
    resolve_ncbi_protein_uid_snapshot,
)
from src.pago_pipeline.storage import sha256_of_lines


class TestNCBIUidSnapshotResolver(unittest.TestCase):

    def _create_valid_latest_snapshot(self, snapshot_root):
        """
        Create a valid temporary latest/ UID snapshot.
        """
        latest_directory = snapshot_root / "latest"
        latest_directory.mkdir()

        uids = ["101", "102", "103", "104"]

        uid_hash = sha256_of_lines(
            text_lines=uids,
            deduplicate_lines_preserving_order=False,
            sort_lines=False,
        )

        manifest = {
            "normalized_protein_uid_count": len(uids),
            "protein_uids_sha256": uid_hash,
        }

        manifest_path = latest_directory / "manifest.json"
        uids_path = latest_directory / "protein_uids.txt"

        manifest_path.write_text(
            json.dumps(manifest),
            encoding="utf-8",
        )

        uids_path.write_text(
            "\n".join(uids) + "\n",
            encoding="utf-8",
        )

        return manifest_path, uids_path

    def test_reuse_or_create_rejects_invalid_sha256(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot_root = Path(temp_dir)

            _, uids_path = self._create_valid_latest_snapshot(
                snapshot_root
            )

            # Change one UID without changing the count.
            uids_path.write_text(
                "101\n102\n999\n104\n",
                encoding="utf-8",
            )

            with patch(
                "src.pago_pipeline.ncbi_snapshot."
                "fetch_ncbi_protein_uid_snapshot"
            ) as mock_fetch:

                with self.assertRaisesRegex(
                    RuntimeError,
                    "SHA-256 mismatch",
                ):
                    resolve_ncbi_protein_uid_snapshot(
                        snapshot_mode=SnapshotMode.reuse_latest_or_create,
                        snapshot_root_directory=snapshot_root,
                    )

                mock_fetch.assert_not_called()

    def test_reuse_latest_rejects_corrupted_snapshot(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot_root = Path(temp_dir)

            manifest_path, _ = self._create_valid_latest_snapshot(
                snapshot_root
            )

            # Replace the expected hash with an incorrect,
            # but syntactically valid, SHA-256 value.
            manifest = json.loads(
                manifest_path.read_text(encoding="utf-8")
            )

            manifest["protein_uids_sha256"] = "0" * 64

            manifest_path.write_text(
                json.dumps(manifest),
                encoding="utf-8",
            )

            with patch(
                "src.pago_pipeline.ncbi_snapshot."
                "fetch_ncbi_protein_uid_snapshot"
            ) as mock_fetch:

                with self.assertRaisesRegex(
                    RuntimeError,
                    "SHA-256 mismatch",
                ):
                    resolve_ncbi_protein_uid_snapshot(
                        snapshot_mode=SnapshotMode.reuse_latest,
                        snapshot_root_directory=snapshot_root,
                    )

                mock_fetch.assert_not_called()

    def test_reuse_or_create_rejects_invalid_uid_count(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot_root = Path(temp_dir)

            _, uids_path = self._create_valid_latest_snapshot(
                snapshot_root
            )

            # Remove one UID without updating the manifest.
            uids_path.write_text(
                "101\n102\n103\n",
                encoding="utf-8",
            )

            with patch(
                "src.pago_pipeline.ncbi_snapshot."
                "fetch_ncbi_protein_uid_snapshot"
            ) as mock_fetch:

                with self.assertRaisesRegex(
                    RuntimeError,
                    "UID count mismatch",
                ):
                    resolve_ncbi_protein_uid_snapshot(
                        snapshot_mode=SnapshotMode.reuse_latest_or_create,
                        snapshot_root_directory=snapshot_root,
                    )

                mock_fetch.assert_not_called()

    def test_reuse_valid_latest_without_ncbi_request(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot_root = Path(temp_dir)
            latest_directory = snapshot_root / "latest"
            latest_directory.mkdir()

            uids = ["101", "102", "103", "104"]

            uid_hash = sha256_of_lines(
                text_lines=uids,
                deduplicate_lines_preserving_order=False,
                sort_lines=False,
            )

            manifest = {
                "normalized_protein_uid_count": len(uids),
                "protein_uids_sha256": uid_hash,
            }

            (latest_directory / "manifest.json").write_text(
                json.dumps(manifest),
                encoding="utf-8",
            )

            (latest_directory / "protein_uids.txt").write_text(
                "\n".join(uids) + "\n",
                encoding="utf-8",
            )

            # A valid snapshot must be reused locally.
            with patch(
                "src.pago_pipeline.ncbi_snapshot."
                "fetch_ncbi_protein_uid_snapshot"
            ) as mock_fetch:

                result = resolve_ncbi_protein_uid_snapshot(
                    snapshot_mode=SnapshotMode.reuse_latest_or_create,
                    snapshot_root_directory=snapshot_root,
                )

                mock_fetch.assert_not_called()

            self.assertEqual(result["protein_uids"], uids)
            self.assertEqual(result["manifest"], manifest)

    def test_incomplete_latest_fails_without_ncbi_request(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot_root = Path(temp_dir)
            latest_directory = snapshot_root / "latest"
            latest_directory.mkdir()

            # The directory exists, but protein_uids.txt is missing.
            (latest_directory / "manifest.json").write_text(
                "{}",
                encoding="utf-8",
            )

            with patch(
                "src.pago_pipeline.ncbi_snapshot."
                "fetch_ncbi_protein_uid_snapshot",
                side_effect=AssertionError(
                    "Unexpected NCBI request"
                ),
            ) as mock_fetch:

                with self.assertRaisesRegex(
                    RuntimeError,
                    "Incomplete latest UID snapshot",
                ):
                    resolve_ncbi_protein_uid_snapshot(
                        snapshot_mode=SnapshotMode.reuse_latest_or_create,
                        snapshot_root_directory=snapshot_root,
                        search_query="PIWI[All Fields] AND Bacteria[Organism]",
                        ncbi_email="test@example.org",
                    )

                mock_fetch.assert_not_called()

    def test_reuse_or_create_fetches_when_latest_absent(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot_root = Path(temp_dir)

            # Confirm that latest/ does not exist.
            self.assertFalse(
                (snapshot_root / "latest").exists()
            )

            query = "PIWI[All Fields] AND Bacteria[Organism]"

            # Simulated output directory for the new snapshot.
            saved_snapshot_directory = (
                snapshot_root / "snapshots" / "mock-snapshot"
            )

            expected_payload = {
                "snapshot_directory": saved_snapshot_directory,
                "protein_uids": ["101", "102"],
            }

            # Mock the three external operations.
            with (
                patch(
                    "src.pago_pipeline.ncbi_snapshot."
                    "fetch_ncbi_protein_uid_snapshot"
                ) as mock_fetch,
                patch(
                    "src.pago_pipeline.ncbi_snapshot."
                    "save_ncbi_protein_uid_snapshot"
                ) as mock_save,
                patch(
                    "src.pago_pipeline.ncbi_snapshot."
                    "load_snapshot_by_directory"
                ) as mock_load,
            ):
                mock_save.return_value = saved_snapshot_directory
                mock_load.return_value = expected_payload

                result = resolve_ncbi_protein_uid_snapshot(
                    snapshot_mode=SnapshotMode.reuse_latest_or_create,
                    snapshot_root_directory=snapshot_root,
                    search_query=query,
                    ncbi_email="test@example.org",
                )

                # The NCBI fetch must be initiated.
                mock_fetch.assert_called_once()

                self.assertEqual(
                    mock_fetch.call_args.kwargs["query"],
                    query,
                )

                self.assertEqual(
                    mock_fetch.call_args.kwargs["ncbi_email"],
                    "test@example.org",
                )

                # The fetched result must be sent to storage.
                mock_save.assert_called_once_with(
                    fetch_result=mock_fetch.return_value,
                    snapshot_root_directory=snapshot_root,
                    update_latest_directory=True,
                )

                # The saved snapshot must be loaded.
                mock_load.assert_called_once_with(
                    snapshot_directory=saved_snapshot_directory,
                )

                # Resolve must return the loaded payload.
                self.assertIs(result, expected_payload)

    def test_create_new_fetches_even_when_latest_is_valid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            snapshot_root = Path(temp_dir)

            # Create an existing valid latest/ snapshot.
            manifest_path, uids_path = (
                self._create_valid_latest_snapshot(snapshot_root)
            )

            self.assertTrue(manifest_path.is_file())
            self.assertTrue(uids_path.is_file())

            query = "PIWI[All Fields] AND Bacteria[Organism]"

            # Simulated directory for a newly created snapshot.
            new_snapshot_directory = (
                snapshot_root / "snapshots" / "new-mock-snapshot"
            )

            expected_payload = {
                "snapshot_directory": new_snapshot_directory,
                "protein_uids": ["201", "202"],
            }

            with (
                patch(
                    "src.pago_pipeline.ncbi_snapshot."
                    "load_latest_snapshot"
                ) as mock_load_latest,
                patch(
                    "src.pago_pipeline.ncbi_snapshot."
                    "fetch_ncbi_protein_uid_snapshot"
                ) as mock_fetch,
                patch(
                    "src.pago_pipeline.ncbi_snapshot."
                    "save_ncbi_protein_uid_snapshot"
                ) as mock_save,
                patch(
                    "src.pago_pipeline.ncbi_snapshot."
                    "load_snapshot_by_directory"
                ) as mock_load,
            ):
                mock_save.return_value = new_snapshot_directory
                mock_load.return_value = expected_payload

                result = resolve_ncbi_protein_uid_snapshot(
                    snapshot_mode=SnapshotMode.create_new,
                    snapshot_root_directory=snapshot_root,
                    search_query=query,
                    ncbi_email="test@example.org",
                    update_latest_directory=True,
                )

                # The existing latest/ must not be reused.
                mock_load_latest.assert_not_called()

                # A new NCBI retrieval must be initiated.
                mock_fetch.assert_called_once()

                self.assertEqual(
                    mock_fetch.call_args.kwargs["query"],
                    query,
                )

                # The new result must be sent to storage.
                mock_save.assert_called_once_with(
                    fetch_result=mock_fetch.return_value,
                    snapshot_root_directory=snapshot_root,
                    update_latest_directory=True,
                )

                # The newly saved snapshot must be loaded.
                mock_load.assert_called_once_with(
                    snapshot_directory=new_snapshot_directory,
                )

                # Resolve must return the new snapshot payload.
                self.assertIs(result, expected_payload)

if __name__ == "__main__":
    unittest.main()
