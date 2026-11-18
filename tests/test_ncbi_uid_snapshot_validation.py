import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

if "Bio" not in sys.modules:
    fake_bio_module = types.ModuleType("Bio")
    fake_bio_module.__version__ = "test"
    fake_entrez_module = types.ModuleType("Bio.Entrez")
    fake_entrez_module.urlopen = lambda *args, **kwargs: None
    fake_seqio_module = types.ModuleType("Bio.SeqIO")
    fake_bio_module.Entrez = fake_entrez_module
    fake_bio_module.SeqIO = fake_seqio_module
    sys.modules["Bio"] = fake_bio_module
    sys.modules["Bio.Entrez"] = fake_entrez_module
    sys.modules["Bio.SeqIO"] = fake_seqio_module

import src.pago_pipeline.ncbi_snapshot as ncbi_snapshot
from src.pago_pipeline.ncbi_api import NCBIProteinUidFetchResult
from src.pago_pipeline.ncbi_snapshot import (
    SnapshotValidationError,
    save_ncbi_protein_uid_snapshot,
)
from src.pago_pipeline.storage import sha256_of_lines

DEFAULT_PROTEIN_UIDS = ["1000", "200", "300"]
PREVIOUS_LATEST_FILES = {
    "protein_uids.txt": b"9\n",
    "manifest.json": b'{"previous": true}',
}


def build_fetch_result(
    *,
    protein_uids: list[str] | None = None,
    **overrides,
) -> NCBIProteinUidFetchResult:
    uids = list(DEFAULT_PROTEIN_UIDS if protein_uids is None else protein_uids)
    fields = dict(
        database_name="protein",
        search_query="PIWI[All Fields] AND Bacteria[Organism]",
        translated_query="PIWI[All Fields]",
        identifier_type="uid",
        retrieved_at_utc="2026-01-01T00:00:00Z",
        ncbi_reported_result_count=len(uids),
        protein_uids=uids,
        raw_protein_uid_count=len(uids),
        normalized_protein_uid_count=len(uids),
        deduplicate_uids=True,
        sort_uids=True,
        protein_uids_sha256=sha256_of_lines(
            text_lines=uids,
            deduplicate_lines_preserving_order=False,
            sort_lines=False,
        ),
        page_size=1000,
        request_delay_seconds=0.1,
        max_retry_attempts=5,
        history_web_env="WEB",
        history_query_key="1",
        python_version="3.14",
        biopython_version="test",
    )
    fields.update(overrides)
    return NCBIProteinUidFetchResult(**fields)


class UidSnapshotValidationTestCase(unittest.TestCase):
    def setUp(self) -> None:
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        self.snapshot_root_directory = Path(temporary_directory.name)
        self.latest_directory = self.snapshot_root_directory / "latest"
        self.snapshots_directory = self.snapshot_root_directory / "snapshots"

    def write_previous_latest(self) -> None:
        self.latest_directory.mkdir(parents=True)
        for file_name, content in PREVIOUS_LATEST_FILES.items():
            (self.latest_directory / file_name).write_bytes(content)

    def read_latest_files(self) -> dict[str, bytes]:
        return {
            path.name: path.read_bytes()
            for path in sorted(self.latest_directory.iterdir())
        }

    def list_snapshot_directories(self) -> list[Path]:
        if not self.snapshots_directory.exists():
            return []
        return sorted(self.snapshots_directory.iterdir())

    def save(self, fetch_result, **kwargs):
        return save_ncbi_protein_uid_snapshot(
            fetch_result=fetch_result,
            snapshot_root_directory=self.snapshot_root_directory,
            **kwargs,
        )

    def assert_rejected_and_nothing_published(
        self,
        *,
        check_name: str,
        action,
    ) -> None:
        with self.assertRaises(SnapshotValidationError) as raised:
            action()

        message = str(raised.exception)
        self.assertIn(f"'{check_name}'", message)
        self.assertIn("observed", message)
        self.assertIn("expected", message)
        self.assertIn("latest/ snapshot was kept unchanged", message)
        self.assertIn("Possible causes", message)
        self.assertIn("Next step", message)
        self.assertEqual(self.read_latest_files(), PREVIOUS_LATEST_FILES)
        self.assertEqual(self.list_snapshot_directories(), [])


class Stage1FetchResultValidationTests(UidSnapshotValidationTestCase):
    def setUp(self) -> None:
        super().setUp()
        self.write_previous_latest()

    def assert_stage_1_rejects(self, fetch_result, check_name: str) -> None:
        self.assert_rejected_and_nothing_published(
            check_name=check_name,
            action=lambda: self.save(fetch_result),
        )

    def test_zero_count_result_is_rejected(self) -> None:
        self.assert_stage_1_rejects(
            build_fetch_result(protein_uids=[]),
            "non_empty_uid_list",
        )

    def test_non_numeric_identifier_is_rejected(self) -> None:
        self.assert_stage_1_rejects(
            build_fetch_result(protein_uids=["1000", "12a"]),
            "protein_uid_format",
        )

    def test_accession_is_rejected(self) -> None:
        self.assert_stage_1_rejects(
            build_fetch_result(protein_uids=["1000", "WP_123456789.1"]),
            "protein_uid_format",
        )

    def test_identifier_with_whitespace_is_rejected(self) -> None:
        self.assert_stage_1_rejects(
            build_fetch_result(protein_uids=["1000", " 200"]),
            "protein_uid_format",
        )

    def test_raw_count_different_from_reported_count_is_rejected(self) -> None:
        self.assert_stage_1_rejects(
            build_fetch_result(ncbi_reported_result_count=10),
            "raw_count_matches_ncbi_reported_count",
        )

    def test_normalized_count_above_raw_count_is_rejected(self) -> None:
        self.assert_stage_1_rejects(
            build_fetch_result(
                raw_protein_uid_count=3,
                ncbi_reported_result_count=3,
                normalized_protein_uid_count=4,
            ),
            "normalized_count_not_above_raw_count",
        )

    def test_normalized_count_different_from_list_length_is_rejected(self) -> None:
        self.assert_stage_1_rejects(
            build_fetch_result(
                raw_protein_uid_count=5,
                ncbi_reported_result_count=5,
                normalized_protein_uid_count=5,
            ),
            "normalized_count_matches_uid_list",
        )

    def test_duplicate_in_normalized_list_is_rejected_when_deduplicating(
        self,
    ) -> None:
        self.assert_stage_1_rejects(
            build_fetch_result(protein_uids=["1000", "1000", "200"]),
            "unique_uids",
        )

    def test_out_of_order_list_is_rejected_when_sorting(self) -> None:
        self.assert_stage_1_rejects(
            build_fetch_result(protein_uids=["300", "200", "1000"]),
            "lexicographic_uid_order",
        )

    def test_hash_that_does_not_match_the_list_is_rejected(self) -> None:
        self.assert_stage_1_rejects(
            build_fetch_result(protein_uids_sha256="0" * 64),
            "uid_list_sha256",
        )


class AcceptedUidResultTests(UidSnapshotValidationTestCase):
    def test_raw_duplicates_removed_by_deduplication_are_accepted(self) -> None:
        fetch_result = build_fetch_result(
            protein_uids=["1000", "200", "300"],
            ncbi_reported_result_count=5,
            raw_protein_uid_count=5,
            normalized_protein_uid_count=3,
        )

        snapshot_directory = self.save(fetch_result)

        persisted_lines = (
            (snapshot_directory / "protein_uids.txt").read_text().splitlines()
        )
        self.assertEqual(persisted_lines, ["1000", "200", "300"])

    def test_lexicographic_order_that_is_not_numeric_order_is_accepted(
        self,
    ) -> None:
        fetch_result = build_fetch_result(protein_uids=["1000", "200", "300"])

        snapshot_directory = self.save(fetch_result)

        self.assertTrue((snapshot_directory / "manifest.json").exists())

    def test_unsorted_list_is_accepted_when_sorting_is_disabled(self) -> None:
        uids = ["300", "200", "1000"]
        fetch_result = build_fetch_result(protein_uids=uids, sort_uids=False)

        snapshot_directory = self.save(fetch_result)

        persisted_lines = (
            (snapshot_directory / "protein_uids.txt").read_text().splitlines()
        )
        self.assertEqual(persisted_lines, uids)


class PromotionTests(UidSnapshotValidationTestCase):
    def test_valid_result_is_promoted_to_latest(self) -> None:
        self.write_previous_latest()

        snapshot_directory = self.save(build_fetch_result())

        self.assertEqual(self.list_snapshot_directories(), [snapshot_directory])
        for file_name in ("protein_uids.txt", "manifest.json"):
            self.assertEqual(
                (self.latest_directory / file_name).read_bytes(),
                (snapshot_directory / file_name).read_bytes(),
            )

    def test_flag_disabled_creates_snapshot_and_leaves_latest_untouched(
        self,
    ) -> None:
        self.write_previous_latest()

        snapshot_directory = self.save(
            build_fetch_result(),
            update_latest_directory=False,
        )

        self.assertTrue((snapshot_directory / "manifest.json").exists())
        self.assertEqual(self.read_latest_files(), PREVIOUS_LATEST_FILES)

    def test_valid_result_is_promoted_when_no_latest_exists(self) -> None:
        snapshot_directory = self.save(build_fetch_result())

        self.assertTrue((self.latest_directory / "manifest.json").exists())
        self.assertTrue(snapshot_directory.exists())


class Stage2PersistedSnapshotValidationTests(UidSnapshotValidationTestCase):
    def setUp(self) -> None:
        super().setUp()
        self.write_previous_latest()

    def assert_stage_2_rejects(self, *, check_name: str, patches) -> None:
        with patches:
            self.assert_rejected_and_nothing_published(
                check_name=check_name,
                action=lambda: self.save(build_fetch_result()),
            )

    @staticmethod
    def write_lines_with(transform):
        def write_text_lines(*, text_lines, output_file_path):
            output_file_path = Path(output_file_path)
            output_file_path.parent.mkdir(parents=True, exist_ok=True)
            output_file_path.write_text(
                "".join(f"{line}\n" for line in transform(list(text_lines))),
                encoding="utf-8",
            )

        return write_text_lines

    @staticmethod
    def build_manifest_with(transform):
        original_build_manifest = ncbi_snapshot._build_snapshot_manifest

        def build_manifest(**kwargs):
            manifest_payload = original_build_manifest(**kwargs)
            transform(manifest_payload)
            return manifest_payload

        return build_manifest

    def test_truncated_uid_file_is_rejected(self) -> None:
        self.assert_stage_2_rejects(
            check_name="line_count_matches_manifest",
            patches=patch.object(
                ncbi_snapshot,
                "_write_text_lines",
                self.write_lines_with(lambda lines: lines[:1]),
            ),
        )

    def test_reordered_uid_file_is_rejected(self) -> None:
        self.assert_stage_2_rejects(
            check_name="lexicographic_uid_order",
            patches=patch.object(
                ncbi_snapshot,
                "_write_text_lines",
                self.write_lines_with(lambda lines: list(reversed(lines))),
            ),
        )

    def test_uid_file_with_a_blank_line_is_rejected(self) -> None:
        self.assert_stage_2_rejects(
            check_name="protein_uid_format",
            patches=patch.object(
                ncbi_snapshot,
                "_write_text_lines",
                self.write_lines_with(lambda lines: lines + [""]),
            ),
        )

    def test_manifest_count_that_disagrees_with_file_is_rejected(self) -> None:
        def raise_count(manifest_payload) -> None:
            manifest_payload["normalized_protein_uid_count"] += 1

        self.assert_stage_2_rejects(
            check_name="line_count_matches_manifest",
            patches=patch.object(
                ncbi_snapshot,
                "_build_snapshot_manifest",
                self.build_manifest_with(raise_count),
            ),
        )

    def test_manifest_hash_that_disagrees_with_file_is_rejected(self) -> None:
        def corrupt_hash(manifest_payload) -> None:
            manifest_payload["protein_uids_sha256"] = "0" * 64

        self.assert_stage_2_rejects(
            check_name="uid_file_sha256",
            patches=patch.object(
                ncbi_snapshot,
                "_build_snapshot_manifest",
                self.build_manifest_with(corrupt_hash),
            ),
        )

    def test_manifest_missing_a_required_key_is_rejected(self) -> None:
        def drop_hash(manifest_payload) -> None:
            manifest_payload.pop("protein_uids_sha256")

        self.assert_stage_2_rejects(
            check_name="manifest_required_keys",
            patches=patch.object(
                ncbi_snapshot,
                "_build_snapshot_manifest",
                self.build_manifest_with(drop_hash),
            ),
        )

    def test_manifest_with_other_format_version_is_rejected(self) -> None:
        def change_version(manifest_payload) -> None:
            manifest_payload["snapshot_format_version"] = "2.0"

        self.assert_stage_2_rejects(
            check_name="snapshot_format_version",
            patches=patch.object(
                ncbi_snapshot,
                "_build_snapshot_manifest",
                self.build_manifest_with(change_version),
            ),
        )

    def test_manifest_raw_count_different_from_reported_count_is_rejected(
        self,
    ) -> None:
        def change_reported_count(manifest_payload) -> None:
            manifest_payload["ncbi_reported_result_count"] = 99

        self.assert_stage_2_rejects(
            check_name="raw_count_matches_ncbi_reported_count",
            patches=patch.object(
                ncbi_snapshot,
                "_build_snapshot_manifest",
                self.build_manifest_with(change_reported_count),
            ),
        )

    def test_failed_stage_2_removes_the_files_it_had_written(self) -> None:
        original_write_text_lines = ncbi_snapshot._write_text_lines
        written_paths: list[Path] = []

        def record_and_truncate(*, text_lines, output_file_path) -> None:
            written_paths.append(Path(output_file_path))
            original_write_text_lines(
                text_lines=list(text_lines)[:1],
                output_file_path=output_file_path,
            )

        with patch.object(ncbi_snapshot, "_write_text_lines", record_and_truncate):
            with self.assertRaises(SnapshotValidationError):
                self.save(build_fetch_result())

        self.assertEqual(len(written_paths), 1)
        self.assertFalse(written_paths[0].exists())


if __name__ == "__main__":
    unittest.main()
