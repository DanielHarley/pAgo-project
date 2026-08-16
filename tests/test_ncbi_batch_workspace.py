import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.error import URLError

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

import src.pago_pipeline.ncbi_api as ncbi_api
from src.pago_pipeline.ncbi_api import (
    NCBIXmlBatchFetchError,
    fetch_ncbi_protein_xml_batches,
)
from src.pago_pipeline.ncbi_batch_workspace import (
    BATCH_PLAN_FILE_NAME,
    NCBIXmlBatchWorkspace,
    build_xml_batch_plan,
)


class FakeFetchHandle:
    def __init__(self, payload: bytes) -> None:
        self.payload = payload
        self.closed = False
        self.url = "https://ncbi.example.test/efetch"
        self.headers = {"content-type": "text/xml"}

    def read(self):
        return self.payload

    def close(self) -> None:
        self.closed = True


def build_gbset_payload(protein_uids: list[str]) -> bytes:
    gbseq_records = b"".join(
        f"<GBSeq><GBSeqid>gi|{protein_uid}</GBSeqid></GBSeq>".encode("utf-8")
        for protein_uid in protein_uids
    )
    return b"<GBSet>" + gbseq_records + b"</GBSet>"


class BatchPlanTests(unittest.TestCase):
    def test_plan_splits_uids_into_index_addressed_batches(self) -> None:
        planned_batches = build_xml_batch_plan(
            protein_uids=["1", "2", "3", "4", "5"],
            batch_size=2,
        )

        self.assertEqual(
            [
                (
                    planned_batch.batch_index,
                    planned_batch.batch_start_index,
                    planned_batch.batch_end_index,
                    planned_batch.protein_uids,
                )
                for planned_batch in planned_batches
            ],
            [
                (1, 0, 1, ["1", "2"]),
                (2, 2, 3, ["3", "4"]),
                (3, 4, 4, ["5"]),
            ],
        )


class BatchWorkspaceReuseTests(unittest.TestCase):
    def setUp(self) -> None:
        with ncbi_api._default_ncbi_xml_circuit_breakers_lock:
            ncbi_api._default_ncbi_xml_circuit_breakers.clear()

        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        self.workspace_directory = Path(temporary_directory.name) / "workspace"

        sleep_patcher = patch("src.pago_pipeline.ncbi_api.time.sleep")
        sleep_patcher.start()
        self.addCleanup(sleep_patcher.stop)

        jitter_patcher = patch(
            "src.pago_pipeline.ncbi_api.random.uniform",
            return_value=0.0,
        )
        jitter_patcher.start()
        self.addCleanup(jitter_patcher.stop)

        self.protein_uids = ["1001", "1002", "1003", "1004"]

    def fetch_batches(self, *, efetch_side_effect):
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=efetch_side_effect,
        ) as efetch_mock:
            fetch_result = fetch_ncbi_protein_xml_batches(
                ncbi_email="test@example.org",
                ncbi_api_key=None,
                protein_uids=self.protein_uids,
                batch_size=2,
                max_retry_attempts=1,
                request_delay_seconds=0.0,
                batch_workspace_directory=self.workspace_directory,
            )

        return fetch_result, efetch_mock

    def run_failing_first_attempt(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[
                FakeFetchHandle(build_gbset_payload(["1001", "1002"])),
                URLError("temporarily unavailable"),
            ],
        ):
            with self.assertRaises(NCBIXmlBatchFetchError):
                fetch_ncbi_protein_xml_batches(
                    ncbi_email="test@example.org",
                    ncbi_api_key=None,
                    protein_uids=self.protein_uids,
                    batch_size=2,
                    max_retry_attempts=1,
                    request_delay_seconds=0.0,
                    batch_workspace_directory=self.workspace_directory,
                )

    def test_failed_run_persists_completed_batches_for_reuse(self) -> None:
        self.run_failing_first_attempt()

        persisted_batch_files = sorted(
            path.name
            for path in (self.workspace_directory / "batches").iterdir()
        )
        self.assertEqual(
            persisted_batch_files,
            ["batch_000001.json", "batch_000001.xml"],
        )
        self.assertTrue(
            (self.workspace_directory / BATCH_PLAN_FILE_NAME).exists()
        )

    def test_second_run_refetches_only_the_missing_batch(self) -> None:
        self.run_failing_first_attempt()

        fetch_result, efetch_mock = self.fetch_batches(
            efetch_side_effect=[
                FakeFetchHandle(build_gbset_payload(["1003", "1004"])),
            ],
        )

        efetch_mock.assert_called_once()
        self.assertEqual(fetch_result.reused_batch_count, 1)
        self.assertEqual(fetch_result.fetched_batch_count, 1)
        self.assertEqual(
            [xml_batch.batch_index for xml_batch in fetch_result.xml_batches],
            [1, 2],
        )
        self.assertTrue(fetch_result.xml_batches[0].reused_from_workspace)
        self.assertFalse(fetch_result.xml_batches[1].reused_from_workspace)

    def test_corrupt_persisted_batch_is_discarded_and_refetched(self) -> None:
        self.run_failing_first_attempt()

        corrupted_batch_file_path = (
            self.workspace_directory / "batches" / "batch_000001.xml"
        )
        corrupted_batch_file_path.write_bytes(b"<GBSet><GBSeq>")

        fetch_result, efetch_mock = self.fetch_batches(
            efetch_side_effect=[
                FakeFetchHandle(build_gbset_payload(["1001", "1002"])),
                FakeFetchHandle(build_gbset_payload(["1003", "1004"])),
            ],
        )

        self.assertEqual(efetch_mock.call_count, 2)
        self.assertEqual(fetch_result.reused_batch_count, 0)

    def test_persisted_batch_with_mismatched_uid_metadata_is_discarded(
        self,
    ) -> None:
        self.run_failing_first_attempt()

        batch_metadata_file_path = (
            self.workspace_directory / "batches" / "batch_000001.json"
        )
        batch_metadata = json.loads(
            batch_metadata_file_path.read_text(encoding="utf-8")
        )
        batch_metadata["protein_uids_sha256"] = "0" * 64
        batch_metadata_file_path.write_text(
            json.dumps(batch_metadata),
            encoding="utf-8",
        )

        fetch_result, efetch_mock = self.fetch_batches(
            efetch_side_effect=[
                FakeFetchHandle(build_gbset_payload(["1001", "1002"])),
                FakeFetchHandle(build_gbset_payload(["1003", "1004"])),
            ],
        )

        self.assertEqual(efetch_mock.call_count, 2)
        self.assertEqual(fetch_result.reused_batch_count, 0)

    def test_incompatible_batch_plan_resets_the_workspace(self) -> None:
        self.run_failing_first_attempt()

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[
                FakeFetchHandle(
                    build_gbset_payload(["1001", "1002", "1003", "1004"])
                ),
            ],
        ) as efetch_mock:
            fetch_result = fetch_ncbi_protein_xml_batches(
                ncbi_email="test@example.org",
                ncbi_api_key=None,
                protein_uids=self.protein_uids,
                batch_size=4,
                max_retry_attempts=1,
                request_delay_seconds=0.0,
                batch_workspace_directory=self.workspace_directory,
            )

        efetch_mock.assert_called_once()
        self.assertEqual(fetch_result.reused_batch_count, 0)
        self.assertFalse(
            (
                self.workspace_directory / "batches" / "batch_000002.xml"
            ).exists()
        )

    def test_resume_can_be_disabled_without_disabling_spill(self) -> None:
        self.run_failing_first_attempt()

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[
                FakeFetchHandle(build_gbset_payload(["1001", "1002"])),
                FakeFetchHandle(build_gbset_payload(["1003", "1004"])),
            ],
        ) as efetch_mock:
            fetch_result = fetch_ncbi_protein_xml_batches(
                ncbi_email="test@example.org",
                ncbi_api_key=None,
                protein_uids=self.protein_uids,
                batch_size=2,
                max_retry_attempts=1,
                request_delay_seconds=0.0,
                batch_workspace_directory=self.workspace_directory,
                enable_batch_resume=False,
            )

        self.assertEqual(efetch_mock.call_count, 2)
        self.assertEqual(fetch_result.reused_batch_count, 0)
        for xml_batch in fetch_result.xml_batches:
            self.assertIsNotNone(xml_batch.xml_payload_file_path)
            self.assertIsNone(xml_batch.xml_payload_bytes)


class EphemeralWorkspaceTests(unittest.TestCase):
    def test_ephemeral_workspace_does_not_write_a_batch_plan(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            workspace = NCBIXmlBatchWorkspace(
                workspace_directory=None,
                plan_identity={"protein_uids_sha256": "abc"},
                ephemeral_parent_directory=temporary_directory_name,
            )
            self.addCleanup(workspace.purge)
            workspace.open(planned_batches=[])

            self.assertTrue(workspace.is_ephemeral)
            self.assertFalse(workspace.resume_enabled)
            self.assertFalse(workspace.batch_plan_file_path.exists())
            self.assertTrue(
                Path(workspace.workspace_directory).is_relative_to(
                    Path(temporary_directory_name)
                )
            )


if __name__ == "__main__":
    unittest.main()
