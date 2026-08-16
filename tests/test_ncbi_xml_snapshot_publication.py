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
import src.pago_pipeline.ncbi_snapshot as ncbi_snapshot
from src.pago_pipeline.ncbi_api import NCBIXmlBatchFetchError
from src.pago_pipeline.ncbi_snapshot import (
    BATCH_WORKSPACE_DIRECTORY_NAME,
    XML_SNAPSHOT_FORMAT_VERSION,
    resolve_ncbi_protein_xml_snapshot,
)
from src.pago_pipeline.ncbi_xml_stream import build_legacy_consolidated_xml_payload
from src.pago_pipeline.storage import sha256_of_file, sha256_of_lines

SEARCH_QUERY = "PIWI[All Fields] AND Bacteria[Organism]"


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
    gbseq_records = "".join(
        f"<GBSeq>\n  <GBSeqid>gi|{protein_uid}</GBSeqid>\n"
        f"  <GBSeq_locus>LOCUS_{protein_uid}</GBSeq_locus>\n</GBSeq>\n"
        for protein_uid in protein_uids
    )
    return (
        '<?xml version="1.0"?>\n<GBSet>\n' + gbseq_records + "</GBSet>\n"
    ).encode("utf-8")


class XmlSnapshotPublicationTests(unittest.TestCase):
    def setUp(self) -> None:
        with ncbi_api._default_ncbi_xml_circuit_breakers_lock:
            ncbi_api._default_ncbi_xml_circuit_breakers.clear()

        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        self.working_directory = Path(temporary_directory.name)
        self.uid_snapshot_root_directory = self.working_directory / "uid_snapshots"
        self.xml_snapshot_root_directory = self.working_directory / "xml_snapshots"

        self.protein_uids = [str(1000 + index) for index in range(6)]
        self._write_source_uid_snapshot()

        sleep_patcher = patch("src.pago_pipeline.ncbi_api.time.sleep")
        sleep_patcher.start()
        self.addCleanup(sleep_patcher.stop)

        jitter_patcher = patch(
            "src.pago_pipeline.ncbi_api.random.uniform",
            return_value=0.0,
        )
        jitter_patcher.start()
        self.addCleanup(jitter_patcher.stop)

    def _write_source_uid_snapshot(self) -> None:
        latest_directory = self.uid_snapshot_root_directory / "latest"
        latest_directory.mkdir(parents=True, exist_ok=True)

        (latest_directory / "protein_uids.txt").write_text(
            "".join(f"{protein_uid}\n" for protein_uid in self.protein_uids),
            encoding="utf-8",
            newline="\n",
        )
        (latest_directory / "manifest.json").write_text(
            json.dumps(
                {
                    "snapshot_format_version": "1.1",
                    "search_query": SEARCH_QUERY,
                    "translated_query": "translated",
                    "protein_uids_sha256": sha256_of_lines(
                        text_lines=self.protein_uids,
                        deduplicate_lines_preserving_order=False,
                        sort_lines=False,
                    ),
                    "normalized_protein_uid_count": len(self.protein_uids),
                    "immutable_snapshot_relative_path": "snapshots/uid-snapshot",
                    "immutable_snapshot_directory_name": "uid-snapshot",
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    def _fake_efetch(self, **parameters):
        return FakeFetchHandle(build_gbset_payload(list(parameters["protein_uids"])))

    def _resolve_snapshot(self, **overrides):
        resolve_parameters = {
            "snapshot_mode": "create_new",
            "snapshot_root_directory": self.xml_snapshot_root_directory,
            "source_uid_snapshot_root_directory": self.uid_snapshot_root_directory,
            "xml_batch_size": 2,
            "xml_request_delay_seconds": 0.0,
            "ncbi_email": "test@example.org",
        }
        resolve_parameters.update(overrides)
        return resolve_ncbi_protein_xml_snapshot(**resolve_parameters)

    def test_published_snapshot_matches_in_memory_consolidation_byte_for_byte(
        self,
    ) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            self._fake_efetch,
        ):
            snapshot_payload = self._resolve_snapshot()

        expected_payload, expected_record_count = (
            build_legacy_consolidated_xml_payload(
                batch_payloads=[
                    build_gbset_payload(self.protein_uids[0:2]),
                    build_gbset_payload(self.protein_uids[2:4]),
                    build_gbset_payload(self.protein_uids[4:6]),
                ]
            )
        )

        xml_file_path = Path(snapshot_payload["xml_file_path"])
        self.assertEqual(xml_file_path.read_bytes(), expected_payload)
        self.assertEqual(
            snapshot_payload["manifest"]["consolidated_record_count"],
            expected_record_count,
        )
        self.assertEqual(
            snapshot_payload["manifest"]["xml_file_sha256"],
            sha256_of_file(input_file_path=xml_file_path),
        )

    def test_publication_writes_the_complete_artifact_set(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            self._fake_efetch,
        ):
            snapshot_payload = self._resolve_snapshot()

        snapshot_directory = Path(snapshot_payload["snapshot_directory"])
        latest_directory = self.xml_snapshot_root_directory / "latest"

        for artifact_directory in (snapshot_directory, latest_directory):
            for artifact_file_name in (
                "manifest.json",
                "protein_uids.txt",
                "protein_records.xml",
            ):
                self.assertTrue(
                    (artifact_directory / artifact_file_name).exists(),
                    f"{artifact_file_name} missing from {artifact_directory}",
                )

        self.assertEqual(
            (latest_directory / "protein_records.xml").read_bytes(),
            (snapshot_directory / "protein_records.xml").read_bytes(),
        )

    def test_manifest_records_format_version_policy_and_telemetry(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            self._fake_efetch,
        ):
            snapshot_payload = self._resolve_snapshot()

        manifest = snapshot_payload["manifest"]
        self.assertEqual(
            manifest["snapshot_format_version"],
            XML_SNAPSHOT_FORMAT_VERSION,
        )
        self.assertEqual(manifest["rettype"], "gp")
        self.assertEqual(manifest["batch_size"], 2)
        self.assertEqual(manifest["batch_count"], 3)
        self.assertEqual(manifest["request_policy"]["max_concurrent_requests"], 1)

        retrieval_telemetry = manifest["retrieval_telemetry"]
        self.assertEqual(
            retrieval_telemetry["stages"]["xml_batches"]["request_count"],
            3,
        )
        self.assertIn(
            "consolidated_xml_write",
            retrieval_telemetry["local_phase_seconds"],
        )
        self.assertIn(
            "consolidated_xml_validation",
            retrieval_telemetry["local_phase_seconds"],
        )
        # The manifest is written before the latest/ copy, so the copy phase can
        # only appear in the payload telemetry, which is built after publication.
        self.assertNotIn(
            "latest_directory_copy",
            retrieval_telemetry["local_phase_seconds"],
        )
        self.assertIn(
            "latest_directory_copy",
            snapshot_payload["retrieval_telemetry"]["local_phase_seconds"],
        )

    def test_publication_validates_the_written_file_exactly_once(self) -> None:
        original_validator = ncbi_snapshot.validate_consolidated_xml_file
        validation_calls: list[dict] = []

        def counting_validator(**parameters):
            validation_calls.append(parameters)
            return original_validator(**parameters)

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            self._fake_efetch,
        ):
            with patch.object(
                ncbi_snapshot,
                "validate_consolidated_xml_file",
                counting_validator,
            ):
                self._resolve_snapshot()

        self.assertEqual(len(validation_calls), 1)
        self.assertEqual(
            validation_calls[0]["expected_protein_uids"],
            self.protein_uids,
        )

    def test_batch_workspace_is_purged_after_a_successful_publication(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            self._fake_efetch,
        ):
            self._resolve_snapshot()

        self.assertFalse(
            (
                self.xml_snapshot_root_directory / BATCH_WORKSPACE_DIRECTORY_NAME
            ).exists()
        )

    def test_failed_run_keeps_the_workspace_so_the_rerun_resumes(self) -> None:
        def fail_on_the_last_batch(**parameters):
            requested_protein_uids = list(parameters["protein_uids"])
            if requested_protein_uids == self.protein_uids[4:6]:
                raise URLError("temporarily unavailable")
            return FakeFetchHandle(build_gbset_payload(requested_protein_uids))

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            fail_on_the_last_batch,
        ):
            with self.assertRaises(NCBIXmlBatchFetchError):
                self._resolve_snapshot(max_retry_attempts=1)

        workspace_directory = (
            self.xml_snapshot_root_directory / BATCH_WORKSPACE_DIRECTORY_NAME
        )
        self.assertTrue(workspace_directory.exists())

        efetch_call_count = {"value": 0}

        def count_and_serve(**parameters):
            efetch_call_count["value"] += 1
            return self._fake_efetch(**parameters)

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            count_and_serve,
        ):
            snapshot_payload = self._resolve_snapshot(max_retry_attempts=1)

        self.assertEqual(efetch_call_count["value"], 1)
        self.assertEqual(snapshot_payload["manifest"]["reused_batch_count"], 2)
        self.assertEqual(snapshot_payload["manifest"]["fetched_batch_count"], 1)
        self.assertFalse(workspace_directory.exists())

    def test_published_snapshot_reloads_and_revalidates(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            self._fake_efetch,
        ):
            created_snapshot_payload = self._resolve_snapshot()

        reused_snapshot_payload = resolve_ncbi_protein_xml_snapshot(
            snapshot_mode="reuse_latest",
            snapshot_root_directory=self.xml_snapshot_root_directory,
            source_uid_snapshot_root_directory=self.uid_snapshot_root_directory,
        )

        self.assertEqual(
            reused_snapshot_payload["manifest"]["xml_file_sha256"],
            created_snapshot_payload["manifest"]["xml_file_sha256"],
        )
        self.assertEqual(
            reused_snapshot_payload["protein_uids"],
            self.protein_uids,
        )

    def test_unsupported_manifest_format_version_is_rejected(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            self._fake_efetch,
        ):
            self._resolve_snapshot()

        latest_manifest_file_path = (
            self.xml_snapshot_root_directory / "latest" / "manifest.json"
        )
        manifest_payload = json.loads(
            latest_manifest_file_path.read_text(encoding="utf-8")
        )
        manifest_payload["snapshot_format_version"] = "9.9"
        latest_manifest_file_path.write_text(
            json.dumps(manifest_payload, indent=2) + "\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(RuntimeError, "format version is not supported"):
            resolve_ncbi_protein_xml_snapshot(
                snapshot_mode="reuse_latest",
                snapshot_root_directory=self.xml_snapshot_root_directory,
                source_uid_snapshot_root_directory=self.uid_snapshot_root_directory,
            )

    def test_previous_manifest_format_version_is_still_accepted(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            self._fake_efetch,
        ):
            self._resolve_snapshot()

        latest_manifest_file_path = (
            self.xml_snapshot_root_directory / "latest" / "manifest.json"
        )
        manifest_payload = json.loads(
            latest_manifest_file_path.read_text(encoding="utf-8")
        )
        manifest_payload["snapshot_format_version"] = "1.0"
        latest_manifest_file_path.write_text(
            json.dumps(manifest_payload, indent=2) + "\n",
            encoding="utf-8",
        )

        reused_snapshot_payload = resolve_ncbi_protein_xml_snapshot(
            snapshot_mode="reuse_latest",
            snapshot_root_directory=self.xml_snapshot_root_directory,
            source_uid_snapshot_root_directory=self.uid_snapshot_root_directory,
        )

        self.assertEqual(
            reused_snapshot_payload["manifest"]["snapshot_format_version"],
            "1.0",
        )

    def test_concurrent_retrieval_publishes_the_same_bytes_as_sequential(
        self,
    ) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            self._fake_efetch,
        ):
            sequential_snapshot_payload = self._resolve_snapshot()

        sequential_xml_sha256 = sequential_snapshot_payload["manifest"][
            "xml_file_sha256"
        ]

        second_snapshot_root_directory = self.working_directory / "xml_snapshots_2"
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            self._fake_efetch,
        ):
            concurrent_snapshot_payload = self._resolve_snapshot(
                snapshot_root_directory=second_snapshot_root_directory,
                max_concurrent_requests=3,
                max_request_starts_per_second=1000.0,
            )

        self.assertEqual(
            concurrent_snapshot_payload["manifest"]["xml_file_sha256"],
            sequential_xml_sha256,
        )


class LatestDirectoryPublicationTests(unittest.TestCase):
    """
    Publishing latest/ must never leave the reader with nothing.

    Deleting the previous directory before the replacement is in place turns a
    briefly locked artifact into a destroyed latest/: the manifest goes, the
    payload stays, and no reader can validate what is left.
    """

    def setUp(self) -> None:
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        self.working_directory = Path(temporary_directory.name)

        self.source_directory = self.working_directory / "source"
        self.source_directory.mkdir()
        (self.source_directory / "manifest.json").write_text(
            '{"version": "new"}',
            encoding="utf-8",
        )
        (self.source_directory / "protein_records.xml").write_text(
            "<GBSet />",
            encoding="utf-8",
        )

        self.latest_directory = self.working_directory / "latest"
        self.files_to_copy = [
            (self.source_directory / "manifest.json", "manifest.json"),
            (self.source_directory / "protein_records.xml", "protein_records.xml"),
        ]

    def _write_existing_latest(self) -> None:
        self.latest_directory.mkdir()
        (self.latest_directory / "manifest.json").write_text(
            '{"version": "old"}',
            encoding="utf-8",
        )
        (self.latest_directory / "protein_records.xml").write_text(
            "<GBSet><!-- old --></GBSet>",
            encoding="utf-8",
        )

    def test_publication_replaces_an_existing_latest_directory(self) -> None:
        self._write_existing_latest()

        ncbi_snapshot._replace_latest_directory(
            latest_directory=self.latest_directory,
            files_to_copy=self.files_to_copy,
        )

        self.assertEqual(
            (self.latest_directory / "manifest.json").read_text(encoding="utf-8"),
            '{"version": "new"}',
        )
        self.assertEqual(
            sorted(path.name for path in self.working_directory.iterdir()),
            ["latest", "source"],
        )

    def test_existing_latest_survives_when_it_cannot_be_moved_aside(self) -> None:
        self._write_existing_latest()

        with patch.object(
            ncbi_snapshot,
            "_rename_directory_with_retries",
            return_value=False,
        ):
            with self.assertRaisesRegex(RuntimeError, "move the existing latest"):
                ncbi_snapshot._replace_latest_directory(
                    latest_directory=self.latest_directory,
                    files_to_copy=self.files_to_copy,
                )

        # The old snapshot is still complete and readable.
        self.assertEqual(
            (self.latest_directory / "manifest.json").read_text(encoding="utf-8"),
            '{"version": "old"}',
        )
        self.assertTrue((self.latest_directory / "protein_records.xml").exists())
        self.assertEqual(
            sorted(path.name for path in self.working_directory.iterdir()),
            ["latest", "source"],
        )

    def test_existing_latest_is_restored_when_publication_cannot_finish(
        self,
    ) -> None:
        self._write_existing_latest()
        rename_call_outcomes = [True, False, True]

        def rename_with_scripted_outcomes(*, source_directory, destination_directory, **_):
            rename_succeeds = rename_call_outcomes.pop(0)
            if rename_succeeds:
                source_directory.replace(destination_directory)

            return rename_succeeds

        with patch.object(
            ncbi_snapshot,
            "_rename_directory_with_retries",
            rename_with_scripted_outcomes,
        ):
            with self.assertRaisesRegex(RuntimeError, "Unable to publish"):
                ncbi_snapshot._replace_latest_directory(
                    latest_directory=self.latest_directory,
                    files_to_copy=self.files_to_copy,
                )

        self.assertEqual(
            (self.latest_directory / "manifest.json").read_text(encoding="utf-8"),
            '{"version": "old"}',
        )

    def test_leftovers_from_an_interrupted_publication_are_swept(self) -> None:
        self._write_existing_latest()
        leftover_directory = (
            self.working_directory / "latest_superseded_abc123"
        )
        leftover_directory.mkdir()
        (leftover_directory / "protein_records.xml").write_text(
            "stale",
            encoding="utf-8",
        )

        ncbi_snapshot._replace_latest_directory(
            latest_directory=self.latest_directory,
            files_to_copy=self.files_to_copy,
        )

        self.assertFalse(leftover_directory.exists())
        self.assertEqual(
            sorted(path.name for path in self.working_directory.iterdir()),
            ["latest", "source"],
        )


if __name__ == "__main__":
    unittest.main()
