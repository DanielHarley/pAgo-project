from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from src.pago_pipeline.ncbi_esearch_preflight import EsearchPreflightResult
from src.pago_pipeline.ncbi_esearch_preflight_snapshot import (
    ARTIFACT_TYPE,
    latest_ncbi_esearch_preflight_snapshot_is_available,
    resolve_ncbi_esearch_preflight_snapshot,
)
from src.pago_pipeline.ncbi_snapshot import SnapshotMode

_QUERY = "(PIWI[All Fields] OR Argonaute[All Fields]) AND (Bacteria[Organism])"


def _fake_result(*, result_count: int, max_uid_count: int) -> EsearchPreflightResult:
    return EsearchPreflightResult(
        search_query=_QUERY,
        translated_query="piwi[All Fields] OR argonaute[All Fields]",
        result_count=result_count,
        history_web_env="WEBENV",
        history_query_key="1",
        retrieved_at_utc="2026-08-30T12:00:00Z",
        max_uid_count=max_uid_count,
        exceeds_max_uid_count=result_count > max_uid_count,
        sample_requested_count=2,
        sample_uid_list=["100", "101"],
        sample_record_count=2,
        sample_records_with_sequence=2,
        sample_records_missing_sequence=0,
        sample_records_with_extractable_uid=2,
        python_version="3.14.0",
        biopython_version="1.85",
    )


class NCBIEsearchPreflightSnapshotTests(unittest.TestCase):
    def test_resolve_materializes_report_and_reuses(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_root = Path(temporary_directory_name) / "preflight"
            with mock.patch(
                "src.pago_pipeline.ncbi_esearch_preflight_snapshot"
                ".run_ncbi_esearch_preflight",
                return_value=_fake_result(result_count=42, max_uid_count=150_000),
            ) as run_mock:
                payload = resolve_ncbi_esearch_preflight_snapshot(
                    snapshot_mode=SnapshotMode.reuse_latest_or_create,
                    snapshot_root_directory=snapshot_root,
                    search_query=_QUERY,
                    ncbi_email="tester@example.com",
                )
                run_mock.assert_called_once()

            manifest = payload["manifest"]
            self.assertEqual(manifest["artifact_type"], ARTIFACT_TYPE)
            self.assertEqual(manifest["result_count"], 42)
            self.assertFalse(manifest["exceeds_max_uid_count"])
            self.assertEqual(manifest["search_query"], _QUERY)
            self.assertEqual(
                Path(payload["sample_uid_file_path"]).read_text().split(),
                ["100", "101"],
            )

            self.assertTrue(
                latest_ncbi_esearch_preflight_snapshot_is_available(
                    snapshot_root_directory=snapshot_root,
                    search_query=_QUERY,
                )
            )
            self.assertFalse(
                latest_ncbi_esearch_preflight_snapshot_is_available(
                    snapshot_root_directory=snapshot_root,
                    search_query="a different query",
                )
            )

            with mock.patch(
                "src.pago_pipeline.ncbi_esearch_preflight_snapshot"
                ".run_ncbi_esearch_preflight"
            ) as run_mock_2:
                resolve_ncbi_esearch_preflight_snapshot(
                    snapshot_mode=SnapshotMode.reuse_latest_or_create,
                    snapshot_root_directory=snapshot_root,
                    search_query=_QUERY,
                    ncbi_email="tester@example.com",
                )
                run_mock_2.assert_not_called()

    def test_resolve_raises_when_result_count_exceeds_max_uid_count(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_root = Path(temporary_directory_name) / "preflight"
            with mock.patch(
                "src.pago_pipeline.ncbi_esearch_preflight_snapshot"
                ".run_ncbi_esearch_preflight",
                return_value=_fake_result(result_count=400_000, max_uid_count=150_000),
            ):
                with self.assertRaisesRegex(RuntimeError, "above the configured"):
                    resolve_ncbi_esearch_preflight_snapshot(
                        snapshot_mode=SnapshotMode.reuse_latest_or_create,
                        snapshot_root_directory=snapshot_root,
                        search_query=_QUERY,
                        ncbi_email="tester@example.com",
                        max_uid_count=150_000,
                    )
            # The report is still materialized for audit despite the raise.
            self.assertTrue(
                latest_ncbi_esearch_preflight_snapshot_is_available(
                    snapshot_root_directory=snapshot_root,
                    search_query=_QUERY,
                )
            )
            # Explicit opt-in returns the payload instead of raising.
            payload = resolve_ncbi_esearch_preflight_snapshot(
                snapshot_mode=SnapshotMode.reuse_latest,
                snapshot_root_directory=snapshot_root,
                allow_exceeds_max_uid_count=True,
            )
            self.assertTrue(payload["manifest"]["exceeds_max_uid_count"])


if __name__ == "__main__":
    unittest.main()
