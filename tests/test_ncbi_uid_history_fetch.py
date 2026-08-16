import contextlib
import io
import sys
import types
import unittest
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
    NCBIXmlResponseValidationError,
    fetch_ncbi_protein_uid_snapshot,
)


class FakeFetchHandle:
    def __init__(self, payload: bytes) -> None:
        self.payload = payload
        self.closed = False
        self.url = "https://ncbi.example.test/efetch"
        self.headers = {"content-type": "text/plain"}

    def read(self):
        return self.payload

    def close(self) -> None:
        self.closed = True


def build_uid_page_payload(protein_uids: list[str]) -> bytes:
    return ("\n".join(protein_uids) + "\n").encode("utf-8")


class ProteinUidHistoryRetrievalTests(unittest.TestCase):
    def setUp(self) -> None:
        self.initial_search_response = {
            "Count": "5",
            "QueryTranslation": "translated query",
            "WebEnv": "MCID_test_web_env",
            "QueryKey": "1",
        }

        entrez_read_patcher = patch.object(
            ncbi_api.Entrez,
            "read",
            return_value=self.initial_search_response,
            create=True,
        )
        entrez_read_patcher.start()
        self.addCleanup(entrez_read_patcher.stop)

        self.esearch_calls: list[dict] = []

        def fake_esearch(**parameters):
            self.esearch_calls.append(parameters)
            return FakeFetchHandle(b"<eSearchResult />")

        esearch_patcher = patch.object(
            ncbi_api.Entrez,
            "esearch",
            fake_esearch,
            create=True,
        )
        esearch_patcher.start()
        self.addCleanup(esearch_patcher.stop)

        sleep_patcher = patch("src.pago_pipeline.ncbi_api.time.sleep")
        self.sleep_mock = sleep_patcher.start()
        self.addCleanup(sleep_patcher.stop)

        jitter_patcher = patch(
            "src.pago_pipeline.ncbi_api.random.uniform",
            return_value=0.0,
        )
        jitter_patcher.start()
        self.addCleanup(jitter_patcher.stop)

    def test_query_is_executed_once_and_pages_read_from_the_history_set(self) -> None:
        efetch_calls: list[dict] = []

        def fake_efetch(**parameters):
            efetch_calls.append(parameters)
            page_start_index = parameters["retstart"]
            all_protein_uids = ["1005", "1004", "1003", "1002", "1001"]
            page_size = parameters["retmax"]
            return FakeFetchHandle(
                build_uid_page_payload(
                    all_protein_uids[
                        page_start_index: page_start_index + page_size
                    ]
                )
            )

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            fake_efetch,
        ):
            fetch_result = fetch_ncbi_protein_uid_snapshot(
                ncbi_email="test@example.org",
                ncbi_api_key="test-key",
                query="PIWI[All Fields]",
                page_size=2,
                request_delay_seconds=0.0,
            )

        self.assertEqual(len(self.esearch_calls), 1)
        self.assertEqual(self.esearch_calls[0]["usehistory"], "y")
        self.assertEqual(self.esearch_calls[0]["retmax"], 0)

        self.assertEqual(len(efetch_calls), 3)
        for efetch_call in efetch_calls:
            self.assertEqual(efetch_call["rettype"], "uilist")
            self.assertEqual(efetch_call["retmode"], "text")
            self.assertEqual(efetch_call["history_web_env"], "MCID_test_web_env")
            self.assertEqual(efetch_call["history_query_key"], "1")
            self.assertNotIn("protein_uids", efetch_call)

        self.assertEqual(
            [efetch_call["retstart"] for efetch_call in efetch_calls],
            [0, 2, 4],
        )
        self.assertEqual(
            fetch_result.protein_uids,
            ["1001", "1002", "1003", "1004", "1005"],
        )
        self.assertEqual(fetch_result.normalized_protein_uid_count, 5)
        self.assertEqual(fetch_result.esearch_request_count, 1)
        self.assertEqual(fetch_result.efetch_request_count, 3)
        self.assertEqual(
            fetch_result.uid_retrieval_strategy,
            "esearch_history_efetch_uilist",
        )
        self.assertEqual(fetch_result.history_web_env, "MCID_test_web_env")

    def test_page_size_is_capped_at_the_efetch_paging_ceiling(self) -> None:
        with self.assertRaisesRegex(ValueError, "EFetch paging ceiling"):
            fetch_ncbi_protein_uid_snapshot(
                ncbi_email="test@example.org",
                ncbi_api_key=None,
                query="PIWI[All Fields]",
                page_size=10_001,
            )

    def test_uid_page_transient_failure_is_retried_under_the_shared_controls(
        self,
    ) -> None:
        page_payload = build_uid_page_payload(["1001", "1002", "1003", "1004", "1005"])

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[
                URLError("temporarily unavailable"),
                FakeFetchHandle(page_payload),
            ],
        ) as efetch_mock:
            fetch_result = fetch_ncbi_protein_uid_snapshot(
                ncbi_email="test@example.org",
                ncbi_api_key=None,
                query="PIWI[All Fields]",
                page_size=10,
                max_retry_attempts=2,
                request_delay_seconds=0.0,
                retry_backoff_initial_seconds=0.25,
            )

        self.assertEqual(efetch_mock.call_count, 2)
        self.assertEqual(fetch_result.normalized_protein_uid_count, 5)
        self.sleep_mock.assert_any_call(0.25)

    def test_uid_page_retry_limit_raises_an_explicit_error(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=URLError("temporarily unavailable"),
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "Failed to extract protein UIDs after 2 attempts",
            ):
                fetch_ncbi_protein_uid_snapshot(
                    ncbi_email="test@example.org",
                    ncbi_api_key=None,
                    query="PIWI[All Fields]",
                    page_size=10,
                    max_retry_attempts=2,
                    request_delay_seconds=0.0,
                    retry_backoff_initial_seconds=0.25,
                )

    def test_non_numeric_uid_page_is_a_permanent_validation_failure(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            return_value=FakeFetchHandle(b"1001\nnot-a-uid\n"),
        ) as efetch_mock:
            with self.assertRaisesRegex(
                RuntimeError,
                "Permanent NCBI protein UID page failure",
            ) as raised_error:
                fetch_ncbi_protein_uid_snapshot(
                    ncbi_email="test@example.org",
                    ncbi_api_key=None,
                    query="PIWI[All Fields]",
                    page_size=10,
                    request_delay_seconds=0.0,
                )

        self.assertIsInstance(
            raised_error.exception.__cause__,
            NCBIXmlResponseValidationError,
        )
        efetch_mock.assert_called_once()

    def test_uid_page_larger_than_requested_is_rejected(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            return_value=FakeFetchHandle(
                build_uid_page_payload(["1001", "1002", "1003"])
            ),
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "more protein UIDs than requested",
            ):
                fetch_ncbi_protein_uid_snapshot(
                    ncbi_email="test@example.org",
                    ncbi_api_key=None,
                    query="PIWI[All Fields]",
                    page_size=2,
                    request_delay_seconds=0.0,
                )

    def test_uid_pages_report_progress_against_the_reported_result_count(
        self,
    ) -> None:
        self.initial_search_response["Count"] = "5"

        def fake_efetch(**parameters):
            page_start_index = parameters["retstart"]
            all_protein_uids = ["1001", "1002", "1003", "1004", "1005"]
            return FakeFetchHandle(
                build_uid_page_payload(
                    all_protein_uids[
                        page_start_index: page_start_index + parameters["retmax"]
                    ]
                )
            )

        captured_output = io.StringIO()
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            fake_efetch,
        ):
            with contextlib.redirect_stdout(captured_output):
                fetch_ncbi_protein_uid_snapshot(
                    ncbi_email="test@example.org",
                    ncbi_api_key=None,
                    query="PIWI[All Fields]",
                    page_size=2,
                    request_delay_seconds=0.0,
                )

        retrieval_log = captured_output.getvalue()

        # Progress advances by identifiers retrieved, not by page count, so the
        # bar tracks the quantity the reader actually cares about.
        self.assertIn("2/5", retrieval_log)
        self.assertIn("4/5", retrieval_log)
        self.assertIn("5/5 100.0%", retrieval_log)
        self.assertIn("page 1", retrieval_log)
        self.assertIn("page 3", retrieval_log)
        self.assertRegex(
            retrieval_log,
            r"Retrieved 5 protein UIDs in [\d.]+s\.",
        )
        self.assertNotIn("Extracted 5 protein UIDs", retrieval_log)

    def test_uid_retrieval_records_stage_telemetry(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            return_value=FakeFetchHandle(
                build_uid_page_payload(
                    ["1001", "1002", "1003", "1004", "1005"]
                )
            ),
        ):
            fetch_result = fetch_ncbi_protein_uid_snapshot(
                ncbi_email="test@example.org",
                ncbi_api_key=None,
                query="PIWI[All Fields]",
                page_size=10,
                request_delay_seconds=0.0,
            )

        telemetry_stages = fetch_result.retrieval_telemetry["stages"]
        self.assertEqual(telemetry_stages["uid_initial_search"]["request_count"], 1)
        self.assertEqual(telemetry_stages["uid_pages"]["request_count"], 1)
        self.assertGreater(
            telemetry_stages["uid_pages"]["response_byte_total"],
            0,
        )


if __name__ == "__main__":
    unittest.main()
