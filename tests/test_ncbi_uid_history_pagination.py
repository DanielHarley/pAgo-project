
import unittest
from unittest.mock import patch

from src.pago_pipeline.ncbi_api import (
    fetch_ncbi_protein_uid_snapshot,
)

class TestNCBIHistoryPagination(unittest.TestCase):

    @patch("src.pago_pipeline.ncbi_api.time.sleep")
    @patch("src.pago_pipeline.ncbi_api.Entrez.read")
    @patch("src.pago_pipeline.ncbi_api.Entrez.esearch")
    def test_all_pages_use_original_history(
        self,
        mock_esearch,
        mock_read,
        mock_sleep,
    ):
        query = "PIWI[All Fields] AND Bacteria[Organism]"

        # Simulate 2,500 distinct protein UIDs.
        expected_uids = [
            str(uid) for uid in range(1, 2501)
        ]

        # One initial search followed by three pages.
        mock_read.side_effect = [
            {
                "Count": "2500",
                "WebEnv": "TEST_WEBENV",
                "QueryKey": "1",
                "QueryTranslation": query,
            },
            {"IdList": expected_uids[:1000]},
            {"IdList": expected_uids[1000:2000]},
            {"IdList": expected_uids[2000:]},
        ]

        result = fetch_ncbi_protein_uid_snapshot(
            ncbi_email="test@example.org",
            ncbi_api_key=None,
            query=query,
            page_size=1000,
            max_retry_attempts=1,
            request_delay_seconds=0.0,
        )

        calls = mock_esearch.call_args_list

        # One initial search + three paginated requests.
        self.assertEqual(len(calls), 4)

        # The initial search uses the original query.
        self.assertEqual(calls[0].kwargs["term"], query)
        self.assertEqual(calls[0].kwargs["retmax"], 0)
        self.assertEqual(calls[0].kwargs["usehistory"], "y")

        # Every page must reference the original history.
        for index, call in enumerate(calls[1:]):
            self.assertEqual(
                call.kwargs.get("term"),
                "#1",
            )
            self.assertEqual(
                call.kwargs.get("WebEnv"),
                "TEST_WEBENV",
            )
            self.assertEqual(
                call.kwargs["retstart"],
                index * 1000,
            )
            self.assertEqual(
                call.kwargs["retmax"],
                1000,
            )

        # Verify the complete result was collected.
        self.assertEqual(
            result.ncbi_reported_result_count,
            2500,
        )
        self.assertEqual(
            result.normalized_protein_uid_count,
            2500,
        )
        self.assertEqual(
            set(result.protein_uids),
            set(expected_uids),
        )

    @patch("src.pago_pipeline.ncbi_api.time.sleep")
    @patch("src.pago_pipeline.ncbi_api.Entrez.read")
    @patch("src.pago_pipeline.ncbi_api.Entrez.esearch")
    def test_history_failure_does_not_repeat_original_query(
        self,
        mock_esearch,
        mock_read,
        mock_sleep,
    ):
        query = "PIWI[All Fields] AND Bacteria[Organism]"

        # Initial search succeeds and creates a History Server result.
        initial_response = {
            "Count": "4",
            "WebEnv": "TEST_WEBENV",
            "QueryKey": "1",
            "QueryTranslation": query,
        }

        # Both attempts to retrieve the first page fail.
        mock_read.side_effect = [
            initial_response,
            RuntimeError("History Server session expired"),
            RuntimeError("History Server session expired"),
        ]

        # The function must report the failure.
        with self.assertRaisesRegex(
            RuntimeError,
            "History Server session expired",
        ):
            fetch_ncbi_protein_uid_snapshot(
                ncbi_email="test@example.org",
                ncbi_api_key=None,
                query=query,
                page_size=2,
                max_retry_attempts=2,
                request_delay_seconds=0.0,
            )

        calls = mock_esearch.call_args_list

        # One initial search and two failed pagination attempts.
        self.assertEqual(len(calls), 3)

        # The initial request uses the original query.
        self.assertEqual(calls[0].kwargs["term"], query)

        # Both retries must reference the original History Server set.
        for call in calls[1:]:
            self.assertEqual(
                call.kwargs["term"],
                "#1",
            )

            self.assertEqual(
                call.kwargs["WebEnv"],
                "TEST_WEBENV",
            )

    @patch("src.pago_pipeline.ncbi_api.time.sleep")
    @patch("src.pago_pipeline.ncbi_api.Entrez.read")
    @patch("src.pago_pipeline.ncbi_api.Entrez.esearch")
    def test_reject_incomplete_normalized_uid_list(
            self,
            mock_esearch,
            mock_read,
            mock_sleep,
    ):
        query = "PIWI[All Fields] AND Bacteria[Organism]"

        # NCBI initially reports four protein UIDs.
        initial_response = {
            "Count": "4",
            "WebEnv": "TEST_WEBENV",
            "QueryKey": "1",
            "QueryTranslation": query,
        }

        # Four entries are retrieved, but one UID is duplicated.
        # Therefore, only three unique UIDs are present.
        mock_read.side_effect = [
            initial_response,
            {"IdList": ["101", "102"]},
            {"IdList": ["102", "104"]},
        ]

        # An incomplete normalized result must be rejected.
        with self.assertRaisesRegex(
                RuntimeError,
                "UID completeness validation failed",
        ):
            fetch_ncbi_protein_uid_snapshot(
                ncbi_email="test@example.org",
                ncbi_api_key=None,
                query=query,
                page_size=2,
                max_retry_attempts=1,
                request_delay_seconds=0.0,
                deduplicate_uids=True,
                sort_uids=True,
            )

    @patch("src.pago_pipeline.ncbi_api.time.sleep")
    @patch("src.pago_pipeline.ncbi_api.Entrez.read")
    @patch("src.pago_pipeline.ncbi_api.Entrez.esearch")
    def test_accept_complete_unique_uid_list(
        self,
        mock_esearch,
        mock_read,
        mock_sleep,
    ):
        query = "PIWI[All Fields] AND Bacteria[Organism]"

        mock_read.side_effect = [
            {
                "Count": "4",
                "WebEnv": "TEST_WEBENV",
                "QueryKey": "1",
                "QueryTranslation": query,
            },
            {"IdList": ["101", "102"]},
            {"IdList": ["103", "104"]},
        ]

        result = fetch_ncbi_protein_uid_snapshot(
            ncbi_email="test@example.org",
            ncbi_api_key=None,
            query=query,
            page_size=2,
            max_retry_attempts=1,
            request_delay_seconds=0.0,
            deduplicate_uids=True,
            sort_uids=True,
        )

        self.assertEqual(
            result.ncbi_reported_result_count, 4
        )
        self.assertEqual(
            result.raw_protein_uid_count, 4
        )
        self.assertEqual(
            result.normalized_protein_uid_count, 4
        )
        self.assertEqual(
            result.protein_uids,
            ["101", "102", "103", "104"],
        )

    @patch("src.pago_pipeline.ncbi_api.time.sleep")
    @patch("src.pago_pipeline.ncbi_api.Entrez.read")
    @patch("src.pago_pipeline.ncbi_api.Entrez.esearch")
    def test_reject_duplicates_when_deduplication_disabled(
        self,
        mock_esearch,
        mock_read,
        mock_sleep,
    ):
        query = "PIWI[All Fields] AND Bacteria[Organism]"

        mock_read.side_effect = [
            {
                "Count": "4",
                "WebEnv": "TEST_WEBENV",
                "QueryKey": "1",
                "QueryTranslation": query,
            },
            {"IdList": ["101", "102"]},
            {"IdList": ["102", "104"]},
        ]

        with self.assertRaisesRegex(
            RuntimeError,
            "UID completeness validation failed",
        ):
            fetch_ncbi_protein_uid_snapshot(
                ncbi_email="test@example.org",
                ncbi_api_key=None,
                query=query,
                page_size=2,
                max_retry_attempts=1,
                request_delay_seconds=0.0,
                deduplicate_uids=False,
                sort_uids=True,
            )

    @patch("src.pago_pipeline.ncbi_api.time.sleep")
    @patch("src.pago_pipeline.ncbi_api.Entrez.read")
    @patch("src.pago_pipeline.ncbi_api.Entrez.esearch")
    def test_reject_raw_uid_count_exceeding_ncbi_count(
                self,
                mock_esearch,
                mock_read,
                mock_sleep,
        ):
            query = "PIWI[All Fields] AND Bacteria[Organism]"

            # NCBI initially reports four protein UIDs.
            initial_response = {
                "Count": "4",
                "WebEnv": "TEST_WEBENV",
                "QueryKey": "1",
                "QueryTranslation": query,
            }

            # Simulate an anomalous second page returning
            # three UIDs despite a requested page size of two.
            mock_read.side_effect = [
                initial_response,
                {"IdList": ["101", "102"]},
                {"IdList": ["103", "104", "105"]},
            ]

            # Five unique UIDs must not be accepted when
            # the original NCBI result count was four.
            with self.assertRaisesRegex(
                    RuntimeError,
                    "UID completeness validation failed",
            ):
                fetch_ncbi_protein_uid_snapshot(
                    ncbi_email="test@example.org",
                    ncbi_api_key=None,
                    query=query,
                    page_size=2,
                    max_retry_attempts=1,
                    request_delay_seconds=0.0,
                    deduplicate_uids=True,
                    sort_uids=True,
                )

if __name__ == "__main__":
    unittest.main()