
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
if __name__ == "__main__":
    unittest.main()