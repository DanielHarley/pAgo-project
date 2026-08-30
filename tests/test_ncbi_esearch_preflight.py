from __future__ import annotations

import unittest
from unittest import mock

from src.pago_pipeline.ncbi_esearch_preflight import (
    build_esearch_preflight_result,
    parse_uilist_text,
    run_ncbi_esearch_preflight,
    summarize_sample_xml,
)

_SAMPLE_GBSET_XML = b"""<?xml version="1.0"?>
<GBSet>
  <GBSeq>
    <GBSeq_other-seqids><GBSeqid>gi|100</GBSeqid></GBSeq_other-seqids>
    <GBSeq_sequence>mkstl</GBSeq_sequence>
  </GBSeq>
  <GBSeq>
    <GBSeq_other-seqids><GBSeqid>gi|101</GBSeqid></GBSeq_other-seqids>
  </GBSeq>
</GBSet>
"""


class EsearchPreflightLogicTests(unittest.TestCase):
    def test_parse_uilist_text(self) -> None:
        self.assertEqual(
            parse_uilist_text(payload_text="100\n101\n\n 102 \n"),
            ["100", "101", "102"],
        )

    def test_summarize_sample_xml(self) -> None:
        summary = summarize_sample_xml(xml_payload_bytes=_SAMPLE_GBSET_XML)
        self.assertEqual(summary["sample_record_count"], 2)
        self.assertEqual(summary["sample_records_with_sequence"], 1)
        self.assertEqual(summary["sample_records_missing_sequence"], 1)
        self.assertEqual(summary["sample_records_with_extractable_uid"], 2)

    def test_build_result_flags_exceeds_max_uid_count(self) -> None:
        below = build_esearch_preflight_result(
            search_query="q",
            retrieved_at_utc="2026-08-30T00:00:00Z",
            result_count=100,
            translated_query="q[All Fields]",
            history_web_env="W",
            history_query_key="1",
            max_uid_count=150,
            sample_requested_count=0,
        )
        self.assertFalse(below.exceeds_max_uid_count)

        above = build_esearch_preflight_result(
            search_query="q",
            retrieved_at_utc="2026-08-30T00:00:00Z",
            result_count=200,
            translated_query="q[All Fields]",
            history_web_env="W",
            history_query_key="1",
            max_uid_count=150,
            sample_requested_count=0,
        )
        self.assertTrue(above.exceeds_max_uid_count)

    def test_run_preflight_uses_history_and_sample(self) -> None:
        esearch_handle = mock.MagicMock()
        esearch_handle.read.return_value = (
            b"<?xml version='1.0'?>"
            b"<eSearchResult><Count>2</Count>"
            b"<QueryTranslation>piwi[All Fields]</QueryTranslation>"
            b"<IdList/><WebEnv>WEBENV</WebEnv><QueryKey>1</QueryKey>"
            b"</eSearchResult>"
        )
        uilist_handle = mock.MagicMock()
        uilist_handle.read.return_value = b"100\n101\n"
        xml_handle = mock.MagicMock()
        xml_handle.read.return_value = _SAMPLE_GBSET_XML

        with mock.patch(
            "src.pago_pipeline.ncbi_esearch_preflight.Entrez"
        ) as entrez_mock, mock.patch(
            "src.pago_pipeline.ncbi_esearch_preflight._configured_ncbi_entrez_urlopen"
        ) as ssl_ctx_mock:
            ssl_ctx_mock.return_value.__enter__.return_value = (None, None)
            entrez_mock.esearch.return_value = esearch_handle
            entrez_mock.efetch.side_effect = [uilist_handle, xml_handle]
            entrez_mock.read.side_effect = lambda handle: {
                "Count": "2",
                "QueryTranslation": "piwi[All Fields]",
                "WebEnv": "WEBENV",
                "QueryKey": "1",
            }

            result = run_ncbi_esearch_preflight(
                search_query="(PIWI[All Fields] OR Argonaute[All Fields])",
                ncbi_email="tester@example.com",
                max_uid_count=1,
                sample_size=5,
            )

        self.assertEqual(result.result_count, 2)
        self.assertEqual(result.translated_query, "piwi[All Fields]")
        self.assertTrue(result.exceeds_max_uid_count)
        self.assertEqual(result.sample_uid_list, ["100", "101"])
        self.assertEqual(result.sample_record_count, 2)
        self.assertEqual(result.sample_records_with_sequence, 1)
        self.assertIsNone(result.sample_fetch_error)


if __name__ == "__main__":
    unittest.main()
