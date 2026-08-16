import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from src.pago_pipeline.ncbi_xml_stream import (
    XmlBatchPayloadSource,
    build_legacy_consolidated_xml_payload,
    validate_consolidated_xml_file,
    validate_xml_batch_payload,
    write_consolidated_xml_document,
)


def build_gbseq_record(
    *,
    protein_uid: str,
    locus: str,
    trailing_text: str = "\n",
) -> str:
    return (
        f"<GBSeq><GBSeqid>gi|{protein_uid}</GBSeqid>"
        f"<GBSeq_locus>{locus}</GBSeq_locus></GBSeq>{trailing_text}"
    )


def build_gbset_batch_payload(
    *,
    records: list[str],
    include_declaration: bool = True,
    root_attributes: str = "",
) -> bytes:
    declaration = '<?xml version="1.0" encoding="UTF-8"?>\n' if include_declaration else ""
    return (
        f"{declaration}<GBSet{root_attributes}>\n"
        + "".join(records)
        + "</GBSet>\n"
    ).encode("utf-8")


class StreamingConsolidationTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        self.working_directory = Path(temporary_directory.name)

    def _assert_streaming_matches_legacy(self, batch_payloads: list[bytes]) -> None:
        legacy_payload, legacy_record_count = (
            build_legacy_consolidated_xml_payload(batch_payloads=batch_payloads)
        )

        output_file_path = self.working_directory / "protein_records.xml"
        write_result = write_consolidated_xml_document(
            batch_payload_sources=[
                XmlBatchPayloadSource(
                    batch_index=batch_index,
                    payload_bytes=batch_payload,
                )
                for batch_index, batch_payload in enumerate(batch_payloads, start=1)
            ],
            output_file_path=output_file_path,
        )

        streamed_payload = output_file_path.read_bytes()
        self.assertEqual(streamed_payload, legacy_payload)
        self.assertEqual(write_result.record_count, legacy_record_count)
        self.assertEqual(write_result.byte_count, len(legacy_payload))

    def test_streaming_output_is_byte_identical_to_in_memory_consolidation(
        self,
    ) -> None:
        self._assert_streaming_matches_legacy(
            [
                build_gbset_batch_payload(
                    records=[
                        build_gbseq_record(protein_uid="1001", locus="AAA"),
                        build_gbseq_record(protein_uid="1002", locus="BBB"),
                    ],
                ),
                build_gbset_batch_payload(
                    records=[
                        build_gbseq_record(protein_uid="1003", locus="CCC"),
                    ],
                    include_declaration=False,
                ),
            ]
        )

    def test_streaming_output_preserves_entities_unicode_and_empty_elements(
        self,
    ) -> None:
        self._assert_streaming_matches_legacy(
            [
                build_gbset_batch_payload(
                    records=[
                        "<GBSeq><GBSeqid>gi|2001</GBSeqid>"
                        "<GBSeq_definition>alpha &amp; beta &lt;core&gt;</GBSeq_definition>"
                        "<GBSeq_comment/>"
                        "<GBSeq_organism>Pseudomonas æruginosa</GBSeq_organism>"
                        '<GBSeq_other value="x&quot;y"/>'
                        "</GBSeq>\n",
                    ],
                ),
                build_gbset_batch_payload(
                    records=[
                        build_gbseq_record(
                            protein_uid="2002",
                            locus="DDD",
                            trailing_text="",
                        ),
                    ],
                ),
            ]
        )

    def test_streaming_output_preserves_root_attributes(self) -> None:
        self._assert_streaming_matches_legacy(
            [
                build_gbset_batch_payload(
                    records=[build_gbseq_record(protein_uid="3001", locus="EEE")],
                    root_attributes=' schema="gbset" release="2026"',
                ),
                build_gbset_batch_payload(
                    records=[build_gbseq_record(protein_uid="3002", locus="FFF")],
                    root_attributes=' schema="gbset" release="2026"',
                ),
            ]
        )

    def test_streaming_consolidation_reads_batches_from_disk(self) -> None:
        batch_payloads = [
            build_gbset_batch_payload(
                records=[build_gbseq_record(protein_uid="4001", locus="GGG")],
            ),
            build_gbset_batch_payload(
                records=[build_gbseq_record(protein_uid="4002", locus="HHH")],
            ),
        ]
        batch_payload_sources = []
        for batch_index, batch_payload in enumerate(batch_payloads, start=1):
            batch_file_path = self.working_directory / f"batch_{batch_index}.xml"
            batch_file_path.write_bytes(batch_payload)
            batch_payload_sources.append(
                XmlBatchPayloadSource(
                    batch_index=batch_index,
                    payload_file_path=batch_file_path,
                )
            )

        output_file_path = self.working_directory / "from_disk.xml"
        write_result = write_consolidated_xml_document(
            batch_payload_sources=batch_payload_sources,
            output_file_path=output_file_path,
        )

        legacy_payload, _ = build_legacy_consolidated_xml_payload(
            batch_payloads=batch_payloads,
        )
        self.assertEqual(output_file_path.read_bytes(), legacy_payload)
        self.assertEqual(write_result.record_count, 2)

    def test_streaming_consolidation_hash_matches_written_bytes(self) -> None:
        import hashlib

        output_file_path = self.working_directory / "hashed.xml"
        write_result = write_consolidated_xml_document(
            batch_payload_sources=[
                XmlBatchPayloadSource(
                    batch_index=1,
                    payload_bytes=build_gbset_batch_payload(
                        records=[
                            build_gbseq_record(protein_uid="5001", locus="III"),
                        ],
                    ),
                )
            ],
            output_file_path=output_file_path,
        )

        self.assertEqual(
            write_result.sha256,
            hashlib.sha256(output_file_path.read_bytes()).hexdigest(),
        )

    def test_streaming_consolidation_rejects_root_tag_mismatch(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "root tag mismatch"):
            write_consolidated_xml_document(
                batch_payload_sources=[
                    XmlBatchPayloadSource(
                        batch_index=1,
                        payload_bytes=build_gbset_batch_payload(
                            records=[
                                build_gbseq_record(
                                    protein_uid="6001",
                                    locus="JJJ",
                                ),
                            ],
                        ),
                    ),
                    XmlBatchPayloadSource(
                        batch_index=2,
                        payload_bytes=b"<NotGBSet><GBSeq /></NotGBSet>",
                    ),
                ],
                output_file_path=self.working_directory / "mismatch.xml",
            )

    def test_streaming_consolidation_removes_partial_output_on_parse_failure(
        self,
    ) -> None:
        output_file_path = self.working_directory / "partial.xml"

        with self.assertRaisesRegex(RuntimeError, "Failed to parse XML batch 1"):
            write_consolidated_xml_document(
                batch_payload_sources=[
                    XmlBatchPayloadSource(
                        batch_index=1,
                        payload_bytes=b"<GBSet><GBSeq>",
                    )
                ],
                output_file_path=output_file_path,
            )

        self.assertFalse(output_file_path.exists())
        self.assertEqual(list(self.working_directory.iterdir()), [])


class StreamingValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        self.working_directory = Path(temporary_directory.name)

    def _write_consolidated_file(self, payload: bytes) -> Path:
        xml_file_path = self.working_directory / "protein_records.xml"
        xml_file_path.write_bytes(payload)
        return xml_file_path

    def test_validation_returns_records_and_uids_in_document_order(self) -> None:
        xml_file_path = self._write_consolidated_file(
            build_gbset_batch_payload(
                records=[
                    build_gbseq_record(protein_uid="1001", locus="AAA"),
                    build_gbseq_record(protein_uid="1002", locus="BBB"),
                ],
            )
        )

        validation_result = validate_consolidated_xml_file(
            xml_file_path=xml_file_path,
            expected_record_count=2,
            expected_protein_uids=["1001", "1002"],
        )

        self.assertEqual(validation_result.record_count, 2)
        self.assertEqual(validation_result.protein_uids, ["1001", "1002"])

    def test_validation_reports_missing_and_unexpected_uids(self) -> None:
        xml_file_path = self._write_consolidated_file(
            build_gbset_batch_payload(
                records=[build_gbseq_record(protein_uid="9999", locus="AAA")],
            )
        )

        with self.assertRaisesRegex(RuntimeError, "record UIDs do not match"):
            validate_consolidated_xml_file(
                xml_file_path=xml_file_path,
                expected_protein_uids=["1001"],
            )

    def test_validation_rejects_malformed_document(self) -> None:
        xml_file_path = self._write_consolidated_file(b"<GBSet><GBSeq>")

        with self.assertRaisesRegex(RuntimeError, "not well-formed"):
            validate_consolidated_xml_file(xml_file_path=xml_file_path)


class XmlBatchPayloadValidationTests(unittest.TestCase):
    def test_batch_validation_accepts_exact_uid_set(self) -> None:
        validation_result = validate_xml_batch_payload(
            xml_payload_bytes=build_gbset_batch_payload(
                records=[
                    build_gbseq_record(protein_uid="1001", locus="AAA"),
                    build_gbseq_record(protein_uid="1002", locus="BBB"),
                ],
            ),
            expected_protein_uids=["1001", "1002"],
        )

        self.assertEqual(validation_result.record_count, 2)
        self.assertEqual(validation_result.root_tag, "GBSet")

    def test_batch_validation_detects_substituted_uid(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "record UIDs do not match"):
            validate_xml_batch_payload(
                xml_payload_bytes=build_gbset_batch_payload(
                    records=[
                        build_gbseq_record(protein_uid="1001", locus="AAA"),
                        build_gbseq_record(protein_uid="9999", locus="BBB"),
                    ],
                ),
                expected_protein_uids=["1001", "1002"],
            )

    def test_batch_validation_can_skip_uid_extraction(self) -> None:
        validation_result = validate_xml_batch_payload(
            xml_payload_bytes=b"<GBSet><GBSeq /><GBSeq /></GBSet>",
            expected_protein_uids=["1001", "1002"],
            validate_protein_uids=False,
        )

        self.assertEqual(validation_result.record_count, 2)
        self.assertEqual(validation_result.protein_uids, [])

    def test_batch_validation_rejects_unexpected_record_tag(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "unexpected child tags"):
            validate_xml_batch_payload(
                xml_payload_bytes=b"<GBSet><NotGBSeq /></GBSet>",
                expected_protein_uids=["1001"],
            )

    def test_batch_validation_propagates_parse_errors(self) -> None:
        with self.assertRaises(ET.ParseError):
            validate_xml_batch_payload(
                xml_payload_bytes=b"<GBSet><GBSeq>",
                expected_protein_uids=["1001"],
            )


if __name__ == "__main__":
    unittest.main()
