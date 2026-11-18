import builtins
import hashlib
import io
import json
import os
import socket
import sys
import tempfile
import threading
import time
import types
import unittest
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import Mock, call, patch
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, urlsplit

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
from src.pago_pipeline.ncbi_snapshot import (
    build_snapshot_directory_name,
    load_latest_xml_snapshot,
    load_snapshot_by_directory,
    load_xml_snapshot_by_directory,
    latest_snapshot_is_available,
    latest_xml_snapshot_is_available,
    get_latest_snapshot_protein_uids_path,
    get_snapshot_protein_uids_path,
    resolve_ncbi_protein_uid_snapshot,
    resolve_ncbi_protein_xml_snapshot,
    save_ncbi_protein_identifier_snapshot,
    save_ncbi_protein_identifier_snapshot_from_file,
    save_ncbi_protein_uid_snapshot,
    save_ncbi_protein_xml_snapshot,
    _validate_saved_consolidated_xml_snapshot,
    _validate_xml_record_identifiers,
    _validate_xml_record_uids,
    _build_xml_snapshot_manifest,
)
from src.pago_pipeline.ncbi_api import (
    NCBICircuitBreakerOpen,
    NCBIProteinIdentifierType,
    NCBIXmlBatchDeadlineExceeded,
    NCBIXmlBatchFetchError,
    NCBIXmlResponseValidationError,
    NCBIXmlCircuitBreaker,
    NCBIProteinUidFetchResult,
    NCBIProteinXmlBatchFetchResult,
    NCBIProteinXmlFetchResult,
    _configured_ncbi_entrez_urlopen,
    _normalize_protein_identifier_list,
    _normalize_protein_uid_list,
    coerce_ncbi_protein_identifier_type,
    _resolve_ncbi_ssl_ca_configuration,
    fetch_ncbi_protein_xml_batches,
    fetch_ncbi_protein_uid_snapshot,
    get_default_ncbi_xml_circuit_breaker,
    validate_ncbi_protein_identifiers,
)
from src.pago_pipeline.storage import sha256_of_file, sha256_of_lines


class FakeFetchHandle:
    def __init__(self, payload: bytes | str, read_callback=None) -> None:
        self.payload = payload
        self.read_callback = read_callback
        self.closed = False
        self.url = "https://ncbi.example.test/efetch"
        self.headers = {"content-type": "text/xml"}

    def read(self):
        if self.read_callback is not None:
            self.read_callback()
        return self.payload

    def close(self) -> None:
        self.closed = True


class ValidateSavedConsolidatedXmlSnapshotTests(unittest.TestCase):
    def _write_xml(self, xml_text: str) -> Path:
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)

        xml_file_path = Path(temporary_directory.name) / "protein_records.xml"
        xml_file_path.write_text(xml_text, encoding="utf-8")
        return xml_file_path

    def test_returns_record_count_for_valid_gbset_snapshot(self) -> None:
        xml_file_path = self._write_xml(
            """<?xml version="1.0" encoding="utf-8"?>
<GBSet>
  <GBSeq><GBSeq_locus>A</GBSeq_locus></GBSeq>
  <GBSeq><GBSeq_locus>B</GBSeq_locus></GBSeq>
</GBSet>
"""
        )

        record_count = _validate_saved_consolidated_xml_snapshot(
            xml_file_path=xml_file_path,
            expected_record_count=2,
        )

        self.assertEqual(record_count, 2)

    def test_raises_for_unexpected_root_tag(self) -> None:
        xml_file_path = self._write_xml(
            """<?xml version="1.0" encoding="utf-8"?>
<NotGBSet>
  <GBSeq><GBSeq_locus>A</GBSeq_locus></GBSeq>
</NotGBSet>
"""
        )

        with self.assertRaisesRegex(RuntimeError, "root tag mismatch"):
            _validate_saved_consolidated_xml_snapshot(
                xml_file_path=xml_file_path,
                expected_record_count=1,
            )

    def test_raises_for_unexpected_child_tag(self) -> None:
        xml_file_path = self._write_xml(
            """<?xml version="1.0" encoding="utf-8"?>
<GBSet>
  <NotGBSeq><value>A</value></NotGBSeq>
</GBSet>
"""
        )

        with self.assertRaisesRegex(RuntimeError, "unexpected child tags"):
            _validate_saved_consolidated_xml_snapshot(
                xml_file_path=xml_file_path,
                expected_record_count=1,
            )

    def test_raises_for_record_count_mismatch(self) -> None:
        xml_file_path = self._write_xml(
            """<?xml version="1.0" encoding="utf-8"?>
<GBSet>
  <GBSeq><GBSeq_locus>A</GBSeq_locus></GBSeq>
</GBSet>
"""
        )

        with self.assertRaisesRegex(RuntimeError, "record count mismatch"):
            _validate_saved_consolidated_xml_snapshot(
                xml_file_path=xml_file_path,
                expected_record_count=2,
            )


class ValidateXmlRecordIdentifiersTests(unittest.TestCase):
    def _write_xml(self, xml_text: str) -> Path:
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        xml_file_path = Path(temporary_directory.name) / "protein_records.xml"
        xml_file_path.write_text(xml_text, encoding="utf-8")
        return xml_file_path

    def test_uid_validation_matches_requested_identifier(self) -> None:
        xml_file_path = self._write_xml(
            """<GBSet>
<GBSeq>
<GBSeqid>gi|1001</GBSeqid>
</GBSeq>
</GBSet>"""
        )
        self.assertEqual(
            _validate_xml_record_identifiers(
                xml_file_path=xml_file_path,
                expected_protein_identifiers=["1001"],
                identifier_type="uid",
            ),
            ["1001"],
        )
        self.assertEqual(
            _validate_xml_record_uids(
                xml_file_path=xml_file_path,
                expected_protein_uids=["1001"],
            ),
            ["1001"],
        )

    def test_uid_validation_rejects_mismatched_identifier(self) -> None:
        xml_file_path = self._write_xml(
            """<GBSet>
<GBSeq>
<GBSeqid>gi|9999</GBSeqid>
</GBSeq>
</GBSet>"""
        )
        mismatch_message = (
            r"record UIDs do not match.*"
            r"Missing: \['1001'\]; Unexpected: \['9999'\]"
        )
        with self.assertRaisesRegex(RuntimeError, mismatch_message):
            _validate_xml_record_identifiers(
                xml_file_path=xml_file_path,
                expected_protein_identifiers=["1001"],
                identifier_type="uid",
            )
        with self.assertRaisesRegex(RuntimeError, mismatch_message):
            _validate_xml_record_uids(
                xml_file_path=xml_file_path,
                expected_protein_uids=["1001"],
            )

    def test_uid_validation_preserves_xml_order_and_duplicate_occurrences(self) -> None:
        xml_file_path = self._write_xml(
            "<GBSet>"
            "<GBSeq><GBSeqid>gi|1002</GBSeqid></GBSeq>"
            "<GBSeq><GBSeqid>gi|1001</GBSeqid></GBSeq>"
            "<GBSeq><GBSeqid>gi|1002</GBSeqid></GBSeq>"
            "</GBSet>"
        )
        expected_identifiers = ["1001", "1002", "1002"]
        for identifier_type in ("uid", NCBIProteinIdentifierType.uid):
            with self.subTest(identifier_type=identifier_type):
                self.assertEqual(
                    _validate_xml_record_identifiers(
                        xml_file_path=xml_file_path,
                        expected_protein_identifiers=expected_identifiers,
                        identifier_type=identifier_type,
                    ),
                    ["1002", "1001", "1002"],
                )
        self.assertEqual(
            _validate_xml_record_uids(
                xml_file_path=xml_file_path,
                expected_protein_uids=expected_identifiers,
            ),
            ["1002", "1001", "1002"],
        )

    def test_uid_validation_rejects_changed_duplicate_counts(self) -> None:
        xml_file_path = self._write_xml(
            "<GBSet>"
            "<GBSeq><GBSeqid>gi|1001</GBSeqid></GBSeq>"
            "<GBSeq><GBSeqid>gi|1002</GBSeqid></GBSeq>"
            "<GBSeq><GBSeqid>gi|1002</GBSeqid></GBSeq>"
            "</GBSet>"
        )
        with self.assertRaisesRegex(
            RuntimeError,
            r"record UIDs do not match.*Missing: \['1001'\]; Unexpected: \['1002'\]",
        ):
            _validate_xml_record_identifiers(
                xml_file_path=xml_file_path,
                expected_protein_identifiers=["1001", "1001", "1002"],
                identifier_type="uid",
            )

    def test_uid_validation_preserves_extraction_errors(self) -> None:
        cases = [
            ("<GBSeqid>ref|WP_012345678.1|</GBSeqid>", "Failed to extract a protein UID"),
            ("<GBSeqid>gi|1001</GBSeqid><GBSeqid>gi|1002</GBSeqid>",
             "multiple conflicting protein UID candidates"),
        ]
        for record_contents, message in cases:
            with self.subTest(record_contents=record_contents):
                xml_file_path = self._write_xml(
                    f"<GBSet><GBSeq>{record_contents}</GBSeq></GBSet>"
                )
                with self.assertRaisesRegex(RuntimeError, message):
                    _validate_xml_record_identifiers(
                        xml_file_path=xml_file_path,
                        expected_protein_identifiers=["1001"],
                        identifier_type="uid",
                    )

    def test_uid_validation_accepts_legacy_gi_forms_in_direct_or_nested_seqids(self) -> None:
        for value in ("gi|12345", "gi|12345|", " \n gi|12345| \t"):
            for nested in (False, True):
                with self.subTest(value=value, nested=nested):
                    record = f"<GBSeqid>{value}</GBSeqid>"
                    if nested:
                        record = f"<GBSeq_other-seqids>{record}</GBSeq_other-seqids>"
                    xml_file_path = self._write_xml(f"<GBSet><GBSeq>{record}</GBSeq></GBSet>")
                    self.assertEqual(
                        _validate_xml_record_uids(
                            xml_file_path=xml_file_path,
                            expected_protein_uids=["12345"],
                        ),
                        ["12345"],
                    )

    def test_uid_validation_ignores_non_uid_fields_and_repeated_candidates(self) -> None:
        xml_file_path = self._write_xml(
            "<GBSet><GBSeq>"
            "<GBSeq_accession-version>WP_012345678.1</GBSeq_accession-version>"
            "<GBSeq_primary-accession>99999</GBSeq_primary-accession>"
            "<GBSeqid>gi|12345</GBSeqid>"
            "<GBSeq_other-seqids>"
            "<GBSeqid>ref|WP_012345678.1|</GBSeqid>"
            "<GBSeqid>ref|99999|</GBSeqid>"
            "<GBSeqid>99999</GBSeqid>"
            "<GBSeqid>gi|99999|ref|WP_012345678.1|</GBSeqid>"
            "<GBSeqid>gi|12345|</GBSeqid>"
            "<GBSeqid>gi|12345</GBSeqid>"
            "</GBSeq_other-seqids>"
            "</GBSeq></GBSet>"
        )
        self.assertEqual(
            _validate_xml_record_identifiers(
                xml_file_path=xml_file_path,
                expected_protein_identifiers=["12345"],
                identifier_type="uid",
            ),
            ["12345"],
        )

    def test_uid_validation_requires_a_complete_gi_candidate_in_each_record(self) -> None:
        invalid_records = [
            "",
            "<GBSeqid />",
            "<GBSeqid> \t </GBSeqid>",
            "<GBSeq_locus>gi|12345</GBSeq_locus>",
            "<GBSeq_primary-accession>12345</GBSeq_primary-accession>",
            "<GBSeq_accession-version>WP_012345678.1</GBSeq_accession-version>",
            *[
                f"<GBSeqid>{value}</GBSeqid>"
                for value in (
                    "12345", "ref|12345|", "ref|WP_012345678.1|", "gi|",
                    "gi|ABC", "gi|-12345", "gi|+12345", "gi|12345.1",
                    "gi|12345abc", "gi|12 345", "gi|12345||", "GI|12345",
                    "prefixgi|12345", "gi|12345|ref|WP_012345678.1|",
                )
            ],
        ]
        for record in invalid_records:
            with self.subTest(record=record):
                # A valid first record must not hide a missing UID in the second.
                xml_file_path = self._write_xml(
                    "<GBSet><GBSeq><GBSeqid>gi|12345</GBSeqid></GBSeq>"
                    f"<GBSeq>{record}</GBSeq></GBSet>"
                )
                with self.assertRaisesRegex(RuntimeError, "Failed to extract a protein UID"):
                    _validate_xml_record_identifiers(
                        xml_file_path=xml_file_path,
                        expected_protein_identifiers=["12345", "12345"],
                        identifier_type="uid",
                    )

    def test_uid_validation_rejects_conflicts_between_direct_and_nested_candidates(self) -> None:
        xml_file_path = self._write_xml(
            "<GBSet><GBSeq>"
            "<GBSeqid>gi|12345</GBSeqid>"
            "<GBSeq_other-seqids>"
            "<GBSeqid>gi|12345|</GBSeqid><GBSeqid>gi|99999</GBSeqid>"
            "</GBSeq_other-seqids>"
            "</GBSeq></GBSet>"
        )
        with self.assertRaisesRegex(RuntimeError, "multiple conflicting protein UID candidates"):
            _validate_xml_record_identifiers(
                xml_file_path=xml_file_path,
                expected_protein_identifiers=["12345"],
                identifier_type="uid",
            )

    def test_uid_validation_rejects_missing_or_extra_record_occurrences(self) -> None:
        for record_count, expected_count, missing, unexpected in (
            (1, 2, ["12345"], []),
            (2, 1, [], ["12345"]),
            (0, 1, ["12345"], []),
            (1, 0, [], ["12345"]),
        ):
            with self.subTest(record_count=record_count, expected_count=expected_count):
                record = "<GBSeq><GBSeqid>gi|12345</GBSeqid></GBSeq>"
                xml_file_path = self._write_xml(f"<GBSet>{record * record_count}</GBSet>")
                with self.assertRaisesRegex(RuntimeError, "record UIDs do not match") as error:
                    _validate_xml_record_identifiers(
                        xml_file_path=xml_file_path,
                        expected_protein_identifiers=["12345"] * expected_count,
                        identifier_type="uid",
                    )
                self.assertIn(f"Missing: {missing}; Unexpected: {unexpected}.", str(error.exception))

    def test_generic_validator_preserves_xml_structure_errors(self) -> None:
        cases = [
            ("<GBSet>", "not well-formed"),
            ("<OtherSet />", "root tag mismatch during UID validation"),
            ("<GBSet><OtherRecord /></GBSet>", "unexpected child tags during UID validation"),
        ]
        for xml_text, message in cases:
            with self.subTest(xml_text=xml_text):
                xml_file_path = self._write_xml(xml_text)
                with self.assertRaisesRegex(RuntimeError, message):
                    _validate_xml_record_identifiers(
                        xml_file_path=xml_file_path,
                        expected_protein_identifiers=["1001"],
                        identifier_type="uid",
                    )

    def test_legacy_wrapper_preserves_custom_tags_and_string_path(self) -> None:
        xml_file_path = self._write_xml(
            "<Records><Protein><GBSeqid>gi|1001</GBSeqid></Protein></Records>"
        )
        self.assertEqual(
            _validate_xml_record_uids(
                xml_file_path=str(xml_file_path),
                expected_protein_uids=["1001"],
                expected_root_tag="Records",
                expected_record_tag="Protein",
            ),
            ["1001"],
        )

    def test_legacy_wrapper_delegates_to_generic_uid_validation(self) -> None:
        identifiers = ["1001"]
        with patch(
            "src.pago_pipeline.ncbi_snapshot._validate_xml_record_identifiers",
            return_value=identifiers,
        ) as validate_mock:
            result = _validate_xml_record_uids(
                xml_file_path="records.xml",
                expected_protein_uids=identifiers,
                expected_root_tag="Records",
                expected_record_tag="Protein",
            )
        self.assertIs(result, identifiers)
        validate_mock.assert_called_once_with(
            xml_file_path="records.xml",
            expected_protein_identifiers=identifiers,
            identifier_type="uid",
            expected_root_tag="Records",
            expected_record_tag="Protein",
        )

    def test_unknown_identifier_types_are_rejected_before_xml_parsing(self) -> None:
        for identifier_type in ("unknown", "protein_accession"):
            with self.subTest(identifier_type=identifier_type):
                with patch("src.pago_pipeline.ncbi_snapshot.ET.parse") as parse_mock:
                    with self.assertRaisesRegex(ValueError, "Invalid NCBI protein identifier_type"):
                        _validate_xml_record_identifiers(
                            xml_file_path="records.xml",
                            expected_protein_identifiers=["WP_012345678.1"],
                            identifier_type=identifier_type,
                        )
                parse_mock.assert_not_called()


class ValidateXmlAccessionVersionTests(unittest.TestCase):
    def _validate_records(
        self,
        records: list[str],
        expected: list[str],
        identifier_type: str = "accession.version",
    ) -> list[str]:
        with tempfile.TemporaryDirectory() as directory:
            xml_file_path = Path(directory) / "protein_records.xml"
            xml_file_path.write_text(
                "<GBSet>" + "".join(f"<GBSeq>{record}</GBSeq>" for record in records) + "</GBSet>",
                encoding="utf-8",
            )
            return _validate_xml_record_identifiers(
                xml_file_path=xml_file_path,
                expected_protein_identifiers=expected,
                identifier_type=identifier_type,
            )

    def test_primary_accession_version_matches_requested_identifier(self) -> None:
        self.assertEqual(
            self._validate_records(
                ["<GBSeq_accession-version>WP_012345678.1</GBSeq_accession-version>"],
                ["WP_012345678.1"],
            ),
            ["WP_012345678.1"],
        )

    def test_primary_accession_version_and_gi_do_not_conflict(self) -> None:
        self.assertEqual(
            self._validate_records(
                [
                    "<GBSeq_accession-version>WP_012345678.1</GBSeq_accession-version>"
                    "<GBSeqid>gi|12345</GBSeqid>"
                ],
                ["WP_012345678.1"],
            ),
            ["WP_012345678.1"],
        )

    def test_direct_seqid_fallback_matches_requested_identifier(self) -> None:
        self.assertEqual(
            self._validate_records(
                ["<GBSeqid>ref|WP_012345678.1|</GBSeqid>"],
                ["WP_012345678.1"],
            ),
            ["WP_012345678.1"],
        )

    def test_primary_accession_alone_cannot_confirm_an_exact_revision(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "Failed to extract an accession.version"):
            self._validate_records(
                ["<GBSeq_primary-accession>WP_012345678</GBSeq_primary-accession>"],
                ["WP_012345678.1"],
            )

    def test_primary_accession_version_rejects_conflicting_seqids(self) -> None:
        for primary_value, fields in (
            ("WP_012345678.1", "<GBSeqid>ref|WP_999999999.1|</GBSeqid>"),
            ("WP_012345678.1", "<GBSeqid>ref|WP_012345678.2|</GBSeqid>"),
            ("WP_012345678.2", "<GBSeqid>ref|WP_012345678.1|</GBSeqid>"),
            (
                "WP_012345678.1",
                "<GBSeq_other-seqids><GBSeqid>ref|WP_999999999.1|</GBSeqid></GBSeq_other-seqids>",
            ),
        ):
            primary = f"<GBSeq_accession-version> {primary_value} </GBSeq_accession-version>"
            with self.subTest(primary=primary_value, fields=fields), self.assertRaisesRegex(
                RuntimeError, "conflicting accession.version candidates.*GBSeq_accession-version.*GBSeqid",
            ):
                self._validate_records([primary + fields], ["WP_012345678.1"])

    def test_primary_accession_version_accepts_matching_seqids_and_non_candidates(self) -> None:
        primary = "<GBSeq_accession-version> WP_012345678.1 </GBSeq_accession-version>"
        other_fields = [
            "<GBSeqid>ref|WP_012345678.1|</GBSeqid>"
            "<GBSeq_other-seqids><GBSeqid>ref|WP_012345678.1|</GBSeqid></GBSeq_other-seqids>",
            "<GBSeqid>ref|WP_999999999|</GBSeqid>",
            "<GBSeq_primary-accession>WP_999999999</GBSeq_primary-accession>"
            "<GBSeq_sequence-version>99</GBSeq_sequence-version>",
            primary,
        ]
        for fields in other_fields:
            with self.subTest(fields=fields):
                self.assertEqual(
                    self._validate_records([primary + fields], ["WP_012345678.1"]),
                    ["WP_012345678.1"],
                )

    def test_versioned_seqids_identify_records_without_a_primary_accession_version(self) -> None:
        records = [
            "<GBSeq_other-seqids><GBSeqid>ref|WP_012345678.1|</GBSeqid></GBSeq_other-seqids>",
            "<GBSeqid>gi|12345</GBSeqid><GBSeqid>ref|WP_012345678.1|</GBSeqid>",
            "<GBSeqid>ref|WP_012345678.1|</GBSeqid>"
            "<GBSeq_other-seqids><GBSeqid>ref|WP_012345678.1|</GBSeqid>"
            "<GBSeqid>gi|99999</GBSeqid></GBSeq_other-seqids>",
            "<GBSeqid>ref|WP_012345678.1|</GBSeqid>"
            "<GBSeq_primary-accession>WP_999999999</GBSeq_primary-accession>"
            "<GBSeq_sequence-version>99</GBSeq_sequence-version>",
            "<GBSeqid>ref|WP_012345678.1|</GBSeqid><GBSeqid>ref|WP_999999999|</GBSeqid>",
        ]
        records.extend(
            f"<GBSeq_other-seqids><GBSeqid> {prefix}|WP_012345678.1| </GBSeqid></GBSeq_other-seqids>"
            for prefix in ("ref", "gb", "emb", "dbj", "tpg", "tpe", "tpd", "sp", "tr", "pir", "prf")
        )
        for record in records:
            with self.subTest(record=record):
                self.assertEqual(
                    self._validate_records([record], ["WP_012345678.1"]),
                    ["WP_012345678.1"],
                )

    def test_empty_or_invalid_primary_does_not_allow_fallback_rescue(self) -> None:
        for value in ("", " \t", "WP_012345678", "WP_012345678.0", "WP_012345678.A", "12345.1"):
            with self.subTest(primary_value=value):
                record = (
                    f"<GBSeq_accession-version>{value}</GBSeq_accession-version>"
                    "<GBSeqid>ref|WP_012345678.1|</GBSeqid>"
                )
                with self.assertRaisesRegex(RuntimeError, "Invalid accession.version in GBSeq_accession-version"):
                    self._validate_records([record], ["WP_012345678.1"])

    def test_conflicting_primary_candidates_are_rejected(self) -> None:
        for conflicting_value in ("WP_999999999.1", "WP_012345678.2"):
            with self.subTest(conflicting_value=conflicting_value):
                record = (
                    "<GBSeq_accession-version>WP_012345678.1</GBSeq_accession-version>"
                    f"<GBSeq_accession-version>{conflicting_value}</GBSeq_accession-version>"
                    "<GBSeqid>ref|WP_012345678.1|</GBSeqid>"
                )
                with self.assertRaisesRegex(RuntimeError, "conflicting accession.version candidates.*GBSeq_accession-version"):
                    self._validate_records([record], ["WP_012345678.1"])

    def test_conflicting_seqid_candidates_are_rejected(self) -> None:
        for second_seqid in (
            "<GBSeqid>ref|WP_999999999.1|</GBSeqid>",
            "<GBSeq_other-seqids><GBSeqid>ref|WP_999999999.1|</GBSeqid></GBSeq_other-seqids>",
            "<GBSeq_other-seqids><GBSeqid>gb|WP_012345678.2|</GBSeqid></GBSeq_other-seqids>",
        ):
            with self.subTest(second_seqid=second_seqid):
                with self.assertRaisesRegex(RuntimeError, "conflicting accession.version candidates.*GBSeqid"):
                    self._validate_records(
                        ["<GBSeqid>ref|WP_012345678.1|</GBSeqid>" + second_seqid],
                        ["WP_012345678.1"],
                    )

    def test_missing_or_unversioned_candidates_cannot_confirm_an_exact_revision(self) -> None:
        invalid_records = [
            "",
            "<GBSeqid>gi|12345</GBSeqid>",
            "<GBSeq_primary-accession>WP_012345678.1</GBSeq_primary-accession>",
            "<GBSeq_sequence-version>WP_012345678.1</GBSeq_sequence-version>",
            "<GBSeq_primary-accession>WP_012345678</GBSeq_primary-accession>"
            "<GBSeq_sequence-version>1</GBSeq_sequence-version>",
            *[
                f"<GBSeq_other-seqids><GBSeqid>{seqid}</GBSeqid></GBSeq_other-seqids>"
                for seqid in (
                    "", " \t", "ref|WP_012345678|", "ref|WP_012345678.0|",
                    "ref|WP_012345678.A|", "ref|12345.1|", "ref|12345|",
                    "WP_012345678.1", "gi|WP_012345678.1|", "lcl|WP_012345678.1|",
                    "ref|WP_012345678.1||", "ref|WP_012345678.1|junk",
                )
            ],
        ]
        for record in invalid_records:
            with self.subTest(record=record):
                with self.assertRaisesRegex(RuntimeError, "Failed to extract an accession.version"):
                    self._validate_records(
                        [record, "<GBSeqid>ref|WP_012345678.1|</GBSeqid>"],
                        ["WP_012345678.1", "WP_012345678.1"],
                    )

    def test_record_accession_or_revision_mismatch_is_rejected(self) -> None:
        for record, unexpected_identifier in (
            (
                "<GBSeq_accession-version>WP_999999999.1</GBSeq_accession-version>",
                "WP_999999999.1",
            ),
            (
                "<GBSeq_accession-version>WP_012345678.2</GBSeq_accession-version>",
                "WP_012345678.2",
            ),
            ("<GBSeqid>ref|WP_012345678.2|</GBSeqid>", "WP_012345678.2"),
        ):
            with self.subTest(record=record):
                with self.assertRaisesRegex(RuntimeError, "accession.version identifiers do not match") as error:
                    self._validate_records([record], ["WP_012345678.1"])
                self.assertIn(
                    "Missing: ['WP_012345678.1']; "
                    f"Unexpected: ['{unexpected_identifier}'].",
                    str(error.exception),
                )

    def test_counter_preserves_duplicates_and_returns_xml_order(self) -> None:
        records = [
            "<GBSeq_accession-version>WP_999999999.1</GBSeq_accession-version>",
            "<GBSeqid>ref|WP_012345678.1|</GBSeqid>",
            "<GBSeq_other-seqids><GBSeqid>ref|WP_999999999.1|</GBSeqid></GBSeq_other-seqids>",
        ]
        self.assertEqual(
            self._validate_records(
                records,
                ["WP_012345678.1", "WP_999999999.1", "WP_999999999.1"],
                identifier_type=NCBIProteinIdentifierType.accession_version,
            ),
            ["WP_999999999.1", "WP_012345678.1", "WP_999999999.1"],
        )
        with self.assertRaisesRegex(RuntimeError, "accession.version identifiers do not match") as error:
            self._validate_records(records, ["WP_012345678.1", "WP_012345678.1", "WP_999999999.1"])
        self.assertIn("Missing: ['WP_012345678.1']; Unexpected: ['WP_999999999.1']", str(error.exception))


class SaveNcbiProteinIdentifierSnapshotTests(unittest.TestCase):
    def test_saves_programmatic_accession_version_identifier_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_root_directory = Path(temporary_directory_name) / "identifiers"

            snapshot_directory = save_ncbi_protein_identifier_snapshot(
                protein_identifiers=[
                    " WP_012345678.1 ",
                    "",
                    "XP_123456789.2",
                    "WP_012345678.1",
                ],
                identifier_type="accession.version",
                snapshot_root_directory=snapshot_root_directory,
                source_description="manual accession.version list",
                retrieved_at_utc="2026-03-21T18:42:03Z",
            )

            snapshot_payload = load_snapshot_by_directory(
                snapshot_directory=snapshot_directory,
            )
            latest_payload = load_snapshot_by_directory(
                snapshot_directory=snapshot_root_directory / "latest",
            )

        self.assertEqual(
            snapshot_payload["protein_uids"],
            ["WP_012345678.1", "XP_123456789.2"],
        )
        self.assertEqual(latest_payload["protein_uids"], snapshot_payload["protein_uids"])
        self.assertEqual(
            snapshot_payload["manifest"]["artifact_type"],
            "ncbi_protein_identifier_snapshot",
        )
        self.assertEqual(
            snapshot_payload["manifest"]["source_type"],
            "provided_identifier_list",
        )
        self.assertEqual(
            snapshot_payload["manifest"]["source_description"],
            "manual accession.version list",
        )
        self.assertNotIn("search_query", snapshot_payload["manifest"])
        self.assertEqual(
            snapshot_payload["manifest"]["identifier_type"],
            "accession.version",
        )
        self.assertEqual(
            snapshot_payload["manifest"]["protein_identifiers_file_name"],
            "protein_identifiers.txt",
        )
        self.assertEqual(
            snapshot_payload["manifest"]["protein_identifiers_sha256"],
            sha256_of_lines(
                text_lines=["WP_012345678.1", "XP_123456789.2"],
                deduplicate_lines_preserving_order=False,
                sort_lines=False,
            ),
        )
        self.assertEqual(
            snapshot_payload["manifest"]["protein_uids_sha256"],
            snapshot_payload["manifest"]["protein_identifiers_sha256"],
        )
        manifest = snapshot_payload["manifest"]
        self.assertEqual(manifest["snapshot_format_version"], "1.1")
        self.assertEqual(manifest["raw_protein_identifier_count"], 4)
        self.assertEqual(manifest["normalized_protein_identifier_count"], 2)
        self.assertEqual(manifest["protein_identifier_count"], 2)
        self.assertEqual(
            manifest["normalized_protein_uid_count"],
            manifest["normalized_protein_identifier_count"],
        )

    def test_saves_identifier_snapshot_from_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            identifier_file_path = temporary_directory / "provided_uids.txt"
            identifier_file_path.write_text(
                "1002\n\n1001\n1002\n",
                encoding="utf-8",
            )
            snapshot_root_directory = temporary_directory / "identifiers"

            snapshot_directory = save_ncbi_protein_identifier_snapshot_from_file(
                identifier_file_path=identifier_file_path,
                identifier_type="uid",
                snapshot_root_directory=snapshot_root_directory,
                source_description="manual uid list",
                retrieved_at_utc="2026-03-21T18:42:03Z",
            )

            snapshot_payload = load_snapshot_by_directory(
                snapshot_directory=snapshot_directory,
            )

        self.assertEqual(snapshot_payload["protein_uids"], ["1001", "1002"])
        self.assertEqual(
            snapshot_payload["manifest"]["source_file_name"],
            "provided_uids.txt",
        )
        self.assertEqual(
            snapshot_payload["manifest"]["source_type"],
            "provided_identifier_list",
        )
        self.assertEqual(
            snapshot_payload["manifest"]["source_description"],
            "manual uid list",
        )
        self.assertNotIn("search_query", snapshot_payload["manifest"])
        self.assertEqual(snapshot_payload["manifest"]["identifier_type"], "uid")

    def test_rejects_invalid_programmatic_identifier_list(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            with self.assertRaisesRegex(
                ValueError,
                "Mixed or invalid identifier lists are not supported",
            ):
                save_ncbi_protein_identifier_snapshot(
                    protein_identifiers=["WP_012345678"],
                    identifier_type="accession.version",
                    snapshot_root_directory=Path(temporary_directory_name),
                )

    def test_legacy_esearch_uid_snapshot_uses_generic_identifier_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_root_directory = Path(temporary_directory_name) / "uid_snapshot"
            protein_uids = ["1001", "1002"]
            protein_uids_sha256 = sha256_of_lines(
                text_lines=protein_uids,
                deduplicate_lines_preserving_order=False,
                sort_lines=False,
            )
            fetch_result = NCBIProteinUidFetchResult(
                database_name="protein",
                search_query="txid9606[Organism]",
                translated_query="txid9606[Organism]",
                identifier_type="uid",
                retrieved_at_utc="2026-03-21T18:42:03Z",
                ncbi_reported_result_count=2,
                protein_uids=protein_uids,
                raw_protein_uid_count=2,
                normalized_protein_uid_count=2,
                deduplicate_uids=True,
                sort_uids=True,
                protein_uids_sha256=protein_uids_sha256,
                page_size=1000,
                request_delay_seconds=0.34,
                max_retry_attempts=5,
                history_web_env="web-env",
                history_query_key="1",
                python_version="3.12.0",
                biopython_version="test",
            )

            snapshot_directory = save_ncbi_protein_uid_snapshot(
                fetch_result=fetch_result,
                snapshot_root_directory=snapshot_root_directory,
            )
            manifest = json.loads(
                (snapshot_directory / "manifest.json").read_text(
                    encoding="utf-8",
                )
            )

        self.assertEqual(
            manifest["artifact_type"],
            "ncbi_protein_identifier_snapshot",
        )
        self.assertEqual(manifest["source_type"], "entrez_search")
        self.assertEqual(manifest["search_query"], "txid9606[Organism]")
        self.assertEqual(manifest["identifier_type"], "uid")
        self.assertEqual(manifest["snapshot_format_version"], "1.1")
        self.assertEqual(manifest["protein_identifiers_file_name"], "protein_identifiers.txt")
        self.assertEqual(manifest["raw_protein_identifier_count"], 2)
        self.assertEqual(manifest["normalized_protein_identifier_count"], 2)
        self.assertEqual(manifest["protein_identifier_count"], 2)
        self.assertEqual(manifest["raw_protein_uid_count"], 2)
        self.assertEqual(manifest["normalized_protein_uid_count"], 2)
        self.assertEqual(manifest["deduplicate_identifiers"], manifest["deduplicate_uids"])
        self.assertEqual(manifest["sort_identifiers"], manifest["sort_uids"])
        self.assertEqual(
            manifest["protein_identifiers_sha256"],
            protein_uids_sha256,
        )

        # Old constructor/serialization contracts remain usable in either domain.
        for identifier_type, identifiers in (
            ("uid", ["1001", "1002"]),
            ("accession.version", ["WP_012345678.1", "XP_123456789.2"]),
        ):
            with self.subTest(identifier_type=identifier_type):
                result = replace(
                    fetch_result,
                    identifier_type=identifier_type,
                    protein_uids=identifiers,
                    raw_protein_uid_count=3,
                    protein_uids_sha256=sha256_of_lines(text_lines=identifiers),
                )
                self.assertIs(result.protein_identifiers, result.protein_uids)
                self.assertEqual(result.raw_protein_identifier_count, 3)
                self.assertEqual(result.normalized_protein_identifier_count, 2)
                self.assertEqual(result.protein_identifier_count, 2)
                self.assertEqual(result.protein_identifiers_sha256, result.protein_uids_sha256)
                self.assertEqual(NCBIProteinUidFetchResult(**asdict(result)), result)
                self.assertNotIn("protein_identifiers", asdict(result))


class NcbiIdentifierManifestFormatTests(unittest.TestCase):
    def test_new_identifier_snapshots_write_1_1_canonical_file_and_legacy_copy(self) -> None:
        for identifier_type, identifiers in (
            ("uid", ["1002", "1001", "1002"]),
            ("accession.version", ["XP_123456789.2", "WP_012345678.1", "XP_123456789.2"]),
        ):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                root = Path(name)
                snapshot = save_ncbi_protein_identifier_snapshot(
                    protein_identifiers=identifiers,
                    identifier_type=identifier_type,
                    snapshot_root_directory=root,
                    deduplicate_identifiers=False,
                    sort_identifiers=False,
                )
                for directory in (snapshot, root / "latest"):
                    payload = load_snapshot_by_directory(snapshot_directory=directory)
                    manifest = payload["manifest"]
                    self.assertEqual(manifest["snapshot_format_version"], "1.1")
                    self.assertEqual(manifest["artifact_type"], "ncbi_protein_identifier_snapshot")
                    self.assertEqual(manifest["source_type"], "provided_identifier_list")
                    self.assertEqual(manifest["identifier_type"], identifier_type)
                    self.assertEqual(manifest["protein_identifiers_file_name"], "protein_identifiers.txt")
                    self.assertEqual(manifest["normalized_protein_identifier_count"], 3)
                    self.assertEqual(payload["protein_uids"], identifiers)
                    self.assertEqual(payload["protein_uids_file_path"], directory / "protein_identifiers.txt")
                    self.assertEqual(
                        (directory / "protein_identifiers.txt").read_bytes(),
                        (directory / "protein_uids.txt").read_bytes(),
                    )
                    self.assertEqual(
                        manifest["protein_identifiers_sha256"],
                        sha256_of_lines(text_lines=identifiers, deduplicate_lines_preserving_order=False),
                    )

    def test_identifier_reuse_and_path_helpers_work_without_legacy_file(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            save_ncbi_protein_identifier_snapshot(
                protein_identifiers=["WP_012345678.1"],
                identifier_type="accession.version",
                snapshot_root_directory=root,
            )
            latest = root / "latest"
            (latest / "protein_uids.txt").unlink()
            self.assertTrue(latest_snapshot_is_available(snapshot_root_directory=root))
            self.assertEqual(get_latest_snapshot_protein_uids_path(snapshot_root_directory=root), latest / "protein_identifiers.txt")
            self.assertEqual(get_snapshot_protein_uids_path(snapshot_directory=latest), latest / "protein_identifiers.txt")
            with patch("src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_uid_snapshot") as fetch_mock:
                for mode in ("reuse_latest", "reuse_latest_or_create"):
                    payload = resolve_ncbi_protein_uid_snapshot(snapshot_mode=mode, snapshot_root_directory=root)
                    self.assertEqual(payload["protein_uids"], ["WP_012345678.1"])
                fetch_mock.assert_not_called()

    def test_loader_uses_legacy_filename_fields_only_when_generic_field_is_absent(self) -> None:
        for filename_field in ("protein_uids_file_name", "snapshot_file_name", None):
            with self.subTest(filename_field=filename_field), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                filename = "provided_ids.txt" if filename_field else "protein_uids.txt"
                (directory / filename).write_text("1001\n", encoding="utf-8")
                manifest = {"snapshot_format_version": "1.0"}
                if filename_field:
                    manifest[filename_field] = filename
                (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
                self.assertEqual(load_snapshot_by_directory(snapshot_directory=directory)["protein_uids"], ["1001"])

    def test_invalid_generic_filename_is_not_replaced_by_valid_legacy_filename(self) -> None:
        for filename in (None, "", " \t", 123):
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as name:
                root = Path(name)
                snapshot = save_ncbi_protein_identifier_snapshot(
                    protein_identifiers=["1001"], identifier_type="uid", snapshot_root_directory=root,
                )
                manifest_file = snapshot / "manifest.json"
                manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
                manifest["protein_identifiers_file_name"] = filename
                manifest_file.write_text(json.dumps(manifest), encoding="utf-8")
                with self.assertRaisesRegex(RuntimeError, "non-empty protein identifier file name"):
                    load_snapshot_by_directory(snapshot_directory=snapshot)

    def test_generic_file_has_priority_over_different_legacy_contents(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            snapshot = save_ncbi_protein_identifier_snapshot(
                protein_identifiers=["WP_012345678.1"], identifier_type="accession.version", snapshot_root_directory=root,
            )
            (snapshot / "protein_uids.txt").write_text("1001\n", encoding="utf-8")
            payload = load_snapshot_by_directory(snapshot_directory=snapshot)
            self.assertEqual(payload["protein_uids"], ["WP_012345678.1"])

    def test_only_unversioned_and_1_0_manifests_default_missing_type_to_uid(self) -> None:
        for version in (None, "1.0", "1.1", "2.0"):
            with self.subTest(version=version), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                (directory / "protein_uids.txt").write_text("1001\n", encoding="utf-8")
                manifest = {"protein_uids_file_name": "protein_uids.txt"}
                if version is not None:
                    manifest["snapshot_format_version"] = version
                (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
                if version in (None, "1.0"):
                    self.assertEqual(load_snapshot_by_directory(snapshot_directory=directory)["protein_uids"], ["1001"])
                else:
                    error = "identifier_type is required" if version == "1.1" else "Unsupported snapshot_format_version"
                    with self.assertRaisesRegex(ValueError, error):
                        load_snapshot_by_directory(snapshot_directory=directory)


class LegacyNcbiManifestCompatibilityTests(unittest.TestCase):
    def _uid_fetch_result(self, identifier_type="uid") -> NCBIProteinUidFetchResult:
        return NCBIProteinUidFetchResult(
            database_name="protein", search_query="test protein query", translated_query=None,
            identifier_type=identifier_type, retrieved_at_utc="2026-10-26T12:00:00Z",
            ncbi_reported_result_count=2, protein_uids=["1001", "1002"],
            raw_protein_uid_count=2, normalized_protein_uid_count=2,
            deduplicate_uids=False, sort_uids=False,
            protein_uids_sha256=sha256_of_lines(text_lines=["1001", "1002"]),
            page_size=1000, request_delay_seconds=0.0, max_retry_attempts=5,
            history_web_env=None, history_query_key=None,
            python_version="3.14.0", biopython_version="test",
        )

    def _write_legacy_source(self, directory: Path) -> tuple[Path, dict]:
        (directory / "protein_uids.txt").write_text("1001\n1002\n", encoding="utf-8")
        manifest = {
            "protein_uids_sha256": sha256_of_lines(text_lines=["1001", "1002"]),
            "normalized_protein_uid_count": 2,
        }
        path = directory / "manifest.json"
        path.write_text(json.dumps(manifest), encoding="utf-8")
        return path, manifest

    def test_esearch_save_does_not_publish_an_inconsistent_identifier_list(self) -> None:
        for identifiers, error in (
            (["1001", "1002", "1003"], "identifier count mismatch"),
            (["1001", "1003"], "identifier hash mismatch"),
        ):
            with self.subTest(identifiers=identifiers), tempfile.TemporaryDirectory() as name:
                root = Path(name)
                result = self._uid_fetch_result()
                original = save_ncbi_protein_uid_snapshot(
                    fetch_result=result, snapshot_root_directory=root,
                )
                latest_manifest = root / "latest" / "manifest.json"
                original_manifest = latest_manifest.read_bytes()
                inconsistent = replace(
                    result, protein_uids=identifiers, retrieved_at_utc="2026-10-26T13:00:00Z",
                )
                with self.assertRaisesRegex(RuntimeError, error):
                    save_ncbi_protein_uid_snapshot(
                        fetch_result=inconsistent, snapshot_root_directory=root,
                    )
                self.assertEqual(latest_manifest.read_bytes(), original_manifest)
                self.assertEqual(list((root / "snapshots").iterdir()), [original])

    def test_reads_minimal_legacy_manifest_without_rewriting_it(self) -> None:
        for version in (None, "1.0"):
            with self.subTest(version=version), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                path, manifest = self._write_legacy_source(directory)
                if version is not None:
                    manifest["snapshot_format_version"] = version
                    path.write_text(json.dumps(manifest), encoding="utf-8")
                original = path.read_bytes()
                loaded = load_snapshot_by_directory(snapshot_directory=directory)
                self.assertEqual(loaded["protein_uids"], ["1001", "1002"])
                self.assertEqual(loaded["manifest"], manifest)
                self.assertNotIn("identifier_type", loaded["manifest"])
                self.assertEqual(path.read_bytes(), original)

    def test_explicit_invalid_type_in_legacy_manifest_does_not_default_to_uid(self) -> None:
        for identifier_type in (None, "", "unknown"):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                path, manifest = self._write_legacy_source(directory)
                manifest["identifier_type"] = identifier_type
                path.write_text(json.dumps(manifest), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "identifier_type"):
                    load_snapshot_by_directory(snapshot_directory=directory)

    def test_legacy_default_does_not_infer_accessions_from_file_contents(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            self._write_legacy_source(directory)
            (directory / "protein_uids.txt").write_text("WP_012345678.1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "must contain only uid identifiers"):
                load_snapshot_by_directory(snapshot_directory=directory)

    def test_new_identifier_snapshots_require_explicit_type_before_writing(self) -> None:
        for workflow in ("provided_list", "esearch_result"):
            for identifier_type in (None, "", "unknown"):
                with self.subTest(workflow=workflow, identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                    root = Path(name) / "new_snapshot"
                    with self.assertRaisesRegex(ValueError, "identifier_type"):
                        if workflow == "provided_list":
                            save_ncbi_protein_identifier_snapshot(
                                protein_identifiers=["1001", "1002"],
                                identifier_type=identifier_type,
                                snapshot_root_directory=root,
                            )
                        else:
                            save_ncbi_protein_uid_snapshot(
                                fetch_result=self._uid_fetch_result(identifier_type),
                                snapshot_root_directory=root,
                            )
                    self.assertFalse(root.exists())

    def test_saving_from_legacy_values_writes_new_explicit_type(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            path, _ = self._write_legacy_source(directory)
            original = path.read_bytes()
            loaded = load_snapshot_by_directory(snapshot_directory=directory)
            snapshot = save_ncbi_protein_uid_snapshot(
                fetch_result=self._uid_fetch_result(),
                snapshot_root_directory=directory / "new_snapshot",
            )
            manifest = json.loads((snapshot / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["snapshot_format_version"], "1.1")
            self.assertEqual(manifest["identifier_type"], "uid")
            self.assertEqual(manifest["normalized_protein_identifier_count"], len(loaded["protein_uids"]))
            self.assertEqual(path.read_bytes(), original)


class SaveNcbiProteinXmlSnapshotTests(unittest.TestCase):
    def _write_source_uid_manifest(
        self,
        directory: Path,
        *,
        search_query: str,
        protein_uids_sha256: str,
        normalized_protein_uid_count: int,
    ) -> tuple[dict[str, object], Path]:
        manifest_payload: dict[str, object] = {
            "search_query": search_query,
            "protein_uids_sha256": protein_uids_sha256,
            "normalized_protein_uid_count": normalized_protein_uid_count,
            "immutable_snapshot_relative_path": "snapshots/source_uid_snapshot",
            "immutable_snapshot_directory_name": "source_uid_snapshot",
        }
        manifest_file_path = directory / "source_manifest.json"
        manifest_file_path.write_text(
            json.dumps(manifest_payload, indent=2) + "\n",
            encoding="utf-8",
        )
        return manifest_payload, manifest_file_path

    def _build_gbseq_xml_record(self, uid: str, locus: str) -> bytes:
        return (
            b"<GBSeq>\n"
            + f"  <GBSeqid>gi|{uid}</GBSeqid>\n".encode("utf-8")
            + f"  <GBSeq_locus>{locus}</GBSeq_locus>\n".encode("utf-8")
            + b"</GBSeq>\n"
        )

    def test_raises_when_consolidated_xml_record_count_is_less_than_expected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            snapshot_root_directory = temporary_directory / "xml_snapshot"
            protein_uids = ["1001", "1002"]
            protein_uids_sha256 = sha256_of_lines(
                text_lines=protein_uids,
                deduplicate_lines_preserving_order=False,
                sort_lines=False,
            )
            source_uid_snapshot_manifest, source_uid_snapshot_manifest_file_path = (
                self._write_source_uid_manifest(
                    temporary_directory,
                    search_query="txid9606[Organism]",
                    protein_uids_sha256=protein_uids_sha256,
                    normalized_protein_uid_count=len(protein_uids),
                )
            )
            fetch_result = NCBIProteinXmlFetchResult(
                database_name="protein",
                identifier_type="uid",
                retrieved_at_utc="2026-03-21T18:42:03Z",
                requested_protein_uid_count=len(protein_uids),
                normalized_protein_uid_count=len(protein_uids),
                protein_uids_sha256=protein_uids_sha256,
                batch_size=100,
                batch_count=1,
                retmode="xml",
                request_delay_seconds=0.34,
                max_retry_attempts=5,
                python_version="3.12.0",
                biopython_version="test",
                xml_batches=[
                    NCBIProteinXmlBatchFetchResult(
                        batch_index=1,
                        batch_start_index=0,
                        batch_end_index=1,
                        protein_uids=protein_uids,
                        protein_uid_count=len(protein_uids),
                        xml_payload_bytes=b"".join(
                            [
                                b'<?xml version="1.0" encoding="utf-8"?>\n',
                                b"<GBSet>\n",
                                self._build_gbseq_xml_record("1001", "ONLY_ONE"),
                                b"</GBSet>\n",
                            ]
                        ),
                        xml_payload_sha256="deadbeef",
                    )
                ],
            )

            with self.assertRaisesRegex(
                RuntimeError,
                "does not match the expected protein UID count",
            ):
                save_ncbi_protein_xml_snapshot(
                    fetch_result=fetch_result,
                    snapshot_root_directory=snapshot_root_directory,
                    source_uid_snapshot_manifest=source_uid_snapshot_manifest,
                    source_uid_snapshot_manifest_file_path=source_uid_snapshot_manifest_file_path,
                    protein_uids=protein_uids,
                )

    def test_raises_when_xml_record_uids_do_not_match_expected_uids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            snapshot_root_directory = temporary_directory / "xml_snapshot"
            protein_uids = ["1001"]
            protein_uids_sha256 = sha256_of_lines(
                text_lines=protein_uids,
                deduplicate_lines_preserving_order=False,
                sort_lines=False,
            )
            source_uid_snapshot_manifest, source_uid_snapshot_manifest_file_path = (
                self._write_source_uid_manifest(
                    temporary_directory,
                    search_query="txid9606[Organism]",
                    protein_uids_sha256=protein_uids_sha256,
                    normalized_protein_uid_count=len(protein_uids),
                )
            )
            fetch_result = NCBIProteinXmlFetchResult(
                database_name="protein",
                identifier_type="uid",
                retrieved_at_utc="2026-03-21T18:42:03Z",
                requested_protein_uid_count=len(protein_uids),
                normalized_protein_uid_count=len(protein_uids),
                protein_uids_sha256=protein_uids_sha256,
                batch_size=100,
                batch_count=1,
                retmode="xml",
                request_delay_seconds=0.34,
                max_retry_attempts=5,
                python_version="3.12.0",
                biopython_version="test",
                xml_batches=[
                    NCBIProteinXmlBatchFetchResult(
                        batch_index=1,
                        batch_start_index=0,
                        batch_end_index=0,
                        protein_uids=protein_uids,
                        protein_uid_count=len(protein_uids),
                        xml_payload_bytes=b"".join(
                            [
                                b'<?xml version="1.0" encoding="utf-8"?>\n',
                                b"<GBSet>\n",
                                self._build_gbseq_xml_record("9999", "WRONG_UID"),
                                b"</GBSet>\n",
                            ]
                        ),
                        xml_payload_sha256="deadbeef",
                    )
                ],
            )

            with self.assertRaisesRegex(RuntimeError, "record UIDs do not match"):
                save_ncbi_protein_xml_snapshot(
                    fetch_result=fetch_result,
                    snapshot_root_directory=snapshot_root_directory,
                    source_uid_snapshot_manifest=source_uid_snapshot_manifest,
                    source_uid_snapshot_manifest_file_path=source_uid_snapshot_manifest_file_path,
                    protein_uids=protein_uids,
                )

    def test_cleans_up_incomplete_snapshot_directory_when_xml_consolidation_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            snapshot_root_directory = temporary_directory / "xml_snapshot"
            protein_uids = ["1001"]
            protein_uids_sha256 = sha256_of_lines(
                text_lines=protein_uids,
                deduplicate_lines_preserving_order=False,
                sort_lines=False,
            )
            retrieved_at_utc = "2026-03-21T18:42:03Z"
            search_query = "txid9606[Organism]"
            source_uid_snapshot_manifest, source_uid_snapshot_manifest_file_path = (
                self._write_source_uid_manifest(
                    temporary_directory,
                    search_query=search_query,
                    protein_uids_sha256=protein_uids_sha256,
                    normalized_protein_uid_count=len(protein_uids),
                )
            )
            fetch_result = NCBIProteinXmlFetchResult(
                database_name="protein",
                identifier_type="uid",
                retrieved_at_utc=retrieved_at_utc,
                requested_protein_uid_count=len(protein_uids),
                normalized_protein_uid_count=len(protein_uids),
                protein_uids_sha256=protein_uids_sha256,
                batch_size=100,
                batch_count=1,
                retmode="xml",
                request_delay_seconds=0.34,
                max_retry_attempts=5,
                python_version="3.12.0",
                biopython_version="test",
                xml_batches=[
                    NCBIProteinXmlBatchFetchResult(
                        batch_index=1,
                        batch_start_index=0,
                        batch_end_index=0,
                        protein_uids=protein_uids,
                        protein_uid_count=len(protein_uids),
                        xml_payload_bytes=b"<GBSet><GBSeq>",
                        xml_payload_sha256="deadbeef",
                    )
                ],
            )
            expected_snapshot_directory = (
                snapshot_root_directory
                / "snapshots"
                / build_snapshot_directory_name(
                    retrieved_at_utc=retrieved_at_utc,
                    search_query=search_query,
                )
            )

            with self.assertRaisesRegex(RuntimeError, "Failed to parse XML batch"):
                save_ncbi_protein_xml_snapshot(
                    fetch_result=fetch_result,
                    snapshot_root_directory=snapshot_root_directory,
                    source_uid_snapshot_manifest=source_uid_snapshot_manifest,
                    source_uid_snapshot_manifest_file_path=source_uid_snapshot_manifest_file_path,
                    protein_uids=protein_uids,
                )

            self.assertFalse(expected_snapshot_directory.exists())

    def test_resolve_reuse_latest_rejects_tampered_xml_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            snapshot_root_directory = temporary_directory / "xml_snapshot"
            protein_uids = ["1001"]
            protein_uids_sha256 = sha256_of_lines(
                text_lines=protein_uids,
                deduplicate_lines_preserving_order=False,
                sort_lines=False,
            )
            source_uid_snapshot_manifest, source_uid_snapshot_manifest_file_path = (
                self._write_source_uid_manifest(
                    temporary_directory,
                    search_query="txid9606[Organism]",
                    protein_uids_sha256=protein_uids_sha256,
                    normalized_protein_uid_count=len(protein_uids),
                )
            )
            fetch_result = NCBIProteinXmlFetchResult(
                database_name="protein",
                identifier_type="uid",
                retrieved_at_utc="2026-03-21T18:42:03Z",
                requested_protein_uid_count=len(protein_uids),
                normalized_protein_uid_count=len(protein_uids),
                protein_uids_sha256=protein_uids_sha256,
                batch_size=100,
                batch_count=1,
                retmode="xml",
                request_delay_seconds=0.34,
                max_retry_attempts=5,
                python_version="3.12.0",
                biopython_version="test",
                xml_batches=[
                    NCBIProteinXmlBatchFetchResult(
                        batch_index=1,
                        batch_start_index=0,
                        batch_end_index=0,
                        protein_uids=protein_uids,
                        protein_uid_count=len(protein_uids),
                        xml_payload_bytes=b"".join(
                            [
                                b'<?xml version="1.0" encoding="utf-8"?>\n',
                                b"<GBSet>\n",
                                self._build_gbseq_xml_record("1001", "ORIGINAL"),
                                b"</GBSet>\n",
                            ]
                        ),
                        xml_payload_sha256="deadbeef",
                    )
                ],
            )

            save_ncbi_protein_xml_snapshot(
                fetch_result=fetch_result,
                snapshot_root_directory=snapshot_root_directory,
                source_uid_snapshot_manifest=source_uid_snapshot_manifest,
                source_uid_snapshot_manifest_file_path=source_uid_snapshot_manifest_file_path,
                protein_uids=protein_uids,
                update_latest_directory=True,
            )

            latest_xml_file_path = snapshot_root_directory / "latest" / "protein_records.xml"
            latest_xml_file_path.write_text(
                """<?xml version="1.0" encoding="utf-8"?>
<GBSet>
  <GBSeq>
    <GBSeqid>gi|9999</GBSeqid>
    <GBSeq_locus>TAMPERED</GBSeq_locus>
  </GBSeq>
</GBSet>
""",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(RuntimeError, "record UIDs do not match"):
                resolve_ncbi_protein_xml_snapshot(
                    snapshot_mode="reuse_latest",
                    snapshot_root_directory=snapshot_root_directory,
                    source_uid_snapshot_root_directory=temporary_directory / "unused_source",
                )

    def test_xml_snapshot_can_use_generic_identifier_source_without_search_query(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            identifier_snapshot_root_directory = temporary_directory / "identifiers"
            identifier_snapshot_directory = save_ncbi_protein_identifier_snapshot(
                protein_identifiers=["1001"],
                identifier_type="uid",
                snapshot_root_directory=identifier_snapshot_root_directory,
                retrieved_at_utc="2026-03-21T18:42:03Z",
                update_latest_directory=False,
            )
            source_identifier_snapshot_payload = load_snapshot_by_directory(
                snapshot_directory=identifier_snapshot_directory,
            )
            source_identifier_snapshot_manifest = source_identifier_snapshot_payload[
                "manifest"
            ]
            source_identifier_snapshot_manifest_file_path = (
                source_identifier_snapshot_payload["manifest_file_path"]
            )
            protein_uids = source_identifier_snapshot_payload["protein_uids"]
            protein_uids_sha256 = source_identifier_snapshot_manifest[
                "protein_uids_sha256"
            ]
            fetch_result = NCBIProteinXmlFetchResult(
                database_name="protein",
                identifier_type="uid",
                retrieved_at_utc="2026-03-21T18:43:03Z",
                requested_protein_uid_count=len(protein_uids),
                normalized_protein_uid_count=len(protein_uids),
                protein_uids_sha256=protein_uids_sha256,
                batch_size=100,
                batch_count=1,
                retmode="xml",
                request_delay_seconds=0.34,
                max_retry_attempts=5,
                python_version="3.12.0",
                biopython_version="test",
                xml_batches=[
                    NCBIProteinXmlBatchFetchResult(
                        batch_index=1,
                        batch_start_index=0,
                        batch_end_index=0,
                        protein_uids=protein_uids,
                        protein_uid_count=len(protein_uids),
                        xml_payload_bytes=b"".join(
                            [
                                b'<?xml version="1.0" encoding="utf-8"?>\n',
                                b"<GBSet>\n",
                                self._build_gbseq_xml_record("1001", "GENERIC_SOURCE"),
                                b"</GBSet>\n",
                            ]
                        ),
                        xml_payload_sha256="deadbeef",
                    )
                ],
            )

            xml_snapshot_directory = save_ncbi_protein_xml_snapshot(
                fetch_result=fetch_result,
                snapshot_root_directory=temporary_directory / "xml_snapshot",
                source_uid_snapshot_manifest=source_identifier_snapshot_manifest,
                source_uid_snapshot_manifest_file_path=(
                    source_identifier_snapshot_manifest_file_path
                ),
                protein_uids=protein_uids,
            )
            xml_manifest = json.loads(
                (xml_snapshot_directory / "manifest.json").read_text(
                    encoding="utf-8",
                )
            )

        self.assertIsNone(xml_manifest["search_query"])
        self.assertEqual(
            xml_manifest["source_identifier_snapshot_artifact_type"],
            "ncbi_protein_identifier_snapshot",
        )
        self.assertEqual(
            xml_manifest["source_identifier_snapshot_source_type"],
            "provided_identifier_list",
        )
        for generic_key, legacy_key in (
            ("requested_protein_identifier_count", "requested_protein_uid_count"),
            ("normalized_protein_identifier_count", "normalized_protein_uid_count"),
            ("protein_identifiers_sha256", "protein_uids_sha256"),
        ):
            self.assertEqual(xml_manifest[generic_key], xml_manifest[legacy_key])
        self.assertEqual(xml_manifest["protein_identifiers_file_name"], "protein_identifiers.txt")
        self.assertEqual(xml_manifest["protein_uids_file_name"], "protein_uids.txt")
        self.assertEqual(xml_manifest["snapshot_format_version"], "1.1")
        self.assertEqual(xml_manifest["source_type"], "identifier_snapshot")
        self.assertEqual(xml_manifest["protein_identifier_count"], 1)
        self.assertEqual(xml_manifest["source_identifier_count"], 1)
        self.assertEqual(xml_manifest["source_identifier_sha256"], protein_uids_sha256)
        self.assertEqual(xml_manifest["batches"][0]["protein_identifier_count"], 1)
        self.assertEqual(xml_manifest["batches"][0]["protein_uid_count"], 1)

    def test_xml_manifest_exposes_generic_metadata_for_both_identifier_domains(self) -> None:
        # Metadata construction remains independent of XML record validation.
        for identifier_type, identifiers in (
            ("uid", ["1001", "1002"]),
            ("accession.version", ["WP_012345678.1", "XP_123456789.2"]),
        ):
            identifier_hash = sha256_of_lines(text_lines=identifiers)
            result = NCBIProteinXmlFetchResult(
                database_name="protein",
                identifier_type=identifier_type,
                retrieved_at_utc="2026-10-09T12:00:00Z",
                requested_protein_uid_count=3,
                normalized_protein_uid_count=2,
                protein_uids_sha256=identifier_hash,
                batch_size=100,
                batch_count=1,
                retmode="xml",
                request_delay_seconds=0.0,
                max_retry_attempts=5,
                python_version="3.14.0",
                biopython_version="test",
                xml_batches=[],
            )
            sources = [
                {"normalized_protein_uid_count": 2, "protein_uids_sha256": identifier_hash},
                {"normalized_protein_identifier_count": 2,
                 "protein_identifiers_sha256": identifier_hash},
                {"protein_identifier_count": 2, "protein_identifiers_sha256": identifier_hash},
            ]
            for source in sources:
                with self.subTest(identifier_type=identifier_type, source=source):
                    manifest = _build_xml_snapshot_manifest(
                        fetch_result=result,
                        immutable_snapshot_directory_name="xml_snapshot",
                        immutable_snapshot_relative_path="snapshots/xml_snapshot",
                        source_uid_snapshot_manifest=source,
                        source_uid_snapshot_manifest_sha256="source-manifest-hash",
                        consolidated_xml_file_name="protein_records.xml",
                        consolidated_xml_file_sha256="xml-hash",
                        consolidated_record_count=2,
                        batch_metadata_records=[],
                    )
                    self.assertEqual(manifest["identifier_type"], identifier_type)
                    self.assertEqual(manifest["requested_protein_identifier_count"], 3)
                    self.assertEqual(manifest["normalized_protein_identifier_count"], 2)
                    self.assertEqual(manifest["protein_identifier_count"], 2)
                    self.assertEqual(manifest["protein_identifiers_sha256"], identifier_hash)
                    self.assertEqual(manifest["protein_identifiers_file_name"], "protein_identifiers.txt")
                    self.assertEqual(manifest["requested_protein_uid_count"], 3)
                    self.assertEqual(manifest["normalized_protein_uid_count"], 2)
                    self.assertEqual(manifest["protein_uids_sha256"], identifier_hash)
                    self.assertEqual(manifest["protein_uids_file_name"], "protein_uids.txt")
                    self.assertEqual(manifest["source_identifier_count"], 2)
                    self.assertEqual(manifest["source_identifier_sha256"], identifier_hash)


class SaveNcbiProteinXmlIdentifierSnapshotTests(unittest.TestCase):
    def _fetch_result(
        self,
        identifiers: list[str],
        identifier_type: str,
        records: list[str] | None = None,
    ) -> NCBIProteinXmlFetchResult:
        if records is None:
            records = [
                f"<GBSeq><GBSeqid>gi|{identifier}</GBSeqid></GBSeq>"
                if identifier_type == "uid" else
                f"<GBSeq><GBSeq_accession-version>{identifier}</GBSeq_accession-version>"
                "<GBSeqid>gi|12345</GBSeqid></GBSeq>"
                for identifier in identifiers
            ]
        batches = []
        for start in range(0, len(records), 2):
            payload = ("<GBSet>" + "".join(records[start:start + 2]) + "</GBSet>").encode("utf-8")
            batches.append(NCBIProteinXmlBatchFetchResult(
                batch_index=len(batches) + 1,
                batch_start_index=start,
                batch_end_index=start + len(records[start:start + 2]) - 1,
                protein_uids=identifiers[start:start + 2],
                protein_uid_count=len(identifiers[start:start + 2]),
                xml_payload_bytes=payload,
                xml_payload_sha256=hashlib.sha256(payload).hexdigest(),
            ))
        return NCBIProteinXmlFetchResult(
            database_name="protein",
            identifier_type=identifier_type,
            retrieved_at_utc="2026-10-16T12:00:00Z",
            requested_protein_uid_count=len(identifiers),
            normalized_protein_uid_count=len(identifiers),
            protein_uids_sha256=sha256_of_lines(
                text_lines=identifiers,
                deduplicate_lines_preserving_order=False,
                sort_lines=False,
            ),
            batch_size=2,
            batch_count=len(batches),
            retmode="xml",
            request_delay_seconds=0.0,
            max_retry_attempts=5,
            python_version="3.14.0",
            biopython_version="test",
            xml_batches=batches,
        )

    def _source_manifest(self, result: NCBIProteinXmlFetchResult) -> dict:
        return {
            "identifier_type": result.identifier_type,
            "source_type": "provided_identifier_list",
            "normalized_protein_identifier_count": result.protein_identifier_count,
            "protein_identifiers_sha256": result.protein_identifiers_sha256,
            "immutable_snapshot_relative_path": "snapshots/source_identifiers",
            "immutable_snapshot_directory_name": "source_identifiers",
        }

    def _save(self, directory: Path, result, identifiers, source) -> Path:
        source_file = directory / "source_manifest.json"
        source_file.write_text(json.dumps(source, indent=2) + "\n", encoding="utf-8")
        return save_ncbi_protein_xml_snapshot(
            fetch_result=result,
            snapshot_root_directory=directory / "xml_snapshot",
            source_uid_snapshot_manifest=source,
            source_uid_snapshot_manifest_file_path=source_file,
            protein_uids=identifiers,
        )

    def test_direct_xml_save_requires_canonical_source_integrity_fields_in_1_1(self) -> None:
        for identifier_type, identifier in (
            ("uid", "1001"), ("accession.version", "WP_012345678.1"),
        ):
            result = self._fetch_result([identifier], identifier_type)
            for field in ("protein_identifiers_sha256", "normalized_protein_identifier_count"):
                with self.subTest(identifier_type=identifier_type, field=field), tempfile.TemporaryDirectory() as name:
                    directory = Path(name)
                    source = self._source_manifest(result)
                    source.update({
                        "snapshot_format_version": "1.1",
                        "protein_uids_sha256": result.protein_identifiers_sha256,
                        "normalized_protein_uid_count": 1,
                    })
                    del source[field]
                    with self.assertRaisesRegex(RuntimeError, f"{field} is required"):
                        self._save(directory, result, [identifier], source)
                    self.assertFalse((directory / "xml_snapshot").exists())

    def test_direct_xml_save_checks_supplied_provenance_hash_against_manifest_bytes(self) -> None:
        for identifier_type, identifier in (
            ("uid", "1001"), ("accession.version", "WP_012345678.1"),
        ):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result([identifier], identifier_type)
                source = self._source_manifest(result)
                source_file = directory / "source_manifest.json"
                arguments = {
                    "fetch_result": result,
                    "snapshot_root_directory": directory / "xml_snapshot",
                    "source_uid_snapshot_manifest": source,
                    "source_uid_snapshot_manifest_file_path": source_file,
                    "source_uid_snapshot_manifest_sha256": "0" * 64,
                    "protein_uids": [identifier],
                }
                with self.assertRaises(FileNotFoundError):
                    save_ncbi_protein_xml_snapshot(**arguments)
                source_file.write_text(json.dumps(source), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "SHA-256 does not match"):
                    save_ncbi_protein_xml_snapshot(**arguments)
                self.assertFalse((directory / "xml_snapshot").exists())

    def test_xml_provenance_is_pinned_before_latest_changes_during_fetch(self) -> None:
        for identifier_type, identifiers, replacement in (
            ("uid", ["1001", "1002"], ["1003", "1004"]),
            ("accession.version", ["WP_012345678.1", "XP_123456789.2"], ["WP_012345678.2", "XP_123456789.3"]),
        ):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                source_root = directory / "identifiers"
                original = save_ncbi_protein_identifier_snapshot(
                    protein_identifiers=identifiers, identifier_type=identifier_type,
                    snapshot_root_directory=source_root, retrieved_at_utc="2026-10-15T12:00:00Z",
                    sort_identifiers=False,
                )
                original_manifest = original / "manifest.json"
                latest_manifest = source_root / "latest/manifest.json"
                expected_hash = hashlib.sha256(original_manifest.read_bytes()).hexdigest()
                original_open = Path.open
                fetch_started = False

                def open_file(path, *args, **kwargs):
                    if fetch_started and path in (original_manifest, latest_manifest):
                        self.fail("Source manifest was read after starting the fetch.")
                    return original_open(path, *args, **kwargs)

                def fetch_from_original_source(**kwargs):
                    nonlocal fetch_started
                    fetch_started = True
                    self.assertEqual(kwargs["protein_uids"], identifiers)
                    self.assertEqual(kwargs["identifier_type"], identifier_type)
                    save_ncbi_protein_identifier_snapshot(
                        protein_identifiers=replacement, identifier_type=identifier_type,
                        snapshot_root_directory=source_root, retrieved_at_utc="2026-10-15T12:00:01Z",
                    )
                    return self._fetch_result(identifiers, identifier_type)

                with patch.object(Path, "open", open_file), patch.object(
                    ncbi_snapshot, "fetch_ncbi_protein_xml_batches", side_effect=fetch_from_original_source,
                ) as fetch_mock:
                    published = resolve_ncbi_protein_xml_snapshot(
                        snapshot_mode="create_new", snapshot_root_directory=directory / "xml_snapshot",
                        source_identifier_snapshot_root_directory=source_root, ncbi_email="test@example.org",
                    )
                fetch_mock.assert_called_once()
                self.assertNotEqual(hashlib.sha256(latest_manifest.read_bytes()).hexdigest(), expected_hash)
                self.assertEqual(published["protein_uids"], identifiers)
                manifest = published["manifest"]
                self.assertEqual(manifest["source_identifier_snapshot_directory_name"], original.name)
                self.assertEqual(manifest["source_identifier_snapshot_manifest_sha256"], expected_hash)
                self.assertEqual(manifest["source_uid_snapshot_manifest_sha256"], expected_hash)

    def test_xml_resolver_rejects_latest_manifest_divergence_from_immutable_source(self) -> None:
        for identifier_type, identifiers in (
            ("uid", ["1001", "1002"]),
            ("accession.version", ["WP_012345678.1", "XP_123456789.2"]),
        ):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                source_root = self._frozen_resolver_source(directory, identifiers, identifier_type)
                latest_manifest = source_root / "latest/manifest.json"
                manifest = json.loads(latest_manifest.read_bytes())
                manifest["source_description"] = "Altered only in latest"
                latest_manifest.write_text(json.dumps(manifest), encoding="utf-8")
                xml_root = directory / "xml_snapshot"
                with patch.object(ncbi_snapshot, "fetch_ncbi_protein_xml_batches") as fetch_mock:
                    with self.assertRaisesRegex(RuntimeError, "differs from its immutable snapshot"):
                        resolve_ncbi_protein_xml_snapshot(
                            snapshot_mode="create_new", snapshot_root_directory=xml_root,
                            source_identifier_snapshot_root_directory=source_root, ncbi_email="test@example.org",
                        )
                    fetch_mock.assert_not_called()
                self.assertFalse(xml_root.exists())

    def test_create_new_xml_can_replace_latest_with_an_invalid_identifier_path(self) -> None:
        for identifier_type, identifiers in (
            ("uid", ["1001"]), ("accession.version", ["WP_012345678.1"]),
        ):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                source_root = self._frozen_resolver_source(directory, identifiers, identifier_type)
                result = self._fetch_result(identifiers, identifier_type)
                self._save(directory, result, identifiers, self._source_manifest(result))
                xml_root = directory / "xml_snapshot"
                latest_manifest = xml_root / "latest/manifest.json"
                manifest = json.loads(latest_manifest.read_bytes())
                manifest["protein_identifiers_file_name"] = "../outside.txt"
                latest_manifest.write_text(json.dumps(manifest), encoding="utf-8")
                result = replace(result, retrieved_at_utc="2026-10-16T12:00:01Z")
                with patch.object(ncbi_snapshot, "fetch_ncbi_protein_xml_batches", return_value=result):
                    published = resolve_ncbi_protein_xml_snapshot(
                        snapshot_mode="create_new", snapshot_root_directory=xml_root,
                        source_identifier_snapshot_root_directory=source_root, ncbi_email="test@example.org",
                    )
                self.assertEqual(published["protein_uids"], identifiers)
                self.assertEqual(self._reuse(directory)["protein_uids"], identifiers)

    def test_invalid_retrieval_timestamps_are_rejected_before_creating_snapshots(self) -> None:
        for identifier_type, identifiers in (
            ("uid", ["1001"]), ("accession.version", ["WP_012345678.1"]),
        ):
            for timestamp in ("../../unexpected", "2026-02-30T12:00:00Z", "2026-10-16T12:00:00", "2026-10-16T12:00:00+00:99"):
                with (
                    self.subTest(identifier_type=identifier_type, timestamp=timestamp),
                    tempfile.TemporaryDirectory() as name,
                ):
                    directory = Path(name)
                    with self.assertRaisesRegex(ValueError, "retrieved_at_utc"):
                        save_ncbi_protein_identifier_snapshot(
                            protein_identifiers=identifiers, identifier_type=identifier_type,
                            snapshot_root_directory=directory / "identifiers", retrieved_at_utc=timestamp,
                        )
                    result = replace(self._fetch_result(identifiers, identifier_type), retrieved_at_utc=timestamp)
                    with self.assertRaisesRegex(ValueError, "retrieved_at_utc"):
                        self._save(directory, result, identifiers, self._source_manifest(result))
                    self.assertFalse((directory / "identifiers").exists())
                    self.assertFalse((directory / "xml_snapshot").exists())
                    self.assertEqual([path.name for path in directory.iterdir()], ["source_manifest.json"])

    def test_retrieval_timestamp_offsets_produce_canonical_utc_directory_names(self) -> None:
        for identifier_type, identifiers in (
            ("uid", ["1001"]), ("accession.version", ["WP_012345678.1"]),
        ):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                timestamp = "2026-10-16T09:00:00.125-03:00"
                identifiers_directory = save_ncbi_protein_identifier_snapshot(
                    protein_identifiers=identifiers, identifier_type=identifier_type,
                    snapshot_root_directory=directory / "identifiers", retrieved_at_utc=timestamp,
                )
                result = replace(self._fetch_result(identifiers, identifier_type), retrieved_at_utc=timestamp)
                xml_directory = self._save(directory, result, identifiers, self._source_manifest(result))
                for snapshot in (identifiers_directory, xml_directory):
                    self.assertTrue(snapshot.name.startswith("2026-10-16T12-00-00.125000Z__q_"))

    def test_xml_resolver_checks_source_file_hash_before_fetch(self) -> None:
        for identifier_type, identifiers in (
            ("uid", ["1001", "1002"]),
            ("accession.version", ["WP_012345678.1", "XP_123456789.2"]),
        ):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                source_root = self._frozen_resolver_source(directory, identifiers, identifier_type)
                (source_root / "latest/protein_identifiers.txt").write_text(
                    "\n".join(reversed(identifiers)) + "\n", encoding="utf-8",
                )
                with patch.object(
                    ncbi_snapshot, "fetch_ncbi_protein_xml_batches",
                    return_value=self._fetch_result(list(reversed(identifiers)), identifier_type),
                ) as fetch_mock:
                    with self.assertRaisesRegex(RuntimeError, "identifier hash mismatch"):
                        resolve_ncbi_protein_xml_snapshot(
                            snapshot_mode="create_new",
                            snapshot_root_directory=directory / "xml_snapshot",
                            source_identifier_snapshot_root_directory=source_root,
                            ncbi_email="test@example.org",
                        )
                    fetch_mock.assert_not_called()
                self.assertFalse((directory / "xml_snapshot").exists())

    def test_xml_resolver_checks_source_count_before_fetch(self) -> None:
        count_fields = (
            "normalized_protein_identifier_count",
            "protein_identifier_count",
            "normalized_protein_uid_count",
        )
        for identifier_type, identifiers in (
            ("uid", ["1001", "1002"]),
            ("accession.version", ["WP_012345678.1", "XP_123456789.2"]),
        ):
            for count_field in count_fields:
                for invalid_count in (999, None, True, "2"):
                    with (
                        self.subTest(identifier_type=identifier_type, field=count_field, count=invalid_count),
                        tempfile.TemporaryDirectory() as name,
                    ):
                        directory = Path(name)
                        source_root = self._frozen_resolver_source(directory, identifiers, identifier_type)
                        manifest_path = source_root / "latest/manifest.json"
                        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                        for preceding_field in count_fields[:count_fields.index(count_field)]:
                            manifest.pop(preceding_field)
                        manifest["snapshot_format_version"] = "1.0"
                        manifest[count_field] = invalid_count
                        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
                        with patch.object(
                            ncbi_snapshot, "fetch_ncbi_protein_xml_batches",
                            return_value=self._fetch_result(identifiers, identifier_type),
                        ) as fetch_mock:
                            with self.assertRaisesRegex(RuntimeError, "identifier count|must be an integer"):
                                resolve_ncbi_protein_xml_snapshot(
                                    snapshot_mode="create_new",
                                    snapshot_root_directory=directory / "xml_snapshot",
                                    source_identifier_snapshot_root_directory=source_root,
                                    ncbi_email="test@example.org",
                                )
                            fetch_mock.assert_not_called()
                        self.assertFalse((directory / "xml_snapshot").exists())

    def test_xml_save_rejects_reordered_list_without_replacing_latest(self) -> None:
        for identifier_type, identifiers in (
            ("uid", ["1001", "1002"]),
            ("accession.version", ["WP_012345678.1", "XP_123456789.2"]),
        ):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(identifiers, identifier_type)
                source = self._source_manifest(result)
                original = self._save(directory, result, identifiers, source)
                latest = directory / "xml_snapshot/latest"
                original_files = {path.name: path.read_bytes() for path in latest.iterdir()}
                next_result = replace(result, retrieved_at_utc="2026-10-16T12:00:01Z")
                with self.assertRaisesRegex(ValueError, "protein_identifiers_sha256 values differ"):
                    self._save(directory, next_result, list(reversed(identifiers)), source)
                self.assertEqual(list((directory / "xml_snapshot/snapshots").iterdir()), [original])
                self.assertEqual(
                    {path.name: path.read_bytes() for path in latest.iterdir()}, original_files,
                )
                self.assertEqual(self._reuse(directory)["protein_uids"], identifiers)

    def test_xml_save_checks_source_count_and_generic_field_precedence(self) -> None:
        for identifier_type, identifiers in (
            ("uid", ["1001", "1002"]),
            ("accession.version", ["WP_012345678.1", "XP_123456789.2"]),
        ):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(identifiers, identifier_type)
                source = self._source_manifest(result)
                source["normalized_protein_identifier_count"] = 999
                source["normalized_protein_uid_count"] = len(identifiers)
                with self.assertRaisesRegex(ValueError, "protein identifier count"):
                    self._save(directory, result, identifiers, source)
                self.assertFalse((directory / "xml_snapshot").exists())
                source["normalized_protein_identifier_count"] = len(identifiers)
                source["normalized_protein_uid_count"] = 999
                self._save(directory, result, identifiers, source)
                self.assertEqual(self._reuse(directory)["protein_uids"], identifiers)

    def test_xml_save_verifies_written_identifier_files_before_publication(self) -> None:
        write_identifier_files = ncbi_snapshot._write_snapshot_identifier_files
        for identifier_type, identifiers in (
            ("uid", ["1001", "1002"]),
            ("accession.version", ["WP_012345678.1", "XP_123456789.2"]),
        ):
            for file_name in ("protein_identifiers.txt", "protein_uids.txt"):
                for changed_identifiers in (list(reversed(identifiers)), identifiers[:1]):
                    with (
                        self.subTest(identifier_type=identifier_type, file=file_name, values=changed_identifiers),
                        tempfile.TemporaryDirectory() as name,
                    ):
                        directory = Path(name)
                        result = self._fetch_result(identifiers, identifier_type)
                        source = self._source_manifest(result)
                        original = self._save(directory, result, identifiers, source)
                        latest = directory / "xml_snapshot/latest"
                        original_files = {path.name: path.read_bytes() for path in latest.iterdir()}

                        def write_changed_identifier_file(**kwargs):
                            files = write_identifier_files(**kwargs)
                            (kwargs["snapshot_directory"] / file_name).write_text(
                                "\n".join(changed_identifiers) + "\n", encoding="utf-8",
                            )
                            return files

                        with patch.object(
                            ncbi_snapshot, "_write_snapshot_identifier_files",
                            side_effect=write_changed_identifier_file,
                        ):
                            with self.assertRaisesRegex(ValueError, "protein identifier|protein_identifiers_sha256"):
                                self._save(
                                    directory, replace(result, retrieved_at_utc="2026-10-16T12:00:01Z"),
                                    identifiers, source,
                                )
                        self.assertEqual(list((directory / "xml_snapshot/snapshots").iterdir()), [original])
                        self.assertEqual(
                            {path.name: path.read_bytes() for path in latest.iterdir()}, original_files,
                        )

    def test_saves_accession_xml_from_generic_identifier_snapshot_and_updates_latest(self) -> None:
        identifiers = ["XP_123456789.2", "WP_012345678.1", "XP_123456789.2"]
        records = [
            "<GBSeq><GBSeqid>ref|XP_123456789.2|</GBSeqid></GBSeq>",
            "<GBSeq><GBSeq_accession-version>WP_012345678.1</GBSeq_accession-version>"
            "<GBSeqid>gi|12345</GBSeqid></GBSeq>",
            "<GBSeq><GBSeq_other-seqids><GBSeqid>ref|XP_123456789.2|</GBSeqid>"
            "</GBSeq_other-seqids></GBSeq>",
        ]
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            source_directory = save_ncbi_protein_identifier_snapshot(
                protein_identifiers=identifiers,
                identifier_type="accession.version",
                snapshot_root_directory=directory / "identifiers",
                deduplicate_identifiers=False,
                sort_identifiers=False,
                update_latest_directory=False,
            )
            source_payload = load_snapshot_by_directory(snapshot_directory=source_directory)
            result = self._fetch_result(identifiers, "accession.version", records)
            snapshot = self._save(directory, result, source_payload["protein_uids"], source_payload["manifest"])
            manifest = json.loads((snapshot / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["identifier_type"], "accession.version")
            self.assertEqual(manifest["normalized_protein_identifier_count"], 3)
            self.assertEqual(manifest["requested_protein_identifier_count"], 3)
            self.assertEqual(manifest["protein_identifier_count"], 3)
            self.assertEqual(manifest["protein_identifiers_sha256"], result.protein_identifiers_sha256)
            self.assertEqual(manifest["protein_identifiers_file_name"], "protein_identifiers.txt")
            self.assertEqual(manifest["snapshot_format_version"], "1.1")
            self.assertEqual(manifest["source_type"], "identifier_snapshot")
            self.assertEqual(manifest["protein_uids_sha256"], manifest["protein_identifiers_sha256"])
            self.assertEqual(manifest["consolidated_record_count"], 3)
            self.assertEqual([batch["protein_identifier_count"] for batch in manifest["batches"]], [2, 1])
            self.assertEqual(manifest["source_identifier_count"], 3)
            self.assertEqual(manifest["source_identifier_sha256"], result.protein_identifiers_sha256)
            for suffix in ("relative_path", "directory_name", "manifest_sha256"):
                self.assertEqual(
                    manifest[f"source_identifier_snapshot_{suffix}"],
                    manifest[f"source_uid_snapshot_{suffix}"],
                )
            self.assertEqual(
                manifest["source_identifier_snapshot_manifest_sha256"],
                sha256_of_file(input_file_path=directory / "source_manifest.json"),
            )
            self.assertEqual(
                manifest["xml_file_sha256"],
                sha256_of_file(input_file_path=snapshot / "protein_records.xml"),
            )
            self.assertEqual((snapshot / "protein_uids.txt").read_text(encoding="utf-8"), "\n".join(identifiers) + "\n")
            for filename in ("manifest.json", "protein_identifiers.txt", "protein_uids.txt", "protein_records.xml"):
                self.assertEqual(
                    (directory / "xml_snapshot/latest" / filename).read_bytes(),
                    (snapshot / filename).read_bytes(),
                )

    def test_saves_uid_xml_with_explicit_or_legacy_source_type(self) -> None:
        result = self._fetch_result(["1001"], "uid")
        for legacy in (False, True):
            with self.subTest(legacy=legacy), tempfile.TemporaryDirectory() as name:
                source = self._source_manifest(result)
                if legacy:
                    source.pop("identifier_type")
                    source["protein_uids_sha256"] = source.pop("protein_identifiers_sha256")
                    source["normalized_protein_uid_count"] = source.pop("normalized_protein_identifier_count")
                snapshot = self._save(Path(name), result, ["1001"], source)
                manifest = json.loads((snapshot / "manifest.json").read_text(encoding="utf-8"))
                self.assertEqual(manifest["identifier_type"], "uid")
                self.assertEqual(manifest["source_identifier_sha256"], result.protein_identifiers_sha256)
                self.assertEqual(manifest["source_identifier_count"], 1)

    def test_rejects_type_mismatch_before_creating_snapshot(self) -> None:
        cases = [
            ("uid", "accession.version", ["1001"]),
            ("accession.version", "uid", ["WP_012345678.1"]),
            ("accession.version", None, ["WP_012345678.1"]),
        ]
        for fetch_type, source_type, identifiers in cases:
            with self.subTest(fetch_type=fetch_type, source_type=source_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(identifiers, fetch_type)
                source = self._source_manifest(result)
                if source_type is None:
                    source.pop("identifier_type")
                else:
                    source["identifier_type"] = source_type
                with self.assertRaisesRegex(ValueError, "identifier_type does not match"):
                    self._save(directory, result, identifiers, source)
                self.assertFalse((directory / "xml_snapshot").exists())

    def test_rejects_explicit_invalid_types_instead_of_defaulting_to_uid(self) -> None:
        for origin in ("fetch", "source"):
            for value in (None, "", "protein_accession"):
                with self.subTest(origin=origin, value=value), tempfile.TemporaryDirectory() as name:
                    directory = Path(name)
                    result = self._fetch_result(["1001"], "uid")
                    source = self._source_manifest(result)
                    if origin == "source":
                        source["identifier_type"] = value
                    else:
                        result = replace(result, identifier_type=value)
                    with self.assertRaisesRegex(ValueError, "identifier_type"):
                        self._save(directory, result, ["1001"], source)
                    self.assertFalse((directory / "xml_snapshot").exists())

    def test_type_mismatch_preserves_existing_latest_xml_snapshot(self) -> None:
        for source_type, fetch_type, identifiers in (
            ("uid", "accession.version", ["1001"]),
            ("accession.version", "uid", ["WP_012345678.1"]),
        ):
            with self.subTest(source_type=source_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(identifiers, source_type)
                source = self._source_manifest(result)
                original = self._save(directory, result, identifiers, source)
                latest = directory / "xml_snapshot/latest"
                original_files = {
                    path.name: path.read_bytes() for path in latest.iterdir() if path.is_file()
                }
                mismatched_result = replace(
                    result,
                    identifier_type=fetch_type,
                    retrieved_at_utc="2026-10-16T12:00:01Z",
                )
                with self.assertRaisesRegex(
                    ValueError, f"Expected {source_type}, got {fetch_type}",
                ):
                    self._save(directory, mismatched_result, identifiers, source)
                self.assertEqual(list((directory / "xml_snapshot/snapshots").iterdir()), [original])
                self.assertEqual(
                    {path.name: path.read_bytes() for path in latest.iterdir() if path.is_file()},
                    original_files,
                )

    def test_generic_source_hash_is_checked_with_precedence_over_legacy_alias(self) -> None:
        result = self._fetch_result(["WP_012345678.1"], "accession.version")
        for valid_generic_hash in (False, True):
            with self.subTest(valid_generic_hash=valid_generic_hash), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                source = self._source_manifest(result)
                if valid_generic_hash:
                    source["protein_uids_sha256"] = "incorrect-legacy-hash"
                    self._save(directory, result, ["WP_012345678.1"], source)
                else:
                    source["protein_uids_sha256"] = result.protein_identifiers_sha256
                    source["protein_identifiers_sha256"] = "0" * 64
                    with self.assertRaisesRegex(ValueError, "protein_identifiers_sha256 values differ"):
                        self._save(directory, result, ["WP_012345678.1"], source)
                    self.assertFalse((directory / "xml_snapshot").exists())

    def test_generic_only_source_hash_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            result = self._fetch_result(["WP_012345678.1"], "accession.version")
            source = self._source_manifest(result)
            source["protein_identifiers_sha256"] = "0" * 64
            with self.assertRaisesRegex(ValueError, "protein_identifiers_sha256 values differ"):
                self._save(directory, result, ["WP_012345678.1"], source)
            self.assertFalse((directory / "xml_snapshot").exists())

    def test_invalid_selected_source_hash_is_rejected_before_creating_xml_snapshot(self) -> None:
        for identifier_type, identifiers in (
            ("uid", ["1001"]),
            ("accession.version", ["WP_012345678.1"]),
        ):
            for hash_field in ("protein_identifiers_sha256", "protein_uids_sha256"):
                for invalid_hash in (None, "", " \t", 123, False, [], {}, "invalid", "0" * 63, "g" * 64):
                    with (
                        self.subTest(identifier_type=identifier_type, hash_field=hash_field, value=invalid_hash),
                        tempfile.TemporaryDirectory() as name,
                    ):
                        directory = Path(name)
                        result = self._fetch_result(identifiers, identifier_type)
                        source = self._source_manifest(result)
                        source["protein_uids_sha256"] = result.protein_identifiers_sha256
                        if hash_field == "protein_uids_sha256":
                            source.pop("protein_identifiers_sha256")
                        source[hash_field] = invalid_hash
                        with self.assertRaisesRegex(ValueError, f"{hash_field} must be a SHA-256 hash"):
                            self._save(directory, result, identifiers, source)
                        self.assertFalse((directory / "xml_snapshot").exists())

    def test_source_manifest_without_either_identifier_hash_is_rejected(self) -> None:
        for identifier_type, identifiers in (
            ("uid", ["1001"]),
            ("accession.version", ["WP_012345678.1"]),
        ):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(identifiers, identifier_type)
                source = self._source_manifest(result)
                source.pop("protein_identifiers_sha256")
                with self.assertRaisesRegex(
                    ValueError, "must define protein_identifiers_sha256 or protein_uids_sha256",
                ):
                    self._save(directory, result, identifiers, source)
                self.assertFalse((directory / "xml_snapshot").exists())

    def test_valid_legacy_source_hash_is_used_when_generic_hash_is_absent(self) -> None:
        for identifier_type, identifiers in (
            ("uid", ["1001"]),
            ("accession.version", ["WP_012345678.1"]),
        ):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(identifiers, identifier_type)
                source = self._source_manifest(result)
                source["protein_uids_sha256"] = source.pop("protein_identifiers_sha256")
                snapshot = self._save(directory, result, identifiers, source)
                manifest = json.loads((snapshot / "manifest.json").read_text(encoding="utf-8"))
                self.assertEqual(manifest["source_identifier_sha256"], result.protein_identifiers_sha256)
                self.assertEqual(self._reuse(directory)["protein_uids"], identifiers)

    def test_null_generic_source_hash_cannot_publish_changed_revision_or_replace_latest(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            identifiers = ["WP_012345678.1"]
            original_result = self._fetch_result(identifiers, "accession.version")
            original = self._save(directory, original_result, identifiers, self._source_manifest(original_result))
            latest = directory / "xml_snapshot/latest"
            original_files = {
                path.name: path.read_bytes() for path in latest.iterdir() if path.is_file()
            }
            source_root = self._frozen_resolver_source(directory, identifiers, "accession.version")
            source_latest = source_root / "latest"
            manifest_path = source_latest / "manifest.json"
            source = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(source["protein_uids_sha256"], original_result.protein_identifiers_sha256)
            source["protein_identifiers_sha256"] = None
            manifest_path.write_text(json.dumps(source), encoding="utf-8")
            changed_identifiers = ["WP_012345678.2"]
            (source_latest / "protein_identifiers.txt").write_text(
                "\n".join(changed_identifiers) + "\n", encoding="utf-8",
            )
            changed_result = replace(
                self._fetch_result(changed_identifiers, "accession.version"),
                retrieved_at_utc="2026-10-16T12:00:01Z",
            )
            with patch(
                "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches", return_value=changed_result,
            ) as fetch_mock:
                with self.assertRaisesRegex(RuntimeError, "protein_identifiers_sha256 must be a non-empty string"):
                    resolve_ncbi_protein_xml_snapshot(
                        snapshot_mode="create_new",
                        snapshot_root_directory=directory / "xml_snapshot",
                        source_identifier_snapshot_root_directory=source_root,
                        ncbi_email="test@example.org",
                    )
            fetch_mock.assert_not_called()
            self.assertEqual(list((directory / "xml_snapshot/snapshots").iterdir()), [original])
            self.assertEqual(
                {path.name: path.read_bytes() for path in latest.iterdir() if path.is_file()},
                original_files,
            )

    def test_invalid_accession_xml_is_not_published_and_previous_latest_is_preserved(self) -> None:
        invalid_records = [
            "<GBSeq><GBSeq_accession-version>WP_012345678.2</GBSeq_accession-version></GBSeq>",
            "<GBSeq><GBSeq_primary-accession>WP_012345678</GBSeq_primary-accession></GBSeq>",
            "<GBSeq><GBSeqid>ref|WP_012345678.1|</GBSeqid>"
            "<GBSeqid>ref|XP_123456789.2|</GBSeqid></GBSeq>",
        ]
        for record in invalid_records:
            with self.subTest(record=record), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                identifiers = ["WP_012345678.1"]
                result = self._fetch_result(identifiers, "accession.version")
                source = self._source_manifest(result)
                original = self._save(directory, result, identifiers, source)
                latest = directory / "xml_snapshot/latest"
                original_files = {filename: (latest / filename).read_bytes() for filename in
                                  ("manifest.json", "protein_uids.txt", "protein_records.xml")}
                invalid_result = replace(
                    self._fetch_result(identifiers, "accession.version", [record]),
                    retrieved_at_utc="2026-10-16T12:00:01Z",
                )
                with self.assertRaises(RuntimeError):
                    self._save(directory, invalid_result, identifiers, source)
                self.assertEqual(list((directory / "xml_snapshot/snapshots").iterdir()), [original])
                for filename, content in original_files.items():
                    self.assertEqual((latest / filename).read_bytes(), content)

    def _reuse(self, directory: Path, mode: str = "reuse_latest") -> dict:
        with patch(
            "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches",
            side_effect=AssertionError("Reuse must not fetch XML from NCBI."),
        ) as fetch_mock, patch(
            "src.pago_pipeline.ncbi_snapshot.resolve_ncbi_protein_uid_snapshot",
            side_effect=AssertionError("Reuse must not load or fetch a source snapshot."),
        ) as source_mock:
            payload = resolve_ncbi_protein_xml_snapshot(
                snapshot_mode=mode,
                snapshot_root_directory=directory / "xml_snapshot",
            )
        fetch_mock.assert_not_called()
        source_mock.assert_not_called()
        return payload

    def test_reuses_accession_xml_from_latest_and_immutable_directory_without_fetch(self) -> None:
        identifiers = ["XP_123456789.2", "WP_012345678.1", "XP_123456789.2"]
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            result = self._fetch_result(identifiers, "accession.version")
            snapshot = self._save(directory, result, identifiers, self._source_manifest(result))
            immutable = load_xml_snapshot_by_directory(snapshot_directory=snapshot)
            self.assertEqual(immutable["protein_uids"], identifiers)
            for mode in ("reuse_latest", "reuse_latest_or_create"):
                with self.subTest(mode=mode):
                    reused = self._reuse(directory, mode)
                    self.assertEqual(reused["protein_uids"], identifiers)
                    self.assertEqual(reused["manifest"]["identifier_type"], "accession.version")
                    self.assertEqual(
                        reused["xml_file_path"],
                        directory / "xml_snapshot/latest/protein_records.xml",
                    )
                    self.assertEqual(reused["xml_file_path"].read_bytes(), immutable["xml_file_path"].read_bytes())

    def test_reuse_prioritizes_generic_metadata_over_conflicting_legacy_fields(self) -> None:
        identifiers = ["WP_012345678.1", "XP_123456789.2"]
        for count_field in ("normalized_protein_identifier_count", "protein_identifier_count"):
            with self.subTest(count_field=count_field), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(identifiers, "accession.version")
                self._save(directory, result, identifiers, self._source_manifest(result))
                manifest_path = directory / "xml_snapshot/latest/manifest.json"
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                if count_field == "protein_identifier_count":
                    manifest["snapshot_format_version"] = "1.0"
                    manifest.pop("normalized_protein_identifier_count")
                else:
                    manifest["protein_identifier_count"] = 999
                manifest["normalized_protein_uid_count"] = 999
                manifest["protein_uids_sha256"] = "wrong-legacy-hash"
                manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
                self.assertEqual(self._reuse(directory)["protein_uids"], identifiers)

    def test_accession_reuse_selects_generic_and_legacy_fields_independently(self) -> None:
        identifiers = ["WP_012345678.1", "WP_999999999.1"]
        for absent_fields in (
            (),
            ("protein_identifiers_file_name",),
            ("protein_identifiers_sha256",),
            ("normalized_protein_identifier_count", "protein_identifier_count"),
            (
                "protein_identifiers_file_name", "protein_identifiers_sha256",
                "normalized_protein_identifier_count", "protein_identifier_count",
            ),
        ):
            with self.subTest(absent_fields=absent_fields), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(identifiers, "accession.version")
                self._save(directory, result, identifiers, self._source_manifest(result))
                latest = directory / "xml_snapshot/latest"
                manifest_path = latest / "manifest.json"
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                manifest["snapshot_format_version"] = "1.0"
                for field in absent_fields:
                    manifest.pop(field)
                if "protein_identifiers_file_name" in manifest:
                    manifest["protein_uids_file_name"] = "unselected_legacy_file.txt"
                    expected_filename = "protein_identifiers.txt"
                else:
                    (latest / "protein_identifiers.txt").unlink()
                    expected_filename = "protein_uids.txt"
                if "protein_identifiers_sha256" in manifest:
                    manifest["protein_uids_sha256"] = "unselected-legacy-hash"
                if "normalized_protein_identifier_count" in manifest:
                    manifest["normalized_protein_uid_count"] = 999
                manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
                original_manifest_bytes = manifest_path.read_bytes()
                for mode in ("reuse_latest", "reuse_latest_or_create"):
                    with self.subTest(mode=mode):
                        reused = self._reuse(directory, mode)
                        self.assertEqual(reused["protein_uids"], identifiers)
                        self.assertEqual(reused["manifest"]["identifier_type"], "accession.version")
                        self.assertEqual(reused["protein_uids_file_path"], latest / expected_filename)
                self.assertEqual(manifest_path.read_bytes(), original_manifest_bytes)

    def test_accession_reuse_checks_legacy_hash_and_count_when_generic_fields_are_absent(self) -> None:
        for absent_fields, legacy_field, invalid_value, message in (
            (
                ("protein_identifiers_sha256",),
                "protein_uids_sha256", "wrong-legacy-hash", "identifier hash mismatch",
            ),
            (
                ("normalized_protein_identifier_count", "protein_identifier_count"),
                "normalized_protein_uid_count", 2, "identifier count mismatch",
            ),
        ):
            with self.subTest(legacy_field=legacy_field), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                identifiers = ["WP_012345678.1"]
                result = self._fetch_result(identifiers, "accession.version")
                self._save(directory, result, identifiers, self._source_manifest(result))
                manifest_path = directory / "xml_snapshot/latest/manifest.json"
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                manifest["snapshot_format_version"] = "1.0"
                for field in absent_fields:
                    manifest.pop(field)
                manifest[legacy_field] = invalid_value
                manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
                original_manifest_bytes = manifest_path.read_bytes()
                for mode in ("reuse_latest", "reuse_latest_or_create"):
                    with self.subTest(mode=mode):
                        with self.assertRaisesRegex(RuntimeError, message):
                            self._reuse(directory, mode)
                self.assertEqual(manifest_path.read_bytes(), original_manifest_bytes)

    def test_reuses_accession_xml_with_generic_only_hash_and_count(self) -> None:
        identifiers = ["WP_012345678.1"]
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            result = self._fetch_result(identifiers, "accession.version")
            self._save(directory, result, identifiers, self._source_manifest(result))
            manifest_path = directory / "xml_snapshot/latest/manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest.pop("normalized_protein_uid_count")
            manifest.pop("protein_uids_sha256")
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            self.assertEqual(self._reuse(directory)["protein_uids"], identifiers)

    def test_reuse_rejects_invalid_generic_fields_instead_of_using_legacy_values(self) -> None:
        cases = [
            ("normalized_protein_identifier_count", 2, "identifier count mismatch"),
            ("normalized_protein_identifier_count", None, "must be an integer"),
            ("normalized_protein_identifier_count", "1", "must be an integer"),
            ("normalized_protein_identifier_count", True, "must be an integer"),
            ("protein_identifier_count", 2, "identifier count mismatch"),
            ("protein_identifier_count", None, "must be an integer"),
            ("protein_identifiers_sha256", "0" * 64, "identifier hash mismatch"),
            ("protein_identifiers_sha256", None, "must be a non-empty string"),
            ("protein_identifiers_sha256", "", "must be a non-empty string"),
            ("protein_identifiers_sha256", 123, "must be a non-empty string"),
        ]
        for field, value, message in cases:
            with self.subTest(field=field, value=value), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(["WP_012345678.1"], "accession.version")
                self._save(directory, result, ["WP_012345678.1"], self._source_manifest(result))
                manifest_path = directory / "xml_snapshot/latest/manifest.json"
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                if field == "protein_identifier_count":
                    manifest["snapshot_format_version"] = "1.0"
                    manifest.pop("normalized_protein_identifier_count")
                manifest[field] = value
                manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
                with self.assertRaisesRegex(RuntimeError, message):
                    self._reuse(directory)

    def test_reuse_defaults_legacy_manifest_to_uid_and_checks_legacy_metadata(self) -> None:
        for corrupt_field in (None, "normalized_protein_uid_count", "protein_uids_sha256"):
            with self.subTest(corrupt_field=corrupt_field), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(["1001"], "uid")
                self._save(directory, result, ["1001"], self._source_manifest(result))
                manifest_path = directory / "xml_snapshot/latest/manifest.json"
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                for field in (
                    "identifier_type", "normalized_protein_identifier_count",
                    "protein_identifier_count", "protein_identifiers_sha256", "protein_identifiers_file_name",
                ):
                    manifest.pop(field)
                manifest["snapshot_format_version"] = "1.0"
                if corrupt_field == "normalized_protein_uid_count":
                    manifest[corrupt_field] = 2
                elif corrupt_field == "protein_uids_sha256":
                    manifest[corrupt_field] = "wrong-legacy-hash"
                manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
                original_manifest_bytes = manifest_path.read_bytes()
                if corrupt_field is None:
                    self.assertEqual(self._reuse(directory)["protein_uids"], ["1001"])
                else:
                    with self.assertRaisesRegex(RuntimeError, "identifier (count|hash) mismatch"):
                        self._reuse(directory)
                self.assertEqual(manifest_path.read_bytes(), original_manifest_bytes)

    def test_reuse_rejects_explicit_invalid_identifier_types(self) -> None:
        for identifier_type in (None, "", "unknown"):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(["1001"], "uid")
                self._save(directory, result, ["1001"], self._source_manifest(result))
                manifest_path = directory / "xml_snapshot/latest/manifest.json"
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                manifest["identifier_type"] = identifier_type
                manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "identifier_type"):
                    self._reuse(directory)

    def test_reuse_rejects_accession_revision_change_even_with_updated_xml_hash(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            result = self._fetch_result(["WP_012345678.1"], "accession.version")
            self._save(directory, result, ["WP_012345678.1"], self._source_manifest(result))
            xml_file = directory / "xml_snapshot/latest/protein_records.xml"
            xml_file.write_bytes(xml_file.read_bytes().replace(b"WP_012345678.1", b"WP_012345678.2"))
            manifest_path = directory / "xml_snapshot/latest/manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["xml_file_sha256"] = sha256_of_file(input_file_path=xml_file)
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            for mode in ("reuse_latest", "reuse_latest_or_create"):
                with self.subTest(mode=mode):
                    with self.assertRaisesRegex(RuntimeError, "accession.version identifiers do not match"):
                        self._reuse(directory, mode)

    def test_reuse_preserves_counter_checks_when_identifier_list_hash_is_updated(self) -> None:
        identifiers = ["WP_012345678.1", "XP_123456789.2", "WP_012345678.1"]
        changed = ["WP_012345678.1", "XP_123456789.2", "XP_123456789.2"]
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            result = self._fetch_result(identifiers, "accession.version")
            self._save(directory, result, identifiers, self._source_manifest(result))
            latest = directory / "xml_snapshot/latest"
            (latest / "protein_identifiers.txt").write_text("\n".join(changed) + "\n", encoding="utf-8")
            manifest_path = latest / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["protein_identifiers_sha256"] = sha256_of_lines(
                text_lines=changed, deduplicate_lines_preserving_order=False, sort_lines=False,
            )
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "accession.version identifiers do not match"):
                self._reuse(directory)

    def test_reuse_retains_xml_record_count_and_file_hash_checks(self) -> None:
        for field, value, message in (
            ("consolidated_record_count", 2, "record count mismatch"),
            ("xml_file_sha256", "0" * 64, "file hash mismatch"),
        ):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(["WP_012345678.1"], "accession.version")
                self._save(directory, result, ["WP_012345678.1"], self._source_manifest(result))
                manifest_path = directory / "xml_snapshot/latest/manifest.json"
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                manifest[field] = value
                manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
                with self.assertRaisesRegex(RuntimeError, message):
                    self._reuse(directory)

    def test_1_1_xml_manifest_reuses_generic_only_files_and_metadata(self) -> None:
        identifiers = ["WP_012345678.1", "XP_123456789.2"]
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            result = self._fetch_result(identifiers, "accession.version")
            snapshot = self._save(directory, result, identifiers, self._source_manifest(result))
            for snapshot_directory in (snapshot, directory / "xml_snapshot/latest"):
                manifest_file = snapshot_directory / "manifest.json"
                manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
                self.assertEqual(manifest["snapshot_format_version"], "1.1")
                self.assertEqual(manifest["source_type"], "identifier_snapshot")
                self.assertEqual(manifest["source_identifier_count"], 2)
                self.assertEqual(manifest["source_identifier_sha256"], result.protein_identifiers_sha256)
                self.assertEqual(manifest["source_identifier_snapshot_relative_path"], "snapshots/source_identifiers")
                self.assertEqual(manifest["source_identifier_snapshot_manifest_sha256"], sha256_of_file(input_file_path=directory / "source_manifest.json"))
                for key in (
                    "protein_uids_file_name", "protein_uids_sha256", "normalized_protein_uid_count",
                    "requested_protein_uid_count", "source_uid_count", "source_uid_sha256",
                    "source_uid_snapshot_relative_path", "source_uid_snapshot_directory_name",
                    "source_uid_snapshot_manifest_sha256",
                ):
                    manifest.pop(key)
                manifest_file.write_text(json.dumps(manifest), encoding="utf-8")
                (snapshot_directory / "protein_uids.txt").unlink()
                self.assertEqual(load_xml_snapshot_by_directory(snapshot_directory=snapshot_directory)["protein_uids"], identifiers)
            self.assertTrue(latest_xml_snapshot_is_available(snapshot_root_directory=directory / "xml_snapshot"))
            for mode in ("reuse_latest", "reuse_latest_or_create"):
                self.assertEqual(self._reuse(directory, mode)["protein_uids"], identifiers)

    def test_1_1_missing_type_is_rejected_in_source_or_reused_xml_manifest(self) -> None:
        for origin in ("source", "xml"):
            with self.subTest(origin=origin), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                result = self._fetch_result(["1001"], "uid")
                source = self._source_manifest(result)
                if origin == "source":
                    source["snapshot_format_version"] = "1.1"
                    source.pop("identifier_type")
                    with self.assertRaisesRegex(ValueError, "identifier_type is required"):
                        self._save(directory, result, ["1001"], source)
                    self.assertFalse((directory / "xml_snapshot").exists())
                else:
                    self._save(directory, result, ["1001"], source)
                    manifest_file = directory / "xml_snapshot/latest/manifest.json"
                    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
                    manifest.pop("identifier_type")
                    manifest_file.write_text(json.dumps(manifest), encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, "identifier_type is required"):
                        self._reuse(directory)

    def test_missing_canonical_identifier_file_does_not_fall_back_to_legacy_copy(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            result = self._fetch_result(["1001"], "uid")
            self._save(directory, result, ["1001"], self._source_manifest(result))
            latest = directory / "xml_snapshot/latest"
            (latest / "protein_identifiers.txt").unlink()
            self.assertTrue((latest / "protein_uids.txt").exists())
            self.assertFalse(latest_xml_snapshot_is_available(snapshot_root_directory=directory / "xml_snapshot"))
            with self.assertRaisesRegex(FileNotFoundError, "protein_identifiers.txt"):
                self._reuse(directory)

    def test_minimal_legacy_source_produces_new_xml_with_explicit_uid_type(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            identifiers = ["1001", "1002"]
            result = self._fetch_result(identifiers, "uid")
            source = {
                "protein_uids_sha256": result.protein_identifiers_sha256,
                "normalized_protein_uid_count": 2,
            }
            snapshot = self._save(directory, result, identifiers, source)
            self.assertEqual(set(source), {"protein_uids_sha256", "normalized_protein_uid_count"})
            manifest = json.loads((snapshot / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["snapshot_format_version"], "1.1")
            self.assertEqual(manifest["identifier_type"], "uid")
            self.assertEqual(manifest["source_identifier_count"], 2)
            self.assertEqual(manifest["source_identifier_sha256"], source["protein_uids_sha256"])
            self.assertEqual(self._reuse(directory)["protein_uids"], identifiers)


    def _frozen_resolver_source(self, directory: Path, identifiers: list[str], identifier_type: str) -> Path:
        root = directory / "identifiers"
        save_ncbi_protein_identifier_snapshot(
            protein_identifiers=identifiers, identifier_type=identifier_type,
            snapshot_root_directory=root, deduplicate_identifiers=False, sort_identifiers=False,
        )
        return root

    def test_accession_snapshot_without_search_query_can_be_created_and_reused(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            identifiers = ["WP_012345678.1", "WP_999999999.1"]
            source_root = self._frozen_resolver_source(directory, identifiers, "accession.version")
            source_manifest_path = source_root / "latest/manifest.json"
            original_source_manifest_bytes = source_manifest_path.read_bytes()
            source_manifest = json.loads(original_source_manifest_bytes)
            self.assertEqual(source_manifest["snapshot_format_version"], "1.1")
            self.assertEqual(source_manifest["identifier_type"], "accession.version")
            self.assertNotIn("search_query", source_manifest)
            records = [
                f"<GBSeq><GBSeq_accession-version>{identifier}</GBSeq_accession-version></GBSeq>"
                for identifier in identifiers
            ]
            result = self._fetch_result(identifiers, "accession.version", records)
            with patch(
                "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches", return_value=result,
            ) as fetch_mock, patch(
                "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_uid_snapshot",
                side_effect=AssertionError("A provided accession list must not trigger ESearch."),
            ) as search_mock:
                created = resolve_ncbi_protein_xml_snapshot(
                    snapshot_mode="create_new",
                    snapshot_root_directory=directory / "xml_snapshot",
                    source_identifier_snapshot_root_directory=source_root,
                    ncbi_email="test@example.org",
                )
            fetch_mock.assert_called_once()
            self.assertEqual(fetch_mock.call_args.kwargs["identifier_type"], "accession.version")
            self.assertEqual(fetch_mock.call_args.kwargs["protein_uids"], identifiers)
            search_mock.assert_not_called()
            self.assertEqual(created["manifest"]["snapshot_format_version"], "1.1")
            self.assertEqual(created["manifest"]["identifier_type"], "accession.version")
            self.assertEqual(created["manifest"]["source_type"], "identifier_snapshot")
            self.assertIsNone(created["manifest"]["search_query"])
            for mode in ("reuse_latest", "reuse_latest_or_create"):
                with self.subTest(mode=mode):
                    reused = self._reuse(directory, mode)
                    self.assertEqual(reused["protein_uids"], identifiers)
                    self.assertEqual(reused["manifest"], created["manifest"])
                    self.assertEqual(reused["xml_file_path"].read_bytes(), created["xml_file_path"].read_bytes())
            self.assertEqual(source_manifest_path.read_bytes(), original_source_manifest_bytes)

    def test_uid_and_accession_workflows_use_generic_sources_without_spreadsheets(self) -> None:
        attempted_spreadsheet_paths = []

        def reject_spreadsheet_open(original_open):
            def guarded_open(file, *args, **kwargs):
                if isinstance(file, (str, bytes, os.PathLike)):
                    path = Path(os.fsdecode(file))
                    if path.suffix.lower() in (".xls", ".xlsx"):
                        attempted_spreadsheet_paths.append(str(path))
                        raise AssertionError("NCBI acquisition must not open a spreadsheet.")
                return original_open(file, *args, **kwargs)

            return guarded_open

        for identifier_type, identifiers in (
            ("uid", ["1001", "1002"]),
            ("accession.version", ["WP_012345678.1", "WP_999999999.1"]),
        ):
            for source_kind in ("list", "text_file"):
                with (
                    self.subTest(identifier_type=identifier_type, source_kind=source_kind),
                    tempfile.TemporaryDirectory() as name,
                ):
                    directory = Path(name)
                    source_root = directory / "identifiers"
                    with patch(
                        "builtins.open", side_effect=reject_spreadsheet_open(builtins.open),
                    ), patch(
                        "io.open", side_effect=reject_spreadsheet_open(io.open),
                    ), patch.dict(sys.modules, {"xlrd": None, "openpyxl": None}):
                        if source_kind == "list":
                            source = save_ncbi_protein_identifier_snapshot(
                                protein_identifiers=identifiers,
                                identifier_type=identifier_type,
                                snapshot_root_directory=source_root,
                            )
                        else:
                            identifier_file = directory / "provided_identifiers.txt"
                            identifier_file.write_text("\n".join(identifiers) + "\n", encoding="utf-8")
                            source = save_ncbi_protein_identifier_snapshot_from_file(
                                identifier_file_path=identifier_file,
                                identifier_type=identifier_type,
                                snapshot_root_directory=source_root,
                            )
                        source_payload = load_snapshot_by_directory(snapshot_directory=source)
                        records = [
                            f"<GBSeq><GBSeqid>gi|{identifier}</GBSeqid></GBSeq>"
                            if identifier_type == "uid" else
                            f"<GBSeq><GBSeq_accession-version>{identifier}</GBSeq_accession-version></GBSeq>"
                            for identifier in identifiers
                        ]
                        xml_payload = ("<GBSet>" + "".join(records) + "</GBSet>").encode("utf-8")
                        with patch(
                            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
                            return_value=FakeFetchHandle(xml_payload),
                        ) as efetch_mock, patch(
                            "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_uid_snapshot",
                            side_effect=AssertionError("A provided identifier snapshot must not trigger ESearch."),
                        ) as search_mock:
                            created = resolve_ncbi_protein_xml_snapshot(
                                snapshot_mode="reuse_latest_or_create",
                                snapshot_root_directory=directory / "xml_snapshot",
                                source_uid_snapshot_root_directory=source_root,
                                ncbi_email="test@example.org",
                                xml_request_delay_seconds=0.0,
                            )
                        reused = self._reuse(directory)
                    efetch_mock.assert_called_once_with(
                        database_name="protein", protein_uids=identifiers, retmode="xml",
                    )
                    search_mock.assert_not_called()
                    self.assertEqual(source_payload["manifest"]["artifact_type"], "ncbi_protein_identifier_snapshot")
                    self.assertEqual(source_payload["manifest"]["source_type"], "provided_identifier_list")
                    self.assertNotIn("search_query", source_payload["manifest"])
                    self.assertEqual(created["manifest"]["identifier_type"], identifier_type)
                    self.assertEqual(created["protein_uids"], identifiers)
                    self.assertEqual(reused["protein_uids"], identifiers)
                    self.assertEqual(reused["manifest"], created["manifest"])
        self.assertEqual(attempted_spreadsheet_paths, [])

    def test_xml_resolver_creates_typed_snapshots_with_either_source_root_name(self) -> None:
        for identifier_type, identifiers in (
            ("uid", ["1002", "1001", "1002"]),
            ("accession.version", ["XP_123456789.2", "WP_012345678.1", "XP_123456789.2"]),
        ):
            for mode in ("create_new", "reuse_latest_or_create"):
                for argument in ("legacy", "generic", "both"):
                    with self.subTest(identifier_type=identifier_type, mode=mode, argument=argument), tempfile.TemporaryDirectory() as name:
                        directory = Path(name)
                        source_root = self._frozen_resolver_source(directory, identifiers, identifier_type)
                        (source_root / "latest/protein_uids.txt").unlink()
                        options = {}
                        if argument in ("legacy", "both"):
                            options["source_uid_snapshot_root_directory"] = source_root
                        if argument in ("generic", "both"):
                            options["source_identifier_snapshot_root_directory"] = str(source_root)
                        result = self._fetch_result(identifiers, identifier_type)
                        with patch(
                            "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches", return_value=result,
                        ) as fetch_mock, patch(
                            "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_uid_snapshot",
                            side_effect=AssertionError("A frozen source must not trigger ESearch."),
                        ):
                            payload = resolve_ncbi_protein_xml_snapshot(
                                snapshot_mode=mode, snapshot_root_directory=directory / "xml_snapshot",
                                ncbi_email="test@example.org", **options,
                            )
                        fetch_mock.assert_called_once()
                        self.assertEqual(fetch_mock.call_args.kwargs["identifier_type"], identifier_type)
                        self.assertEqual(fetch_mock.call_args.kwargs["protein_uids"], identifiers)
                        self.assertEqual(payload["protein_uids"], identifiers)
                        self.assertEqual(payload["manifest"]["identifier_type"], identifier_type)
                        self.assertEqual(payload["manifest"]["source_type"], "identifier_snapshot")
                        self.assertEqual(payload["manifest"]["source_identifier_sha256"], result.protein_identifiers_sha256)

    def test_xml_resolver_reuses_without_source_roots_or_credentials(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            result = self._fetch_result(["WP_012345678.1"], "accession.version")
            self._save(directory, result, ["WP_012345678.1"], self._source_manifest(result))
            conflicting_sources = {
                "source_uid_snapshot_root_directory": directory / "unused_old",
                "source_identifier_snapshot_root_directory": directory / "unused_new",
            }
            for mode, options in (
                ("reuse_latest", {}),
                ("reuse_latest", conflicting_sources),
                ("reuse_latest_or_create", {}),
            ):
                with self.subTest(mode=mode, options=options), patch(
                    "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches",
                    side_effect=AssertionError("Reusing XML must not fetch."),
                ), patch(
                    "src.pago_pipeline.ncbi_snapshot.resolve_ncbi_protein_uid_snapshot",
                    side_effect=AssertionError("Reusing XML must not load a source."),
                ):
                    payload = resolve_ncbi_protein_xml_snapshot(
                        snapshot_mode=mode, snapshot_root_directory=directory / "xml_snapshot", **options,
                    )
                    self.assertEqual(payload["protein_uids"], ["WP_012345678.1"])
            with self.assertRaisesRegex(ValueError, "must refer to the same source directory"):
                resolve_ncbi_protein_xml_snapshot(
                    snapshot_mode="reuse_latest_or_create",
                    snapshot_root_directory=directory / "xml_snapshot",
                    **conflicting_sources,
                )

    def test_automatic_xml_reuse_validates_xml_only_for_a_matching_source(self) -> None:
        for identifier_type, original, replacement in (
            ("uid", "1001", "1002"),
            ("accession.version", "WP_012345678.1", "WP_012345678.2"),
        ):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                source_root = directory / "identifiers"
                xml_root = directory / "xml_snapshot"
                arguments = {
                    "snapshot_mode": "reuse_latest_or_create",
                    "snapshot_root_directory": xml_root,
                    "source_identifier_snapshot_root_directory": source_root,
                }
                save_ncbi_protein_identifier_snapshot(
                    protein_identifiers=[original], identifier_type=identifier_type,
                    snapshot_root_directory=source_root, retrieved_at_utc="2026-10-15T12:00:00Z",
                )
                first_result = self._fetch_result([original], identifier_type)
                second_result = replace(
                    self._fetch_result([replacement], identifier_type),
                    retrieved_at_utc="2026-10-16T13:00:00Z",
                )
                with patch(
                    "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches",
                    side_effect=[first_result, second_result],
                ) as fetch:
                    initial = resolve_ncbi_protein_xml_snapshot(**arguments, ncbi_email="test@example.com")
                    cached = resolve_ncbi_protein_xml_snapshot(**arguments)
                    self.assertEqual(cached["manifest"], initial["manifest"])
                    fetch.assert_called_once()

                    latest_xml = cached["xml_file_path"]
                    pointer = (
                        "version https://git-lfs.github.com/spec/v1\n"
                        f"oid sha256:{cached['manifest']['xml_file_sha256']}\n"
                        f"size {latest_xml.stat().st_size}\n"
                    ).encode("ascii")
                    latest_xml.write_bytes(pointer)
                    with self.assertRaisesRegex(RuntimeError, "not well-formed"):
                        resolve_ncbi_protein_xml_snapshot(**arguments, ncbi_email="test@example.com")
                    fetch.assert_called_once()

                    updated_source = save_ncbi_protein_identifier_snapshot(
                        protein_identifiers=[replacement], identifier_type=identifier_type,
                        snapshot_root_directory=source_root, retrieved_at_utc="2026-10-15T13:00:00Z",
                    )
                    with self.assertRaisesRegex(ValueError, "ncbi_email is required"):
                        resolve_ncbi_protein_xml_snapshot(**arguments)
                    fetch.assert_called_once()
                    self.assertEqual(latest_xml.read_bytes(), pointer)

                    refreshed = resolve_ncbi_protein_xml_snapshot(**arguments, ncbi_email="test@example.com")
                    self.assertEqual(fetch.call_count, 2)
                    self.assertEqual(fetch.call_args.kwargs["protein_uids"], [replacement])
                    self.assertEqual(refreshed["protein_uids"], [replacement])
                    self.assertEqual(
                        load_latest_xml_snapshot(snapshot_root_directory=xml_root)["protein_uids"], [replacement],
                    )
                    self.assertEqual(
                        refreshed["manifest"]["source_identifier_snapshot_manifest_sha256"],
                        hashlib.sha256((updated_source / "manifest.json").read_bytes()).hexdigest(),
                    )

    def test_automatic_xml_reuse_accepts_line_ending_conversion_of_legacy_source(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            source_root = directory / "identifiers"
            source_latest = source_root / "latest"
            source_latest.mkdir(parents=True)
            result = self._fetch_result(["1001"], "uid")
            source = {
                "protein_uids_sha256": result.protein_identifiers_sha256,
                "normalized_protein_uid_count": 1,
            }
            (source_latest / "protein_uids.txt").write_bytes(b"1001\n")
            manifest_path = source_latest / "manifest.json"
            manifest_path.write_bytes((json.dumps(source, indent=2) + "\n").encode("utf-8"))
            arguments = {
                "snapshot_mode": "reuse_latest_or_create",
                "snapshot_root_directory": directory / "xml_snapshot",
                "source_uid_snapshot_root_directory": source_root,
            }
            with patch(
                "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches",
                return_value=result,
            ) as fetch:
                initial = resolve_ncbi_protein_xml_snapshot(**arguments, ncbi_email="test@example.com")
                manifest_path.write_bytes(manifest_path.read_bytes().replace(b"\n", b"\r\n"))
                cached = resolve_ncbi_protein_xml_snapshot(**arguments)
                self.assertEqual(cached["manifest"], initial["manifest"])
                source["source_description"] = "Different source metadata"
                manifest_path.write_text(json.dumps(source), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "ncbi_email is required"):
                    resolve_ncbi_protein_xml_snapshot(**arguments)
                fetch.assert_called_once()

    def test_xml_resolver_rejects_missing_or_conflicting_creation_source_arguments(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            for options, message in (
                ({}, "is required when creating an XML snapshot"),
                ({"source_uid_snapshot_root_directory": directory / "old",
                  "source_identifier_snapshot_root_directory": directory / "new"},
                 "must refer to the same source directory"),
            ):
                with self.subTest(options=options), patch(
                    "src.pago_pipeline.ncbi_snapshot.resolve_ncbi_protein_uid_snapshot",
                ) as source_mock, patch(
                    "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches",
                ) as fetch_mock, patch(
                    "src.pago_pipeline.ncbi_snapshot.save_ncbi_protein_xml_snapshot",
                ) as save_mock:
                    with self.assertRaisesRegex(ValueError, message):
                        resolve_ncbi_protein_xml_snapshot(
                            snapshot_mode="create_new", snapshot_root_directory=directory / "xml_snapshot",
                            ncbi_email="test@example.org", **options,
                        )
                    source_mock.assert_not_called()
                    fetch_mock.assert_not_called()
                    save_mock.assert_not_called()
                    self.assertFalse((directory / "xml_snapshot").exists())

    def test_xml_resolver_requires_a_frozen_source_before_fetch(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            with patch("src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches") as fetch_mock:
                with self.assertRaisesRegex(FileNotFoundError, "source identifier snapshot"):
                    resolve_ncbi_protein_xml_snapshot(
                        snapshot_mode="create_new", snapshot_root_directory=directory / "xml_snapshot",
                        source_identifier_snapshot_root_directory=directory / "missing_source",
                        ncbi_email="test@example.org",
                    )
                fetch_mock.assert_not_called()
            self.assertFalse((directory / "missing_source").exists())
            self.assertFalse((directory / "xml_snapshot").exists())

    def test_xml_resolver_rejects_fetch_type_mismatch_before_calling_save(self) -> None:
        for source_type, fetch_type, identifiers in (
            ("uid", "accession.version", ["1001"]),
            ("accession.version", "uid", ["WP_012345678.1"]),
            ("uid", None, ["1001"]),
            ("uid", "unknown", ["1001"]),
        ):
            with self.subTest(source_type=source_type, fetch_type=fetch_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                source_root = self._frozen_resolver_source(directory, identifiers, source_type)
                result = replace(self._fetch_result(identifiers, source_type), identifier_type=fetch_type)
                with patch(
                    "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches", return_value=result,
                ) as fetch_mock, patch(
                    "src.pago_pipeline.ncbi_snapshot.save_ncbi_protein_xml_snapshot",
                ) as save_mock:
                    with self.assertRaisesRegex(ValueError, "identifier_type"):
                        resolve_ncbi_protein_xml_snapshot(
                            snapshot_mode="create_new", snapshot_root_directory=directory / "xml_snapshot",
                            source_identifier_snapshot_root_directory=source_root, ncbi_email="test@example.org",
                        )
                    self.assertEqual(fetch_mock.call_args.kwargs["identifier_type"], source_type)
                    save_mock.assert_not_called()
                self.assertFalse((directory / "xml_snapshot").exists())

    def test_xml_resolver_rejects_invalid_source_type_before_fetch(self) -> None:
        for identifier_type in ("absent", None, "", "unknown"):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                source_root = self._frozen_resolver_source(directory, ["1001"], "uid")
                manifest_file = source_root / "latest/manifest.json"
                manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
                if identifier_type == "absent":
                    manifest.pop("identifier_type")
                else:
                    manifest["identifier_type"] = identifier_type
                manifest_file.write_text(json.dumps(manifest), encoding="utf-8")
                with patch("src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches") as fetch_mock:
                    with self.assertRaisesRegex(ValueError, "identifier_type"):
                        resolve_ncbi_protein_xml_snapshot(
                            snapshot_mode="create_new", snapshot_root_directory=directory / "xml_snapshot",
                            source_identifier_snapshot_root_directory=source_root, ncbi_email="test@example.org",
                        )
                    fetch_mock.assert_not_called()

    def test_xml_resolver_creates_explicitly_typed_xml_from_legacy_uid_source(self) -> None:
        for version in (None, "1.0"):
            with self.subTest(version=version), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                source_root = directory / "identifiers"
                latest = source_root / "latest"
                latest.mkdir(parents=True)
                identifiers = ["1001", "1002"]
                result = self._fetch_result(identifiers, "uid")
                (latest / "protein_uids.txt").write_text("1001\n1002\n", encoding="utf-8")
                manifest = {"protein_uids_sha256": result.protein_identifiers_sha256,
                            "normalized_protein_uid_count": 2}
                if version is not None:
                    manifest["snapshot_format_version"] = version
                manifest_file = latest / "manifest.json"
                manifest_file.write_text(json.dumps(manifest), encoding="utf-8")
                original = manifest_file.read_bytes()
                with patch(
                    "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches", return_value=result,
                ) as fetch_mock:
                    payload = resolve_ncbi_protein_xml_snapshot(
                        snapshot_mode="create_new", snapshot_root_directory=directory / "xml_snapshot",
                        source_uid_snapshot_root_directory=source_root, ncbi_email="test@example.org",
                    )
                self.assertEqual(fetch_mock.call_args.kwargs["identifier_type"], "uid")
                self.assertEqual(payload["manifest"]["identifier_type"], "uid")
                self.assertEqual(payload["manifest"]["snapshot_format_version"], "1.1")
                self.assertEqual(payload["manifest"]["source_identifier_sha256"], result.protein_identifiers_sha256)
                self.assertEqual(
                    payload["manifest"]["source_identifier_snapshot_manifest_sha256"],
                    hashlib.sha256(original).hexdigest(),
                )
                self.assertEqual(manifest_file.read_bytes(), original)


class NCBIXmlFetchInfrastructureTests(unittest.TestCase):
    def setUp(self) -> None:
        with ncbi_api._default_ncbi_xml_circuit_breakers_lock:
            ncbi_api._default_ncbi_xml_circuit_breakers.clear()

    def test_identifier_type_coercion_accepts_supported_values(self) -> None:
        self.assertEqual(
            coerce_ncbi_protein_identifier_type(None),
            NCBIProteinIdentifierType.uid,
        )
        self.assertEqual(
            coerce_ncbi_protein_identifier_type("uid"),
            NCBIProteinIdentifierType.uid,
        )
        self.assertEqual(
            coerce_ncbi_protein_identifier_type("accession.version"),
            NCBIProteinIdentifierType.accession_version,
        )

    def test_identifier_type_coercion_rejects_unknown_value(self) -> None:
        with self.assertRaisesRegex(
            ValueError,
            "Invalid NCBI protein identifier_type",
        ):
            coerce_ncbi_protein_identifier_type("accession")

    def test_identifier_validation_rejects_mixed_uid_and_accession_list(self) -> None:
        with self.assertRaisesRegex(
            ValueError,
            "Mixed or invalid identifier lists are not supported",
        ):
            validate_ncbi_protein_identifiers(
                protein_identifiers=["1001", "WP_012345678.1"],
                identifier_type="uid",
            )

    def test_generic_identifier_normalizer_accepts_valid_uids_with_whitespace(self) -> None:
        normalized_identifiers = _normalize_protein_identifier_list(
            protein_identifiers=["1001", " 1002 "],
            identifier_type="uid",
            deduplicate_identifiers=False,
            sort_identifiers=False,
        )
        self.assertEqual(normalized_identifiers, ["1001", "1002"])

    def test_generic_identifier_normalizer_rejects_accession_in_uid_mode(self) -> None:
        with self.assertRaisesRegex(ValueError, "must contain only uid identifiers"):
            _normalize_protein_identifier_list(
                protein_identifiers=["WP_012345678.1"],
                identifier_type="uid",
                deduplicate_identifiers=False,
                sort_identifiers=False,
            )

    def test_generic_identifier_normalizer_rejects_unknown_identifier_type(self) -> None:
        for identifiers in (["1001"], ["WP_012345678.1"], []):
            with self.subTest(identifiers=identifiers):
                with self.assertRaisesRegex(ValueError, "Invalid NCBI protein identifier_type"):
                    _normalize_protein_identifier_list(
                        protein_identifiers=identifiers,
                        identifier_type="unknown",
                        deduplicate_identifiers=False,
                        sort_identifiers=False,
                    )

    def test_generic_identifier_normalizer_strips_deduplicates_sorts_and_validates_uid(
        self,
    ) -> None:
        normalized_identifiers = _normalize_protein_identifier_list(
            protein_identifiers=[" 1002 ", "", 1001, "1002"],
            identifier_type="uid",
            deduplicate_identifiers=True,
            sort_identifiers=True,
        )

        self.assertEqual(normalized_identifiers, ["1001", "1002"])

    def test_legacy_uid_normalizer_delegates_to_uid_identifier_contract(self) -> None:
        normalized_uids = _normalize_protein_uid_list(
            protein_uids=[" 1002 ", "", "1001", "1002"],
            deduplicate_uids=True,
            sort_uids=True,
        )

        self.assertEqual(normalized_uids, ["1001", "1002"])

        with self.assertRaisesRegex(
            ValueError,
            "Mixed or invalid identifier lists are not supported",
        ):
            _normalize_protein_uid_list(
                protein_uids=["1001", "WP_012345678.1"],
                deduplicate_uids=True,
                sort_uids=True,
            )

    def test_generic_identifier_normalizer_accepts_accession_version(self) -> None:
        for identifiers in (
            ["WP_012345678.1", " XP_123456789.2 "],
            [" WP_012345678.1 ", "", "XP_123456789.2"],
        ):
            with self.subTest(identifiers=identifiers):
                normalized_identifiers = _normalize_protein_identifier_list(
                    protein_identifiers=identifiers,
                    identifier_type="accession.version",
                    deduplicate_identifiers=False,
                    sort_identifiers=False,
                )
                self.assertEqual(
                    normalized_identifiers,
                    ["WP_012345678.1", "XP_123456789.2"],
                )

    def test_generic_identifier_normalizer_rejects_invalid_accession_versions(
        self,
    ) -> None:
        invalid_accession_versions = [
            "WP_012345678",
            "WP_012345678.A",
            "WP_012345678.0",
            "123456789",
            "123456789.1",
        ]

        for protein_identifier in invalid_accession_versions:
            with self.subTest(protein_identifier=protein_identifier):
                with self.assertRaisesRegex(
                    ValueError,
                    "Mixed or invalid identifier lists are not supported",
                ):
                    _normalize_protein_identifier_list(
                        protein_identifiers=[protein_identifier],
                        identifier_type="accession.version",
                        deduplicate_identifiers=False,
                        sort_identifiers=False,
                    )

    def test_xml_fetch_retries_transient_failure_and_returns_batch(self) -> None:
        fetch_error = URLError("timed out")
        successful_payload = b"<GBSet><GBSeq /><GBSeq /></GBSet>"
        successful_handle = FakeFetchHandle(successful_payload)

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[fetch_error, successful_handle],
        ) as efetch_mock:
            with patch(
                "src.pago_pipeline.ncbi_api.random.uniform",
                return_value=0.0,
            ):
                with patch("src.pago_pipeline.ncbi_api.time.sleep") as sleep_mock:
                    fetch_result = fetch_ncbi_protein_xml_batches(
                        ncbi_email="test@example.org",
                        ncbi_api_key=None,
                        protein_uids=["1001", "1002"],
                        batch_size=2,
                        max_retry_attempts=2,
                        request_delay_seconds=0.01,
                        fetch_timeout_seconds=5.0,
                        batch_deadline_seconds=30.0,
                        retry_backoff_initial_seconds=0.25,
                    )

        self.assertEqual(efetch_mock.call_count, 2)
        self.assertEqual(fetch_result.batch_count, 1)
        self.assertEqual(fetch_result.identifier_type, "uid")
        self.assertEqual(
            fetch_result.xml_batches[0].xml_payload_bytes,
            successful_payload,
        )
        sleep_mock.assert_any_call(0.25)
        sleep_mock.assert_any_call(0.01)

    def test_xml_fetch_accepts_explicit_accession_version_identifiers(self) -> None:
        successful_payload = b"<GBSet><GBSeq /></GBSet>"

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            return_value=FakeFetchHandle(successful_payload),
        ) as efetch_mock:
            with patch(
                "src.pago_pipeline.ncbi_api.random.uniform",
                return_value=0.0,
            ):
                with patch("src.pago_pipeline.ncbi_api.time.sleep"):
                    fetch_result = fetch_ncbi_protein_xml_batches(
                        ncbi_email="test@example.org",
                        ncbi_api_key=None,
                        protein_uids=["WP_012345678.1"],
                        identifier_type="accession.version",
                        batch_size=1,
                        request_delay_seconds=0.0,
                    )

        self.assertEqual(fetch_result.identifier_type, "accession.version")
        self.assertEqual(
            efetch_mock.call_args.kwargs["protein_uids"],
            ["WP_012345678.1"],
        )

    def test_xml_fetch_rejects_unknown_identifier_type_before_network(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
        ) as efetch_mock:
            with self.assertRaisesRegex(
                ValueError,
                "Invalid NCBI protein identifier_type",
            ):
                fetch_ncbi_protein_xml_batches(
                    ncbi_email="test@example.org",
                    ncbi_api_key=None,
                    protein_uids=["1001"],
                    identifier_type="protein_accession",
                )

        efetch_mock.assert_not_called()

    def test_xml_fetch_normalizes_identifiers_preserving_order_and_duplicates(self) -> None:
        cases = [
            ({}, [1002, " 1001 ", "", "1002", " \t"], ["1002", "1001", "1002"]),
            (
                {"identifier_type": "accession.version"},
                [" XP_123456789.2 ", "WP_012345678.1", "", "XP_123456789.2", " \t"],
                ["XP_123456789.2", "WP_012345678.1", "XP_123456789.2"],
            ),
        ]
        for options, raw_identifiers, expected_identifiers in cases:
            with self.subTest(options=options):
                with patch(
                    "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
                    side_effect=[
                        FakeFetchHandle(b"<GBSet><GBSeq /><GBSeq /></GBSet>"),
                        FakeFetchHandle(b"<GBSet><GBSeq /></GBSet>"),
                    ],
                ) as efetch_mock, patch("src.pago_pipeline.ncbi_api.time.sleep"):
                    fetch_result = fetch_ncbi_protein_xml_batches(
                        ncbi_email="test@example.org",
                        ncbi_api_key=None,
                        protein_uids=raw_identifiers,
                        batch_size=2,
                        request_delay_seconds=0.0,
                        **options,
                    )

                self.assertEqual(
                    fetch_result.identifier_type,
                    options.get("identifier_type", "uid"),
                )
                self.assertEqual(fetch_result.requested_protein_uid_count, 5)
                self.assertEqual(fetch_result.normalized_protein_uid_count, 3)
                self.assertEqual(fetch_result.requested_protein_identifier_count, 5)
                self.assertEqual(fetch_result.normalized_protein_identifier_count, 3)
                self.assertEqual(fetch_result.protein_identifier_count, 3)
                self.assertEqual(
                    fetch_result.protein_identifiers_sha256,
                    fetch_result.protein_uids_sha256,
                )
                serialized = asdict(fetch_result)
                self.assertEqual(serialized["requested_protein_uid_count"], 5)
                self.assertNotIn("protein_identifiers_sha256", serialized)
                for batch in fetch_result.xml_batches:
                    self.assertIs(batch.protein_identifiers, batch.protein_uids)
                    self.assertEqual(batch.protein_identifier_count, batch.protein_uid_count)
                    self.assertEqual(NCBIProteinXmlBatchFetchResult(**asdict(batch)), batch)
                    self.assertNotIn("protein_identifiers", asdict(batch))
                self.assertEqual(fetch_result.batch_count, 2)
                self.assertEqual(
                    fetch_result.protein_uids_sha256,
                    sha256_of_lines(
                        text_lines=expected_identifiers,
                        deduplicate_lines_preserving_order=False,
                        sort_lines=False,
                    ),
                )
                self.assertEqual(
                    [request.kwargs["protein_uids"] for request in efetch_mock.call_args_list],
                    [expected_identifiers[:2], expected_identifiers[2:]],
                )
                self.assertEqual(
                    [batch.protein_uids for batch in fetch_result.xml_batches],
                    [expected_identifiers[:2], expected_identifiers[2:]],
                )
                self.assertEqual(
                    [
                        (batch.batch_start_index, batch.batch_end_index, batch.protein_uid_count)
                        for batch in fetch_result.xml_batches
                    ],
                    [(0, 1, 2), (2, 2, 1)],
                )

    def test_xml_fetch_sends_both_identifier_domains_in_efetch_id_parameter(self) -> None:
        for identifier_type in ("uid", "accession.version"):
            for count, method in ((2, "GET"), (200, "POST")):
                with self.subTest(identifier_type=identifier_type, method=method):
                    identifiers = [
                        str(1000 + index) if identifier_type == "uid"
                        else f"WP_{index:09d}.1"
                        for index in range(count)
                    ]
                    payload = b"<GBSet>" + b"<GBSeq />" * count + b"</GBSet>"
                    with patch.object(
                        ncbi_api.Entrez, "urlopen", return_value=FakeFetchHandle(payload),
                    ) as urlopen_mock, patch("src.pago_pipeline.ncbi_api.time.sleep"):
                        fetch_result = fetch_ncbi_protein_xml_batches(
                            ncbi_email="test@example.org",
                            ncbi_api_key=None,
                            protein_uids=identifiers,
                            identifier_type=identifier_type,
                            batch_size=count,
                            request_delay_seconds=0.0,
                        )

                    urlopen_mock.assert_called_once()
                    request = urlopen_mock.call_args.args[0]
                    self.assertEqual(request.get_method(), method)
                    parameters = parse_qs(
                        request.data.decode("utf-8") if request.data is not None
                        else urlsplit(request.full_url).query
                    )
                    self.assertEqual(parameters["id"], [",".join(identifiers)])
                    self.assertEqual(parameters["db"], ["protein"])
                    self.assertEqual(parameters["retmode"], ["xml"])
                    self.assertEqual(fetch_result.identifier_type, identifier_type)

    def test_xml_fetch_rejects_mixed_identifier_list_before_network(self) -> None:
        cases = [{}, {"identifier_type": "uid"}, {"identifier_type": "accession.version"}]
        for options in cases:
            with self.subTest(options=options):
                with patch(
                    "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
                ) as efetch_mock:
                    with self.assertRaisesRegex(
                        ValueError,
                        "Mixed or invalid identifier lists are not supported",
                    ):
                        fetch_ncbi_protein_xml_batches(
                            ncbi_email="test@example.org",
                            ncbi_api_key=None,
                            protein_uids=["1001", "WP_012345678.1"],
                            batch_size=1,
                            **options,
                        )

                efetch_mock.assert_not_called()

    def test_xml_fetch_rejects_invalid_or_empty_identifiers_before_network(self) -> None:
        cases = [
            ({}, ["WP_012345678.1"]),
            ({"identifier_type": "accession.version"}, ["1001"]),
            *[
                ({"identifier_type": "accession.version"}, [identifier])
                for identifier in ("WP_012345678", "WP_012345678.A", "WP_012345678.0")
            ],
            *[
                ({"identifier_type": identifier_type}, identifiers)
                for identifier_type in ("uid", "accession.version")
                for identifiers in ([], ["", " \t"])
            ],
        ]
        for options, identifiers in cases:
            with self.subTest(options=options, identifiers=identifiers):
                with patch(
                    "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
                ) as efetch_mock:
                    with self.assertRaises(ValueError):
                        fetch_ncbi_protein_xml_batches(
                            ncbi_email="test@example.org",
                            ncbi_api_key=None,
                            protein_uids=identifiers,
                            **options,
                        )
                efetch_mock.assert_not_called()

    def test_xml_fetch_retries_socket_timeout_raised_during_response_read(
        self,
    ) -> None:
        timed_out_handle = FakeFetchHandle(
            b"",
            read_callback=lambda: (_ for _ in ()).throw(
                socket.timeout("response read timed out")
            ),
        )
        successful_payload = b"<GBSet><GBSeq /></GBSet>"

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[
                timed_out_handle,
                FakeFetchHandle(successful_payload),
            ],
        ) as efetch_mock:
            with patch(
                "src.pago_pipeline.ncbi_api.random.uniform",
                return_value=0.0,
            ):
                with patch("src.pago_pipeline.ncbi_api.time.sleep") as sleep_mock:
                    fetch_result = fetch_ncbi_protein_xml_batches(
                        ncbi_email="test@example.org",
                        ncbi_api_key=None,
                        protein_uids=["1001"],
                        batch_size=1,
                        max_retry_attempts=2,
                        request_delay_seconds=0.0,
                        fetch_timeout_seconds=0.25,
                        batch_deadline_seconds=5.0,
                        retry_backoff_initial_seconds=0.1,
                    )

        self.assertTrue(timed_out_handle.closed)
        self.assertEqual(efetch_mock.call_count, 2)
        self.assertEqual(
            fetch_result.xml_batches[0].xml_payload_bytes,
            successful_payload,
        )
        sleep_mock.assert_any_call(0.1)

    def test_xml_fetch_applies_exponential_backoff_with_configured_cap(
        self,
    ) -> None:
        successful_payload = b"<GBSet><GBSeq /></GBSet>"

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[
                URLError("temporary failure 1"),
                URLError("temporary failure 2"),
                URLError("temporary failure 3"),
                FakeFetchHandle(successful_payload),
            ],
        ) as efetch_mock:
            with patch(
                "src.pago_pipeline.ncbi_api.random.uniform",
                return_value=0.0,
            ):
                with patch("src.pago_pipeline.ncbi_api.time.sleep") as sleep_mock:
                    fetch_ncbi_protein_xml_batches(
                        ncbi_email="test@example.org",
                        ncbi_api_key=None,
                        protein_uids=["1001"],
                        batch_size=1,
                        max_retry_attempts=4,
                        request_delay_seconds=0.0,
                        fetch_timeout_seconds=1.0,
                        batch_deadline_seconds=10.0,
                        retry_backoff_initial_seconds=0.25,
                        retry_backoff_multiplier=2.0,
                        retry_backoff_max_seconds=0.6,
                    )

        self.assertEqual(efetch_mock.call_count, 4)
        self.assertEqual(
            sleep_mock.call_args_list[:3],
            [call(0.25), call(0.5), call(0.6)],
        )

    def test_xml_fetch_retries_transient_http_429_and_503(self) -> None:
        http_429 = HTTPError(
            "https://ncbi.example.test/efetch",
            429,
            "Too Many Requests",
            None,
            None,
        )
        http_503 = HTTPError(
            "https://ncbi.example.test/efetch",
            503,
            "Service Unavailable",
            None,
            None,
        )
        self.addCleanup(http_429.close)
        self.addCleanup(http_503.close)

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[
                http_429,
                http_503,
                FakeFetchHandle(b"<GBSet><GBSeq /></GBSet>"),
            ],
        ) as efetch_mock:
            with patch(
                "src.pago_pipeline.ncbi_api.random.uniform",
                return_value=0.0,
            ):
                with patch("src.pago_pipeline.ncbi_api.time.sleep"):
                    fetch_ncbi_protein_xml_batches(
                        ncbi_email="test@example.org",
                        ncbi_api_key=None,
                        protein_uids=["1001"],
                        max_retry_attempts=3,
                        request_delay_seconds=0.0,
                        retry_backoff_initial_seconds=0.1,
                    )

        self.assertEqual(efetch_mock.call_count, 3)

    def test_xml_fetch_does_not_retry_permanent_http_error(self) -> None:
        http_400 = HTTPError(
            "https://ncbi.example.test/efetch",
            400,
            "Bad Request",
            None,
            None,
        )
        self.addCleanup(http_400.close)

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=http_400,
        ) as efetch_mock:
            with patch("src.pago_pipeline.ncbi_api.time.sleep") as sleep_mock:
                with self.assertRaisesRegex(
                    RuntimeError,
                    "Permanent NCBI XML batch failure; not retrying",
                ):
                    fetch_ncbi_protein_xml_batches(
                        ncbi_email="test@example.org",
                        ncbi_api_key=None,
                        protein_uids=["1001"],
                        max_retry_attempts=3,
                        request_delay_seconds=0.0,
                    )

        efetch_mock.assert_called_once()
        sleep_mock.assert_not_called()

    def test_xml_fetch_passes_network_timeout_to_entrez_urlopen(self) -> None:
        captured_urlopen_timeout: dict[str, float | None] = {"timeout": None}

        def fake_urlopen(*args, **kwargs):
            captured_urlopen_timeout["timeout"] = kwargs.get("timeout")
            return FakeFetchHandle(b"<GBSet><GBSeq /></GBSet>")

        with patch.object(ncbi_api.Entrez, "urlopen", fake_urlopen, create=True):
            with patch(
                "src.pago_pipeline.ncbi_api.random.uniform",
                return_value=0.0,
            ):
                with patch("src.pago_pipeline.ncbi_api.time.sleep"):
                    fetch_ncbi_protein_xml_batches(
                        ncbi_email="test@example.org",
                        ncbi_api_key=None,
                        protein_uids=["1001"],
                        batch_size=1,
                        request_delay_seconds=0.0,
                        fetch_timeout_seconds=12.5,
                        batch_deadline_seconds=30.0,
                    )

        self.assertEqual(captured_urlopen_timeout["timeout"], 12.5)

    def test_xml_fetch_does_not_mutate_global_biopython_retry_policy(
        self,
    ) -> None:
        captured_entrez_retry_policy: dict[str, float | int | None] = {}

        def fake_urlopen(*args, **kwargs):
            captured_entrez_retry_policy["max_tries"] = getattr(
                ncbi_api.Entrez,
                "max_tries",
                None,
            )
            captured_entrez_retry_policy["sleep_between_tries"] = getattr(
                ncbi_api.Entrez,
                "sleep_between_tries",
                None,
            )
            return FakeFetchHandle(b"<GBSet><GBSeq /></GBSet>")

        with patch.object(ncbi_api.Entrez, "max_tries", 3, create=True):
            with patch.object(
                ncbi_api.Entrez,
                "sleep_between_tries",
                15,
                create=True,
            ):
                with patch.object(
                    ncbi_api.Entrez,
                    "urlopen",
                    fake_urlopen,
                    create=True,
                ):
                    with patch(
                        "src.pago_pipeline.ncbi_api.random.uniform",
                        return_value=0.0,
                    ):
                        with patch("src.pago_pipeline.ncbi_api.time.sleep"):
                            fetch_ncbi_protein_xml_batches(
                                ncbi_email="test@example.org",
                                ncbi_api_key=None,
                                protein_uids=["1001"],
                                batch_size=1,
                                request_delay_seconds=0.0,
                                fetch_timeout_seconds=12.5,
                                batch_deadline_seconds=30.0,
                            )

                self.assertEqual(ncbi_api.Entrez.max_tries, 3)
                self.assertEqual(ncbi_api.Entrez.sleep_between_tries, 15)

        self.assertEqual(captured_entrez_retry_policy["max_tries"], 3)
        self.assertEqual(captured_entrez_retry_policy["sleep_between_tries"], 15)

    def test_xml_fetch_raises_explicit_error_after_retry_limit(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=URLError("temporarily unavailable"),
        ):
            with patch("src.pago_pipeline.ncbi_api.time.sleep") as sleep_mock:
                with self.assertRaisesRegex(
                    NCBIXmlBatchFetchError,
                    "Failed to fetch NCBI XML batch after 2 attempts",
                ):
                    fetch_ncbi_protein_xml_batches(
                        ncbi_email="test@example.org",
                        ncbi_api_key=None,
                        protein_uids=["1001", "1002"],
                        batch_size=2,
                        max_retry_attempts=2,
                        request_delay_seconds=0.01,
                        fetch_timeout_seconds=5.0,
                        batch_deadline_seconds=30.0,
                        retry_backoff_initial_seconds=0.5,
                    )

        sleep_mock.assert_called_once_with(0.5)

    def test_xml_fetch_deadline_aborts_slow_read(self) -> None:
        read_started = threading.Event()
        release_blocked_read = threading.Event()
        blocked_handle_closed = threading.Event()

        def block_response_read() -> None:
            read_started.set()
            release_blocked_read.wait(timeout=2.0)

        blocked_handle = FakeFetchHandle(
            b"<GBSet><GBSeq /></GBSet>",
            read_callback=block_response_read,
        )

        def close_blocked_handle() -> None:
            blocked_handle.closed = True
            blocked_handle_closed.set()

        blocked_handle.close = close_blocked_handle
        deadline_test_started_at = time.perf_counter()
        try:
            with patch(
                "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
                return_value=blocked_handle,
            ):
                with self.assertRaisesRegex(
                    NCBIXmlBatchDeadlineExceeded,
                    "deadline exhausted",
                ):
                    fetch_ncbi_protein_xml_batches(
                        ncbi_email="test@example.org",
                        ncbi_api_key=None,
                        protein_uids=["1001"],
                        batch_size=1,
                        max_retry_attempts=1,
                        request_delay_seconds=0.0,
                        fetch_timeout_seconds=1.0,
                        batch_deadline_seconds=0.05,
                    )
        finally:
            release_blocked_read.set()

        deadline_test_elapsed_seconds = (
            time.perf_counter() - deadline_test_started_at
        )
        self.assertTrue(read_started.is_set())
        self.assertTrue(blocked_handle_closed.wait(timeout=0.5))
        self.assertLess(deadline_test_elapsed_seconds, 1.0)

    def test_xml_fetch_raises_deadline_when_next_backoff_will_not_fit(
        self,
    ) -> None:
        fake_time = {"value": 0.0}
        circuit_breaker = NCBIXmlCircuitBreaker(
            failure_threshold=3,
            cooldown_seconds=60.0,
        )

        def fail_after_consuming_deadline(*args, **kwargs):
            fake_time["value"] = 0.75
            raise URLError("temporarily unavailable")

        with patch(
            "src.pago_pipeline.ncbi_api.time.monotonic",
            side_effect=lambda: fake_time["value"],
        ):
            with patch(
                "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
                side_effect=fail_after_consuming_deadline,
            ) as efetch_mock:
                with patch(
                    "src.pago_pipeline.ncbi_api.time.sleep",
                ) as sleep_mock:
                    with self.assertRaisesRegex(
                        NCBIXmlBatchDeadlineExceeded,
                        "deadline exhausted",
                    ):
                        fetch_ncbi_protein_xml_batches(
                            ncbi_email="test@example.org",
                            ncbi_api_key=None,
                            protein_uids=["1001"],
                            max_retry_attempts=2,
                            request_delay_seconds=0.0,
                            fetch_timeout_seconds=1.0,
                            batch_deadline_seconds=1.0,
                            retry_backoff_initial_seconds=0.5,
                            circuit_breaker=circuit_breaker,
                        )

        efetch_mock.assert_called_once()
        sleep_mock.assert_not_called()
        self.assertEqual(circuit_breaker.consecutive_failure_count, 1)
        self.assertIsNone(circuit_breaker.opened_at_seconds)

    def test_default_circuit_breaker_accumulates_failures_and_recovers(
        self,
    ) -> None:
        clock = {"value": 0.0}
        successful_payload = b"<GBSet><GBSeq /></GBSet>"

        with patch(
            "src.pago_pipeline.ncbi_api.time.monotonic",
            side_effect=lambda: clock["value"],
        ):
            with patch(
                "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
                side_effect=[
                    URLError("first degraded batch"),
                    URLError("second degraded batch"),
                    FakeFetchHandle(successful_payload),
                    URLError("failure after recovery"),
                ],
            ) as efetch_mock:
                with patch(
                    "src.pago_pipeline.ncbi_api.random.uniform",
                    return_value=0.0,
                ):
                    with patch("src.pago_pipeline.ncbi_api.time.sleep"):
                        for expected_failure_count in (1, 2):
                            with self.assertRaises(NCBIXmlBatchFetchError):
                                fetch_ncbi_protein_xml_batches(
                                    ncbi_email="test@example.org",
                                    ncbi_api_key=None,
                                    protein_uids=["1001"],
                                    max_retry_attempts=1,
                                    request_delay_seconds=0.0,
                                    circuit_breaker_failure_threshold=2,
                                    circuit_breaker_cooldown_seconds=30.0,
                                )
                            circuit_breaker = (
                                get_default_ncbi_xml_circuit_breaker(
                                    failure_threshold=2,
                                    cooldown_seconds=30.0,
                                )
                            )
                            self.assertEqual(
                                circuit_breaker.consecutive_failure_count,
                                expected_failure_count,
                            )

                        clock["value"] = 1.0
                        request_count_before_open_call = efetch_mock.call_count
                        with self.assertRaisesRegex(
                            NCBICircuitBreakerOpen,
                            "circuit breaker is open",
                        ):
                            fetch_ncbi_protein_xml_batches(
                                ncbi_email="test@example.org",
                                ncbi_api_key=None,
                                protein_uids=["1001"],
                                max_retry_attempts=1,
                                request_delay_seconds=0.0,
                                circuit_breaker_failure_threshold=2,
                                circuit_breaker_cooldown_seconds=30.0,
                            )
                        self.assertEqual(
                            efetch_mock.call_count,
                            request_count_before_open_call,
                        )

                        clock["value"] = 31.0
                        fetch_ncbi_protein_xml_batches(
                            ncbi_email="test@example.org",
                            ncbi_api_key=None,
                            protein_uids=["1001"],
                            max_retry_attempts=1,
                            request_delay_seconds=0.0,
                            circuit_breaker_failure_threshold=2,
                            circuit_breaker_cooldown_seconds=30.0,
                        )
                        self.assertEqual(
                            circuit_breaker.consecutive_failure_count,
                            0,
                        )
                        self.assertIsNone(circuit_breaker.opened_at_seconds)

                        clock["value"] = 32.0
                        with self.assertRaises(NCBIXmlBatchFetchError):
                            fetch_ncbi_protein_xml_batches(
                                ncbi_email="test@example.org",
                                ncbi_api_key=None,
                                protein_uids=["1001"],
                                max_retry_attempts=1,
                                request_delay_seconds=0.0,
                                circuit_breaker_failure_threshold=2,
                                circuit_breaker_cooldown_seconds=30.0,
                            )

        self.assertEqual(efetch_mock.call_count, 4)
        self.assertEqual(circuit_breaker.consecutive_failure_count, 1)
        self.assertIsNone(circuit_breaker.opened_at_seconds)

    def test_xml_fetch_rejects_partial_response_without_retrying(self) -> None:
        partial_payload = b"<GBSet><GBSeq><GBSeqid>1001</GBSeqid></GBSeq></GBSet>"

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            return_value=FakeFetchHandle(partial_payload),
        ) as efetch_mock:
            with patch("src.pago_pipeline.ncbi_api.time.sleep") as sleep_mock:
                with self.assertRaisesRegex(
                    RuntimeError,
                    "partial or mismatched XML batch",
                ) as raised_error:
                    fetch_ncbi_protein_xml_batches(
                        ncbi_email="test@example.org",
                        ncbi_api_key=None,
                        protein_uids=["1001", "1002"],
                        batch_size=2,
                        max_retry_attempts=3,
                        request_delay_seconds=0.0,
                    )

        self.assertIsInstance(
            raised_error.exception.__cause__,
            NCBIXmlResponseValidationError,
        )
        efetch_mock.assert_called_once()
        sleep_mock.assert_not_called()

    def test_xml_fetch_configuration_rejects_invalid_controls(self) -> None:
        invalid_control_cases = [
            ("fetch_timeout_seconds", {"fetch_timeout_seconds": 0.0}),
            ("batch_deadline_seconds", {"batch_deadline_seconds": 0.0}),
            ("request_delay_seconds", {"request_delay_seconds": -0.01}),
            (
                "retry_backoff_initial_seconds",
                {"retry_backoff_initial_seconds": 0.0},
            ),
            ("retry_backoff_multiplier", {"retry_backoff_multiplier": 0.5}),
            ("retry_backoff_max_seconds", {"retry_backoff_max_seconds": 0.0}),
            (
                "circuit_breaker_failure_threshold",
                {"circuit_breaker_failure_threshold": 0},
            ),
            (
                "circuit_breaker_cooldown_seconds",
                {"circuit_breaker_cooldown_seconds": -1.0},
            ),
        ]

        for parameter_name, overrides in invalid_control_cases:
            with self.subTest(parameter_name=parameter_name):
                with self.assertRaisesRegex(ValueError, parameter_name):
                    fetch_ncbi_protein_xml_batches(
                        ncbi_email="test@example.org",
                        ncbi_api_key=None,
                        protein_uids=["1001"],
                        **overrides,
                    )


class NCBIApiSslConfigurationTests(unittest.TestCase):
    def test_resolve_ssl_ca_configuration_returns_none_when_no_inputs_are_provided(
        self,
    ) -> None:
        with patch.dict(os.environ, {}, clear=True):
            resolved_ssl_ca_file, resolved_ssl_ca_directory = (
                _resolve_ncbi_ssl_ca_configuration(
                    ssl_ca_file=None,
                    ssl_ca_directory=None,
                )
            )

        self.assertIsNone(resolved_ssl_ca_file)
        self.assertIsNone(resolved_ssl_ca_directory)

    def test_resolve_ssl_ca_configuration_uses_environment_variables(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            ssl_ca_file_path = temporary_directory / "corporate-ca.pem"
            ssl_ca_file_path.write_text("dummy-ca", encoding="utf-8")
            ssl_ca_directory_path = temporary_directory / "certs"
            ssl_ca_directory_path.mkdir()

            with patch.dict(
                os.environ,
                {
                    "NCBI_SSL_CERT_FILE": str(ssl_ca_file_path),
                    "NCBI_SSL_CERT_DIR": str(ssl_ca_directory_path),
                },
                clear=False,
            ):
                resolved_ssl_ca_file, resolved_ssl_ca_directory = (
                    _resolve_ncbi_ssl_ca_configuration(
                        ssl_ca_file=None,
                        ssl_ca_directory=None,
                    )
                )

        self.assertEqual(resolved_ssl_ca_file, str(ssl_ca_file_path))
        self.assertEqual(resolved_ssl_ca_directory, str(ssl_ca_directory_path))

    def test_ncbi_ssl_environment_variables_override_generic_ssl_environment(
        self,
    ) -> None:
        with patch.dict(
            os.environ,
            {
                "NCBI_SSL_CERT_FILE": "ncbi-file.pem",
                "NCBI_SSL_CERT_DIR": "ncbi-cert-dir",
                "SSL_CERT_FILE": "generic-file.pem",
                "SSL_CERT_DIR": "generic-cert-dir",
            },
            clear=False,
        ):
            resolved_ssl_ca_file, resolved_ssl_ca_directory = (
                _resolve_ncbi_ssl_ca_configuration(
                    ssl_ca_file=None,
                    ssl_ca_directory=None,
                )
            )

        self.assertEqual(resolved_ssl_ca_file, "ncbi-file.pem")
        self.assertEqual(resolved_ssl_ca_directory, "ncbi-cert-dir")

    def test_explicit_ssl_ca_configuration_overrides_environment(self) -> None:
        with patch.dict(
            os.environ,
            {
                "NCBI_SSL_CERT_FILE": "env-file.pem",
                "NCBI_SSL_CERT_DIR": "env-cert-dir",
            },
            clear=False,
        ):
            resolved_ssl_ca_file, resolved_ssl_ca_directory = (
                _resolve_ncbi_ssl_ca_configuration(
                    ssl_ca_file="explicit-file.pem",
                    ssl_ca_directory="explicit-cert-dir",
                )
            )

        self.assertEqual(resolved_ssl_ca_file, "explicit-file.pem")
        self.assertEqual(resolved_ssl_ca_directory, "explicit-cert-dir")

    def test_configured_entrez_urlopen_preserves_default_behavior_without_ca_configuration(
        self,
    ) -> None:
        captured_call: dict[str, object] = {}

        def fake_urlopen(*args, **kwargs):
            captured_call["args"] = args
            captured_call["context"] = kwargs.get("context")
            return "ok"

        with patch.object(ncbi_api.Entrez, "urlopen", fake_urlopen, create=True):
            with patch.dict(os.environ, {}, clear=True):
                with _configured_ncbi_entrez_urlopen(
                    ssl_ca_file=None,
                    ssl_ca_directory=None,
                ):
                    result = ncbi_api.Entrez.urlopen("request-object")

        self.assertEqual(result, "ok")
        self.assertEqual(captured_call["args"], ("request-object",))
        self.assertIsNone(captured_call["context"])

    def test_configured_entrez_urlopen_injects_ssl_context(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            ssl_ca_file_path = temporary_directory / "corporate-ca.pem"
            ssl_ca_file_path.write_text("dummy-ca", encoding="utf-8")
            ssl_ca_directory_path = temporary_directory / "certs"
            ssl_ca_directory_path.mkdir()

            captured_call: dict[str, object] = {}

            def fake_urlopen(*args, **kwargs):
                captured_call["args"] = args
                captured_call["context"] = kwargs.get("context")
                return "ok"

            ssl_context = object()

            with patch.object(ncbi_api.Entrez, "urlopen", fake_urlopen, create=True):
                with patch(
                    "src.pago_pipeline.ncbi_api.ssl.create_default_context",
                    return_value=ssl_context,
                ) as create_default_context:
                    with _configured_ncbi_entrez_urlopen(
                        ssl_ca_file=str(ssl_ca_file_path),
                        ssl_ca_directory=str(ssl_ca_directory_path),
                    ):
                        result = ncbi_api.Entrez.urlopen("request-object")

        self.assertEqual(result, "ok")
        self.assertEqual(captured_call["args"], ("request-object",))
        self.assertIs(captured_call["context"], ssl_context)
        create_default_context.assert_called_once_with(
            cafile=str(ssl_ca_file_path),
            capath=str(ssl_ca_directory_path),
        )

    def test_configured_entrez_urlopen_resets_active_context_after_exception(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            ssl_ca_file_path = temporary_directory / "corporate-ca.pem"
            ssl_ca_file_path.write_text("dummy-ca", encoding="utf-8")
            ssl_ca_directory_path = temporary_directory / "certs"
            ssl_ca_directory_path.mkdir()

            captured_contexts: list[object | None] = []

            def fake_urlopen(*args, **kwargs):
                captured_contexts.append(kwargs.get("context"))
                return "ok"

            ssl_context = object()

            with patch.object(ncbi_api.Entrez, "urlopen", fake_urlopen, create=True):
                with patch(
                    "src.pago_pipeline.ncbi_api.ssl.create_default_context",
                    return_value=ssl_context,
                ):
                    with self.assertRaisesRegex(RuntimeError, "forced failure"):
                        with _configured_ncbi_entrez_urlopen(
                            ssl_ca_file=str(ssl_ca_file_path),
                            ssl_ca_directory=str(ssl_ca_directory_path),
                        ):
                            ncbi_api.Entrez.urlopen("during-error")
                            raise RuntimeError("forced failure")

                    ncbi_api.Entrez.urlopen("after-error")

        self.assertEqual(captured_contexts, [ssl_context, None])

    def test_fetch_uid_snapshot_raises_clear_error_for_invalid_ssl_ca_file(self) -> None:
        with self.assertRaisesRegex(
            ValueError,
            "ssl_ca_file does not exist",
        ):
            fetch_ncbi_protein_uid_snapshot(
                ncbi_email="test@example.org",
                ncbi_api_key=None,
                query="txid9606[Organism]",
                ssl_ca_file="C:\\definitely-missing-corporate-ca.pem",
            )

    def test_fetch_uid_snapshot_raises_clear_error_for_invalid_ssl_ca_directory(
        self,
    ) -> None:
        with self.assertRaisesRegex(
            ValueError,
            "ssl_ca_directory does not exist",
        ):
            fetch_ncbi_protein_uid_snapshot(
                ncbi_email="test@example.org",
                ncbi_api_key=None,
                query="txid9606[Organism]",
                ssl_ca_directory="C:\\definitely-missing-corporate-certs",
            )

    def test_fetch_uid_snapshot_raises_clear_error_when_ssl_ca_file_is_directory(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            with self.assertRaisesRegex(
                ValueError,
                "ssl_ca_file is not a file",
            ):
                fetch_ncbi_protein_uid_snapshot(
                    ncbi_email="test@example.org",
                    ncbi_api_key=None,
                    query="txid9606[Organism]",
                    ssl_ca_file=temporary_directory_name,
                )

    def test_fetch_uid_snapshot_raises_clear_error_when_ssl_ca_directory_is_file(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            ssl_ca_file_path = Path(temporary_directory_name) / "corporate-ca.pem"
            ssl_ca_file_path.write_text("dummy-ca", encoding="utf-8")

            with self.assertRaisesRegex(
                ValueError,
                "ssl_ca_directory is not a directory",
            ):
                fetch_ncbi_protein_uid_snapshot(
                    ncbi_email="test@example.org",
                    ncbi_api_key=None,
                    query="txid9606[Organism]",
                    ssl_ca_directory=str(ssl_ca_file_path),
                )

    def test_fetch_uid_snapshot_raises_clear_error_for_certificate_verification_failure(
        self,
    ) -> None:
        with patch.object(
            ncbi_api.Entrez,
            "esearch",
            side_effect=URLError(
                "CERTIFICATE_VERIFY_FAILED: self-signed certificate in certificate chain"
            ),
            create=True,
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "NCBI_SSL_CERT_FILE / NCBI_SSL_CERT_DIR",
            ):
                fetch_ncbi_protein_uid_snapshot(
                    ncbi_email="test@example.org",
                    ncbi_api_key=None,
                    query="txid9606[Organism]",
                )


class ResolveNcbiSnapshotSslPropagationTests(unittest.TestCase):
    def test_uid_snapshot_resolve_propagates_ssl_ca_arguments(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_root_directory = Path(temporary_directory_name) / "uid_snapshot"

            with patch(
                "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_uid_snapshot",
                return_value=Mock(),
            ) as fetch_mock:
                with patch(
                    "src.pago_pipeline.ncbi_snapshot.save_ncbi_protein_uid_snapshot",
                    return_value=snapshot_root_directory / "snapshots" / "snapshot-1",
                ):
                    with patch(
                        "src.pago_pipeline.ncbi_snapshot.load_snapshot_by_directory",
                        return_value={"manifest": {}, "protein_uids": []},
                    ):
                        resolve_ncbi_protein_uid_snapshot(
                            snapshot_mode="create_new",
                            snapshot_root_directory=snapshot_root_directory,
                            search_query="txid9606[Organism]",
                            ssl_ca_file="corporate-ca.pem",
                            ssl_ca_directory="corporate-certs",
                            ncbi_email="test@example.org",
                        )

        self.assertEqual(
            fetch_mock.call_args.kwargs["ssl_ca_file"],
            "corporate-ca.pem",
        )
        self.assertEqual(
            fetch_mock.call_args.kwargs["ssl_ca_directory"],
            "corporate-certs",
        )

    def test_xml_snapshot_resolve_propagates_fetch_policy_and_shared_breaker(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            snapshot_root_directory = temporary_directory / "xml_snapshot"
            source_uid_snapshot_root_directory = temporary_directory / "uid_snapshot"

            source_uid_snapshot_payload = {
                "protein_uids": ["1001", "1002"],
                "manifest": {
                    "protein_uids_sha256": sha256_of_lines(text_lines=["1001", "1002"]),
                    "normalized_protein_uid_count": 2,
                },
                "manifest_file_path": source_uid_snapshot_root_directory / "manifest.json",
            }
            manifest_bytes = json.dumps(source_uid_snapshot_payload["manifest"]).encode("utf-8")
            source_uid_snapshot_payload["manifest_file_bytes"] = manifest_bytes
            source_uid_snapshot_payload["manifest_file_sha256"] = hashlib.sha256(manifest_bytes).hexdigest()

            with patch(
                "src.pago_pipeline.ncbi_snapshot.resolve_ncbi_protein_uid_snapshot",
                return_value=source_uid_snapshot_payload,
            ):
                with patch(
                    "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches",
                    return_value=Mock(identifier_type="uid"),
                ) as fetch_mock:
                    with patch(
                        "src.pago_pipeline.ncbi_snapshot.save_ncbi_protein_xml_snapshot",
                        return_value=snapshot_root_directory / "snapshots" / "snapshot-1",
                    ):
                        with patch(
                            "src.pago_pipeline.ncbi_snapshot.load_xml_snapshot_by_directory",
                            return_value={"manifest": {}, "protein_uids": ["1001", "1002"]},
                        ):
                            resolve_ncbi_protein_xml_snapshot(
                                snapshot_mode="create_new",
                                snapshot_root_directory=snapshot_root_directory,
                                source_uid_snapshot_root_directory=(
                                    source_uid_snapshot_root_directory
                                ),
                                ssl_ca_file="corporate-ca.pem",
                                ssl_ca_directory="corporate-certs",
                                fetch_timeout_seconds=12.0,
                                batch_deadline_seconds=45.0,
                                retry_backoff_initial_seconds=0.25,
                                retry_backoff_multiplier=3.0,
                                retry_backoff_max_seconds=4.0,
                                circuit_breaker_failure_threshold=2,
                                circuit_breaker_cooldown_seconds=30.0,
                                ncbi_email="test@example.org",
                            )
                            resolve_ncbi_protein_xml_snapshot(
                                snapshot_mode="create_new",
                                snapshot_root_directory=snapshot_root_directory,
                                source_uid_snapshot_root_directory=(
                                    source_uid_snapshot_root_directory
                                ),
                                circuit_breaker_failure_threshold=2,
                                circuit_breaker_cooldown_seconds=30.0,
                                ncbi_email="test@example.org",
                            )

        self.assertEqual(
            fetch_mock.call_args_list[0].kwargs["ssl_ca_file"],
            "corporate-ca.pem",
        )
        self.assertEqual(
            fetch_mock.call_args_list[0].kwargs["ssl_ca_directory"],
            "corporate-certs",
        )
        self.assertEqual(
            fetch_mock.call_args_list[0].kwargs["identifier_type"],
            NCBIProteinIdentifierType.uid,
        )
        self.assertEqual(
            fetch_mock.call_args_list[0].kwargs["fetch_timeout_seconds"],
            12.0,
        )
        self.assertEqual(
            fetch_mock.call_args_list[0].kwargs["batch_deadline_seconds"],
            45.0,
        )
        self.assertEqual(
            fetch_mock.call_args_list[0].kwargs["retry_backoff_initial_seconds"],
            0.25,
        )
        self.assertEqual(
            fetch_mock.call_args_list[0].kwargs["retry_backoff_multiplier"],
            3.0,
        )
        self.assertEqual(
            fetch_mock.call_args_list[0].kwargs["retry_backoff_max_seconds"],
            4.0,
        )
        self.assertIs(
            fetch_mock.call_args_list[0].kwargs["circuit_breaker"],
            fetch_mock.call_args_list[1].kwargs["circuit_breaker"],
        )

    def test_xml_snapshot_resolve_propagates_source_identifier_type(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            snapshot_root_directory = temporary_directory / "xml_snapshot"
            source_uid_snapshot_root_directory = temporary_directory / "uid_snapshot"

            source_uid_snapshot_payload = {
                "protein_uids": ["WP_012345678.1"],
                "manifest": {
                    "identifier_type": "accession.version",
                    "protein_identifiers_sha256": sha256_of_lines(text_lines=["WP_012345678.1"]),
                    "normalized_protein_identifier_count": 1,
                },
                "manifest_file_path": source_uid_snapshot_root_directory / "manifest.json",
            }
            manifest_bytes = json.dumps(source_uid_snapshot_payload["manifest"]).encode("utf-8")
            source_uid_snapshot_payload["manifest_file_bytes"] = manifest_bytes
            source_uid_snapshot_payload["manifest_file_sha256"] = hashlib.sha256(manifest_bytes).hexdigest()

            with patch(
                "src.pago_pipeline.ncbi_snapshot.resolve_ncbi_protein_uid_snapshot",
                return_value=source_uid_snapshot_payload,
            ):
                with patch(
                    "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches",
                    return_value=Mock(identifier_type="accession.version"),
                ) as fetch_mock:
                    with patch(
                        "src.pago_pipeline.ncbi_snapshot.save_ncbi_protein_xml_snapshot",
                        return_value=snapshot_root_directory / "snapshots" / "snapshot-1",
                    ):
                        with patch(
                            "src.pago_pipeline.ncbi_snapshot.load_xml_snapshot_by_directory",
                            return_value={
                                "manifest": {},
                                "protein_uids": ["WP_012345678.1"],
                            },
                        ):
                            resolve_ncbi_protein_xml_snapshot(
                                snapshot_mode="create_new",
                                snapshot_root_directory=snapshot_root_directory,
                                source_uid_snapshot_root_directory=(
                                    source_uid_snapshot_root_directory
                                ),
                                ncbi_email="test@example.org",
                            )

        self.assertEqual(
            fetch_mock.call_args.kwargs["identifier_type"],
            NCBIProteinIdentifierType.accession_version,
        )


if __name__ == "__main__":
    unittest.main()
