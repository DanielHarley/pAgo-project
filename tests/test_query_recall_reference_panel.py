from __future__ import annotations

import csv
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from src.pago_pipeline.query_recall_reference_panel import (
    DEFAULT_PANEL_DIRECTORY,
    validate_query_recall_reference_panel,
)

_LOCK_FILE_NAME = "panel_lock.json"
_PANEL_FILE_NAME = "query_recall_reference_panel.csv"


def _copy_committed_panel(destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for name in ("panel_lock.json", "query_recall_reference_panel.csv", "curation_notes.md"):
        shutil.copy2(DEFAULT_PANEL_DIRECTORY / name, destination / name)


def _read_rows(panel_directory: Path) -> list[dict[str, str]]:
    with (panel_directory / _PANEL_FILE_NAME).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        return list(csv.DictReader(handle))


def _write_rows(panel_directory: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = list(rows[0].keys())
    with (panel_directory / _PANEL_FILE_NAME).open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _refresh_panel_hash_in_lock(panel_directory: Path) -> None:
    import hashlib

    lock_path = panel_directory / _LOCK_FILE_NAME
    payload = json.loads(lock_path.read_text(encoding="utf-8"))
    panel_bytes = (panel_directory / _PANEL_FILE_NAME).read_bytes()
    payload["panel_file_sha256"] = hashlib.sha256(panel_bytes).hexdigest()
    lock_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _mutate_lock(panel_directory: Path, **overrides: object) -> None:
    lock_path = panel_directory / _LOCK_FILE_NAME
    payload = json.loads(lock_path.read_text(encoding="utf-8"))
    payload.update(overrides)
    lock_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


class QueryRecallReferencePanelCommittedResourceTests(unittest.TestCase):
    """Validate the actually-committed panel resources, not a fixture copy."""

    def test_committed_panel_validates_offline(self) -> None:
        result = validate_query_recall_reference_panel()
        self.assertEqual(len(result.rows), 21)

    def test_lock_integrity(self) -> None:
        result = validate_query_recall_reference_panel()
        self.assertEqual(result.lock_payload["artifact_type"], "query_recall_reference_panel_lock")
        self.assertEqual(result.lock_payload["format_version"], "1.0")
        self.assertEqual(result.lock_payload["reference_count"], 21)

    def test_schema(self) -> None:
        result = validate_query_recall_reference_panel()
        for row in result.rows:
            self.assertTrue(row.accession)
            self.assertIsInstance(row.sequence_length, int)

    def test_reference_counts(self) -> None:
        result = validate_query_recall_reference_panel()
        families = [row.ago_family for row in result.rows]
        self.assertEqual(families.count("PAGO"), 14)
        self.assertEqual(families.count("PIWI_RE"), 7)

    def test_piwi_re_clade_semantics(self) -> None:
        result = validate_query_recall_reference_panel()
        for row in result.rows:
            if row.ago_family == "PIWI_RE":
                self.assertEqual(row.pago_clade, "UNRESOLVED")

    def test_evidence_tier_values(self) -> None:
        result = validate_query_recall_reference_panel()
        allowed = {
            "EXPERIMENTAL",
            "LITERATURE_PHYLOGENETIC",
            "CURATED_COMPUTATIONAL",
            "DATABASE_ANNOTATION",
        }
        for row in result.rows:
            self.assertIn(row.reference_label_evidence, allowed)

    def test_sequence_hash_format(self) -> None:
        result = validate_query_recall_reference_panel()
        for row in result.rows:
            self.assertRegex(row.sequence_sha256, r"^[0-9a-f]{64}$")

    def test_sequence_length_validation(self) -> None:
        result = validate_query_recall_reference_panel()
        for row in result.rows:
            self.assertGreater(row.sequence_length, 0)

    def test_use_restriction_is_machine_readable(self) -> None:
        result = validate_query_recall_reference_panel()
        self.assertEqual(result.lock_payload["intended_use"], "NCBI_QUERY_COVERAGE_AUDIT")
        self.assertFalse(result.lock_payload["detector_validation_eligible"])
        forbidden = set(result.lock_payload["forbidden_uses"])
        self.assertIn("DETECTOR_VALIDATION", forbidden)
        self.assertIn("DETECTOR_FINAL_HOLDOUT", forbidden)
        self.assertIn("DETECTOR_CALIBRATION", forbidden)
        self.assertIn("CLASSIFIER_GOLD_STANDARD", forbidden)

    def test_afago_reconciliation(self) -> None:
        result = validate_query_recall_reference_panel()
        row = next(r for r in result.rows if r.accession == "WP_010878815.1")
        self.assertEqual(row.protein_short_name, "AfAgo")
        self.assertEqual(row.ago_family, "PAGO")
        self.assertEqual(row.pago_clade, "UNRESOLVED")
        self.assertEqual(row.historical_query_recall_clade, "SHORT")

    def test_siago_reconciliation(self) -> None:
        result = validate_query_recall_reference_panel()
        row = next(r for r in result.rows if r.accession == "WP_012735993.1")
        self.assertEqual(row.protein_short_name, "SiAgo")
        self.assertEqual(row.ago_family, "PAGO")
        self.assertEqual(row.pago_clade, "UNRESOLVED")
        self.assertEqual(row.historical_query_recall_clade, "SHORT")

    def test_ngago_reconciliation(self) -> None:
        result = validate_query_recall_reference_panel()
        row = next(r for r in result.rows if r.accession == "WP_005580376.1")
        self.assertEqual(row.protein_short_name, "NgAgo")
        self.assertEqual(row.ago_family, "PAGO")
        self.assertEqual(row.pago_clade, "LONG_A")
        self.assertEqual(row.historical_query_recall_clade, "LONG_B")

    def test_ttago_verification_status_is_unchanged_from_historical(self) -> None:
        # No authoritative decision changes TtAgo's evidence-strength field;
        # it must stay at its historical value, independent of the accession-
        # identity question tracked below.
        result = validate_query_recall_reference_panel()
        row = next(r for r in result.rows if r.accession == "WP_011174533.1")
        self.assertEqual(row.protein_short_name, "TtAgo")
        self.assertEqual(row.pago_clade, "LONG_A")
        self.assertEqual(row.verification_status, "verified")

    def test_ttago_accession_curation_status_is_unresolved(self) -> None:
        # The open accession-identity question (SCIENTIFIC_STATE.md, Open
        # scientific question 2) is tracked independently of verification_status.
        result = validate_query_recall_reference_panel()
        row = next(r for r in result.rows if r.accession == "WP_011174533.1")
        self.assertEqual(row.accession_curation_status, "UNRESOLVED")
        # Independence: verification_status is unaffected by this open question.
        self.assertEqual(row.verification_status, "verified")

    def test_accession_curation_status_values(self) -> None:
        result = validate_query_recall_reference_panel()
        allowed = {"RESOLVED", "UNRESOLVED"}
        for row in result.rows:
            self.assertIn(row.accession_curation_status, allowed)
        unresolved = [
            row.accession
            for row in result.rows
            if row.accession_curation_status == "UNRESOLVED"
        ]
        self.assertEqual(unresolved, ["WP_011174533.1"])

    def test_no_piwi_re_pago_clade_value_exists_anywhere(self) -> None:
        result = validate_query_recall_reference_panel()
        for row in result.rows:
            self.assertNotEqual(row.pago_clade, "PIWI_RE")

    def test_declared_lock_descriptions_match_the_csv(self) -> None:
        # The lock's descriptive fields describe the committed CSV and are
        # recomputable offline, so a stale count, schema, or allowed-value set
        # must not survive validation. DATABASE_ANNOTATION is declared as an
        # explicit zero and is checked like any other value.
        result = validate_query_recall_reference_panel()
        lock = result.lock_payload
        self.assertEqual(lock["counts_by_ago_family"], {"PAGO": 14, "PIWI_RE": 7})
        self.assertEqual(
            lock["counts_by_evidence_tier"],
            {
                "CURATED_COMPUTATIONAL": 6,
                "DATABASE_ANNOTATION": 0,
                "EXPERIMENTAL": 12,
                "LITERATURE_PHYLOGENETIC": 3,
            },
        )
        self.assertEqual(
            lock["counts_by_verification_status"], {"provisional": 5, "verified": 16}
        )
        self.assertEqual(
            lock["counts_by_accession_curation_status"],
            {"RESOLVED": 20, "UNRESOLVED": 1},
        )
        self.assertEqual(
            tuple(lock["schema_fields"]),
            (
                "accession",
                "protein_short_name",
                "organism",
                "ago_family",
                "pago_clade",
                "historical_query_recall_clade",
                "sequence_sha256",
                "sequence_length",
                "uniprot_accession",
                "reference_label_source",
                "reference_label_evidence",
                "verification_status",
                "accession_curation_status",
                "notes",
            ),
        )
        self.assertEqual(set(lock["allowed_ago_family_values"]), {"PAGO", "PIWI_RE"})
        self.assertEqual(
            set(lock["allowed_pago_clade_values"]),
            {"LONG_A", "LONG_B", "SHORT", "UNRESOLVED"},
        )
        self.assertEqual(
            set(lock["allowed_evidence_tier_values"]),
            {
                "EXPERIMENTAL",
                "LITERATURE_PHYLOGENETIC",
                "CURATED_COMPUTATIONAL",
                "DATABASE_ANNOTATION",
            },
        )
        self.assertEqual(
            set(lock["allowed_verification_status_values"]),
            {"verified", "provisional"},
        )
        self.assertEqual(
            set(lock["allowed_accession_curation_status_values"]),
            {"RESOLVED", "UNRESOLVED"},
        )


class QueryRecallReferencePanelRejectionTests(unittest.TestCase):
    """Negative tests: corrupted or reverted inputs must be rejected."""

    def test_rejects_duplicate_accession(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            rows = _read_rows(panel_directory)
            rows[1]["accession"] = rows[0]["accession"]
            _write_rows(panel_directory, rows)
            _refresh_panel_hash_in_lock(panel_directory)
            with self.assertRaisesRegex(RuntimeError, "duplicate accessions"):
                validate_query_recall_reference_panel(panel_directory=panel_directory)

    def test_rejects_bad_sequence_hash(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            rows = _read_rows(panel_directory)
            rows[0]["sequence_sha256"] = "not-a-valid-hash"
            _write_rows(panel_directory, rows)
            _refresh_panel_hash_in_lock(panel_directory)
            with self.assertRaisesRegex(RuntimeError, "sequence_sha256"):
                validate_query_recall_reference_panel(panel_directory=panel_directory)

    def test_rejects_piwi_re_with_non_unresolved_pago_clade(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            rows = _read_rows(panel_directory)
            for row in rows:
                if row["ago_family"] == "PIWI_RE":
                    row["pago_clade"] = "LONG_A"
                    break
            _write_rows(panel_directory, rows)
            _refresh_panel_hash_in_lock(panel_directory)
            with self.assertRaisesRegex(RuntimeError, "PIWI_RE rows must carry"):
                validate_query_recall_reference_panel(panel_directory=panel_directory)

    def test_rejects_unknown_evidence_tier(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            rows = _read_rows(panel_directory)
            rows[0]["reference_label_evidence"] = "MADE_UP_TIER"
            _write_rows(panel_directory, rows)
            _refresh_panel_hash_in_lock(panel_directory)
            with self.assertRaisesRegex(RuntimeError, "reference_label_evidence"):
                validate_query_recall_reference_panel(panel_directory=panel_directory)

    def test_rejects_missing_use_restriction(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            _mutate_lock(panel_directory, detector_validation_eligible=True)
            with self.assertRaisesRegex(RuntimeError, "detector_validation_eligible"):
                validate_query_recall_reference_panel(panel_directory=panel_directory)

    def test_rejects_ngago_reverting_to_historical_long_b(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            rows = _read_rows(panel_directory)
            for row in rows:
                if row["accession"] == "WP_005580376.1":
                    row["pago_clade"] = "LONG_B"
                    break
            _write_rows(panel_directory, rows)
            _refresh_panel_hash_in_lock(panel_directory)
            with self.assertRaisesRegex(RuntimeError, "WP_005580376.1"):
                validate_query_recall_reference_panel(panel_directory=panel_directory)

    def test_rejects_afago_reverting_to_historical_short(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            rows = _read_rows(panel_directory)
            for row in rows:
                if row["accession"] == "WP_010878815.1":
                    row["pago_clade"] = "SHORT"
                    break
            _write_rows(panel_directory, rows)
            _refresh_panel_hash_in_lock(panel_directory)
            with self.assertRaisesRegex(RuntimeError, "WP_010878815.1"):
                validate_query_recall_reference_panel(panel_directory=panel_directory)

    def test_rejects_siago_reverting_to_historical_short(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            rows = _read_rows(panel_directory)
            for row in rows:
                if row["accession"] == "WP_012735993.1":
                    row["pago_clade"] = "SHORT"
                    break
            _write_rows(panel_directory, rows)
            _refresh_panel_hash_in_lock(panel_directory)
            with self.assertRaisesRegex(RuntimeError, "WP_012735993.1"):
                validate_query_recall_reference_panel(panel_directory=panel_directory)

    def test_rejects_ttago_accession_curation_status_marked_resolved(self) -> None:
        # This must be caught independently of verification_status, which is
        # deliberately left untouched by this mutation.
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            rows = _read_rows(panel_directory)
            for row in rows:
                if row["accession"] == "WP_011174533.1":
                    row["accession_curation_status"] = "RESOLVED"
                    self.assertEqual(row["verification_status"], "verified")
                    break
            _write_rows(panel_directory, rows)
            _refresh_panel_hash_in_lock(panel_directory)
            with self.assertRaisesRegex(RuntimeError, "WP_011174533.1"):
                validate_query_recall_reference_panel(panel_directory=panel_directory)

    def test_rejects_unknown_accession_curation_status(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            rows = _read_rows(panel_directory)
            rows[0]["accession_curation_status"] = "MADE_UP_STATUS"
            _write_rows(panel_directory, rows)
            _refresh_panel_hash_in_lock(panel_directory)
            with self.assertRaisesRegex(RuntimeError, "accession_curation_status"):
                validate_query_recall_reference_panel(panel_directory=panel_directory)

    def test_rejects_lock_count_diverging_from_csv(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            _mutate_lock(
                panel_directory,
                counts_by_evidence_tier={
                    "EXPERIMENTAL": 11,
                    "LITERATURE_PHYLOGENETIC": 3,
                    "CURATED_COMPUTATIONAL": 6,
                    "DATABASE_ANNOTATION": 0,
                },
            )
            with self.assertRaisesRegex(RuntimeError, "counts_by_evidence_tier"):
                validate_query_recall_reference_panel(panel_directory=panel_directory)

    def test_rejects_lock_schema_fields_diverging_from_csv(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            _mutate_lock(
                panel_directory,
                schema_fields=[
                    "accession",
                    "protein_short_name",
                    "organism",
                    "ago_family",
                    "pago_clade",
                    "historical_query_recall_clade",
                    "sequence_sha256",
                    "sequence_length",
                    "uniprot_accession",
                    "reference_label_source",
                    "reference_label_evidence",
                    "verification_status",
                    "notes",
                ],
            )
            with self.assertRaisesRegex(RuntimeError, "schema_fields"):
                validate_query_recall_reference_panel(panel_directory=panel_directory)

    def test_rejects_lock_allowed_values_diverging_from_loader(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            panel_directory = Path(name)
            _copy_committed_panel(panel_directory)
            _mutate_lock(
                panel_directory,
                allowed_verification_status_values=[
                    "verified",
                    "provisional",
                    "MADE_UP_STATUS",
                ],
            )
            with self.assertRaisesRegex(
                RuntimeError, "allowed_verification_status_values"
            ):
                validate_query_recall_reference_panel(panel_directory=panel_directory)


if __name__ == "__main__":
    unittest.main()
