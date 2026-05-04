import csv
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path

if "Bio" not in sys.modules:
    fake_bio_module = types.ModuleType("Bio")
    fake_bio_module.__version__ = "test"
    fake_entrez_module = types.ModuleType("Bio.Entrez")
    fake_seqio_module = types.ModuleType("Bio.SeqIO")
    fake_bio_module.Entrez = fake_entrez_module
    fake_bio_module.SeqIO = fake_seqio_module
    sys.modules["Bio"] = fake_bio_module
    sys.modules["Bio.Entrez"] = fake_entrez_module
    sys.modules["Bio.SeqIO"] = fake_seqio_module

from src.pago_pipeline.ncbi_metadata_csv import export_ncbi_protein_metadata_csv
from src.pago_pipeline.ncbi_metadata_qc import run_ncbi_protein_metadata_csv_qc


class RunNcbiProteinMetadataCsvQcTests(unittest.TestCase):
    def test_qc_report_matches_manifest_and_source_xml(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            xml_file_path = temporary_directory / "protein_records.xml"
            csv_file_path = temporary_directory / "protein_metadata.csv"
            metadata_manifest_file_path = temporary_directory / "metadata_manifest.json"
            qc_report_file_path = temporary_directory / "qc_report.json"

            xml_file_path.write_text(
                """<?xml version="1.0" encoding="utf-8"?>
<GBSet>
  <GBSeq>
    <GBSeqid>gi|1001|</GBSeqid>
    <GBSeq_locus>SEQ_A</GBSeq_locus>
    <GBSeq_taxonomy>Bacteria; Proteobacteria</GBSeq_taxonomy>
    <GBSeq_feature-table>
      <GBFeature>
        <GBFeature_key>Region</GBFeature_key>
        <GBFeature_quals>
          <GBQualifier>
            <GBQualifier_name>region_name</GBQualifier_name>
            <GBQualifier_value>PAZ</GBQualifier_value>
          </GBQualifier>
        </GBFeature_quals>
      </GBFeature>
      <GBFeature>
        <GBFeature_key>Region</GBFeature_key>
        <GBFeature_quals>
          <GBQualifier>
            <GBQualifier_name>note</GBQualifier_name>
            <GBQualifier_value>Second region</GBQualifier_value>
          </GBQualifier>
        </GBFeature_quals>
      </GBFeature>
    </GBSeq_feature-table>
  </GBSeq>
  <GBSeq>
    <GBSeqid>gi|1002|</GBSeqid>
    <GBSeq_locus>SEQ_B</GBSeq_locus>
    <GBSeq_taxonomy>Bacteria; Firmicutes</GBSeq_taxonomy>
    <GBSeq_references>
      <GBReference>
        <GBReference_reference>1</GBReference_reference>
        <GBReference_title>Reference title</GBReference_title>
      </GBReference>
    </GBSeq_references>
  </GBSeq>
</GBSet>
""",
                encoding="utf-8",
            )

            export_ncbi_protein_metadata_csv(
                xml_file_path=xml_file_path,
                output_csv_file_path=csv_file_path,
                output_manifest_file_path=metadata_manifest_file_path,
            )

            qc_result = run_ncbi_protein_metadata_csv_qc(
                csv_file_path=csv_file_path,
                metadata_manifest_file_path=metadata_manifest_file_path,
                source_xml_file_path=xml_file_path,
                source_xml_manifest_payload={
                    "consolidated_record_count": 2,
                    "immutable_snapshot_relative_path": "snapshots/example",
                },
                output_report_file_path=qc_report_file_path,
            )

            self.assertTrue(qc_result.schema_match_with_manifest)
            self.assertTrue(qc_result.row_count_match_with_manifest)
            self.assertTrue(qc_result.row_count_match_with_source_xml)
            self.assertEqual(qc_result.duplicate_protein_uid_count, 0)
            self.assertEqual(qc_result.empty_protein_uid_count, 0)

            qc_report_payload = json.loads(
                qc_report_file_path.read_text(encoding="utf-8")
            )
            self.assertEqual(
                qc_report_payload["checks"]["protein_uid_has_no_duplicates"],
                True,
            )
            self.assertEqual(
                qc_report_payload["checks"]["row_count_matches_source_xml"],
                True,
            )
            self.assertIn(
                "region",
                qc_report_payload["xml_flattening_risks"]["repeated_feature_keys"],
            )

    def test_qc_detects_duplicate_ids_and_normalization_collisions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            xml_file_path = temporary_directory / "protein_records.xml"
            csv_file_path = temporary_directory / "protein_metadata.csv"
            metadata_manifest_file_path = temporary_directory / "metadata_manifest.json"
            qc_report_file_path = temporary_directory / "qc_report.json"

            xml_file_path.write_text(
                """<?xml version="1.0" encoding="utf-8"?>
<GBSet>
  <GBSeq>
    <GBSeqid>gi|1001|</GBSeqid>
    <GBSeq_locus>SEQ_A</GBSeq_locus>
    <GBSeq_feature-table>
      <GBFeature>
        <GBFeature_key>source</GBFeature_key>
        <GBFeature_quals>
          <GBQualifier>
            <GBQualifier_name>db-xref</GBQualifier_name>
            <GBQualifier_value>taxon:1</GBQualifier_value>
          </GBQualifier>
          <GBQualifier>
            <GBQualifier_name>db_xref</GBQualifier_name>
            <GBQualifier_value>taxon:2</GBQualifier_value>
          </GBQualifier>
        </GBFeature_quals>
      </GBFeature>
    </GBSeq_feature-table>
  </GBSeq>
</GBSet>
""",
                encoding="utf-8",
            )

            export_ncbi_protein_metadata_csv(
                xml_file_path=xml_file_path,
                output_csv_file_path=csv_file_path,
                output_manifest_file_path=metadata_manifest_file_path,
            )

            with csv_file_path.open("r", encoding="utf-8", newline="") as csv_file:
                rows = list(csv.DictReader(csv_file))
                fieldnames = list(rows[0].keys())

            duplicated_rows = [rows[0], dict(rows[0])]
            with csv_file_path.open("w", encoding="utf-8", newline="") as csv_file:
                csv_writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                csv_writer.writeheader()
                csv_writer.writerows(duplicated_rows)

            qc_result = run_ncbi_protein_metadata_csv_qc(
                csv_file_path=csv_file_path,
                metadata_manifest_file_path=metadata_manifest_file_path,
                source_xml_file_path=xml_file_path,
                source_xml_manifest_payload={"consolidated_record_count": 1},
                output_report_file_path=qc_report_file_path,
            )

            self.assertFalse(qc_result.row_count_match_with_manifest)
            self.assertFalse(qc_result.row_count_match_with_source_xml)
            self.assertEqual(qc_result.duplicate_protein_uid_count, 1)
            self.assertEqual(qc_result.normalization_collision_count, 1)

            qc_report_payload = json.loads(
                qc_report_file_path.read_text(encoding="utf-8")
            )
            self.assertEqual(
                qc_report_payload["checks"]["protein_uid_has_no_duplicates"],
                False,
            )
            self.assertIn(
                "db_xref",
                qc_report_payload["xml_flattening_risks"]["normalization_collisions"][
                    "feature_qualifier_names"
                ],
            )

    def test_qc_reports_schema_drift_against_reference_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            reference_xml_file_path = temporary_directory / "reference_records.xml"
            current_xml_file_path = temporary_directory / "current_records.xml"
            reference_csv_file_path = temporary_directory / "reference_metadata.csv"
            current_csv_file_path = temporary_directory / "current_metadata.csv"
            reference_manifest_file_path = temporary_directory / "reference_manifest.json"
            current_manifest_file_path = temporary_directory / "current_manifest.json"

            reference_xml_file_path.write_text(
                """<?xml version="1.0" encoding="utf-8"?>
<GBSet>
  <GBSeq>
    <GBSeqid>gi|1001|</GBSeqid>
    <GBSeq_locus>SEQ_A</GBSeq_locus>
    <GBSeq_feature-table>
      <GBFeature>
        <GBFeature_key>source</GBFeature_key>
      </GBFeature>
    </GBSeq_feature-table>
  </GBSeq>
</GBSet>
""",
                encoding="utf-8",
            )
            current_xml_file_path.write_text(
                """<?xml version="1.0" encoding="utf-8"?>
<GBSet>
  <GBSeq>
    <GBSeqid>gi|1001|</GBSeqid>
    <GBSeq_locus>SEQ_A</GBSeq_locus>
    <GBSeq_taxonomy>Bacteria; Proteobacteria; Gammaproteobacteria</GBSeq_taxonomy>
    <GBSeq_feature-table>
      <GBFeature>
        <GBFeature_key>source</GBFeature_key>
      </GBFeature>
      <GBFeature>
        <GBFeature_key>Region</GBFeature_key>
        <GBFeature_quals>
          <GBQualifier>
            <GBQualifier_name>region_name</GBQualifier_name>
            <GBQualifier_value>PIWI</GBQualifier_value>
          </GBQualifier>
        </GBFeature_quals>
      </GBFeature>
    </GBSeq_feature-table>
  </GBSeq>
</GBSet>
""",
                encoding="utf-8",
            )

            export_ncbi_protein_metadata_csv(
                xml_file_path=reference_xml_file_path,
                output_csv_file_path=reference_csv_file_path,
                output_manifest_file_path=reference_manifest_file_path,
            )
            export_ncbi_protein_metadata_csv(
                xml_file_path=current_xml_file_path,
                output_csv_file_path=current_csv_file_path,
                output_manifest_file_path=current_manifest_file_path,
            )

            qc_report = run_ncbi_protein_metadata_csv_qc(
                csv_file_path=current_csv_file_path,
                metadata_manifest_file_path=current_manifest_file_path,
                reference_metadata_manifest_file_path=reference_manifest_file_path,
            )

            self.assertTrue(qc_report.schema_match_with_manifest)

            current_qc_payload = run_ncbi_protein_metadata_csv_qc(
                csv_file_path=current_csv_file_path,
                metadata_manifest_file_path=current_manifest_file_path,
                reference_metadata_manifest_file_path=reference_manifest_file_path,
                output_report_file_path=temporary_directory / "schema_drift_qc.json",
            )
            self.assertTrue(current_qc_payload.schema_match_with_manifest)

            qc_report_payload = json.loads(
                (temporary_directory / "schema_drift_qc.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertIn("taxonomy__03", qc_report_payload["schema_drift"]["added_columns"])
            self.assertIn(
                "region",
                qc_report_payload["schema_drift"]["added_feature_keys"],
            )


if __name__ == "__main__":
    unittest.main()
