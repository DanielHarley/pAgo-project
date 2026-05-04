import csv
import json
import sys
import tempfile
import types
import unittest
import xml.etree.ElementTree as ET
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

from src.pago_pipeline.ncbi_metadata_csv import (
    export_ncbi_protein_metadata_csv,
    flatten_gbseq_metadata,
    inspect_ncbi_protein_metadata_schema,
)


class FlattenGbseqMetadataTests(unittest.TestCase):
    def test_collects_dynamic_metadata_from_gbseq_record(self) -> None:
        gbseq_element = ET.fromstring(
            """
<GBSeq>
  <GBSeq_locus>SEQ_A</GBSeq_locus>
  <GBSeq_length>321</GBSeq_length>
  <GBSeq_taxonomy>Bacteria; Proteobacteria</GBSeq_taxonomy>
  <GBSeq_keywords>
    <GBKeyword>argonaute</GBKeyword>
    <GBKeyword>piwi</GBKeyword>
  </GBSeq_keywords>
  <GBSeq_other-seqids>
    <GBSeqid>gi|1001|</GBSeqid>
    <GBSeqid>ref|WP_000001.1|</GBSeqid>
  </GBSeq_other-seqids>
  <GBSeq_references>
    <GBReference>
      <GBReference_reference>1</GBReference_reference>
      <GBReference_title>First title</GBReference_title>
      <GBReference_authors>
        <GBAuthor>Alice</GBAuthor>
        <GBAuthor>Bob</GBAuthor>
      </GBReference_authors>
    </GBReference>
    <GBReference>
      <GBReference_reference>2</GBReference_reference>
      <GBReference_journal>Journal A</GBReference_journal>
    </GBReference>
  </GBSeq_references>
  <GBSeq_feature-table>
    <GBFeature>
      <GBFeature_key>source</GBFeature_key>
      <GBFeature_location>1..321</GBFeature_location>
      <GBFeature_intervals>
        <GBInterval>
          <GBInterval_from>1</GBInterval_from>
          <GBInterval_to>321</GBInterval_to>
        </GBInterval>
      </GBFeature_intervals>
      <GBFeature_quals>
        <GBQualifier>
          <GBQualifier_name>organism</GBQualifier_name>
          <GBQualifier_value>Testus organismus</GBQualifier_value>
        </GBQualifier>
        <GBQualifier>
          <GBQualifier_name>db_xref</GBQualifier_name>
          <GBQualifier_value>taxon:123</GBQualifier_value>
        </GBQualifier>
      </GBFeature_quals>
    </GBFeature>
    <GBFeature>
      <GBFeature_key>CDS</GBFeature_key>
      <GBFeature_location>5..300</GBFeature_location>
      <GBFeature_intervals>
        <GBInterval>
          <GBInterval_from>5</GBInterval_from>
          <GBInterval_to>300</GBInterval_to>
          <GBInterval_accession>NC_000001</GBInterval_accession>
        </GBInterval>
      </GBFeature_intervals>
      <GBFeature_quals>
        <GBQualifier>
          <GBQualifier_name>coded_by</GBQualifier_name>
          <GBQualifier_value>NC_000001:5..300</GBQualifier_value>
        </GBQualifier>
      </GBFeature_quals>
    </GBFeature>
  </GBSeq_feature-table>
</GBSeq>
"""
        )

        row = flatten_gbseq_metadata(gbseq_element=gbseq_element)

        self.assertEqual(row["protein_uid"], "1001")
        self.assertEqual(row["gbseq__locus"], "SEQ_A")
        self.assertEqual(row["gbseq__length"], "321")
        self.assertEqual(row["taxonomy__raw"], "Bacteria; Proteobacteria")
        self.assertEqual(row["taxonomy__01"], "Bacteria")
        self.assertEqual(row["taxonomy__02"], "Proteobacteria")
        self.assertEqual(row["gbseq__keywords__keyword"], "argonaute;piwi")
        self.assertEqual(
            row["gbseq__other_seqids"],
            "gi|1001|;ref|WP_000001.1|",
        )
        self.assertEqual(row["reference__count"], "2")
        self.assertEqual(row["reference__reference"], "1;2")
        self.assertEqual(row["reference__authors"], "Alice;Bob")
        self.assertEqual(row["feature__keys_present"], "source;cds")
        self.assertEqual(row["feature__source__location"], "1..321")
        self.assertEqual(
            row["feature__source__qual__organism"],
            "Testus organismus",
        )
        self.assertEqual(row["feature__source__qual__db_xref"], "taxon:123")
        self.assertEqual(row["feature__cds__location"], "5..300")
        self.assertEqual(row["feature__cds__interval__from"], "5")
        self.assertEqual(row["feature__cds__interval__to"], "300")
        self.assertEqual(
            row["feature__cds__interval__accession"],
            "NC_000001",
        )
        self.assertEqual(
            row["feature__cds__qual__coded_by"],
            "NC_000001:5..300",
        )


class InspectNcbiProteinMetadataSchemaTests(unittest.TestCase):
    def _write_xml_file(self, xml_text: str) -> Path:
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)

        xml_file_path = Path(temporary_directory.name) / "protein_records.xml"
        xml_file_path.write_text(xml_text, encoding="utf-8")
        return xml_file_path

    def test_discovers_union_of_dynamic_columns_and_taxonomy_depth(self) -> None:
        xml_file_path = self._write_xml_file(
            """<?xml version="1.0" encoding="utf-8"?>
<GBSet>
  <GBSeq>
    <GBSeqid>gi|1001|</GBSeqid>
    <GBSeq_locus>SEQ_A</GBSeq_locus>
    <GBSeq_taxonomy>Bacteria; Proteobacteria</GBSeq_taxonomy>
    <GBSeq_feature-table>
      <GBFeature>
        <GBFeature_key>source</GBFeature_key>
        <GBFeature_quals>
          <GBQualifier>
            <GBQualifier_name>organism</GBQualifier_name>
            <GBQualifier_value>Species A</GBQualifier_value>
          </GBQualifier>
        </GBFeature_quals>
      </GBFeature>
    </GBSeq_feature-table>
  </GBSeq>
  <GBSeq>
    <GBSeqid>gi|1002|</GBSeqid>
    <GBSeq_locus>SEQ_B</GBSeq_locus>
    <GBSeq_taxonomy>Bacteria; Firmicutes; Bacilli</GBSeq_taxonomy>
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
    </GBSeq_feature-table>
  </GBSeq>
</GBSet>
"""
        )

        schema_metadata = inspect_ncbi_protein_metadata_schema(
            xml_file_path=xml_file_path,
        )

        self.assertEqual(schema_metadata.row_count, 2)
        self.assertEqual(schema_metadata.max_taxonomy_depth, 3)
        self.assertIn("protein_uid", schema_metadata.columns)
        self.assertIn("taxonomy__03", schema_metadata.columns)
        self.assertIn(
            "feature__source__qual__organism",
            schema_metadata.columns,
        )
        self.assertIn(
            "feature__region__qual__region_name",
            schema_metadata.columns,
        )
        self.assertEqual(
            schema_metadata.observed_feature_keys,
            ["region", "source"],
        )
        self.assertEqual(
            schema_metadata.observed_feature_qualifiers,
            {
                "region": ["region_name"],
                "source": ["organism"],
            },
        )


class ExportNcbiProteinMetadataCsvTests(unittest.TestCase):
    def test_writes_csv_and_manifest_with_dynamic_columns(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            xml_file_path = temporary_directory / "protein_records.xml"
            csv_file_path = temporary_directory / "protein_metadata.csv"
            manifest_file_path = temporary_directory / "manifest.json"

            xml_file_path.write_text(
                """<?xml version="1.0" encoding="utf-8"?>
<GBSet>
  <GBSeq>
    <GBSeqid>gi|1001|</GBSeqid>
    <GBSeq_locus>SEQ_A</GBSeq_locus>
    <GBSeq_taxonomy>Bacteria; Proteobacteria</GBSeq_taxonomy>
    <GBSeq_feature-table>
      <GBFeature>
        <GBFeature_key>source</GBFeature_key>
        <GBFeature_quals>
          <GBQualifier>
            <GBQualifier_name>organism</GBQualifier_name>
            <GBQualifier_value>Species A</GBQualifier_value>
          </GBQualifier>
        </GBFeature_quals>
      </GBFeature>
    </GBSeq_feature-table>
  </GBSeq>
  <GBSeq>
    <GBSeqid>gi|1002|</GBSeqid>
    <GBSeq_locus>SEQ_B</GBSeq_locus>
    <GBSeq_taxonomy>Bacteria; Firmicutes; Bacilli</GBSeq_taxonomy>
    <GBSeq_references>
      <GBReference>
        <GBReference_reference>1</GBReference_reference>
        <GBReference_title>Reference title</GBReference_title>
      </GBReference>
    </GBSeq_references>
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
    </GBSeq_feature-table>
  </GBSeq>
</GBSet>
""",
                encoding="utf-8",
            )

            export_result = export_ncbi_protein_metadata_csv(
                xml_file_path=xml_file_path,
                output_csv_file_path=csv_file_path,
                output_manifest_file_path=manifest_file_path,
            )

            self.assertEqual(export_result.row_count, 2)
            self.assertTrue(csv_file_path.exists())
            self.assertTrue(manifest_file_path.exists())

            with csv_file_path.open("r", encoding="utf-8", newline="") as csv_file:
                rows = list(csv.DictReader(csv_file))

            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["protein_uid"], "1001")
            self.assertEqual(rows[0]["taxonomy__03"], "")
            self.assertEqual(rows[0]["feature__source__qual__organism"], "Species A")
            self.assertEqual(rows[0]["feature__region__qual__region_name"], "")
            self.assertEqual(rows[1]["protein_uid"], "1002")
            self.assertEqual(rows[1]["taxonomy__03"], "Bacilli")
            self.assertEqual(rows[1]["reference__count"], "1")
            self.assertEqual(rows[1]["reference__title"], "Reference title")
            self.assertEqual(rows[1]["feature__region__qual__region_name"], "PAZ")

            manifest_payload = json.loads(
                manifest_file_path.read_text(encoding="utf-8")
            )
            self.assertEqual(
                manifest_payload["artifact_type"],
                "ncbi_protein_metadata_csv",
            )
            self.assertEqual(manifest_payload["row_count"], 2)
            self.assertEqual(manifest_payload["column_count"], len(rows[0].keys()))
            self.assertEqual(manifest_payload["max_taxonomy_depth"], 3)
            self.assertIn("feature__region__qual__region_name", manifest_payload["columns"])


if __name__ == "__main__":
    unittest.main()
