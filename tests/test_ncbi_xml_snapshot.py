import sys
import tempfile
import types
import unittest
from pathlib import Path

if "Bio" not in sys.modules:
    fake_bio_module = types.ModuleType("Bio")
    fake_bio_module.__version__ = "test"
    fake_entrez_module = types.ModuleType("Bio.Entrez")
    fake_bio_module.Entrez = fake_entrez_module
    sys.modules["Bio"] = fake_bio_module
    sys.modules["Bio.Entrez"] = fake_entrez_module

from src.pago_pipeline.ncbi_snapshot import (
    _validate_saved_consolidated_xml_snapshot,
)


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


if __name__ == "__main__":
    unittest.main()
