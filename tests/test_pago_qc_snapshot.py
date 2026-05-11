from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from src.pago_pipeline.pago_qc_snapshot import (
    DEFAULT_EVIDENCE_COUNTS_FILE_NAME,
    DEFAULT_EVIDENCE_FLAGS_FILE_NAME,
    DEFAULT_SUSPICIOUS_TERMS_FILE_NAME,
    DEFAULT_TOP_PRODUCTS_FILE_NAME,
    DEFAULT_TOP_REGION_NAMES_FILE_NAME,
    _validate_loaded_pago_qc_evidence_inventory_payload,
)
from src.pago_pipeline.storage import sha256_of_file


class PagoQcSnapshotValidationTests(unittest.TestCase):
    def test_pago_qc_inventory_validation_requires_output_file_sha256(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            output_file_names_by_key = {
                "evidence_flags_file": DEFAULT_EVIDENCE_FLAGS_FILE_NAME,
                "evidence_counts_file": DEFAULT_EVIDENCE_COUNTS_FILE_NAME,
                "top_region_names_file": DEFAULT_TOP_REGION_NAMES_FILE_NAME,
                "top_products_file": DEFAULT_TOP_PRODUCTS_FILE_NAME,
                "suspicious_terms_file": DEFAULT_SUSPICIOUS_TERMS_FILE_NAME,
            }
            for output_file_name in output_file_names_by_key.values():
                (temporary_directory / output_file_name).write_text(
                    "header\n",
                    encoding="utf-8",
                )

            manifest_payload = {
                "artifact_type": "pago_qc_evidence_inventory",
                "output_files": {
                    output_file_key: {
                        "file_name": output_file_name,
                        "sha256": sha256_of_file(
                            input_file_path=temporary_directory / output_file_name
                        ),
                    }
                    for output_file_key, output_file_name in (
                        output_file_names_by_key.items()
                    )
                },
            }
            del manifest_payload["output_files"]["evidence_counts_file"]["sha256"]

            with self.assertRaisesRegex(RuntimeError, "non-empty sha256"):
                _validate_loaded_pago_qc_evidence_inventory_payload(
                    snapshot_directory=temporary_directory,
                    manifest_payload=manifest_payload,
                )
