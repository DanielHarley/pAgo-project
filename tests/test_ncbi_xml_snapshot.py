import os
import sys
import tempfile
import types
import unittest
import json
from pathlib import Path
from unittest.mock import Mock, patch
from urllib.error import URLError

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
from src.pago_pipeline.ncbi_snapshot import (
    build_snapshot_directory_name,
    resolve_ncbi_protein_uid_snapshot,
    resolve_ncbi_protein_xml_snapshot,
    save_ncbi_protein_xml_snapshot,
    _validate_saved_consolidated_xml_snapshot,
)
from src.pago_pipeline.ncbi_api import (
    NCBIProteinXmlBatchFetchResult,
    NCBIProteinXmlFetchResult,
    _configured_ncbi_entrez_urlopen,
    _resolve_ncbi_ssl_ca_configuration,
    fetch_ncbi_protein_uid_snapshot,
)
from src.pago_pipeline.storage import sha256_of_lines


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

    def test_xml_snapshot_resolve_propagates_ssl_ca_arguments(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            snapshot_root_directory = temporary_directory / "xml_snapshot"
            source_uid_snapshot_root_directory = temporary_directory / "uid_snapshot"

            source_uid_snapshot_payload = {
                "protein_uids": ["1001", "1002"],
                "manifest": {},
                "manifest_file_path": source_uid_snapshot_root_directory / "manifest.json",
            }

            with patch(
                "src.pago_pipeline.ncbi_snapshot.resolve_ncbi_protein_uid_snapshot",
                return_value=source_uid_snapshot_payload,
            ):
                with patch(
                    "src.pago_pipeline.ncbi_snapshot.fetch_ncbi_protein_xml_batches",
                    return_value=Mock(),
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
                                source_uid_snapshot_root_directory=source_uid_snapshot_root_directory,
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


if __name__ == "__main__":
    unittest.main()
