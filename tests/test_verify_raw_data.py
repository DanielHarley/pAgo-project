from __future__ import annotations

import hashlib
import io
import json
import os
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from scripts.verify_raw_data import _verify_uid_file, _verify_xml_file, verify_raw_data
from src.pago_pipeline.ncbi_snapshot import (
    load_snapshot_by_directory,
    load_xml_snapshot_by_directory,
    resolve_ncbi_protein_uid_snapshot,
    resolve_ncbi_protein_xml_snapshot,
)
from src.pago_pipeline.storage import sha256_of_lines


class VerifyRawDataTests(unittest.TestCase):
    def _write_integrity_snapshot(self, root: Path, identifier_type: str, include_xml: bool) -> tuple[Path, dict]:
        directory = root / "latest"
        directory.mkdir()
        identifiers = ["1001", "1002", "1001"]
        if identifier_type == "accession.version":
            identifiers = ["WP_012345678.1", "XP_123456789.2", "WP_012345678.1"]
        identifier_text = "\n".join(identifiers) + "\n"
        for filename in ("protein_identifiers.txt", "protein_uids.txt"):
            (directory / filename).write_text(identifier_text, encoding="utf-8")
        identifier_hash = sha256_of_lines(
            text_lines=identifiers, deduplicate_lines_preserving_order=False, sort_lines=False,
        )
        manifest = {
            "snapshot_format_version": "1.1",
            "artifact_type": "ncbi_protein_identifier_snapshot",
            "identifier_type": identifier_type,
            "protein_identifiers_file_name": "protein_identifiers.txt",
            "protein_identifiers_sha256": identifier_hash,
            "normalized_protein_identifier_count": len(identifiers),
            "protein_identifier_count": len(identifiers),
            "protein_uids_file_name": "protein_uids.txt",
            "protein_uids_sha256": identifier_hash,
            "normalized_protein_uid_count": len(identifiers),
        }
        if include_xml:
            records = []
            for identifier in identifiers[1:] + identifiers[:1]:
                identity = f"<GBSeqid>gi|{identifier}</GBSeqid>"
                if identifier_type == "accession.version":
                    identity = f"<GBSeq_accession-version>{identifier}</GBSeq_accession-version>"
                records.append(f"<GBSeq>{identity}<GBSeq_sequence>AAAA</GBSeq_sequence></GBSeq>")
            xml_bytes = ("<GBSet>" + "".join(records) + "</GBSet>").encode("utf-8")
            (directory / "protein_records.xml").write_bytes(xml_bytes)
            manifest.update({
                "artifact_type": "ncbi_protein_xml_snapshot",
                "xml_file_name": "protein_records.xml",
                "xml_file_sha256": hashlib.sha256(xml_bytes).hexdigest(),
                "consolidated_record_count": len(identifiers),
            })
        self._write_manifest(directory, manifest)
        return directory, manifest

    def test_identifier_paths_cannot_escape_snapshot_in_reuse_or_verification(self) -> None:
        for identifier_type in ("uid", "accession.version"):
            for include_xml in (False, True):
                with (
                    self.subTest(domain=identifier_type, xml=include_xml),
                    tempfile.TemporaryDirectory() as name,
                ):
                    root = Path(name)
                    directory, manifest = self._write_integrity_snapshot(root, identifier_type, include_xml)
                    outside_file = root / "outside.txt"
                    outside_file.write_bytes((directory / "protein_identifiers.txt").read_bytes())
                    load = load_xml_snapshot_by_directory if include_xml else load_snapshot_by_directory
                    resolve = resolve_ncbi_protein_xml_snapshot if include_xml else resolve_ncbi_protein_uid_snapshot
                    for file_name in (
                        "../outside.txt", r"..\outside.txt", str(outside_file),
                        r"C:outside.txt", r"\\server\share\outside.txt",
                    ):
                        with self.subTest(file_name=file_name):
                            manifest["protein_identifiers_file_name"] = file_name
                            self._write_manifest(directory, manifest)
                            with self.assertRaisesRegex(RuntimeError, "inside the snapshot directory"):
                                load(snapshot_directory=directory)
                            with self.assertRaisesRegex(RuntimeError, "inside the snapshot directory"):
                                resolve(snapshot_mode="reuse_latest_or_create", snapshot_root_directory=root)
                            status, _, errors = self._run_verifier(root)
                            self.assertEqual(status, 1)
                            self.assertIn("inside the snapshot directory", errors)

    def test_xml_paths_cannot_escape_snapshot_in_reuse_or_verification(self) -> None:
        for identifier_type in ("uid", "accession.version"):
            with self.subTest(domain=identifier_type), tempfile.TemporaryDirectory() as name:
                root = Path(name)
                directory, manifest = self._write_integrity_snapshot(root, identifier_type, True)
                outside_file = root / "outside.xml"
                (directory / "protein_records.xml").rename(outside_file)
                for file_name in (
                    "../outside.xml", r"..\outside.xml", str(outside_file),
                    r"C:outside.xml", r"\\server\share\outside.xml",
                ):
                    with self.subTest(file_name=file_name):
                        manifest["xml_file_name"] = file_name
                        self._write_manifest(directory, manifest)
                        with self.assertRaisesRegex(RuntimeError, "inside the snapshot directory"):
                            load_xml_snapshot_by_directory(snapshot_directory=directory)
                        for mode in ("reuse_latest", "reuse_latest_or_create"):
                            with self.assertRaisesRegex(RuntimeError, "inside the snapshot directory"):
                                resolve_ncbi_protein_xml_snapshot(
                                    snapshot_mode=mode, snapshot_root_directory=root,
                                )
                        status, _, errors = self._run_verifier(root)
                        self.assertEqual(status, 1)
                        self.assertIn("inside the snapshot directory", errors)

    def test_identifier_file_may_use_a_contained_relative_subdirectory(self) -> None:
        for identifier_type in ("uid", "accession.version"):
            with self.subTest(domain=identifier_type), tempfile.TemporaryDirectory() as name:
                root = Path(name)
                directory, manifest = self._write_integrity_snapshot(root, identifier_type, True)
                subdirectory = directory / "lists"
                subdirectory.mkdir()
                (directory / "protein_identifiers.txt").rename(subdirectory / "identifiers.txt")
                manifest["protein_identifiers_file_name"] = r"lists\identifiers.txt"
                self._write_manifest(directory, manifest)
                payload = load_xml_snapshot_by_directory(snapshot_directory=directory)
                self.assertEqual(payload["protein_uids_file_path"], subdirectory / "identifiers.txt")
                self.assertEqual(self._run_verifier(root)[0], 0)

    def test_identifier_symlink_cannot_target_a_file_outside_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            directory, manifest = self._write_integrity_snapshot(root, "uid", False)
            outside_file = root / "outside.txt"
            outside_file.write_bytes((directory / "protein_identifiers.txt").read_bytes())
            link = directory / "linked_identifiers.txt"
            try:
                link.symlink_to(outside_file)
            except OSError as error:
                if os.name == "nt" and error.winerror in (5, 1314):
                    self.skipTest("Windows does not grant symbolic link creation in this environment.")
                raise
            manifest["protein_identifiers_file_name"] = link.name
            self._write_manifest(directory, manifest)
            with self.assertRaisesRegex(RuntimeError, "inside the snapshot directory"):
                load_snapshot_by_directory(snapshot_directory=directory)
            status, _, errors = self._run_verifier(root)
            self.assertEqual(status, 1)
            self.assertIn("inside the snapshot directory", errors)

    def test_1_1_reuse_and_verifier_reject_removed_integrity_fields(self) -> None:
        for identifier_type in ("uid", "accession.version"):
            for include_xml in (False, True):
                fields = ["protein_identifiers_sha256", "normalized_protein_identifier_count"]
                if include_xml:
                    fields.extend(("xml_file_sha256", "consolidated_record_count"))
                for field in fields:
                    with (
                        self.subTest(domain=identifier_type, xml=include_xml, field=field),
                        tempfile.TemporaryDirectory() as name,
                    ):
                        root = Path(name)
                        directory, manifest = self._write_integrity_snapshot(root, identifier_type, include_xml)
                        manifest.pop(field)
                        self._write_manifest(directory, manifest)
                        original_manifest = (directory / "manifest.json").read_bytes()
                        load = load_xml_snapshot_by_directory if include_xml else load_snapshot_by_directory
                        with self.assertRaisesRegex(RuntimeError, field):
                            load(snapshot_directory=directory)
                        resolve = resolve_ncbi_protein_xml_snapshot if include_xml else resolve_ncbi_protein_uid_snapshot
                        for mode in ("reuse_latest", "reuse_latest_or_create"):
                            with self.assertRaisesRegex(RuntimeError, field):
                                resolve(snapshot_mode=mode, snapshot_root_directory=root)
                        status, output, errors = self._run_verifier(root)
                        self.assertEqual(status, 1)
                        self.assertIn(field, errors)
                        self.assertNotIn("Verified", output)
                        self.assertEqual((directory / "manifest.json").read_bytes(), original_manifest)

    def test_1_1_reuse_and_verifier_reject_invalid_hashes_and_counts(self) -> None:
        for identifier_type in ("uid", "accession.version"):
            for field, values in (
                ("protein_identifiers_sha256", (None, "", "0" * 63, "g" * 64, 123)),
                ("xml_file_sha256", (None, "", "0" * 63, "g" * 64, 123)),
                ("normalized_protein_identifier_count", (None, True, 3.0, "3", -1, 999)),
                ("consolidated_record_count", (None, True, 3.0, "3", -1, 999)),
            ):
                for value in values:
                    with (
                        self.subTest(domain=identifier_type, field=field, value=value),
                        tempfile.TemporaryDirectory() as name,
                    ):
                        root = Path(name)
                        directory, manifest = self._write_integrity_snapshot(root, identifier_type, True)
                        manifest[field] = value
                        self._write_manifest(directory, manifest)
                        with self.assertRaises(RuntimeError):
                            load_xml_snapshot_by_directory(snapshot_directory=directory)
                        status, output, errors = self._run_verifier(root)
                        self.assertEqual(status, 1)
                        self.assertTrue(errors)
                        self.assertNotIn("Verified", output)

    def test_1_1_verifier_checks_xml_identity_structure_and_actual_counts(self) -> None:
        for identifier_type in ("uid", "accession.version"):
            for change in ("identity", "structure", "record_count", "sequence_without_hash"):
                with (
                    self.subTest(domain=identifier_type, change=change),
                    tempfile.TemporaryDirectory() as name,
                ):
                    root = Path(name)
                    directory, manifest = self._write_integrity_snapshot(root, identifier_type, True)
                    xml_file = directory / "protein_records.xml"
                    xml_bytes = xml_file.read_bytes()
                    if change == "identity":
                        xml_bytes = xml_bytes.replace(b"1002", b"1009").replace(b"XP_123456789.2", b"XP_123456789.9")
                    elif change == "structure":
                        xml_bytes = xml_bytes.replace(b"GBSet", b"UnexpectedRoot")
                    elif change == "record_count":
                        first_record_end = xml_bytes.index(b"</GBSeq>") + len(b"</GBSeq>")
                        xml_bytes = b"<GBSet>" + xml_bytes[first_record_end:]
                    else:
                        xml_bytes = xml_bytes.replace(b"AAAA", b"TTTT")
                    xml_file.write_bytes(xml_bytes)
                    manifest["xml_file_sha256"] = hashlib.sha256(xml_bytes).hexdigest()
                    if change == "sequence_without_hash":
                        manifest.pop("xml_file_sha256")
                    self._write_manifest(directory, manifest)
                    with self.assertRaises(RuntimeError):
                        load_xml_snapshot_by_directory(snapshot_directory=directory)
                    self.assertEqual(self._run_verifier(root)[0], 1)

    def test_1_1_loader_and_verifier_reject_conflicting_xml_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            directory, manifest = self._write_integrity_snapshot(root, "accession.version", True)
            xml_file = directory / "protein_records.xml"
            xml_bytes = xml_file.read_bytes().replace(
                b"</GBSeq_accession-version>",
                b"</GBSeq_accession-version><GBSeqid>ref|WP_999999999.1|</GBSeqid>",
                1,
            )
            xml_file.write_bytes(xml_bytes)
            manifest["xml_file_sha256"] = hashlib.sha256(xml_bytes).hexdigest()
            self._write_manifest(directory, manifest)

            with self.assertRaisesRegex(RuntimeError, "conflicting accession.version candidates"):
                load_xml_snapshot_by_directory(snapshot_directory=directory)
            status, _, errors = self._run_verifier(root)
            self.assertEqual(status, 1)
            self.assertIn("conflicting accession.version candidates", errors)

    def test_unsupported_versions_do_not_enable_legacy_integrity_rules(self) -> None:
        for identifier_type in ("uid", "accession.version"):
            for include_xml in (False, True):
                with (
                    self.subTest(domain=identifier_type, xml=include_xml),
                    tempfile.TemporaryDirectory() as name,
                ):
                    root = Path(name)
                    directory, manifest = self._write_integrity_snapshot(root, identifier_type, include_xml)
                    del manifest["protein_identifiers_sha256"]
                    del manifest["normalized_protein_identifier_count"]
                    if include_xml:
                        del manifest["consolidated_record_count"]
                    load = load_xml_snapshot_by_directory if include_xml else load_snapshot_by_directory
                    for version in (None, "1.1 ", "1.2", False, 1.1):
                        with self.subTest(version=version):
                            manifest["snapshot_format_version"] = version
                            self._write_manifest(directory, manifest)
                            with self.assertRaisesRegex(ValueError, "Unsupported snapshot_format_version"):
                                load(snapshot_directory=directory)
                            status, _, errors = self._run_verifier(root)
                            self.assertEqual(status, 1)
                            self.assertIn("Unsupported snapshot_format_version", errors)

    def test_complete_1_1_and_legacy_snapshots_pass_reuse_and_verification(self) -> None:
        for identifier_type in ("uid", "accession.version"):
            for include_xml in (False, True):
                for version in (None, "1.0", "1.1"):
                    with (
                        self.subTest(domain=identifier_type, xml=include_xml, version=version),
                        tempfile.TemporaryDirectory() as name,
                    ):
                        root = Path(name)
                        directory, manifest = self._write_integrity_snapshot(root, identifier_type, include_xml)
                        if version != "1.1":
                            manifest.pop("protein_identifiers_sha256")
                            manifest.pop("normalized_protein_identifier_count")
                            manifest.pop("protein_identifier_count")
                            manifest.pop("protein_identifiers_file_name")
                            if identifier_type == "uid":
                                manifest.pop("identifier_type")
                        if version is None:
                            manifest.pop("snapshot_format_version")
                        else:
                            manifest["snapshot_format_version"] = version
                        self._write_manifest(directory, manifest)
                        original_manifest = (directory / "manifest.json").read_bytes()
                        load = load_xml_snapshot_by_directory if include_xml else load_snapshot_by_directory
                        payload = load(snapshot_directory=directory)
                        self.assertEqual(len(payload["protein_uids"]), 3)
                        self.assertEqual(self._run_verifier(root)[0], 0)
                        self.assertEqual((directory / "manifest.json").read_bytes(), original_manifest)

    def test_legacy_xml_keeps_optional_counts_and_uid_default(self) -> None:
        for version in (None, "1.0"):
            with self.subTest(version=version), tempfile.TemporaryDirectory() as name:
                root = Path(name)
                directory, manifest = self._write_integrity_snapshot(root, "uid", True)
                for field in (
                    "identifier_type", "protein_identifiers_sha256",
                    "normalized_protein_identifier_count", "protein_identifier_count",
                    "normalized_protein_uid_count", "consolidated_record_count",
                ):
                    manifest.pop(field)
                if version is None:
                    manifest.pop("snapshot_format_version")
                else:
                    manifest["snapshot_format_version"] = version
                self._write_manifest(directory, manifest)
                original_manifest = (directory / "manifest.json").read_bytes()
                self.assertEqual(
                    load_xml_snapshot_by_directory(snapshot_directory=directory)["protein_uids"],
                    ["1001", "1002", "1001"],
                )
                self.assertEqual(self._run_verifier(root)[0], 0)
                self.assertEqual((directory / "manifest.json").read_bytes(), original_manifest)

    def test_1_1_lfs_verification_requires_complete_integrity_metadata(self) -> None:
        for identifier_type in ("uid", "accession.version"):
            with self.subTest(domain=identifier_type), tempfile.TemporaryDirectory() as name:
                root = Path(name)
                directory, manifest = self._write_integrity_snapshot(root, identifier_type, True)
                xml_file = directory / "protein_records.xml"
                xml_file.write_text(
                    "version https://git-lfs.github.com/spec/v1\n"
                    f"oid sha256:{manifest['xml_file_sha256']}\nsize {xml_file.stat().st_size}\n",
                    encoding="utf-8",
                )
                self.assertEqual(self._run_verifier(root)[0], 0)
                for field in (
                    "protein_identifiers_sha256", "normalized_protein_identifier_count",
                    "xml_file_sha256", "consolidated_record_count",
                ):
                    incomplete_manifest = dict(manifest)
                    incomplete_manifest.pop(field)
                    self._write_manifest(directory, incomplete_manifest)
                    self.assertEqual(self._run_verifier(root)[0], 1)

    def _write_manifest(self, directory: Path, manifest: dict) -> None:
        (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    def _run_verifier(self, root: Path) -> tuple[int, str, str]:
        output, errors = io.StringIO(), io.StringIO()
        with redirect_stdout(output), redirect_stderr(errors):
            status = verify_raw_data(root)
        return status, output.getvalue(), errors.getvalue()

    def test_generic_only_manifests_are_verified_and_counted_for_both_domains(self) -> None:
        for identifier_type, identifier in (("uid", "1001"), ("accession.version", "WP_012345678.1")):
            with self.subTest(identifier_type=identifier_type), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                (directory / "protein_identifiers.txt").write_text(identifier + "\n", encoding="utf-8")
                self._write_manifest(directory, {
                    "snapshot_format_version": "1.1", "artifact_type": "ncbi_protein_identifier_snapshot",
                    "identifier_type": identifier_type,
                    "normalized_protein_identifier_count": 1,
                    "protein_identifiers_file_name": "protein_identifiers.txt",
                    "protein_identifiers_sha256": sha256_of_lines(text_lines=[identifier]),
                })
                status, output, errors = self._run_verifier(directory)
                self.assertEqual(status, 0)
                self.assertIn("Verified 1 raw data identities from 1 manifests.", output)
                self.assertEqual(errors, "")

    def test_generic_hash_and_filename_override_conflicting_legacy_fields_without_double_counting(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            (directory / "protein_identifiers.txt").write_text("WP_012345678.1\n", encoding="utf-8")
            (directory / "protein_uids.txt").write_text("1001\n", encoding="utf-8")
            self._write_manifest(directory, {
                "snapshot_format_version": "1.1", "identifier_type": "accession.version",
                "protein_identifiers_file_name": "protein_identifiers.txt",
                "normalized_protein_identifier_count": 1,
                "protein_uids_file_name": "protein_uids.txt",
                "protein_identifiers_sha256": sha256_of_lines(text_lines=["WP_012345678.1"]),
                "protein_uids_sha256": "incorrect-legacy-hash",
            })
            status, output, errors = self._run_verifier(directory)
            self.assertEqual(status, 0)
            self.assertIn("Verified 1 raw data identities", output)
            self.assertEqual(errors, "")

    def test_hash_and_filename_fallbacks_are_selected_independently(self) -> None:
        for hash_field in ("protein_identifiers_sha256", "protein_uids_sha256"):
            for file_field in ("protein_identifiers_file_name", "protein_uids_file_name"):
                with self.subTest(hash_field=hash_field, file_field=file_field), tempfile.TemporaryDirectory() as name:
                    directory = Path(name)
                    (directory / "provided_ids.txt").write_text("WP_012345678.1\n", encoding="utf-8")
                    self._write_manifest(directory, {
                        "snapshot_format_version": "1.0", "identifier_type": "accession.version",
                        file_field: "provided_ids.txt",
                        hash_field: sha256_of_lines(text_lines=["WP_012345678.1"]),
                    })
                    status, output, errors = self._run_verifier(directory)
                    self.assertEqual(status, 0)
                    self.assertIn("Verified 1 raw data identities", output)
                    self.assertEqual(errors, "")

    def test_missing_or_invalid_generic_hash_is_an_error_instead_of_legacy_fallback_or_skip(self) -> None:
        for value in (None, "", " \t", 123, [], "absent"):
            with self.subTest(value=value), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                (directory / "protein_uids.txt").write_text("1001\n", encoding="utf-8")
                manifest = {
                    "snapshot_format_version": "1.1", "artifact_type": "ncbi_protein_identifier_snapshot",
                    "identifier_type": "uid", "protein_uids_file_name": "protein_uids.txt",
                    "protein_identifiers_sha256": value,
                    "protein_uids_sha256": sha256_of_lines(text_lines=["1001"]),
                }
                if value == "absent":
                    manifest.pop("protein_identifiers_sha256")
                    manifest.pop("protein_uids_sha256")
                self._write_manifest(directory, manifest)
                status, output, errors = self._run_verifier(directory)
                self.assertEqual(status, 1)
                self.assertIn("protein_identifiers_sha256", errors)
                self.assertNotIn("Verified", output)

    def test_invalid_or_missing_generic_file_is_not_replaced_by_legacy_copy(self) -> None:
        for filename in (None, "", 123, "missing_identifiers.txt"):
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                (directory / "protein_uids.txt").write_text("1001\n", encoding="utf-8")
                self._write_manifest(directory, {
                    "identifier_type": "uid",
                    "protein_identifiers_file_name": filename,
                    "protein_uids_file_name": "protein_uids.txt",
                    "protein_identifiers_sha256": sha256_of_lines(text_lines=["1001"]),
                })
                status, output, errors = self._run_verifier(directory)
                self.assertEqual(status, 1)
                self.assertTrue("Invalid identifier file name" in errors or "Missing protein identifier file" in errors)
                self.assertNotIn("Verified", output)

    def test_generic_hash_mismatch_rejects_the_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            (directory / "protein_identifiers.txt").write_text("WP_012345678.2\n", encoding="utf-8")
            self._write_manifest(directory, {
                "snapshot_format_version": "1.1", "identifier_type": "accession.version",
                "protein_identifiers_file_name": "protein_identifiers.txt",
                "normalized_protein_identifier_count": 1,
                "protein_identifiers_sha256": sha256_of_lines(text_lines=["WP_012345678.1"]),
            })
            status, output, errors = self._run_verifier(directory)
            self.assertEqual(status, 1)
            self.assertIn("identifier hash mismatch", errors)
            self.assertNotIn("Verified", output)

    def test_legacy_only_manifests_are_still_verified_and_counted(self) -> None:
        for file_field in (None, "protein_uids_file_name", "snapshot_file_name"):
            with self.subTest(file_field=file_field), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                (directory / "protein_uids.txt").write_text("1001\n", encoding="utf-8")
                manifest = {
                    "protein_uids_sha256": sha256_of_lines(text_lines=["1001"]),
                    "normalized_protein_uid_count": 1,
                }
                if file_field:
                    manifest[file_field] = "protein_uids.txt"
                self._write_manifest(directory, manifest)
                status, output, errors = self._run_verifier(directory)
                self.assertEqual(status, 0)
                self.assertIn("Verified 1 raw data identities", output)
                self.assertEqual(errors, "")

    def test_generic_xml_snapshots_count_identifier_and_xml_checks_including_lfs(self) -> None:
        xml_bytes = b"<GBSet><GBSeq><GBSeq_accession-version>WP_012345678.1</GBSeq_accession-version></GBSeq></GBSet>"
        xml_sha256 = hashlib.sha256(xml_bytes).hexdigest()
        for lfs_pointer in (False, True):
            with self.subTest(lfs_pointer=lfs_pointer), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                (directory / "protein_identifiers.txt").write_text("WP_012345678.1\n", encoding="utf-8")
                if lfs_pointer:
                    (directory / "protein_records.xml").write_text(
                        "version https://git-lfs.github.com/spec/v1\n"
                        f"oid sha256:{xml_sha256}\nsize {len(xml_bytes)}\n", encoding="utf-8",
                    )
                else:
                    (directory / "protein_records.xml").write_bytes(xml_bytes)
                self._write_manifest(directory, {
                    "snapshot_format_version": "1.1", "artifact_type": "ncbi_protein_xml_snapshot",
                    "source_type": "identifier_snapshot", "identifier_type": "accession.version",
                    "normalized_protein_identifier_count": 1, "consolidated_record_count": 1,
                    "protein_identifiers_file_name": "protein_identifiers.txt",
                    "protein_identifiers_sha256": sha256_of_lines(text_lines=["WP_012345678.1"]),
                    "xml_file_name": "protein_records.xml", "xml_file_sha256": xml_sha256,
                })
                status, output, errors = self._run_verifier(directory)
                self.assertEqual(status, 0)
                self.assertIn("Verified 2 raw data identities from 1 manifests.", output)
                self.assertEqual(errors, "")

    def test_partial_xml_metadata_cannot_be_silently_skipped(self) -> None:
        for manifest in (
            {"xml_file_sha256": "a" * 64},
            {"xml_file_name": "protein_records.xml"},
            {"xml_file_name": "protein_records.xml", "xml_file_sha256": None},
            {"xml_file_name": "", "xml_file_sha256": "a" * 64},
        ):
            with self.subTest(manifest=manifest), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                self._write_manifest(directory, manifest)
                status, output, errors = self._run_verifier(directory)
                self.assertEqual(status, 1)
                self.assertIn("Invalid XML", errors)
                self.assertNotIn("Verified", output)

    def test_mixed_manifest_collection_counts_generic_and_legacy_identities(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            for label, manifest, filename in (
                ("legacy", {"protein_uids_sha256": sha256_of_lines(text_lines=["1001"])}, "protein_uids.txt"),
                ("generic", {"snapshot_format_version": "1.1", "identifier_type": "uid",
                             "normalized_protein_identifier_count": 1,
                             "protein_identifiers_sha256": sha256_of_lines(text_lines=["1001"])}, "protein_identifiers.txt"),
                ("xml_only", {"xml_file_name": "protein_records.xml", "xml_file_sha256": hashlib.sha256(b"<GBSet />").hexdigest()}, None),
            ):
                directory = root / label
                directory.mkdir()
                if filename:
                    (directory / filename).write_text("1001\n", encoding="utf-8")
                else:
                    (directory / "protein_records.xml").write_bytes(b"<GBSet />")
                self._write_manifest(directory, manifest)
            status, output, errors = self._run_verifier(root)
            self.assertEqual(status, 0)
            self.assertIn("Verified 3 raw data identities from 3 manifests.", output)
            self.assertEqual(errors, "")

    def test_minimal_legacy_manifest_is_verified_using_default_uid_file(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            path = directory / "protein_uids.txt"
            path.write_text("1001\n", encoding="utf-8")
            manifest = {
                "protein_uids_sha256": sha256_of_lines(text_lines=["1001"]),
                "normalized_protein_uid_count": 1,
            }
            self.assertEqual(_verify_uid_file(directory, manifest), [])
            path.write_text("1002\n", encoding="utf-8")
            errors = _verify_uid_file(directory, manifest)
            self.assertEqual(len(errors), 1)
            self.assertIn("identifier hash mismatch", errors[0])

    def test_default_uid_is_restricted_to_legacy_manifest_versions(self) -> None:
        for version in (None, "1.0", "1.1", "2.0"):
            with self.subTest(version=version), tempfile.TemporaryDirectory() as name:
                directory = Path(name)
                (directory / "protein_uids.txt").write_text("1001\n", encoding="utf-8")
                manifest = {
                    "protein_uids_sha256": sha256_of_lines(text_lines=["1001"]),
                    "protein_uids_file_name": "protein_uids.txt",
                }
                if version is not None:
                    manifest["snapshot_format_version"] = version
                errors = _verify_uid_file(directory, manifest)
                if version in (None, "1.0"):
                    self.assertEqual(errors, [])
                else:
                    self.assertEqual(len(errors), 1)
                    expected = "identifier_type is required" if version == "1.1" else "Unsupported snapshot_format_version"
                    self.assertIn(expected, errors[0])

    def test_explicit_null_type_in_legacy_manifest_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            (directory / "protein_uids.txt").write_text("1001\n", encoding="utf-8")
            manifest = {
                "snapshot_format_version": "1.0",
                "identifier_type": None,
                "protein_uids_sha256": sha256_of_lines(text_lines=["1001"]),
            }
            errors = _verify_uid_file(directory, manifest)
            self.assertEqual(len(errors), 1)
            self.assertIn("identifier_type must not be null", errors[0])

    def test_uid_file_without_identifier_type_defaults_to_uid(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_directory = Path(temporary_directory_name)
            uid_file_path = snapshot_directory / "protein_uids.txt"
            protein_uids = ["1001", "1002"]
            uid_file_path.write_text("\n".join(protein_uids) + "\n", encoding="utf-8")
            manifest = {
                "snapshot_file_name": uid_file_path.name,
                "protein_uids_sha256": sha256_of_lines(
                    text_lines=protein_uids,
                    deduplicate_lines_preserving_order=False,
                    sort_lines=False,
                ),
            }

            self.assertEqual(_verify_uid_file(snapshot_directory, manifest), [])

    def test_identifier_file_uses_generic_manifest_fields(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_directory = Path(temporary_directory_name)
            identifier_file_path = snapshot_directory / "protein_identifiers.txt"
            protein_identifiers = ["WP_012345678.1"]
            identifier_file_path.write_text("WP_012345678.1\n", encoding="utf-8")
            manifest = {
                "protein_identifiers_file_name": identifier_file_path.name,
                "protein_identifiers_sha256": sha256_of_lines(
                    text_lines=protein_identifiers,
                    deduplicate_lines_preserving_order=False,
                    sort_lines=False,
                ),
                "identifier_type": "accession.version",
            }

            self.assertEqual(_verify_uid_file(snapshot_directory, manifest), [])

    def test_uid_file_rejects_unknown_identifier_type(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_directory = Path(temporary_directory_name)
            uid_file_path = snapshot_directory / "protein_uids.txt"
            protein_uids = ["1001"]
            uid_file_path.write_text("1001\n", encoding="utf-8")
            manifest = {
                "snapshot_file_name": uid_file_path.name,
                "protein_uids_sha256": sha256_of_lines(
                    text_lines=protein_uids,
                    deduplicate_lines_preserving_order=False,
                    sort_lines=False,
                ),
                "identifier_type": "accession",
            }

            errors = _verify_uid_file(snapshot_directory, manifest)

            self.assertEqual(len(errors), 1)
            self.assertIn("Invalid NCBI protein identifier_type", errors[0])

    def test_uid_file_rejects_mixed_identifier_list(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_directory = Path(temporary_directory_name)
            uid_file_path = snapshot_directory / "protein_uids.txt"
            protein_uids = ["1001", "WP_012345678.1"]
            uid_file_path.write_text("\n".join(protein_uids) + "\n", encoding="utf-8")
            manifest = {
                "snapshot_file_name": uid_file_path.name,
                "protein_uids_sha256": sha256_of_lines(
                    text_lines=protein_uids,
                    deduplicate_lines_preserving_order=False,
                    sort_lines=False,
                ),
                "identifier_type": "uid",
            }

            errors = _verify_uid_file(snapshot_directory, manifest)

            self.assertEqual(len(errors), 1)
            self.assertIn("Mixed or invalid identifier lists", errors[0])

    def test_materialized_xml_is_verified_from_file_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_directory = Path(temporary_directory_name)
            xml_file_path = snapshot_directory / "protein_records.xml"
            xml_bytes = b"<GBSet><GBSeq /></GBSet>\n"
            xml_file_path.write_bytes(xml_bytes)

            expected_sha256 = hashlib.sha256(xml_bytes).hexdigest()
            manifest = {
                "xml_file_name": xml_file_path.name,
                "xml_file_sha256": expected_sha256,
            }

            self.assertEqual(_verify_xml_file(snapshot_directory, manifest), [])

    def test_matching_lfs_pointer_is_verified_without_materialized_xml(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_directory = Path(temporary_directory_name)
            xml_file_path = snapshot_directory / "protein_records.xml"
            expected_sha256 = (
                "ada0563932d20f54f68c614df3c40d4fa3520038ad438a41cc471a1eacb970bf"
            )
            xml_file_path.write_text(
                "version https://git-lfs.github.com/spec/v1\n"
                f"oid sha256:{expected_sha256}\n"
                "size 363034184\n",
                encoding="utf-8",
            )
            manifest = {
                "xml_file_name": xml_file_path.name,
                "xml_file_sha256": expected_sha256,
            }

            self.assertEqual(_verify_xml_file(snapshot_directory, manifest), [])

    def test_mismatching_lfs_pointer_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            snapshot_directory = Path(temporary_directory_name)
            xml_file_path = snapshot_directory / "protein_records.xml"
            expected_sha256 = "a" * 64
            actual_sha256 = "b" * 64
            xml_file_path.write_text(
                "version https://git-lfs.github.com/spec/v1\n"
                f"oid sha256:{actual_sha256}\n"
                "size 363034184\n",
                encoding="utf-8",
            )
            manifest = {
                "xml_file_name": xml_file_path.name,
                "xml_file_sha256": expected_sha256,
            }

            errors = _verify_xml_file(snapshot_directory, manifest)

            self.assertEqual(len(errors), 1)
            self.assertIn("XML LFS pointer SHA-256 mismatch", errors[0])
            self.assertIn(expected_sha256, errors[0])
            self.assertIn(actual_sha256, errors[0])


if __name__ == "__main__":
    unittest.main()
