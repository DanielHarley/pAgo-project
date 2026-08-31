from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from src.pago_pipeline.pfam_hmm_bundle import (
    build_pfam_hmm_bundle,
    load_pfam_hmm_bundle_registry,
    validate_pfam_hmm_reference_data,
)
from src.pago_pipeline.storage import sha256_of_file
from tests.pfam_hmm_test_support import (
    DEFAULT_TEST_MODELS,
    build_fake_pyhmmer,
    build_test_hmm_payload,
    read_registry,
    write_registry,
    write_test_reference_bundle,
)


class PfamHmmBundleTests(unittest.TestCase):
    def test_validation_and_bundle_build_preserve_canonical_order_deterministically(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            registry_path, hmm_directory = write_test_reference_bundle(
                reference_directory=temporary_directory / "reference"
            )
            fake_pyhmmer = build_fake_pyhmmer()

            validation = validate_pfam_hmm_reference_data(
                registry_file_path=registry_path,
                hmm_directory=hmm_directory,
                pyhmmer_module=fake_pyhmmer,
            )
            self.assertEqual(validation.registry.source_release, "test-release-1")
            self.assertEqual(validation.registry.model_count, 2)
            self.assertEqual(
                [model.inspected_model.accession for model in validation.validated_models],
                ["PF02171.22", "PF16487.10"],
            )

            first_result = build_pfam_hmm_bundle(
                output_directory=temporary_directory / "first",
                registry_file_path=registry_path,
                hmm_directory=hmm_directory,
                pyhmmer_module=fake_pyhmmer,
            )
            second_result = build_pfam_hmm_bundle(
                output_directory=temporary_directory / "second",
                registry_file_path=registry_path,
                hmm_directory=hmm_directory,
                pyhmmer_module=fake_pyhmmer,
            )

            self.assertEqual(
                first_result.bundle_file_path.read_bytes(),
                second_result.bundle_file_path.read_bytes(),
            )
            bundle_text = first_result.bundle_file_path.read_text(encoding="utf-8")
            self.assertLess(bundle_text.index("NAME  Piwi"), bundle_text.index("NAME  ArgoMid"))
            self.assertEqual(len(first_result.pressed_file_paths), 4)
            self.assertTrue(
                all(path.stat().st_size > 0 for path in first_result.pressed_file_paths)
            )
            self.assertEqual(
                [path.read_bytes() for path in first_result.pressed_file_paths],
                [path.read_bytes() for path in second_result.pressed_file_paths],
            )

    def test_validation_rejects_inventory_and_source_hash_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            registry_path, hmm_directory = write_test_reference_bundle(
                reference_directory=temporary_directory / "reference"
            )
            fake_pyhmmer = build_fake_pyhmmer()

            unexpected_file_path = hmm_directory / "Unexpected__PF99999.hmm"
            unexpected_file_path.write_bytes(
                build_test_hmm_payload(name="Unexpected", accession="PF99999.1")
            )
            with self.assertRaisesRegex(RuntimeError, "inventory"):
                validate_pfam_hmm_reference_data(
                    registry_file_path=registry_path,
                    hmm_directory=hmm_directory,
                    pyhmmer_module=fake_pyhmmer,
                )
            unexpected_file_path.unlink()

            source_file_path = hmm_directory / "Piwi__PF02171.hmm"
            temporarily_missing_file_path = hmm_directory / "Piwi__PF02171.saved"
            source_file_path.replace(temporarily_missing_file_path)
            with self.assertRaisesRegex(RuntimeError, "inventory"):
                validate_pfam_hmm_reference_data(
                    registry_file_path=registry_path,
                    hmm_directory=hmm_directory,
                    pyhmmer_module=fake_pyhmmer,
                )
            temporarily_missing_file_path.replace(source_file_path)

            source_file_path.write_bytes(source_file_path.read_bytes() + b"\n")
            with self.assertRaisesRegex(RuntimeError, "SHA-256 mismatch"):
                validate_pfam_hmm_reference_data(
                    registry_file_path=registry_path,
                    hmm_directory=hmm_directory,
                    pyhmmer_module=fake_pyhmmer,
                )

    def test_registry_rejects_noncanonical_order_and_missing_release(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            registry_path, _ = write_test_reference_bundle(
                reference_directory=temporary_directory / "reference"
            )
            registry = read_registry(registry_path)
            registry["models"] = list(reversed(registry["models"]))
            write_registry(registry_path, registry)
            with self.assertRaisesRegex(RuntimeError, "canonical order"):
                load_pfam_hmm_bundle_registry(registry_file_path=registry_path)

            registry["models"] = list(reversed(registry["models"]))
            registry.pop("source_release")
            write_registry(registry_path, registry)
            with self.assertRaisesRegex(RuntimeError, "source_release"):
                load_pfam_hmm_bundle_registry(registry_file_path=registry_path)

    def test_registry_requires_version_pinned_pfam_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            registry_path, _ = write_test_reference_bundle(
                reference_directory=temporary_directory / "reference"
            )

            registry = read_registry(registry_path)
            self.assertEqual(registry["models"][0]["accession"], "PF02171.22")
            self.assertIsNotNone(
                load_pfam_hmm_bundle_registry(registry_file_path=registry_path)
            )

            registry["models"][0]["accession"] = "PF02171"
            write_registry(registry_path, registry)
            with self.assertRaisesRegex(RuntimeError, "version-pinned Pfam accession"):
                load_pfam_hmm_bundle_registry(registry_file_path=registry_path)

    def test_validation_rejects_model_identity_count_alphabet_and_missing_ga(self) -> None:
        cases = (
            ("name", "WrongName", "name mismatch"),
            ("alphabet", "dna", "alphabet mismatch"),
            ("missing_ga", None, "missing required gathering cutoffs"),
            ("second_model", None, "model count mismatch"),
        )
        for case_name, replacement, expected_error in cases:
            with self.subTest(case_name=case_name):
                with tempfile.TemporaryDirectory() as temporary_directory_name:
                    temporary_directory = Path(temporary_directory_name)
                    registry_path, hmm_directory = write_test_reference_bundle(
                        reference_directory=temporary_directory / "reference"
                    )
                    registry = read_registry(registry_path)
                    first_model = registry["models"][0]
                    source_file_path = hmm_directory / first_model["file_name"]

                    if case_name == "name":
                        first_model[case_name] = replacement
                    elif case_name == "alphabet":
                        source_file_path.write_bytes(
                            build_test_hmm_payload(
                                name="Piwi",
                                accession="PF02171.22",
                                alphabet="dna",
                            )
                        )
                        first_model["sha256"] = sha256_of_file(
                            input_file_path=source_file_path
                        )
                    elif case_name == "missing_ga":
                        source_file_path.write_bytes(
                            build_test_hmm_payload(
                                name="Piwi",
                                accession="PF02171.22",
                                gathering=None,
                            )
                        )
                        first_model["sha256"] = sha256_of_file(
                            input_file_path=source_file_path
                        )
                    else:
                        source_file_path.write_bytes(
                            source_file_path.read_bytes()
                            + build_test_hmm_payload(
                                name="Extra",
                                accession="PF99998.1",
                            )
                        )
                        first_model["sha256"] = sha256_of_file(
                            input_file_path=source_file_path
                        )
                    write_registry(registry_path, registry)

                    with self.assertRaisesRegex(RuntimeError, expected_error):
                        validate_pfam_hmm_reference_data(
                            registry_file_path=registry_path,
                            hmm_directory=hmm_directory,
                            pyhmmer_module=build_fake_pyhmmer(),
                        )

    def test_validation_only_api_does_not_require_hmmpress_or_write_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            registry_path, hmm_directory = write_test_reference_bundle(
                reference_directory=temporary_directory / "reference",
                models=DEFAULT_TEST_MODELS,
            )
            validation = validate_pfam_hmm_reference_data(
                registry_file_path=registry_path,
                hmm_directory=hmm_directory,
                pyhmmer_module=build_fake_pyhmmer(hmmpress_function=None),
            )
            self.assertEqual(validation.registry.model_count, 2)
            build_result = build_pfam_hmm_bundle(
                output_directory=temporary_directory / "unpressed_fixture_output",
                registry_file_path=registry_path,
                hmm_directory=hmm_directory,
                pyhmmer_module=build_fake_pyhmmer(hmmpress_function=None),
                press_bundle=False,
            )
            self.assertTrue(build_result.bundle_file_path.is_file())
            self.assertEqual(build_result.pressed_file_paths, ())
            self.assertEqual(
                sorted(path.name for path in hmm_directory.iterdir()),
                [
                    "ArgoMid__PF16487.hmm",
                    "Piwi__PF02171.hmm",
                    "pfam_hmm_bundle_lock.json",
                ],
            )

    def test_build_refuses_to_overwrite_any_canonical_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            registry_path, hmm_directory = write_test_reference_bundle(
                reference_directory=temporary_directory / "reference"
            )
            output_directory = temporary_directory / "output"
            output_directory.mkdir()
            preexisting_pressed_path = output_directory / "pago_domains.hmm.h3i"
            preexisting_pressed_path.write_bytes(b"preserve-me")

            with self.assertRaisesRegex(FileExistsError, "Refusing to overwrite"):
                build_pfam_hmm_bundle(
                    output_directory=output_directory,
                    registry_file_path=registry_path,
                    hmm_directory=hmm_directory,
                    pyhmmer_module=build_fake_pyhmmer(),
                )
            self.assertEqual(preexisting_pressed_path.read_bytes(), b"preserve-me")
            self.assertEqual(
                sorted(path.name for path in output_directory.iterdir()),
                ["pago_domains.hmm.h3i"],
            )


if __name__ == "__main__":
    unittest.main()
