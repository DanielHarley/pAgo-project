from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from src.pago_pipeline.apaz_hmm_build import (
    APAZ_MODEL_SPECS,
    _verify_hmm_artifact_on_disk,
    build_apaz_hmms,
    validate_apaz_hmm_build_inputs,
)
from src.pago_pipeline.storage import write_json_atomic
from tests.apaz_hmm_test_support import (
    build_real_amino_hmm_for_tests,
    deterministic_fake_hmm_builder,
    read_json,
    refresh_apaz_locked_seed,
    refresh_apaz_partitions_lock_hash,
    write_apaz_reference_fixture,
    write_stockholm,
    write_valid_amino_hmm,
)


class ApazHmmBuildTests(unittest.TestCase):
    def test_apaz_build_validates_global_union_and_builds_six_real_hmms(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            reference_directory = temporary_directory / "apaz_reference"
            output_directory = temporary_directory / "apaz_hmms"
            write_apaz_reference_fixture(reference_directory=reference_directory)

            validation = validate_apaz_hmm_build_inputs(
                reference_directory=reference_directory
            )
            self.assertEqual(len(validation.seed_artifacts), 6)
            self.assertEqual(
                validation.seed_artifacts[0].alignment_summary.sequence_count, 5
            )
            self.assertEqual(
                set(validation.partition_by_accession.values()),
                {"BUILD", "CALIBRATION", "FINAL_HOLDOUT"},
            )

            result = build_apaz_hmms(
                reference_directory=reference_directory,
                output_directory=output_directory,
                random_seed=42,
            )
            repeated_result = build_apaz_hmms(
                reference_directory=reference_directory,
                output_directory=temporary_directory / "apaz_hmms_repeated",
                random_seed=42,
            )
            self.assertEqual(len(result.hmm_artifacts), len(APAZ_MODEL_SPECS))
            self.assertEqual(result.builder_identity, "pyhmmer/0.12.3")
            self.assertTrue(
                all(
                    artifact.hmm_file_path.stat().st_size > 0
                    and artifact.match_state_count >= 1
                    for artifact in result.hmm_artifacts
                )
            )
            self.assertEqual(
                [artifact.hmm_sha256 for artifact in result.hmm_artifacts],
                [artifact.hmm_sha256 for artifact in repeated_result.hmm_artifacts],
            )

    def test_apaz_real_hmms_carry_no_volatile_provenance_lines(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            reference_directory = directory / "apaz_reference"
            write_apaz_reference_fixture(reference_directory=reference_directory)
            first = build_apaz_hmms(
                reference_directory=reference_directory,
                output_directory=directory / "a",
                random_seed=42,
            )
            second = build_apaz_hmms(
                reference_directory=reference_directory,
                output_directory=directory / "b",
                random_seed=42,
            )
            for artifact in first.hmm_artifacts:
                payload = artifact.hmm_file_path.read_bytes()
                self.assertNotIn(b"\nCOM ", payload)
                self.assertNotIn(b"\nDATE ", payload)
            self.assertEqual(
                [a.hmm_sha256 for a in first.hmm_artifacts],
                [a.hmm_sha256 for a in second.hmm_artifacts],
            )

    def test_apaz_rejects_global_subgroup_drift_and_subgroup_overlap(self) -> None:
        for case_name in ("global_missing", "subgroup_overlap"):
            with self.subTest(case_name=case_name):
                with tempfile.TemporaryDirectory() as temporary_directory_name:
                    reference_directory = (
                        Path(temporary_directory_name) / "apaz_reference"
                    )
                    write_apaz_reference_fixture(
                        reference_directory=reference_directory
                    )
                    lock_file_path = reference_directory / "seeds_lock.json"
                    if case_name == "global_missing":
                        seed_file_path = reference_directory / "apaz_global.sto"
                        write_stockholm(
                            file_path=seed_file_path,
                            sequence_by_accession={
                                f"APAZ{index}.1": "ACDEFG"
                                for index in range(1, 5)
                            },
                        )
                        refresh_apaz_locked_seed(
                            lock_file_path=lock_file_path,
                            seed_file_path=seed_file_path,
                            sequence_count=4,
                        )
                        expected_error = "exact union"
                    else:
                        seed_file_path = reference_directory / "apaz_Ib.sto"
                        write_stockholm(
                            file_path=seed_file_path,
                            sequence_by_accession={"APAZ1.1": "ACDEFG"},
                        )
                        refresh_apaz_locked_seed(
                            lock_file_path=lock_file_path,
                            seed_file_path=seed_file_path,
                        )
                        expected_error = "accession-disjoint"

                    with self.assertRaisesRegex(RuntimeError, expected_error):
                        validate_apaz_hmm_build_inputs(
                            reference_directory=reference_directory
                        )

    def test_apaz_builder_failure_removes_every_partial_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory_name:
            temporary_directory = Path(temporary_directory_name)
            reference_directory = temporary_directory / "apaz_reference"
            output_directory = temporary_directory / "output"
            write_apaz_reference_fixture(reference_directory=reference_directory)
            calls = 0

            def failing_builder(**kwargs) -> None:
                nonlocal calls
                calls += 1
                deterministic_fake_hmm_builder(**kwargs)
                if calls == 3:
                    raise RuntimeError("simulated HMM build failure")

            with self.assertRaisesRegex(RuntimeError, "simulated HMM build failure"):
                build_apaz_hmms(
                    reference_directory=reference_directory,
                    output_directory=output_directory,
                    hmm_builder=failing_builder,
                    builder_identity="fixture-builder/1",
                )
            self.assertEqual(list(output_directory.iterdir()), [])

    def test_apaz_build_rejects_nondeterministic_or_coerced_random_seeds(
        self,
    ) -> None:
        for random_seed in (0, 42.5, True):
            with self.subTest(random_seed=random_seed):
                with tempfile.TemporaryDirectory() as temporary_directory_name:
                    reference_directory = (
                        Path(temporary_directory_name) / "apaz_reference"
                    )
                    write_apaz_reference_fixture(
                        reference_directory=reference_directory
                    )
                    with self.assertRaisesRegex(ValueError, "random_seed"):
                        build_apaz_hmms(
                            reference_directory=reference_directory,
                            output_directory=(
                                Path(temporary_directory_name) / "output"
                            ),
                            hmm_builder=deterministic_fake_hmm_builder,
                            builder_identity="fixture-builder/1",
                            random_seed=random_seed,
                        )


class ApazHmmStructuralValidationTests(unittest.TestCase):
    """B3.2: every written HMM is reopened from disk and structurally proven."""

    def test_accepts_a_single_well_formed_amino_model(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            hmm_path = Path(name) / "good.hmm"
            write_valid_amino_hmm(output_path=hmm_path, model_name="APAZ_IA")
            report = _verify_hmm_artifact_on_disk(
                hmm_file_path=hmm_path, expected_model_name="APAZ_IA"
            )
            self.assertGreaterEqual(report.match_state_count, 1)
            self.assertEqual(report.file_size_bytes, hmm_path.stat().st_size)

    def test_rejects_a_corrupt_hmm_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            hmm_path = Path(name) / "corrupt.hmm"
            write_valid_amino_hmm(output_path=hmm_path, model_name="APAZ_IA")
            payload = hmm_path.read_bytes()
            hmm_path.write_bytes(payload[:20] + b"\nGARBLED BODY\n" + payload[40:])
            with self.assertRaisesRegex(RuntimeError, "parse"):
                _verify_hmm_artifact_on_disk(
                    hmm_file_path=hmm_path, expected_model_name="APAZ_IA"
                )

    def test_rejects_an_empty_hmm_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            hmm_path = Path(name) / "empty.hmm"
            hmm_path.write_bytes(b"")
            with self.assertRaisesRegex(RuntimeError, "empty"):
                _verify_hmm_artifact_on_disk(
                    hmm_file_path=hmm_path, expected_model_name="APAZ_IA"
                )

    def test_rejects_more_than_one_model_in_one_file(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            hmm_path = Path(name) / "multi.hmm"
            write_valid_amino_hmm(output_path=hmm_path, model_name="APAZ_IA")
            payload = hmm_path.read_bytes()
            hmm_path.write_bytes(payload + payload)
            with self.assertRaisesRegex(RuntimeError, "more than one model"):
                _verify_hmm_artifact_on_disk(
                    hmm_file_path=hmm_path, expected_model_name="APAZ_IA"
                )

    def test_rejects_an_unexpected_model_name(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            hmm_path = Path(name) / "named.hmm"
            write_valid_amino_hmm(output_path=hmm_path, model_name="APAZ_IB")
            with self.assertRaisesRegex(RuntimeError, "name mismatch"):
                _verify_hmm_artifact_on_disk(
                    hmm_file_path=hmm_path, expected_model_name="APAZ_IA"
                )

    def test_rejects_a_non_amino_alphabet(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            seed_path = directory / "dna.sto"
            write_stockholm(
                file_path=seed_path,
                sequence_by_accession={
                    "AAA.1": "ACGTACGTACGT",
                    "BBB.1": "ACGTACGAACGT",
                    "CCC.1": "ACGTTCGTACGT",
                },
            )
            hmm_path = directory / "dna.hmm"
            build_real_amino_hmm_for_tests(
                seed_alignment_file_path=seed_path,
                output_hmm_file_path=hmm_path,
                model_name="APAZ_IA",
                alphabet_kind="dna",
            )
            with self.assertRaisesRegex(RuntimeError, "amino"):
                _verify_hmm_artifact_on_disk(
                    hmm_file_path=hmm_path, expected_model_name="APAZ_IA"
                )

    def test_build_aborts_and_cleans_up_when_an_artifact_is_unparseable(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            reference_directory = directory / "apaz_reference"
            output_directory = directory / "output"
            write_apaz_reference_fixture(reference_directory=reference_directory)

            def corrupting_builder(
                *, output_hmm_file_path: Path, **_: object
            ) -> None:
                output_hmm_file_path.write_bytes(b"HMMER3/f [bad]\nnot a model\n")

            with self.assertRaisesRegex(RuntimeError, "parse|model"):
                build_apaz_hmms(
                    reference_directory=reference_directory,
                    output_directory=output_directory,
                    hmm_builder=corrupting_builder,
                    builder_identity="fixture-builder/1",
                )
            self.assertEqual(list(output_directory.iterdir()), [])


class ApazHmmBuildInputRejectionTests(unittest.TestCase):
    """B3.4: the build refuses inputs that break the frozen B2 v2 contract."""

    def test_rejects_a_seed_accession_outside_the_build_partition(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            reference_directory = Path(name) / "apaz_reference"
            write_apaz_reference_fixture(reference_directory=reference_directory)
            partitions_path = reference_directory / "apaz_partitions.csv"
            with partitions_path.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            for row in rows:
                if row["accession"] == "APAZ1.1":
                    row["partition"] = "CALIBRATION"
            with partitions_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            refresh_apaz_partitions_lock_hash(reference_directory=reference_directory)
            with self.assertRaisesRegex(RuntimeError, "non-BUILD"):
                validate_apaz_hmm_build_inputs(
                    reference_directory=reference_directory
                )

    def test_rejects_a_partition_table_without_the_v2_split_group_column(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as name:
            reference_directory = Path(name) / "apaz_reference"
            write_apaz_reference_fixture(reference_directory=reference_directory)
            partitions_path = reference_directory / "apaz_partitions.csv"
            with partitions_path.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            fieldnames = [c for c in rows[0] if c != "split_group_id"]
            with partitions_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fieldnames, extrasaction="ignore"
                )
                writer.writeheader()
                writer.writerows(rows)
            refresh_apaz_partitions_lock_hash(reference_directory=reference_directory)
            with self.assertRaisesRegex(RuntimeError, "split_group_id"):
                validate_apaz_hmm_build_inputs(
                    reference_directory=reference_directory
                )

    def test_rejects_an_incompatible_seed_lock_format_version(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            reference_directory = Path(name) / "apaz_reference"
            write_apaz_reference_fixture(reference_directory=reference_directory)
            lock_path = reference_directory / "seeds_lock.json"
            payload = read_json(lock_path)
            payload["lock_format_version"] = "1.0"
            write_json_atomic(payload=payload, output_file_path=lock_path)
            with self.assertRaisesRegex(RuntimeError, "lock_format_version"):
                validate_apaz_hmm_build_inputs(
                    reference_directory=reference_directory
                )

    def test_rejects_a_seed_whose_bytes_diverge_from_the_lock_hash(self) -> None:
        with tempfile.TemporaryDirectory() as name:
            reference_directory = Path(name) / "apaz_reference"
            write_apaz_reference_fixture(reference_directory=reference_directory)
            write_stockholm(
                file_path=reference_directory / "apaz_Ia.sto",
                sequence_by_accession={"APAZ1.1": "GGGGGG"},
            )
            with self.assertRaisesRegex(RuntimeError, "hash mismatch"):
                validate_apaz_hmm_build_inputs(
                    reference_directory=reference_directory
                )


if __name__ == "__main__":
    unittest.main()
