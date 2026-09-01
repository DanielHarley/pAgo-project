"""Test support for the APAZ HMM builder (Phase B, sub-step B3).

Everything here depends only on the APAZ reference layer (``apaz_hmm_build``)
and the shared storage helpers. It deliberately imports no clade, MID-PIWI,
or downstream (B4+) module, so the APAZ HMM tests run from a checkout that
contains only A + B0 + B1 + B2 + B3.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

from src.pago_pipeline.apaz_hmm_build import APAZ_MODEL_SPECS
from src.pago_pipeline.storage import (
    sha256_of_file,
    write_bytes_atomic,
    write_json_atomic,
)


def sequence_sha256(token: str) -> str:
    return hashlib.sha256(token.encode("utf-8")).hexdigest()


def read_json(file_path: Path) -> dict[str, object]:
    return json.loads(Path(file_path).read_text(encoding="utf-8"))


def write_stockholm(
    *,
    file_path: Path,
    sequence_by_accession: dict[str, str],
) -> None:
    lines = ["# STOCKHOLM 1.0"]
    lines.extend(
        f"{accession} {sequence}"
        for accession, sequence in sequence_by_accession.items()
    )
    lines.append("//")
    file_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_real_amino_hmm_for_tests(
    *,
    seed_alignment_file_path: Path,
    output_hmm_file_path: Path,
    model_name: str,
    random_seed: int = 42,
    alphabet_kind: str = "amino",
) -> None:
    """Construct one genuine, structurally valid HMM for tests.

    Post-build validation in ``apaz_hmm_build`` reopens every artifact from
    disk, so test doubles must emit real HMMER3 profiles rather than
    placeholders.
    """
    import pyhmmer

    alphabet = (
        pyhmmer.easel.Alphabet.amino()
        if alphabet_kind == "amino"
        else pyhmmer.easel.Alphabet.dna()
    )
    with pyhmmer.easel.MSAFile(
        str(seed_alignment_file_path), digital=True, alphabet=alphabet
    ) as alignment_file:
        alignment = alignment_file.read()
    if alignment is None:
        raise RuntimeError(
            f"Could not read a test MSA from {seed_alignment_file_path}."
        )
    alignment.name = model_name.encode("ascii")
    builder = pyhmmer.plan7.Builder(alphabet, seed=int(random_seed))
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(alignment, background)
    hmm.name = model_name.encode("ascii")
    for volatile_attribute in ("creation_time", "command_line"):
        try:
            setattr(hmm, volatile_attribute, None)
        except (AttributeError, TypeError):
            pass
    serialized_hmm = io.BytesIO()
    hmm.write(serialized_hmm)
    write_bytes_atomic(
        binary_payload=serialized_hmm.getvalue(),
        output_file_path=output_hmm_file_path,
    )


def deterministic_fake_hmm_builder(
    *,
    seed_alignment_file_path: Path,
    output_hmm_file_path: Path,
    model_name: str,
    random_seed: int,
) -> None:
    """Injected builder double: a real HMMER3 profile built without touching
    the production pyhmmer version pin or identity-resolution path."""
    build_real_amino_hmm_for_tests(
        seed_alignment_file_path=seed_alignment_file_path,
        output_hmm_file_path=output_hmm_file_path,
        model_name=model_name,
        random_seed=random_seed,
    )


def write_valid_amino_hmm(*, output_path: Path, model_name: str) -> None:
    """Write one genuine amino HMM (three-sequence toy seed) to ``output_path``."""
    seed_path = output_path.with_suffix(".seed.sto")
    write_stockholm(
        file_path=seed_path,
        sequence_by_accession={
            "AAA.1": "ACDEFGHIKL",
            "BBB.1": "ACDEFGHIKM",
            "CCC.1": "ACDEFGHLKL",
        },
    )
    build_real_amino_hmm_for_tests(
        seed_alignment_file_path=seed_path,
        output_hmm_file_path=output_path,
        model_name=model_name,
    )


def write_apaz_reference_fixture(*, reference_directory: Path) -> None:
    """Write a minimal, self-consistent APAZ v2 reference directory.

    Ships no ``split_groups/`` tree, so ``validate_apaz_hmm_build_inputs``
    exercises its own structural checks rather than delegating to the frozen
    B2 canonical validators.
    """
    reference_directory.mkdir(parents=True, exist_ok=True)
    partition_rows: list[dict[str, str]] = []
    lock_files = []
    subgroup_specs = [
        spec for spec in APAZ_MODEL_SPECS if spec.model_id != "global"
    ]
    sequence_by_subgroup_accession = {
        f"APAZ{index}.1": "ACDEFG"
        for index, _ in enumerate(subgroup_specs, start=1)
    }
    for index, spec in enumerate(subgroup_specs, start=1):
        accession = f"APAZ{index}.1"
        seed_path = reference_directory / spec.seed_file_name
        write_stockholm(
            file_path=seed_path,
            sequence_by_accession={accession: "ACDEFG"},
        )
        partition_rows.append(
            {
                "accession": accession,
                "homology_cluster_id": f"apaz_build_group_{index}",
                "split_group_id": f"apaz_build_group_{index}",
                "partition": "BUILD",
                "sequence_sha256": sequence_sha256("ACDEFG_build_" + str(index)),
            }
        )
        lock_files.append(
            {
                "path": spec.seed_file_name,
                "sha256": sha256_of_file(input_file_path=seed_path),
                "sequence_count": 1,
                "alignment_length": 6,
                "alphabet": "amino",
                "subgroup": spec.lock_subgroup,
            }
        )

    global_spec = next(
        spec for spec in APAZ_MODEL_SPECS if spec.model_id == "global"
    )
    global_seed_path = reference_directory / global_spec.seed_file_name
    write_stockholm(
        file_path=global_seed_path,
        sequence_by_accession=sequence_by_subgroup_accession,
    )
    lock_files.insert(
        0,
        {
            "path": global_spec.seed_file_name,
            "sha256": sha256_of_file(input_file_path=global_seed_path),
            "sequence_count": len(sequence_by_subgroup_accession),
            "alignment_length": 6,
            "alphabet": "amino",
            "subgroup": global_spec.lock_subgroup,
        },
    )
    partition_rows.extend(
        [
            {
                "accession": "APAZ_CAL.1",
                "homology_cluster_id": "apaz_calibration_group",
                "split_group_id": "apaz_calibration_group",
                "partition": "CALIBRATION",
                "sequence_sha256": sequence_sha256("ACDEFG_calibration"),
            },
            {
                "accession": "APAZ_HOLD.1",
                "homology_cluster_id": "apaz_holdout_group",
                "split_group_id": "apaz_holdout_group",
                "partition": "FINAL_HOLDOUT",
                "sequence_sha256": sequence_sha256("ACDEFG_holdout"),
            },
        ]
    )
    partitions_path = reference_directory / "apaz_partitions.csv"
    with partitions_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "accession", "homology_cluster_id", "split_group_id",
                "partition", "sequence_sha256",
            ],
        )
        writer.writeheader()
        writer.writerows(partition_rows)
    write_json_atomic(
        payload={
            "artifact_type": "apaz_seed_lock",
            "lock_format_version": "2.0",
            "curation_protocol_version": "fixture-2",
            "partition_protocol": "apaz_partition_v2_mmseqs90_80",
            "partition_salt": "fixture-salt",
            "files": lock_files,
            "partitions_file": {
                "path": "apaz_partitions.csv",
                "sha256": sha256_of_file(input_file_path=partitions_path),
            },
        },
        output_file_path=reference_directory / "seeds_lock.json",
    )


def refresh_apaz_partitions_lock_hash(*, reference_directory: Path) -> None:
    lock_path = reference_directory / "seeds_lock.json"
    payload = read_json(lock_path)
    partitions_entry = payload["partitions_file"]
    if not isinstance(partitions_entry, dict):
        raise AssertionError("Fixture partitions_file must be an object.")
    partitions_entry["sha256"] = sha256_of_file(
        input_file_path=reference_directory / "apaz_partitions.csv"
    )
    write_json_atomic(payload=payload, output_file_path=lock_path)


def refresh_apaz_locked_seed(
    *,
    lock_file_path: Path,
    seed_file_path: Path,
    sequence_count: int | None = None,
) -> None:
    payload = read_json(lock_file_path)
    files = payload["files"]
    if not isinstance(files, list):
        raise AssertionError("Fixture lock files must be a list.")
    entry = next(
        item
        for item in files
        if isinstance(item, dict) and item.get("path") == seed_file_path.name
    )
    entry["sha256"] = sha256_of_file(input_file_path=seed_file_path)
    if sequence_count is not None:
        entry["sequence_count"] = sequence_count
    write_json_atomic(payload=payload, output_file_path=lock_file_path)
