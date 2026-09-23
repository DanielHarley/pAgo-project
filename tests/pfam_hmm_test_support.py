from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Iterable

from src.pago_pipeline.pfam_hmm_bundle import (
    REGISTRY_ARTIFACT_TYPE,
    REGISTRY_FORMAT_VERSION,
)
from src.pago_pipeline.storage import sha256_of_file, write_json_atomic


DEFAULT_TEST_MODELS: tuple[dict[str, object], ...] = (
    {
        "canonical_order": 1,
        "file_name": "Piwi__PF02171.hmm",
        "accession": "PF02171.22",
        "name": "Piwi",
        "alphabet": "amino",
        "ga": (21.0, 21.0),
    },
    {
        "canonical_order": 2,
        "file_name": "ArgoMid__PF16487.hmm",
        "accession": "PF16487.10",
        "name": "ArgoMid",
        "alphabet": "amino",
        "ga": (18.5, 18.5),
    },
)


class FakeAlphabet:
    def __init__(self, name: str) -> None:
        self.name = name

    def is_amino(self) -> bool:
        return self.name == "amino"


class FakeCutoffs:
    def __init__(self, gathering: tuple[float, float] | None) -> None:
        self.gathering = gathering

    def gathering_available(self) -> bool:
        return self.gathering is not None


class FakeHmm:
    def __init__(
        self,
        *,
        name: str,
        accession: str,
        alphabet: str,
        gathering: tuple[float, float] | None,
    ) -> None:
        self.name = name.encode("utf-8")
        self.accession = accession.encode("utf-8")
        self.alphabet = FakeAlphabet(alphabet)
        self.cutoffs = FakeCutoffs(gathering)


def build_test_hmm_payload(
    *,
    name: str,
    accession: str,
    alphabet: str = "amino",
    gathering: tuple[float, float] | None = (20.0, 20.0),
) -> bytes:
    lines = [
        "HMMER3/f [test fixture]",
        f"NAME  {name}",
        f"ACC   {accession}",
        f"ALPH  {alphabet}",
    ]
    if gathering is not None:
        lines.append(f"GA    {gathering[0]:.1f} {gathering[1]:.1f};")
    lines.append("//")
    return ("\n".join(lines) + "\n").encode("utf-8")


def _parse_test_hmm_file(file_path: str | Path) -> list[FakeHmm]:
    payload = Path(file_path).read_text(encoding="utf-8")
    parsed_models: list[FakeHmm] = []
    for raw_block in payload.split("//"):
        block = raw_block.strip()
        if not block:
            continue
        fields: dict[str, str] = {}
        for line in block.splitlines():
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                fields[parts[0]] = parts[1].strip().rstrip(";")
        if "NAME" not in fields:
            continue
        gathering: tuple[float, float] | None = None
        if "GA" in fields:
            values = fields["GA"].split()
            gathering = (float(values[0]), float(values[1]))
        parsed_models.append(
            FakeHmm(
                name=fields["NAME"],
                accession=fields.get("ACC", ""),
                alphabet=fields.get("ALPH", ""),
                gathering=gathering,
            )
        )
    return parsed_models


class FakeHmmFile:
    def __init__(self, file_path: str | Path) -> None:
        self.models = _parse_test_hmm_file(file_path)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        return None

    def __iter__(self):
        return iter(self.models)


def _fake_hmmpress(hmms: Iterable[FakeHmm], output: str | Path) -> None:
    model_accessions = [hmm.accession for hmm in hmms]
    canonical_payload = b"|".join(model_accessions)
    for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
        Path(str(output) + suffix).write_bytes(
            suffix.encode("ascii") + b":" + canonical_payload
        )


def build_fake_pyhmmer(
    *,
    version: str = "0.test",
    hmmpress_function=_fake_hmmpress,
):
    hmmer = SimpleNamespace()
    if hmmpress_function is not None:
        hmmer.hmmpress = hmmpress_function
    return SimpleNamespace(
        __version__=version,
        plan7=SimpleNamespace(HMMFile=FakeHmmFile),
        hmmer=hmmer,
    )


def write_test_reference_bundle(
    *,
    reference_directory: Path,
    models: tuple[dict[str, object], ...] = DEFAULT_TEST_MODELS,
    source_release: str = "test-release-1",
) -> tuple[Path, Path]:
    reference_directory.mkdir(parents=True, exist_ok=True)
    registry_models: list[dict[str, object]] = []
    for model in models:
        file_path = reference_directory / str(model["file_name"])
        file_path.write_bytes(
            build_test_hmm_payload(
                name=str(model["name"]),
                accession=str(model["accession"]),
                alphabet=str(model["alphabet"]),
                gathering=model.get("ga"),
            )
        )
        registry_models.append(
            {
                "canonical_order": model["canonical_order"],
                "file_name": model["file_name"],
                "accession": model["accession"],
                "name": model["name"],
                "alphabet": model["alphabet"],
                "model_count": 1,
                "sha256": sha256_of_file(input_file_path=file_path),
                "gathering_cutoffs_required": True,
            }
        )

    registry_file_path = reference_directory / "pfam_hmm_bundle_lock.json"
    write_json_atomic(
        payload={
            "artifact_type": REGISTRY_ARTIFACT_TYPE,
            "lock_format_version": REGISTRY_FORMAT_VERSION,
            "source_database": "Pfam",
            "source_release": source_release,
            "model_count": len(registry_models),
            "models": registry_models,
        },
        output_file_path=registry_file_path,
    )
    return registry_file_path, reference_directory


def read_registry(registry_file_path: Path) -> dict[str, object]:
    return json.loads(registry_file_path.read_text(encoding="utf-8"))


def write_registry(registry_file_path: Path, payload: dict[str, object]) -> None:
    write_json_atomic(payload=payload, output_file_path=registry_file_path)
