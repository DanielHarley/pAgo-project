from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import requests


INTERPRO_RELEASE: Final = "109.0"
PFAM_RELEASE: Final = "38.2"
LOCK_FORMAT_VERSION: Final = "1.0"
ARTIFACT_TYPE: Final = "pfam_hmm_bundle_lock"
INTERPRO_ENTRY_API: Final = "https://www.ebi.ac.uk/interpro/api/entry/pfam"
INTERPRO_RELEASE_NOTES_URL: Final = (
    "https://www.ebi.ac.uk/interpro/release_notes/109.0/"
)
INTERPRO_LICENSE_URL: Final = "https://www.ebi.ac.uk/interpro/about/license/"


@dataclass(frozen=True)
class RequestedProfile:
    accession: str
    expected_name: str
    file_stem: str
    scientific_role: str


REQUESTED_PROFILES: Final = (
    RequestedProfile("PF02171", "Piwi", "Piwi__PF02171", "Argonaute core PIWI domain"),
    RequestedProfile("PF02170", "PAZ", "PAZ__PF02170", "canonical PAZ domain"),
    RequestedProfile("PF16486", "ArgoN", "ArgoN__PF16486", "prokaryotic Argonaute N domain"),
    RequestedProfile("PF08699", "ArgoL1", "ArgoL1__PF08699", "Argonaute linker 1 domain"),
    RequestedProfile("PF16488", "ArgoL2", "ArgoL2__PF16488", "Argonaute linker 2 domain"),
    RequestedProfile("PF16487", "ArgoMid", "ArgoMid__PF16487", "Argonaute MID domain"),
    RequestedProfile("PF02146", "SIR2", "SIR2__PF02146", "SIR2-associated effector domain"),
    RequestedProfile("PF13676", "TIR_2", "TIR_2__PF13676", "TIR-associated effector domain"),
    RequestedProfile("PF01582", "TIR", "TIR__PF01582", "TIR-associated effector domain"),
    RequestedProfile("PF04471", "Mrr_cat", "Mrr_cat__PF04471", "Mrr nuclease domain"),
)


_HEADER_PATTERN: Final = re.compile(r"^(NAME|ACC|ALPH|GA)\s+(.+?)\s*$", re.MULTILINE)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _write_bytes_atomic(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        temporary_path.replace(path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


def _write_json_atomic(path: Path, payload: object) -> None:
    serialized = json.dumps(
        payload,
        ensure_ascii=False,
        indent=2,
        sort_keys=True,
    ).encode("utf-8") + b"\n"
    _write_bytes_atomic(path, serialized)


def _download_profile(
    *,
    session: requests.Session,
    profile: RequestedProfile,
    timeout_seconds: float,
) -> tuple[bytes, str]:
    url = f"{INTERPRO_ENTRY_API}/{profile.accession}/?annotation=hmm"
    response = session.get(url, timeout=timeout_seconds)
    response.raise_for_status()
    if response.headers.get("Content-Type", "").split(";", 1)[0] != "application/gzip":
        raise RuntimeError(
            f"InterPro returned an unexpected content type for {profile.accession}."
        )
    try:
        uncompressed = gzip.decompress(response.content)
    except (OSError, EOFError) as error:
        raise RuntimeError(
            f"InterPro returned an invalid gzip payload for {profile.accession}."
        ) from error
    if not uncompressed.endswith(b"//\n"):
        raise RuntimeError(
            f"Downloaded HMM for {profile.accession} has no canonical terminator."
        )
    return uncompressed, url


def _parse_header(hmm_payload: bytes) -> dict[str, str]:
    try:
        text = hmm_payload.decode("ascii")
    except UnicodeDecodeError as error:
        raise RuntimeError("Pfam HMM payload is not ASCII HMMER text.") from error
    header: dict[str, str] = {}
    for key, value in _HEADER_PATTERN.findall(text):
        header.setdefault(key, value)
    missing = sorted({"NAME", "ACC", "ALPH", "GA"} - set(header))
    if missing:
        raise RuntimeError(f"Pfam HMM header is missing fields: {missing}.")
    return header


def _validate_with_pyhmmer(
    *,
    hmm_path: Path,
    expected_accession: str,
    expected_name: str,
) -> None:
    try:
        import pyhmmer
    except ImportError as error:
        raise RuntimeError(
            "PyHMMER is required to materialize and validate the Pfam bundle."
        ) from error

    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hmm_file:
        models = list(hmm_file)
    if len(models) != 1:
        raise RuntimeError(f"{hmm_path.name} must contain exactly one HMM.")
    model = models[0]
    raw_accession = model.accession or ""
    raw_name = model.name or ""
    accession = (
        raw_accession.decode("ascii")
        if isinstance(raw_accession, bytes)
        else str(raw_accession)
    )
    name = raw_name.decode("ascii") if isinstance(raw_name, bytes) else str(raw_name)
    if accession != expected_accession or name != expected_name:
        raise RuntimeError(
            f"PyHMMER identity mismatch for {hmm_path.name}: {accession}, {name}."
        )
    if not model.alphabet.is_amino():
        raise RuntimeError(f"{hmm_path.name} is not an amino-acid HMM.")
    if not model.cutoffs.gathering_available():
        raise RuntimeError(f"{hmm_path.name} has no Pfam gathering cutoff.")


def materialize(*, output_directory: Path, timeout_seconds: float) -> Path:
    output_directory.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": (
                "pAgo-project reference materializer "
                "(offline scientific reference preparation)"
            )
        }
    )

    metadata_response = session.get(
        f"{INTERPRO_ENTRY_API}/{REQUESTED_PROFILES[0].accession}/",
        timeout=timeout_seconds,
    )
    metadata_response.raise_for_status()
    observed_interpro_release = metadata_response.headers.get("InterPro-Version")
    if observed_interpro_release != INTERPRO_RELEASE:
        raise RuntimeError(
            "InterPro release drift detected. Expected "
            f"{INTERPRO_RELEASE}, got {observed_interpro_release!r}. Review and "
            "update the pinned catalog before materializing new reference data."
        )

    model_entries: list[dict[str, object]] = []
    staged_payloads: list[tuple[Path, bytes, str, str]] = []
    for canonical_order, profile in enumerate(REQUESTED_PROFILES, start=1):
        hmm_payload, source_url = _download_profile(
            session=session,
            profile=profile,
            timeout_seconds=timeout_seconds,
        )
        header = _parse_header(hmm_payload)
        versioned_accession = header["ACC"]
        if versioned_accession.split(".", 1)[0] != profile.accession:
            raise RuntimeError(
                f"Downloaded accession mismatch for {profile.accession}: "
                f"{versioned_accession}."
            )
        if header["NAME"] != profile.expected_name:
            raise RuntimeError(
                f"Downloaded name mismatch for {profile.accession}: "
                f"{header['NAME']!r}."
            )
        if header["ALPH"].lower() not in {"amino", "amino acid"}:
            raise RuntimeError(
                f"Downloaded alphabet mismatch for {profile.accession}: "
                f"{header['ALPH']!r}."
            )

        file_name = f"{profile.file_stem}.hmm"
        final_path = output_directory / file_name
        staged_payloads.append((final_path, hmm_payload, versioned_accession, header["NAME"]))
        model_entries.append(
            {
                "accession": versioned_accession,
                "alphabet": "amino",
                "canonical_order": canonical_order,
                "file_name": file_name,
                "gathering_cutoffs_required": True,
                "model_count": 1,
                "name": header["NAME"],
                "scientific_role": profile.scientific_role,
                "sha256": _sha256(hmm_payload),
                "source_url": source_url,
            }
        )

    expected_hmm_names = {entry[0].name for entry in staged_payloads}
    unexpected_hmms = sorted(
        path.name
        for path in output_directory.glob("*.hmm")
        if path.name not in expected_hmm_names
    )
    if unexpected_hmms:
        raise RuntimeError(
            "Refusing to materialize beside undeclared HMM files: "
            + ", ".join(unexpected_hmms)
        )

    for final_path, payload, versioned_accession, expected_name in staged_payloads:
        _write_bytes_atomic(final_path, payload)
        _validate_with_pyhmmer(
            hmm_path=final_path,
            expected_accession=versioned_accession,
            expected_name=expected_name,
        )

    lock_payload = {
        "artifact_type": ARTIFACT_TYPE,
        "curation_scope": (
            "Argonaute core domains and accessory SIR2, TIR, and Mrr domains "
            "specified for Phase B"
        ),
        "interpro_release": INTERPRO_RELEASE,
        "license": "CC0-1.0",
        "license_url": INTERPRO_LICENSE_URL,
        "lock_format_version": LOCK_FORMAT_VERSION,
        "model_count": len(model_entries),
        "models": model_entries,
        "pfam_release": PFAM_RELEASE,
        "release_notes_url": INTERPRO_RELEASE_NOTES_URL,
        "source_database": "Pfam",
        "source_release": PFAM_RELEASE,
    }
    lock_path = output_directory / "pfam_hmm_bundle_lock.json"
    _write_json_atomic(lock_path, lock_payload)
    return lock_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Materialize the commit-pinned Pfam 38.2 HMM bundle used by Phase B. "
            "This is a maintenance command. Pipeline execution remains offline."
        )
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=(
            Path(__file__).resolve().parents[1]
            / "src"
            / "pago_pipeline"
            / "resources"
            / "pfam_hmm"
        ),
    )
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    arguments = parser.parse_args()
    lock_path = materialize(
        output_directory=arguments.output_directory,
        timeout_seconds=arguments.timeout_seconds,
    )
    print(lock_path)


if __name__ == "__main__":
    main()
