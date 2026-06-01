from __future__ import annotations

import gzip
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TypeAlias
from urllib.request import Request, urlopen

import pandas as pd
import pyhmmer

from src.pago_pipeline.storage import (
    sha256_of_file,
    write_json_atomic,
)

PathLike: TypeAlias = str | Path

DEFAULT_INTERPRO_PFAM_HMM_URL_TEMPLATE = (
    "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{accession}?annotation=hmm"
)

DEFAULT_PAGO_PFAM_ACCESSIONS: tuple[str, ...] = (
    "PF16486",  # ArgoN
    "PF08699",  # ArgoL1
    "PF02170",  # PAZ
    "PF16488",  # ArgoL2
    "PF16487",  # ArgoMid
    "PF02171",  # Piwi
    "PF18157",  # MID_pPIWI_RE
    "PF18155",  # pPIWI_RE_Z
    "PF18156",  # pPIWI_RE_Y
    "PF18154",  # pPIWI_RE_REase
)

CLASSIC_AGO_PFAM_ACCESSIONS: frozenset[str] = frozenset(
    {
        "PF16486",
        "PF08699",
        "PF02170",
        "PF16488",
        "PF16487",
        "PF02171",
    }
)
PPIWI_RE_MODULE_PFAM_ACCESSIONS: frozenset[str] = frozenset(
    {
        "PF18157",
    }
)
PPIWI_RE_ASSOCIATED_PFAM_ACCESSIONS: frozenset[str] = frozenset(
    {
        "PF18155",
        "PF18156",
        "PF18154",
    }
)

DEFAULT_DOMAIN_HITS_FILE_NAME = "pyhmmer_pago_domain_hits.csv"
DEFAULT_SEQUENCE_SUMMARY_FILE_NAME = "pyhmmer_pago_sequence_summary.csv"
DEFAULT_ANNOTATED_METADATA_FILE_NAME = "sweep_genes_metadata_with_pyhmmer_pago.csv"
DEFAULT_PROFILE_HMM_FILE_NAME = "pago_pfam_profiles.hmm"
DEFAULT_MANIFEST_FILE_NAME = "manifest.json"


@dataclass(frozen=True)
class PyhmmerPagoResult:
    output_directory: Path
    profile_hmm_file_path: Path
    domain_hits_file_path: Path
    sequence_summary_file_path: Path
    annotated_metadata_file_path: Path
    manifest_file_path: Path
    domain_hits: pd.DataFrame
    sequence_summary: pd.DataFrame
    annotated_metadata: pd.DataFrame
    manifest: dict[str, Any]


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _current_utc_timestamp() -> str:
    return (
        datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat()
        .replace("+00:00", "Z")
    )


def _strip_pfam_version(accession: str | None) -> str:
    if not accession:
        return ""
    return str(accession).split(".", maxsplit=1)[0]


def _download_pfam_hmm_bytes(
    *,
    accession: str,
    url_template: str = DEFAULT_INTERPRO_PFAM_HMM_URL_TEMPLATE,
    timeout_seconds: int = 120,
) -> bytes:
    url = url_template.format(accession=accession)
    request = Request(url, headers={"User-Agent": "pAgo-project/0.1"})

    with urlopen(request, timeout=timeout_seconds) as response:
        response_payload = response.read()

    if response_payload.startswith(b"\x1f\x8b"):
        return gzip.decompress(response_payload)

    return response_payload


def download_pfam_hmm_profiles(
    *,
    accessions: Sequence[str],
    output_hmm_file_path: PathLike,
    force: bool = False,
) -> Path:
    resolved_output_hmm_file_path = _as_path(output_hmm_file_path)
    if resolved_output_hmm_file_path.exists() and not force:
        return resolved_output_hmm_file_path

    resolved_output_hmm_file_path.parent.mkdir(parents=True, exist_ok=True)
    normalized_accessions = [str(accession).strip() for accession in accessions]
    normalized_accessions = [accession for accession in normalized_accessions if accession]

    with resolved_output_hmm_file_path.open("wb") as output_handle:
        for accession in normalized_accessions:
            hmm_bytes = _download_pfam_hmm_bytes(accession=accession)
            output_handle.write(hmm_bytes.rstrip() + b"\n")

    return resolved_output_hmm_file_path


def _load_sequence_block(
    *,
    fasta_file_path: PathLike,
    target_record_ids: Iterable[str] | None = None,
) -> pyhmmer.easel.DigitalSequenceBlock:
    resolved_fasta_file_path = _as_path(fasta_file_path)
    alphabet = pyhmmer.easel.Alphabet.amino()

    with pyhmmer.easel.SequenceFile(
        str(resolved_fasta_file_path),
        digital=True,
        alphabet=alphabet,
    ) as sequence_file:
        sequence_block = sequence_file.read_block()

    if target_record_ids is None:
        return sequence_block

    target_record_id_set = {
        str(record_id).strip()
        for record_id in target_record_ids
        if str(record_id).strip()
    }
    selected_sequences = [
        sequence
        for sequence in sequence_block
        if str(sequence.name) in target_record_id_set
    ]

    if not selected_sequences:
        raise ValueError("None of the requested target_record_ids were found in the FASTA.")

    return pyhmmer.easel.DigitalSequenceBlock(alphabet, selected_sequences)


def run_pago_hmmsearch(
    *,
    fasta_file_path: PathLike,
    profile_hmm_file_path: PathLike,
    target_record_ids: Iterable[str] | None = None,
    cpus: int = 0,
    bit_cutoffs: str | None = "trusted",
    evalue_threshold: float = 1e-3,
    domain_evalue_threshold: float = 1e-3,
    include_only_reported_domains: bool = True,
    include_only_included_domains: bool = True,
) -> pd.DataFrame:
    resolved_profile_hmm_file_path = _as_path(profile_hmm_file_path)
    sequence_block = _load_sequence_block(
        fasta_file_path=fasta_file_path,
        target_record_ids=target_record_ids,
    )

    search_options: dict[str, Any] = {}
    if bit_cutoffs:
        search_options["bit_cutoffs"] = bit_cutoffs
    else:
        search_options["E"] = evalue_threshold
        search_options["domE"] = domain_evalue_threshold

    domain_hit_rows: list[dict[str, Any]] = []

    with pyhmmer.plan7.HMMFile(str(resolved_profile_hmm_file_path)) as hmm_file:
        for top_hits in pyhmmer.hmmer.hmmsearch(
            hmm_file,
            sequence_block,
            cpus=cpus,
            **search_options,
        ):
            query_hmm = top_hits.query
            profile_accession = _strip_pfam_version(query_hmm.accession)
            profile_accession_versioned = str(query_hmm.accession or "")
            profile_name = str(query_hmm.name or "")
            profile_description = str(query_hmm.description or "")

            for hit in top_hits:
                for domain_number, domain in enumerate(hit.domains, start=1):
                    if include_only_reported_domains and not domain.reported:
                        continue
                    if include_only_included_domains and not domain.included:
                        continue

                    alignment = domain.alignment
                    domain_hit_rows.append(
                        {
                            "record_id": str(hit.name),
                            "target_accession": (
                                None if hit.accession is None else str(hit.accession)
                            ),
                            "target_description": (
                                None
                                if hit.description is None
                                else str(hit.description)
                            ),
                            "target_length": int(hit.length),
                            "profile_name": profile_name,
                            "profile_accession": profile_accession,
                            "profile_accession_versioned": profile_accession_versioned,
                            "profile_description": profile_description,
                            "profile_length": int(query_hmm.M),
                            "full_sequence_evalue": float(hit.evalue),
                            "full_sequence_score": float(hit.score),
                            "full_sequence_bias": float(hit.bias),
                            "domain_number": int(domain_number),
                            "domain_i_evalue": float(domain.i_evalue),
                            "domain_c_evalue": float(domain.c_evalue),
                            "domain_score": float(domain.score),
                            "domain_bias": float(domain.bias),
                            "domain_reported": bool(domain.reported),
                            "domain_included": bool(domain.included),
                            "target_env_from": int(domain.env_from),
                            "target_env_to": int(domain.env_to),
                            "target_alignment_from": int(alignment.target_from),
                            "target_alignment_to": int(alignment.target_to),
                            "hmm_alignment_from": int(alignment.hmm_from),
                            "hmm_alignment_to": int(alignment.hmm_to),
                        }
                    )

    return pd.DataFrame(domain_hit_rows)


def summarize_pago_hmmer_hits(
    *,
    domain_hits: pd.DataFrame,
) -> pd.DataFrame:
    if domain_hits.empty:
        return pd.DataFrame(
            columns=[
                "record_id",
                "pyhmmer_pago_domain_count",
                "pyhmmer_pago_profile_count",
                "pyhmmer_pago_profiles",
                "pyhmmer_pago_best_domain_i_evalue",
                "pyhmmer_pago_best_domain_score",
                "pyhmmer_has_selected_ago_related_profile",
                "pyhmmer_has_any_classic_ago_profile",
                "pyhmmer_has_piwi_domain",
                "pyhmmer_has_paz_domain",
                "pyhmmer_has_ppiwi_re_module",
                "pyhmmer_has_ppiwi_re_accessory_profile",
                "pyhmmer_is_classic_pago_candidate",
                "pyhmmer_is_ppiwi_re_candidate",
                "pyhmmer_is_ppiwi_re_accessory_candidate",
                "pyhmmer_requires_manual_review",
                "pyhmmer_evidence_class",
            ]
        )

    summary_rows: list[dict[str, Any]] = []
    for record_id, record_hits in domain_hits.groupby("record_id", sort=True):
        profile_accessions = sorted(
            {
                str(accession)
                for accession in record_hits["profile_accession"].dropna()
                if str(accession)
            }
        )
        profile_names = sorted(
            {
                str(profile_name)
                for profile_name in record_hits["profile_name"].dropna()
                if str(profile_name)
            }
        )
        profile_accession_set = set(profile_accessions)
        has_any_classic_ago_profile = bool(
            profile_accession_set & CLASSIC_AGO_PFAM_ACCESSIONS
        )
        has_piwi_domain = "PF02171" in profile_accession_set
        has_paz_domain = "PF02170" in profile_accession_set
        has_ppiwi_re_module = bool(
            profile_accession_set & PPIWI_RE_MODULE_PFAM_ACCESSIONS
        )
        has_ppiwi_re_accessory_profile = bool(
            profile_accession_set & PPIWI_RE_ASSOCIATED_PFAM_ACCESSIONS
        )
        is_classic_pago_candidate = has_piwi_domain
        is_ppiwi_re_candidate = has_ppiwi_re_module
        is_ppiwi_re_accessory_candidate = (
            has_ppiwi_re_accessory_profile
            and not is_classic_pago_candidate
            and not is_ppiwi_re_candidate
            and not has_any_classic_ago_profile
        )
        requires_manual_review = (
            has_any_classic_ago_profile
            and not is_classic_pago_candidate
            and not is_ppiwi_re_candidate
        )

        if is_classic_pago_candidate:
            evidence_class = "classic_piwi_profile_candidate"
        elif is_ppiwi_re_candidate:
            evidence_class = "ppiwi_re_module_candidate"
        elif is_ppiwi_re_accessory_candidate:
            evidence_class = "ppiwi_re_accessory_only"
        elif requires_manual_review:
            evidence_class = "isolated_non_piwi_classic_ago_profile"
        else:
            evidence_class = "selected_profile_hit_unclassified"

        summary_rows.append(
            {
                "record_id": record_id,
                "pyhmmer_pago_domain_count": int(len(record_hits)),
                "pyhmmer_pago_profile_count": int(len(profile_accessions)),
                "pyhmmer_pago_profiles": ";".join(profile_names),
                "pyhmmer_pago_profile_accessions": ";".join(profile_accessions),
                "pyhmmer_pago_best_domain_i_evalue": float(
                    record_hits["domain_i_evalue"].min()
                ),
                "pyhmmer_pago_best_domain_score": float(
                    record_hits["domain_score"].max()
                ),
                "pyhmmer_has_selected_ago_related_profile": True,
                "pyhmmer_has_any_classic_ago_profile": has_any_classic_ago_profile,
                "pyhmmer_has_piwi_domain": has_piwi_domain,
                "pyhmmer_has_paz_domain": has_paz_domain,
                "pyhmmer_has_ppiwi_re_module": has_ppiwi_re_module,
                "pyhmmer_has_ppiwi_re_accessory_profile": (
                    has_ppiwi_re_accessory_profile
                ),
                "pyhmmer_is_classic_pago_candidate": is_classic_pago_candidate,
                "pyhmmer_is_ppiwi_re_candidate": is_ppiwi_re_candidate,
                "pyhmmer_is_ppiwi_re_accessory_candidate": (
                    is_ppiwi_re_accessory_candidate
                ),
                "pyhmmer_requires_manual_review": requires_manual_review,
                "pyhmmer_evidence_class": evidence_class,
            }
        )

    return pd.DataFrame(summary_rows)


def annotate_metadata_with_pyhmmer(
    *,
    sequence_metadata: pd.DataFrame,
    sequence_summary: pd.DataFrame,
) -> pd.DataFrame:
    annotated_metadata = sequence_metadata.merge(
        sequence_summary,
        on="record_id",
        how="left",
    )

    boolean_columns = [
        "pyhmmer_has_selected_ago_related_profile",
        "pyhmmer_has_any_classic_ago_profile",
        "pyhmmer_has_piwi_domain",
        "pyhmmer_has_paz_domain",
        "pyhmmer_has_ppiwi_re_module",
        "pyhmmer_has_ppiwi_re_accessory_profile",
        "pyhmmer_is_classic_pago_candidate",
        "pyhmmer_is_ppiwi_re_candidate",
        "pyhmmer_is_ppiwi_re_accessory_candidate",
        "pyhmmer_requires_manual_review",
    ]
    integer_columns = [
        "pyhmmer_pago_domain_count",
        "pyhmmer_pago_profile_count",
    ]

    for column in boolean_columns:
        if column in annotated_metadata:
            annotated_metadata[column] = annotated_metadata[column].fillna(False)

    for column in integer_columns:
        if column in annotated_metadata:
            annotated_metadata[column] = (
                annotated_metadata[column].fillna(0).astype("int64")
            )

    for column in [
        "pyhmmer_pago_profiles",
        "pyhmmer_pago_profile_accessions",
    ]:
        if column in annotated_metadata:
            annotated_metadata[column] = annotated_metadata[column].fillna("")

    if "pyhmmer_evidence_class" in annotated_metadata:
        annotated_metadata["pyhmmer_evidence_class"] = annotated_metadata[
            "pyhmmer_evidence_class"
        ].fillna("no_selected_profile_hit")

    return annotated_metadata


def build_pyhmmer_pago_annotations(
    *,
    protein_fasta_file_path: PathLike,
    sequence_metadata: pd.DataFrame,
    output_directory: PathLike,
    profile_accessions: Sequence[str] = DEFAULT_PAGO_PFAM_ACCESSIONS,
    target_record_ids: Iterable[str] | None = None,
    cpus: int = 0,
    bit_cutoffs: str | None = "trusted",
    force_profile_download: bool = False,
) -> PyhmmerPagoResult:
    resolved_output_directory = _as_path(output_directory)
    resolved_output_directory.mkdir(parents=True, exist_ok=True)

    resolved_profile_hmm_file_path = (
        resolved_output_directory / DEFAULT_PROFILE_HMM_FILE_NAME
    )
    domain_hits_file_path = resolved_output_directory / DEFAULT_DOMAIN_HITS_FILE_NAME
    sequence_summary_file_path = (
        resolved_output_directory / DEFAULT_SEQUENCE_SUMMARY_FILE_NAME
    )
    annotated_metadata_file_path = (
        resolved_output_directory / DEFAULT_ANNOTATED_METADATA_FILE_NAME
    )
    manifest_file_path = resolved_output_directory / DEFAULT_MANIFEST_FILE_NAME

    profile_hmm_file_path = download_pfam_hmm_profiles(
        accessions=profile_accessions,
        output_hmm_file_path=resolved_profile_hmm_file_path,
        force=force_profile_download,
    )
    effective_target_record_ids = (
        target_record_ids
        if target_record_ids is not None
        else sequence_metadata["record_id"].astype(str).tolist()
    )
    domain_hits = run_pago_hmmsearch(
        fasta_file_path=protein_fasta_file_path,
        profile_hmm_file_path=profile_hmm_file_path,
        target_record_ids=effective_target_record_ids,
        cpus=cpus,
        bit_cutoffs=bit_cutoffs,
    )
    sequence_summary = summarize_pago_hmmer_hits(domain_hits=domain_hits)
    annotated_metadata = annotate_metadata_with_pyhmmer(
        sequence_metadata=sequence_metadata,
        sequence_summary=sequence_summary,
    )

    domain_hits.to_csv(domain_hits_file_path, index=False)
    sequence_summary.to_csv(sequence_summary_file_path, index=False)
    annotated_metadata.to_csv(annotated_metadata_file_path, index=False)

    classic_piwi_candidate_count = int(
        annotated_metadata["pyhmmer_is_classic_pago_candidate"].sum()
    )
    ppiwi_re_module_candidate_count = int(
        annotated_metadata["pyhmmer_is_ppiwi_re_candidate"].sum()
    )
    ppiwi_re_accessory_only_count = int(
        annotated_metadata["pyhmmer_is_ppiwi_re_accessory_candidate"].sum()
    )
    isolated_non_piwi_classic_ago_profile_count = int(
        (
            annotated_metadata["pyhmmer_evidence_class"]
            == "isolated_non_piwi_classic_ago_profile"
        ).sum()
    )

    manifest: dict[str, Any] = {
        "artifact_type": "pyhmmer_pago_annotation",
        "snapshot_created_at_utc": _current_utc_timestamp(),
        "protein_fasta_file_path": str(_as_path(protein_fasta_file_path)),
        "protein_fasta_file_sha256": sha256_of_file(
            input_file_path=protein_fasta_file_path,
        ),
        "profile_hmm_file_name": profile_hmm_file_path.name,
        "profile_hmm_file_sha256": sha256_of_file(
            input_file_path=profile_hmm_file_path,
        ),
        "profile_accessions": list(profile_accessions),
        "pyhmmer_version": pyhmmer.__version__,
        "cpus": cpus,
        "bit_cutoffs": bit_cutoffs,
        "input_sequence_count": int(len(sequence_metadata)),
        "domain_hit_count": int(len(domain_hits)),
        "sequence_with_selected_profile_hit_count": int(len(sequence_summary)),
        "sequence_with_piwi_domain_count": int(
            annotated_metadata["pyhmmer_has_piwi_domain"].sum()
        ),
        "sequence_with_ppiwi_re_module_count": ppiwi_re_module_candidate_count,
        "sequence_with_ppiwi_re_accessory_only_count": (
            ppiwi_re_accessory_only_count
        ),
        "sequence_with_classic_piwi_profile_candidate_count": (
            classic_piwi_candidate_count
        ),
        "sequence_with_ppiwi_re_module_candidate_count": (
            ppiwi_re_module_candidate_count
        ),
        "sequence_with_isolated_non_piwi_classic_ago_profile_count": (
            isolated_non_piwi_classic_ago_profile_count
        ),
        "domain_hits_file_name": domain_hits_file_path.name,
        "sequence_summary_file_name": sequence_summary_file_path.name,
        "annotated_metadata_file_name": annotated_metadata_file_path.name,
    }
    write_json_atomic(payload=manifest, output_file_path=manifest_file_path)

    return PyhmmerPagoResult(
        output_directory=resolved_output_directory,
        profile_hmm_file_path=profile_hmm_file_path,
        domain_hits_file_path=domain_hits_file_path,
        sequence_summary_file_path=sequence_summary_file_path,
        annotated_metadata_file_path=annotated_metadata_file_path,
        manifest_file_path=manifest_file_path,
        domain_hits=domain_hits,
        sequence_summary=sequence_summary,
        annotated_metadata=annotated_metadata,
        manifest=manifest,
    )
