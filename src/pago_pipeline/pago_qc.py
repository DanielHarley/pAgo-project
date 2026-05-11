from __future__ import annotations

from collections.abc import Sequence
from enum import Enum

import pandas as pd

REGION_TEXT_COLUMNS: tuple[str, ...] = (
    "feature__region__qual__region_name",
    "feature__region__qual__note",
    "feature__region__qual__db_xref",
)

SITE_TEXT_COLUMNS: tuple[str, ...] = (
    "feature__site__qual__site_type",
    "feature__site__qual__note",
    "feature__site__qual__db_xref",
)

ANNOTATION_TEXT_COLUMNS: tuple[str, ...] = (
    "feature__region__qual__region_name",
    "feature__region__qual__note",
    "feature__region__qual__db_xref",
    "feature__protein__qual__product",
    "feature__protein__qual__name",
    "feature__protein__qual__note",
    "gbseq__definition",
    "gbseq__comment",
    "gbseq__keywords__keyword",
)

CORE_OUTPUT_COLUMNS: tuple[str, ...] = (
    "protein_uid",
    "gbseq__accession_version",
    "gbseq__definition",
    "gbseq__organism",
    "taxonomy__raw",
    "taxonomy__03",
    "taxonomy__04",
    "feature__protein__qual__product",
    "feature__region__qual__region_name",
    "feature__region__qual__db_xref",
    "feature__site__qual__site_type",
    "feature__site__qual__note",
)

CLASSIC_PIWI_REGION_PATTERNS: tuple[str, ...] = (
    r"\bpiwi\b",
    r"\bpiwi-like\b",
    r"\bpiwi_piwi-like_proark\b",
    r"\bargonaute\b",
)

CLASSIC_PIWI_TEXT_PATTERNS: tuple[str, ...] = (
    r"\bpiwi\b",
    r"\bargonaute\b",
    r"\bpago\b",
)

PPIWI_RE_REGION_PATTERNS: tuple[str, ...] = (
    r"\bppiwi[_\-\s]?re(?:[_\-\s]?[xyz])?\b",
    r"\bmid[_\-\s]?ppiwi[_\-\s]?re\b",
    r"\brnaseh[_\-\s]?ppiwi[_\-\s]?re\b",
    r"\bpiwi[_\-\s]?re\b",
)

PPIWI_RE_TEXT_PATTERNS: tuple[str, ...] = (
    r"\bppiwi[_\-\s]?re(?:[_\-\s]?[xyz])?\b",
    r"\bpiwi[_\-\s]?re\b",
    r"\bppiwi[_\-\s]?re\s+module\b",
)

CLASSIC_MID_REGION_PATTERNS: tuple[str, ...] = (r"\bmid\b",)

PPIWI_RE_MID_REGION_PATTERNS: tuple[str, ...] = (
    r"\bmid[_\-\s]?ppiwi[_\-\s]?re\b",
)

PAZ_REGION_PATTERNS: tuple[str, ...] = (
    r"\bpaz\b",
    r"\bpaz_[0-9]+\b",
)

N_TERMINAL_REGION_PATTERNS: tuple[str, ...] = (
    r"\bn-terminal\b",
    r"\bargonaute n\b",
    r"\bago n\b",
)

UBIG_PATTERNS: tuple[str, ...] = (r"\bubig\b",)

SAM_METHYLTRANSFERASE_PATTERNS: tuple[str, ...] = (
    r"\bsam-dependent methyltransferase\b",
    r"\bclass i sam-dependent methyltransferase\b",
    r"\bmethyltransferase domain-containing protein\b",
)

PARTIAL_OR_FRAGMENT_PATTERNS: tuple[str, ...] = (
    r"\bpartial\b",
    r"\bfragment\b",
    r"\btruncated\b",
    r"\bincomplete\b",
)

ACTIVE_OR_FUNCTIONAL_SITE_PATTERNS: tuple[str, ...] = (
    r"\bactive\b",
    r"\bfunctional\b",
    r"\bbinding\b",
    r"\banchor(?:ing)?\b",
    r"\bcatalytic\b",
    r"\bmetal\b",
)

DEPRECATED_ALIAS_BY_COLUMN: dict[str, str] = {
    "has_piwi_region": "has_any_piwi_related_region",
    "has_piwi_text_anywhere": "has_classic_piwi_text_anywhere",
    "has_active_site_annotation": "has_active_or_functional_site_annotation",
}

FLAG_COLUMNS: tuple[str, ...] = (
    "has_classic_piwi_text_anywhere",
    "has_ppiwi_re_text_anywhere",
    "has_any_piwi_related_text_anywhere",
    "has_piwi_text_anywhere",
    "has_any_piwi_related_region",
    "has_any_piwi_related_evidence",
    "has_piwi_region",
    "has_classic_piwi_region",
    "has_ppiwi_re_region",
    "has_classic_mid_region",
    "has_ppiwi_re_mid_region",
    "has_mid_region",
    "has_paz_region",
    "has_n_terminal_region",
    "has_active_or_functional_site_annotation",
    "has_active_site_annotation",
    "has_cdd_region",
    "has_ubig_term",
    "has_sam_methyltransferase_term",
    "is_probable_methyltransferase_noise",
    "has_conflicting_piwi_and_methyltransferase_evidence",
    "is_methyltransferase_noise_or_conflict",
    "has_partial_or_fragment_keyword",
    "is_short_lt_300",
    "is_possible_partial_300_599",
    "is_compatible_long_pago_600_900",
    "is_large_or_fused_901_1300",
    "is_long_outlier_gt_1300",
)

REQUIRED_LABEL_INPUT_COLUMNS: tuple[str, ...] = (
    "is_methyltransferase_noise_or_conflict",
    "has_ppiwi_re_region",
    "has_ppiwi_re_text_anywhere",
    "has_classic_piwi_region",
    "has_classic_piwi_text_anywhere",
    "has_partial_or_fragment_keyword",
    "length_bin",
    "is_compatible_long_pago_600_900",
    "has_active_or_functional_site_annotation",
    "is_large_or_fused_901_1300",
    "is_possible_partial_300_599",
    "has_ppiwi_re_mid_region",
    "has_ubig_term",
    "has_sam_methyltransferase_term",
    "has_conflicting_piwi_and_methyltransferase_evidence",
    "is_short_lt_300",
    "is_long_outlier_gt_1300",
)


class PagoQcPrimaryLabel(str, Enum):
    CLASSIC_PIWI_CANDIDATE = "classic_piwi_candidate"
    PIWI_RE_CANDIDATE = "piwi_re_candidate"
    METHYLTRANSFERASE_NOISE_OR_UNRESOLVED_CONFLICT = (
        "methyltransferase_noise_or_unresolved_conflict"
    )
    LOW_EVIDENCE_OR_UNRELATED = "low_evidence_or_unrelated"


class PagoQcDecision(str, Enum):
    INCLUDE = "include"
    EXCLUDE = "exclude"
    REVIEW = "review"
    SEPARATE_DATASET = "separate_dataset"


def combine_text_columns(
    *,
    dataframe: pd.DataFrame,
    column_names: Sequence[str],
) -> pd.Series:
    combined = pd.Series("", index=dataframe.index, dtype="string")
    for column_name in column_names:
        if column_name not in dataframe.columns:
            continue
        column_values = dataframe[column_name].fillna("").astype(str)
        combined = combined.str.cat(column_values, sep=" ")
    return combined.str.lower()


def assign_length_bin(sequence_length: object) -> str:
    numeric_length = pd.to_numeric(sequence_length, errors="coerce")
    if pd.isna(numeric_length):
        return "unknown"
    if numeric_length < 300:
        return "lt_300"
    if 300 <= numeric_length <= 599:
        return "300_599"
    if 600 <= numeric_length <= 900:
        return "600_900"
    if 901 <= numeric_length <= 1300:
        return "901_1300"
    return "gt_1300"


def contains_any_pattern(
    *,
    text_series: pd.Series,
    regex_patterns: Sequence[str],
) -> pd.Series:
    result_mask = pd.Series(False, index=text_series.index)
    for regex_pattern in regex_patterns:
        result_mask = result_mask | text_series.str.contains(
            regex_pattern,
            regex=True,
            case=False,
            na=False,
        )
    return result_mask


def build_pago_qc_evidence_flags(
    *,
    metadata_dataframe: pd.DataFrame,
) -> pd.DataFrame:
    if "protein_uid" not in metadata_dataframe.columns:
        raise RuntimeError("metadata_dataframe must contain protein_uid.")

    metadata_dataframe_copy = metadata_dataframe.copy()
    metadata_dataframe_copy["protein_uid"] = (
        metadata_dataframe_copy["protein_uid"].astype(str).str.strip()
    )

    existing_core_columns = [
        column_name
        for column_name in CORE_OUTPUT_COLUMNS
        if column_name in metadata_dataframe_copy.columns
    ]
    evidence_dataframe = metadata_dataframe_copy[existing_core_columns].copy()

    sequence_length = pd.to_numeric(
        metadata_dataframe_copy.get("gbseq__length"),
        errors="coerce",
    )
    region_text = combine_text_columns(
        dataframe=metadata_dataframe_copy,
        column_names=REGION_TEXT_COLUMNS,
    )
    site_text = combine_text_columns(
        dataframe=metadata_dataframe_copy,
        column_names=SITE_TEXT_COLUMNS,
    )
    annotation_text = combine_text_columns(
        dataframe=metadata_dataframe_copy,
        column_names=ANNOTATION_TEXT_COLUMNS,
    )

    evidence_dataframe["sequence_length"] = sequence_length
    evidence_dataframe["length_bin"] = sequence_length.apply(assign_length_bin)
    evidence_dataframe["has_classic_piwi_text_anywhere"] = contains_any_pattern(
        text_series=annotation_text,
        regex_patterns=CLASSIC_PIWI_TEXT_PATTERNS,
    )
    evidence_dataframe["has_ppiwi_re_text_anywhere"] = contains_any_pattern(
        text_series=annotation_text,
        regex_patterns=PPIWI_RE_TEXT_PATTERNS,
    )
    evidence_dataframe["has_any_piwi_related_text_anywhere"] = (
        evidence_dataframe["has_classic_piwi_text_anywhere"]
        | evidence_dataframe["has_ppiwi_re_text_anywhere"]
    )
    evidence_dataframe["has_piwi_text_anywhere"] = evidence_dataframe[
        "has_classic_piwi_text_anywhere"
    ]
    evidence_dataframe["has_classic_piwi_region"] = contains_any_pattern(
        text_series=region_text,
        regex_patterns=CLASSIC_PIWI_REGION_PATTERNS,
    )
    evidence_dataframe["has_ppiwi_re_region"] = contains_any_pattern(
        text_series=region_text,
        regex_patterns=PPIWI_RE_REGION_PATTERNS,
    )
    evidence_dataframe["has_any_piwi_related_region"] = (
        evidence_dataframe["has_classic_piwi_region"]
        | evidence_dataframe["has_ppiwi_re_region"]
    )
    evidence_dataframe["has_any_piwi_related_evidence"] = (
        evidence_dataframe["has_any_piwi_related_region"]
        | evidence_dataframe["has_any_piwi_related_text_anywhere"]
    )
    evidence_dataframe["has_piwi_region"] = evidence_dataframe[
        "has_any_piwi_related_region"
    ]
    evidence_dataframe["has_ppiwi_re_mid_region"] = contains_any_pattern(
        text_series=region_text,
        regex_patterns=PPIWI_RE_MID_REGION_PATTERNS,
    )
    evidence_dataframe["has_classic_mid_region"] = contains_any_pattern(
        text_series=region_text,
        regex_patterns=CLASSIC_MID_REGION_PATTERNS,
    ) & ~evidence_dataframe["has_ppiwi_re_mid_region"]
    evidence_dataframe["has_mid_region"] = (
        evidence_dataframe["has_classic_mid_region"]
        | evidence_dataframe["has_ppiwi_re_mid_region"]
    )
    evidence_dataframe["has_paz_region"] = contains_any_pattern(
        text_series=region_text,
        regex_patterns=PAZ_REGION_PATTERNS,
    )
    evidence_dataframe["has_n_terminal_region"] = contains_any_pattern(
        text_series=region_text,
        regex_patterns=N_TERMINAL_REGION_PATTERNS,
    )
    evidence_dataframe["has_active_or_functional_site_annotation"] = (
        contains_any_pattern(
            text_series=site_text,
            regex_patterns=ACTIVE_OR_FUNCTIONAL_SITE_PATTERNS,
        )
    )
    evidence_dataframe["has_active_site_annotation"] = evidence_dataframe[
        "has_active_or_functional_site_annotation"
    ]

    region_db_xref_column = "feature__region__qual__db_xref"
    if region_db_xref_column in metadata_dataframe_copy.columns:
        evidence_dataframe["has_cdd_region"] = (
            metadata_dataframe_copy[region_db_xref_column]
            .fillna("")
            .astype(str)
            .str.contains(
                r"CDD:",
                regex=True,
                case=False,
                na=False,
            )
        )
    else:
        evidence_dataframe["has_cdd_region"] = False

    evidence_dataframe["has_ubig_term"] = contains_any_pattern(
        text_series=annotation_text,
        regex_patterns=UBIG_PATTERNS,
    )
    evidence_dataframe["has_sam_methyltransferase_term"] = contains_any_pattern(
        text_series=annotation_text,
        regex_patterns=SAM_METHYLTRANSFERASE_PATTERNS,
    )
    evidence_dataframe["has_partial_or_fragment_keyword"] = contains_any_pattern(
        text_series=annotation_text,
        regex_patterns=PARTIAL_OR_FRAGMENT_PATTERNS,
    )
    evidence_dataframe["is_short_lt_300"] = sequence_length < 300
    evidence_dataframe["is_possible_partial_300_599"] = sequence_length.between(
        300,
        599,
        inclusive="both",
    )
    evidence_dataframe["is_compatible_long_pago_600_900"] = sequence_length.between(
        600,
        900,
        inclusive="both",
    )
    evidence_dataframe["is_large_or_fused_901_1300"] = sequence_length.between(
        901,
        1300,
        inclusive="both",
    )
    evidence_dataframe["is_long_outlier_gt_1300"] = sequence_length > 1300
    has_methyltransferase_evidence = (
        evidence_dataframe["has_ubig_term"]
        | evidence_dataframe["has_sam_methyltransferase_term"]
    )
    evidence_dataframe["is_probable_methyltransferase_noise"] = (
        has_methyltransferase_evidence
        & ~evidence_dataframe["has_any_piwi_related_evidence"]
    )
    evidence_dataframe["has_conflicting_piwi_and_methyltransferase_evidence"] = (
        has_methyltransferase_evidence
        & evidence_dataframe["has_any_piwi_related_evidence"]
    )
    evidence_dataframe["is_methyltransferase_noise_or_conflict"] = (
        evidence_dataframe["is_probable_methyltransferase_noise"]
        | evidence_dataframe["has_conflicting_piwi_and_methyltransferase_evidence"]
    )

    return evidence_dataframe


def assign_pago_qc_primary_label(row: pd.Series) -> str:
    if row["is_methyltransferase_noise_or_conflict"]:
        return PagoQcPrimaryLabel.METHYLTRANSFERASE_NOISE_OR_UNRESOLVED_CONFLICT.value
    if row["has_ppiwi_re_region"] or row["has_ppiwi_re_text_anywhere"]:
        return PagoQcPrimaryLabel.PIWI_RE_CANDIDATE.value
    if row["has_classic_piwi_region"] or row["has_classic_piwi_text_anywhere"]:
        return PagoQcPrimaryLabel.CLASSIC_PIWI_CANDIDATE.value
    return PagoQcPrimaryLabel.LOW_EVIDENCE_OR_UNRELATED.value


def assign_pago_qc_decision(row: pd.Series) -> str:
    primary_label = row["primary_label"]
    length_bin = row["length_bin"]

    if (
        primary_label
        == PagoQcPrimaryLabel.METHYLTRANSFERASE_NOISE_OR_UNRESOLVED_CONFLICT.value
    ):
        return PagoQcDecision.EXCLUDE.value
    if primary_label == PagoQcPrimaryLabel.PIWI_RE_CANDIDATE.value:
        return PagoQcDecision.SEPARATE_DATASET.value
    if primary_label == PagoQcPrimaryLabel.CLASSIC_PIWI_CANDIDATE.value:
        if row["has_partial_or_fragment_keyword"]:
            return PagoQcDecision.REVIEW.value
        if length_bin == "600_900":
            return PagoQcDecision.INCLUDE.value
        if length_bin in {"300_599", "901_1300"}:
            return PagoQcDecision.REVIEW.value
        return PagoQcDecision.EXCLUDE.value
    return PagoQcDecision.EXCLUDE.value


def assign_pago_qc_confidence_score(row: pd.Series) -> int:
    score = 0
    primary_label = row["primary_label"]

    if primary_label == PagoQcPrimaryLabel.CLASSIC_PIWI_CANDIDATE.value:
        if row["has_classic_piwi_region"]:
            score += 3
        if row["has_classic_piwi_text_anywhere"]:
            score += 1
        if row["is_compatible_long_pago_600_900"]:
            score += 2
        if row["has_active_or_functional_site_annotation"]:
            score += 1
        if row["has_partial_or_fragment_keyword"]:
            score -= 2
        if row["is_large_or_fused_901_1300"]:
            score -= 1
        if row["is_possible_partial_300_599"]:
            score -= 1
    elif primary_label == PagoQcPrimaryLabel.PIWI_RE_CANDIDATE.value:
        if row["has_ppiwi_re_region"]:
            score += 3
        if row["has_ppiwi_re_text_anywhere"]:
            score += 1
        if row["has_ppiwi_re_mid_region"]:
            score += 1
        if row["has_partial_or_fragment_keyword"]:
            score -= 1
    elif (
        primary_label
        == PagoQcPrimaryLabel.METHYLTRANSFERASE_NOISE_OR_UNRESOLVED_CONFLICT.value
    ):
        if row["has_ubig_term"]:
            score += 2
        if row["has_sam_methyltransferase_term"]:
            score += 2
        if row["has_conflicting_piwi_and_methyltransferase_evidence"]:
            score += 1

    return max(0, int(score))


def build_pago_qc_rationale(row: pd.Series) -> str:
    reasons: list[str] = []

    if row["has_classic_piwi_region"]:
        reasons.append("classic PIWI/Argonaute region evidence")
    if row["has_classic_piwi_text_anywhere"]:
        reasons.append("classic PIWI/Argonaute text evidence")
    if row["has_ppiwi_re_region"]:
        reasons.append("PIWI-RE region evidence")
    if row["has_ppiwi_re_text_anywhere"]:
        reasons.append("PIWI-RE text evidence")
    if row["has_ubig_term"]:
        reasons.append("UbiG methyltransferase evidence")
    if row["has_sam_methyltransferase_term"]:
        reasons.append("SAM-dependent methyltransferase evidence")
    if row["has_conflicting_piwi_and_methyltransferase_evidence"]:
        reasons.append(
            "conflicting PIWI-like and methyltransferase evidence; "
            "excluded from conservative classic pAgo positive dataset "
            "pending review"
        )
    if row["has_partial_or_fragment_keyword"]:
        reasons.append("partial/fragment/truncated keyword")
    if row["is_compatible_long_pago_600_900"]:
        reasons.append("length compatible with long pAgo candidate")
    if row["is_possible_partial_300_599"]:
        reasons.append("length 300-599 aa, possible partial/truncated candidate")
    if row["is_large_or_fused_901_1300"]:
        reasons.append("length 901-1300 aa, possible large/fused candidate")
    if row["is_short_lt_300"]:
        reasons.append("length <300 aa")
    if row["is_long_outlier_gt_1300"]:
        reasons.append("length >1300 aa")

    return "; ".join(reasons) if reasons else "no strong pAgo-related evidence"


def validate_label_input_columns(
    *,
    evidence_dataframe: pd.DataFrame,
) -> None:
    missing_columns = [
        column_name
        for column_name in REQUIRED_LABEL_INPUT_COLUMNS
        if column_name not in evidence_dataframe.columns
    ]
    if missing_columns:
        raise RuntimeError(
            "evidence_dataframe is missing required label input columns: "
            f"{missing_columns}."
        )


def build_pago_qc_labelled_records(
    *,
    evidence_dataframe: pd.DataFrame,
) -> pd.DataFrame:
    validate_label_input_columns(evidence_dataframe=evidence_dataframe)
    labelled_dataframe = evidence_dataframe.copy()
    labelled_dataframe["primary_label"] = labelled_dataframe.apply(
        assign_pago_qc_primary_label,
        axis=1,
    )
    labelled_dataframe["qc_decision"] = labelled_dataframe.apply(
        assign_pago_qc_decision,
        axis=1,
    )
    labelled_dataframe["confidence_score"] = labelled_dataframe.apply(
        assign_pago_qc_confidence_score,
        axis=1,
    )
    labelled_dataframe["rationale"] = labelled_dataframe.apply(
        build_pago_qc_rationale,
        axis=1,
    )
    return labelled_dataframe


def build_value_counts_dataframe(
    *,
    dataframe: pd.DataFrame,
    column_name: str,
) -> pd.DataFrame:
    value_counts = dataframe[column_name].value_counts(dropna=False)
    total_count = int(value_counts.sum())
    return pd.DataFrame(
        [
            {
                column_name: value,
                "count": int(count),
                "fraction": float(count / total_count) if total_count else 0.0,
            }
            for value, count in value_counts.items()
        ]
    )


def build_label_counts_dataframe(
    *,
    labelled_dataframe: pd.DataFrame,
) -> pd.DataFrame:
    return build_value_counts_dataframe(
        dataframe=labelled_dataframe,
        column_name="primary_label",
    )


def build_filter_decision_counts_dataframe(
    *,
    labelled_dataframe: pd.DataFrame,
) -> pd.DataFrame:
    return build_value_counts_dataframe(
        dataframe=labelled_dataframe,
        column_name="qc_decision",
    )


def build_evidence_counts_dataframe(
    *,
    evidence_dataframe: pd.DataFrame,
    flag_columns: Sequence[str] = FLAG_COLUMNS,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "flag": flag_column,
                "count": int(evidence_dataframe[flag_column].sum()),
                "fraction": float(evidence_dataframe[flag_column].mean()),
            }
            for flag_column in flag_columns
        ]
    )


def build_top_value_count_dataframe(
    *,
    dataframe: pd.DataFrame,
    column_name: str,
) -> pd.DataFrame:
    if column_name not in dataframe.columns:
        return pd.DataFrame({column_name: [], "count": []})

    return (
        dataframe[column_name]
        .fillna("")
        .astype(str)
        .value_counts(dropna=False)
        .rename_axis(column_name)
        .reset_index(name="count")
    )


def build_suspicious_terms_dataframe(
    *,
    evidence_dataframe: pd.DataFrame,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "term_family": "ubig",
                "pattern": "|".join(UBIG_PATTERNS),
                "count": int(evidence_dataframe["has_ubig_term"].sum()),
            },
            {
                "term_family": "sam_methyltransferase",
                "pattern": "|".join(SAM_METHYLTRANSFERASE_PATTERNS),
                "count": int(
                    evidence_dataframe["has_sam_methyltransferase_term"].sum()
                ),
            },
            {
                "term_family": "probable_methyltransferase_noise",
                "pattern": (
                    "methyltransferase evidence without PIWI-related evidence"
                ),
                "count": int(
                    evidence_dataframe["is_probable_methyltransferase_noise"].sum()
                ),
            },
            {
                "term_family": (
                    "conflicting_piwi_and_methyltransferase_evidence"
                ),
                "pattern": "methyltransferase evidence with PIWI-related evidence",
                "count": int(
                    evidence_dataframe[
                        "has_conflicting_piwi_and_methyltransferase_evidence"
                    ].sum()
                ),
            },
        ]
    )
