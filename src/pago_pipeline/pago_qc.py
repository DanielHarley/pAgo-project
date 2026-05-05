from __future__ import annotations

from collections.abc import Sequence

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

PPIWI_RE_REGION_PATTERNS: tuple[str, ...] = (
    r"\bppiwi[_-]?re\b",
    r"\bppiwi_re_[xyz]\b",
    r"\bmid[_-]?ppiwi[_-]?re\b",
    r"\brnaseh[_-]?ppiwi[_-]?re\b",
    r"\bpiwi[_-]?re\b",
)

MID_REGION_PATTERNS: tuple[str, ...] = (
    r"\bmid\b",
    r"\bmid[_-]?ppiwi[_-]?re\b",
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

FLAG_COLUMNS: tuple[str, ...] = (
    "has_piwi_text_anywhere",
    "has_piwi_region",
    "has_classic_piwi_region",
    "has_ppiwi_re_region",
    "has_mid_region",
    "has_paz_region",
    "has_n_terminal_region",
    "has_active_site_annotation",
    "has_cdd_region",
    "has_ubig_term",
    "has_sam_methyltransferase_term",
    "is_probable_methyltransferase_noise",
    "has_conflicting_piwi_and_methyltransferase_evidence",
    "has_partial_or_fragment_keyword",
    "is_short_lt_300",
    "is_possible_partial_300_599",
    "is_compatible_long_pago_600_900",
    "is_large_or_fused_901_1300",
    "is_long_outlier_gt_1300",
)


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
    evidence_dataframe["has_piwi_text_anywhere"] = annotation_text.str.contains(
        r"\bpiwi\b|\bargonaute\b|\bpago\b",
        regex=True,
        case=False,
        na=False,
    )
    evidence_dataframe["has_classic_piwi_region"] = contains_any_pattern(
        text_series=region_text,
        regex_patterns=CLASSIC_PIWI_REGION_PATTERNS,
    )
    evidence_dataframe["has_ppiwi_re_region"] = contains_any_pattern(
        text_series=region_text,
        regex_patterns=PPIWI_RE_REGION_PATTERNS,
    )
    evidence_dataframe["has_piwi_region"] = (
        evidence_dataframe["has_classic_piwi_region"]
        | evidence_dataframe["has_ppiwi_re_region"]
    )
    evidence_dataframe["has_mid_region"] = contains_any_pattern(
        text_series=region_text,
        regex_patterns=MID_REGION_PATTERNS,
    )
    evidence_dataframe["has_paz_region"] = contains_any_pattern(
        text_series=region_text,
        regex_patterns=PAZ_REGION_PATTERNS,
    )
    evidence_dataframe["has_n_terminal_region"] = contains_any_pattern(
        text_series=region_text,
        regex_patterns=N_TERMINAL_REGION_PATTERNS,
    )
    evidence_dataframe["has_active_site_annotation"] = site_text.str.contains(
        r"\bactive\b",
        regex=True,
        case=False,
        na=False,
    )

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
    evidence_dataframe["is_probable_methyltransferase_noise"] = (
        evidence_dataframe["has_ubig_term"]
        | evidence_dataframe["has_sam_methyltransferase_term"]
    ) & ~evidence_dataframe["has_piwi_region"]
    evidence_dataframe["has_conflicting_piwi_and_methyltransferase_evidence"] = (
        evidence_dataframe["has_ubig_term"]
        | evidence_dataframe["has_sam_methyltransferase_term"]
    ) & evidence_dataframe["has_piwi_region"]

    return evidence_dataframe


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
                "pattern": "methyltransferase evidence without PIWI region evidence",
                "count": int(
                    evidence_dataframe["is_probable_methyltransferase_noise"].sum()
                ),
            },
            {
                "term_family": (
                    "conflicting_piwi_and_methyltransferase_evidence"
                ),
                "pattern": "methyltransferase evidence with PIWI region evidence",
                "count": int(
                    evidence_dataframe[
                        "has_conflicting_piwi_and_methyltransferase_evidence"
                    ].sum()
                ),
            },
        ]
    )
