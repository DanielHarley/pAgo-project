from __future__ import annotations

import pandas as pd

from src.pago_pipeline.pago_qc import (
    build_pago_qc_evidence_flags,
    build_pago_qc_labelled_records,
)


def test_build_pago_qc_evidence_flags_marks_core_evidence_classes() -> None:
    metadata_dataframe = pd.DataFrame(
        [
            {
                "protein_uid": "1",
                "gbseq__length": "750",
                "feature__region__qual__region_name": "PIWI",
                "feature__region__qual__note": "PIWI domain",
                "feature__region__qual__db_xref": "CDD:441040",
                "feature__protein__qual__product": "Piwi domain-containing protein",
                "feature__site__qual__site_type": "active",
                "feature__site__qual__note": "5' RNA guide strand anchoring site",
            },
            {
                "protein_uid": "2",
                "gbseq__length": "982",
                "feature__region__qual__region_name": (
                    "pPIWI_RE_X;MID_pPIWI_RE;RNaseH_pPIWI_RE"
                ),
                "feature__region__qual__note": "restriction module",
                "feature__region__qual__db_xref": "CDD:123;CDD:456",
                "feature__protein__qual__product": (
                    "pPIWI_RE module domain-containing protein"
                ),
                "feature__site__qual__site_type": "",
                "feature__site__qual__note": "",
            },
            {
                "protein_uid": "3",
                "gbseq__length": "210",
                "feature__region__qual__region_name": "UbiG",
                "feature__region__qual__note": "methyltransferase family",
                "feature__region__qual__db_xref": "CDD:999",
                "feature__protein__qual__product": (
                    "class I SAM-dependent methyltransferase"
                ),
                "feature__site__qual__site_type": "",
                "feature__site__qual__note": "",
            },
        ]
    )

    evidence_dataframe = build_pago_qc_evidence_flags(
        metadata_dataframe=metadata_dataframe
    )

    classic_row = evidence_dataframe.iloc[0]
    assert classic_row["has_classic_piwi_region"]
    assert classic_row["has_any_piwi_related_region"]
    assert classic_row["has_any_piwi_related_evidence"]
    assert classic_row["has_piwi_region"]
    assert classic_row["has_active_or_functional_site_annotation"]
    assert classic_row["has_active_site_annotation"]
    assert classic_row["length_bin"] == "600_900"
    assert not classic_row["is_probable_methyltransferase_noise"]

    ppiwi_re_row = evidence_dataframe.iloc[1]
    assert ppiwi_re_row["has_ppiwi_re_region"]
    assert ppiwi_re_row["has_ppiwi_re_text_anywhere"]
    assert ppiwi_re_row["has_any_piwi_related_text_anywhere"]
    assert ppiwi_re_row["has_any_piwi_related_evidence"]
    assert ppiwi_re_row["has_ppiwi_re_mid_region"]
    assert ppiwi_re_row["has_mid_region"]
    assert ppiwi_re_row["has_piwi_region"]
    assert not ppiwi_re_row["is_probable_methyltransferase_noise"]

    noise_row = evidence_dataframe.iloc[2]
    assert noise_row["has_ubig_term"]
    assert noise_row["has_sam_methyltransferase_term"]
    assert noise_row["is_short_lt_300"]
    assert noise_row["is_probable_methyltransferase_noise"]
    assert noise_row["is_methyltransferase_noise_or_conflict"]


def test_build_pago_qc_labelled_records_assigns_core_decisions() -> None:
    metadata_dataframe = pd.DataFrame(
        [
            {
                "protein_uid": "1",
                "gbseq__length": "742",
                "feature__region__qual__region_name": "Piwi_piwi-like_ProArk",
                "feature__protein__qual__product": "argonaute/piwi family protein",
            },
            {
                "protein_uid": "2",
                "gbseq__length": "1032",
                "feature__region__qual__region_name": "",
                "feature__protein__qual__product": (
                    "pPIWI_RE_Z domain-containing protein"
                ),
            },
            {
                "protein_uid": "3",
                "gbseq__length": "244",
                "feature__region__qual__region_name": "Methyltransf_24;Piwi-like",
                "feature__protein__qual__product": "SAM-dependent methyltransferase",
            },
            {
                "protein_uid": "4",
                "gbseq__length": "301",
                "feature__region__qual__region_name": "",
                "feature__protein__qual__product": (
                    "SAM-dependent methyltransferase with PIWI RE module"
                ),
            },
            {
                "protein_uid": "5",
                "gbseq__length": "520",
                "feature__region__qual__region_name": "",
                "feature__protein__qual__product": "Piwi domain-containing protein",
                "gbseq__definition": "fragment",
            },
        ]
    )

    evidence_dataframe = build_pago_qc_evidence_flags(
        metadata_dataframe=metadata_dataframe
    )
    labelled_dataframe = build_pago_qc_labelled_records(
        evidence_dataframe=evidence_dataframe
    )

    classic_row = labelled_dataframe.iloc[0]
    assert classic_row["primary_label"] == "classic_piwi_candidate"
    assert classic_row["qc_decision"] == "include"
    assert classic_row["confidence_score"] >= 4

    ppiwi_re_row = labelled_dataframe.iloc[1]
    assert ppiwi_re_row["has_ppiwi_re_text_anywhere"]
    assert ppiwi_re_row["primary_label"] == "piwi_re_candidate"
    assert ppiwi_re_row["qc_decision"] == "separate_dataset"

    conflict_row = labelled_dataframe.iloc[2]
    assert conflict_row["has_conflicting_piwi_and_methyltransferase_evidence"]
    assert conflict_row["is_methyltransferase_noise_or_conflict"]
    assert (
        conflict_row["primary_label"]
        == "methyltransferase_noise_or_unresolved_conflict"
    )
    assert conflict_row["qc_decision"] == "exclude"
    assert "excluded from conservative classic pAgo positive dataset" in (
        conflict_row["rationale"]
    )

    text_conflict_row = labelled_dataframe.iloc[3]
    assert text_conflict_row["has_ppiwi_re_text_anywhere"]
    assert text_conflict_row["has_conflicting_piwi_and_methyltransferase_evidence"]
    assert not text_conflict_row["is_probable_methyltransferase_noise"]
    assert (
        text_conflict_row["primary_label"]
        == "methyltransferase_noise_or_unresolved_conflict"
    )

    weak_fragment_row = labelled_dataframe.iloc[4]
    assert weak_fragment_row["primary_label"] == "classic_piwi_candidate"
    assert weak_fragment_row["qc_decision"] == "review"
    assert weak_fragment_row["confidence_score"] == 0


def test_build_pago_qc_labelled_records_reports_missing_input_columns() -> None:
    try:
        build_pago_qc_labelled_records(evidence_dataframe=pd.DataFrame())
    except RuntimeError as exc:
        assert "missing required label input columns" in str(exc)
    else:
        raise AssertionError("Expected RuntimeError for missing label input columns.")
