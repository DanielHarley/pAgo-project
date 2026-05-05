from __future__ import annotations

import pandas as pd

from src.pago_pipeline.pago_qc import build_pago_qc_evidence_flags


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
    assert classic_row["has_piwi_region"]
    assert classic_row["has_active_site_annotation"]
    assert not classic_row["is_probable_methyltransferase_noise"]

    ppiwi_re_row = evidence_dataframe.iloc[1]
    assert ppiwi_re_row["has_ppiwi_re_region"]
    assert ppiwi_re_row["has_mid_region"]
    assert ppiwi_re_row["has_piwi_region"]
    assert not ppiwi_re_row["is_probable_methyltransferase_noise"]

    noise_row = evidence_dataframe.iloc[2]
    assert noise_row["has_ubig_term"]
    assert noise_row["has_sam_methyltransferase_term"]
    assert noise_row["is_short_lt_300"]
    assert noise_row["is_probable_methyltransferase_noise"]
