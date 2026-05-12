from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from src.pago_pipeline.pago_qc import PagoQcDecision, PagoQcPrimaryLabel


@dataclass(frozen=True)
class PagoQcFilteredDatasetDefinition:
    dataset_name: str
    file_name: str
    description: str
    qc_decision: str
    primary_label: str | None = None

    @property
    def selection_rule(self) -> str:
        if self.primary_label is None:
            return f"qc_decision == {self.qc_decision!r}"
        return (
            f"primary_label == {self.primary_label!r} "
            f"and qc_decision == {self.qc_decision!r}"
        )


CLASSIC_PAGO_HIGH_PRECISION_DATASET_NAME = "classic_pago_high_precision"
CLASSIC_PAGO_REVIEW_DATASET_NAME = "classic_pago_review"
PIWI_RE_DATASET_NAME = "piwi_re"
EXCLUDED_DATASET_NAME = "excluded"

FILTERED_DATASET_DEFINITIONS: tuple[PagoQcFilteredDatasetDefinition, ...] = (
    PagoQcFilteredDatasetDefinition(
        dataset_name=CLASSIC_PAGO_HIGH_PRECISION_DATASET_NAME,
        file_name="classic_pago_high_precision_records.csv",
        description="Conservative classic pAgo positive dataset.",
        qc_decision=PagoQcDecision.INCLUDE.value,
        primary_label=PagoQcPrimaryLabel.CLASSIC_PIWI_CANDIDATE.value,
    ),
    PagoQcFilteredDatasetDefinition(
        dataset_name=CLASSIC_PAGO_REVIEW_DATASET_NAME,
        file_name="classic_pago_review_records.csv",
        description="Classic pAgo candidates held for review before inclusion.",
        qc_decision=PagoQcDecision.REVIEW.value,
        primary_label=PagoQcPrimaryLabel.CLASSIC_PIWI_CANDIDATE.value,
    ),
    PagoQcFilteredDatasetDefinition(
        dataset_name=PIWI_RE_DATASET_NAME,
        file_name="piwi_re_records.csv",
        description="PIWI-RE candidates materialized as a separate dataset.",
        qc_decision=PagoQcDecision.SEPARATE_DATASET.value,
        primary_label=PagoQcPrimaryLabel.PIWI_RE_CANDIDATE.value,
    ),
    PagoQcFilteredDatasetDefinition(
        dataset_name=EXCLUDED_DATASET_NAME,
        file_name="excluded_records.csv",
        description=(
            "Records excluded from the conservative classic pAgo positive "
            "dataset. This subset is not a biologically validated negative class."
        ),
        qc_decision=PagoQcDecision.EXCLUDE.value,
    ),
)

REQUIRED_FILTER_INPUT_COLUMNS: tuple[str, ...] = (
    "protein_uid",
    "primary_label",
    "qc_decision",
)

FILTERED_DATASET_COUNT_COLUMNS: tuple[str, ...] = (
    "dataset_name",
    "count",
    "fraction",
    "selection_rule",
    "description",
)


def build_filter_policy_payload() -> list[dict[str, object]]:
    return [
        {
            "dataset_name": dataset_definition.dataset_name,
            "file_name": dataset_definition.file_name,
            "primary_label": dataset_definition.primary_label,
            "qc_decision": dataset_definition.qc_decision,
            "selection_rule": dataset_definition.selection_rule,
            "description": dataset_definition.description,
        }
        for dataset_definition in FILTERED_DATASET_DEFINITIONS
    ]


def build_filter_policy_sha256() -> str:
    serialized_payload = json.dumps(
        build_filter_policy_payload(),
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(serialized_payload.encode("utf-8")).hexdigest()


def _sum_boolean_masks(
    *,
    boolean_mask_by_name: Mapping[str, pd.Series],
) -> pd.Series:
    membership_count: pd.Series | None = None
    for boolean_mask in boolean_mask_by_name.values():
        numeric_mask = boolean_mask.astype(int)
        if membership_count is None:
            membership_count = numeric_mask
        else:
            membership_count = membership_count + numeric_mask

    if membership_count is None:
        return pd.Series(dtype=int)
    return membership_count


def validate_filter_input_columns(
    *,
    labelled_records_dataframe: pd.DataFrame,
) -> None:
    missing_columns = [
        column_name
        for column_name in REQUIRED_FILTER_INPUT_COLUMNS
        if column_name not in labelled_records_dataframe.columns
    ]
    if missing_columns:
        raise RuntimeError(
            "labelled_records_dataframe is missing required filter input "
            f"columns: {missing_columns}."
        )


def _copy_matching_records(
    *,
    labelled_records_dataframe: pd.DataFrame,
    mask: pd.Series,
) -> pd.DataFrame:
    return labelled_records_dataframe.loc[mask].copy().reset_index(drop=True)


def _build_mask_for_dataset_definition(
    *,
    labelled_records_dataframe: pd.DataFrame,
    dataset_definition: PagoQcFilteredDatasetDefinition,
) -> pd.Series:
    mask = labelled_records_dataframe["qc_decision"] == dataset_definition.qc_decision
    if dataset_definition.primary_label is not None:
        mask = (
            mask
            & (
                labelled_records_dataframe["primary_label"]
                == dataset_definition.primary_label
            )
        )
    return mask


def build_pago_qc_filtered_datasets(
    *,
    labelled_records_dataframe: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    validate_filter_input_columns(
        labelled_records_dataframe=labelled_records_dataframe,
    )

    mask_by_dataset_name = {
        dataset_definition.dataset_name: _build_mask_for_dataset_definition(
            labelled_records_dataframe=labelled_records_dataframe,
            dataset_definition=dataset_definition,
        )
        for dataset_definition in FILTERED_DATASET_DEFINITIONS
    }
    membership_count = _sum_boolean_masks(boolean_mask_by_name=mask_by_dataset_name)
    unassigned_record_count = int((membership_count == 0).sum())
    multiply_assigned_record_count = int((membership_count > 1).sum())
    if unassigned_record_count or multiply_assigned_record_count:
        raise RuntimeError(
            "Filtered dataset rules must assign each labelled record to exactly "
            "one dataset. "
            f"Unassigned records: {unassigned_record_count}; "
            f"multiply assigned records: {multiply_assigned_record_count}."
        )

    return {
        dataset_name: _copy_matching_records(
            labelled_records_dataframe=labelled_records_dataframe,
            mask=mask,
        )
        for dataset_name, mask in mask_by_dataset_name.items()
    }


def build_filtered_dataset_counts_dataframe(
    *,
    filtered_dataset_by_name: Mapping[str, pd.DataFrame],
    total_record_count: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for dataset_definition in FILTERED_DATASET_DEFINITIONS:
        if dataset_definition.dataset_name not in filtered_dataset_by_name:
            raise RuntimeError(
                "filtered_dataset_by_name is missing required dataset "
                f"{dataset_definition.dataset_name!r}."
            )
        dataset_dataframe = filtered_dataset_by_name[dataset_definition.dataset_name]
        dataset_count = int(len(dataset_dataframe))
        rows.append(
            {
                "dataset_name": dataset_definition.dataset_name,
                "count": dataset_count,
                "fraction": (
                    float(dataset_count / total_record_count)
                    if total_record_count
                    else 0.0
                ),
                "selection_rule": dataset_definition.selection_rule,
                "description": dataset_definition.description,
            }
        )

    assigned_record_count = int(sum(row["count"] for row in rows))
    if assigned_record_count != total_record_count:
        raise RuntimeError(
            "Filtered dataset counts must sum to total_record_count. "
            f"Expected {total_record_count}, got {assigned_record_count}."
        )

    return pd.DataFrame(rows, columns=list(FILTERED_DATASET_COUNT_COLUMNS))
