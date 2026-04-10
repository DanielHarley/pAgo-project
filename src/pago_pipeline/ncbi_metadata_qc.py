from __future__ import annotations

import csv
from collections import Counter, defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, TypeAlias

from src.pago_pipeline.ncbi_metadata_csv import (
    _clean_text,
    _normalize_name,
    iter_gbseq_elements,
)
from src.pago_pipeline.storage import read_json_file, sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path

CRITICAL_METADATA_COLUMNS: tuple[str, ...] = (
    "protein_uid",
    "gbseq__locus",
    "taxonomy__raw",
    "feature__keys_present",
    "reference__count",
)


@dataclass(frozen=True)
class NCBIProteinMetadataCsvQcResult:
    """
    Summary returned after auditing one metadata CSV export.
    """

    csv_file_path: Path
    metadata_manifest_file_path: Optional[Path]
    report_file_path: Optional[Path]
    row_count: int
    column_count: int
    empty_protein_uid_count: int
    duplicate_protein_uid_count: int
    fully_empty_column_count: int
    schema_match_with_manifest: Optional[bool]
    row_count_match_with_manifest: Optional[bool]
    row_count_match_with_source_xml: Optional[bool]
    normalization_collision_count: int


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _fraction(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return 0.0

    return numerator / denominator


def _load_metadata_manifest(
    *,
    metadata_manifest_file_path: Optional[PathLike],
) -> tuple[Optional[Path], Optional[dict[str, object]]]:
    if metadata_manifest_file_path is None:
        return None, None

    resolved_manifest_file_path = _as_path(metadata_manifest_file_path)
    manifest_payload = read_json_file(input_file_path=resolved_manifest_file_path)
    if not isinstance(manifest_payload, dict):
        raise RuntimeError("Metadata manifest must deserialize into a JSON object.")

    return resolved_manifest_file_path, manifest_payload


def _read_metadata_csv_summary(
    *,
    csv_file_path: PathLike,
) -> dict[str, object]:
    resolved_csv_file_path = _as_path(csv_file_path)

    with resolved_csv_file_path.open("r", encoding="utf-8", newline="") as csv_file:
        csv_reader = csv.DictReader(csv_file)
        csv_columns = list(csv_reader.fieldnames or [])

        nonempty_count_by_column = {column_name: 0 for column_name in csv_columns}
        protein_uid_counter: Counter[str] = Counter()
        row_count = 0
        empty_protein_uid_count = 0

        for row in csv_reader:
            row_count += 1

            protein_uid = _clean_text(row.get("protein_uid"))
            if protein_uid:
                protein_uid_counter[protein_uid] += 1
            else:
                empty_protein_uid_count += 1

            for column_name in csv_columns:
                cell_text = _clean_text(row.get(column_name))
                if cell_text:
                    nonempty_count_by_column[column_name] += 1

    duplicate_protein_uids = sorted(
        protein_uid
        for protein_uid, occurrence_count in protein_uid_counter.items()
        if occurrence_count > 1
    )

    completeness_by_column: dict[str, dict[str, object]] = {}
    fully_empty_columns: list[str] = []
    partially_missing_columns: list[str] = []

    for column_name in csv_columns:
        nonempty_count = nonempty_count_by_column[column_name]
        missing_count = row_count - nonempty_count
        nonempty_fraction = _fraction(nonempty_count, row_count)

        completeness_by_column[column_name] = {
            "nonempty_count": nonempty_count,
            "missing_count": missing_count,
            "nonempty_fraction": nonempty_fraction,
        }

        if nonempty_count == 0:
            fully_empty_columns.append(column_name)
        elif missing_count > 0:
            partially_missing_columns.append(column_name)

    critical_column_completeness = {
        column_name: completeness_by_column[column_name]
        for column_name in CRITICAL_METADATA_COLUMNS
        if column_name in completeness_by_column
    }

    return {
        "csv_file_path": resolved_csv_file_path,
        "row_count": row_count,
        "columns": csv_columns,
        "column_count": len(csv_columns),
        "empty_protein_uid_count": empty_protein_uid_count,
        "duplicate_protein_uid_count": len(duplicate_protein_uids),
        "duplicate_protein_uids": duplicate_protein_uids,
        "completeness_by_column": completeness_by_column,
        "critical_column_completeness": critical_column_completeness,
        "fully_empty_columns": fully_empty_columns,
        "partially_missing_columns": partially_missing_columns,
    }


def _audit_xml_flattening_risks(
    *,
    source_xml_file_path: PathLike,
) -> dict[str, object]:
    resolved_source_xml_file_path = _as_path(source_xml_file_path)

    raw_name_groups: dict[str, defaultdict[str, set[str]]] = {
        "gbseq_tags": defaultdict(set),
        "reference_tags": defaultdict(set),
        "feature_keys": defaultdict(set),
        "feature_qualifier_names": defaultdict(set),
        "interval_tags": defaultdict(set),
    }
    feature_repeat_stats: defaultdict[str, dict[str, int]] = defaultdict(
        lambda: {
            "rows_with_multiple_occurrences": 0,
            "max_occurrences_in_one_row": 0,
        }
    )

    for gbseq_element in iter_gbseq_elements(xml_file_path=resolved_source_xml_file_path):
        per_row_feature_counter: Counter[str] = Counter()

        for child_element in list(gbseq_element):
            if child_element.tag == "GBSeq_feature-table":
                for feature_element in child_element.findall("GBFeature"):
                    raw_feature_key = (
                        _clean_text(feature_element.findtext("GBFeature_key"))
                        or "unknown"
                    )
                    normalized_feature_key = _normalize_name(raw_feature_key)
                    raw_name_groups["feature_keys"][normalized_feature_key].add(
                        raw_feature_key
                    )
                    per_row_feature_counter[normalized_feature_key] += 1

                    for qualifier_element in feature_element.findall(
                        "GBFeature_quals/GBQualifier"
                    ):
                        raw_qualifier_name = _clean_text(
                            qualifier_element.findtext("GBQualifier_name")
                        )
                        if raw_qualifier_name:
                            normalized_qualifier_name = _normalize_name(
                                raw_qualifier_name
                            )
                            raw_name_groups["feature_qualifier_names"][
                                normalized_qualifier_name
                            ].add(raw_qualifier_name)

                    for interval_element in feature_element.findall(
                        "GBFeature_intervals/GBInterval"
                    ):
                        for interval_child_element in list(interval_element):
                            raw_interval_tag = interval_child_element.tag
                            normalized_interval_tag = _normalize_name(raw_interval_tag)
                            raw_name_groups["interval_tags"][
                                normalized_interval_tag
                            ].add(raw_interval_tag)

                continue

            if child_element.tag == "GBSeq_references":
                for reference_element in child_element.findall("GBReference"):
                    for reference_child_element in list(reference_element):
                        raw_reference_tag = reference_child_element.tag
                        normalized_reference_tag = _normalize_name(raw_reference_tag)
                        raw_name_groups["reference_tags"][
                            normalized_reference_tag
                        ].add(raw_reference_tag)
                continue

            raw_gbseq_tag = child_element.tag
            normalized_gbseq_tag = _normalize_name(raw_gbseq_tag)
            raw_name_groups["gbseq_tags"][normalized_gbseq_tag].add(raw_gbseq_tag)

        for feature_key, occurrence_count in per_row_feature_counter.items():
            if occurrence_count > 1:
                feature_repeat_stats[feature_key]["rows_with_multiple_occurrences"] += 1

            if (
                occurrence_count
                > feature_repeat_stats[feature_key]["max_occurrences_in_one_row"]
            ):
                feature_repeat_stats[feature_key][
                    "max_occurrences_in_one_row"
                ] = occurrence_count

    normalization_collisions = {
        category_name: {
            normalized_name: sorted(raw_names)
            for normalized_name, raw_names in normalized_name_groups.items()
            if len(raw_names) > 1
        }
        for category_name, normalized_name_groups in raw_name_groups.items()
    }
    normalization_collisions = {
        category_name: category_collisions
        for category_name, category_collisions in normalization_collisions.items()
        if category_collisions
    }

    repeated_feature_keys = {
        feature_key: repeat_stats
        for feature_key, repeat_stats in sorted(feature_repeat_stats.items())
        if repeat_stats["rows_with_multiple_occurrences"] > 0
    }

    normalization_collision_count = sum(
        len(category_collisions)
        for category_collisions in normalization_collisions.values()
    )

    return {
        "source_xml_file_path": resolved_source_xml_file_path,
        "normalization_collisions": normalization_collisions,
        "normalization_collision_count": normalization_collision_count,
        "repeated_feature_keys": repeated_feature_keys,
    }


def _build_schema_drift_summary(
    *,
    metadata_manifest_payload: Optional[Mapping[str, object]],
    reference_metadata_manifest_file_path: Optional[PathLike],
) -> Optional[dict[str, object]]:
    if reference_metadata_manifest_file_path is None:
        return None

    if metadata_manifest_payload is None:
        raise ValueError(
            "reference_metadata_manifest_file_path requires a current metadata manifest."
        )

    reference_manifest_payload = read_json_file(
        input_file_path=_as_path(reference_metadata_manifest_file_path)
    )
    if not isinstance(reference_manifest_payload, dict):
        raise RuntimeError(
            "Reference metadata manifest must deserialize into a JSON object."
        )

    current_columns = metadata_manifest_payload.get("columns")
    reference_columns = reference_manifest_payload.get("columns")
    if not isinstance(current_columns, list) or not isinstance(reference_columns, list):
        raise RuntimeError("Both metadata manifests must expose a list of columns.")

    current_column_set = {str(column_name) for column_name in current_columns}
    reference_column_set = {str(column_name) for column_name in reference_columns}

    current_feature_keys = metadata_manifest_payload.get("observed_feature_keys", [])
    reference_feature_keys = reference_manifest_payload.get("observed_feature_keys", [])
    if not isinstance(current_feature_keys, list) or not isinstance(
        reference_feature_keys,
        list,
    ):
        raise RuntimeError(
            "Both metadata manifests must expose observed_feature_keys as lists."
        )

    return {
        "reference_metadata_manifest_file_path": str(
            _as_path(reference_metadata_manifest_file_path)
        ),
        "added_columns": sorted(current_column_set - reference_column_set),
        "removed_columns": sorted(reference_column_set - current_column_set),
        "current_max_taxonomy_depth": metadata_manifest_payload.get(
            "max_taxonomy_depth"
        ),
        "reference_max_taxonomy_depth": reference_manifest_payload.get(
            "max_taxonomy_depth"
        ),
        "added_feature_keys": sorted(
            set(map(str, current_feature_keys)) - set(map(str, reference_feature_keys))
        ),
        "removed_feature_keys": sorted(
            set(map(str, reference_feature_keys)) - set(map(str, current_feature_keys))
        ),
    }


def run_ncbi_protein_metadata_csv_qc(
    *,
    csv_file_path: PathLike,
    metadata_manifest_file_path: Optional[PathLike] = None,
    source_xml_file_path: Optional[PathLike] = None,
    source_xml_manifest_payload: Optional[Mapping[str, object]] = None,
    reference_metadata_manifest_file_path: Optional[PathLike] = None,
    output_report_file_path: Optional[PathLike] = None,
) -> NCBIProteinMetadataCsvQcResult:
    """
    Audit one exported metadata CSV and optionally compare it against its source XML.
    """

    csv_summary = _read_metadata_csv_summary(csv_file_path=csv_file_path)
    resolved_csv_file_path = _as_path(csv_file_path)
    resolved_metadata_manifest_file_path, metadata_manifest_payload = (
        _load_metadata_manifest(
            metadata_manifest_file_path=metadata_manifest_file_path,
        )
    )

    row_count_match_with_manifest: Optional[bool] = None
    schema_match_with_manifest: Optional[bool] = None

    if metadata_manifest_payload is not None:
        manifest_row_count = metadata_manifest_payload.get("row_count")
        manifest_columns = metadata_manifest_payload.get("columns")

        if isinstance(manifest_row_count, int):
            row_count_match_with_manifest = (
                csv_summary["row_count"] == manifest_row_count
            )

        if isinstance(manifest_columns, list):
            schema_match_with_manifest = (
                csv_summary["columns"] == manifest_columns
            )

    row_count_match_with_source_xml: Optional[bool] = None
    xml_qc_summary: Optional[dict[str, object]] = None
    if source_xml_file_path is not None:
        xml_qc_summary = _audit_xml_flattening_risks(
            source_xml_file_path=source_xml_file_path,
        )

        if source_xml_manifest_payload is not None:
            source_xml_record_count = source_xml_manifest_payload.get(
                "consolidated_record_count"
            )
            if isinstance(source_xml_record_count, int):
                row_count_match_with_source_xml = (
                    csv_summary["row_count"] == source_xml_record_count
                )

    schema_drift_summary = _build_schema_drift_summary(
        metadata_manifest_payload=metadata_manifest_payload,
        reference_metadata_manifest_file_path=reference_metadata_manifest_file_path,
    )

    report_payload: dict[str, object] = {
        "artifact_type": "ncbi_protein_metadata_csv_qc",
        "snapshot_format_version": "1.0",
        "csv_file_path": str(resolved_csv_file_path),
        "csv_file_sha256": sha256_of_file(input_file_path=resolved_csv_file_path),
        "row_count": csv_summary["row_count"],
        "column_count": csv_summary["column_count"],
        "columns": csv_summary["columns"],
        "empty_protein_uid_count": csv_summary["empty_protein_uid_count"],
        "duplicate_protein_uid_count": csv_summary["duplicate_protein_uid_count"],
        "duplicate_protein_uids": csv_summary["duplicate_protein_uids"],
        "critical_column_completeness": csv_summary["critical_column_completeness"],
        "fully_empty_columns": csv_summary["fully_empty_columns"],
        "partially_missing_columns": csv_summary["partially_missing_columns"],
        "checks": {
            "row_count_matches_metadata_manifest": row_count_match_with_manifest,
            "schema_matches_metadata_manifest": schema_match_with_manifest,
            "row_count_matches_source_xml": row_count_match_with_source_xml,
            "protein_uid_has_no_empty_values": (
                csv_summary["empty_protein_uid_count"] == 0
            ),
            "protein_uid_has_no_duplicates": (
                csv_summary["duplicate_protein_uid_count"] == 0
            ),
        },
        "completeness_by_column": csv_summary["completeness_by_column"],
    }

    if resolved_metadata_manifest_file_path is not None:
        report_payload["metadata_manifest_file_path"] = str(
            resolved_metadata_manifest_file_path
        )
        report_payload["metadata_manifest_file_sha256"] = sha256_of_file(
            input_file_path=resolved_metadata_manifest_file_path
        )

    if metadata_manifest_payload is not None:
        report_payload["metadata_manifest_summary"] = {
            "row_count": metadata_manifest_payload.get("row_count"),
            "column_count": metadata_manifest_payload.get("column_count"),
            "max_taxonomy_depth": metadata_manifest_payload.get("max_taxonomy_depth"),
            "observed_feature_keys": metadata_manifest_payload.get(
                "observed_feature_keys"
            ),
        }

    if source_xml_file_path is not None:
        resolved_source_xml_file_path = _as_path(source_xml_file_path)
        report_payload["source_xml_file_path"] = str(resolved_source_xml_file_path)
        report_payload["source_xml_file_sha256"] = sha256_of_file(
            input_file_path=resolved_source_xml_file_path
        )

    if source_xml_manifest_payload is not None:
        report_payload["source_xml_manifest_summary"] = {
            "retrieved_at_utc": source_xml_manifest_payload.get("retrieved_at_utc"),
            "consolidated_record_count": source_xml_manifest_payload.get(
                "consolidated_record_count"
            ),
            "immutable_snapshot_relative_path": source_xml_manifest_payload.get(
                "immutable_snapshot_relative_path"
            ),
        }

    if xml_qc_summary is not None:
        report_payload["xml_flattening_risks"] = {
            "normalization_collisions": xml_qc_summary["normalization_collisions"],
            "normalization_collision_count": xml_qc_summary[
                "normalization_collision_count"
            ],
            "repeated_feature_keys": xml_qc_summary["repeated_feature_keys"],
        }

    if schema_drift_summary is not None:
        report_payload["schema_drift"] = schema_drift_summary

    resolved_output_report_file_path: Optional[Path] = None
    if output_report_file_path is not None:
        resolved_output_report_file_path = _as_path(output_report_file_path)
        write_json_atomic(
            payload=report_payload,
            output_file_path=resolved_output_report_file_path,
        )

    normalization_collision_count = 0
    if xml_qc_summary is not None:
        normalization_collision_count = int(
            xml_qc_summary["normalization_collision_count"]
        )

    return NCBIProteinMetadataCsvQcResult(
        csv_file_path=resolved_csv_file_path,
        metadata_manifest_file_path=resolved_metadata_manifest_file_path,
        report_file_path=resolved_output_report_file_path,
        row_count=int(csv_summary["row_count"]),
        column_count=int(csv_summary["column_count"]),
        empty_protein_uid_count=int(csv_summary["empty_protein_uid_count"]),
        duplicate_protein_uid_count=int(csv_summary["duplicate_protein_uid_count"]),
        fully_empty_column_count=len(csv_summary["fully_empty_columns"]),
        schema_match_with_manifest=schema_match_with_manifest,
        row_count_match_with_manifest=row_count_match_with_manifest,
        row_count_match_with_source_xml=row_count_match_with_source_xml,
        normalization_collision_count=normalization_collision_count,
    )
