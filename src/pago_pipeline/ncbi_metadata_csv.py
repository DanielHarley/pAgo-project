from __future__ import annotations

import csv
import re
import tempfile
from xml.etree import ElementTree
from collections import defaultdict
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Optional, TypeAlias, cast

from src.pago_pipeline.ncbi_snapshot import load_xml_snapshot_by_directory
from src.pago_pipeline.storage import sha256_of_file, write_json_atomic

PathLike: TypeAlias = str | Path
FeatureQualifierMap: TypeAlias = dict[str, list[str]]

_GB_PREFIXES_TO_STRIP = (
    "GBSeq_",
    "GBReference_",
    "GBFeature_",
    "GBQualifier_",
    "GBInterval_",
)


@dataclass(frozen=True)
class NCBIProteinMetadataCsvExportResult:
    """
    Summary returned after exporting a flattened CSV from one XML snapshot.
    """

    xml_file_path: Path
    csv_file_path: Path
    manifest_file_path: Optional[Path]
    row_count: int
    column_count: int
    columns: list[str]
    max_taxonomy_depth: int
    observed_feature_keys: list[str]
    observed_feature_qualifiers: FeatureQualifierMap
    xml_file_sha256: str
    csv_file_sha256: str


@dataclass(frozen=True)
class NCBIProteinMetadataSchemaInspectionResult:
    """
    Dynamic schema discovered from a consolidated XML snapshot.
    """

    xml_file_path: Path
    row_count: int
    columns: list[str]
    max_taxonomy_depth: int
    observed_feature_keys: list[str]
    observed_feature_qualifiers: FeatureQualifierMap


def _as_path(path_like: PathLike) -> Path:
    return Path(path_like)


def _clean_text(text: str | None) -> str:
    return (text or "").strip()


def _normalize_name(text: str) -> str:
    normalized_text = text.strip()

    for prefix in _GB_PREFIXES_TO_STRIP:
        if normalized_text.startswith(prefix):
            normalized_text = normalized_text[len(prefix) :]
            break

    # Generic fallback for tag families such as GBAuthor / GBKeyword.
    if (
        normalized_text.startswith("GB")
        and len(normalized_text) > 2
        and normalized_text[2].isupper()
    ):
        normalized_text = normalized_text[2:]

    normalized_text = normalized_text.replace("-", "_")
    if re.fullmatch(r"[A-Z0-9_]+", normalized_text):
        normalized_text = normalized_text.lower()
    else:
        normalized_text = re.sub(r"(?<!^)(?=[A-Z])", "_", normalized_text)
        normalized_text = normalized_text.lower()
    normalized_text = re.sub(r"[^a-z0-9_]+", "_", normalized_text)
    normalized_text = re.sub(r"_+", "_", normalized_text).strip("_")

    return normalized_text or "value"


def _join_values(values: list[str], delimiter: str = ";") -> str:
    return delimiter.join(value for value in values if value)


def _new_row_lists() -> DefaultDict[str, list[str]]:
    return defaultdict(list)


def _new_observed_feature_qualifiers() -> DefaultDict[str, set[str]]:
    return defaultdict(set)


def _extract_ncbi_protein_uid_from_gbseq_element(
    *,
    gbseq_element: ElementTree.Element,
) -> str:
    uid_candidates: set[str] = set()

    for gbseqid_element in gbseq_element.findall(".//GBSeqid"):
        gbseqid_text = _clean_text(gbseqid_element.text)
        if not gbseqid_text:
            continue

        uid_match = re.fullmatch(r"gi\|(\d+)\|?", gbseqid_text)
        if uid_match is not None:
            uid_candidates.add(uid_match.group(1))
            continue

        if gbseqid_text.isdigit():
            uid_candidates.add(gbseqid_text)

    if not uid_candidates:
        raise RuntimeError(
            "Failed to extract a protein UID from GBSeqid fields in one XML record."
        )

    if len(uid_candidates) != 1:
        raise RuntimeError(
            "Found multiple conflicting protein UID candidates in one XML record: "
            f"{sorted(uid_candidates)}."
        )

    return next(iter(uid_candidates))


def iter_gbseq_elements(
    *,
    xml_file_path: PathLike,
) -> Iterator[ElementTree.Element]:
    """
    Yield one GBSeq element at a time using iterparse to limit memory growth.
    """

    resolved_xml_file_path = _as_path(xml_file_path)
    context = ElementTree.iterparse(resolved_xml_file_path, events=("start", "end"))

    _, root_element = next(context)
    if root_element.tag != "GBSet":
        raise RuntimeError(
            "XML metadata export expects a GBSet root element. "
            f"Found {root_element.tag!r}."
        )

    for event, element in context:
        if event == "end" and element.tag == "GBSeq":
            yield element
            element.clear()
            root_element.clear()


def _collect_leaf_text_values(
    *,
    element: ElementTree.Element,
    prefix_parts: list[str],
    row_lists: DefaultDict[str, list[str]],
) -> None:
    child_elements = list(element)

    if not child_elements:
        text_value = _clean_text(element.text)
        if text_value:
            row_lists["__".join(prefix_parts)].append(text_value)
        return

    for child_element in child_elements:
        _collect_leaf_text_values(
            element=child_element,
            prefix_parts=prefix_parts + [_normalize_name(child_element.tag)],
            row_lists=row_lists,
        )


def _append_reference_values(
    *,
    reference_elements: list[ElementTree.Element],
    row_lists: DefaultDict[str, list[str]],
) -> None:
    row_lists["reference__count"].append(str(len(reference_elements)))

    for reference_element in reference_elements:
        for child_element in list(reference_element):
            if child_element.tag == "GBReference_authors":
                author_values = [
                    _clean_text(author_element.text)
                    for author_element in child_element.findall("GBAuthor")
                ]
                for author_value in author_values:
                    if author_value:
                        row_lists["reference__authors"].append(author_value)
                continue

            child_tag_name = _normalize_name(child_element.tag)
            if list(child_element):
                _collect_leaf_text_values(
                    element=child_element,
                    prefix_parts=["reference", child_tag_name],
                    row_lists=row_lists,
                )
                continue

            child_text = _clean_text(child_element.text)
            if child_text:
                row_lists[f"reference__{child_tag_name}"].append(child_text)


def _append_feature_values(
    *,
    feature_elements: list[ElementTree.Element],
    row_lists: DefaultDict[str, list[str]],
    observed_feature_keys: set[str],
    observed_feature_qualifiers: DefaultDict[str, set[str]],
) -> None:
    for feature_element in feature_elements:
        feature_key_text = _clean_text(feature_element.findtext("GBFeature_key"))
        feature_key_name = _normalize_name(feature_key_text or "unknown")

        observed_feature_keys.add(feature_key_name)
        row_lists["feature__keys_present"].append(feature_key_name)

        feature_location_text = _clean_text(
            feature_element.findtext("GBFeature_location")
        )
        if feature_location_text:
            row_lists[f"feature__{feature_key_name}__location"].append(
                feature_location_text
            )

        for interval_element in feature_element.findall("GBFeature_intervals/GBInterval"):
            for interval_child_element in list(interval_element):
                interval_child_text = _clean_text(interval_child_element.text)
                if not interval_child_text:
                    continue

                interval_child_name = _normalize_name(interval_child_element.tag)
                row_lists[
                    f"feature__{feature_key_name}__interval__{interval_child_name}"
                ].append(interval_child_text)

        for qualifier_element in feature_element.findall("GBFeature_quals/GBQualifier"):
            qualifier_name_text = _clean_text(
                qualifier_element.findtext("GBQualifier_name")
            )
            qualifier_value_text = _clean_text(
                qualifier_element.findtext("GBQualifier_value")
            )

            if not qualifier_name_text or not qualifier_value_text:
                continue

            qualifier_name = _normalize_name(qualifier_name_text)
            observed_feature_qualifiers[feature_key_name].add(qualifier_name)
            row_lists[
                f"feature__{feature_key_name}__qual__{qualifier_name}"
            ].append(qualifier_value_text)


def flatten_gbseq_metadata(
    *,
    gbseq_element: ElementTree.Element,
    observed_feature_keys: Optional[set[str]] = None,
    observed_feature_qualifiers: Optional[DefaultDict[str, set[str]]] = None,
) -> dict[str, str]:
    """
    Flatten one GBSeq record into a CSV-ready row with dynamic columns.
    """

    row_lists = _new_row_lists()
    row_values: dict[str, str] = {
        "protein_uid": _extract_ncbi_protein_uid_from_gbseq_element(
            gbseq_element=gbseq_element
        )
    }

    resolved_observed_feature_keys = (
        observed_feature_keys if observed_feature_keys is not None else set()
    )
    resolved_observed_feature_qualifiers = (
        observed_feature_qualifiers
        if observed_feature_qualifiers is not None
        else _new_observed_feature_qualifiers()
    )

    for child_element in list(gbseq_element):
        if child_element.tag == "GBSeq_feature-table":
            _append_feature_values(
                feature_elements=child_element.findall("GBFeature"),
                row_lists=row_lists,
                observed_feature_keys=resolved_observed_feature_keys,
                observed_feature_qualifiers=resolved_observed_feature_qualifiers,
            )
            continue

        if child_element.tag == "GBSeq_references":
            _append_reference_values(
                reference_elements=child_element.findall("GBReference"),
                row_lists=row_lists,
            )
            continue

        if child_element.tag == "GBSeq_taxonomy":
            taxonomy_text = _clean_text(child_element.text)
            if taxonomy_text:
                row_values["taxonomy__raw"] = taxonomy_text
                taxonomy_tokens = [
                    token.strip()
                    for token in taxonomy_text.split(";")
                    if token.strip()
                ]
                for taxonomy_index, taxonomy_token in enumerate(
                    taxonomy_tokens,
                    start=1,
                ):
                    row_values[f"taxonomy__{taxonomy_index:02d}"] = taxonomy_token
            continue

        if child_element.tag == "GBSeq_other-seqids":
            for other_seqid_element in child_element.findall("GBSeqid"):
                other_seqid_text = _clean_text(other_seqid_element.text)
                if other_seqid_text:
                    row_lists["gbseq__other_seqids"].append(other_seqid_text)
            continue

        child_name = _normalize_name(child_element.tag)
        if list(child_element):
            _collect_leaf_text_values(
                element=child_element,
                prefix_parts=["gbseq", child_name],
                row_lists=row_lists,
            )
            continue

        child_text = _clean_text(child_element.text)
        if child_text:
            row_values[f"gbseq__{child_name}"] = child_text

    for column_name, column_values in row_lists.items():
        joined_text = _join_values(column_values)
        if joined_text:
            row_values[column_name] = joined_text

    return row_values


def _column_sort_key(column_name: str) -> tuple[int, int, str]:
    if column_name == "protein_uid":
        return 0, 0, column_name

    if column_name.startswith("gbseq__"):
        return 1, 0, column_name

    if column_name == "taxonomy__raw":
        return 2, 0, column_name

    taxonomy_match = re.fullmatch(r"taxonomy__(\d+)", column_name)
    if taxonomy_match is not None:
        return 2, int(taxonomy_match.group(1)), column_name

    if column_name == "reference__count":
        return 3, 0, column_name

    if column_name.startswith("reference__"):
        return 3, 1, column_name

    if column_name == "feature__keys_present":
        return 4, 0, column_name

    if column_name.startswith("feature__"):
        return 4, 1, column_name

    return 5, 0, column_name


def inspect_ncbi_protein_metadata_schema(
    *,
    xml_file_path: PathLike,
) -> NCBIProteinMetadataSchemaInspectionResult:
    """
    Inspect one consolidated XML file and return the dynamic CSV schema metadata.
    """

    resolved_xml_file_path = _as_path(xml_file_path)
    discovered_columns: set[str] = {"protein_uid"}
    row_count = 0
    max_taxonomy_depth = 0
    observed_feature_keys: set[str] = set()
    observed_feature_qualifiers = _new_observed_feature_qualifiers()

    for gbseq_element in iter_gbseq_elements(xml_file_path=resolved_xml_file_path):
        row_values = flatten_gbseq_metadata(
            gbseq_element=gbseq_element,
            observed_feature_keys=observed_feature_keys,
            observed_feature_qualifiers=observed_feature_qualifiers,
        )
        discovered_columns.update(row_values.keys())
        row_count += 1

        taxonomy_depth = len(
            [column_name for column_name in row_values if re.fullmatch(r"taxonomy__\d+", column_name)]
        )
        if taxonomy_depth > max_taxonomy_depth:
            max_taxonomy_depth = taxonomy_depth

    sorted_columns = sorted(discovered_columns, key=_column_sort_key)

    return NCBIProteinMetadataSchemaInspectionResult(
        xml_file_path=resolved_xml_file_path,
        row_count=row_count,
        columns=sorted_columns,
        max_taxonomy_depth=max_taxonomy_depth,
        observed_feature_keys=sorted(observed_feature_keys),
        observed_feature_qualifiers={
            feature_key: sorted(qualifier_names)
            for feature_key, qualifier_names in sorted(
                observed_feature_qualifiers.items()
            )
        },
    )


def _write_csv_atomic(
    *,
    output_csv_file_path: PathLike,
    columns: list[str],
    xml_file_path: PathLike,
) -> Path:
    resolved_output_csv_file_path = _as_path(output_csv_file_path)
    resolved_output_csv_file_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        dir=resolved_output_csv_file_path.parent,
        encoding="utf-8",
        newline="",
    ) as temporary_file:
        resolved_temporary_file_path = Path(temporary_file.name)
        csv_writer = csv.DictWriter(
            temporary_file,
            fieldnames=columns,
            extrasaction="ignore",
        )
        csv_writer.writeheader()

        for gbseq_element in iter_gbseq_elements(xml_file_path=xml_file_path):
            csv_writer.writerow(flatten_gbseq_metadata(gbseq_element=gbseq_element))

    resolved_temporary_file_path.replace(resolved_output_csv_file_path)
    return resolved_output_csv_file_path


def export_ncbi_protein_metadata_csv(
    *,
    xml_file_path: PathLike,
    output_csv_file_path: PathLike,
    output_manifest_file_path: Optional[PathLike] = None,
    source_snapshot_payload: Optional[Mapping[str, object]] = None,
) -> NCBIProteinMetadataCsvExportResult:
    """
    Export one wide CSV where each row corresponds to one GBSeq record.
    """

    resolved_xml_file_path = _as_path(xml_file_path)
    schema_metadata = inspect_ncbi_protein_metadata_schema(
        xml_file_path=resolved_xml_file_path,
    )
    columns = schema_metadata.columns

    resolved_output_csv_file_path = _write_csv_atomic(
        output_csv_file_path=output_csv_file_path,
        columns=columns,
        xml_file_path=resolved_xml_file_path,
    )

    xml_file_sha256 = sha256_of_file(input_file_path=resolved_xml_file_path)
    csv_file_sha256 = sha256_of_file(input_file_path=resolved_output_csv_file_path)

    resolved_output_manifest_file_path: Optional[Path] = None
    if output_manifest_file_path is not None:
        manifest_file_path = _as_path(output_manifest_file_path)
        resolved_output_manifest_file_path = manifest_file_path
        manifest_payload: dict[str, object] = {
            "artifact_type": "ncbi_protein_metadata_csv",
            "snapshot_format_version": "1.0",
            "xml_file_name": resolved_xml_file_path.name,
            "xml_file_path": str(resolved_xml_file_path),
            "xml_file_sha256": xml_file_sha256,
            "csv_file_name": resolved_output_csv_file_path.name,
            "csv_file_path": str(resolved_output_csv_file_path),
            "csv_file_sha256": csv_file_sha256,
            "row_count": schema_metadata.row_count,
            "column_count": len(columns),
            "columns": columns,
            "max_taxonomy_depth": schema_metadata.max_taxonomy_depth,
            "observed_feature_keys": schema_metadata.observed_feature_keys,
            "observed_feature_qualifiers": schema_metadata.observed_feature_qualifiers,
        }

        if source_snapshot_payload is not None:
            source_manifest = source_snapshot_payload.get("manifest")
            source_manifest_file_path = source_snapshot_payload.get("manifest_file_path")

            if isinstance(source_manifest, Mapping):
                manifest_payload["source_xml_snapshot_relative_path"] = source_manifest.get(
                    "immutable_snapshot_relative_path"
                )
                manifest_payload["source_xml_snapshot_directory_name"] = source_manifest.get(
                    "immutable_snapshot_directory_name"
                )

            if isinstance(source_manifest_file_path, (str, Path)):
                manifest_payload["source_xml_snapshot_manifest_sha256"] = sha256_of_file(
                    input_file_path=source_manifest_file_path
                )

        write_json_atomic(
            payload=manifest_payload,
            output_file_path=manifest_file_path,
        )

    return NCBIProteinMetadataCsvExportResult(
        xml_file_path=resolved_xml_file_path,
        csv_file_path=resolved_output_csv_file_path,
        manifest_file_path=resolved_output_manifest_file_path,
        row_count=schema_metadata.row_count,
        column_count=len(columns),
        columns=columns,
        max_taxonomy_depth=schema_metadata.max_taxonomy_depth,
        observed_feature_keys=schema_metadata.observed_feature_keys,
        observed_feature_qualifiers=schema_metadata.observed_feature_qualifiers,
        xml_file_sha256=xml_file_sha256,
        csv_file_sha256=csv_file_sha256,
    )


def export_ncbi_protein_metadata_csv_from_xml_snapshot(
    *,
    snapshot_directory: PathLike,
    output_csv_file_path: PathLike,
    output_manifest_file_path: Optional[PathLike] = None,
) -> NCBIProteinMetadataCsvExportResult:
    """
    Resolve one validated XML snapshot directory and export its flattened CSV.
    """

    snapshot_payload = load_xml_snapshot_by_directory(
        snapshot_directory=snapshot_directory,
    )
    xml_file_path = snapshot_payload.get("xml_file_path")
    if not isinstance(xml_file_path, Path):
        raise RuntimeError(
            "Resolved XML snapshot payload is missing a valid xml_file_path."
        )

    return export_ncbi_protein_metadata_csv(
        xml_file_path=xml_file_path,
        output_csv_file_path=output_csv_file_path,
        output_manifest_file_path=output_manifest_file_path,
        source_snapshot_payload=cast(Mapping[str, object], snapshot_payload),
    )
