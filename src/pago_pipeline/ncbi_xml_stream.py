from __future__ import annotations

import hashlib
import io
import re
import tempfile
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Union

PathLike = Union[str, Path]

DEFAULT_CONSOLIDATED_ROOT_TAG = "GBSet"
DEFAULT_CONSOLIDATED_RECORD_TAG = "GBSeq"

_GI_PROTEIN_UID_PATTERN = re.compile(r"gi\|(\d+)\|?")

_XML_DECLARATION_BYTES = b"<?xml version='1.0' encoding='utf-8'?>\n"


def _as_path(path_like: PathLike) -> Path:
    """
    Convert PathLike input into pathlib.Path explicitly.
    """
    return Path(path_like)


@dataclass(frozen=True)
class XmlBatchPayloadSource:
    """
    One ordered consolidation input, held either in memory or on disk.

    Streaming persistence keeps fetched batches on disk, so consolidation must
    be able to consume a batch without the caller first materializing all
    payloads in memory. In-memory payloads remain supported for small runs and
    for tests.
    """
    batch_index: int
    payload_bytes: Optional[bytes] = None
    payload_file_path: Optional[PathLike] = None

    def open_payload_stream(self):
        if self.payload_bytes is not None:
            return io.BytesIO(self.payload_bytes)

        if self.payload_file_path is None:
            raise ValueError(
                "XmlBatchPayloadSource requires payload_bytes or "
                f"payload_file_path for batch {self.batch_index}."
            )

        return _as_path(self.payload_file_path).open("rb")


@dataclass(frozen=True)
class ConsolidatedXmlWriteResult:
    """
    Outcome of one streaming consolidation write.
    """
    output_file_path: Path
    root_tag: str
    record_count: int
    byte_count: int
    sha256: str


@dataclass(frozen=True)
class XmlBatchValidationResult:
    """
    Outcome of validating one fetched XML batch payload in a single pass.
    """
    root_tag: str
    record_count: int
    protein_uids: List[str]


@dataclass
class ConsolidatedXmlValidationResult:
    """
    Outcome of one streaming validation pass over a consolidated XML file.
    """
    root_tag: str
    record_count: int
    protein_uids: List[str] = field(default_factory=list)


def extract_protein_uid_from_gbseq_element(
    *,
    gbseq_element: ET.Element,
) -> str:
    """
    Extract one protein UID from a GBSeq record using its GBSeqid fields.
    """
    uid_candidates: set[str] = set()

    for gbseqid_element in gbseq_element.findall(".//GBSeqid"):
        gbseqid_text = (gbseqid_element.text or "").strip()
        if not gbseqid_text:
            continue

        uid_match = _GI_PROTEIN_UID_PATTERN.fullmatch(gbseqid_text)
        if uid_match is not None:
            uid_candidates.add(uid_match.group(1))

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


def _stream_root_children(
    *,
    payload_stream,
    on_root_element: Optional[Callable[[str, Dict[str, str]], None]],
    on_record_element: Callable[[ET.Element], None],
    on_serializable_record: Optional[Callable[[ET.Element], None]] = None,
) -> str:
    """
    Stream the direct children of an XML root without holding the document.

    ``on_record_element`` observes each record as soon as its end tag is read.
    ``on_serializable_record`` observes the same record one event later, which
    is the earliest point at which ElementTree has assigned the record's tail
    text. Byte-exact re-serialization requires that tail, so consolidation uses
    the delayed callback while validation uses the immediate one.
    """
    parser_events = ET.iterparse(payload_stream, events=("start", "end"))
    element_depth = 0
    root_tag: Optional[str] = None
    root_element: Optional[ET.Element] = None
    pending_record_element: Optional[ET.Element] = None

    def flush_pending_record() -> None:
        nonlocal pending_record_element

        if pending_record_element is None:
            return

        if on_serializable_record is not None:
            on_serializable_record(pending_record_element)

        if root_element is not None:
            try:
                root_element.remove(pending_record_element)
            except ValueError:
                pass

        pending_record_element.clear()
        pending_record_element = None

    for parser_event, element in parser_events:
        if parser_event == "start":
            element_depth += 1
            if element_depth == 1:
                root_element = element
                root_tag = element.tag
                if on_root_element is not None:
                    on_root_element(element.tag, dict(element.attrib))
            elif element_depth == 2:
                flush_pending_record()
            continue

        element_depth -= 1
        if element_depth == 1:
            on_record_element(element)
            if on_serializable_record is None:
                if root_element is not None:
                    try:
                        root_element.remove(element)
                    except ValueError:
                        pass
                element.clear()
            else:
                pending_record_element = element
        elif element_depth == 0:
            flush_pending_record()

    if root_tag is None:
        raise ET.ParseError("no element found")

    return root_tag


def validate_xml_batch_payload(
    *,
    xml_payload_bytes: bytes,
    expected_protein_uids: Optional[Sequence[str]] = None,
    expected_record_tag: str = DEFAULT_CONSOLIDATED_RECORD_TAG,
    validate_protein_uids: bool = True,
) -> XmlBatchValidationResult:
    """
    Validate one fetched XML batch payload structurally in a single pass.

    The batch is checked for well-formedness, record tag consistency, record
    count, and, when requested, the exact set of protein UIDs it carries. The
    exact-UID check is what makes an omitted or substituted record diagnosable
    at the batch that produced it rather than at final consolidation.
    """
    extracted_protein_uids: List[str] = []
    record_count = 0
    invalid_record_tags: set[str] = set()

    def observe_record(record_element: ET.Element) -> None:
        nonlocal record_count

        record_count += 1
        if record_element.tag != expected_record_tag:
            invalid_record_tags.add(record_element.tag)
            return

        if validate_protein_uids:
            extracted_protein_uids.append(
                extract_protein_uid_from_gbseq_element(
                    gbseq_element=record_element,
                )
            )

    root_tag = _stream_root_children(
        payload_stream=io.BytesIO(xml_payload_bytes),
        on_root_element=None,
        on_record_element=observe_record,
    )

    if invalid_record_tags:
        raise RuntimeError(
            "NCBI XML batch contains unexpected child tags under "
            f"{root_tag}: {sorted(invalid_record_tags)}."
        )

    if expected_protein_uids is not None and validate_protein_uids:
        expected_uid_counter = Counter(expected_protein_uids)
        extracted_uid_counter = Counter(extracted_protein_uids)
        if expected_uid_counter != extracted_uid_counter:
            missing_uids = sorted(
                (expected_uid_counter - extracted_uid_counter).elements()
            )
            unexpected_uids = sorted(
                (extracted_uid_counter - expected_uid_counter).elements()
            )
            raise RuntimeError(
                "NCBI XML batch record UIDs do not match the requested protein "
                f"UIDs. Missing: {missing_uids[:5]}; "
                f"Unexpected: {unexpected_uids[:5]}."
            )

    return XmlBatchValidationResult(
        root_tag=root_tag,
        record_count=record_count,
        protein_uids=extracted_protein_uids,
    )


def _build_consolidated_root_open_tag(
    *,
    root_tag: str,
    root_attributes: Dict[str, str],
) -> bytes:
    """
    Build the consolidated root open tag exactly as ElementTree would write it.

    The consolidated snapshot hash is part of the artifact contract, so the
    streaming writer must reproduce ElementTree's serialization byte for byte,
    including attribute escaping. Serializing an empty element and converting
    its self-closing suffix keeps that responsibility inside ElementTree.
    """
    empty_root_element = ET.Element(root_tag, root_attributes)
    serialized_empty_root = ET.tostring(empty_root_element, encoding="utf-8")

    if serialized_empty_root.endswith(b" />"):
        return serialized_empty_root[: -len(b" />")] + b">"

    closing_tag = f"</{root_tag}>".encode("utf-8")
    if serialized_empty_root.endswith(closing_tag):
        return serialized_empty_root[: -len(closing_tag)]

    raise RuntimeError(
        "Unable to derive the consolidated XML root open tag for "
        f"{root_tag!r}."
    )


def write_consolidated_xml_document(
    *,
    batch_payload_sources: Sequence[XmlBatchPayloadSource],
    output_file_path: PathLike,
    write_chunk_size: int = 1024 * 1024,
) -> ConsolidatedXmlWriteResult:
    """
    Merge ordered XML batches into one document without holding it in memory.

    Each batch is streamed with ``iterparse`` and each record is re-serialized
    straight into the destination file, so peak memory is bounded by the largest
    single record instead of by the whole consolidated document. The SHA-256 is
    computed from the same bytes as they are written, which removes the separate
    hashing pass over the finished file.

    The output is byte-identical to building one in-memory ElementTree and
    calling ``ET.tostring(root, encoding="utf-8", xml_declaration=True)``.
    """
    if not batch_payload_sources:
        raise ValueError("batch_payload_sources must contain at least one batch.")

    resolved_output_file_path = _as_path(output_file_path)
    resolved_output_file_path.parent.mkdir(parents=True, exist_ok=True)

    consolidated_hasher = hashlib.sha256()
    written_byte_count = 0
    total_record_count = 0
    consolidated_root_tag: Optional[str] = None
    consolidated_root_attributes: Optional[Dict[str, str]] = None

    with tempfile.NamedTemporaryFile(
        mode="wb",
        delete=False,
        dir=resolved_output_file_path.parent,
        prefix=f"{resolved_output_file_path.name}_tmp_",
    ) as temporary_file:
        resolved_temporary_file_path = Path(temporary_file.name)

        def write_consolidated_bytes(payload: bytes) -> None:
            nonlocal written_byte_count

            temporary_file.write(payload)
            consolidated_hasher.update(payload)
            written_byte_count += len(payload)

        try:
            write_consolidated_bytes(_XML_DECLARATION_BYTES)

            for batch_payload_source in batch_payload_sources:
                current_batch_index = batch_payload_source.batch_index

                def observe_root_element(
                    root_tag: str,
                    root_attributes: Dict[str, str],
                ) -> None:
                    nonlocal consolidated_root_tag
                    nonlocal consolidated_root_attributes

                    if consolidated_root_tag is None:
                        consolidated_root_tag = root_tag
                        consolidated_root_attributes = root_attributes
                        write_consolidated_bytes(
                            _build_consolidated_root_open_tag(
                                root_tag=root_tag,
                                root_attributes=root_attributes,
                            )
                        )
                    elif root_tag != consolidated_root_tag:
                        raise RuntimeError(
                            "XML batch root tag mismatch during consolidation. "
                            f"Expected {consolidated_root_tag}, got {root_tag} "
                            f"in batch {current_batch_index}."
                        )

                def count_record_element(record_element: ET.Element) -> None:
                    nonlocal total_record_count

                    total_record_count += 1

                def write_record_element(record_element: ET.Element) -> None:
                    write_consolidated_bytes(
                        ET.tostring(record_element, encoding="utf-8")
                    )

                payload_stream = batch_payload_source.open_payload_stream()
                try:
                    _stream_root_children(
                        payload_stream=payload_stream,
                        on_root_element=observe_root_element,
                        on_record_element=count_record_element,
                        on_serializable_record=write_record_element,
                    )
                except ET.ParseError as error:
                    raise RuntimeError(
                        f"Failed to parse XML batch {current_batch_index}: {error}"
                    ) from error
                finally:
                    payload_stream.close()

            if consolidated_root_tag is None:
                raise RuntimeError("Failed to create consolidated XML root.")

            write_consolidated_bytes(
                f"</{consolidated_root_tag}>".encode("utf-8")
            )
            temporary_file.flush()
        except BaseException:
            temporary_file.close()
            resolved_temporary_file_path.unlink(missing_ok=True)
            raise

    resolved_temporary_file_path.replace(resolved_output_file_path)

    return ConsolidatedXmlWriteResult(
        output_file_path=resolved_output_file_path,
        root_tag=consolidated_root_tag,
        record_count=total_record_count,
        byte_count=written_byte_count,
        sha256=consolidated_hasher.hexdigest(),
    )


def validate_consolidated_xml_file(
    *,
    xml_file_path: PathLike,
    expected_root_tag: str = DEFAULT_CONSOLIDATED_ROOT_TAG,
    expected_record_tag: str = DEFAULT_CONSOLIDATED_RECORD_TAG,
    expected_record_count: Optional[int] = None,
    expected_protein_uids: Optional[Sequence[str]] = None,
    validation_context: str = "",
) -> ConsolidatedXmlValidationResult:
    """
    Validate a saved consolidated XML snapshot in one streaming pass.

    Structure, record tags, record count, and the exact protein UID multiset are
    all checked while the file is read once. Collapsing these checks removes the
    repeated full parses that a consolidated snapshot previously paid on every
    creation and on every reload.
    """
    resolved_xml_file_path = _as_path(xml_file_path)
    context_suffix = f" {validation_context}" if validation_context else ""

    extracted_protein_uids: List[str] = []
    invalid_record_tags: set[str] = set()
    record_count = 0
    collect_protein_uids = expected_protein_uids is not None

    def observe_record(record_element: ET.Element) -> None:
        nonlocal record_count

        record_count += 1
        if record_element.tag != expected_record_tag:
            invalid_record_tags.add(record_element.tag)
            return

        if collect_protein_uids:
            extracted_protein_uids.append(
                extract_protein_uid_from_gbseq_element(
                    gbseq_element=record_element,
                )
            )

    with resolved_xml_file_path.open("rb") as xml_file_handle:
        try:
            root_tag = _stream_root_children(
                payload_stream=xml_file_handle,
                on_root_element=None,
                on_record_element=observe_record,
            )
        except ET.ParseError as error:
            raise RuntimeError(
                f"Saved consolidated XML snapshot is not well-formed: {error}"
            ) from error

    if root_tag != expected_root_tag:
        raise RuntimeError(
            "Saved consolidated XML snapshot root tag mismatch"
            f"{context_suffix}. Expected {expected_root_tag}, got {root_tag}."
        )

    if invalid_record_tags:
        raise RuntimeError(
            "Saved consolidated XML snapshot contains unexpected child tags "
            f"under {expected_root_tag}{context_suffix}: "
            f"{sorted(invalid_record_tags)}."
        )

    if expected_record_count is not None and record_count != expected_record_count:
        raise RuntimeError(
            "Saved consolidated XML snapshot record count mismatch. "
            f"Expected {expected_record_count}, got {record_count}."
        )

    if expected_protein_uids is not None:
        expected_uid_counter = Counter(expected_protein_uids)
        extracted_uid_counter = Counter(extracted_protein_uids)
        if extracted_uid_counter != expected_uid_counter:
            missing_uids = sorted(
                (expected_uid_counter - extracted_uid_counter).elements()
            )
            unexpected_uids = sorted(
                (extracted_uid_counter - expected_uid_counter).elements()
            )
            raise RuntimeError(
                "Saved consolidated XML snapshot record UIDs do not match the "
                "expected protein UIDs. "
                f"Missing: {missing_uids[:5]}; Unexpected: {unexpected_uids[:5]}."
            )

    return ConsolidatedXmlValidationResult(
        root_tag=root_tag,
        record_count=record_count,
        protein_uids=extracted_protein_uids,
    )


def build_legacy_consolidated_xml_payload(
    *,
    batch_payloads: Iterable[bytes],
) -> tuple[bytes, int]:
    """
    Build a consolidated payload with the original in-memory algorithm.

    Retained as the reference implementation that the streaming writer is
    verified against. It is not used by the retrieval workflow because it holds
    the whole consolidated document in memory.
    """
    consolidated_root: Optional[ET.Element] = None
    expected_root_tag: Optional[str] = None
    total_record_count = 0

    for batch_index, batch_payload in enumerate(batch_payloads, start=1):
        try:
            batch_root = ET.fromstring(batch_payload)
        except ET.ParseError as error:
            raise RuntimeError(
                f"Failed to parse XML batch {batch_index}: {error}"
            ) from error

        if consolidated_root is None:
            consolidated_root = ET.Element(batch_root.tag, batch_root.attrib)
            expected_root_tag = batch_root.tag
        elif batch_root.tag != expected_root_tag:
            raise RuntimeError(
                "XML batch root tag mismatch during consolidation. "
                f"Expected {expected_root_tag}, got {batch_root.tag} "
                f"in batch {batch_index}."
            )

        batch_children = list(batch_root)
        for child_element in batch_children:
            consolidated_root.append(child_element)

        total_record_count += len(batch_children)

    if consolidated_root is None:
        raise RuntimeError("Failed to create consolidated XML root.")

    return (
        ET.tostring(consolidated_root, encoding="utf-8", xml_declaration=True),
        total_record_count,
    )
