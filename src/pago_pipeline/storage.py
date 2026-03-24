from __future__ import annotations

import hashlib
import json
import tempfile
from pathlib import Path
from typing import Any, Iterable, List, Union

PathLike = Union[str, Path]


def _normalize_text_lines(
    *,
    text_lines: Iterable[str],
    deduplicate_lines_preserving_order: bool = True,
    sort_lines: bool = False,
) -> List[str]:
    """
    Normalize an iterable of text lines by:
    - casting to str
    - stripping surrounding whitespace
    - removing empty lines
    - optionally deduplicating while preserving first occurrence order
    - optionally sorting canonically
    """
    normalized_nonempty_lines = [
        str(item).strip()
        for item in text_lines
        if str(item).strip()
    ]

    if deduplicate_lines_preserving_order:
        seen = set()
        normalized_nonempty_lines = [
            item
            for item in normalized_nonempty_lines
            if not (item in seen or seen.add(item))
        ]

    if sort_lines:
        normalized_nonempty_lines = sorted(normalized_nonempty_lines)

    return normalized_nonempty_lines


def write_text_lines_to_file(
    *,
    text_lines: Iterable[str],
    output_file_path: PathLike,
    deduplicate_lines_preserving_order: bool = True,
    sort_lines: bool = False,
    encoding: str = "utf-8",
    line_ending: str = "\n",
) -> Path:
    """
    Persist iterable[str] as a text file with one entry per line.

    Guarantees:
    - parent directory creation
    - optional deduplication (stable order)
    - optional sorting
    - atomic write (temp -> replace) to avoid partial artifacts
    """
    output_file_path = Path(output_file_path)
    output_file_path.parent.mkdir(parents=True, exist_ok=True)

    normalized_nonempty_lines = _normalize_text_lines(
        text_lines=text_lines,
        deduplicate_lines_preserving_order=deduplicate_lines_preserving_order,
        sort_lines=sort_lines,
    )

    with tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        dir=output_file_path.parent,
        encoding=encoding,
        newline="\n",
    ) as temporary_file:
        temporary_file_path = Path(temporary_file.name)

        for line in normalized_nonempty_lines:
            temporary_file.write(line + line_ending)

    temporary_file_path.replace(output_file_path)
    return output_file_path


def read_text_lines_from_file(
    *,
    input_file_path: PathLike,
    encoding: str = "utf-8",
    strip_lines: bool = True,
    skip_empty_lines: bool = True,
) -> List[str]:
    """
    Read a text file as a list of lines.

    Useful when you want to reuse a frozen snapshot instead of calling NCBI again.
    """
    input_file_path = Path(input_file_path)

    with input_file_path.open("r", encoding=encoding, newline=None) as file_handle:
        lines = file_handle.readlines()

    if strip_lines:
        lines = [line.strip() for line in lines]
    else:
        lines = [line.rstrip("\r\n") for line in lines]

    if skip_empty_lines:
        lines = [line for line in lines if line]

    return lines


def write_json_atomic(
    *,
    payload: Any,
    output_file_path: PathLike,
    encoding: str = "utf-8",
    indent: int = 2,
    sort_keys: bool = True,
    ensure_ascii: bool = False,
) -> Path:
    """
    Atomically write a JSON file.

    Appropriate for manifests, metadata, provenance, query parameters, hashes, etc.
    """
    output_file_path = Path(output_file_path)
    output_file_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        dir=output_file_path.parent,
        encoding=encoding,
        newline="\n",
    ) as temporary_file:
        temporary_file_path = Path(temporary_file.name)

        json.dump(
            payload,
            temporary_file,
            indent=indent,
            sort_keys=sort_keys,
            ensure_ascii=ensure_ascii,
        )
        temporary_file.write("\n")

    temporary_file_path.replace(output_file_path)
    return output_file_path


def read_json_file(
    *,
    input_file_path: PathLike,
    encoding: str = "utf-8",
) -> Any:
    """
    Read and deserialize a JSON file.
    """
    input_file_path = Path(input_file_path)

    with input_file_path.open("r", encoding=encoding) as file_handle:
        return json.load(file_handle)


def sha256_of_lines(
    *,
    text_lines: Iterable[str],
    deduplicate_lines_preserving_order: bool = True,
    sort_lines: bool = False,
    encoding: str = "utf-8",
    line_ending: str = "\n",
) -> str:
    """
    Compute SHA-256 from a normalized sequence of text lines.

    Important:
    The hash depends on the normalization policy. For reproducible dataset identity,
    use the same deduplication/sorting policy here that you use when saving the file.

    Example for canonical dataset identity:
        deduplicate_lines_preserving_order=True
        sort_lines=True
    """
    normalized_nonempty_lines = _normalize_text_lines(
        text_lines=text_lines,
        deduplicate_lines_preserving_order=deduplicate_lines_preserving_order,
        sort_lines=sort_lines,
    )

    normalized_text = line_ending.join(normalized_nonempty_lines)
    return hashlib.sha256(normalized_text.encode(encoding)).hexdigest()


def sha256_of_file(
    *,
    input_file_path: PathLike,
    chunk_size: int = 1024 * 1024,
) -> str:
    """
    Compute SHA-256 from the raw bytes of an existing file.
    """
    input_file_path = Path(input_file_path)
    hasher = hashlib.sha256()

    with input_file_path.open("rb") as file_handle:
        while True:
            chunk = file_handle.read(chunk_size)
            if not chunk:
                break
            hasher.update(chunk)

    return hasher.hexdigest()


def save_ncbi_protein_ids_as_txt(
    *,
    ncbi_protein_id_list: List[str],
    output_txt_file_path: PathLike,
    deduplicate_ids: bool = True,
    sort_ids: bool = False,
) -> Path:
    """
    Domain wrapper: save NCBI protein IDs (one per line) as .txt artifact.
    """
    return write_text_lines_to_file(
        text_lines=ncbi_protein_id_list,
        output_file_path=output_txt_file_path,
        deduplicate_lines_preserving_order=deduplicate_ids,
        sort_lines=sort_ids,
    )