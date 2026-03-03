from __future__ import annotations
from pathlib import Path
from typing import Iterable, Union, List
import tempfile


PathLike = Union[str, Path]


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

    normalized_nonempty_lines = [
        str(item).strip()
        for item in text_lines
        if str(item).strip()
    ]

    if deduplicate_lines_preserving_order:
        seen = set()
        normalized_nonempty_lines = [
            item for item in normalized_nonempty_lines
            if not (item in seen or seen.add(item))
        ]

    if sort_lines:
        normalized_nonempty_lines = sorted(normalized_nonempty_lines)

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