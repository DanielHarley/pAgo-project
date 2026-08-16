"""
Check that requesting ``rettype=gp`` does not change the EFetch XML payload.

NBK25499 documents ``rettype=gp`` with ``retmode=xml`` as the GBSeq XML
combination for ``db=protein``. The pipeline historically sent ``retmode=xml``
with no ``rettype`` at all, which is not a documented combination even though
NCBI currently resolves it to GBSeq XML. Pinning the documented value removes
that dependency on undocumented behavior, but only if it is a no-op today.

The discriminating test is a paired fetch of the same UIDs in the same run.
Comparing a fresh payload against a months-old snapshot hash cannot distinguish
a format change from ordinary record revision; comparing two payloads fetched
seconds apart can.

Cost: two NCBI EFetch requests.

Usage:
    python scripts/ncbi_rettype_invariance_check.py
    python scripts/ncbi_rettype_invariance_check.py --protein-uid-count 100
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path

from dotenv import find_dotenv, load_dotenv

REPOSITORY_ROOT_DIRECTORY = Path(__file__).resolve().parent.parent
if str(REPOSITORY_ROOT_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT_DIRECTORY))

from src.pago_pipeline.ncbi_api import (  # noqa: E402
    NCBI_PROTEIN_GBSEQ_XML_RETTYPE,
    fetch_ncbi_protein_xml_batches,
    read_ncbi_xml_batch_payload_bytes,
)
from src.pago_pipeline.storage import read_json_file, read_text_lines_from_file  # noqa: E402

DEFAULT_PROTEIN_UIDS_FILE_PATH = (
    Path("data") / "01-raw" / "protein_uid_snapshots" / "latest" / "protein_uids.txt"
)
DEFAULT_XML_MANIFEST_FILE_PATH = (
    Path("data") / "01-raw" / "protein_xml_snapshots" / "latest" / "manifest.json"
)


def parse_command_line_arguments() -> argparse.Namespace:
    argument_parser = argparse.ArgumentParser(
        description=(
            "Fetch the same protein UIDs with and without rettype=gp and "
            "compare the raw payloads."
        ),
    )
    argument_parser.add_argument(
        "--protein-uid-count",
        type=int,
        default=100,
        help="Number of leading UIDs from the frozen snapshot to fetch.",
    )
    argument_parser.add_argument(
        "--protein-uids-file",
        type=Path,
        default=REPOSITORY_ROOT_DIRECTORY / DEFAULT_PROTEIN_UIDS_FILE_PATH,
        help="Frozen protein_uids.txt used as the request input.",
    )
    argument_parser.add_argument(
        "--xml-manifest-file",
        type=Path,
        default=REPOSITORY_ROOT_DIRECTORY / DEFAULT_XML_MANIFEST_FILE_PATH,
        help=(
            "Frozen XML snapshot manifest used for the secondary comparison "
            "against the first recorded batch hash."
        ),
    )
    return argument_parser.parse_args()


def fetch_single_batch_payload(
    *,
    ncbi_email: str,
    ncbi_api_key: str | None,
    protein_uids: list[str],
    rettype: str | None,
) -> bytes:
    fetch_result = fetch_ncbi_protein_xml_batches(
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        protein_uids=protein_uids,
        batch_size=len(protein_uids),
        rettype=rettype,
        spill_batch_payloads_to_disk=False,
        max_retry_attempts=3,
    )
    return read_ncbi_xml_batch_payload_bytes(
        xml_batch=fetch_result.xml_batches[0],
    )


def main() -> int:
    command_line_arguments = parse_command_line_arguments()

    dotenv_file_path = find_dotenv(usecwd=False)
    if dotenv_file_path:
        load_dotenv(dotenv_path=dotenv_file_path, override=True)

    ncbi_email = os.getenv("NCBI_EMAIL")
    ncbi_api_key = os.getenv("NCBI_API_KEY")
    if not ncbi_email:
        print(
            "NCBI_EMAIL is not configured. Add it to .env before running this "
            "check."
        )
        return 2

    protein_uids = read_text_lines_from_file(
        input_file_path=command_line_arguments.protein_uids_file,
    )[: command_line_arguments.protein_uid_count]

    if not protein_uids:
        print(
            "No protein UIDs were read from "
            f"{command_line_arguments.protein_uids_file}."
        )
        return 2

    print(
        f"Comparing EFetch payloads for {len(protein_uids)} protein UIDs "
        f"({protein_uids[0]}..{protein_uids[-1]})."
    )
    print("API key configured:", bool(ncbi_api_key))

    payload_without_rettype = fetch_single_batch_payload(
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        protein_uids=protein_uids,
        rettype=None,
    )
    payload_with_gp_rettype = fetch_single_batch_payload(
        ncbi_email=ncbi_email,
        ncbi_api_key=ncbi_api_key,
        protein_uids=protein_uids,
        rettype=NCBI_PROTEIN_GBSEQ_XML_RETTYPE,
    )

    sha256_without_rettype = hashlib.sha256(payload_without_rettype).hexdigest()
    sha256_with_gp_rettype = hashlib.sha256(payload_with_gp_rettype).hexdigest()
    payloads_are_identical = sha256_without_rettype == sha256_with_gp_rettype

    print()
    print("=== PRIMARY CRITERION: rettype invariance ===")
    print(f"A  retmode=xml, no rettype : {sha256_without_rettype}")
    print(f"   bytes                   : {len(payload_without_rettype)}")
    print(f"B  retmode=xml, rettype=gp : {sha256_with_gp_rettype}")
    print(f"   bytes                   : {len(payload_with_gp_rettype)}")
    print(f"sha256(A) == sha256(B)     : {payloads_are_identical}")

    print()
    print("=== SECONDARY SIGNAL: comparison against the frozen snapshot ===")
    frozen_batch_sha256 = None
    if command_line_arguments.xml_manifest_file.exists():
        frozen_manifest = read_json_file(
            input_file_path=command_line_arguments.xml_manifest_file,
        )
        frozen_batches = frozen_manifest.get("batches") or []
        if frozen_batches:
            frozen_batch_sha256 = frozen_batches[0].get("xml_payload_sha256")

    if frozen_batch_sha256 is None:
        print("No frozen batch hash was available for comparison.")
    else:
        print(f"Frozen batches[0].xml_payload_sha256: {frozen_batch_sha256}")
        print(f"A matches frozen batch: {sha256_without_rettype == frozen_batch_sha256}")
        print(f"B matches frozen batch: {sha256_with_gp_rettype == frozen_batch_sha256}")
        if payloads_are_identical and sha256_without_rettype != frozen_batch_sha256:
            print(
                "A and B agree but differ from the frozen batch: this indicates "
                "record revision since the snapshot was taken, not a format "
                "change. Record it and continue."
            )

    print()
    if payloads_are_identical:
        print("PASS: rettype=gp is a no-op today and can be adopted safely.")
        return 0

    print(
        "FAIL: rettype=gp changed the payload. Do not adopt it without "
        "re-examining the downstream GBSeq schema assumptions."
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
