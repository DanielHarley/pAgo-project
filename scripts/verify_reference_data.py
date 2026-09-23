from __future__ import annotations

import argparse
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Per-layer validators are imported lazily inside each verifier so that a single
# approved scope (for example ``--scope pfam``) runs without importing any
# reference layer that has not been reintegrated to `master` yet. Only the
# layers actually present in this repository are registered in VERIFIERS
# below; a future reference-layer PR adds its own verifier function and
# registry entry rather than restructuring this dispatcher.


@dataclass(frozen=True)
class VerificationSummary:
    reference_layer: str
    details: tuple[str, ...]


Verifier = Callable[[], VerificationSummary]


def verify_pfam_reference_data() -> VerificationSummary:
    from src.pago_pipeline.pfam_hmm_bundle import validate_pfam_hmm_reference_data

    result = validate_pfam_hmm_reference_data()
    return VerificationSummary(
        reference_layer="pfam_hmm",
        details=(
            f"source={result.registry.source_database} {result.registry.source_release}",
            f"models={len(result.validated_models)}",
            f"registry_sha256={result.registry.registry_sha256}",
            f"pyhmmer={result.pyhmmer_version}",
        ),
    )


VERIFIERS: dict[str, Verifier] = {
    "pfam": verify_pfam_reference_data,
}


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate committed pAgo reference-layer locks, SHA-256 hashes, "
            "and inventories without modifying any artifact or touching the "
            "network."
        )
    )
    parser.add_argument(
        "--scope",
        action="append",
        choices=tuple(VERIFIERS),
        help=(
            "Reference layer to validate. Repeat for multiple layers. "
            "Omission validates every registered layer."
        ),
    )
    return parser


def verify_selected_reference_layers(
    *,
    scopes: Sequence[str] | None = None,
) -> tuple[VerificationSummary, ...]:
    selected_scopes = tuple(scopes) if scopes else tuple(VERIFIERS)
    if len(set(selected_scopes)) != len(selected_scopes):
        raise ValueError("Each verification scope may be selected only once.")
    return tuple(VERIFIERS[scope]() for scope in selected_scopes)


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_argument_parser().parse_args(argv)
    summaries = verify_selected_reference_layers(scopes=arguments.scope)
    for summary in summaries:
        print(f"{summary.reference_layer}: OK")
        for detail in summary.details:
            print(f"  {detail}")
    print(f"Verified reference layers: {len(summaries)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
