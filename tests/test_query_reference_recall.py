from __future__ import annotations

import re
import unittest
from pathlib import Path

import pandas as pd

from src.pago_pipeline.query_reference_recall import (
    ReferenceMatchMethod,
    build_matching_strategy_sha256,
    compute_query_reference_recall,
    normalize_protein_sequence,
    protein_sequence_sha256,
)

_FIXTURE_CSV = (
    Path(__file__).resolve().parent / "fixtures" / "query_recall_reference_set.csv"
)

# Pinned so any change to the matching methodology forces this test (and the
# snapshot-reuse invariant) to be revisited deliberately.
_EXPECTED_MATCHING_STRATEGY_SHA256 = (
    "3460b048fc6de363ddf9282c2943a44c284e51d2c48092ca14426600b2871a08"
)

_RSAGO_ACCESSION = "ABP72561.1"
_RSAGO_SEQUENCE_SHA256 = (
    "cbdb6bb64718c9e8ca78a34ac8445eff1556cb87b5ad687026373ed401c5fb36"
)

_HEX64 = re.compile(r"^[0-9a-f]{64}$")

# Reference sequences used by the synthetic logic tests. Only their sha256 goes
# into the reference frame - the retrieved frame carries the raw sequence.
_REF_SEQUENCES = {
    "WP_001.1": "MAAAAKKKK",
    "WP_002.2": "MBBBBKKKK",
    "WP_003.1": "MCCCCKKKK",
    "WP_004.1": "MDDDDKKKK",
    "WP_005.1": "MEEEEKKKK",
}


def _reference_dataframe() -> pd.DataFrame:
    rows = [
        ("WP_001.1", "PAGO", "LONG_A", "EXPERIMENTAL"),
        ("WP_002.2", "PAGO", "LONG_A", "EXPERIMENTAL"),
        ("WP_003.1", "PAGO", "LONG_B", "EXPERIMENTAL"),
        ("WP_004.1", "PAGO", "SHORT", "EXPERIMENTAL"),
        ("WP_005.1", "PIWI_RE", "UNRESOLVED", "LITERATURE_PHYLOGENETIC"),
    ]
    return pd.DataFrame(
        [
            {
                "accession": accession,
                "ago_family": ago_family,
                "clade": clade,
                "sequence_sha256": protein_sequence_sha256(
                    _REF_SEQUENCES[accession]
                ),
                "reference_label_source": "src",
                "reference_label_evidence": evidence,
            }
            for accession, ago_family, clade, evidence in rows
        ]
    )


def _detail_row(result, accession: str) -> pd.Series:
    frame = result.detail_dataframe
    return frame[frame["accession"] == accession].iloc[0]


class QueryReferenceRecallLogicTests(unittest.TestCase):
    def test_matching_strategy_sha256_is_pinned(self) -> None:
        self.assertEqual(
            build_matching_strategy_sha256(), _EXPECTED_MATCHING_STRATEGY_SHA256
        )

    def test_sequence_normalization_strips_whitespace_and_uppercases(self) -> None:
        self.assertEqual(normalize_protein_sequence("  ma k\n l "), "MAKL")
        self.assertEqual(
            protein_sequence_sha256("MA KL"), protein_sequence_sha256("makl")
        )

    def test_exact_accession_version_match(self) -> None:
        retrieved = pd.DataFrame(
            {
                "protein_uid": ["11"],
                "gbseq__accession_version": ["WP_001.1"],
                "gbseq__sequence": ["TOTALLYDIFFERENT"],
            }
        )
        result = compute_query_reference_recall(
            reference_dataframe=_reference_dataframe(),
            retrieved_metadata_dataframe=retrieved,
        )
        row = _detail_row(result, "WP_001.1")
        self.assertEqual(
            row["match_method"], ReferenceMatchMethod.EXACT_ACCESSION_VERSION.value
        )
        self.assertEqual(row["matched_accession"], "WP_001.1")
        self.assertEqual(row["matched_protein_uid"], "11")
        self.assertTrue(row["recovered"])
        self.assertTrue(row["recovered_exact_accession"])

    def test_same_base_accession_different_version(self) -> None:
        retrieved = pd.DataFrame(
            {
                "protein_uid": ["12"],
                "gbseq__accession_version": ["WP_002.9"],
                "gbseq__sequence": ["SOMETHING"],
            }
        )
        result = compute_query_reference_recall(
            reference_dataframe=_reference_dataframe(),
            retrieved_metadata_dataframe=retrieved,
        )
        row = _detail_row(result, "WP_002.2")
        self.assertEqual(
            row["match_method"], ReferenceMatchMethod.SAME_BASE_ACCESSION.value
        )
        self.assertEqual(row["matched_accession"], "WP_002.9")
        self.assertTrue(row["recovered"])
        self.assertTrue(row["recovered_exact_accession"])

    def test_different_accession_identical_sequence_is_equivalent_only(self) -> None:
        # RsAgo-shaped case: the reference sequence is retrieved under a wholly
        # unrelated accession. Recovered for the retrieval-equivalent reading,
        # NOT for the strict-accession reading.
        retrieved = pd.DataFrame(
            {
                "protein_uid": ["77"],
                "gbseq__accession_version": ["A4WYU7.1"],
                "gbseq__sequence": [_REF_SEQUENCES["WP_003.1"]],
            }
        )
        result = compute_query_reference_recall(
            reference_dataframe=_reference_dataframe(),
            retrieved_metadata_dataframe=retrieved,
        )
        row = _detail_row(result, "WP_003.1")
        self.assertEqual(
            row["match_method"], ReferenceMatchMethod.SEQUENCE_SHA256.value
        )
        self.assertEqual(row["matched_accession"], "A4WYU7.1")
        self.assertEqual(row["matched_protein_uid"], "77")
        self.assertEqual(row["sequence_match_count"], 1)
        self.assertTrue(row["recovered"])
        self.assertFalse(row["recovered_exact_accession"])

        self.assertEqual(
            result.stratum_exact_recall["long_b_reference_recall"], 0.0
        )
        self.assertEqual(
            result.stratum_equivalent_recall["long_b_reference_recall"], 1.0
        )
        self.assertEqual(result.exact_recovered_count, 0)
        self.assertEqual(result.equivalent_recovered_count, 1)

    def test_different_sequence_is_not_recovered(self) -> None:
        retrieved = pd.DataFrame(
            {
                "protein_uid": ["99"],
                "gbseq__accession_version": ["WP_888888.1"],
                "gbseq__sequence": ["MZZZZQQQQ"],
            }
        )
        result = compute_query_reference_recall(
            reference_dataframe=_reference_dataframe(),
            retrieved_metadata_dataframe=retrieved,
        )
        for accession in _REF_SEQUENCES:
            row = _detail_row(result, accession)
            self.assertEqual(
                row["match_method"], ReferenceMatchMethod.NONE.value
            )
            self.assertFalse(row["recovered"])
        self.assertEqual(result.exact_recovered_count, 0)
        self.assertEqual(result.equivalent_recovered_count, 0)

    def test_exact_accession_takes_priority_over_sequence_hash(self) -> None:
        # WP_001.1 is retrieved verbatim but with a mutated sequence; a second,
        # unrelated record carries the true WP_001.1 sequence. The exact
        # accession hit must win, and it must be scored as strict-accession.
        retrieved = pd.DataFrame(
            {
                "protein_uid": ["1", "2"],
                "gbseq__accession_version": ["WP_001.1", "DECOY_1.1"],
                "gbseq__sequence": [
                    "MUTATEDSEQUENCE",
                    _REF_SEQUENCES["WP_001.1"],
                ],
            }
        )
        result = compute_query_reference_recall(
            reference_dataframe=_reference_dataframe(),
            retrieved_metadata_dataframe=retrieved,
        )
        row = _detail_row(result, "WP_001.1")
        self.assertEqual(
            row["match_method"], ReferenceMatchMethod.EXACT_ACCESSION_VERSION.value
        )
        self.assertEqual(row["matched_accession"], "WP_001.1")
        self.assertTrue(row["recovered_exact_accession"])

    def test_same_base_accession_takes_priority_over_sequence_hash(self) -> None:
        retrieved = pd.DataFrame(
            {
                "protein_uid": ["1", "2"],
                "gbseq__accession_version": ["WP_002.5", "DECOY_2.1"],
                "gbseq__sequence": [
                    "MUTATEDSEQUENCE",
                    _REF_SEQUENCES["WP_002.2"],
                ],
            }
        )
        result = compute_query_reference_recall(
            reference_dataframe=_reference_dataframe(),
            retrieved_metadata_dataframe=retrieved,
        )
        row = _detail_row(result, "WP_002.2")
        self.assertEqual(
            row["match_method"], ReferenceMatchMethod.SAME_BASE_ACCESSION.value
        )
        self.assertEqual(row["matched_accession"], "WP_002.5")

    def test_multiple_accessions_same_hash_resolved_deterministically(self) -> None:
        retrieved = pd.DataFrame(
            {
                "protein_uid": ["30", "31", "32"],
                "gbseq__accession_version": ["ZZ_9.1", "AA_1.1", "MM_5.1"],
                "gbseq__sequence": [_REF_SEQUENCES["WP_005.1"]] * 3,
            }
        )
        first = compute_query_reference_recall(
            reference_dataframe=_reference_dataframe(),
            retrieved_metadata_dataframe=retrieved,
        )
        second = compute_query_reference_recall(
            reference_dataframe=_reference_dataframe(),
            retrieved_metadata_dataframe=retrieved.sample(frac=1, random_state=0),
        )
        for result in (first, second):
            row = _detail_row(result, "WP_005.1")
            self.assertEqual(
                row["match_method"], ReferenceMatchMethod.SEQUENCE_SHA256.value
            )
            # Lexicographic minimum of the colliding accessions.
            self.assertEqual(row["matched_accession"], "AA_1.1")
            self.assertEqual(row["matched_protein_uid"], "31")
            self.assertEqual(row["sequence_match_count"], 3)

    def test_stratified_recall_and_both_readings(self) -> None:
        retrieved = pd.DataFrame(
            {
                "protein_uid": ["11", "12", "13"],
                # WP_001.1 exact; WP_002.9 same base as WP_002.2; a decoy that
                # carries WP_003.1's sequence.
                "gbseq__accession_version": ["WP_001.1", "WP_002.9", "DECOY.1"],
                "gbseq__sequence": [
                    _REF_SEQUENCES["WP_001.1"],
                    _REF_SEQUENCES["WP_002.2"],
                    _REF_SEQUENCES["WP_003.1"],
                ],
            }
        )
        result = compute_query_reference_recall(
            reference_dataframe=_reference_dataframe(),
            retrieved_metadata_dataframe=retrieved,
        )
        self.assertEqual(result.reference_count, 5)
        self.assertEqual(result.exact_recovered_count, 2)
        self.assertEqual(result.equivalent_recovered_count, 3)
        self.assertAlmostEqual(
            result.stratum_exact_recall["overall_reference_recall"], 2 / 5
        )
        self.assertAlmostEqual(
            result.stratum_equivalent_recall["overall_reference_recall"], 3 / 5
        )
        self.assertAlmostEqual(
            result.stratum_exact_recall["long_a_reference_recall"], 1.0
        )
        self.assertAlmostEqual(
            result.stratum_equivalent_recall["long_b_reference_recall"], 1.0
        )
        self.assertAlmostEqual(
            result.stratum_exact_recall["long_b_reference_recall"], 0.0
        )
        for metric_name in (
            "overall_reference_recall",
            "long_a_reference_recall",
            "long_b_reference_recall",
            "short_reference_recall",
            "piwi_re_reference_recall",
        ):
            self.assertEqual(
                result.stratum_recall_status[metric_name], "EVALUABLE"
            )

    def test_empty_stratum_is_not_evaluable_not_zero(self) -> None:
        reference_dataframe = _reference_dataframe()
        reference_dataframe = reference_dataframe[
            reference_dataframe["clade"] != "SHORT"
        ]
        result = compute_query_reference_recall(
            reference_dataframe=reference_dataframe,
            retrieved_metadata_dataframe=pd.DataFrame(
                {
                    "protein_uid": ["1"],
                    "gbseq__accession_version": ["WP_001.1"],
                    "gbseq__sequence": ["X"],
                }
            ),
        )
        self.assertIsNone(
            result.stratum_exact_recall["short_reference_recall"]
        )
        self.assertIsNone(
            result.stratum_equivalent_recall["short_reference_recall"]
        )
        self.assertEqual(
            result.stratum_recall_status["short_reference_recall"],
            "NOT_EVALUABLE",
        )
        short_row = result.summary_dataframe[
            result.summary_dataframe["stratum"] == "SHORT"
        ].iloc[0]
        self.assertEqual(short_row["reference_count"], 0)
        self.assertTrue(pd.isna(short_row["exact_accession_recall"]))
        self.assertTrue(pd.isna(short_row["retrieval_equivalent_recall"]))
        self.assertEqual(short_row["recall_status"], "NOT_EVALUABLE")
        # A present-but-unrecovered stratum is still EVALUABLE at 0.0.
        self.assertEqual(
            result.stratum_recall_status["long_b_reference_recall"], "EVALUABLE"
        )
        self.assertEqual(
            result.stratum_exact_recall["long_b_reference_recall"], 0.0
        )

    def test_missing_reference_columns_raise(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "missing required columns"):
            compute_query_reference_recall(
                reference_dataframe=pd.DataFrame([{"accession": "x"}]),
                retrieved_metadata_dataframe=pd.DataFrame(
                    {"gbseq__accession_version": ["x"]}
                ),
            )


class CommittedReferenceSetTests(unittest.TestCase):
    def setUp(self) -> None:
        self.reference_dataframe = pd.read_csv(_FIXTURE_CSV, dtype=str).fillna("")

    def test_committed_reference_set_is_well_formed(self) -> None:
        for column_name in (
            "accession",
            "ago_family",
            "clade",
            "sequence_sha256",
            "sequence_length",
            "reference_label_source",
            "reference_label_evidence",
            "verification_status",
        ):
            self.assertIn(column_name, self.reference_dataframe.columns)
        self.assertGreaterEqual(len(self.reference_dataframe), 21)
        self.assertTrue(
            (self.reference_dataframe["accession"].str.len() > 0).all()
        )
        self.assertTrue(
            self.reference_dataframe["accession"].str.contains(r"\.\d+$").all()
        )
        clades = set(self.reference_dataframe["clade"].str.upper())
        # PIWI_RE is an ago_family, not a pAgo clade, so it is not an allowed
        # value in the clade column.
        self.assertTrue(
            clades <= {"LONG_A", "LONG_B", "SHORT", "UNRESOLVED", "NA", ""}
        )
        for required_clade in ("LONG_A", "LONG_B", "SHORT"):
            self.assertIn(required_clade, clades)

        piwi_re_rows = self.reference_dataframe[
            self.reference_dataframe["ago_family"].str.upper() == "PIWI_RE"
        ]
        self.assertGreaterEqual(len(piwi_re_rows), 1)
        # Every PIWI-RE row leaves the pAgo clade unresolved.
        self.assertEqual(
            set(piwi_re_rows["clade"].str.upper()), {"UNRESOLVED"}
        )

        self.assertTrue(
            set(self.reference_dataframe["reference_label_evidence"].str.upper())
            <= {
                "EXPERIMENTAL",
                "LITERATURE_PHYLOGENETIC",
                "CURATED_COMPUTATIONAL",
                "DATABASE_ANNOTATION",
            }
        )
        self.assertTrue(
            set(self.reference_dataframe["verification_status"].str.lower())
            <= {"verified", "provisional"}
        )

    def test_stratum_sizes_match_curation(self) -> None:
        clade = self.reference_dataframe["clade"].str.upper()
        family = self.reference_dataframe["ago_family"].str.upper()
        self.assertEqual(int((clade == "LONG_A").sum()), 8)
        self.assertEqual(int((clade == "LONG_B").sum()), 2)
        self.assertEqual(int((clade == "SHORT").sum()), 4)
        self.assertEqual(int((family == "PIWI_RE").sum()), 7)

    def test_reference_sequence_hashes_are_well_formed(self) -> None:
        digests = self.reference_dataframe["sequence_sha256"].str.strip()
        self.assertEqual(len(digests), 21)
        self.assertTrue((digests != "").all())
        self.assertTrue(digests.map(lambda value: bool(_HEX64.match(value))).all())
        # Every reference sequence is distinct in this panel.
        self.assertEqual(digests.nunique(), len(digests))

        lengths = self.reference_dataframe["sequence_length"].str.strip()
        self.assertTrue((lengths.str.match(r"^\d+$")).all())
        self.assertTrue((lengths.astype(int) > 0).all())

    def test_rsago_row_carries_the_investigated_sequence_hash(self) -> None:
        rsago = self.reference_dataframe[
            self.reference_dataframe["accession"] == _RSAGO_ACCESSION
        ]
        self.assertEqual(len(rsago), 1)
        rsago_row = rsago.iloc[0]
        self.assertEqual(rsago_row["clade"].upper(), "LONG_B")
        self.assertEqual(
            rsago_row["sequence_sha256"].strip().lower(), _RSAGO_SEQUENCE_SHA256
        )
        self.assertEqual(rsago_row["sequence_length"].strip(), "777")

    def test_rsago_is_recovered_by_sequence_identity_under_alias_accession(
        self,
    ) -> None:
        # The RsAgo scenario from the first real Phase A run: ABP72561.1 is not
        # among the retrieved accessions, but its byte-identical sequence is,
        # under the Swiss-Prot alias A4WYU7.1 (same NCBI Identical Protein Group).
        # The committed sha256 is derived offline; here we drive an equivalent
        # collision with a synthetic sequence so the test stays hermetic.
        synthetic_rsago_sequence = "M" + "RSAGO" * 40
        rsago_row = self.reference_dataframe[
            self.reference_dataframe["accession"] == _RSAGO_ACCESSION
        ].copy()
        rsago_row["sequence_sha256"] = protein_sequence_sha256(
            synthetic_rsago_sequence
        )
        retrieved = pd.DataFrame(
            {
                "protein_uid": ["2500461169"],
                "gbseq__accession_version": ["A4WYU7.1"],
                "gbseq__sequence": [synthetic_rsago_sequence],
            }
        )
        result = compute_query_reference_recall(
            reference_dataframe=rsago_row,
            retrieved_metadata_dataframe=retrieved,
        )
        detail = result.detail_dataframe.iloc[0]
        self.assertEqual(detail["accession"], _RSAGO_ACCESSION)
        self.assertEqual(
            detail["match_method"], ReferenceMatchMethod.SEQUENCE_SHA256.value
        )
        self.assertEqual(detail["matched_accession"], "A4WYU7.1")
        self.assertEqual(detail["matched_protein_uid"], "2500461169")
        self.assertTrue(detail["recovered"])
        self.assertFalse(detail["recovered_exact_accession"])
        self.assertEqual(
            result.stratum_exact_recall["long_b_reference_recall"], 0.0
        )
        self.assertEqual(
            result.stratum_equivalent_recall["long_b_reference_recall"], 1.0
        )


if __name__ == "__main__":
    unittest.main()
