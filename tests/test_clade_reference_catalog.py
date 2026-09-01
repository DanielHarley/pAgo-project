from __future__ import annotations

import unittest

from src.pago_pipeline.clade_reference_catalog import (
    AGO_TYPE_TO_CLADE,
    MIDPIWI_COORDINATE_CONVENTION,
    coordinate_convention_evidence,
    parse_coordinate,
    parse_table_s1_rows,
    prove_coordinate_convention,
    slice_midpiwi,
)

# A tiny synthetic protein whose MID and PIWI residues are individually
# identifiable, so an off-by-one in the slice shows up immediately.
#   positions 1..4  : header  "HEAD"
#   positions 5..12 : MID     "mmmmmmmm"
#   positions 13..20: PIWI    "pppppppp"
#   positions 21..24: tail    "TAIL"
SYNTHETIC_SEQUENCE = "HEAD" + "M" * 8 + "P" * 8 + "TAIL"
MID_START, MID_END = 5, 12
PIWI_START, PIWI_END = 13, 20


class SliceMidpiwiTests(unittest.TestCase):
    def test_one_based_inclusive_slice_is_exact(self) -> None:
        region = slice_midpiwi(
            SYNTHETIC_SEQUENCE, mid_start=MID_START, piwi_end=PIWI_END
        )
        self.assertEqual(region, "M" * 8 + "P" * 8)
        self.assertEqual(len(region), PIWI_END - MID_START + 1)

    def test_off_by_one_low_would_capture_a_head_residue(self) -> None:
        # mid_start - 1 (a classic 1-vs-0 base confusion) must NOT be what the
        # function does: it would prepend the last HEAD residue.
        wrong = SYNTHETIC_SEQUENCE[(MID_START - 1) - 1:PIWI_END]
        self.assertTrue(wrong.startswith("D"))  # last residue of "HEAD"
        self.assertNotEqual(
            slice_midpiwi(SYNTHETIC_SEQUENCE, mid_start=MID_START, piwi_end=PIWI_END),
            wrong,
        )

    def test_off_by_one_high_would_drop_the_last_piwi_residue(self) -> None:
        # Treating piwi_end as an exclusive bound drops one PIWI residue.
        wrong = SYNTHETIC_SEQUENCE[MID_START - 1:PIWI_END - 1]
        self.assertEqual(len(wrong), PIWI_END - MID_START)
        self.assertNotEqual(
            slice_midpiwi(SYNTHETIC_SEQUENCE, mid_start=MID_START, piwi_end=PIWI_END),
            wrong,
        )

    def test_piwi_running_to_the_c_terminus_is_valid(self) -> None:
        region = slice_midpiwi(SYNTHETIC_SEQUENCE, mid_start=MID_START,
                               piwi_end=len(SYNTHETIC_SEQUENCE))
        self.assertTrue(region.endswith("TAIL"))

    def test_rejects_coordinates_outside_the_sequence(self) -> None:
        with self.assertRaises(ValueError):
            slice_midpiwi(SYNTHETIC_SEQUENCE, mid_start=0, piwi_end=10)
        with self.assertRaises(ValueError):
            slice_midpiwi(SYNTHETIC_SEQUENCE, mid_start=5, piwi_end=999)
        with self.assertRaises(ValueError):
            slice_midpiwi(SYNTHETIC_SEQUENCE, mid_start=15, piwi_end=10)


class CoordinateConventionTests(unittest.TestCase):
    def _record(self, *, length, paz, mid, piwi):
        return {
            "source_length": length,
            "source_paz_start": paz[0], "source_paz_end": paz[1],
            "source_mid_start": mid[0], "source_mid_end": mid[1],
            "source_piwi_start": piwi[0], "source_piwi_end": piwi[1],
        }

    def test_one_based_inclusive_contiguous_records_prove_the_convention(self) -> None:
        records = [
            self._record(length=706, paz=(165, 264), mid=(306, 485), piwi=(486, 692)),
            self._record(length=517, paz=(None, None), mid=(37, 276), piwi=(277, 503)),
            self._record(length=685, paz=(175, 270), mid=(342, 462), piwi=(463, 685)),
        ]
        evidence = prove_coordinate_convention(records)
        self.assertEqual(evidence["convention"], MIDPIWI_COORDINATE_CONVENTION)
        self.assertEqual(evidence["piwi_end_gt_length"], 0)
        self.assertEqual(evidence["mid_end_plus_one_eq_piwi_start"], 3)
        self.assertEqual(evidence["mid_end_eq_piwi_start"], 0)

    def test_a_zero_based_half_open_signature_is_rejected(self) -> None:
        records = [
            self._record(length=706, paz=(165, 305), mid=(305, 485), piwi=(485, 692)),
        ]
        with self.assertRaisesRegex(RuntimeError, "1-based"):
            prove_coordinate_convention(records)

    def test_a_piwi_end_past_the_protein_is_rejected(self) -> None:
        records = [
            self._record(length=400, paz=(10, 90), mid=(100, 200), piwi=(201, 401)),
        ]
        with self.assertRaisesRegex(RuntimeError, "1-based"):
            prove_coordinate_convention(records)


class ParseTests(unittest.TestCase):
    def test_parse_coordinate_handles_range_and_null(self) -> None:
        self.assertEqual(parse_coordinate("306-485"), (306, 485))
        self.assertEqual(parse_coordinate("NULL"), (None, None))
        self.assertEqual(parse_coordinate(""), (None, None))

    def test_parse_table_s1_rows_maps_ago_type_to_clade(self) -> None:
        rows = [
            ["WP_1.1", "desc", "Species one", 111.0, 700.0, "longA",
             "10-90", "100-300", "301-690", "normal", "YK", "NA", "NA"],
            ["WP_2.1", "desc", "Species two", 222.0, 500.0, "short",
             "NULL", "30-260", "261-495", "NULL", "HK", "NA", "NA"],
            ["WP_3.1", "desc", "Species three", 333.0, 610.0, "unkn",
             "NULL", "250-390", "391-603", "NULL", "NA", "NA", "NA"],
            ["WP_4.1", "desc", "Species four", 444.0, 470.0, "longB_trun",
             "NULL", "40-170", "171-460", "NULL", "YK", "NA", "NA"],
        ]
        records = parse_table_s1_rows(rows)
        self.assertEqual(
            [r["source_phylogenetic_clade"] for r in records],
            ["LONG_A", "SHORT", "UNRESOLVED", "LONG_B"],
        )
        self.assertEqual([r["truncated_flag"] for r in records],
                         [False, False, False, True])
        self.assertTrue(all(r["ago_family"] == "PAGO" for r in records))
        self.assertEqual(AGO_TYPE_TO_CLADE["longA_trun"], "LONG_A")

    def test_parse_table_s1_rows_rejects_unknown_ago_type(self) -> None:
        rows = [["WP_x.1", "d", "s", 1.0, 100.0, "medium",
                 "NULL", "10-40", "41-99", "NULL", "NA", "NA", "NA"]]
        with self.assertRaisesRegex(RuntimeError, "Ago_type"):
            parse_table_s1_rows(rows)


if __name__ == "__main__":
    unittest.main()
