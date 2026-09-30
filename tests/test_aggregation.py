"""集計規則（臨床段階、重複、欠測、細胞分類）を検証する。"""

import math
import unittest
from collections.abc import Callable

from autoimmune_atlas.aggregation import (
    cell_catalog,
    expression_metadata,
    expression_state,
    filtered_records,
    location_classes,
    summarize,
    target_class_of,
    target_class_options,
    target_matches,
)
from autoimmune_atlas.models import (
    AggregationSnapshot,
    CellDefinition,
    ExpressionRow,
    ExpressionStateInput,
    TargetAnnotation,
)
from tests.core_fixture import core_snapshot


class AggregationTests(unittest.TestCase):
    def test_all_source_cells_count_each_target_and_drug_once(self) -> None:
        row = summarize(core_snapshot(), "all", 0.5, level="all", disease_ids=["D1"])[0]
        self.assertEqual(row["cell_id"], "all")
        self.assertEqual(set(row["member_cell_ids"]), {"T4", "T8", "B1", "X"})
        self.assertEqual((row["count"], row["denominator"], row["unknown"]), (3, 3, 0))
        self.assertEqual((row["drug_count"], row["drug_denominator"]), (2, 2))

    def test_mixed_cells_preserve_group_and_source_counts_with_distinct_ids(
        self,
    ) -> None:
        snapshot = core_snapshot()
        catalog = cell_catalog(snapshot, "mixed")
        ids = [cell["id"] for cell in catalog]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertIn("group:X", ids)
        self.assertIn("X", ids)
        selected = ["group:CL_0000084", "T4", "X", "group:X"]
        rows = summarize(
            snapshot, "all", 0.5, level="mixed", cell_ids=selected, disease_ids=["D1"]
        )
        self.assertEqual({row["cell_id"] for row in rows}, set(selected))
        for row in rows:
            original = next(
                item
                for item in summarize(snapshot, "all", 0.5, level=row["cell_level"])
                if item["disease_id"] == "D1" and item["cell_id"] == row["ontology_id"]
            )
            for key in (
                "count",
                "percent",
                "drug_count",
                "drug_percent",
                "denominator",
                "unknown",
                "member_cell_ids",
            ):
                self.assertEqual(row[key], original[key], key)
        self.assertEqual(
            summarize(snapshot, "all", 0.5, level="mixed", cell_ids=[]), []
        )

    def test_stage_is_maximized_within_disease_and_canonical_drug(self) -> None:
        rows = filtered_records(core_snapshot(), "antibody", "phase3")
        self.assertEqual(
            {(r["disease_id"], r["drug_id"], r["target_id"]) for r in rows},
            {
                ("D1", "A-SALT", "G1"),
                ("D1", "A", "G2"),
                ("D1", "B", "G3"),
                ("D1", "U", ""),
            },
        )
        self.assertEqual({r["canonical_stage"] for r in rows}, {"PHASE_3"})

    def test_canonical_stage_is_chosen_before_modality_filter(self) -> None:
        snapshot = core_snapshot()
        snapshot["records"] = snapshot["records"][:2]
        snapshot["records"][1]["drug_type"] = "Protein"
        rows = filtered_records(snapshot, "antibody", "phase3")
        self.assertEqual(
            [(row["drug_id"], row["canonical_stage"]) for row in rows],
            [("A-SALT", "PHASE_3")],
        )

    def test_canonical_drug_and_multiple_targets_are_counted_once(self) -> None:
        row = next(
            r
            for r in summarize(core_snapshot(), "all", 0.5, level="group")
            if r["disease_id"] == "D1" and r["cell_id"] == "CL_0000084"
        )
        self.assertEqual((row["count"], row["denominator"]), (2, 3))
        self.assertEqual((row["drug_count"], row["drug_denominator"]), (1, 2))
        self.assertEqual(
            (row["total_drugs"], row["mapped_drugs"], row["unmapped_drugs"]), (3, 2, 1)
        )
        self.assertEqual({r["drug_id"] for r in row["records"]}, {"A", "A-SALT"})

    def test_target_and_drug_unknown_states_are_independent(self) -> None:
        row = next(
            r
            for r in summarize(core_snapshot(), "all", 0.5, level="cell")
            if r["disease_id"] == "D1" and r["cell_id"] == "T4"
        )
        self.assertEqual((row["count"], row["unknown"]), (1, 1))
        self.assertEqual((row["drug_count"], row["unknown_drugs"]), (1, 1))
        percent = row["percent"]
        assert percent is not None
        self.assertAlmostEqual(percent, 100 / 3)
        self.assertEqual(row["drug_percent"], 50)

    def test_unmapped_drug_does_not_change_known_percentages(self) -> None:
        snapshot = core_snapshot()
        snapshot["records"] = [snapshot["records"][1], snapshot["records"][3]]
        row = next(
            r
            for r in summarize(snapshot, "all", 3.0, level="cell")
            if r["cell_id"] == "T4" and r["disease_id"] == "D1"
        )
        self.assertEqual(row["count"], 0)
        self.assertEqual(row["drug_count"], 0)
        self.assertEqual(row["percent"], 0)
        self.assertEqual(row["drug_percent"], 0)

    def test_unknown_expression_keeps_count_zero_and_percent_unknown(self) -> None:
        snapshot = core_snapshot()
        snapshot["records"] = [snapshot["records"][2]]
        row = next(
            r
            for r in summarize(snapshot, "all", 0.5, level="cell")
            if r["cell_id"] == "T4" and r["disease_id"] == "D1"
        )
        self.assertEqual((row["count"], row["drug_count"]), (0, 0))
        self.assertEqual((row["unknown"], row["unknown_drugs"]), (1, 1))
        self.assertIsNone(row["percent"])
        self.assertIsNone(row["drug_percent"])

    def test_relative_median_uses_all_reference_cells_and_missing_is_unknown(
        self,
    ) -> None:
        snapshot = core_snapshot()
        b_cell = next(
            r
            for r in summarize(snapshot, "all", 0.5, method="relative", level="cell")
            if r["disease_id"] == "D1" and r["cell_id"] == "B1"
        )
        self.assertEqual(b_cell["count"], 1)
        t4 = next(
            r
            for r in summarize(snapshot, "all", 0.5, method="relative", level="cell")
            if r["disease_id"] == "D1" and r["cell_id"] == "T4"
        )
        self.assertEqual(t4["unknown"], 1)
        self.assertEqual(
            expression_metadata(snapshot)["G1", "B1"]["target_median"], 0.55
        )

    def test_stage_filter_boundaries_follow_mixed_phase_policy(self) -> None:
        snapshot = core_snapshot()
        stages = (
            "PHASE_1",
            "PHASE_1_2",
            "PHASE_2",
            "PHASE_2_3",
            "PHASE_3",
            "PREAPPROVAL",
            "APPROVAL",
            "PHASE_4",
        )
        snapshot["records"] = [
            {
                **snapshot["records"][1],
                "drug_id": value,
                "drug": value,
                "canonical_drug_id": value,
                "canonical_drug": value,
                "stage": value,
            }
            for value in stages
        ]
        expected = {
            "phase1": set(stages),
            "phase2": set(stages[2:]),
            "phase3": set(stages[4:]),
            "approved": {"APPROVAL", "PHASE_4"},
        }
        for stage, wanted in expected.items():
            with self.subTest(stage=stage):
                self.assertEqual(
                    {row["drug_id"] for row in filtered_records(snapshot, stage=stage)},
                    wanted,
                )

    def test_cpm_and_specificity_boundaries_are_inclusive(self) -> None:
        row = next(
            r
            for r in summarize(
                core_snapshot(),
                "all",
                0.5,
                method="specificity",
                specificity_threshold=0.75,
                level="cell",
            )
            if r["disease_id"] == "D1" and r["cell_id"] == "T4"
        )
        self.assertEqual(row["count"], 1)
        self.assertEqual(row["records"][0]["specificity_score"], 0.75)
        boundary = core_snapshot()
        boundary["expression"]["G1"][0]["median"] = 0.5
        row = next(
            r
            for r in summarize(boundary, "all", 0.5, method="specificity", level="cell")
            if r["disease_id"] == "D1" and r["cell_id"] == "T4"
        )
        self.assertIn("G1", {record["target_id"] for record in row["records"]})
        at_threshold: ExpressionStateInput = {
            "median": 0.5,
            "target_median": 0.5,
            "specificity_score": 0.75,
        }
        for method in ("fixed", "relative", "specificity"):
            self.assertTrue(expression_state(at_threshold, 0.5, method, 0.75))
        self.assertFalse(expression_state({"median": 0.49}, 0.5, "fixed"))
        metadata = expression_metadata(core_snapshot())["G1", "T4"]
        self.assertTrue(expression_state(metadata, 0.5, "specificity", 0.75))
        self.assertTrue(
            expression_state(
                {"median": 1, "specificity_score": 0.5}, 0.5, "specificity"
            )
        )
        with self.assertRaises(ValueError):
            expression_state(metadata, 0.5, "guess")

    def test_group_catalog_unions_t_cells_and_keeps_unknown_parent(self) -> None:
        groups = cell_catalog(core_snapshot(), "group")
        self.assertIn(
            {"id": "CL_0000084", "name": "T cell", "members": ["T4", "T8"]}, groups
        )
        self.assertIn({"id": "X", "name": "Novel cell", "members": ["X"]}, groups)
        self.assertIn(
            {"id": "T4", "name": "CD4 T cell", "members": ["T4"]},
            cell_catalog(core_snapshot(), "cell"),
        )

    def test_catalog_uses_lineage_group_priority_and_keeps_every_cell(self) -> None:
        definitions: list[tuple[str, str, str, str, list[str]]] = [
            ("CD8_EFFECTOR", "A CD8 effector", "OTHER", "Other", ["CL_0000084", "CD8"]),
            ("UNKNOWN_Z", "Zeta cell", "UNKNOWN_Z", "Zeta group", []),
            ("CD4_MEMORY", "A CD4 memory", "OTHER", "Other", ["CL_0000084", "CD4"]),
            ("B", "B cell", "CL_0000945", "B lineage", ["CL_0000945"]),
            ("CD8", "CD8 T cell", "OTHER", "Other", ["CL_0000084"]),
            ("UNKNOWN_A", "Alpha cell", "UNKNOWN_A", "Alpha group", []),
            ("CD4", "CD4 T cell", "OTHER", "Other", ["CL_0000084"]),
        ]
        cell_table: dict[str, CellDefinition] = {
            cell_id: {
                "name": name,
                "parent_id": parent_id,
                "parent": parent,
                "ancestor_ids": ancestors,
            }
            for cell_id, name, parent_id, parent, ancestors in definitions
        }
        rows: list[ExpressionRow] = [
            {"cell_id": cell_id, "median": 1, "specificity_score": None}
            for cell_id, *_ in definitions
        ]
        snapshot: AggregationSnapshot = {
            "schema": 4,
            "diseases": [],
            "records": [],
            "cells": cell_table,
            "expression": {"G": rows},
        }

        groups = cell_catalog(snapshot, "group")
        cells = cell_catalog(snapshot, "cell")
        self.assertEqual(
            [group["id"] for group in groups],
            ["CL_0000084", "CL_0000945", "UNKNOWN_A", "UNKNOWN_Z"],
        )
        self.assertEqual(
            groups[0]["members"], ["CD4", "CD4_MEMORY", "CD8", "CD8_EFFECTOR"]
        )
        self.assertEqual(
            [cell["id"] for cell in cells],
            ["CD4", "CD4_MEMORY", "CD8", "CD8_EFFECTOR", "B", "UNKNOWN_A", "UNKNOWN_Z"],
        )
        self.assertEqual(
            {cell["id"] for cell in cells}, {row["cell_id"] for row in rows}
        )

        snapshot["cells"] = dict(reversed(list(cell_table.items())))
        self.assertEqual(cell_catalog(snapshot, "group"), groups)
        self.assertEqual(cell_catalog(snapshot, "cell"), cells)

    def test_catalog_uses_curated_cell_state_order_with_lineage_fallback(self) -> None:
        definitions: list[tuple[str, str, str, str, list[str]]] = [
            (
                "CL_0000896",
                "activated CD4-positive, alpha-beta T cell",
                "CL_0000624",
                "CD4-positive, alpha-beta T cell",
                ["CL_0000084", "CL_0000624"],
            ),
            (
                "CL_0000232",
                "erythrocyte",
                "CL_0000764",
                "erythroid lineage cell",
                ["CL_0000764"],
            ),
            (
                "CL_0000018",
                "spermatid",
                "CL_0000015",
                "male germ cell",
                ["CL_0000015"],
            ),
            (
                "CL_0000624",
                "CD4-positive, alpha-beta T cell",
                "CL_0000084",
                "T cell",
                ["CL_0000084"],
            ),
            (
                "CL_9999999",
                "future CD4 subtype",
                "CL_0000624",
                "CD4-positive, alpha-beta T cell",
                ["CL_0000084", "CL_0000624"],
            ),
            (
                "CL_0000020",
                "spermatogonium",
                "CL_0000015",
                "male germ cell",
                ["CL_0000015"],
            ),
            ("CL_0000084", "T cell", "CL_0000084", "T cell", []),
            (
                "CL_0000038",
                "erythroid progenitor cell",
                "CL_0000764",
                "erythroid lineage cell",
                ["CL_0000764"],
            ),
            (
                "CL_0000895",
                "naive thymus-derived CD4-positive, alpha-beta T cell",
                "CL_0000624",
                "CD4-positive, alpha-beta T cell",
                ["CL_0000084", "CL_0000624"],
            ),
            (
                "CL_0000015",
                "male germ cell",
                "CL_0000015",
                "male germ cell",
                [],
            ),
            (
                "CL_0000017",
                "spermatocyte",
                "CL_0000015",
                "male germ cell",
                ["CL_0000015"],
            ),
        ]
        cells: dict[str, CellDefinition] = {
            cell_id: {
                "name": name,
                "parent_id": parent_id,
                "parent": parent,
                "ancestor_ids": ancestors,
            }
            for cell_id, name, parent_id, parent, ancestors in definitions
        }
        rows: list[ExpressionRow] = [
            {"cell_id": cell_id, "median": 1, "specificity_score": None}
            for cell_id, *_ in definitions
        ]
        snapshot: AggregationSnapshot = {
            "schema": 4,
            "diseases": [],
            "records": [],
            "cells": cells,
            "expression": {"G": rows},
        }

        groups = {group["id"]: group["members"] for group in cell_catalog(snapshot)}

        self.assertEqual(
            groups["CL_0000084"],
            [
                "CL_0000084",
                "CL_0000624",
                "CL_0000895",
                "CL_0000896",
                "CL_9999999",
            ],
        )
        self.assertEqual(groups["CL_0000764"], ["CL_0000038", "CL_0000232"])
        self.assertEqual(
            groups["CL_0000015"],
            ["CL_0000015", "CL_0000020", "CL_0000017", "CL_0000018"],
        )
        self.assertCountEqual(
            [cell_id for members in groups.values() for cell_id in members],
            [cell_id for cell_id, *_ in definitions],
        )

    def test_expression_cell_missing_from_cells_table_fails(self) -> None:
        snapshot = core_snapshot()
        snapshot["expression"]["G2"].append(
            {"cell_id": "GHOST", "median": 1.0, "specificity_score": None}
        )
        with self.assertRaisesRegex(ValueError, "cells に無い"):
            expression_metadata(snapshot)

    def test_metadata_rows_carry_cell_name_from_cells_table(self) -> None:
        metadata = expression_metadata(core_snapshot())
        self.assertEqual(metadata["G1", "T4"]["cell"], "CD4 T cell")
        self.assertEqual(metadata["G3", "X"]["cell"], "Novel cell")

    def test_invalid_modes_and_nonfinite_thresholds_fail(self) -> None:
        invalid_calls: tuple[tuple[str, Callable[[], object]], ...] = (
            ("stage", lambda: summarize(core_snapshot(), "all", 0.5, stage="late")),
            ("method", lambda: summarize(core_snapshot(), "all", 0.5, method="guess")),
            ("level", lambda: summarize(core_snapshot(), "all", 0.5, level="organ")),
            (
                "specificity_threshold",
                lambda: summarize(
                    core_snapshot(), "all", 0.5, specificity_threshold=math.inf
                ),
            ),
        )
        for name, call in invalid_calls:
            with self.subTest(parameter=name), self.assertRaises(ValueError):
                call()
        with self.assertRaises(ValueError):
            summarize(core_snapshot(), "all", math.nan)
        with self.assertRaises(ValueError):
            summarize(core_snapshot(), "all", -0.1)
        with self.assertRaises(ValueError):
            summarize(core_snapshot(), "all", 0.5, specificity_threshold=1.1)
        with self.assertRaises(ValueError):
            filtered_records(core_snapshot(), "other", "phase3")


def _annotation(*locations: str, target_class: str | None = None) -> TargetAnnotation:
    return {
        "target_class": target_class,
        "locations": [
            {"location": location, "source": "uniprot"} for location in locations
        ],
    }


class TargetFilterTests(unittest.TestCase):
    def test_location_classes_follow_leading_term_in_class_order(self) -> None:
        cases: tuple[tuple[TargetAnnotation | None, list[str]], ...] = (
            (None, ["unknown"]),
            (_annotation(), ["unknown"]),
            (_annotation("Cytoplasm"), ["intracellular"]),
            (_annotation("Nucleus, nucleoplasm"), ["intracellular"]),
            (
                _annotation("Nucleus", "Secreted", "Cell membrane ; Single-pass"),
                ["secreted", "cell_surface", "intracellular"],
            ),
            (_annotation("Extracellular space"), ["secreted"]),
            (_annotation("Predicted to be secreted"), ["secreted"]),
            (_annotation("Apical cell membrane"), ["cell_surface"]),
            (_annotation("Plasma membrane"), ["cell_surface"]),
            (_annotation("Cell junctions"), ["cell_surface"]),
            (_annotation("Golgi apparatus membrane"), ["intracellular"]),
            (_annotation("Late endosome membrane"), ["intracellular"]),
            (_annotation("Rough endoplasmic reticulum"), ["intracellular"]),
            (_annotation("Membrane ; Multi-pass membrane protein"), ["unknown"]),
            (_annotation("Note=Found in the cytoplasm"), ["unknown"]),
            (_annotation("Cell projection, cilium"), ["unknown"]),
        )
        for annotation, expected in cases:
            with self.subTest(annotation=annotation):
                self.assertEqual(location_classes(annotation), expected)

    def test_target_class_of_uses_unknown_when_missing(self) -> None:
        self.assertEqual(target_class_of(_annotation(target_class="Enzyme")), "Enzyme")
        self.assertEqual(target_class_of(_annotation()), "unknown")
        self.assertEqual(target_class_of(None), "unknown")

    def test_target_matches_requires_both_class_and_location(self) -> None:
        enzyme = _annotation("Cytoplasm", target_class="Enzyme")
        self.assertTrue(target_matches(enzyme, "all", "all"))
        self.assertTrue(target_matches(enzyme, "Enzyme", "intracellular"))
        self.assertFalse(target_matches(enzyme, "Enzyme", "secreted"))
        self.assertFalse(target_matches(enzyme, "Kinase", "intracellular"))
        self.assertTrue(target_matches(None, "unknown", "unknown"))

    def test_class_options_are_sorted_and_exclude_unknown(self) -> None:
        self.assertEqual(
            target_class_options(core_snapshot()), ["Enzyme", "Membrane receptor"]
        )

    def test_filtered_records_keep_matching_targets_and_carry_classes(self) -> None:
        rows = filtered_records(core_snapshot(), "all", "phase3")
        by_target = {r["target_id"]: r for r in rows}
        self.assertEqual(by_target["G1"].get("target_class"), "Enzyme")
        self.assertEqual(by_target["G1"].get("location_classes"), ["intracellular"])
        self.assertEqual(
            by_target["G2"].get("location_classes"), ["secreted", "cell_surface"]
        )
        self.assertEqual(by_target["G3"].get("target_class"), "unknown")
        self.assertEqual(by_target[""].get("location_classes"), ["unknown"])
        enzyme = filtered_records(core_snapshot(), "all", "phase3", "Enzyme")
        # 標的で絞っても、canonical stage は絞る前の全元記録から求める。
        self.assertEqual(
            [(r["drug_id"], r["canonical_stage"]) for r in enzyme],
            [("A-SALT", "PHASE_3")],
        )
        surface = filtered_records(
            core_snapshot(), "all", "phase3", location="cell_surface"
        )
        self.assertEqual([r["target_id"] for r in surface], ["G2"])
        unknown = filtered_records(core_snapshot(), "all", "phase3", "unknown")
        self.assertEqual({r["target_id"] for r in unknown}, {"G3", ""})
        row = summarize(
            core_snapshot(),
            "all",
            0.5,
            level="all",
            disease_ids=["D1"],
            target_class="Membrane receptor",
        )[0]
        self.assertEqual(row["denominator"], 1)

    def test_invalid_target_filters_fail(self) -> None:
        with self.assertRaises(ValueError):
            filtered_records(core_snapshot(), "all", "phase3", "Kinase")
        with self.assertRaises(ValueError):
            filtered_records(core_snapshot(), "all", "phase3", location="organ")
        with self.assertRaises(ValueError):
            summarize(core_snapshot(), "all", 0.5, location="organ")
