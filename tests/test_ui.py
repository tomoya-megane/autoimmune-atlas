"""schema 2 の Dash UI と表示用ロジックを検証する。"""

from __future__ import annotations

import json
import math
import os
import runpy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from dash import html

from backend import aggregation as atlas
from backend import refresh, snapshot
from backend.ui import (
    application,
    callbacks,
    components,
    config,
    figures,
    layout,
)

ASSETS_PATH = Path(__file__).resolve().parents[1] / "assets"


def fixture() -> dict:
    """段階、成分重複、欠測、group union を含む最小 fixture を返す。"""
    common = {"disease_id": "MONDO_RA_TEST", "disease": "rheumatoid arthritis"}
    return {
        "schema": 2,
        "root": {"id": "MONDO_AUTOIMMUNE", "name": "autoimmune disease"},
        "data_version": {"year": 26, "month": 6, "iteration": 1},
        "retrieved_at": "2026-09-22T12:00:00Z",
        "source": "https://api.platform.opentargets.org/",
        "diseases": [
            {"id": "MONDO_RA_TEST", "name": "rheumatoid arthritis", "status": "ready"}
        ],
        "records": [
            {
                **common,
                "drug_id": "DRUG_FORM_1",
                "drug": "=alpha salt",
                "canonical_drug_id": "DRUG_ALPHA",
                "canonical_drug": "alpha",
                "drug_type": "Protein",
                "stage": "PHASE_3",
                "target_id": "ENSG_TARGET_1",
                "target": "TARGET1",
                "mechanism": "inhibitor",
                "action_types": ["INHIBITOR"],
                "references": [
                    {
                        "source": "PubMed",
                        "ids": ["123"],
                        "urls": ["https://pubmed.ncbi.nlm.nih.gov/123/"],
                    }
                ],
            },
            {
                **common,
                "drug_id": "DRUG_FORM_2",
                "drug": "alpha hydrate",
                "canonical_drug_id": "DRUG_ALPHA",
                "canonical_drug": "alpha",
                "drug_type": "Protein",
                "stage": "PHASE_2",
                "target_id": "ENSG_TARGET_2",
                "target": "TARGET2",
                "mechanism": "modulator",
                "action_types": ["MODULATOR"],
                "references": [],
            },
            {
                **common,
                "drug_id": "DRUG_BETA",
                "drug": "beta",
                "canonical_drug_id": "DRUG_BETA",
                "canonical_drug": "beta",
                "drug_type": "Antibody",
                "stage": "APPROVAL",
                "target_id": None,
                "target": None,
                "mechanism": None,
                "action_types": [],
                "references": [
                    {"source": "ClinicalTrials", "ids": ["NCT1"], "urls": []}
                ],
            },
            {
                **common,
                "drug_id": "DRUG_GAMMA",
                "drug": "gamma",
                "canonical_drug_id": "DRUG_GAMMA",
                "canonical_drug": "gamma",
                "drug_type": "Small molecule",
                "stage": "PHASE_3",
                "target_id": "ENSG_TARGET_2",
                "target": "TARGET2",
                "mechanism": "agonist",
                "action_types": ["AGONIST"],
                "references": [],
            },
        ],
        "expression": {
            "ENSG_TARGET_1": [
                {
                    "cell_id": "CL_B_ONE",
                    "cell": "memory B cell",
                    "median": 2.0,
                    "specificity_score": 0.8,
                    "parent_id": "CL_B_GROUP",
                    "parent": "B cell",
                    "ancestor_ids": ["CL_B_GROUP"],
                },
                {
                    "cell_id": "CL_B_TWO",
                    "cell": "naive B cell",
                    "median": 0.1,
                    "specificity_score": 0.2,
                    "parent_id": "CL_B_GROUP",
                    "parent": "B cell",
                    "ancestor_ids": ["CL_B_GROUP"],
                },
                {
                    "cell_id": "CL_T_ONE",
                    "cell": "CD8-positive T cell",
                    "median": 0.2,
                    "specificity_score": 0.1,
                    "parent_id": "CL_T_PARENT",
                    "parent": "T lymphocyte",
                    "ancestor_ids": [atlas.T_CELL_ID],
                },
            ],
            "ENSG_TARGET_2": [
                {
                    "cell_id": "CL_B_ONE",
                    "cell": "memory B cell",
                    "median": None,
                    "specificity_score": None,
                    "parent_id": "CL_B_GROUP",
                    "parent": "B cell",
                    "ancestor_ids": ["CL_B_GROUP"],
                },
                {
                    "cell_id": "CL_B_TWO",
                    "cell": "naive B cell",
                    "median": 0.6,
                    "specificity_score": 0.75,
                    "parent_id": "CL_B_GROUP",
                    "parent": "B cell",
                    "ancestor_ids": ["CL_B_GROUP"],
                },
                {
                    "cell_id": "CL_T_ONE",
                    "cell": "CD8-positive T cell",
                    "median": 1.0,
                    "specificity_score": 0.7,
                    "parent_id": "CL_T_PARENT",
                    "parent": "T lymphocyte",
                    "ancestor_ids": [atlas.T_CELL_ID],
                },
            ],
        },
    }


def summary(**overrides) -> dict:
    row = {
        "disease_id": "D1",
        "disease": "Disease one",
        "cell_id": "C1",
        "cell": "B cell",
        "count": 1,
        "percent": 20.0,
        "denominator": 5,
        "unknown": 0,
        "drug_count": 1,
        "drug_percent": 25.0,
        "drug_denominator": 4,
        "unknown_drugs": 0,
        "mapped_drugs": 4,
        "total_drugs": 4,
        "unmapped_drugs": 0,
        "status": "complete",
        "member_cell_ids": ["C1"],
        "records": [],
    }
    row.update(overrides)
    return row


class DiseaseTreeTests(unittest.TestCase):
    def test_sections_follow_parents_and_place_each_term_once(self) -> None:
        family = {
            "id": "ROOT",
            "label": "root disease",
            "diseases": [
                {"id": "ROOT", "name": "root disease", "parent_ids": []},
                {"id": "CHILD", "name": "child", "parent_ids": ["ROOT"]},
                {"id": "LEAF", "name": "leaf", "parent_ids": ["ROOT"]},
                {"id": "GRAND", "name": "grandchild", "parent_ids": ["CHILD"]},
                {
                    "id": "SHARED",
                    "name": "shared term",
                    "parent_ids": ["ROOT", "GRAND"],
                },
                {"id": "LOOP_A", "name": "loop a", "parent_ids": ["LOOP_B"]},
                {"id": "LOOP_B", "name": "loop b", "parent_ids": ["LOOP_A"]},
            ],
        }
        sections = components._disease_checklist_sections(family)
        ids = {s["id"]: [d["id"] for d in s["diseases"]] for s in sections}
        self.assertEqual(
            ids,
            {
                "disease-family-ROOT": ["ROOT"],
                "disease-details-ROOT": ["LEAF"],
                "disease-family-ROOT-CHILD": ["CHILD"],
                "disease-family-ROOT-GRAND": ["GRAND"],
                "disease-details-ROOT-GRAND": ["SHARED"],
                "disease-family-ROOT-LOOP_A": ["LOOP_A"],
                "disease-details-ROOT-LOOP_A": ["LOOP_B"],
            },
        )
        placed = [d for s in sections for d in s["diseases"]]
        self.assertEqual(len(placed), len(family["diseases"]))
        layout = components.disease_selector(
            {"diseases": [{**d, "status": "ready"} for d in family["diseases"]]},
            [],
        )
        summaries = []

        def walk(component):
            if isinstance(component, html.Summary):
                summaries.append(component.children)
            for child in getattr(component, "children", None) or []:
                if not isinstance(child, str):
                    walk(child)

        walk(layout)
        self.assertIn("child (3 terms)", summaries)
        self.assertIn("child details (2 terms)", summaries)
        self.assertIn("grandchild (2 terms)", summaries)

    def test_family_without_root_keeps_flat_top_level(self) -> None:
        family = {
            "id": "other",
            "label": "Other terms",
            "diseases": [
                {"id": "B", "name": "b term", "parent_ids": []},
                {"id": "A", "name": "a term", "parent_ids": ["MISSING"]},
            ],
        }
        sections = components._disease_checklist_sections(family)
        self.assertEqual(
            [(s["id"], [d["id"] for d in s["diseases"]]) for s in sections],
            [("disease-family-other", ["A", "B"])],
        )


class FigureTests(unittest.TestCase):
    """2種類の指標で 0、下限、欠測を独立に扱う。"""

    def test_target_and_drug_lower_bounds_are_independent(self) -> None:
        row = summary(
            count=0, percent=0.0, drug_count=1, drug_percent=25.0, unknown_drugs=1
        )
        self.assertEqual(figures._display_value(row, "count"), "0")
        self.assertEqual(figures._display_value(row, "count", "drug"), "≥1")
        self.assertEqual(figures._display_value(row, "percent"), "0.0%")
        self.assertEqual(figures._display_value(row, "percent", "drug"), "≥25.0%")

    def test_unmapped_drug_does_not_change_percent_lower_bound(self) -> None:
        row = summary(
            count=2,
            percent=40.0,
            drug_count=1,
            drug_percent=25.0,
            unmapped_drugs=1,
            total_drugs=5,
            status="partial",
        )
        self.assertEqual(figures._display_value(row, "count"), "≥2")
        self.assertEqual(figures._display_value(row, "percent"), "40.0%")
        self.assertEqual(figures._display_value(row, "count", "drug"), "≥1")
        self.assertEqual(figures._display_value(row, "percent", "drug"), "25.0%")

    def test_long_disease_ticks_wrap_without_changing_data(self) -> None:
        name = "anti-neutrophil cytoplasmic antibody-associated vasculitis"
        rows = [summary(disease=name), summary(disease_id="D2", disease="Short name")]
        for kind in ("target", "drug"):
            figure = figures.build_figure(rows, ["D1", "D2"], ["C1"], "count", kind)
            self.assertEqual(
                list(figure.layout.xaxis.ticktext),
                [
                    "anti-neutrophil cytoplasmic<br>antibody-associated vasculitis",
                    "Short name",
                ],
            )
            self.assertEqual(list(figure.data[0].x), [name, "Short name"])
            self.assertEqual(figure.data[0].customdata[0][0], ["D1", "C1"])
            self.assertEqual(figure.layout.xaxis.tickangle, -45)

    def test_heatmap_cells_omit_lower_bound_mark_but_keep_hover_and_values(
        self,
    ) -> None:
        row = summary(
            count=2,
            percent=40.0,
            drug_count=1,
            drug_percent=25.0,
            unknown=1,
            unknown_drugs=1,
        )
        for kind, measure, label, value in (
            ("target", "count", "2", 2),
            ("target", "percent", "40.0%", 40.0),
            ("drug", "count", "1", 1),
            ("drug", "percent", "25.0%", 25.0),
        ):
            figure = figures.build_figure([row], ["D1"], ["C1"], measure, kind)
            self.assertEqual(figure.data[0].text[0][0], label)
            self.assertEqual(figure.data[0].z[0][0], value)
            self.assertIn("≥" + label, figure.data[0].hovertext[0][0])
            self.assertIn("≥ is a lower bound", figure.data[0].hovertext[0][0])

    def test_color_ranges_follow_each_visible_measure(self) -> None:
        rows = [
            summary(),
            summary(cell_id="C2", cell="T cell", percent=60, drug_percent=40),
        ]
        target = figures.build_figure(rows, ["D1"], ["C1", "C2"], "percent", "target")
        drug = figures.build_figure(rows, ["D1"], ["C1", "C2"], "percent", "drug")
        self.assertEqual((target.data[0].zmin, target.data[0].zmax), (20, 60))
        self.assertEqual((drug.data[0].zmin, drug.data[0].zmax), (25, 40))
        flat = figures.build_figure([summary(percent=0)], ["D1"], ["C1"], "percent")
        self.assertEqual((flat.data[0].zmin, flat.data[0].zmax), (0, 1))
        missing = figures.build_figure(
            [summary(percent=None)], ["D1"], ["C1"], "percent"
        )
        self.assertEqual((missing.data[0].zmin, missing.data[0].zmax), (0, 1))

    def test_missing_percent_uses_same_fill_and_text_as_zero(self):
        rows = [
            summary(percent=0, drug_percent=0),
            summary(cell_id="C2", cell="T cell", percent=None, drug_percent=None),
        ]
        for kind in ("target", "drug"):
            figure = figures.build_figure(rows, ["D1"], ["C1", "C2"], "percent", kind)
            self.assertEqual(figure.data[0].z[0][0], figure.data[0].z[1][0])
            self.assertEqual(figure.data[0].text[0][0], figure.data[0].text[1][0])
        self.assertIsNone(rows[1]["percent"])
        self.assertIsNone(rows[1]["drug_percent"])

    def test_two_figures_have_aligned_axes_and_zero_for_na(self) -> None:
        rows = [
            summary(),
            summary(
                cell_id="C2",
                cell="T cell",
                count=None,
                percent=None,
                drug_count=None,
                drug_percent=None,
                unknown=1,
                unknown_drugs=1,
                status="unavailable",
            ),
        ]
        target = figures.build_figure(rows, ["D1"], ["C1", "C2"], "count", "target")
        drug = figures.build_figure(rows, ["D1"], ["C1", "C2"], "count", "drug")
        for figure in (target, drug):
            self.assertIsNone(figure.layout.width)
            self.assertEqual(figure.data[0].colorbar.lenmode, "pixels")
            self.assertEqual(figure.data[0].colorbar.len, 240)
            self.assertEqual(figure.data[0].colorbar.y, 1)
            self.assertEqual(figure.data[0].colorbar.yanchor, "top")
        for figure in (target, drug):
            self.assertFalse(figure.layout.yaxis.autorange)
            self.assertEqual(list(figure.layout.yaxis.range), [1.5, -0.5])
        self.assertEqual(list(target.data[0].x), list(drug.data[0].x))
        self.assertEqual(list(target.data[0].y), list(drug.data[0].y))
        self.assertEqual(len(target.data), 1)
        self.assertEqual(target.data[0].z[1][0], 0)
        self.assertEqual(target.data[0].text[1][0], "0")
        self.assertIsNone(rows[1]["count"])
        hover = target.data[0].hovertext[1][0]
        self.assertIn("Denominator: 5", hover)
        self.assertIn("Unresolved in denominator: 1", hover)
        self.assertIn("Drugs with known targets / all drugs: 4 / 4", hover)
        self.assertIn("not a measured zero", hover)

    def test_count_zero_with_unknown_evidence_has_no_cross(self) -> None:
        row = summary(
            count=0,
            drug_count=0,
            percent=None,
            drug_percent=None,
            unknown=1,
            unknown_drugs=1,
            status="partial",
        )
        for kind in ("target", "drug"):
            figure = figures.build_figure([row], ["D1"], ["C1"], "count", kind)
            self.assertEqual(figure.data[0].z[0][0], 0)
            self.assertEqual(figure.data[0].text[0][0], "0")
            self.assertEqual(len(figure.data), 1)
            self.assertIn("≥0", figure.data[0].hovertext[0][0])
        percent = figures.build_figure([row], ["D1"], ["C1"], "percent")
        self.assertEqual(percent.data[0].z[0][0], 0)
        self.assertEqual(percent.data[0].text[0][0], "0.0%")
        self.assertIsNone(row["percent"])

    def test_no_value_hover_distinguishes_unavailable_from_empty_denominator(
        self,
    ) -> None:
        unavailable = summary(
            status="unavailable", count=None, percent=None, denominator=0
        )
        empty = summary(percent=None, denominator=0)
        self.assertIn(
            "Disease data unavailable",
            figures._hover_text(unavailable, "percent", "target"),
        )
        self.assertNotIn(
            "No eligible items", figures._hover_text(unavailable, "percent", "target")
        )
        self.assertIn(
            "No eligible items in the percentage denominator",
            figures._hover_text(empty, "percent", "target"),
        )


class EvidenceTests(unittest.TestCase):
    """詳細が陽性以外の元記録も保持する。"""

    def test_source_links_only_show_open_targets_pages(self) -> None:
        links = str(components._reference_links(fixture()["records"][0]))
        self.assertIn("Open Targets target", links)
        self.assertIn("Open Targets drug", links)
        self.assertIn("Open Targets canonical drug", links)
        self.assertIn("Open Targets disease", links)
        self.assertNotIn("PubMed", links)

    def setUp(self) -> None:
        self.snapshot = fixture()
        self.rows = callbacks.visible_rows(
            self.snapshot, "all", 0.5, ["MONDO_RA_TEST"], ["CL_B_GROUP"]
        )

    def test_source_records_use_only_selected_disease_without_pair_table(self) -> None:
        other = {"id": "MONDO_OTHER", "name": "other disease", "status": "ready"}
        self.snapshot["diseases"].append(other)
        self.snapshot["records"].append(
            {
                **self.snapshot["records"][0],
                "disease_id": other["id"],
                "disease": other["name"],
            }
        )
        all_rows = atlas.summarize(self.snapshot, "all", 0.5, level="all")
        panel = layout.detail_panel(all_rows, "MONDO_RA_TEST", self.snapshot)
        self.assertNotIn("Drug–target pairs", str(panel))
        self.assertNotIn("Targets meeting rule in any source cell", str(panel))
        self.assertNotIn("rheumatoid arthritis · all source cell types", str(panel))
        self.assertIn("Drug–target records for rheumatoid arthritis", str(panel))
        self.assertIn("Drug–target records", str(panel.children[0]))
        self.assertIn("Relative target expression", str(panel.children[2]))
        source_children = panel.children[1].children
        self.assertEqual(source_children[-1].className, "source-pagination")
        self.assertEqual(source_children[1].children.id, "source-records-page")
        context = source_children[0].data
        records = [
            record
            for record in atlas.filtered_records(
                self.snapshot, context["modality"], context["stage"]
            )
            if record["disease_id"] == context["disease_id"]
        ]
        self.assertEqual(
            {record["disease_id"] for record in records}, {"MONDO_RA_TEST"}
        )
        self.assertEqual(len(records), 4)

    def test_continuous_expression_uses_target_wise_z_score_and_raw_cpm_hover(
        self,
    ) -> None:
        evidence = figures.evidence_rows(
            self.snapshot,
            self.rows[0],
            modality="all",
            stage="phase3",
            threshold=0.5,
            method="fixed",
            specificity=0.75,
        )
        figure = figures.expression_figure(self.snapshot, evidence)
        targets = list(figure.data[0].x)
        cells = list(figure.data[0].y)
        target_index = next(i for i, label in enumerate(targets) if "TARGET1" in label)
        cell_index = cells.index("memory B cell")
        self.assertGreater(figure.data[0].z[cell_index][target_index], 0)
        target2_index = next(i for i, label in enumerate(targets) if "TARGET2" in label)
        self.assertAlmostEqual(
            figure.data[0].z[cells.index("CD8-positive T cell")][target2_index], 1
        )
        self.assertAlmostEqual(
            figure.data[0].z[cells.index("naive B cell")][target2_index], -1
        )
        self.assertIsNone(figure.data[0].z[cell_index][target2_index])
        self.assertIn(
            "Target-wise z-score:", figure.data[0].hovertext[cell_index][target_index]
        )
        self.assertIn(
            "Median CPM: 2", figure.data[0].hovertext[cell_index][target_index]
        )
        self.assertAlmostEqual(
            figure.data[0].zmax,
            max(
                abs(value)
                for row in figure.data[0].z
                for value in row
                if value is not None
            ),
        )
        self.assertEqual(figure.data[0].zmin, -figure.data[0].zmax)
        self.assertEqual(figure.data[0].colorbar.lenmode, "pixels")
        self.assertEqual(figure.data[0].colorbar.len, 240)
        self.assertEqual(figure.data[0].colorbar.y, 1)
        self.assertEqual(figure.data[0].colorbar.yanchor, "top")
        missing = next(
            trace for trace in figure.data if trace.name == "Missing expression"
        )
        self.assertEqual(missing.mode, "text")
        self.assertEqual(missing.text, "0")
        self.assertEqual(list(missing.x), ["TARGET2 (ENSG_TARGET_2)"])
        self.assertEqual(list(missing.y), ["memory B cell"])
        self.assertEqual(figure.layout.xaxis.side, "top")
        self.assertEqual(list(figure.layout.yaxis.range), [len(cells) - 0.5, -0.5])

    def test_group_expression_averages_cpm_before_transform_and_keeps_scale(self):
        figure = figures.expression_figure(
            self.snapshot, self.snapshot["records"], grouped=True
        )
        heatmap = figure.data[0]
        row = list(heatmap.y).index("B cell (group)")
        col = list(heatmap.x).index("TARGET1 (ENSG_TARGET_1)")
        values = [math.log2(1 + x) for x in (2, 0.1, 0.2)]
        center = sum(values) / len(values)
        spread = (sum((x - center) ** 2 for x in values) / len(values)) ** 0.5
        self.assertAlmostEqual(
            heatmap.z[row][col], (math.log2(1 + 1.05) - center) / spread
        )
        self.assertIn("Mean CPM: 1.05", heatmap.hovertext[row][col])
        self.assertIn("Observed source cell types: 2 / 2", heatmap.hovertext[row][col])
        col2 = list(heatmap.x).index("TARGET2 (ENSG_TARGET_2)")
        self.assertIn("Mean CPM: 0.6", heatmap.hovertext[row][col2])
        self.assertIn("Observed source cell types: 1 / 2", heatmap.hovertext[row][col2])
        original = figures.expression_figure(self.snapshot, self.snapshot["records"])
        self.assertEqual(original.data[0].x, heatmap.x)
        for i, name in enumerate(original.data[0].y):
            self.assertEqual(
                original.data[0].z[i], heatmap.z[list(heatmap.y).index(name)]
            )
        metadata = atlas.expression_metadata(self.snapshot)
        figures.expression_figure(
            self.snapshot, self.snapshot["records"], metadata, grouped=True
        )
        self.assertEqual(metadata, atlas.expression_metadata(self.snapshot))
        catalog = atlas.cell_catalog(self.snapshot, "mixed")
        collapsed = figures.expression_view(figure, catalog, [])
        expanded = figures.expression_view(figure, catalog, ["group:CL_B_GROUP"])
        self.assertEqual(len(collapsed.data[0].y), 2)
        self.assertEqual(len(expanded.data[0].y), 4)
        self.assertEqual(collapsed.data[0].zmax, expanded.data[0].zmax)
        self.assertEqual(collapsed.data[0].x, expanded.data[0].x)
        self.assertEqual(len(figure.data[0].y), 5)
        for trace in expanded.data[1:]:
            self.assertTrue(set(trace.y) <= set(expanded.data[0].y))
            if trace.name == "Meets expression rule":
                self.assertNotIn("B cell (group)", trace.y)

    def test_group_expression_is_missing_when_all_members_are_missing(self):
        for item in self.snapshot["expression"]["ENSG_TARGET_2"]:
            item["median"] = None
        figure = figures.expression_figure(
            self.snapshot, self.snapshot["records"], grouped=True
        )
        column = list(figure.data[0].x).index("TARGET2 (ENSG_TARGET_2)")
        self.assertTrue(all(row[column] is None for row in figure.data[0].z))
        row = list(figure.data[0].y).index("B cell (group)")
        self.assertIn(
            "Observed source cell types: 0 / 2", figure.data[0].hovertext[row][column]
        )

    def test_constant_target_expression_has_zero_z_score(self) -> None:
        for item in self.snapshot["expression"]["ENSG_TARGET_1"]:
            item["median"] = 2
        evidence = figures.evidence_rows(
            self.snapshot,
            self.rows[0],
            modality="all",
            stage="phase3",
            threshold=0.5,
            method="fixed",
            specificity=0.75,
        )
        figure = figures.expression_figure(self.snapshot, evidence)
        target_index = next(
            i for i, label in enumerate(figure.data[0].x) if "TARGET1" in label
        )
        self.assertTrue(all(row[target_index] == 0 for row in figure.data[0].z))
        constant_only = figures.expression_figure(
            self.snapshot,
            [record for record in evidence if record["target_id"] == "ENSG_TARGET_1"],
        )
        self.assertEqual(
            (constant_only.data[0].zmin, constant_only.data[0].zmax), (-1, 1)
        )

    def test_expression_markers_follow_current_rule(self) -> None:
        records = self.snapshot["records"]
        fixed = figures.expression_figure(
            self.snapshot, records, threshold=0.5, method="fixed"
        )
        specific = figures.expression_figure(
            self.snapshot, records, threshold=0.5, method="specificity"
        )

        def marked_cells(figure):
            marks = next(
                trace for trace in figure.data if trace.name == "Meets expression rule"
            )
            return set(zip(marks.x, marks.y))

        self.assertEqual(
            marked_cells(fixed),
            {
                ("TARGET1 (ENSG_TARGET_1)", "memory B cell"),
                ("TARGET2 (ENSG_TARGET_2)", "naive B cell"),
                ("TARGET2 (ENSG_TARGET_2)", "CD8-positive T cell"),
            },
        )
        self.assertEqual(
            marked_cells(specific),
            {
                ("TARGET1 (ENSG_TARGET_1)", "memory B cell"),
                ("TARGET2 (ENSG_TARGET_2)", "naive B cell"),
                ("TARGET2 (ENSG_TARGET_2)", "CD8-positive T cell"),
            },
        )

    def test_detail_expression_markers_use_applied_rule(self) -> None:
        panel = layout.detail_panel(
            self.rows, "MONDO_RA_TEST", self.snapshot, method="specificity"
        )
        self.assertIn("○ Source cell meets expression rule", str(panel.children[2]))
        graph = panel.children[3].children.children[0]
        self.assertEqual(graph.id, "expression-heatmap")
        self.assertTrue(graph.config["responsive"])
        self.assertEqual(panel.children[1].children[0].data["method"], "specificity")

    def test_expression_columns_group_similar_profiles_and_put_missing_last(
        self,
    ) -> None:
        snapshot = fixture()
        cells = ["CL_B_ONE", "CL_B_TWO", "CL_T_ONE"]
        profiles = {
            "A": [10, 1, 1],
            "B": [1, 10, 10],
            "C": [10, 1, 1],
            "D": [1, 10, 10],
        }
        metadata = {
            (target, cell): {"median": value}
            for target, values in profiles.items()
            for cell, value in zip(cells, values)
        }
        records = [
            {"target_id": target, "target": target}
            for target in ["A", "B", "C", "D", "E"]
        ]
        labels = list(figures.expression_figure(snapshot, records, metadata).data[0].x)
        self.assertEqual(labels[-1], "E (E)")
        self.assertEqual(abs(labels.index("A (A)") - labels.index("C (C)")), 1)
        self.assertEqual(abs(labels.index("B (B)") - labels.index("D (D)")), 1)

    def test_continuous_expression_shows_all_cells(
        self,
    ) -> None:
        row = {
            **self.rows[0],
            "member_cell_ids": list(reversed(self.rows[0]["member_cell_ids"])),
        }
        evidence = figures.evidence_rows(
            self.snapshot,
            row,
            modality="all",
            stage="phase3",
            threshold=0.5,
            method="fixed",
            specificity=0.75,
        )
        figure = figures.expression_figure(self.snapshot, evidence)
        self.assertEqual(
            list(figure.data[0].y),
            ["CD8-positive T cell", "memory B cell", "naive B cell"],
        )

    def test_selection_rejects_stale_click_and_accepts_either_heatmap(self) -> None:
        stale = {"points": [{"customdata": ["OLD", "CELL"]}]}
        valid = {"points": [{"customdata": ["MONDO_RA_TEST", "CL_B_GROUP"]}]}
        visible = ["MONDO_RA_TEST"]
        self.assertIsNone(
            callbacks.resolve_disease_selection("target-heatmap", stale, None, visible)
        )
        self.assertEqual(
            callbacks.resolve_disease_selection("target-heatmap", valid, None, visible),
            "MONDO_RA_TEST",
        )
        self.assertEqual(
            callbacks.resolve_disease_selection("drug-heatmap", valid, None, visible),
            "MONDO_RA_TEST",
        )
        self.assertEqual(
            callbacks.resolve_disease_selection(
                "detail-disease", None, "MONDO_RA_TEST", visible
            ),
            "MONDO_RA_TEST",
        )


class InputTests(unittest.TestCase):
    """空欄と不正値で表示値と計算値がずれない。"""

    def test_empty_and_invalid_thresholds_use_visible_defaults(self) -> None:
        self.assertEqual(callbacks.effective_filters(None, ""), (0.5, 0.5, []))
        minimum, specificity, errors = callbacks.effective_filters("NaN", 1.5)
        self.assertEqual((minimum, specificity), (0.5, 0.5))
        self.assertEqual(len(errors), 2)

    def test_schema1_requires_refresh(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            path.write_text(json.dumps({"schema": 1}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Refresh"):
                snapshot.load_snapshot(path)


class ApplicationPathTests(unittest.TestCase):
    def test_data_paths_keep_the_repository_data_directory(self):
        expected = ASSETS_PATH.parent / "data" / "snapshot.json"
        self.assertEqual(snapshot.DATA_PATH, expected)
        self.assertEqual(refresh.DATA_PATH, expected)

    def test_entrypoint_and_assets_do_not_depend_on_working_directory(self):
        entrypoint = ASSETS_PATH.parent / "app.py"
        previous = Path.cwd()
        try:
            with tempfile.TemporaryDirectory() as directory:
                os.chdir(directory)
                with patch("backend.snapshot.load_snapshot", return_value=None):
                    namespace = runpy.run_path(
                        str(entrypoint), run_name="app_path_test"
                    )
            dash_app = namespace["app"]
            self.assertEqual(Path(dash_app.config.assets_folder), ASSETS_PATH)
            client = dash_app.server.test_client()
            for asset in ("/assets/style.css", "/assets/icons/info.svg"):
                response = client.get(asset)
                try:
                    self.assertEqual(response.status_code, 200)
                finally:
                    response.close()
        finally:
            os.chdir(previous)


class CallbackTests(unittest.TestCase):
    """HTTP 経由で 階層選択、2図、詳細の callback を確認する。"""

    @classmethod
    def setUpClass(cls) -> None:
        cls.snapshot = fixture()
        cls.application = application.create_app(
            cls.snapshot, assets_folder=ASSETS_PATH
        )
        cls.client = cls.application.server.test_client()
        cls.components = {}

        def collect(node):
            if isinstance(node, dict):
                props = node.get("props", {})
                if "id" in props:
                    cls.components[props["id"]] = props
                for value in node.values():
                    collect(value)
            elif isinstance(node, list):
                for value in node:
                    collect(value)

        collect(cls.client.get("/_dash-layout").get_json())

    def _callback_key(self, output_id: str) -> str:
        return next(key for key in self.application.callback_map if output_id in key)

    def _post(self, output_id: str, values: dict, changed: str):
        key = self._callback_key(output_id)
        callback = self.application.callback_map[key]
        self.assertIn(
            changed, [f"{item['id']}.{item['property']}" for item in callback["inputs"]]
        )
        outputs = (
            callback["output"]
            if isinstance(callback["output"], list)
            else [callback["output"]]
        )
        response = self.client.post(
            "/_dash-update-component",
            json={
                "output": key,
                "outputs": [
                    {"id": item.component_id, "property": item.component_property}
                    for item in outputs
                ]
                if isinstance(callback["output"], list)
                else {
                    "id": outputs[0].component_id,
                    "property": outputs[0].component_property,
                },
                "inputs": [
                    {**item, "value": values[item["id"], item["property"]]}
                    for item in callback["inputs"]
                ],
                "state": [
                    {**item, "value": values[item["id"], item["property"]]}
                    for item in callback["state"]
                ],
                "changedPropIds": [changed],
            },
        )
        self.assertEqual(response.status_code, 200, response.data)
        return response.get_json()["response"]

    def _apply(self, values):
        result = self._post("applied-parameters.data", values, "update-button.n_clicks")
        values[("applied-parameters", "data")] = result["applied-parameters"]["data"]

    def _post_applied(self, output_id, values):
        self._apply(values)
        return self._post(output_id, values, "applied-parameters.data")

    def _click_details(self, values, graph):
        selected = self._post("detail-disease.value", values, f"{graph}.clickData")
        values = dict(values)
        values[("detail-disease", "value")] = selected["detail-disease"]["value"]
        return self._post("details.children", values, "detail-disease.value")

    def _values(self) -> dict:
        values = {
            ("measure", "value"): "percent",
            ("modality", "value"): "all",
            ("stage", "value"): "phase3",
            ("method", "value"): "fixed",
            ("threshold", "value"): None,
            ("specificity", "value"): None,
            ("diseases", "value"): ["MONDO_RA_TEST"],
            ("detail-disease", "value"): "MONDO_RA_TEST",
            ("target-heatmap", "clickData"): None,
            ("drug-heatmap", "clickData"): None,
            ("update-button", "n_clicks"): 1,
            ("heatmap-view", "value"): "target",
            ("expanded-cell-groups", "data"): [],
        }
        values[("applied-parameters", "data")] = {
            key: values[(key, "value")] for key in config.PARAMETER_IDS
        }
        return values

    def test_expand_cells_preserves_scale(self):
        values = self._values()
        before = self._post("target-heatmap.figure", values, "applied-parameters.data")
        values[("expanded-cell-groups", "data")] = ["group:CL_B_GROUP"]
        expanded = self._post(
            "target-heatmap.figure", values, "expanded-cell-groups.data"
        )
        for graph in ("target-heatmap", "drug-heatmap"):
            self.assertEqual(
                expanded[graph]["figure"]["layout"]["xaxis"]["tickangle"], -45
            )
            old = before[graph]["figure"]["data"][0]
            new = expanded[graph]["figure"]["data"][0]
            self.assertGreater(len(new["y"]), len(old["y"]))
            self.assertEqual((old["zmin"], old["zmax"]), (new["zmin"], new["zmax"]))
            self.assertEqual(old["z"][0], new["z"][0])
        values[("expanded-cell-groups", "data")] = []
        collapsed = self._post(
            "target-heatmap.figure", values, "expanded-cell-groups.data"
        )
        self.assertEqual(
            before["target-heatmap"]["figure"], collapsed["target-heatmap"]["figure"]
        )
        self.assertEqual(
            figures.heatmap_cell_ids(
                self.snapshot, ["group:CL_B_GROUP", "CL_B_ONE"], []
            ),
            ["group:CL_B_GROUP", "CL_B_ONE"],
        )

    def test_heatmap_button_toggles_and_ignores_rendered_buttons(self):
        button = {
            "type": "heatmap-cell-toggle",
            "kind": "target",
            "cell": "group:CL_B_GROUP",
        }
        for clicks, expanded, expected in (
            (1, [], ["group:CL_B_GROUP"]),
            (1, ["group:CL_B_GROUP"], []),
            (0, [], None),
        ):
            response = self.client.post(
                "/_dash-update-component",
                json={
                    "output": "expanded-cell-groups.data",
                    "outputs": {"id": "expanded-cell-groups", "property": "data"},
                    "inputs": [
                        [{"id": button, "property": "n_clicks", "value": clicks}]
                    ],
                    "state": [
                        {
                            "id": "expanded-cell-groups",
                            "property": "data",
                            "value": expanded,
                        }
                    ],
                    "changedPropIds": [
                        json.dumps(button, sort_keys=True, separators=(",", ":"))
                        + ".n_clicks"
                    ],
                },
            )
            if expected is None:
                self.assertNotIn(
                    "expanded-cell-groups", response.get_json().get("response", {})
                )
            else:
                self.assertEqual(response.status_code, 200)
                self.assertEqual(
                    response.get_json()["response"]["expanded-cell-groups"]["data"],
                    expected,
                )

    def test_expansion_reuses_summary_and_catalog_but_new_filters_recompute(self):
        values = self._values()
        self._post("target-heatmap.figure", values, "applied-parameters.data")
        with (
            patch.object(atlas, "summarize", wraps=atlas.summarize) as summarize,
            patch.object(atlas, "cell_catalog", wraps=atlas.cell_catalog) as catalog,
        ):
            values[("expanded-cell-groups", "data")] = ["group:CL_B_GROUP"]
            self._post("target-heatmap.figure", values, "expanded-cell-groups.data")
            values[("expanded-cell-groups", "data")] = []
            self._post("target-heatmap.figure", values, "expanded-cell-groups.data")
            summarize.assert_not_called()
            catalog.assert_not_called()
            values[("applied-parameters", "data")] = {
                **values[("applied-parameters", "data")],
                "threshold": 123,
            }
            result = self._post(
                "target-heatmap.figure", values, "applied-parameters.data"
            )
            summarize.assert_called_once()
            self.assertEqual(summarize.call_args.args[2], 123)
            self.assertEqual(
                result["target-heatmap"]["figure"]["data"][0]["z"][1][0], 0
            )

    def test_update_applies_parameters_together(
        self,
    ) -> None:
        values = self._values()
        applied = self._post(
            "applied-parameters.data", values, "update-button.n_clicks"
        )
        values[("applied-parameters", "data")] = applied["applied-parameters"]["data"]
        before = self._post("target-heatmap.figure", values, "applied-parameters.data")
        values[("threshold", "value")] = 2
        values[("measure", "value")] = "count"
        values[("heatmap-view", "value")] = "drug"
        status = self._post("update-status.children", values, "threshold.value")
        self.assertIn("not applied", status["update-status"]["children"])
        unchanged = self._post(
            "target-heatmap.figure", values, "applied-parameters.data"
        )
        self.assertEqual(before, unchanged)
        applied = self._post(
            "applied-parameters.data", values, "update-button.n_clicks"
        )
        values[("applied-parameters", "data")] = applied["applied-parameters"]["data"]
        after = self._post("target-heatmap.figure", values, "applied-parameters.data")
        self.assertNotEqual(before, after)
        self.assertEqual(
            after["target-heatmap"]["figure"]["data"][0]["y"],
            ["T cell (group)", "B cell (group)"],
        )
        self.assertIn("CPM ≥ 2", json.dumps(after["matrix-note"], ensure_ascii=False))
        panels = self._post(
            "target-heatmap-panel.hidden", values, "applied-parameters.data"
        )
        self.assertTrue(panels["target-heatmap-panel"]["hidden"])
        selectors = self._post(
            "detail-disease.options", values, "applied-parameters.data"
        )
        self.assertEqual(selectors["detail-disease"]["value"], "MONDO_RA_TEST")
        details = self._post("details.children", values, "detail-disease.value")
        self.assertIn(
            "Drug–target records for rheumatoid arthritis",
            json.dumps(details, ensure_ascii=False),
        )
        status = self._post("update-status.children", values, "applied-parameters.data")
        self.assertEqual(status["update-status"]["children"], "Settings applied.")
        parameter_ids = {
            "threshold",
            "specificity",
            "modality",
            "stage",
            "method",
            "diseases",
            "measure",
            "heatmap-view",
        }
        for output in (
            "target-heatmap.figure",
            "detail-disease.options",
            "details.children",
            "target-heatmap-panel.hidden",
        ):
            callback = self.application.callback_map[self._callback_key(output)]
            self.assertFalse(
                parameter_ids
                & {item["id"] for item in callback["inputs"] + callback["state"]}
            )

    def test_disease_tree_and_search_sync_without_selecting_descendants(self) -> None:
        snapshot = fixture()
        for item_id, name, parents in (
            ("MONDO_0007915", "systemic lupus erythematosus", []),
            ("MONDO_0008383", "rheumatoid arthritis", []),
            ("EFO_0009459", "ACPA-positive rheumatoid arthritis", ["MONDO_0008383"]),
        ):
            snapshot["diseases"].append(
                {"id": item_id, "name": name, "status": "ready", "parent_ids": parents}
            )
            snapshot["records"].append(
                {**snapshot["records"][0], "disease_id": item_id, "disease": name}
            )
        self.application = application.create_app(snapshot, assets_folder=ASSETS_PATH)
        self.client = self.application.server.test_client()
        values = self._values()
        values[("diseases", "value")] = [
            "EFO_0009459",
            "MONDO_0008383",
            "MONDO_0007915",
            "MONDO_0008383",
        ]
        for family in ("MONDO_0007915", "MONDO_0008383", "other"):
            values[(f"disease-family-{family}", "value")] = []
        values[("disease-details-MONDO_0008383", "value")] = []
        synced = self._post("diseases.value", values, "diseases.value")
        expected = ["MONDO_0008383", "EFO_0009459", "MONDO_0007915"]
        self.assertEqual(synced["diseases"]["value"], expected)
        self.assertEqual(
            synced["disease-family-MONDO_0008383"]["value"], ["MONDO_0008383"]
        )
        self.assertEqual(
            synced["disease-details-MONDO_0008383"]["value"], ["EFO_0009459"]
        )
        values[("diseases", "value")] = synced["diseases"]["value"]
        values[("disease-family-MONDO_0008383", "value")] = []
        synced = self._post(
            "diseases.value", values, "disease-family-MONDO_0008383.value"
        )
        self.assertEqual(synced["diseases"]["value"], ["EFO_0009459", "MONDO_0007915"])
        values[("diseases", "value")] = synced["diseases"]["value"]
        values[("disease-family-MONDO_0008383", "value")] = ["MONDO_0008383"]
        synced = self._post(
            "diseases.value", values, "disease-family-MONDO_0008383.value"
        )
        self.assertEqual(synced["diseases"]["value"], expected)
        values[("diseases", "value")] = synced["diseases"]["value"]
        synced = self._post(
            "diseases.value", values, "disease-details-MONDO_0008383.value"
        )
        self.assertEqual(
            synced["diseases"]["value"], ["MONDO_0008383", "MONDO_0007915"]
        )
        values[("diseases", "value")] = list(reversed(synced["diseases"]["value"]))
        figures = self._post_applied("target-heatmap.figure", values)
        for graph in ("target-heatmap", "drug-heatmap"):
            self.assertEqual(
                figures[graph]["figure"]["data"][0]["x"],
                ["rheumatoid arthritis", "systemic lupus erythematosus"],
            )
        values[("diseases", "value")] = []
        cleared = self._post("diseases.value", values, "diseases.value")
        self.assertTrue(all(result["value"] == [] for result in cleared.values()))

    def test_valid_details_survive_filter_changes_without_eager_table(self) -> None:
        values = self._values()
        for changed in (
            "threshold.value",
            "measure.value",
            "modality.value",
        ):
            with self.subTest(changed=changed):
                result = self._post_applied("details.children", values)
                rendered = json.dumps(result)
                self.assertIn(
                    "Drug–target records for rheumatoid arthritis",
                    json.dumps(result, ensure_ascii=False),
                )
                self.assertNotIn('"children": "Original drug"', rendered)
        result = self._post_applied("details.children", values)
        self.assertIn(
            "Drug–target records for rheumatoid arthritis",
            json.dumps(result, ensure_ascii=False),
        )

    def test_heatmap_selection_keeps_disease_wide_details(self) -> None:
        values = self._values()
        self._apply(values)
        for graph in ("target-heatmap", "drug-heatmap"):
            with self.subTest(graph=graph):
                values[(graph, "clickData")] = {
                    "points": [
                        {"customdata": ["MONDO_RA_TEST", "group:" + atlas.T_CELL_ID]}
                    ]
                }
                selected = self._post(
                    "detail-disease.value", values, f"{graph}.clickData"
                )
                self.assertEqual(selected["detail-disease"]["value"], "MONDO_RA_TEST")
                result = self._post_applied("details.children", values)
                self.assertIn(
                    "Drug–target records for rheumatoid arthritis",
                    json.dumps(result, ensure_ascii=False),
                )

    def test_expression_expansion_uses_cached_values_and_applied_rule(self):
        values = {
            ("source-context", "data"): {
                "disease_id": "MONDO_RA_TEST",
                "modality": "all",
                "stage": "phase3",
                "threshold": 0.5,
                "method": "specificity",
                "specificity": 0.75,
            },
            ("expanded-expression-groups", "data"): [],
        }
        collapsed = self._post(
            "expression-heatmap.figure", values, "source-context.data"
        )
        values[("expanded-expression-groups", "data")] = ["group:CL_B_GROUP"]
        with patch.object(
            callbacks,
            "expression_figure",
            side_effect=AssertionError("Expansion recomputed expression"),
        ):
            expanded = self._post(
                "expression-heatmap.figure", values, "expanded-expression-groups.data"
            )
        for result in (collapsed, expanded):
            graph = result["expression-heatmap"]
            self.assertEqual(
                graph["style"]["height"], f"{graph['figure']['layout']['height']}px"
            )
        before = collapsed["expression-heatmap"]["figure"]["data"][0]
        after = expanded["expression-heatmap"]["figure"]["data"][0]
        self.assertEqual(len(before["y"]), 2)
        self.assertEqual(len(after["y"]), 4)
        self.assertEqual(before["zmax"], after["zmax"])
        marks = next(
            trace
            for trace in expanded["expression-heatmap"]["figure"]["data"]
            if trace.get("name") == "Meets expression rule"
        )
        self.assertEqual(set(marks["y"]), {"memory B cell", "naive B cell"})
        self.assertEqual(len(expanded["expression-row-controls"]["children"]), 4)

    def test_source_records_are_visible_and_paged_without_losing_rows(self) -> None:
        self.assertIn("source-records-page.children", self.application.callback_map)
        context = {"disease_id": "MONDO_RA_TEST", "modality": "all", "stage": "phase3"}
        values = {("source-page", "value"): 1, ("source-context", "data"): context}
        records = atlas.filtered_records(self.snapshot, "all", "phase3")
        seen = []
        with patch.object(callbacks, "SOURCE_PAGE_SIZE", 2):
            for page in range(1, math.ceil(len(records) / 2) + 1):
                values[("source-page", "value")] = page
                result = self._post(
                    "source-records-page.children", values, "source-page.value"
                )
                children = result["source-records-page"]["children"]["props"][
                    "children"
                ]
                table_rows = children[1]["props"]["children"]["props"]["children"][1][
                    "props"
                ]["children"]
                self.assertLessEqual(len(table_rows), 2)
                seen.extend(
                    row["props"]["children"][1]["props"]["children"]
                    for row in table_rows
                )
                headers = children[1]["props"]["children"]["props"]["children"][0][
                    "props"
                ]["children"]["props"]["children"]
                labels_in_table = [
                    header["props"]["children"][0]
                    if isinstance(header["props"]["children"], list)
                    else header["props"]["children"]
                    for header in headers
                ]
                self.assertEqual(
                    labels_in_table,
                    [
                        "Canonical drug",
                        "Original drug",
                        "Modality",
                        "Original / canonical stage",
                        "Target",
                        "Action / mechanism",
                        "Open Targets links",
                    ],
                )
        self.assertEqual(
            seen, [f"{record['drug']} ({record['drug_id']})" for record in records]
        )

    def test_layout_has_all_modalities_and_defaults(self) -> None:
        self.assertEqual(self.client.get("/").status_code, 200)
        self.assertEqual(self.components["stage"]["value"], "phase3")
        self.assertEqual(self.components["measure"]["value"], "percent")
        self.assertEqual(
            [option["value"] for option in self.components["measure"]["options"]],
            ["percent", "count"],
        )
        self.assertEqual(
            self.components["applied-parameters"]["data"]["measure"],
            "percent",
        )
        self.assertEqual(self.components["method"]["value"], "specificity")
        self.assertEqual(
            self.components["applied-parameters"]["data"]["method"],
            "specificity",
        )
        self.assertEqual(self.components["specificity"]["value"], 0.5)
        self.assertEqual(
            self.components["applied-parameters"]["data"]["specificity"], 0.5
        )
        self.assertNotIn("level", self.components)
        self.assertNotIn("detail-cell", self.components)
        self.assertNotIn("download-button", self.components)
        self.assertNotIn("download", self.components)
        self.assertFalse(
            any("download" in key for key in self.application.callback_map)
        )
        self.assertNotIn("cells", self.components)
        self.assertNotIn("cells", config.PARAMETER_IDS)
        self.assertEqual(
            [option["value"] for option in self.components["modality"]["options"]],
            ["all", *[value for _, value in atlas.DRUG_TYPE_MODALITIES]],
        )
        layout = self.client.get("/_dash-layout").get_json()
        controls = next(
            node
            for node in layout["props"]["children"]
            if node.get("props", {}).get("className") == "panel controls"
        )
        self.assertEqual(
            controls["props"]["children"][0]["props"]["children"][0], "Settings panel"
        )
        groups = controls["props"]["children"][1]
        scope = groups["props"]["children"][0]
        scope_controls = scope["props"]["children"][1]
        self.assertEqual(scope_controls["type"], "Div")
        self.assertEqual(
            [
                node["props"]["children"][1]["props"]["id"]
                for node in scope_controls["props"]["children"]
            ],
            ["diseases"],
        )
        self.assertEqual(
            [
                group["props"]["children"][0]["props"]["children"]
                for group in groups["props"]["children"]
                if "filter-group" in group["props"].get("className", "").split()
            ],
            ["Comparison scope", "Drug evidence", "Expression criteria", "Display"],
        )
        browser = scope_controls["props"]["children"][0]["props"]["children"][2]
        self.assertEqual(browser["type"], "Div")
        self.assertEqual(browser["props"]["children"][0]["type"], "H4")
        self.assertFalse(self.components["specificity"]["disabled"])

    def test_specificity_input_follows_rule_and_keeps_value(self) -> None:
        values = self._values()
        values[("specificity", "value")] = 0.8
        for method, disabled in (
            ("fixed", True),
            ("specificity", False),
            ("relative", True),
            ("specificity", False),
        ):
            values[("method", "value")] = method
            response = self._post("specificity.disabled", values, "method.value")
            self.assertEqual(response["specificity"], {"disabled": disabled})
        figures = self._post_applied("target-heatmap.figure", values)
        self.assertIn(
            "specificity ≥ 0.8",
            json.dumps(figures["matrix-note"]["children"], ensure_ascii=False),
        )

    def test_help_buttons_describe_unique_tooltips_and_errors_stay_visible(
        self,
    ) -> None:
        values = self._values()
        values[("threshold", "value")] = -1
        figures = self._post_applied("target-heatmap.figure", values)
        note = figures["matrix-note"]["children"]["props"]["children"]
        self.assertEqual(note[0]["props"]["children"], "2 disease–cell combinations")
        self.assertEqual(note[2]["props"]["role"], "alert")
        self.assertIn("must be finite", note[2]["props"]["children"])

        def nodes(item):
            if isinstance(item, dict):
                if "props" in item:
                    yield item
                for value in item.values():
                    yield from nodes(value)
            elif isinstance(item, list):
                for value in item:
                    yield from nodes(value)

        layout = self.client.get("/_dash-layout").get_json()
        selected = self._values()
        selected[("target-heatmap", "clickData")] = {
            "points": [{"customdata": ["MONDO_RA_TEST", "group:CL_B_GROUP"]}]
        }
        details = self._click_details(selected, "target-heatmap")
        surfaces = [
            layout,
            figures["matrix-note"]["children"],
            details["details"]["children"],
        ]
        tips = [
            node["props"]["id"]
            for node in nodes(surfaces)
            if node["props"].get("role") == "tooltip"
        ]
        buttons = [
            node["props"]
            for node in nodes(surfaces)
            if node.get("type") == "Button"
            and node["props"].get("className") == "info-button"
        ]
        self.assertEqual(len(tips), len(set(tips)))
        self.assertEqual(set(tips), {button["aria-describedby"] for button in buttons})
        self.assertTrue(
            all(
                button["type"] == "button" and button["aria-label"]
                for button in buttons
            )
        )

    def test_both_heatmaps_callbacks(self) -> None:
        values = self._values()
        figures = self._post_applied("target-heatmap.figure", values)
        for graph in ("target-heatmap", "drug-heatmap"):
            self.assertEqual(figures[graph]["style"]["minWidth"], "600px")
        self.assertEqual(
            figures["target-heatmap"]["figure"]["data"][0]["x"],
            figures["drug-heatmap"]["figure"]["data"][0]["x"],
        )
        self.assertIn(
            "CPM ≥ 0.5",
            json.dumps(figures["matrix-note"]["children"], ensure_ascii=False),
        )
        self.assertNotIn(
            "specificity ≥",
            json.dumps(figures["matrix-note"]["children"], ensure_ascii=False),
        )

        values[("method", "value")] = "relative"
        relative = self._post_applied("target-heatmap.figure", values)
        self.assertIn(
            "full-reference target median",
            json.dumps(relative["matrix-note"]["children"], ensure_ascii=False),
        )

    def test_heatmap_view_switch_keeps_both_graphs_and_details(self) -> None:
        self.assertEqual(self.components["heatmap-view"]["value"], "target")
        self.assertFalse(self.components["target-heatmap-panel"].get("hidden", False))
        self.assertTrue(self.components["drug-heatmap-panel"]["hidden"])
        values = self._values()
        figures = self._post_applied("target-heatmap.figure", values)
        self.assertIn("figure", figures["target-heatmap"])
        self.assertIn("figure", figures["drug-heatmap"])
        click = {"points": [{"customdata": ["MONDO_RA_TEST", "group:CL_B_GROUP"]}]}
        for view, hidden, graph in (
            ("target", (False, True), "target-heatmap"),
            ("drug", (True, False), "drug-heatmap"),
        ):
            values[("heatmap-view", "value")] = view
            panels = self._post_applied("target-heatmap-panel.hidden", values)
            self.assertEqual(
                (
                    panels["target-heatmap-panel"]["hidden"],
                    panels["drug-heatmap-panel"]["hidden"],
                ),
                hidden,
            )
            case = dict(values)
            case[(graph, "clickData")] = click
            details = self._click_details(case, graph)
            self.assertIn(
                "Relative target expression by cell type",
                str(details["details"]["children"]),
            )

    def test_figures_keep_cell_lineage_order(self) -> None:
        values = self._values()
        figures = self._post_applied("target-heatmap.figure", values)
        self.assertEqual(
            figures["target-heatmap"]["figure"]["data"][0]["y"],
            ["T cell (group)", "B cell (group)"],
        )
        self.assertEqual(
            figures["drug-heatmap"]["figure"]["data"][0]["y"],
            ["T cell (group)", "B cell (group)"],
        )
        selectors = self._post_applied("detail-disease.options", values)
        self.assertEqual(
            [item["value"] for item in selectors["detail-disease"]["options"]],
            ["MONDO_RA_TEST"],
        )

    def test_details_from_both_clicks(self) -> None:
        values = self._values()
        click = {"points": [{"customdata": ["MONDO_RA_TEST", "group:CL_B_GROUP"]}]}
        for heatmap in ("target-heatmap", "drug-heatmap"):
            case = dict(values)
            case[(heatmap, "clickData")] = click
            details = self._click_details(case, heatmap)
            rendered = str(details["details"]["children"])
            self.assertIn("Relative target expression by cell type", rendered)
            self.assertIn("Drug–target records for rheumatoid arthritis", rendered)


if __name__ == "__main__":
    unittest.main()
