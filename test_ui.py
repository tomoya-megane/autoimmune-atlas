"""schema 2 の Dash UI と表示用ロジックを検証する。"""

from __future__ import annotations

import csv
import io
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import app
import atlas


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


class FigureTests(unittest.TestCase):
    """2種類の指標で 0、下限、欠測を独立に扱う。"""

    def test_target_and_drug_lower_bounds_are_independent(self) -> None:
        row = summary(
            count=0, percent=0.0, drug_count=1, drug_percent=25.0, unknown_drugs=1
        )
        self.assertEqual(app._display_value(row, "count"), "0")
        self.assertEqual(app._display_value(row, "count", "drug"), "≥1")
        self.assertEqual(app._display_value(row, "percent"), "0.0%")
        self.assertEqual(app._display_value(row, "percent", "drug"), "≥25.0%")

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
        self.assertEqual(app._display_value(row, "count"), "≥2")
        self.assertEqual(app._display_value(row, "percent"), "40.0%")
        self.assertEqual(app._display_value(row, "count", "drug"), "≥1")
        self.assertEqual(app._display_value(row, "percent", "drug"), "25.0%")

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
            figure = app.build_figure([row], ["D1"], ["C1"], measure, kind)
            self.assertEqual(figure.data[0].text[0][0], label)
            self.assertEqual(figure.data[0].z[0][0], value)
            self.assertIn("≥" + label, figure.data[0].hovertext[0][0])
            self.assertIn("≥ is a lower bound", figure.data[0].hovertext[0][0])

    def test_two_figures_have_aligned_axes_and_cross_for_na(self) -> None:
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
                status="partial",
            ),
        ]
        target = app.build_figure(rows, ["D1"], ["C1", "C2"], "count", "target")
        drug = app.build_figure(rows, ["D1"], ["C1", "C2"], "count", "drug")
        self.assertEqual(list(target.data[0].x), list(drug.data[0].x))
        self.assertEqual(list(target.data[0].y), list(drug.data[0].y))
        self.assertIn("No value (see hover)", [trace.name for trace in target.data])
        missing = next(
            trace for trace in target.data if trace.name == "No value (see hover)"
        )
        self.assertIn("Denominator: 5", missing.hovertext[0])
        self.assertIn("Unknown in denominator: 1", missing.hovertext[0])
        self.assertIn("Mapped / all canonical drugs: 4 / 4", missing.hovertext[0])
        self.assertIn("Status: partial", missing.hovertext[0])
        self.assertFalse(
            any(
                getattr(trace.marker, "symbol", None) == "square-open"
                for trace in target.data
                if trace.type == "scatter"
            )
        )

    def test_no_value_hover_distinguishes_unavailable_from_empty_denominator(
        self,
    ) -> None:
        unavailable = summary(
            status="unavailable", count=None, percent=None, denominator=0
        )
        empty = summary(percent=None, denominator=0)
        self.assertIn(
            "Disease data unavailable",
            app._hover_text(unavailable, "percent", "target"),
        )
        self.assertNotIn(
            "No eligible items", app._hover_text(unavailable, "percent", "target")
        )
        self.assertIn(
            "No eligible items in the percentage denominator",
            app._hover_text(empty, "percent", "target"),
        )


class EvidenceTests(unittest.TestCase):
    """詳細と CSV が陽性以外の元記録も保持する。"""

    def setUp(self) -> None:
        self.snapshot = fixture()
        self.rows = app.visible_rows(
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
        panel = app.detail_panel(all_rows, "MONDO_RA_TEST", self.snapshot)
        self.assertNotIn("Drug–target pairs", str(panel))
        self.assertNotIn("Targets meeting rule in any source cell", str(panel))
        self.assertNotIn("rheumatoid arthritis · all source cell types", str(panel))
        self.assertIn("Filtered drug records for rheumatoid arthritis", str(panel))
        self.assertIn("Filtered drug records", str(panel.children[0]))
        self.assertIn("Relative expression", str(panel.children[2]))
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
        evidence = app.evidence_rows(
            self.snapshot,
            self.rows[0],
            modality="all",
            stage="phase3",
            threshold=0.5,
            method="fixed",
            specificity=0.75,
        )
        figure = app.expression_figure(self.snapshot, evidence)
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
        self.assertEqual((figure.data[0].zmin, figure.data[0].zmax), (-3, 3))
        self.assertEqual(figure.data[0].colorbar.lenmode, "pixels")
        self.assertEqual(figure.data[0].colorbar.len, 240)
        self.assertEqual(figure.data[0].colorbar.y, 1)
        self.assertEqual(figure.data[0].colorbar.yanchor, "top")
        missing = next(
            trace for trace in figure.data if trace.name == "Missing expression"
        )
        self.assertEqual(list(missing.x), ["TARGET2 (ENSG_TARGET_2)"])
        self.assertEqual(list(missing.y), ["memory B cell"])
        self.assertEqual(figure.layout.xaxis.side, "top")
        self.assertEqual(list(figure.layout.yaxis.range), [len(cells) - 0.5, -0.5])

    def test_constant_target_expression_has_zero_z_score(self) -> None:
        for item in self.snapshot["expression"]["ENSG_TARGET_1"]:
            item["median"] = 2
        evidence = app.evidence_rows(
            self.snapshot,
            self.rows[0],
            modality="all",
            stage="phase3",
            threshold=0.5,
            method="fixed",
            specificity=0.75,
        )
        figure = app.expression_figure(self.snapshot, evidence)
        target_index = next(
            i for i, label in enumerate(figure.data[0].x) if "TARGET1" in label
        )
        self.assertTrue(all(row[target_index] == 0 for row in figure.data[0].z))

    def test_expression_markers_follow_current_rule(self) -> None:
        records = self.snapshot["records"]
        fixed = app.expression_figure(
            self.snapshot, records, threshold=0.5, method="fixed"
        )
        specific = app.expression_figure(
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
            },
        )

    def test_detail_expression_markers_use_applied_rule(self) -> None:
        panel = app.detail_panel(
            self.rows, "MONDO_RA_TEST", self.snapshot, method="specificity"
        )
        self.assertIn("○ Meets expression rule", str(panel.children[2]))
        figure = panel.children[3].children.figure
        marks = next(
            trace for trace in figure.data if trace.name == "Meets expression rule"
        )
        self.assertEqual(len(marks.x), 2)

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
        labels = list(app.expression_figure(snapshot, records, metadata).data[0].x)
        self.assertEqual(labels[-1], "E (E)")
        self.assertEqual(abs(labels.index("A (A)") - labels.index("C (C)")), 1)
        self.assertEqual(abs(labels.index("B (B)") - labels.index("D (D)")), 1)

    def test_continuous_expression_shows_all_cells_and_csv_keeps_member_order(
        self,
    ) -> None:
        row = {
            **self.rows[0],
            "member_cell_ids": list(reversed(self.rows[0]["member_cell_ids"])),
        }
        evidence = app.evidence_rows(
            self.snapshot,
            row,
            modality="all",
            stage="phase3",
            threshold=0.5,
            method="fixed",
            specificity=0.75,
        )
        figure = app.expression_figure(self.snapshot, evidence)
        self.assertEqual(
            list(figure.data[0].y),
            ["CD8-positive T cell", "memory B cell", "naive B cell"],
        )
        exported = app.export_rows([row], "count", snapshot=self.snapshot)
        self.assertEqual(
            [
                item["evidence_cell_id"]
                for item in exported
                if item["row_type"] == "expression_evidence"
            ][:2],
            ["CL_B_ONE", "CL_B_TWO"],
        )

    def test_all_states_and_source_cells_are_exported(self) -> None:
        exported = app.export_rows(self.rows, "count", snapshot=self.snapshot)
        self.assertTrue(
            {"positive", "negative", "unknown", "unmapped"}.issubset(
                {row["support_state"] for row in exported}
            )
        )
        self.assertEqual(
            {row["row_type"] for row in exported},
            {"summary", "source_record", "expression_evidence"},
        )
        self.assertTrue(
            all(
                row["group_cell_id"] == "CL_B_GROUP"
                for row in exported
                if row["row_type"] != "source_record"
            )
        )
        self.assertIn("CL_B_ONE", {row["evidence_cell_id"] for row in exported})
        self.assertEqual({row["minimum_cpm"] for row in exported}, {0.5})
        self.assertTrue(
            all("target_count" in row and "drug_count" in row for row in exported)
        )
        self.assertEqual(
            {row["expression_source"] for row in exported}, {"Tabula Sapiens"}
        )
        content = atlas.to_csv(exported)
        parsed = list(csv.DictReader(io.StringIO(content)))
        risky = next(row for row in parsed if row["drug_id"] == "DRUG_FORM_1")
        self.assertEqual(risky["drug"], "'=alpha salt")

    def test_selection_rejects_stale_click_and_accepts_either_heatmap(self) -> None:
        stale = {"points": [{"customdata": ["OLD", "CELL"]}]}
        valid = {"points": [{"customdata": ["MONDO_RA_TEST", "CL_B_GROUP"]}]}
        visible = ["MONDO_RA_TEST"]
        self.assertIsNone(
            app.resolve_disease_selection("target-heatmap", stale, None, visible)
        )
        self.assertEqual(
            app.resolve_disease_selection("target-heatmap", valid, None, visible),
            "MONDO_RA_TEST",
        )
        self.assertEqual(
            app.resolve_disease_selection("drug-heatmap", valid, None, visible),
            "MONDO_RA_TEST",
        )
        self.assertEqual(
            app.resolve_disease_selection(
                "detail-disease", None, "MONDO_RA_TEST", visible
            ),
            "MONDO_RA_TEST",
        )


class InputTests(unittest.TestCase):
    """空欄と不正値で表示値と計算値がずれない。"""

    def test_empty_and_invalid_thresholds_use_visible_defaults(self) -> None:
        self.assertEqual(app.effective_filters(None, ""), (0.5, 0.75, []))
        minimum, specificity, errors = app.effective_filters("NaN", 1.5)
        self.assertEqual((minimum, specificity), (0.5, 0.75))
        self.assertEqual(len(errors), 2)

    def test_schema1_requires_refresh(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            path.write_text(json.dumps({"schema": 1}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Refresh"):
                app.load_snapshot(path)


class CallbackTests(unittest.TestCase):
    """HTTP 経由で 階層選択、2図、詳細、CSV の callback を確認する。"""

    @classmethod
    def setUpClass(cls) -> None:
        cls.snapshot = fixture()
        cls.application = app.create_app(cls.snapshot)
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
            ("cells", "value"): ["group:CL_B_GROUP"],
            ("detail-disease", "value"): "MONDO_RA_TEST",
            ("target-heatmap", "clickData"): None,
            ("drug-heatmap", "clickData"): None,
            ("download-button", "n_clicks"): 1,
            ("update-button", "n_clicks"): 1,
            ("heatmap-view", "value"): "target",
        }
        values[("applied-parameters", "data")] = {
            key: values[(key, "value")] for key in app.PARAMETER_IDS
        }
        return values

    def test_update_applies_parameters_together_and_csv_uses_displayed_settings(
        self,
    ) -> None:
        values = self._values()
        applied = self._post(
            "applied-parameters.data", values, "update-button.n_clicks"
        )
        values[("applied-parameters", "data")] = applied["applied-parameters"]["data"]
        before = self._post("target-heatmap.figure", values, "applied-parameters.data")
        values[("threshold", "value")] = 2
        values[("cells", "value")] = ["CL_B_ONE"]
        values[("measure", "value")] = "count"
        values[("heatmap-view", "value")] = "drug"
        status = self._post("update-status.children", values, "threshold.value")
        self.assertIn("not applied", status["update-status"]["children"])
        unchanged = self._post(
            "target-heatmap.figure", values, "applied-parameters.data"
        )
        self.assertEqual(before, unchanged)
        exported = self._post("download.data", values, "download-button.n_clicks")
        rows = list(
            csv.DictReader(
                io.StringIO(exported["download"]["data"]["content"].lstrip("\ufeff"))
            )
        )
        self.assertEqual({row["minimum_cpm"] for row in rows}, {"0.5"})
        self.assertEqual(
            {row["group_cell_id"] for row in rows if row["row_type"] == "summary"},
            {"group:CL_B_GROUP"},
        )
        applied = self._post(
            "applied-parameters.data", values, "update-button.n_clicks"
        )
        values[("applied-parameters", "data")] = applied["applied-parameters"]["data"]
        after = self._post("target-heatmap.figure", values, "applied-parameters.data")
        self.assertNotEqual(before, after)
        self.assertEqual(
            after["target-heatmap"]["figure"]["data"][0]["y"], ["memory B cell"]
        )
        self.assertIn("CPM > 2", json.dumps(after["matrix-note"], ensure_ascii=False))
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
            "Filtered drug records for rheumatoid arthritis",
            json.dumps(details, ensure_ascii=False),
        )
        exported = self._post("download.data", values, "download-button.n_clicks")
        rows = list(
            csv.DictReader(
                io.StringIO(exported["download"]["data"]["content"].lstrip("\ufeff"))
            )
        )
        self.assertEqual({row["minimum_cpm"] for row in rows}, {"2.0"})
        self.assertEqual({row["measure"] for row in rows}, {"count"})
        status = self._post("update-status.children", values, "applied-parameters.data")
        self.assertEqual(status["update-status"]["children"], "Settings applied.")
        parameter_ids = {
            "threshold",
            "specificity",
            "modality",
            "stage",
            "method",
            "diseases",
            "cells",
            "measure",
            "heatmap-view",
        }
        for output in (
            "target-heatmap.figure",
            "detail-disease.options",
            "details.children",
            "target-heatmap-panel.hidden",
            "download.data",
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
        self.application = app.create_app(snapshot)
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
            "cells.value",
        ):
            with self.subTest(changed=changed):
                result = self._post_applied("details.children", values)
                rendered = json.dumps(result)
                self.assertIn(
                    "Filtered drug records for rheumatoid arthritis",
                    json.dumps(result, ensure_ascii=False),
                )
                self.assertNotIn('"children": "Original drug"', rendered)
        values[("cells", "value")] = []
        result = self._post_applied("details.children", values)
        self.assertIn(
            "Filtered drug records for rheumatoid arthritis",
            json.dumps(result, ensure_ascii=False),
        )

    def test_heatmap_selection_keeps_disease_wide_details(self) -> None:
        values = self._values()
        values[("cells", "value")] = ["group:" + atlas.T_CELL_ID, "group:CL_B_GROUP"]
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
                    "Filtered drug records for rheumatoid arthritis",
                    json.dumps(result, ensure_ascii=False),
                )

    def test_source_records_are_visible_and_paged_without_losing_rows(self) -> None:
        self.assertIn("source-records-page.children", self.application.callback_map)
        context = {"disease_id": "MONDO_RA_TEST", "modality": "all", "stage": "phase3"}
        values = {("source-page", "value"): 1, ("source-context", "data"): context}
        records = atlas.filtered_records(self.snapshot, "all", "phase3")
        seen = []
        with patch.object(app, "SOURCE_PAGE_SIZE", 2):
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
                labels_in_table = [header["props"]["children"] for header in headers]
                self.assertEqual(
                    labels_in_table,
                    [
                        "Canonical drug",
                        "Original drug",
                        "Modality",
                        "Record stage / highest disease-specific drug stage",
                        "Target",
                        "Action / mechanism",
                        "Sources",
                    ],
                )
        self.assertEqual(
            seen, [f"{record['drug']} ({record['drug_id']})" for record in records]
        )

    def test_layout_has_all_modalities_and_defaults(self) -> None:
        self.assertEqual(self.client.get("/").status_code, 200)
        self.assertEqual(self.components["stage"]["value"], "phase3")
        self.assertEqual(self.components["method"]["value"], "fixed")
        self.assertNotIn("level", self.components)
        self.assertNotIn("detail-cell", self.components)
        self.assertTrue(
            all(
                value.startswith("group:")
                for value in self.components["cells"]["value"]
            )
        )
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
        groups = controls["props"]["children"][0]
        scope = groups["props"]["children"][2]
        scope_controls = scope["props"]["children"][1]
        self.assertEqual(scope_controls["type"], "Div")
        self.assertEqual(
            [
                node["props"]["children"][1]["props"]["id"]
                for node in scope_controls["props"]["children"]
            ],
            ["diseases", "cells"],
        )
        self.assertEqual(
            [
                group["props"]["children"][0]["props"]["children"]
                for group in groups["props"]["children"]
            ],
            ["Drug evidence", "Expression criteria", "Comparison scope", "Display"],
        )
        self.assertTrue(self.components["specificity"]["disabled"])

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
        self.assertEqual(note[0]["props"]["children"], "1 disease–cell combinations")
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

    def test_cell_tree_and_search_sync_independent_group_and_source_selection(
        self,
    ) -> None:
        values = self._values()
        for group in app._cell_selection_catalog(self.snapshot):
            for section in group["sections"]:
                values[(section["id"], "value")] = []
        values[("cells", "value")] = ["CL_B_ONE", "group:CL_B_GROUP"]
        synced = self._post("cells.value", values, "cells.value")
        self.assertEqual(synced["cells"]["value"], ["group:CL_B_GROUP", "CL_B_ONE"])
        self.assertEqual(synced["cell-group-CL_B_GROUP"]["value"], ["group:CL_B_GROUP"])
        self.assertEqual(synced["cell-details-CL_B_GROUP"]["value"], ["CL_B_ONE"])
        values[("cells", "value")] = synced["cells"]["value"]
        removed = self._post("cells.value", values, "cell-group-CL_B_GROUP.value")
        self.assertEqual(removed["cells"]["value"], ["CL_B_ONE"])
        values[("cells", "value")] = []
        cleared = self._post("cells.value", values, "cells.value")
        self.assertTrue(all(result["value"] == [] for result in cleared.values()))

    def test_both_heatmaps_callbacks(self) -> None:
        values = self._values()
        figures = self._post_applied("target-heatmap.figure", values)
        self.assertEqual(
            figures["target-heatmap"]["figure"]["data"][0]["x"],
            figures["drug-heatmap"]["figure"]["data"][0]["x"],
        )
        self.assertIn(
            "CPM > 0.5",
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
                "Relative expression of targets across all source cell types",
                str(details["details"]["children"]),
            )

    def test_figures_ignore_reversed_selection_order(self) -> None:
        values = self._values()
        values[("cells", "value")] = ["group:CL_B_GROUP", "group:" + atlas.T_CELL_ID]
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

    def test_mixed_selection_csv_preserves_row_ids_and_group_evidence(self) -> None:
        values = self._values()
        values[("cells", "value")] = ["group:CL_B_GROUP", "CL_B_ONE"]
        figures = self._post_applied("target-heatmap.figure", values)
        self.assertEqual(
            figures["target-heatmap"]["figure"]["data"][0]["y"],
            ["B cell (group)", "memory B cell"],
        )
        downloaded = self._post("download.data", values, "download-button.n_clicks")
        exported = list(
            csv.DictReader(
                io.StringIO(downloaded["download"]["data"]["content"].lstrip("\ufeff"))
            )
        )
        summaries = {
            row["group_cell_id"]: row
            for row in exported
            if row["row_type"] == "summary"
        }
        self.assertEqual(
            {
                key: (row["ontology_cell_id"], row["level"])
                for key, row in summaries.items()
            },
            {
                "group:CL_B_GROUP": ("CL_B_GROUP", "group"),
                "CL_B_ONE": ("CL_B_ONE", "cell"),
            },
        )
        evidence = [row for row in exported if row["row_type"] == "expression_evidence"]
        self.assertEqual(
            {row["evidence_cell_id"] for row in evidence if row["level"] == "group"},
            {"CL_B_ONE", "CL_B_TWO"},
        )
        self.assertEqual(
            {row["evidence_cell_id"] for row in evidence if row["level"] == "cell"},
            {"CL_B_ONE"},
        )
        sources = [row for row in exported if row["row_type"] == "source_record"]
        self.assertEqual(
            len(sources), len(atlas.filtered_records(self.snapshot, "all", "phase3"))
        )

    def test_details_from_both_clicks_and_csv(self) -> None:
        values = self._values()
        click = {"points": [{"customdata": ["MONDO_RA_TEST", "group:CL_B_GROUP"]}]}
        for heatmap in ("target-heatmap", "drug-heatmap"):
            case = dict(values)
            case[(heatmap, "clickData")] = click
            details = self._click_details(case, heatmap)
            rendered = str(details["details"]["children"])
            self.assertIn(
                "Relative expression of targets across all source cell types", rendered
            )
            self.assertIn("Filtered drug records for rheumatoid arthritis", rendered)

        downloaded = self._post("download.data", values, "download-button.n_clicks")
        exported = list(
            csv.DictReader(
                io.StringIO(downloaded["download"]["data"]["content"].lstrip("\ufeff"))
            )
        )
        self.assertTrue(exported)
        self.assertEqual({row["minimum_cpm"] for row in exported}, {"0.5"})
        self.assertEqual({row["specificity_threshold"] for row in exported}, {"0.75"})
        self.assertEqual({row["stage_filter"] for row in exported}, {"phase3"})
        self.assertTrue(
            {"positive", "negative", "unknown", "unmapped"}.issubset(
                {row["support_state"] for row in exported}
            )
        )
        self.assertEqual(
            {row["row_type"] for row in exported},
            {"summary", "source_record", "expression_evidence"},
        )


if __name__ == "__main__":
    unittest.main()
