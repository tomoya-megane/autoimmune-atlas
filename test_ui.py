"""schema 2 の Dash UI と表示用ロジックを検証する。"""

from __future__ import annotations

import csv
import io
import json
import math
import tempfile
import unittest
from pathlib import Path

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
        "diseases": [{"id": "MONDO_RA_TEST", "name": "rheumatoid arthritis", "status": "ready"}],
        "records": [
            {**common, "drug_id": "DRUG_FORM_1", "drug": "=alpha salt", "canonical_drug_id": "DRUG_ALPHA", "canonical_drug": "alpha", "drug_type": "Protein", "stage": "PHASE_3", "target_id": "ENSG_TARGET_1", "target": "TARGET1", "mechanism": "inhibitor", "action_types": ["INHIBITOR"], "references": [{"source": "PubMed", "ids": ["123"], "urls": ["https://pubmed.ncbi.nlm.nih.gov/123/"]}]},
            {**common, "drug_id": "DRUG_FORM_2", "drug": "alpha hydrate", "canonical_drug_id": "DRUG_ALPHA", "canonical_drug": "alpha", "drug_type": "Protein", "stage": "PHASE_2", "target_id": "ENSG_TARGET_2", "target": "TARGET2", "mechanism": "modulator", "action_types": ["MODULATOR"], "references": []},
            {**common, "drug_id": "DRUG_BETA", "drug": "beta", "canonical_drug_id": "DRUG_BETA", "canonical_drug": "beta", "drug_type": "Antibody", "stage": "APPROVAL", "target_id": None, "target": None, "mechanism": None, "action_types": [], "references": [{"source": "ClinicalTrials", "ids": ["NCT1"], "urls": []}]},
            {**common, "drug_id": "DRUG_GAMMA", "drug": "gamma", "canonical_drug_id": "DRUG_GAMMA", "canonical_drug": "gamma", "drug_type": "Small molecule", "stage": "PHASE_3", "target_id": "ENSG_TARGET_2", "target": "TARGET2", "mechanism": "agonist", "action_types": ["AGONIST"], "references": []},
        ],
        "expression": {
            "ENSG_TARGET_1": [
                {"cell_id": "CL_B_ONE", "cell": "memory B cell", "median": 2.0, "specificity_score": .8, "parent_id": "CL_B_GROUP", "parent": "B cell", "ancestor_ids": ["CL_B_GROUP"]},
                {"cell_id": "CL_B_TWO", "cell": "naive B cell", "median": .1, "specificity_score": .2, "parent_id": "CL_B_GROUP", "parent": "B cell", "ancestor_ids": ["CL_B_GROUP"]},
                {"cell_id": "CL_T_ONE", "cell": "CD8-positive T cell", "median": .2, "specificity_score": .1, "parent_id": "CL_T_PARENT", "parent": "T lymphocyte", "ancestor_ids": [atlas.T_CELL_ID]},
            ],
            "ENSG_TARGET_2": [
                {"cell_id": "CL_B_ONE", "cell": "memory B cell", "median": None, "specificity_score": None, "parent_id": "CL_B_GROUP", "parent": "B cell", "ancestor_ids": ["CL_B_GROUP"]},
                {"cell_id": "CL_B_TWO", "cell": "naive B cell", "median": .6, "specificity_score": .75, "parent_id": "CL_B_GROUP", "parent": "B cell", "ancestor_ids": ["CL_B_GROUP"]},
                {"cell_id": "CL_T_ONE", "cell": "CD8-positive T cell", "median": 1.0, "specificity_score": .7, "parent_id": "CL_T_PARENT", "parent": "T lymphocyte", "ancestor_ids": [atlas.T_CELL_ID]},
            ],
        },
    }


def summary(**overrides) -> dict:
    row = {
        "disease_id": "D1", "disease": "Disease one", "cell_id": "C1", "cell": "B cell",
        "count": 1, "percent": 20.0, "denominator": 5, "unknown": 0,
        "drug_count": 1, "drug_percent": 25.0, "drug_denominator": 4, "unknown_drugs": 0,
        "mapped_drugs": 4, "total_drugs": 4, "unmapped_drugs": 0, "status": "complete",
        "member_cell_ids": ["C1"], "records": [],
    }
    row.update(overrides)
    return row


class FigureTests(unittest.TestCase):
    """2種類の指標で 0、下限、欠測を独立に扱う。"""

    def test_target_and_drug_lower_bounds_are_independent(self) -> None:
        row = summary(count=0, percent=0.0, drug_count=1, drug_percent=25.0, unknown_drugs=1)
        self.assertEqual(app._display_value(row, "count"), "0")
        self.assertEqual(app._display_value(row, "count", "drug"), "≥1")
        self.assertEqual(app._display_value(row, "percent"), "0.0%")
        self.assertEqual(app._display_value(row, "percent", "drug"), "≥25.0%")

    def test_unmapped_drug_does_not_change_percent_lower_bound(self) -> None:
        row = summary(count=2, percent=40.0, drug_count=1, drug_percent=25.0, unmapped_drugs=1, total_drugs=5, status="partial")
        self.assertEqual(app._display_value(row, "count"), "≥2")
        self.assertEqual(app._display_value(row, "percent"), "40.0%")
        self.assertEqual(app._display_value(row, "count", "drug"), "≥1")
        self.assertEqual(app._display_value(row, "percent", "drug"), "25.0%")

    def test_two_figures_have_aligned_axes_and_cross_for_na(self) -> None:
        rows = [summary(), summary(cell_id="C2", cell="T cell", count=None, percent=None, drug_count=None, drug_percent=None, unknown=1, unknown_drugs=1, status="partial")]
        target = app.build_figure(rows, ["D1"], ["C1", "C2"], "count", "target")
        drug = app.build_figure(rows, ["D1"], ["C1", "C2"], "count", "drug")
        self.assertEqual(list(target.data[0].x), list(drug.data[0].x))
        self.assertEqual(list(target.data[0].y), list(drug.data[0].y))
        self.assertIn("Unknown / missing", [trace.name for trace in target.data])
        missing = next(trace for trace in target.data if trace.name == "Unknown / missing")
        self.assertIn("Denominator: 5", missing.hovertext[0])
        self.assertIn("Unknown in denominator: 1", missing.hovertext[0])
        self.assertIn("Mapped / all canonical drugs: 4 / 4", missing.hovertext[0])
        self.assertIn("Status: partial", missing.hovertext[0])
        self.assertFalse(any(getattr(trace.marker, "symbol", None) == "square-open" for trace in target.data if trace.type == "scatter"))


class EvidenceTests(unittest.TestCase):
    """詳細と CSV が陽性以外の元記録も保持する。"""

    def setUp(self) -> None:
        self.snapshot = fixture()
        self.rows = app.visible_rows(self.snapshot, "all", .5, ["MONDO_RA_TEST"], ["CL_B_GROUP"])

    def test_continuous_expression_uses_log2_color_and_raw_cpm_hover(self) -> None:
        evidence = app.evidence_rows(self.snapshot, self.rows[0], modality="all", stage="phase3", threshold=.5, method="fixed", specificity=.75)
        figure = app.expression_figure(self.snapshot, self.rows[0], evidence)
        targets = list(figure.data[0].y)
        cells = list(figure.data[0].x)
        target_index = next(i for i, label in enumerate(targets) if "TARGET1" in label)
        cell_index = cells.index("memory B cell")
        self.assertAlmostEqual(figure.data[0].z[target_index][cell_index], math.log2(3))
        self.assertIn("Raw CPM: 2", figure.data[0].hovertext[target_index][cell_index])
        self.assertIn("Missing expression", [trace.name for trace in figure.data])

    def test_all_states_and_source_cells_are_exported(self) -> None:
        exported = app.export_rows(self.rows, "count", snapshot=self.snapshot)
        self.assertTrue({"positive", "negative", "unknown", "unmapped"}.issubset({row["support_state"] for row in exported}))
        self.assertEqual({row["row_type"] for row in exported}, {"summary", "source_record", "expression_evidence"})
        self.assertTrue(all(row["group_cell_id"] == "CL_B_GROUP" for row in exported if row["row_type"] != "source_record"))
        self.assertIn("CL_B_ONE", {row["evidence_cell_id"] for row in exported})
        self.assertEqual({row["minimum_cpm"] for row in exported}, {.5})
        self.assertTrue(all("target_count" in row and "drug_count" in row for row in exported))
        self.assertEqual({row["expression_source"] for row in exported}, {"Tabula Sapiens"})
        content = atlas.to_csv(exported)
        parsed = list(csv.DictReader(io.StringIO(content)))
        risky = next(row for row in parsed if row["drug_id"] == "DRUG_FORM_1")
        self.assertEqual(risky["drug"], "'=alpha salt")

    def test_selection_rejects_stale_click_and_accepts_either_heatmap(self) -> None:
        stale = {"points": [{"customdata": ["OLD", "CELL"]}]}
        valid = {"points": [{"customdata": ["MONDO_RA_TEST", "CL_B_GROUP"]}]}
        self.assertIsNone(app.resolve_selection("target-heatmap", stale, None, None, self.rows))
        self.assertEqual(app.resolve_selection("target-heatmap", valid, None, None, self.rows), ("MONDO_RA_TEST", "CL_B_GROUP"))
        self.assertEqual(app.resolve_selection("drug-heatmap", valid, None, None, self.rows), ("MONDO_RA_TEST", "CL_B_GROUP"))
        self.assertEqual(app.resolve_selection("detail-cell", None, "MONDO_RA_TEST", "CL_B_GROUP", self.rows), ("MONDO_RA_TEST", "CL_B_GROUP"))


class InputTests(unittest.TestCase):
    """空欄と不正値で表示値と計算値がずれない。"""

    def test_empty_and_invalid_thresholds_use_visible_defaults(self) -> None:
        self.assertEqual(app.effective_filters(None, ""), (.5, .75, []))
        minimum, specificity, errors = app.effective_filters("NaN", 1.5)
        self.assertEqual((minimum, specificity), (.5, .75))
        self.assertEqual(len(errors), 2)

    def test_schema1_requires_refresh(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            path.write_text(json.dumps({"schema": 1}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Refresh"):
                app.load_snapshot(path)


class CallbackTests(unittest.TestCase):
    """HTTP 経由で level、2図、詳細、CSV の callback を確認する。"""

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
        outputs = callback["output"] if isinstance(callback["output"], list) else [callback["output"]]
        response = self.client.post("/_dash-update-component", json={
            "output": key,
            "outputs": [{"id": item.component_id, "property": item.component_property} for item in outputs] if isinstance(callback["output"], list) else {"id": outputs[0].component_id, "property": outputs[0].component_property},
            "inputs": [{**item, "value": values[item["id"], item["property"]]} for item in callback["inputs"]],
            "state": [{**item, "value": values[item["id"], item["property"]]} for item in callback["state"]],
            "changedPropIds": [changed],
        })
        self.assertEqual(response.status_code, 200, response.data)
        return response.get_json()["response"]

    def _values(self) -> dict:
        return {
            ("measure", "value"): "percent", ("modality", "value"): "all", ("stage", "value"): "phase3",
            ("method", "value"): "fixed", ("threshold", "value"): None, ("specificity", "value"): None,
            ("level", "value"): "group", ("diseases", "value"): ["MONDO_RA_TEST"], ("cells", "value"): ["CL_B_GROUP"],
            ("detail-disease", "value"): "MONDO_RA_TEST", ("detail-cell", "value"): "CL_B_GROUP",
            ("target-heatmap", "clickData"): None, ("drug-heatmap", "clickData"): None,
            ("download-button", "n_clicks"): 1,
        }

    def test_layout_has_all_modalities_and_defaults(self) -> None:
        self.assertEqual(self.client.get("/").status_code, 200)
        self.assertEqual(self.components["stage"]["value"], "phase3")
        self.assertEqual(self.components["method"]["value"], "fixed")
        self.assertEqual(self.components["level"]["value"], "group")
        self.assertEqual([option["value"] for option in self.components["modality"]["options"]], ["all", *[value for _, value in atlas.DRUG_TYPE_MODALITIES]])

    def test_level_and_both_heatmaps_callbacks(self) -> None:
        values = self._values()
        values[("level", "value")] = "cell"
        level_response = self._post("cells.options", values, "level.value")
        self.assertIn("CL_B_ONE", {item["value"] for item in level_response["cells"]["options"]})
        self.assertEqual(set(level_response["cells"]["value"]), {"CL_B_ONE", "CL_B_TWO"})

        values[("cells", "value")] = ["CL_B_ONE"]
        values[("level", "value")] = "group"
        group_response = self._post("cells.options", values, "level.value")
        self.assertEqual(group_response["cells"]["value"], ["CL_B_GROUP"])

        values = self._values()
        figures = self._post("target-heatmap.figure", values, "measure.value")
        self.assertEqual(figures["target-heatmap"]["figure"]["data"][0]["x"], figures["drug-heatmap"]["figure"]["data"][0]["x"])
        self.assertIn("CPM > 0.5", figures["matrix-note"]["children"])
        self.assertNotIn("specificity ≥", figures["matrix-note"]["children"])

        values[("method", "value")] = "relative"
        relative = self._post("target-heatmap.figure", values, "method.value")
        self.assertIn("full-reference target median", relative["matrix-note"]["children"])

    def test_details_from_both_clicks_and_csv(self) -> None:
        values = self._values()
        click = {"points": [{"customdata": ["MONDO_RA_TEST", "CL_B_GROUP"]}]}
        for heatmap in ("target-heatmap", "drug-heatmap"):
            case = dict(values); case[(heatmap, "clickData")] = click
            details = self._post("details.children", case, f"{heatmap}.clickData")
            rendered = str(details["details"]["children"])
            self.assertIn("Continuous expression", rendered)
            self.assertIn("All filtered source records", rendered)

        downloaded = self._post("download.data", values, "download-button.n_clicks")
        exported = list(csv.DictReader(io.StringIO(downloaded["download"]["data"]["content"].lstrip("\ufeff"))))
        self.assertTrue(exported)
        self.assertEqual({row["minimum_cpm"] for row in exported}, {"0.5"})
        self.assertEqual({row["specificity_threshold"] for row in exported}, {"0.75"})
        self.assertEqual({row["stage_filter"] for row in exported}, {"phase3"})
        self.assertTrue({"positive", "negative", "unknown", "unmapped"}.issubset({row["support_state"] for row in exported}))
        self.assertEqual({row["row_type"] for row in exported}, {"summary", "source_record", "expression_evidence"})


if __name__ == "__main__":
    unittest.main()
