"""Dash UI の表示用ロジックを検証する。"""

from __future__ import annotations

import csv
import io
import unittest

import app
import atlas


class FigureTests(unittest.TestCase):
    """0、下限値、欠測が同じ色や文字にならないことを検証する。"""

    def test_zero_partial_and_missing_are_distinct(self) -> None:
        rows = [
            _summary("C1", "B cell", 0, 0.0, "complete"),
            _summary("C2", "T cell", 2, 40.0, "partial", unknown=1),
            _summary("C3", "Monocyte", None, None, "partial", unmapped=1),
        ]

        figure = app.build_figure(rows, ["D1"], ["C1", "C2", "C3"], "count")

        heatmap = figure.data[0]
        self.assertEqual(list(heatmap.z[0]), [0])
        self.assertEqual(list(heatmap.z[1]), [2])
        self.assertEqual(list(heatmap.text[0]), ["0"])
        self.assertEqual(list(heatmap.text[1]), ["≥2"])
        self.assertIsNone(heatmap.z[2][0])
        self.assertIn("Unknown / missing", [trace.name for trace in figure.data])
        self.assertIn("Drugs without mapped targets: 1", heatmap.hovertext[2][0])
        missing = next(trace for trace in figure.data if trace.name == "Unknown / missing")
        self.assertEqual(list(missing.customdata[0]), ["D1", "C3"])
        self.assertFalse(any(trace.type == "scatter" and trace.marker.symbol == "square-open" for trace in figure.data))

    def test_unmapped_drugs_do_not_make_known_target_percentage_a_lower_bound(self) -> None:
        row = _summary("C1", "B cell", 2, 40.0, "partial", unmapped=1)
        self.assertEqual(app._display_value(row, "count"), "≥2")
        self.assertEqual(app._display_value(row, "percent"), "40.0%")
        self.assertEqual(app.export_rows([row], "percent")[0]["display_value"], "40.0%")
        row["unknown"] = 1
        self.assertEqual(app._display_value(row, "percent"), "≥40.0%")

    def test_csv_keeps_missing_blank_and_status(self) -> None:
        rows = [_summary("C1", "B cell", None, None, "partial", unknown=5)]

        content = atlas.to_csv(app.export_rows(rows, "count"))
        exported = next(csv.DictReader(io.StringIO(content)))

        self.assertEqual(exported["count"], "")
        self.assertEqual(exported["display_value"], "")
        self.assertEqual(exported["status"], "partial")
        self.assertEqual(exported["denominator"], "5")

    def test_missing_detail_does_not_claim_absence(self) -> None:
        panel = app.detail_panel(
            [_summary("C1", "B cell", None, None, "partial", unknown=5)],
            ("D1", "C1"),
        )

        self.assertIn("insufficient", str(panel))


class SelectionTests(unittest.TestCase):
    """フィルター変更後に以前のセルを参照しないことを検証する。"""

    def test_stale_click_is_cleared(self) -> None:
        visible = [_summary("C2", "T cell", 1, 50.0, "complete", disease_id="D2")]
        stale_click = {"points": [{"customdata": ["D1", "C1"]}]}

        self.assertIsNone(
            app.resolve_selection("heatmap", stale_click, None, None, visible),
        )

    def test_accessible_selectors_choose_visible_cell(self) -> None:
        visible = [_summary("C2", "T cell", 1, 50.0, "complete", disease_id="D2")]

        self.assertEqual(
            app.resolve_selection("detail-cell", None, "D2", "C2", visible),
            ("D2", "C2"),
        )


class AppTests(unittest.TestCase):
    """Dash がブラウザーへ初期ページを返すことを検証する。"""

    def test_index_renders(self) -> None:
        client = app.app.server.test_client()
        for path in ("/", "/_dash-layout", "/_dash-dependencies"):
            with self.subTest(path=path):
                response = client.get(path)
                self.assertEqual(response.status_code, 200)
        self.assertIn(b"Autoimmune", client.get("/").data)

    def test_expression_callbacks_and_csv(self) -> None:
        """実際のコールバック経由で発現表示、詳細、CSV の連携を確認する。"""
        client = app.app.server.test_client()
        layout = client.get("/_dash-layout").get_json()
        components = {}

        def collect(node):
            if isinstance(node, dict):
                props = node.get("props", {})
                if "id" in props:
                    components[props["id"]] = props
                for value in node.values():
                    collect(value)
            elif isinstance(node, list):
                for value in node:
                    collect(value)

        collect(layout)
        self.assertNotIn("mode", components)
        self.assertNotIn("annotation-note", components)
        self.assertEqual(
            [option["value"] for option in components["modality"]["options"]],
            ["all", "small_molecule", "antibody", "protein", "cell", "gene", "enzyme", "oligonucleotide", "unknown"],
        )
        values = {
            ("measure", "value"): "percent",
            ("modality", "value"): "protein",
            ("threshold", "value"): None,
            ("diseases", "value"): ["EFO_1001466"],
            ("cells", "value"): ["CL_0000015"],
            ("detail-disease", "value"): "EFO_1001466",
            ("detail-cell", "value"): "CL_0000015",
            ("heatmap", "clickData"): None,
            ("download-button", "n_clicks"): 1,
        }
        responses = {}
        for key, callback in app.app.callback_map.items():
            outputs = callback["output"]
            outputs = outputs if isinstance(outputs, list) else [outputs]
            output_data = [{"id": item.component_id, "property": item.component_property} for item in outputs]
            for item in callback["inputs"] + callback["state"]:
                self.assertIn(item["id"], components)
            response = client.post("/_dash-update-component", json={
                "output": key,
                "outputs": output_data if isinstance(callback["output"], list) else output_data[0],
                "inputs": [{**item, "value": values[item["id"], item["property"]]} for item in callback["inputs"]],
                "state": [{**item, "value": values[item["id"], item["property"]]} for item in callback["state"]],
                "changedPropIds": ["detail-cell.value"],
            })
            self.assertEqual(response.status_code, 200, response.data)
            responses.update(response.get_json()["response"])
        self.assertEqual(components["threshold"]["value"], 0.5)
        self.assertEqual(app.effective_threshold(0), 0)
        expected = app.visible_rows(app.SNAPSHOT, "protein", None, ["EFO_1001466"], ["CL_0000015"])[0]
        self.assertEqual(responses["heatmap"]["figure"]["data"][0]["z"][0][0], expected["percent"])
        self.assertIn("Median CPM > 0.5", responses["matrix-note"]["children"])
        self.assertIn("Tabula Sapiens", str(responses["details"]))
        exported = list(csv.DictReader(io.StringIO(responses["download"]["data"]["content"].lstrip("\ufeff"))))
        self.assertTrue(exported)
        self.assertEqual({row["mode"] for row in exported}, {"expression"})
        self.assertEqual({row["modality_filter"] for row in exported}, {"protein"})
        self.assertEqual({row["modality"] for row in exported}, {"protein"})
        self.assertEqual({row["expression_threshold"] for row in exported}, {"0.5"})
        self.assertTrue(all(row["data_version"] and row["retrieved_at"] for row in exported))


def _summary(
    cell_id: str,
    cell: str,
    count: int | None,
    percent: float | None,
    status: str,
    *,
    disease_id: str = "D1",
    unknown: int = 0,
    unmapped: int = 0,
) -> dict:
    return {
        "disease_id": disease_id,
        "disease": "Disease " + disease_id,
        "cell_id": cell_id,
        "cell": cell,
        "count": count,
        "percent": percent,
        "denominator": 5,
        "unknown": unknown,
        "unmapped_drugs": unmapped,
        "status": status,
        "records": [],
    }


if __name__ == "__main__":
    unittest.main()
