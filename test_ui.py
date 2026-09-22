"""Dash UI の表示用ロジックを検証する。"""

from __future__ import annotations

import csv
import io
import tempfile
import unittest
from pathlib import Path

import app
import atlas


class AnnotationTests(unittest.TestCase):
    """注釈 CSV の入力境界を検証する。"""

    def test_missing_annotation_file_is_empty(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(app.load_annotations(Path(directory) / "missing.csv"), [])

    def test_invalid_source_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "annotations.csv"
            path.write_text(
                "disease_id,drug_id,target_id,cell_id,cell,status,source,note\n"
                "D1,R1,T1,C1,B cell,yes,ftp://example.org,evidence\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "HTTP"):
                app.load_annotations(path)

    def test_conflicting_duplicate_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "annotations.csv"
            path.write_text(
                "disease_id,drug_id,target_id,cell_id,cell,status,source,note\n"
                "D1,R1,T1,C1,B cell,yes,https://example.org,a\n"
                "D1,R1,T1,C1,B cell,no,https://example.org,b\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "duplicate"):
                app.load_annotations(path)


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
