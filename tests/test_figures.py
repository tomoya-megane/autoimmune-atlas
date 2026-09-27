"""ヒートマップとドットプロットの図の値、軸、色、大きさを検証する。"""

from __future__ import annotations

import math
import unittest
from collections.abc import Sequence
from typing import Literal, cast

from autoimmune_atlas.ui import (
    figures,
)
from tests.ui_fixture import (
    ScatterTrace,
    _figure,
    _heatmap,
    _scatter,
    summary,
)


class FigureTests(unittest.TestCase):
    """2種類の指標で 0、下限、欠測を独立に扱う。"""

    def test_target_and_drug_lower_bounds_are_independent(self) -> None:
        row = summary(
            count=0, percent=0.0, drug_count=1, drug_percent=25.0, unknown_drugs=1
        )
        self.assertEqual(figures.display_value(row, "count"), "0")
        self.assertEqual(figures.display_value(row, "count", "drug"), "≥1")
        self.assertEqual(figures.display_value(row, "percent"), "0.0%")
        self.assertEqual(figures.display_value(row, "percent", "drug"), "≥25.0%")

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
        self.assertEqual(figures.display_value(row, "count"), "≥2")
        self.assertEqual(figures.display_value(row, "percent"), "40.0%")
        self.assertEqual(figures.display_value(row, "count", "drug"), "≥1")
        self.assertEqual(figures.display_value(row, "percent", "drug"), "25.0%")

    def test_long_disease_ticks_wrap_without_changing_data(self) -> None:
        name = "anti-neutrophil cytoplasmic antibody-associated vasculitis"
        rows = [summary(disease=name), summary(disease_id="D2", disease="Short name")]
        for kind in ("target", "drug"):
            figure = figures.build_figure(rows, ["D1", "D2"], ["C1"], "count", kind)
            heatmap = _heatmap(figure)
            ticktext = _figure(figure).layout.xaxis.ticktext
            assert ticktext is not None
            self.assertEqual(
                list(ticktext),
                [
                    "anti-neutrophil cytoplasmic<br>antibody-associated vasculitis",
                    "Short name",
                ],
            )
            self.assertEqual(list(heatmap.x), [name, "Short name"])
            self.assertEqual(heatmap.customdata[0][0], ["D1", "C1"])
            self.assertEqual(_figure(figure).layout.xaxis.tickangle, -45)

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
        cases: tuple[
            tuple[
                Literal["target", "drug"],
                Literal["count", "percent"],
                str,
                int | float,
            ],
            ...,
        ] = (
            ("target", "count", "2", 2),
            ("target", "percent", "40.0%", 40.0),
            ("drug", "count", "1", 1),
            ("drug", "percent", "25.0%", 25.0),
        )
        for kind, measure, label, value in cases:
            figure = figures.build_figure([row], ["D1"], ["C1"], measure, kind)
            heatmap = _heatmap(figure)
            self.assertEqual(heatmap.text[0][0], label)
            self.assertEqual(heatmap.z[0][0], value)
            self.assertIn("≥" + label, heatmap.hovertext[0][0])
            self.assertIn("≥ is a lower bound", heatmap.hovertext[0][0])

    def test_color_ranges_follow_each_visible_measure(self) -> None:
        rows = [
            summary(),
            summary(cell_id="C2", cell="T cell", percent=60, drug_percent=40),
        ]
        target = figures.build_figure(rows, ["D1"], ["C1", "C2"], "percent", "target")
        drug = figures.build_figure(rows, ["D1"], ["C1", "C2"], "percent", "drug")
        self.assertEqual((_heatmap(target).zmin, _heatmap(target).zmax), (20, 60))
        self.assertEqual((_heatmap(drug).zmin, _heatmap(drug).zmax), (25, 40))
        flat = figures.build_figure([summary(percent=0)], ["D1"], ["C1"], "percent")
        self.assertEqual((_heatmap(flat).zmin, _heatmap(flat).zmax), (0, 1))
        missing = figures.build_figure(
            [summary(percent=None)], ["D1"], ["C1"], "percent"
        )
        self.assertEqual((_heatmap(missing).zmin, _heatmap(missing).zmax), (0, 1))

    def test_missing_percent_uses_same_fill_and_text_as_zero(self):
        rows = [
            summary(percent=0, drug_percent=0),
            summary(cell_id="C2", cell="T cell", percent=None, drug_percent=None),
        ]
        for kind in ("target", "drug"):
            figure = figures.build_figure(rows, ["D1"], ["C1", "C2"], "percent", kind)
            heatmap = _heatmap(figure)
            self.assertEqual(heatmap.z[0][0], heatmap.z[1][0])
            self.assertEqual(heatmap.text[0][0], heatmap.text[1][0])
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
            figure_data = _figure(figure)
            heatmap = _heatmap(figure)
            self.assertIsNone(figure_data.layout.width)
            self.assertEqual(heatmap.colorbar.lenmode, "pixels")
            self.assertEqual(heatmap.colorbar.len, 240)
            self.assertEqual(heatmap.colorbar.y, 1)
            self.assertEqual(heatmap.colorbar.yanchor, "top")
        for figure in (target, drug):
            self.assertFalse(_figure(figure).layout.yaxis.autorange)
            axis_range = _figure(figure).layout.yaxis.range
            assert axis_range is not None
            self.assertEqual(list(axis_range), [1.5, -0.5])
        target_heatmap = _heatmap(target)
        self.assertEqual(list(target_heatmap.x), list(_heatmap(drug).x))
        self.assertEqual(list(target_heatmap.y), list(_heatmap(drug).y))
        self.assertEqual(len(_figure(target).data), 1)
        self.assertEqual(target_heatmap.z[1][0], 0)
        self.assertEqual(target_heatmap.text[1][0], "0")
        self.assertIsNone(rows[1]["count"])
        hover = target_heatmap.hovertext[1][0]
        self.assertIn("Denominator: 5", hover)
        self.assertIn("Unresolved in denominator: 1", hover)
        self.assertIn(
            "Canonical drugs with targets / all canonical drugs: 4 / 4", hover
        )
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
            heatmap = _heatmap(figure)
            self.assertEqual(heatmap.z[0][0], 0)
            self.assertEqual(heatmap.text[0][0], "0")
            self.assertEqual(len(_figure(figure).data), 1)
            self.assertIn("≥0", heatmap.hovertext[0][0])
        percent = figures.build_figure([row], ["D1"], ["C1"], "percent")
        self.assertEqual(_heatmap(percent).z[0][0], 0)
        self.assertEqual(_heatmap(percent).text[0][0], "0.0%")
        self.assertIsNone(row["percent"])

    def test_no_value_hover_distinguishes_unavailable_from_empty_denominator(
        self,
    ) -> None:
        unavailable = summary(
            status="unavailable", count=None, percent=None, denominator=0
        )
        empty = summary(percent=None, denominator=0)
        self.assertIn(
            "Disease data not loaded",
            figures.hover_text(unavailable, "percent", "target"),
        )
        self.assertNotIn(
            "No eligible items", figures.hover_text(unavailable, "percent", "target")
        )
        self.assertIn(
            "No eligible items in the percentage denominator",
            figures.hover_text(empty, "percent", "target"),
        )

    def test_dot_area_and_color_map_count_and_percent_for_each_kind(self) -> None:
        rows = [
            summary(count=1, percent=20, drug_count=2, drug_percent=25),
            summary(
                cell_id="C2",
                cell="T cell",
                count=4,
                percent=80,
                drug_count=8,
                drug_percent=75,
            ),
        ]
        for kind, counts, colors, maximum in (
            ("target", [1, 4], [20, 80], 4),
            ("drug", [2, 8], [25, 75], 8),
        ):
            with self.subTest(kind=kind):
                figure = figures.build_dot_figure(
                    rows,
                    ["D1"],
                    ["C1", "C2"],
                    cast(figures.Kind, kind),
                    scale_rows=rows,
                )
                trace = _scatter(figure, "Values")
                marker = trace.marker
                marker_sizes = marker.size
                marker_colors = marker.color
                customdata = trace.customdata
                hovertext = trace.hovertext
                assert isinstance(marker_sizes, Sequence)
                assert isinstance(marker_colors, Sequence) and not isinstance(
                    marker_colors, str
                )
                assert customdata is not None and hovertext is not None
                expected_sizes = [count**0.75 for count in counts]
                self.assertEqual(list(marker_sizes), expected_sizes)
                self.assertEqual(list(marker_colors), colors)
                self.assertEqual(marker.sizemode, "area")
                self.assertEqual(marker.sizeref, 2 * maximum**0.75 / 22**2)
                self.assertEqual((marker.cmin, marker.cmax), (min(colors), max(colors)))
                self.assertTrue(marker.showscale)
                self.assertEqual(customdata[0], ["D1", "C1"])
                self.assertIn("Count:", hovertext[0])
                self.assertIn("Percent:", hovertext[0])
                _, _, legend = figures.dot_size_scale(rows, cast(figures.Kind, kind))
                self.assertEqual(legend[-1], (maximum, 22))
                self.assertEqual(legend[0][1], 22 * (counts[0] / maximum) ** 0.375)
                self.assertGreater(legend[0][1], 22 * math.sqrt(counts[0] / maximum))
                self.assertEqual(
                    [diameter for _, diameter in legend],
                    sorted(diameter for _, diameter in legend),
                )
                assert marker.sizeref is not None
                for count, (_, legend_diameter) in zip(counts, legend, strict=True):
                    self.assertAlmostEqual(
                        math.sqrt(2 * math.pow(count, 0.75) / marker.sizeref),
                        legend_diameter,
                    )
                self.assertNotIn(
                    str(maximum),
                    [
                        cast(ScatterTrace, trace).name
                        for trace in _figure(figure).data
                        if hasattr(trace, "name")
                    ],
                )

    def test_dot_leaves_zero_and_missing_blank_without_changing_categories(
        self,
    ) -> None:
        rows = [
            summary(count=0, percent=0),
            summary(
                cell_id="C2",
                cell="Missing first",
                count=None,
                percent=None,
                unknown=1,
                status="unavailable",
            ),
            summary(cell_id="C3", cell="Positive last", count=3, percent=60),
        ]
        figure = figures.build_dot_figure(
            rows, ["D1"], ["C1", "C2", "C3"], "target", scale_rows=rows
        )
        positive = _scatter(figure, "Values")
        with self.assertRaisesRegex(AssertionError, "trace not found: Zero"):
            _scatter(figure, "Zero")
        with self.assertRaisesRegex(AssertionError, "trace not found: Not available"):
            _scatter(figure, "Not available")
        plotted_cells = [
            y
            for trace in _figure(figure).data
            if hasattr(trace, "y")
            for y in (cast(ScatterTrace, trace).y or ())
        ]
        self.assertNotIn("B cell", plotted_cells)
        self.assertNotIn("Missing first", plotted_cells)
        self.assertEqual(list(positive.y), ["Positive last"])
        self.assertEqual(_figure(figure).layout.xaxis.type, "category")
        self.assertEqual(_figure(figure).layout.xaxis.categoryorder, "array")
        x_categories = _figure(figure).layout.xaxis.categoryarray
        x_range = _figure(figure).layout.xaxis.range
        assert x_categories is not None and x_range is not None
        self.assertEqual(list(x_categories), ["Disease one"])
        self.assertEqual(list(x_range), [-0.5, 0.5])
        self.assertFalse(_figure(figure).layout.xaxis.autorange)
        self.assertEqual(_figure(figure).layout.yaxis.categoryorder, "array")
        y_categories = _figure(figure).layout.yaxis.categoryarray
        y_range = _figure(figure).layout.yaxis.range
        assert y_categories is not None and y_range is not None
        self.assertEqual(
            list(y_categories), ["B cell", "Missing first", "Positive last"]
        )
        self.assertEqual(list(y_range), [2.5, -0.5])

    def test_dot_keeps_percent_colorbar_when_no_positive_dot_exists(self) -> None:
        for row in (
            summary(count=0, percent=0),
            summary(count=None, percent=None, status="unavailable"),
        ):
            with self.subTest(status=row["status"], percent=row["percent"]):
                figure = figures.build_dot_figure([row], ["D1"], ["C1"])
                values = _scatter(figure, "Values")
                with self.assertRaisesRegex(AssertionError, "trace not found: Zero"):
                    _scatter(figure, "Zero")
                with self.assertRaisesRegex(
                    AssertionError, "trace not found: Not available"
                ):
                    _scatter(figure, "Not available")
                self.assertEqual(list(values.x), [None, None])
                self.assertEqual(
                    list(cast(Sequence[float], values.marker.size)), [0, 0]
                )
                self.assertTrue(values.marker.showscale)
                self.assertEqual((values.marker.cmin, values.marker.cmax), (0, 1))

    def test_dot_scale_includes_collapsed_source_cells(self) -> None:
        group = summary(count=1, percent=20)
        child = summary(cell_id="C2", cell="Source cell", count=9, percent=90)
        collapsed = figures.build_dot_figure(
            [group], ["D1"], ["C1"], "target", scale_rows=[group, child]
        )
        expanded = figures.build_dot_figure(
            [group, child], ["D1"], ["C1", "C2"], "target", scale_rows=[group, child]
        )
        for attribute in ("sizeref", "cmin", "cmax"):
            self.assertEqual(
                getattr(_scatter(collapsed, "Values").marker, attribute),
                getattr(_scatter(expanded, "Values").marker, attribute),
            )
