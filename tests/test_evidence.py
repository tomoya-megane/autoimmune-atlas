"""根拠の行、詳細パネル、発現の表示を検証する。"""

from __future__ import annotations

import math
import unittest
from collections.abc import Mapping, Sequence
from typing import Protocol, cast, override

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas.models import (
    Disease,
    ExpressionStateInput,
    FilteredRecord,
    Snapshot,
    SummaryRow,
    TargetRecord,
)
from autoimmune_atlas.ui import (
    callbacks,
    components,
    figures,
    layout,
)
from tests.ui_fixture import (
    JsonValue,
    ScatterTrace,
    _component,
    _figure,
    _heatmap,
    _scatter,
    fixture,
)


class GraphComponent(Protocol):
    config: Mapping[str, object] | None
    id: str | dict[str, JsonValue] | None


class StoreComponent(Protocol):
    data: JsonValue


def _scatters(figure: object) -> list[ScatterTrace]:
    traces = [trace for trace in _figure(figure).data[1:] if not hasattr(trace, "z")]
    assert all(
        all(hasattr(trace, field) for field in ("name", "x", "y", "mode", "text"))
        for trace in traces
    )
    return [cast(ScatterTrace, trace) for trace in traces]


def _component_list(component: object) -> list[object]:
    children = _component(component).children
    assert isinstance(children, list)
    return cast(list[object], children)


class EvidenceTests(unittest.TestCase):
    """詳細が陽性以外の元記録も保持する。"""

    snapshot: Snapshot = fixture()
    rows: list[SummaryRow] = []

    def test_source_links_only_show_open_targets_pages(self) -> None:
        record: FilteredRecord = {
            **fixture()["records"][0],
            "canonical_stage": "PHASE_3",
        }
        links = str(components.reference_links(record))
        self.assertIn("Open Targets target", links)
        self.assertIn("Open Targets drug", links)
        self.assertIn("Open Targets disease", links)
        self.assertNotIn("PubMed", links)

    @override
    def setUp(self) -> None:
        self.snapshot = fixture()
        self.rows = callbacks.visible_rows(
            self.snapshot, "all", 0.5, ["MONDO_RA_TEST"], ["CL_B_GROUP"]
        )

    def test_source_records_use_only_selected_disease_without_pair_table(self) -> None:
        other: Disease = {
            "id": "MONDO_OTHER",
            "name": "other disease",
            "status": "ready",
        }
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
        panel_children = _component_list(panel)
        self.assertIn("Drug–target records", str(panel_children[0]))
        self.assertIn("Relative target expression", str(panel_children[2]))
        source_children = _component_list(panel_children[1])
        grid = source_children[-1]
        assert hasattr(grid, "id")
        self.assertEqual(cast(GraphComponent, grid).id, "source-records-grid")
        context_store = source_children[0]
        assert hasattr(context_store, "data")
        context = cast(StoreComponent, context_store).data
        assert isinstance(context, dict)
        modality = context["modality"]
        stage = context["stage"]
        disease_id = context["disease_id"]
        assert isinstance(modality, str)
        assert isinstance(stage, str)
        assert isinstance(disease_id, str)
        records = [
            record
            for record in atlas.filtered_records(self.snapshot, modality, stage)
            if record["disease_id"] == disease_id
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
        heatmap = _heatmap(figure)
        targets = list(heatmap.x)
        cells = list(heatmap.y)
        target_index = next(i for i, label in enumerate(targets) if "TARGET1" in label)
        cell_index = cells.index("memory B cell")
        value = heatmap.z[cell_index][target_index]
        assert value is not None
        self.assertGreater(value, 0)
        target2_index = next(i for i, label in enumerate(targets) if "TARGET2" in label)
        cd8_value = heatmap.z[cells.index("CD8-positive T cell")][target2_index]
        naive_b_value = heatmap.z[cells.index("naive B cell")][target2_index]
        assert cd8_value is not None and naive_b_value is not None
        self.assertAlmostEqual(cd8_value, 1)
        self.assertAlmostEqual(naive_b_value, -1)
        self.assertIsNone(heatmap.z[cell_index][target2_index])
        self.assertIn(
            "Target-wise z-score:", heatmap.hovertext[cell_index][target_index]
        )
        self.assertIn("Median CPM: 2", heatmap.hovertext[cell_index][target_index])
        maximum = heatmap.zmax
        assert maximum is not None
        self.assertAlmostEqual(
            maximum,
            max(abs(value) for row in heatmap.z for value in row if value is not None),
        )
        self.assertEqual(heatmap.zmin, -maximum)
        self.assertEqual(heatmap.colorbar.lenmode, "pixels")
        self.assertEqual(heatmap.colorbar.len, 240)
        self.assertEqual(heatmap.colorbar.y, 1)
        self.assertEqual(heatmap.colorbar.yanchor, "top")
        missing = _heatmap(figure, "Missing expression")
        self.assertEqual(list(missing.x), targets)
        self.assertEqual(list(missing.y), cells)
        self.assertEqual(
            {
                (cells[row], targets[column])
                for row, values in enumerate(missing.z)
                for column, value in enumerate(values)
                if value is not None
            },
            {("memory B cell", "TARGET2 (ENSG_TARGET_2)")},
        )
        self.assertIn("Shown in gray", missing.hovertext[cell_index][target2_index])
        dots = _scatter(figure, "CELLEX specificity")
        self.assertEqual(dots.mode, "markers")
        self.assertEqual(dots.marker.sizemode, "diameter")
        self.assertEqual(dots.marker.sizeref, 1)
        self.assertEqual(figures.expression_dot_diameter(0), 3)
        self.assertEqual(figures.expression_dot_diameter(0.0001), 3)
        self.assertEqual(figures.expression_dot_diameter(0.01), 3)
        self.assertGreater(figures.expression_dot_diameter(0.1), 3)
        self.assertLess(
            figures.expression_dot_diameter(0.1),
            figures.expression_dot_diameter(0.8),
        )
        self.assertEqual(
            set(zip(dots.x, dots.y, strict=True)),
            {
                ("TARGET1 (ENSG_TARGET_1)", "memory B cell"),
                ("TARGET1 (ENSG_TARGET_1)", "naive B cell"),
                ("TARGET1 (ENSG_TARGET_1)", "CD8-positive T cell"),
                ("TARGET2 (ENSG_TARGET_2)", "naive B cell"),
                ("TARGET2 (ENSG_TARGET_2)", "CD8-positive T cell"),
            },
        )
        dot_sizes = dict(
            zip(
                zip(dots.x, dots.y, strict=True),
                cast(Sequence[float], dots.marker.size),
                strict=True,
            )
        )
        self.assertEqual(
            dot_sizes,
            {
                (
                    "TARGET1 (ENSG_TARGET_1)",
                    "memory B cell",
                ): figures.expression_dot_diameter(0.8),
                (
                    "TARGET1 (ENSG_TARGET_1)",
                    "naive B cell",
                ): figures.expression_dot_diameter(0.2),
                (
                    "TARGET1 (ENSG_TARGET_1)",
                    "CD8-positive T cell",
                ): figures.expression_dot_diameter(0.1),
                (
                    "TARGET2 (ENSG_TARGET_2)",
                    "naive B cell",
                ): figures.expression_dot_diameter(0.75),
                (
                    "TARGET2 (ENSG_TARGET_2)",
                    "CD8-positive T cell",
                ): figures.expression_dot_diameter(0.7),
            },
        )
        self.assertEqual(len(dots.customdata or ()), len(dots.x))
        assert dots.hovertext is not None
        self.assertIn("Expression rule:", dots.hovertext[0])
        self.assertEqual(_figure(figure).layout.xaxis.side, "top")
        axis_range = _figure(figure).layout.yaxis.range
        assert axis_range is not None
        self.assertEqual(list(axis_range), [len(cells) - 0.5, -0.5])

    def test_group_expression_averages_cpm_before_transform_and_keeps_scale(self):
        figure = figures.expression_figure(
            self.snapshot, self.snapshot["records"], grouped=True
        )
        heatmap = _heatmap(figure)
        row = list(heatmap.y).index("B cell (group)")
        col = list(heatmap.x).index("TARGET1 (ENSG_TARGET_1)")
        values = [math.log2(1 + x) for x in (2, 0.1, 0.2)]
        center = sum(values) / len(values)
        spread = math.sqrt(sum((x - center) ** 2 for x in values) / len(values))
        grouped_value = heatmap.z[row][col]
        assert grouped_value is not None
        self.assertAlmostEqual(grouped_value, (math.log2(1 + 1.05) - center) / spread)
        self.assertIn("Mean CPM: 1.05", heatmap.hovertext[row][col])
        self.assertIn("Observed source cell types: 2 / 2", heatmap.hovertext[row][col])
        dots = _scatter(figure, "CELLEX specificity")
        assert dots.hovertext is not None
        dot_index = list(zip(dots.x, dots.y, strict=True)).index(
            ("TARGET1 (ENSG_TARGET_1)", "B cell (group)")
        )
        self.assertEqual(
            cast(Sequence[float], dots.marker.size)[dot_index],
            figures.expression_dot_diameter(0.5),
        )
        self.assertIn("Mean CELLEX specificity: 0.5", dots.hovertext[dot_index])
        self.assertIn(
            "CELLEX-observed source cell types: 2 / 2", dots.hovertext[dot_index]
        )
        target2_dot = list(zip(dots.x, dots.y, strict=True)).index(
            ("TARGET2 (ENSG_TARGET_2)", "B cell (group)")
        )
        self.assertEqual(
            cast(Sequence[float], dots.marker.size)[target2_dot],
            figures.expression_dot_diameter(0.75),
        )
        self.assertIn(
            "CELLEX-observed source cell types: 1 / 2",
            dots.hovertext[target2_dot],
        )
        col2 = list(heatmap.x).index("TARGET2 (ENSG_TARGET_2)")
        self.assertIn("Mean CPM: 0.6", heatmap.hovertext[row][col2])
        self.assertIn("Observed source cell types: 1 / 2", heatmap.hovertext[row][col2])
        original = figures.expression_figure(self.snapshot, self.snapshot["records"])
        original_heatmap = _heatmap(original)
        self.assertEqual(original_heatmap.x, heatmap.x)
        for i, name in enumerate(original_heatmap.y):
            self.assertEqual(
                original_heatmap.z[i], heatmap.z[list(heatmap.y).index(name)]
            )
        metadata = atlas.expression_metadata(self.snapshot)
        figures.expression_figure(
            self.snapshot, self.snapshot["records"], metadata, grouped=True
        )
        self.assertEqual(metadata, atlas.expression_metadata(self.snapshot))
        catalog = atlas.cell_catalog(self.snapshot, "mixed")
        collapsed = figures.expression_view(figure, catalog, [])
        expanded = figures.expression_view(figure, catalog, ["group:CL_B_GROUP"])
        collapsed_heatmap = _heatmap(collapsed)
        expanded_heatmap = _heatmap(expanded)
        self.assertEqual(len(collapsed_heatmap.y), 2)
        self.assertEqual(len(expanded_heatmap.y), 4)
        self.assertEqual(collapsed_heatmap.zmax, expanded_heatmap.zmax)
        self.assertEqual(collapsed_heatmap.x, expanded_heatmap.x)
        self.assertEqual(len(heatmap.y), 5)
        for trace in _scatters(expanded):
            self.assertTrue(set(trace.y) <= set(expanded_heatmap.y))
            if trace.name == "Meets expression rule":
                self.assertNotIn("B cell (group)", trace.y)

    def test_expression_view_fixes_x_range_to_heatmap_cell_edges(self):
        figure = figures.expression_figure(
            self.snapshot, self.snapshot["records"], grouped=True
        )
        catalog = atlas.cell_catalog(self.snapshot, "mixed")

        for expanded in ([], ["group:CL_B_GROUP"]):
            view = figures.expression_view(figure, catalog, expanded)
            xaxis = _figure(view).layout.xaxis
            axis_range = xaxis.range
            assert axis_range is not None
            self.assertEqual(list(axis_range), [-0.5, 1.5])
            self.assertFalse(xaxis.autorange)

            dot = figures.expression_view(figure, catalog, expanded, "dot")
            self.assertEqual(len(_figure(dot).data), 1)
            dots = _scatter(dot, "CELLEX specificity")
            x_categories = _figure(dot).layout.xaxis.categoryarray
            y_categories = _figure(dot).layout.yaxis.categoryarray
            view_heatmap = _heatmap(view)
            assert x_categories is not None and y_categories is not None
            assert view_heatmap.x is not None and view_heatmap.y is not None
            self.assertEqual(list(x_categories), list(view_heatmap.x))
            self.assertEqual(list(y_categories), list(view_heatmap.y))
            self.assertTrue(set(dots.y) <= set(view_heatmap.y))
            self.assertEqual(
                {
                    len(dots.x),
                    len(dots.y),
                    len(dots.hovertext or ()),
                    len(dots.customdata or ()),
                    len(cast(Sequence[float], dots.marker.size)),
                    len(cast(Sequence[float], dots.marker.color)),
                },
                {len(dots.x)},
            )

        hidden_base = figures.expression_figure(
            self.snapshot, self.snapshot["records"], grouped=True
        )
        hidden_dots = _scatter(hidden_base, "CELLEX specificity")
        hidden_indices = [
            index for index, name in enumerate(hidden_dots.y) if "(group)" not in name
        ]
        hidden_dots.x = [hidden_dots.x[index] for index in hidden_indices]
        hidden_dots.y = [hidden_dots.y[index] for index in hidden_indices]
        hidden_dots.hovertext = [
            (hidden_dots.hovertext or ())[index] for index in hidden_indices
        ]
        hidden_dots.customdata = [
            (hidden_dots.customdata or ())[index] for index in hidden_indices
        ]
        hidden_dots.marker.size = [
            cast(Sequence[float], hidden_dots.marker.size)[index]
            for index in hidden_indices
        ]
        hidden_dots.marker.color = [
            cast(Sequence[float], hidden_dots.marker.color)[index]
            for index in hidden_indices
        ]
        hidden_view = figures.expression_view(hidden_base, catalog, [], "dot")
        self.assertEqual(
            list(_scatter(hidden_view, "CELLEX specificity").x), [None, None]
        )

    def test_group_expression_is_missing_when_all_members_are_missing(self):
        for item in self.snapshot["expression"]["ENSG_TARGET_2"]:
            item["median"] = None
        figure = figures.expression_figure(
            self.snapshot, self.snapshot["records"], grouped=True
        )
        heatmap = _heatmap(figure)
        column = list(heatmap.x).index("TARGET2 (ENSG_TARGET_2)")
        self.assertTrue(all(row[column] is None for row in heatmap.z))
        row = list(heatmap.y).index("B cell (group)")
        self.assertIn(
            "Observed source cell types: 0 / 2", heatmap.hovertext[row][column]
        )

    def test_expression_dot_shows_zero_and_blanks_unavailable_values(self) -> None:
        zero_group_snapshot = fixture()
        zero_group_snapshot["expression"]["ENSG_TARGET_1"][0]["specificity_score"] = 0
        zero_group_snapshot["expression"]["ENSG_TARGET_1"][1]["specificity_score"] = 0
        zero_group_figure = figures.expression_figure(
            zero_group_snapshot, zero_group_snapshot["records"], grouped=True
        )
        zero_group_dots = _scatter(zero_group_figure, "CELLEX specificity")
        zero_group_sizes = dict(
            zip(
                zip(zero_group_dots.x, zero_group_dots.y, strict=True),
                cast(Sequence[float], zero_group_dots.marker.size),
                strict=True,
            )
        )
        self.assertEqual(
            zero_group_sizes[("TARGET1 (ENSG_TARGET_1)", "B cell (group)")], 3
        )

        snapshot = fixture()
        snapshot["expression"]["ENSG_TARGET_1"][0]["specificity_score"] = None
        snapshot["expression"]["ENSG_TARGET_1"][1]["specificity_score"] = None
        snapshot["expression"]["ENSG_TARGET_1"][2]["specificity_score"] = 0
        snapshot["expression"]["ENSG_TARGET_2"][0]["specificity_score"] = 0.9
        figure = figures.expression_figure(snapshot, snapshot["records"], grouped=True)
        dots = _scatter(figure, "CELLEX specificity")
        coordinates = set(zip(dots.x, dots.y, strict=True))
        self.assertNotIn(("TARGET1 (ENSG_TARGET_1)", "B cell (group)"), coordinates)
        self.assertIn(("TARGET1 (ENSG_TARGET_1)", "CD8-positive T cell"), coordinates)
        zero_index = list(zip(dots.x, dots.y, strict=True)).index(
            ("TARGET1 (ENSG_TARGET_1)", "CD8-positive T cell")
        )
        self.assertEqual(cast(Sequence[float], dots.marker.size)[zero_index], 3)
        self.assertNotIn(("TARGET2 (ENSG_TARGET_2)", "memory B cell"), coordinates)
        heatmap = _heatmap(figure)
        row = list(heatmap.y).index("B cell (group)")
        column = list(heatmap.x).index("TARGET1 (ENSG_TARGET_1)")
        self.assertIn(
            "Mean CELLEX specificity: missing", heatmap.hovertext[row][column]
        )
        self.assertIn(
            "CELLEX-observed source cell types: 0 / 2",
            heatmap.hovertext[row][column],
        )

        for items in snapshot["expression"].values():
            for item in items:
                item["specificity_score"] = None
        empty_base = figures.expression_figure(
            snapshot, snapshot["records"], grouped=True
        )
        empty = figures.expression_view(
            empty_base, atlas.cell_catalog(snapshot, "mixed"), [], "dot"
        )
        carrier = _scatter(empty, "CELLEX specificity")
        self.assertEqual(list(carrier.x), [None, None])
        self.assertTrue(carrier.marker.showscale)
        empty_categories = _figure(empty).layout.xaxis.categoryarray
        empty_heatmap_x = _heatmap(empty_base).x
        assert empty_categories is not None and empty_heatmap_x is not None
        self.assertEqual(
            list(empty_categories),
            list(empty_heatmap_x),
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
        heatmap = _heatmap(figure)
        target_index = next(
            i for i, label in enumerate(heatmap.x) if "TARGET1" in label
        )
        self.assertTrue(all(row[target_index] == 0 for row in heatmap.z))
        constant_only = figures.expression_figure(
            self.snapshot,
            [record for record in evidence if record["target_id"] == "ENSG_TARGET_1"],
        )
        self.assertEqual(
            (_heatmap(constant_only).zmin, _heatmap(constant_only).zmax), (-1, 1)
        )
        constant_dots = _scatter(constant_only, "CELLEX specificity")
        self.assertEqual(
            list(cast(Sequence[float], constant_dots.marker.color)), [0, 0, 0]
        )

    def test_expression_markers_follow_current_rule(self) -> None:
        records = self.snapshot["records"]
        fixed = figures.expression_figure(
            self.snapshot, records, threshold=0.5, method="fixed"
        )
        specific = figures.expression_figure(
            self.snapshot, records, threshold=0.5, method="specificity"
        )

        def marked_cells(figure: object) -> set[tuple[str, str]]:
            marks = _scatter(figure, "Meets expression rule")
            return set(zip(marks.x, marks.y, strict=True))

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
        panel_children = _component_list(panel)
        self.assertIn("expression-chart-type", str(panel_children[2]))
        self.assertIn("○ Cell type meets expression rule", str(panel_children[4]))
        self.assertIn("expression-dot-key", str(panel_children[5]))
        graph_container = _component(_component(panel_children[3]).children)
        graph_children = graph_container.children
        assert isinstance(graph_children, list) and graph_children
        graph = cast(list[object], graph_children)[0]
        assert hasattr(graph, "id") and hasattr(graph, "config")
        typed_graph = cast(GraphComponent, graph)
        self.assertEqual(typed_graph.id, "expression-heatmap")
        assert typed_graph.config is not None
        self.assertTrue(typed_graph.config["responsive"])
        stores = _component_list(panel_children[1])
        assert stores and hasattr(stores[0], "data")
        store_data = cast(StoreComponent, stores[0]).data
        assert isinstance(store_data, dict)
        self.assertEqual(store_data["method"], "specificity")

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
        metadata: dict[tuple[str, str], ExpressionStateInput] = {
            (target, cell): {"median": value}
            for target, values in profiles.items()
            for cell, value in zip(cells, values, strict=True)
        }
        records: list[TargetRecord] = [
            {"target_id": target, "target": target}
            for target in ["A", "B", "C", "D", "E"]
        ]
        labels = list(
            _heatmap(figures.expression_figure(snapshot, records, metadata)).x
        )
        self.assertEqual(labels[-1], "E (E)")
        self.assertEqual(abs(labels.index("A (A)") - labels.index("C (C)")), 1)
        self.assertEqual(abs(labels.index("B (B)") - labels.index("D (D)")), 1)

    def test_continuous_expression_shows_all_cells(
        self,
    ) -> None:
        row: SummaryRow = {
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
            list(_heatmap(figure).y),
            ["CD8-positive T cell", "memory B cell", "naive B cell"],
        )
