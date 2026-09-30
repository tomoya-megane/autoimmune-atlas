"""Dash の callback とレイアウトの初期値を検証する。"""

from __future__ import annotations

import json
import unittest
from collections.abc import Iterator
from typing import ClassVar, Protocol, cast, override
from unittest.mock import patch

from flask.testing import FlaskClient
from werkzeug.test import TestResponse

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas.models import (
    FilteredRecord,
    Snapshot,
)
from autoimmune_atlas.ui import (
    application,
    components,
    config,
    drugs_callbacks,
    figures,
)
from tests.ui_fixture import (
    ASSETS_PATH,
    DashApplication,
    JsonValue,
    fixture,
)

type JsonObject = dict[str, JsonValue]


type CallbackValues = dict[tuple[str, str], JsonValue]


class GridComponent(Protocol):
    id: str
    rowData: list[dict[str, str]]
    columnDefs: list[dict[str, object]]
    defaultColDef: dict[str, object]
    dashGridOptions: dict[str, object]


def _json_value(value: object) -> JsonValue:
    if value is None or isinstance(value, bool | float | int | str):
        return value
    if isinstance(value, list):
        return [_json_value(item) for item in cast(list[object], value)]
    if isinstance(value, dict):
        result: JsonObject = {}
        for key, item in cast(dict[object, object], value).items():
            assert isinstance(key, str)
            result[key] = _json_value(item)
        return result
    raise AssertionError(f"not JSON-compatible: {type(value).__name__}")


def _json_object(value: object) -> JsonObject:
    parsed = _json_value(value)
    assert isinstance(parsed, dict)
    return parsed


def _json_array(value: JsonValue) -> list[JsonValue]:
    assert isinstance(value, list)
    return value


def _json_string(value: JsonValue) -> str:
    assert isinstance(value, str)
    return value


def _at(value: JsonValue, *path: str | int) -> JsonValue:
    current = value
    for key in path:
        if isinstance(key, str):
            current = _json_object(current)[key]
        else:
            current = _json_array(current)[key]
    return current


def _response_json(response: TestResponse) -> JsonObject:
    raw = cast(object, json.loads(response.get_data(as_text=True)))
    return _json_object(raw)


class CallbackTests(unittest.TestCase):
    """HTTP 経由で 階層選択、2図、詳細の callback を確認する。"""

    snapshot: ClassVar[Snapshot]
    application: ClassVar[DashApplication | None] = None
    client: ClassVar[FlaskClient | None] = None
    components: ClassVar[dict[str, JsonObject]]

    @classmethod
    @override
    def setUpClass(cls) -> None:
        cls.snapshot = fixture()
        app = application.create_app(cls.snapshot, assets_folder=ASSETS_PATH)
        assert hasattr(app, "callback_map") and hasattr(app, "server")
        typed_app = cast(DashApplication, cast(object, app))
        cls.application = typed_app
        client = typed_app.server.test_client()
        cls.client = client
        cls.components = {}

        def collect(node: JsonValue) -> None:
            if isinstance(node, dict):
                props_value = node.get("props", {})
                assert isinstance(props_value, dict)
                props = props_value
                component_id = props.get("id")
                if isinstance(component_id, str):
                    cls.components[component_id] = props
                for value in node.values():
                    collect(value)
            elif isinstance(node, list):
                for value in node:
                    collect(value)

        collect(_response_json(client.get("/_dash-layout")))

    def _application(self) -> DashApplication:
        assert self.application is not None
        return self.application

    def _client(self) -> FlaskClient:
        assert self.client is not None
        return self.client

    def _callback_key(self, output_id: str) -> str:
        # 遺伝子ページの genetics- の ID と取り違えないよう、出力の ID を完全一致で比べる。
        return next(
            key
            for key in self._application().callback_map
            if output_id in key.strip(".").split("...")
        )

    def _post(self, output_id: str, values: CallbackValues, changed: str) -> JsonObject:
        key = self._callback_key(output_id)
        callback = self._application().callback_map[key]
        self.assertIn(
            changed, [f"{item['id']}.{item['property']}" for item in callback["inputs"]]
        )
        outputs = (
            callback["output"]
            if isinstance(callback["output"], list)
            else [callback["output"]]
        )
        response = self._client().post(
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
        self.assertEqual(response.status_code, 200, response.get_data())
        body = _response_json(response)["response"]
        assert isinstance(body, dict)
        return body

    def _apply(self, values: CallbackValues) -> None:
        result = self._post("applied-parameters.data", values, "update-button.n_clicks")
        applied = result["applied-parameters"]
        assert isinstance(applied, dict)
        values[("applied-parameters", "data")] = applied["data"]

    def _post_applied(self, output_id: str, values: CallbackValues) -> JsonObject:
        self._apply(values)
        return self._post(output_id, values, "applied-parameters.data")

    def _select_details(
        self, values: CallbackValues, disease: str = "MONDO_RA_TEST"
    ) -> JsonObject:
        values = dict(values)
        values[("detail-disease", "value")] = disease
        return self._post("details.children", values, "detail-disease.value")

    def _values(self) -> CallbackValues:
        values: CallbackValues = {
            ("measure", "value"): "percent",
            ("modality", "value"): "all",
            ("target-class", "value"): "all",
            ("target-location", "value"): "all",
            ("stage", "value"): "phase3",
            ("method", "value"): "fixed",
            ("threshold", "value"): None,
            ("specificity", "value"): None,
            ("diseases", "value"): ["MONDO_RA_TEST"],
            ("detail-disease", "value"): "MONDO_RA_TEST",
            ("update-button", "n_clicks"): 1,
            ("heatmap-view", "value"): "target",
            ("chart-type", "value"): "heatmap",
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
                _at(expanded, graph, "figure", "layout", "xaxis", "tickangle"),
                -45,
            )
            old = _at(before, graph, "figure", "data", 0)
            new = _at(expanded, graph, "figure", "data", 0)
            self.assertGreater(
                len(_json_array(_at(new, "y"))),
                len(_json_array(_at(old, "y"))),
            )
            self.assertEqual(
                (_at(old, "zmin"), _at(old, "zmax")),
                (_at(new, "zmin"), _at(new, "zmax")),
            )
            self.assertEqual(_at(old, "z", 0), _at(new, "z", 0))
        values[("expanded-cell-groups", "data")] = []
        collapsed = self._post(
            "target-heatmap.figure", values, "expanded-cell-groups.data"
        )
        self.assertEqual(
            _at(before, "target-heatmap", "figure"),
            _at(collapsed, "target-heatmap", "figure"),
        )
        self.assertEqual(
            figures.heatmap_cell_ids(
                self.snapshot, ["group:CL_B_GROUP", "CL_B_ONE"], []
            ),
            ["group:CL_B_GROUP", "CL_B_ONE"],
        )

    def test_dot_expansion_preserves_area_and_color_scales(self) -> None:
        values = self._values()
        values[("chart-type", "value")] = "dot"
        before = self._post("target-heatmap.figure", values, "chart-type.value")
        values[("expanded-cell-groups", "data")] = ["group:CL_B_GROUP"]
        expanded = self._post(
            "target-heatmap.figure", values, "expanded-cell-groups.data"
        )
        for graph in ("target-heatmap", "drug-heatmap"):
            old_marker = _at(before, graph, "figure", "data", 0, "marker")
            new_marker = _at(expanded, graph, "figure", "data", 0, "marker")
            self.assertEqual(
                (
                    _at(old_marker, "sizeref"),
                    _at(old_marker, "cmin"),
                    _at(old_marker, "cmax"),
                ),
                (
                    _at(new_marker, "sizeref"),
                    _at(new_marker, "cmin"),
                    _at(new_marker, "cmax"),
                ),
            )
            self.assertGreater(
                len(
                    _json_array(
                        _at(
                            expanded,
                            graph,
                            "figure",
                            "layout",
                            "yaxis",
                            "categoryarray",
                        )
                    )
                ),
                len(
                    _json_array(
                        _at(before, graph, "figure", "layout", "yaxis", "categoryarray")
                    )
                ),
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
            response = self._client().post(
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
            response_json = _response_json(response)
            if expected is None:
                self.assertNotIn(
                    "expanded-cell-groups",
                    _json_object(response_json.get("response", {})),
                )
            else:
                self.assertEqual(response.status_code, 200)
                self.assertEqual(
                    _at(
                        response_json,
                        "response",
                        "expanded-cell-groups",
                        "data",
                    ),
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
                **_json_object(values[("applied-parameters", "data")]),
                "threshold": 123,
            }
            result = self._post(
                "target-heatmap.figure", values, "applied-parameters.data"
            )
            summarize.assert_called_once()
            self.assertEqual(summarize.call_args.args[2], 123)
            self.assertEqual(
                _at(result, "target-heatmap", "figure", "data", 0, "z", 1, 0), 0
            )

    def test_update_applies_parameters_together(
        self,
    ) -> None:
        values = self._values()
        applied = self._post(
            "applied-parameters.data", values, "update-button.n_clicks"
        )
        values[("applied-parameters", "data")] = _at(
            applied, "applied-parameters", "data"
        )
        before = self._post("target-heatmap.figure", values, "applied-parameters.data")
        values[("threshold", "value")] = 2
        values[("measure", "value")] = "count"
        values[("heatmap-view", "value")] = "drug"
        status = self._post("update-status.children", values, "threshold.value")
        self.assertIn(
            "not applied", _json_string(_at(status, "update-status", "children"))
        )
        unchanged = self._post(
            "target-heatmap.figure", values, "applied-parameters.data"
        )
        self.assertEqual(before, unchanged)
        applied = self._post(
            "applied-parameters.data", values, "update-button.n_clicks"
        )
        values[("applied-parameters", "data")] = _at(
            applied, "applied-parameters", "data"
        )
        after = self._post("target-heatmap.figure", values, "applied-parameters.data")
        self.assertNotEqual(before, after)
        self.assertEqual(
            _at(after, "target-heatmap", "figure", "data", 0, "y"),
            ["T cell (group)", "B cell (group)"],
        )
        self.assertIn("CPM ≥ 2", json.dumps(after["matrix-note"], ensure_ascii=False))
        panels = self._post(
            "target-heatmap-panel.hidden", values, "applied-parameters.data"
        )
        hidden = _at(panels, "target-heatmap-panel", "hidden")
        assert isinstance(hidden, bool)
        self.assertTrue(hidden)
        selectors = self._post(
            "detail-disease.options", values, "applied-parameters.data"
        )
        self.assertEqual(_at(selectors, "detail-disease", "value"), "MONDO_RA_TEST")
        details = self._post("details.children", values, "detail-disease.value")
        self.assertIn(
            "Drug–target records for rheumatoid arthritis",
            json.dumps(details, ensure_ascii=False),
        )
        status = self._post("update-status.children", values, "applied-parameters.data")
        self.assertEqual(_at(status, "update-status", "children"), "Settings applied")
        parameter_ids = {
            "threshold",
            "specificity",
            "modality",
            "target-class",
            "target-location",
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
            callback = self._application().callback_map[self._callback_key(output)]
            self.assertFalse(
                parameter_ids
                & {item["id"] for item in callback["inputs"] + callback["state"]}
            )

    def test_target_filters_narrow_records_and_appear_in_note(self) -> None:
        values = self._values()
        values[("target-class", "value")] = "Membrane receptor"
        values[("target-location", "value")] = "cell_surface"
        after = self._post_applied("target-heatmap.figure", values)
        note = json.dumps(after["matrix-note"], ensure_ascii=False)
        self.assertIn("Target class: Membrane receptor", note)
        self.assertIn("Target location: Cell surface", note)
        details = json.dumps(self._select_details(values), ensure_ascii=False)
        self.assertIn("TARGET2", details)
        self.assertNotIn("TARGET1", details)
        self.assertIn("Cell surface", details)

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
        app = application.create_app(snapshot, assets_folder=ASSETS_PATH)
        assert hasattr(app, "callback_map") and hasattr(app, "server")
        typed_app = cast(DashApplication, cast(object, app))
        type(self).application = typed_app
        type(self).client = typed_app.server.test_client()
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
        self.assertEqual(_at(synced, "diseases", "value"), expected)
        self.assertEqual(
            _at(synced, "disease-family-MONDO_0008383", "value"),
            ["MONDO_0008383"],
        )
        self.assertEqual(
            _at(synced, "disease-details-MONDO_0008383", "value"),
            ["EFO_0009459"],
        )
        values[("diseases", "value")] = _at(synced, "diseases", "value")
        values[("disease-family-MONDO_0008383", "value")] = []
        synced = self._post(
            "diseases.value", values, "disease-family-MONDO_0008383.value"
        )
        self.assertEqual(
            _at(synced, "diseases", "value"),
            ["EFO_0009459", "MONDO_0007915"],
        )
        values[("diseases", "value")] = _at(synced, "diseases", "value")
        values[("disease-family-MONDO_0008383", "value")] = ["MONDO_0008383"]
        synced = self._post(
            "diseases.value", values, "disease-family-MONDO_0008383.value"
        )
        self.assertEqual(_at(synced, "diseases", "value"), expected)
        values[("diseases", "value")] = _at(synced, "diseases", "value")
        synced = self._post(
            "diseases.value", values, "disease-details-MONDO_0008383.value"
        )
        self.assertEqual(
            _at(synced, "diseases", "value"),
            ["MONDO_0008383", "MONDO_0007915"],
        )
        selected_diseases = _json_array(_at(synced, "diseases", "value"))
        values[("diseases", "value")] = list(reversed(selected_diseases))
        figures = self._post_applied("target-heatmap.figure", values)
        for graph in ("target-heatmap", "drug-heatmap"):
            self.assertEqual(
                _at(figures, graph, "figure", "data", 0, "x"),
                ["rheumatoid arthritis", "systemic lupus erythematosus"],
            )
        values[("diseases", "value")] = []
        cleared = self._post("diseases.value", values, "diseases.value")
        self.assertTrue(all(_at(result, "value") == [] for result in cleared.values()))

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

    def test_expression_expansion_uses_cached_values_and_applied_rule(self):
        values: CallbackValues = {
            ("source-context", "data"): {
                "disease_id": "MONDO_RA_TEST",
                "modality": "all",
                "target_class": "all",
                "location": "all",
                "stage": "phase3",
                "threshold": 0.5,
                "method": "specificity",
                "specificity": 0.75,
            },
            ("expanded-expression-groups", "data"): [],
            ("expression-chart-type", "value"): "heatmap",
        }
        collapsed = self._post(
            "expression-heatmap.figure", values, "source-context.data"
        )
        values[("expanded-expression-groups", "data")] = ["group:CL_B_GROUP"]
        with patch.object(
            drugs_callbacks,
            "expression_figure",
            side_effect=AssertionError("Expansion recomputed expression"),
        ):
            expanded = self._post(
                "expression-heatmap.figure", values, "expanded-expression-groups.data"
            )
        for result in (collapsed, expanded):
            graph = _at(result, "expression-heatmap")
            self.assertEqual(
                _at(graph, "style", "height"),
                f"{_at(graph, 'figure', 'layout', 'height')}px",
            )
        before = _at(collapsed, "expression-heatmap", "figure", "data", 0)
        after = _at(expanded, "expression-heatmap", "figure", "data", 0)
        self.assertEqual(len(_json_array(_at(before, "y"))), 2)
        self.assertEqual(len(_json_array(_at(after, "y"))), 4)
        self.assertEqual(_at(before, "zmax"), _at(after, "zmax"))
        marks = next(
            _json_object(trace)
            for trace in _json_array(
                _at(expanded, "expression-heatmap", "figure", "data")
            )
            if _json_object(trace).get("name") == "Meets expression rule"
        )
        self.assertEqual(
            set(_json_string(item) for item in _json_array(marks["y"])),
            {"memory B cell", "naive B cell"},
        )
        self.assertEqual(
            len(_json_array(_at(expanded, "expression-row-controls", "children"))),
            4,
        )
        values[("expression-chart-type", "value")] = "dot"
        with patch.object(
            drugs_callbacks,
            "expression_figure",
            side_effect=AssertionError("Chart switch recomputed expression"),
        ):
            dot = self._post(
                "expression-heatmap.figure", values, "expression-chart-type.value"
            )
        self.assertEqual(
            _at(dot, "expression-heatmap", "figure", "data", 0, "name"),
            "CELLEX specificity",
        )
        self.assertEqual(_at(dot, "expression-dot-key", "style"), {"display": "flex"})
        dot_key = json.dumps(dot["expression-dot-key"], ensure_ascii=False)
        self.assertIn("CELLEX specificity (area)", dot_key)
        self.assertIn('"width": "3px"', dot_key)
        self.assertIn('"width": "11px"', dot_key)
        self.assertIn('"width": "15.5563px"', dot_key)
        self.assertIn('"width": "22px"', dot_key)
        self.assertIn("Empty cell: missing CPM or CELLEX", dot_key)
        self.assertEqual(
            _at(dot, "expression-heatmap-key", "style"), {"display": "none"}
        )

    def test_source_records_grid_keeps_all_rows_and_enables_sort_and_filter(
        self,
    ) -> None:
        records = [
            record
            for record in atlas.filtered_records(self.snapshot, "all", "phase3")
            if record["disease_id"] == "MONDO_RA_TEST"
        ]
        salt: FilteredRecord = {
            **records[0],
            "drug_id": "CHEMBL_SALT",
            "drug": "salt form",
            "drug_type": "Salt",
        }
        records = [*records, salt]
        rows = components.evidence_rows(records)
        grid = cast(GridComponent, cast(object, components.evidence_grid(rows)))
        self.assertEqual(grid.id, "source-records-grid")
        expected_pairs = list(
            dict.fromkeys(
                (record["canonical_drug_id"], record["target_id"]) for record in records
            )
        )
        self.assertEqual(
            [
                (row["canonical_drug"].rsplit("(", 1)[1][:-1], row["target"])
                for row in grid.rowData
            ],
            [
                (
                    canonical_id,
                    f"{next(r['target'] for r in records if r['target_id'] == target_id) or '—'} ({target_id or '—'})",
                )
                for canonical_id, target_id in expected_pairs
            ],
        )
        self.assertEqual(len(rows), len(records) - 1)
        self.assertIn("Salt", rows[0]["modality"])
        self.assertNotIn("salt form", str(rows))
        self.assertEqual(
            [column["headerName"] for column in grid.columnDefs],
            [
                "Canonical drug",
                "Modality",
                "Canonical stage",
                "Target",
                "Target class",
                "Target location",
                "Action / mechanism",
                "Open Targets links",
            ],
        )
        self.assertTrue(grid.defaultColDef["sortable"])
        self.assertTrue(grid.defaultColDef["floatingFilter"])
        self.assertEqual(grid.defaultColDef["filter"], "agTextColumnFilter")
        links = grid.columnDefs[-1]
        self.assertEqual(links["cellRenderer"], "markdown")
        self.assertEqual(links["linkTarget"], "_blank")
        self.assertFalse(links["sortable"])
        self.assertTrue(grid.dashGridOptions["pagination"])
        self.assertEqual(
            grid.dashGridOptions["paginationPageSize"], config.SOURCE_PAGE_SIZE
        )
        self.assertIn("[Open Targets drug](", grid.rowData[0]["links"])
        self.assertNotIn(
            "source-records-page.children", self._application().callback_map
        )

    def test_layout_has_all_modalities_and_defaults(self) -> None:
        self.assertEqual(self._client().get("/").status_code, 200)
        self.assertEqual(_at(self.components["stage"], "value"), "phase3")
        self.assertEqual(_at(self.components["measure"], "value"), "percent")
        self.assertEqual(
            [
                _at(option, "value")
                for option in _json_array(_at(self.components["measure"], "options"))
            ],
            ["percent", "count"],
        )
        self.assertEqual(
            _at(self.components["applied-parameters"], "data", "measure"),
            "percent",
        )
        self.assertEqual(_at(self.components["method"], "value"), "specificity")
        self.assertEqual(
            _at(self.components["applied-parameters"], "data", "method"),
            "specificity",
        )
        self.assertEqual(_at(self.components["specificity"], "value"), 0.5)
        self.assertEqual(
            _at(self.components["applied-parameters"], "data", "specificity"), 0.5
        )
        self.assertNotIn("level", self.components)
        self.assertNotIn("detail-cell", self.components)
        self.assertNotIn("download-button", self.components)
        self.assertNotIn("download", self.components)
        self.assertFalse(
            any("download" in key for key in self._application().callback_map)
        )
        self.assertNotIn("cells", self.components)
        self.assertNotIn("cells", config.PARAMETER_IDS)
        self.assertEqual(
            [
                _at(option, "value")
                for option in _json_array(_at(self.components["modality"], "options"))
            ],
            ["all", *[value for _, value in atlas.DRUG_TYPE_MODALITIES]],
        )
        self.assertEqual(
            [
                _at(option, "value")
                for option in _json_array(
                    _at(self.components["target-class"], "options")
                )
            ],
            ["all", "Enzyme", "Membrane receptor", "unknown"],
        )
        self.assertEqual(
            [
                _at(option, "value")
                for option in _json_array(
                    _at(self.components["target-location"], "options")
                )
            ],
            ["all", *[value for value, _ in atlas.LOCATION_CLASSES]],
        )
        self.assertEqual(
            _at(self.components["applied-parameters"], "data", "target-location"),
            "all",
        )
        layout = _response_json(self._client().get("/_dash-layout"))
        drug_page = next(
            node
            for node in _json_array(_at(layout, "props", "children"))
            if _at(node, "props", "id") == "drug-page"
        )
        controls = next(
            _json_object(node)
            for node in _json_array(_at(drug_page, "props", "children"))
            if _at(node, "props", "className") == "panel controls"
        )
        self.assertEqual(
            _at(controls, "props", "children", 0, "props", "children", 0),
            "Settings",
        )
        groups = _at(controls, "props", "children", 1)
        scope = _at(groups, "props", "children", 0)
        scope_controls = _at(scope, "props", "children", 1)
        self.assertEqual(_at(scope_controls, "type"), "Div")
        self.assertEqual(
            [
                _at(node, "props", "children", 1, "props", "id")
                for node in _json_array(_at(scope_controls, "props", "children"))
            ],
            ["diseases"],
        )
        self.assertEqual(
            [
                _at(group, "props", "children", 0, "props", "children")
                for group in _json_array(_at(groups, "props", "children"))
                if "filter-group"
                in _json_string(
                    _json_object(_at(group, "props")).get("className", "")
                ).split()
            ],
            ["Comparison scope", "Drug evidence", "Expression criteria", "Display"],
        )
        browser = _at(scope_controls, "props", "children", 0, "props", "children", 2)
        self.assertEqual(_at(browser, "type"), "Div")
        self.assertEqual(_at(browser, "props", "children", 0, "type"), "H4")
        disabled = _at(self.components["specificity"], "disabled")
        assert isinstance(disabled, bool)
        self.assertFalse(disabled)

    def test_chart_type_switches_immediately_and_preserves_applied_filters(
        self,
    ) -> None:
        self.assertEqual(_at(self.components["chart-type"], "value"), "heatmap")
        self.assertEqual(
            [
                _at(option, "value")
                for option in _json_array(_at(self.components["chart-type"], "options"))
            ],
            ["heatmap", "dot"],
        )
        self.assertNotIn("chart-type", config.PARAMETER_IDS)
        values = self._values()
        values[("applied-parameters", "data")] = {
            **_json_object(values[("applied-parameters", "data")]),
            "threshold": 0.71,
        }
        self._post("target-heatmap.figure", values, "applied-parameters.data")
        values[("chart-type", "value")] = "dot"
        with patch.object(atlas, "summarize", wraps=atlas.summarize) as summarize:
            result = self._post("target-heatmap.figure", values, "chart-type.value")
        summarize.assert_not_called()
        self.assertEqual(
            _at(result, "target-heatmap", "figure", "data", 0, "type"), "scatter"
        )
        self.assertEqual(
            _at(result, "target-heatmap", "figure", "data", 0, "customdata", 0),
            ["MONDO_RA_TEST", "group:" + atlas.T_CELL_ID],
        )
        self.assertEqual(_at(result, "target-dot-key", "style"), {"display": "flex"})
        target_key = json.dumps(result["target-dot-key"], ensure_ascii=False)
        self.assertIn("Count (area)", target_key)
        self.assertIn('"width": "22px"', target_key)
        self.assertIn("Empty cell: zero or missing", target_key)
        details = self._select_details(values)
        self.assertIn("Drug–target records for rheumatoid arthritis", str(details))

    def test_dot_plot_disables_measure_without_changing_its_value(self) -> None:
        values = self._values()
        for chart_type, disabled in (("dot", True), ("heatmap", False)):
            values[("chart-type", "value")] = chart_type
            result = self._post("measure.options", values, "chart-type.value")
            self.assertEqual(
                result["measure"],
                {
                    "options": [
                        {"label": "Percent", "value": "percent", "disabled": disabled},
                        {"label": "Count", "value": "count", "disabled": disabled},
                    ]
                },
            )
            self.assertEqual(values[("measure", "value")], "percent")

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
            json.dumps(_at(figures, "matrix-note", "children"), ensure_ascii=False),
        )

    def test_help_buttons_describe_unique_tooltips_and_errors_stay_visible(
        self,
    ) -> None:
        values = self._values()
        values[("threshold", "value")] = -1
        figures = self._post_applied("target-heatmap.figure", values)
        note = _json_array(_at(figures, "matrix-note", "children", "props", "children"))
        self.assertEqual(
            _at(note[0], "props", "children"), "2 disease–cell combinations"
        )
        self.assertEqual(_at(note[2], "props", "role"), "alert")
        self.assertIn("must be finite", _json_string(_at(note[2], "props", "children")))

        def nodes(item: JsonValue) -> Iterator[JsonObject]:
            if isinstance(item, dict):
                if "props" in item:
                    yield item
                for value in item.values():
                    yield from nodes(value)
            elif isinstance(item, list):
                for value in item:
                    yield from nodes(value)

        layout = _response_json(self._client().get("/_dash-layout"))
        details = self._select_details(self._values())
        surfaces: list[JsonValue] = [
            layout,
            _at(figures, "matrix-note", "children"),
            _at(details, "details", "children"),
        ]
        tips: list[str] = []
        buttons: list[JsonObject] = []
        for node in nodes(surfaces):
            props = _json_object(node["props"])
            if props.get("role") == "tooltip":
                tips.append(_json_string(props["id"]))
            if node.get("type") == "Button" and props.get("className") == "info-button":
                buttons.append(props)
        descriptions = {_json_string(button["aria-describedby"]) for button in buttons}
        self.assertEqual(len(tips), len(set(tips)))
        self.assertEqual(set(tips), descriptions)
        self.assertTrue(
            all(
                button["type"] == "button" and bool(_json_string(button["aria-label"]))
                for button in buttons
            )
        )

    def test_both_heatmaps_callbacks(self) -> None:
        values = self._values()
        figures = self._post_applied("target-heatmap.figure", values)
        for graph in ("target-heatmap", "drug-heatmap"):
            self.assertEqual(_at(figures, graph, "style", "minWidth"), "600px")
        self.assertEqual(
            _at(figures, "target-heatmap", "figure", "data", 0, "x"),
            _at(figures, "drug-heatmap", "figure", "data", 0, "x"),
        )
        self.assertIn(
            "CPM ≥ 0.5",
            json.dumps(_at(figures, "matrix-note", "children"), ensure_ascii=False),
        )
        self.assertNotIn(
            "specificity ≥",
            json.dumps(_at(figures, "matrix-note", "children"), ensure_ascii=False),
        )

        values[("method", "value")] = "relative"
        relative = self._post_applied("target-heatmap.figure", values)
        self.assertIn(
            "target-relative median",
            json.dumps(_at(relative, "matrix-note", "children"), ensure_ascii=False),
        )

    def test_heatmap_view_switch_keeps_both_graphs_and_details(self) -> None:
        self.assertEqual(_at(self.components["heatmap-view"], "value"), "target")
        self.assertFalse(self.components["target-heatmap-panel"].get("hidden", False))
        initial_hidden = self.components["drug-heatmap-panel"]["hidden"]
        assert isinstance(initial_hidden, bool)
        self.assertTrue(initial_hidden)
        values = self._values()
        figures = self._post_applied("target-heatmap.figure", values)
        self.assertIn("figure", _json_object(figures["target-heatmap"]))
        self.assertIn("figure", _json_object(figures["drug-heatmap"]))
        for view, hidden in (
            ("target", (False, True)),
            ("drug", (True, False)),
        ):
            values[("heatmap-view", "value")] = view
            panels = self._post_applied("target-heatmap-panel.hidden", values)
            self.assertEqual(
                (
                    _at(panels, "target-heatmap-panel", "hidden"),
                    _at(panels, "drug-heatmap-panel", "hidden"),
                ),
                hidden,
            )
            details = self._select_details(values)
            self.assertIn(
                "Relative target expression by cell type",
                str(_at(details, "details", "children")),
            )

    def test_figures_keep_cell_lineage_order(self) -> None:
        values = self._values()
        figures = self._post_applied("target-heatmap.figure", values)
        self.assertEqual(
            _at(figures, "target-heatmap", "figure", "data", 0, "y"),
            ["T cell (group)", "B cell (group)"],
        )
        self.assertEqual(
            _at(figures, "drug-heatmap", "figure", "data", 0, "y"),
            ["T cell (group)", "B cell (group)"],
        )
        selectors = self._post_applied("detail-disease.options", values)
        self.assertEqual(
            [
                _at(item, "value")
                for item in _json_array(_at(selectors, "detail-disease", "options"))
            ],
            ["MONDO_RA_TEST"],
        )

    def test_details_follow_disease_selector(self) -> None:
        details = self._select_details(self._values())
        rendered = str(_at(details, "details", "children"))
        self.assertIn("Relative target expression by cell type", rendered)
        self.assertIn("Drug–target records for rheumatoid arthritis", rendered)

    def test_detail_selector_ignores_heatmap_clicks(self) -> None:
        callback = self._application().callback_map[
            self._callback_key("detail-disease.value")
        ]
        self.assertNotIn("clickData", [item["property"] for item in callback["inputs"]])
