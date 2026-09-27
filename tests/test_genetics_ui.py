"""遺伝子ページの図、部品、レイアウト、callback を検証する。"""

from __future__ import annotations

import json
import unittest
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import ClassVar, Protocol, TypedDict, cast, override

from flask.testing import FlaskClient

from autoimmune_atlas import genetics
from autoimmune_atlas.models import DiseaseFamily, SummaryRow
from autoimmune_atlas.ui import components, figures, genetics_layout, layout
from autoimmune_atlas.ui.application import create_app
from autoimmune_atlas.ui.genetics_callbacks import effective_score
from tests.genetics_fixture import genetics_snapshot
from tests.test_genetics import snapshot as core_ui_snapshot

ASSETS_PATH = Path(__file__).resolve().parents[1] / "assets"


class _Title(Protocol):
    text: str


class _ColorBar(Protocol):
    title: _Title


class _Heatmap(Protocol):
    colorbar: _ColorBar


class _Figure(Protocol):
    data: Sequence[_Heatmap]


class _Button(Protocol):
    id: dict[str, str]


class _NumberInput(Protocol):
    min: float
    max: float
    value: float


class _Grid(Protocol):
    columnDefs: list[dict[str, object]]


class _Classed(Protocol):
    className: str


class _Identified(Protocol):
    id: str


class _Linked(Protocol):
    href: str


def gene_row(**overrides: object) -> SummaryRow:
    base: dict[str, object] = {
        "disease_id": "D1",
        "disease": "Disease one",
        "cell_id": "group:LYMPH",
        "cell": "Lymphocyte (group)",
        "ontology_id": "LYMPH",
        "cell_level": "group",
        "count": 3,
        "percent": 75.0,
        "denominator": 4,
        "unknown": 0,
        "unmapped_drugs": 0,
        "status": "complete",
        "records": [],
        "drug_count": None,
        "drug_percent": None,
        "drug_denominator": 0,
        "unknown_drugs": 0,
        "total_drugs": 0,
        "mapped_drugs": 0,
        "member_cell_ids": ["T4", "T8"],
    }
    base.update(overrides)
    return cast(SummaryRow, cast(object, base))


class GeneFigureTests(unittest.TestCase):
    def test_gene_hover_names_genes_and_omits_drug_coverage(self) -> None:
        text = figures.hover_text(gene_row(), "percent", "gene")
        self.assertIn("Genes: 75.0%", text)
        self.assertNotIn("Canonical drugs", text)
        dot = figures.dot_hover_text(gene_row(), "gene")
        self.assertIn("Count: 3", dot)
        self.assertNotIn("Canonical drugs", dot)

    def test_zero_genes_above_threshold_has_its_own_message(self) -> None:
        row = gene_row(count=0, percent=None, denominator=0)
        text = figures.hover_text(row, "count", "gene")
        self.assertIn("No genes at or above the score threshold", text)
        text = figures.hover_text(row, "percent", "gene")
        self.assertIn("No genes at or above the score threshold", text)

    def test_gene_figure_title(self) -> None:
        figure = figures.build_figure(
            [gene_row()], ["D1"], ["group:LYMPH"], "count", "gene"
        )
        heatmap = cast(_Figure, cast(object, figure)).data[0]
        self.assertEqual(heatmap.colorbar.title.text, "Genes (count)")


class PrefixedComponentTests(unittest.TestCase):
    def test_checklist_sections_take_prefix(self) -> None:
        family: DiseaseFamily = {
            "id": "ROOT",
            "label": "Root",
            "diseases": [
                {"id": "ROOT", "name": "Root"},
                {"id": "CHILD", "name": "Child", "parent_ids": ["ROOT"]},
            ],
        }
        ids = [
            s["id"] for s in components.disease_checklist_sections(family, "genetics-")
        ]
        self.assertEqual(
            ids, ["genetics-disease-family-ROOT", "genetics-disease-details-ROOT"]
        )
        plain = [s["id"] for s in components.disease_checklist_sections(family)]
        self.assertEqual(plain, ["disease-family-ROOT", "disease-details-ROOT"])

    def test_row_controls_use_kind_specific_toggle_type(self) -> None:
        names = {"group:LYMPH": "Lymphocyte (group)"}
        button = components.heatmap_row_controls(
            names, ["group:LYMPH"], [], "genetics"
        )[0]
        self.assertEqual(
            cast(_Button, cast(object, button)).id["type"], "genetics-cell-toggle"
        )
        button = components.heatmap_row_controls(names, ["group:LYMPH"], [], "target")[
            0
        ]
        self.assertEqual(
            cast(_Button, cast(object, button)).id["type"], "heatmap-cell-toggle"
        )


def _walk(node: object) -> Iterator[object]:
    yield node
    children = getattr(node, "children", None)
    if isinstance(children, list):
        for child in cast(list[object], children):
            yield from _walk(child)
    elif children is not None:
        yield from _walk(cast(object, children))


class GeneticsLayoutTests(unittest.TestCase):
    def test_page_has_prefixed_controls_and_nav_marks_current(self) -> None:
        page = genetics_layout.genetics_page(core_ui_snapshot(), genetics_snapshot())
        found = {
            cast(_Identified, c).id
            for c in _walk(page)
            if isinstance(getattr(c, "id", None), str)
        }
        for component_id in (
            "genetics-page",
            "genetics-diseases",
            "genetics-score",
            "genetics-method",
            "genetics-threshold",
            "genetics-specificity",
            "genetics-measure",
            "genetics-update-button",
            "genetics-applied-parameters",
            "genetics-heatmap",
            "genetics-chart-type",
            "genetics-detail-disease",
            "genetics-details",
        ):
            self.assertIn(component_id, found)
        self.assertNotIn("diseases", found)
        self.assertNotIn("update-button", found)
        self.assertTrue({"genetics-comparison-tip", "genetics-selection-tip"} <= found)
        self.assertEqual([c for c in found if not c.startswith("genetics-")], [])
        current = [c for c in _walk(page) if getattr(c, "aria-current", None) == "page"]
        self.assertEqual([cast(_Linked, c).href for c in current], ["/genetics"])

    def test_score_input_range_and_default(self) -> None:
        page = genetics_layout.genetics_page(core_ui_snapshot(), genetics_snapshot())
        score = cast(
            _NumberInput,
            next(c for c in _walk(page) if getattr(c, "id", None) == "genetics-score"),
        )
        self.assertEqual(score.min, 0.1)
        self.assertEqual(score.max, 1)
        self.assertEqual(score.value, 0.5)

    def test_gene_rows_have_scores_datasource_columns_and_evidence_link(self) -> None:
        data = genetics_snapshot()
        genes = genetics.genes_for_disease(data, "D1", 0.5)
        rows = genetics_layout.gene_rows(genes, "D1", data["datasources"])
        self.assertEqual(rows[0]["gene"], "Gene 1 (G1)")
        self.assertEqual(rows[0]["score"], "0.900")
        self.assertEqual(rows[0]["gwas_credible_sets"], "0.900")
        self.assertEqual(rows[0]["eva"], "—")
        self.assertIn(
            "https://platform.opentargets.org/evidence/G1/D1", rows[0]["links"]
        )
        grid = cast(
            _Grid, cast(object, genetics_layout.gene_grid(rows, data["datasources"]))
        )
        fields = [c["field"] for c in grid.columnDefs]
        self.assertEqual(
            fields, ["gene", "score", "gwas_credible_sets", "eva", "links"]
        )

    def test_detail_panel_without_selection_and_unavailable_page(self) -> None:
        panel = genetics_layout.genetics_detail_panel(
            [],
            None,
            genetics_snapshot(),
            score_threshold=0.5,
            threshold=0.5,
            method="fixed",
            specificity=0.5,
        )
        self.assertEqual(cast(_Classed, cast(object, panel)).className, "empty-note")
        page = genetics_layout.genetics_unavailable_page("versions differ")
        self.assertIn("refresh-genetics", str(page))


class _Dependency(Protocol):
    component_id: str
    component_property: str


class _CallbackInput(TypedDict):
    id: str
    property: str


class _CallbackDefinition(TypedDict):
    inputs: list[_CallbackInput]
    output: _Dependency | list[_Dependency]
    state: list[_CallbackInput]


class _Server(Protocol):
    def test_client(self) -> FlaskClient: ...


class _Application(Protocol):
    callback_map: dict[str, _CallbackDefinition]
    server: _Server


def _post(
    test: unittest.TestCase,
    application: _Application,
    client: FlaskClient,
    output_id: str,
    values: dict[str, object],
) -> dict[str, object]:
    """test_ui.py の CallbackTests._post と同じ形で callback を 1 回呼ぶ。"""
    key = next(k for k in application.callback_map if output_id in k)
    callback = application.callback_map[key]
    outputs = callback["output"]
    response = client.post(
        "/_dash-update-component",
        json={
            "output": key,
            "outputs": [
                {"id": o.component_id, "property": o.component_property}
                for o in outputs
            ]
            if isinstance(outputs, list)
            else {
                "id": outputs.component_id,
                "property": outputs.component_property,
            },
            "inputs": [
                {**item, "value": values.get(f"{item['id']}.{item['property']}")}
                for item in callback["inputs"]
            ],
            "state": [
                {**item, "value": values.get(f"{item['id']}.{item['property']}")}
                for item in callback["state"]
            ],
            "changedPropIds": [
                f"{item['id']}.{item['property']}" for item in callback["inputs"]
            ][:1],
        },
    )
    try:
        test.assertEqual(response.status_code, 200, response.get_data())
        body = cast(dict[str, object], json.loads(response.get_data(as_text=True)))
        return cast(dict[str, object], body["response"])
    finally:
        response.close()


def _layout_json(app: object) -> object:
    """/_dash-layout の JSON を返す。"""
    response = cast(_Application, app).server.test_client().get("/_dash-layout")
    try:
        return cast(object, json.loads(response.get_data(as_text=True)))
    finally:
        response.close()


def _collect_ids(node: object) -> set[str]:
    """レイアウトの JSON から、文字列の id を全部集める。"""
    found: set[str] = set()
    if isinstance(node, dict):
        mapping = cast(dict[str, object], node)
        component_id = mapping.get("id")
        if isinstance(component_id, str):
            found.add(component_id)
        for value in mapping.values():
            found |= _collect_ids(value)
    elif isinstance(node, list):
        for item in cast(list[object], node):
            found |= _collect_ids(item)
    return found


class GeneticsCallbackTests(unittest.TestCase):
    """HTTP 経由で遺伝子ページの比較図と詳細の callback を確かめる。"""

    application: ClassVar[_Application]
    client: ClassVar[FlaskClient]

    @classmethod
    @override
    def setUpClass(cls) -> None:
        app = create_app(
            core_ui_snapshot(), assets_folder=ASSETS_PATH, genetics=genetics_snapshot()
        )
        cls.application = cast(_Application, cast(object, app))
        cls.client = cls.application.server.test_client()

    def _post(self, output_id: str, values: dict[str, object]) -> dict[str, object]:
        return _post(self, self.application, self.client, output_id, values)

    def test_score_input_falls_back_below_floor_and_on_empty(self) -> None:
        self.assertEqual(
            effective_score(0.05),
            (0.5, "Score threshold must be finite and between 0.1 and 1; using 0.5."),
        )
        self.assertEqual(effective_score(""), (0.5, None))
        self.assertEqual(effective_score(0.3), (0.3, None))

    def test_update_applies_parameters_and_details_follow(self) -> None:
        applied = {
            "genetics-measure": "count",
            "genetics-score": 0.5,
            "genetics-method": "fixed",
            "genetics-threshold": 0.5,
            "genetics-specificity": 0.5,
            "genetics-diseases": ["D1", "D2"],
        }
        response = self._post(
            "genetics-heatmap.figure",
            {
                "genetics-applied-parameters.data": applied,
                "genetics-expanded-cell-groups.data": [],
                "genetics-chart-type.value": "heatmap",
            },
        )
        heatmap = cast(dict[str, dict[str, object]], response["genetics-heatmap"])
        data = cast(list[dict[str, object]], heatmap["figure"]["data"])
        self.assertEqual(data[0]["x"], ["Disease one", "Disease two"])
        note = cast(dict[str, object], response["genetics-matrix-note"])["children"]
        self.assertIn("Score threshold", json.dumps(note))
        details = self._post(
            "genetics-details.children",
            {
                "genetics-detail-disease.value": "D1",
                "genetics-applied-parameters.data": applied,
            },
        )
        self.assertIn(
            "Genetically associated genes for Disease one", json.dumps(details)
        )
        empty = self._post(
            "genetics-details.children",
            {
                "genetics-detail-disease.value": "D2",
                "genetics-applied-parameters.data": applied,
            },
        )
        self.assertIn("No genes at or above the score threshold", json.dumps(empty))


class RouterTests(unittest.TestCase):
    """URL で 2 ページを切り替える shell を確かめる。"""

    def test_layout_holds_both_pages_and_url_toggles_hidden(self) -> None:
        app = cast(
            _Application,
            cast(
                object,
                create_app(
                    core_ui_snapshot(),
                    assets_folder=ASSETS_PATH,
                    genetics=genetics_snapshot(),
                ),
            ),
        )
        client = app.server.test_client()
        ids = _collect_ids(_layout_json(app))
        for component_id in (
            "url",
            "drug-page",
            "genetics-page",
            "diseases",
            "genetics-diseases",
        ):
            self.assertIn(component_id, ids)
        for path, expected in (
            ("/genetics", [True, False]),
            ("/", [False, True]),
            (None, [False, True]),
        ):
            response = cast(
                dict[str, dict[str, object]],
                _post(self, app, client, "drug-page.hidden", {"url.pathname": path}),
            )
            self.assertEqual(
                [response["drug-page"]["hidden"], response["genetics-page"]["hidden"]],
                expected,
            )
        for route in ("/", "/genetics"):
            page = client.get(route)
            try:
                self.assertEqual(page.status_code, 200)
            finally:
                page.close()

    def test_version_mismatch_and_missing_genetics_show_notice_but_keep_drug_page(
        self,
    ) -> None:
        stale = genetics_snapshot()
        stale["data_version"] = {"year": "26", "month": "06", "iteration": None}
        for data in (stale, None):
            app = create_app(
                core_ui_snapshot(), assets_folder=ASSETS_PATH, genetics=data
            )
            layout_data = _layout_json(app)
            ids = _collect_ids(layout_data)
            self.assertIn("drug-page", ids)
            self.assertIn("genetics-page", ids)
            self.assertNotIn("genetics-heatmap", ids)
            self.assertIn("refresh-genetics", json.dumps(layout_data))

    def test_missing_snapshot_keeps_both_pages(self) -> None:
        app = create_app(None, assets_folder=ASSETS_PATH, genetics=genetics_snapshot())
        ids = _collect_ids(_layout_json(app))
        self.assertTrue({"url", "drug-page", "genetics-page"} <= ids)
        self.assertNotIn("genetics-heatmap", ids)

    def test_nav_marks_current_page(self) -> None:
        drug = layout.dashboard_layout(core_ui_snapshot())
        links = [
            c for c in _walk(drug) if getattr(c, "className", None) == "app-page-link"
        ]
        self.assertEqual([cast(_Linked, c).href for c in links], ["/", "/genetics"])
        current = [c for c in links if getattr(c, "aria-current", None) == "page"]
        self.assertEqual([cast(_Linked, c).href for c in current], ["/"])
        self.assertEqual(getattr(drug, "id", None), "drug-page")
