"""schema 2 の Dash UI と表示用ロジックを検証する。"""

from __future__ import annotations

import json
import math
import os
import runpy
import tempfile
import unittest
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import ClassVar, Literal, Protocol, TypedDict, Unpack, cast, override
from unittest.mock import patch

from flask.testing import FlaskClient
from werkzeug.test import TestResponse

from backend import aggregation as atlas
from backend import refresh, snapshot
from backend.models import (
    Disease,
    DiseaseFamily,
    DrugRecord,
    ExpressionStateInput,
    FilteredRecord,
    RootIdentity,
    Snapshot,
    SummaryRow,
    TargetRecord,
)
from backend.ui import (
    application,
    callbacks,
    components,
    config,
    figures,
    layout,
)

ASSETS_PATH = Path(__file__).resolve().parents[1] / "assets"


type JsonScalar = bool | float | int | str | None
type JsonValue = JsonScalar | list[JsonValue] | dict[str, JsonValue]
type JsonObject = dict[str, JsonValue]
type CallbackValues = dict[tuple[str, str], JsonValue]


class SummaryOverrides(TypedDict, total=False):
    disease_id: str
    disease: str
    cell_id: str
    cell: str
    count: int | None
    percent: float | None
    denominator: int
    unknown: int
    drug_count: int | None
    drug_percent: float | None
    drug_denominator: int
    unknown_drugs: int
    mapped_drugs: int
    total_drugs: int
    unmapped_drugs: int
    status: str
    member_cell_ids: list[str]


class FigureAxis(Protocol):
    autorange: bool | str | None
    categoryarray: Sequence[str] | None
    categoryorder: str | None
    range: Sequence[float | str] | None
    side: str | None
    tickangle: float | None
    ticktext: Sequence[str] | None
    type: str | None


class FigureLayout(Protocol):
    width: int | None
    xaxis: FigureAxis
    yaxis: FigureAxis


class FigureData(Protocol):
    data: tuple[object, ...]
    layout: FigureLayout


class ColorBar(Protocol):
    lenmode: str | None
    len: float | None
    y: float | None
    yanchor: str | None


class HeatmapTrace(Protocol):
    colorbar: ColorBar
    customdata: Sequence[Sequence[object]]
    hovertext: Sequence[Sequence[str]]
    name: str | None
    text: Sequence[Sequence[str]]
    x: Sequence[str]
    y: Sequence[str]
    z: Sequence[Sequence[float | None]]
    zmax: float | None
    zmin: float | None


class DotMarker(Protocol):
    cmax: float | None
    cmin: float | None
    color: Sequence[float] | str | None
    showscale: bool | None
    size: Sequence[float] | float | None
    sizemode: str | None
    sizeref: float | None
    symbol: str | None


class ScatterTrace(Protocol):
    customdata: Sequence[Sequence[str]] | None
    hovertext: Sequence[str] | None
    marker: DotMarker
    mode: str | None
    name: str | None
    text: str | Sequence[str] | None
    x: Sequence[str]
    y: Sequence[str]


class ChildrenComponent(Protocol):
    children: object


class ClassComponent(Protocol):
    className: str | None


class GraphComponent(Protocol):
    config: Mapping[str, object] | None
    id: str | dict[str, JsonValue] | None


class StoreComponent(Protocol):
    data: JsonValue


class CallbackDependency(Protocol):
    component_id: str
    component_property: str


class CallbackInput(TypedDict):
    id: str
    property: str


class CallbackDefinition(TypedDict):
    inputs: list[CallbackInput]
    output: CallbackDependency | list[CallbackDependency]
    state: list[CallbackInput]


class DashServer(Protocol):
    def test_client(self) -> FlaskClient: ...


class DashApplication(Protocol):
    callback_map: dict[str, CallbackDefinition]
    config: Mapping[str, object]
    server: DashServer


def _figure(figure: object) -> FigureData:
    assert hasattr(figure, "data") and hasattr(figure, "layout")
    return cast(FigureData, figure)


def _heatmap(figure: object) -> HeatmapTrace:
    traces = _figure(figure).data
    fields = (
        "colorbar",
        "customdata",
        "hovertext",
        "name",
        "text",
        "x",
        "y",
        "z",
        "zmax",
        "zmin",
    )
    assert traces and all(hasattr(traces[0], name) for name in fields)
    return cast(HeatmapTrace, traces[0])


def _scatter(figure: object, name: str) -> ScatterTrace:
    for trace in _figure(figure).data:
        if hasattr(trace, "name") and cast(ScatterTrace, trace).name == name:
            assert all(hasattr(trace, field) for field in ("x", "y", "mode", "text"))
            return cast(ScatterTrace, trace)
    raise AssertionError(f"trace not found: {name}")


def _scatters(figure: object) -> list[ScatterTrace]:
    traces = _figure(figure).data[1:]
    assert all(
        all(hasattr(trace, field) for field in ("name", "x", "y", "mode", "text"))
        for trace in traces
    )
    return [cast(ScatterTrace, trace) for trace in traces]


def _component(component: object) -> ChildrenComponent:
    assert hasattr(component, "children")
    return cast(ChildrenComponent, component)


def _component_list(component: object) -> list[object]:
    children = _component(component).children
    assert isinstance(children, list)
    return cast(list[object], children)


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


def fixture() -> Snapshot:
    """段階、成分重複、欠測、group union を含む最小 fixture を返す。"""
    root: RootIdentity = {
        "id": "MONDO_AUTOIMMUNE",
        "name": "autoimmune disease",
    }
    records: list[DrugRecord] = [
        {
            "disease_id": "MONDO_RA_TEST",
            "disease": "rheumatoid arthritis",
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
            "disease_id": "MONDO_RA_TEST",
            "disease": "rheumatoid arthritis",
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
            "disease_id": "MONDO_RA_TEST",
            "disease": "rheumatoid arthritis",
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
            "references": [{"source": "ClinicalTrials", "ids": ["NCT1"], "urls": []}],
        },
        {
            "disease_id": "MONDO_RA_TEST",
            "disease": "rheumatoid arthritis",
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
    ]
    return {
        "schema": 2,
        "root": root,
        "data_version": {"year": "26", "month": "6", "iteration": None},
        "retrieved_at": "2026-09-22T12:00:00Z",
        "source": "https://api.platform.opentargets.org/",
        "diseases": [
            {"id": "MONDO_RA_TEST", "name": "rheumatoid arthritis", "status": "ready"}
        ],
        "records": records,
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


def summary(**overrides: Unpack[SummaryOverrides]) -> SummaryRow:
    row: SummaryRow = {
        "disease_id": "D1",
        "disease": "Disease one",
        "cell_id": "C1",
        "cell": "B cell",
        "ontology_id": "C1",
        "cell_level": "group",
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
    merged = dict(row)
    merged.update(overrides)
    return cast(SummaryRow, cast(object, merged))


class DiseaseTreeTests(unittest.TestCase):
    def test_sections_follow_parents_and_place_each_term_once(self) -> None:
        family: DiseaseFamily = {
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
        sections = components.disease_checklist_sections(family)
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
        selector_snapshot = fixture()
        selector_snapshot["diseases"] = [
            {**d, "status": "ready"} for d in family["diseases"]
        ]
        layout = components.disease_selector(selector_snapshot, [])
        summaries: list[object] = []

        def walk(component: object) -> None:
            if not hasattr(component, "children"):
                return
            if type(component).__name__ == "Summary":
                summaries.append(_component(component).children)
            children = _component(component).children
            if children is None:
                return
            if not isinstance(children, list):
                children = [children]
            for child in cast(list[object], children):
                if not isinstance(child, str):
                    walk(child)

        walk(layout)
        self.assertIn("child (3 terms)", summaries)
        self.assertIn("child details (2 terms)", summaries)
        self.assertIn("grandchild (2 terms)", summaries)

    def test_family_without_root_keeps_flat_top_level(self) -> None:
        family: DiseaseFamily = {
            "id": "other",
            "label": "Other terms",
            "diseases": [
                {"id": "B", "name": "b term", "parent_ids": []},
                {"id": "A", "name": "a term", "parent_ids": ["MISSING"]},
            ],
        }
        sections = components.disease_checklist_sections(family)
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
            "Disease data unavailable",
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
        self.assertIn("Open Targets canonical drug", links)
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
        pagination = source_children[-1]
        assert hasattr(pagination, "className")
        self.assertEqual(
            cast(ClassComponent, pagination).className, "source-pagination"
        )
        table = _component(source_children[1]).children
        assert hasattr(table, "id")
        self.assertEqual(cast(GraphComponent, table).id, "source-records-page")
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
        missing = _scatter(figure, "Missing expression")
        self.assertEqual(missing.mode, "text")
        self.assertEqual(missing.text, "0")
        self.assertEqual(list(missing.x), ["TARGET2 (ENSG_TARGET_2)"])
        self.assertEqual(list(missing.y), ["memory B cell"])
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
        self.assertIn("○ Source cell meets expression rule", str(panel_children[2]))
        self.assertIn("expression-chart-type", str(panel_children[2]))
        self.assertIn("expression-dot-key", str(panel_children[2]))
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

    def test_selection_rejects_stale_click_and_accepts_either_heatmap(self) -> None:
        stale: callbacks.ClickData = {"points": [{"customdata": ["OLD", "CELL"]}]}
        valid: callbacks.ClickData = {
            "points": [{"customdata": ["MONDO_RA_TEST", "CL_B_GROUP"]}]
        }
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
                    namespace = cast(
                        dict[str, object],
                        runpy.run_path(str(entrypoint), run_name="app_path_test"),
                    )
            dash_app_value = namespace["app"]
            assert hasattr(dash_app_value, "server") and hasattr(
                dash_app_value, "config"
            )
            dash_app = cast(DashApplication, dash_app_value)
            assets_folder = dash_app.config["assets_folder"]
            assert isinstance(assets_folder, str)
            self.assertEqual(Path(assets_folder), ASSETS_PATH)
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
        return next(key for key in self._application().callback_map if output_id in key)

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

    def _click_details(self, values: CallbackValues, graph: str) -> JsonObject:
        selected = self._post("detail-disease.value", values, f"{graph}.clickData")
        values = dict(values)
        detail = selected["detail-disease"]
        assert isinstance(detail, dict)
        values[("detail-disease", "value")] = detail["value"]
        return self._post("details.children", values, "detail-disease.value")

    def _values(self) -> CallbackValues:
        values: CallbackValues = {
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
        self.assertEqual(_at(status, "update-status", "children"), "Settings applied.")
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
            callback = self._application().callback_map[self._callback_key(output)]
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

    def test_heatmap_selection_keeps_disease_wide_details(self) -> None:
        values = self._values()
        self._apply(values)
        for graph in ("target-heatmap", "drug-heatmap"):
            with self.subTest(graph=graph):
                values[(graph, "clickData")] = _json_value(
                    {
                        "points": [
                            {
                                "customdata": [
                                    "MONDO_RA_TEST",
                                    "group:" + atlas.T_CELL_ID,
                                ]
                            }
                        ]
                    }
                )
                selected = self._post(
                    "detail-disease.value", values, f"{graph}.clickData"
                )
                self.assertEqual(
                    _at(selected, "detail-disease", "value"), "MONDO_RA_TEST"
                )
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
            callbacks,
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
            callbacks,
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
        self.assertIn("Blank: missing CPM/CELLEX", dot_key)
        self.assertEqual(
            _at(dot, "expression-heatmap-key", "style"), {"display": "none"}
        )

    def test_source_records_are_visible_and_paged_without_losing_rows(self) -> None:
        self.assertIn("source-records-page.children", self._application().callback_map)
        context = {"disease_id": "MONDO_RA_TEST", "modality": "all", "stage": "phase3"}
        values: CallbackValues = {
            ("source-page", "value"): 1,
            ("source-context", "data"): _json_value(context),
        }
        records = atlas.filtered_records(self.snapshot, "all", "phase3")
        seen: list[JsonValue] = []
        with patch.object(callbacks, "SOURCE_PAGE_SIZE", 2):
            for page in range(1, math.ceil(len(records) / 2) + 1):
                values[("source-page", "value")] = page
                result = self._post(
                    "source-records-page.children", values, "source-page.value"
                )
                children = _json_array(
                    _at(result, "source-records-page", "children", "props", "children")
                )
                table_rows = _json_array(
                    _at(
                        children[1],
                        "props",
                        "children",
                        "props",
                        "children",
                        1,
                        "props",
                        "children",
                    )
                )
                self.assertLessEqual(len(table_rows), 2)
                seen.extend(
                    _at(row, "props", "children", 1, "props", "children")
                    for row in table_rows
                )
                headers = _json_array(
                    _at(
                        children[1],
                        "props",
                        "children",
                        "props",
                        "children",
                        0,
                        "props",
                        "children",
                        "props",
                        "children",
                    )
                )
                labels_in_table = [
                    _json_array(_at(header, "props", "children"))[0]
                    if isinstance(_at(header, "props", "children"), list)
                    else _at(header, "props", "children")
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
        layout = _response_json(self._client().get("/_dash-layout"))
        controls = next(
            _json_object(node)
            for node in _json_array(_at(layout, "props", "children"))
            if _at(node, "props", "className") == "panel controls"
        )
        self.assertEqual(
            _at(controls, "props", "children", 0, "props", "children", 0),
            "Settings panel",
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
        self.assertIn("Count (scaled size)", target_key)
        self.assertIn('"width": "22px"', target_key)
        self.assertIn("Blank: zero or unavailable", target_key)
        click = _json_value(
            {"points": [{"customdata": ["MONDO_RA_TEST", "group:" + atlas.T_CELL_ID]}]}
        )
        values[("target-heatmap", "clickData")] = click
        details = self._click_details(values, "target-heatmap")
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
        selected = self._values()
        selected[("target-heatmap", "clickData")] = _json_value(
            {"points": [{"customdata": ["MONDO_RA_TEST", "group:CL_B_GROUP"]}]}
        )
        details = self._click_details(selected, "target-heatmap")
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
            "full-reference target median",
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
        click = _json_value(
            {"points": [{"customdata": ["MONDO_RA_TEST", "group:CL_B_GROUP"]}]}
        )
        for view, hidden, graph in (
            ("target", (False, True), "target-heatmap"),
            ("drug", (True, False), "drug-heatmap"),
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
            case = dict(values)
            case[(graph, "clickData")] = click
            details = self._click_details(case, graph)
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

    def test_details_from_both_clicks(self) -> None:
        values = self._values()
        click = _json_value(
            {"points": [{"customdata": ["MONDO_RA_TEST", "group:CL_B_GROUP"]}]}
        )
        for heatmap in ("target-heatmap", "drug-heatmap"):
            case = dict(values)
            case[(heatmap, "clickData")] = click
            details = self._click_details(case, heatmap)
            rendered = str(_at(details, "details", "children"))
            self.assertIn("Relative target expression by cell type", rendered)
            self.assertIn("Drug–target records for rheumatoid arthritis", rendered)


if __name__ == "__main__":
    unittest.main()
