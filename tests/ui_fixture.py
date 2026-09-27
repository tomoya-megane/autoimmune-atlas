"""UI のテストが共有する型、スナップショット、補助関数。"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Protocol, TypedDict, Unpack, cast

from flask.testing import FlaskClient

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas.models import (
    DrugRecord,
    RootIdentity,
    Snapshot,
    SummaryRow,
)

__all__ = [
    "ASSETS_PATH",
    "JsonScalar",
    "JsonValue",
    "SummaryOverrides",
    "FigureAxis",
    "FigureLayout",
    "FigureData",
    "ColorBar",
    "HeatmapTrace",
    "DotMarker",
    "ScatterTrace",
    "ChildrenComponent",
    "ClassComponent",
    "CallbackDependency",
    "CallbackInput",
    "CallbackDefinition",
    "DashServer",
    "DashApplication",
    "_figure",
    "_heatmap",
    "_scatter",
    "_component",
    "fixture",
    "summary",
]

ASSETS_PATH = Path(__file__).resolve().parents[1] / "assets"


type JsonScalar = bool | float | int | str | None


type JsonValue = JsonScalar | list[JsonValue] | dict[str, JsonValue]


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


def _heatmap(figure: object, name: str | None = None) -> HeatmapTrace:
    traces = _figure(figure).data
    if name is not None:
        traces = [
            trace
            for trace in traces
            if hasattr(trace, "z") and cast(HeatmapTrace, trace).name == name
        ]
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


def _component(component: object) -> ChildrenComponent:
    assert hasattr(component, "children")
    return cast(ChildrenComponent, component)


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
        "schema": 3,
        "root": root,
        "data_version": {"year": "26", "month": "6", "iteration": None},
        "retrieved_at": "2026-09-22T12:00:00Z",
        "source": "https://api.platform.opentargets.org/",
        "diseases": [
            {"id": "MONDO_RA_TEST", "name": "rheumatoid arthritis", "status": "ready"}
        ],
        "records": records,
        "datasources": [],
        "associations": {"MONDO_RA_TEST": []},
        "cells": {
            "CL_B_ONE": {
                "name": "memory B cell",
                "parent_id": "CL_B_GROUP",
                "parent": "B cell",
                "ancestor_ids": ["CL_B_GROUP"],
            },
            "CL_B_TWO": {
                "name": "naive B cell",
                "parent_id": "CL_B_GROUP",
                "parent": "B cell",
                "ancestor_ids": ["CL_B_GROUP"],
            },
            "CL_T_ONE": {
                "name": "CD8-positive T cell",
                "parent_id": "CL_T_PARENT",
                "parent": "T lymphocyte",
                "ancestor_ids": [atlas.T_CELL_ID],
            },
        },
        "expression": {
            "ENSG_TARGET_1": [
                {"cell_id": "CL_B_ONE", "median": 2.0, "specificity_score": 0.8},
                {"cell_id": "CL_B_TWO", "median": 0.1, "specificity_score": 0.2},
                {"cell_id": "CL_T_ONE", "median": 0.2, "specificity_score": 0.1},
            ],
            "ENSG_TARGET_2": [
                {"cell_id": "CL_B_ONE", "median": None, "specificity_score": None},
                {"cell_id": "CL_B_TWO", "median": 0.6, "specificity_score": 0.75},
                {"cell_id": "CL_T_ONE", "median": 1.0, "specificity_score": 0.7},
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
