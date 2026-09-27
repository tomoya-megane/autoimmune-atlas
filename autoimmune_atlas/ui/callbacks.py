"""Dash アプリのコールバック。"""

import math
from collections.abc import Callable, Sequence
from functools import lru_cache
from typing import TypedDict, cast

import plotly.graph_objects as go  # pyright: ignore[reportMissingTypeStubs] - Plotly に型スタブがない。
from dash import ALL, Dash, Input, Output, State, ctx, html, no_update

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas.disease_catalog import disease_catalog, ordered_disease_ids
from autoimmune_atlas.models import Snapshot, SummaryRow
from autoimmune_atlas.ui.components import (
    disease_checklist_sections,
    heatmap_row_controls,
    info_tip,
)
from autoimmune_atlas.ui.config import (
    DEFAULT_EXPRESSION_THRESHOLD,
    DEFAULT_SPECIFICITY_THRESHOLD,
    METHOD_LABELS,
    PARAMETER_IDS,
    STAGE_LABELS,
)
from autoimmune_atlas.ui.figures import (
    Measure,
    build_dot_figure,
    build_figure,
    disease_label_lines,
    dot_size_scale,
    expression_dot_size_scale,
    expression_figure,
    expression_view,
    heatmap_cell_ids,
    measure_fields,
)
from autoimmune_atlas.ui.layout import detail_panel

type NumberInput = int | float | str | None
type ParameterValue = str | int | float | list[str] | None

AppliedParameters = TypedDict(
    "AppliedParameters",
    {
        "measure": Measure,
        "modality": str,
        "stage": str,
        "method": str,
        "threshold": NumberInput,
        "specificity": NumberInput,
        "diseases": list[str],
        "heatmap-view": str,
    },
)


class SourceContext(TypedDict):
    disease_id: str
    modality: str
    stage: str
    threshold: float
    method: str
    specificity: float


class CellToggleId(TypedDict):
    type: str
    kind: str
    cell: str


def _applied_parameters(values: tuple[ParameterValue, ...]) -> AppliedParameters:
    """Dash の設定値をコールバック間で共有する形にする。"""
    # Dash の各 control が値の型を固定するが、callback デコレーターからはその型を取得できない。
    return cast(
        AppliedParameters,
        cast(object, dict(zip(PARAMETER_IDS, values, strict=True))),
    )


def effective_number(
    value: NumberInput,
    default: float,
    label: str,
    maximum: float | None = None,
    *,
    minimum: float | None = None,
) -> tuple[float, str | None]:
    """空欄には初期値を使い、不正値は理由を示して初期値へ戻す。"""
    lower = minimum if minimum is not None else 0
    if value is None or value == "":
        return default, None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default, f"Invalid {label}; using {default:g}."
    if (
        isinstance(value, bool)
        or not math.isfinite(number)
        or number < lower
        or (maximum is not None and number > maximum)
    ):
        bound = (
            f" between {lower:g} and {maximum:g}"
            if maximum is not None
            else f" at or above {lower:g}"
        )
        return (
            default,
            f"{label[0].upper() + label[1:]} must be finite and{bound}; using {default:g}.",
        )
    return number, None


def effective_filters(
    threshold: NumberInput, specificity: NumberInput
) -> tuple[float, float, list[str]]:
    """実際に集計へ渡す閾値と入力エラーを返す。"""
    minimum, minimum_error = effective_number(
        threshold, DEFAULT_EXPRESSION_THRESHOLD, "minimum CPM"
    )
    specificity_value, specificity_error = effective_number(
        specificity, DEFAULT_SPECIFICITY_THRESHOLD, "CELLEX specificity", 1
    )
    return (
        minimum,
        specificity_value,
        [message for message in (minimum_error, specificity_error) if message],
    )


def effective_threshold(threshold: NumberInput) -> float:
    """後方互換用に、適用される最低 CPM だけを返す。"""
    return effective_filters(threshold, DEFAULT_SPECIFICITY_THRESHOLD)[0]


def visible_rows(
    snapshot: Snapshot,
    modality: str,
    threshold: NumberInput,
    disease_ids: Sequence[str] | None,
    cell_ids: Sequence[str] | None,
    *,
    stage: str = "phase3",
    method: str = "fixed",
    specificity: NumberInput = DEFAULT_SPECIFICITY_THRESHOLD,
    level: str = "group",
) -> list[SummaryRow]:
    """集計結果を現在の表示範囲へ絞る。"""
    minimum, specificity_value, _ = effective_filters(threshold, specificity)
    return atlas.summarize(
        snapshot,
        modality,
        minimum,
        stage=stage,
        method=method,
        specificity_threshold=specificity_value,
        level=level,
        cell_ids=list(cell_ids or ()),
        disease_ids=list(disease_ids or ()),
    )


def make_toggle(
    valid_groups: set[str],
) -> Callable[[list[int | None], list[str] | None], object]:
    """行見出しのクリックで、大分類の展開を切り替える callback の本体を作る。"""

    def toggle_cell_group(
        _clicks: list[int | None], expanded: list[str] | None
    ) -> object:
        # Newly rendered buttons have zero clicks; only user clicks toggle a group.
        triggered = cast(CellToggleId | str | None, ctx.triggered_id)
        inputs_list = cast(list[list[dict[str, object]]], ctx.inputs_list)
        if not isinstance(triggered, dict) or not any(
            item["id"] == triggered and item.get("value") for item in inputs_list[0]
        ):
            return no_update
        cell_id = triggered["cell"]
        if cell_id not in valid_groups:
            return no_update
        expanded_set = set(expanded or ()) & valid_groups
        expanded_set.symmetric_difference_update([cell_id])
        return sorted(expanded_set)

    return toggle_cell_group


def register_callbacks(application: Dash, snapshot: Snapshot) -> None:
    """schema 2 dashboard のコールバックを登録する。"""
    heatmap_catalog = atlas.cell_catalog(snapshot, "mixed")
    heatmap_names = {cell["id"]: cell["name"] for cell in heatmap_catalog}
    heatmap_groups = {
        cell["id"] for cell in heatmap_catalog if cell.get("cell_level") == "group"
    }

    @lru_cache(maxsize=1)
    def comparison_rows(
        modality: str,
        minimum: float,
        disease_ids: tuple[str, ...],
        cell_ids: tuple[str, ...],
        stage: str,
        method: str,
        specificity: float,
    ) -> list[SummaryRow]:
        # ponytail: retain one filter combination per app; enlarge only for concurrent users.
        return visible_rows(
            snapshot,
            modality,
            minimum,
            disease_ids,
            cell_ids,
            stage=stage,
            method=method,
            specificity=specificity,
            level="mixed",
        )

    toggle_cell_group = make_toggle(heatmap_groups)
    application.callback(  # pyright: ignore[reportUnknownMemberType] - Dash の callback メソッドの型が不完全である。
        Output("expanded-cell-groups", "data"),
        Input({"type": "heatmap-cell-toggle", "kind": ALL, "cell": ALL}, "n_clicks"),
        State("expanded-cell-groups", "data"),
        prevent_initial_call=True,
    )(toggle_cell_group)

    application.callback(  # pyright: ignore[reportUnknownMemberType] - Dash の callback メソッドの型が不完全である。
        Output("expanded-expression-groups", "data"),
        Input({"type": "expression-cell-toggle", "kind": ALL, "cell": ALL}, "n_clicks"),
        State("expanded-expression-groups", "data"),
        prevent_initial_call=True,
    )(toggle_cell_group)

    @lru_cache(maxsize=1)
    def expression_base(
        disease: str,
        modality: str,
        stage: str,
        threshold: float,
        method: str,
        specificity: float,
    ) -> go.Figure:
        # ponytail: retain one selected disease; enlarge only for concurrent users.
        records = [
            record
            for record in atlas.filtered_records(snapshot, modality, stage)
            if record["disease_id"] == disease
        ]
        return expression_figure(
            snapshot,
            records,
            threshold=threshold,
            method=method,
            specificity=specificity,
            grouped=True,
        )

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("expression-heatmap", "figure"),
        Output("expression-heatmap", "style"),
        Output("expression-row-controls", "children"),
        Output("expression-heatmap-key", "style"),
        Output("expression-dot-key", "children"),
        Output("expression-dot-key", "style"),
        Input("source-context", "data"),
        Input("expanded-expression-groups", "data"),
        Input("expression-chart-type", "value"),
    )
    def update_expression(
        context: SourceContext | None,
        expanded: list[str] | None,
        chart_type: str,
    ) -> tuple[object, object, object, object, object, object]:
        if not context:
            return no_update, no_update, no_update, no_update, no_update, no_update
        base = expression_base(
            *(
                context[key]
                for key in (
                    "disease_id",
                    "modality",
                    "stage",
                    "threshold",
                    "method",
                    "specificity",
                )
            )
        )
        figure = expression_view(base, heatmap_catalog, expanded or [], chart_type)
        base_heatmap = cast(go.Heatmap, base.data[0])
        x_values = cast(tuple[object, ...] | None, base_heatmap.x)
        cells = heatmap_cell_ids(
            snapshot, tuple(heatmap_groups), expanded, catalog=heatmap_catalog
        )
        _, size_values = expression_dot_size_scale()
        dot_key = [
            html.Span("CELLEX specificity (area)", className="dot-key-title"),
            *[
                html.Span(
                    [
                        html.Span(
                            className="dot-key-circle",
                            style={
                                "width": f"{diameter:g}px",
                                "height": f"{diameter:g}px",
                            },
                        ),
                        f"{value:g}",
                    ],
                    className="dot-key-item",
                )
                for value, diameter in size_values
            ],
            html.Span(
                "Empty cell: missing CPM or CELLEX",
                className="dot-key-item",
            ),
        ]
        return (
            figure,
            {
                "minWidth": f"{max(640, 320 + 28 * len(x_values if x_values is not None else ()))}px",
                "height": f"{figure.layout.height}px",  # pyright: ignore[reportUnknownMemberType] - Plotly に型スタブがない。
            },
            heatmap_row_controls(heatmap_names, cells, expanded or [], "expression"),
            {"display": "inline-flex" if chart_type == "heatmap" else "none"},
            dot_key,
            {"display": "flex" if chart_type == "dot" else "none"},
        )

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("applied-parameters", "data"),
        Input("update-button", "n_clicks"),
        *[State(item, "value") for item in PARAMETER_IDS],
        prevent_initial_call=True,
    )
    def apply_parameters(
        _clicks: int | None, *values: ParameterValue
    ) -> AppliedParameters:
        return _applied_parameters(values)

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("update-status", "children"),
        Input("applied-parameters", "data"),
        *[Input(item, "value") for item in PARAMETER_IDS],
    )
    def parameter_status(applied: AppliedParameters, *values: ParameterValue) -> str:
        return (
            "Settings applied"
            if _applied_parameters(values) == applied
            else "Changes not applied; click Update"
        )

    families = [
        section
        for group in disease_catalog(snapshot)
        for family in group["families"]
        for section in disease_checklist_sections(family)
    ]
    family_ids = [section["id"] for section in families]

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("diseases", "value"),
        *[Output(item, "value") for item in family_ids],
        Input("diseases", "value"),
        *[Input(item, "value") for item in family_ids],
    )
    def sync_disease_selection(
        selected: list[str] | None, *family_values: list[str] | None
    ) -> list[list[str]]:
        selected_set = set(selected or ())
        triggered_id = cast(str | dict[str, str] | None, ctx.triggered_id)
        if isinstance(triggered_id, str) and triggered_id in family_ids:
            index = family_ids.index(triggered_id)
            members = {d["id"] for d in families[index]["diseases"]}
            selected_set = (selected_set - members) | (
                set(family_values[index] or []) & members
            )
        ordered = ordered_disease_ids(snapshot, list(selected_set))
        return [
            ordered,
            *[
                [d["id"] for d in family["diseases"] if d["id"] in ordered]
                for family in families
            ],
        ]

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("specificity", "disabled"), Input("method", "value")
    )
    def toggle_specificity(
        method: str | None,
    ) -> bool:
        return method != "specificity"

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("measure", "options"), Input("chart-type", "value")
    )
    def toggle_measure(chart_type: str | None) -> list[dict[str, str | bool]]:
        disabled = chart_type == "dot"
        return [
            {"label": label, "value": value, "disabled": disabled}
            for label, value in (("Percent", "percent"), ("Count", "count"))
        ]

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("target-heatmap-panel", "hidden"),
        Output("drug-heatmap-panel", "hidden"),
        Input("applied-parameters", "data"),
    )
    def select_heatmap(
        applied: AppliedParameters,
    ) -> tuple[bool, bool]:
        view = applied["heatmap-view"]
        return view == "drug", view != "drug"

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("target-heatmap", "figure"),
        Output("drug-heatmap", "figure"),
        Output("target-heatmap", "style"),
        Output("drug-heatmap", "style"),
        Output("matrix-note", "children"),
        Output("target-row-controls", "children"),
        Output("drug-row-controls", "children"),
        Output("target-row-controls", "style"),
        Output("drug-row-controls", "style"),
        Output("target-dot-key", "children"),
        Output("drug-dot-key", "children"),
        Output("target-dot-key", "style"),
        Output("drug-dot-key", "style"),
        Input("applied-parameters", "data"),
        Input("expanded-cell-groups", "data"),
        Input("chart-type", "value"),
    )
    def update_figures(
        applied: AppliedParameters, expanded: list[str] | None, chart_type: str | None
    ) -> tuple[
        go.Figure,
        go.Figure,
        dict[str, str],
        dict[str, str],
        html.Div,
        list[html.Button | html.Div],
        list[html.Button | html.Div],
        dict[str, str],
        dict[str, str],
        list[html.Span],
        list[html.Span],
        dict[str, str],
        dict[str, str],
    ]:
        measure = applied["measure"]
        modality = applied["modality"]
        stage = applied["stage"]
        method = applied["method"]
        threshold = applied["threshold"]
        specificity = applied["specificity"]
        disease_ids = applied["diseases"]
        cell_ids = heatmap_groups
        minimum, specificity_value, errors = effective_filters(threshold, specificity)
        disease_ids = ordered_disease_ids(snapshot, disease_ids)
        scale_cell_ids = heatmap_cell_ids(
            snapshot,
            tuple(cell_ids),
            tuple(heatmap_groups),
            catalog=heatmap_catalog,
        )
        cell_ids = heatmap_cell_ids(
            snapshot, tuple(cell_ids), expanded, catalog=heatmap_catalog
        )
        scale_rows = comparison_rows(
            modality,
            minimum,
            tuple(disease_ids),
            tuple(scale_cell_ids),
            stage,
            method,
            specificity_value,
        )
        rows = [row for row in scale_rows if row["cell_id"] in cell_ids]
        names = {row["id"]: row["name"] for row in snapshot["diseases"]}
        label_extents = [
            max(map(len, lines)) * 4.5 + (len(lines) - 1) * 12
            for d in disease_ids
            for lines in [disease_label_lines(names.get(d, d))]
        ]
        top_margin = max(130, math.ceil(max(label_extents, default=0)) + 35)
        right_margin = max(90, math.ceil(label_extents[-1]) if label_extents else 0)
        target_figure, drug_figure = (
            build_dot_figure(
                rows,
                disease_ids or [],
                cell_ids or [],
                kind,
                scale_rows=scale_rows,
            )
            if chart_type == "dot"
            else build_figure(
                rows,
                disease_ids or [],
                cell_ids or [],
                measure,
                kind,
                scale_rows=scale_rows,
            )
            for kind in ("target", "drug")
        )
        for figure in (target_figure, drug_figure):
            figure.update_layout(  # pyright: ignore[reportUnknownMemberType] - Plotly に型スタブがない。
                height=top_margin + 70 + 28 * max(1, len(cell_ids)),
                margin={"l": 240, "r": right_margin, "t": top_margin, "b": 70},
            )
            figure.update_xaxes(  # pyright: ignore[reportUnknownMemberType] - Plotly に型スタブがない。
                automargin=False, tickangle=-45
            )
            figure.update_yaxes(  # pyright: ignore[reportUnknownMemberType] - Plotly に型スタブがない。
                title=None, automargin=False, fixedrange=True
            )
            if chart_type == "dot":
                figure.update_traces(  # pyright: ignore[reportUnknownMemberType] - Plotly に型スタブがない。
                    marker_colorbar_len=min(240, max(28, 28 * len(cell_ids))),
                    selector={"type": "scatter", "name": "Values"},
                )
            else:
                figure.update_traces(  # pyright: ignore[reportUnknownMemberType] - Plotly に型スタブがない。
                    colorbar_len=min(240, max(28, 28 * len(cell_ids))),
                    selector={"type": "heatmap"},
                )
        if chart_type == "dot":
            target_missing = sum(
                row[measure_fields("target", "count")[0]] is None
                or row[measure_fields("target", "percent")[0]] is None
                for row in rows
            )
            drug_missing = sum(
                row[measure_fields("drug", "count")[0]] is None
                or row[measure_fields("drug", "percent")[0]] is None
                for row in rows
            )
        else:
            target_missing = sum(
                row[measure_fields("target", measure)[0]] is None for row in rows
            )
            drug_missing = sum(
                row[measure_fields("drug", measure)[0]] is None for row in rows
            )
        condition = f"median CPM ≥ {minimum:g}"
        if method == "relative":
            condition += " and CPM ≥ the target-relative median"
        elif method == "specificity":
            condition += f" and CELLEX specificity ≥ {specificity_value:g}"
        modality_label = (
            "All modalities"
            if modality == "all"
            else next(
                (
                    label
                    for label, value in atlas.DRUG_TYPE_MODALITIES
                    if value == modality
                ),
                modality,
            )
        )
        display = (
            f"Dot plot: dot area uses a compressed Count scale, color shows Percent, zero and missing entries are empty (missing: targets {target_missing}, canonical drugs {drug_missing})"
            if chart_type == "dot"
            else f"Heatmap measure: {measure.capitalize()}, missing entries shown as zero (targets {target_missing}, canonical drugs {drug_missing})"
        )
        note = f"Applied · Clinical stage: {STAGE_LABELS[stage]} · Drug modality: {modality_label} · Expression rule: {METHOD_LABELS[method]} ({condition}) · Cells: all groups and expanded cell types · {display} · {len(rows)} disease–cell combinations"
        status = html.Span(f"{len(rows)} disease–cell combinations")
        error_note = (
            html.Span(" ".join(errors), className="filter-errors", role="alert")
            if errors
            else None
        )
        graph_style = {
            "minWidth": f"{max(600, 240 + right_margin + 50 * len(disease_ids))}px"
        }

        def dot_key(kind: str) -> list[html.Span]:
            _, _, values = dot_size_scale(
                scale_rows, "drug" if kind == "drug" else "target"
            )
            return [
                html.Span("Count (area)", className="dot-key-title"),
                *[
                    html.Span(
                        [
                            html.Span(
                                className="dot-key-circle",
                                style={
                                    "width": f"{diameter:g}px",
                                    "height": f"{diameter:g}px",
                                },
                            ),
                            f"{value:g}",
                        ],
                        className="dot-key-item",
                    )
                    for value, diameter in values
                ],
                html.Span("Empty cell: zero or missing", className="dot-key-item"),
            ]

        key_style = {"display": "flex" if chart_type == "dot" else "none"}
        return (
            target_figure,
            drug_figure,
            graph_style,
            graph_style,
            html.Div(
                [status, info_tip("applied", "applied filters", note), error_note],
                className="matrix-status",
            ),
            heatmap_row_controls(heatmap_names, cell_ids, expanded or [], "target"),
            heatmap_row_controls(heatmap_names, cell_ids, expanded or [], "drug"),
            {"top": f"{top_margin}px"},
            {"top": f"{top_margin}px"},
            dot_key("target"),
            dot_key("drug"),
            key_style,
            key_style,
        )

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("detail-disease", "options"),
        Output("detail-disease", "value"),
        Input("applied-parameters", "data"),
        State("detail-disease", "value"),
    )
    def update_detail_selector(
        applied: AppliedParameters, current_disease: str | None
    ) -> tuple[list[dict[str, str]], str | None]:
        disease_ids = ordered_disease_ids(snapshot, applied["diseases"])
        disease_names = {row["id"]: row["name"] for row in snapshot["diseases"]}
        return (
            [
                {"label": disease_names[item], "value": item}
                for item in disease_ids
                if item in disease_names
            ],
            current_disease
            if current_disease in disease_ids
            else next(iter(disease_ids), None),
        )

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("details", "children"),
        Input("detail-disease", "value"),
        Input("applied-parameters", "data"),
    )
    def update_details(
        detail_disease: str | None, applied: AppliedParameters
    ) -> html.Div:
        modality = applied["modality"]
        stage = applied["stage"]
        method = applied["method"]
        threshold = applied["threshold"]
        specificity = applied["specificity"]
        disease_ids = applied["diseases"]
        minimum, specificity_value, _ = effective_filters(threshold, specificity)
        selected = detail_disease if detail_disease in disease_ids else None
        rows = atlas.summarize(
            snapshot,
            modality,
            minimum,
            stage=stage,
            method=method,
            specificity_threshold=specificity_value,
            level="all",
            disease_ids=[selected] if selected else [],
        )
        return detail_panel(
            rows,
            selected,
            snapshot,
            modality=modality,
            stage=stage,
            threshold=minimum,
            method=method,
            specificity=specificity_value,
        )
