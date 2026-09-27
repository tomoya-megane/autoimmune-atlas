"""遺伝子ページのコールバック。薬剤ページの callback を genetics- の ID で写す。"""

import math
from functools import lru_cache
from typing import TypedDict, cast

import plotly.graph_objects as go  # pyright: ignore[reportMissingTypeStubs] - Plotly に型スタブがない。
from dash import ALL, Dash, Input, Output, State, ctx, html, no_update

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas import genetics as gene_data
from autoimmune_atlas.disease_catalog import disease_catalog, ordered_disease_ids
from autoimmune_atlas.models import (
    CellCatalogEntry,
    GeneticsSnapshot,
    Snapshot,
    SummaryRow,
    TargetRecord,
)
from autoimmune_atlas.ui.callbacks import (
    NumberInput,
    ParameterValue,
    effective_filters,
    effective_number,
    make_toggle,
)
from autoimmune_atlas.ui.components import (
    disease_checklist_sections,
    heatmap_row_controls,
    info_tip,
)
from autoimmune_atlas.ui.config import (
    DEFAULT_SCORE_THRESHOLD,
    GENETICS_PARAMETER_IDS,
    METHOD_LABELS,
    SCORE_FLOOR,
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
from autoimmune_atlas.ui.genetics_layout import genetics_detail_panel

GeneticsParameters = TypedDict(
    "GeneticsParameters",
    {
        "genetics-measure": Measure,
        "genetics-score": NumberInput,
        "genetics-method": str,
        "genetics-threshold": NumberInput,
        "genetics-specificity": NumberInput,
        "genetics-diseases": list[str],
    },
)


class GeneticsSourceContext(TypedDict):
    disease_id: str
    score: float
    threshold: float
    method: str
    specificity: float


def _applied_parameters(values: tuple[ParameterValue, ...]) -> GeneticsParameters:
    """Dash の設定値をコールバック間で共有する形にする。"""
    return cast(
        GeneticsParameters,
        cast(object, dict(zip(GENETICS_PARAMETER_IDS, values, strict=True))),
    )


def effective_score(value: NumberInput) -> tuple[float, str | None]:
    """空欄は初期値、0.1 未満と 1 超と不正値は理由を示して初期値へ戻す。"""
    return effective_number(
        value, DEFAULT_SCORE_THRESHOLD, "score threshold", 1, minimum=SCORE_FLOOR
    )


def register_genetics_callbacks(
    application: Dash, snapshot: Snapshot, genetics: GeneticsSnapshot
) -> None:
    """遺伝子ページのコールバックを登録する。snapshot は merged_snapshot を通したもの。"""
    # 遺伝子は薬剤の標的より桁違いに多く、発現の表と細胞の一覧を毎回作ると 1 回の更新に数秒かかる。
    # 登録時に 1 度だけ作り、各 callback で使い回す。
    metadata = atlas.expression_metadata(snapshot)
    heatmap_catalog = atlas.cell_catalog(snapshot, "mixed")
    all_catalog: list[CellCatalogEntry] = [
        {
            "id": "all",
            "name": "All source cell types",
            "members": [cell["id"] for cell in atlas.cell_catalog(snapshot, "cell")],
        }
    ]
    heatmap_names = {cell["id"]: cell["name"] for cell in heatmap_catalog}
    heatmap_groups = {
        cell["id"] for cell in heatmap_catalog if cell.get("cell_level") == "group"
    }

    @lru_cache(maxsize=1)
    def comparison_rows(
        score: float,
        minimum: float,
        disease_ids: tuple[str, ...],
        cell_ids: tuple[str, ...],
        method: str,
        specificity: float,
    ) -> list[SummaryRow]:
        # ponytail: retain one filter combination per app; enlarge only for concurrent users.
        return gene_data.summarize_genes(
            snapshot,
            genetics,
            score,
            minimum,
            method=method,
            specificity_threshold=specificity,
            level="mixed",
            cell_ids=list(cell_ids),
            disease_ids=list(disease_ids),
            metadata=metadata,
            catalog=heatmap_catalog,
        )

    toggle_cell_group = make_toggle(heatmap_groups)
    application.callback(  # pyright: ignore[reportUnknownMemberType] - Dash の callback メソッドの型が不完全である。
        Output("genetics-expanded-cell-groups", "data"),
        Input({"type": "genetics-cell-toggle", "kind": ALL, "cell": ALL}, "n_clicks"),
        State("genetics-expanded-cell-groups", "data"),
        prevent_initial_call=True,
    )(toggle_cell_group)

    application.callback(  # pyright: ignore[reportUnknownMemberType] - Dash の callback メソッドの型が不完全である。
        Output("genetics-expanded-expression-groups", "data"),
        Input(
            {"type": "genetics-expression-cell-toggle", "kind": ALL, "cell": ALL},
            "n_clicks",
        ),
        State("genetics-expanded-expression-groups", "data"),
        prevent_initial_call=True,
    )(toggle_cell_group)

    @lru_cache(maxsize=1)
    def expression_base(
        disease: str,
        score: float,
        threshold: float,
        method: str,
        specificity: float,
    ) -> go.Figure:
        # ponytail: retain one selected disease; enlarge only for concurrent users.
        records: list[TargetRecord] = [
            {"target_id": gene["target_id"], "target": gene["target"]}
            for gene in gene_data.genes_for_disease(genetics, disease, score)
        ]
        return expression_figure(
            snapshot,
            records,
            metadata,
            threshold=threshold,
            method=method,
            specificity=specificity,
            grouped=True,
        )

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("genetics-expression-heatmap", "figure"),
        Output("genetics-expression-heatmap", "style"),
        Output("genetics-expression-row-controls", "children"),
        Output("genetics-expression-heatmap-key", "style"),
        Output("genetics-expression-dot-key", "children"),
        Output("genetics-expression-dot-key", "style"),
        Input("genetics-source-context", "data"),
        Input("genetics-expanded-expression-groups", "data"),
        Input("genetics-expression-chart-type", "value"),
    )
    def update_expression(
        context: GeneticsSourceContext | None,
        expanded: list[str] | None,
        chart_type: str,
    ) -> tuple[object, object, object, object, object, object]:
        if not context:
            return no_update, no_update, no_update, no_update, no_update, no_update
        base = expression_base(
            context["disease_id"],
            context["score"],
            context["threshold"],
            context["method"],
            context["specificity"],
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
            heatmap_row_controls(
                heatmap_names, cells, expanded or [], "genetics-expression"
            ),
            {"display": "inline-flex" if chart_type == "heatmap" else "none"},
            dot_key,
            {"display": "flex" if chart_type == "dot" else "none"},
        )

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("genetics-applied-parameters", "data"),
        Input("genetics-update-button", "n_clicks"),
        *[State(item, "value") for item in GENETICS_PARAMETER_IDS],
        prevent_initial_call=True,
    )
    def apply_parameters(
        _clicks: int | None, *values: ParameterValue
    ) -> GeneticsParameters:
        return _applied_parameters(values)

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("genetics-update-status", "children"),
        Input("genetics-applied-parameters", "data"),
        *[Input(item, "value") for item in GENETICS_PARAMETER_IDS],
    )
    def parameter_status(applied: GeneticsParameters, *values: ParameterValue) -> str:
        return (
            "Settings applied"
            if _applied_parameters(values) == applied
            else "Changes not applied; click Update"
        )

    families = [
        section
        for group in disease_catalog(snapshot)
        for family in group["families"]
        for section in disease_checklist_sections(family, "genetics-")
    ]
    family_ids = [section["id"] for section in families]

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("genetics-diseases", "value"),
        *[Output(item, "value") for item in family_ids],
        Input("genetics-diseases", "value"),
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
        Output("genetics-specificity", "disabled"), Input("genetics-method", "value")
    )
    def toggle_specificity(method: str | None) -> bool:
        return method != "specificity"

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("genetics-measure", "options"), Input("genetics-chart-type", "value")
    )
    def toggle_measure(chart_type: str | None) -> list[dict[str, str | bool]]:
        disabled = chart_type == "dot"
        return [
            {"label": label, "value": value, "disabled": disabled}
            for label, value in (("Percent", "percent"), ("Count", "count"))
        ]

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("genetics-heatmap", "figure"),
        Output("genetics-heatmap", "style"),
        Output("genetics-matrix-note", "children"),
        Output("genetics-row-controls", "children"),
        Output("genetics-row-controls", "style"),
        Output("genetics-dot-key", "children"),
        Output("genetics-dot-key", "style"),
        Input("genetics-applied-parameters", "data"),
        Input("genetics-expanded-cell-groups", "data"),
        Input("genetics-chart-type", "value"),
    )
    def update_figures(
        applied: GeneticsParameters,
        expanded: list[str] | None,
        chart_type: str | None,
    ) -> tuple[
        go.Figure,
        dict[str, str],
        html.Div,
        list[html.Button | html.Div],
        dict[str, str],
        list[html.Span],
        dict[str, str],
    ]:
        measure = applied["genetics-measure"]
        method = applied["genetics-method"]
        score, score_error = effective_score(applied["genetics-score"])
        minimum, specificity_value, errors = effective_filters(
            applied["genetics-threshold"], applied["genetics-specificity"]
        )
        if score_error:
            errors = [score_error, *errors]
        disease_ids = ordered_disease_ids(snapshot, applied["genetics-diseases"])
        scale_cell_ids = heatmap_cell_ids(
            snapshot,
            tuple(heatmap_groups),
            tuple(heatmap_groups),
            catalog=heatmap_catalog,
        )
        cell_ids = heatmap_cell_ids(
            snapshot, tuple(heatmap_groups), expanded, catalog=heatmap_catalog
        )
        scale_rows = comparison_rows(
            score,
            minimum,
            tuple(disease_ids),
            tuple(scale_cell_ids),
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
        figure = (
            build_dot_figure(
                rows, disease_ids or [], cell_ids or [], "gene", scale_rows=scale_rows
            )
            if chart_type == "dot"
            else build_figure(
                rows,
                disease_ids or [],
                cell_ids or [],
                measure,
                "gene",
                scale_rows=scale_rows,
            )
        )
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
            gene_missing = sum(
                row[measure_fields("gene", "count")[0]] is None
                or row[measure_fields("gene", "percent")[0]] is None
                for row in rows
            )
        else:
            figure.update_traces(  # pyright: ignore[reportUnknownMemberType] - Plotly に型スタブがない。
                colorbar_len=min(240, max(28, 28 * len(cell_ids))),
                selector={"type": "heatmap"},
            )
            gene_missing = sum(
                row[measure_fields("gene", measure)[0]] is None for row in rows
            )
        condition = f"median CPM ≥ {minimum:g}"
        if method == "relative":
            condition += " and CPM ≥ the target-relative median"
        elif method == "specificity":
            condition += f" and CELLEX specificity ≥ {specificity_value:g}"
        display = (
            f"Dot plot: dot area uses a compressed Count scale, color shows Percent, zero and missing entries are empty (missing: genes {gene_missing})"
            if chart_type == "dot"
            else f"Heatmap measure: {measure.capitalize()}, missing entries shown as zero (genes {gene_missing})"
        )
        note = f"Applied · Score threshold: ≥ {score:g} · Expression rule: {METHOD_LABELS[method]} ({condition}) · Cells: all groups and expanded cell types · {display} · {len(rows)} disease–cell combinations"
        status = html.Span(f"{len(rows)} disease–cell combinations")
        error_note = (
            html.Span(" ".join(errors), className="filter-errors", role="alert")
            if errors
            else None
        )
        _, _, values = dot_size_scale(scale_rows, "gene")
        dot_key = [
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
        return (
            figure,
            {"minWidth": f"{max(600, 240 + right_margin + 50 * len(disease_ids))}px"},
            html.Div(
                [
                    status,
                    info_tip("genetics-applied", "applied filters", note),
                    error_note,
                ],
                className="matrix-status",
            ),
            heatmap_row_controls(heatmap_names, cell_ids, expanded or [], "genetics"),
            {"top": f"{top_margin}px"},
            dot_key,
            {"display": "flex" if chart_type == "dot" else "none"},
        )

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("genetics-detail-disease", "options"),
        Output("genetics-detail-disease", "value"),
        Input("genetics-applied-parameters", "data"),
        State("genetics-detail-disease", "value"),
    )
    def update_detail_selector(
        applied: GeneticsParameters,
        current_disease: str | None,
    ) -> tuple[list[dict[str, str]], str | None]:
        disease_ids = ordered_disease_ids(snapshot, applied["genetics-diseases"])
        disease_names = {row["id"]: row["name"] for row in snapshot["diseases"]}
        # 比較図のマスのクリックでは切り替えない。表示中なら今の選択を残し、なければ先頭にする。
        selected = current_disease if current_disease in disease_ids else None
        return (
            [
                {"label": disease_names[item], "value": item}
                for item in disease_ids
                if item in disease_names
            ],
            selected or next(iter(disease_ids), None),
        )

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("genetics-details", "children"),
        Input("genetics-detail-disease", "value"),
        Input("genetics-applied-parameters", "data"),
    )
    def update_details(
        detail_disease: str | None, applied: GeneticsParameters
    ) -> html.Div:
        method = applied["genetics-method"]
        score, _ = effective_score(applied["genetics-score"])
        minimum, specificity_value, _ = effective_filters(
            applied["genetics-threshold"], applied["genetics-specificity"]
        )
        selected = (
            detail_disease if detail_disease in applied["genetics-diseases"] else None
        )
        rows = gene_data.summarize_genes(
            snapshot,
            genetics,
            score,
            minimum,
            method=method,
            specificity_threshold=specificity_value,
            level="all",
            disease_ids=[selected] if selected else [],
            metadata=metadata,
            catalog=all_catalog,
        )
        return genetics_detail_panel(
            rows,
            selected,
            genetics,
            score_threshold=score,
            threshold=minimum,
            method=method,
            specificity=specificity_value,
        )
