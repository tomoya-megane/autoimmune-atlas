"""Dash アプリのコールバック。"""

import math
from functools import lru_cache

from dash import ALL, Dash, Input, Output, State, ctx, html, no_update

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas.disease_catalog import disease_catalog, ordered_disease_ids
from autoimmune_atlas.ui.components import (
    _disease_checklist_sections,
    _evidence_table,
    heatmap_row_controls,
    info_tip,
)
from autoimmune_atlas.ui.config import (
    DEFAULT_EXPRESSION_THRESHOLD,
    DEFAULT_SPECIFICITY_THRESHOLD,
    METHOD_LABELS,
    PARAMETER_IDS,
    SOURCE_PAGE_SIZE,
    STAGE_LABELS,
)
from autoimmune_atlas.ui.figures import (
    _measure_fields,
    build_figure,
    disease_label_lines,
    expression_figure,
    expression_view,
    heatmap_cell_ids,
)
from autoimmune_atlas.ui.layout import detail_panel


def _effective_number(value, default, label, maximum=None):
    """空欄には初期値を使い、不正値は理由を示して初期値へ戻す。"""
    if value is None or value == "":
        return default, None
    try:
        number = float(value)
    except TypeError, ValueError:
        return default, f"Invalid {label}; using {default:g}."
    if (
        isinstance(value, bool)
        or not math.isfinite(number)
        or number < 0
        or (maximum is not None and number > maximum)
    ):
        bound = (
            f" between 0 and {maximum:g}" if maximum is not None else " at or above 0"
        )
        return (
            default,
            f"{label.capitalize()} must be finite and{bound}; using {default:g}.",
        )
    return number, None


def effective_filters(threshold, specificity):
    """実際に集計へ渡す閾値と入力エラーを返す。"""
    minimum, minimum_error = _effective_number(
        threshold, DEFAULT_EXPRESSION_THRESHOLD, "minimum CPM"
    )
    specificity_value, specificity_error = _effective_number(
        specificity, DEFAULT_SPECIFICITY_THRESHOLD, "specificity threshold", 1
    )
    return (
        minimum,
        specificity_value,
        [message for message in (minimum_error, specificity_error) if message],
    )


def effective_threshold(threshold) -> float:
    """後方互換用に、適用される最低 CPM だけを返す。"""
    return effective_filters(threshold, DEFAULT_SPECIFICITY_THRESHOLD)[0]


def visible_rows(
    snapshot,
    modality,
    threshold,
    disease_ids,
    cell_ids,
    *,
    stage="phase3",
    method="fixed",
    specificity=DEFAULT_SPECIFICITY_THRESHOLD,
    level="group",
):
    """backend の集計結果を現在の表示範囲へ絞る。"""
    minimum, specificity_value, _ = effective_filters(threshold, specificity)
    return atlas.summarize(
        snapshot,
        modality,
        minimum,
        stage=stage,
        method=method,
        specificity_threshold=specificity_value,
        level=level,
        cell_ids=cell_ids or [],
        disease_ids=disease_ids or [],
    )


def resolve_disease_selection(triggered_id, click_data, disease_id, visible_ids):
    """表示中の疾患に限って、クリックまたは選択欄を受け付ける。"""
    selection = disease_id
    if (
        triggered_id in {"target-heatmap", "drug-heatmap", "heatmap"}
        and click_data
        and click_data.get("points")
    ):
        custom = click_data["points"][0].get("customdata")
        if isinstance(custom, (list, tuple)) and custom:
            selection = custom[0]
    return selection if selection in visible_ids else None


def register_callbacks(application: Dash, snapshot: dict) -> None:
    """schema 2 dashboard のコールバックを登録する。"""
    heatmap_catalog = atlas.cell_catalog(snapshot, "mixed")
    heatmap_names = {cell["id"]: cell["name"] for cell in heatmap_catalog}
    heatmap_groups = {
        cell["id"] for cell in heatmap_catalog if cell["cell_level"] == "group"
    }

    @lru_cache(maxsize=1)
    def comparison_rows(
        modality, minimum, disease_ids, cell_ids, stage, method, specificity
    ):
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

    @application.callback(
        Output("expanded-cell-groups", "data"),
        Input({"type": "heatmap-cell-toggle", "kind": ALL, "cell": ALL}, "n_clicks"),
        State("expanded-cell-groups", "data"),
        prevent_initial_call=True,
    )
    def toggle_cell_group(_clicks, expanded):
        # Newly rendered buttons have zero clicks; only user clicks toggle a group.
        triggered = ctx.triggered_id
        if not isinstance(triggered, dict) or not any(
            item["id"] == triggered and item.get("value") for item in ctx.inputs_list[0]
        ):
            return no_update
        cell_id = triggered["cell"]
        valid = heatmap_groups
        if cell_id not in valid:
            return no_update
        expanded = set(expanded or []) & valid
        expanded.symmetric_difference_update([cell_id])
        return sorted(expanded)

    application.callback(
        Output("expanded-expression-groups", "data"),
        Input({"type": "expression-cell-toggle", "kind": ALL, "cell": ALL}, "n_clicks"),
        State("expanded-expression-groups", "data"),
        prevent_initial_call=True,
    )(toggle_cell_group)

    @lru_cache(maxsize=1)
    def expression_base(disease, modality, stage, threshold, method, specificity):
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

    @application.callback(
        Output("expression-heatmap", "figure"),
        Output("expression-heatmap", "style"),
        Output("expression-row-controls", "children"),
        Input("source-context", "data"),
        Input("expanded-expression-groups", "data"),
    )
    def update_expression(context, expanded):
        if not context:
            return no_update, no_update, no_update
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
        figure = expression_view(base, heatmap_catalog, expanded or [])
        cells = heatmap_cell_ids(
            snapshot, heatmap_groups, expanded, catalog=heatmap_catalog
        )
        return (
            figure,
            {
                "minWidth": f"{max(640, 320 + 28 * len(figure.data[0].x))}px",
                "height": f"{figure.layout.height}px",
            },
            heatmap_row_controls(heatmap_names, cells, expanded or [], "expression"),
        )

    @application.callback(
        Output("applied-parameters", "data"),
        Input("update-button", "n_clicks"),
        *[State(item, "value") for item in PARAMETER_IDS],
        prevent_initial_call=True,
    )
    def apply_parameters(_clicks, *values):
        return dict(zip(PARAMETER_IDS, values))

    @application.callback(
        Output("update-status", "children"),
        Input("applied-parameters", "data"),
        *[Input(item, "value") for item in PARAMETER_IDS],
    )
    def parameter_status(applied, *values):
        return (
            "Settings applied."
            if dict(zip(PARAMETER_IDS, values)) == applied
            else "Changes not applied. Click Update."
        )

    families = [
        section
        for group in disease_catalog(snapshot)
        for family in group["families"]
        for section in _disease_checklist_sections(family)
    ]
    family_ids = [section["id"] for section in families]

    @application.callback(
        Output("diseases", "value"),
        *[Output(item, "value") for item in family_ids],
        Input("diseases", "value"),
        *[Input(item, "value") for item in family_ids],
    )
    def sync_disease_selection(selected, *family_values):
        selected = set(selected or [])
        if ctx.triggered_id in family_ids:
            index = family_ids.index(ctx.triggered_id)
            members = {d["id"] for d in families[index]["diseases"]}
            selected = (selected - members) | (
                set(family_values[index] or []) & members
            )
        ordered = ordered_disease_ids(snapshot, selected)
        return [
            ordered,
            *[
                [d["id"] for d in family["diseases"] if d["id"] in ordered]
                for family in families
            ],
        ]

    @application.callback(Output("specificity", "disabled"), Input("method", "value"))
    def toggle_specificity(method):
        return method != "specificity"

    @application.callback(
        Output("target-heatmap-panel", "hidden"),
        Output("drug-heatmap-panel", "hidden"),
        Input("applied-parameters", "data"),
    )
    def select_heatmap(applied):
        view = applied["heatmap-view"]
        return view == "drug", view != "drug"

    @application.callback(
        Output("target-heatmap", "figure"),
        Output("drug-heatmap", "figure"),
        Output("target-heatmap", "style"),
        Output("drug-heatmap", "style"),
        Output("matrix-note", "children"),
        Output("target-row-controls", "children"),
        Output("drug-row-controls", "children"),
        Output("target-row-controls", "style"),
        Output("drug-row-controls", "style"),
        Input("applied-parameters", "data"),
        Input("expanded-cell-groups", "data"),
    )
    def update_figures(applied, expanded):
        (
            measure,
            modality,
            stage,
            method,
            threshold,
            specificity,
            disease_ids,
        ) = (applied[key] for key in PARAMETER_IDS[:-1])
        cell_ids = heatmap_groups
        minimum, specificity_value, errors = effective_filters(threshold, specificity)
        disease_ids = ordered_disease_ids(snapshot, disease_ids)
        scale_cell_ids = heatmap_cell_ids(
            snapshot, cell_ids, heatmap_groups, catalog=heatmap_catalog
        )
        cell_ids = heatmap_cell_ids(
            snapshot, cell_ids, expanded, catalog=heatmap_catalog
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
            build_figure(
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
            figure.update_layout(
                height=top_margin + 70 + 28 * max(1, len(cell_ids)),
                margin={"l": 240, "r": right_margin, "t": top_margin, "b": 70},
            )
            figure.update_xaxes(automargin=False, tickangle=-45)
            figure.update_yaxes(title=None, automargin=False, fixedrange=True)
            figure.update_traces(
                colorbar_len=min(240, max(28, 28 * len(cell_ids))),
                selector={"type": "heatmap"},
            )
        target_missing = sum(
            row[_measure_fields("target", measure)[0]] is None for row in rows
        )
        drug_missing = sum(
            row[_measure_fields("drug", measure)[0]] is None for row in rows
        )
        condition = f"median CPM ≥ {minimum:g}"
        if method == "relative":
            condition += " and CPM ≥ the full-reference target median"
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
        note = f"Applied: Clinical stage: {STAGE_LABELS[stage]}; Drug modality: {modality_label}; Expression rule: {METHOD_LABELS[method]} ({condition}); Cells: all groups and expanded source cells; Measure: {measure.capitalize()}. {len(rows)} disease–cell combinations. Unavailable entries shown as zero: targets {target_missing}, drugs {drug_missing}."
        status = html.Span(f"{len(rows)} disease–cell combinations")
        error_note = (
            html.Span(" ".join(errors), className="filter-errors", role="alert")
            if errors
            else None
        )
        graph_style = {
            "minWidth": f"{max(600, 240 + right_margin + 50 * len(disease_ids))}px"
        }
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
        )

    @application.callback(
        Output("detail-disease", "options"),
        Output("detail-disease", "value"),
        Input("applied-parameters", "data"),
        Input("target-heatmap", "clickData"),
        Input("drug-heatmap", "clickData"),
        State("detail-disease", "value"),
    )
    def update_detail_selector(applied, target_click, drug_click, current_disease):
        disease_ids = ordered_disease_ids(snapshot, applied["diseases"])
        disease_names = {row["id"]: row["name"] for row in snapshot["diseases"]}
        triggered = ctx.triggered_id
        click_data = (
            target_click
            if triggered == "target-heatmap"
            else drug_click
            if triggered == "drug-heatmap"
            else None
        )
        selected = resolve_disease_selection(
            triggered, click_data, current_disease, disease_ids
        )
        return (
            [
                {"label": disease_names[item], "value": item}
                for item in disease_ids
                if item in disease_names
            ],
            selected or next(iter(disease_ids), None),
        )

    @application.callback(
        Output("details", "children"),
        Input("detail-disease", "value"),
        Input("applied-parameters", "data"),
    )
    def update_details(detail_disease, applied):
        (
            _measure,
            modality,
            stage,
            method,
            threshold,
            specificity,
            disease_ids,
        ) = (applied[key] for key in PARAMETER_IDS[:-1])
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

    @application.callback(
        Output("source-records-page", "children"),
        Input("source-page", "value"),
        Input("source-context", "data"),
    )
    def show_source_records(page, context):
        if not context:
            return None
        records = [
            record
            for record in atlas.filtered_records(
                snapshot, context["modality"], context["stage"]
            )
            if record["disease_id"] == context["disease_id"]
        ]
        pages = max(1, math.ceil(len(records) / SOURCE_PAGE_SIZE))
        page = (
            min(max(page, 1), pages)
            if isinstance(page, int) and not isinstance(page, bool)
            else 1
        )
        start = (page - 1) * SOURCE_PAGE_SIZE
        end = min(start + SOURCE_PAGE_SIZE, len(records))
        return html.Div(
            [
                html.P(
                    f"Rows {start + 1 if records else 0}–{end} of {len(records)}",
                    role="status",
                    className="matrix-note",
                ),
                _evidence_table(records[start:end]),
            ]
        )
