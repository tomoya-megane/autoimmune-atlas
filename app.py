"""自己免疫疾患の薬剤標的と細胞型を閲覧する Dash アプリ。"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import plotly.graph_objects as go
from dash import Dash, Input, Output, State, ctx, dcc, html

import atlas

BASE_DIR = Path(__file__).parent
SNAPSHOT_PATH = BASE_DIR / "data" / "snapshot.json"
FILTER_IDS = {"measure", "modality", "threshold", "diseases", "cells"}
DEFAULT_EXPRESSION_THRESHOLD = 0.5


def load_snapshot(path: Path) -> dict | None:
    """スナップショットを読み、UI が前提とする構造を確認する。

    Parameters
    ----------
    path : Path
        JSON ファイルのパス。

    Returns
    -------
    dict | None
        読み込んだスナップショット。ファイルが無い場合は None。

    Raises
    ------
    ValueError
        schema または必須の要素が不正な場合。

    """
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as handle:
        snapshot = json.load(handle)
    if not isinstance(snapshot, dict) or snapshot.get("schema") != 1:
        raise ValueError("snapshot.json must use schema 1")
    required = {"root", "retrieved_at", "source", "diseases", "records", "expression"}
    missing = required - snapshot.keys()
    if missing:
        raise ValueError("Missing required fields in snapshot.json: " + ", ".join(sorted(missing)))
    if not isinstance(snapshot["diseases"], list) or not isinstance(snapshot["records"], list):
        raise ValueError("diseases and records must be arrays")
    if not isinstance(snapshot["expression"], dict):
        raise ValueError("expression must be an object")
    return snapshot


def cell_catalog(snapshot: dict) -> list[tuple[str, str]]:
    """発現データに現れる細胞型を名前順で返す。"""
    cells = {
        row["cell_id"]: row["cell"]
        for rows in snapshot.get("expression", {}).values()
        for row in rows
    }
    return sorted(cells.items(), key=lambda item: item[1].casefold())


def choose_defaults(
    items: list[tuple[str, str]],
    terms: tuple[str, ...],
    limit: int,
    preferred_ids: set[str] | None = None,
) -> list[str]:
    """よく使う名前を優先し、実在する ID だけを初期選択する。"""
    selected = []
    names = {item_id: name.casefold() for item_id, name in items}
    for term in terms:
        exact = next((item_id for item_id, name in names.items() if name == term.casefold()), None)
        if exact is not None and exact not in selected:
            selected.append(exact)
    for term in terms:
        for item_id, name in items:
            if term.casefold() in name.casefold() and item_id not in selected:
                selected.append(item_id)
                break
    preferred = preferred_ids or set()
    for item_id, _ in sorted(items, key=lambda item: item[0] not in preferred):
        if len(selected) >= limit:
            break
        if item_id not in selected:
            selected.append(item_id)
    return selected[:limit]


def format_data_version(value: object) -> str:
    """Open Targets の版を画面と CSV で使える短い文字列にする。"""
    if not isinstance(value, dict):
        return str(value or "")
    parts = [value.get("year"), value.get("month"), value.get("iteration")]
    return ".".join(str(part) for part in parts if part is not None)


def visible_rows(
    snapshot: dict,
    modality: str,
    threshold: float | None,
    disease_ids: list[str] | None,
    cell_ids: list[str] | None,
) -> list[dict]:
    """集計結果を現在の疾患と細胞型へ絞る。"""
    selected_diseases = set(disease_ids or [])
    selected_cells = set(cell_ids or [])
    rows = atlas.summarize(snapshot, modality, effective_threshold(threshold))
    return [
        row
        for row in rows
        if row["disease_id"] in selected_diseases and row["cell_id"] in selected_cells
    ]


def effective_threshold(threshold: float | None) -> float:
    """入力が空なら初期値を使い、負の値は 0 に丸める。"""
    return DEFAULT_EXPRESSION_THRESHOLD if threshold is None else max(0.0, threshold)


def _is_lower_bound(row: dict, measure: str) -> bool:
    """割合は既知標的内の未判定、標的数は標的不明の薬剤も考慮する。"""
    value = row["percent"] if measure == "percent" else row["count"]
    return value is not None and row["status"] == "partial" and (
        row["unknown"] > 0 or (measure != "percent" and row["unmapped_drugs"] > 0)
    )


def _display_value(row: dict, measure: str) -> str:
    value = row["percent"] if measure == "percent" else row["count"]
    if value is None:
        return ""
    rendered = f"{value:.1f}%" if measure == "percent" else str(value)
    return "≥" + rendered if _is_lower_bound(row, measure) else rendered


def _hover_text(row: dict, measure: str) -> str:
    value = _display_value(row, measure) or "Unknown / missing"
    return "<br>".join(
        (
            f"<b>{row['disease']} × {row['cell']}</b>",
            f"Value: {value}",
            "≥ indicates a lower bound; unresolved targets may increase this value" if _is_lower_bound(row, measure) else "All targets within the displayed scope are assessed" if value != "Unknown / missing" else "Insufficient data to determine a value",
            f"Known targets: {row['denominator']}",
            f"Unassessed targets: {row['unknown']}",
            f"Drugs without mapped targets: {row['unmapped_drugs']}",
            f"Status: {row['status']}",
        ),
    )


def build_figure(
    rows: list[dict],
    disease_ids: list[str],
    cell_ids: list[str],
    measure: str,
) -> go.Figure:
    """0、下限値、欠測を区別した heatmap を作る。"""
    lookup = {(row["disease_id"], row["cell_id"]): row for row in rows}
    disease_names = {
        row["disease_id"]: row["disease"] for row in rows
    }
    cell_names = {row["cell_id"]: row["cell"] for row in rows}
    z, text_grid, hover_grid, custom_grid = [], [], [], []
    missing_x, missing_y, missing_custom = [], [], []
    for cell_id in cell_ids:
        z_row, text_row, hover_row, custom_row = [], [], [], []
        for disease_id in disease_ids:
            row = lookup.get((disease_id, cell_id))
            value = None if row is None else row["percent"] if measure == "percent" else row["count"]
            z_row.append(value)
            text_row.append("" if row is None else _display_value(row, measure))
            hover_row.append("No data" if row is None else _hover_text(row, measure))
            custom_row.append([disease_id, cell_id])
            if value is None:
                missing_x.append(disease_names.get(disease_id, disease_id))
                missing_y.append(cell_names.get(cell_id, cell_id))
                missing_custom.append([disease_id, cell_id])
        z.append(z_row)
        text_grid.append(text_row)
        hover_grid.append(hover_row)
        custom_grid.append(custom_row)
    values = [value for z_row in z for value in z_row if value is not None]
    zmax = 100 if measure == "percent" else max(values, default=1) or 1
    figure = go.Figure(
        go.Heatmap(
            x=[disease_names.get(item, item) for item in disease_ids],
            y=[cell_names.get(item, item) for item in cell_ids],
            z=z,
            zmin=0,
            zmax=zmax,
            colorscale=((0.0, "#f5f7fa"), (0.35, "#9fc4e0"), (1.0, "#1769aa")),
            colorbar={"title": "Known targets (%)" if measure == "percent" else "Target count"},
            text=text_grid,
            texttemplate="%{text}",
            hovertext=hover_grid,
            hovertemplate="%{hovertext}<extra></extra>",
            customdata=custom_grid,
            xgap=2,
            ygap=2,
        ),
    )
    if missing_x:
        figure.add_trace(
            go.Scatter(
                x=missing_x,
                y=missing_y,
                customdata=missing_custom,
                mode="markers",
                marker={"symbol": "x", "size": 10, "color": "#8a99a6"},
                name="Unknown / missing",
                hovertemplate="Unknown / missing<extra></extra>",
            ),
        )
    if not disease_ids or not cell_ids:
        figure.add_annotation(text="Select diseases and cell types", showarrow=False)
    figure.update_layout(
        template="plotly_white",
        font={"family": "Arial, sans-serif", "size": 12, "color": "#263238"},
        width=max(760, min(2800, 260 + 125 * len(disease_ids))),
        height=max(440, min(1400, 190 + 28 * len(cell_ids))),
        margin={"l": 180, "r": 40, "t": 35, "b": 150},
        legend={"orientation": "h", "y": 1.08, "x": 0},
        hoverlabel={"align": "left"},
    )
    figure.update_xaxes(tickangle=-35, side="bottom", title="Disease")
    figure.update_yaxes(autorange="reversed", title="Cell type")
    return figure


def resolve_selection(
    triggered_id: str | None,
    click_data: dict | None,
    disease_id: str | None,
    cell_id: str | None,
    rows: list[dict],
) -> tuple[str, str] | None:
    """表示中の組合せに限って詳細の選択を受け付ける。"""
    visible = {(row["disease_id"], row["cell_id"]) for row in rows}
    selection = None
    if triggered_id == "heatmap" and click_data and click_data.get("points"):
        custom = click_data["points"][0].get("customdata")
        if isinstance(custom, (list, tuple)) and len(custom) >= 2:
            selection = custom[0], custom[1]
    elif disease_id and cell_id:
        selection = disease_id, cell_id
    return selection if selection in visible else None


def detail_panel(rows: list[dict], selection: tuple[str, str] | None) -> html.Div:
    """選択した疾患と細胞型の根拠を native table で表示する。"""
    if selection is None:
        return html.Div("Select a disease and cell type using the heatmap or the selectors above.", className="empty-note")
    row = next(item for item in rows if (item["disease_id"], item["cell_id"]) == selection)
    state = {
        "complete": "Fully assessed",
        "partial": "Unassessed targets or drugs without mapped targets",
        "unavailable": "Unavailable",
    }[row["status"]]
    summary = html.Div(
        [
            html.H3(f"{row['disease']} × {row['cell']}"),
            html.P(
                f"Status: {state} | Target count: {_display_value(row, 'count') or 'Unknown'} | "
                f"Share of known targets in this disease: {_display_value(row, 'percent') or 'Unknown'} | "
                f"Known targets: {row['denominator']} | Unassessed targets: {row['unknown']} | "
                f"Drugs without mapped targets: {row['unmapped_drugs']}",
                className="detail-summary",
            ),
        ],
    )
    if not row["records"]:
        message = (
            "Available data are insufficient to determine an expressed-target count for this selection."
            if row["count"] is None
            else "No mapped drug targets have measured expression above the threshold for this selection."
        )
        return html.Div([summary, html.P(message, className="empty-note")])
    body = []
    for record in row["records"]:
        source = (
            html.A("Source", href=record["evidence"], target="_blank", rel="noreferrer")
            if record.get("evidence")
            else "—"
        )
        body.append(
            html.Tr(
                [
                    html.Td(record["drug"]),
                    html.Td(record["drug_type"]),
                    html.Td(record["stage"]),
                    html.Td(record["target"]),
                    html.Td(record["mechanism"] or "—"),
                    html.Td(record.get("note") or "—"),
                    html.Td(source),
                ],
            ),
        )
    return html.Div(
        [
            summary,
            html.Div(
                html.Table(
                    [
                        html.Thead(html.Tr([html.Th(item) for item in ("Drug", "Modality", "Stage", "Target", "Mechanism of action", "Expression evidence", "Source")])),
                        html.Tbody(body),
                    ],
                ),
                className="table-scroll",
            ),
        ],
    )


def export_rows(
    rows: list[dict],
    measure: str,
    *,
    modality_filter: str = "",
    expression_threshold: float | None = None,
    snapshot: dict | None = None,
) -> list[dict]:
    """matrix と根拠を 1 行形式へ展開する。"""
    output = []
    detail_fields = ("drug_id", "drug", "modality", "drug_type", "stage", "target_id", "target", "mechanism", "evidence", "note")
    snapshot = snapshot or {}
    for row in rows:
        details = row["records"] or [{}]
        for detail in details:
            output.append(
                {
                    "disease_id": row["disease_id"],
                    "disease": row["disease"],
                    "cell_id": row["cell_id"],
                    "cell": row["cell"],
                    "mode": "expression",
                    "measure": measure,
                    "modality_filter": modality_filter,
                    "expression_threshold": expression_threshold,
                    "data_version": format_data_version(snapshot.get("data_version")),
                    "retrieved_at": snapshot.get("retrieved_at", ""),
                    "display_value": _display_value(row, measure),
                    "count": row["count"],
                    "percent": row["percent"],
                    "status": row["status"],
                    "denominator": row["denominator"],
                    "unknown_targets": row["unknown"],
                    "unmapped_drugs": row["unmapped_drugs"],
                    **{field: detail.get(field, "") for field in detail_fields},
                },
            )
    return output


def unavailable_layout(error: Exception | None = None) -> html.Main:
    """データ未取得または読み込み失敗を説明する。"""
    if error is None:
        title = "No data loaded yet"
        message = "Create data/snapshot.json to display the heatmap."
        detail = "pixi run --as-is python fetch_data.py"
    else:
        title = "Unable to load data"
        message = "Correct the input files and restart the app."
        detail = f"{type(error).__name__}: {error}"
    return html.Main(
        [
            html.Section(
                [
                    html.P("AUTOIMMUNE DRUG–CELL ATLAS", className="eyebrow"),
                    html.H1(title),
                    html.P(message),
                    html.Code(detail),
                ],
                className="unavailable-card",
            ),
        ],
        className="shell unavailable",
    )


def dashboard_layout(snapshot: dict) -> html.Main:
    """読み込んだデータから dashboard の初期画面を作る。"""
    diseases = [(row["id"], row["name"]) for row in snapshot["diseases"]]
    cells = cell_catalog(snapshot)
    default_diseases = choose_defaults(
        diseases,
        (
            "rheumatoid arthritis",
            "systemic lupus erythematosus",
            "multiple sclerosis",
            "systemic sclerosis",
            "Sjogren syndrome",
            "myasthenia gravis",
            "psoriatic arthritis",
            "type 1 diabetes mellitus",
        ),
        8,
        {row["disease_id"] for row in snapshot["records"]},
    )
    default_cells = choose_defaults(
        cells,
        ("b cell", "t cell", "monocyte", "macrophage", "dendritic", "natural killer", "neutrophil", "fibroblast", "endothelial", "epithelial", "plasma cell"),
        16,
    )
    ready = sum(row["status"] == "ready" for row in snapshot["diseases"])
    unavailable = len(snapshot["diseases"]) - ready
    unclassified = sum(row.get("unclassified_stages", 0) for row in snapshot["diseases"])
    data_version = format_data_version(snapshot.get("data_version"))
    source = "https://platform.opentargets.org/"
    return html.Main(
        [
            html.Nav(
                [
                    html.A("Autoimmune Atlas", href="#", className="app-brand"),
                    html.Div(
                        [
                            html.A("Comparison", href="#comparison"),
                            html.A("Evidence", href="#evidence"),
                            html.A("Open Targets ↗", href=source, target="_blank", rel="noreferrer"),
                        ],
                        className="app-nav",
                    ),
                ],
                className="app-bar",
                **{"aria-label": "Main navigation"},
            ),
            html.Header(
                [
                    html.Div(
                        [
                            html.P("AUTOIMMUNE DISEASE / DRUG TARGETS", className="eyebrow"),
                            html.H1("Drug target expression by cell type"),
                            html.P("Explore which cell types express the targets of drugs reaching Phase III or later in each autoimmune disease.", className="lede"),
                        ],
                    ),
                    html.Div(
                        [
                            html.Div([html.Strong(str(len(snapshot["diseases"]))), html.Span("Disease terms")], className="stat"),
                            html.Div([html.Strong(str(len(snapshot["records"]))), html.Span("Drug–target records")], className="stat"),
                            html.Div([html.Strong(str(len(snapshot["expression"]))), html.Span("Targets queried")], className="stat"),
                        ],
                        className="stats",
                    ),
                ],
                className="hero",
            ),
            html.Section(
                [
                    html.Div([html.Span("Data status", className="meta-label"), html.Strong(f"Ready {ready} / unavailable {unavailable}")]),
                    html.Div([html.Span("Retrieved at", className="meta-label"), html.Strong(snapshot["retrieved_at"])]),
                    html.Div([html.Span("Unknown / withdrawn stage records", className="meta-label"), html.Strong(str(unclassified))]),
                    html.Div([html.Span("Data source", className="meta-label"), html.A(f"Open Targets {data_version}".strip(), href=source, target="_blank", rel="noreferrer")]),
                ],
                className="source-bar",
            ),
            html.Section(
                [
                    html.Div(
                        [
                            html.Div([html.Label("Measure", htmlFor="measure"), dcc.RadioItems(id="measure", options=[{"label": "Target count", "value": "count"}, {"label": "Share of known targets within disease (%)", "value": "percent"}], value="count", inline=True)], className="control"),
                            html.Div([html.Label("Drug modality", htmlFor="modality"), dcc.Dropdown(id="modality", options=[{"label": label, "value": value} for label, value in (("All", "all"),) + atlas.DRUG_TYPE_MODALITIES], value="all", clearable=False)], className="control"),
                            html.Div([html.Label("Expression threshold (median CPM > value)", htmlFor="threshold"), dcc.Input(id="threshold", type="number", min=0, step=0.1, value=DEFAULT_EXPRESSION_THRESHOLD)], className="control"),
                        ],
                        className="control-grid compact",
                    ),
                    html.Details(
                        [
                            html.Summary("Select diseases and cell types"),
                            html.Div(
                                [
                                    html.Div([html.Label("Diseases", htmlFor="diseases"), dcc.Dropdown(id="diseases", options=[{"label": name, "value": item_id} for item_id, name in diseases], value=default_diseases, multi=True, searchable=True)], className="control"),
                                    html.Div([html.Label("Cell types", htmlFor="cells"), dcc.Dropdown(id="cells", options=[{"label": name, "value": item_id} for item_id, name in cells], value=default_cells, multi=True, searchable=True)], className="control"),
                                ],
                                className="control-grid selectors",
                            ),
                        ],
                        className="filter-details",
                    ),
                    html.Details(
                        [
                            html.Summary("How to read this view"),
                            html.P("This view uses healthy-donor Tabula Sapiens pseudobulk data and selects targets whose median across donors exceeds the threshold. It does not establish expression in diseased tissue or therapeutic efficacy."),
                            html.P("Differences between diseases reflect their drug targets, not disease-specific changes in expression."),
                            html.P("Example: ≥42 means at least 42 targets; unresolved targets may increase the count. A value of 42 is fully assessed within the selected scope. Percentages carry ≥ only when known targets in the denominator remain unassessed."),
                        ],
                        className="methods-note",
                    ),
                ],
                className="panel controls",
            ),
            html.Section(
                [
                    html.Div([html.H2("Expressed drug targets by disease and cell type"), html.Button("Download CSV", id="download-button", n_clicks=0), dcc.Download(id="download")], className="section-heading"),
                    html.P("Color shows the number or share of targets above the expression threshold, not expression levels or drug efficacy.", className="matrix-note"),
                    html.P(id="matrix-note", className="matrix-note"),
                    html.P("≥ Lower bound · × Unknown / missing · Click a cell to inspect its evidence", className="matrix-note"),
                    html.Div(dcc.Graph(id="heatmap", config={"displaylogo": False, "responsive": False}), className="graph-scroll"),
                ],
                className="panel matrix-panel",
                id="comparison",
            ),
            html.Section(
                [
                    html.Div([html.H2("Selection details"), html.P("Use the selectors below to open the same details without clicking the heatmap.")], className="section-heading detail-heading"),
                    html.Div(
                        [
                            html.Div([html.Label("Disease", htmlFor="detail-disease"), dcc.Dropdown(id="detail-disease", clearable=False)], className="control"),
                            html.Div([html.Label("Cell type", htmlFor="detail-cell"), dcc.Dropdown(id="detail-cell", clearable=False)], className="control"),
                        ],
                        className="control-grid selectors",
                    ),
                    html.Div(id="details", className="details"),
                ],
                className="panel",
                id="evidence",
            ),
            html.Footer("Counts represent distinct targets, not drugs. Percentage denominators include known targets only; drugs without mapped targets are excluded."),
        ],
        className="shell",
    )


try:
    SNAPSHOT = load_snapshot(SNAPSHOT_PATH)
    LOAD_ERROR: Exception | None = None
except (OSError, json.JSONDecodeError, ValueError, TypeError, KeyError) as error:
    SNAPSHOT = None
    LOAD_ERROR = error

app = Dash(__name__, title="Autoimmune Target Expression Atlas")
app.layout = unavailable_layout(LOAD_ERROR) if SNAPSHOT is None else dashboard_layout(SNAPSHOT)


if SNAPSHOT is not None:

    @app.callback(
        Output("heatmap", "figure"),
        Output("matrix-note", "children"),
        Input("measure", "value"),
        Input("modality", "value"),
        Input("threshold", "value"),
        Input("diseases", "value"),
        Input("cells", "value"),
    )
    def update_figure(measure, modality, threshold, disease_ids, cell_ids):
        """現在の条件で heatmap と欠測の説明を更新する。"""
        rows = visible_rows(SNAPSHOT, modality, threshold, disease_ids, cell_ids)
        figure = build_figure(rows, disease_ids or [], cell_ids or [], measure)
        missing = sum((row["percent"] if measure == "percent" else row["count"]) is None for row in rows)
        partial = sum(_is_lower_bound(row, measure) for row in rows)
        note = f"Median CPM > {effective_threshold(threshold):g} | {len(rows)} combinations | Lower bounds {partial} | Unknown / missing {missing}"
        if measure == "percent":
            unmapped = sum(
                max((row["unmapped_drugs"] for row in rows if row["disease_id"] == disease_id), default=0)
                for disease_id in set(row["disease_id"] for row in rows)
            )
            note += f" | Denominator: known targets only | Disease–drug records without mapped targets: {unmapped}"
        return figure, note

    @app.callback(
        Output("detail-disease", "options"),
        Output("detail-disease", "value"),
        Output("detail-cell", "options"),
        Output("detail-cell", "value"),
        Input("diseases", "value"),
        Input("cells", "value"),
        State("detail-disease", "value"),
        State("detail-cell", "value"),
    )
    def update_detail_selectors(disease_ids, cell_ids, current_disease, current_cell):
        """詳細欄の選択肢を heatmap の表示範囲へ揃える。"""
        disease_ids = disease_ids or []
        cell_ids = cell_ids or []
        disease_names = {row["id"]: row["name"] for row in SNAPSHOT["diseases"]}
        cell_names = dict(cell_catalog(SNAPSHOT))
        return (
            [{"label": disease_names[item], "value": item} for item in disease_ids],
            current_disease if current_disease in disease_ids else next(iter(disease_ids), None),
            [{"label": cell_names[item], "value": item} for item in cell_ids],
            current_cell if current_cell in cell_ids else next(iter(cell_ids), None),
        )

    @app.callback(
        Output("details", "children"),
        Input("heatmap", "clickData"),
        Input("detail-disease", "value"),
        Input("detail-cell", "value"),
        Input("measure", "value"),
        Input("modality", "value"),
        Input("threshold", "value"),
        Input("diseases", "value"),
        Input("cells", "value"),
    )
    def update_details(click_data, detail_disease, detail_cell, _measure, modality, threshold, disease_ids, cell_ids):
        """クリックまたは選択欄から、現在の条件に合う詳細だけを表示する。"""
        rows = visible_rows(SNAPSHOT, modality, threshold, disease_ids, cell_ids)
        triggered = ctx.triggered_id
        if triggered in FILTER_IDS:
            selection = None
        else:
            selection = resolve_selection(triggered, click_data, detail_disease, detail_cell, rows)
        return detail_panel(rows, selection)

    @app.callback(
        Output("download", "data"),
        Input("download-button", "n_clicks"),
        State("measure", "value"),
        State("modality", "value"),
        State("threshold", "value"),
        State("diseases", "value"),
        State("cells", "value"),
        prevent_initial_call=True,
    )
    def download_csv(_clicks, measure, modality, threshold, disease_ids, cell_ids):
        """表示中の matrix と対応する根拠を CSV にする。"""
        rows = visible_rows(SNAPSHOT, modality, threshold, disease_ids, cell_ids)
        content = "\ufeff" + atlas.to_csv(
            export_rows(
                rows,
                measure,
                modality_filter=modality,
                expression_threshold=effective_threshold(threshold),
                snapshot=SNAPSHOT,
            ),
        )
        return {"content": content, "filename": "autoimmune-drug-cell-matrix.csv", "type": "text/csv;charset=utf-8"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Autoimmune Atlas")
    parser.add_argument("--dev", action="store_true", help="Automatically reload when code changes")
    args = parser.parse_args()
    app.run(host="127.0.0.1", port=8050, debug=args.dev)
