"""自己免疫疾患の薬剤標的と細胞型を閲覧する Dash アプリ。"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timedelta, timezone
from pathlib import Path

import plotly.graph_objects as go
from dash import Dash, Input, Output, State, ctx, dcc, html

import atlas
from disease_catalog import disease_catalog, ordered_disease_ids

BASE_DIR = Path(__file__).parent
SNAPSHOT_PATH = BASE_DIR / "data" / "snapshot.json"
DEFAULT_EXPRESSION_THRESHOLD = 0.5
DEFAULT_SPECIFICITY_THRESHOLD = 0.75
SOURCE_PAGE_SIZE = 50
PARAMETER_IDS = ("measure", "modality", "stage", "method", "threshold", "specificity", "diseases", "cells", "heatmap-view")
STAGE_LABELS = {"phase1": "Phase I or later", "phase2": "Phase II or later", "phase3": "Phase III or later", "approved": "Approval reached"}
METHOD_LABELS = {"fixed": "Fixed CPM", "relative": "Fixed CPM + Target-relative median", "specificity": "Fixed CPM + CELLEX specificity"}


def info_tip(key: str, label: str, description: str):
    """操作対象の近くに、キーボードでも読める説明を置く。"""
    tip_id = f"{key}-tip"
    return html.Span([
        html.Button("ⓘ", type="button", className="info-button", **{"aria-label": f"About {label}", "aria-describedby": tip_id}),
        html.Span(description, id=tip_id, role="tooltip", className="info-content"),
    ], className="info-tip")


def load_snapshot(path: Path) -> dict | None:
    """schema 2 の保存済みスナップショットを読む。"""
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as handle:
        snapshot = json.load(handle)
    if not isinstance(snapshot, dict):
        raise ValueError("snapshot.json must contain a JSON object")
    if snapshot.get("schema") != 2:
        raise ValueError(
            "Unsupported snapshot format. Refresh it with: pixi run refresh"
        )
    required = {"root", "retrieved_at", "source", "diseases", "records", "expression"}
    missing = required - snapshot.keys()
    if missing:
        raise ValueError("Missing required fields in snapshot.json: " + ", ".join(sorted(missing)))
    if not isinstance(snapshot["diseases"], list) or not isinstance(snapshot["records"], list):
        raise ValueError("diseases and records must be arrays")
    if not isinstance(snapshot["expression"], dict):
        raise ValueError("expression must be an object")
    return snapshot


def cell_catalog(snapshot: dict, level: str = "group") -> list[tuple[str, str]]:
    """公開関数が返す細胞カタログを選択欄向けに整える。"""
    return [(row["id"], row["name"]) for row in atlas.cell_catalog(snapshot, level)]


def _ordered_cell_ids(snapshot: dict, level: str, cell_ids) -> list[str]:
    selected = set(cell_ids or [])
    return [cell_id for cell_id, _ in cell_catalog(snapshot, level) if cell_id in selected]


def _disease_checklist_sections(family):
    """疾患本体と詳細の選択欄を、描画と同期で同じ範囲に分ける。"""
    main = [d for d in family["diseases"] if d["id"] == family["id"]]
    details = [d for d in family["diseases"] if d["id"] != family["id"]]
    sections = [{"id": f"disease-family-{family['id']}", "diseases": main or details}]
    if main and details:
        sections.append({"id": f"disease-details-{family['id']}", "diseases": details})
    return sections


def disease_selector(snapshot, selected):
    """検索欄と、群から開けるチェック欄を同じ選択へ結び付ける。"""
    catalog = disease_catalog(snapshot)
    options, groups = [], []
    for group in catalog:
        families = []
        for family in group["families"]:
            choices = [{"label": d["name"], "value": d["id"]} for d in family["diseases"]]
            options.extend(choices)
            sections = _disease_checklist_sections(family)
            checklists = [dcc.Checklist(id=section["id"], options=[{"label": d["name"], "value": d["id"]} for d in section["diseases"]], value=[d["id"] for d in section["diseases"] if d["id"] in selected], className="disease-checklist") for section in sections]
            contents = checklists[0]
            if len(checklists) > 1:
                contents = html.Div([checklists[0], html.Details([html.Summary(f"{family['label']} details ({len(sections[1]['diseases'])} terms)"), checklists[1]])], className="disease-families")
            families.append(html.Details([html.Summary(f"{family['label']} ({len(choices)} terms)"), contents]) if len(choices) > 1 else contents)
        groups.append(html.Details([html.Summary(group["label"]), html.Div(families, className="disease-families")]))
    return html.Div([
        html.Div([html.Label("Diseases", htmlFor="diseases"), info_tip("diseases", "disease groups", "Browse groups to select individual diseases and related terms, or search by name. Each checkbox selects only that term; parent and child terms are never combined. Groups organize browsing, not diagnostic classification.")], className="label-help"),
        dcc.Dropdown(id="diseases", options=options, value=ordered_disease_ids(snapshot, selected), multi=True, searchable=True),
        html.Details([html.Summary("Browse disease groups"), html.Div(groups, className="disease-tree")], className="disease-browser"),
    ], className="control")


def _cell_selection_catalog(snapshot):
    """大分類と元細胞を、同じ ID の場合も別の選択肢として定義する。"""
    catalog = {cell["id"]: cell for cell in atlas.cell_catalog(snapshot, "mixed")}
    return [{"name": group["name"], "sections": [
        {"id": f"cell-group-{group['id']}", "cells": [catalog["group:" + group["id"]]]},
        {"id": f"cell-details-{group['id']}", "cells": [catalog[member] for member in group["members"]]},
    ]} for group in atlas.cell_catalog(snapshot, "group")]


def cell_selector(snapshot, selected):
    """疾患と同じ開閉操作で、大分類と細分類を独立して選べるようにする。"""
    groups = []
    for group in _cell_selection_catalog(snapshot):
        checklists = [dcc.Checklist(id=section["id"], options=[{"label": cell["name"], "value": cell["id"]} for cell in section["cells"]], value=[cell["id"] for cell in section["cells"] if cell["id"] in selected], className="disease-checklist") for section in group["sections"]]
        groups.append(html.Details([
            html.Summary(group["name"]),
            html.Div([checklists[0], html.Details([html.Summary(f"{group['name']} details ({len(group['sections'][1]['cells'])} cells)"), checklists[1]])], className="disease-families"),
        ]))
    return html.Div([
        html.Div([html.Label("Cells", htmlFor="cells"), info_tip("cells", "cells", "Select groups and source cells independently. A group counts the union of targets or drugs meeting the rule in any member; CPM values are never added or averaged. A source cell uses only its own expression. Selecting a group does not select its members. Rows follow lineage order, not support strength.")], className="label-help"),
        dcc.Dropdown(id="cells", options=[{"label": name, "value": cell_id} for cell_id, name in cell_catalog(snapshot, "mixed")], value=selected, multi=True, searchable=True),
        html.Details([html.Summary("Browse cell groups"), html.Div(groups, className="disease-tree")], className="disease-browser"),
    ], className="control")


def choose_defaults(items, terms, limit, preferred_ids=None) -> list[str]:
    """名前と実データの有無から、存在する項目だけを初期選択する。"""
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


def _default_cells(snapshot: dict) -> list[str]:
    return ["group:" + cell["id"] for cell in atlas.cell_catalog(snapshot, "group")]


def format_data_version(value: object) -> str:
    """Open Targets の版を画面と CSV で使える短い文字列にする。"""
    if not isinstance(value, dict):
        return str(value or "")
    return ".".join(str(part) for part in (value.get("year"), value.get("month"), value.get("iteration")) if part is not None)


def _effective_number(value, default, label, maximum=None):
    """空欄には初期値を使い、不正値は理由を示して初期値へ戻す。"""
    if value is None or value == "":
        return default, None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default, f"Invalid {label}; using {default:g}."
    if isinstance(value, bool) or not math.isfinite(number) or number < 0 or (maximum is not None and number > maximum):
        bound = f" between 0 and {maximum:g}" if maximum is not None else " at or above 0"
        return default, f"{label.capitalize()} must be finite and{bound}; using {default:g}."
    return number, None


def effective_filters(threshold, specificity):
    """実際に集計へ渡す閾値と入力エラーを返す。"""
    minimum, minimum_error = _effective_number(threshold, DEFAULT_EXPRESSION_THRESHOLD, "minimum CPM")
    specificity_value, specificity_error = _effective_number(specificity, DEFAULT_SPECIFICITY_THRESHOLD, "specificity threshold", 1)
    return minimum, specificity_value, [message for message in (minimum_error, specificity_error) if message]


def effective_threshold(threshold) -> float:
    """後方互換用に、適用される最低 CPM だけを返す。"""
    return effective_filters(threshold, DEFAULT_SPECIFICITY_THRESHOLD)[0]


def visible_rows(snapshot, modality, threshold, disease_ids, cell_ids, *, stage="phase3", method="fixed", specificity=DEFAULT_SPECIFICITY_THRESHOLD, level="group"):
    """backend の集計結果を現在の表示範囲へ絞る。"""
    minimum, specificity_value, _ = effective_filters(threshold, specificity)
    return atlas.summarize(snapshot, modality, minimum, stage=stage, method=method, specificity_threshold=specificity_value, level=level, cell_ids=cell_ids or [], disease_ids=disease_ids or [])


def _measure_fields(kind: str, measure: str):
    if kind == "drug":
        return "drug_percent" if measure == "percent" else "drug_count", "unknown_drugs", "drug_denominator"
    return "percent" if measure == "percent" else "count", "unknown", "denominator"


def _is_lower_bound(row: dict, measure: str, kind: str = "target") -> bool:
    """実数と割合で異なる下限判定を、その指標自身の未判定から求める。"""
    value_key, unknown_key, _ = _measure_fields(kind, measure)
    return row.get(value_key) is not None and (row.get(unknown_key, 0) > 0 or (measure == "count" and row.get("unmapped_drugs", 0) > 0))


def _display_value(row: dict, measure: str, kind: str = "target") -> str:
    value = row.get(_measure_fields(kind, measure)[0])
    if value is None:
        return ""
    rendered = f"{value:.1f}%" if measure == "percent" else str(value)
    return "≥" + rendered if _is_lower_bound(row, measure, kind) else rendered


def _hover_text(row: dict, measure: str, kind: str) -> str:
    value_key, unknown_key, denominator_key = _measure_fields(kind, measure)
    value = _display_value(row, measure, kind) or "No value"
    label = "Drugs" if kind == "drug" else "Targets"
    if row["status"] == "unavailable":
        assessment = "Disease data unavailable"
    elif measure == "percent" and row.get(denominator_key, 0) == 0:
        assessment = "No eligible items in the percentage denominator"
    elif _is_lower_bound(row, measure, kind):
        assessment = "≥ is a lower bound; unresolved evidence may increase this value"
    elif row.get(value_key) is None:
        assessment = "Support cannot be determined from available target and expression data."
    elif row[value_key] == 0:
        assessment = f"No qualifying {label.lower()} under this rule"
    else:
        assessment = f"{label} meeting the expression rule"
    return "<br>".join((
        f"<b>{row['disease']} × {row['cell']}</b>", f"{label}: {value}", assessment,
        f"Denominator: {row.get(denominator_key, 0)}", f"Unknown in denominator: {row.get(unknown_key, 0)}",
        f"Mapped / all canonical drugs: {row.get('mapped_drugs', 0)} / {row.get('total_drugs', 0)}",
        f"Drugs without mapped targets: {row.get('unmapped_drugs', 0)}", f"Status: {row['status']}",
    ))


def build_figure(rows, disease_ids, cell_ids, measure, kind="target") -> go.Figure:
    """0、下限値、欠測を区別した target または drug heatmap を作る。"""
    value_key = _measure_fields(kind, measure)[0]
    lookup = {(row["disease_id"], row["cell_id"]): row for row in rows}
    disease_names = {row["disease_id"]: row["disease"] for row in rows}
    cell_names = {row["cell_id"]: row["cell"] for row in rows}
    z, texts, hovers, customs = [], [], [], []
    missing_x, missing_y, missing_custom, missing_hover = [], [], [], []
    for cell_id in cell_ids:
        z_row, text_row, hover_row, custom_row = [], [], [], []
        for disease_id in disease_ids:
            row = lookup.get((disease_id, cell_id))
            value = None if row is None else row.get(value_key)
            z_row.append(value)
            text_row.append("" if row is None else _display_value(row, measure, kind).removeprefix("≥"))
            hover_row.append("No summary for this combination" if row is None else _hover_text(row, measure, kind))
            custom_row.append([disease_id, cell_id])
            if value is None:
                missing_x.append(disease_names.get(disease_id, disease_id)); missing_y.append(cell_names.get(cell_id, cell_id)); missing_custom.append([disease_id, cell_id]); missing_hover.append("No summary for this combination" if row is None else _hover_text(row, measure, kind))
        z.append(z_row); texts.append(text_row); hovers.append(hover_row); customs.append(custom_row)
    values = [value for z_row in z for value in z_row if value is not None]
    title = "Drug" if kind == "drug" else "Target"
    figure = go.Figure(go.Heatmap(
        x=[disease_names.get(item, item) for item in disease_ids], y=[cell_names.get(item, item) for item in cell_ids], z=z,
        zmin=0, zmax=100 if measure == "percent" else max(values, default=1) or 1,
        colorscale=((0, "#f5f7fa"), (.35, "#9fc4e0"), (1, "#1769aa")), colorbar={"title": f"{title} {'share (%)' if measure == 'percent' else 'count'}"},
        text=texts, texttemplate="%{text}", hovertext=hovers, hovertemplate="%{hovertext}<extra></extra>", customdata=customs, xgap=2, ygap=2,
    ))
    if missing_x:
        figure.add_trace(go.Scatter(x=missing_x, y=missing_y, customdata=missing_custom, hovertext=missing_hover, mode="markers", marker={"symbol": "x", "size": 10, "color": "#8a99a6"}, name="No value (see hover)", hovertemplate="%{hovertext}<extra></extra>"))
    if not disease_ids or not cell_ids:
        figure.add_annotation(text="Select diseases and cell types", showarrow=False)
    figure.update_layout(template="plotly_white", font={"family": "Arial, sans-serif", "size": 12, "color": "#263238"}, width=max(760, min(2800, 260 + 125 * len(disease_ids))), height=max(420, min(1400, 180 + 28 * len(cell_ids))), margin={"l": 180, "r": 40, "t": 130, "b": 70}, legend={"orientation": "h", "y": -.08, "yanchor": "top", "x": 0}, hoverlabel={"align": "left"})
    figure.update_xaxes(tickangle=-45, side="top", title="Disease", automargin=True)
    figure.update_yaxes(autorange="reversed", title="Cell type")
    return figure


def resolve_selection(triggered_id, click_data, disease_id, cell_id, rows):
    """表示中の組合せに限って、クリックまたは選択欄を受け付ける。"""
    visible = {(row["disease_id"], row["cell_id"]) for row in rows}
    selection = None
    if triggered_id in {"target-heatmap", "drug-heatmap", "heatmap"} and click_data and click_data.get("points"):
        custom = click_data["points"][0].get("customdata")
        if isinstance(custom, (list, tuple)) and len(custom) >= 2:
            selection = custom[0], custom[1]
    elif disease_id and cell_id:
        selection = disease_id, cell_id
    return selection if selection in visible else None


def evidence_rows(snapshot, row, *, modality, stage, threshold, method, specificity, metadata=None, records=None, cell_names=None):
    """選択範囲の全元記録へ、公開判定関数による細胞別状態を添える。"""
    metadata = metadata if metadata is not None else atlas.expression_metadata(snapshot)
    records = records if records is not None else [record for record in atlas.filtered_records(snapshot, modality, stage) if record["disease_id"] == row["disease_id"]]
    cell_names = cell_names if cell_names is not None else {item["id"]: item["name"] for item in atlas.cell_catalog(snapshot, "cell")}
    members = _ordered_cell_ids(snapshot, "cell", row["member_cell_ids"])
    output = []
    for record in records:
        base = {**record, "group_cell_id": row["cell_id"], "group_cell": row["cell"], "references_json": json.dumps(record.get("references") or [], ensure_ascii=False, sort_keys=True)}
        if not record.get("target_id"):
            output.append({**base, "evidence_cell_id": "", "evidence_cell": "", "cpm": None, "specificity_score": None, "target_median": None, "support_state": "unmapped", "contributing": False})
            continue
        for member in members:
            item = metadata.get((record["target_id"], member))
            state = atlas.expression_state(item, threshold, method, specificity)
            output.append({
                **base, "evidence_cell_id": member, "evidence_cell": item.get("cell", cell_names.get(member, member)) if item else cell_names.get(member, member),
                "cpm": item.get("median") if item else None, "specificity_score": item.get("specificity_score") if item else None,
                "target_median": item.get("target_median") if item else None,
                "support_state": "positive" if state is True else "negative" if state is False else "unknown", "contributing": state is True,
            })
    return output


def expression_figure(snapshot, records, metadata=None) -> go.Figure:
    """全元細胞について、現在の薬剤条件に含まれる標的の連続発現量を示す。"""
    metadata = metadata if metadata is not None else atlas.expression_metadata(snapshot)
    target_names = {}
    for record in records:
        if record.get("target_id"):
            target_names.setdefault(record["target_id"], record.get("target") or record["target_id"])
    cells = atlas.cell_catalog(snapshot, "cell")
    members = [cell["id"] for cell in cells]
    member_names = {cell["id"]: cell["name"] for cell in cells}
    targets = sorted(target_names, key=lambda item: (target_names[item].casefold(), item))
    z, hover, missing_x, missing_y = [], [], [], []
    for member in members:
        z_row, hover_row = [], []
        for target_id in targets:
            label = f"{target_names[target_id]} ({target_id})"
            item = metadata.get((target_id, member)); cpm = item.get("median") if item else None
            z_row.append(math.log2(1 + cpm) if cpm is not None else None)
            hover_row.append("<br>".join((
                f"<b>{label}</b>", f"Cell: {member_names[member]} ({member})", f"Median CPM: {cpm:g}" if cpm is not None else "Median CPM: missing",
                f"CELLEX specificity: {item.get('specificity_score'):g}" if item and item.get("specificity_score") is not None else "CELLEX specificity: missing",
                f"Target-relative median: {item.get('target_median'):g}" if item and item.get("target_median") is not None else "Target-relative median: missing",
            )))
            if cpm is None:
                missing_x.append(label); missing_y.append(member_names[member])
        z.append(z_row); hover.append(hover_row)
    target_labels = [f"{target_names[t]} ({t})" for t in targets]
    figure = go.Figure(go.Heatmap(x=target_labels, y=[member_names[m] for m in members], z=z, zmin=0, colorscale="Blues", colorbar={"title": "log2(1 + CPM)"}, hovertext=hover, hovertemplate="%{hovertext}<extra></extra>", xgap=2, ygap=2))
    if missing_x:
        figure.add_trace(go.Scatter(x=missing_x, y=missing_y, mode="markers", marker={"symbol": "x", "size": 9, "color": "#8a99a6"}, name="Missing expression", hovertemplate="Expression missing<extra></extra>"))
    if not targets:
        figure.add_annotation(text="No known targets in the selected scope", showarrow=False)
    figure.update_layout(template="plotly_white", width=max(900, 420 + 32 * len(targets)), height=max(360, 170 + 25 * len(members)), margin={"l": 280, "r": 40, "t": 110, "b": 60}, font={"family": "Arial, sans-serif", "size": 11, "color": "#263238"}, hoverlabel={"align": "left"})
    figure.update_xaxes(tickangle=-45, side="top", title="Known target", tickvals=target_labels, ticktext=[target_names[t] for t in targets], automargin=True)
    figure.update_yaxes(autorange="reversed", title="Source cell type", automargin=True)
    return figure


def _reference_links(record):
    links = []
    for reference in record.get("references") or []:
        source, urls, ids = reference.get("source") or "Source", reference.get("urls") or [], reference.get("ids") or []
        if urls:
            links.extend(html.A(f"{source}: {ids[i] if i < len(ids) else i + 1}", href=url, target="_blank", rel="noreferrer") for i, url in enumerate(urls))
        else:
            links.append(html.Span(f"{source}: {', '.join(map(str, ids))}" if ids else source))
    if record.get("target_id"):
        links.append(html.A("Open Targets target", href="https://platform.opentargets.org/target/" + record["target_id"], target="_blank", rel="noreferrer"))
    if record.get("drug_id"):
        links.append(html.A("Open Targets drug", href="https://platform.opentargets.org/drug/" + record["drug_id"], target="_blank", rel="noreferrer"))
    if record.get("canonical_drug_id") and record.get("canonical_drug_id") != record.get("drug_id"):
        links.append(html.A("Open Targets canonical drug", href="https://platform.opentargets.org/drug/" + record["canonical_drug_id"], target="_blank", rel="noreferrer"))
    if record.get("disease_id"):
        links.append(html.A("Open Targets disease", href="https://platform.opentargets.org/disease/" + record["disease_id"], target="_blank", rel="noreferrer"))
    children = []
    for index, item in enumerate(links):
        if index:
            children.append(html.Br())
        children.append(item)
    return html.Div(children) if children else "—"


def _evidence_table(records):
    headers = ("Expression result", "Canonical drug", "Original drug", "Modality", "Record stage / highest disease-specific drug stage", "Target", "Action / mechanism", "Source cell", "Median CPM", "CELLEX specificity", "Reference-cell median CPM", "Sources")
    state_labels = {"positive": "Meets rule", "negative": "Does not meet rule", "unknown": "Not assessed", "unmapped": "No mapped target"}
    body = []
    for record in records:
        actions, mechanism = ", ".join(record.get("action_types") or []), record.get("mechanism") or "—"
        body.append(html.Tr([
            html.Td(state_labels[record["support_state"]]), html.Td(f"{record.get('canonical_drug') or '—'} ({record.get('canonical_drug_id') or '—'})"),
            html.Td(f"{record.get('drug') or '—'} ({record.get('drug_id') or '—'})"), html.Td(record.get("drug_type") or "—"),
            html.Td(f"{record.get('stage') or '—'} / {record.get('canonical_stage') or '—'}"), html.Td(f"{record.get('target') or '—'} ({record.get('target_id') or '—'})"),
            html.Td(f"{actions or '—'} / {mechanism}"), html.Td(f"{record.get('evidence_cell') or '—'} ({record.get('evidence_cell_id') or '—'})"),
            html.Td("—" if record.get("cpm") is None else f"{record['cpm']:g}"), html.Td("—" if record.get("specificity_score") is None else f"{record['specificity_score']:g}"),
            html.Td("—" if record.get("target_median") is None else f"{record['target_median']:g}"), html.Td(_reference_links(record)),
        ]))
    return html.Div(html.Table([html.Thead(html.Tr([html.Th(item) for item in headers])), html.Tbody(body)]), className="table-scroll")


def detail_panel(rows, selection, snapshot=None, *, modality="all", stage="phase3", threshold=DEFAULT_EXPRESSION_THRESHOLD, method="fixed", specificity=DEFAULT_SPECIFICITY_THRESHOLD):
    """選択した組合せの集計、連続発現、全元記録を表示する。"""
    if selection is None:
        return html.Div("Select a disease and cell type using the displayed heatmap or the selectors above.", className="empty-note")
    row = next((item for item in rows if (item["disease_id"], item["cell_id"]) == selection), None)
    if row is None:
        return html.Div("The previous selection is outside the current filters. Select a visible cell.", className="empty-note")
    summary = html.Div([html.H3(f"{row['disease']} × {row['cell']}"), html.P(
        f"Targets meeting rule: {_display_value(row, 'count') or 'NA'}; {_display_value(row, 'percent') or 'NA'} of {row['denominator']} known targets ({row['unknown']} not assessed). "
        f"Canonical drugs meeting rule: {_display_value(row, 'count', 'drug') or 'NA'}; {_display_value(row, 'percent', 'drug') or 'NA'} of {row['drug_denominator']} drugs with known targets ({row['unknown_drugs']} not assessed). "
        f"Drug mapping: {row['mapped_drugs']} of {row['total_drugs']} canonical drugs have known targets; {row['unmapped_drugs']} have none.", className="detail-summary")])
    if snapshot is None:
        return html.Div([summary, html.P("Underlying schema 2 data are required for evidence details.", className="empty-note")])
    metadata = atlas.expression_metadata(snapshot)
    filtered = [record for record in atlas.filtered_records(snapshot, modality, stage) if record["disease_id"] == row["disease_id"]]
    cell_names = {item["id"]: item["name"] for item in atlas.cell_catalog(snapshot, "cell")}
    records = evidence_rows(snapshot, row, modality=modality, stage=stage, threshold=threshold, method=method, specificity=specificity, metadata=metadata, records=filtered, cell_names=cell_names)
    positive = [record for record in records if record["contributing"]]
    positive_pairs = {}
    for record in positive:
        key = record.get("canonical_drug_id"), record.get("target_id")
        item = positive_pairs.setdefault(key, {"drug": record.get("canonical_drug") or record.get("drug"), "target": record.get("target") or record.get("target_id"), "cells": {}})
        item["cells"][record["evidence_cell_id"]] = record.get("evidence_cell") or record["evidence_cell_id"]
    pair_table = html.Div(html.Table([
        html.Thead(html.Tr([html.Th(label, scope="col") for label in ("Drug", "Target", "Source cells meeting rule")])),
        html.Tbody([html.Tr([
            html.Td(item["drug"]), html.Td(item["target"]),
            html.Td(html.Details([
                html.Summary(f"{len(item['cells'])} {'cell' if len(item['cells']) == 1 else 'cells'}"),
                html.Ul([html.Li(cell) for cell in item["cells"].values()]),
            ], open=False)),
        ]) for item in positive_pairs.values()]),
    ], **{"aria-label": "Drug–target pairs meeting the expression rule"}), className="table-scroll pair-table") if positive_pairs else html.P("No drug–target pairs meet the applied rule.", className="empty-note")
    source_context = {key: row[key] for key in ("disease_id", "cell_id", "cell", "member_cell_ids")}
    source_context.update(modality=modality, stage=stage, threshold=threshold, method=method, specificity=specificity)
    pages = max(1, math.ceil(len(records) / SOURCE_PAGE_SIZE))
    return html.Div([
        summary, html.H3(["Expression of targets across all source cell types", info_tip("expression", "target expression", "Targets come from the selected disease and drug filters. All healthy reference source cell types are shown, regardless of the selected cell above. Color is log2(1 + median CPM); values are not from disease samples.")]),
        html.Div(dcc.Graph(figure=expression_figure(snapshot, records, metadata), config={"displaylogo": False, "responsive": False}), className="graph-scroll expression-graph"),
        html.H3("Drug–target pairs meeting the expression rule"),
        html.P("Targets meeting the expression rule in the selected cells, and drugs acting on those targets. Select a cell count to show cell names.", className="detail-summary"),
        pair_table,
        html.H3(["Filtered drug records by source cell type", info_tip("source-records", "source records", "Rows repeat for each source cell type. Records that do not meet the rule, cannot be assessed, or lack a mapped target remain visible.")]),
        html.Details([
            html.Summary("Show source records", id="source-toggle", n_clicks=0),
            dcc.Store(id="source-context", data=source_context),
            html.Div([
                html.Label("Page", htmlFor="source-page"),
                dcc.Dropdown(id="source-page", options=[{"label": f"{page} / {pages}", "value": page} for page in range(1, pages + 1)], value=1, clearable=False, searchable=False),
            ], className="source-pagination"),
            dcc.Loading(html.Div(id="source-records-page")),
        ], open=False, className="source-records"),
    ])


def export_rows(rows, measure, *, modality_filter="all", stage_filter="phase3", method="fixed", expression_threshold=DEFAULT_EXPRESSION_THRESHOLD, specificity_threshold=DEFAULT_SPECIFICITY_THRESHOLD, level="group", snapshot=None):
    """集計、元薬剤記録、発現根拠を重複させず一つの CSV に並べる。"""
    snapshot = snapshot or {}
    base = {
        "row_type": "", "disease_id": "", "disease": "", "group_cell_id": "", "group_cell": "", "ontology_cell_id": "", "evidence_cell_id": "", "evidence_cell": "",
        "measure": measure, "modality_filter": modality_filter, "stage_filter": stage_filter, "method": method,
        "minimum_cpm": expression_threshold, "specificity_threshold": specificity_threshold, "level": level,
        "data_version": format_data_version(snapshot.get("data_version")), "retrieved_at": snapshot.get("retrieved_at", ""), "source": snapshot.get("source", ""),
        "expression_source": "Tabula Sapiens", "expression_unit": "CPM(pseudobulk sum[counts])", "expression_summary": "donor median",
        "target_display": "", "target_count": "", "target_percent": "", "target_denominator": "", "unknown_targets": "",
        "drug_display": "", "drug_count": "", "drug_percent": "", "drug_denominator": "", "unknown_drugs": "",
        "mapped_drugs": "", "total_drugs": "", "unmapped_drugs": "", "status": "", "support_state": "", "contributing": "",
        "canonical_drug_id": "", "canonical_drug": "", "drug_id": "", "drug": "", "drug_type": "", "original_stage": "", "canonical_stage": "",
        "target_id": "", "target": "", "mechanism": "", "action_types": "", "raw_cpm": "", "specificity_score": "", "target_relative_median": "", "references": "",
    }
    output = []
    for row in rows:
        output.append({
            **base, "row_type": "summary", "disease_id": row["disease_id"], "disease": row["disease"], "group_cell_id": row["cell_id"], "group_cell": row["cell"],
            "ontology_cell_id": row.get("ontology_id", row["cell_id"]), "level": row.get("cell_level", level),
            "target_display": _display_value(row, measure), "target_count": row["count"], "target_percent": row["percent"], "target_denominator": row["denominator"], "unknown_targets": row["unknown"],
            "drug_display": _display_value(row, measure, "drug"), "drug_count": row["drug_count"], "drug_percent": row["drug_percent"], "drug_denominator": row["drug_denominator"], "unknown_drugs": row["unknown_drugs"],
            "mapped_drugs": row["mapped_drugs"], "total_drugs": row["total_drugs"], "unmapped_drugs": row["unmapped_drugs"], "status": row["status"], "support_state": "summary",
        })
    if snapshot.get("schema") != 2:
        return output
    selected_diseases = {row["disease_id"] for row in rows}
    records = [record for record in atlas.filtered_records(snapshot, modality_filter, stage_filter) if record["disease_id"] in selected_diseases]
    for record in records:
        output.append({
            **base, "row_type": "source_record", "disease_id": record["disease_id"], "disease": record["disease"],
            "support_state": "unmapped" if not record.get("target_id") else "source_record",
            "canonical_drug_id": record.get("canonical_drug_id", ""), "canonical_drug": record.get("canonical_drug", ""),
            "drug_id": record.get("drug_id", ""), "drug": record.get("drug", ""), "drug_type": record.get("drug_type", ""),
            "original_stage": record.get("stage", ""), "canonical_stage": record.get("canonical_stage", ""),
            "target_id": record.get("target_id", ""), "target": record.get("target", ""), "mechanism": record.get("mechanism", ""),
            "action_types": "|".join(record.get("action_types") or []), "references": json.dumps(record.get("references") or [], ensure_ascii=False, sort_keys=True),
        })
    metadata = atlas.expression_metadata(snapshot)
    cell_names = {item["id"]: item["name"] for item in atlas.cell_catalog(snapshot, "cell")}
    cell_rank = {cell_id: rank for rank, cell_id in enumerate(cell_names)}
    targets_by_disease = {}
    for record in records:
        if record.get("target_id"):
            targets_by_disease.setdefault(record["disease_id"], {})[record["target_id"]] = record.get("target") or record["target_id"]
    for row in rows:
        for target_id, target in targets_by_disease.get(row["disease_id"], {}).items():
            for member in sorted(row["member_cell_ids"], key=cell_rank.__getitem__):
                item = metadata.get((target_id, member))
                state = atlas.expression_state(item, expression_threshold, method, specificity_threshold)
                output.append({
                    **base, "row_type": "expression_evidence", "disease_id": row["disease_id"], "disease": row["disease"],
                    "group_cell_id": row["cell_id"], "group_cell": row["cell"], "evidence_cell_id": member,
                    "ontology_cell_id": row.get("ontology_id", row["cell_id"]), "level": row.get("cell_level", level),
                    "evidence_cell": item.get("cell", cell_names.get(member, member)) if item else cell_names.get(member, member),
                    "support_state": "positive" if state is True else "negative" if state is False else "unknown", "contributing": state is True,
                    "target_id": target_id, "target": target, "raw_cpm": item.get("median", "") if item else "",
                    "specificity_score": item.get("specificity_score", "") if item else "", "target_relative_median": item.get("target_median", "") if item else "",
                })
    return output


def unavailable_layout(error=None) -> html.Main:
    """データ未取得または再取得が必要な状態を説明する。"""
    title = "No data loaded yet" if error is None else "Data refresh required"
    return html.Main([html.Section([html.P("AUTOIMMUNE DRUG–CELL ATLAS", className="eyebrow"), html.H1(title), html.P("Refresh the snapshot to load the comparison."), html.Code("pixi run refresh"), *([html.P(str(error), className="empty-note")] if error else [])], className="unavailable-card")], className="shell unavailable")


def dashboard_layout(snapshot) -> html.Main:
    """schema 2 のデータから dashboard の初期画面を作る。"""
    diseases = [(row["id"], row["name"]) for row in snapshot["diseases"]]
    default_diseases = choose_defaults(diseases, (
        "systemic lupus erythematosus", "systemic sclerosis", "Sjogren syndrome",
        "rheumatoid arthritis", "multiple sclerosis", "myasthenia gravis",
        "type 1 diabetes mellitus", "Graves disease", "pemphigus", "autoimmune hepatitis",
    ), 10, {row["disease_id"] for row in snapshot["records"]})
    default_cells = _default_cells(snapshot)
    retrieved_at = datetime.fromisoformat(snapshot["retrieved_at"]).astimezone(timezone(timedelta(hours=9))).strftime("%Y-%m-%d %H:%M JST")
    ready = sum(row["status"] == "ready" for row in snapshot["diseases"]); data_version = format_data_version(snapshot.get("data_version")); source = "https://platform.opentargets.org/"
    modality_options = [{"label": "All", "value": "all"}] + [{"label": label, "value": value} for label, value in atlas.DRUG_TYPE_MODALITIES]
    control = lambda label, component: html.Div([html.Label(label, htmlFor=component.id), component], className="control")
    applied = dict(zip(PARAMETER_IDS, ("count", "all", "phase3", "fixed", DEFAULT_EXPRESSION_THRESHOLD, DEFAULT_SPECIFICITY_THRESHOLD, ordered_disease_ids(snapshot, default_diseases), default_cells, "target")))
    return html.Main([
        html.Nav([html.A("Autoimmune Atlas", href="#", className="app-brand"), html.Div([html.A("Comparison", href="#comparison"), html.A("Evidence", href="#evidence"), html.A("Open Targets ↗", href=source, target="_blank", rel="noreferrer")], className="app-nav")], className="app-bar", **{"aria-label": "Main navigation"}),
        html.Header([html.Div([html.P("AUTOIMMUNE DISEASE / DRUG TARGETS", className="eyebrow"), html.H1(["Drug targets by cell type", info_tip("overview", "the atlas", "Compare which drug targets meet an expression rule in healthy reference cells across diseases. Disease differences reflect eligible drug and target sets, not disease-specific expression. This does not establish treatment efficacy; disease records may include symptom or comorbidity treatment.")])]),
        html.Section([
            html.Div([html.Span("Disease terms", className="meta-label"), html.Strong(f"{ready} loaded / {len(snapshot['diseases'])} total")]),
            html.Div([html.Span("Drug records", className="meta-label"), html.Strong(str(len(snapshot["records"])))]),
            html.Div([html.Span("Targets", className="meta-label"), html.Strong(str(len(snapshot["expression"])))]),
            html.Div([html.Span("Retrieved at", className="meta-label"), html.Strong(retrieved_at, title=snapshot["retrieved_at"])]),
            html.Div([html.Span("Expression reference", className="meta-label"), html.Strong("Tabula Sapiens")]),
            html.Div([html.Span("Data source", className="meta-label"), html.A(f"Open Targets {data_version}".strip(), href=source, target="_blank", rel="noreferrer")]),
        ], className="source-bar", **{"aria-label": "Snapshot data"}),
        ], className="hero"),
        html.Section([
            html.Div([
                html.Section([html.H3("Drug evidence"),
                    html.Div([html.Div([html.Label("Clinical stage", htmlFor="stage"), info_tip("stage", "clinical stage", "Highest recorded stage of a canonical drug in this disease. Approval reached may include withdrawn drugs.")], className="label-help"), dcc.Dropdown(id="stage", options=[{"label": label, "value": value} for label, value in (("Phase I or later", "phase1"), ("Phase II or later", "phase2"), ("Phase III or later", "phase3"), ("Approval reached", "approved"))], value="phase3", clearable=False)], className="control"),
                    control("Drug modality", dcc.Dropdown(id="modality", options=modality_options, value="all", clearable=False)),
                ], className="filter-group"),
                html.Section([html.H3("Expression criteria"),
                    control("Expression rule", dcc.Dropdown(id="method", options=[{"label": label, "value": value} for value, label in METHOD_LABELS.items()], value="fixed", clearable=False)),
                    # Dash 4.4.1 は max を省略すると増減時に NaN になるため、上限なしを明示する。
                    html.Div([html.Div([html.Label("Minimum CPM (strict >)", htmlFor="threshold"), info_tip("threshold", "minimum CPM", f"Every rule requires median CPM above this value. Blank uses {DEFAULT_EXPRESSION_THRESHOLD:g} CPM.")], className="label-help"), dcc.Input(id="threshold", type="number", min=0, max=None, step=.1, value=DEFAULT_EXPRESSION_THRESHOLD)], className="control"),
                    html.Div([html.Div([html.Label("CELLEX specificity (≥)", htmlFor="specificity"), info_tip("specificity", "CELLEX specificity", f"Only used with Fixed CPM + CELLEX specificity. Blank uses {DEFAULT_SPECIFICITY_THRESHOLD:g}.")], className="label-help"), dcc.Input(id="specificity", type="number", min=0, max=1, step=.05, value=DEFAULT_SPECIFICITY_THRESHOLD, disabled=True)], className="control"),
                ], className="filter-group"),
                html.Div([
                    html.H3("Comparison scope"),
                    html.Div([
                        disease_selector(snapshot, default_diseases),
                        cell_selector(snapshot, default_cells),
                    ], className="control-grid scope-controls"),
                ], className="filter-group scope-group"),
                html.Section([html.H3("Display"),
                    html.Fieldset([html.Legend(["View", info_tip("view", "heatmap view", "Targets are distinct genes meeting the rule. Drug forms mapped to the same active ingredient count once if any known target meets the rule.")]), dcc.RadioItems(id="heatmap-view", options=[{"label": "Distinct targets", "value": "target"}, {"label": "Canonical drugs", "value": "drug"}], value="target", inline=True)], className="control radio-control"),
                    html.Fieldset([html.Legend(["Measure", info_tip("measure", "measure", "Percent divides by known targets or canonical drugs with known targets for the current disease and drug filters. Each denominator stays fixed across cells; percent is not a share of cells.")]), dcc.RadioItems(id="measure", options=[{"label": "Count", "value": "count"}, {"label": "Percent", "value": "percent"}], value="count", inline=True)], className="control radio-control"),
                ], className="filter-group"),
            ], className="filter-groups"),
            html.Div([
                html.Button("Update", id="update-button", type="button", n_clicks=0),
                html.Span("Settings applied.", id="update-status", role="status"),
            ], className="update-actions"),
            dcc.Store(id="applied-parameters", data=applied),
        ], className="panel controls"),
        html.Section([html.Div([html.H2(["Cell-type comparison", info_tip("comparison", "cell-type comparison", "Heatmap cells show confirmed support without a ≥ mark; hover, details, and CSV values mark lower bounds. × can mean unavailable disease data, unresolved target or expression data, or no eligible percentage denominator. A zero count means no qualifying targets or drugs under the rule.")]), html.Div([html.Button("Download CSV", id="download-button", n_clicks=0), info_tip("csv", "CSV download", "CSV uses the settings last applied with Update, matching the displayed results. It includes summaries, expression values and source drug records."), dcc.Download(id="download")], className="download-actions")], className="section-heading"), html.Div(id="matrix-note", className="matrix-note", role="status"), html.Div([html.Article([html.H3("Distinct targets"), html.Div(dcc.Graph(id="target-heatmap", config={"displaylogo": False, "responsive": False}), className="graph-scroll")], id="target-heatmap-panel"), html.Article([html.H3("Canonical drugs"), html.Div(dcc.Graph(id="drug-heatmap", config={"displaylogo": False, "responsive": False}), className="graph-scroll")], id="drug-heatmap-panel", hidden=True)], className="heatmap-stack"), html.P("× No value · Hover for values · Click for details", className="matrix-note")], className="panel matrix-panel", id="comparison"),
        html.Section([html.Div([html.H2(["Selection details", info_tip("selection", "selection details", "The selectors resolve the disease–cell summary and evidence tables. The expression heatmap shows all healthy reference source cell types for targets in the selected disease and drug filters.")])], className="section-heading detail-heading"), html.Div([control("Disease", dcc.Dropdown(id="detail-disease", clearable=False)), control("Cell", dcc.Dropdown(id="detail-cell", clearable=False))], className="control-grid selectors"), html.Div(id="details", className="details")], className="panel", id="evidence"),
        html.Footer("Drug records: Open Targets. Healthy reference expression: Tabula Sapiens."),
    ], className="shell")


def register_callbacks(application: Dash, snapshot: dict) -> None:
    """schema 2 dashboard のコールバックを登録する。"""
    @application.callback(
        Output("applied-parameters", "data"), Input("update-button", "n_clicks"),
        *[State(item, "value") for item in PARAMETER_IDS], prevent_initial_call=True,
    )
    def apply_parameters(_clicks, *values):
        return dict(zip(PARAMETER_IDS, values))

    @application.callback(
        Output("update-status", "children"), Input("applied-parameters", "data"),
        *[Input(item, "value") for item in PARAMETER_IDS],
    )
    def parameter_status(applied, *values):
        return "Settings applied." if dict(zip(PARAMETER_IDS, values)) == applied else "Changes not applied. Click Update."

    families = [section for group in disease_catalog(snapshot) for family in group["families"] for section in _disease_checklist_sections(family)]
    family_ids = [section["id"] for section in families]

    @application.callback(
        Output("diseases", "value"), *[Output(item, "value") for item in family_ids],
        Input("diseases", "value"), *[Input(item, "value") for item in family_ids],
    )
    def sync_disease_selection(selected, *family_values):
        selected = set(selected or [])
        if ctx.triggered_id in family_ids:
            index = family_ids.index(ctx.triggered_id)
            members = {d["id"] for d in families[index]["diseases"]}
            selected = (selected - members) | (set(family_values[index] or []) & members)
        ordered = ordered_disease_ids(snapshot, selected)
        return [ordered, *[[d["id"] for d in family["diseases"] if d["id"] in ordered] for family in families]]

    @application.callback(Output("specificity", "disabled"), Input("method", "value"))
    def toggle_specificity(method):
        return method != "specificity"

    @application.callback(Output("target-heatmap-panel", "hidden"), Output("drug-heatmap-panel", "hidden"), Input("applied-parameters", "data"))
    def select_heatmap(applied):
        view = applied["heatmap-view"]
        return view == "drug", view != "drug"

    cell_sections = [section for group in _cell_selection_catalog(snapshot) for section in group["sections"]]
    cell_section_ids = [section["id"] for section in cell_sections]

    @application.callback(
        Output("cells", "value"), *[Output(item, "value") for item in cell_section_ids],
        Input("cells", "value"), *[Input(item, "value") for item in cell_section_ids],
    )
    def sync_cell_selection(selected, *section_values):
        selected = set(selected or [])
        if ctx.triggered_id in cell_section_ids:
            index = cell_section_ids.index(ctx.triggered_id)
            members = {cell["id"] for cell in cell_sections[index]["cells"]}
            selected = (selected - members) | (set(section_values[index] or []) & members)
        ordered = _ordered_cell_ids(snapshot, "mixed", selected)
        return [ordered, *[[cell["id"] for cell in section["cells"] if cell["id"] in ordered] for section in cell_sections]]

    @application.callback(Output("target-heatmap", "figure"), Output("drug-heatmap", "figure"), Output("matrix-note", "children"), Input("applied-parameters", "data"))
    def update_figures(applied):
        measure, modality, stage, method, threshold, specificity, disease_ids, cell_ids = (applied[key] for key in PARAMETER_IDS[:-1])
        minimum, specificity_value, errors = effective_filters(threshold, specificity)
        disease_ids = ordered_disease_ids(snapshot, disease_ids)
        cell_ids = _ordered_cell_ids(snapshot, "mixed", cell_ids)
        rows = visible_rows(snapshot, modality, minimum, disease_ids, cell_ids, stage=stage, method=method, specificity=specificity_value, level="mixed")
        target_figure, drug_figure = build_figure(rows, disease_ids or [], cell_ids or [], measure, "target"), build_figure(rows, disease_ids or [], cell_ids or [], measure, "drug")
        target_missing = sum(row[_measure_fields("target", measure)[0]] is None for row in rows); drug_missing = sum(row[_measure_fields("drug", measure)[0]] is None for row in rows)
        condition = f"median CPM > {minimum:g}"
        if method == "relative":
            condition += " and CPM ≥ the full-reference target median"
        elif method == "specificity":
            condition += f" and CELLEX specificity ≥ {specificity_value:g}"
        modality_label = "All modalities" if modality == "all" else next((label for label, value in atlas.DRUG_TYPE_MODALITIES if value == modality), modality)
        note = f"Applied: Clinical stage: {STAGE_LABELS[stage]}; Drug modality: {modality_label}; Expression rule: {METHOD_LABELS[method]} ({condition}); Cells: selected groups and source cells; Measure: {measure.capitalize()}. {len(rows)} disease–cell combinations. Heatmap entries without a value: targets {target_missing}, drugs {drug_missing}."
        status = html.Span(f"{len(rows)} disease–cell combinations")
        error_note = html.Span(" ".join(errors), className="filter-errors", role="alert") if errors else None
        return target_figure, drug_figure, html.Div([status, info_tip("applied", "applied filters", note), error_note], className="matrix-status")

    @application.callback(Output("detail-disease", "options"), Output("detail-disease", "value"), Output("detail-cell", "options"), Output("detail-cell", "value"), Input("applied-parameters", "data"), Input("target-heatmap", "clickData"), Input("drug-heatmap", "clickData"), State("detail-disease", "value"), State("detail-cell", "value"))
    def update_detail_selectors(applied, target_click, drug_click, current_disease, current_cell):
        disease_ids, cell_ids = applied["diseases"], applied["cells"]
        disease_ids, cell_ids = ordered_disease_ids(snapshot, disease_ids), _ordered_cell_ids(snapshot, "mixed", cell_ids); disease_names = {row["id"]: row["name"] for row in snapshot["diseases"]}; cell_names = dict(cell_catalog(snapshot, "mixed"))
        triggered = ctx.triggered_id
        click_data = target_click if triggered == "target-heatmap" else drug_click if triggered == "drug-heatmap" else None
        visible = [{"disease_id": disease, "cell_id": cell} for disease in disease_ids for cell in cell_ids]
        selected = resolve_selection(triggered, click_data, current_disease, current_cell, visible)
        if selected:
            current_disease, current_cell = selected
        return ([{"label": disease_names[item], "value": item} for item in disease_ids if item in disease_names], current_disease if current_disease in disease_ids else next(iter(disease_ids), None), [{"label": cell_names[item], "value": item} for item in cell_ids if item in cell_names], current_cell if current_cell in cell_ids else next(iter(cell_ids), None))

    @application.callback(Output("details", "children"), Input("detail-disease", "value"), Input("detail-cell", "value"), Input("applied-parameters", "data"))
    def update_details(detail_disease, detail_cell, applied):
        _measure, modality, stage, method, threshold, specificity, disease_ids, cell_ids = (applied[key] for key in PARAMETER_IDS[:-1])
        minimum, specificity_value, _ = effective_filters(threshold, specificity)
        rows = visible_rows(snapshot, modality, minimum, disease_ids, cell_ids, stage=stage, method=method, specificity=specificity_value, level="mixed")
        selection = resolve_selection("detail-cell", None, detail_disease, detail_cell, rows)
        return detail_panel(rows, selection, snapshot, modality=modality, stage=stage, threshold=minimum, method=method, specificity=specificity_value)

    @application.callback(Output("source-records-page", "children"), Input("source-toggle", "n_clicks"), Input("source-page", "value"), Input("source-context", "data"))
    def show_source_records(clicks, page, context):
        # Summary のクリックは、キーボード操作でも開閉ごとに一度発生する。
        if not clicks or clicks % 2 == 0 or not context:
            return None
        records = evidence_rows(snapshot, context, **{key: context[key] for key in ("modality", "stage", "threshold", "method", "specificity")})
        pages = max(1, math.ceil(len(records) / SOURCE_PAGE_SIZE))
        page = min(max(page, 1), pages) if isinstance(page, int) and not isinstance(page, bool) else 1
        start = (page - 1) * SOURCE_PAGE_SIZE
        end = min(start + SOURCE_PAGE_SIZE, len(records))
        return html.Div([
            html.P(f"Rows {start + 1 if records else 0}–{end} of {len(records)}", role="status", className="matrix-note"),
            _evidence_table(records[start:end]),
        ])

    @application.callback(Output("download", "data"), Input("download-button", "n_clicks"), State("applied-parameters", "data"), prevent_initial_call=True)
    def download_csv(_clicks, applied):
        measure, modality, stage, method, threshold, specificity, disease_ids, cell_ids = (applied[key] for key in PARAMETER_IDS[:-1])
        minimum, specificity_value, _ = effective_filters(threshold, specificity)
        rows = visible_rows(snapshot, modality, minimum, disease_ids, cell_ids, stage=stage, method=method, specificity=specificity_value, level="mixed")
        content = "\ufeff" + atlas.to_csv(export_rows(rows, measure, modality_filter=modality, stage_filter=stage, method=method, expression_threshold=minimum, specificity_threshold=specificity_value, level="mixed", snapshot=snapshot))
        return {"content": content, "filename": "autoimmune-drug-cell-matrix.csv", "type": "text/csv;charset=utf-8"}


def create_app(snapshot: dict | None, error: Exception | None = None) -> Dash:
    """保存済みデータまたは明示的な fixture から Dash アプリを作る。"""
    # 元記録の操作部は、選択した詳細を描画するときに追加する。
    application = Dash(__name__, title="Autoimmune Target Expression Atlas", suppress_callback_exceptions=True)
    application.layout = unavailable_layout(error) if snapshot is None else dashboard_layout(snapshot)
    if snapshot is not None:
        register_callbacks(application, snapshot)
    return application


try:
    SNAPSHOT = load_snapshot(SNAPSHOT_PATH); LOAD_ERROR = None
except (OSError, json.JSONDecodeError, ValueError, TypeError, KeyError) as error:
    SNAPSHOT = None; LOAD_ERROR = error

app = create_app(SNAPSHOT, LOAD_ERROR)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Autoimmune Atlas"); parser.add_argument("--dev", action="store_true", help="Automatically reload when code changes"); args = parser.parse_args()
    app.run(host="127.0.0.1", port=8050, debug=args.dev)
