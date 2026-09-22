"""自己免疫疾患の薬剤標的と細胞型を閲覧する Dash アプリ。"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import plotly.graph_objects as go
from dash import Dash, Input, Output, State, ctx, dcc, html

import atlas

BASE_DIR = Path(__file__).parent
SNAPSHOT_PATH = BASE_DIR / "data" / "snapshot.json"
DEFAULT_EXPRESSION_THRESHOLD = 0.5
DEFAULT_SPECIFICITY_THRESHOLD = 0.75
FILTER_IDS = {"measure", "modality", "stage", "method", "threshold", "specificity", "level", "diseases", "cells"}
STAGE_LABELS = {"phase1": "Phase I or later", "phase2": "Phase II or later", "phase3": "Phase III or later", "approved": "Approval reached"}
METHOD_LABELS = {"fixed": "Fixed CPM", "relative": "Target-relative median", "specificity": "CELLEX specificity"}


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
            "This snapshot uses schema 1 and cannot support clinical-stage or specificity filters. "
            "Refresh it with: pixi run --as-is python fetch_data.py"
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


def _default_cells(snapshot: dict, level: str) -> list[str]:
    return choose_defaults(
        cell_catalog(snapshot, level),
        ("b cell", "t cell", "monocyte", "macrophage", "dendritic", "natural killer", "neutrophil", "fibroblast", "endothelial", "epithelial", "plasma cell"),
        16,
    )


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
    selected_diseases, selected_cells = set(disease_ids or []), set(cell_ids or [])
    return [
        row
        for row in atlas.summarize(snapshot, modality, minimum, stage=stage, method=method, specificity_threshold=specificity_value, level=level)
        if row["disease_id"] in selected_diseases and row["cell_id"] in selected_cells
    ]


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
    value = _display_value(row, measure, kind) or "Unknown / missing"
    label = "Drugs" if kind == "drug" else "Targets"
    assessment = "≥ is a lower bound; unresolved evidence may increase this value" if _is_lower_bound(row, measure, kind) else "Insufficient data to determine a value" if row.get(value_key) is None else f"All {label.lower()} in this denominator are assessed"
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
            text_row.append("" if row is None else _display_value(row, measure, kind))
            hover_row.append("No data" if row is None else _hover_text(row, measure, kind))
            custom_row.append([disease_id, cell_id])
            if value is None:
                missing_x.append(disease_names.get(disease_id, disease_id)); missing_y.append(cell_names.get(cell_id, cell_id)); missing_custom.append([disease_id, cell_id]); missing_hover.append("No data" if row is None else _hover_text(row, measure, kind))
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
        figure.add_trace(go.Scatter(x=missing_x, y=missing_y, customdata=missing_custom, hovertext=missing_hover, mode="markers", marker={"symbol": "x", "size": 10, "color": "#8a99a6"}, name="Unknown / missing", hovertemplate="%{hovertext}<extra></extra>"))
    if not disease_ids or not cell_ids:
        figure.add_annotation(text="Select diseases and cell types", showarrow=False)
    figure.update_layout(template="plotly_white", font={"family": "Arial, sans-serif", "size": 12, "color": "#263238"}, width=max(760, min(2800, 260 + 125 * len(disease_ids))), height=max(420, min(1400, 180 + 28 * len(cell_ids))), margin={"l": 180, "r": 40, "t": 20, "b": 150}, legend={"orientation": "h", "y": 1.08, "x": 0}, hoverlabel={"align": "left"})
    figure.update_xaxes(tickangle=-35, side="bottom", title="Disease")
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
    output = []
    for record in records:
        base = {**record, "group_cell_id": row["cell_id"], "group_cell": row["cell"], "references_json": json.dumps(record.get("references") or [], ensure_ascii=False, sort_keys=True)}
        if not record.get("target_id"):
            output.append({**base, "evidence_cell_id": "", "evidence_cell": "", "cpm": None, "specificity_score": None, "target_median": None, "support_state": "unmapped", "contributing": False})
            continue
        for member in row["member_cell_ids"]:
            item = metadata.get((record["target_id"], member))
            state = atlas.expression_state(item, threshold, method, specificity)
            output.append({
                **base, "evidence_cell_id": member, "evidence_cell": item.get("cell", cell_names.get(member, member)) if item else cell_names.get(member, member),
                "cpm": item.get("median") if item else None, "specificity_score": item.get("specificity_score") if item else None,
                "target_median": item.get("target_median") if item else None,
                "support_state": "positive" if state is True else "negative" if state is False else "unknown", "contributing": state is True,
            })
    return output


def expression_figure(snapshot, row, records, metadata=None) -> go.Figure:
    """選択した元細胞について、全既知標的の連続発現量を示す。"""
    metadata = metadata if metadata is not None else atlas.expression_metadata(snapshot)
    target_names = {}
    for record in records:
        if record.get("target_id"):
            target_names.setdefault(record["target_id"], record.get("target") or record["target_id"])
    members = list(row["member_cell_ids"])
    member_names = {member: next((metadata[target, member]["cell"] for target in target_names if (target, member) in metadata), member) for member in members}
    targets = sorted(target_names, key=lambda item: (target_names[item].casefold(), item))
    z, hover, missing_x, missing_y = [], [], [], []
    for target_id in targets:
        z_row, hover_row = [], []
        label = f"{target_names[target_id]} ({target_id})"
        for member in members:
            item = metadata.get((target_id, member)); cpm = item.get("median") if item else None
            z_row.append(math.log2(1 + cpm) if cpm is not None else None)
            hover_row.append("<br>".join((
                f"<b>{label}</b>", f"Cell: {member_names[member]} ({member})", f"Raw CPM: {cpm:g}" if cpm is not None else "Raw CPM: missing",
                f"CELLEX specificity: {item.get('specificity_score'):g}" if item and item.get("specificity_score") is not None else "CELLEX specificity: missing",
                f"Target-relative median: {item.get('target_median'):g}" if item and item.get("target_median") is not None else "Target-relative median: missing",
            )))
            if cpm is None:
                missing_x.append(member_names[member]); missing_y.append(label)
        z.append(z_row); hover.append(hover_row)
    figure = go.Figure(go.Heatmap(x=[member_names[m] for m in members], y=[f"{target_names[t]} ({t})" for t in targets], z=z, zmin=0, colorscale="Blues", colorbar={"title": "log2(1 + CPM)"}, hovertext=hover, hovertemplate="%{hovertext}<extra></extra>", xgap=2, ygap=2))
    if missing_x:
        figure.add_trace(go.Scatter(x=missing_x, y=missing_y, mode="markers", marker={"symbol": "x", "size": 9, "color": "#8a99a6"}, name="Missing expression", hovertemplate="Expression missing<extra></extra>"))
    if not targets:
        figure.add_annotation(text="No known targets in the selected scope", showarrow=False)
    figure.update_layout(template="plotly_white", height=max(360, min(1200, 170 + 25 * len(targets))), margin={"l": 190, "r": 40, "t": 20, "b": 110}, font={"family": "Arial, sans-serif", "size": 11, "color": "#263238"}, hoverlabel={"align": "left"})
    figure.update_xaxes(tickangle=-30, title="Source cell type"); figure.update_yaxes(autorange="reversed", title="Known target")
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
    headers = ("State", "Canonical drug", "Original drug", "Modality", "Original / canonical stage", "Target", "Action / mechanism", "Source cell", "CPM", "Specificity", "Target median", "Sources")
    body = []
    for record in records:
        actions, mechanism = ", ".join(record.get("action_types") or []), record.get("mechanism") or "—"
        body.append(html.Tr([
            html.Td(record["support_state"]), html.Td(f"{record.get('canonical_drug') or '—'} ({record.get('canonical_drug_id') or '—'})"),
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
        return html.Div("Select a disease and cell type using either heatmap or the selectors above.", className="empty-note")
    row = next((item for item in rows if (item["disease_id"], item["cell_id"]) == selection), None)
    if row is None:
        return html.Div("The previous selection is outside the current filters. Select a visible cell.", className="empty-note")
    summary = html.Div([html.H3(f"{row['disease']} × {row['cell']}"), html.P(
        f"Targets: {_display_value(row, 'count') or 'NA'}; {_display_value(row, 'percent') or 'NA'} of {row['denominator']} known targets ({row['unknown']} unknown). "
        f"Drugs: {_display_value(row, 'count', 'drug') or 'NA'}; {_display_value(row, 'percent', 'drug') or 'NA'} of {row['drug_denominator']} mapped canonical drugs ({row['unknown_drugs']} unknown). "
        f"Coverage: {row['mapped_drugs']} mapped / {row['total_drugs']} all canonical drugs; {row['unmapped_drugs']} without mapped targets.", className="detail-summary")])
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
        item = positive_pairs.setdefault(key, {"drug": record.get("canonical_drug") or record.get("drug"), "target": record.get("target") or record.get("target_id"), "cells": []})
        cell = record.get("evidence_cell") or record.get("evidence_cell_id")
        if cell not in item["cells"]:
            item["cells"].append(cell)
    positive_list = html.Ul([html.Li(f"{item['drug']} → {item['target']} in {', '.join(item['cells'])}") for item in positive_pairs.values()], className="evidence-list") if positive_pairs else html.P("No positive evidence under the applied rule.", className="empty-note")
    return html.Div([
        summary, html.H3("Continuous expression for all known targets"),
        html.P("Color is log2(1 + CPM); hover shows raw CPM. Group rows display member cell types separately and never combine their CPM values.", className="matrix-note"),
        html.Div(dcc.Graph(figure=expression_figure(snapshot, row, records, metadata), config={"displaylogo": False}), className="graph-scroll expression-graph"),
        html.H3("Positive evidence"), positive_list, html.Details([html.Summary("All filtered source records"), html.P("Negative, unknown, and unmapped records remain here so the aggregate can be traced.", className="matrix-note"), _evidence_table(records)], open=False, className="source-records"),
    ])


def export_rows(rows, measure, *, modality_filter="all", stage_filter="phase3", method="fixed", expression_threshold=DEFAULT_EXPRESSION_THRESHOLD, specificity_threshold=DEFAULT_SPECIFICITY_THRESHOLD, level="group", snapshot=None):
    """集計、元薬剤記録、発現根拠を重複させず一つの CSV に並べる。"""
    snapshot = snapshot or {}
    base = {
        "row_type": "", "disease_id": "", "disease": "", "group_cell_id": "", "group_cell": "", "evidence_cell_id": "", "evidence_cell": "",
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
    targets_by_disease = {}
    for record in records:
        if record.get("target_id"):
            targets_by_disease.setdefault(record["disease_id"], {})[record["target_id"]] = record.get("target") or record["target_id"]
    for row in rows:
        for target_id, target in targets_by_disease.get(row["disease_id"], {}).items():
            for member in row["member_cell_ids"]:
                item = metadata.get((target_id, member))
                state = atlas.expression_state(item, expression_threshold, method, specificity_threshold)
                output.append({
                    **base, "row_type": "expression_evidence", "disease_id": row["disease_id"], "disease": row["disease"],
                    "group_cell_id": row["cell_id"], "group_cell": row["cell"], "evidence_cell_id": member,
                    "evidence_cell": item.get("cell", cell_names.get(member, member)) if item else cell_names.get(member, member),
                    "support_state": "positive" if state is True else "negative" if state is False else "unknown", "contributing": state is True,
                    "target_id": target_id, "target": target, "raw_cpm": item.get("median", "") if item else "",
                    "specificity_score": item.get("specificity_score", "") if item else "", "target_relative_median": item.get("target_median", "") if item else "",
                })
    return output


def unavailable_layout(error=None) -> html.Main:
    """データ未取得または再取得が必要な状態を説明する。"""
    title = "No data loaded yet" if error is None else "Data refresh required"
    detail = "pixi run --as-is python fetch_data.py" if error is None else str(error)
    return html.Main([html.Section([html.P("AUTOIMMUNE DRUG–CELL ATLAS", className="eyebrow"), html.H1(title), html.P("Create the schema 2 snapshot to display stage, specificity, and cell-level results."), html.Code(detail)], className="unavailable-card")], className="shell unavailable")


def dashboard_layout(snapshot) -> html.Main:
    """schema 2 のデータから dashboard の初期画面を作る。"""
    diseases = [(row["id"], row["name"]) for row in snapshot["diseases"]]
    default_diseases = choose_defaults(diseases, ("rheumatoid arthritis", "systemic lupus erythematosus", "multiple sclerosis", "systemic sclerosis", "Sjogren syndrome", "myasthenia gravis", "psoriatic arthritis", "type 1 diabetes mellitus"), 8, {row["disease_id"] for row in snapshot["records"]})
    group_cells, default_cells = cell_catalog(snapshot, "group"), _default_cells(snapshot, "group")
    ready = sum(row["status"] == "ready" for row in snapshot["diseases"]); data_version = format_data_version(snapshot.get("data_version")); source = "https://platform.opentargets.org/"
    modality_options = [{"label": "All", "value": "all"}] + [{"label": label, "value": value} for label, value in atlas.DRUG_TYPE_MODALITIES]
    control = lambda label, component: html.Div([html.Label(label, htmlFor=component.id), component], className="control")
    return html.Main([
        html.Nav([html.A("Autoimmune Atlas", href="#", className="app-brand"), html.Div([html.A("Comparison", href="#comparison"), html.A("Evidence", href="#evidence"), html.A("Open Targets ↗", href=source, target="_blank", rel="noreferrer")], className="app-nav")], className="app-bar", **{"aria-label": "Main navigation"}),
        html.Header([html.Div([html.P("AUTOIMMUNE DISEASE / DRUG TARGETS", className="eyebrow"), html.H1("Drug target support by cell type"), html.P("Compare target breadth and canonical-drug coverage under one shared expression rule.", className="lede")]), html.Div([html.Div([html.Strong(str(len(snapshot["diseases"]))), html.Span("Disease terms")], className="stat"), html.Div([html.Strong(str(len(snapshot["records"]))), html.Span("Source records")], className="stat"), html.Div([html.Strong(str(len(snapshot["expression"]))), html.Span("Targets queried")], className="stat")], className="stats")], className="hero"),
        html.Section([html.Div([html.Span("Data status", className="meta-label"), html.Strong(f"Ready {ready} / unavailable {len(snapshot['diseases']) - ready}")]), html.Div([html.Span("Retrieved at", className="meta-label"), html.Strong(snapshot["retrieved_at"])]), html.Div([html.Span("Expression reference", className="meta-label"), html.Strong("Tabula Sapiens")]), html.Div([html.Span("Data source", className="meta-label"), html.A(f"Open Targets {data_version}".strip(), href=source, target="_blank", rel="noreferrer")])], className="source-bar"),
        html.Section([
            html.Div([
                control("Measure", dcc.RadioItems(id="measure", options=[{"label": "Count", "value": "count"}, {"label": "Percent", "value": "percent"}], value="count", inline=True)),
                control("Clinical stage", dcc.Dropdown(id="stage", options=[{"label": label, "value": value} for label, value in (("Phase I or later", "phase1"), ("Phase II or later", "phase2"), ("Phase III or later", "phase3"), ("Approval reached", "approved"))], value="phase3", clearable=False)),
                control("Expression rule", dcc.Dropdown(id="method", options=[{"label": label, "value": value} for label, value in (("Fixed CPM", "fixed"), ("Target-relative median", "relative"), ("CELLEX specificity", "specificity"))], value="fixed", clearable=False)),
                control("Minimum CPM (strict >)", dcc.Input(id="threshold", type="number", min=0, step=.1, value=DEFAULT_EXPRESSION_THRESHOLD)),
                control("Specificity (≥)", dcc.Input(id="specificity", type="number", min=0, max=1, step=.05, value=DEFAULT_SPECIFICITY_THRESHOLD)),
                control("Cell level", dcc.RadioItems(id="level", options=[{"label": "Groups", "value": "group"}, {"label": "Source cells", "value": "cell"}], value="group", inline=True)),
                control("Drug modality", dcc.Dropdown(id="modality", options=modality_options, value="all", clearable=False)),
            ], className="control-grid compact"),
            html.Details([html.Summary("Select diseases and cells"), html.Div([control("Diseases", dcc.Dropdown(id="diseases", options=[{"label": name, "value": item_id} for item_id, name in diseases], value=default_diseases, multi=True, searchable=True)), control("Cells", dcc.Dropdown(id="cells", options=[{"label": name, "value": item_id} for item_id, name in group_cells], value=default_cells, multi=True, searchable=True))], className="control-grid selectors")], className="filter-details"),
            html.Details([html.Summary("How to read this view"), html.P("Healthy-atlas expression is not disease-specific expression or evidence of efficacy. Historical approvals may be withdrawn, and disease records can include symptom or comorbidity treatment."), html.P("Groups are unions of source cell types. Their values can increase with the number and granularity of member subtypes; member CPM values are never added or averaged."), html.P("≥ marks a lower bound. × marks an unknown value. Percentages use their own known-item denominator and can be 0 even when a count remains unknown because of an unmapped drug.")], className="methods-note"),
        ], className="panel controls"),
        html.Section([html.Div([html.H2("Shared-condition comparison"), html.Button("Download CSV", id="download-button", n_clicks=0), dcc.Download(id="download")], className="section-heading"), html.P("CSV rows are labeled summary, source_record, or expression_evidence. Join source and expression rows by disease_id + target_id; group_cell_id and evidence_cell_id identify the displayed group and measured source cell.", className="matrix-note"), html.P(id="matrix-note", className="matrix-note", role="status"), html.Div([html.Article([html.H3("Distinct targets"), html.Div(dcc.Graph(id="target-heatmap", config={"displaylogo": False, "responsive": False}), className="graph-scroll")]), html.Article([html.H3("Canonical drugs"), html.Div(dcc.Graph(id="drug-heatmap", config={"displaylogo": False, "responsive": False}), className="graph-scroll")])], className="heatmap-stack"), html.P("≥ Lower bound · × Unknown / missing · Click either heatmap to inspect evidence", className="matrix-note")], className="panel matrix-panel", id="comparison"),
        html.Section([html.Div([html.H2("Selection details"), html.P("The selectors resolve the same disease–cell evidence as a heatmap click.")], className="section-heading detail-heading"), html.Div([control("Disease", dcc.Dropdown(id="detail-disease", clearable=False)), control("Cell", dcc.Dropdown(id="detail-cell", clearable=False))], className="control-grid selectors"), html.Div(id="details", className="details")], className="panel", id="evidence"),
        html.Footer("Target and canonical-drug summaries use separate denominators and unknown counts."),
    ], className="shell")


def register_callbacks(application: Dash, snapshot: dict) -> None:
    """schema 2 dashboard のコールバックを登録する。"""
    @application.callback(Output("cells", "options"), Output("cells", "value"), Input("level", "value"), State("cells", "value"))
    def update_cell_options(level, current_cells):
        catalog = cell_catalog(snapshot, level)
        groups = atlas.cell_catalog(snapshot, "group")
        group_members = {item["id"]: item["members"] for item in groups}
        member_group = {member: item["id"] for item in groups for member in item["members"]}
        selected = []
        for item in current_cells or []:
            mapped = group_members.get(item, [item]) if level == "cell" else [item if item in group_members else member_group.get(item)]
            selected.extend(value for value in mapped if value and value not in selected)
        available = {item_id for item_id, _ in catalog}
        selected = [item for item in selected if item in available]
        return [{"label": name, "value": item_id} for item_id, name in catalog], selected or _default_cells(snapshot, level)

    @application.callback(Output("target-heatmap", "figure"), Output("drug-heatmap", "figure"), Output("matrix-note", "children"), Input("measure", "value"), Input("modality", "value"), Input("stage", "value"), Input("method", "value"), Input("threshold", "value"), Input("specificity", "value"), Input("level", "value"), Input("diseases", "value"), Input("cells", "value"))
    def update_figures(measure, modality, stage, method, threshold, specificity, level, disease_ids, cell_ids):
        minimum, specificity_value, errors = effective_filters(threshold, specificity)
        rows = visible_rows(snapshot, modality, minimum, disease_ids, cell_ids, stage=stage, method=method, specificity=specificity_value, level=level)
        target_figure, drug_figure = build_figure(rows, disease_ids or [], cell_ids or [], measure, "target"), build_figure(rows, disease_ids or [], cell_ids or [], measure, "drug")
        target_missing = sum(row[_measure_fields("target", measure)[0]] is None for row in rows); drug_missing = sum(row[_measure_fields("drug", measure)[0]] is None for row in rows)
        condition = f"median CPM > {minimum:g}"
        if method == "relative":
            condition += " and CPM ≥ the full-reference target median"
        elif method == "specificity":
            condition += f" and CELLEX specificity ≥ {specificity_value:g}"
        note = f"Applied: {STAGE_LABELS[stage]}; {METHOD_LABELS[method]} ({condition}); {level}; {measure}; {modality}. {len(rows)} combinations; target NA {target_missing}; drug NA {drug_missing}."
        return target_figure, drug_figure, note + ((" " + " ".join(errors)) if errors else "")

    @application.callback(Output("detail-disease", "options"), Output("detail-disease", "value"), Output("detail-cell", "options"), Output("detail-cell", "value"), Input("diseases", "value"), Input("cells", "value"), Input("level", "value"), State("detail-disease", "value"), State("detail-cell", "value"))
    def update_detail_selectors(disease_ids, cell_ids, level, current_disease, current_cell):
        disease_ids, cell_ids = disease_ids or [], cell_ids or []; disease_names = {row["id"]: row["name"] for row in snapshot["diseases"]}; cell_names = dict(cell_catalog(snapshot, level))
        return ([{"label": disease_names[item], "value": item} for item in disease_ids if item in disease_names], current_disease if current_disease in disease_ids else next(iter(disease_ids), None), [{"label": cell_names[item], "value": item} for item in cell_ids if item in cell_names], current_cell if current_cell in cell_ids else next(iter(cell_ids), None))

    @application.callback(Output("details", "children"), Input("target-heatmap", "clickData"), Input("drug-heatmap", "clickData"), Input("detail-disease", "value"), Input("detail-cell", "value"), Input("measure", "value"), Input("modality", "value"), Input("stage", "value"), Input("method", "value"), Input("threshold", "value"), Input("specificity", "value"), Input("level", "value"), Input("diseases", "value"), Input("cells", "value"))
    def update_details(target_click, drug_click, detail_disease, detail_cell, _measure, modality, stage, method, threshold, specificity, level, disease_ids, cell_ids):
        minimum, specificity_value, _ = effective_filters(threshold, specificity)
        rows = visible_rows(snapshot, modality, minimum, disease_ids, cell_ids, stage=stage, method=method, specificity=specificity_value, level=level); triggered = ctx.triggered_id
        if triggered in FILTER_IDS:
            selection = None
        else:
            click_data = target_click if triggered == "target-heatmap" else drug_click if triggered == "drug-heatmap" else None
            selection = resolve_selection(triggered, click_data, detail_disease, detail_cell, rows)
        return detail_panel(rows, selection, snapshot, modality=modality, stage=stage, threshold=minimum, method=method, specificity=specificity_value)

    @application.callback(Output("download", "data"), Input("download-button", "n_clicks"), State("measure", "value"), State("modality", "value"), State("stage", "value"), State("method", "value"), State("threshold", "value"), State("specificity", "value"), State("level", "value"), State("diseases", "value"), State("cells", "value"), prevent_initial_call=True)
    def download_csv(_clicks, measure, modality, stage, method, threshold, specificity, level, disease_ids, cell_ids):
        minimum, specificity_value, _ = effective_filters(threshold, specificity)
        rows = visible_rows(snapshot, modality, minimum, disease_ids, cell_ids, stage=stage, method=method, specificity=specificity_value, level=level)
        content = "\ufeff" + atlas.to_csv(export_rows(rows, measure, modality_filter=modality, stage_filter=stage, method=method, expression_threshold=minimum, specificity_threshold=specificity_value, level=level, snapshot=snapshot))
        return {"content": content, "filename": "autoimmune-drug-cell-matrix.csv", "type": "text/csv;charset=utf-8"}


def create_app(snapshot: dict | None, error: Exception | None = None) -> Dash:
    """保存済みデータまたは明示的な fixture から Dash アプリを作る。"""
    application = Dash(__name__, title="Autoimmune Target Expression Atlas")
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
