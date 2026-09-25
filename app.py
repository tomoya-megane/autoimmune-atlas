"""自己免疫疾患の薬剤標的と細胞型を閲覧する Dash アプリ。"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timedelta, timezone
from functools import lru_cache
from pathlib import Path
from statistics import mean, pstdev
from textwrap import wrap

import plotly.graph_objects as go
from dash import ALL, Dash, Input, Output, State, ctx, dcc, html, no_update
from scipy.cluster.hierarchy import leaves_list, linkage

import atlas
from disease_catalog import disease_catalog, ordered_disease_ids

BASE_DIR = Path(__file__).parent
SNAPSHOT_PATH = BASE_DIR / "data" / "snapshot.json"
DEFAULT_EXPRESSION_THRESHOLD = 0.5
DEFAULT_SPECIFICITY_THRESHOLD = 0.5
SOURCE_PAGE_SIZE = 10
PARAMETER_IDS = (
    "measure",
    "modality",
    "stage",
    "method",
    "threshold",
    "specificity",
    "diseases",
    "heatmap-view",
)
STAGE_LABELS = {
    "phase1": "Phase I or later",
    "phase2": "Phase II or later",
    "phase3": "Phase III or later",
    "approved": "Approval reached",
}
METHOD_LABELS = {
    "fixed": "Fixed CPM",
    "relative": "Fixed CPM + Target-relative median",
    "specificity": "Fixed CPM + CELLEX specificity",
}


def info_tip(key: str, label: str, description: str):
    """操作対象の近くに、キーボードでも読める説明を置く。"""
    tip_id = f"{key}-tip"
    return html.Span(
        [
            html.Button(
                "ⓘ",
                type="button",
                className="info-button",
                **{"aria-label": f"About {label}", "aria-describedby": tip_id},
            ),
            html.Span(description, id=tip_id, role="tooltip", className="info-content"),
        ],
        className="info-tip",
    )


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
        raise ValueError(
            "Missing required fields in snapshot.json: " + ", ".join(sorted(missing))
        )
    if not isinstance(snapshot["diseases"], list) or not isinstance(
        snapshot["records"], list
    ):
        raise ValueError("diseases and records must be arrays")
    if not isinstance(snapshot["expression"], dict):
        raise ValueError("expression must be an object")
    return snapshot


def cell_catalog(snapshot: dict, level: str = "group") -> list[tuple[str, str]]:
    """公開関数が返す細胞カタログを選択欄向けに整える。"""
    return [(row["id"], row["name"]) for row in atlas.cell_catalog(snapshot, level)]


def _ordered_cell_ids(snapshot: dict, level: str, cell_ids) -> list[str]:
    selected = set(cell_ids or [])
    return [
        cell_id for cell_id, _ in cell_catalog(snapshot, level) if cell_id in selected
    ]


def _disease_tree(family):
    """ファミリー内の親子関係を、各語が 1 つの親の下にだけ現れる木にする。

    複数の親を持つ語は、より深い親の下に置き、同じ深さなら名前順で決める。
    親子が循環して根から届かない語は、最上位に置いて失わない。
    戻り値は (最上位の語の一覧, 語 ID → 子の一覧)。
    """
    members = {d["id"]: d for d in family["diseases"]}
    parents_of = {
        did: [p for p in d.get("parent_ids", []) if p in members and p != did]
        for did, d in members.items()
    }
    depth = {}

    def node_depth(did, trail=()):
        if did not in depth:
            parents = [p for p in parents_of[did] if p not in trail]
            depth[did] = (
                1 + max(node_depth(p, trail + (did,)) for p in parents)
                if parents
                else 0
            )
        return depth[did]

    def by_name(did):
        return members[did]["name"].casefold(), did

    parent = {}
    for did in sorted(members, key=by_name):
        if parents_of[did]:
            parent[did] = min(
                parents_of[did], key=lambda p: (-node_depth(p), *by_name(p))
            )

    def reaches_top(did, trail=()):
        return did not in parent or (
            did not in trail and reaches_top(parent[did], trail + (did,))
        )

    children = {did: [] for did in members}
    top = []
    for did in sorted(members, key=by_name):
        if did in parent and reaches_top(did):
            children[parent[did]].append(did)
        else:
            parent.pop(did, None)
            top.append(did)
    top.sort(key=lambda did: (did != family["id"], *by_name(did)))
    return top, children


def _disease_checklist_sections(family):
    """疾患本体と各階層の選択欄を、描画と同期で同じ範囲に分ける。

    子を持つ語は、ファミリーと同じ形にする。
    本体だけの section（kind が self）と、子のうち葉だけを入れた section（kind が
    children）を持ち、子を持つ子は同じ形で入れ子になる。
    先頭の section はファミリーの本体と、親を持たない葉である。
    """
    top, children = _disease_tree(family)
    members = {d["id"]: d for d in family["diseases"]}
    root = family["id"] if family["id"] in members else None

    def section_id(prefix, node):
        suffix = "" if node == root else f"-{node}"
        return f"{prefix}-{family['id']}{suffix}"

    sections = [
        {
            "id": f"disease-family-{family['id']}",
            "diseases": [members[d] for d in top if d == root or not children[d]],
            "node": None,
            "kind": "self",
        }
    ]
    queue = [d for d in top if children[d]]
    while queue:
        node = queue.pop(0)
        if node != root:
            sections.append(
                {
                    "id": section_id("disease-family", node),
                    "diseases": [members[node]],
                    "node": node,
                    "kind": "self",
                }
            )
        leaves = [d for d in children[node] if not children[d]]
        if leaves:
            sections.append(
                {
                    "id": section_id("disease-details", node),
                    "diseases": [members[d] for d in leaves],
                    "node": node,
                    "kind": "children",
                }
            )
        queue.extend(d for d in children[node] if children[d])
    return sections


def disease_selector(snapshot, selected):
    """検索欄と、群から開けるチェック欄を同じ選択へ結び付ける。"""
    catalog = disease_catalog(snapshot)
    options, groups = [], []
    for group in catalog:
        families = []
        for family in group["families"]:
            choices = [
                {"label": d["name"], "value": d["id"]} for d in family["diseases"]
            ]
            options.extend(choices)
            sections = _disease_checklist_sections(family)
            top, children = _disease_tree(family)
            root = (
                family["id"]
                if family["id"] in {d["id"] for d in family["diseases"]}
                else None
            )
            by_kind = {(s["kind"], s["node"]): s for s in sections}
            names = {d["id"]: d["name"] for d in family["diseases"]}

            def descendant_count(node):
                return sum(1 + descendant_count(c) for c in children[node])

            def checklist(section):
                return dcc.Checklist(
                    id=section["id"],
                    options=[
                        {"label": d["name"], "value": d["id"]}
                        for d in section["diseases"]
                    ],
                    value=[d["id"] for d in section["diseases"] if d["id"] in selected],
                    className="disease-checklist",
                )

            def render_details(node):
                """子を持つ語の details。葉のチェック欄と、子を持つ子の入れ子を並べる。"""
                inner = []
                if ("children", node) in by_kind:
                    inner.append(checklist(by_kind[("children", node)]))
                inner.extend(render_node(c) for c in children[node] if children[c])
                return html.Details(
                    [
                        html.Summary(
                            f"{names[node]} details ({descendant_count(node)} terms)"
                        ),
                        html.Div(inner, className="disease-families"),
                    ]
                )

            def render_node(node):
                """本体のチェック欄と details を、ファミリーと同じ形で包む。"""
                return html.Details(
                    [
                        html.Summary(
                            f"{names[node]} ({1 + descendant_count(node)} terms)"
                        ),
                        html.Div(
                            [checklist(by_kind[("self", node)]), render_details(node)],
                            className="disease-families",
                        ),
                    ]
                )

            top_items = [checklist(sections[0])]
            if root and children[root]:
                top_items.append(render_details(root))
            top_items.extend(render_node(d) for d in top if children[d] and d != root)
            contents = (
                top_items[0]
                if len(top_items) == 1
                else html.Div(top_items, className="disease-families")
            )
            families.append(
                html.Details(
                    [
                        html.Summary(f"{family['label']} ({len(choices)} terms)"),
                        contents,
                    ]
                )
                if len(choices) > 1
                else contents
            )
        groups.append(
            html.Details(
                [
                    html.Summary(group["label"]),
                    html.Div(families, className="disease-families"),
                ]
            )
        )
    return html.Div(
        [
            html.Div(
                [
                    html.Label("Diseases", htmlFor="diseases"),
                    info_tip(
                        "diseases",
                        "disease groups",
                        "Browse groups to select individual diseases and related terms, or search by name. Each checkbox selects only that term; parent and child terms are never combined. Groups organize browsing, not diagnostic classification.",
                    ),
                ],
                className="label-help",
            ),
            dcc.Dropdown(
                id="diseases",
                options=options,
                value=ordered_disease_ids(snapshot, selected),
                multi=True,
                searchable=True,
            ),
            html.Div(
                [
                    html.H4("Browse disease groups"),
                    html.Div(groups, className="disease-tree"),
                ],
                className="disease-browser",
            ),
        ],
        className="control",
    )


def choose_defaults(items, terms, limit, preferred_ids=None) -> list[str]:
    """名前と実データの有無から、存在する項目だけを初期選択する。"""
    selected = []
    names = {item_id: name.casefold() for item_id, name in items}
    for term in terms:
        exact = next(
            (item_id for item_id, name in names.items() if name == term.casefold()),
            None,
        )
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


def heatmap_cell_ids(snapshot, selected, expanded, *, catalog=None):
    """選択済みの行に、展開中の大分類の元細胞を追加する。"""
    visible = set(selected or [])
    if catalog is None:
        catalog = atlas.cell_catalog(snapshot, "mixed")
    for cell in catalog:
        if (
            cell["cell_level"] == "group"
            and cell["id"] in visible
            and cell["id"] in (expanded or [])
        ):
            visible.update(cell["members"])
    return [cell["id"] for cell in catalog if cell["id"] in visible]


def heatmap_row_controls(names, cell_ids, expanded, kind):
    """図の各行に揃えた、キーボードでも操作できる行見出し。"""
    return [
        html.Button(
            ("▼ " if cell_id in expanded else "▶ ") + names[cell_id],
            id={
                "type": "expression-cell-toggle"
                if kind == "expression"
                else "heatmap-cell-toggle",
                "kind": kind,
                "cell": cell_id,
            },
            n_clicks=0,
            title=("Collapse " if cell_id in expanded else "Expand ") + names[cell_id],
            **{"aria-expanded": "true" if cell_id in expanded else "false"},
        )
        if cell_id.startswith("group:")
        else html.Div(
            names[cell_id], className="heatmap-child-label", title=names[cell_id]
        )
        for cell_id in cell_ids
    ]


def format_data_version(value: object) -> str:
    """Open Targets の版を画面で使える短い文字列にする。"""
    if not isinstance(value, dict):
        return str(value or "")
    return ".".join(
        str(part)
        for part in (value.get("year"), value.get("month"), value.get("iteration"))
        if part is not None
    )


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


def _measure_fields(kind: str, measure: str):
    if kind == "drug":
        return (
            "drug_percent" if measure == "percent" else "drug_count",
            "unknown_drugs",
            "drug_denominator",
        )
    return "percent" if measure == "percent" else "count", "unknown", "denominator"


def _is_lower_bound(row: dict, measure: str, kind: str = "target") -> bool:
    """実数と割合で異なる下限判定を、その指標自身の未判定から求める。"""
    value_key, unknown_key, _ = _measure_fields(kind, measure)
    return row.get(value_key) is not None and (
        row.get(unknown_key, 0) > 0
        or (measure == "count" and row.get("unmapped_drugs", 0) > 0)
    )


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
        assessment = (
            "Support cannot be determined from available target and expression data."
        )
    elif row[value_key] == 0:
        assessment = f"No qualifying {label.lower()} under this rule"
    else:
        assessment = f"{label} meeting the expression rule"
    return "<br>".join(
        (
            f"<b>{row['disease']} × {row['cell']}</b>",
            f"{label}: {value}",
            assessment,
            f"Denominator: {row.get(denominator_key, 0)}",
            f"Unknown in denominator: {row.get(unknown_key, 0)}",
            f"Mapped / all canonical drugs: {row.get('mapped_drugs', 0)} / {row.get('total_drugs', 0)}",
            f"Drugs without mapped targets: {row.get('unmapped_drugs', 0)}",
            f"Status: {row['status']}",
        )
    )


def disease_label_lines(name):
    """長い疾患名を単語とハイフン付きの語を保って折り返す。"""
    return wrap(name, width=32, break_long_words=False, break_on_hyphens=False) or [
        name
    ]


def build_figure(
    rows, disease_ids, cell_ids, measure, kind="target", *, scale_rows=None
) -> go.Figure:
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
            text_row.append(
                ""
                if row is None
                else _display_value(row, measure, kind).removeprefix("≥")
            )
            hover_row.append(
                "No summary for this combination"
                if row is None
                else _hover_text(row, measure, kind)
            )
            custom_row.append([disease_id, cell_id])
            if value is None:
                missing_x.append(disease_names.get(disease_id, disease_id))
                missing_y.append(cell_names.get(cell_id, cell_id))
                missing_custom.append([disease_id, cell_id])
                missing_hover.append(
                    "No summary for this combination"
                    if row is None
                    else _hover_text(row, measure, kind)
                )
        z.append(z_row)
        texts.append(text_row)
        hovers.append(hover_row)
        customs.append(custom_row)
    values = (
        [row[value_key] for row in scale_rows if row.get(value_key) is not None]
        if scale_rows is not None
        else [value for z_row in z for value in z_row if value is not None]
    )
    color_min = min(values, default=0)
    color_max = max(values, default=1)
    if color_min == color_max:
        color_min, color_max = 0, max(1, color_max)
    title = "Drug" if kind == "drug" else "Target"
    figure = go.Figure(
        go.Heatmap(
            x=[disease_names.get(item, item) for item in disease_ids],
            y=[cell_names.get(item, item) for item in cell_ids],
            z=z,
            zmin=color_min,
            zmax=color_max,
            colorscale="Greens",
            colorbar={
                "title": f"{title} {'share (%)' if measure == 'percent' else 'count'}",
                "lenmode": "pixels",
                "len": 240,
                "y": 1,
                "yanchor": "top",
            },
            text=texts,
            texttemplate="%{text}",
            hovertext=hovers,
            hovertemplate="%{hovertext}<extra></extra>",
            customdata=customs,
            xgap=2,
            ygap=2,
        )
    )
    if missing_x:
        figure.add_trace(
            go.Scatter(
                x=missing_x,
                y=missing_y,
                customdata=missing_custom,
                hovertext=missing_hover,
                mode="markers",
                marker={"symbol": "x", "size": 10, "color": "#8a99a6"},
                name="No value (see hover)",
                hovertemplate="%{hovertext}<extra></extra>",
            )
        )
    if not disease_ids or not cell_ids:
        figure.add_annotation(text="Select diseases and cell types", showarrow=False)
    figure.update_layout(
        template="plotly_white",
        font={"family": "Arial, sans-serif", "size": 12, "color": "#263238"},
        height=max(420, min(1400, 180 + 28 * len(cell_ids))),
        margin={"l": 180, "r": 40, "t": 130, "b": 70},
        legend={"orientation": "h", "y": -0.08, "yanchor": "top", "x": 0},
        hoverlabel={"align": "left"},
    )
    labels = [disease_names.get(item, item) for item in disease_ids]
    figure.update_xaxes(
        tickangle=-45,
        side="top",
        title="Disease",
        automargin=True,
        tickmode="array",
        tickvals=labels,
        ticktext=["<br>".join(disease_label_lines(name)) for name in labels],
    )
    figure.update_yaxes(autorange="reversed", title="Cell type")
    return figure


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


def evidence_rows(
    snapshot,
    row,
    *,
    modality,
    stage,
    threshold,
    method,
    specificity,
    metadata=None,
    records=None,
    cell_names=None,
):
    """選択範囲の全元記録へ、公開判定関数による細胞別状態を添える。"""
    metadata = metadata if metadata is not None else atlas.expression_metadata(snapshot)
    records = (
        records
        if records is not None
        else [
            record
            for record in atlas.filtered_records(snapshot, modality, stage)
            if record["disease_id"] == row["disease_id"]
        ]
    )
    cell_names = (
        cell_names
        if cell_names is not None
        else {item["id"]: item["name"] for item in atlas.cell_catalog(snapshot, "cell")}
    )
    members = _ordered_cell_ids(snapshot, "cell", row["member_cell_ids"])
    output = []
    for record in records:
        base = {
            **record,
            "group_cell_id": row["cell_id"],
            "group_cell": row["cell"],
            "references_json": json.dumps(
                record.get("references") or [], ensure_ascii=False, sort_keys=True
            ),
        }
        if not record.get("target_id"):
            output.append(
                {
                    **base,
                    "evidence_cell_id": "",
                    "evidence_cell": "",
                    "cpm": None,
                    "specificity_score": None,
                    "target_median": None,
                    "support_state": "unmapped",
                    "contributing": False,
                }
            )
            continue
        for member in members:
            item = metadata.get((record["target_id"], member))
            state = atlas.expression_state(item, threshold, method, specificity)
            output.append(
                {
                    **base,
                    "evidence_cell_id": member,
                    "evidence_cell": item.get("cell", cell_names.get(member, member))
                    if item
                    else cell_names.get(member, member),
                    "cpm": item.get("median") if item else None,
                    "specificity_score": item.get("specificity_score")
                    if item
                    else None,
                    "target_median": item.get("target_median") if item else None,
                    "support_state": "positive"
                    if state is True
                    else "negative"
                    if state is False
                    else "unknown",
                    "contributing": state is True,
                }
            )
    return output


def _clustered_targets(targets, profiles):
    """全細胞の発現パターンが近い標的を隣接させ、欠測を末尾に置く。"""
    complete = [target for target in targets if target in profiles]
    if len(complete) > 1:
        tree = linkage(
            [profiles[target] for target in complete],
            method="average",
            metric="euclidean",
            optimal_ordering=True,
        )
        complete = [complete[index] for index in leaves_list(tree)]
    return complete + [target for target in targets if target not in profiles]


def expression_figure(
    snapshot,
    records,
    metadata=None,
    *,
    threshold=DEFAULT_EXPRESSION_THRESHOLD,
    method="fixed",
    specificity=DEFAULT_SPECIFICITY_THRESHOLD,
    grouped=False,
) -> go.Figure:
    """全元細胞について、現在の薬剤条件に含まれる標的の連続発現量を示す。"""
    metadata = metadata if metadata is not None else atlas.expression_metadata(snapshot)
    target_names = {}
    for record in records:
        if record.get("target_id"):
            target_names.setdefault(
                record["target_id"], record.get("target") or record["target_id"]
            )
    cells = atlas.cell_catalog(snapshot, "cell")
    members = [cell["id"] for cell in cells]
    targets = sorted(
        target_names, key=lambda item: (target_names[item].casefold(), item)
    )
    # 表示用の大分類平均は、標準化とクラスタリングの基準に含めない。
    target_stats = {}
    profiles = {}
    for target_id in targets:
        values = [
            math.log2(1 + metadata[target_id, member]["median"])
            for member in members
            if (target_id, member) in metadata
            and metadata[target_id, member].get("median") is not None
        ]
        target_stats[target_id] = (mean(values), pstdev(values)) if values else (0, 0)
        if members and len(values) == len(members):
            center, spread = target_stats[target_id]
            profiles[target_id] = tuple(
                (value - center) / spread if spread else 0 for value in values
            )
    targets = _clustered_targets(targets, profiles)
    z, hover, missing_x, missing_y = [], [], [], []
    positive_x, positive_y, positive_hover = [], [], []
    display_cells = atlas.cell_catalog(snapshot, "mixed") if grouped else cells
    member_names = {cell["id"]: cell["name"] for cell in display_cells}
    for cell in display_cells:
        member = cell["id"]
        is_group = cell.get("cell_level") == "group"
        z_row, hover_row = [], []
        for target_id in targets:
            label = f"{target_names[target_id]} ({target_id})"
            item = metadata.get((target_id, member))
            cpm = item.get("median") if item else None
            observed = []
            if is_group:
                observed = [
                    metadata[target_id, child]["median"]
                    for child in cell["members"]
                    if (target_id, child) in metadata
                    and metadata[target_id, child].get("median") is not None
                ]
                cpm = mean(observed) if observed else None
            center, spread = target_stats[target_id]
            score = (
                (math.log2(1 + cpm) - center) / spread
                if cpm is not None and spread
                else (0 if cpm is not None else None)
            )
            z_row.append(score)
            hover_row.append(
                "<br>".join(
                    (
                        f"<b>{label}</b>",
                        f"Cell: {member_names[member]} ({member})",
                        f"Target-wise z-score: {score:.2f}"
                        if score is not None
                        else "Target-wise z-score: missing",
                        f"Median CPM: {cpm:g}"
                        if cpm is not None
                        else "Median CPM: missing",
                        f"CELLEX specificity: {item.get('specificity_score'):g}"
                        if item and item.get("specificity_score") is not None
                        else "CELLEX specificity: missing",
                        f"Target-relative median: {item.get('target_median'):g}"
                        if item and item.get("target_median") is not None
                        else "Target-relative median: missing",
                    )
                )
            )
            if is_group:
                hover_row[-1] = "<br>".join(
                    (
                        f"<b>{label}</b>",
                        f"Cell group: {cell['name']}",
                        f"Mean CPM: {cpm:g}"
                        if cpm is not None
                        else "Mean CPM: missing",
                        f"Observed source cells: {len(observed)} / {len(cell['members'])}",
                        f"Target-wise z-score: {score:.2f}"
                        if score is not None
                        else "Target-wise z-score: missing",
                        "Arithmetic mean of available source-cell median CPM values; equal weight per cell type.",
                    )
                )
            if cpm is None:
                missing_x.append(label)
                missing_y.append(member_names[member])
            elif (
                not is_group
                and atlas.expression_state(item, threshold, method, specificity) is True
            ):
                positive_x.append(label)
                positive_y.append(member_names[member])
                positive_hover.append(hover_row[-1] + "<br>Expression rule: met")
        z.append(z_row)
        hover.append(hover_row)
    target_labels = [f"{target_names[t]} ({t})" for t in targets]
    color_limit = (
        max(
            (abs(value) for z_row in z for value in z_row if value is not None),
            default=0,
        )
        or 1
    )
    figure = go.Figure(
        go.Heatmap(
            x=target_labels,
            y=[cell["name"] for cell in display_cells],
            z=z,
            zmin=-color_limit,
            zmax=color_limit,
            colorscale="RdBu_r",
            colorbar={
                "title": "Target-wise z-score",
                "lenmode": "pixels",
                "len": 240,
                "y": 1,
                "yanchor": "top",
            },
            hovertext=hover,
            hovertemplate="%{hovertext}<extra></extra>",
            xgap=2,
            ygap=2,
        )
    )
    if missing_x:
        figure.add_trace(
            go.Scatter(
                x=missing_x,
                y=missing_y,
                mode="markers",
                marker={"symbol": "x", "size": 9, "color": "#8a99a6"},
                name="Missing expression",
                hovertemplate="Expression missing<extra></extra>",
            )
        )
    if positive_x:
        figure.add_trace(
            go.Scattergl(
                x=positive_x,
                y=positive_y,
                mode="markers",
                marker={
                    "symbol": "circle",
                    "size": 7,
                    "color": "white",
                    "line": {"color": "#263238", "width": 1},
                },
                name="Meets expression rule",
                showlegend=False,
                hovertext=positive_hover,
                hovertemplate="%{hovertext}<extra></extra>",
            )
        )
    if not targets:
        figure.add_annotation(
            text="No known targets in the selected scope", showarrow=False
        )
    figure.update_layout(
        template="plotly_white",
        height=max(360, 170 + 25 * len(display_cells)),
        margin={"l": 280, "r": 40, "t": 110, "b": 60},
        font={"family": "Arial, sans-serif", "size": 11, "color": "#263238"},
        hoverlabel={"align": "left"},
    )
    figure.update_xaxes(
        tickangle=-45,
        side="top",
        title="Known target",
        tickvals=target_labels,
        ticktext=[target_names[t] for t in targets],
        automargin=True,
    )
    figure.update_yaxes(autorange="reversed", title="Source cell type", automargin=True)
    if members:
        figure.update_yaxes(range=[len(display_cells) - 0.5, -0.5], autorange=False)
    return figure


def expression_view(base, catalog, expanded):
    """集計済みの図から表示行だけを選び、色範囲と標的順を保つ。"""
    groups = [cell["id"] for cell in catalog if cell["cell_level"] == "group"]
    visible = set(heatmap_cell_ids(None, groups, expanded, catalog=catalog))
    names = {cell["name"] for cell in catalog if cell["id"] in visible}
    figure = go.Figure(base)
    heatmap = figure.data[0]
    indices = [i for i, name in enumerate(heatmap.y) if name in names]
    for attr in ("y", "z", "hovertext"):
        values = getattr(heatmap, attr)
        setattr(heatmap, attr, [values[i] for i in indices])
    for trace in figure.data[1:]:
        keep = [i for i, name in enumerate(trace.y) if name in names]
        for attr in ("x", "y", "hovertext"):
            values = getattr(trace, attr)
            if values is not None:
                setattr(trace, attr, [values[i] for i in keep])
    count = max(1, len(indices))
    figure.update_layout(
        height=110 + 60 + 28 * count, margin=dict(l=280, r=40, t=110, b=60)
    )
    figure.update_xaxes(automargin=False)
    figure.update_yaxes(
        range=[count - 0.5, -0.5],
        autorange=False,
        automargin=False,
        fixedrange=True,
        showticklabels=False,
        title=None,
    )
    heatmap.colorbar.len = min(240, 28 * count)
    return figure


def _reference_links(record):
    links = []
    if record.get("target_id"):
        links.append(
            html.A(
                "Open Targets target",
                href="https://platform.opentargets.org/target/" + record["target_id"],
                target="_blank",
                rel="noreferrer",
            )
        )
    if record.get("drug_id"):
        links.append(
            html.A(
                "Open Targets drug",
                href="https://platform.opentargets.org/drug/" + record["drug_id"],
                target="_blank",
                rel="noreferrer",
            )
        )
    if record.get("canonical_drug_id") and record.get(
        "canonical_drug_id"
    ) != record.get("drug_id"):
        links.append(
            html.A(
                "Open Targets canonical drug",
                href="https://platform.opentargets.org/drug/"
                + record["canonical_drug_id"],
                target="_blank",
                rel="noreferrer",
            )
        )
    if record.get("disease_id"):
        links.append(
            html.A(
                "Open Targets disease",
                href="https://platform.opentargets.org/disease/" + record["disease_id"],
                target="_blank",
                rel="noreferrer",
            )
        )
    children = []
    for index, item in enumerate(links):
        if index:
            children.append(html.Br())
        children.append(item)
    return html.Div(children) if children else "—"


def _evidence_table(records):
    headers = (
        "Canonical drug",
        "Original drug",
        "Modality",
        "Record stage / highest disease-specific drug stage",
        "Target",
        "Action / mechanism",
        "Open Targets links",
    )
    body = []
    for record in records:
        actions, mechanism = (
            ", ".join(record.get("action_types") or []),
            record.get("mechanism") or "—",
        )
        body.append(
            html.Tr(
                [
                    html.Td(
                        f"{record.get('canonical_drug') or '—'} ({record.get('canonical_drug_id') or '—'})"
                    ),
                    html.Td(
                        f"{record.get('drug') or '—'} ({record.get('drug_id') or '—'})"
                    ),
                    html.Td(record.get("drug_type") or "—"),
                    html.Td(
                        f"{record.get('stage') or '—'} / {record.get('canonical_stage') or '—'}"
                    ),
                    html.Td(
                        f"{record.get('target') or '—'} ({record.get('target_id') or '—'})"
                    ),
                    html.Td(f"{actions or '—'} / {mechanism}"),
                    html.Td(_reference_links(record)),
                ]
            )
        )
    return html.Div(
        html.Table(
            [html.Thead(html.Tr([html.Th(item) for item in headers])), html.Tbody(body)]
        ),
        className="table-scroll",
    )


def detail_panel(
    rows,
    selection,
    snapshot=None,
    *,
    modality="all",
    stage="phase3",
    threshold=DEFAULT_EXPRESSION_THRESHOLD,
    method="fixed",
    specificity=DEFAULT_SPECIFICITY_THRESHOLD,
):
    """選択した疾患の連続発現と元記録を表示する。"""
    if selection is None:
        return html.Div(
            "Select a disease using the displayed heatmap or the selector above.",
            className="empty-note",
        )
    row = next((item for item in rows if item["disease_id"] == selection), None)
    if row is None:
        return html.Div(
            "The previous disease is outside the current filters. Select a visible disease.",
            className="empty-note",
        )
    if snapshot is None:
        return html.Div(
            [
                html.P(
                    "Underlying schema 2 data are required for evidence details.",
                    className="empty-note",
                ),
            ]
        )
    filtered = [
        record
        for record in atlas.filtered_records(snapshot, modality, stage)
        if record["disease_id"] == row["disease_id"]
    ]
    source_context = {
        "disease_id": row["disease_id"],
        "modality": modality,
        "stage": stage,
        "threshold": threshold,
        "method": method,
        "specificity": specificity,
    }
    pages = max(1, math.ceil(len(filtered) / SOURCE_PAGE_SIZE))
    return html.Div(
        [
            html.H3(
                [
                    f"Filtered drug records for {row['disease']}",
                    info_tip(
                        "source-records",
                        "source records",
                        "Only records for the selected disease and the applied clinical stage and modality are shown. Each original drug record appears once.",
                    ),
                ]
            ),
            html.Div(
                [
                    dcc.Store(id="source-context", data=source_context),
                    dcc.Loading(html.Div(id="source-records-page")),
                    html.Div(
                        [
                            html.Label("Page", htmlFor="source-page"),
                            dcc.Dropdown(
                                id="source-page",
                                options=[
                                    {"label": f"{page} / {pages}", "value": page}
                                    for page in range(1, pages + 1)
                                ],
                                value=1,
                                clearable=False,
                                searchable=False,
                            ),
                        ],
                        className="source-pagination",
                    ),
                ],
                className="source-records",
            ),
            html.H3(
                [
                    "Relative expression of targets across all source cell types",
                    info_tip(
                        "expression",
                        "target expression",
                        "Click a cell group to expand or collapse its source cells. Group values are arithmetic means of available source-cell median CPM values, with equal weight per cell type; hover shows coverage. The group mean is log2(1 + CPM) transformed and standardized using the same target-wise mean and standard deviation as all healthy reference source cells. A white dot with a dark outline marks an individual source-cell pair meeting the applied expression rule; group averages have no dot. Targets with similar z-score patterns are grouped together; targets with missing expression are shown last. Group averages never enter target-relative medians, z-score reference statistics, clustering, or expression-rule calculations. Expansion does not change the color range or target order. The color range follows the largest absolute z-score and stays centered on zero; hover shows the z-score and raw median CPM. This is not disease-sample expression.",
                    ),
                    html.Small(
                        "○ Meets expression rule", className="expression-marker-key"
                    ),
                ],
                className="detail-expression-heading",
            ),
            html.Div(
                html.Div(
                    [
                        dcc.Graph(
                            id="expression-heatmap",
                            config={"displaylogo": False, "responsive": True},
                        ),
                        html.Div(
                            id="expression-row-controls",
                            className="heatmap-row-controls",
                            style={"top": "110px", "width": "280px"},
                        ),
                    ],
                    className="expandable-heatmap",
                ),
                className="graph-scroll expression-graph",
            ),
        ]
    )


def unavailable_layout(error=None) -> html.Main:
    """データ未取得または再取得が必要な状態を説明する。"""
    title = "No data loaded yet" if error is None else "Data refresh required"
    return html.Main(
        [
            html.Section(
                [
                    html.P("AUTOIMMUNE DRUG–CELL ATLAS", className="eyebrow"),
                    html.H1(title),
                    html.P("Refresh the snapshot to load the comparison."),
                    html.Code("pixi run refresh"),
                    *([html.P(str(error), className="empty-note")] if error else []),
                ],
                className="unavailable-card",
            )
        ],
        className="shell unavailable",
    )


def dashboard_layout(snapshot) -> html.Main:
    """schema 2 のデータから dashboard の初期画面を作る。"""
    diseases = [(row["id"], row["name"]) for row in snapshot["diseases"]]
    default_diseases = choose_defaults(
        diseases,
        (
            "systemic lupus erythematosus",
            "systemic sclerosis",
            "Sjogren syndrome",
            "rheumatoid arthritis",
            "multiple sclerosis",
            "dermatomyositis",
            "type 1 diabetes mellitus",
            "anti-neutrophil cytoplasmic antibody-associated vasculitis",
            "pemphigus",
            "autoimmune hepatitis",
        ),
        10,
        {row["disease_id"] for row in snapshot["records"]},
    )
    retrieved_at = (
        datetime.fromisoformat(snapshot["retrieved_at"])
        .astimezone(timezone(timedelta(hours=9)))
        .strftime("%Y-%m-%d %H:%M JST")
    )
    ready = sum(row["status"] == "ready" for row in snapshot["diseases"])
    data_version = format_data_version(snapshot.get("data_version"))
    source = "https://platform.opentargets.org/"
    modality_options = [{"label": "All", "value": "all"}] + [
        {"label": label, "value": value} for label, value in atlas.DRUG_TYPE_MODALITIES
    ]

    def control(label, component):
        return html.Div(
            [html.Label(label, htmlFor=component.id), component], className="control"
        )

    applied = dict(
        zip(
            PARAMETER_IDS,
            (
                "percent",
                "all",
                "phase3",
                "specificity",
                DEFAULT_EXPRESSION_THRESHOLD,
                DEFAULT_SPECIFICITY_THRESHOLD,
                ordered_disease_ids(snapshot, default_diseases),
                "target",
            ),
        )
    )
    return html.Main(
        [
            html.Nav(
                [
                    html.A("Autoimmune Atlas", href="#", className="app-brand"),
                    html.Div(
                        [
                            html.A("Comparison", href="#comparison"),
                            html.A("Evidence", href="#evidence"),
                            html.A(
                                "Open Targets ↗",
                                href=source,
                                target="_blank",
                                rel="noreferrer",
                            ),
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
                            html.P(
                                "AUTOIMMUNE DISEASE / DRUG TARGETS", className="eyebrow"
                            ),
                            html.H1(
                                [
                                    "Drug targets by cell type",
                                    info_tip(
                                        "overview",
                                        "the atlas",
                                        "Compare which drug targets meet an expression rule in healthy reference cells across diseases. Disease differences reflect eligible drug and target sets, not disease-specific expression. This does not establish treatment efficacy; disease records may include symptom or comorbidity treatment.",
                                    ),
                                ]
                            ),
                        ]
                    ),
                    html.Section(
                        [
                            html.Div(
                                [
                                    html.Span("Disease terms", className="meta-label"),
                                    html.Strong(
                                        f"{ready} loaded / {len(snapshot['diseases'])} total"
                                    ),
                                ]
                            ),
                            html.Div(
                                [
                                    html.Span("Drug records", className="meta-label"),
                                    html.Strong(str(len(snapshot["records"]))),
                                ]
                            ),
                            html.Div(
                                [
                                    html.Span("Targets", className="meta-label"),
                                    html.Strong(str(len(snapshot["expression"]))),
                                ]
                            ),
                            html.Div(
                                [
                                    html.Span("Retrieved at", className="meta-label"),
                                    html.Strong(
                                        retrieved_at, title=snapshot["retrieved_at"]
                                    ),
                                ]
                            ),
                            html.Div(
                                [
                                    html.Span(
                                        "Expression reference", className="meta-label"
                                    ),
                                    html.Strong("Tabula Sapiens"),
                                ]
                            ),
                            html.Div(
                                [
                                    html.Span("Data source", className="meta-label"),
                                    html.A(
                                        f"Open Targets {data_version}".strip(),
                                        href=source,
                                        target="_blank",
                                        rel="noreferrer",
                                    ),
                                ]
                            ),
                        ],
                        className="source-bar",
                        **{"aria-label": "Snapshot data"},
                    ),
                ],
                className="hero",
            ),
            html.Section(
                [
                    html.H2("Settings panel"),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H3("Comparison scope"),
                                    html.Div(
                                        [
                                            disease_selector(
                                                snapshot, default_diseases
                                            ),
                                        ],
                                        className="control-grid scope-controls",
                                    ),
                                ],
                                className="filter-group scope-group",
                            ),
                            html.Section(
                                [
                                    html.H3("Drug evidence"),
                                    html.Div(
                                        [
                                            html.Div(
                                                [
                                                    html.Label(
                                                        "Clinical stage",
                                                        htmlFor="stage",
                                                    ),
                                                    info_tip(
                                                        "stage",
                                                        "clinical stage",
                                                        "Highest recorded stage of a canonical drug in this disease. Approval reached may include withdrawn drugs.",
                                                    ),
                                                ],
                                                className="label-help",
                                            ),
                                            dcc.Dropdown(
                                                id="stage",
                                                options=[
                                                    {"label": label, "value": value}
                                                    for label, value in (
                                                        ("Phase I or later", "phase1"),
                                                        ("Phase II or later", "phase2"),
                                                        (
                                                            "Phase III or later",
                                                            "phase3",
                                                        ),
                                                        (
                                                            "Approval reached",
                                                            "approved",
                                                        ),
                                                    )
                                                ],
                                                value="phase3",
                                                clearable=False,
                                            ),
                                        ],
                                        className="control",
                                    ),
                                    control(
                                        "Drug modality",
                                        dcc.Dropdown(
                                            id="modality",
                                            options=modality_options,
                                            value="all",
                                            clearable=False,
                                        ),
                                    ),
                                ],
                                className="filter-group",
                            ),
                            html.Section(
                                [
                                    html.H3("Expression criteria"),
                                    control(
                                        "Expression rule",
                                        dcc.Dropdown(
                                            id="method",
                                            options=[
                                                {"label": label, "value": value}
                                                for value, label in METHOD_LABELS.items()
                                            ],
                                            value="specificity",
                                            clearable=False,
                                        ),
                                    ),
                                    # Dash 4.4.1 は max を省略すると増減時に NaN になるため、上限なしを明示する。
                                    html.Div(
                                        [
                                            html.Div(
                                                [
                                                    html.Label(
                                                        "Minimum CPM (≥)",
                                                        htmlFor="threshold",
                                                    ),
                                                    info_tip(
                                                        "threshold",
                                                        "minimum CPM",
                                                        f"Every rule requires median CPM at least this value. Blank uses {DEFAULT_EXPRESSION_THRESHOLD:g} CPM.",
                                                    ),
                                                ],
                                                className="label-help",
                                            ),
                                            dcc.Input(
                                                id="threshold",
                                                type="number",
                                                min=0,
                                                max=None,
                                                step=0.1,
                                                value=DEFAULT_EXPRESSION_THRESHOLD,
                                            ),
                                        ],
                                        className="control",
                                    ),
                                    html.Div(
                                        [
                                            html.Div(
                                                [
                                                    html.Label(
                                                        "CELLEX specificity (≥)",
                                                        htmlFor="specificity",
                                                    ),
                                                    info_tip(
                                                        "specificity",
                                                        "CELLEX specificity",
                                                        f"Only used with Fixed CPM + CELLEX specificity. Blank uses {DEFAULT_SPECIFICITY_THRESHOLD:g}.",
                                                    ),
                                                ],
                                                className="label-help",
                                            ),
                                            dcc.Input(
                                                id="specificity",
                                                type="number",
                                                min=0,
                                                max=1,
                                                step=0.05,
                                                value=DEFAULT_SPECIFICITY_THRESHOLD,
                                                disabled=False,
                                            ),
                                        ],
                                        className="control",
                                    ),
                                ],
                                className="filter-group",
                            ),
                            html.Section(
                                [
                                    html.H3("Display"),
                                    html.Fieldset(
                                        [
                                            html.Legend(
                                                [
                                                    "View",
                                                    info_tip(
                                                        "view",
                                                        "heatmap view",
                                                        "Targets are distinct genes meeting the rule. Drug forms mapped to the same active ingredient count once if any known target meets the rule.",
                                                    ),
                                                ]
                                            ),
                                            dcc.RadioItems(
                                                id="heatmap-view",
                                                options=[
                                                    {
                                                        "label": "Distinct targets",
                                                        "value": "target",
                                                    },
                                                    {
                                                        "label": "Canonical drugs",
                                                        "value": "drug",
                                                    },
                                                ],
                                                value="target",
                                                inline=True,
                                            ),
                                        ],
                                        className="control radio-control",
                                    ),
                                    html.Fieldset(
                                        [
                                            html.Legend(
                                                [
                                                    "Measure",
                                                    info_tip(
                                                        "measure",
                                                        "measure",
                                                        "Percent divides by known targets or canonical drugs with known targets for the current disease and drug filters. Each denominator stays fixed across cells; percent is not a share of cells.",
                                                    ),
                                                ]
                                            ),
                                            dcc.RadioItems(
                                                id="measure",
                                                options=[
                                                    {
                                                        "label": "Percent",
                                                        "value": "percent",
                                                    },
                                                    {
                                                        "label": "Count",
                                                        "value": "count",
                                                    },
                                                ],
                                                value="percent",
                                                inline=True,
                                            ),
                                        ],
                                        className="control radio-control",
                                    ),
                                ],
                                className="filter-group",
                            ),
                            html.Div(
                                [
                                    html.Button(
                                        "Update",
                                        id="update-button",
                                        type="button",
                                        n_clicks=0,
                                    ),
                                    html.Span(
                                        "Settings applied.",
                                        id="update-status",
                                        role="status",
                                    ),
                                ],
                                className="update-actions",
                            ),
                        ],
                        className="filter-groups",
                    ),
                    dcc.Store(id="applied-parameters", data=applied),
                    dcc.Store(id="expanded-cell-groups", data=[]),
                    dcc.Store(id="expanded-expression-groups", data=[]),
                ],
                className="panel controls",
            ),
            html.Section(
                [
                    html.Div(
                        [
                            html.H2(
                                [
                                    "Cell-type comparison",
                                    info_tip(
                                        "comparison",
                                        "cell-type comparison",
                                        "Click a cell group to expand or collapse its source cells. A group counts the union of qualifying targets or drugs across its members; CPM values are never added or averaged. Heatmap cells show confirmed support without a ≥ mark; hover values mark lower bounds. × can mean unavailable disease data, unresolved target or expression data, or no eligible percentage denominator. A zero count means no qualifying targets or drugs under the rule.",
                                    ),
                                ]
                            ),
                        ],
                        className="section-heading",
                    ),
                    html.Div(id="matrix-note", className="matrix-note", role="status"),
                    html.Div(
                        [
                            html.Article(
                                [
                                    html.H3("Distinct targets"),
                                    html.Div(
                                        html.Div(
                                            [
                                                dcc.Graph(
                                                    id="target-heatmap",
                                                    style={"minWidth": "600px"},
                                                    config={
                                                        "displaylogo": False,
                                                        "responsive": True,
                                                    },
                                                ),
                                                html.Div(
                                                    id="target-row-controls",
                                                    className="heatmap-row-controls",
                                                ),
                                            ],
                                            className="expandable-heatmap",
                                        ),
                                        className="graph-scroll",
                                    ),
                                ],
                                id="target-heatmap-panel",
                            ),
                            html.Article(
                                [
                                    html.H3("Canonical drugs"),
                                    html.Div(
                                        html.Div(
                                            [
                                                dcc.Graph(
                                                    id="drug-heatmap",
                                                    style={"minWidth": "600px"},
                                                    config={
                                                        "displaylogo": False,
                                                        "responsive": True,
                                                    },
                                                ),
                                                html.Div(
                                                    id="drug-row-controls",
                                                    className="heatmap-row-controls",
                                                ),
                                            ],
                                            className="expandable-heatmap",
                                        ),
                                        className="graph-scroll",
                                    ),
                                ],
                                id="drug-heatmap-panel",
                                hidden=True,
                            ),
                        ],
                        className="heatmap-stack",
                    ),
                    html.P(
                        "Hover for values · Click for details",
                        className="matrix-note",
                    ),
                ],
                className="panel matrix-panel",
                id="comparison",
            ),
            html.Section(
                [
                    html.Div(
                        [
                            html.H2(
                                [
                                    "Selection details",
                                    info_tip(
                                        "selection",
                                        "selection details",
                                        "Select a disease to see target expression across healthy reference source cells and its filtered drug records. Click a comparison heatmap column to select its disease.",
                                    ),
                                ]
                            )
                        ],
                        className="section-heading detail-heading",
                    ),
                    html.Div(
                        [
                            control(
                                "Disease",
                                dcc.Dropdown(id="detail-disease", clearable=False),
                            )
                        ],
                        className="control-grid selectors",
                    ),
                    html.Div(id="details", className="details"),
                ],
                className="panel",
                id="evidence",
            ),
        ],
        className="shell",
    )


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
        note = f"Applied: Clinical stage: {STAGE_LABELS[stage]}; Drug modality: {modality_label}; Expression rule: {METHOD_LABELS[method]} ({condition}); Cells: all groups and expanded source cells; Measure: {measure.capitalize()}. {len(rows)} disease–cell combinations. Heatmap entries without a value: targets {target_missing}, drugs {drug_missing}."
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


def create_app(snapshot: dict | None, error: Exception | None = None) -> Dash:
    """保存済みデータまたは明示的な fixture から Dash アプリを作る。"""
    # 元記録の操作部は、選択した詳細を描画するときに追加する。
    application = Dash(
        __name__,
        title="Autoimmune Target Expression Atlas",
        suppress_callback_exceptions=True,
    )
    application.layout = (
        unavailable_layout(error) if snapshot is None else dashboard_layout(snapshot)
    )
    if snapshot is not None:
        register_callbacks(application, snapshot)
    return application


try:
    SNAPSHOT = load_snapshot(SNAPSHOT_PATH)
    LOAD_ERROR = None
except (OSError, json.JSONDecodeError, ValueError, TypeError, KeyError) as error:
    SNAPSHOT = None
    LOAD_ERROR = error

app = create_app(SNAPSHOT, LOAD_ERROR)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Autoimmune Atlas")
    parser.add_argument(
        "--dev", action="store_true", help="Automatically reload when code changes"
    )
    args = parser.parse_args()
    app.run(host="127.0.0.1", port=8050, debug=args.dev)
