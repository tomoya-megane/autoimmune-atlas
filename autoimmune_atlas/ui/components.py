"""Dash の再利用可能な画面部品。"""

from collections.abc import Collection, Iterable, Mapping, Sequence
from typing import Literal, TypedDict

import dash_ag_grid as dag  # pyright: ignore[reportMissingTypeStubs] - dash-ag-grid に型スタブがない。
from dash import dcc, html
from dash.development.base_component import Component

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas.disease_catalog import disease_catalog, ordered_disease_ids
from autoimmune_atlas.models import (
    CatalogDisease,
    DiseaseCatalogGroup,
    DiseaseFamily,
    FilteredRecord,
    Snapshot,
)
from autoimmune_atlas.ui.config import SOURCE_PAGE_SIZE


class DiseaseSection(TypedDict):
    """同じ選択欄に表示する疾患の区分。"""

    id: str
    diseases: list[CatalogDisease]
    node: str | None
    kind: Literal["self", "children"]


def page_nav(current: Literal["drugs", "genetics"]) -> html.Nav:
    """2 ページへのリンクと Open Targets へのリンクを持つ上部バー。"""

    def link(label: str, href: str, key: str) -> dcc.Link:
        # dcc.Link は pushState で切り替えるので、再読み込みせず、もう一方のページの状態が残る。
        # dcc.Link は aria-current を受け取らないので、表示中のページは current クラスで示す。
        if key == current:
            return dcc.Link(
                label,
                href=href,
                className="app-page-link current",
                title="Current page",
            )
        return dcc.Link(label, href=href, className="app-page-link")

    return html.Nav(
        [
            # dcc.Link なので、ロゴから薬剤ページへ戻るときも再読み込みしない。
            dcc.Link("Autoimmune Atlas", href="/", className="app-brand"),
            html.Div(
                [
                    link("Drug targets", "/", "drugs"),
                    link("Genetic associations", "/genetics", "genetics"),
                    html.A(
                        "Open Targets",
                        href="https://platform.opentargets.org/",
                        target="_blank",
                        rel="noreferrer",
                    ),
                ],
                className="app-nav",
            ),
        ],
        className="app-bar",
        **{  # pyright: ignore[reportArgumentType] - Dash の型定義に ARIA kwargs がない。
            "aria-label": "Main navigation"
        },
    )


def info_tip(key: str, label: str, description: str) -> html.Span:
    """操作対象の近くに、キーボードでも読める説明を置く。"""
    tip_id = f"{key}-tip"
    return html.Span(
        [
            html.Button(
                html.Span(
                    className="ui-icon icon-info",
                    **{  # pyright: ignore[reportArgumentType] - Dash の型定義に ARIA kwargs がない。
                        "aria-hidden": "true"
                    },
                ),
                type="button",
                className="info-button",
                **{  # pyright: ignore[reportArgumentType] - Dash の型定義に ARIA kwargs がない。
                    "aria-label": f"About {label}",
                    "aria-describedby": tip_id,
                },
            ),
            html.Span(description, id=tip_id, role="tooltip", className="info-content"),
        ],
        className="info-tip",
    )


def cell_catalog(snapshot: Snapshot, level: str = "group") -> list[tuple[str, str]]:
    """公開関数が返す細胞カタログを選択欄向けに整える。"""
    return [(row["id"], row["name"]) for row in atlas.cell_catalog(snapshot, level)]


def ordered_cell_ids(
    snapshot: Snapshot, level: str, cell_ids: Iterable[str] | None
) -> list[str]:
    selected = set(cell_ids or ())
    return [
        cell_id for cell_id, _ in cell_catalog(snapshot, level) if cell_id in selected
    ]


def _disease_tree(
    family: DiseaseFamily,
) -> tuple[list[str], dict[str, list[str]]]:
    """ファミリー内の親子関係を、各語が 1 つの親の下にだけ現れる木にする。

    複数の親を持つ語は、より深い親の下に置き、同じ深さなら名前順で決める。
    親子が循環して根から届かない語は、最上位に置いて失わない。
    戻り値は (最上位の語の一覧, 語 ID → 子の一覧)。
    """
    members: dict[str, CatalogDisease] = {
        disease["id"]: disease for disease in family["diseases"]
    }
    parents_of: dict[str, list[str]] = {
        did: [p for p in d.get("parent_ids", []) if p in members and p != did]
        for did, d in members.items()
    }
    depth: dict[str, int] = {}

    def node_depth(did: str, trail: tuple[str, ...] = ()) -> int:
        if did not in depth:
            parents = [p for p in parents_of[did] if p not in trail]
            depth[did] = (
                1 + max(node_depth(p, trail + (did,)) for p in parents)
                if parents
                else 0
            )
        return depth[did]

    def by_name(did: str) -> tuple[str, str]:
        return members[did]["name"].casefold(), did

    parent: dict[str, str] = {}
    for did in sorted(members, key=by_name):
        if parents_of[did]:
            parent[did] = min(
                parents_of[did], key=lambda p: (-node_depth(p), *by_name(p))
            )

    def reaches_top(did: str, trail: tuple[str, ...] = ()) -> bool:
        return did not in parent or (
            did not in trail and reaches_top(parent[did], trail + (did,))
        )

    children: dict[str, list[str]] = {did: [] for did in members}
    top: list[str] = []
    for did in sorted(members, key=by_name):
        if did in parent and reaches_top(did):
            children[parent[did]].append(did)
        else:
            parent.pop(did, None)
            top.append(did)
    top.sort(key=lambda did: (did != family["id"], *by_name(did)))
    return top, children


def disease_checklist_sections(
    family: DiseaseFamily, prefix: str = ""
) -> list[DiseaseSection]:
    """疾患本体と各階層の選択欄を、描画と同期で同じ範囲に分ける。

    子を持つ語は、ファミリーと同じ形にする。
    本体だけの section（kind が self）と、子のうち葉だけを入れた section（kind が
    children）を持ち、子を持つ子は同じ形で入れ子になる。
    先頭の section はファミリーの本体と、親を持たない葉である。
    """
    top, children = _disease_tree(family)
    members = {d["id"]: d for d in family["diseases"]}
    root = family["id"] if family["id"] in members else None

    def section_id(kind: str, node: str) -> str:
        suffix = "" if node == root else f"-{node}"
        return f"{prefix}{kind}-{family['id']}{suffix}"

    sections: list[DiseaseSection] = [
        {
            "id": f"{prefix}disease-family-{family['id']}",
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


def _disease_family_selector(
    family: DiseaseFamily, selected_ids: set[str], prefix: str = ""
) -> tuple[list[dict[str, str]], Component]:
    """1つの疾患ファミリーの選択肢と階層表示を組み立てる。"""
    choices = [
        {"label": disease["name"], "value": disease["id"]}
        for disease in family["diseases"]
    ]
    sections = disease_checklist_sections(family, prefix)
    top, children = _disease_tree(family)
    root = (
        family["id"]
        if family["id"] in {disease["id"] for disease in family["diseases"]}
        else None
    )
    by_kind = {(section["kind"], section["node"]): section for section in sections}
    names = {disease["id"]: disease["name"] for disease in family["diseases"]}

    def descendant_count(node: str) -> int:
        return sum(1 + descendant_count(child) for child in children[node])

    def checklist(section: DiseaseSection) -> dcc.Checklist:
        return dcc.Checklist(
            id=section["id"],
            options=[
                {"label": disease["name"], "value": disease["id"]}
                for disease in section["diseases"]
            ],
            value=[
                disease["id"]
                for disease in section["diseases"]
                if disease["id"] in selected_ids
            ],
            className="disease-checklist",
        )

    def render_details(node: str) -> html.Details:
        """子を持つ語の details。葉のチェック欄と、子を持つ子の入れ子を並べる。"""
        inner: list[Component] = []
        if ("children", node) in by_kind:
            inner.append(checklist(by_kind[("children", node)]))
        inner.extend(render_node(child) for child in children[node] if children[child])
        return html.Details(
            [
                html.Summary(f"{names[node]} details ({descendant_count(node)} terms)"),
                html.Div(inner, className="disease-families"),
            ]
        )

    def render_node(node: str) -> html.Details:
        """本体のチェック欄と details を、ファミリーと同じ形で包む。"""
        return html.Details(
            [
                html.Summary(f"{names[node]} ({1 + descendant_count(node)} terms)"),
                html.Div(
                    [checklist(by_kind[("self", node)]), render_details(node)],
                    className="disease-families",
                ),
            ]
        )

    top_items: list[Component] = [checklist(sections[0])]
    if root and children[root]:
        top_items.append(render_details(root))
    top_items.extend(
        render_node(node) for node in top if children[node] and node != root
    )
    contents = (
        top_items[0]
        if len(top_items) == 1
        else html.Div(top_items, className="disease-families")
    )
    component = (
        html.Details(
            [
                html.Summary(f"{family['label']} ({len(choices)} terms)"),
                contents,
            ]
        )
        if len(choices) > 1
        else contents
    )
    return choices, component


def disease_selector(
    snapshot: Snapshot, selected: Iterable[str] | None, prefix: str = ""
) -> html.Div:
    """検索欄と、群から開けるチェック欄を同じ選択へ結び付ける。"""
    catalog: list[DiseaseCatalogGroup] = disease_catalog(snapshot)
    selected_ids = set(selected or ())
    options: list[dict[str, str]] = []
    groups: list[Component] = []
    for group in catalog:
        families: list[Component] = []
        for family in group["families"]:
            choices, component = _disease_family_selector(family, selected_ids, prefix)
            options.extend(choices)
            families.append(component)
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
                    html.Label("Diseases", htmlFor=f"{prefix}diseases"),
                    info_tip(
                        f"{prefix}diseases",
                        "diseases",
                        "Diseases sets which diseases the comparison includes. Each checkbox selects one disease; selecting a parent does not select its children. The disease groups are navigation aids, not a diagnostic classification. Search by name or browse the disease groups, then click Update.",
                    ),
                ],
                className="label-help",
            ),
            dcc.Dropdown(
                id=f"{prefix}diseases",
                options=options,
                value=ordered_disease_ids(snapshot, list(selected_ids)),
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


def heatmap_row_controls(
    names: Mapping[str, str],
    cell_ids: Iterable[str],
    expanded: Collection[str],
    kind: str,
) -> list[html.Button | html.Div]:
    """図の各行に揃えた、キーボードでも操作できる行見出し。"""
    return [
        html.Button(
            names[cell_id],
            id={
                "type": "heatmap-cell-toggle"
                if kind in ("target", "drug")
                else f"{kind}-cell-toggle",
                "kind": kind,
                "cell": cell_id,
            },
            n_clicks=0,
            title=("Collapse " if cell_id in expanded else "Expand ") + names[cell_id],
            **{  # pyright: ignore[reportArgumentType] - Dash の型定義に ARIA kwargs がない。
                "aria-expanded": "true" if cell_id in expanded else "false"
            },
        )
        if cell_id.startswith("group:")
        else html.Div(
            names[cell_id], className="heatmap-child-label", title=names[cell_id]
        )
        for cell_id in cell_ids
    ]


def reference_links(record: FilteredRecord) -> str:
    """Open Targets へのリンクを、AgGrid の markdown セル用にまとめる。"""
    links: list[tuple[str, str]] = []
    target_id = record.get("target_id")
    if target_id:
        links.append(
            (
                "Open Targets target",
                "https://platform.opentargets.org/target/" + target_id,
            )
        )
    if record.get("canonical_drug_id"):
        links.append(
            (
                "Open Targets drug",
                "https://platform.opentargets.org/drug/" + record["canonical_drug_id"],
            )
        )
    if record.get("disease_id"):
        links.append(
            (
                "Open Targets disease",
                "https://platform.opentargets.org/disease/" + record["disease_id"],
            )
        )
    return "  \n".join(f"[{label}]({url})" for label, url in links) or "—"


EVIDENCE_COLUMNS: tuple[tuple[str, str, str | None], ...] = (
    (
        "canonical_drug",
        "Canonical drug",
        "The canonical drug (Open Targets parentMolecule) is the unit used for counting. Original drug forms that share a parentMolecule are merged into one row; without a parentMolecule, the original drug form is the canonical drug.",
    ),
    ("modality", "Modality", None),
    (
        "stage",
        "Canonical stage",
        "The highest stage across the original drug forms of the same canonical drug in this disease, used for filtering.",
    ),
    ("target", "Target", None),
    ("action_mechanism", "Action / mechanism", None),
    ("links", "Open Targets links", None),
)


def _joined(values: Iterable[str | None]) -> str:
    unique = list(dict.fromkeys(value for value in values if value))
    return " · ".join(unique) or "—"


def evidence_rows(records: Sequence[FilteredRecord]) -> list[dict[str, str]]:
    """元記録を canonical drug と標的の組ごとに 1 行へ畳む。"""
    groups: dict[tuple[str, str | None], list[FilteredRecord]] = {}
    for record in records:
        groups.setdefault(
            (record["canonical_drug_id"], record.get("target_id")), []
        ).append(record)
    rows: list[dict[str, str]] = []
    for members in groups.values():
        first = members[0]
        rows.append(
            {
                "canonical_drug": f"{first.get('canonical_drug') or '—'} ({first['canonical_drug_id']})",
                "modality": _joined(member.get("drug_type") for member in members),
                "stage": first.get("canonical_stage") or "—",
                "target": f"{first.get('target') or '—'} ({first.get('target_id') or '—'})",
                "action_mechanism": _joined(
                    f"{', '.join(member.get('action_types') or []) or '—'} / {member.get('mechanism') or '—'}"
                    for member in members
                ),
                "links": reference_links(first),
            }
        )
    return rows


def evidence_grid(rows: Sequence[dict[str, str]]) -> dag.AgGrid:
    """列ごとにソートと絞り込みができる、薬剤と標的の記録の表。"""
    column_defs: list[dict[str, object]] = []
    for field, header, help_text in EVIDENCE_COLUMNS:
        column: dict[str, object] = {"field": field, "headerName": header}
        if help_text:
            column["headerTooltip"] = help_text
        if field == "links":
            column.update(
                {
                    "cellRenderer": "markdown",
                    "linkTarget": "_blank",
                    "sortable": False,
                    "filter": False,
                    "floatingFilter": False,
                }
            )
        column_defs.append(column)
    return dag.AgGrid(
        id="source-records-grid",
        columnDefs=column_defs,
        rowData=list(rows),
        defaultColDef={
            "sortable": True,
            "filter": "agTextColumnFilter",
            "floatingFilter": True,
            "resizable": True,
            "wrapText": True,
            "autoHeight": True,
            "minWidth": 120,
        },
        columnSize="responsiveSizeToFit",
        dashGridOptions={
            "pagination": True,
            "paginationPageSize": SOURCE_PAGE_SIZE,
            "paginationPageSizeSelector": False,
            "animateRows": False,
            "domLayout": "autoHeight",
            "tooltipShowDelay": 300,
            "suppressCellFocus": True,
        },
        className="ag-theme-alpine source-records-grid",
    )
