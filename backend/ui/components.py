"""Dash の再利用可能な画面部品。"""

from collections.abc import Collection, Iterable, Mapping, Sequence
from typing import Literal, TypedDict

from dash import dcc, html
from dash.development.base_component import Component

from backend import aggregation as atlas
from backend.disease_catalog import disease_catalog, ordered_disease_ids
from backend.models import (
    CatalogDisease,
    DiseaseCatalogGroup,
    DiseaseFamily,
    FilteredRecord,
    Snapshot,
)


class DiseaseSection(TypedDict):
    """同じ選択欄に表示する疾患の区分。"""

    id: str
    diseases: list[CatalogDisease]
    node: str | None
    kind: Literal["self", "children"]


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


def disease_checklist_sections(family: DiseaseFamily) -> list[DiseaseSection]:
    """疾患本体と各階層の選択欄を、描画と同期で同じ範囲に分ける。

    子を持つ語は、ファミリーと同じ形にする。
    本体だけの section（kind が self）と、子のうち葉だけを入れた section（kind が
    children）を持ち、子を持つ子は同じ形で入れ子になる。
    先頭の section はファミリーの本体と、親を持たない葉である。
    """
    top, children = _disease_tree(family)
    members = {d["id"]: d for d in family["diseases"]}
    root = family["id"] if family["id"] in members else None

    def section_id(prefix: str, node: str) -> str:
        suffix = "" if node == root else f"-{node}"
        return f"{prefix}-{family['id']}{suffix}"

    sections: list[DiseaseSection] = [
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


def disease_selector(snapshot: Snapshot, selected: Iterable[str] | None) -> html.Div:
    """検索欄と、群から開けるチェック欄を同じ選択へ結び付ける。"""
    catalog: list[DiseaseCatalogGroup] = disease_catalog(snapshot)
    selected_ids = set(selected or ())
    options: list[dict[str, str]] = []
    groups: list[Component] = []
    for group in catalog:
        families: list[Component] = []
        for family in group["families"]:
            choices = [
                {"label": d["name"], "value": d["id"]} for d in family["diseases"]
            ]
            options.extend(choices)
            sections = disease_checklist_sections(family)
            top, children = _disease_tree(family)
            root = (
                family["id"]
                if family["id"] in {d["id"] for d in family["diseases"]}
                else None
            )
            by_kind = {(s["kind"], s["node"]): s for s in sections}
            names = {d["id"]: d["name"] for d in family["diseases"]}

            def descendant_count(node: str) -> int:
                return sum(1 + descendant_count(c) for c in children[node])

            def checklist(section: DiseaseSection) -> dcc.Checklist:
                return dcc.Checklist(
                    id=section["id"],
                    options=[
                        {"label": d["name"], "value": d["id"]}
                        for d in section["diseases"]
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
                inner.extend(render_node(c) for c in children[node] if children[c])
                return html.Details(
                    [
                        html.Summary(
                            f"{names[node]} details ({descendant_count(node)} terms)"
                        ),
                        html.Div(inner, className="disease-families"),
                    ]
                )

            def render_node(node: str) -> html.Details:
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

            top_items: list[Component] = [checklist(sections[0])]
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
                        "Search by name or browse disease groups. Each checkbox selects one disease term; selecting a parent does not select its children. Click Update to apply your selection. These groups are navigation aids, not a diagnostic classification.",
                    ),
                ],
                className="label-help",
            ),
            dcc.Dropdown(
                id="diseases",
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
                "type": "expression-cell-toggle"
                if kind == "expression"
                else "heatmap-cell-toggle",
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


def reference_links(record: FilteredRecord) -> html.Div | str:
    links: list[html.A] = []
    target_id = record.get("target_id")
    if target_id:
        links.append(
            html.A(
                "Open Targets target",
                href="https://platform.opentargets.org/target/" + target_id,
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
    children: list[Component] = []
    for index, item in enumerate(links):
        if index:
            children.append(html.Br())
        children.append(item)
    return html.Div(children) if children else "—"


def evidence_table(records: Sequence[FilteredRecord]) -> html.Div:
    headers = (
        "Canonical drug",
        "Original drug",
        "Modality",
        "Original / canonical stage",
        "Target",
        "Action / mechanism",
        "Open Targets links",
    )
    column_help = {
        "Canonical drug": (
            "canonical-drug",
            "Drug used for counting: original forms sharing Open Targets' parentMolecule are counted once. If no parent is recorded, the original drug is used.",
        ),
        "Original drug": (
            "original-drug",
            "Drug name and ID in the source record, which may describe a salt or another form. Original does not mean originator brand. Forms of the same canonical drug can have different targets, stages or references.",
        ),
        "Original / canonical stage": (
            "record-stage",
            "Left: stage of this original drug in this disease. Right: highest stage across original forms of the same canonical drug in this disease, used for filtering.",
        ),
    }
    body: list[html.Tr] = []
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
                    html.Td(reference_links(record)),
                ]
            )
        )
    return html.Div(
        html.Table(
            [
                html.Thead(
                    html.Tr(
                        [
                            html.Th(
                                [
                                    item,
                                    info_tip(
                                        column_help[item][0],
                                        item.lower(),
                                        column_help[item][1],
                                    ),
                                ]
                            )
                            if item in column_help
                            else html.Th(item)
                            for item in headers
                        ]
                    )
                ),
                html.Tbody(body),
            ]
        ),
        className="table-scroll",
    )
