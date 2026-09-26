"""Dash アプリの画面構成。"""

import math
from collections.abc import Collection, Iterable, Sequence
from datetime import datetime, timedelta, timezone
from typing import cast

from dash import dcc, html
from dash.development.base_component import Component

from backend import aggregation as atlas
from backend.disease_catalog import ordered_disease_ids
from backend.models import Snapshot, SummaryRow
from backend.ui.components import disease_selector, info_tip
from backend.ui.config import (
    DEFAULT_EXPRESSION_THRESHOLD,
    DEFAULT_SPECIFICITY_THRESHOLD,
    METHOD_LABELS,
    PARAMETER_IDS,
    SOURCE_PAGE_SIZE,
)


def choose_defaults(
    items: Sequence[tuple[str, str]],
    terms: Iterable[str],
    limit: int,
    preferred_ids: Collection[str] | None = None,
) -> list[str]:
    """名前と実データの有無から、存在する項目だけを初期選択する。"""
    selected: list[str] = []
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


def format_data_version(value: object) -> str:
    """Open Targets の版を画面で使える短い文字列にする。"""
    if not isinstance(value, dict):
        return str(value or "")
    version = cast(dict[str, object], value)
    return ".".join(
        str(part)
        for part in (
            version.get("year"),
            version.get("month"),
            version.get("iteration"),
        )
        if part is not None
    )


def detail_panel(
    rows: Sequence[SummaryRow],
    selection: str | None,
    snapshot: Snapshot | None = None,
    *,
    modality: str = "all",
    stage: str = "phase3",
    threshold: float = DEFAULT_EXPRESSION_THRESHOLD,
    method: str = "fixed",
    specificity: float = DEFAULT_SPECIFICITY_THRESHOLD,
) -> html.Div:
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
                    "Refresh the data with pixi run refresh to view evidence details.",
                    className="empty-note",
                ),
            ]
        )
    filtered = [
        record
        for record in atlas.filtered_records(snapshot, modality, stage)
        if record["disease_id"] == row["disease_id"]
    ]
    source_context: dict[str, str | float] = {
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
                    f"Drug–target records for {row['disease']}",
                    info_tip(
                        "source-records",
                        "source records",
                        "Records match the selected disease and applied drug filters. Each row is an original drug–target pair, so a drug with several targets appears on several rows. Expression thresholds do not filter this table. Ten rows are shown per page; — means unavailable information.",
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
                    "Relative target expression by cell type",
                    info_tip(
                        "expression",
                        "target expression",
                        "Healthy reference expression for known targets of the selected disease's filtered drugs. Click a group name to show its source cell types. The same reference data are used for every disease; differences between diseases reflect their drug and target sets, not expression in patients.",
                    ),
                    html.Small(
                        "○ Source cell meets expression rule",
                        className="expression-marker-key",
                    ),
                    html.Small(
                        [
                            "Group means",
                            info_tip(
                                "expression-groups",
                                "group mean expression",
                                "Each group shows the arithmetic mean of available source-cell donor-median CPM values, with equal weight per cell type. Hover shows how many cell types have data. All missing means no group value. Group means are display-only: they never enter reference medians, z-score reference statistics, clustering, expression rules or comparison counts. White dots apply only to individual source cells.",
                            ),
                            "Color & zeros",
                            info_tip(
                                "expression-scale",
                                "expression colors and zeros",
                                "Red is higher and blue lower expression relative to the same target across all source cell types, not relative to other targets. Values are log2(1 + CPM), standardized using source-cell values only; group means use that same transformation. Targets are ordered by source-cell patterns. Expanding groups changes neither order nor color scale. A printed 0 marks missing expression; hover identifies it. Neutral color at z = 0 means reference-average log expression, not absence of expression.",
                            ),
                        ],
                        className="expression-marker-key",
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


def unavailable_layout(error: Exception | None = None) -> html.Main:
    """データ未取得または再取得が必要な状態を説明する。"""
    title = "No data loaded yet" if error is None else "Data refresh required"
    return html.Main(
        [
            html.Section(
                [
                    html.P("AUTOIMMUNE DRUG–CELL ATLAS", className="eyebrow"),
                    html.H1(title),
                    html.P(
                        "Run the command below to download data, then restart the app."
                    ),
                    html.Code("pixi run refresh"),
                    *([html.P(str(error), className="empty-note")] if error else []),
                ],
                className="unavailable-card",
            )
        ],
        className="shell unavailable",
    )


def dashboard_layout(snapshot: Snapshot) -> html.Main:
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

    def control(
        label: str,
        component_id: str,
        component: dcc.Dropdown,
        description: str | None = None,
    ) -> html.Div:
        title: Component = html.Label(label, htmlFor=component_id)
        if description:
            title = html.Div(
                [title, info_tip(component_id, label.lower(), description)],
                className="label-help",
            )
        return html.Div([title, component], className="control")

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
            strict=True,
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
                                "Open Targets",
                                href=source,
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
                                    "Drug targets across diseases and cell types",
                                    info_tip(
                                        "overview",
                                        "the atlas",
                                        "Snapshot totals above the settings describe all saved data, not the filtered results. Compare which drug targets meet an expression rule in healthy reference cells across diseases. Disease differences reflect eligible drug and target sets, not disease-specific expression. This does not establish treatment efficacy; disease records may include symptom or comorbidity treatment.",
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
                                    html.Span(
                                        "Drug–target records", className="meta-label"
                                    ),
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
                        **{  # pyright: ignore[reportArgumentType] - Dash の型定義に ARIA kwargs がない。
                            "aria-label": "Snapshot data"
                        },
                    ),
                ],
                className="hero",
            ),
            html.Section(
                [
                    html.H2(
                        [
                            "Settings panel",
                            info_tip(
                                "settings",
                                "applying settings",
                                "Edit the settings, then click Update to apply them together. Until then, both heatmaps and disease details keep the previous settings. Group expansion, detail disease selection and table pages update immediately.",
                            ),
                        ]
                    ),
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
                                                        "Filters by the highest recorded stage across original forms of a canonical drug within each disease. An earlier-stage original record may remain in the table. Approval reached can include withdrawn drugs and does not guarantee current approval.",
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
                                        "modality",
                                        dcc.Dropdown(
                                            id="modality",
                                            options=modality_options,
                                            value="all",
                                            clearable=False,
                                        ),
                                        "Filters original drug records by the modality recorded in Open Targets. All includes every modality; Unknown means no mapped modality.",
                                    ),
                                ],
                                className="filter-group",
                            ),
                            html.Section(
                                [
                                    html.H3("Expression criteria"),
                                    control(
                                        "Expression rule",
                                        "method",
                                        dcc.Dropdown(
                                            id="method",
                                            options=[
                                                {"label": label, "value": value}
                                                for value, label in METHOD_LABELS.items()
                                            ],
                                            value="specificity",
                                            clearable=False,
                                        ),
                                        "Fixed CPM uses Minimum CPM alone. Target-relative median also requires expression at or above that target's median across all source cell types; incomplete reference data make that median unavailable. CELLEX also requires the selected specificity threshold. Group averages never enter these rules.",
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
                                                        f"Minimum donor-median CPM for an individual source cell type. Applies to every rule; group averages are not tested. Blank uses {DEFAULT_EXPRESSION_THRESHOLD:g} CPM.",
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
                                                        f"Requires the source-cell CELLEX score to be at least this value, in addition to Minimum CPM. Used only with the CELLEX rule; missing scores remain unresolved. Blank uses {DEFAULT_SPECIFICITY_THRESHOLD:g}.",
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
                                                        "Distinct targets counts unique genes meeting the expression rule. Canonical drugs counts each parent drug once when any known target meets the rule. A group counts a target or drug once if any member cell type qualifies.",
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
                                                        "Count shows qualifying targets or drugs. Percent divides by all known targets, or drugs with known targets, for that disease and the applied drug filters. The denominator is fixed across cell types. These are not cell proportions or probabilities of treatment benefit.",
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
                                    "Compare diseases by cell type",
                                    info_tip(
                                        "comparison",
                                        "cell-type comparison",
                                        "Compare qualifying targets or drugs across diseases. Group rows count distinct targets or drugs across their members, not average expression. Zero-colored cells can mean no qualifying evidence or an unavailable value; hover explains which. A ≥ in the hover means the result is a lower bound. Color scales are separate for targets and drugs and can change when filters change.",
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
                        "Click a group name to expand · Click a heatmap cell for disease details · Hover to distinguish zero from unavailable",
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
                                    "Drug records and target expression",
                                    info_tip(
                                        "selection",
                                        "selection details",
                                        "Click a comparison heatmap cell or choose a disease below. Details update immediately within the last applied settings. Drug records show the underlying evidence; the expression map includes all known targets of those drugs, even when they do not meet the expression rule.",
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
                                "detail-disease",
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
