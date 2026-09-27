"""遺伝子ページの画面構成。"""

from collections.abc import Sequence
from datetime import datetime, timedelta, timezone

import dash_ag_grid as dag  # pyright: ignore[reportMissingTypeStubs] - dash-ag-grid に型スタブがない。
from dash import dcc, html
from dash.development.base_component import Component

from autoimmune_atlas import genetics as gene_data
from autoimmune_atlas.disease_catalog import ordered_disease_ids
from autoimmune_atlas.models import (
    GeneAssociation,
    GeneticsSnapshot,
    Snapshot,
    SummaryRow,
)
from autoimmune_atlas.ui.components import disease_selector, info_tip, page_nav
from autoimmune_atlas.ui.config import (
    DEFAULT_DISEASE_TERMS,
    DEFAULT_EXPRESSION_THRESHOLD,
    DEFAULT_SCORE_THRESHOLD,
    DEFAULT_SPECIFICITY_THRESHOLD,
    GENETICS_PARAMETER_IDS,
    METHOD_LABELS,
    SCORE_FLOOR,
    SOURCE_PAGE_SIZE,
)
from autoimmune_atlas.ui.layout import choose_defaults, format_data_version

OPEN_TARGETS = "https://platform.opentargets.org/"

GENE_COLUMNS: tuple[tuple[str, str, str | None], ...] = (
    ("gene", "Gene", None),
    (
        "score",
        "Genetic association score",
        "The Open Targets datatype score for genetic association, combined across the genetic datasources. It is not a probability of causality or of treatment benefit.",
    ),
    ("links", "Open Targets", None),
)


def gene_rows(
    genes: Sequence[GeneAssociation], disease_id: str, datasources: Sequence[str]
) -> list[dict[str, str]]:
    """閾値以上の遺伝子 1 件を 1 行にし、datasource ごとのスコアを列にする。"""
    rows: list[dict[str, str]] = []
    for gene in genes:
        row = {
            "gene": f"{gene['target']} ({gene['target_id']})",
            "score": f"{gene['score']:.3f}",
            "links": f"[Evidence](https://platform.opentargets.org/evidence/{gene['target_id']}/{disease_id})",
        }
        for datasource in datasources:
            value = gene["datasource_scores"].get(datasource)
            row[datasource] = "—" if value is None else f"{value:.3f}"
        rows.append(row)
    return rows


def gene_grid(rows: Sequence[dict[str, str]], datasources: Sequence[str]) -> dag.AgGrid:
    """列ごとにソートと絞り込みができる、関連遺伝子の表。"""
    column_defs: list[dict[str, object]] = []
    for field, header, help_text in GENE_COLUMNS[:2]:
        column: dict[str, object] = {"field": field, "headerName": header}
        if help_text:
            column["headerTooltip"] = help_text
        if field == "score":
            column["sort"] = "desc"
        column_defs.append(column)
    for datasource in datasources:
        column_defs.append(
            {
                "field": datasource,
                "headerName": datasource.replace("_", " "),
                "headerTooltip": f"Datasource score from {datasource}. — means no evidence from this datasource.",
            }
        )
    column_defs.append(
        {
            "field": "links",
            "headerName": "Open Targets",
            "cellRenderer": "markdown",
            "linkTarget": "_blank",
            "sortable": False,
            "filter": False,
            "floatingFilter": False,
        }
    )
    return dag.AgGrid(
        id="genetics-gene-grid",
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


def genetics_detail_panel(
    rows: Sequence[SummaryRow],
    selection: str | None,
    genetics: GeneticsSnapshot,
    *,
    score_threshold: float,
    threshold: float,
    method: str,
    specificity: float,
) -> html.Div:
    """選択した疾患の関連遺伝子の表と連続発現を表示する。"""
    if selection is None:
        return html.Div(
            "Select a disease in the comparison charts or the selector above.",
            className="empty-note",
        )
    row = next((item for item in rows if item["disease_id"] == selection), None)
    if row is None:
        return html.Div(
            "The previous disease is outside the current filters. Select a visible disease.",
            className="empty-note",
        )
    genes = gene_data.genes_for_disease(genetics, row["disease_id"], score_threshold)
    source_context: dict[str, str | float] = {
        "disease_id": row["disease_id"],
        "score": score_threshold,
        "threshold": threshold,
        "method": method,
        "specificity": specificity,
    }
    # 遺伝子が無くても発現欄は残し、発現の callback が空の旨を表示する。
    table: Component = (
        gene_grid(
            gene_rows(genes, row["disease_id"], genetics["datasources"]),
            genetics["datasources"],
        )
        if genes
        else html.P(
            "No genes at or above the score threshold for this disease.",
            className="empty-note",
        )
    )
    return html.Div(
        [
            html.H3(
                [
                    f"Genetically associated genes for {row['disease']}",
                    info_tip(
                        "genetics-source-records",
                        "genetically associated genes",
                        "This table lists the genes whose genetic association score for the selected disease is at or above the applied threshold. Each row is one gene. Datasource columns show the score from each genetic datasource; — means no evidence from that datasource. Expression rules do not filter this table. Ten rows are shown per page. Click a column header to sort, or type in the boxes under the headers to filter. Evidence opens the Open Targets evidence page for the gene and disease.",
                    ),
                ]
            ),
            html.Div(
                [
                    dcc.Store(id="genetics-source-context", data=source_context),
                    table,
                ],
                className="source-records",
            ),
            html.H3(
                [
                    "Relative gene expression by cell type",
                    info_tip(
                        "genetics-expression",
                        "gene expression",
                        "This chart shows healthy reference expression for every gene at or above the applied genetic association score. The heatmap and the dot plot show the same data. In the dot plot, dot area represents CELLEX specificity and color is the target-wise z-score; zero and very small specificity values use 3 px dots, and an empty cell means missing CPM or CELLEX. The same reference data are used for every disease, so differences between diseases reflect their associated gene sets, not expression in patients. Click a cell group name to show its source cell types.",
                    ),
                    html.Fieldset(
                        [
                            html.Legend("Chart"),
                            dcc.RadioItems(
                                id="genetics-expression-chart-type",
                                options=[
                                    {"label": "Heatmap", "value": "heatmap"},
                                    {"label": "Dot plot", "value": "dot"},
                                ],
                                value="heatmap",
                                inline=True,
                                persistence=True,
                                persistence_type="session",
                            ),
                        ],
                        className="control radio-control",
                    ),
                    html.Small(
                        [
                            "Group means",
                            info_tip(
                                "genetics-expression-groups",
                                "group means",
                                "A group mean is the equal-weight arithmetic mean across the source cell types of a cell group that have values: donor-median CPM, the median CPM across donors, in both charts and CELLEX specificity for dot size. The hover shows how many source cell types have values for each measure. A cell group whose values are all missing has no group mean. Group means are display-only: they never enter reference medians, z-score reference statistics, clustering, expression rules or comparison counts.",
                            ),
                        ],
                        className="expression-marker-key",
                    ),
                    html.Small(
                        [
                            "Color",
                            info_tip(
                                "genetics-expression-scale",
                                "expression colors",
                                "The color is the target-wise z-score, shared by both charts: red is higher and blue lower expression relative to the same target across all source cell types, not relative to other targets. Values are log2(1 + CPM), standardized over source cell types only; group means use the same transformation. Targets are ordered by their patterns across source cell types. Expanding a cell group changes neither the order nor the color scale. The neutral color at z = 0 means reference-average log expression, not absence of expression.",
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
                            id="genetics-expression-heatmap",
                            config={"displaylogo": False, "responsive": True},
                        ),
                        html.Div(
                            id="genetics-expression-row-controls",
                            className="heatmap-row-controls",
                            style={"top": "110px", "width": "280px"},
                        ),
                    ],
                    className="expandable-heatmap",
                ),
                className="graph-scroll expression-graph",
            ),
            html.Small(
                "○ Cell type meets expression rule · Gray: missing expression",
                id="genetics-expression-heatmap-key",
                className="expression-marker-key",
                style={"display": "inline-flex"},
            ),
            html.Div(
                id="genetics-expression-dot-key",
                className="dot-key expression-dot-key",
                style={"display": "none"},
            ),
        ]
    )


def genetics_unavailable_page(reason: str) -> html.Main:
    """genetics.json が無いか版が合わないときの案内。"""
    return html.Main(
        [
            page_nav("genetics"),
            html.Section(
                [
                    html.H1("Genetic associations are not available"),
                    html.P(reason),
                    html.P(
                        "Refresh the genetic association data with pixi run refresh-genetics after pixi run refresh, then restart the app."
                    ),
                ],
                className="panel",
            ),
        ],
        className="shell",
        id="genetics-page",
    )


def _labelled(
    label: str, component_id: str, component: Component, description: str
) -> html.Div:
    """ラベルと ⓘ を付けた設定欄。"""
    return html.Div(
        [
            html.Div(
                [
                    html.Label(label, htmlFor=component_id),
                    info_tip(component_id, label.lower(), description),
                ],
                className="label-help",
            ),
            component,
        ],
        className="control",
    )


def _meta(label: str, value: Component | str) -> html.Div:
    """概要欄の 1 項目。"""
    return html.Div([html.Span(label, className="meta-label"), value])


def genetics_page(snapshot: Snapshot, genetics: GeneticsSnapshot) -> html.Main:
    """遺伝子ページの初期画面を作る。"""
    default_diseases = choose_defaults(
        [(d["id"], d["name"]) for d in snapshot["diseases"]],
        DEFAULT_DISEASE_TERMS,
        10,
        set(genetics["associations"]),
    )
    retrieved_at = (
        datetime.fromisoformat(genetics["retrieved_at"])
        .astimezone(timezone(timedelta(hours=9)))
        .strftime("%Y-%m-%d %H:%M JST")
    )
    loaded = sum(d["id"] in genetics["associations"] for d in snapshot["diseases"])
    applied = dict(
        zip(
            GENETICS_PARAMETER_IDS,
            (
                "percent",
                DEFAULT_SCORE_THRESHOLD,
                "specificity",
                DEFAULT_EXPRESSION_THRESHOLD,
                DEFAULT_SPECIFICITY_THRESHOLD,
                ordered_disease_ids(snapshot, default_diseases),
            ),
            strict=True,
        )
    )
    return html.Main(
        [
            page_nav("genetics"),
            html.Header(
                [
                    html.Div(
                        [
                            html.P(
                                "GENETIC ASSOCIATIONS · HEALTHY REFERENCE EXPRESSION",
                                className="eyebrow",
                            ),
                            html.H1(
                                [
                                    "Genetically associated genes across diseases and cell types",
                                    info_tip(
                                        "genetics-overview",
                                        "the atlas",
                                        "This atlas compares which genetically associated genes meet an expression rule in healthy reference cell types across diseases. Differences between diseases reflect their associated gene sets, not disease-specific expression. The genetic association score is combined by Open Targets from several genetic datasources and does not establish causality or treatment efficacy. The snapshot totals above the settings describe all saved data, not the filtered results.",
                                    ),
                                ]
                            ),
                        ]
                    ),
                    html.Section(
                        [
                            _meta(
                                "Disease terms",
                                html.Strong(
                                    f"{loaded} loaded / {len(snapshot['diseases'])} total"
                                ),
                            ),
                            _meta(
                                "Gene–disease associations (score ≥ 0.1)",
                                html.Strong(
                                    str(
                                        sum(
                                            len(g)
                                            for g in genetics["associations"].values()
                                        )
                                    )
                                ),
                            ),
                            _meta(
                                "Genes with expression",
                                html.Strong(
                                    str(
                                        len(
                                            set(snapshot["expression"])
                                            | set(genetics["expression"])
                                        )
                                    )
                                ),
                            ),
                            _meta(
                                "Retrieved at",
                                html.Strong(
                                    retrieved_at, title=genetics["retrieved_at"]
                                ),
                            ),
                            _meta(
                                "Expression reference", html.Strong("Tabula Sapiens")
                            ),
                            _meta(
                                "Data source",
                                html.A(
                                    f"Open Targets {format_data_version(genetics['data_version'])}".strip(),
                                    href=OPEN_TARGETS,
                                    target="_blank",
                                    rel="noreferrer",
                                ),
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
                            "Settings",
                            info_tip(
                                "genetics-settings",
                                "settings",
                                "These settings define the comparison and the disease details. Until you click Update, both charts and the disease details keep the previous settings. Cell group expansion, the detail disease and table pages update immediately. Edit the settings, then click Update to apply them together.",
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
                                                snapshot,
                                                default_diseases,
                                                prefix="genetics-",
                                            ),
                                        ],
                                        className="control-grid scope-controls",
                                    ),
                                ],
                                className="filter-group scope-group",
                            ),
                            html.Section(
                                [
                                    html.H3("Genetic evidence"),
                                    _labelled(
                                        "Genetic association score (≥)",
                                        "genetics-score",
                                        dcc.Input(
                                            id="genetics-score",
                                            type="number",
                                            min=SCORE_FLOOR,
                                            max=1,
                                            step=0.05,
                                            value=DEFAULT_SCORE_THRESHOLD,
                                        ),
                                        "The score threshold keeps genes whose Open Targets genetic association score for the disease is at or above this value. Scores below 0.1 are not stored, so the threshold cannot go below 0.1. Empty or out-of-range values use 0.5.",
                                    ),
                                ],
                                className="filter-group",
                            ),
                            html.Section(
                                [
                                    html.H3("Expression criteria"),
                                    _labelled(
                                        "Expression rule",
                                        "genetics-method",
                                        dcc.Dropdown(
                                            id="genetics-method",
                                            options=[
                                                {"label": label, "value": value}
                                                for value, label in METHOD_LABELS.items()
                                            ],
                                            value="specificity",
                                            clearable=False,
                                        ),
                                        "The expression rule decides whether a target counts as expressed in a cell type. Fixed CPM uses Minimum CPM alone. Target-relative median also requires expression at or above that target's median across all cell types in the reference; when reference data are incomplete, that median is missing and the rule is unresolved. CELLEX also requires at least the selected CELLEX specificity. Group means never enter these rules.",
                                    ),
                                    # Dash 4.4.1 は max を省略すると増減時に NaN になるため、上限なしを明示する。
                                    _labelled(
                                        "Minimum CPM (≥)",
                                        "genetics-threshold",
                                        dcc.Input(
                                            id="genetics-threshold",
                                            type="number",
                                            min=0,
                                            max=None,
                                            step=0.1,
                                            value=DEFAULT_EXPRESSION_THRESHOLD,
                                        ),
                                        f"Minimum CPM is the lowest donor-median CPM, the median CPM across donors, that a cell type must reach. It applies to every expression rule; group means are not tested. Empty uses {DEFAULT_EXPRESSION_THRESHOLD:g} CPM.",
                                    ),
                                    _labelled(
                                        "CELLEX specificity (≥)",
                                        "genetics-specificity",
                                        dcc.Input(
                                            id="genetics-specificity",
                                            type="number",
                                            min=0,
                                            max=1,
                                            step=0.05,
                                            value=DEFAULT_SPECIFICITY_THRESHOLD,
                                            disabled=False,
                                        ),
                                        f"CELLEX specificity is the lowest CELLEX score that a cell type must reach in addition to Minimum CPM. It is used only with the CELLEX expression rule, and a missing score leaves the rule unresolved. Empty uses {DEFAULT_SPECIFICITY_THRESHOLD:g}.",
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
                                                    "Measure (heatmap)",
                                                    info_tip(
                                                        "genetics-measure",
                                                        "measure",
                                                        "Count is the number of qualifying genes. Percent divides that number by all genes at or above the score threshold for the disease; the denominator is the same for every cell type.",
                                                    ),
                                                ]
                                            ),
                                            dcc.RadioItems(
                                                id="genetics-measure",
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
                                        id="genetics-update-button",
                                        type="button",
                                        n_clicks=0,
                                    ),
                                    html.Span(
                                        "Settings applied",
                                        id="genetics-update-status",
                                        role="status",
                                    ),
                                ],
                                className="update-actions",
                            ),
                        ],
                        className="filter-groups",
                    ),
                    dcc.Store(id="genetics-applied-parameters", data=applied),
                    dcc.Store(id="genetics-expanded-cell-groups", data=[]),
                    dcc.Store(id="genetics-expanded-expression-groups", data=[]),
                ],
                className="panel controls",
            ),
            html.Section(
                [
                    html.Div(
                        [
                            html.H2(
                                [
                                    "Comparison by cell type",
                                    info_tip(
                                        "genetics-comparison",
                                        "comparison by cell type",
                                        "The comparison counts genes at or above the score threshold that meet the expression rule for each disease and cell type. A cell group row counts distinct genes across its cell types, not mean expression. The heatmap shows the selected measure. In the dot plot, dot area uses a compressed Count scale and color shows Percent; an empty cell means zero or missing. A ≥ in the hover marks a lower bound. The scale can change when settings change. Click a cell group name to expand it, click a cell or mark for disease details, or hover for counts and evidence status.",
                                    ),
                                ]
                            ),
                            html.Fieldset(
                                [
                                    html.Legend("Chart"),
                                    dcc.RadioItems(
                                        id="genetics-chart-type",
                                        options=[
                                            {"label": "Heatmap", "value": "heatmap"},
                                            {"label": "Dot plot", "value": "dot"},
                                        ],
                                        value="heatmap",
                                        inline=True,
                                    ),
                                ],
                                className="control radio-control",
                            ),
                            html.Div(
                                id="genetics-matrix-note",
                                className="matrix-note",
                                role="status",
                            ),
                        ],
                        className="section-heading",
                    ),
                    html.Div(
                        [
                            html.Article(
                                [
                                    html.H3("Genes"),
                                    html.Div(
                                        html.Div(
                                            [
                                                dcc.Graph(
                                                    id="genetics-heatmap",
                                                    style={"minWidth": "600px"},
                                                    config={
                                                        "displaylogo": False,
                                                        "responsive": True,
                                                    },
                                                ),
                                                html.Div(
                                                    id="genetics-row-controls",
                                                    className="heatmap-row-controls",
                                                ),
                                            ],
                                            className="expandable-heatmap",
                                        ),
                                        className="graph-scroll",
                                    ),
                                    html.Div(
                                        id="genetics-dot-key",
                                        className="dot-key",
                                        style={"display": "none"},
                                    ),
                                ],
                                id="genetics-heatmap-panel",
                            ),
                        ],
                        className="heatmap-stack",
                    ),
                ],
                className="panel matrix-panel",
                id="genetics-comparison",
            ),
            html.Section(
                [
                    html.Div(
                        [
                            html.H2(
                                [
                                    "Disease details",
                                    info_tip(
                                        "genetics-selection",
                                        "disease details",
                                        "Disease details show the genetically associated genes and their expression for one disease. They update immediately within the last applied settings. The table lists every gene at or above the score threshold, and the expression chart includes all of them, even when they do not meet the expression rule. Click a cell or mark in the comparison chart, or choose a disease below.",
                                    ),
                                ]
                            )
                        ],
                        className="section-heading detail-heading",
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Label(
                                        "Disease", htmlFor="genetics-detail-disease"
                                    ),
                                    dcc.Dropdown(
                                        id="genetics-detail-disease", clearable=False
                                    ),
                                ],
                                className="control",
                            )
                        ],
                        className="control-grid selectors",
                    ),
                    html.Div(id="genetics-details", className="details"),
                ],
                className="panel",
                id="genetics-evidence",
            ),
        ],
        className="shell",
        id="genetics-page",
    )
