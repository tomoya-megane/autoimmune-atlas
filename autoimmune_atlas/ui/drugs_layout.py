"""Dash アプリの画面構成。"""

from collections.abc import Sequence
from datetime import datetime, timedelta, timezone

from dash import dcc, html
from dash.development.base_component import Component

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas.disease_catalog import ordered_disease_ids
from autoimmune_atlas.models import Snapshot, SummaryRow
from autoimmune_atlas.ui.components import (
    choose_defaults,
    disease_selector,
    evidence_grid,
    evidence_rows,
    format_data_version,
    info_tip,
    page_nav,
)
from autoimmune_atlas.ui.config import (
    DEFAULT_DISEASE_TERMS,
    DEFAULT_EXPRESSION_THRESHOLD,
    DEFAULT_SPECIFICITY_THRESHOLD,
    METHOD_LABELS,
    PARAMETER_IDS,
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
            "Choose a disease in the selector above.",
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
    pairs = evidence_rows(filtered)
    source_context: dict[str, str | float] = {
        "disease_id": row["disease_id"],
        "modality": modality,
        "stage": stage,
        "threshold": threshold,
        "method": method,
        "specificity": specificity,
    }
    return html.Div(
        [
            html.H3(
                [
                    f"Drug–target records for {row['disease']}",
                    info_tip(
                        "source-records",
                        "drug–target records",
                        "This table lists drug–target records for the selected disease under the applied drug filters. Each row is a canonical drug–target pair; original drug forms such as salts are merged, and a canonical drug with several targets appears on several rows. Expression rules do not filter this table. Ten rows are shown per page. — means missing. Click a column header to sort, or type in the boxes under the headers to filter.",
                    ),
                ]
            ),
            html.Div(
                [
                    dcc.Store(id="source-context", data=source_context),
                    evidence_grid(pairs),
                ],
                className="source-records",
            ),
            html.H3(
                [
                    "Relative target expression by cell type",
                    info_tip(
                        "expression",
                        "target expression",
                        "This chart shows healthy reference expression for every target of the selected disease's filtered canonical drugs. The heatmap and the dot plot show the same data. In the dot plot, dot area represents CELLEX specificity and color is the target-wise z-score; zero and very small specificity values use 3 px dots, and an empty cell means missing CPM or CELLEX. The same reference data are used for every disease, so differences between diseases reflect their drug and target sets, not expression in patients. Click a cell group name to show its source cell types.",
                    ),
                    html.Fieldset(
                        [
                            html.Legend("Chart"),
                            dcc.RadioItems(
                                id="expression-chart-type",
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
                                "expression-groups",
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
                                "expression-scale",
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
            html.Small(
                "○ Cell type meets expression rule · Gray: missing expression",
                id="expression-heatmap-key",
                className="expression-marker-key",
                style={"display": "inline-flex"},
            ),
            html.Div(
                id="expression-dot-key",
                className="dot-key expression-dot-key",
                style={"display": "none"},
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
                    html.P(
                        "DRUG TARGETS · HEALTHY REFERENCE EXPRESSION",
                        className="eyebrow",
                    ),
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
        id="drug-page",
    )


def dashboard_layout(snapshot: Snapshot) -> html.Main:
    """schema 2 のデータから dashboard の初期画面を作る。"""
    diseases = [(row["id"], row["name"]) for row in snapshot["diseases"]]
    default_diseases = choose_defaults(
        diseases,
        DEFAULT_DISEASE_TERMS,
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
            page_nav("drugs"),
            html.Header(
                [
                    html.Div(
                        [
                            html.P(
                                "DRUG TARGETS · HEALTHY REFERENCE EXPRESSION",
                                className="eyebrow",
                            ),
                            html.H1(
                                [
                                    "Drug targets across diseases and cell types",
                                    info_tip(
                                        "overview",
                                        "the atlas",
                                        "This atlas compares which drug targets meet an expression rule in healthy reference cell types across diseases. Differences between diseases reflect their eligible drug and target sets, not disease-specific expression. The snapshot totals above the settings describe all saved data, not the filtered results. The atlas does not establish treatment efficacy, and disease records may include treatments for symptoms or comorbidities.",
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
                            "Settings",
                            info_tip(
                                "settings",
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
                                                        "The clinical stage filter uses the highest recorded stage across the original drug forms of a canonical drug within each disease. A record of an original drug form at an earlier stage may remain in the table. Approval reached can include withdrawn drugs and does not guarantee current approval.",
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
                                        "The drug modality filter uses the modality that Open Targets records for each original drug form. All includes every modality; Unknown means no mapped modality.",
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
                                        "The expression rule decides whether a target counts as expressed in a cell type. Fixed CPM uses Minimum CPM alone. Target-relative median also requires expression at or above that target's median across all cell types in the reference; when reference data are incomplete, that median is missing and the rule is unresolved. CELLEX also requires at least the selected CELLEX specificity. Group means never enter these rules.",
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
                                                        f"Minimum CPM is the lowest donor-median CPM, the median CPM across donors, that a cell type must reach. It applies to every expression rule; group means are not tested. Empty uses {DEFAULT_EXPRESSION_THRESHOLD:g} CPM.",
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
                                                        f"CELLEX specificity is the lowest CELLEX score that a cell type must reach in addition to Minimum CPM. It is used only with the CELLEX expression rule, and a missing score leaves the rule unresolved. Empty uses {DEFAULT_SPECIFICITY_THRESHOLD:g}.",
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
                                                        "view",
                                                        "The view chooses what the comparison counts. Distinct targets counts each target that meets the expression rule once. Canonical drugs counts each canonical drug once when any of its targets meets the rule. A cell group counts a target or canonical drug once if any of its cell types qualifies.",
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
                                                    "Measure (heatmap)",
                                                    info_tip(
                                                        "measure",
                                                        "measure",
                                                        "The measure chooses the value that the heatmap shows. Count is the number of qualifying targets or canonical drugs. Percent divides that number by all targets, or by all canonical drugs with targets, for the disease under the applied drug filters; the denominator is the same for every cell type. In the dot plot, dot area uses a compressed Count scale and color shows Percent. These values are not cell proportions or probabilities of treatment benefit.",
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
                                        "Settings applied",
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
                                    "Comparison by cell type",
                                    info_tip(
                                        "comparison",
                                        "comparison by cell type",
                                        "The comparison counts qualifying targets or canonical drugs for each disease and cell type. A cell group row counts distinct targets or canonical drugs across its cell types, not mean expression. The heatmap shows the selected measure. In the dot plot, dot area uses a compressed Count scale and color shows Percent; an empty cell means zero or missing. A ≥ in the hover marks a lower bound. Scales are separate for targets and canonical drugs and can change when filters change. Click a cell group name to expand it, or hover for counts and evidence status.",
                                    ),
                                ]
                            ),
                            html.Fieldset(
                                [
                                    html.Legend("Chart"),
                                    dcc.RadioItems(
                                        id="chart-type",
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
                                id="matrix-note", className="matrix-note", role="status"
                            ),
                        ],
                        className="section-heading",
                    ),
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
                                    html.Div(
                                        id="target-dot-key",
                                        className="dot-key",
                                        style={"display": "none"},
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
                                    html.Div(
                                        id="drug-dot-key",
                                        className="dot-key",
                                        style={"display": "none"},
                                    ),
                                ],
                                id="drug-heatmap-panel",
                                hidden=True,
                            ),
                        ],
                        className="heatmap-stack",
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
                                    "Disease details",
                                    info_tip(
                                        "selection",
                                        "disease details",
                                        "Disease details show the drug–target records and target expression for one disease. They update immediately within the last applied settings. The records show the underlying evidence, and the expression chart includes every target of those canonical drugs, even when it does not meet the expression rule. Choose a disease below.",
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
        id="drug-page",
    )
