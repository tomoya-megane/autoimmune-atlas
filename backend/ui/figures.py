"""比較と発現の図を組み立てる。"""

import json
import math
from collections.abc import Mapping, Sequence
from statistics import mean, pstdev
from textwrap import wrap
from typing import Literal, Protocol, cast

import plotly.graph_objects as go  # pyright: ignore[reportMissingTypeStubs] -- Plotly は型スタブを配布していない。
from scipy.cluster.hierarchy import (  # pyright: ignore[reportMissingTypeStubs] -- SciPy は型スタブを配布していない。
    leaves_list,  # pyright: ignore[reportUnknownVariableType] -- SciPy は型を公開していない。
    linkage,  # pyright: ignore[reportUnknownVariableType] -- SciPy は型を公開していない。
)

from backend import aggregation as atlas
from backend.models import (
    CellCatalogEntry,
    DrugRecord,
    EvidenceBase,
    EvidenceTableRow,
    ExpressionMetadata,
    ExpressionStateInput,
    FilteredRecord,
    Snapshot,
    SummaryRow,
    TargetRecord,
)
from backend.ui.components import ordered_cell_ids
from backend.ui.config import (
    DEFAULT_EXPRESSION_THRESHOLD,
    DEFAULT_SPECIFICITY_THRESHOLD,
)

type MetricValueKey = Literal["drug_percent", "drug_count", "percent", "count"]
type UnknownKey = Literal["unknown_drugs", "unknown"]
type DenominatorKey = Literal["drug_denominator", "denominator"]
type Measure = Literal["percent", "count"]
type Kind = Literal["target", "drug"]
type ExpressionInput = ExpressionMetadata | ExpressionStateInput
type ExpressionRecord = DrugRecord | FilteredRecord | EvidenceTableRow | TargetRecord


class _Linkage(Protocol):
    def __call__(
        self,
        observations: Sequence[Sequence[float]],
        *,
        method: str,
        metric: str,
        optimal_ordering: bool,
    ) -> object: ...


class _LeavesList(Protocol):
    def __call__(self, tree: object) -> Sequence[int]: ...


class _FigureOps(Protocol):
    data: tuple[object, ...]

    def add_annotation(self, **kwargs: object) -> object: ...

    def add_trace(self, trace: object) -> object: ...

    def update_layout(self, **kwargs: object) -> object: ...

    def update_xaxes(self, **kwargs: object) -> object: ...

    def update_yaxes(self, **kwargs: object) -> object: ...


class _HeatmapData(Protocol):
    y: Sequence[str] | None
    z: Sequence[Sequence[int | float | None]] | None
    hovertext: Sequence[Sequence[str]] | None

    def update(self, **kwargs: object) -> object: ...


class _ScatterData(Protocol):
    x: Sequence[str] | None
    y: Sequence[str] | None
    hovertext: Sequence[str] | None


def heatmap_cell_ids(
    snapshot: Snapshot | None,
    selected: Sequence[str] | None,
    expanded: Sequence[str] | None,
    *,
    catalog: Sequence[CellCatalogEntry] | None = None,
) -> list[str]:
    """選択済みの行に、展開中の大分類の元細胞を追加する。"""
    visible = set(selected or [])
    if catalog is None:
        if snapshot is None:
            raise ValueError("snapshot is required when catalog is omitted")
        catalog = atlas.cell_catalog(snapshot, "mixed")
    for cell in catalog:
        if (
            cell.get("cell_level") == "group"
            and cell["id"] in visible
            and cell["id"] in (expanded or [])
        ):
            visible.update(cell["members"])
    return [cell["id"] for cell in catalog if cell["id"] in visible]


def measure_fields(
    kind: str, measure: str
) -> tuple[
    MetricValueKey,
    UnknownKey,
    DenominatorKey,
]:
    if kind == "drug":
        return (
            "drug_percent" if measure == "percent" else "drug_count",
            "unknown_drugs",
            "drug_denominator",
        )
    return "percent" if measure == "percent" else "count", "unknown", "denominator"


def _is_lower_bound(row: SummaryRow, measure: Measure, kind: Kind = "target") -> bool:
    """実数と割合で異なる下限判定を、その指標自身の未判定から求める。"""
    value_key, unknown_key, _ = measure_fields(kind, measure)
    return row.get(value_key) is not None and (
        row.get(unknown_key, 0) > 0
        or (measure == "count" and row.get("unmapped_drugs", 0) > 0)
    )


def display_value(row: SummaryRow, measure: Measure, kind: Kind = "target") -> str:
    value_key = measure_fields(kind, measure)[0]
    if value_key == "drug_percent":
        value = row["drug_percent"]
    elif value_key == "drug_count":
        value = row["drug_count"]
    elif value_key == "percent":
        value = row["percent"]
    else:
        value = row["count"]
    if value is None:
        return ""
    rendered = f"{value:.1f}%" if measure == "percent" else str(value)
    return "≥" + rendered if _is_lower_bound(row, measure, kind) else rendered


def hover_text(row: SummaryRow, measure: Measure, kind: Kind) -> str:
    value_key, unknown_key, denominator_key = measure_fields(kind, measure)
    value = display_value(row, measure, kind) or "Not available"
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
            f"<b>{row['disease']} / {row['cell']}</b>",
            f"{label}: {value}",
            assessment,
            *(
                ["Displayed as zero for readability; not a measured zero."]
                if row.get(value_key) is None
                else []
            ),
            f"Denominator: {row.get(denominator_key, 0)}",
            f"Unresolved in denominator: {row.get(unknown_key, 0)}",
            f"Drugs with known targets / all drugs: {row.get('mapped_drugs', 0)} / {row.get('total_drugs', 0)}",
        )
    )


def disease_label_lines(name: str) -> list[str]:
    """長い疾患名を単語とハイフン付きの語を保って折り返す。"""
    return wrap(name, width=32, break_long_words=False, break_on_hyphens=False) or [
        name
    ]


def build_figure(
    rows: list[SummaryRow],
    disease_ids: Sequence[str],
    cell_ids: Sequence[str],
    measure: Measure,
    kind: Kind = "target",
    *,
    scale_rows: list[SummaryRow] | None = None,
) -> go.Figure:
    """0、下限値、欠測を区別した target または drug heatmap を作る。"""
    value_key = measure_fields(kind, measure)[0]
    lookup = {(row["disease_id"], row["cell_id"]): row for row in rows}
    disease_names = {row["disease_id"]: row["disease"] for row in rows}
    cell_names = {row["cell_id"]: row["cell"] for row in rows}
    z: list[list[int | float | None]] = []
    texts: list[list[str]] = []
    hovers: list[list[str]] = []
    customs: list[list[list[str]]] = []
    for cell_id in cell_ids:
        z_row: list[int | float | None] = []
        text_row: list[str] = []
        hover_row: list[str] = []
        custom_row: list[list[str]] = []
        for disease_id in disease_ids:
            row = lookup.get((disease_id, cell_id))
            value = None if row is None else row.get(value_key)
            z_row.append(value)
            text_row.append(
                display_value(row, measure, kind).removeprefix("≥")
                if row is not None and value is not None
                else ("0.0%" if measure == "percent" else "0")
            )
            hover_row.append(
                "No summary for this combination"
                if row is None
                else hover_text(row, measure, kind)
            )
            custom_row.append([disease_id, cell_id])
        z.append(z_row)
        texts.append(text_row)
        hovers.append(hover_row)
        customs.append(custom_row)
    values: list[int | float] = (
        [value for row in scale_rows if (value := row.get(value_key)) is not None]
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
            # 欠測を0にするのは描画用の配列だけで、集計と色範囲の計算には加えない。
            z=[[0 if value is None else value for value in row] for row in z],
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
    figure_ops = cast(_FigureOps, cast(object, figure))
    if not disease_ids or not cell_ids:
        figure_ops.add_annotation(
            text="Select diseases and cell types", showarrow=False
        )
    figure_ops.update_layout(
        template="plotly_white",
        font={"family": "Arial, sans-serif", "size": 12, "color": "#263238"},
        height=max(420, min(1400, 180 + 28 * len(cell_ids))),
        margin={"l": 180, "r": 40, "t": 130, "b": 70},
        legend={"orientation": "h", "y": -0.08, "yanchor": "top", "x": 0},
        hoverlabel={"align": "left"},
    )
    labels = [disease_names.get(item, item) for item in disease_ids]
    figure_ops.update_xaxes(
        tickangle=-45,
        side="top",
        title="Disease",
        automargin=True,
        tickmode="array",
        tickvals=labels,
        ticktext=["<br>".join(disease_label_lines(name)) for name in labels],
    )
    # 欠測表示用の文字が追加されても、行見出しとマスの中心を揃える。
    figure_ops.update_yaxes(
        range=[max(1, len(cell_ids)) - 0.5, -0.5],
        autorange=False,
        title="Cell type",
    )
    return figure


def evidence_rows(
    snapshot: Snapshot,
    row: SummaryRow,
    *,
    modality: str,
    stage: str,
    threshold: object,
    method: str,
    specificity: object,
    metadata: Mapping[tuple[str, str], ExpressionMetadata] | None = None,
    records: Sequence[FilteredRecord] | None = None,
    cell_names: Mapping[str, str] | None = None,
) -> list[EvidenceTableRow]:
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
    members = ordered_cell_ids(snapshot, "cell", row["member_cell_ids"])
    output: list[EvidenceTableRow] = []
    for record in records:
        base: EvidenceBase = {
            **record,
            "group_cell_id": row["cell_id"],
            "group_cell": row["cell"],
            "references_json": json.dumps(
                record.get("references") or [], ensure_ascii=False, sort_keys=True
            ),
        }
        target_id = record.get("target_id")
        if not target_id:
            evidence: EvidenceTableRow = {
                **base,
                "evidence_cell_id": "",
                "evidence_cell": "",
                "cpm": None,
                "specificity_score": None,
                "target_median": None,
                "support_state": "unmapped",
                "contributing": False,
            }
            output.append(evidence)
            continue
        for member in members:
            item = metadata.get((target_id, member))
            state = atlas.expression_state(item, threshold, method, specificity)
            evidence = {
                **base,
                "evidence_cell_id": member,
                "evidence_cell": item.get("cell", cell_names.get(member, member))
                if item
                else cell_names.get(member, member),
                "cpm": item.get("median") if item else None,
                "specificity_score": item.get("specificity_score") if item else None,
                "target_median": item.get("target_median") if item else None,
                "support_state": "positive"
                if state is True
                else "negative"
                if state is False
                else "unknown",
                "contributing": state is True,
            }
            output.append(evidence)
    return output


def _clustered_targets(
    targets: Sequence[str], profiles: Mapping[str, Sequence[float]]
) -> list[str]:
    """全細胞の発現パターンが近い標的を隣接させ、欠測を末尾に置く。"""
    complete = [target for target in targets if target in profiles]
    if len(complete) > 1:
        tree = cast(_Linkage, linkage)(
            [profiles[target] for target in complete],
            method="average",
            metric="euclidean",
            optimal_ordering=True,
        )
        complete = [complete[index] for index in cast(_LeavesList, leaves_list)(tree)]
    return complete + [target for target in targets if target not in profiles]


def expression_figure(
    snapshot: Snapshot,
    records: Sequence[ExpressionRecord],
    metadata: Mapping[tuple[str, str], ExpressionInput] | None = None,
    *,
    threshold: object = DEFAULT_EXPRESSION_THRESHOLD,
    method: str = "fixed",
    specificity: object = DEFAULT_SPECIFICITY_THRESHOLD,
    grouped: bool = False,
) -> go.Figure:
    """全元細胞について、現在の薬剤条件に含まれる標的の連続発現量を示す。"""
    metadata = metadata if metadata is not None else atlas.expression_metadata(snapshot)
    target_names: dict[str, str] = {}
    for record in records:
        target_id = record.get("target_id")
        if target_id:
            target_names.setdefault(target_id, record.get("target") or target_id)
    cells = atlas.cell_catalog(snapshot, "cell")
    members = [cell["id"] for cell in cells]
    targets = sorted(
        target_names, key=lambda item: (target_names[item].casefold(), item)
    )
    # 表示用の大分類平均は、標準化とクラスタリングの基準に含めない。
    target_stats: dict[str, tuple[float, float]] = {}
    profiles: dict[str, tuple[float, ...]] = {}
    for target_id in targets:
        values: list[float] = []
        for member in members:
            item = metadata.get((target_id, member))
            median_value = item.get("median") if item else None
            if median_value is not None:
                values.append(math.log2(1 + median_value))
        target_stats[target_id] = (mean(values), pstdev(values)) if values else (0, 0)
        if members and len(values) == len(members):
            center, spread = target_stats[target_id]
            profiles[target_id] = tuple(
                (value - center) / spread if spread else 0 for value in values
            )
    targets = _clustered_targets(targets, profiles)
    z: list[list[float | None]] = []
    hover: list[list[str]] = []
    missing_x: list[str] = []
    missing_y: list[str] = []
    missing_hover: list[str] = []
    positive_x: list[str] = []
    positive_y: list[str] = []
    positive_hover: list[str] = []
    display_cells = atlas.cell_catalog(snapshot, "mixed") if grouped else cells
    member_names = {cell["id"]: cell["name"] for cell in display_cells}
    for cell in display_cells:
        member = cell["id"]
        is_group = cell.get("cell_level") == "group"
        z_row: list[float | None] = []
        hover_row: list[str] = []
        for target_id in targets:
            label = f"{target_names[target_id]} ({target_id})"
            item = metadata.get((target_id, member))
            cpm = item.get("median") if item else None
            observed: list[float] = []
            if is_group:
                for child in cell["members"]:
                    child_item = metadata.get((target_id, child))
                    median_value = child_item.get("median") if child_item else None
                    if median_value is not None:
                        observed.append(median_value)
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
                        f"Observed source cell types: {len(observed)} / {len(cell['members'])}",
                        f"Target-wise z-score: {score:.2f}"
                        if score is not None
                        else "Target-wise z-score: missing",
                        "Equal-weight mean across cell types with data; missing values are excluded.",
                        "Display only: excluded from reference medians, standardization, clustering and rule evaluation.",
                    )
                )
            if cpm is None:
                missing_x.append(label)
                missing_y.append(member_names[member])
                missing_hover.append(
                    hover_row[-1]
                    + "<br>Displayed as 0; expression is missing, not measured as zero."
                )
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
    figure_ops = cast(_FigureOps, cast(object, figure))
    if missing_x:
        figure_ops.add_trace(
            go.Scatter(
                x=missing_x,
                y=missing_y,
                mode="text",
                text="0",
                textfont={"size": 12, "color": "#263238"},
                name="Missing expression",
                hovertext=missing_hover,
                hovertemplate="%{hovertext}<extra></extra>",
            )
        )
    if positive_x:
        figure_ops.add_trace(
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
        figure_ops.add_annotation(
            text="No known targets in the selected scope", showarrow=False
        )
    figure_ops.update_layout(
        template="plotly_white",
        height=max(360, 170 + 25 * len(display_cells)),
        margin={"l": 280, "r": 40, "t": 110, "b": 60},
        font={"family": "Arial, sans-serif", "size": 11, "color": "#263238"},
        hoverlabel={"align": "left"},
    )
    figure_ops.update_xaxes(
        tickangle=-45,
        side="top",
        title="Known target",
        tickvals=target_labels,
        ticktext=[target_names[t] for t in targets],
        automargin=True,
    )
    figure_ops.update_yaxes(
        autorange="reversed", title="Source cell type", automargin=True
    )
    if members:
        figure_ops.update_yaxes(range=[len(display_cells) - 0.5, -0.5], autorange=False)
    return figure


def expression_view(
    base: go.Figure,
    catalog: Sequence[CellCatalogEntry],
    expanded: Sequence[str],
) -> go.Figure:
    """集計済みの図から表示行だけを選び、色範囲と標的順を保つ。"""
    groups = [cell["id"] for cell in catalog if cell.get("cell_level") == "group"]
    visible = set(heatmap_cell_ids(None, groups, expanded, catalog=catalog))
    names = {cell["name"] for cell in catalog if cell["id"] in visible}
    figure = go.Figure(base)
    figure_ops = cast(_FigureOps, cast(object, figure))
    heatmap = cast(_HeatmapData, figure_ops.data[0])
    heatmap_y = heatmap.y if heatmap.y is not None else ()
    indices = [i for i, name in enumerate(heatmap_y) if name in names]
    heatmap.y = [heatmap_y[i] for i in indices]
    if heatmap.z is not None:
        heatmap.z = [heatmap.z[i] for i in indices]
    if heatmap.hovertext is not None:
        heatmap.hovertext = [heatmap.hovertext[i] for i in indices]
    for trace_value in figure_ops.data[1:]:
        trace = cast(_ScatterData, trace_value)
        trace_y = trace.y if trace.y is not None else ()
        keep = [i for i, name in enumerate(trace_y) if name in names]
        if trace.x is not None:
            trace.x = [trace.x[i] for i in keep]
        if trace.y is not None:
            trace.y = [trace.y[i] for i in keep]
        if trace.hovertext is not None:
            trace.hovertext = [trace.hovertext[i] for i in keep]
    count = max(1, len(indices))
    figure_ops.update_layout(
        height=110 + 60 + 28 * count, margin=dict(l=280, r=40, t=110, b=60)
    )
    figure_ops.update_xaxes(automargin=False)
    figure_ops.update_yaxes(
        range=[count - 0.5, -0.5],
        autorange=False,
        automargin=False,
        fixedrange=True,
        showticklabels=False,
        title=None,
    )
    heatmap.update(colorbar_len=min(240, 28 * count))
    return figure
