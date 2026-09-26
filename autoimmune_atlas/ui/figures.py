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

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas.models import (
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
from autoimmune_atlas.ui.components import ordered_cell_ids
from autoimmune_atlas.ui.config import (
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
    x: Sequence[str] | None
    y: Sequence[str] | None
    z: Sequence[Sequence[int | float | None]] | None
    hovertext: Sequence[Sequence[str]] | None

    def update(self, **kwargs: object) -> object: ...


class _ColorbarData(Protocol):
    len: int | None


class _MarkerData(Protocol):
    cmax: float | None
    cmin: float | None
    color: Sequence[float] | None
    colorbar: _ColorbarData
    size: Sequence[float] | None


class _ScatterData(Protocol):
    customdata: Sequence[object] | None
    marker: _MarkerData
    name: str | None
    visible: bool | None
    x: Sequence[str | None] | None
    y: Sequence[str | None] | None
    hovertext: Sequence[str] | None

    def update(self, **kwargs: object) -> object: ...


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
    value = display_value(row, measure, kind) or "missing"
    label = "Canonical drugs" if kind == "drug" else "Targets"
    if row["status"] == "unavailable":
        assessment = "Disease data not loaded"
    elif measure == "percent" and row.get(denominator_key, 0) == 0:
        assessment = "No eligible items in the percentage denominator"
    elif _is_lower_bound(row, measure, kind):
        assessment = "≥ is a lower bound; unresolved evidence may increase this value"
    elif row.get(value_key) is None:
        assessment = "The expression rule is unresolved because target or expression data are missing."
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
            f"Canonical drugs with targets / all canonical drugs: {row.get('mapped_drugs', 0)} / {row.get('total_drugs', 0)}",
        )
    )


def dot_hover_text(row: SummaryRow, kind: Kind) -> str:
    """Dot plot の面積と色を、欠測と下限を含めて一緒に説明する。"""
    _, unknown_key, denominator_key = measure_fields(kind, "percent")
    count = display_value(row, "count", kind) or "missing"
    percent = display_value(row, "percent", kind) or "missing"
    label = "Canonical drugs" if kind == "drug" else "Targets"
    if row["status"] == "unavailable":
        assessment = "Disease data not loaded"
    elif row.get(denominator_key, 0) == 0:
        assessment = "No eligible items in the percentage denominator"
    elif _is_lower_bound(row, "count", kind) or _is_lower_bound(row, "percent", kind):
        assessment = "≥ is a lower bound; unresolved evidence may increase this value"
    elif (
        row.get(measure_fields(kind, "count")[0]) is None
        or row.get(measure_fields(kind, "percent")[0]) is None
    ):
        assessment = "The expression rule is unresolved because target or expression data are missing."
    elif row[measure_fields(kind, "count")[0]] == 0:
        assessment = f"No qualifying {label.lower()} under this rule"
    else:
        assessment = f"{label} meeting the expression rule"
    return "<br>".join(
        (
            f"<b>{row['disease']} / {row['cell']}</b>",
            f"Count: {count}",
            f"Percent: {percent}",
            assessment,
            f"Denominator: {row.get(denominator_key, 0)}",
            f"Unresolved in denominator: {row.get(unknown_key, 0)}",
            f"Canonical drugs with targets / all canonical drugs: {row.get('mapped_drugs', 0)} / {row.get('total_drugs', 0)}",
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
    title = "Canonical drugs" if kind == "drug" else "Targets"
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
                "title": f"{title} {'(%)' if measure == 'percent' else '(count)'}",
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
        font={"family": "Arial, sans-serif", "size": 12, "color": "#2d2d2d"},
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


def dot_size_scale(
    rows: Sequence[SummaryRow], kind: Kind
) -> tuple[int | float, float, list[tuple[int | float, float]]]:
    """Dot plot と面積キーで共有する件数の尺度を返す。"""
    count_key = measure_fields(kind, "count")[0]
    counts: list[int | float] = []
    for row in rows:
        count = cast(int | float | None, row.get(count_key))
        if count is not None and count > 0:
            counts.append(count)
    size_max = max(counts, default=1)
    unique_counts = sorted(set(counts))
    legend_counts: list[int | float] = []
    if unique_counts:
        for index in (0, len(unique_counts) // 2, len(unique_counts) - 1):
            if unique_counts[index] not in legend_counts:
                legend_counts.append(unique_counts[index])
    return (
        size_max,
        2 * size_max**0.75 / 22**2,
        [(value, 22 * (value / size_max) ** 0.375) for value in legend_counts],
    )


def expression_dot_diameter(value: int | float) -> float:
    """CELLEX 値を、0 も見える最大径 22 px の円へ変換する。"""
    return max(3, 22 * math.sqrt(value))


def expression_dot_size_scale() -> tuple[float, tuple[tuple[float, float], ...]]:
    """発現 Dot plot と凡例で共有する固定尺度を返す。"""
    values = (0.0, 0.25, 0.5, 1.0)
    return 1, tuple((value, expression_dot_diameter(value)) for value in values)


def build_dot_figure(
    rows: list[SummaryRow],
    disease_ids: Sequence[str],
    cell_ids: Sequence[str],
    kind: Kind = "target",
    *,
    scale_rows: list[SummaryRow] | None = None,
) -> go.Figure:
    """面積を件数の 0.75 乗、色を割合に固定した dot plot を作る。"""
    count_key = measure_fields(kind, "count")[0]
    percent_key = measure_fields(kind, "percent")[0]
    lookup = {(row["disease_id"], row["cell_id"]): row for row in rows}
    disease_names = {row["disease_id"]: row["disease"] for row in rows}
    cell_names = {row["cell_id"]: row["cell"] for row in rows}
    disease_labels = [disease_names.get(item, item) for item in disease_ids]
    cell_labels = [cell_names.get(item, item) for item in cell_ids]
    x: list[str] = []
    y: list[str] = []
    hover: list[str] = []
    custom: list[list[str]] = []
    sizes: list[int | float] = []
    colors: list[int | float] = []
    for cell_id, cell_label in zip(cell_ids, cell_labels, strict=True):
        for disease_id, disease_label in zip(disease_ids, disease_labels, strict=True):
            row = lookup.get((disease_id, cell_id))
            count = cast(
                int | float | None, None if row is None else row.get(count_key)
            )
            percent = cast(
                int | float | None, None if row is None else row.get(percent_key)
            )
            if row is None or count is None or count <= 0 or percent is None:
                continue
            x.append(disease_label)
            y.append(cell_label)
            custom.append([disease_id, cell_id])
            hover.append(dot_hover_text(row, kind))
            sizes.append(math.pow(count, 0.75))
            colors.append(percent)

    scale_rows = rows if scale_rows is None else scale_rows
    scale_percents: list[int | float] = []
    for row in scale_rows:
        percent = cast(int | float | None, row.get(percent_key))
        if percent is not None:
            scale_percents.append(percent)
    _, size_ref, _ = dot_size_scale(scale_rows, kind)
    color_min = min(scale_percents, default=0)
    color_max = max(scale_percents, default=1)
    if color_min == color_max:
        color_min, color_max = 0, max(1, color_max)
    title = "Canonical drugs" if kind == "drug" else "Targets"
    has_values = bool(x)
    value_x = x if has_values else [None, None]
    value_y = y if has_values else [None, None]
    value_sizes = sizes if has_values else [0, 0]
    value_colors = colors if has_values else [color_min, color_max]
    value_hover = hover if has_values else ["", ""]
    value_custom = custom if has_values else [None, None]
    figure = go.Figure()
    figure_ops = cast(_FigureOps, cast(object, figure))
    figure_ops.add_trace(
        go.Scatter(
            x=value_x,
            y=value_y,
            mode="markers",
            name="Values",
            showlegend=False,
            marker={
                "size": value_sizes,
                "color": value_colors,
                "sizemode": "area",
                "sizeref": size_ref,
                "colorscale": "Greens",
                "cmin": color_min,
                "cmax": color_max,
                "showscale": True,
                "colorbar": {
                    "title": f"{title} (%)",
                    "lenmode": "pixels",
                    "len": 240,
                    "y": 1,
                    "yanchor": "top",
                },
                "line": {"color": "#35543a", "width": 0.6},
            },
            hovertext=value_hover,
            hovertemplate="%{hovertext}<extra></extra>",
            customdata=value_custom,
        )
    )
    if not disease_ids or not cell_ids:
        figure_ops.add_annotation(
            text="Select diseases and cell types", showarrow=False
        )
    figure_ops.update_layout(
        template="plotly_white",
        font={"family": "Arial, sans-serif", "size": 12, "color": "#2d2d2d"},
        height=max(420, min(1400, 180 + 28 * len(cell_ids))),
        margin={"l": 180, "r": 40, "t": 130, "b": 70},
        showlegend=False,
        hoverlabel={"align": "left"},
    )
    figure_ops.update_xaxes(
        type="category",
        categoryorder="array",
        categoryarray=disease_labels,
        range=[-0.5, max(0.5, len(disease_labels) - 0.5)],
        autorange=False,
        tickangle=-45,
        side="top",
        title="Disease",
        automargin=True,
        tickmode="array",
        tickvals=disease_labels,
        ticktext=["<br>".join(disease_label_lines(name)) for name in disease_labels],
    )
    figure_ops.update_yaxes(
        type="category",
        categoryorder="array",
        categoryarray=cell_labels,
        range=[max(1, len(cell_labels)) - 0.5, -0.5],
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
    missing_z: list[list[int | None]] = []
    missing_hover: list[list[str]] = []
    positive_x: list[str] = []
    positive_y: list[str] = []
    positive_hover: list[str] = []
    dot_x: list[str] = []
    dot_y: list[str] = []
    dot_sizes: list[float] = []
    dot_colors: list[float] = []
    dot_hover: list[str] = []
    dot_custom: list[list[str]] = []
    display_cells = atlas.cell_catalog(snapshot, "mixed") if grouped else cells
    member_names = {cell["id"]: cell["name"] for cell in display_cells}
    for cell in display_cells:
        member = cell["id"]
        is_group = cell.get("cell_level") == "group"
        z_row: list[float | None] = []
        hover_row: list[str] = []
        missing_row: list[int | None] = []
        missing_hover_row: list[str] = []
        for target_id in targets:
            label = f"{target_names[target_id]} ({target_id})"
            item = metadata.get((target_id, member))
            cpm = item.get("median") if item else None
            observed: list[float] = []
            specificity_score = item.get("specificity_score") if item else None
            specificity_observed: list[float] = []
            if is_group:
                for child in cell["members"]:
                    child_item = metadata.get((target_id, child))
                    median_value = child_item.get("median") if child_item else None
                    if median_value is not None:
                        observed.append(median_value)
                    specificity_value = (
                        child_item.get("specificity_score") if child_item else None
                    )
                    if specificity_value is not None:
                        specificity_observed.append(specificity_value)
                cpm = mean(observed) if observed else None
                specificity_score = (
                    mean(specificity_observed) if specificity_observed else None
                )
            center, spread = target_stats[target_id]
            score = (
                (math.log2(1 + cpm) - center) / spread
                if cpm is not None and spread
                else (0 if cpm is not None else None)
            )
            z_row.append(score)
            state = (
                None
                if is_group
                else atlas.expression_state(item, threshold, method, specificity)
            )
            rule_state = (
                "display only for groups"
                if is_group
                else "met"
                if state is True
                else "not met"
                if state is False
                else "unresolved"
            )
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
                        f"CELLEX specificity: {specificity_score:g}"
                        if specificity_score is not None
                        else "CELLEX specificity: missing",
                        f"Target-relative median: {item.get('target_median'):g}"
                        if item and item.get("target_median") is not None
                        else "Target-relative median: missing",
                        f"Expression rule: {rule_state}",
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
                        f"Mean CELLEX specificity: {specificity_score:g}"
                        if specificity_score is not None
                        else "Mean CELLEX specificity: missing",
                        f"CELLEX-observed source cell types: {len(specificity_observed)} / {len(cell['members'])}",
                        "Equal-weight mean across source cell types with values; missing values are excluded.",
                        "Display only: excluded from reference medians, standardization, clustering and rule evaluation.",
                    )
                )
            missing_row.append(1 if cpm is None else None)
            missing_hover_row.append(
                hover_row[-1]
                + "<br>Shown in gray; expression is missing, not measured as zero."
                if cpm is None
                else ""
            )
            if cpm is not None and not is_group and state is True:
                positive_x.append(label)
                positive_y.append(member_names[member])
                positive_hover.append(hover_row[-1])
            if (
                cpm is not None
                and specificity_score is not None
                and specificity_score >= 0
            ):
                dot_x.append(label)
                dot_y.append(member_names[member])
                dot_sizes.append(expression_dot_diameter(specificity_score))
                dot_colors.append(score if score is not None else 0)
                dot_hover.append(hover_row[-1])
                dot_custom.append([target_id, member])
        z.append(z_row)
        hover.append(hover_row)
        missing_z.append(missing_row)
        missing_hover.append(missing_hover_row)
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
            hoverongaps=False,
            xgap=2,
            ygap=2,
        )
    )
    figure_ops = cast(_FigureOps, cast(object, figure))
    if any(value is not None for row in missing_z for value in row):
        figure_ops.add_trace(
            go.Heatmap(
                x=target_labels,
                y=[cell["name"] for cell in display_cells],
                z=missing_z,
                colorscale=[[0, "#d9d9d9"], [1, "#d9d9d9"]],
                showscale=False,
                name="Missing expression",
                hovertext=missing_hover,
                hovertemplate="%{hovertext}<extra></extra>",
                hoverongaps=False,
                xgap=2,
                ygap=2,
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
                    "line": {"color": "#2d2d2d", "width": 1},
                },
                name="Meets expression rule",
                showlegend=False,
                hovertext=positive_hover,
                hovertemplate="%{hovertext}<extra></extra>",
            )
        )
    dot_size_ref, _ = expression_dot_size_scale()
    has_dots = bool(dot_x)
    figure_ops.add_trace(
        go.Scatter(
            x=dot_x if has_dots else [None, None],
            y=dot_y if has_dots else [None, None],
            mode="markers",
            marker={
                "size": dot_sizes if has_dots else [0, 0],
                "color": dot_colors if has_dots else [-color_limit, color_limit],
                "sizemode": "diameter",
                "sizeref": dot_size_ref,
                "colorscale": "RdBu_r",
                "cmin": -color_limit,
                "cmax": color_limit,
                "showscale": True,
                "colorbar": {
                    "title": "Target-wise z-score",
                    "lenmode": "pixels",
                    "len": 240,
                    "y": 1,
                    "yanchor": "top",
                },
                "line": {"color": "#4a4a4a", "width": 0.6},
            },
            name="CELLEX specificity",
            showlegend=False,
            hovertext=dot_hover if has_dots else ["", ""],
            hovertemplate="%{hovertext}<extra></extra>",
            customdata=dot_custom if has_dots else [None, None],
            visible=False,
        )
    )
    if not targets:
        figure_ops.add_annotation(
            text="No targets in the selected scope", showarrow=False
        )
    figure_ops.update_layout(
        template="plotly_white",
        height=max(360, 170 + 25 * len(display_cells)),
        margin={"l": 280, "r": 40, "t": 110, "b": 60},
        font={"family": "Arial, sans-serif", "size": 12, "color": "#2d2d2d"},
        hoverlabel={"align": "left"},
    )
    figure_ops.update_xaxes(
        tickangle=-45,
        side="top",
        title="Target",
        tickvals=target_labels,
        ticktext=[target_names[t] for t in targets],
        automargin=True,
    )
    figure_ops.update_yaxes(autorange="reversed", title="Cell type", automargin=True)
    if members:
        figure_ops.update_yaxes(range=[len(display_cells) - 0.5, -0.5], autorange=False)
    return figure


def expression_view(
    base: go.Figure,
    catalog: Sequence[CellCatalogEntry],
    expanded: Sequence[str],
    chart_type: str = "heatmap",
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
        if isinstance(trace_value, go.Heatmap):
            # 欠損セルを灰色で塗る heatmap は、主 heatmap と同じ行だけ残す。
            overlay = cast(_HeatmapData, cast(object, trace_value))
            overlay.y = [heatmap_y[i] for i in indices]
            if overlay.z is not None:
                overlay.z = [overlay.z[i] for i in indices]
            if overlay.hovertext is not None:
                overlay.hovertext = [overlay.hovertext[i] for i in indices]
            continue
        trace = cast(_ScatterData, trace_value)
        trace_y = trace.y if trace.y is not None else ()
        keep = [i for i, name in enumerate(trace_y) if name is None or name in names]
        if trace.x is not None:
            trace.x = [trace.x[i] for i in keep]
        if trace.y is not None:
            trace.y = [trace.y[i] for i in keep]
        if trace.hovertext is not None:
            trace.hovertext = [trace.hovertext[i] for i in keep]
        if trace.customdata is not None:
            trace.customdata = [trace.customdata[i] for i in keep]
        if trace.name == "CELLEX specificity":
            marker_sizes = trace.marker.size or ()
            marker_colors = trace.marker.color or ()
            trace.marker.size = [marker_sizes[i] for i in keep]
            trace.marker.color = [marker_colors[i] for i in keep]
            if not trace.x:
                trace.x = [None, None]
                trace.y = [None, None]
                trace.hovertext = ["", ""]
                trace.customdata = [None, None]
                trace.marker.size = [0, 0]
                trace.marker.color = [trace.marker.cmin or 0, trace.marker.cmax or 1]
    count = max(1, len(indices))
    target_count = max(1, len(heatmap.x or ()))
    figure_ops.update_layout(
        height=110 + 60 + 28 * count, margin=dict(l=280, r=40, t=110, b=60)
    )
    target_labels = list(heatmap.x or ())
    cell_labels = list(heatmap.y or ())
    if chart_type == "dot":
        dots = next(
            trace
            for trace in figure_ops.data[1:]
            if cast(_ScatterData, trace).name == "CELLEX specificity"
        )
        dot_trace = cast(_ScatterData, dots)
        dot_trace.update(visible=True)
        dot_trace.marker.colorbar.len = min(240, 28 * count)
        figure_ops.data = (dots,)
    else:
        figure_ops.data = tuple(
            trace
            for trace in figure_ops.data
            if cast(_ScatterData, trace).name != "CELLEX specificity"
        )
    figure_ops.update_xaxes(
        type="category",
        categoryorder="array",
        categoryarray=target_labels,
        range=[-0.5, target_count - 0.5],
        autorange=False,
        automargin=False,
    )
    figure_ops.update_yaxes(
        type="category",
        categoryorder="array",
        categoryarray=cell_labels,
        range=[count - 0.5, -0.5],
        autorange=False,
        automargin=False,
        fixedrange=True,
        showticklabels=False,
        title=None,
    )
    heatmap.update(colorbar_len=min(240, 28 * count))
    return figure
