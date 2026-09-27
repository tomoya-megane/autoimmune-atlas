"""genetic association の保存データを読み、遺伝子ページ用に集計する。"""

import json
import math
from pathlib import Path
from typing import cast

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas.models import (
    CellCatalogEntry,
    GeneAssociation,
    GeneticsSnapshot,
    Snapshot,
    SummaryRow,
)

GENETICS_PATH = Path(__file__).resolve().parent.parent / "data" / "genetics.json"
GENETICS_SCHEMA = 1


def load_genetics(path: Path = GENETICS_PATH) -> GeneticsSnapshot | None:
    """schema 1 の genetics.json を読む。無ければ None、形が違えば ValueError。"""
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as handle:
        raw = cast(object, json.load(handle))
    if not isinstance(raw, dict):
        raise ValueError("genetics.json must contain a JSON object")
    data = cast(dict[str, object], raw)
    if data.get("schema") != GENETICS_SCHEMA:
        raise ValueError(
            "Unsupported genetics format. Refresh it with: pixi run refresh-genetics"
        )
    required = {
        "data_version",
        "retrieved_at",
        "source",
        "score_floor",
        "datasources",
        "associations",
        "expression",
    }
    missing = required - data.keys()
    if missing:
        raise ValueError(
            "Missing required fields in genetics.json: " + ", ".join(sorted(missing))
        )
    if not isinstance(data["associations"], dict) or not isinstance(
        data["expression"], dict
    ):
        raise ValueError("associations and expression must be objects")
    return cast(GeneticsSnapshot, cast(object, data))


def version_matches(snapshot: Snapshot, genetics: GeneticsSnapshot) -> bool:
    """2 つのファイルが同じ Open Targets の版から取得されたかを返す。"""
    return snapshot.get("data_version") == genetics["data_version"]


def merged_snapshot(snapshot: Snapshot, genetics: GeneticsSnapshot) -> Snapshot:
    """発現データを結合した Snapshot を返す。同じ遺伝子は snapshot.json を優先する。"""
    expression = {**genetics["expression"], **snapshot["expression"]}
    return cast(Snapshot, cast(object, {**snapshot, "expression": expression}))


def genes_for_disease(
    genetics: GeneticsSnapshot, disease_id: str, score_threshold: float
) -> list[GeneAssociation]:
    """閾値以上（等号を含む）の関連遺伝子をスコアの降順で返す。"""
    genes = [
        gene
        for gene in genetics["associations"].get(disease_id, [])
        if gene["score"] >= score_threshold
    ]
    return sorted(genes, key=lambda gene: (-gene["score"], gene["target_id"]))


def summarize_genes(
    snapshot: Snapshot,
    genetics: GeneticsSnapshot,
    score_threshold: object,
    threshold: object,
    *,
    method: str = "fixed",
    specificity_threshold: object = 0.5,
    level: str = "group",
    cell_ids: list[str] | None = None,
    disease_ids: list[str] | None = None,
) -> list[SummaryRow]:
    """疾患・細胞ごとに、閾値以上の遺伝子のうち発現判定が陽性の数と割合を返す。

    薬剤ページの標的と同じ三値判定、大分類の集約、割合の NA の規則を使う。
    薬剤に関する列は空の値で埋め、図の共有に使う。
    引数の検証は、遺伝子が 0 件の疾患でも行われるよう、発現の判定より前に行う。
    """
    if method not in {"fixed", "relative", "specificity"}:
        raise ValueError(f"未対応の発現判定方法: {method}")
    if (
        isinstance(threshold, bool)
        or not isinstance(threshold, (int, float))
        or not math.isfinite(threshold)
        or threshold < 0
    ):
        raise ValueError("最低 CPM は 0 以上の有限の数値にしてください")
    if (
        isinstance(specificity_threshold, bool)
        or not isinstance(specificity_threshold, (int, float))
        or not math.isfinite(specificity_threshold)
        or not 0 <= specificity_threshold <= 1
    ):
        raise ValueError("特異性閾値は 0 以上 1 以下の有限の数値にしてください")
    if (
        isinstance(score_threshold, bool)
        or not isinstance(score_threshold, (int, float))
        or not math.isfinite(score_threshold)
        or not 0 <= score_threshold <= 1
    ):
        raise ValueError("スコアの閾値は 0 以上 1 以下の有限の数値にしてください")
    metadata = atlas.expression_metadata(snapshot)
    catalog: list[CellCatalogEntry] = (
        [
            {
                "id": "all",
                "name": "All source cell types",
                "members": [
                    cell["id"] for cell in atlas.cell_catalog(snapshot, "cell")
                ],
            }
        ]
        if level == "all"
        else atlas.cell_catalog(snapshot, level)
    )
    if cell_ids is not None:
        selected_cells = set(cell_ids)
        catalog = [cell for cell in catalog if cell["id"] in selected_cells]
    selected_diseases = None if disease_ids is None else set(disease_ids)
    output: list[SummaryRow] = []
    for disease in snapshot.get("diseases", []):
        if selected_diseases is not None and disease["id"] not in selected_diseases:
            continue
        genes = genes_for_disease(genetics, disease["id"], score_threshold)
        targets = [gene["target_id"] for gene in genes]
        loaded = disease["id"] in genetics["associations"]
        for cell in catalog:
            states = {
                target: atlas.aggregate_states(
                    [
                        atlas.expression_state(
                            metadata.get((target, member)),
                            threshold,
                            method,
                            specificity_threshold,
                        )
                        for member in cell["members"]
                    ]
                )
                for target in targets
            }
            positive = {t for t, s in states.items() if s is True}
            unknown = {t for t, s in states.items() if s is None}
            output.append(
                {
                    "disease_id": disease["id"],
                    "disease": disease["name"],
                    "cell_id": cell["id"],
                    "cell": cell["name"],
                    "ontology_id": cell.get("ontology_id", cell["id"]),
                    "cell_level": cell.get("cell_level", level),
                    "count": len(positive) if loaded else None,
                    "percent": atlas.percent_of(positive, unknown, len(targets))
                    if loaded
                    else None,
                    "denominator": len(targets),
                    "unknown": len(unknown),
                    "unmapped_drugs": 0,
                    "status": "unavailable"
                    if not loaded
                    else "partial"
                    if unknown
                    else "complete",
                    "records": [],
                    "drug_count": None,
                    "drug_percent": None,
                    "drug_denominator": 0,
                    "unknown_drugs": 0,
                    "total_drugs": 0,
                    "mapped_drugs": 0,
                    "member_cell_ids": cell["members"],
                }
            )
    return output
