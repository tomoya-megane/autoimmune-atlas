"""薬剤、標的、細胞の対応から表示用の集計を作る。"""

from __future__ import annotations

import csv
import io
import math
from collections import defaultdict
from statistics import median

DRUG_TYPE_MODALITIES = (
    ("Small molecule", "small_molecule"),
    ("Antibody", "antibody"),
    ("Protein", "protein"),
    ("Cell", "cell"),
    ("Gene", "gene"),
    ("Enzyme", "enzyme"),
    ("Oligonucleotide", "oligonucleotide"),
    ("Antibody drug conjugate", "antibody_drug_conjugate"),
    ("Vaccine component", "vaccine_component"),
    ("Oligosaccharide", "oligosaccharide"),
    ("Unknown", "unknown"),
)
DRUG_TYPE_TO_MODALITY = dict(DRUG_TYPE_MODALITIES)
STAGE_FILTERS = {
    "phase1": {"PHASE_1", "PHASE_1_2", "PHASE_2", "PHASE_2_3", "PHASE_3", "PREAPPROVAL", "APPROVAL", "PHASE_4"},
    "phase2": {"PHASE_2", "PHASE_2_3", "PHASE_3", "PREAPPROVAL", "APPROVAL", "PHASE_4"},
    "phase3": {"PHASE_3", "PREAPPROVAL", "APPROVAL", "PHASE_4"},
    "approved": {"APPROVAL", "PHASE_4"},
}
STAGE_ORDER = {stage: rank for rank, stages in enumerate((
    ("PHASE_1",), ("PHASE_1_2",), ("PHASE_2",), ("PHASE_2_3",),
    ("PHASE_3",), ("PREAPPROVAL",), ("APPROVAL", "PHASE_4"),
)) for stage in stages}
T_CELL_ID = "CL_0000084"
CELL_GROUP_ORDER = (
    # 免疫・造血・幹細胞：T、B、自然リンパ球、骨髄系、樹状、顆粒球、造血、赤血球、幹細胞
    "CL_0000084", "CL_0000945", "CL_0001065", "CL_0000766", "CL_0000451",
    "CL_0000094", "CL_0000988", "CL_0000764", "CL_0000034",
    # 間質・血管・上皮系：間質、線維芽、結合組織、収縮、内皮、上皮、分泌
    "CL_0000499", "CL_0000057", "CL_0002320", "CL_0000183", "CL_0000115",
    "CL_0000066", "CL_0000151",
    # 神経・色素・生殖系：神経、グリア、メラノサイト、雌性生殖、雄性生殖
    "CL_0000540", "CL_0000125", "CL_0000148", "CL_0000021", "CL_0000015",
)
CELL_GROUP_RANK = {cell_id: rank for rank, cell_id in enumerate(CELL_GROUP_ORDER)}


def _validate_snapshot(snapshot: dict) -> None:
    if snapshot.get("schema") != 2:
        raise ValueError("snapshot must use schema 2")


def _record_modality(row: dict) -> str:
    try:
        return DRUG_TYPE_TO_MODALITY[row["drug_type"]]
    except KeyError as error:
        raise ValueError(f"未対応の薬剤型: {row.get('drug_type')}") from error


def filtered_records(snapshot: dict, modality: str = "all", stage: str = "phase3") -> list[dict]:
    """疾患内の有効成分ごとの最高段階を使って元記録を絞る。"""
    _validate_snapshot(snapshot)
    if modality != "all" and modality not in DRUG_TYPE_TO_MODALITY.values():
        raise ValueError(f"未対応のモダリティ: {modality}")
    if stage not in STAGE_FILTERS:
        raise ValueError(f"未対応の臨床段階フィルター: {stage}")
    candidates = []
    maximum = {}
    for row in snapshot.get("records", []):
        row_modality = _record_modality(row)
        clinical_stage = row.get("stage")
        if clinical_stage not in STAGE_ORDER:
            raise ValueError(f"集計できない臨床段階: {clinical_stage}")
        canonical_id = row.get("canonical_drug_id")
        canonical_name = row.get("canonical_drug")
        if not canonical_id or not canonical_name:
            raise ValueError("schema 2 の有効成分 ID または名称がありません")
        key = row["disease_id"], canonical_id
        if key not in maximum or STAGE_ORDER[clinical_stage] > STAGE_ORDER[maximum[key]]:
            maximum[key] = clinical_stage
        if modality == "all" or row_modality == modality:
            candidates.append({**row, "modality": row_modality})
    return [
        {**row, "canonical_stage": maximum[row["disease_id"], row["canonical_drug_id"]]}
        for row in candidates
        if maximum[row["disease_id"], row["canonical_drug_id"]] in STAGE_FILTERS[stage]
    ]


def _cell_memberships(snapshot: dict) -> dict[str, dict]:
    cells = {}
    for rows in snapshot.get("expression", {}).values():
        for row in rows:
            cell_id = row["cell_id"]
            ancestors = tuple(row.get("ancestor_ids") or ())
            if cell_id == T_CELL_ID or T_CELL_ID in ancestors:
                group_id, group_name = T_CELL_ID, "T cell"
            else:
                group_id = row.get("parent_id") or cell_id
                group_name = row.get("parent") or row["cell"]
            definition = {
                "id": cell_id, "name": row["cell"], "group_id": group_id, "group_name": group_name,
                "ancestor_ids": tuple(sorted(set(ancestors))),
            }
            if cell_id in cells and cells[cell_id] != definition:
                raise ValueError(f"細胞型の親分類が標的間で一貫しません: {cell_id}")
            cells[cell_id] = definition
    return cells


def _lineage_order(cells: dict[str, dict], members: list[str]) -> list[str]:
    """観測された祖先だけを使い、親を子より先に一度ずつ並べる。"""
    member_set = set(members)
    key = lambda cell_id: (cells[cell_id]["name"].casefold(), cell_id)
    parents = {}
    for cell_id in members:
        candidates = member_set.intersection(cells[cell_id]["ancestor_ids"])
        closest = [
            candidate for candidate in candidates
            if not any(candidate in cells[other]["ancestor_ids"] for other in candidates if other != candidate)
        ]
        if closest:
            parents[cell_id] = min(closest, key=key)
    children = defaultdict(list)
    for child, parent in parents.items():
        children[parent].append(child)
    ordered, seen = [], set()

    def visit(cell_id):
        if cell_id in seen:
            return
        seen.add(cell_id)
        ordered.append(cell_id)
        for child in sorted(children[cell_id], key=key):
            visit(child)

    for cell_id in sorted((cell_id for cell_id in members if cell_id not in parents), key=key):
        visit(cell_id)
    for cell_id in sorted(members, key=key):
        visit(cell_id)
    return ordered


def cell_catalog(snapshot: dict, level: str = "group") -> list[dict]:
    """細胞型または大分類と、所属する元細胞 ID を返す。"""
    _validate_snapshot(snapshot)
    if level not in {"group", "cell", "mixed"}:
        raise ValueError(f"未対応の細胞分類レベル: {level}")
    cells = _cell_memberships(snapshot)
    groups = {}
    for item in cells.values():
        group = groups.setdefault(item["group_id"], {"id": item["group_id"], "name": item["group_name"], "members": []})
        if group["name"] != item["group_name"]:
            raise ValueError(f"親分類の名称が一貫しません: {item['group_id']}")
        group["members"].append(item["id"])
    ordered_groups = sorted(
        groups.values(),
        key=lambda item: (CELL_GROUP_RANK.get(item["id"], len(CELL_GROUP_RANK)), item["name"].casefold(), item["id"]),
    )
    for group in ordered_groups:
        group["members"] = _lineage_order(cells, group["members"])
    if level == "group":
        return ordered_groups
    if level == "mixed":
        mixed = []
        for group in ordered_groups:
            mixed.append({**group, "id": "group:" + group["id"], "name": group["name"] + " (group)", "ontology_id": group["id"], "cell_level": "group"})
            for cell_id in group["members"]:
                name = cells[cell_id]["name"] + (" (source)" if cell_id == group["id"] else "")
                mixed.append({"id": cell_id, "name": name, "members": [cell_id], "ontology_id": cell_id, "cell_level": "cell"})
        return mixed
    return [
        {"id": cell_id, "name": cells[cell_id]["name"], "members": [cell_id]}
        for group in ordered_groups for cell_id in group["members"]
    ]


def expression_metadata(snapshot: dict) -> dict[tuple[str, str], dict]:
    """各標的・細胞の発現値と、全参照細胞から求めた標的内中央値を返す。"""
    _validate_snapshot(snapshot)
    all_cells = set(_cell_memberships(snapshot))
    result = {}
    for target_id, rows in snapshot.get("expression", {}).items():
        by_cell = {row["cell_id"]: row for row in rows}
        values = [row.get("median") for row in rows]
        target_median = median(values) if set(by_cell) == all_cells and values and all(value is not None for value in values) else None
        for row in rows:
            result[target_id, row["cell_id"]] = {**row, "target_median": target_median}
    return result


def expression_state(metadata: dict | None, threshold: float, method: str, specificity_threshold: float = 0.75) -> bool | None:
    """一つの標的・細胞について陽性、陰性、未知を返す。"""
    if method not in {"fixed", "relative", "specificity"}:
        raise ValueError(f"未対応の発現判定方法: {method}")
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not math.isfinite(threshold) or threshold < 0:
        raise ValueError("最低 CPM は 0 以上の有限の数値にしてください")
    if isinstance(specificity_threshold, bool) or not isinstance(specificity_threshold, (int, float)) or not math.isfinite(specificity_threshold) or not 0 <= specificity_threshold <= 1:
        raise ValueError("特異性閾値は 0 以上 1 以下の有限の数値にしてください")
    return _expression_state(metadata, threshold, method, specificity_threshold)


def _expression_state(metadata: dict | None, threshold: float, method: str, specificity_threshold: float) -> bool | None:
    if metadata is None or metadata.get("median") is None:
        return None
    if metadata["median"] <= threshold:
        return False
    if method == "fixed":
        return True
    if method == "relative":
        return None if metadata["target_median"] is None else metadata["median"] >= metadata["target_median"]
    score = metadata.get("specificity_score")
    return None if score is None else score >= specificity_threshold


def _aggregate(states: list[bool | None]) -> bool | None:
    if any(state is True for state in states):
        return True
    if states and all(state is False for state in states):
        return False
    return None


def _value(positives: set, unknown: set, extra_unknown: bool = False) -> int | None:
    if positives:
        return len(positives)
    return None if unknown or extra_unknown else 0


def _percent(positives: set, unknown: set, denominator: int) -> float | None:
    if not denominator:
        return None
    if positives:
        return 100 * len(positives) / denominator
    return None if unknown else 0


def summarize(
    snapshot: dict,
    modality: str,
    threshold: float,
    *,
    stage: str = "phase3",
    method: str = "fixed",
    specificity_threshold: float = 0.75,
    level: str = "group",
    cell_ids: list[str] | None = None,
    disease_ids: list[str] | None = None,
) -> list[dict]:
    """疾患・細胞または全元細胞の標的数と有効成分数を三値判定で集計する。"""
    if method not in {"fixed", "relative", "specificity"}:
        raise ValueError(f"未対応の発現判定方法: {method}")
    if level not in {"group", "cell", "mixed", "all"}:
        raise ValueError(f"未対応の細胞分類レベル: {level}")
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not math.isfinite(threshold) or threshold < 0:
        raise ValueError("最低 CPM は 0 以上の有限の数値にしてください")
    if isinstance(specificity_threshold, bool) or not isinstance(specificity_threshold, (int, float)) or not math.isfinite(specificity_threshold) or not 0 <= specificity_threshold <= 1:
        raise ValueError("特異性閾値は 0 以上 1 以下の有限の数値にしてください")
    records_by_disease = defaultdict(list)
    for row in filtered_records(snapshot, modality, stage):
        records_by_disease[row["disease_id"]].append(row)
    metadata = expression_metadata(snapshot)
    catalog = ([{"id": "all", "name": "All source cell types", "members": [cell["id"] for cell in cell_catalog(snapshot, "cell")]}]
               if level == "all" else cell_catalog(snapshot, level))
    if cell_ids is not None:
        selected_cells = set(cell_ids)
        catalog = [cell for cell in catalog if cell["id"] in selected_cells]
    selected_diseases = None if disease_ids is None else set(disease_ids)
    output = []
    for disease in snapshot.get("diseases", []):
        if selected_diseases is not None and disease["id"] not in selected_diseases:
            continue
        records = records_by_disease[disease["id"]]
        targets = {row["target_id"] for row in records if row["target_id"]}
        drug_targets = defaultdict(set)
        drugs = set()
        for row in records:
            drug_id = row["canonical_drug_id"]
            drugs.add(drug_id)
            if row["target_id"]:
                drug_targets[drug_id].add(row["target_id"])
        mapped_drugs = set(drug_targets)
        unmapped_drugs = drugs - mapped_drugs
        for cell in catalog:
            target_states = {
                target: _aggregate([_expression_state(metadata.get((target, member)), threshold, method, specificity_threshold) for member in cell["members"]])
                for target in targets
            }
            positive_targets = {target for target, state in target_states.items() if state is True}
            unknown_targets = {target for target, state in target_states.items() if state is None}
            drug_states = {drug: _aggregate([target_states[target] for target in drug_targets[drug]]) for drug in mapped_drugs}
            positive_drugs = {drug for drug, state in drug_states.items() if state is True}
            unknown_drugs = {drug for drug, state in drug_states.items() if state is None}
            selected = []
            for row in records:
                target = row["target_id"]
                if target not in positive_targets:
                    continue
                for member in cell["members"]:
                    item = metadata.get((target, member))
                    if _expression_state(item, threshold, method, specificity_threshold) is True:
                        selected.append({
                            **row, "cell_id": member, "cell": item["cell"], "cpm": item["median"],
                            "specificity_score": item.get("specificity_score"), "target_median": item["target_median"],
                            "evidence": "https://platform.opentargets.org/target/" + target,
                            "note": f"Tabula Sapiens / median {item['median']:g} CPM",
                        })
            complete = disease.get("status") == "ready"
            incomplete = bool(unknown_targets or unknown_drugs or unmapped_drugs)
            target_count = _value(positive_targets, unknown_targets, bool(unmapped_drugs)) if complete else None
            drug_count = _value(positive_drugs, unknown_drugs, bool(unmapped_drugs)) if complete else None
            output.append({
                "disease_id": disease["id"], "disease": disease["name"], "cell_id": cell["id"], "cell": cell["name"],
                "ontology_id": cell.get("ontology_id", cell["id"]), "cell_level": cell.get("cell_level", level),
                "count": target_count, "percent": _percent(positive_targets, unknown_targets, len(targets)) if complete else None,
                "denominator": len(targets), "unknown": len(unknown_targets), "unmapped_drugs": len(unmapped_drugs),
                "status": "unavailable" if not complete else "partial" if incomplete else "complete", "records": selected,
                "drug_count": drug_count, "drug_percent": _percent(positive_drugs, unknown_drugs, len(mapped_drugs)) if complete else None,
                "drug_denominator": len(mapped_drugs), "unknown_drugs": len(unknown_drugs), "total_drugs": len(drugs),
                "mapped_drugs": len(mapped_drugs), "member_cell_ids": cell["members"],
            })
    return output


def to_csv(rows: list[dict]) -> str:
    """欠測を空欄にし、表計算ソフトで数式と解釈される文字列を保護する。"""
    if not rows:
        return ""
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=list(rows[0]))
    writer.writeheader()
    for row in rows:
        writer.writerow({key: "'" + value if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")) else value for key, value in row.items()})
    return output.getvalue()
