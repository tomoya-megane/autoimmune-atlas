"""薬剤、標的、細胞の対応から表示用の集計を作る。"""

from __future__ import annotations

import csv
import io
from collections import defaultdict


def summarize(snapshot: dict, annotations: list[dict], mode: str, modality: str, threshold: float) -> list[dict]:
    """同じ標的を一度だけ数え、未判定を下限値または欠測として残す。

    Parameters
    ----------
    snapshot : dict
        取得した疾患、薬剤と発現量。
    annotations : list[dict]
        疾患と薬剤と標的を指定した、出典付きの薬効細胞注釈。
    mode : str
        expression または mechanism。
    modality : str
        all、small_molecule、antibody、other。
    threshold : float
        発現陽性とする中央値の下限。この値を超えると陽性。

    Returns
    -------
    list[dict]
        疾患と細胞の組合せごとの標的数、割合、欠測数と該当薬剤。

    """
    cells = {}
    expression = {}
    for target, rows in snapshot.get("expression", {}).items():
        for row in rows:
            cells[row["cell_id"]] = row["cell"]
            expression[target, row["cell_id"]] = row["median"]
    evidence = {}
    for row in annotations:
        cells[row["cell_id"]] = row["cell"]
        if "drug_id" in row:
            evidence[row["disease_id"], row["drug_id"], row["target_id"], row["cell_id"]] = row
    records = defaultdict(list)
    for row in snapshot.get("records", []):
        if modality == "all" or row["modality"] == modality:
            records[row["disease_id"]].append(row)
    output = []
    for disease in snapshot["diseases"]:
        drugs = records[disease["id"]]
        targets = {r["target_id"] for r in drugs if r["target_id"]}
        unmapped = {r["drug_id"] for r in drugs if not r["target_id"]}
        for cell_id, cell in sorted(cells.items(), key=lambda item: item[1]):
            positive, negative = set(), set()
            selected = []
            for row in drugs:
                target = row["target_id"]
                if not target:
                    continue
                value = expression.get((target, cell_id))
                annotation = evidence.get((disease["id"], row["drug_id"], target, cell_id))
                if mode == "expression":
                    state = None if value is None else value > threshold
                    source = "https://platform.opentargets.org/target/" + target
                    note = f"Tabula Sapiens / median {value:g} CPM" if value is not None else ""
                else:
                    state = None if annotation is None else annotation["status"] == "yes"
                    source = annotation["source"] if annotation else ""
                    note = annotation["note"] if annotation else ""
                if state is True:
                    positive.add(target)
                    selected.append({**row, "cell_id": cell_id, "cell": cell, "evidence": source, "note": note})
                elif state is False:
                    negative.add((row["drug_id"], target))
            # 薬効モードの陰性は、その標的を持つ全薬剤で確認できたときに限る。
            ruled_out = {
                target for target in targets
                if all((r["drug_id"], target) in negative for r in drugs if r["target_id"] == target)
            }
            unknown = len(targets - positive - ruled_out)
            complete = disease["status"] == "ready"
            incomplete = unknown > 0 or bool(unmapped)
            count = len(positive) if complete and (positive or not incomplete) else None
            # 分母は既知の標的だけ。標的未判明の薬剤数は別に返す。
            percent = 100 * count / len(targets) if count is not None and targets else None
            status = "unavailable" if not complete else "partial" if incomplete else "complete"
            output.append({
                "disease_id": disease["id"], "disease": disease["name"],
                "cell_id": cell_id, "cell": cell, "count": count, "percent": percent,
                "denominator": len(targets), "unknown": unknown,
                "unmapped_drugs": len(unmapped), "status": status, "records": selected,
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
        writer.writerow({
            key: "'" + value if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")) else value
            for key, value in row.items()
        })
    return output.getvalue()
