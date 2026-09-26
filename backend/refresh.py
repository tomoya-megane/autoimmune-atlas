"""Open Targets の公開 API からアプリ用のスナップショットを作る。"""

from __future__ import annotations

import base64
import http.client
import json
import math
import os
import tempfile
import time
from collections.abc import Callable, Mapping
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from typing import NotRequired, TypedDict, cast
from urllib.parse import urlparse

from backend.aggregation import (
    DRUG_TYPE_TO_MODALITY,
    STAGE_FILTERS,
    cell_catalog,
)
from backend.disease_catalog import SCOPE_ROOTS
from backend.models import (
    DataVersion,
    Disease,
    DrugRecord,
    ExpressionRow,
    Reference,
    Snapshot,
    SnapshotRoot,
)


class _DiseaseScope(TypedDict):
    id: str
    name: str
    descendants: list[str]


class _DiseaseScopeResponse(TypedDict):
    disease: _DiseaseScope | None


class _Drug(TypedDict):
    id: str
    name: str
    drugType: str
    parentMolecule: _ParentDrug | None


class _ParentDrug(TypedDict):
    id: str
    name: str


class ClinicalCandidate(TypedDict):
    maxClinicalStage: str
    drug: _Drug | None


class _NormalizedDrug(TypedDict):
    drug_id: str
    drug: str
    canonical_drug_id: str
    canonical_drug: str
    modality: str
    drug_type: str
    stage: str


class _Biosample(TypedDict):
    biosampleId: str
    biosampleName: str
    ancestors: NotRequired[list[str] | None]


class ExpressionApiRow(TypedDict):
    datasourceId: str
    unit: str
    median: object
    specificity_score: NotRequired[object]
    tissueBiosample: object | None
    celltypeBiosample: _Biosample | None
    celltypeBiosampleParent: _Biosample | None


class _Target(TypedDict):
    id: str
    approvedSymbol: str


class _ApiReference(TypedDict):
    source: str | None
    ids: list[str] | None
    urls: list[str] | None


class MechanismApiRow(TypedDict):
    mechanismOfAction: str
    actionType: str | None
    targets: list[_Target]
    references: list[_ApiReference] | None


class _MechanismAccumulator(TypedDict):
    target_id: str
    target: str
    mechanisms: set[str]
    action_types: set[str]
    references: dict[tuple[str, tuple[str, ...], tuple[str, ...]], Reference]


class _MechanismRecord(TypedDict):
    target_id: str
    target: str
    mechanism: str
    action_types: list[str]
    references: list[Reference]


class _ClinicalDrug(_NormalizedDrug):
    disease_id: str
    disease: str


class _DiseaseParent(TypedDict):
    id: str


class _DiseaseClinical(TypedDict):
    id: str
    name: str
    parents: list[_DiseaseParent]
    drugAndClinicalCandidates: _ClinicalBlock


class _ClinicalBlock(TypedDict):
    count: int
    rows: list[ClinicalCandidate]


class _ClinicalResponse(TypedDict):
    diseases: list[_DiseaseClinical]


class _DrugMechanisms(TypedDict):
    rows: list[MechanismApiRow]


class _DrugResponseRow(TypedDict):
    id: str
    mechanismsOfAction: _DrugMechanisms | None


class _DrugsResponse(TypedDict):
    drugs: list[_DrugResponseRow]


class _ExpressionBlock(TypedDict):
    count: int
    rows: list[ExpressionApiRow]


class _ExpressionTarget(TypedDict):
    baselineExpression: _ExpressionBlock


class _ExpressionResponse(TypedDict):
    target: _ExpressionTarget | None


class _Meta(TypedDict):
    dataVersion: DataVersion


class _VersionResponse(TypedDict):
    meta: _Meta


_DiseaseQuery = Callable[[str, dict[str, str]], dict[str, object]]

API_HOST = "api.platform.opentargets.org"
API_PATH = "/api/v4/graphql"
ROOT_ID = SCOPE_ROOTS[0][0]
KNOWN_STAGES = STAGE_FILTERS["phase1"] | {
    "UNKNOWN",
    "WITHDRAWAL",
    "EARLY_PHASE_1",
    "IND",
    "PRECLINICAL",
    "PHASE_0",
}
DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "snapshot.json"


def query_api(
    query: str, variables: Mapping[str, object] | None = None
) -> dict[str, object]:
    """proxy を含む既存の通信環境で照会し、不完全な応答を失敗として扱う。"""
    proxy = urlparse(
        os.environ.get("https_proxy") or os.environ.get("HTTPS_PROXY") or ""
    )
    for attempt in range(3):
        connection = http.client.HTTPSConnection(
            proxy.hostname or API_HOST,
            proxy.port if proxy.hostname else 443,
            timeout=45,
        )
        if proxy.hostname:
            headers: dict[str, str] = {}
            if proxy.username:
                token = base64.b64encode(
                    f"{proxy.username}:{proxy.password or ''}".encode()
                ).decode()
                headers["Proxy-Authorization"] = "Basic " + token
            connection.set_tunnel(API_HOST, 443, headers)
        try:
            connection.request(
                "POST",
                API_PATH,
                json.dumps({"query": query, "variables": variables or {}}),
                {"Content-Type": "application/json", "Accept": "application/json"},
            )
            response = connection.getresponse()
            body = response.read()
            if response.status != 200:
                raise RuntimeError(f"Open Targets: HTTP {response.status}")
            raw = cast(object, json.loads(body))
            if not isinstance(raw, dict):
                raise RuntimeError("Open Targets: response must be an object")
            payload = cast(dict[str, object], raw)
            if payload.get("errors") or "data" not in payload:
                raise RuntimeError(f"Open Targets: {payload.get('errors', payload)}")
            data = payload["data"]
            if not isinstance(data, dict):
                raise RuntimeError("Open Targets: data must be an object")
            return cast(dict[str, object], data)
        except (OSError, http.client.HTTPException, ValueError, RuntimeError):
            if attempt == 2:
                raise
            time.sleep(attempt + 1)
        finally:
            connection.close()
    raise RuntimeError("Open Targets の照会が失敗しました")


def resolve_disease_ids(
    query: _DiseaseQuery = query_api,
) -> tuple[list[SnapshotRoot], list[str]]:
    """起点ごとに下位語を取り、和集合を返す。起点が 1 つでも無ければ止める。"""
    roots: list[SnapshotRoot] = []
    ids: set[str] = set()
    for root_id, include_descendants in SCOPE_ROOTS:
        response = cast(
            _DiseaseScopeResponse,
            cast(
                object,
                query(
                    "query($id:String!){disease(efoId:$id){id name descendants}}",
                    {"id": root_id},
                ),
            ),
        )
        disease = response["disease"]
        if disease is None:
            raise ValueError(f"起点の疾患が見つかりません: {root_id}")
        members = {root_id}
        if include_descendants:
            members |= set(disease["descendants"])
        roots.append(
            {
                "id": root_id,
                "name": disease["name"],
                "include_descendants": include_descendants,
                "count": len(members),
            }
        )
        ids |= members
    return roots, sorted(ids)


def normalize_drugs(rows: list[ClinicalCandidate]) -> list[_NormalizedDrug]:
    """薬剤全体の段階を使わず、指定疾患での Phase I 以降を選ぶ。"""
    drugs: list[_NormalizedDrug] = []
    seen: set[str] = set()
    for row in rows:
        stage = row["maxClinicalStage"]
        if stage not in KNOWN_STAGES:
            raise ValueError(f"未対応の臨床段階: {stage}")
        if stage not in STAGE_FILTERS["phase1"]:
            continue
        drug = row.get("drug")
        if not drug:
            raise ValueError("臨床段階を満たす行に薬剤情報がありません")
        if drug["id"] in seen:
            raise ValueError(f"疾患内の薬剤行が重複しています: {drug['id']}")
        seen.add(drug["id"])
        kind = drug["drugType"]
        try:
            modality = DRUG_TYPE_TO_MODALITY[kind]
        except KeyError as error:
            raise ValueError(f"未対応の薬剤型: {kind}") from error
        parent = drug.get("parentMolecule") or drug
        drugs.append(
            {
                "drug_id": drug["id"],
                "drug": drug["name"].lower(),
                "canonical_drug_id": parent["id"],
                "canonical_drug": parent["name"].lower(),
                "modality": modality,
                "drug_type": kind,
                "stage": stage,
            }
        )
    return drugs


def extract_expression(rows: list[ExpressionApiRow]) -> list[ExpressionRow]:
    """Tabula Sapiens の細胞型別 pseudobulk の中央値を取り出す。"""
    result: dict[str, ExpressionRow] = {}
    for row in rows:
        cell = row.get("celltypeBiosample")
        if (
            row["datasourceId"] != "tabula_sapiens"
            or not cell
            or row.get("tissueBiosample")
        ):
            continue
        if row["unit"] != "CPM(pseudobulk sum[counts])":
            raise ValueError(f"発現単位が変わりました: {row['unit']}")
        value = row["median"]
        if value is not None and (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0
        ):
            raise ValueError("発現中央値が不正です")
        specificity = row.get("specificity_score")
        if specificity is not None and (
            isinstance(specificity, bool)
            or not isinstance(specificity, (int, float))
            or not math.isfinite(specificity)
            or not 0 <= specificity <= 1
        ):
            raise ValueError("特異性スコアが不正です")
        cell_id = cell["biosampleId"]
        if cell_id in result:
            raise ValueError(f"細胞型別の発現行が重複しています: {cell_id}")
        parent = row.get("celltypeBiosampleParent")
        result[cell_id] = {
            "cell_id": cell_id,
            "cell": cell["biosampleName"],
            "median": value,
            "specificity_score": specificity,
            "parent_id": parent["biosampleId"] if parent else None,
            "parent": parent["biosampleName"] if parent else None,
            "ancestor_ids": list(cell.get("ancestors") or []),
        }
    return list(result.values())


def extract_mechanisms(rows: list[MechanismApiRow]) -> list[_MechanismRecord]:
    """標的ごとに作用機序と最低限の出典を重複なくまとめる。"""
    targets: dict[str, _MechanismAccumulator] = {}
    for mechanism in rows:
        for target in mechanism.get("targets", []):
            item = targets.setdefault(
                target["id"],
                {
                    "target_id": target["id"],
                    "target": target["approvedSymbol"],
                    "mechanisms": set(),
                    "action_types": set(),
                    "references": {},
                },
            )
            item["mechanisms"].add(mechanism["mechanismOfAction"])
            action_type = mechanism.get("actionType")
            if action_type:
                item["action_types"].add(action_type)
            for reference in mechanism.get("references") or []:
                normalized: Reference = {
                    "source": reference.get("source") or "",
                    "ids": list(reference.get("ids") or []),
                    "urls": list(reference.get("urls") or []),
                }
                key = (
                    normalized["source"],
                    tuple(normalized["ids"]),
                    tuple(normalized["urls"]),
                )
                item["references"][key] = normalized
    return [
        {
            "target_id": item["target_id"],
            "target": item["target"],
            "mechanism": "; ".join(sorted(item["mechanisms"])),
            "action_types": sorted(item["action_types"]),
            "references": list(item["references"].values()),
        }
        for item in sorted(targets.values(), key=lambda item: item["target_id"])
    ]


def fetch_expression(target_id: str) -> tuple[str, list[ExpressionRow]]:
    """発現データを全ページ取得してから、細胞型別の値を選ぶ。"""
    query = """query($id:String!, $page:Int!) {
      target(ensemblId:$id) { baselineExpression(page:{index:$page,size:3000}) {
        count rows { datasourceId unit median
        specificity_score tissueBiosample {biosampleId biosampleName}
          celltypeBiosampleParent {biosampleId biosampleName}
          celltypeBiosample {biosampleId biosampleName ancestors} }
      } }
    }"""
    rows: list[ExpressionApiRow] = []
    page = 0
    while True:
        response = cast(
            _ExpressionResponse,
            cast(object, query_api(query, {"id": target_id, "page": page})),
        )
        target = response["target"]
        if target is None:
            raise ValueError(f"標的が見つかりません: {target_id}")
        block = target["baselineExpression"]
        rows.extend(block["rows"])
        if len(rows) >= block["count"]:
            if len(rows) != block["count"]:
                raise ValueError("発現行の件数が一致しません")
            break
        if not block["rows"]:
            raise ValueError("発現データの途中のページが空です")
        page += 1
    return target_id, extract_expression(rows)


def save_snapshot(path: Path, snapshot: object) -> None:
    """全取得が成功してからファイルを置き換え、失敗時は前回の内容を保つ。"""
    content = json.dumps(
        snapshot, ensure_ascii=False, allow_nan=False, separators=(",", ":")
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent, delete=False
        ) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        temporary.replace(path)
    finally:
        if temporary and temporary.exists():
            temporary.unlink()


def main() -> None:
    """疾患、薬剤の標的、発現量を順番に取得する。"""
    version_query = "query{meta{dataVersion{year month iteration}}}"
    version = cast(_VersionResponse, cast(object, query_api(version_query)))["meta"][
        "dataVersion"
    ]
    roots, ids = resolve_disease_ids()
    for root in roots:
        print(f"起点: {root['name']} {root['count']} 疾患", flush=True)
    print(f"対象: {len(roots)} 起点の和集合 {len(ids)} 疾患", flush=True)
    diseases: list[Disease] = []
    clinical: list[_ClinicalDrug] = []
    clinical_query = """query($ids:[String!]!){diseases(efoIds:$ids){id name parents{id}
      drugAndClinicalCandidates {count rows {maxClinicalStage drug{id name drugType parentMolecule{id name}}}}
    }}"""
    for start in range(0, len(ids), 5):
        chunk = ids[start : start + 5]
        result = cast(
            _ClinicalResponse,
            cast(object, query_api(clinical_query, {"ids": chunk})),
        )["diseases"]
        if {d["id"] for d in result} != set(chunk):
            raise ValueError("下位疾患の取得結果に不足があります")
        for disease in result:
            block = disease["drugAndClinicalCandidates"]
            rows = block.get("rows", [])
            if len(rows) != block["count"]:
                raise ValueError(f"薬剤の取得に不足があります: {disease['id']}")
            diseases.append(
                {
                    "id": disease["id"],
                    "name": disease["name"],
                    "parent_ids": [parent["id"] for parent in disease["parents"]],
                    "status": "ready",
                    "unclassified_stages": sum(
                        r["maxClinicalStage"] in {"UNKNOWN", "WITHDRAWAL"} for r in rows
                    ),
                }
            )
            for drug in normalize_drugs(rows):
                clinical_drug: _ClinicalDrug = {
                    "disease_id": disease["id"],
                    "disease": disease["name"],
                    **drug,
                }
                clinical.append(clinical_drug)
        print(f"疾患: {min(start + 5, len(ids))}/{len(ids)}", flush=True)
    drug_ids = sorted({d["drug_id"] for d in clinical})
    mechanisms: dict[str, list[_MechanismRecord]] = {}
    for start in range(0, len(drug_ids), 20):
        chunk = drug_ids[start : start + 20]
        result = cast(
            _DrugsResponse,
            cast(
                object,
                query_api(
                    """query($ids:[String!]!){drugs(chemblIds:$ids){id mechanismsOfAction {rows{
          mechanismOfAction actionType targets{id approvedSymbol} references{source ids urls}
        }}}}""",
                    {"ids": chunk},
                ),
            ),
        )["drugs"]
        if {d["id"] for d in result} != set(chunk):
            raise ValueError("薬剤の標的情報に不足があります")
        for drug in result:
            mechanisms[drug["id"]] = extract_mechanisms(
                (drug.get("mechanismsOfAction") or {}).get("rows", [])
            )
        print(
            f"薬剤の標的: {min(start + 20, len(drug_ids))}/{len(drug_ids)}", flush=True
        )
    records: dict[tuple[str, str, str], DrugRecord] = {}
    for row in clinical:
        for target in mechanisms[row["drug_id"]] or [
            {
                "target_id": "",
                "target": "Unknown",
                "mechanism": "",
                "action_types": [],
                "references": [],
            }
        ]:
            key = row["disease_id"], row["drug_id"], target["target_id"]
            record: DrugRecord = {**row, **target}
            records[key] = record
    target_ids = sorted(
        {target_id for record in records.values() if (target_id := record["target_id"])}
    )
    expression: dict[str, list[ExpressionRow]] = {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        for index, (target, rows) in enumerate(
            pool.map(fetch_expression, target_ids), 1
        ):
            expression[target] = rows
            print(f"細胞型別発現: {index}/{len(target_ids)}", flush=True)
    final_version = cast(_VersionResponse, cast(object, query_api(version_query)))[
        "meta"
    ]["dataVersion"]
    if final_version != version:
        raise ValueError("取得中にデータの版が変わりました。再取得してください")
    snapshot: Snapshot = {
        "schema": 2,
        "root": ROOT_ID,
        "roots": roots,
        "data_version": version,
        "retrieved_at": datetime.now(UTC).isoformat(),
        "source": f"https://{API_HOST}{API_PATH}",
        "diseases": sorted(diseases, key=lambda d: d["name"]),
        "records": list(records.values()),
        "expression": expression,
    }
    cell_catalog(snapshot)
    save_snapshot(DATA_PATH, snapshot)
    print(
        f"保存: {DATA_PATH}\n疾患 {len(diseases)} / 薬剤 {len(drug_ids)} / 標的 {len(target_ids)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
