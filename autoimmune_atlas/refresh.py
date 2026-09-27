"""Open Targets の公開 API からアプリ用のスナップショットを作る。"""

from __future__ import annotations

import json
import math
import tempfile
import time
from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from typing import NotRequired, TypedDict, cast

import httpx
import msgspec

from autoimmune_atlas.aggregation import (
    DRUG_TYPE_TO_MODALITY,
    STAGE_FILTERS,
    cell_catalog,
)
from autoimmune_atlas.disease_catalog import SCOPE_ROOTS
from autoimmune_atlas.models import (
    CellDefinition,
    DataVersion,
    Disease,
    DrugRecord,
    ExpressionRow,
    GeneAssociation,
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


class _Score(TypedDict):
    id: str
    score: float


class _AssociationTarget(TypedDict):
    id: str
    approvedSymbol: str


class _AssociationRow(TypedDict):
    target: _AssociationTarget
    datatypeScores: list[_Score]
    datasourceScores: list[_Score]


class _AssociationBlock(TypedDict):
    count: int
    rows: list[_AssociationRow]


class _AssociationDisease(TypedDict):
    associatedTargets: _AssociationBlock


class _AssociationResponse(TypedDict):
    disease: _AssociationDisease | None


class _EvidenceRow(TypedDict):
    datasourceId: str
    datatypeId: str


class _EvidenceBlock(TypedDict):
    rows: list[_EvidenceRow]


class _EvidenceDisease(TypedDict):
    evidences: _EvidenceBlock


class _EvidenceResponse(TypedDict):
    disease: _EvidenceDisease | None


class _StoredRow(TypedDict):
    """再利用のために読む発現の行。schema 2 と genetics.json は cell などを持ち、schema 3 は持たない。"""

    cell_id: str
    median: float | None
    specificity_score: NotRequired[float | None]
    cell: NotRequired[str]
    parent_id: NotRequired[str | None]
    parent: NotRequired[str | None]
    ancestor_ids: NotRequired[list[str] | None]


class _StoredFile(TypedDict):
    data_version: NotRequired[DataVersion]
    cells: NotRequired[dict[str, CellDefinition]]
    expression: NotRequired[dict[str, list[_StoredRow]]]


_DiseaseQuery = Callable[[str, dict[str, str]], dict[str, object]]
_Query = Callable[[str, Mapping[str, object]], dict[str, object]]

API_HOST = "api.platform.opentargets.org"
API_PATH = "/api/v4/graphql"
API_URL = f"https://{API_HOST}{API_PATH}"
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
LEGACY_GENETICS_PATH = DATA_PATH.parent / "genetics.json"
PAGE_SIZE = 500
GENETIC_DATATYPE = "genetic_association"

ASSOCIATION_QUERY = """query($id:String!,$index:Int!,$size:Int!){
  disease(efoId:$id){ associatedTargets(
    page:{index:$index,size:$size}, enableIndirect:false, orderByScore:"genetic_association"
  ){ count rows{ target{id approvedSymbol} datatypeScores{id score} datasourceScores{id score} } } }
}"""

EVIDENCE_QUERY = """query($id:String!,$gene:String!,$datasource:String!){
  disease(efoId:$id){ evidences(ensemblIds:[$gene], datasourceIds:[$datasource], size:1){
    rows{ datasourceId datatypeId } } }
}"""

VERSION_QUERY = "query{meta{dataVersion{year month iteration}}}"


def query_api(
    query: str,
    variables: Mapping[str, object] | None = None,
    *,
    transport: httpx.BaseTransport | None = None,
) -> dict[str, object]:
    """GraphQL を照会し、不完全な応答を失敗として扱う。3 回まで再試行する。

    proxy は httpx が環境変数（HTTPS_PROXY、https_proxy）から読む。
    transport はテストが偽の応答を差し込むためにある。
    """
    for attempt in range(3):
        try:
            with httpx.Client(timeout=45, transport=transport) as client:
                response = client.post(
                    API_URL, json={"query": query, "variables": variables or {}}
                )
            if response.status_code != 200:
                raise RuntimeError(f"Open Targets: HTTP {response.status_code}")
            raw = cast(object, response.json())
            if not isinstance(raw, dict):
                raise RuntimeError("Open Targets: response must be an object")
            payload = cast(dict[str, object], raw)
            if payload.get("errors") or "data" not in payload:
                raise RuntimeError(f"Open Targets: {payload.get('errors', payload)}")
            data = payload["data"]
            if not isinstance(data, dict):
                raise RuntimeError("Open Targets: data must be an object")
            return cast(dict[str, object], data)
        except (httpx.HTTPError, ValueError, RuntimeError):
            if attempt == 2:
                raise
            time.sleep(attempt + 1)
    raise RuntimeError("Open Targets の照会が失敗しました")


def resolve_disease_ids(
    query: _DiseaseQuery = query_api,
) -> tuple[list[SnapshotRoot], list[str]]:
    """起点ごとに下位語を取り、和集合を返す。起点が 1 つでも無ければ止める。"""
    roots: list[SnapshotRoot] = []
    ids: set[str] = set()
    for root_id, include_descendants in SCOPE_ROOTS:
        response = msgspec.convert(
            query(
                "query($id:String!){disease(efoId:$id){id name descendants}}",
                {"id": root_id},
            ),
            _DiseaseScopeResponse,
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


def extract_expression(
    rows: list[ExpressionApiRow],
) -> tuple[list[ExpressionRow], dict[str, CellDefinition]]:
    """Tabula Sapiens の細胞型別 pseudobulk の中央値と、細胞の定数を取り出す。"""
    result: dict[str, ExpressionRow] = {}
    cells: dict[str, CellDefinition] = {}
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
            "median": value,
            "specificity_score": specificity,
        }
        cells[cell_id] = {
            "name": cell["biosampleName"],
            "parent_id": parent["biosampleId"] if parent else None,
            "parent": parent["biosampleName"] if parent else None,
            "ancestor_ids": list(cell.get("ancestors") or []),
        }
    return list(result.values()), cells


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


def fetch_expression(
    target_id: str,
) -> tuple[str, list[ExpressionRow], dict[str, CellDefinition]]:
    """発現データを全ページ取得してから、細胞型別の値と細胞の定数を選ぶ。"""
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
        response = msgspec.convert(
            query_api(query, {"id": target_id, "page": page}), _ExpressionResponse
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
    kept, cells = extract_expression(rows)
    return target_id, kept, cells


def merge_cells(
    cells: dict[str, CellDefinition], found: Mapping[str, CellDefinition]
) -> None:
    """細胞の定数を表へ足し、既にある定義と食い違えば失敗にする。"""
    for cell_id, definition in found.items():
        if cell_id in cells and cells[cell_id] != definition:
            raise ValueError(f"細胞型の名前か親分類が標的間で一貫しません: {cell_id}")
        cells[cell_id] = definition


def collect_associations(
    query: _Query, disease_ids: list[str]
) -> dict[str, list[GeneAssociation]]:
    """疾患ごとに genetic スコアの降順で取り、スコアを持たない遺伝子に当たった時点で止める。"""
    output: dict[str, list[GeneAssociation]] = {}
    for disease_id in disease_ids:
        genes: list[GeneAssociation] = []
        index = 0
        while True:
            response = msgspec.convert(
                query(
                    ASSOCIATION_QUERY,
                    {"id": disease_id, "index": index, "size": PAGE_SIZE},
                ),
                _AssociationResponse,
            )
            disease = response["disease"]
            if disease is None:
                raise ValueError(f"疾患が見つかりません: {disease_id}")
            block = disease["associatedTargets"]
            stop = False
            for row in block["rows"]:
                genetic = next(
                    (
                        s["score"]
                        for s in row["datatypeScores"]
                        if s["id"] == GENETIC_DATATYPE
                    ),
                    None,
                )
                if genetic is None:
                    stop = True
                    break
                genes.append(
                    {
                        "target_id": row["target"]["id"],
                        "target": row["target"]["approvedSymbol"],
                        "score": genetic,
                        "datasource_scores": {
                            s["id"]: s["score"] for s in row["datasourceScores"]
                        },
                    }
                )
            index += 1
            if stop or index * PAGE_SIZE >= block["count"]:
                break
            if not block["rows"]:
                # HTTP と GraphQL が成功しても、ページが途中で切れることがある。保存させない。
                raise ValueError(f"関連遺伝子の途中のページが空です: {disease_id}")
        # API は降順で返すが、保存の形として並びを保証する
        output[disease_id] = sorted(genes, key=lambda g: (-g["score"], g["target_id"]))
    return output


def classify_datasources(
    query: _Query, associations: Mapping[str, list[GeneAssociation]]
) -> list[str]:
    """現れた datasource ごとに根拠を 1 件取り、genetic_association のものだけ残す。"""
    sample: dict[str, tuple[str, str]] = {}
    for disease_id, genes in associations.items():
        for gene in genes:
            for datasource in gene["datasource_scores"]:
                sample.setdefault(datasource, (disease_id, gene["target_id"]))
    genetic: list[str] = []
    for datasource, (disease_id, gene_id) in sorted(sample.items()):
        response = msgspec.convert(
            query(
                EVIDENCE_QUERY,
                {"id": disease_id, "gene": gene_id, "datasource": datasource},
            ),
            _EvidenceResponse,
        )
        disease = response["disease"]
        rows = disease["evidences"]["rows"] if disease else []
        if not rows:
            raise ValueError(f"datasource の根拠が取れません: {datasource}")
        if rows[0]["datatypeId"] == GENETIC_DATATYPE:
            genetic.append(datasource)
    return genetic


def restrict_datasources(
    associations: Mapping[str, list[GeneAssociation]], datasources: list[str]
) -> dict[str, list[GeneAssociation]]:
    """datasource ごとのスコアを genetic association のものに絞る。"""
    keep = set(datasources)
    return {
        disease_id: [
            {
                **gene,
                "datasource_scores": {
                    key: value
                    for key, value in gene["datasource_scores"].items()
                    if key in keep
                },
            }
            for gene in genes
        ]
        for disease_id, genes in associations.items()
    }


def reusable_expression(
    paths: Sequence[Path], version: DataVersion
) -> tuple[dict[str, list[ExpressionRow]], dict[str, CellDefinition]]:
    """同じ版の保存データから、発現の行と細胞の定数を取り出す。

    schema 2 の snapshot.json、schema 1 の genetics.json、schema 3 の snapshot.json を読む。
    発現の値は先に読んだファイルの遺伝子を優先する。
    細胞の定数は、後のファイルにある遺伝子の行でも全部確かめ、食い違えば失敗にする。
    無い、読めない、版が違うファイルは飛ばす。
    load_snapshot() は schema 3 しか受け付けないので、ここでは緩い型で読む。
    """
    expression: dict[str, list[ExpressionRow]] = {}
    cells: dict[str, CellDefinition] = {}
    for path in paths:
        try:
            data = msgspec.json.decode(path.read_bytes(), type=_StoredFile)
        except (OSError, msgspec.DecodeError, msgspec.ValidationError):
            continue
        if data.get("data_version") != version:
            continue
        merge_cells(cells, data.get("cells") or {})
        for target_id, rows in (data.get("expression") or {}).items():
            kept: list[ExpressionRow] = []
            for row in rows:
                if "cell" in row:
                    # schema 2 と genetics.json の行。細胞の定数を表へ移す。
                    merge_cells(
                        cells,
                        {
                            row["cell_id"]: {
                                "name": row["cell"],
                                "parent_id": row.get("parent_id"),
                                "parent": row.get("parent"),
                                "ancestor_ids": list(row.get("ancestor_ids") or []),
                            }
                        },
                    )
                kept.append(
                    {
                        "cell_id": row["cell_id"],
                        "median": row["median"],
                        "specificity_score": row.get("specificity_score"),
                    }
                )
            if target_id not in expression:
                expression[target_id] = kept
    return expression, cells


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
    """疾患、薬剤の標的、関連遺伝子、発現量を順番に取得し、schema 3 で保存する。"""
    version = msgspec.convert(query_api(VERSION_QUERY), _VersionResponse)["meta"][
        "dataVersion"
    ]
    roots, ids = resolve_disease_ids(query_api)
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
        result = msgspec.convert(
            query_api(clinical_query, {"ids": chunk}), _ClinicalResponse
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
        result = msgspec.convert(
            query_api(
                """query($ids:[String!]!){drugs(chemblIds:$ids){id mechanismsOfAction {rows{
          mechanismOfAction actionType targets{id approvedSymbol} references{source ids urls}
        }}}}""",
                {"ids": chunk},
            ),
            _DrugsResponse,
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
    disease_ids = [disease["id"] for disease in diseases]
    associations = collect_associations(query_api, disease_ids)
    print(
        f"疾患: {len(disease_ids)} / 関連遺伝子を持つ疾患: {sum(bool(g) for g in associations.values())}",
        flush=True,
    )
    datasources = classify_datasources(query_api, associations)
    print(f"genetic association の datasource: {', '.join(datasources)}", flush=True)
    associations = restrict_datasources(associations, datasources)
    drug_targets = {
        target_id for record in records.values() if (target_id := record["target_id"])
    }
    gene_targets = {g["target_id"] for genes in associations.values() for g in genes}
    target_ids = sorted(drug_targets | gene_targets)
    reused, cells = reusable_expression([DATA_PATH, LEGACY_GENETICS_PATH], version)
    expression: dict[str, list[ExpressionRow]] = {
        target_id: reused[target_id] for target_id in target_ids if target_id in reused
    }
    needed = [target_id for target_id in target_ids if target_id not in expression]
    print(
        f"細胞型別発現: 再利用 {len(expression)} / 取得 {len(needed)} / 全体 {len(target_ids)}",
        flush=True,
    )
    with ThreadPoolExecutor(max_workers=2) as pool:
        for index, (target, rows, found) in enumerate(
            pool.map(fetch_expression, needed), 1
        ):
            expression[target] = rows
            merge_cells(cells, found)
            print(f"細胞型別発現: {index}/{len(needed)}", flush=True)
    unknown_cells = {
        row["cell_id"] for rows in expression.values() for row in rows
    } - set(cells)
    if unknown_cells:
        raise ValueError(
            f"cells に無い細胞 ID が発現にあります: {sorted(unknown_cells)[:5]}"
        )
    final_version = msgspec.convert(query_api(VERSION_QUERY), _VersionResponse)["meta"][
        "dataVersion"
    ]
    if final_version != version:
        raise ValueError("取得中にデータの版が変わりました。再取得してください")
    snapshot: Snapshot = {
        "schema": 3,
        "root": ROOT_ID,
        "roots": roots,
        "data_version": version,
        "retrieved_at": datetime.now(UTC).isoformat(),
        "source": API_URL,
        "diseases": sorted(diseases, key=lambda d: d["name"]),
        "records": list(records.values()),
        "associations": associations,
        "datasources": datasources,
        "cells": cells,
        "expression": expression,
    }
    cell_catalog(snapshot)
    save_snapshot(DATA_PATH, snapshot)
    print(
        f"保存: {DATA_PATH}\n疾患 {len(diseases)} / 薬剤 {len(drug_ids)} / 標的 {len(drug_targets)} / 関連遺伝子 {len(gene_targets)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
