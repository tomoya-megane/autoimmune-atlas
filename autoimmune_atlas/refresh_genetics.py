"""Open Targets の genetic association から遺伝子ページ用のデータを作る。"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from typing import TypedDict, cast

from autoimmune_atlas.genetics import GENETICS_PATH, GENETICS_SCHEMA, SCORE_FLOOR
from autoimmune_atlas.models import (
    DataVersion,
    ExpressionRow,
    GeneAssociation,
    GeneticsSnapshot,
    Snapshot,
)
from autoimmune_atlas.refresh import (
    API_HOST,
    API_PATH,
    fetch_expression,
    query_api,
    save_snapshot,
)
from autoimmune_atlas.snapshot import load_snapshot

PAGE_SIZE = 500
GENETIC_DATATYPE = "genetic_association"

Query = Callable[[str, Mapping[str, object]], dict[str, object]]


class _Score(TypedDict):
    id: str
    score: float


class _Target(TypedDict):
    id: str
    approvedSymbol: str


class _AssociationRow(TypedDict):
    target: _Target
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


class _Meta(TypedDict):
    dataVersion: DataVersion


class _VersionResponse(TypedDict):
    meta: _Meta


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


def collect_associations(
    query: Query, disease_ids: list[str], floor: float = SCORE_FLOOR
) -> dict[str, list[GeneAssociation]]:
    """疾患ごとに genetic スコアの降順で取り、floor を下回った時点で止める。"""
    output: dict[str, list[GeneAssociation]] = {}
    for disease_id in disease_ids:
        genes: list[GeneAssociation] = []
        index = 0
        while True:
            response = cast(
                _AssociationResponse,
                cast(
                    object,
                    query(
                        ASSOCIATION_QUERY,
                        {"id": disease_id, "index": index, "size": PAGE_SIZE},
                    ),
                ),
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
                if genetic is None or genetic < floor:
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
            if stop or not block["rows"] or index * PAGE_SIZE >= block["count"]:
                break
        output[disease_id] = genes
    return output


def classify_datasources(
    query: Query, associations: Mapping[str, list[GeneAssociation]]
) -> list[str]:
    """現れた datasource ごとに根拠を 1 件取り、genetic_association のものだけ残す。"""
    sample: dict[str, tuple[str, str]] = {}
    for disease_id, genes in associations.items():
        for gene in genes:
            for datasource in gene["datasource_scores"]:
                sample.setdefault(datasource, (disease_id, gene["target_id"]))
    genetic: list[str] = []
    for datasource, (disease_id, gene_id) in sorted(sample.items()):
        response = cast(
            _EvidenceResponse,
            cast(
                object,
                query(
                    EVIDENCE_QUERY,
                    {"id": disease_id, "gene": gene_id, "datasource": datasource},
                ),
            ),
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


def main() -> None:
    """snapshot.json の疾患について関連遺伝子と不足分の発現を取り、genetics.json を書く。"""
    snapshot: Snapshot | None = load_snapshot()
    if snapshot is None:
        raise SystemExit(
            "snapshot.json がありません。先に pixi run refresh を実行してください"
        )
    version = cast(_VersionResponse, cast(object, query_api(VERSION_QUERY)))["meta"][
        "dataVersion"
    ]
    if snapshot.get("data_version") != version:
        raise SystemExit(
            "snapshot.json の版が API と違います。先に pixi run refresh を実行してください"
        )
    disease_ids = [disease["id"] for disease in snapshot["diseases"]]
    associations = collect_associations(query_api, disease_ids)
    print(
        f"疾患: {len(disease_ids)} / 関連遺伝子を持つ疾患: {sum(bool(g) for g in associations.values())}",
        flush=True,
    )
    datasources = classify_datasources(query_api, associations)
    print(f"genetic association の datasource: {', '.join(datasources)}", flush=True)
    associations = restrict_datasources(associations, datasources)
    needed = sorted(
        {gene["target_id"] for genes in associations.values() for gene in genes}
        - set(snapshot["expression"])
    )
    expression: dict[str, list[ExpressionRow]] = {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        for index, (target, rows) in enumerate(pool.map(fetch_expression, needed), 1):
            expression[target] = rows
            print(f"細胞型別発現: {index}/{len(needed)}", flush=True)
    final_version = cast(_VersionResponse, cast(object, query_api(VERSION_QUERY)))[
        "meta"
    ]["dataVersion"]
    if final_version != version:
        raise ValueError("取得中にデータの版が変わりました。再取得してください")
    genetics: GeneticsSnapshot = {
        "schema": GENETICS_SCHEMA,
        "data_version": version,
        "retrieved_at": datetime.now(UTC).isoformat(),
        "source": f"https://{API_HOST}{API_PATH}",
        "score_floor": SCORE_FLOOR,
        "datasources": datasources,
        "associations": associations,
        "expression": expression,
    }
    save_snapshot(GENETICS_PATH, genetics)
    total_genes = len({g["target_id"] for gs in associations.values() for g in gs})
    print(
        f"保存: {GENETICS_PATH}\n疾患 {len(disease_ids)} / 遺伝子 {total_genes}",
        flush=True,
    )


if __name__ == "__main__":
    main()
