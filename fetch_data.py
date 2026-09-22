"""Open Targets の公開 API からアプリ用のスナップショットを作る。"""

from __future__ import annotations

import base64
import http.client
import json
import math
import os
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

API_HOST = "api.platform.opentargets.org"
API_PATH = "/api/v4/graphql"
ROOT_ID = "MONDO_0007179"
LATE_STAGES = {"PHASE_3", "PREAPPROVAL", "PHASE_4", "APPROVAL"}
KNOWN_STAGES = LATE_STAGES | {"UNKNOWN", "WITHDRAWAL", "PHASE_2_3", "PHASE_2", "PHASE_1_2", "PHASE_1", "EARLY_PHASE_1", "IND", "PRECLINICAL", "PHASE_0"}
DATA_PATH = Path(__file__).parent / "data" / "snapshot.json"


def query_api(query: str, variables: dict | None = None) -> dict:
    """proxy を含む既存の通信環境で照会し、不完全な応答を失敗として扱う。"""
    proxy = urlparse(os.environ.get("https_proxy") or os.environ.get("HTTPS_PROXY") or "")
    for attempt in range(3):
        connection = http.client.HTTPSConnection(proxy.hostname or API_HOST, proxy.port if proxy.hostname else 443, timeout=45)
        if proxy.hostname:
            headers = {}
            if proxy.username:
                token = base64.b64encode(f"{proxy.username}:{proxy.password or ''}".encode()).decode()
                headers["Proxy-Authorization"] = "Basic " + token
            connection.set_tunnel(API_HOST, 443, headers)
        try:
            connection.request("POST", API_PATH, json.dumps({"query": query, "variables": variables or {}}), {"Content-Type": "application/json", "Accept": "application/json"})
            response = connection.getresponse()
            body = response.read()
            if response.status != 200:
                raise RuntimeError(f"Open Targets: HTTP {response.status}")
            payload = json.loads(body)
            if payload.get("errors") or "data" not in payload:
                raise RuntimeError(f"Open Targets: {payload.get('errors', payload)}")
            return payload["data"]
        except (OSError, http.client.HTTPException, ValueError, RuntimeError):
            if attempt == 2:
                raise
            time.sleep(attempt + 1)
        finally:
            connection.close()
    raise RuntimeError("Open Targets の照会が失敗しました")


def normalize_drugs(rows: list[dict]) -> list[dict]:
    """薬剤全体の段階を使わず、指定疾患での Phase III 以降を選ぶ。"""
    drugs = []
    for row in rows:
        stage = row["maxClinicalStage"]
        if stage not in KNOWN_STAGES:
            raise ValueError(f"未対応の臨床段階: {stage}")
        if stage not in LATE_STAGES:
            continue
        drug = row.get("drug")
        if not drug:
            raise ValueError("臨床段階を満たす行に薬剤情報がありません")
        kind = drug["drugType"]
        modality = {"Small molecule": "small_molecule", "Antibody": "antibody", "Unknown": "unknown"}.get(kind, "other")
        drugs.append({"drug_id": drug["id"], "drug": drug["name"].lower(), "modality": modality, "drug_type": kind, "stage": stage})
    return drugs


def extract_expression(rows: list[dict]) -> list[dict]:
    """Tabula Sapiens の細胞型別 pseudobulk の中央値を取り出す。"""
    result = {}
    for row in rows:
        cell = row.get("celltypeBiosample")
        if row["datasourceId"] != "tabula_sapiens" or not cell or row.get("tissueBiosample"):
            continue
        if row["unit"] != "CPM(pseudobulk sum[counts])":
            raise ValueError(f"発現単位が変わりました: {row['unit']}")
        value = row["median"]
        if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0):
            raise ValueError("発現中央値が不正です")
        cell_id = cell["biosampleId"]
        if cell_id in result:
            raise ValueError(f"細胞型別の発現行が重複しています: {cell_id}")
        result[cell_id] = {"cell_id": cell_id, "cell": cell["biosampleName"], "median": value}
    return list(result.values())


def fetch_expression(target_id: str) -> tuple[str, list[dict]]:
    """発現データを全ページ取得してから、細胞型別の値を選ぶ。"""
    query = """query($id:String!, $page:Int!) {
      target(ensemblId:$id) { baselineExpression(page:{index:$page,size:3000}) {
        count rows { datasourceId unit median
          tissueBiosample {biosampleId biosampleName}
          celltypeBiosample {biosampleId biosampleName} }
      } }
    }"""
    rows = []
    page = 0
    while True:
        target = query_api(query, {"id": target_id, "page": page})["target"]
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


def save_snapshot(path: Path, snapshot: dict) -> None:
    """全取得が成功してからファイルを置き換え、失敗時は前回の内容を保つ。"""
    content = json.dumps(snapshot, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        temporary.replace(path)
    finally:
        if temporary and temporary.exists():
            temporary.unlink()


def main() -> None:
    """疾患、薬剤の標的、発現量を順番に取得する。"""
    version_query = 'query{meta{dataVersion{year month iteration}}}'
    version = query_api(version_query)["meta"]["dataVersion"]
    root = query_api('query($id:String!){disease(efoId:$id){id name descendants}}', {"id": ROOT_ID})["disease"]
    ids = sorted(set(root["descendants"]))
    print(f"対象: {root['name']} の下位 {len(ids)} 疾患", flush=True)
    diseases, clinical = [], []
    clinical_query = """query($ids:[String!]!){diseases(efoIds:$ids){id name
      drugAndClinicalCandidates {count rows {maxClinicalStage drug{id name drugType}}}
    }}"""
    for start in range(0, len(ids), 5):
        chunk = ids[start:start + 5]
        result = query_api(clinical_query, {"ids": chunk})["diseases"]
        if {d["id"] for d in result} != set(chunk):
            raise ValueError("下位疾患の取得結果に不足があります")
        for disease in result:
            block = disease["drugAndClinicalCandidates"]
            rows = block.get("rows", [])
            if len(rows) != block["count"]:
                raise ValueError(f"薬剤の取得に不足があります: {disease['id']}")
            diseases.append({"id": disease["id"], "name": disease["name"], "status": "ready", "unclassified_stages": sum(r["maxClinicalStage"] in {"UNKNOWN", "WITHDRAWAL"} for r in rows)})
            for drug in normalize_drugs(rows):
                clinical.append({"disease_id": disease["id"], "disease": disease["name"], **drug})
        print(f"疾患: {min(start + 5, len(ids))}/{len(ids)}", flush=True)
    drug_ids = sorted({d["drug_id"] for d in clinical})
    mechanisms = {}
    for start in range(0, len(drug_ids), 20):
        chunk = drug_ids[start:start + 20]
        result = query_api("""query($ids:[String!]!){drugs(chemblIds:$ids){id mechanismsOfAction {rows{mechanismOfAction targets{id approvedSymbol}}}}}""", {"ids": chunk})["drugs"]
        if {d["id"] for d in result} != set(chunk):
            raise ValueError("薬剤の標的情報に不足があります")
        for drug in result:
            targets = {}
            for mechanism in (drug.get("mechanismsOfAction") or {}).get("rows", []):
                for target in mechanism.get("targets", []):
                    item = targets.setdefault(target["id"], {"target_id": target["id"], "target": target["approvedSymbol"], "mechanisms": set()})
                    item["mechanisms"].add(mechanism["mechanismOfAction"])
            mechanisms[drug["id"]] = [{"target_id": t["target_id"], "target": t["target"], "mechanism": "; ".join(sorted(t["mechanisms"]))} for t in targets.values()]
        print(f"薬剤の標的: {min(start + 20, len(drug_ids))}/{len(drug_ids)}", flush=True)
    records = {}
    for row in clinical:
        for target in mechanisms[row["drug_id"]] or [{"target_id": "", "target": "未判明", "mechanism": ""}]:
            key = row["disease_id"], row["drug_id"], target["target_id"]
            records[key] = {**row, **target}
    target_ids = sorted({r["target_id"] for r in records.values()} - {""})
    expression = {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        for index, (target, rows) in enumerate(pool.map(fetch_expression, target_ids), 1):
            expression[target] = rows
            print(f"細胞型別発現: {index}/{len(target_ids)}", flush=True)
    if query_api(version_query)["meta"]["dataVersion"] != version:
        raise ValueError("取得中にデータの版が変わりました。再取得してください")
    snapshot = {"schema": 1, "root": ROOT_ID, "data_version": version, "retrieved_at": datetime.now(timezone.utc).isoformat(), "source": f"https://{API_HOST}{API_PATH}", "diseases": sorted(diseases, key=lambda d: d["name"]), "records": list(records.values()), "expression": expression}
    save_snapshot(DATA_PATH, snapshot)
    print(f"保存: {DATA_PATH}\n疾患 {len(diseases)} / 薬剤 {len(drug_ids)} / 標的 {len(target_ids)}", flush=True)


if __name__ == "__main__":
    main()
