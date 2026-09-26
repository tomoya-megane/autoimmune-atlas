"""genetic association の保存データを読み、遺伝子ページ用に集計する。"""

import json
from pathlib import Path
from typing import cast

from autoimmune_atlas.models import GeneticsSnapshot, Snapshot

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
