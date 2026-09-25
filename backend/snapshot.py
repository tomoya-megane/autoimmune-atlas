"""保存済みスナップショットを検証して読み込む。"""

import json
from pathlib import Path
from typing import cast

from backend.models import Snapshot

DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "snapshot.json"


def load_snapshot(path: Path = DATA_PATH) -> Snapshot | None:
    """schema 2 の保存済みスナップショットを読む。"""
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as handle:
        raw = cast(object, json.load(handle))
    if not isinstance(raw, dict):
        raise ValueError("snapshot.json must contain a JSON object")
    # JSON object と確認済み。以降は既存の schema 2 契約を検証する。
    snapshot = cast(dict[str, object], raw)
    if snapshot.get("schema") != 2:
        raise ValueError(
            "Unsupported snapshot format. Refresh it with: pixi run refresh"
        )
    required = {"root", "retrieved_at", "source", "diseases", "records", "expression"}
    missing = required - snapshot.keys()
    if missing:
        raise ValueError(
            "Missing required fields in snapshot.json: " + ", ".join(sorted(missing))
        )
    if not isinstance(snapshot["diseases"], list) or not isinstance(
        snapshot["records"], list
    ):
        raise ValueError("diseases and records must be arrays")
    if not isinstance(snapshot["expression"], dict):
        raise ValueError("expression must be an object")
    # このファイルは refresh が生成する。上の既存検証後は schema 2 として扱う。
    return cast(Snapshot, cast(object, snapshot))
