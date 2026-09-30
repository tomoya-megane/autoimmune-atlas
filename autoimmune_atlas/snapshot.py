"""保存済みスナップショットを検証して読み込む。"""

from pathlib import Path

import msgspec

from autoimmune_atlas.models import Snapshot

DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "snapshot.json"
REFRESH_HINT = "Refresh it with: pixi run refresh"


def load_snapshot(path: Path = DATA_PATH) -> Snapshot | None:
    """schema 4 の保存済みスナップショットを読む。無ければ None、形が違えば ValueError。"""
    if not path.exists():
        return None
    try:
        snapshot = msgspec.json.decode(path.read_bytes(), type=Snapshot)
    except msgspec.ValidationError as error:
        # schema 2 のファイルは cells が、schema 3 のファイルは targets が無いのでここに来る。
        raise ValueError(
            f"Unsupported snapshot format ({error}). {REFRESH_HINT}"
        ) from error
    except msgspec.DecodeError as error:
        raise ValueError(f"snapshot.json is not valid JSON: {error}") from error
    if snapshot["schema"] != 4:
        raise ValueError(
            f"Unsupported snapshot schema {snapshot['schema']}. {REFRESH_HINT}"
        )
    return snapshot
