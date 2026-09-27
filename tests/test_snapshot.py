"""スナップショットの読み込みを検証する。"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from autoimmune_atlas import snapshot


class SnapshotLoadTests(unittest.TestCase):
    """古い schema のスナップショットは読み込まず、更新を求める。"""

    def test_schema1_and_schema2_require_refresh(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            for old in (1, 2):
                path.write_text(json.dumps({"schema": old}), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "Refresh"):
                    snapshot.load_snapshot(path)
