"""スナップショットの読み込みを検証する。"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from autoimmune_atlas import snapshot
from tests.ui_fixture import fixture


class SnapshotLoadTests(unittest.TestCase):
    """古い schema のスナップショットは読み込まず、更新を求める。"""

    def test_schema1_and_schema2_require_refresh(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            for old in (1, 2):
                path.write_text(json.dumps({"schema": old}), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "Refresh"):
                    snapshot.load_snapshot(path)

    def test_schema4_loads_and_schema3_requires_refresh(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            data = fixture()
            path.write_text(json.dumps(data), encoding="utf-8")
            loaded = snapshot.load_snapshot(path)
            assert loaded is not None
            self.assertEqual(loaded["targets"], data["targets"])
            path.write_text(json.dumps({**data, "schema": 3}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "schema 3"):
                snapshot.load_snapshot(path)
            # schema 3 の実ファイルは targets を持たない。
            old = {key: value for key, value in data.items() if key != "targets"}
            path.write_text(json.dumps({**old, "schema": 3}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Refresh"):
                snapshot.load_snapshot(path)
