"""genetics.json の読み込みと遺伝子の集計を検証する。"""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import cast

from autoimmune_atlas import genetics
from autoimmune_atlas.models import Snapshot
from tests.core_fixture import core_snapshot
from tests.genetics_fixture import genetics_snapshot


def snapshot() -> Snapshot:
    core = core_snapshot()
    return cast(
        Snapshot,
        cast(
            object,
            {
                **core,
                "root": "MONDO_0007179",
                "data_version": {"year": "26", "month": "09", "iteration": None},
                "retrieved_at": "2026-09-27T00:00:00+00:00",
                "source": "https://api.platform.opentargets.org/api/v4/graphql",
            },
        ),
    )


class LoadTests(unittest.TestCase):
    def test_missing_file_returns_none_and_bad_schema_raises(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "genetics.json"
            self.assertIsNone(genetics.load_genetics(path))
            path.write_text(json.dumps({"schema": 2}), encoding="utf-8")
            with self.assertRaises(ValueError):
                genetics.load_genetics(path)
            path.write_text(json.dumps(genetics_snapshot()), encoding="utf-8")
            loaded = genetics.load_genetics(path)
            assert loaded is not None
            self.assertEqual(loaded["score_floor"], 0.1)

    def test_version_match_and_merged_expression_prefer_snapshot(self) -> None:
        base = snapshot()
        data = genetics_snapshot()
        self.assertTrue(genetics.version_matches(base, data))
        data["data_version"] = {"year": "26", "month": "06", "iteration": None}
        self.assertFalse(genetics.version_matches(base, data))
        data = genetics_snapshot()
        data["expression"]["G1"] = [
            {"cell_id": "T4", "cell": "CD4 T cell", "median": 99.0}
        ]
        merged = genetics.merged_snapshot(base, data)
        self.assertEqual(merged["expression"]["G1"][0]["median"], 2.0)
        self.assertIn("G9", merged["expression"])
        self.assertNotIn("G9", base["expression"])
