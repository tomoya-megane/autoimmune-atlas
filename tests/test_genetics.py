"""genetics.json の読み込みと遺伝子の集計を検証する。"""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import cast, override

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas import genetics
from autoimmune_atlas.models import GeneticsSnapshot, Snapshot
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
            {"cell_id": "T4", "median": 99.0, "specificity_score": None}
        ]
        merged = genetics.merged_snapshot(base, data)
        self.assertEqual(merged["expression"]["G1"][0]["median"], 2.0)
        self.assertIn("G9", merged["expression"])
        self.assertNotIn("G9", base["expression"])


class SummarizeTests(unittest.TestCase):
    base: Snapshot = snapshot()
    data: GeneticsSnapshot = genetics_snapshot()

    @override
    def setUp(self) -> None:
        self.base = genetics.merged_snapshot(snapshot(), genetics_snapshot())
        self.data = genetics_snapshot()

    def test_precomputed_metadata_and_catalog_give_same_rows(self) -> None:
        catalog = atlas.cell_catalog(self.base, "group")
        cell_ids = [cell["id"] for cell in catalog][:2]
        diseases = ["D1", "D2"]
        plain = genetics.summarize_genes(
            self.base,
            self.data,
            0.5,
            0.5,
            level="group",
            cell_ids=cell_ids,
            disease_ids=diseases,
        )
        cached = genetics.summarize_genes(
            self.base,
            self.data,
            0.5,
            0.5,
            metadata=atlas.expression_metadata(self.base),
            level="group",
            cell_ids=cell_ids,
            disease_ids=diseases,
            catalog=catalog,
        )
        self.assertTrue(plain)
        self.assertEqual(cached, plain)
        # 渡した catalog は cell_ids の絞り込みで書き換えない。
        self.assertEqual(len(catalog), len(atlas.cell_catalog(self.base, "group")))

    def test_threshold_is_inclusive_and_sorted_by_score(self) -> None:
        genes = genetics.genes_for_disease(self.data, "D1", 0.5)
        self.assertEqual([g["target_id"] for g in genes], ["G1", "G2", "G9"])
        self.assertEqual(genetics.genes_for_disease(self.data, "D2", 0.5), [])
        self.assertEqual(genetics.genes_for_disease(self.data, "MISSING", 0.5), [])

    def test_counts_share_expression_rules_and_zero_genes_give_na_percent(self):
        rows = genetics.summarize_genes(
            self.base, self.data, 0.5, 0.5, method="fixed", level="group"
        )
        t_cell = next(
            r for r in rows if r["disease_id"] == "D1" and r["cell_id"] == "CL_0000084"
        )
        # T4 と T8 は ancestor_ids により T cell（CL_0000084）の大分類にまとまる。
        # G1 (T4 2.0), G2 (T8 3.0), G9 (T4 5.0) はその大分類のどこかで 0.5 以上
        self.assertEqual(t_cell["count"], 3)
        self.assertEqual(t_cell["denominator"], 3)
        self.assertEqual(t_cell["percent"], 100.0)
        self.assertEqual(t_cell["status"], "complete")
        self.assertEqual(t_cell["drug_denominator"], 0)
        self.assertEqual(t_cell["records"], [])
        empty = next(
            r for r in rows if r["disease_id"] == "D2" and r["cell_id"] == "CL_0000084"
        )
        self.assertEqual(empty["count"], 0)
        self.assertIsNone(empty["percent"])
        self.assertEqual(empty["denominator"], 0)
        self.assertEqual(empty["status"], "complete")

    def test_gene_without_expression_is_unknown_not_error(self) -> None:
        data = genetics_snapshot()
        data["associations"]["D1"].append(
            {
                "target_id": "G_NOEXPR",
                "target": "No expr",
                "score": 0.7,
                "datasource_scores": {},
            }
        )
        rows = genetics.summarize_genes(self.base, data, 0.5, 0.5, level="group")
        t_cell = next(
            r for r in rows if r["disease_id"] == "D1" and r["cell_id"] == "CL_0000084"
        )
        self.assertEqual(t_cell["count"], 3)
        self.assertEqual(t_cell["unknown"], 1)
        self.assertEqual(t_cell["denominator"], 4)
        self.assertEqual(t_cell["status"], "partial")

    def test_disease_and_cell_filters_and_all_level(self) -> None:
        rows = genetics.summarize_genes(
            self.base, self.data, 0.5, 0.5, level="all", disease_ids=["D1"]
        )
        self.assertEqual(
            [(r["disease_id"], r["cell_id"]) for r in rows], [("D1", "all")]
        )
        rows = genetics.summarize_genes(
            self.base, self.data, 0.5, 0.5, level="mixed", cell_ids=["group:CL_0000084"]
        )
        self.assertEqual({r["cell_id"] for r in rows}, {"group:CL_0000084"})

    def test_invalid_arguments_raise_before_expression_lookup(self) -> None:
        # D2 は遺伝子を持たないので、検証が発現の判定より前に行われることを確かめられる
        with self.assertRaises(ValueError):
            genetics.summarize_genes(
                self.base, self.data, 0.5, 0.5, method="bogus", disease_ids=["D2"]
            )
        with self.assertRaises(ValueError):
            genetics.summarize_genes(self.base, self.data, 1.5, 0.5, disease_ids=["D2"])

    def test_disease_missing_from_associations_is_unavailable(self) -> None:
        data = genetics_snapshot()
        del data["associations"]["D2"]
        rows = genetics.summarize_genes(
            self.base, data, 0.5, 0.5, level="group", disease_ids=["D1", "D2"]
        )
        missing = next(
            r for r in rows if r["disease_id"] == "D2" and r["cell_id"] == "CL_0000084"
        )
        self.assertIsNone(missing["count"])
        self.assertIsNone(missing["percent"])
        self.assertEqual(missing["denominator"], 0)
        self.assertEqual(missing["status"], "unavailable")
        t_cell = next(
            r for r in rows if r["disease_id"] == "D1" and r["cell_id"] == "CL_0000084"
        )
        self.assertEqual(t_cell["count"], 3)
        self.assertEqual(t_cell["status"], "complete")
