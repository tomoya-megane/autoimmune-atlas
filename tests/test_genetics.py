"""snapshot の関連遺伝子の集計を検証する。"""

import unittest
from typing import cast

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas import genetics
from autoimmune_atlas.models import Snapshot
from tests.core_fixture import core_snapshot


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


class SummarizeTests(unittest.TestCase):
    base: Snapshot = snapshot()

    def test_zero_and_small_thresholds_keep_every_scored_gene(self) -> None:
        for threshold in (0, 0.02):
            genes = genetics.genes_for_disease(self.base, "D1", threshold)
            self.assertEqual([g["target_id"] for g in genes], ["G1", "G2", "G9", "G3"])
        rows = genetics.summarize_genes(
            self.base, 0, 0.5, level="group", disease_ids=["D1"]
        )
        t_cell = next(r for r in rows if r["cell_id"] == "CL_0000084")
        self.assertEqual(t_cell["denominator"], 4)

    def test_precomputed_metadata_and_catalog_give_same_rows(self) -> None:
        catalog = atlas.cell_catalog(self.base, "group")
        cell_ids = [cell["id"] for cell in catalog][:2]
        diseases = ["D1", "D2"]
        plain = genetics.summarize_genes(
            self.base,
            0.5,
            0.5,
            level="group",
            cell_ids=cell_ids,
            disease_ids=diseases,
        )
        cached = genetics.summarize_genes(
            self.base,
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
        genes = genetics.genes_for_disease(self.base, "D1", 0.5)
        self.assertEqual([g["target_id"] for g in genes], ["G1", "G2", "G9"])
        self.assertEqual(genetics.genes_for_disease(self.base, "D2", 0.5), [])
        self.assertEqual(genetics.genes_for_disease(self.base, "MISSING", 0.5), [])

    def test_counts_share_expression_rules_and_zero_genes_give_na_percent(self):
        rows = genetics.summarize_genes(
            self.base, 0.5, 0.5, method="fixed", level="group"
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
        base = snapshot()
        base["associations"]["D1"].append(
            {
                "target_id": "G_NOEXPR",
                "target": "No expr",
                "score": 0.7,
                "datasource_scores": {},
            }
        )
        rows = genetics.summarize_genes(base, 0.5, 0.5, level="group")
        t_cell = next(
            r for r in rows if r["disease_id"] == "D1" and r["cell_id"] == "CL_0000084"
        )
        self.assertEqual(t_cell["count"], 3)
        self.assertEqual(t_cell["unknown"], 1)
        self.assertEqual(t_cell["denominator"], 4)
        self.assertEqual(t_cell["status"], "partial")

    def test_disease_and_cell_filters_and_all_level(self) -> None:
        rows = genetics.summarize_genes(
            self.base, 0.5, 0.5, level="all", disease_ids=["D1"]
        )
        self.assertEqual(
            [(r["disease_id"], r["cell_id"]) for r in rows], [("D1", "all")]
        )
        rows = genetics.summarize_genes(
            self.base, 0.5, 0.5, level="mixed", cell_ids=["group:CL_0000084"]
        )
        self.assertEqual({r["cell_id"] for r in rows}, {"group:CL_0000084"})

    def test_invalid_arguments_raise_before_expression_lookup(self) -> None:
        # D2 は遺伝子を持たないので、検証が発現の判定より前に行われることを確かめられる
        with self.assertRaises(ValueError):
            genetics.summarize_genes(
                self.base, 0.5, 0.5, method="bogus", disease_ids=["D2"]
            )
        with self.assertRaises(ValueError):
            genetics.summarize_genes(self.base, 1.5, 0.5, disease_ids=["D2"])

    def test_disease_missing_from_associations_is_unavailable(self) -> None:
        base = snapshot()
        del base["associations"]["D2"]
        rows = genetics.summarize_genes(
            base, 0.5, 0.5, level="group", disease_ids=["D1", "D2"]
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
