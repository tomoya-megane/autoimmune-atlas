"""疾患の閲覧用分類が、選択対象と集計単位を変えないことを確かめる。"""

import unittest

from backend.disease_catalog import (
    DISEASE_GROUPS,
    SCOPE_ROOTS,
    disease_catalog,
    ordered_disease_ids,
)
from backend.models import CatalogDisease, DiseaseCatalogInput


class ScopeRootTests(unittest.TestCase):
    def test_roots_are_unique_and_start_with_autoimmune_disease(self):
        ids = [root_id for root_id, _ in SCOPE_ROOTS]
        self.assertEqual(ids[0], "MONDO_0007179")
        self.assertEqual(len(ids), len(set(ids)))

    def test_every_added_root_has_a_browsing_family(self):
        family_ids = {root for _, _, ids in DISEASE_GROUPS for root in ids}
        missing = [
            root_id for root_id, _ in SCOPE_ROOTS[1:] if root_id not in family_ids
        ]
        self.assertEqual(missing, [])


class DiseaseCatalogTests(unittest.TestCase):
    def test_groups_and_families_follow_names_with_other_last(self):
        diseases: list[CatalogDisease] = [
            {
                "id": "MONDO_0007915",
                "name": "systemic lupus erythematosus",
                "parent_ids": [],
            },
            {"id": "MONDO_0008383", "name": "rheumatoid arthritis", "parent_ids": []},
            {"id": "MONDO_0011849", "name": "psoriatic arthritis", "parent_ids": []},
            {
                "id": "MONDO_0005147",
                "name": "type 1 diabetes mellitus",
                "parent_ids": [],
            },
            {"id": "OTHER", "name": "unclassified term", "parent_ids": []},
        ]
        snapshot: DiseaseCatalogInput = {"diseases": diseases}
        catalog = disease_catalog(snapshot)
        self.assertEqual(
            [group["id"] for group in catalog],
            ["endocrine", "joints", "systemic", "other"],
        )
        self.assertEqual(
            [family["label"] for family in catalog[1]["families"]],
            ["psoriatic arthritis", "rheumatoid arthritis"],
        )
        self.assertEqual(
            ordered_disease_ids(snapshot, [d["id"] for d in diseases]),
            [
                "MONDO_0005147",
                "MONDO_0011849",
                "MONDO_0008383",
                "MONDO_0007915",
                "OTHER",
            ],
        )

    def test_nearest_family_keeps_psoriatic_arthritis_separate_and_all_terms_available(
        self,
    ):
        diseases: list[CatalogDisease] = [
            {
                "id": "EFO_0009459",
                "name": "ACPA-positive rheumatoid arthritis",
                "parent_ids": ["MONDO_0008383"],
            },
            {
                "id": "MONDO_0011849",
                "name": "psoriatic arthritis",
                "parent_ids": ["MONDO_0008383"],
            },
            {"id": "MONDO_0008383", "name": "rheumatoid arthritis", "parent_ids": []},
            {"id": "NEW_TERM", "name": "New term", "parent_ids": ["MISSING_PARENT"]},
            {"id": "CYCLE", "name": "Cyclic term", "parent_ids": ["CYCLE"]},
        ]
        catalog = disease_catalog({"diseases": diseases})
        families = {f["id"]: f for g in catalog for f in g["families"]}
        self.assertEqual(
            [d["id"] for d in families["MONDO_0008383"]["diseases"]],
            ["MONDO_0008383", "EFO_0009459"],
        )
        self.assertEqual(
            [d["id"] for d in families["MONDO_0011849"]["diseases"]], ["MONDO_0011849"]
        )
        all_ids = [
            d["id"] for g in catalog for f in g["families"] for d in f["diseases"]
        ]
        self.assertEqual(set(all_ids), {d["id"] for d in diseases})
        self.assertEqual(len(all_ids), len(diseases))
        self.assertEqual(
            ordered_disease_ids(
                {"diseases": diseases},
                ["EFO_0009459", "MONDO_0008383", "MONDO_0008383", "INVALID"],
            ),
            ["MONDO_0008383", "EFO_0009459"],
        )
        self.assertEqual(
            ordered_disease_ids({"diseases": diseases}, ["MONDO_0008383"]),
            ["MONDO_0008383"],
        )


if __name__ == "__main__":
    unittest.main()
