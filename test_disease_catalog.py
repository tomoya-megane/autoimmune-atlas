"""疾患の閲覧用分類が、選択対象と集計単位を変えないことを確かめる。"""

import unittest

from disease_catalog import disease_catalog, ordered_disease_ids


class DiseaseCatalogTests(unittest.TestCase):
    def test_nearest_family_keeps_psoriatic_arthritis_separate_and_all_terms_available(self):
        diseases = [
            {"id": "EFO_0009459", "name": "ACPA-positive rheumatoid arthritis", "parent_ids": ["MONDO_0008383"]},
            {"id": "MONDO_0011849", "name": "psoriatic arthritis", "parent_ids": ["MONDO_0008383"]},
            {"id": "MONDO_0008383", "name": "rheumatoid arthritis", "parent_ids": []},
            {"id": "NEW_TERM", "name": "New term", "parent_ids": ["MISSING_PARENT"]},
            {"id": "CYCLE", "name": "Cyclic term", "parent_ids": ["CYCLE"]},
        ]
        catalog = disease_catalog({"diseases": diseases})
        families = {f["id"]: f for g in catalog for f in g["families"]}
        self.assertEqual([d["id"] for d in families["MONDO_0008383"]["diseases"]], ["MONDO_0008383", "EFO_0009459"])
        self.assertEqual([d["id"] for d in families["MONDO_0011849"]["diseases"]], ["MONDO_0011849"])
        all_ids = [d["id"] for g in catalog for f in g["families"] for d in f["diseases"]]
        self.assertEqual(set(all_ids), {d["id"] for d in diseases})
        self.assertEqual(len(all_ids), len(diseases))
        self.assertEqual(ordered_disease_ids({"diseases": diseases}, ["EFO_0009459", "MONDO_0008383", "MONDO_0008383", "INVALID"]), ["MONDO_0008383", "EFO_0009459"])
        self.assertEqual(ordered_disease_ids({"diseases": diseases}, ["MONDO_0008383"]), ["MONDO_0008383"])


if __name__ == "__main__":
    unittest.main()
