"""genetic association の取得と datasource の判定を検証する。"""

import unittest
from collections.abc import Mapping
from unittest.mock import patch

from autoimmune_atlas.refresh_genetics import (
    classify_datasources,
    collect_associations,
    restrict_datasources,
)


def _row(
    target: str, genetic: float | None, extra: dict[str, float]
) -> dict[str, object]:
    datatype = (
        [] if genetic is None else [{"id": "genetic_association", "score": genetic}]
    )
    datasource = [{"id": key, "score": value} for key, value in extra.items()]
    return {
        "target": {"id": target, "approvedSymbol": target.lower()},
        "datatypeScores": datatype + [{"id": "literature", "score": 0.3}],
        "datasourceScores": datasource,
    }


PAGES: dict[tuple[str, int], list[dict[str, object]]] = {
    ("D1", 0): [
        _row("G1", 0.9, {"gwas_credible_sets": 0.9, "europepmc": 0.3}),
        _row("G2", 0.5, {"eva": 0.5}),
    ],
    ("D1", 1): [
        _row("G3", 0.05, {"eva": 0.05}),
        _row("G4", None, {"europepmc": 0.3}),
    ],
    ("D2", 0): [_row("G5", None, {"europepmc": 0.3})],
}


def fake_query(query: str, variables: Mapping[str, object]) -> dict[str, object]:
    if "associatedTargets" in query:
        disease = str(variables["id"])
        index = int(str(variables["index"]))
        total = sum(len(v) for (d, _), v in PAGES.items() if d == disease)
        return {
            "disease": {
                "associatedTargets": {
                    "count": total,
                    "rows": PAGES.get((disease, index), []),
                }
            }
        }
    if "evidences" in query:
        datasource = str(variables["datasource"])
        datatype = "literature" if datasource == "europepmc" else "genetic_association"
        return {
            "disease": {
                "evidences": {
                    "rows": [{"datasourceId": datasource, "datatypeId": datatype}]
                }
            }
        }
    raise AssertionError(query)


class RefreshGeneticsTests(unittest.TestCase):
    def test_paging_stops_below_floor_and_keeps_scored_genes(self) -> None:
        associations = collect_associations(fake_query, ["D1", "D2"], 0.1)
        self.assertEqual([g["target_id"] for g in associations["D1"]], ["G1", "G2"])
        self.assertEqual(associations["D2"], [])
        self.assertEqual(associations["D1"][0]["target"], "g1")

    def test_page_size_is_two_in_this_test_only_via_patch(self) -> None:
        # PAGE_SIZE を 2 にして 2 ページ目まで読むことを確かめる。
        with patch("autoimmune_atlas.refresh_genetics.PAGE_SIZE", 2):
            associations = collect_associations(fake_query, ["D1"], 0.1)
        self.assertEqual(len(associations["D1"]), 2)

    def test_datasources_are_classified_by_one_evidence_each(self) -> None:
        associations = collect_associations(fake_query, ["D1"], 0.1)
        datasources = classify_datasources(fake_query, associations)
        self.assertEqual(datasources, ["eva", "gwas_credible_sets"])
        restricted = restrict_datasources(associations, datasources)
        self.assertEqual(
            restricted["D1"][0]["datasource_scores"], {"gwas_credible_sets": 0.9}
        )

    def test_missing_disease_raises(self) -> None:
        def missing(_query: str, _variables: Mapping[str, object]) -> dict[str, object]:
            return {"disease": None}

        with self.assertRaises(ValueError):
            collect_associations(missing, ["D1"], 0.1)
