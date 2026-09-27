"""genetic association の取得と datasource の判定を検証する。"""

import contextlib
import io
import unittest
from collections.abc import Mapping
from typing import cast
from unittest.mock import MagicMock, patch

from autoimmune_atlas.genetics import GENETICS_SCHEMA, SCORE_FLOOR
from autoimmune_atlas.models import GeneticsSnapshot
from autoimmune_atlas.refresh_genetics import (
    VERSION_QUERY,
    classify_datasources,
    collect_associations,
    main,
    restrict_datasources,
)

MODULE = "autoimmune_atlas.refresh_genetics"
VERSION = {"year": "26", "month": "09", "iteration": None}


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


PAGES: dict[tuple[str, int], list[Mapping[str, object]]] = {
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
    def test_first_page_keeps_scored_genes_and_stops_at_count(self) -> None:
        associations = collect_associations(fake_query, ["D1", "D2"], 0.1)
        self.assertEqual([g["target_id"] for g in associations["D1"]], ["G1", "G2"])
        self.assertEqual(associations["D2"], [])
        self.assertEqual(associations["D1"][0]["target"], "g1")

    def test_paging_stops_at_first_gene_below_floor(self) -> None:
        # PAGE_SIZE を 2 にして 2 ページ目を読ませ、0.1 未満の G3 で止まることを確かめる。
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


def _snapshot() -> dict[str, object]:
    """G1 の発現だけを持つ snapshot.json の代わり。"""
    return {
        "data_version": VERSION,
        "diseases": [{"id": "D1"}],
        "expression": {"G1": []},
    }


def _query_with_versions(versions: list[Mapping[str, object]]) -> object:
    """版の問い合わせには versions を順に返し、残りは fake_query に渡す。"""

    def query(
        text: str, variables: Mapping[str, object] | None = None
    ) -> dict[str, object]:
        if text == VERSION_QUERY:
            return {"meta": {"dataVersion": versions.pop(0)}}
        return fake_query(text, variables or {})

    return query


def _fetch(target: str) -> tuple[str, list[Mapping[str, object]], dict[str, object]]:
    return target, [{"target_id": target}], {}


class RefreshGeneticsMainTests(unittest.TestCase):
    """main() が版を確かめ、不足分の発現だけを取って保存することを確かめる。"""

    def _run(
        self, versions: list[Mapping[str, object]], fetch: MagicMock, save: MagicMock
    ) -> None:
        with (
            patch(f"{MODULE}.load_snapshot", return_value=_snapshot()),
            patch(f"{MODULE}.query_api", side_effect=_query_with_versions(versions)),
            patch(f"{MODULE}.fetch_expression", fetch),
            patch(f"{MODULE}.save_snapshot", save),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            main()

    def test_start_version_mismatch_stops_without_saving(self) -> None:
        other = {"year": "26", "month": "06", "iteration": None}
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with self.assertRaises(SystemExit):
            self._run([other], fetch, save)
        fetch.assert_not_called()
        save.assert_not_called()

    def test_fetches_only_missing_expression_and_saves_payload(self) -> None:
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        self._run([VERSION, VERSION], fetch, save)
        self.assertEqual([c.args[0] for c in fetch.call_args_list], ["G2"])
        save.assert_called_once()
        payload = cast(GeneticsSnapshot, save.call_args.args[1])
        self.assertEqual(payload["schema"], GENETICS_SCHEMA)
        self.assertEqual(payload["score_floor"], SCORE_FLOOR)
        self.assertEqual(payload["data_version"], VERSION)
        self.assertEqual(payload["datasources"], ["eva", "gwas_credible_sets"])
        self.assertEqual(
            [
                (g["target_id"], g["datasource_scores"])
                for g in payload["associations"]["D1"]
            ],
            [("G1", {"gwas_credible_sets": 0.9}), ("G2", {"eva": 0.5})],
        )
        self.assertEqual(list(payload["expression"]), ["G2"])

    def test_version_change_during_refresh_raises_and_saves_nothing(self) -> None:
        later = {"year": "26", "month": "12", "iteration": None}
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with self.assertRaises(ValueError):
            self._run([VERSION, later], fetch, save)
        save.assert_not_called()
