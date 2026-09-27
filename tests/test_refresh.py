"""公開 API からの取得、正規化、保存を検証する。"""

import contextlib
import io
import json
import unittest
from collections.abc import Callable, Mapping
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import cast
from unittest.mock import MagicMock, patch

import httpx

from autoimmune_atlas.models import CellDefinition, DataVersion, ExpressionRow
from autoimmune_atlas.refresh import (
    VERSION_QUERY,
    ClinicalCandidate,
    ExpressionApiRow,
    MechanismApiRow,
    classify_datasources,
    collect_associations,
    extract_expression,
    extract_mechanisms,
    main,
    merge_cells,
    normalize_drugs,
    query_api,
    resolve_disease_ids,
    restrict_datasources,
    reusable_expression,
    save_snapshot,
)


class RefreshTests(unittest.TestCase):
    def test_scope_is_the_union_of_roots_and_body_only_roots_add_no_descendants(
        self,
    ) -> None:
        tree = {
            "MONDO_0007179": ["A", "B"],
            "MONDO_0003346": ["AIDS_ARTERITIS", "B"],
        }

        def fake_query(_query: str, variables: dict[str, str]) -> dict[str, object]:
            root_id = variables["id"]
            if root_id not in tree:
                return (
                    {"disease": None}
                    if root_id == "MISSING"
                    else {
                        "disease": {"id": root_id, "name": root_id, "descendants": []}
                    }
                )
            return {
                "disease": {
                    "id": root_id,
                    "name": root_id,
                    "descendants": tree[root_id],
                }
            }

        roots, ids = resolve_disease_ids(fake_query)
        self.assertEqual(roots[0]["id"], "MONDO_0007179")
        self.assertEqual(roots[0]["count"], 3)
        body_only = next(r for r in roots if r["id"] == "MONDO_0003346")
        self.assertFalse(body_only["include_descendants"])
        self.assertEqual(body_only["count"], 1)
        self.assertIn("A", ids)
        self.assertIn("MONDO_0003346", ids)
        self.assertNotIn("AIDS_ARTERITIS", ids)
        self.assertEqual(ids.count("B"), 1)
        with (
            patch("autoimmune_atlas.refresh.SCOPE_ROOTS", (("MISSING", True),)),
            self.assertRaises(ValueError),
        ):
            resolve_disease_ids(fake_query)

    def test_phase_one_drugs_and_parent_molecule_are_preserved(self) -> None:
        kinds = (
            "Small molecule",
            "Antibody",
            "Protein",
            "Cell",
            "Gene",
            "Enzyme",
            "Oligonucleotide",
            "Antibody drug conjugate",
            "Vaccine component",
            "Oligosaccharide",
            "Unknown",
        )
        modalities = (
            "small_molecule",
            "antibody",
            "protein",
            "cell",
            "gene",
            "enzyme",
            "oligonucleotide",
            "antibody_drug_conjugate",
            "vaccine_component",
            "oligosaccharide",
            "unknown",
        )
        rows: list[ClinicalCandidate] = [
            {
                "maxClinicalStage": "PHASE_1",
                "drug": {
                    "id": str(index),
                    "name": kind,
                    "drugType": kind,
                    "parentMolecule": (
                        {"id": "P", "name": "Parent"} if index == 0 else None
                    ),
                },
            }
            for index, kind in enumerate(kinds)
        ]
        drugs = normalize_drugs(rows)
        self.assertEqual([drug["drug_type"] for drug in drugs], list(kinds))
        self.assertEqual([drug["modality"] for drug in drugs], list(modalities))
        self.assertEqual(drugs[0]["canonical_drug_id"], "P")
        self.assertEqual(drugs[1]["canonical_drug_id"], "1")

    def test_withdrawal_is_not_stored_but_approval_and_phase_four_are(self) -> None:
        rows: list[ClinicalCandidate] = [
            {
                "maxClinicalStage": stage,
                "drug": {
                    "id": stage,
                    "name": stage,
                    "drugType": "Antibody",
                    "parentMolecule": None,
                },
            }
            for stage in ("WITHDRAWAL", "APPROVAL", "PHASE_4")
        ]
        self.assertEqual(
            [row["stage"] for row in normalize_drugs(rows)], ["APPROVAL", "PHASE_4"]
        )

    def test_duplicate_drug_rows_are_rejected(self) -> None:
        row: ClinicalCandidate = {
            "maxClinicalStage": "PHASE_1",
            "drug": {
                "id": "A",
                "name": "Alpha",
                "drugType": "Antibody",
                "parentMolecule": None,
            },
        }
        with self.assertRaisesRegex(ValueError, "重複"):
            normalize_drugs([row, row])

    def test_expression_splits_values_from_cell_definitions(self) -> None:
        base: ExpressionApiRow = {
            "datasourceId": "tabula_sapiens",
            "unit": "CPM(pseudobulk sum[counts])",
            "median": 2.0,
            "specificity_score": 0.8,
            "celltypeBiosample": {
                "biosampleId": "C1",
                "biosampleName": "B cell",
                "ancestors": ["CL_1", "CL_2"],
            },
            "celltypeBiosampleParent": {
                "biosampleId": "P1",
                "biosampleName": "Lymphocyte",
            },
            "tissueBiosample": None,
        }
        rows, cells = extract_expression([base])
        self.assertEqual(
            rows, [{"cell_id": "C1", "median": 2.0, "specificity_score": 0.8}]
        )
        self.assertEqual(
            cells,
            {
                "C1": {
                    "name": "B cell",
                    "parent_id": "P1",
                    "parent": "Lymphocyte",
                    "ancestor_ids": ["CL_1", "CL_2"],
                }
            },
        )
        self.assertEqual(
            extract_expression(
                [
                    {**base, "datasourceId": "DICE"},
                    {**base, "tissueBiosample": {"biosampleId": "T1"}},
                ]
            ),
            ([], {}),
        )
        invalid_rows: tuple[tuple[str, ExpressionApiRow], ...] = (
            ("median", {**base, "median": -1}),
            ("specificity_score", {**base, "specificity_score": 1.1}),
        )
        for key, invalid_row in invalid_rows:
            with self.subTest(key=key), self.assertRaises(ValueError):
                extract_expression([invalid_row])
        invalid_unit: ExpressionApiRow = {**base, "unit": "TPM"}
        with self.assertRaisesRegex(ValueError, "単位"):
            extract_expression([invalid_unit])

    def test_merge_cells_rejects_inconsistent_definitions(self) -> None:
        cells: dict[str, CellDefinition] = {}
        first: CellDefinition = {
            "name": "B cell",
            "parent_id": "P",
            "parent": "Lymphocyte",
            "ancestor_ids": ["CL_1"],
        }
        merge_cells(cells, {"C1": first})
        merge_cells(cells, {"C1": first, "C2": {**first, "name": "T cell"}})
        self.assertEqual(set(cells), {"C1", "C2"})
        with self.assertRaisesRegex(ValueError, "一貫"):
            merge_cells(cells, {"C1": {**first, "ancestor_ids": []}})

    def test_mechanism_references_are_deduplicated(self) -> None:
        row: MechanismApiRow = {
            "mechanismOfAction": "blocks",
            "actionType": "BLOCKER",
            "targets": [{"id": "G1", "approvedSymbol": "GENE1"}],
            "references": [
                {"source": "FDA", "ids": ["1"], "urls": ["https://example.test"]},
                {"source": "FDA", "ids": ["1"], "urls": ["https://example.test"]},
            ],
        }
        self.assertEqual(
            extract_mechanisms([row]),
            [
                {
                    "target_id": "G1",
                    "target": "GENE1",
                    "mechanism": "blocks",
                    "action_types": ["BLOCKER"],
                    "references": [
                        {
                            "source": "FDA",
                            "ids": ["1"],
                            "urls": ["https://example.test"],
                        }
                    ],
                }
            ],
        )

    def test_failed_serialization_preserves_previous_snapshot(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            path.write_text('{"previous": true}')
            with self.assertRaises(ValueError):
                save_snapshot(path, {"value": float("nan")})
            self.assertEqual(path.read_text(), '{"previous": true}')


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


ROWS: dict[str, list[Mapping[str, object]]] = {
    "D1": [
        _row("G1", 0.9, {"gwas_credible_sets": 0.9, "europepmc": 0.3}),
        _row("G2", 0.5, {"eva": 0.5}),
        _row("G3", 0.05, {"eva": 0.05}),
        _row("G4", None, {"europepmc": 0.3}),
        _row("G5", None, {"europepmc": 0.2}),
    ],
    "D2": [_row("G6", None, {"europepmc": 0.3})],
}


def fake_query(query: str, variables: Mapping[str, object]) -> dict[str, object]:
    if "associatedTargets" in query:
        disease = str(variables["id"])
        index = int(str(variables["index"]))
        size = int(str(variables["size"]))
        rows = ROWS.get(disease, [])
        return {
            "disease": {
                "associatedTargets": {
                    "count": len(rows),
                    "rows": rows[index * size : (index + 1) * size],
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


class AssociationTests(unittest.TestCase):
    def test_every_scored_gene_is_kept_and_unscored_genes_stop_paging(self) -> None:
        associations = collect_associations(fake_query, ["D1", "D2"])
        # G3 は 0.05 でも残る。スコアを持たない G4 以降は入らない。
        self.assertEqual(
            [g["target_id"] for g in associations["D1"]], ["G1", "G2", "G3"]
        )
        self.assertEqual(associations["D2"], [])
        self.assertEqual(associations["D1"][0]["target"], "g1")

    def test_paging_stops_on_the_page_that_holds_the_first_unscored_gene(self) -> None:
        # PAGE_SIZE を 2 にすると D1 は 3 ページになる。G4 のある 2 ページ目で止まり、3 ページ目は読まない。
        calls: list[int] = []

        def counting(query: str, variables: Mapping[str, object]) -> dict[str, object]:
            if "associatedTargets" in query:
                calls.append(int(str(variables["index"])))
            return fake_query(query, variables)

        with patch("autoimmune_atlas.refresh.PAGE_SIZE", 2):
            associations = collect_associations(counting, ["D1"])
        self.assertEqual(
            [g["target_id"] for g in associations["D1"]], ["G1", "G2", "G3"]
        )
        self.assertEqual(calls, [0, 1])

    def test_empty_page_before_count_is_reached_fails(self) -> None:
        def truncated(query: str, variables: Mapping[str, object]) -> dict[str, object]:
            response = fake_query(query, variables)
            if "associatedTargets" in query:
                block = cast(
                    dict[str, object],
                    cast(dict[str, object], response["disease"])["associatedTargets"],
                )
                block["count"] = 10
                if int(str(variables["index"])) > 0:
                    block["rows"] = []
            return response

        # PAGE_SIZE=3 なら 1 ページ目は G1、G2、G3 で全部スコアがあり、
        # 2 ページ目が空で count 10 に達していないので失敗する。
        with (
            patch("autoimmune_atlas.refresh.PAGE_SIZE", 3),
            self.assertRaisesRegex(ValueError, "途中のページが空"),
        ):
            collect_associations(truncated, ["D1"])

    def test_empty_first_page_below_page_size_fails(self) -> None:
        # count が 1 ページに収まる疾患で 1 ページ目が空で返っても、空の配列を保存しない。
        def empty(query: str, variables: Mapping[str, object]) -> dict[str, object]:
            if "associatedTargets" in query:
                return {"disease": {"associatedTargets": {"count": 10, "rows": []}}}
            return fake_query(query, variables)

        with self.assertRaisesRegex(ValueError, "途中のページが空"):
            collect_associations(empty, ["D1"])

    def test_datasources_are_classified_by_one_evidence_each(self) -> None:
        associations = collect_associations(fake_query, ["D1"])
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
            collect_associations(missing, ["D1"])


VERSION: DataVersion = {"year": "26", "month": "09", "iteration": None}
MODULE = "autoimmune_atlas.refresh"
CELL: CellDefinition = {
    "name": "B cell",
    "parent_id": None,
    "parent": None,
    "ancestor_ids": [],
}
SCHEMA2_ROW: dict[str, object] = {
    "cell_id": "C1",
    "cell": "B cell",
    "median": 2.0,
    "specificity_score": 0.8,
    "parent_id": None,
    "parent": None,
    "ancestor_ids": [],
}


class ReuseTests(unittest.TestCase):
    def test_reuses_same_version_rows_from_schema2_and_genetics_files(self) -> None:
        with TemporaryDirectory() as directory:
            snapshot_path = Path(directory) / "snapshot.json"
            genetics_path = Path(directory) / "genetics.json"
            snapshot_path.write_text(
                json.dumps(
                    {
                        "schema": 2,
                        "data_version": VERSION,
                        "expression": {"G1": [SCHEMA2_ROW]},
                    }
                ),
                encoding="utf-8",
            )
            genetics_path.write_text(
                json.dumps(
                    {
                        "schema": 1,
                        "data_version": VERSION,
                        "expression": {
                            "G1": [{**SCHEMA2_ROW, "median": 99.0}],
                            "G2": [{**SCHEMA2_ROW, "median": 0.5}],
                        },
                    }
                ),
                encoding="utf-8",
            )
            expression, cells = reusable_expression(
                [snapshot_path, genetics_path], VERSION
            )
        # snapshot.json を先に読むので G1 は 2.0 のまま。G2 は genetics.json から来る。
        self.assertEqual(
            expression,
            {
                "G1": [{"cell_id": "C1", "median": 2.0, "specificity_score": 0.8}],
                "G2": [{"cell_id": "C1", "median": 0.5, "specificity_score": 0.8}],
            },
        )
        self.assertEqual(cells, {"C1": CELL})

    def test_cell_definitions_are_checked_even_for_genes_already_reused(self) -> None:
        with TemporaryDirectory() as directory:
            snapshot_path = Path(directory) / "snapshot.json"
            genetics_path = Path(directory) / "genetics.json"
            snapshot_path.write_text(
                json.dumps(
                    {
                        "schema": 2,
                        "data_version": VERSION,
                        "expression": {"G1": [SCHEMA2_ROW]},
                    }
                ),
                encoding="utf-8",
            )
            genetics_path.write_text(
                json.dumps(
                    {
                        "schema": 1,
                        "data_version": VERSION,
                        "expression": {"G1": [{**SCHEMA2_ROW, "cell": "T cell"}]},
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "一貫"):
                reusable_expression([snapshot_path, genetics_path], VERSION)

    def test_reuses_schema3_cells_table(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            path.write_text(
                json.dumps(
                    {
                        "schema": 3,
                        "data_version": VERSION,
                        "cells": {"C1": CELL},
                        "expression": {
                            "G1": [
                                {
                                    "cell_id": "C1",
                                    "median": 1.0,
                                    "specificity_score": None,
                                }
                            ]
                        },
                    }
                ),
                encoding="utf-8",
            )
            expression, cells = reusable_expression([path], VERSION)
        self.assertEqual(list(expression), ["G1"])
        self.assertEqual(cells, {"C1": CELL})

    def test_other_version_missing_and_broken_files_are_skipped(self) -> None:
        with TemporaryDirectory() as directory:
            other = Path(directory) / "other.json"
            other.write_text(
                json.dumps(
                    {
                        "schema": 3,
                        "data_version": {
                            "year": "26",
                            "month": "06",
                            "iteration": None,
                        },
                        "cells": {},
                        "expression": {"G1": []},
                    }
                ),
                encoding="utf-8",
            )
            broken = Path(directory) / "broken.json"
            broken.write_text("{not json", encoding="utf-8")
            missing = Path(directory) / "missing.json"
            self.assertEqual(
                reusable_expression([other, broken, missing], VERSION), ({}, {})
            )


def _api(versions: list[Mapping[str, object]]) -> Callable[..., dict[str, object]]:
    """main() が打つ 6 種類の照会に応答する。版の照会には versions を順に返す。"""

    def query(
        text: str, variables: Mapping[str, object] | None = None
    ) -> dict[str, object]:
        variables = variables or {}
        if text == VERSION_QUERY:
            return {"meta": {"dataVersion": versions.pop(0)}}
        if "descendants" in text:
            return {"disease": {"id": "R", "name": "root", "descendants": ["D1", "D2"]}}
        if "drugAndClinicalCandidates" in text:
            ids = cast(list[str], variables["ids"])
            return {
                "diseases": [
                    {
                        "id": disease_id,
                        "name": disease_id.lower(),
                        "parents": [],
                        "drugAndClinicalCandidates": {
                            "count": 1 if disease_id == "D1" else 0,
                            "rows": [
                                {
                                    "maxClinicalStage": "PHASE_3",
                                    "drug": {
                                        "id": "A",
                                        "name": "Alpha",
                                        "drugType": "Antibody",
                                        "parentMolecule": None,
                                    },
                                }
                            ]
                            if disease_id == "D1"
                            else [],
                        },
                    }
                    for disease_id in ids
                ]
            }
        if "mechanismsOfAction" in text:
            return {
                "drugs": [
                    {
                        "id": "A",
                        "mechanismsOfAction": {
                            "rows": [
                                {
                                    "mechanismOfAction": "blocks",
                                    "actionType": "BLOCKER",
                                    "targets": [{"id": "G1", "approvedSymbol": "g1"}],
                                    "references": [],
                                }
                            ]
                        },
                    }
                ]
            }
        return fake_query(text, variables)

    return query


def _fetch(target: str) -> tuple[str, list[ExpressionRow], dict[str, CellDefinition]]:
    return (
        target,
        [{"cell_id": "C1", "median": 1.0, "specificity_score": None}],
        {"C1": CELL},
    )


class MainTests(unittest.TestCase):
    """main() が関連遺伝子まで取り、同じ版の発現を再利用し、全件成功したときだけ保存することを確かめる。"""

    def _run(
        self,
        versions: list[Mapping[str, object]],
        fetch: MagicMock,
        save: MagicMock,
        directory: str,
    ) -> None:
        with (
            patch(f"{MODULE}.SCOPE_ROOTS", (("R", True),)),
            patch(f"{MODULE}.DATA_PATH", Path(directory) / "snapshot.json"),
            patch(f"{MODULE}.LEGACY_GENETICS_PATH", Path(directory) / "genetics.json"),
            patch(f"{MODULE}.query_api", side_effect=_api(versions)),
            patch(f"{MODULE}.fetch_expression", fetch),
            patch(f"{MODULE}.save_snapshot", save),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            main()

    def test_fetches_union_of_targets_and_genes_and_saves_schema3(self) -> None:
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with TemporaryDirectory() as directory:
            self._run([VERSION, VERSION], fetch, save, directory)
        # 薬剤の標的 G1 と関連遺伝子 G1、G2、G3 の和集合
        self.assertEqual(
            sorted(c.args[0] for c in fetch.call_args_list), ["G1", "G2", "G3"]
        )
        save.assert_called_once()
        payload = cast(dict[str, object], save.call_args.args[1])
        self.assertEqual(payload["schema"], 3)
        self.assertEqual(payload["data_version"], VERSION)
        self.assertEqual(payload["datasources"], ["eva", "gwas_credible_sets"])
        self.assertEqual(payload["cells"], {"C1": CELL})
        associations = cast(dict[str, list[dict[str, object]]], payload["associations"])
        self.assertEqual(
            [(g["target_id"], g["datasource_scores"]) for g in associations["D1"]],
            [
                ("G1", {"gwas_credible_sets": 0.9}),
                ("G2", {"eva": 0.5}),
                ("G3", {"eva": 0.05}),
            ],
        )
        self.assertEqual(associations["D2"], [])
        self.assertEqual(
            sorted(cast(dict[str, object], payload["expression"])), ["G1", "G2", "G3"]
        )

    def test_reuses_same_version_expression_and_skips_fetch(self) -> None:
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with TemporaryDirectory() as directory:
            (Path(directory) / "snapshot.json").write_text(
                json.dumps(
                    {
                        "schema": 2,
                        "data_version": VERSION,
                        "expression": {
                            "G1": [{**SCHEMA2_ROW, "median": 7.0}],
                            "G2": [SCHEMA2_ROW],
                            "G3": [SCHEMA2_ROW],
                        },
                    }
                ),
                encoding="utf-8",
            )
            self._run([VERSION, VERSION], fetch, save, directory)
        fetch.assert_not_called()
        payload = cast(dict[str, object], save.call_args.args[1])
        self.assertEqual(payload["cells"], {"C1": CELL})
        expression = cast(dict[str, list[dict[str, object]]], payload["expression"])
        self.assertEqual(expression["G1"][0]["median"], 7.0)

    def test_reuse_and_fetch_mix_and_cells_match(self) -> None:
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with TemporaryDirectory() as directory:
            (Path(directory) / "snapshot.json").write_text(
                json.dumps(
                    {
                        "schema": 2,
                        "data_version": VERSION,
                        "expression": {"G1": [SCHEMA2_ROW]},
                    }
                ),
                encoding="utf-8",
            )
            self._run([VERSION, VERSION], fetch, save, directory)
        self.assertEqual([c.args[0] for c in fetch.call_args_list], ["G2", "G3"])
        payload = cast(dict[str, object], save.call_args.args[1])
        self.assertEqual(
            sorted(cast(dict[str, object], payload["expression"])), ["G1", "G2", "G3"]
        )
        self.assertEqual(payload["cells"], {"C1": CELL})

    def test_saved_cells_are_limited_to_saved_expression(self) -> None:
        row = {"cell_id": "C1", "median": 1.0, "specificity_score": None}
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with TemporaryDirectory() as directory:
            (Path(directory) / "snapshot.json").write_text(
                json.dumps(
                    {
                        "schema": 3,
                        "data_version": VERSION,
                        "cells": {"C1": CELL, "C9": {**CELL, "name": "extra"}},
                        "expression": {
                            "G1": [row],
                            "G2": [row],
                            "G3": [row],
                            "GX": [{**row, "cell_id": "C9"}],
                        },
                    }
                ),
                encoding="utf-8",
            )
            self._run([VERSION, VERSION], fetch, save, directory)
        fetch.assert_not_called()
        payload = cast(dict[str, object], save.call_args.args[1])
        self.assertEqual(payload["cells"], {"C1": CELL})
        self.assertNotIn("GX", cast(dict[str, object], payload["expression"]))

    def test_other_version_is_refetched(self) -> None:
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with TemporaryDirectory() as directory:
            (Path(directory) / "snapshot.json").write_text(
                json.dumps(
                    {
                        "schema": 3,
                        "data_version": {
                            "year": "26",
                            "month": "06",
                            "iteration": None,
                        },
                        "cells": {"C1": CELL},
                        "expression": {"G1": [], "G2": [], "G3": []},
                    }
                ),
                encoding="utf-8",
            )
            self._run([VERSION, VERSION], fetch, save, directory)
        self.assertEqual(len(fetch.call_args_list), 3)

    def test_version_change_during_refresh_raises_and_saves_nothing(self) -> None:
        later = {"year": "26", "month": "12", "iteration": None}
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with TemporaryDirectory() as directory, self.assertRaises(ValueError):
            self._run([VERSION, later], fetch, save, directory)
        save.assert_not_called()

    def test_error_before_fetch_skips_fetch_and_save(self) -> None:
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with (
            TemporaryDirectory() as directory,
            patch(
                f"{MODULE}.collect_associations",
                side_effect=ValueError("途中のページが空"),
            ),
            self.assertRaisesRegex(ValueError, "途中のページが空"),
        ):
            self._run([VERSION, VERSION], fetch, save, directory)
        fetch.assert_not_called()
        save.assert_not_called()

    def test_inconsistent_cells_across_genes_raise(self) -> None:
        def inconsistent(
            target: str,
        ) -> tuple[str, list[ExpressionRow], dict[str, CellDefinition]]:
            cell: CellDefinition = {**CELL, "name": target}
            return (
                target,
                [{"cell_id": "C1", "median": 1.0, "specificity_score": None}],
                {"C1": cell},
            )

        fetch, save = MagicMock(side_effect=inconsistent), MagicMock()
        with (
            TemporaryDirectory() as directory,
            self.assertRaisesRegex(ValueError, "一貫"),
        ):
            self._run([VERSION, VERSION], fetch, save, directory)
        save.assert_not_called()


def _graphql(bodies: list[httpx.Response]) -> httpx.MockTransport:
    """応答を順に返す。要求の本文が GraphQL の形であることも確かめる。"""

    def handle(request: httpx.Request) -> httpx.Response:
        payload = cast(dict[str, object], json.loads(request.content))
        assert set(payload) == {"query", "variables"}, payload
        return bodies.pop(0)

    return httpx.MockTransport(handle)


class QueryApiTests(unittest.TestCase):
    def test_returns_data_and_retries_transient_failures(self) -> None:
        transport = _graphql(
            [
                httpx.Response(502),
                httpx.Response(200, json={"data": {"meta": {"ok": True}}}),
            ]
        )
        with patch("autoimmune_atlas.refresh.time.sleep"):
            self.assertEqual(
                query_api("query{meta}", transport=transport), {"meta": {"ok": True}}
            )

    def test_graphql_errors_and_missing_data_raise_after_retries(self) -> None:
        for body in (
            {"errors": [{"message": "bad"}]},
            {"something": 1},
            [1, 2],
        ):
            transport = _graphql([httpx.Response(200, json=body)] * 3)
            with (
                self.subTest(body=body),
                patch("autoimmune_atlas.refresh.time.sleep"),
                self.assertRaises(RuntimeError),
            ):
                query_api("query{meta}", transport=transport)


if __name__ == "__main__":
    unittest.main()
