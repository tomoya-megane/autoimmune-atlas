"""取得と集計の臨床段階、重複、欠測、細胞分類を検証する。"""

import math
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from autoimmune_atlas.aggregation import (
    cell_catalog,
    expression_metadata,
    expression_state,
    filtered_records,
    summarize,
)
from autoimmune_atlas.refresh import (
    _extract_mechanisms,
    extract_expression,
    normalize_drugs,
    resolve_disease_ids,
    save_snapshot,
)


def _snapshot() -> dict:
    return {
        "schema": 2,
        "diseases": [
            {"id": "D1", "name": "Disease one", "status": "ready"},
            {"id": "D2", "name": "Disease two", "status": "ready"},
        ],
        "records": [
            {
                "disease_id": "D1",
                "drug_id": "A-SALT",
                "drug": "alpha salt",
                "canonical_drug_id": "A",
                "canonical_drug": "alpha",
                "drug_type": "Antibody",
                "stage": "PHASE_1",
                "target_id": "G1",
                "target": "Gene 1",
                "mechanism": "blocks",
                "references": [
                    {"source": "FDA", "ids": ["1"], "urls": ["https://example.test/1"]}
                ],
            },
            {
                "disease_id": "D1",
                "drug_id": "A",
                "drug": "alpha",
                "canonical_drug_id": "A",
                "canonical_drug": "alpha",
                "drug_type": "Antibody",
                "stage": "PHASE_3",
                "target_id": "G2",
                "target": "Gene 2",
                "mechanism": "blocks",
                "references": [],
            },
            {
                "disease_id": "D1",
                "drug_id": "B",
                "drug": "beta",
                "canonical_drug_id": "B",
                "canonical_drug": "beta",
                "drug_type": "Antibody",
                "stage": "PHASE_3",
                "target_id": "G3",
                "target": "Gene 3",
                "mechanism": "inhibits",
                "references": [],
            },
            {
                "disease_id": "D1",
                "drug_id": "U",
                "drug": "unknown",
                "canonical_drug_id": "U",
                "canonical_drug": "unknown",
                "drug_type": "Antibody",
                "stage": "PHASE_3",
                "target_id": "",
                "target": "Unknown",
                "mechanism": "",
                "references": [],
            },
            {
                "disease_id": "D2",
                "drug_id": "A",
                "drug": "alpha",
                "canonical_drug_id": "A",
                "canonical_drug": "alpha",
                "drug_type": "Antibody",
                "stage": "PHASE_1",
                "target_id": "G1",
                "target": "Gene 1",
                "mechanism": "blocks",
                "references": [],
            },
        ],
        "expression": {
            "G1": [
                {
                    "cell_id": "T4",
                    "cell": "CD4 T cell",
                    "median": 2.0,
                    "specificity_score": 0.75,
                    "parent_id": "LYMPH",
                    "parent": "Lymphocyte",
                    "ancestor_ids": ["CL_0000084"],
                },
                {
                    "cell_id": "T8",
                    "cell": "CD8 T cell",
                    "median": 0.1,
                    "specificity_score": 0.2,
                    "parent_id": "LYMPH",
                    "parent": "Lymphocyte",
                    "ancestor_ids": ["CL_0000084"],
                },
                {
                    "cell_id": "B1",
                    "cell": "B cell",
                    "median": 1.0,
                    "specificity_score": 0.8,
                    "parent_id": "B-GROUP",
                    "parent": "B lineage",
                    "ancestor_ids": [],
                },
                {
                    "cell_id": "X",
                    "cell": "Novel cell",
                    "median": 0.1,
                    "specificity_score": 0.1,
                    "parent_id": None,
                    "parent": None,
                    "ancestor_ids": [],
                },
            ],
            "G2": [
                {
                    "cell_id": "T4",
                    "cell": "CD4 T cell",
                    "median": 0.2,
                    "specificity_score": 0.1,
                    "parent_id": "LYMPH",
                    "parent": "Lymphocyte",
                    "ancestor_ids": ["CL_0000084"],
                },
                {
                    "cell_id": "T8",
                    "cell": "CD8 T cell",
                    "median": 3.0,
                    "specificity_score": 0.9,
                    "parent_id": "LYMPH",
                    "parent": "Lymphocyte",
                    "ancestor_ids": ["CL_0000084"],
                },
                {
                    "cell_id": "B1",
                    "cell": "B cell",
                    "median": 0.1,
                    "specificity_score": 0.1,
                    "parent_id": "B-GROUP",
                    "parent": "B lineage",
                    "ancestor_ids": [],
                },
                {
                    "cell_id": "X",
                    "cell": "Novel cell",
                    "median": 0.1,
                    "specificity_score": 0.1,
                    "parent_id": None,
                    "parent": None,
                    "ancestor_ids": [],
                },
            ],
            "G3": [
                {
                    "cell_id": "T4",
                    "cell": "CD4 T cell",
                    "median": None,
                    "specificity_score": None,
                    "parent_id": "LYMPH",
                    "parent": "Lymphocyte",
                    "ancestor_ids": ["CL_0000084"],
                },
                {
                    "cell_id": "T8",
                    "cell": "CD8 T cell",
                    "median": 0.1,
                    "specificity_score": 0.9,
                    "parent_id": "LYMPH",
                    "parent": "Lymphocyte",
                    "ancestor_ids": ["CL_0000084"],
                },
                {
                    "cell_id": "X",
                    "cell": "Novel cell",
                    "median": 2.0,
                    "specificity_score": 0.9,
                    "parent_id": None,
                    "parent": None,
                    "ancestor_ids": [],
                },
            ],
        },
    }


class AtlasTests(unittest.TestCase):
    def test_all_source_cells_count_each_target_and_drug_once(self) -> None:
        row = summarize(_snapshot(), "all", 0.5, level="all", disease_ids=["D1"])[0]
        self.assertEqual(row["cell_id"], "all")
        self.assertEqual(set(row["member_cell_ids"]), {"T4", "T8", "B1", "X"})
        self.assertEqual((row["count"], row["denominator"], row["unknown"]), (3, 3, 0))
        self.assertEqual((row["drug_count"], row["drug_denominator"]), (2, 2))

    def test_mixed_cells_preserve_group_and_source_counts_with_distinct_ids(
        self,
    ) -> None:
        snapshot = _snapshot()
        catalog = cell_catalog(snapshot, "mixed")
        ids = [cell["id"] for cell in catalog]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertIn("group:X", ids)
        self.assertIn("X", ids)
        selected = ["group:CL_0000084", "T4", "X", "group:X"]
        rows = summarize(
            snapshot, "all", 0.5, level="mixed", cell_ids=selected, disease_ids=["D1"]
        )
        self.assertEqual({row["cell_id"] for row in rows}, set(selected))
        for row in rows:
            original = next(
                item
                for item in summarize(snapshot, "all", 0.5, level=row["cell_level"])
                if item["disease_id"] == "D1" and item["cell_id"] == row["ontology_id"]
            )
            for key in (
                "count",
                "percent",
                "drug_count",
                "drug_percent",
                "denominator",
                "unknown",
                "member_cell_ids",
            ):
                self.assertEqual(row[key], original[key], key)
        self.assertEqual(
            summarize(snapshot, "all", 0.5, level="mixed", cell_ids=[]), []
        )

    def test_stage_is_maximized_within_disease_and_canonical_drug(self) -> None:
        rows = filtered_records(_snapshot(), "antibody", "phase3")
        self.assertEqual(
            {(r["disease_id"], r["drug_id"], r["target_id"]) for r in rows},
            {
                ("D1", "A-SALT", "G1"),
                ("D1", "A", "G2"),
                ("D1", "B", "G3"),
                ("D1", "U", ""),
            },
        )
        self.assertEqual({r["canonical_stage"] for r in rows}, {"PHASE_3"})

    def test_canonical_stage_is_chosen_before_modality_filter(self) -> None:
        snapshot = _snapshot()
        snapshot["records"] = snapshot["records"][:2]
        snapshot["records"][1]["drug_type"] = "Protein"
        rows = filtered_records(snapshot, "antibody", "phase3")
        self.assertEqual(
            [(row["drug_id"], row["canonical_stage"]) for row in rows],
            [("A-SALT", "PHASE_3")],
        )

    def test_canonical_drug_and_multiple_targets_are_counted_once(self) -> None:
        row = next(
            r
            for r in summarize(_snapshot(), "all", 0.5, level="group")
            if r["disease_id"] == "D1" and r["cell_id"] == "CL_0000084"
        )
        self.assertEqual((row["count"], row["denominator"]), (2, 3))
        self.assertEqual((row["drug_count"], row["drug_denominator"]), (1, 2))
        self.assertEqual(
            (row["total_drugs"], row["mapped_drugs"], row["unmapped_drugs"]), (3, 2, 1)
        )
        self.assertEqual({r["drug_id"] for r in row["records"]}, {"A", "A-SALT"})

    def test_target_and_drug_unknown_states_are_independent(self) -> None:
        row = next(
            r
            for r in summarize(_snapshot(), "all", 0.5, level="cell")
            if r["disease_id"] == "D1" and r["cell_id"] == "T4"
        )
        self.assertEqual((row["count"], row["unknown"]), (1, 1))
        self.assertEqual((row["drug_count"], row["unknown_drugs"]), (1, 1))
        self.assertAlmostEqual(row["percent"], 100 / 3)
        self.assertEqual(row["drug_percent"], 50)

    def test_unmapped_drug_does_not_change_known_percentages(self) -> None:
        snapshot = _snapshot()
        snapshot["records"] = [snapshot["records"][1], snapshot["records"][3]]
        row = next(
            r
            for r in summarize(snapshot, "all", 3.0, level="cell")
            if r["cell_id"] == "T4" and r["disease_id"] == "D1"
        )
        self.assertEqual(row["count"], 0)
        self.assertEqual(row["drug_count"], 0)
        self.assertEqual(row["percent"], 0)
        self.assertEqual(row["drug_percent"], 0)

    def test_unknown_expression_keeps_count_zero_and_percent_unknown(self) -> None:
        snapshot = _snapshot()
        snapshot["records"] = [snapshot["records"][2]]
        row = next(
            r
            for r in summarize(snapshot, "all", 0.5, level="cell")
            if r["cell_id"] == "T4" and r["disease_id"] == "D1"
        )
        self.assertEqual((row["count"], row["drug_count"]), (0, 0))
        self.assertEqual((row["unknown"], row["unknown_drugs"]), (1, 1))
        self.assertIsNone(row["percent"])
        self.assertIsNone(row["drug_percent"])

    def test_relative_median_uses_all_reference_cells_and_missing_is_unknown(
        self,
    ) -> None:
        snapshot = _snapshot()
        b_cell = next(
            r
            for r in summarize(snapshot, "all", 0.5, method="relative", level="cell")
            if r["disease_id"] == "D1" and r["cell_id"] == "B1"
        )
        self.assertEqual(b_cell["count"], 1)
        t4 = next(
            r
            for r in summarize(snapshot, "all", 0.5, method="relative", level="cell")
            if r["disease_id"] == "D1" and r["cell_id"] == "T4"
        )
        self.assertEqual(t4["unknown"], 1)
        self.assertEqual(
            expression_metadata(snapshot)["G1", "B1"]["target_median"], 0.55
        )

    def test_stage_filter_boundaries_follow_mixed_phase_policy(self) -> None:
        snapshot = _snapshot()
        stages = (
            "PHASE_1",
            "PHASE_1_2",
            "PHASE_2",
            "PHASE_2_3",
            "PHASE_3",
            "PREAPPROVAL",
            "APPROVAL",
            "PHASE_4",
        )
        snapshot["records"] = [
            {
                **snapshot["records"][1],
                "drug_id": value,
                "drug": value,
                "canonical_drug_id": value,
                "canonical_drug": value,
                "stage": value,
            }
            for value in stages
        ]
        expected = {
            "phase1": set(stages),
            "phase2": set(stages[2:]),
            "phase3": set(stages[4:]),
            "approved": {"APPROVAL", "PHASE_4"},
        }
        for stage, wanted in expected.items():
            with self.subTest(stage=stage):
                self.assertEqual(
                    {row["drug_id"] for row in filtered_records(snapshot, stage=stage)},
                    wanted,
                )

    def test_cpm_and_specificity_boundaries_are_inclusive(self) -> None:
        row = next(
            r
            for r in summarize(
                _snapshot(),
                "all",
                0.5,
                method="specificity",
                specificity_threshold=0.75,
                level="cell",
            )
            if r["disease_id"] == "D1" and r["cell_id"] == "T4"
        )
        self.assertEqual(row["count"], 1)
        self.assertEqual(row["records"][0]["specificity_score"], 0.75)
        boundary = _snapshot()
        boundary["expression"]["G1"][0]["median"] = 0.5
        row = next(
            r
            for r in summarize(boundary, "all", 0.5, method="specificity", level="cell")
            if r["disease_id"] == "D1" and r["cell_id"] == "T4"
        )
        self.assertIn("G1", {record["target_id"] for record in row["records"]})
        at_threshold = {"median": 0.5, "target_median": 0.5, "specificity_score": 0.75}
        for method in ("fixed", "relative", "specificity"):
            self.assertTrue(expression_state(at_threshold, 0.5, method, 0.75))
        self.assertFalse(expression_state({"median": 0.49}, 0.5, "fixed"))
        metadata = expression_metadata(_snapshot())["G1", "T4"]
        self.assertTrue(expression_state(metadata, 0.5, "specificity", 0.75))
        self.assertTrue(
            expression_state(
                {"median": 1, "specificity_score": 0.5}, 0.5, "specificity"
            )
        )
        with self.assertRaises(ValueError):
            expression_state(metadata, 0.5, "guess")

    def test_group_catalog_unions_t_cells_and_keeps_unknown_parent(self) -> None:
        groups = cell_catalog(_snapshot(), "group")
        self.assertIn(
            {"id": "CL_0000084", "name": "T cell", "members": ["T4", "T8"]}, groups
        )
        self.assertIn({"id": "X", "name": "Novel cell", "members": ["X"]}, groups)
        self.assertIn(
            {"id": "T4", "name": "CD4 T cell", "members": ["T4"]},
            cell_catalog(_snapshot(), "cell"),
        )

    def test_catalog_uses_lineage_group_priority_and_keeps_every_cell(self) -> None:
        rows = [
            {
                "cell_id": "CD8_EFFECTOR",
                "cell": "A CD8 effector",
                "median": 1,
                "parent_id": "OTHER",
                "parent": "Other",
                "ancestor_ids": ["CL_0000084", "CD8"],
            },
            {
                "cell_id": "UNKNOWN_Z",
                "cell": "Zeta cell",
                "median": 1,
                "parent_id": "UNKNOWN_Z",
                "parent": "Zeta group",
                "ancestor_ids": [],
            },
            {
                "cell_id": "CD4_MEMORY",
                "cell": "A CD4 memory",
                "median": 1,
                "parent_id": "OTHER",
                "parent": "Other",
                "ancestor_ids": ["CL_0000084", "CD4"],
            },
            {
                "cell_id": "B",
                "cell": "B cell",
                "median": 1,
                "parent_id": "CL_0000945",
                "parent": "B lineage",
                "ancestor_ids": ["CL_0000945"],
            },
            {
                "cell_id": "CD8",
                "cell": "CD8 T cell",
                "median": 1,
                "parent_id": "OTHER",
                "parent": "Other",
                "ancestor_ids": ["CL_0000084"],
            },
            {
                "cell_id": "UNKNOWN_A",
                "cell": "Alpha cell",
                "median": 1,
                "parent_id": "UNKNOWN_A",
                "parent": "Alpha group",
                "ancestor_ids": [],
            },
            {
                "cell_id": "CD4",
                "cell": "CD4 T cell",
                "median": 1,
                "parent_id": "OTHER",
                "parent": "Other",
                "ancestor_ids": ["CL_0000084"],
            },
        ]
        snapshot = {"schema": 2, "expression": {"G": rows}}

        groups = cell_catalog(snapshot, "group")
        cells = cell_catalog(snapshot, "cell")
        self.assertEqual(
            [group["id"] for group in groups],
            ["CL_0000084", "CL_0000945", "UNKNOWN_A", "UNKNOWN_Z"],
        )
        self.assertEqual(
            groups[0]["members"], ["CD4", "CD4_MEMORY", "CD8", "CD8_EFFECTOR"]
        )
        self.assertEqual(
            [cell["id"] for cell in cells],
            ["CD4", "CD4_MEMORY", "CD8", "CD8_EFFECTOR", "B", "UNKNOWN_A", "UNKNOWN_Z"],
        )
        self.assertEqual(
            {cell["id"] for cell in cells}, {row["cell_id"] for row in rows}
        )

        snapshot["expression"]["G"] = list(reversed(rows))
        self.assertEqual(cell_catalog(snapshot, "group"), groups)
        self.assertEqual(cell_catalog(snapshot, "cell"), cells)

    def test_group_membership_must_be_consistent_across_targets(self) -> None:
        snapshot = _snapshot()
        snapshot["expression"]["G2"][0]["ancestor_ids"] = []
        with self.assertRaisesRegex(ValueError, "一貫"):
            cell_catalog(snapshot)

    def test_invalid_modes_and_nonfinite_thresholds_fail(self) -> None:
        for kwargs in (
            {"stage": "late"},
            {"method": "guess"},
            {"level": "organ"},
            {"specificity_threshold": math.inf},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                summarize(snapshot=_snapshot(), modality="all", threshold=0.5, **kwargs)
        with self.assertRaises(ValueError):
            summarize(_snapshot(), "all", math.nan)
        with self.assertRaises(ValueError):
            summarize(_snapshot(), "all", -0.1)
        with self.assertRaises(ValueError):
            summarize(_snapshot(), "all", 0.5, specificity_threshold=1.1)
        with self.assertRaises(ValueError):
            filtered_records(_snapshot(), "other", "phase3")


class ImportTests(unittest.TestCase):
    def test_scope_is_the_union_of_roots_and_body_only_roots_add_no_descendants(
        self,
    ) -> None:
        tree = {
            "MONDO_0007179": ["A", "B"],
            "MONDO_0003346": ["AIDS_ARTERITIS", "B"],
        }

        def fake_query(query: str, variables: dict) -> dict:
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
        with patch("autoimmune_atlas.refresh.SCOPE_ROOTS", (("MISSING", True),)):
            with self.assertRaises(ValueError):
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
        rows = [
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
        rows = [
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
        row = {
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

    def test_expression_keeps_specificity_parent_and_ancestors(self) -> None:
        base = {
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
        self.assertEqual(
            extract_expression([base]),
            [
                {
                    "cell_id": "C1",
                    "cell": "B cell",
                    "median": 2.0,
                    "specificity_score": 0.8,
                    "parent_id": "P1",
                    "parent": "Lymphocyte",
                    "ancestor_ids": ["CL_1", "CL_2"],
                }
            ],
        )
        self.assertEqual(
            extract_expression(
                [
                    {**base, "datasourceId": "DICE"},
                    {**base, "tissueBiosample": {"biosampleId": "T1"}},
                ]
            ),
            [],
        )
        for key, value in (("median", -1), ("specificity_score", 1.1)):
            with self.subTest(key=key), self.assertRaises(ValueError):
                extract_expression([{**base, key: value}])
        with self.assertRaisesRegex(ValueError, "単位"):
            extract_expression([{**base, "unit": "TPM"}])

    def test_mechanism_references_are_deduplicated(self) -> None:
        row = {
            "mechanismOfAction": "blocks",
            "actionType": "BLOCKER",
            "targets": [{"id": "G1", "approvedSymbol": "GENE1"}],
            "references": [
                {"source": "FDA", "ids": ["1"], "urls": ["https://example.test"]},
                {"source": "FDA", "ids": ["1"], "urls": ["https://example.test"]},
            ],
        }
        self.assertEqual(
            _extract_mechanisms([row]),
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


if __name__ == "__main__":
    unittest.main()
