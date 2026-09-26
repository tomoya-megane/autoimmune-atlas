"""公開 API からの取得、正規化、保存を検証する。"""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from autoimmune_atlas.refresh import (
    ClinicalCandidate,
    ExpressionApiRow,
    MechanismApiRow,
    extract_expression,
    extract_mechanisms,
    normalize_drugs,
    resolve_disease_ids,
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

    def test_expression_keeps_specificity_parent_and_ancestors(self) -> None:
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


if __name__ == "__main__":
    unittest.main()
