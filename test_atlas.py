"""集計の重複、分母、欠測の扱いを確認する。"""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from atlas import DRUG_TYPE_MODALITIES, summarize, to_csv
from fetch_data import extract_expression, normalize_drugs, save_snapshot


class AtlasTests(unittest.TestCase):
    def setUp(self) -> None:
        self.snapshot = {
            "diseases": [{"id": "D1", "name": "Disease", "status": "ready"}],
            "records": [
                {"disease_id": "D1", "drug_id": "A", "drug": "Alpha", "modality": "antibody", "drug_type": "Antibody", "target_id": "G1", "target": "Gene1"},
                {"disease_id": "D1", "drug_id": "B", "drug": "Beta", "modality": "antibody", "drug_type": "Antibody", "target_id": "G1", "target": "Gene1"},
                {"disease_id": "D1", "drug_id": "C", "drug": "Gamma", "modality": "small_molecule", "drug_type": "Small molecule", "target_id": "G2", "target": "Gene2"},
            ],
            "expression": {
                "G1": [{"cell_id": "C1", "cell": "B cell", "median": 2.0}],
                "G2": [{"cell_id": "C1", "cell": "B cell", "median": 0.1}],
            },
        }

    def test_unique_targets_and_filtered_denominator(self) -> None:
        all_rows = summarize(self.snapshot, "all", 0.5)
        self.assertEqual(all_rows[0]["count"], 1)
        self.assertEqual(all_rows[0]["denominator"], 2)
        self.assertEqual(all_rows[0]["percent"], 50)
        antibodies = summarize(self.snapshot, "antibody", 0.5)
        self.assertEqual(antibodies[0]["percent"], 100)
        self.assertEqual(len(antibodies[0]["records"]), 2)

    def test_missing_expression_is_not_zero(self) -> None:
        del self.snapshot["expression"]["G1"]
        row = summarize(self.snapshot, "all", 0.5)[0]
        self.assertIsNone(row["count"])
        self.assertIsNone(row["percent"])
        self.assertEqual(row["unknown"], 1)

    def test_unmapped_drugs_are_reported_outside_known_target_denominator(self) -> None:
        self.snapshot["records"].append({"disease_id": "D1", "drug_id": "D", "drug": "Delta", "modality": "other", "drug_type": "Cell", "target_id": "", "target": "未判明"})
        row = summarize(self.snapshot, "all", 0.5)[0]
        self.assertEqual(row["percent"], 50)
        self.assertEqual(row["denominator"], 2)
        self.assertEqual(row["unmapped_drugs"], 1)
        self.assertEqual(row["status"], "partial")

    def test_original_drug_type_replaces_legacy_coarse_modality(self) -> None:
        self.snapshot["records"].extend([
            {"disease_id": "D1", "drug_id": "D", "drug": "Delta", "modality": "other", "drug_type": "Protein", "target_id": "G3", "target": "Gene3"},
            {"disease_id": "D1", "drug_id": "E", "drug": "Epsilon", "modality": "unknown", "drug_type": "Unknown", "target_id": "G4", "target": "Gene4"},
            {"disease_id": "D1", "drug_id": "F", "drug": "Zeta", "modality": "other", "drug_type": "Cell", "target_id": "", "target": "未判明"},
        ])
        self.snapshot["expression"].update({
            "G3": [{"cell_id": "C1", "cell": "B cell", "median": 2.0}],
            "G4": [{"cell_id": "C1", "cell": "B cell", "median": 2.0}],
        })

        protein = summarize(self.snapshot, "protein", 0.5)[0]
        unknown = summarize(self.snapshot, "unknown", 0.5)[0]
        cell = summarize(self.snapshot, "cell", 0.5)[0]

        self.assertEqual((protein["denominator"], protein["records"][0]["modality"]), (1, "protein"))
        self.assertEqual((unknown["denominator"], unknown["records"][0]["modality"]), (1, "unknown"))
        self.assertEqual((cell["count"], cell["denominator"], cell["unmapped_drugs"]), (None, 0, 1))

    def test_partial_positive_is_a_lower_bound(self) -> None:
        del self.snapshot["expression"]["G2"]
        row = summarize(self.snapshot, "all", 0.5)[0]
        self.assertEqual(row["count"], 1)
        self.assertEqual(row["status"], "partial")
        self.assertEqual(row["unknown"], 1)

    def test_expression_threshold_is_strict(self) -> None:
        row = summarize(self.snapshot, "all", 2.0)[0]
        self.assertEqual(row["count"], 0)
        self.assertEqual(row["unknown"], 0)

    def test_fetch_error_and_empty_denominator(self) -> None:
        self.snapshot["diseases"][0]["status"] = "error"
        self.assertIsNone(summarize(self.snapshot, "all", 0.5)[0]["count"])
        self.snapshot["diseases"][0]["status"] = "ready"
        self.snapshot["records"] = []
        rows = summarize(self.snapshot, "all", 0.5)
        self.assertEqual(rows[0]["count"], 0)
        self.assertIsNone(rows[0]["percent"])

    def test_csv_preserves_unicode_and_neutralizes_formulas(self) -> None:
        csv = to_csv([{"薬剤": "=1+1", "細胞": "B細胞", "数": None}])
        self.assertIn("B細胞", csv)
        self.assertIn("'=1+1", csv)


class ImportTests(unittest.TestCase):
    def test_disease_specific_phase_and_all_drug_types(self) -> None:
        kinds = ("Small molecule", "Antibody", "Protein", "Cell", "Gene", "Enzyme", "Oligonucleotide", "Unknown")
        rows = [{"maxClinicalStage": "PHASE_2", "drug": {"id": "A", "name": "A", "drugType": "Antibody", "maximumClinicalStage": "APPROVAL"}}]
        rows.extend({"maxClinicalStage": "PHASE_3", "drug": {"id": str(index), "name": kind, "drugType": kind}} for index, kind in enumerate(kinds))
        drugs = normalize_drugs(rows)
        self.assertEqual([drug["drug_type"] for drug in drugs], list(kinds))
        self.assertEqual([drug["modality"] for drug in drugs], [value for _, value in DRUG_TYPE_MODALITIES])

    def test_cell_only_expression_excludes_tissue_and_other_sources(self) -> None:
        base = {"datasourceId": "tabula_sapiens", "unit": "CPM(pseudobulk sum[counts])", "median": 2.0, "celltypeBiosample": {"biosampleId": "C1", "biosampleName": "B cell"}, "tissueBiosample": None}
        rows = [base, {**base, "tissueBiosample": {"biosampleId": "T1"}}, {**base, "datasourceId": "DICE"}]
        self.assertEqual(len(extract_expression(rows)), 1)
        self.assertEqual(extract_expression(rows)[0]["median"], 2)
        with self.assertRaises(ValueError):
            extract_expression([{**base, "median": -1}])

    def test_failed_serialization_preserves_previous_snapshot(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            path.write_text('{"previous": true}')
            with self.assertRaises(ValueError):
                save_snapshot(path, {"value": float("nan")})
            self.assertEqual(path.read_text(), '{"previous": true}')


if __name__ == "__main__":
    unittest.main()
