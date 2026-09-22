"""集計の重複、分母、欠測の扱いを確認する。"""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from atlas import summarize, to_csv
from fetch_data import extract_expression, normalize_drugs, save_snapshot


class AtlasTests(unittest.TestCase):
    def setUp(self) -> None:
        self.snapshot = {
            "diseases": [{"id": "D1", "name": "Disease", "status": "ready"}],
            "records": [
                {"disease_id": "D1", "drug_id": "A", "drug": "Alpha", "modality": "antibody", "target_id": "G1", "target": "Gene1"},
                {"disease_id": "D1", "drug_id": "B", "drug": "Beta", "modality": "antibody", "target_id": "G1", "target": "Gene1"},
                {"disease_id": "D1", "drug_id": "C", "drug": "Gamma", "modality": "small_molecule", "target_id": "G2", "target": "Gene2"},
            ],
            "expression": {
                "G1": [{"cell_id": "C1", "cell": "B cell", "median": 2.0}],
                "G2": [{"cell_id": "C1", "cell": "B cell", "median": 0.1}],
            },
        }

    def test_unique_targets_and_filtered_denominator(self) -> None:
        all_rows = summarize(self.snapshot, [], "expression", "all", 0.5)
        self.assertEqual(all_rows[0]["count"], 1)
        self.assertEqual(all_rows[0]["denominator"], 2)
        self.assertEqual(all_rows[0]["percent"], 50)
        antibodies = summarize(self.snapshot, [], "expression", "antibody", 0.5)
        self.assertEqual(antibodies[0]["percent"], 100)
        self.assertEqual(len(antibodies[0]["records"]), 2)

    def test_missing_expression_is_not_zero(self) -> None:
        del self.snapshot["expression"]["G1"]
        row = summarize(self.snapshot, [], "expression", "all", 0.5)[0]
        self.assertIsNone(row["count"])
        self.assertIsNone(row["percent"])
        self.assertEqual(row["unknown"], 1)

    def test_unmapped_drugs_are_reported_outside_known_target_denominator(self) -> None:
        self.snapshot["records"].append({"disease_id": "D1", "drug_id": "D", "drug": "Delta", "modality": "other", "target_id": "", "target": "未判明"})
        row = summarize(self.snapshot, [], "expression", "all", 0.5)[0]
        self.assertEqual(row["percent"], 50)
        self.assertEqual(row["denominator"], 2)
        self.assertEqual(row["unmapped_drugs"], 1)
        self.assertEqual(row["status"], "partial")

    def test_partial_positive_is_a_lower_bound(self) -> None:
        del self.snapshot["expression"]["G2"]
        row = summarize(self.snapshot, [], "expression", "all", 0.5)[0]
        self.assertEqual(row["count"], 1)
        self.assertEqual(row["status"], "partial")
        self.assertEqual(row["unknown"], 1)

    def test_mechanism_requires_matching_drug_and_disease(self) -> None:
        annotations = [{"disease_id": "D1", "drug_id": "A", "target_id": "G1", "cell_id": "C1", "cell": "B cell", "status": "yes", "source": "https://example.org/evidence", "note": "test"}]
        row = summarize(self.snapshot, annotations, "mechanism", "all", 0.5)[0]
        self.assertEqual([r["drug_id"] for r in row["records"]], ["A"])
        self.assertEqual(row["count"], 1)
        self.assertEqual(row["unknown"], 1)
        annotations[0]["disease_id"] = "D2"
        row = summarize(self.snapshot, annotations, "mechanism", "all", 0.5)[0]
        self.assertIsNone(row["count"])

    def test_fetch_error_and_empty_denominator(self) -> None:
        self.snapshot["diseases"][0]["status"] = "error"
        self.assertIsNone(summarize(self.snapshot, [], "expression", "all", 0.5)[0]["count"])
        self.snapshot["diseases"][0]["status"] = "ready"
        self.snapshot["records"] = []
        rows = summarize(self.snapshot, [{"cell_id": "C1", "cell": "B cell"}], "expression", "all", 0.5)
        self.assertEqual(rows[0]["count"], 0)
        self.assertIsNone(rows[0]["percent"])

    def test_csv_preserves_unicode_and_neutralizes_formulas(self) -> None:
        csv = to_csv([{"薬剤": "=1+1", "細胞": "B細胞", "数": None}])
        self.assertIn("B細胞", csv)
        self.assertIn("'=1+1", csv)


class ImportTests(unittest.TestCase):
    def test_disease_specific_phase_and_missing_targets(self) -> None:
        rows = [
            {"maxClinicalStage": "PHASE_2", "drug": {"id": "A", "name": "A", "drugType": "Antibody", "maximumClinicalStage": "APPROVAL"}},
            {"maxClinicalStage": "PHASE_3", "drug": {"id": "B", "name": "B", "drugType": "Protein"}},
            {"maxClinicalStage": "APPROVAL", "drug": {"id": "C", "name": "C", "drugType": "Unknown"}},
        ]
        drugs = normalize_drugs(rows)
        self.assertEqual([d["drug_id"] for d in drugs], ["B", "C"])
        self.assertEqual(drugs[1]["modality"], "unknown")

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
