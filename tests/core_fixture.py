"""集計と取得のテストが共有する最小のスナップショット。"""

from autoimmune_atlas.models import CoreSnapshot


def core_snapshot() -> CoreSnapshot:
    return {
        "schema": 3,
        "diseases": [
            {"id": "D1", "name": "Disease one", "status": "ready"},
            {"id": "D2", "name": "Disease two", "status": "ready"},
        ],
        "datasources": ["gwas_credible_sets", "eva"],
        "associations": {
            "D1": [
                {
                    "target_id": "G1",
                    "target": "Gene 1",
                    "score": 0.9,
                    "datasource_scores": {"gwas_credible_sets": 0.9},
                },
                {
                    "target_id": "G2",
                    "target": "Gene 2",
                    "score": 0.6,
                    "datasource_scores": {"gwas_credible_sets": 0.4, "eva": 0.6},
                },
                {
                    "target_id": "G9",
                    "target": "Gene 9",
                    "score": 0.5,
                    "datasource_scores": {"eva": 0.5},
                },
                {
                    "target_id": "G3",
                    "target": "Gene 3",
                    "score": 0.2,
                    "datasource_scores": {"eva": 0.2},
                },
            ],
            "D2": [],
        },
        "records": [
            {
                "disease_id": "D1",
                "disease": "Disease one",
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
                "disease": "Disease one",
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
                "disease": "Disease one",
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
                "disease": "Disease one",
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
                "disease": "Disease two",
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
        "cells": {
            "T4": {
                "name": "CD4 T cell",
                "parent_id": "LYMPH",
                "parent": "Lymphocyte",
                "ancestor_ids": ["CL_0000084"],
            },
            "T8": {
                "name": "CD8 T cell",
                "parent_id": "LYMPH",
                "parent": "Lymphocyte",
                "ancestor_ids": ["CL_0000084"],
            },
            "B1": {
                "name": "B cell",
                "parent_id": "B-GROUP",
                "parent": "B lineage",
                "ancestor_ids": [],
            },
            "X": {
                "name": "Novel cell",
                "parent_id": None,
                "parent": None,
                "ancestor_ids": [],
            },
        },
        "expression": {
            "G1": [
                {"cell_id": "T4", "median": 2.0, "specificity_score": 0.75},
                {"cell_id": "T8", "median": 0.1, "specificity_score": 0.2},
                {"cell_id": "B1", "median": 1.0, "specificity_score": 0.8},
                {"cell_id": "X", "median": 0.1, "specificity_score": 0.1},
            ],
            "G2": [
                {"cell_id": "T4", "median": 0.2, "specificity_score": 0.1},
                {"cell_id": "T8", "median": 3.0, "specificity_score": 0.9},
                {"cell_id": "B1", "median": 0.1, "specificity_score": 0.1},
                {"cell_id": "X", "median": 0.1, "specificity_score": 0.1},
            ],
            "G3": [
                {"cell_id": "T4", "median": None, "specificity_score": None},
                {"cell_id": "T8", "median": 0.1, "specificity_score": 0.9},
                {"cell_id": "X", "median": 2.0, "specificity_score": 0.9},
            ],
            "G9": [
                {"cell_id": "T4", "median": 5.0, "specificity_score": 0.9},
                {"cell_id": "T8", "median": 0.0, "specificity_score": 0.0},
                {"cell_id": "B1", "median": 0.0, "specificity_score": 0.0},
                {"cell_id": "X", "median": 0.0, "specificity_score": 0.0},
            ],
        },
    }
