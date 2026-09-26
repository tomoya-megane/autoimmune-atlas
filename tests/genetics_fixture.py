"""遺伝子ページのテストが共有する最小の genetics スナップショット。"""

from autoimmune_atlas.models import GeneticsSnapshot


def genetics_snapshot() -> GeneticsSnapshot:
    return {
        "schema": 1,
        "data_version": {"year": "26", "month": "09", "iteration": None},
        "retrieved_at": "2026-09-27T00:00:00+00:00",
        "source": "https://api.platform.opentargets.org/api/v4/graphql",
        "score_floor": 0.1,
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
        "expression": {
            "G9": [
                {
                    "cell_id": "T4",
                    "cell": "CD4 T cell",
                    "median": 5.0,
                    "specificity_score": 0.9,
                    "parent_id": "LYMPH",
                    "parent": "Lymphocyte",
                    "ancestor_ids": ["CL_0000084"],
                },
                {
                    "cell_id": "T8",
                    "cell": "CD8 T cell",
                    "median": 0.0,
                    "specificity_score": 0.0,
                    "parent_id": "LYMPH",
                    "parent": "Lymphocyte",
                    "ancestor_ids": ["CL_0000084"],
                },
                {
                    "cell_id": "B1",
                    "cell": "B cell",
                    "median": 0.0,
                    "specificity_score": 0.0,
                    "parent_id": "B-GROUP",
                    "parent": "B lineage",
                    "ancestor_ids": [],
                },
                {
                    "cell_id": "X",
                    "cell": "Novel cell",
                    "median": 0.0,
                    "specificity_score": 0.0,
                    "parent_id": None,
                    "parent": None,
                    "ancestor_ids": [],
                },
            ]
        },
    }
