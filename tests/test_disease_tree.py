"""疾患の選択ツリーの節分けと配置を検証する。"""

from __future__ import annotations

import unittest
from typing import cast

from autoimmune_atlas.models import (
    DiseaseFamily,
)
from autoimmune_atlas.ui import (
    components,
)
from tests.ui_fixture import (
    _component,
    fixture,
)


class DiseaseTreeTests(unittest.TestCase):
    def test_sections_follow_parents_and_place_each_term_once(self) -> None:
        family: DiseaseFamily = {
            "id": "ROOT",
            "label": "root disease",
            "diseases": [
                {"id": "ROOT", "name": "root disease", "parent_ids": []},
                {"id": "CHILD", "name": "child", "parent_ids": ["ROOT"]},
                {"id": "LEAF", "name": "leaf", "parent_ids": ["ROOT"]},
                {"id": "GRAND", "name": "grandchild", "parent_ids": ["CHILD"]},
                {
                    "id": "SHARED",
                    "name": "shared term",
                    "parent_ids": ["ROOT", "GRAND"],
                },
                {"id": "LOOP_A", "name": "loop a", "parent_ids": ["LOOP_B"]},
                {"id": "LOOP_B", "name": "loop b", "parent_ids": ["LOOP_A"]},
            ],
        }
        sections = components.disease_checklist_sections(family)
        ids = {s["id"]: [d["id"] for d in s["diseases"]] for s in sections}
        self.assertEqual(
            ids,
            {
                "disease-family-ROOT": ["ROOT"],
                "disease-details-ROOT": ["LEAF"],
                "disease-family-ROOT-CHILD": ["CHILD"],
                "disease-family-ROOT-GRAND": ["GRAND"],
                "disease-details-ROOT-GRAND": ["SHARED"],
                "disease-family-ROOT-LOOP_A": ["LOOP_A"],
                "disease-details-ROOT-LOOP_A": ["LOOP_B"],
            },
        )
        placed = [d for s in sections for d in s["diseases"]]
        self.assertEqual(len(placed), len(family["diseases"]))
        selector_snapshot = fixture()
        selector_snapshot["diseases"] = [
            {**d, "status": "ready"} for d in family["diseases"]
        ]
        layout = components.disease_selector(selector_snapshot, [])
        summaries: list[object] = []

        def walk(component: object) -> None:
            if not hasattr(component, "children"):
                return
            if type(component).__name__ == "Summary":
                summaries.append(_component(component).children)
            children = _component(component).children
            if children is None:
                return
            if not isinstance(children, list):
                children = [children]
            for child in cast(list[object], children):
                if not isinstance(child, str):
                    walk(child)

        walk(layout)
        self.assertIn("child (3 terms)", summaries)
        self.assertIn("child details (2 terms)", summaries)
        self.assertIn("grandchild (2 terms)", summaries)

    def test_family_without_root_keeps_flat_top_level(self) -> None:
        family: DiseaseFamily = {
            "id": "other",
            "label": "Other terms",
            "diseases": [
                {"id": "B", "name": "b term", "parent_ids": []},
                {"id": "A", "name": "a term", "parent_ids": ["MISSING"]},
            ],
        }
        sections = components.disease_checklist_sections(family)
        self.assertEqual(
            [(s["id"], [d["id"] for d in s["diseases"]]) for s in sections],
            [("disease-family-other", ["A", "B"])],
        )
