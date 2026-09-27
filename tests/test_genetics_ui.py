"""遺伝子ページの図、部品、レイアウト、callback を検証する。"""

from __future__ import annotations

import unittest
from collections.abc import Sequence
from typing import Protocol, cast

from autoimmune_atlas.models import DiseaseFamily, SummaryRow
from autoimmune_atlas.ui import components, figures


class _Title(Protocol):
    text: str


class _ColorBar(Protocol):
    title: _Title


class _Heatmap(Protocol):
    colorbar: _ColorBar


class _Figure(Protocol):
    data: Sequence[_Heatmap]


class _Button(Protocol):
    id: dict[str, str]


def gene_row(**overrides: object) -> SummaryRow:
    base: dict[str, object] = {
        "disease_id": "D1",
        "disease": "Disease one",
        "cell_id": "group:LYMPH",
        "cell": "Lymphocyte (group)",
        "ontology_id": "LYMPH",
        "cell_level": "group",
        "count": 3,
        "percent": 75.0,
        "denominator": 4,
        "unknown": 0,
        "unmapped_drugs": 0,
        "status": "complete",
        "records": [],
        "drug_count": None,
        "drug_percent": None,
        "drug_denominator": 0,
        "unknown_drugs": 0,
        "total_drugs": 0,
        "mapped_drugs": 0,
        "member_cell_ids": ["T4", "T8"],
    }
    base.update(overrides)
    return cast(SummaryRow, cast(object, base))


class GeneFigureTests(unittest.TestCase):
    def test_gene_hover_names_genes_and_omits_drug_coverage(self) -> None:
        text = figures.hover_text(gene_row(), "percent", "gene")
        self.assertIn("Genes: 75.0%", text)
        self.assertNotIn("Canonical drugs", text)
        dot = figures.dot_hover_text(gene_row(), "gene")
        self.assertIn("Count: 3", dot)
        self.assertNotIn("Canonical drugs", dot)

    def test_zero_genes_above_threshold_has_its_own_message(self) -> None:
        row = gene_row(count=0, percent=None, denominator=0)
        text = figures.hover_text(row, "count", "gene")
        self.assertIn("No genes at or above the score threshold", text)
        text = figures.hover_text(row, "percent", "gene")
        self.assertIn("No genes at or above the score threshold", text)

    def test_gene_figure_title(self) -> None:
        figure = figures.build_figure(
            [gene_row()], ["D1"], ["group:LYMPH"], "count", "gene"
        )
        heatmap = cast(_Figure, cast(object, figure)).data[0]
        self.assertEqual(heatmap.colorbar.title.text, "Genes (count)")


class PrefixedComponentTests(unittest.TestCase):
    def test_checklist_sections_take_prefix(self) -> None:
        family: DiseaseFamily = {
            "id": "ROOT",
            "label": "Root",
            "diseases": [
                {"id": "ROOT", "name": "Root"},
                {"id": "CHILD", "name": "Child", "parent_ids": ["ROOT"]},
            ],
        }
        ids = [
            s["id"] for s in components.disease_checklist_sections(family, "genetics-")
        ]
        self.assertEqual(
            ids, ["genetics-disease-family-ROOT", "genetics-disease-details-ROOT"]
        )
        plain = [s["id"] for s in components.disease_checklist_sections(family)]
        self.assertEqual(plain, ["disease-family-ROOT", "disease-details-ROOT"])

    def test_row_controls_use_kind_specific_toggle_type(self) -> None:
        names = {"group:LYMPH": "Lymphocyte (group)"}
        button = components.heatmap_row_controls(
            names, ["group:LYMPH"], [], "genetics"
        )[0]
        self.assertEqual(
            cast(_Button, cast(object, button)).id["type"], "genetics-cell-toggle"
        )
        button = components.heatmap_row_controls(names, ["group:LYMPH"], [], "target")[
            0
        ]
        self.assertEqual(
            cast(_Button, cast(object, button)).id["type"], "heatmap-cell-toggle"
        )
