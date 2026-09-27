"""遺伝子ページの図、部品、レイアウト、callback を検証する。"""

from __future__ import annotations

import unittest
from collections.abc import Iterator, Sequence
from typing import Protocol, cast

from autoimmune_atlas import genetics
from autoimmune_atlas.models import DiseaseFamily, SummaryRow
from autoimmune_atlas.ui import components, figures, genetics_layout
from tests.genetics_fixture import genetics_snapshot
from tests.test_genetics import snapshot as core_ui_snapshot


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


class _NumberInput(Protocol):
    min: float
    max: float
    value: float


class _Grid(Protocol):
    columnDefs: list[dict[str, object]]


class _Classed(Protocol):
    className: str


class _Identified(Protocol):
    id: str


class _Linked(Protocol):
    href: str


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


def _walk(node: object) -> Iterator[object]:
    yield node
    children = getattr(node, "children", None)
    if isinstance(children, list):
        for child in cast(list[object], children):
            yield from _walk(child)
    elif children is not None:
        yield from _walk(cast(object, children))


class GeneticsLayoutTests(unittest.TestCase):
    def test_page_has_prefixed_controls_and_nav_marks_current(self) -> None:
        page = genetics_layout.genetics_page(core_ui_snapshot(), genetics_snapshot())
        found = {
            cast(_Identified, c).id
            for c in _walk(page)
            if isinstance(getattr(c, "id", None), str)
        }
        for component_id in (
            "genetics-page",
            "genetics-diseases",
            "genetics-score",
            "genetics-method",
            "genetics-threshold",
            "genetics-specificity",
            "genetics-measure",
            "genetics-update-button",
            "genetics-applied-parameters",
            "genetics-heatmap",
            "genetics-chart-type",
            "genetics-detail-disease",
            "genetics-details",
        ):
            self.assertIn(component_id, found)
        self.assertNotIn("diseases", found)
        self.assertNotIn("update-button", found)
        self.assertTrue({"genetics-comparison-tip", "genetics-selection-tip"} <= found)
        self.assertEqual([c for c in found if not c.startswith("genetics-")], [])
        current = [c for c in _walk(page) if getattr(c, "aria-current", None) == "page"]
        self.assertEqual([cast(_Linked, c).href for c in current], ["/genetics"])

    def test_score_input_range_and_default(self) -> None:
        page = genetics_layout.genetics_page(core_ui_snapshot(), genetics_snapshot())
        score = cast(
            _NumberInput,
            next(c for c in _walk(page) if getattr(c, "id", None) == "genetics-score"),
        )
        self.assertEqual(score.min, 0.1)
        self.assertEqual(score.max, 1)
        self.assertEqual(score.value, 0.5)

    def test_gene_rows_have_scores_datasource_columns_and_evidence_link(self) -> None:
        data = genetics_snapshot()
        genes = genetics.genes_for_disease(data, "D1", 0.5)
        rows = genetics_layout.gene_rows(genes, "D1", data["datasources"])
        self.assertEqual(rows[0]["gene"], "Gene 1 (G1)")
        self.assertEqual(rows[0]["score"], "0.900")
        self.assertEqual(rows[0]["gwas_credible_sets"], "0.900")
        self.assertEqual(rows[0]["eva"], "—")
        self.assertIn(
            "https://platform.opentargets.org/evidence/G1/D1", rows[0]["links"]
        )
        grid = cast(
            _Grid, cast(object, genetics_layout.gene_grid(rows, data["datasources"]))
        )
        fields = [c["field"] for c in grid.columnDefs]
        self.assertEqual(
            fields, ["gene", "score", "gwas_credible_sets", "eva", "links"]
        )

    def test_detail_panel_without_selection_and_unavailable_page(self) -> None:
        panel = genetics_layout.genetics_detail_panel(
            [],
            None,
            genetics_snapshot(),
            score_threshold=0.5,
            threshold=0.5,
            method="fixed",
            specificity=0.5,
        )
        self.assertEqual(cast(_Classed, cast(object, panel)).className, "empty-note")
        page = genetics_layout.genetics_unavailable_page("versions differ")
        self.assertIn("refresh-genetics", str(page))
