"""閾値の入力と初期選択の既定値を検証する。"""

from __future__ import annotations

import unittest
from typing import ClassVar

from autoimmune_atlas.ui import (
    callbacks,
    layout,
)


class InputTests(unittest.TestCase):
    """空欄と不正値で表示値と計算値がずれない。"""

    def test_empty_and_invalid_thresholds_use_visible_defaults(self) -> None:
        self.assertEqual(callbacks.effective_filters(None, ""), (0.5, 0.5, []))
        minimum, specificity, errors = callbacks.effective_filters("NaN", 1.5)
        self.assertEqual((minimum, specificity), (0.5, 0.5))
        self.assertEqual(len(errors), 2)


class DefaultsTests(unittest.TestCase):
    """初期選択は、完全一致した語で部分一致の用語を足さない。"""

    ITEMS: ClassVar[list[tuple[str, str]]] = [
        ("A", "systemic lupus erythematosus"),
        ("B", "autosomal systemic lupus erythematosus type 16"),
        ("C", "relapsing-remitting multiple sclerosis"),
    ]

    def test_exact_match_skips_substring_match(self) -> None:
        self.assertEqual(
            layout.choose_defaults(
                self.ITEMS, ["systemic lupus erythematosus"], 2, {"A", "C"}
            ),
            ["A", "C"],
        )

    def test_term_without_exact_match_falls_back_to_substring(self) -> None:
        self.assertEqual(
            layout.choose_defaults(
                self.ITEMS,
                ["systemic lupus erythematosus", "multiple sclerosis"],
                2,
            ),
            ["A", "C"],
        )
