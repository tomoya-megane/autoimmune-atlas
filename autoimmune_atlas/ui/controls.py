"""両ページの callback が共有する、入力値の検証と行の切り替え。"""

import math
from collections.abc import Callable
from typing import TypedDict, cast

from dash import ctx, no_update

from autoimmune_atlas.ui.config import (
    DEFAULT_EXPRESSION_THRESHOLD,
    DEFAULT_SPECIFICITY_THRESHOLD,
)

type NumberInput = int | float | str | None
type ParameterValue = str | int | float | list[str] | None


class CellToggleId(TypedDict):
    type: str
    kind: str
    cell: str


def effective_number(
    value: NumberInput,
    default: float,
    label: str,
    maximum: float | None = None,
    *,
    minimum: float | None = None,
) -> tuple[float, str | None]:
    """空欄には初期値を使い、不正値は理由を示して初期値へ戻す。"""
    lower = minimum if minimum is not None else 0
    if value is None or value == "":
        return default, None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default, f"Invalid {label}; using {default:g}."
    if (
        isinstance(value, bool)
        or not math.isfinite(number)
        or number < lower
        or (maximum is not None and number > maximum)
    ):
        bound = (
            f" between {lower:g} and {maximum:g}"
            if maximum is not None
            else f" at or above {lower:g}"
        )
        return (
            default,
            f"{label[0].upper() + label[1:]} must be finite and{bound}; using {default:g}.",
        )
    return number, None


def effective_filters(
    threshold: NumberInput, specificity: NumberInput
) -> tuple[float, float, list[str]]:
    """実際に集計へ渡す閾値と入力エラーを返す。"""
    minimum, minimum_error = effective_number(
        threshold, DEFAULT_EXPRESSION_THRESHOLD, "minimum CPM"
    )
    specificity_value, specificity_error = effective_number(
        specificity, DEFAULT_SPECIFICITY_THRESHOLD, "CELLEX specificity", 1
    )
    return (
        minimum,
        specificity_value,
        [message for message in (minimum_error, specificity_error) if message],
    )


def make_toggle(
    valid_groups: set[str],
) -> Callable[[list[int | None], list[str] | None], object]:
    """行見出しのクリックで、大分類の展開を切り替える callback の本体を作る。"""

    def toggle_cell_group(
        _clicks: list[int | None], expanded: list[str] | None
    ) -> object:
        # Newly rendered buttons have zero clicks; only user clicks toggle a group.
        triggered = cast(CellToggleId | str | None, ctx.triggered_id)
        inputs_list = cast(list[list[dict[str, object]]], ctx.inputs_list)
        if not isinstance(triggered, dict) or not any(
            item["id"] == triggered and item.get("value") for item in inputs_list[0]
        ):
            return no_update
        cell_id = triggered["cell"]
        if cell_id not in valid_groups:
            return no_update
        expanded_set = set(expanded or ()) & valid_groups
        expanded_set.symmetric_difference_update([cell_id])
        return sorted(expanded_set)

    return toggle_cell_group
