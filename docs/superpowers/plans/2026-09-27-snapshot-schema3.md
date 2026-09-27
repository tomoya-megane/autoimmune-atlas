# 保存データの統合（schema 3）の実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 薬剤の記録、genetic association の関連遺伝子、細胞型別の発現を 1 つの `data/snapshot.json`（schema 3）にまとめ、`pixi run refresh` 1 回で取得する。

**Architecture:** `models.py` の `ExpressionRow` を `cell_id`、`median`、`specificity_score` の 3 キーに減らし、細胞の定数（名前、親、祖先）は `snapshot["cells"]` の表に 1 回だけ置く。`aggregation.py` は `cells` から大分類を作り、`expression_metadata()` が実行時の行に `cell` 名を付ける。`genetics.py` は `snapshot["associations"]` を読み、`genetics.json` と `refresh_genetics.py` は消す。`refresh.py` の `main()` が関連遺伝子の取得も行い、同じ版の発現は既存のファイルから再利用する。保存データと API の応答の検証は msgspec に任せ、手書きの `isinstance` と二重の `cast` を消す。API の呼び出しは httpx にし、proxy の手書きの処理を消す。

**Tech Stack:** Python 3.13、Dash 4、Plotly 7、dash-ag-grid、msgspec 0.21、httpx 0.28、unittest、Ruff、basedpyright。msgspec と httpx は Tomoya が `pixi add` 済み。それ以外の依存パッケージは足さない。

**Spec:** `docs/snapshot-design.md`（保存の形と取得の手順）。表示の規則は `docs/genetics-design.md` と `docs/design.md`。

## Global Constraints

- msgspec と httpx 以外の依存パッケージを足さない（`pixi add` を打たない）
- 画面の文言は英語、コードのコメント・docstring・設計ノートは日本語
- `load_snapshot()` は schema 3 だけを受け付ける。schema 2 と `genetics.json` の読み込みは残さない。schema 2 のファイルは「Refresh it with: pixi run refresh」を含む ValueError にする
- 保存する関連遺伝子に下限を設けない。画面の初期値は 0.5、入力欄の範囲は 0 から 1、`step` は 0.01
- 発現の行は名前付きの dict のまま（配列にしない）。`{cell_id, median, specificity_score}` の 3 キー
- T 細胞の扱い（`CL_0000084` を祖先に持つ細胞を T 細胞にまとめる）は変えない
- 発現の取得の並列数は 2 のまま
- 集計と表示の規則は変えない。既存のテストは fixture を schema 3 に書き換えるだけで通す
- 保存は `json.dumps(allow_nan=False)` のまま。msgspec は NaN をエラーにせず書き出すので、保存には使わない
- 各タスクの最後に `pixi run --locked --no-install check` が通ること。コミットごとにテストと型チェックが通る状態を保つ。`Snapshot` の型に必須キーを足すタスクでは、`Snapshot` を作る全部の箇所（`refresh.main()`、`tests/core_fixture.py`、`tests/test_ui.py` の fixture、`tests/test_genetics.py` の `snapshot()`）を同じタスクで直す
- コミットメッセージは `<型>(<スコープ>): <説明>`。`Co-Authored-By` と `Claude-Session` の署名を付けない。`git commit` の前に Tomoya へ確認を取る（計画の承認は一括の許可ではない）
- Tomoya の dev サーバーが動いているあいだは `git checkout <commit>` を打たない
- 補助ファイル（試作のスクリプト、計測の出力）はリポジトリの外の一時領域に置く

## Review Focus

- `expression` の行の `cell_id` が `cells` に無い。`expression_metadata()` が ValueError にし、エラーも出ないまま名前の無い細胞が図に出ることを防ぐ（Task 1）
- 同じ版の既存ファイルから全遺伝子を再利用できるとき、`fetch_expression` が 1 回も呼ばれず、それでも `cells` が空にならない。先に読んだファイルにある遺伝子でも、後のファイルの行の細胞の定義は確かめる（Task 3）
- 関連遺伝子のページングで、総件数に達する前に空のページが返る。正常終了と見なさず失敗にし、保存しない（Task 3）
- 再利用の候補のファイルが壊れた JSON か、版が違う。読み飛ばして全部取り直し、失敗にしない（Task 3）
- 閾値の入力欄に 0、0.02、-0.1、空欄。0 と 0.02 はそのまま使い、-0.1 と空欄は 0.5 に戻して理由を表示する（Task 2）

---

## タスクの分け方

コミットごとにテストが通る状態を保つため、次の順で進める。

1. **Task 1** は保存の形（`ExpressionRow` の 3 キー化と `cells` の表）だけを変える。`genetics.py` の `merged_snapshot()` と `GeneticsSnapshot` は残し、`genetics_fixture.py` の発現の行も 3 キーにする。`refresh.main()` は `cells` を持つ schema 3 を書くようにする。読み込みは msgspec にする
2. **Task 2** で `associations` と `datasources` を `snapshot` に統合し、`GeneticsSnapshot`、`genetics.json` の読み込み、`create_app()` の `genetics` 引数を消す。`refresh.main()` は空の `associations` と `datasources` を書く（Task 3 で埋める）
3. **Task 3** で `refresh.py` を 1 本にし、`refresh_genetics.py` を消す。API の応答の検証を msgspec にする
4. **Task 4** で `query_api()` を httpx にする。取得の動作を変えるので、schema 3 の切り替えとは別のコミットにする
5. **Task 5** は文書と `pixi.toml`
6. **Task 6** は実データでの検証。`pixi run refresh` は Tomoya が打つ

Task 1 と Task 2 のコミットの時点では、`pixi run refresh` が書く `snapshot.json` は `associations` が空で、遺伝子ページに遺伝子が出ない。Task 3 まで通してから打つ。

## ファイル構成

| ファイル | 役割 |
| --- | --- |
| `autoimmune_atlas/models.py` | `ExpressionRow` を 3 キーに、`CellDefinition` を追加、`Snapshot` に `cells`、`associations`、`datasources` を追加、`GeneticsSnapshot` を削除 |
| `autoimmune_atlas/aggregation.py` | schema 3 の検証、`cells` からの大分類、`expression_metadata()` で `cell` 名を付ける |
| `autoimmune_atlas/snapshot.py` | msgspec による schema 3 の読み込みと検証 |
| `autoimmune_atlas/genetics.py` | `snapshot["associations"]` を読む集計だけを残す |
| `autoimmune_atlas/refresh.py` | 取得を 1 本にする。関連遺伝子、datasource の判定、`cells` の一貫性の検査、同じ版の発現の再利用、msgspec による応答の検証、httpx による呼び出し |
| `autoimmune_atlas/refresh_genetics.py` | 削除 |
| `autoimmune_atlas/ui/config.py` | `SCORE_FLOOR` の再輸出を削除 |
| `autoimmune_atlas/ui/genetics_layout.py` | `genetics_page(snapshot)`、閾値の入力欄の `min` と `step`、概要の文言 |
| `autoimmune_atlas/ui/genetics_callbacks.py` | `register_genetics_callbacks(application, snapshot)`、閾値の下限 0 |
| `autoimmune_atlas/ui/application.py` | `create_app()` から `genetics` と `genetics_error` を削除 |
| `app.py` | `genetics.json` の読み込みを削除 |
| `pixi.toml` | `refresh-genetics` タスクを削除 |
| `README.md`、`docs/design.md`、`docs/genetics-design.md`、`docs/snapshot-design.md` | 保存の形と取得の手順の記述を更新、検証の実測を記録 |
| `tests/core_fixture.py` | schema 3 の fixture。`cells`、`associations`、`datasources` を持つ |
| `tests/genetics_fixture.py` | Task 1 で発現の行を 3 キーにし、Task 2 で削除 |
| `tests/test_aggregation.py`、`tests/test_ui.py`、`tests/test_genetics.py`、`tests/test_genetics_ui.py` | fixture を schema 3 に書き換え、閾値 0 のテストを追加 |
| `tests/test_refresh.py` | 関連遺伝子の取得、`cells` の一貫性、再利用、`main()`、httpx の再試行のテストを追加 |
| `tests/test_refresh_genetics.py` | Task 1 で `_fetch` を 3 要素にし、Task 3 で削除（テストは `test_refresh.py` へ移す） |

---

### Task 1: 発現の行を 3 キーにし、細胞の定数を `cells` の表に移す

**Files:**
- Modify: `autoimmune_atlas/models.py:78-90`（`ExpressionRow`、`ExpressionMetadata`）、`:114-137`（`Snapshot`、`AggregationSnapshot`、`CoreSnapshot`）
- Modify: `autoimmune_atlas/aggregation.py:297-299`（`_validate_snapshot`）、`:330`（エラー文）、`:346-369`（`_cell_memberships`）、`:490-510`（`expression_metadata`）
- Modify: `autoimmune_atlas/snapshot.py`
- Modify: `autoimmune_atlas/refresh.py:320-362`（`extract_expression`）、`:408-437`（`fetch_expression`）、`:545-572`（`main()` の発現と保存）
- Modify: `autoimmune_atlas/refresh_genetics.py:232-236`
- Modify: `tests/core_fixture.py`、`tests/genetics_fixture.py`、`tests/test_aggregation.py`、`tests/test_ui.py`、`tests/test_genetics.py`、`tests/test_refresh.py`、`tests/test_refresh_genetics.py:129-130`
- Test: `tests/test_aggregation.py`、`tests/test_ui.py`、`tests/test_refresh.py`

**Interfaces:**
- Consumes: なし（最初のタスク）
- Produces:
    - `models.ExpressionRow(TypedDict)`: `cell_id: str`、`median: float | None`、`specificity_score: float | None`
    - `models.CellDefinition(TypedDict)`: `name: str`、`parent_id: str | None`、`parent: str | None`、`ancestor_ids: list[str]`
    - `models.ExpressionMetadata(ExpressionRow)`: `cell: str`、`target_median: float | None`
    - `models.Snapshot`、`models.CoreSnapshot` に `cells: dict[str, CellDefinition]`。`models.AggregationSnapshot` に `cells: NotRequired[dict[str, CellDefinition]]`
    - `aggregation._validate_snapshot()` は `schema == 3` を要求する
    - `aggregation.expression_metadata(snapshot)` の返す行は `cell` 名を持つ。`cells` に無い `cell_id` は ValueError
    - `snapshot.load_snapshot()` は msgspec で `Snapshot` として読み、`schema == 3` を要求する
    - `refresh.extract_expression(rows) -> tuple[list[ExpressionRow], dict[str, CellDefinition]]`
    - `refresh.fetch_expression(target_id) -> tuple[str, list[ExpressionRow], dict[str, CellDefinition]]`
    - `refresh.merge_cells(cells: dict[str, CellDefinition], found: Mapping[str, CellDefinition]) -> None`（食い違えば ValueError）

- [ ] **Step 1: `models.py` を書き換える**

`ExpressionRow` と `ExpressionMetadata`（78 行目から 90 行目）を次に置き換える。

```python
class ExpressionRow(TypedDict):
    """保存する発現の行。細胞の定数は Snapshot の cells にある。"""

    cell_id: str
    median: float | None
    specificity_score: float | None


class CellDefinition(TypedDict):
    """細胞 1 つの定数。遺伝子が違っても同じなので 1 回だけ保存する。"""

    name: str
    parent_id: str | None
    parent: str | None
    ancestor_ids: list[str]


class ExpressionMetadata(ExpressionRow):
    """実行時の行。cells から引いた cell 名と、標的内の中央値を持つ。"""

    cell: str
    target_median: float | None
```

`Snapshot`、`AggregationSnapshot`、`CoreSnapshot` に `cells` を足す。

```python
class Snapshot(TypedDict):
    schema: int
    root: str | RootIdentity
    roots: NotRequired[list[SnapshotRoot]]
    data_version: NotRequired[DataVersion]
    retrieved_at: str
    source: str
    diseases: list[Disease]
    records: list[DrugRecord]
    cells: dict[str, CellDefinition]
    expression: dict[str, list[ExpressionRow]]


class AggregationSnapshot(TypedDict):
    schema: int
    diseases: NotRequired[list[Disease]]
    records: NotRequired[list[DrugRecord]]
    cells: NotRequired[dict[str, CellDefinition]]
    expression: NotRequired[dict[str, list[ExpressionRow]]]


class CoreSnapshot(TypedDict):
    schema: int
    diseases: list[Disease]
    records: list[DrugRecord]
    cells: dict[str, CellDefinition]
    expression: dict[str, list[ExpressionRow]]
```

`GeneticsSnapshot` はこのタスクでは触らない。

- [ ] **Step 2: `aggregation.py` を書き換える**

`_validate_snapshot` を schema 3 にする。

```python
def _validate_snapshot(snapshot: SnapshotInput) -> None:
    if snapshot.get("schema") != 3:
        raise ValueError("snapshot must use schema 3")
```

`filtered_records` の 330 行目のエラー文を `"有効成分 ID または名称がありません"` にする。

`_cell_memberships` を `cells` の表から作る形にする。

```python
def _cell_memberships(
    snapshot: SnapshotInput,
) -> dict[str, CellMembership]:
    """cells の表から、細胞ごとの名前、大分類、祖先を作る。"""
    cells: dict[str, CellMembership] = {}
    for cell_id, definition in snapshot.get("cells", {}).items():
        ancestors = tuple(definition.get("ancestor_ids") or ())
        if cell_id == T_CELL_ID or T_CELL_ID in ancestors:
            group_id, group_name = T_CELL_ID, "T cell"
        else:
            group_id = definition.get("parent_id") or cell_id
            group_name = definition.get("parent") or definition["name"]
        cells[cell_id] = {
            "id": cell_id,
            "name": definition["name"],
            "group_id": group_id,
            "group_name": group_name,
            "ancestor_ids": tuple(sorted(set(ancestors))),
        }
    return cells
```

`expression_metadata` で `cell` 名を付け、`cells` に無い細胞を失敗にする。

```python
def expression_metadata(
    snapshot: SnapshotInput,
) -> dict[tuple[str, str], ExpressionMetadata]:
    """各標的・細胞の発現値に cell 名と、全参照細胞から求めた標的内中央値を付けて返す。"""
    _validate_snapshot(snapshot)
    cells = _cell_memberships(snapshot)
    all_cells = set(cells)
    result: dict[tuple[str, str], ExpressionMetadata] = {}
    for target_id, rows in snapshot.get("expression", {}).items():
        by_cell = {row["cell_id"]: row for row in rows}
        values = [row["median"] for row in rows]
        numeric_values = [value for value in values if value is not None]
        target_median = (
            median(numeric_values)
            if set(by_cell) == all_cells
            and values
            and len(numeric_values) == len(values)
            else None
        )
        for row in rows:
            cell_id = row["cell_id"]
            if cell_id not in cells:
                raise ValueError(f"cells に無い細胞 ID が発現にあります: {cell_id}")
            result[target_id, cell_id] = {
                **row,
                "cell": cells[cell_id]["name"],
                "target_median": target_median,
            }
    return result
```

- [ ] **Step 3: `snapshot.py` を msgspec で書き直す**

ファイル全体を次にする。手書きの `isinstance` の検証は消え、発現の行 1 つ 1 つまで型が確かめられる。

```python
"""保存済みスナップショットを検証して読み込む。"""

from pathlib import Path

import msgspec

from autoimmune_atlas.models import Snapshot

DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "snapshot.json"
REFRESH_HINT = "Refresh it with: pixi run refresh"


def load_snapshot(path: Path = DATA_PATH) -> Snapshot | None:
    """schema 3 の保存済みスナップショットを読む。無ければ None、形が違えば ValueError。"""
    if not path.exists():
        return None
    try:
        snapshot = msgspec.json.decode(path.read_bytes(), type=Snapshot)
    except msgspec.ValidationError as error:
        # schema 2 のファイルは cells が無いのでここに来る。
        raise ValueError(
            f"Unsupported snapshot format ({error}). {REFRESH_HINT}"
        ) from error
    except msgspec.DecodeError as error:
        raise ValueError(f"snapshot.json is not valid JSON: {error}") from error
    if snapshot["schema"] != 3:
        raise ValueError(
            f"Unsupported snapshot schema {snapshot['schema']}. {REFRESH_HINT}"
        )
    return snapshot
```

`app.py` の `except (OSError, json.JSONDecodeError, ValueError, TypeError, KeyError)` はそのままでよい（`ValueError` に包んでいる）。

`msgspec.json.decode(..., type=Snapshot)` は、TypedDict の `NotRequired`、`str | RootIdentity` のような union、`float | None` を受け付け、余分なキーは無視する。この版（0.21.1）で確かめてある。

- [ ] **Step 4: `refresh.py` の `extract_expression`、`fetch_expression`、`merge_cells`、`main()` の保存を直す**

`extract_expression` を、3 キーの行と細胞の定数の 2 つを返す形にする。

```python
def extract_expression(
    rows: list[ExpressionApiRow],
) -> tuple[list[ExpressionRow], dict[str, CellDefinition]]:
    """Tabula Sapiens の細胞型別 pseudobulk の中央値と、細胞の定数を取り出す。"""
    result: dict[str, ExpressionRow] = {}
    cells: dict[str, CellDefinition] = {}
    for row in rows:
        cell = row.get("celltypeBiosample")
        if (
            row["datasourceId"] != "tabula_sapiens"
            or not cell
            or row.get("tissueBiosample")
        ):
            continue
        if row["unit"] != "CPM(pseudobulk sum[counts])":
            raise ValueError(f"発現単位が変わりました: {row['unit']}")
        value = row["median"]
        if value is not None and (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0
        ):
            raise ValueError("発現中央値が不正です")
        specificity = row.get("specificity_score")
        if specificity is not None and (
            isinstance(specificity, bool)
            or not isinstance(specificity, (int, float))
            or not math.isfinite(specificity)
            or not 0 <= specificity <= 1
        ):
            raise ValueError("特異性スコアが不正です")
        cell_id = cell["biosampleId"]
        if cell_id in result:
            raise ValueError(f"細胞型別の発現行が重複しています: {cell_id}")
        parent = row.get("celltypeBiosampleParent")
        result[cell_id] = {
            "cell_id": cell_id,
            "median": value,
            "specificity_score": specificity,
        }
        cells[cell_id] = {
            "name": cell["biosampleName"],
            "parent_id": parent["biosampleId"] if parent else None,
            "parent": parent["biosampleName"] if parent else None,
            "ancestor_ids": list(cell.get("ancestors") or []),
        }
    return list(result.values()), cells
```

`fetch_expression` の戻り値を 3 要素にする。ページングは変えない。

```python
def fetch_expression(
    target_id: str,
) -> tuple[str, list[ExpressionRow], dict[str, CellDefinition]]:
    """発現データを全ページ取得してから、細胞型別の値と細胞の定数を選ぶ。"""
    ...（ページングは今のまま）...
    kept, cells = extract_expression(rows)
    return target_id, kept, cells
```

`merge_cells` を `save_snapshot` の前に足す。

```python
def merge_cells(
    cells: dict[str, CellDefinition], found: Mapping[str, CellDefinition]
) -> None:
    """細胞の定数を表へ足し、既にある定義と食い違えば失敗にする。"""
    for cell_id, definition in found.items():
        if cell_id in cells and cells[cell_id] != definition:
            raise ValueError(f"細胞型の名前か親分類が標的間で一貫しません: {cell_id}")
        cells[cell_id] = definition
```

`main()` の発現の取得と保存（545 行目から 572 行目）を次にする。`Snapshot` に `cells` が必須になったので、辞書に入れないと型チェックが失敗する。

```python
expression: dict[str, list[ExpressionRow]] = {}
cells: dict[str, CellDefinition] = {}
with ThreadPoolExecutor(max_workers=2) as pool:
    for index, (target, rows, found) in enumerate(
        pool.map(fetch_expression, target_ids), 1
    ):
        expression[target] = rows
        merge_cells(cells, found)
        print(f"細胞型別発現: {index}/{len(target_ids)}", flush=True)
final_version = cast(_VersionResponse, cast(object, query_api(version_query)))["meta"][
    "dataVersion"
]
if final_version != version:
    raise ValueError("取得中にデータの版が変わりました。再取得してください")
snapshot: Snapshot = {
    "schema": 3,
    "root": ROOT_ID,
    "roots": roots,
    "data_version": version,
    "retrieved_at": datetime.now(UTC).isoformat(),
    "source": f"https://{API_HOST}{API_PATH}",
    "diseases": sorted(diseases, key=lambda d: d["name"]),
    "records": list(records.values()),
    "cells": cells,
    "expression": expression,
}
cell_catalog(snapshot)
save_snapshot(DATA_PATH, snapshot)
```

`import` に `CellDefinition` を足し、`from collections.abc import Callable, Mapping` はそのまま使う。

`refresh_genetics.py:234` の `pool.map(fetch_expression, needed)` の受け側を `for index, (target, rows, _cells) in enumerate(...)` と 3 つにする。`refresh_genetics.main()` は `load_snapshot()` を呼ぶので、Task 1 のコミットの時点では schema 3 の `snapshot.json` が無いと止まる。Task 3 で消すので直さない。

`tests/test_refresh_genetics.py:129-130` の `_fetch` を 3 要素にする。

```python
def _fetch(target: str) -> tuple[str, list[Mapping[str, object]], dict[str, object]]:
    return target, [{"target_id": target}], {}
```

- [ ] **Step 5: `tests/core_fixture.py` と `tests/genetics_fixture.py` を schema 3 にする**

`core_fixture.py` は `"schema": 3` にし、`"records"` の後ろに `cells` を足し、`expression` の各行を 3 キーにする。

```python
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
        },
```

`genetics_fixture.py` の `expression["G9"]` も 3 キーにする。

```python
        "expression": {
            "G9": [
                {"cell_id": "T4", "median": 5.0, "specificity_score": 0.9},
                {"cell_id": "T8", "median": 0.0, "specificity_score": 0.0},
                {"cell_id": "B1", "median": 0.0, "specificity_score": 0.0},
                {"cell_id": "X", "median": 0.0, "specificity_score": 0.0},
            ]
        },
```

`tests/test_genetics.py:53-55` の行を `{"cell_id": "T4", "median": 99.0, "specificity_score": None}` にする。

- [ ] **Step 6: `tests/test_ui.py` と `tests/test_aggregation.py` の inline の snapshot を schema 3 にする**

行番号ではなく `grep` で見つける。

```bash
grep -n '"ancestor_ids"\|"parent_id"\|"schema": 2' tests/
```

`tests/test_ui.py` の `ui_snapshot`（360 行目付近）は、`"schema": 3` にし、`"records": records,` の後ろに `cells` を足し、行を 3 キーにする。

```python
        "cells": {
            "CL_B_ONE": {
                "name": "memory B cell",
                "parent_id": "CL_B_GROUP",
                "parent": "B cell",
                "ancestor_ids": ["CL_B_GROUP"],
            },
            "CL_B_TWO": {
                "name": "naive B cell",
                "parent_id": "CL_B_GROUP",
                "parent": "B cell",
                "ancestor_ids": ["CL_B_GROUP"],
            },
            "CL_T_ONE": {
                "name": "CD8-positive T cell",
                "parent_id": "CL_T_PARENT",
                "parent": "T lymphocyte",
                "ancestor_ids": [atlas.T_CELL_ID],
            },
        },
        "expression": {
            "ENSG_TARGET_1": [
                {"cell_id": "CL_B_ONE", "median": 2.0, "specificity_score": 0.8},
                {"cell_id": "CL_B_TWO", "median": 0.1, "specificity_score": 0.2},
                {"cell_id": "CL_T_ONE", "median": 0.2, "specificity_score": 0.1},
            ],
            "ENSG_TARGET_2": [
                {"cell_id": "CL_B_ONE", "median": None, "specificity_score": None},
                {"cell_id": "CL_B_TWO", "median": 0.6, "specificity_score": 0.75},
                {"cell_id": "CL_T_ONE", "median": 1.0, "specificity_score": 0.7},
            ],
        },
```

`tests/test_aggregation.py` の 2 つの inline の `AggregationSnapshot`（`test_catalog_orders_groups_and_cells`、`test_catalog_uses_curated_cell_state_order_with_lineage_fallback` の中）は、`rows` を作る箇所を `cells` と `rows` の 2 つに分ける。2 つ目のテストは次の形になる。

```python
        cells: dict[str, CellDefinition] = {
            cell_id: {
                "name": name,
                "parent_id": parent_id,
                "parent": parent,
                "ancestor_ids": ancestors,
            }
            for cell_id, name, parent_id, parent, ancestors in definitions
        }
        rows: list[ExpressionRow] = [
            {"cell_id": cell_id, "median": 1, "specificity_score": None}
            for cell_id, *_ in definitions
        ]
        snapshot: AggregationSnapshot = {
            "schema": 3,
            "diseases": [],
            "records": [],
            "cells": cells,
            "expression": {"G": rows},
        }
```

1 つ目のテストも同じ形にする。行の順序を逆にして同じ結果になることを確かめる箇所（`snapshot["expression"]["G"] = list(reversed(rows))`）は、`cells` の順序を逆にする形（`snapshot["cells"] = dict(reversed(list(cells.items())))`）に変える。大分類の順序は `cells` の並びに依存しないことを確かめるのが目的である。

`test_group_membership_must_be_consistent_across_targets` は、行が `ancestor_ids` を持たなくなったので成り立たない。次に置き換える。

```python
def test_expression_cell_missing_from_cells_table_fails(self) -> None:
    snapshot = core_snapshot()
    snapshot["expression"]["G2"].append(
        {"cell_id": "GHOST", "median": 1.0, "specificity_score": None}
    )
    with self.assertRaisesRegex(ValueError, "cells に無い"):
        expression_metadata(snapshot)


def test_metadata_rows_carry_cell_name_from_cells_table(self) -> None:
    metadata = expression_metadata(core_snapshot())
    self.assertEqual(metadata["G1", "T4"]["cell"], "CD4 T cell")
    self.assertEqual(metadata["G3", "X"]["cell"], "Novel cell")
```

`tests/test_ui.py` の `test_schema1_requires_refresh` に schema 2 の場合を足す。

```python
    def test_schema1_and_schema2_require_refresh(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            for old in (1, 2):
                path.write_text(json.dumps({"schema": old}), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "Refresh"):
                    snapshot.load_snapshot(path)
```

`tests/test_refresh.py` の `test_expression_keeps_specificity_parent_and_ancestors` を次に置き換え、`merge_cells` のテストを足す。

```python
def test_expression_splits_values_from_cell_definitions(self) -> None:
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
    rows, cells = extract_expression([base])
    self.assertEqual(rows, [{"cell_id": "C1", "median": 2.0, "specificity_score": 0.8}])
    self.assertEqual(
        cells,
        {
            "C1": {
                "name": "B cell",
                "parent_id": "P1",
                "parent": "Lymphocyte",
                "ancestor_ids": ["CL_1", "CL_2"],
            }
        },
    )
    self.assertEqual(
        extract_expression(
            [
                {**base, "datasourceId": "DICE"},
                {**base, "tissueBiosample": {"biosampleId": "T1"}},
            ]
        ),
        ([], {}),
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


def test_merge_cells_rejects_inconsistent_definitions(self) -> None:
    cells: dict[str, CellDefinition] = {}
    first: CellDefinition = {
        "name": "B cell",
        "parent_id": "P",
        "parent": "Lymphocyte",
        "ancestor_ids": ["CL_1"],
    }
    merge_cells(cells, {"C1": first})
    merge_cells(cells, {"C1": first, "C2": {**first, "name": "T cell"}})
    self.assertEqual(set(cells), {"C1", "C2"})
    with self.assertRaisesRegex(ValueError, "一貫"):
        merge_cells(cells, {"C1": {**first, "ancestor_ids": []}})
```

`CellDefinition` を `autoimmune_atlas.models` から、`merge_cells` を `autoimmune_atlas.refresh` から import する。

書き換えたら、`grep` が空になることを確かめる。

```bash
grep -n '"ancestor_ids"\|"parent_id"\|"schema": 2' tests/ autoimmune_atlas/
```

`cells` の辞書の中と、`refresh.py` の `extract_expression`、`tests/test_genetics.py:113` のコメント以外に残っていなければよい。

- [ ] **Step 7: 点検を全部通す**

Run: `pixi run --locked --no-install check`
Expected: format-check、lint、typecheck、test、test-help のすべてが PASS。

typecheck で `tests/test_ui.py` や `tests/test_aggregation.py` に `CellDefinition`、`ExpressionRow` の import が要ると言われたら足す。

- [ ] **Step 8: コミットする（Tomoya に確認してから）**

```bash
git add autoimmune_atlas/models.py autoimmune_atlas/aggregation.py autoimmune_atlas/snapshot.py autoimmune_atlas/refresh.py autoimmune_atlas/refresh_genetics.py tests/ pixi.toml pixi.lock
git commit -m "refactor(snapshot): 発現の行を 3 キーにし細胞の定数を cells の表へ移動"
```

`pixi.toml` と `pixi.lock` は Tomoya が msgspec と httpx を足した差分で、このコミットに含める。

---

### Task 2: 関連遺伝子を `snapshot` に統合し、遺伝子ページを `snapshot` だけで動かす

**Files:**
- Modify: `autoimmune_atlas/models.py`（`Snapshot`、`CoreSnapshot`、`AggregationSnapshot` に `associations` と `datasources`。`GeneticsSnapshot` を削除）
- Modify: `autoimmune_atlas/genetics.py`（`load_genetics`、`version_matches`、`merged_snapshot`、`GENETICS_PATH`、`GENETICS_SCHEMA`、`SCORE_FLOOR` を削除）
- Modify: `autoimmune_atlas/refresh.py`（`main()` の辞書に空の `associations` と `datasources`）
- Modify: `autoimmune_atlas/ui/config.py:3-4`、`autoimmune_atlas/ui/genetics_layout.py`、`autoimmune_atlas/ui/genetics_callbacks.py`、`autoimmune_atlas/ui/application.py`、`app.py`
- Modify: `autoimmune_atlas/refresh_genetics.py`（`genetics.py` から消える名前の import を止め、定数を自分で持つ。Task 3 で削除するまでの一時的な処置）
- Modify: `tests/core_fixture.py`、`tests/test_ui.py`（fixture に `associations` と `datasources`）、`tests/test_genetics.py`、`tests/test_genetics_ui.py`、`tests/test_refresh_genetics.py`
- Delete: `tests/genetics_fixture.py`
- Test: `tests/test_genetics.py`、`tests/test_genetics_ui.py`

**Interfaces:**
- Consumes: Task 1 の `ExpressionRow`、`CellDefinition`、`expression_metadata()`
- Produces:
    - `models.Snapshot`、`models.CoreSnapshot` に `associations: dict[str, list[GeneAssociation]]`、`datasources: list[str]`。`models.AggregationSnapshot` には `NotRequired`
    - `genetics.genes_for_disease(snapshot: SnapshotInput, disease_id: str, score_threshold: float) -> list[GeneAssociation]`
    - `genetics.summarize_genes(snapshot: SnapshotInput, score_threshold: object, threshold: object, *, method, specificity_threshold, level, cell_ids, disease_ids, metadata, catalog) -> list[SummaryRow]`（`genetics` の引数を消す）
    - `ui.genetics_layout.genetics_page(snapshot: Snapshot) -> html.Main`
    - `ui.genetics_layout.genetics_detail_panel(rows, selection, snapshot: Snapshot, *, score_threshold, threshold, method, specificity) -> html.Div`
    - `ui.genetics_layout.genetics_unavailable_page(reason: str) -> html.Main`（文言だけ変わる）
    - `ui.genetics_callbacks.register_genetics_callbacks(application: Dash, snapshot: Snapshot) -> None`
    - `ui.genetics_callbacks.effective_score(value) -> tuple[float, str | None]`（下限 0）
    - `ui.application.create_app(snapshot: Snapshot | None, error: Exception | None = None, *, assets_folder: str | Path) -> Dash`

- [ ] **Step 1: fixture を先に書く**

`tests/core_fixture.py` の `"records"` の前に、`genetics_fixture.py` から `associations` と `datasources` を移す。`core_snapshot()` は `CoreSnapshot` を返し続ける。

```python
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
```

`expression` に `G9` の 4 行（`genetics_fixture.py` の Task 1 で 3 キーにしたもの）を足す。`G9` は薬剤の標的ではなく関連遺伝子だけの遺伝子で、発現データを持つ。

`tests/test_ui.py` の `ui_snapshot` にも、`"records": records,` の後ろに足す。`create_app()` が薬剤ページと遺伝子ページの両方を組み立てるようになるので、無いと `genetics_page()` が KeyError になる。

```python
        "datasources": [],
        "associations": {"MONDO_RA_TEST": []},
```

`tests/genetics_fixture.py` を削除する。

- [ ] **Step 2: テストを書き換える**

`tests/test_genetics.py`:

- docstring を `"""snapshot の関連遺伝子の集計を検証する。"""` にする
- `snapshot()` のヘルパーは、`core_snapshot()` に `root`、`data_version`、`retrieved_at`、`source` を足す形のまま。`json`、`Path`、`TemporaryDirectory`、`GeneticsSnapshot`、`genetics_fixture` の import を消す
- `LoadTests` クラスを丸ごと消す（`load_genetics`、`version_matches`、`merged_snapshot` は無くなる）
- `SummarizeTests` の `data` と `setUp` を消し、`self.base = snapshot()` にする。`summarize_genes(self.base, self.data, ...)` を `summarize_genes(self.base, ...)` に、`genes_for_disease(self.data, ...)` を `genes_for_disease(self.base, ...)` にする
- `test_gene_without_expression_is_unknown_not_error` と `test_disease_missing_from_associations_is_unavailable` は、`data = genetics_snapshot()` の代わりに `base = snapshot()` を作って `base["associations"]` を変える
- 閾値 0 と 0.02 のテストを足す

```python
def test_zero_and_small_thresholds_keep_every_scored_gene(self) -> None:
    for threshold in (0, 0.02):
        genes = genetics.genes_for_disease(self.base, "D1", threshold)
        self.assertEqual([g["target_id"] for g in genes], ["G1", "G2", "G9", "G3"])
    rows = genetics.summarize_genes(
        self.base, 0, 0.5, level="group", disease_ids=["D1"]
    )
    t_cell = next(r for r in rows if r["cell_id"] == "CL_0000084")
    self.assertEqual(t_cell["denominator"], 4)
```

`tests/test_genetics_ui.py`:

- `genetics_fixture` の import を消す。`genetics_snapshot()` を使っていた箇所は `core_ui_snapshot()` を使う
- `genetics_layout.genetics_page(core_ui_snapshot(), genetics_snapshot())` → `genetics_layout.genetics_page(core_ui_snapshot())`
- `genetics_detail_panel([], None, genetics_snapshot(), ...)` → `genetics_detail_panel([], None, core_ui_snapshot(), ...)`
- `create_app(core_ui_snapshot(), assets_folder=ASSETS_PATH, genetics=genetics_snapshot())` → `create_app(core_ui_snapshot(), assets_folder=ASSETS_PATH)`
- `test_score_input_range_and_default` は `score.min == 0`、`score.step == 0.01` を確かめる。`_NumberInput` の Protocol に `step: float` を足す
- `test_score_input_minimum_follows_stored_floor` を消す
- `test_gene_rows_have_scores_datasource_columns_and_evidence_link` は `data = core_ui_snapshot()` にし、`data["datasources"]` はそのまま使える
- `test_detail_panel_without_selection_and_unavailable_page` の `assertIn("refresh-genetics", str(page))` を `assertIn("pixi run refresh", str(page))` にする
- `test_score_input_falls_back_below_floor_and_on_empty` を次に置き換える

```python
    def test_score_input_accepts_zero_and_falls_back_on_negative_or_empty(self) -> None:
        self.assertEqual(effective_score(0), (0, None))
        self.assertEqual(effective_score(0.02), (0.02, None))
        self.assertEqual(
            effective_score(-0.1),
            (0.5, "Score threshold must be finite and between 0 and 1; using 0.5."),
        )
        self.assertEqual(effective_score(""), (0.5, None))
        self.assertEqual(effective_score(0.3), (0.3, None))
```

- `test_version_mismatch_and_missing_genetics_show_notice_but_keep_drug_page` と `test_genetics_read_error_is_shown_in_the_notice` を消す
- `test_missing_snapshot_keeps_both_pages` は `create_app(None, assets_folder=ASSETS_PATH)` にし、`assertIn("snapshot is not available", ...)` にする

`tests/test_refresh_genetics.py` は Task 3 で消すが、このタスクでは `from autoimmune_atlas.genetics import GENETICS_SCHEMA, SCORE_FLOOR` と `GeneticsSnapshot` の import が失敗する。`from autoimmune_atlas.refresh_genetics import GENETICS_SCHEMA, SCORE_FLOOR` に変え、`payload = cast(GeneticsSnapshot, ...)` を `payload = cast(dict[str, object], save.call_args.args[1])` にする。`payload["associations"]["D1"]` の添字は `cast(dict[str, list[dict[str, object]]], payload["associations"])["D1"]` を挟む。

- [ ] **Step 3: テストが失敗することを確かめる**

Run: `pixi run --locked --no-install test`
Expected: `test_genetics.py`、`test_genetics_ui.py` が import か引数の数で失敗する。

- [ ] **Step 4: `models.py`、`snapshot.py`、`refresh.py` を書き換える**

`models.py` の `Snapshot` と `CoreSnapshot` の `cells` の前に足す。

```python
    associations: dict[str, list[GeneAssociation]]
    datasources: list[str]
```

`AggregationSnapshot` には `NotRequired` で足す。`GeneAssociation` は末尾で定義されているので、`Snapshot` より前へ移す。`GeneticsSnapshot` を削除する。

`snapshot.py` は変えない。msgspec が `Snapshot` の型から `associations` と `datasources` を必須として検証する。

`refresh.py` の `main()` の辞書に、`"cells": cells,` の前に次を足す。

```python
        # Task 3 で関連遺伝子の取得を足すまでの仮の値。
        "associations": {},
        "datasources": [],
```

- [ ] **Step 5: `genetics.py` を書き換える**

ファイル全体を次にする。`summarize_genes` の本体（検証、metadata、catalog、集計）は今のものをそのまま使い、`genetics` の引数と `GeneticsSnapshot` の参照だけを変える。

```python
"""snapshot の genetic association を読み、遺伝子ページ用に集計する。"""

import math
from collections.abc import Mapping

from autoimmune_atlas import aggregation as atlas
from autoimmune_atlas.models import (
    CellCatalogEntry,
    ExpressionMetadata,
    GeneAssociation,
    SummaryRow,
)


def genes_for_disease(
    snapshot: atlas.SnapshotInput, disease_id: str, score_threshold: float
) -> list[GeneAssociation]:
    """閾値以上（等号を含む）の関連遺伝子をスコアの降順で返す。"""
    genes = [
        gene
        for gene in snapshot.get("associations", {}).get(disease_id, [])
        if gene["score"] >= score_threshold
    ]
    return sorted(genes, key=lambda gene: (-gene["score"], gene["target_id"]))


def summarize_genes(
    snapshot: atlas.SnapshotInput,
    score_threshold: object,
    threshold: object,
    *,
    method: str = "fixed",
    specificity_threshold: object = 0.5,
    level: str = "group",
    cell_ids: list[str] | None = None,
    disease_ids: list[str] | None = None,
    metadata: Mapping[tuple[str, str], ExpressionMetadata] | None = None,
    catalog: list[CellCatalogEntry] | None = None,
) -> list[SummaryRow]:
    """疾患・細胞ごとに、閾値以上の遺伝子のうち発現判定が陽性の数と割合を返す。

    薬剤ページの標的と同じ三値判定、大分類の集約、割合の NA の規則を使う。
    薬剤に関する列は空の値で埋め、図の共有に使う。
    引数の検証は、遺伝子が 0 件の疾患でも行われるよう、発現の判定より前に行う。
    metadata と catalog は、呼び出し元が先に計算したものを渡すと再計算を省ける。
    catalog は level に対応したものを渡す。cell_ids の絞り込みは渡した catalog にも行う。
    """
    ...（検証、metadata、catalog、selected_diseases は今のまま）...
    associations = snapshot.get("associations", {})
    output: list[SummaryRow] = []
    for disease in snapshot.get("diseases", []):
        if selected_diseases is not None and disease["id"] not in selected_diseases:
            continue
        genes = genes_for_disease(snapshot, disease["id"], score_threshold)
        targets = [gene["target_id"] for gene in genes]
        loaded = disease["id"] in associations
        ...（以降は今のまま）...
```

`atlas.SnapshotInput` は `aggregation.py` で `type SnapshotInput = ...` と定義されている。`AggregationSnapshot` の `associations` が `NotRequired` なので `.get("associations", {})` で読む。

- [ ] **Step 6: UI を書き換える**

`autoimmune_atlas/ui/config.py`: 3 行目と 4 行目（`SCORE_FLOOR` の再輸出とそのコメント）を消す。

`autoimmune_atlas/ui/genetics_callbacks.py`:

- import から `GeneticsSnapshot` と `SCORE_FLOOR` を消す
- `effective_score` を次にする

```python
def effective_score(value: NumberInput) -> tuple[float, str | None]:
    """空欄は初期値、0 未満と 1 超と不正値は理由を示して初期値へ戻す。"""
    return effective_number(value, DEFAULT_SCORE_THRESHOLD, "score threshold", 1)
```

- `register_genetics_callbacks(application: Dash, snapshot: Snapshot) -> None` にし、docstring を `"""遺伝子ページのコールバックを登録する。"""` にする。関数の中で `genetics` を渡していた箇所（`gene_data.summarize_genes(snapshot, genetics, ...)`、`gene_data.genes_for_disease(genetics, ...)`、`genetics_detail_panel(..., genetics, ...)`）は `snapshot` にする。`grep -n genetics autoimmune_atlas/ui/genetics_callbacks.py` で、変数 `genetics` を指すものを全部直す（`genetics-` の ID と `gene_data` は残す）

`autoimmune_atlas/ui/genetics_layout.py`:

- import から `GeneticsSnapshot` と `SCORE_FLOOR` を消す
- `genetics_detail_panel(rows, selection, snapshot: Snapshot, *, ...)` にし、`genetics["datasources"]` を `snapshot["datasources"]`、`gene_data.genes_for_disease(genetics, ...)` を `gene_data.genes_for_disease(snapshot, ...)` にする
- `genetics_page(snapshot: Snapshot) -> html.Main` にする。`genetics["associations"]` を `snapshot["associations"]`、`genetics["retrieved_at"]` を `snapshot["retrieved_at"]` にし、`score_floor = ...` の行とそのコメントを消す。概要の `"Gene–disease associations (score ≥ 0.1)"` を `"Gene–disease associations"` に、`len(set(snapshot["expression"]) | set(genetics["expression"]))` を `len(snapshot["expression"])` にする
- 閾値の入力欄を次にする

```python
(
    dcc.Input(
        id="genetics-score",
        type="number",
        min=0,
        max=1,
        step=0.01,
        value=DEFAULT_SCORE_THRESHOLD,
    ),
)
(
    "The score threshold keeps genes whose Open Targets genetic association score for the disease is at or above this value. Empty or out-of-range values use 0.5.",
)
```

- `genetics_unavailable_page` の docstring を `"""snapshot.json が無いか読めないときの案内。"""` に、案内の文を `"Refresh the data with pixi run refresh, then restart the app."` にする

`autoimmune_atlas/ui/application.py`:

```python
"""Dash アプリを組み立てる。"""

from pathlib import Path

from dash import Dash, Input, Output, dcc, html

from autoimmune_atlas.models import Snapshot
from autoimmune_atlas.ui.callbacks import register_callbacks
from autoimmune_atlas.ui.genetics_callbacks import register_genetics_callbacks
from autoimmune_atlas.ui.genetics_layout import genetics_page, genetics_unavailable_page
from autoimmune_atlas.ui.layout import dashboard_layout, unavailable_layout

NO_SNAPSHOT = "The snapshot is not available, so the genetics page cannot be shown."


def create_app(
    snapshot: Snapshot | None,
    error: Exception | None = None,
    *,
    assets_folder: str | Path,
) -> Dash:
    """保存済みデータまたは明示的な fixture から、2 ページの Dash アプリを作る。"""
    application = Dash(
        __name__,
        title="Autoimmune Atlas",
        suppress_callback_exceptions=True,
        assets_folder=str(assets_folder),
    )
    if snapshot is None:
        drug_page = unavailable_layout(error)
        gene_page = genetics_unavailable_page(NO_SNAPSHOT)
    else:
        drug_page = dashboard_layout(snapshot)
        gene_page = genetics_page(snapshot)
        register_genetics_callbacks(application, snapshot)
    # 初期表示は薬剤ページにし、URL の callback が正しいほうを出す。
    gene_page.hidden = True  # pyright: ignore[reportAttributeAccessIssue] - Dash の型定義に hidden がない。
    application.layout = html.Div(
        [dcc.Location(id="url", refresh=False), drug_page, gene_page]
    )
    if snapshot is not None:
        register_callbacks(application, snapshot)

    @application.callback(  # pyright: ignore[reportAny, reportUnknownMemberType] - Dash の callback デコレーターに型情報がない。
        Output("drug-page", "hidden"),
        Output("genetics-page", "hidden"),
        Input("url", "pathname"),
    )
    def switch_page(pathname: str | None) -> tuple[bool, bool]:
        """/genetics なら遺伝子ページを、それ以外は薬剤ページを出す。"""
        on_genetics = pathname == "/genetics"
        return on_genetics, not on_genetics

    return application
```

`app.py`: `load_genetics` の import と 21 行目から 28 行目の try ブロックを消し、`create_app(snapshot_data, load_error, assets_folder=BASE_DIR / "assets")` にする。`sys` の import は使わなくなるので消す。

`autoimmune_atlas/refresh_genetics.py`（Task 3 で消すまでの一時的な処置）: 10 行目の import を消し、次を自分で持つ。

```python
GENETICS_PATH = Path(__file__).resolve().parent.parent / "data" / "genetics.json"
GENETICS_SCHEMA = 1
SCORE_FLOOR = 0.1
```

`from pathlib import Path` を足す。`GeneticsSnapshot` の import も消し、`genetics: GeneticsSnapshot = {...}` を `genetics: dict[str, object] = {...}` にする。

- [ ] **Step 7: 点検を全部通す**

Run: `pixi run --locked --no-install check`
Expected: すべて PASS。

- [ ] **Step 8: コミットする（Tomoya に確認してから）**

```bash
git add autoimmune_atlas/models.py autoimmune_atlas/genetics.py autoimmune_atlas/refresh.py autoimmune_atlas/refresh_genetics.py autoimmune_atlas/ui/ app.py tests/
git commit -m "refactor(genetics): 関連遺伝子を snapshot.json へ統合し genetics.json の読み込みを削除"
```

---

### Task 3: 取得を 1 本にし、同じ版の発現を再利用する

**Files:**
- Modify: `autoimmune_atlas/refresh.py`
- Delete: `autoimmune_atlas/refresh_genetics.py`、`tests/test_refresh_genetics.py`
- Modify: `pixi.toml:11`（`refresh-genetics` の行を消す）
- Test: `tests/test_refresh.py`

**Interfaces:**
- Consumes: Task 1 の `extract_expression()`、`fetch_expression()`、`merge_cells()`、Task 2 の `Snapshot`（`associations`、`datasources`、`cells`）
- Produces:
    - `refresh.LEGACY_GENETICS_PATH: Path`（`data/genetics.json`。再利用のためだけに読む）
    - `refresh.PAGE_SIZE = 500`、`refresh.GENETIC_DATATYPE = "genetic_association"`、`refresh.ASSOCIATION_QUERY`、`refresh.EVIDENCE_QUERY`、`refresh.VERSION_QUERY`
    - `refresh.collect_associations(query, disease_ids: list[str]) -> dict[str, list[GeneAssociation]]`（下限なし。スコアを持たない遺伝子で止める。総件数に達する前の空ページは ValueError）
    - `refresh.classify_datasources(query, associations) -> list[str]`、`refresh.restrict_datasources(associations, datasources) -> dict[str, list[GeneAssociation]]`（`refresh_genetics.py` から移す。中身は同じ）
    - `refresh.reusable_expression(paths: Sequence[Path], version: DataVersion) -> tuple[dict[str, list[ExpressionRow]], dict[str, CellDefinition]]`
    - `refresh.main() -> None`。テストは `autoimmune_atlas.refresh.query_api`、`fetch_expression`、`save_snapshot`、`DATA_PATH`、`LEGACY_GENETICS_PATH`、`SCOPE_ROOTS`、`collect_associations` を `patch` する

- [ ] **Step 1: `refresh_genetics.py` のテストを `test_refresh.py` へ移し、偽の API をページサイズに従う形にする**

`tests/test_refresh_genetics.py` の `_row`、`fake_query`、`RefreshGeneticsTests` を `tests/test_refresh.py` に移す。import は `from autoimmune_atlas.refresh import (..., VERSION_QUERY, classify_datasources, collect_associations, main, restrict_datasources, reusable_expression)` にする。

`PAGES` は、ページごとの辞書ではなく、疾患ごとの全行の一覧にし、偽の API が `index` と `size` で切り出す。今の形（2 件ずつの固定ページ）だと、通常の `PAGE_SIZE`（500）では 1 ページ目で `count` に達して 2 ページ目を読まず、G3 が入らない。

```python
ROWS: dict[str, list[Mapping[str, object]]] = {
    "D1": [
        _row("G1", 0.9, {"gwas_credible_sets": 0.9, "europepmc": 0.3}),
        _row("G2", 0.5, {"eva": 0.5}),
        _row("G3", 0.05, {"eva": 0.05}),
        _row("G4", None, {"europepmc": 0.3}),
        _row("G5", None, {"europepmc": 0.2}),
    ],
    "D2": [_row("G6", None, {"europepmc": 0.3})],
}


def fake_query(query: str, variables: Mapping[str, object]) -> dict[str, object]:
    if "associatedTargets" in query:
        disease = str(variables["id"])
        index = int(str(variables["index"]))
        size = int(str(variables["size"]))
        rows = ROWS.get(disease, [])
        return {
            "disease": {
                "associatedTargets": {
                    "count": len(rows),
                    "rows": rows[index * size : (index + 1) * size],
                }
            }
        }
    if "evidences" in query:
        datasource = str(variables["datasource"])
        datatype = "literature" if datasource == "europepmc" else "genetic_association"
        return {
            "disease": {
                "evidences": {
                    "rows": [{"datasourceId": datasource, "datatypeId": datatype}]
                }
            }
        }
    raise AssertionError(query)
```

`RefreshGeneticsTests` を `AssociationTests` に改名し、次のようにする。

```python
class AssociationTests(unittest.TestCase):
    def test_every_scored_gene_is_kept_and_unscored_genes_stop_paging(self) -> None:
        associations = collect_associations(fake_query, ["D1", "D2"])
        # G3 は 0.05 でも残る。スコアを持たない G4 以降は入らない。
        self.assertEqual(
            [g["target_id"] for g in associations["D1"]], ["G1", "G2", "G3"]
        )
        self.assertEqual(associations["D2"], [])
        self.assertEqual(associations["D1"][0]["target"], "g1")

    def test_paging_stops_on_the_page_that_holds_the_first_unscored_gene(self) -> None:
        # PAGE_SIZE を 2 にすると D1 は 3 ページになる。G4 のある 2 ページ目で止まり、3 ページ目は読まない。
        calls: list[int] = []

        def counting(query: str, variables: Mapping[str, object]) -> dict[str, object]:
            if "associatedTargets" in query:
                calls.append(int(str(variables["index"])))
            return fake_query(query, variables)

        with patch("autoimmune_atlas.refresh.PAGE_SIZE", 2):
            associations = collect_associations(counting, ["D1"])
        self.assertEqual(
            [g["target_id"] for g in associations["D1"]], ["G1", "G2", "G3"]
        )
        self.assertEqual(calls, [0, 1])

    def test_empty_page_before_count_is_reached_fails(self) -> None:
        def truncated(query: str, variables: Mapping[str, object]) -> dict[str, object]:
            response = fake_query(query, variables)
            if "associatedTargets" in query:
                block = cast(
                    dict[str, object],
                    cast(dict[str, object], response["disease"])["associatedTargets"],
                )
                block["count"] = 10
                if int(str(variables["index"])) > 0:
                    block["rows"] = []
            return response

        # PAGE_SIZE=3 なら 1 ページ目は G1、G2、G3 で全部スコアがあり、
        # 2 ページ目が空で count 10 に達していないので失敗する。
        with (
            patch("autoimmune_atlas.refresh.PAGE_SIZE", 3),
            self.assertRaisesRegex(ValueError, "途中のページが空"),
        ):
            collect_associations(truncated, ["D1"])

    def test_datasources_are_classified_by_one_evidence_each(self) -> None:
        associations = collect_associations(fake_query, ["D1"])
        datasources = classify_datasources(fake_query, associations)
        self.assertEqual(datasources, ["eva", "gwas_credible_sets"])
        restricted = restrict_datasources(associations, datasources)
        self.assertEqual(
            restricted["D1"][0]["datasource_scores"], {"gwas_credible_sets": 0.9}
        )

    def test_missing_disease_raises(self) -> None:
        def missing(_query: str, _variables: Mapping[str, object]) -> dict[str, object]:
            return {"disease": None}

        with self.assertRaises(ValueError):
            collect_associations(missing, ["D1"])
```

- [ ] **Step 2: 再利用と `main()` のテストを書く**

```python
VERSION: Mapping[str, object] = {"year": "26", "month": "09", "iteration": None}
MODULE = "autoimmune_atlas.refresh"
CELL: CellDefinition = {
    "name": "B cell",
    "parent_id": None,
    "parent": None,
    "ancestor_ids": [],
}
SCHEMA2_ROW: dict[str, object] = {
    "cell_id": "C1",
    "cell": "B cell",
    "median": 2.0,
    "specificity_score": 0.8,
    "parent_id": None,
    "parent": None,
    "ancestor_ids": [],
}


class ReuseTests(unittest.TestCase):
    def test_reuses_same_version_rows_from_schema2_and_genetics_files(self) -> None:
        with TemporaryDirectory() as directory:
            snapshot_path = Path(directory) / "snapshot.json"
            genetics_path = Path(directory) / "genetics.json"
            snapshot_path.write_text(
                json.dumps(
                    {
                        "schema": 2,
                        "data_version": VERSION,
                        "expression": {"G1": [SCHEMA2_ROW]},
                    }
                ),
                encoding="utf-8",
            )
            genetics_path.write_text(
                json.dumps(
                    {
                        "schema": 1,
                        "data_version": VERSION,
                        "expression": {
                            "G1": [{**SCHEMA2_ROW, "median": 99.0}],
                            "G2": [{**SCHEMA2_ROW, "median": 0.5}],
                        },
                    }
                ),
                encoding="utf-8",
            )
            expression, cells = reusable_expression(
                [snapshot_path, genetics_path], VERSION
            )
        # snapshot.json を先に読むので G1 は 2.0 のまま。G2 は genetics.json から来る。
        self.assertEqual(
            expression,
            {
                "G1": [{"cell_id": "C1", "median": 2.0, "specificity_score": 0.8}],
                "G2": [{"cell_id": "C1", "median": 0.5, "specificity_score": 0.8}],
            },
        )
        self.assertEqual(cells, {"C1": CELL})

    def test_cell_definitions_are_checked_even_for_genes_already_reused(self) -> None:
        with TemporaryDirectory() as directory:
            snapshot_path = Path(directory) / "snapshot.json"
            genetics_path = Path(directory) / "genetics.json"
            snapshot_path.write_text(
                json.dumps(
                    {
                        "schema": 2,
                        "data_version": VERSION,
                        "expression": {"G1": [SCHEMA2_ROW]},
                    }
                ),
                encoding="utf-8",
            )
            genetics_path.write_text(
                json.dumps(
                    {
                        "schema": 1,
                        "data_version": VERSION,
                        "expression": {"G1": [{**SCHEMA2_ROW, "cell": "T cell"}]},
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "一貫"):
                reusable_expression([snapshot_path, genetics_path], VERSION)

    def test_reuses_schema3_cells_table(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "snapshot.json"
            path.write_text(
                json.dumps(
                    {
                        "schema": 3,
                        "data_version": VERSION,
                        "cells": {"C1": CELL},
                        "expression": {
                            "G1": [
                                {
                                    "cell_id": "C1",
                                    "median": 1.0,
                                    "specificity_score": None,
                                }
                            ]
                        },
                    }
                ),
                encoding="utf-8",
            )
            expression, cells = reusable_expression([path], VERSION)
        self.assertEqual(list(expression), ["G1"])
        self.assertEqual(cells, {"C1": CELL})

    def test_other_version_missing_and_broken_files_are_skipped(self) -> None:
        with TemporaryDirectory() as directory:
            other = Path(directory) / "other.json"
            other.write_text(
                json.dumps(
                    {
                        "schema": 3,
                        "data_version": {
                            "year": "26",
                            "month": "06",
                            "iteration": None,
                        },
                        "cells": {},
                        "expression": {"G1": []},
                    }
                ),
                encoding="utf-8",
            )
            broken = Path(directory) / "broken.json"
            broken.write_text("{not json", encoding="utf-8")
            missing = Path(directory) / "missing.json"
            self.assertEqual(
                reusable_expression([other, broken, missing], VERSION), ({}, {})
            )
```

`main()` のテスト。偽の API は、版、起点、疾患の薬剤、薬剤の標的の 4 種類に応答し、関連遺伝子と根拠は `fake_query` に渡す。

```python
def _api(versions: list[Mapping[str, object]]) -> Callable[..., dict[str, object]]:
    """main() が打つ 6 種類の照会に応答する。版の照会には versions を順に返す。"""

    def query(
        text: str, variables: Mapping[str, object] | None = None
    ) -> dict[str, object]:
        variables = variables or {}
        if text == VERSION_QUERY:
            return {"meta": {"dataVersion": versions.pop(0)}}
        if "descendants" in text:
            return {"disease": {"id": "R", "name": "root", "descendants": ["D1", "D2"]}}
        if "drugAndClinicalCandidates" in text:
            ids = cast(list[str], variables["ids"])
            return {
                "diseases": [
                    {
                        "id": disease_id,
                        "name": disease_id.lower(),
                        "parents": [],
                        "drugAndClinicalCandidates": {
                            "count": 1 if disease_id == "D1" else 0,
                            "rows": [
                                {
                                    "maxClinicalStage": "PHASE_3",
                                    "drug": {
                                        "id": "A",
                                        "name": "Alpha",
                                        "drugType": "Antibody",
                                        "parentMolecule": None,
                                    },
                                }
                            ]
                            if disease_id == "D1"
                            else [],
                        },
                    }
                    for disease_id in ids
                ]
            }
        if "mechanismsOfAction" in text:
            return {
                "drugs": [
                    {
                        "id": "A",
                        "mechanismsOfAction": {
                            "rows": [
                                {
                                    "mechanismOfAction": "blocks",
                                    "actionType": "BLOCKER",
                                    "targets": [{"id": "G1", "approvedSymbol": "g1"}],
                                    "references": [],
                                }
                            ]
                        },
                    }
                ]
            }
        return fake_query(text, variables)

    return query


def _fetch(target: str) -> tuple[str, list[ExpressionRow], dict[str, CellDefinition]]:
    return (
        target,
        [{"cell_id": "C1", "median": 1.0, "specificity_score": None}],
        {"C1": CELL},
    )


class MainTests(unittest.TestCase):
    """main() が関連遺伝子まで取り、同じ版の発現を再利用し、全件成功したときだけ保存することを確かめる。"""

    def _run(
        self,
        versions: list[Mapping[str, object]],
        fetch: MagicMock,
        save: MagicMock,
        directory: str,
    ) -> None:
        with (
            patch(f"{MODULE}.SCOPE_ROOTS", (("R", True),)),
            patch(f"{MODULE}.DATA_PATH", Path(directory) / "snapshot.json"),
            patch(f"{MODULE}.LEGACY_GENETICS_PATH", Path(directory) / "genetics.json"),
            patch(f"{MODULE}.query_api", side_effect=_api(versions)),
            patch(f"{MODULE}.fetch_expression", fetch),
            patch(f"{MODULE}.save_snapshot", save),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            main()

    def test_fetches_union_of_targets_and_genes_and_saves_schema3(self) -> None:
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with TemporaryDirectory() as directory:
            self._run([VERSION, VERSION], fetch, save, directory)
        # 薬剤の標的 G1 と関連遺伝子 G1、G2、G3 の和集合
        self.assertEqual(
            sorted(c.args[0] for c in fetch.call_args_list), ["G1", "G2", "G3"]
        )
        save.assert_called_once()
        payload = cast(dict[str, object], save.call_args.args[1])
        self.assertEqual(payload["schema"], 3)
        self.assertEqual(payload["data_version"], VERSION)
        self.assertEqual(payload["datasources"], ["eva", "gwas_credible_sets"])
        self.assertEqual(payload["cells"], {"C1": CELL})
        associations = cast(dict[str, list[dict[str, object]]], payload["associations"])
        self.assertEqual(
            [(g["target_id"], g["datasource_scores"]) for g in associations["D1"]],
            [
                ("G1", {"gwas_credible_sets": 0.9}),
                ("G2", {"eva": 0.5}),
                ("G3", {"eva": 0.05}),
            ],
        )
        self.assertEqual(associations["D2"], [])
        self.assertEqual(
            sorted(cast(dict[str, object], payload["expression"])), ["G1", "G2", "G3"]
        )

    def test_reuses_same_version_expression_and_skips_fetch(self) -> None:
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with TemporaryDirectory() as directory:
            (Path(directory) / "snapshot.json").write_text(
                json.dumps(
                    {
                        "schema": 2,
                        "data_version": VERSION,
                        "expression": {
                            "G1": [{**SCHEMA2_ROW, "median": 7.0}],
                            "G2": [SCHEMA2_ROW],
                            "G3": [SCHEMA2_ROW],
                        },
                    }
                ),
                encoding="utf-8",
            )
            self._run([VERSION, VERSION], fetch, save, directory)
        fetch.assert_not_called()
        payload = cast(dict[str, object], save.call_args.args[1])
        self.assertEqual(payload["cells"], {"C1": CELL})
        expression = cast(dict[str, list[dict[str, object]]], payload["expression"])
        self.assertEqual(expression["G1"][0]["median"], 7.0)

    def test_other_version_is_refetched(self) -> None:
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with TemporaryDirectory() as directory:
            (Path(directory) / "snapshot.json").write_text(
                json.dumps(
                    {
                        "schema": 3,
                        "data_version": {
                            "year": "26",
                            "month": "06",
                            "iteration": None,
                        },
                        "cells": {"C1": CELL},
                        "expression": {"G1": [], "G2": [], "G3": []},
                    }
                ),
                encoding="utf-8",
            )
            self._run([VERSION, VERSION], fetch, save, directory)
        self.assertEqual(len(fetch.call_args_list), 3)

    def test_version_change_during_refresh_raises_and_saves_nothing(self) -> None:
        later = {"year": "26", "month": "12", "iteration": None}
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with TemporaryDirectory() as directory, self.assertRaises(ValueError):
            self._run([VERSION, later], fetch, save, directory)
        save.assert_not_called()

    def test_failed_association_paging_saves_nothing(self) -> None:
        fetch, save = MagicMock(side_effect=_fetch), MagicMock()
        with (
            TemporaryDirectory() as directory,
            patch(
                f"{MODULE}.collect_associations",
                side_effect=ValueError("途中のページが空"),
            ),
            self.assertRaisesRegex(ValueError, "途中のページが空"),
        ):
            self._run([VERSION, VERSION], fetch, save, directory)
        fetch.assert_not_called()
        save.assert_not_called()

    def test_inconsistent_cells_across_genes_raise(self) -> None:
        def inconsistent(
            target: str,
        ) -> tuple[str, list[ExpressionRow], dict[str, CellDefinition]]:
            cell: CellDefinition = {**CELL, "name": target}
            return (
                target,
                [{"cell_id": "C1", "median": 1.0, "specificity_score": None}],
                {"C1": cell},
            )

        fetch, save = MagicMock(side_effect=inconsistent), MagicMock()
        with (
            TemporaryDirectory() as directory,
            self.assertRaisesRegex(ValueError, "一貫"),
        ):
            self._run([VERSION, VERSION], fetch, save, directory)
        save.assert_not_called()
```

`contextlib`、`io`、`json`、`Callable`、`Mapping`、`MagicMock`、`cast` を import する。`Path` と `TemporaryDirectory` は既に import されている。`SCOPE_ROOTS` を `(("R", True),)` にするので、疾患の一覧は `D1`、`D2`、`R` の 3 つになり、`associations["R"]` は空の配列になる（テストは `D1` と `D2` だけを見る）。

`tests/test_refresh_genetics.py` を削除する。

- [ ] **Step 3: テストが失敗することを確かめる**

Run: `pixi run --locked --no-install test`
Expected: `test_refresh.py` が import エラーで失敗する（`collect_associations` などが `refresh` に無い）。

- [ ] **Step 4: `refresh.py` に関連遺伝子の取得を移し、応答の検証を msgspec にする**

`refresh_genetics.py` から、次を `refresh.py` へ移す。`_Score`、`_Target`（`refresh.py` に同名のクラスがあるので `_AssociationTarget` に改名）、`_AssociationRow`、`_AssociationBlock`、`_AssociationDisease`、`_AssociationResponse`、`_EvidenceRow`、`_EvidenceBlock`、`_EvidenceDisease`、`_EvidenceResponse`、`PAGE_SIZE`、`GENETIC_DATATYPE`、`ASSOCIATION_QUERY`、`EVIDENCE_QUERY`、`VERSION_QUERY`、`collect_associations`、`classify_datasources`、`restrict_datasources`。`refresh_genetics.py` の `Query`（`Callable[[str, Mapping[str, object]], dict[str, object]]`）は `_Query` として移す。`resolve_disease_ids` の `_DiseaseQuery` はそのまま残し、2 つの型を統合しない。統合すると、`tests/test_refresh.py` の既存の `fake_query`（第 2 引数が `dict[str, str]`）が `Mapping[str, object]` を受ける Callable に代入できず、typecheck で失敗する。

API の応答を TypedDict に当てる二重の `cast`（`cast(_X, cast(object, query(...)))` の形。`resolve_disease_ids`、`fetch_expression`、`main()`、移してくる `collect_associations` と `classify_datasources` に合わせて 8 か所）は、`msgspec.convert(query(...), _X)` にする。`msgspec.convert` は TypedDict の必須キーと型を確かめ、余分なキーは無視する。応答の形が違えば `msgspec.ValidationError` になる。`main()` の中では `ValueError` として扱えるように、`msgspec.ValidationError` は `ValueError` の派生であることを利用する（そのまま伝播させてよい）。

`collect_associations` は下限を消し、総件数に達する前の空ページを失敗にする。

```python
def collect_associations(
    query: _Query, disease_ids: list[str]
) -> dict[str, list[GeneAssociation]]:
    """疾患ごとに genetic スコアの降順で取り、スコアを持たない遺伝子に当たった時点で止める。"""
    output: dict[str, list[GeneAssociation]] = {}
    for disease_id in disease_ids:
        genes: list[GeneAssociation] = []
        index = 0
        while True:
            response = msgspec.convert(
                query(
                    ASSOCIATION_QUERY,
                    {"id": disease_id, "index": index, "size": PAGE_SIZE},
                ),
                _AssociationResponse,
            )
            disease = response["disease"]
            if disease is None:
                raise ValueError(f"疾患が見つかりません: {disease_id}")
            block = disease["associatedTargets"]
            stop = False
            for row in block["rows"]:
                genetic = next(
                    (
                        s["score"]
                        for s in row["datatypeScores"]
                        if s["id"] == GENETIC_DATATYPE
                    ),
                    None,
                )
                if genetic is None:
                    stop = True
                    break
                genes.append(
                    {
                        "target_id": row["target"]["id"],
                        "target": row["target"]["approvedSymbol"],
                        "score": genetic,
                        "datasource_scores": {
                            s["id"]: s["score"] for s in row["datasourceScores"]
                        },
                    }
                )
            index += 1
            if stop or index * PAGE_SIZE >= block["count"]:
                break
            if not block["rows"]:
                # HTTP と GraphQL が成功しても、ページが途中で切れることがある。保存させない。
                raise ValueError(f"関連遺伝子の途中のページが空です: {disease_id}")
        # API は降順で返すが、保存の形として並びを保証する
        output[disease_id] = sorted(genes, key=lambda g: (-g["score"], g["target_id"]))
    return output
```

`main()` にあった `version_query` のローカル変数は `VERSION_QUERY` に置き換える。

- [ ] **Step 5: `reusable_expression` を書く**

読む側の型を TypedDict で定め、msgspec に検証させる。schema 2 と `genetics.json` の行は `cell` などを持ち、schema 3 の行は持たないので、どちらも `NotRequired` にする。

```python
class _StoredRow(TypedDict):
    """再利用のために読む発現の行。schema 2 と genetics.json は cell などを持ち、schema 3 は持たない。"""

    cell_id: str
    median: float | None
    specificity_score: NotRequired[float | None]
    cell: NotRequired[str]
    parent_id: NotRequired[str | None]
    parent: NotRequired[str | None]
    ancestor_ids: NotRequired[list[str] | None]


class _StoredFile(TypedDict):
    data_version: NotRequired[DataVersion]
    cells: NotRequired[dict[str, CellDefinition]]
    expression: NotRequired[dict[str, list[_StoredRow]]]


LEGACY_GENETICS_PATH = DATA_PATH.parent / "genetics.json"


def reusable_expression(
    paths: Sequence[Path], version: DataVersion
) -> tuple[dict[str, list[ExpressionRow]], dict[str, CellDefinition]]:
    """同じ版の保存データから、発現の行と細胞の定数を取り出す。

    schema 2 の snapshot.json、schema 1 の genetics.json、schema 3 の snapshot.json を読む。
    発現の値は先に読んだファイルの遺伝子を優先する。
    細胞の定数は、後のファイルにある遺伝子の行でも全部確かめ、食い違えば失敗にする。
    無い、読めない、版が違うファイルは飛ばす。
    load_snapshot() は schema 3 しか受け付けないので、ここでは緩い型で読む。
    """
    expression: dict[str, list[ExpressionRow]] = {}
    cells: dict[str, CellDefinition] = {}
    for path in paths:
        try:
            data = msgspec.json.decode(path.read_bytes(), type=_StoredFile)
        except (OSError, msgspec.DecodeError, msgspec.ValidationError):
            continue
        if data.get("data_version") != version:
            continue
        merge_cells(cells, data.get("cells") or {})
        for target_id, rows in (data.get("expression") or {}).items():
            kept: list[ExpressionRow] = []
            for row in rows:
                if "cell" in row:
                    # schema 2 と genetics.json の行。細胞の定数を表へ移す。
                    merge_cells(
                        cells,
                        {
                            row["cell_id"]: {
                                "name": row["cell"],
                                "parent_id": row.get("parent_id"),
                                "parent": row.get("parent"),
                                "ancestor_ids": list(row.get("ancestor_ids") or []),
                            }
                        },
                    )
                kept.append(
                    {
                        "cell_id": row["cell_id"],
                        "median": row["median"],
                        "specificity_score": row.get("specificity_score"),
                    }
                )
            if target_id not in expression:
                expression[target_id] = kept
    return expression, cells
```

`from collections.abc import Callable, Mapping, Sequence` と `import msgspec` を足す。

- [ ] **Step 6: `main()` を 1 本にする**

`main()` の関連遺伝子より前（版、起点、疾患の薬剤、薬剤の標的、`records`）は今のままにし、`target_ids` を作るところから次にする。`resolve_disease_ids()` は `resolve_disease_ids(query_api)` と明示して呼ぶ。引数を省くと def 時の `query_api` に束縛され、テストの `patch` が効かない。

```python
disease_ids = [disease["id"] for disease in diseases]
associations = collect_associations(query_api, disease_ids)
print(
    f"疾患: {len(disease_ids)} / 関連遺伝子を持つ疾患: {sum(bool(g) for g in associations.values())}",
    flush=True,
)
datasources = classify_datasources(query_api, associations)
print(f"genetic association の datasource: {', '.join(datasources)}", flush=True)
associations = restrict_datasources(associations, datasources)
drug_targets = {
    target_id for record in records.values() if (target_id := record["target_id"])
}
gene_targets = {g["target_id"] for genes in associations.values() for g in genes}
target_ids = sorted(drug_targets | gene_targets)
reused, cells = reusable_expression([DATA_PATH, LEGACY_GENETICS_PATH], version)
expression: dict[str, list[ExpressionRow]] = {
    target_id: reused[target_id] for target_id in target_ids if target_id in reused
}
needed = [target_id for target_id in target_ids if target_id not in expression]
print(
    f"細胞型別発現: 再利用 {len(expression)} / 取得 {len(needed)} / 全体 {len(target_ids)}",
    flush=True,
)
with ThreadPoolExecutor(max_workers=2) as pool:
    for index, (target, rows, found) in enumerate(
        pool.map(fetch_expression, needed), 1
    ):
        expression[target] = rows
        merge_cells(cells, found)
        print(f"細胞型別発現: {index}/{len(needed)}", flush=True)
unknown_cells = {row["cell_id"] for rows in expression.values() for row in rows} - set(
    cells
)
if unknown_cells:
    raise ValueError(
        f"cells に無い細胞 ID が発現にあります: {sorted(unknown_cells)[:5]}"
    )
final_version = msgspec.convert(query_api(VERSION_QUERY), _VersionResponse)["meta"][
    "dataVersion"
]
if final_version != version:
    raise ValueError("取得中にデータの版が変わりました。再取得してください")
snapshot: Snapshot = {
    "schema": 3,
    "root": ROOT_ID,
    "roots": roots,
    "data_version": version,
    "retrieved_at": datetime.now(UTC).isoformat(),
    "source": f"https://{API_HOST}{API_PATH}",
    "diseases": sorted(diseases, key=lambda d: d["name"]),
    "records": list(records.values()),
    "associations": associations,
    "datasources": datasources,
    "cells": cells,
    "expression": expression,
}
cell_catalog(snapshot)
save_snapshot(DATA_PATH, snapshot)
print(
    f"保存: {DATA_PATH}\n疾患 {len(diseases)} / 薬剤 {len(drug_ids)} / 標的 {len(drug_targets)} / 関連遺伝子 {len(gene_targets)}",
    flush=True,
)
```

`main()` の docstring を `"""疾患、薬剤の標的、関連遺伝子、発現量を順番に取得し、schema 3 で保存する。"""` にする。`import` に `GeneAssociation` を足す。

`autoimmune_atlas/refresh_genetics.py` を削除する。

`pixi.toml` の `refresh-genetics = "python -m autoimmune_atlas.refresh_genetics"` の行を消す。

- [ ] **Step 7: 点検を全部通す**

Run: `pixi run --locked --no-install check`
Expected: すべて PASS。

`grep -rn "refresh_genetics\|refresh-genetics\|genetics.json\|GeneticsSnapshot\|SCORE_FLOOR\|cast(object" autoimmune_atlas/ tests/ pixi.toml app.py` が、`refresh.py` の `LEGACY_GENETICS_PATH` の定義と docstring、`tests/test_refresh.py` の再利用のテスト以外に何も返さないことを確かめる。

- [ ] **Step 8: コミットする（Tomoya に確認してから）**

```bash
git add autoimmune_atlas/refresh.py autoimmune_atlas/refresh_genetics.py tests/test_refresh.py tests/test_refresh_genetics.py pixi.toml
git commit -m "feat(refresh): 関連遺伝子の取得を refresh に統合し同じ版の発現を再利用"
```

---

### Task 4: `query_api()` を httpx にする

**Files:**
- Modify: `autoimmune_atlas/refresh.py:5-17`（import）、`:198-246`（`query_api`）
- Test: `tests/test_refresh.py`

**Interfaces:**
- Consumes: なし
- Produces: `refresh.query_api(query: str, variables: Mapping[str, object] | None = None, *, transport: httpx.BaseTransport | None = None) -> dict[str, object]`。呼び出し側の形は変えない

- [ ] **Step 1: テストを書く**

`httpx.MockTransport` で、再試行と失敗の扱いを確かめる。

```python
import httpx


def _graphql(bodies: list[httpx.Response]) -> httpx.MockTransport:
    """応答を順に返す。要求の本文が GraphQL の形であることも確かめる。"""

    def handle(request: httpx.Request) -> httpx.Response:
        payload = json.loads(request.content)
        assert set(payload) == {"query", "variables"}, payload
        return bodies.pop(0)

    return httpx.MockTransport(handle)


class QueryApiTests(unittest.TestCase):
    def test_returns_data_and_retries_transient_failures(self) -> None:
        transport = _graphql(
            [
                httpx.Response(502),
                httpx.Response(200, json={"data": {"meta": {"ok": True}}}),
            ]
        )
        with patch("autoimmune_atlas.refresh.time.sleep"):
            self.assertEqual(
                query_api("query{meta}", transport=transport), {"meta": {"ok": True}}
            )

    def test_graphql_errors_and_missing_data_raise_after_retries(self) -> None:
        for body in (
            {"errors": [{"message": "bad"}]},
            {"something": 1},
            [1, 2],
        ):
            transport = _graphql([httpx.Response(200, json=body)] * 3)
            with (
                self.subTest(body=body),
                patch("autoimmune_atlas.refresh.time.sleep"),
                self.assertRaises(RuntimeError),
            ):
                query_api("query{meta}", transport=transport)
```

`query_api` を import に足す。

- [ ] **Step 2: テストが失敗することを確かめる**

Run: `pixi run --locked --no-install test tests.test_refresh`
Expected: `query_api()` が `transport` を受け取らず TypeError で失敗する。

- [ ] **Step 3: `query_api` を書き直す**

`base64`、`http.client`、`os`、`urlparse` の import を消し、`import httpx` を足す。`time` は再試行の待ちで使うので残す。

```python
API_URL = f"https://{API_HOST}{API_PATH}"


def query_api(
    query: str,
    variables: Mapping[str, object] | None = None,
    *,
    transport: httpx.BaseTransport | None = None,
) -> dict[str, object]:
    """GraphQL を照会し、不完全な応答を失敗として扱う。3 回まで再試行する。

    proxy は httpx が環境変数（HTTPS_PROXY、https_proxy）から読む。
    transport はテストが偽の応答を差し込むためにある。
    """
    for attempt in range(3):
        try:
            with httpx.Client(timeout=45, transport=transport) as client:
                response = client.post(
                    API_URL, json={"query": query, "variables": variables or {}}
                )
            if response.status_code != 200:
                raise RuntimeError(f"Open Targets: HTTP {response.status_code}")
            raw = cast(object, response.json())
            if not isinstance(raw, dict):
                raise RuntimeError("Open Targets: response must be an object")
            payload = cast(dict[str, object], raw)
            if payload.get("errors") or "data" not in payload:
                raise RuntimeError(f"Open Targets: {payload.get('errors', payload)}")
            data = payload["data"]
            if not isinstance(data, dict):
                raise RuntimeError("Open Targets: data must be an object")
            return cast(dict[str, object], data)
        except (httpx.HTTPError, ValueError, RuntimeError):
            if attempt == 2:
                raise
            time.sleep(attempt + 1)
    raise RuntimeError("Open Targets の照会が失敗しました")
```

`main()` の `"source": f"https://{API_HOST}{API_PATH}"` は `API_URL` にする。

- [ ] **Step 4: 点検を全部通す**

Run: `pixi run --locked --no-install check`
Expected: すべて PASS。

sandbox の外で 1 回だけ実際の API を打ち、proxy を通ることを確かめる。Tomoya が打つ。

```bash
pixi run -- python -c "from autoimmune_atlas.refresh import query_api, VERSION_QUERY; print(query_api(VERSION_QUERY))"
```

- [ ] **Step 5: コミットする（Tomoya に確認してから）**

```bash
git add autoimmune_atlas/refresh.py tests/test_refresh.py
git commit -m "refactor(refresh): API の照会を httpx にし proxy の手書きの処理を削除"
```

---

### Task 5: 文書を保存の形と取得の手順に合わせる

**Files:**
- Modify: `README.md:17`、`README.md` のコード構成の木（`refresh_genetics.py` と `test_refresh_genetics.py` の行）、`README.md:230-255`（データ更新と検証）
- Modify: `docs/design.md:290`、`:302`
- Modify: `docs/genetics-design.md:25-27`、`:36-90`（「データは別ファイルに取得し、版が一致するときだけ表示する」の節）、`:95`、`:149`、`:157`、`:177`、`:192`、`:202`

**Interfaces:**
- Consumes: Task 3 と Task 4 の `pixi run refresh` の振る舞い
- Produces: なし

- [ ] **Step 1: `README.md` を直す**

17 行目を次にする。

```text
データ更新は `pixi run refresh` が `autoimmune_atlas.refresh` をモジュールとして実行する。
```

コード構成の木から `refresh_genetics.py` と `test_refresh_genetics.py` の行を消し、`genetics.py` の説明を `遺伝学的関連遺伝子の集計` にする。

「データ更新と検証」の節（230 行目付近から `pixi run test` の前まで）を次にする。

````markdown
## データ更新と検証

公開 API から取り直す。
初回は、発現データの取得に 1.5 時間から 2 時間かかる。

```bash
pixi run refresh
```

疾患、薬剤の標的、genetic association の関連遺伝子、細胞型別の発現を順に取り、全取得と件数の検証が成功したときだけ `data/snapshot.json` を置き換える。
途中で失敗した場合は前回のデータを保つ。
取得日時と Open Targets のデータ版はスナップショットに保存する。
取得時点の公開データを使うため、環境間で取得日時や版が異なると結果も変わり得る。
更新後はアプリを再起動する。

同じ版の発現データは、既存の `data/snapshot.json` と `data/genetics.json` から再利用し、取り直さない。
版が変わっていれば全部取り直す。
schema 2 の `snapshot.json` と `genetics.json` を使っている場合も、この更新で schema 3 の 1 ファイルに切り替わる。
切り替えが済んだら `data/genetics.json` は消してよい。
保存の形と取得の手順は[保存データの統合と一括取得の設計](docs/snapshot-design.md)にある。
````

- [ ] **Step 2: `docs/design.md` を直す**

290 行目の段落を次にする。

```markdown
遺伝学的関連遺伝子のページは、`snapshot.json` の `associations` を読む別の集計（`autoimmune_atlas/genetics.py`）で `/genetics` に置き、薬剤ページの発現判定と細胞分類と図を共有する。
```

302 行目の根拠の項目の末尾 `` `autoimmune_atlas/refresh_genetics.py` が取得し、`autoimmune_atlas/genetics.py` が集計する。`` を `` `autoimmune_atlas/refresh.py` が取得し、`autoimmune_atlas/genetics.py` が集計する。`` にする。

- [ ] **Step 3: `docs/genetics-design.md` を直す**

25 行目から 27 行目の段落（「保存する下限スコアは 0.1 とし」から「0.1 から 1 の範囲に限る。」まで）を次にする。

```markdown
保存する関連遺伝子に下限は設けず、スコアを持つ遺伝子を全部保存する。
画面の閾値は 0 から 1 の範囲で変えられる。
初期版では下限 0.1 を置いていたが、弱い関連が見えないことが分かったので無くした。
```

「データは別ファイルに取得し、版が一致するときだけ表示する」の節（見出しから、次の見出し「集計は発現判定と細胞分類だけを共有する」の直前まで）を、次の節に置き換える。旧版の大きさとメモリの記述は、この節の末尾に経緯として残す。

```markdown
## データは snapshot.json に保存し、下限を設けない

関連遺伝子は `pixi run refresh` が薬剤の記録と一緒に取得し、`data/snapshot.json` の `associations` と `datasources` に保存する。
保存の形と取得の手順は[保存データの統合と一括取得の設計](snapshot-design.md)にある。
初期版では `data/genetics.json` を別に置き、スコア 0.1 以上だけを保存していた。
2 つのファイルの版を突き合わせる手間と、下限のせいで弱い関連が見えないことが分かったので、1 ファイルに戻し、下限を無くした。

取得の規則は次のとおり。

- 疾患ごとに `Disease.associatedTargets(enableIndirect: false, orderByScore: "genetic_association")` を 500 件ずつ取得し、genetic association のスコアを持たない遺伝子に当たった時点でその疾患のページングを止める。総件数に達する前に空のページが返れば失敗にし、保存しない
- 遺伝子ごとに、`datatypeScores` の genetic association のスコアと、`datasourceScores` のうち genetic association に属する datasource のスコアを保存する。どの datasource が genetic association に属するかは、取得した結果から決める。全疾患の `datasourceScores` に現れた datasource ID ごとに、それを持つ疾患と遺伝子の組を 1 つ選び、`Disease.evidences(ensemblIds: [遺伝子], datasourceIds: [datasource], size: 1)` で根拠を 1 件取り、その `datatypeId` が `genetic_association` のものだけを残す。datasource の一覧をコードに書かないのは、版が上がって datasource が増えても取りこぼさないためである。API の `associationDatasources` は版 26.09 で空の配列を返したので、使わない

`orderByScore: "genetic_association"` の並びは、genetic association のスコアを持つ遺伝子を降順に置き、スコアを持たない遺伝子をその後ろに置く。
関節リウマチで確かめたところ、全件を取って数えた 0.1 以上の件数と、並べて 0.1 を下回るまで数えた件数が一致した。

`snapshot.json` が無いか読めないときは、遺伝子ページに「Refresh the data with pixi run refresh」と表示し、比較図を出さない。
薬剤ページも同じ理由で出ない。

遺伝子ごとの代表的な根拠（GWAS の study、L2G スコア、変異）は初期版では保存しない。
根拠の行は疾患と遺伝子の組ごとに数十件あり、疾患全体では数万件になる。
代表 1 件を選ぶと、それが疾患関連スコアの唯一の根拠と読まれる。
詳細の表からは Open Targets の evidence ページ（`/evidence/<遺伝子 ID>/<疾患 ID>`）へリンクし、根拠はそこで読む。

初期版の `genetics.json` は `snapshot.json` のおよそ 5 倍の大きさで、2 つのファイルを読むとプロセスのメモリは約 3 GB まで増え、画面の構成を作り始めるまでに約 9 秒かかった。
発現の行ごとに、細胞で決まる値（`ancestor_ids`、`parent`、`cell`）を繰り返し持っていたためである。
schema 3 で細胞の表を 1 回だけ持つようにした結果は、[保存データの統合と一括取得の設計](snapshot-design.md)の「検証」にある。
```

残りの行を次のように直す。

- 95 行目：`summarize_genes(snapshot, genetics, score_threshold, threshold, *, ...)` から 2 番目の引数 `genetics` を消し、`summarize_genes(snapshot, score_threshold, threshold, *, ...)` にする
- 149 行目の表の Genetic evidence の行：`0.1 から 1 の数値入力` を `0 から 1 の数値入力` にし、末尾の `既存の数値の検証は 0 未満しか弾かないので、下限 0.1 の検証を足す` を消す
- 157 行目：`` `genetics.json` 全体の値として `` を `` `snapshot.json` 全体の値として `` に、`閾値 0.1 以上の関連の数` を `関連の数` にする
- 177 行目と 192 行目：`` `genetics.json` `` を `` `snapshot.json` `` にする
- 202 行目：`tests/test_refresh_genetics.py` の項目を次にする

```markdown
- `tests/test_refresh.py`：偽の query 関数を渡し、ページングがスコアを持たない遺伝子で止まること、途中の空ページで失敗すること、同じ版の発現を再利用すること、全件成功時だけ保存することを確かめる
```

- [ ] **Step 4: Markdown の点検を通す**

Run: `pixi run --locked --no-install lint-markdown`
Expected: PASS。

`grep -n "refresh-genetics\|refresh_genetics\|test_refresh_genetics" README.md docs/design.md docs/genetics-design.md` が何も返さないことを確かめる。`docs/snapshot-design.md` と、この計画には、経緯として残る。

- [ ] **Step 5: コミットする（Tomoya に確認してから）**

```bash
git add README.md docs/design.md docs/genetics-design.md
git commit -m "docs(snapshot): 保存データの統合に合わせて README と設計ノートを更新"
```

---

### Task 6: 実データで検証し、実測を設計書に書く

**Files:**
- Modify: `docs/snapshot-design.md`（「検証」の節の末尾に実測を足す）

**Interfaces:**
- Consumes: Task 3 と Task 4 の `pixi run refresh`、Task 2 の遺伝子ページ
- Produces: なし

- [ ] **Step 1: Tomoya が `pixi run refresh` を打つ**

エージェントは打たない。
既存の `data/snapshot.json`（schema 2）と `data/genetics.json` が版 26.09 なら、約 2,400 遺伝子ぶんの取得で 30 分程度。
Open Targets の版が上がっていれば全部取り直しになり、1.5 時間から 2 時間かかる。

```bash
pixi run refresh
```

出力の「細胞型別発現: 再利用 N / 取得 M / 全体 K」の行を控える。
再利用が 0 なら版が変わっている。

- [ ] **Step 2: 大きさ、メモリ、起動時間を測る**

計測のスクリプトは一時領域（scratchpad）に置き、`pixi run --locked --no-install -- python <パス>` で打つ。

```python
import resource
import time

start = time.perf_counter()
from autoimmune_atlas.snapshot import load_snapshot

snapshot = load_snapshot()
assert snapshot is not None
loaded = time.perf_counter()
from autoimmune_atlas.ui.application import create_app

create_app(snapshot, assets_folder="assets")
built = time.perf_counter()
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30
print(f"read {loaded - start:.1f} s / build {built - loaded:.1f} s / rss {rss:.2f} GB")
print(
    f"genes {len(snapshot['expression'])} / cells {len(snapshot['cells'])} / "
    f"diseases with genes {sum(bool(g) for g in snapshot['associations'].values())}"
)
```

ファイルの大きさは `ls -lh data/snapshot.json` で見る。

- [ ] **Step 3: 実ブラウザーで関節リウマチの分母を確かめる**

Tomoya が `pixi run dev` でアプリを起動する。
エージェントはブラウザーで `http://127.0.0.1:8050/genetics` を開き、Genetic association score の入力欄を 0 にして Update を押し、詳細の疾患に rheumatoid arthritis を選ぶ。
比較図のホバーか詳細の表で、分母（関連遺伝子の数）が版 26.09 で 697 になることを確かめる。
版が違えば、`associations` の rheumatoid arthritis の配列の長さと同じであることを確かめる。

```bash
pixi run --locked --no-install -- python -c "
from autoimmune_atlas.snapshot import load_snapshot
s = load_snapshot()
ra = next(d for d in s['diseases'] if d['name'] == 'rheumatoid arthritis')
print(ra['id'], len(s['associations'][ra['id']]))
"
```

入力欄に 0.02 を入れても弾かれないことも確かめる。

- [ ] **Step 4: 実測を設計書に書く**

`docs/snapshot-design.md` の「検証」の節の末尾に、次の形で段落を足す。値は Step 1 から 3 の実測に置き換える。

```markdown
版 26.09 で測った結果は次のとおりである。
`snapshot.json` は NNN MB で、schema 2 の 2 ファイル（合わせて約 590 MB）から X 分の 1 になった。
読み込みに N 秒、画面の構成に N 秒かかり、プロセスのメモリは N GB だった。
再利用が働き、取得したのは N 遺伝子で、N 分かかった。
関節リウマチの閾値 0 の分母は 697 で、`associations` の配列の長さと一致した。
```

- [ ] **Step 5: 点検とコミット（Tomoya に確認してから）**

Run: `pixi run --locked --no-install lint-markdown`
Expected: PASS。

```bash
git add docs/snapshot-design.md
git commit -m "docs(snapshot): schema 3 の大きさとメモリと起動時間の実測を記録"
```

---

## 自己点検

- **仕様の網羅**：下限なし（Task 2、3）、`cells` の表（Task 1、3）、集計は `cells` を読む（Task 1）、`ExpressionRow` と `ExpressionMetadata` の分離（Task 1）、`load_snapshot` は schema 3 だけ（Task 1、2）、fixture の統合（Task 1、2）、取得 1 本と再利用（Task 3）、`refresh_genetics` の削除（Task 3）、`create_app` の引数の削除（Task 2）、画面の変更 4 点（Task 2）、全件成功時だけ保存（Task 3 の空ページの失敗と版の確認）、検証 5 点（`cells` の一貫性は Task 1 の `merge_cells` と Task 3 の再利用のテスト、`test_aggregation` の同値は Task 1 の fixture の書き換えで既存のテストが通ることで確かめる、実データの計測は Task 6）
- **コミットごとの整合**：Task 1 で `Snapshot.cells` を必須にすると同時に `refresh.main()` と全 fixture に `cells` を足す。Task 2 で `associations` と `datasources` を必須にすると同時に `refresh.main()`（仮の空の値）、`core_fixture.py`、`test_ui.py` の fixture に足す。`test_refresh_genetics.py` の `_fetch` は Task 1 で 3 要素にする
- **プレースホルダー**：Task 6 の実測の値だけが空欄で、Step 1 から 3 で埋める。それ以外に「後で」は無い
- **型の一貫**：`fetch_expression` の 3 要素の戻り値は Task 1 で決め、Task 3 の `main()` と `_fetch` が同じ形で受ける。`merge_cells` は Task 1 で定義し、Task 3 の `reusable_expression` と `main()` が使う。`reusable_expression` は `Sequence[Path]` を受け、`main()` は 2 要素の list を渡す。`genes_for_disease` と `summarize_genes` の第 1 引数は Task 2 で `atlas.SnapshotInput` にし、UI からは `Snapshot` を渡す。`_DiseaseQuery` と `_Query` は別のまま残す。`query_api` の `transport` は Task 4 で足し、`main()` は使わない
- **Codex のセカンドオピニオン（2026-09-27）で直した点**：Task 1 と Task 2 で `Snapshot` を作る全箇所を同じタスクで直す。`reusable_expression` は再利用済みの遺伝子の行でも細胞の定義を確かめる。偽の API は `index` と `size` で切り出す。総件数に達する前の空ページを失敗にする。文書の残りの参照（`genetics-design.md` の 25、95、149、202 行目、`design.md` の 302 行目、`README.md` の 47 行目）を Task 5 に入れ、Task 6 は追記だけにする
- **Review Focus**：5 項目それぞれに、Task 1（`test_expression_cell_missing_from_cells_table_fails`）、Task 3（`test_reuses_same_version_expression_and_skips_fetch`、`test_cell_definitions_are_checked_even_for_genes_already_reused`、`test_empty_page_before_count_is_reached_fails`、`test_other_version_missing_and_broken_files_are_skipped`）、Task 2（`test_score_input_accepts_zero_and_falls_back_on_negative_or_empty`）のテストがある
