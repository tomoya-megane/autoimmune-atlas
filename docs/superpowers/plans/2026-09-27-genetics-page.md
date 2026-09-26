# 遺伝学的関連遺伝子ページの実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Open Targets の genetic association で疾患と結びついた遺伝子について、既存の薬剤ページと同じ細胞型別の発現比較を `/genetics` に表示する。

**Architecture:** `refresh_genetics.py` が `data/genetics.json` を別に作り、`genetics.py` が既存の発現判定と細胞分類を呼んで遺伝子用の集計を返す。画面は `dcc.Location` で `/` と `/genetics` を切り替え、両ページのレイアウトを起動時に組み立てて `hidden` を入れ替える。比較図は既存の `build_figure` と `build_dot_figure` を `kind="gene"` で共有し、詳細は遺伝子用の表と既存の発現図で作る。

**Tech Stack:** Python 3.13、Dash 4、Plotly 7、dash-ag-grid、unittest、Ruff、basedpyright。新しい依存パッケージは入れない。

**Spec:** `docs/genetics-design.md`

## Global Constraints

- 新しい依存パッケージを導入しない（`pixi add` を打たない）
- 画面の文言は英語、コードのコメント・docstring・設計ノートは日本語
- 既存の `data/snapshot.json` のスキーマ、`autoimmune_atlas/refresh.py` の保存部分は変えない
- `create_app()` の既存の呼び方（`create_app(snapshot, error, assets_folder=...)`）と `tests/test_ui.py` はそのまま通る
- 遺伝子集合は datatype `genetic_association` のスコアが閾値以上（等号を含む）。保存する下限は 0.1、画面の初期値は 0.5、入力欄の範囲は 0.1 から 1
- `genetic_literature` は含めない。`enableIndirect: false`
- 遺伝子ページの部品の ID には `genetics-` を前に付ける
- 各タスクの最後に `pixi run --locked --no-install check` が通ること
- コミットメッセージは `<型>(<スコープ>): <説明>`。署名を付けない。`git commit` の前に Tomoya へ確認を取る（計画の承認は一括の許可ではない）
- `git checkout <commit>` を Tomoya の dev サーバーが動いているあいだに打たない
- Tomoya が読める言語（Python、Bash、JavaScript は最小限）で書く

## Review Focus

- `genetics.json` はあるが `data_version` が `snapshot.json` と違う。遺伝子ページは案内文を出し、薬剤ページは動く（Task 1、Task 7 で確かめる）
- 閾値以上の遺伝子が 0 件の疾患。`count` は 0、`percent` は `None`、ホバーは "No genes at or above the score threshold"（Task 2、Task 4）
- 閾値の入力欄に 0.05 や空欄。0.5 に戻し、適用した値を画面に記録する（Task 6）
- `genetics.json` の遺伝子が `snapshot.json` の発現データにも `genetics.json` の発現データにも無い。三値判定で未判定として数え、例外にしない（Task 2）
- `/genetics` を直接開いたとき、Plotly の図が幅 0 で描かれる。`resize` を送る JavaScript と、実ブラウザーの確認（Task 7）

---

## ファイル構成

| ファイル | 役割 |
| --- | --- |
| `autoimmune_atlas/models.py` | `GeneAssociation`、`GeneticsSnapshot` の型を足す |
| `autoimmune_atlas/genetics.py` | 新規。`genetics.json` の読み込みと検証、発現の結合、遺伝子集合の選択、`summarize_genes` |
| `autoimmune_atlas/refresh_genetics.py` | 新規。取得、datasource の判定、保存。`python -m` で実行 |
| `autoimmune_atlas/ui/figures.py` | `Kind` に `gene` を足し、ホバーと題を遺伝子用に切り替える |
| `autoimmune_atlas/ui/components.py` | `disease_selector` と `disease_checklist_sections` に ID の接頭辞、`heatmap_row_controls` にトグルの種類を足す |
| `autoimmune_atlas/ui/genetics_layout.py` | 新規。遺伝子ページのレイアウト、詳細欄、遺伝子の表 |
| `autoimmune_atlas/ui/genetics_callbacks.py` | 新規。遺伝子ページの callback |
| `autoimmune_atlas/ui/application.py` | 2 ページを包む shell と `dcc.Location` の切り替え |
| `app.py` | `genetics.json` を読み、`create_app` に渡す |
| `assets/router.js` | 新規。ページ切り替え時に `resize` を送る |
| `pixi.toml` | `refresh-genetics` タスク |
| `docs/design.md`、`README.md`、`docs/README.md` | 範囲とコード構成の更新 |
| `tests/test_genetics.py`、`tests/test_refresh_genetics.py`、`tests/test_genetics_ui.py` | 新規のテスト |
| `tests/genetics_fixture.py` | 新規。テストが共有する `GeneticsSnapshot` |

---

### Task 1: `genetics.json` の型、読み込み、発現の結合

**Files:**
- Modify: `autoimmune_atlas/models.py:192`（末尾に追加）
- Create: `autoimmune_atlas/genetics.py`
- Create: `tests/genetics_fixture.py`
- Test: `tests/test_genetics.py`

**Interfaces:**
- Consumes: `models.Snapshot`、`models.ExpressionRow`、`models.DataVersion`
- Produces:
    - `models.GeneAssociation(TypedDict)`: `target_id: str`、`target: str`、`score: float`、`datasource_scores: dict[str, float]`
    - `models.GeneticsSnapshot(TypedDict)`: `schema: int`、`data_version: DataVersion`、`retrieved_at: str`、`source: str`、`score_floor: float`、`datasources: list[str]`、`associations: dict[str, list[GeneAssociation]]`、`expression: dict[str, list[ExpressionRow]]`
    - `genetics.GENETICS_PATH: Path`（`data/genetics.json`）
    - `genetics.load_genetics(path: Path = GENETICS_PATH) -> GeneticsSnapshot | None`
    - `genetics.version_matches(snapshot: Snapshot, genetics: GeneticsSnapshot) -> bool`
    - `genetics.merged_snapshot(snapshot: Snapshot, genetics: GeneticsSnapshot) -> Snapshot`（`expression` を結合した新しい dict。同じ遺伝子は `snapshot` を優先）

- [ ] **Step 1: fixture を書く**

`tests/genetics_fixture.py`:

```python
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
```

`G1`、`G2`、`G3` の発現は `tests/core_fixture.py` の `core_snapshot()` にある。`G9` はこちらにだけある。`D2` は取得済みだが 0 件の疾患である。

- [ ] **Step 2: 失敗するテストを書く**

`tests/test_genetics.py`:

```python
"""genetics.json の読み込みと遺伝子の集計を検証する。"""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import cast

from autoimmune_atlas import genetics
from autoimmune_atlas.models import Snapshot
from tests.core_fixture import core_snapshot
from tests.genetics_fixture import genetics_snapshot


def snapshot() -> Snapshot:
    core = core_snapshot()
    return cast(
        Snapshot,
        cast(
            object,
            {
                **core,
                "root": "MONDO_0007179",
                "data_version": {"year": "26", "month": "09", "iteration": None},
                "retrieved_at": "2026-09-27T00:00:00+00:00",
                "source": "https://api.platform.opentargets.org/api/v4/graphql",
            },
        ),
    )


class LoadTests(unittest.TestCase):
    def test_missing_file_returns_none_and_bad_schema_raises(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "genetics.json"
            self.assertIsNone(genetics.load_genetics(path))
            path.write_text(json.dumps({"schema": 2}), encoding="utf-8")
            with self.assertRaises(ValueError):
                genetics.load_genetics(path)
            path.write_text(json.dumps(genetics_snapshot()), encoding="utf-8")
            loaded = genetics.load_genetics(path)
            assert loaded is not None
            self.assertEqual(loaded["score_floor"], 0.1)

    def test_version_match_and_merged_expression_prefer_snapshot(self) -> None:
        base = snapshot()
        data = genetics_snapshot()
        self.assertTrue(genetics.version_matches(base, data))
        data["data_version"] = {"year": "26", "month": "06", "iteration": None}
        self.assertFalse(genetics.version_matches(base, data))
        data = genetics_snapshot()
        data["expression"]["G1"] = [
            {"cell_id": "T4", "cell": "CD4 T cell", "median": 99.0}
        ]
        merged = genetics.merged_snapshot(base, data)
        self.assertEqual(merged["expression"]["G1"][0]["median"], 2.0)
        self.assertIn("G9", merged["expression"])
        self.assertNotIn("G9", base["expression"])
```

- [ ] **Step 3: 失敗を確かめる**

Run: `pixi run --locked --no-install python -m unittest tests.test_genetics -v`
Expected: FAIL（`ImportError: cannot import name 'genetics'` または `GeneticsSnapshot`）

- [ ] **Step 4: 型を足す**

`autoimmune_atlas/models.py` の末尾に追加:

```python
class GeneAssociation(TypedDict):
    target_id: str
    target: str
    score: float
    datasource_scores: dict[str, float]


class GeneticsSnapshot(TypedDict):
    schema: int
    data_version: DataVersion
    retrieved_at: str
    source: str
    score_floor: float
    datasources: list[str]
    associations: dict[str, list[GeneAssociation]]
    expression: dict[str, list[ExpressionRow]]
```

- [ ] **Step 5: `genetics.py` を書く（読み込みと結合まで）**

`autoimmune_atlas/genetics.py`:

```python
"""genetic association の保存データを読み、遺伝子ページ用に集計する。"""

import json
from pathlib import Path
from typing import cast

from autoimmune_atlas.models import GeneticsSnapshot, Snapshot

GENETICS_PATH = Path(__file__).resolve().parent.parent / "data" / "genetics.json"
GENETICS_SCHEMA = 1


def load_genetics(path: Path = GENETICS_PATH) -> GeneticsSnapshot | None:
    """schema 1 の genetics.json を読む。無ければ None、形が違えば ValueError。"""
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as handle:
        raw = cast(object, json.load(handle))
    if not isinstance(raw, dict):
        raise ValueError("genetics.json must contain a JSON object")
    data = cast(dict[str, object], raw)
    if data.get("schema") != GENETICS_SCHEMA:
        raise ValueError(
            "Unsupported genetics format. Refresh it with: pixi run refresh-genetics"
        )
    required = {
        "data_version",
        "retrieved_at",
        "source",
        "score_floor",
        "datasources",
        "associations",
        "expression",
    }
    missing = required - data.keys()
    if missing:
        raise ValueError(
            "Missing required fields in genetics.json: " + ", ".join(sorted(missing))
        )
    if not isinstance(data["associations"], dict) or not isinstance(
        data["expression"], dict
    ):
        raise ValueError("associations and expression must be objects")
    return cast(GeneticsSnapshot, cast(object, data))


def version_matches(snapshot: Snapshot, genetics: GeneticsSnapshot) -> bool:
    """2 つのファイルが同じ Open Targets の版から取得されたかを返す。"""
    return snapshot.get("data_version") == genetics["data_version"]


def merged_snapshot(snapshot: Snapshot, genetics: GeneticsSnapshot) -> Snapshot:
    """発現データを結合した Snapshot を返す。同じ遺伝子は snapshot.json を優先する。"""
    expression = {**genetics["expression"], **snapshot["expression"]}
    return cast(Snapshot, cast(object, {**snapshot, "expression": expression}))
```

- [ ] **Step 6: テストを通す**

Run: `pixi run --locked --no-install python -m unittest tests.test_genetics -v`
Expected: PASS（2 tests）

- [ ] **Step 7: 検査とコミット**

Run: `pixi run --locked --no-install check`
Expected: すべて成功

Tomoya の許可を得てから:

```bash
git add autoimmune_atlas/models.py autoimmune_atlas/genetics.py tests/genetics_fixture.py tests/test_genetics.py
git commit -m "feat(genetics): genetics.json の型と読み込みと発現の結合を追加"
```

---

### Task 2: 遺伝子集合の選択と `summarize_genes`

**Files:**
- Modify: `autoimmune_atlas/genetics.py`
- Test: `tests/test_genetics.py`

**Interfaces:**
- Consumes: `aggregation.expression_metadata`、`aggregation.expression_state`、`aggregation.cell_catalog`、`aggregation._aggregate`（公開版を足す。下記）、`aggregation._percent`
- Produces:
    - `genetics.genes_for_disease(genetics: GeneticsSnapshot, disease_id: str, score_threshold: float) -> list[GeneAssociation]`（スコア降順、閾値は等号を含む）
    - `genetics.summarize_genes(snapshot: Snapshot, genetics: GeneticsSnapshot, score_threshold: float, threshold: float, *, method: str = "fixed", specificity_threshold: float = 0.5, level: str = "group", cell_ids: list[str] | None = None, disease_ids: list[str] | None = None) -> list[SummaryRow]`
    - 返す `SummaryRow` は薬剤の列を次に固定する: `records=[]`、`unmapped_drugs=0`、`drug_count=None`、`drug_percent=None`、`drug_denominator=0`、`unknown_drugs=0`、`total_drugs=0`、`mapped_drugs=0`。`snapshot` は `merged_snapshot` を通したものを渡す

`aggregation.py` の `_aggregate` と `_percent` は先頭に `_` が付いている。`genetics.py` から使うために、`aggregation.py` に公開の別名を足す:

```python
aggregate_states = _aggregate
percent_of = _percent
```

（`aggregation.py:566` の `_aggregate` と `:574` の `_percent` の定義の直後に置く）

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_genetics.py` に追加:

```python
class SummarizeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.base = genetics.merged_snapshot(snapshot(), genetics_snapshot())
        self.data = genetics_snapshot()

    def test_threshold_is_inclusive_and_sorted_by_score(self) -> None:
        genes = genetics.genes_for_disease(self.data, "D1", 0.5)
        self.assertEqual([g["target_id"] for g in genes], ["G1", "G2", "G9"])
        self.assertEqual(genetics.genes_for_disease(self.data, "D2", 0.5), [])
        self.assertEqual(genetics.genes_for_disease(self.data, "MISSING", 0.5), [])

    def test_counts_share_expression_rules_and_zero_genes_give_na_percent(self):
        rows = genetics.summarize_genes(
            self.base, self.data, 0.5, 0.5, method="fixed", level="group"
        )
        lymph = next(
            r for r in rows if r["disease_id"] == "D1" and r["cell_id"] == "LYMPH"
        )
        # G1 (T4 2.0), G2 (T8 3.0), G9 (T4 5.0) は LYMPH のどこかで 0.5 以上
        self.assertEqual(lymph["count"], 3)
        self.assertEqual(lymph["denominator"], 3)
        self.assertEqual(lymph["percent"], 100.0)
        self.assertEqual(lymph["status"], "complete")
        self.assertEqual(lymph["drug_denominator"], 0)
        self.assertEqual(lymph["records"], [])
        empty = next(
            r for r in rows if r["disease_id"] == "D2" and r["cell_id"] == "LYMPH"
        )
        self.assertEqual(empty["count"], 0)
        self.assertIsNone(empty["percent"])
        self.assertEqual(empty["denominator"], 0)
        self.assertEqual(empty["status"], "complete")

    def test_gene_without_expression_is_unknown_not_error(self) -> None:
        data = genetics_snapshot()
        data["associations"]["D1"].append(
            {
                "target_id": "G_NOEXPR",
                "target": "No expr",
                "score": 0.7,
                "datasource_scores": {},
            }
        )
        rows = genetics.summarize_genes(self.base, data, 0.5, 0.5, level="group")
        lymph = next(
            r for r in rows if r["disease_id"] == "D1" and r["cell_id"] == "LYMPH"
        )
        self.assertEqual(lymph["count"], 3)
        self.assertEqual(lymph["unknown"], 1)
        self.assertEqual(lymph["denominator"], 4)
        self.assertEqual(lymph["status"], "partial")

    def test_disease_and_cell_filters_and_all_level(self) -> None:
        rows = genetics.summarize_genes(
            self.base, self.data, 0.5, 0.5, level="all", disease_ids=["D1"]
        )
        self.assertEqual(
            [(r["disease_id"], r["cell_id"]) for r in rows], [("D1", "all")]
        )
        rows = genetics.summarize_genes(
            self.base, self.data, 0.5, 0.5, level="mixed", cell_ids=["group:LYMPH"]
        )
        self.assertEqual({r["cell_id"] for r in rows}, {"group:LYMPH"})
```

- [ ] **Step 2: 失敗を確かめる**

Run: `pixi run --locked --no-install python -m unittest tests.test_genetics -v`
Expected: FAIL（`AttributeError: module 'autoimmune_atlas.genetics' has no attribute 'genes_for_disease'`）

- [ ] **Step 3: 実装する**

`autoimmune_atlas/aggregation.py` の `_percent` の定義の直後に:

```python
# genetics.py が同じ三値の集約と割合の規則を使うための公開名。
aggregate_states = _aggregate
percent_of = _percent
```

`autoimmune_atlas/genetics.py` に追加（import を `from autoimmune_atlas import aggregation as atlas` と `from autoimmune_atlas.models import CellCatalogEntry, GeneAssociation, GeneticsSnapshot, Snapshot, SummaryRow` に増やす）:

```python
def genes_for_disease(
    genetics: GeneticsSnapshot, disease_id: str, score_threshold: float
) -> list[GeneAssociation]:
    """閾値以上（等号を含む）の関連遺伝子をスコアの降順で返す。"""
    genes = [
        gene
        for gene in genetics["associations"].get(disease_id, [])
        if gene["score"] >= score_threshold
    ]
    return sorted(genes, key=lambda gene: (-gene["score"], gene["target_id"]))


def summarize_genes(
    snapshot: Snapshot,
    genetics: GeneticsSnapshot,
    score_threshold: float,
    threshold: float,
    *,
    method: str = "fixed",
    specificity_threshold: float = 0.5,
    level: str = "group",
    cell_ids: list[str] | None = None,
    disease_ids: list[str] | None = None,
) -> list[SummaryRow]:
    """疾患・細胞ごとに、閾値以上の遺伝子のうち発現判定が陽性の数と割合を返す。

    薬剤ページの標的と同じ三値判定、大分類の集約、割合の NA の規則を使う。
    薬剤に関する列は空の値で埋め、図の共有に使う。
    """
    metadata = atlas.expression_metadata(snapshot)
    catalog: list[CellCatalogEntry] = (
        [
            {
                "id": "all",
                "name": "All source cell types",
                "members": [
                    cell["id"] for cell in atlas.cell_catalog(snapshot, "cell")
                ],
            }
        ]
        if level == "all"
        else atlas.cell_catalog(snapshot, level)
    )
    if cell_ids is not None:
        selected_cells = set(cell_ids)
        catalog = [cell for cell in catalog if cell["id"] in selected_cells]
    selected_diseases = None if disease_ids is None else set(disease_ids)
    output: list[SummaryRow] = []
    for disease in snapshot.get("diseases", []):
        if selected_diseases is not None and disease["id"] not in selected_diseases:
            continue
        genes = genes_for_disease(genetics, disease["id"], score_threshold)
        targets = [gene["target_id"] for gene in genes]
        loaded = disease["id"] in genetics["associations"]
        for cell in catalog:
            states = {
                target: atlas.aggregate_states(
                    [
                        atlas.expression_state(
                            metadata.get((target, member)),
                            threshold,
                            method,
                            specificity_threshold,
                        )
                        for member in cell["members"]
                    ]
                )
                for target in targets
            }
            positive = {t for t, s in states.items() if s is True}
            unknown = {t for t, s in states.items() if s is None}
            output.append(
                {
                    "disease_id": disease["id"],
                    "disease": disease["name"],
                    "cell_id": cell["id"],
                    "cell": cell["name"],
                    "ontology_id": cell.get("ontology_id", cell["id"]),
                    "cell_level": cell.get("cell_level", level),
                    "count": len(positive) if loaded else None,
                    "percent": atlas.percent_of(positive, unknown, len(targets))
                    if loaded
                    else None,
                    "denominator": len(targets),
                    "unknown": len(unknown),
                    "unmapped_drugs": 0,
                    "status": "unavailable"
                    if not loaded
                    else "partial"
                    if unknown
                    else "complete",
                    "records": [],
                    "drug_count": None,
                    "drug_percent": None,
                    "drug_denominator": 0,
                    "unknown_drugs": 0,
                    "total_drugs": 0,
                    "mapped_drugs": 0,
                    "member_cell_ids": cell["members"],
                }
            )
    return output
```

`aggregation.expression_state` は公開関数で、閾値の検証をしてから `_expression_state` を呼ぶ。1 セルごとに検証が走るが、対象は数百遺伝子 × 170 細胞なので許容する。遅ければ `_expression_state` を公開名にする。

- [ ] **Step 4: テストを通す**

Run: `pixi run --locked --no-install python -m unittest tests.test_genetics -v`
Expected: PASS（6 tests）

- [ ] **Step 5: 検査とコミット**

Run: `pixi run --locked --no-install check`

Tomoya の許可を得てから:

```bash
git add autoimmune_atlas/aggregation.py autoimmune_atlas/genetics.py tests/test_genetics.py
git commit -m "feat(genetics): 閾値以上の関連遺伝子を三値判定で集計する summarize_genes を追加"
```

---

### Task 3: `refresh_genetics.py` と `pixi run refresh-genetics`

**Files:**
- Create: `autoimmune_atlas/refresh_genetics.py`
- Modify: `pixi.toml:11`（`refresh` の次の行）
- Test: `tests/test_refresh_genetics.py`

**Interfaces:**
- Consumes: `refresh.query_api`、`refresh.fetch_expression`、`refresh.save_snapshot`、`refresh.API_HOST`、`refresh.API_PATH`、`snapshot.load_snapshot`、`genetics.GENETICS_PATH`、`genetics.version_matches`
- Produces:
    - `refresh_genetics.SCORE_FLOOR = 0.1`、`PAGE_SIZE = 500`
    - `refresh_genetics.collect_associations(query, disease_ids: list[str], floor: float = SCORE_FLOOR) -> dict[str, list[RawAssociation]]`（`RawAssociation` は `target_id`、`target`、`score`、`datasource_scores` を持ち、`datasource_scores` はまだ全 datasource を含む）
    - `refresh_genetics.classify_datasources(query, associations) -> list[str]`（genetic_association に属する datasource ID、名前順）
    - `refresh_genetics.restrict_datasources(associations, datasources) -> dict[str, list[GeneAssociation]]`
    - `refresh_genetics.main() -> None`

型は `refresh.py` と同じく `TypedDict` で書く。`query` の型は `Callable[[str, Mapping[str, object]], dict[str, object]]`。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_refresh_genetics.py`:

```python
"""genetic association の取得と datasource の判定を検証する。"""

import unittest
from collections.abc import Mapping

from autoimmune_atlas.refresh_genetics import (
    classify_datasources,
    collect_associations,
    restrict_datasources,
)


def _row(
    target: str, genetic: float | None, extra: dict[str, float]
) -> dict[str, object]:
    datatype = (
        [] if genetic is None else [{"id": "genetic_association", "score": genetic}]
    )
    datasource = [{"id": key, "score": value} for key, value in extra.items()]
    return {
        "target": {"id": target, "approvedSymbol": target.lower()},
        "datatypeScores": datatype + [{"id": "literature", "score": 0.3}],
        "datasourceScores": datasource,
    }


PAGES: dict[tuple[str, int], list[dict[str, object]]] = {
    ("D1", 0): [
        _row("G1", 0.9, {"gwas_credible_sets": 0.9, "europepmc": 0.3}),
        _row("G2", 0.5, {"eva": 0.5}),
    ],
    ("D1", 1): [
        _row("G3", 0.05, {"eva": 0.05}),
        _row("G4", None, {"europepmc": 0.3}),
    ],
    ("D2", 0): [_row("G5", None, {"europepmc": 0.3})],
}


def fake_query(query: str, variables: Mapping[str, object]) -> dict[str, object]:
    if "associatedTargets" in query:
        disease = str(variables["id"])
        index = int(str(variables["index"]))
        total = sum(len(v) for (d, _), v in PAGES.items() if d == disease)
        return {
            "disease": {
                "associatedTargets": {
                    "count": total,
                    "rows": PAGES.get((disease, index), []),
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


class RefreshGeneticsTests(unittest.TestCase):
    def test_paging_stops_below_floor_and_keeps_scored_genes(self) -> None:
        associations = collect_associations(fake_query, ["D1", "D2"], 0.1)
        self.assertEqual([g["target_id"] for g in associations["D1"]], ["G1", "G2"])
        self.assertEqual(associations["D2"], [])
        self.assertEqual(associations["D1"][0]["target"], "g1")

    def test_page_size_is_two_in_this_test_only_via_patch(self) -> None:
        # PAGE_SIZE を 2 にして 2 ページ目まで読むことを確かめる。
        from unittest.mock import patch

        with patch("autoimmune_atlas.refresh_genetics.PAGE_SIZE", 2):
            associations = collect_associations(fake_query, ["D1"], 0.1)
        self.assertEqual(len(associations["D1"]), 2)

    def test_datasources_are_classified_by_one_evidence_each(self) -> None:
        associations = collect_associations(fake_query, ["D1"], 0.1)
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
            collect_associations(missing, ["D1"], 0.1)
```

`fake_query` の 1 ページ目に 2 行、2 ページ目に 2 行を置いてあるので、`PAGE_SIZE` が 500 のままだと `count` 4 に対して 1 ページ目で `2 < 4` となり 2 ページ目を読む。テスト 1 はそれで通る。テスト 2 は `PAGE_SIZE` を 2 にしても同じ結果になることを確かめるもので、ページングの終了条件（`index * PAGE_SIZE >= count`）に依存しないことを見る。

- [ ] **Step 2: 失敗を確かめる**

Run: `pixi run --locked --no-install python -m unittest tests.test_refresh_genetics -v`
Expected: FAIL（`ModuleNotFoundError: No module named 'autoimmune_atlas.refresh_genetics'`）

- [ ] **Step 3: 実装する**

`autoimmune_atlas/refresh_genetics.py`:

```python
"""Open Targets の genetic association から遺伝子ページ用のデータを作る。"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from typing import TypedDict, cast

from autoimmune_atlas.genetics import GENETICS_PATH, GENETICS_SCHEMA
from autoimmune_atlas.models import (
    DataVersion,
    ExpressionRow,
    GeneAssociation,
    GeneticsSnapshot,
    Snapshot,
)
from autoimmune_atlas.refresh import (
    API_HOST,
    API_PATH,
    fetch_expression,
    query_api,
    save_snapshot,
)
from autoimmune_atlas.snapshot import load_snapshot

SCORE_FLOOR = 0.1
PAGE_SIZE = 500
GENETIC_DATATYPE = "genetic_association"

Query = Callable[[str, Mapping[str, object]], dict[str, object]]


class _Score(TypedDict):
    id: str
    score: float


class _Target(TypedDict):
    id: str
    approvedSymbol: str


class _AssociationRow(TypedDict):
    target: _Target
    datatypeScores: list[_Score]
    datasourceScores: list[_Score]


class _AssociationBlock(TypedDict):
    count: int
    rows: list[_AssociationRow]


class _AssociationDisease(TypedDict):
    associatedTargets: _AssociationBlock


class _AssociationResponse(TypedDict):
    disease: _AssociationDisease | None


class _EvidenceRow(TypedDict):
    datasourceId: str
    datatypeId: str


class _EvidenceBlock(TypedDict):
    rows: list[_EvidenceRow]


class _EvidenceDisease(TypedDict):
    evidences: _EvidenceBlock


class _EvidenceResponse(TypedDict):
    disease: _EvidenceDisease | None


class _Meta(TypedDict):
    dataVersion: DataVersion


class _VersionResponse(TypedDict):
    meta: _Meta


ASSOCIATION_QUERY = """query($id:String!,$index:Int!,$size:Int!){
  disease(efoId:$id){ associatedTargets(
    page:{index:$index,size:$size}, enableIndirect:false, orderByScore:"genetic_association"
  ){ count rows{ target{id approvedSymbol} datatypeScores{id score} datasourceScores{id score} } } }
}"""

EVIDENCE_QUERY = """query($id:String!,$gene:String!,$datasource:String!){
  disease(efoId:$id){ evidences(ensemblIds:[$gene], datasourceIds:[$datasource], size:1){
    rows{ datasourceId datatypeId } } }
}"""

VERSION_QUERY = "query{meta{dataVersion{year month iteration}}}"


def collect_associations(
    query: Query, disease_ids: list[str], floor: float = SCORE_FLOOR
) -> dict[str, list[GeneAssociation]]:
    """疾患ごとに genetic スコアの降順で取り、floor を下回った時点で止める。"""
    output: dict[str, list[GeneAssociation]] = {}
    for disease_id in disease_ids:
        genes: list[GeneAssociation] = []
        index = 0
        while True:
            response = cast(
                _AssociationResponse,
                cast(
                    object,
                    query(
                        ASSOCIATION_QUERY,
                        {"id": disease_id, "index": index, "size": PAGE_SIZE},
                    ),
                ),
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
                if genetic is None or genetic < floor:
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
            if stop or not block["rows"] or index * PAGE_SIZE >= block["count"]:
                break
        output[disease_id] = genes
    return output


def classify_datasources(
    query: Query, associations: Mapping[str, list[GeneAssociation]]
) -> list[str]:
    """現れた datasource ごとに根拠を 1 件取り、genetic_association のものだけ残す。"""
    sample: dict[str, tuple[str, str]] = {}
    for disease_id, genes in associations.items():
        for gene in genes:
            for datasource in gene["datasource_scores"]:
                sample.setdefault(datasource, (disease_id, gene["target_id"]))
    genetic: list[str] = []
    for datasource, (disease_id, gene_id) in sorted(sample.items()):
        response = cast(
            _EvidenceResponse,
            cast(
                object,
                query(
                    EVIDENCE_QUERY,
                    {"id": disease_id, "gene": gene_id, "datasource": datasource},
                ),
            ),
        )
        disease = response["disease"]
        rows = disease["evidences"]["rows"] if disease else []
        if not rows:
            raise ValueError(f"datasource の根拠が取れません: {datasource}")
        if rows[0]["datatypeId"] == GENETIC_DATATYPE:
            genetic.append(datasource)
    return genetic


def restrict_datasources(
    associations: Mapping[str, list[GeneAssociation]], datasources: list[str]
) -> dict[str, list[GeneAssociation]]:
    """datasource ごとのスコアを genetic association のものに絞る。"""
    keep = set(datasources)
    return {
        disease_id: [
            {
                **gene,
                "datasource_scores": {
                    key: value
                    for key, value in gene["datasource_scores"].items()
                    if key in keep
                },
            }
            for gene in genes
        ]
        for disease_id, genes in associations.items()
    }


def main() -> None:
    """snapshot.json の疾患について関連遺伝子と不足分の発現を取り、genetics.json を書く。"""
    snapshot: Snapshot | None = load_snapshot()
    if snapshot is None:
        raise SystemExit(
            "snapshot.json がありません。先に pixi run refresh を実行してください"
        )
    version = cast(_VersionResponse, cast(object, query_api(VERSION_QUERY)))["meta"][
        "dataVersion"
    ]
    if snapshot.get("data_version") != version:
        raise SystemExit(
            "snapshot.json の版が API と違います。先に pixi run refresh を実行してください"
        )
    disease_ids = [disease["id"] for disease in snapshot["diseases"]]
    associations = collect_associations(query_api, disease_ids)
    print(
        f"疾患: {len(disease_ids)} / 関連遺伝子を持つ疾患: {sum(bool(g) for g in associations.values())}",
        flush=True,
    )
    datasources = classify_datasources(query_api, associations)
    print(f"genetic association の datasource: {', '.join(datasources)}", flush=True)
    associations = restrict_datasources(associations, datasources)
    needed = sorted(
        {gene["target_id"] for genes in associations.values() for gene in genes}
        - set(snapshot["expression"])
    )
    expression: dict[str, list[ExpressionRow]] = {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        for index, (target, rows) in enumerate(pool.map(fetch_expression, needed), 1):
            expression[target] = rows
            print(f"細胞型別発現: {index}/{len(needed)}", flush=True)
    final_version = cast(_VersionResponse, cast(object, query_api(VERSION_QUERY)))[
        "meta"
    ]["dataVersion"]
    if final_version != version:
        raise ValueError("取得中にデータの版が変わりました。再取得してください")
    genetics: GeneticsSnapshot = {
        "schema": GENETICS_SCHEMA,
        "data_version": version,
        "retrieved_at": datetime.now(UTC).isoformat(),
        "source": f"https://{API_HOST}{API_PATH}",
        "score_floor": SCORE_FLOOR,
        "datasources": datasources,
        "associations": associations,
        "expression": expression,
    }
    save_snapshot(GENETICS_PATH, genetics)
    print(
        f"保存: {GENETICS_PATH}\n疾患 {len(disease_ids)} / 遺伝子 {len(needed) + len(set(snapshot['expression']) & {g['target_id'] for gs in associations.values() for g in gs})}",
        flush=True,
    )


if __name__ == "__main__":
    main()
```

`pixi.toml` の `refresh = ...` の次の行に:

```toml
refresh-genetics = "python -m autoimmune_atlas.refresh_genetics"
```

- [ ] **Step 4: テストを通す**

Run: `pixi run --locked --no-install python -m unittest tests.test_refresh_genetics -v`
Expected: PASS（4 tests）

- [ ] **Step 5: 検査とコミット**

Run: `pixi run --locked --no-install check`

Tomoya の許可を得てから:

```bash
git add autoimmune_atlas/refresh_genetics.py pixi.toml tests/test_refresh_genetics.py
git commit -m "feat(genetics): genetic association と不足分の発現を取得する refresh-genetics を追加"
```

- [ ] **Step 6: Tomoya に実データの取得を頼む**

`pixi run refresh-genetics` は Tomoya が自分の端末で打つ。数千件の発現の取得があるので、`refresh` の発現の取得の 5 倍程度を目安に伝える。完了したら `data/genetics.json` の大きさと、出力の「genetic association の datasource」の行を報告してもらう。

---

### Task 4: 図と部品を遺伝子ページで共有できる形にする

**Files:**
- Modify: `autoimmune_atlas/ui/figures.py:39`（`Kind`）、`:170-234`（ホバー）、`:295`、`:439`（題）
- Modify: `autoimmune_atlas/ui/components.py:125-172`（`disease_checklist_sections`）、`:260-309`（`disease_selector`）、`:312-340`（`heatmap_row_controls`）
- Test: `tests/test_genetics_ui.py`

**Interfaces:**
- Produces:
    - `figures.Kind = Literal["target", "drug", "gene"]`。`kind == "gene"` のとき、ホバーの見出しは "Genes"、薬剤の網羅性の行を出さない、`count == 0` かつ `denominator == 0` のとき "No genes at or above the score threshold"、色の題は "Genes"
    - `components.disease_checklist_sections(family, prefix: str = "")`。`prefix` を `disease-family-` と `disease-details-` の前に付ける
    - `components.disease_selector(snapshot, selected, prefix: str = "")`。`Dropdown` の ID を `f"{prefix}diseases"`、`Label` の `htmlFor` も同じ、`info_tip` のキーを `f"{prefix}diseases"`
    - `components.heatmap_row_controls(names, cell_ids, expanded, kind)`。`kind` が `"target"`、`"drug"` なら type `heatmap-cell-toggle`、`"expression"` なら `expression-cell-toggle`、それ以外は `f"{kind}-cell-toggle"`

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_genetics_ui.py`:

```python
"""遺伝子ページの図、部品、レイアウト、callback を検証する。"""

from __future__ import annotations

import unittest
from typing import cast

from autoimmune_atlas.models import DiseaseFamily, SummaryRow
from autoimmune_atlas.ui import components, figures


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
        heatmap = figure.data[0]
        self.assertEqual(heatmap.colorbar.title.text, "Genes (count)")  # pyright: ignore[reportAttributeAccessIssue, reportUnknownMemberType]


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
        self.assertEqual(button.id["type"], "genetics-cell-toggle")  # pyright: ignore[reportAttributeAccessIssue, reportIndexIssue, reportUnknownMemberType]
        button = components.heatmap_row_controls(names, ["group:LYMPH"], [], "target")[
            0
        ]
        self.assertEqual(button.id["type"], "heatmap-cell-toggle")  # pyright: ignore[reportAttributeAccessIssue, reportIndexIssue, reportUnknownMemberType]
```

- [ ] **Step 2: 失敗を確かめる**

Run: `pixi run --locked --no-install python -m unittest tests.test_genetics_ui -v`
Expected: FAIL（`hover_text` に "Canonical drugs" が含まれる、`disease_checklist_sections` が `prefix` を受け付けない）

- [ ] **Step 3: `figures.py` を直す**

`figures.py:39`:

```python
type Kind = Literal["target", "drug", "gene"]
```

`hover_text`（`:170`）を次に置き換える:

```python
def _kind_label(kind: Kind) -> str:
    return {"drug": "Canonical drugs", "gene": "Genes"}.get(kind, "Targets")


def _coverage_lines(row: SummaryRow, kind: Kind) -> list[str]:
    """薬剤の網羅性の行。遺伝子には網羅性の概念が無いので出さない。"""
    if kind == "gene":
        return []
    return [
        f"Canonical drugs with targets / all canonical drugs: {row.get('mapped_drugs', 0)} / {row.get('total_drugs', 0)}"
    ]


def hover_text(row: SummaryRow, measure: Measure, kind: Kind) -> str:
    value_key, unknown_key, denominator_key = measure_fields(kind, measure)
    value = display_value(row, measure, kind) or "missing"
    label = _kind_label(kind)
    if row["status"] == "unavailable":
        assessment = "Disease data not loaded"
    elif kind == "gene" and row.get(denominator_key, 0) == 0:
        assessment = "No genes at or above the score threshold"
    elif measure == "percent" and row.get(denominator_key, 0) == 0:
        assessment = "No eligible items in the percentage denominator"
    elif _is_lower_bound(row, measure, kind):
        assessment = "≥ is a lower bound; unresolved evidence may increase this value"
    elif row.get(value_key) is None:
        assessment = "The expression rule is unresolved because target or expression data are missing."
    elif row[value_key] == 0:
        assessment = f"No qualifying {label.lower()} under this rule"
    else:
        assessment = f"{label} meeting the expression rule"
    return "<br>".join(
        (
            f"<b>{row['disease']} / {row['cell']}</b>",
            f"{label}: {value}",
            assessment,
            *(
                ["Displayed as zero for readability; not a measured zero."]
                if row.get(value_key) is None
                else []
            ),
            f"Denominator: {row.get(denominator_key, 0)}",
            f"Unresolved in denominator: {row.get(unknown_key, 0)}",
            *_coverage_lines(row, kind),
        )
    )
```

`dot_hover_text`（`:203`）も同じ形に: `label = _kind_label(kind)`、`row["status"] == "unavailable"` の次に `elif kind == "gene" and row.get(denominator_key, 0) == 0: assessment = "No genes at or above the score threshold"` を挟み、末尾の網羅性の行を `*_coverage_lines(row, kind)` に置き換える。

`build_figure`（`:295`）と `build_dot_figure`（`:439`）の `title = "Canonical drugs" if kind == "drug" else "Targets"` を `title = _kind_label(kind)` に置き換える。

`measure_fields` は `kind == "drug"` 以外を標的の列として扱うので、そのままでよい。

- [ ] **Step 4: `components.py` を直す**

`disease_checklist_sections(family: DiseaseFamily, prefix: str = "")`。`section_id` と先頭の section の ID に `prefix` を前置する:

```python
    def section_id(kind: str, node: str) -> str:
        suffix = "" if node == root else f"-{node}"
        return f"{prefix}{kind}-{family['id']}{suffix}"

    sections: list[DiseaseSection] = [
        {
            "id": f"{prefix}disease-family-{family['id']}",
            ...
```

`_disease_family_selector(family, selected_ids, prefix: str = "")` に `prefix` を通し、`sections = disease_checklist_sections(family, prefix)` にする。

`disease_selector(snapshot, selected, prefix: str = "")`:

```python
            choices, component = _disease_family_selector(family, selected_ids, prefix)
    ...
                    html.Label("Diseases", htmlFor=f"{prefix}diseases"),
                    info_tip(
                        f"{prefix}diseases",
                        "diseases",
                        "Diseases sets which diseases ...",  # 文言は既存のまま
                    ),
    ...
            dcc.Dropdown(
                id=f"{prefix}diseases",
```

`heatmap_row_controls` の `"type"`:

```python
                "type": "heatmap-cell-toggle"
                if kind in ("target", "drug")
                else f"{kind}-cell-toggle",
```

`kind == "expression"` は `expression-cell-toggle` になり、既存と同じである。

- [ ] **Step 5: テストを通す**

Run: `pixi run --locked --no-install python -m unittest tests.test_genetics_ui tests.test_ui -v`
Expected: PASS。`test_ui.py` の既存テストも変わらず通る

- [ ] **Step 6: 検査とコミット**

Run: `pixi run --locked --no-install check`

Tomoya の許可を得てから:

```bash
git add autoimmune_atlas/ui/figures.py autoimmune_atlas/ui/components.py tests/test_genetics_ui.py
git commit -m "refactor(ui): 図の kind に gene を足し疾患選択欄と行見出しに ID の接頭辞を追加"
```

---

### Task 5: 遺伝子ページのレイアウト、詳細欄、遺伝子の表

**Files:**
- Create: `autoimmune_atlas/ui/genetics_layout.py`
- Modify: `autoimmune_atlas/ui/config.py`（定数を足す）
- Test: `tests/test_genetics_ui.py`

**Interfaces:**
- Consumes: `layout.choose_defaults`、`layout.format_data_version`、`components.disease_selector(prefix="genetics-")`、`components.info_tip`、`config.METHOD_LABELS`、`genetics.genes_for_disease`、`genetics.GeneticsSnapshot`
- Produces:
    - `config.DEFAULT_SCORE_THRESHOLD = 0.5`、`config.SCORE_FLOOR = 0.1`、`config.GENETICS_PARAMETER_IDS = ("genetics-measure", "genetics-score", "genetics-method", "genetics-threshold", "genetics-specificity", "genetics-diseases")`
    - `genetics_layout.GENE_COLUMNS`、`genetics_layout.gene_rows(genes: Sequence[GeneAssociation], disease_id: str, datasources: Sequence[str]) -> list[dict[str, str]]`、`genetics_layout.gene_grid(rows, datasources) -> dag.AgGrid`（ID `genetics-gene-grid`）
    - `genetics_layout.genetics_detail_panel(rows: Sequence[SummaryRow], selection: str | None, genetics: GeneticsSnapshot, *, score_threshold: float, threshold: float, method: str, specificity: float) -> html.Div`。`dcc.Store(id="genetics-source-context")` に `disease_id`、`score`、`threshold`、`method`、`specificity` を入れる。発現図の `dcc.Graph` は `genetics-expression-heatmap`、行見出しは `genetics-expression-row-controls`、図の切り替えは `genetics-expression-chart-type`、凡例は `genetics-expression-heatmap-key` と `genetics-expression-dot-key`
    - `genetics_layout.genetics_page(snapshot: Snapshot, genetics: GeneticsSnapshot) -> html.Main`（ID `genetics-page`、クラス `shell`）
    - `genetics_layout.genetics_unavailable_page(reason: str) -> html.Main`（ID `genetics-page`。案内文と `pixi run refresh-genetics` を出す）
    - `genetics_layout.page_nav(current: Literal["drugs", "genetics"]) -> html.Nav`。`Drug targets`（`/`）と `Genetic associations`（`/genetics`）のリンクを持ち、`current` に `aria-current="page"` を付ける。Task 7 で薬剤ページの `Nav` もこれに置き換える

上部パネルの部品 ID は薬剤ページの ID に `genetics-` を付けたもの: `genetics-method`、`genetics-threshold`、`genetics-specificity`、`genetics-measure`、`genetics-update-button`、`genetics-update-status`、`genetics-applied-parameters`、`genetics-expanded-cell-groups`、`genetics-expanded-expression-groups`、`genetics-chart-type`、`genetics-matrix-note`、`genetics-heatmap`、`genetics-row-controls`、`genetics-dot-key`、`genetics-detail-disease`、`genetics-details`。新しい設定は `genetics-score`（`dcc.Input`、`type="number"`、`min=SCORE_FLOOR`、`max=1`、`step=0.05`、`value=DEFAULT_SCORE_THRESHOLD`）。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_genetics_ui.py` に追加（import に `from autoimmune_atlas import genetics` と `from autoimmune_atlas.ui import genetics_layout`、`from tests.genetics_fixture import genetics_snapshot`、`from tests.test_genetics import snapshot as core_ui_snapshot` を足す）:

```python
class GeneticsLayoutTests(unittest.TestCase):
    def setUp(self) -> None:
        self.snapshot = core_ui_snapshot()
        self.genetics = genetics_snapshot()

    def _ids(self, node: object, found: set[str]) -> None:
        component_id = getattr(node, "id", None)
        if isinstance(component_id, str):
            found.add(component_id)
        for child in getattr(node, "children", None) or []:
            self._ids(child, found)
        children = getattr(node, "children", None)
        if children is not None and not isinstance(children, list):
            self._ids(children, found)

    def test_page_has_prefixed_controls_and_nav_marks_current(self) -> None:
        page = genetics_layout.genetics_page(self.snapshot, self.genetics)
        found: set[str] = set()
        self._ids(page, found)
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

    def test_score_input_range_and_default(self) -> None:
        page = genetics_layout.genetics_page(self.snapshot, self.genetics)
        score = next(
            c for c in _walk(page) if getattr(c, "id", None) == "genetics-score"
        )
        self.assertEqual(score.min, 0.1)  # pyright: ignore[reportAttributeAccessIssue, reportUnknownMemberType]
        self.assertEqual(score.max, 1)  # pyright: ignore[reportAttributeAccessIssue, reportUnknownMemberType]
        self.assertEqual(score.value, 0.5)  # pyright: ignore[reportAttributeAccessIssue, reportUnknownMemberType]

    def test_gene_rows_have_scores_datasource_columns_and_evidence_link(self) -> None:
        genes = genetics.genes_for_disease(self.genetics, "D1", 0.5)
        rows = genetics_layout.gene_rows(genes, "D1", self.genetics["datasources"])
        self.assertEqual(rows[0]["gene"], "Gene 1 (G1)")
        self.assertEqual(rows[0]["score"], "0.900")
        self.assertEqual(rows[0]["gwas_credible_sets"], "0.900")
        self.assertEqual(rows[0]["eva"], "—")
        self.assertIn(
            "https://platform.opentargets.org/evidence/G1/D1", rows[0]["links"]
        )
        grid = genetics_layout.gene_grid(rows, self.genetics["datasources"])
        fields = [c["field"] for c in grid.columnDefs]  # pyright: ignore[reportAttributeAccessIssue, reportUnknownMemberType, reportUnknownVariableType]
        self.assertEqual(
            fields, ["gene", "score", "gwas_credible_sets", "eva", "links"]
        )

    def test_detail_panel_without_selection_and_unavailable_page(self) -> None:
        panel = genetics_layout.genetics_detail_panel(
            [],
            None,
            self.genetics,
            score_threshold=0.5,
            threshold=0.5,
            method="fixed",
            specificity=0.5,
        )
        self.assertEqual(panel.className, "empty-note")  # pyright: ignore[reportAttributeAccessIssue, reportUnknownMemberType]
        page = genetics_layout.genetics_unavailable_page("versions differ")
        self.assertIn("refresh-genetics", str(page))


def _walk(node: object):
    yield node
    children = getattr(node, "children", None)
    if isinstance(children, list):
        for child in children:  # pyright: ignore[reportUnknownVariableType]
            yield from _walk(child)
    elif children is not None:
        yield from _walk(children)
```

`GeneticsLayoutTests._ids` は `_walk` を使って書き直してよい（`found = {c.id for c in _walk(page) if isinstance(getattr(c, "id", None), str)}`）。どちらかに揃える。

- [ ] **Step 2: 失敗を確かめる**

Run: `pixi run --locked --no-install python -m unittest tests.test_genetics_ui -v`
Expected: FAIL（`ImportError: cannot import name 'genetics_layout'`）

- [ ] **Step 3: `config.py` に定数を足す**

```python
DEFAULT_SCORE_THRESHOLD = 0.5
SCORE_FLOOR = 0.1
GENETICS_PARAMETER_IDS = (
    "genetics-measure",
    "genetics-score",
    "genetics-method",
    "genetics-threshold",
    "genetics-specificity",
    "genetics-diseases",
)
```

- [ ] **Step 4: `genetics_layout.py` を書く**

薬剤ページの `layout.py:245-805` を手本にし、次の違いで組む。

```python
"""遺伝子ページの画面構成。"""

from collections.abc import Sequence
from datetime import datetime, timedelta, timezone
from typing import Literal

import dash_ag_grid as dag  # pyright: ignore[reportMissingTypeStubs] - dash-ag-grid に型スタブがない。
from dash import dcc, html

from autoimmune_atlas import genetics as gene_data
from autoimmune_atlas.disease_catalog import ordered_disease_ids
from autoimmune_atlas.models import (
    GeneAssociation,
    GeneticsSnapshot,
    Snapshot,
    SummaryRow,
)
from autoimmune_atlas.ui.components import disease_selector, info_tip
from autoimmune_atlas.ui.config import (
    DEFAULT_EXPRESSION_THRESHOLD,
    DEFAULT_SCORE_THRESHOLD,
    DEFAULT_SPECIFICITY_THRESHOLD,
    GENETICS_PARAMETER_IDS,
    METHOD_LABELS,
    SCORE_FLOOR,
    SOURCE_PAGE_SIZE,
)
from autoimmune_atlas.ui.layout import choose_defaults, format_data_version

OPEN_TARGETS = "https://platform.opentargets.org/"
DEFAULT_DISEASE_TERMS = (
    "systemic lupus erythematosus",
    "systemic sclerosis",
    "Sjogren syndrome",
    "rheumatoid arthritis",
    "multiple sclerosis",
    "dermatomyositis",
    "type 1 diabetes mellitus",
    "anti-neutrophil cytoplasmic antibody-associated vasculitis",
    "pemphigus",
    "autoimmune hepatitis",
)


def page_nav(current: Literal["drugs", "genetics"]) -> html.Nav:
    """2 ページへのリンクと Open Targets へのリンクを持つ上部バー。"""

    def link(label: str, href: str, key: str) -> html.A:
        extra = {"aria-current": "page"} if key == current else {}
        return html.A(
            label,
            href=href,
            className="app-page-link",
            **extra,  # pyright: ignore[reportArgumentType] - Dash の型定義に ARIA kwargs がない。
        )

    return html.Nav(
        [
            html.A("Autoimmune Atlas", href="/", className="app-brand"),
            html.Div(
                [
                    link("Drug targets", "/", "drugs"),
                    link("Genetic associations", "/genetics", "genetics"),
                    html.A(
                        "Open Targets",
                        href=OPEN_TARGETS,
                        target="_blank",
                        rel="noreferrer",
                    ),
                ],
                className="app-nav",
            ),
        ],
        className="app-bar",
        **{  # pyright: ignore[reportArgumentType] - Dash の型定義に ARIA kwargs がない。
            "aria-label": "Main navigation"
        },
    )
```

`DEFAULT_DISEASE_TERMS` は `layout.py:250-261` の 10 疾患をここへ移し、`layout.py` からも参照する（重複を避ける）。

`gene_rows` と `gene_grid`:

```python
GENE_COLUMNS: tuple[tuple[str, str, str | None], ...] = (
    ("gene", "Gene", None),
    (
        "score",
        "Genetic association score",
        "The Open Targets datatype score for genetic association, combined across the genetic datasources. It is not a probability of causality or of treatment benefit.",
    ),
    ("links", "Open Targets", None),
)


def gene_rows(
    genes: Sequence[GeneAssociation], disease_id: str, datasources: Sequence[str]
) -> list[dict[str, str]]:
    """閾値以上の遺伝子 1 件を 1 行にし、datasource ごとのスコアを列にする。"""
    rows: list[dict[str, str]] = []
    for gene in genes:
        row = {
            "gene": f"{gene['target']} ({gene['target_id']})",
            "score": f"{gene['score']:.3f}",
            "links": f"[Evidence](https://platform.opentargets.org/evidence/{gene['target_id']}/{disease_id})",
        }
        for datasource in datasources:
            value = gene["datasource_scores"].get(datasource)
            row[datasource] = "—" if value is None else f"{value:.3f}"
        rows.append(row)
    return rows


def gene_grid(rows: Sequence[dict[str, str]], datasources: Sequence[str]) -> dag.AgGrid:
    """列ごとにソートと絞り込みができる、関連遺伝子の表。"""
    column_defs: list[dict[str, object]] = []
    for field, header, help_text in GENE_COLUMNS[:2]:
        column: dict[str, object] = {"field": field, "headerName": header}
        if help_text:
            column["headerTooltip"] = help_text
        if field == "score":
            column["sort"] = "desc"
        column_defs.append(column)
    for datasource in datasources:
        column_defs.append(
            {
                "field": datasource,
                "headerName": datasource.replace("_", " "),
                "headerTooltip": f"Datasource score from {datasource}. — means no evidence from this datasource.",
            }
        )
    column_defs.append(
        {
            "field": "links",
            "headerName": "Open Targets",
            "cellRenderer": "markdown",
            "linkTarget": "_blank",
            "sortable": False,
            "filter": False,
            "floatingFilter": False,
        }
    )
    return dag.AgGrid(
        id="genetics-gene-grid",
        columnDefs=column_defs,
        rowData=list(rows),
        defaultColDef={
            "sortable": True,
            "filter": "agTextColumnFilter",
            "floatingFilter": True,
            "resizable": True,
            "wrapText": True,
            "autoHeight": True,
            "minWidth": 120,
        },
        columnSize="responsiveSizeToFit",
        dashGridOptions={
            "pagination": True,
            "paginationPageSize": SOURCE_PAGE_SIZE,
            "paginationPageSizeSelector": False,
            "animateRows": False,
            "domLayout": "autoHeight",
            "tooltipShowDelay": 300,
            "suppressCellFocus": True,
        },
        className="ag-theme-alpine source-records-grid",
    )
```

`genetics_detail_panel` は `layout.detail_panel`（`layout.py:73-217`）と同じ構造で、次を変える。

- `selection is None` と、`rows` に疾患が無いときの `empty-note` は同じ文言
- `genes = gene_data.genes_for_disease(genetics, row["disease_id"], score_threshold)`
- 表の見出しは `f"Genetically associated genes for {row['disease']}"`、ⓘ は "This table lists the genes whose genetic association score for the selected disease is at or above the applied threshold. Each row is one gene. Datasource columns show the score from each genetic datasource; — means no evidence from that datasource. Expression rules do not filter this table. Ten rows are shown per page. Click a column header to sort, or type in the boxes under the headers to filter. Evidence opens the Open Targets evidence page for the gene and disease."
- `dcc.Store(id="genetics-source-context", data={"disease_id": ..., "score": score_threshold, "threshold": threshold, "method": method, "specificity": specificity})`
- `genes` が空なら表の代わりに `html.P("No genes at or above the score threshold for this disease.", className="empty-note")`
- 発現の見出しは "Relative gene expression by cell type"、ⓘ の文言は薬剤ページの "every target of the selected disease's filtered canonical drugs" を "every gene at or above the applied genetic association score" に置き換え、残りは同じ
- 発現図の ID は `genetics-expression-heatmap`、`genetics-expression-row-controls`、`genetics-expression-chart-type`（`persistence=True`、`persistence_type="session"`）、`genetics-expression-heatmap-key`、`genetics-expression-dot-key`

`genetics_page(snapshot, genetics)` は `layout.dashboard_layout` と同じ骨格で、次を変える。

- `page_nav("genetics")`
- eyebrow "GENETIC ASSOCIATIONS · HEALTHY REFERENCE EXPRESSION"、H1 "Genetically associated genes across diseases and cell types"、ⓘ "This atlas compares which genetically associated genes meet an expression rule in healthy reference cell types across diseases. Differences between diseases reflect their associated gene sets, not disease-specific expression. The genetic association score is combined by Open Targets from several genetic datasources and does not establish causality or treatment efficacy. The snapshot totals above the settings describe all saved data, not the filtered results."
- 概要欄: "Disease terms" `f"{loaded} loaded / {len(snapshot['diseases'])} total"`（`loaded = sum(d["id"] in genetics["associations"] for d in snapshot["diseases"])`）、"Gene–disease associations (score ≥ 0.1)" `sum(len(g) for g in genetics["associations"].values())`、"Genes with expression" `len(set(snapshot["expression"]) | set(genetics["expression"]))`（`snapshot` は結合前でも後でも同じ値）、"Retrieved at"（`genetics["retrieved_at"]` を JST に）、"Expression reference" "Tabula Sapiens"、"Data source" `f"Open Targets {format_data_version(genetics['data_version'])}"`
- Settings の区分: Comparison scope（`disease_selector(snapshot, default_diseases, prefix="genetics-")`）、Genetic evidence（`genetics-score` の `dcc.Input`。ⓘ "The score threshold keeps genes whose Open Targets genetic association score for the disease is at or above this value. Scores below 0.1 are not stored, so the threshold cannot go below 0.1. Empty or out-of-range values use 0.5."）、Expression criteria（薬剤ページと同じ 3 つを `genetics-` 付きで）、Display（`genetics-measure` の Percent / Count だけ。ⓘ "Count is the number of qualifying genes. Percent divides that number by all genes at or above the score threshold for the disease; the denominator is the same for every cell type."）
- `dcc.Store(id="genetics-applied-parameters", data=dict(zip(GENETICS_PARAMETER_IDS, ("percent", DEFAULT_SCORE_THRESHOLD, "specificity", DEFAULT_EXPRESSION_THRESHOLD, DEFAULT_SPECIFICITY_THRESHOLD, ordered_disease_ids(snapshot, default_diseases)), strict=True)))`、`genetics-expanded-cell-groups`、`genetics-expanded-expression-groups`
- 比較図は 1 つだけ: `html.Article([html.H3("Genes"), ... dcc.Graph(id="genetics-heatmap"), html.Div(id="genetics-row-controls", className="heatmap-row-controls") ..., html.Div(id="genetics-dot-key", className="dot-key", style={"display": "none"})], id="genetics-heatmap-panel")`。`genetics-chart-type` の RadioItems と `genetics-matrix-note`
- 詳細: `genetics-detail-disease` の Dropdown と `html.Div(id="genetics-details", className="details")`
- `default_diseases = choose_defaults([(d["id"], d["name"]) for d in snapshot["diseases"]], DEFAULT_DISEASE_TERMS, 10, set(genetics["associations"]))`
- 戻り値 `html.Main([...], className="shell", id="genetics-page")`

`genetics_unavailable_page(reason)`:

```python
def genetics_unavailable_page(reason: str) -> html.Main:
    """genetics.json が無いか版が合わないときの案内。"""
    return html.Main(
        [
            page_nav("genetics"),
            html.Section(
                [
                    html.H1("Genetic associations are not available"),
                    html.P(reason),
                    html.P(
                        "Refresh the genetic association data with pixi run refresh-genetics after pixi run refresh, then restart the app."
                    ),
                ],
                className="panel",
            ),
        ],
        className="shell",
        id="genetics-page",
    )
```

- [ ] **Step 5: `layout.py` の 10 疾患を共有する**

`layout.py:250-261` のタプルを `DEFAULT_DISEASE_TERMS` の参照に置き換える。循環 import を避けるため、`DEFAULT_DISEASE_TERMS` は `config.py` に置き、`layout.py` と `genetics_layout.py` の両方が `config` から読む（上の `genetics_layout.py` の定義は `config.py` へ移す）。

- [ ] **Step 6: テストを通す**

Run: `pixi run --locked --no-install python -m unittest tests.test_genetics_ui tests.test_ui -v`
Expected: PASS

- [ ] **Step 7: 検査とコミット**

Run: `pixi run --locked --no-install check`

Tomoya の許可を得てから:

```bash
git add autoimmune_atlas/ui/genetics_layout.py autoimmune_atlas/ui/config.py autoimmune_atlas/ui/layout.py tests/test_genetics_ui.py
git commit -m "feat(genetics): 遺伝子ページのレイアウトと関連遺伝子の表と詳細欄を追加"
```

---

### Task 6: 遺伝子ページの callback

**Files:**
- Create: `autoimmune_atlas/ui/genetics_callbacks.py`
- Test: `tests/test_genetics_ui.py`

**Interfaces:**
- Consumes: `genetics.summarize_genes`、`genetics.genes_for_disease`、`genetics_layout.genetics_detail_panel`、`figures.build_figure`、`figures.build_dot_figure`、`figures.expression_figure`、`figures.expression_view`、`figures.heatmap_cell_ids`、`figures.dot_size_scale`、`figures.expression_dot_size_scale`、`figures.disease_label_lines`、`components.heatmap_row_controls`、`components.disease_checklist_sections(family, "genetics-")`、`callbacks.effective_filters`、`callbacks._effective_number`（公開名 `effective_number` を足す）、`callbacks.resolve_disease_selection`
- Produces:
    - `genetics_callbacks.effective_score(value: NumberInput) -> tuple[float, str | None]`（空欄と範囲外は `DEFAULT_SCORE_THRESHOLD`。`SCORE_FLOOR` 未満も範囲外）
    - `genetics_callbacks.register_genetics_callbacks(application: Dash, snapshot: Snapshot, genetics: GeneticsSnapshot) -> None`。`snapshot` は `merged_snapshot` を通したもの

`callbacks.py:90` の `_effective_number` に `minimum: float | None = None` を足し、`number < 0` を `number < (minimum if minimum is not None else 0)` にする。エラー文の `" at or above 0"` も `minimum` に合わせる。公開名 `effective_number = _effective_number` を定義の直後に置く。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_genetics_ui.py` に追加。`tests/test_ui.py:1501-1560` の `CallbackTests` と同じく Flask の test client で callback を呼ぶ。`_callback_key`、`_response_json`、`CallbackDefinition` などの補助は `test_ui.py` から import せず、必要な最小のものをこのファイルに写す（`test_ui.py` は 2400 行あり、import すると fixture の読み込みまで走る）。

```python
class GeneticsCallbackTests(unittest.TestCase):
    application: ClassVar[Dash]
    client: ClassVar[FlaskClient]

    @classmethod
    def setUpClass(cls) -> None:
        cls.application = create_app(
            core_ui_snapshot(), assets_folder=ASSETS_PATH, genetics=genetics_snapshot()
        )
        cls.client = cls.application.server.test_client()

    def _callback(
        self,
        output: str,
        inputs: dict[str, object],
        states: dict[str, object] | None = None,
    ) -> object:
        """Dash の /_dash-update-component を 1 回呼び、指定した output の値を返す。"""
        callback_map = cast(dict[str, object], self.application.callback_map)
        key = next(k for k in callback_map if output in k)
        payload = {
            "output": key,
            "outputs": [
                {"id": i.split(".")[0], "property": i.split(".")[1]}
                for i in key.strip(".").split("...")
            ],
            "inputs": [
                {"id": k.split(".")[0], "property": k.split(".")[1], "value": v}
                for k, v in inputs.items()
            ],
            "state": [
                {"id": k.split(".")[0], "property": k.split(".")[1], "value": v}
                for k, v in (states or {}).items()
            ],
            "changedPropIds": list(inputs),
        }
        response = self.client.post("/_dash-update-component", json=payload)
        try:
            self.assertEqual(response.status_code, 200)
            return cast(dict[str, object], response.get_json())["response"]
        finally:
            response.close()

    def test_score_input_falls_back_below_floor_and_on_empty(self) -> None:
        self.assertEqual(
            effective_score(0.05),
            (0.5, "Score threshold must be finite and between 0.1 and 1; using 0.5."),
        )
        self.assertEqual(effective_score(""), (0.5, None))
        self.assertEqual(effective_score(0.3), (0.3, None))

    def test_update_applies_parameters_and_details_follow(self) -> None:
        applied = {
            "genetics-measure": "count",
            "genetics-score": 0.5,
            "genetics-method": "fixed",
            "genetics-threshold": 0.5,
            "genetics-specificity": 0.5,
            "genetics-diseases": ["D1", "D2"],
        }
        response = self._callback(
            "genetics-heatmap.figure",
            {
                "genetics-applied-parameters.data": applied,
                "genetics-expanded-cell-groups.data": [],
                "genetics-chart-type.value": "heatmap",
            },
        )
        figure = cast(dict[str, dict[str, object]], response)["genetics-heatmap"][
            "figure"
        ]
        data = cast(list[dict[str, object]], cast(dict[str, object], figure)["data"])
        self.assertEqual(
            cast(dict[str, object], data[0])["x"], ["Disease one", "Disease two"]
        )
        note = cast(dict[str, dict[str, object]], response)["genetics-matrix-note"][
            "children"
        ]
        self.assertIn("Score threshold", str(note))
        details = self._callback(
            "genetics-details.children",
            {
                "genetics-detail-disease.value": "D1",
                "genetics-applied-parameters.data": applied,
            },
        )
        self.assertIn("Genetically associated genes for Disease one", str(details))
        empty = self._callback(
            "genetics-details.children",
            {
                "genetics-detail-disease.value": "D2",
                "genetics-applied-parameters.data": applied,
            },
        )
        self.assertIn("No genes at or above the score threshold", str(empty))
```

`_callback` の `key` の組み立ては Dash の内部形式に依存する。`tests/test_ui.py:1545-1600` の `_callback_key` と `_call` の実装を読んで同じ形にする（`callback_map` のキーは単一 output なら `"id.prop"`、複数なら `"..id1.prop1...id2.prop2.."`）。上のコードはその形を前提にした概略なので、実装時に `test_ui.py` の補助関数の形に揃える。

- [ ] **Step 2: 失敗を確かめる**

Run: `pixi run --locked --no-install python -m unittest tests.test_genetics_ui -v`
Expected: FAIL（`create_app()` が `genetics` を受け付けない、`genetics_callbacks` が無い）

`create_app` の `genetics` 引数は Task 7 で足す。このタスクでは先に `application.py` に引数だけ足し、渡されたら `register_genetics_callbacks` を呼ぶ最小の変更を入れる（レイアウトの切り替えは Task 7）。

- [ ] **Step 3: `genetics_callbacks.py` を書く**

`callbacks.register_callbacks`（`callbacks.py:187-711`）を手本に、次の callback を登録する。ID はすべて `genetics-` 付き。

1. `genetics-expanded-cell-groups` の toggle。`Input({"type": "genetics-cell-toggle", "kind": ALL, "cell": ALL}, "n_clicks")`。本体は `callbacks.toggle_cell_group` と同じ判定（`ctx.triggered_id` と `ctx.inputs_list` を見る）。`callbacks.py:224-240` の関数を `heatmap_groups` を引数に取る形で `callbacks.py` に切り出し（`make_toggle(heatmap_groups)`）、両方から使う
2. `genetics-expanded-expression-groups` の toggle。`Input({"type": "genetics-expression-cell-toggle", ...})`
3. `genetics-applied-parameters` ← `genetics-update-button` の `n_clicks` と `GENETICS_PARAMETER_IDS` の `State`
4. `genetics-update-status` ← `genetics-applied-parameters` と `GENETICS_PARAMETER_IDS` の `Input`
5. `genetics-diseases` と `genetics-disease-family-*` の同期。`disease_checklist_sections(family, "genetics-")` の ID を使い、本体は `callbacks.py:383-401` と同じ
6. `genetics-specificity` の `disabled` ← `genetics-method`
7. `genetics-measure` の `options` ← `genetics-chart-type`
8. 比較図: Output は `genetics-heatmap` の `figure` と `style`、`genetics-matrix-note` の `children`、`genetics-row-controls` の `children` と `style`、`genetics-dot-key` の `children` と `style`。Input は `genetics-applied-parameters`、`genetics-expanded-cell-groups`、`genetics-chart-type`。`callbacks.py:450-638` の `update_figures` から `kind="gene"` の 1 図分だけを写す。集計は `@lru_cache(maxsize=1)` の `comparison_rows(score, minimum, disease_ids, cell_ids, method, specificity)` → `gene_data.summarize_genes(snapshot, genetics, score, minimum, method=method, specificity_threshold=specificity, level="mixed", cell_ids=list(cell_ids), disease_ids=list(disease_ids))`。note は `f"Applied · Score threshold: ≥ {score:g} · Expression rule: {METHOD_LABELS[method]} ({condition}) · Cells: all groups and expanded cell types · {display} · {len(rows)} disease–cell combinations"`。`display` の文言から "canonical drugs" を除く（`(missing: genes {gene_missing})`）。`heatmap_row_controls(heatmap_names, cell_ids, expanded or [], "genetics")`
9. `genetics-detail-disease` の `options` と `value` ← `genetics-applied-parameters`、`genetics-heatmap` の `clickData`、`State("genetics-detail-disease", "value")`。`resolve_disease_selection` は `triggered_id in {"target-heatmap", "drug-heatmap", "heatmap"}` を見るので、`"genetics-heatmap"` を集合に足す（`callbacks.py:177`）
10. `genetics-details` ← `genetics-detail-disease`、`genetics-applied-parameters`。`summarize_genes(..., level="all", disease_ids=[selected])` の結果を `genetics_detail_panel` に渡す
11. 発現図: `callbacks.py:250-344` の `expression_base` と `update_expression` を写し、`records = [{"target_id": g["target_id"], "target": g["target"]} for g in gene_data.genes_for_disease(genetics, disease, score)]` を `expression_figure(snapshot, records, threshold=..., method=..., specificity=..., grouped=True)` に渡す。Input は `genetics-source-context`、`genetics-expanded-expression-groups`、`genetics-expression-chart-type`。行見出しは `heatmap_row_controls(heatmap_names, cells, expanded or [], "genetics-expression")`

`effective_score`:

```python
def effective_score(value: NumberInput) -> tuple[float, str | None]:
    """空欄は初期値、0.1 未満と 1 超と不正値は理由を示して初期値へ戻す。"""
    return effective_number(
        value, DEFAULT_SCORE_THRESHOLD, "score threshold", 1, minimum=SCORE_FLOOR
    )
```

`_effective_number` のエラー文は `f"{label[0].upper() + label[1:]} must be finite and{bound}; using {default:g}."` で、`bound` は `minimum` と `maximum` があるとき `f" between {minimum:g} and {maximum:g}"`。上のテストの期待文はこれに合わせてある。

- [ ] **Step 4: `application.py` に `genetics` 引数を足す（最小）**

```python
def create_app(
    snapshot: Snapshot | None,
    error: Exception | None = None,
    *,
    assets_folder: str | Path,
    genetics: GeneticsSnapshot | None = None,
) -> Dash:
```

`snapshot is not None and genetics is not None` のとき、`merged = merged_snapshot(snapshot, genetics)` を作り、`register_genetics_callbacks(application, merged, genetics)` を呼ぶ。レイアウトはこのタスクでは変えない（`suppress_callback_exceptions=True` なので、レイアウトに無い ID への callback 登録は失敗しない）。

- [ ] **Step 5: テストを通す**

Run: `pixi run --locked --no-install python -m unittest tests.test_genetics_ui tests.test_ui -v`
Expected: PASS

- [ ] **Step 6: 検査とコミット**

Run: `pixi run --locked --no-install check`

Tomoya の許可を得てから:

```bash
git add autoimmune_atlas/ui/genetics_callbacks.py autoimmune_atlas/ui/callbacks.py autoimmune_atlas/ui/application.py tests/test_genetics_ui.py
git commit -m "feat(genetics): 遺伝子ページの比較図と詳細と発現図の callback を追加"
```

---

### Task 7: ページの切り替え、上部バー、起動

**Files:**
- Modify: `autoimmune_atlas/ui/application.py`
- Modify: `autoimmune_atlas/ui/layout.py:309-332`（`Nav` を `page_nav("drugs")` に）、`:307`（`html.Main` に `id="drug-page"`）
- Modify: `app.py:13-24`
- Create: `assets/router.js`
- Modify: `assets/style.css`（`.app-page-link[aria-current="page"]` の下線）
- Test: `tests/test_genetics_ui.py`

**Interfaces:**
- Produces:
    - `create_app(...)` の `layout` は `html.Div([dcc.Location(id="url", refresh=False), drug_page, genetics_page])`。`drug_page` は `dashboard_layout(snapshot)`（ID `drug-page`）か `unavailable_layout(error)`。`genetics_page` は `genetics_page(merged, genetics)`、`genetics is None` なら `genetics_unavailable_page("The genetic association data have not been retrieved.")`、版が違えば `genetics_unavailable_page("The genetic association data were retrieved from a different Open Targets release than the drug data.")`。`snapshot is None` なら `genetics_unavailable_page(...)` のみで callback は登録しない
    - callback: `Output("drug-page", "hidden")`、`Output("genetics-page", "hidden")` ← `Input("url", "pathname")`。`pathname == "/genetics"` なら `(True, False)`、それ以外は `(False, True)`
    - `assets/router.js`: `dcc.Location` の変化は `popstate` では拾えないので、`MutationObserver` で `#drug-page` と `#genetics-page` の `hidden` 属性の変化を監視し、変わったら `window.dispatchEvent(new Event("resize"))` を送る

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_genetics_ui.py` に追加:

```python
class RouterTests(unittest.TestCase):
    def test_layout_holds_both_pages_and_url_toggles_hidden(self) -> None:
        app = create_app(
            core_ui_snapshot(), assets_folder=ASSETS_PATH, genetics=genetics_snapshot()
        )
        client = app.server.test_client()
        layout = _response_json(client.get("/_dash-layout"))
        ids = _collect_ids(layout)
        self.assertIn("url", ids)
        self.assertIn("drug-page", ids)
        self.assertIn("genetics-page", ids)
        self.assertIn("diseases", ids)
        self.assertIn("genetics-diseases", ids)
        for path, expected in (
            ("/genetics", [True, False]),
            ("/", [False, True]),
            (None, [False, True]),
        ):
            response = _call(client, app, "drug-page.hidden", {"url.pathname": path})
            self.assertEqual(
                [response["drug-page"]["hidden"], response["genetics-page"]["hidden"]],
                expected,
            )
        for route in ("/", "/genetics"):
            page = client.get(route)
            try:
                self.assertEqual(page.status_code, 200)
            finally:
                page.close()

    def test_version_mismatch_and_missing_genetics_show_notice_but_keep_drug_page(
        self,
    ) -> None:
        stale = genetics_snapshot()
        stale["data_version"] = {"year": "26", "month": "06", "iteration": None}
        for data in (stale, None):
            app = create_app(
                core_ui_snapshot(), assets_folder=ASSETS_PATH, genetics=data
            )
            layout = _response_json(app.server.test_client().get("/_dash-layout"))
            ids = _collect_ids(layout)
            self.assertIn("drug-page", ids)
            self.assertIn("genetics-page", ids)
            self.assertNotIn("genetics-heatmap", ids)
            self.assertIn("refresh-genetics", json.dumps(layout))

    def test_nav_marks_current_page(self) -> None:
        drug = layout.dashboard_layout(core_ui_snapshot())
        links = [
            c for c in _walk(drug) if getattr(c, "className", None) == "app-page-link"
        ]
        self.assertEqual([getattr(c, "href") for c in links], ["/", "/genetics"])
        current = [c for c in links if getattr(c, "aria-current", None) == "page"]
        self.assertEqual([getattr(c, "href") for c in current], ["/"])
```

`_response_json`、`_collect_ids`、`_call` は `test_ui.py:1520-1600` の同名の補助と同じ形でこのファイルに写す。`getattr(c, "aria-current", None)` は Dash が ARIA 属性をそのまま属性名で持つことに依存するので、動かなければ `c.to_plotly_json()["props"].get("aria-current")` を使う。

- [ ] **Step 2: 失敗を確かめる**

Run: `pixi run --locked --no-install python -m unittest tests.test_genetics_ui.RouterTests -v`
Expected: FAIL（`url` が無い）

- [ ] **Step 3: `application.py` を書き換える**

```python
"""Dash アプリを組み立てる。"""

from pathlib import Path

from dash import Dash, Input, Output, dcc, html

from autoimmune_atlas.genetics import merged_snapshot, version_matches
from autoimmune_atlas.models import GeneticsSnapshot, Snapshot
from autoimmune_atlas.ui.callbacks import register_callbacks
from autoimmune_atlas.ui.genetics_callbacks import register_genetics_callbacks
from autoimmune_atlas.ui.genetics_layout import genetics_page, genetics_unavailable_page
from autoimmune_atlas.ui.layout import dashboard_layout, unavailable_layout

NOT_RETRIEVED = "The genetic association data have not been retrieved."
VERSION_MISMATCH = (
    "The genetic association data were retrieved from a different Open Targets "
    "release than the drug data."
)


def create_app(
    snapshot: Snapshot | None,
    error: Exception | None = None,
    *,
    assets_folder: str | Path,
    genetics: GeneticsSnapshot | None = None,
) -> Dash:
    """保存済みデータまたは明示的な fixture から、2 ページの Dash アプリを作る。"""
    application = Dash(
        __name__,
        title="Autoimmune Atlas",
        suppress_callback_exceptions=True,
        assets_folder=str(assets_folder),
    )
    drug_page = (
        unavailable_layout(error) if snapshot is None else dashboard_layout(snapshot)
    )
    if snapshot is None or genetics is None:
        gene_page = genetics_unavailable_page(NOT_RETRIEVED)
    elif not version_matches(snapshot, genetics):
        gene_page = genetics_unavailable_page(VERSION_MISMATCH)
    else:
        merged = merged_snapshot(snapshot, genetics)
        gene_page = genetics_page(merged, genetics)
        register_genetics_callbacks(application, merged, genetics)
    gene_page.hidden = True
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
        on_genetics = pathname == "/genetics"
        return on_genetics, not on_genetics

    return application
```

`unavailable_layout(error)` が返す部品にも `id="drug-page"` を付ける（`layout.py` の該当箇所を確かめて足す）。

- [ ] **Step 4: `layout.py` の `Nav` を置き換える**

`layout.py:309-332` の `html.Nav([...])` を `page_nav("drugs")` に置き換え、`html.Main(...)` に `id="drug-page"` を足す。`page_nav` は `genetics_layout.py` から import する。`genetics_layout.py` は `layout.py` から `choose_defaults` と `format_data_version` を import しているので、循環 import になる。`page_nav` は `components.py` に置き、両方が `components` から import する形に変える。

- [ ] **Step 5: `app.py` を直す**

```python
from autoimmune_atlas.genetics import load_genetics

...
try:
    snapshot_data = load_snapshot()
    genetics_data = load_genetics()
    load_error = None
except (OSError, json.JSONDecodeError, ValueError, TypeError, KeyError) as error:
    snapshot_data = None
    genetics_data = None
    load_error = error

app = create_app(
    snapshot_data,
    load_error,
    assets_folder=BASE_DIR / "assets",
    genetics=genetics_data,
)
```

`genetics.json` の読み込みだけが失敗したときに薬剤ページまで止めないよう、`load_genetics()` は別の `try` で包み、失敗したら `genetics_data = None` にして標準エラーに理由を出す。

- [ ] **Step 6: `assets/router.js` と CSS**

`assets/router.js`:

```javascript
// ページの hidden が切り替わったら、Plotly に幅を測り直させる。
// hidden の中で描かれた図は幅 0 のままで、窓の大きさが変わるまで描き直されない。
function watchPages() {
    const pages = ["drug-page", "genetics-page"].map((id) => document.getElementById(id)).filter(Boolean);
    if (pages.length < 2) return false;
    const observer = new MutationObserver(() => window.dispatchEvent(new Event("resize")));
    pages.forEach((page) => observer.observe(page, {attributes: true, attributeFilter: ["hidden"]}));
    return true;
}

const bootstrap = new MutationObserver(() => {
    if (watchPages()) bootstrap.disconnect();
});
bootstrap.observe(document.body, {childList: true, subtree: true});
```

`assets/style.css` の `.app-nav a` の近くに:

```css
.app-page-link[aria-current="page"] { text-decoration: underline; text-underline-offset: 4px; }
```

- [ ] **Step 7: テストを通す**

Run: `pixi run --locked --no-install python -m unittest tests.test_genetics_ui tests.test_ui -v`
Expected: PASS。`test_ui.py` の `CallbackTests` は `/_dash-layout` から ID を集めるので、shell が 1 段増えても通る。`test_layout_has_all_modalities_and_defaults` などが `layout.dashboard_layout(...)` の直下の構造を見ているなら、`Nav` の置き換えで壊れていないかを確かめる

- [ ] **Step 8: 実ブラウザーで確かめる**

Tomoya に `pixi run --locked --no-install start` を頼み、Chrome で次を見る。

1. `http://127.0.0.1:8050/` で薬剤ページが出て、上部バーに 2 つのリンクがあり、`Drug targets` に下線がある
2. `Genetic associations` を押すと `/genetics` になり、比較図が幅いっぱいに描かれる（幅 0 のままなら `router.js` を直す）
3. `/genetics` で閾値を 0.3 にして `Update`、詳細で関節リウマチを選び、表と発現図が出る。遺伝子数が多い疾患の発現図が横スクロールで読めるかを見る
4. `Drug targets` へ戻ると、薬剤ページの設定と展開状態が残っている
5. `http://127.0.0.1:8050/genetics` を直接開いても 2 と同じになる
6. pemphigus の列が 0 で、ホバーに "No genes at or above the score threshold" が出る

見つかった問題は、このタスクの中で直す。

- [ ] **Step 9: 検査とコミット**

Run: `pixi run --locked --no-install check`

Tomoya の許可を得てから:

```bash
git add autoimmune_atlas/ui/application.py autoimmune_atlas/ui/layout.py autoimmune_atlas/ui/components.py autoimmune_atlas/ui/genetics_layout.py app.py assets/router.js assets/style.css tests/test_genetics_ui.py
git commit -m "feat(ui): URL で薬剤ページと遺伝子ページを切り替える上部バーと router を追加"
```

---

### Task 8: 設計ノートと README の更新

**Files:**
- Modify: `docs/design.md:274-284`（「初期版の範囲を保つ」）、`:286-291`（根拠）
- Modify: `docs/genetics-design.md`（実装で変えた点を反映）
- Modify: `README.md`（コード構成の木と `refresh-genetics`）
- Modify: `docs/README.md`

- [ ] **Step 1: `docs/design.md` を直す**

「初期版の範囲を保つ」の `分子 X の入力と横並び比較、組織別の比較、治療目的の分類、薬効細胞の手作業注釈は今回実装しない。` の前に 1 段落:

```markdown
遺伝学的関連遺伝子のページは、別の取得データ（`data/genetics.json`）と別の集計（`autoimmune_atlas/genetics.py`）で `/genetics` に置き、薬剤ページの発現判定と細胞分類と図を共有する。
ページの切り替えは `dcc.Location` で行い、両ページのレイアウトを起動時に組み立てて `hidden` を入れ替える。
遺伝子集合の定め方、保存の形、表示の規則は[遺伝学的関連遺伝子のページの設計](genetics-design.md)にある。
```

「根拠と仕様の対応を残す」に 1 行:

```markdown
- [Open Targets: association scores](https://platform-docs.opentargets.org/associations)：genetic association の datatype スコアと datasource スコアの合成。`autoimmune_atlas/refresh_genetics.py` が取得し、`autoimmune_atlas/genetics.py` が集計する。
```

- [ ] **Step 2: `docs/genetics-design.md` を実装に合わせる**

- 「集計は発現判定と細胞分類だけを共有する」の表: `genes` の行を消し、「返す行は薬剤ページの `SummaryRow` と同じ型で、薬剤の列は空の値（0、`None`、空の配列）で埋める。比較図を `kind="gene"` で共有するためである」と書く
- 「ページは URL で切り替え」: `page_nav` を `components.py` に置いたこと、`router.js` が `MutationObserver` で `hidden` を監視すること
- 「遺伝子ページの操作と表示」: 詳細の表の列を実装どおり（Gene、Genetic association score、datasource ごと、Open Targets）に揃える。datasource の表示名は ID のアンダースコアを空白に置き換えたもの
- 実データで確かめた結果（発現図の可読性、pemphigus の表示）を「検証」に 1 行ずつ足す

- [ ] **Step 3: `README.md` と `docs/README.md`**

`README.md` のコード構成の木に `genetics.py`、`refresh_genetics.py`、`ui/genetics_layout.py`、`ui/genetics_callbacks.py`、`assets/router.js`、`tests/test_genetics.py`、`tests/test_refresh_genetics.py`、`tests/test_genetics_ui.py` を足す。冒頭の説明に「遺伝学的関連遺伝子のページ（`/genetics`）も持つ」と 1 文足す。データ更新の説明に `pixi run refresh-genetics` を足し、`pixi run refresh` のあとに実行することを書く。

`docs/README.md` はすでに `genetics-design.md` を載せている。変更なし。

- [ ] **Step 4: 点検とコミット**

Run: `pixi run --locked --no-install lint-markdown && pixi run --locked --no-install check`

Tomoya の許可を得てから:

```bash
git add docs/design.md docs/genetics-design.md README.md
git commit -m "docs(genetics): 遺伝子ページの実装に合わせて設計ノートと README を更新"
```

---

## 自己点検の記録

- 仕様の各節と担当タスク: 遺伝子集合の定義 → Task 2、6。データの取得と保存 → Task 1、3。集計の共有 → Task 2、4。ページの切り替え → Task 7。操作と表示 → Task 5、6。説明と用語 → Task 5。検証 → 各タスクのテストと Task 7 Step 8
- 仕様の `summarize_genes` の `genes` キーは、図を共有するために `SummaryRow` と同じ型に変えた。Task 8 で仕様を直す
- 型の一貫: `Kind` に `gene`、`effective_number(value, default, label, maximum=None, *, minimum=None)`、`disease_selector(snapshot, selected, prefix="")`、`heatmap_row_controls(..., kind)` の kind は `"genetics"` と `"genetics-expression"`、`page_nav` は `components.py`
- Review Focus の 5 項目は、それぞれ Task 1、2、4、6、7 のテストか手動確認に入っている
