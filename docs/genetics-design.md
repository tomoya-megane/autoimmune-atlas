# 遺伝学的関連遺伝子の発現から細胞への支持を比較する

薬剤の標的の代わりに、Open Targets の genetic association で疾患と結びついた遺伝子を使って、同じ細胞型別の比較を別ページに表示する。
この文書は、実装前に決めた仕様である。
実装を終えたら、集計規則と表示規則を[現行の設計](design.md)へ移し、この文書は経緯の記録として残す。

薬剤ページは「疾患 → 臨床開発された薬剤 → 標的遺伝子 → 健常参照での細胞型別発現」を比較する。
遺伝子ページはこの鎖の前半を「疾患 → 遺伝学的に関連する遺伝子」に置き換え、後半の発現データ、3 つの発現基準、細胞の大分類と展開、実数と割合の表示規則をそのまま使う。
薬剤という単位が無いので、比較図は標的数に相当する 1 種類だけになる。
ここでいう支持も、健常時の発現から推定した相対的な対応であり、疾患の原因細胞や患者組織での発現変化を表さない。

## 疾患ごとの遺伝子集合をスコアの閾値で定める

Open Targets の疾患と標的の関連スコアのうち、datatype `genetic_association` のスコアが閾値以上の遺伝子を、その疾患の関連遺伝子とする。
初期値は 0.5 で、画面から変更できる。
薬剤ページの臨床段階フィルターに当たる設定であり、比較図の分子と分母の両方がこの集合で決まる。

datatype `genetic_literature`（UniProt の文献由来の変異記述）は含めない。
Open Targets の datatype の区分に合わせ、独自にスコアを合成しない。

疾患の下位語の根拠は含めない（`enableIndirect: false`）。
薬剤ページと同じ単位で疾患を比べるためである。
ACPA 陽性関節リウマチのようなサブタイプは、疾患の選択欄で親と並んで個別に選べる。

保存する下限スコアは 0.1 とし、画面ではこの値より下に閾値を下げられない。
再取得せずに閾値を変えられる範囲を、保存の時点で決めるためである。
0.1 未満の関連は保存しないので、閾値の入力欄は 0.1 から 1 の範囲に限る。

閾値以上の遺伝子が 0 件の疾患は、取得済みだが該当なしとして扱う。
比較図では実数を 0、割合を NA とし、ホバーで「no genes at or above the score threshold」と示す。
未取得（not loaded）とは理由を分ける。
確かめた範囲では、pemphigus と dermatomyositis がこれにあたる。

## データは別ファイルに取得し、版が一致するときだけ表示する

`pixi run refresh-genetics` が `autoimmune_atlas.refresh_genetics` を実行し、`data/genetics.json` を書く。
既存の `data/snapshot.json`、そのスキーマ、`autoimmune_atlas/refresh.py` の保存部分は変えない。
薬剤だけ、遺伝子だけを取り直せるようにするためである。

取得の手順は次のとおり。

1. `snapshot.json` を読み、疾患の一覧と発現データを持つ標的の一覧を得る。取得元の疾患の範囲は `snapshot.json` に従い、起点からの解決をやり直さない
2. データの版を `meta.dataVersion` で取得し、`snapshot.json` の `data_version` と一致しなければ、その場で止めて `pixi run refresh` を先に実行するよう案内する
3. 疾患ごとに `Disease.associatedTargets(enableIndirect: false, orderByScore: "genetic_association")` を 500 件ずつ取得し、genetic association のスコアが 0.1 を下回った時点でその疾患のページングを止める
4. 遺伝子ごとに、`datatypeScores` の genetic association のスコアと、`datasourceScores` のうち genetic association に属する datasource のスコアを保存する。どの datasource が genetic association に属するかは、取得した結果から決める。全疾患の `datasourceScores` に現れた datasource ID ごとに、それを持つ疾患と遺伝子の組を 1 つ選び、`Disease.evidences(ensemblIds: [遺伝子], datasourceIds: [datasource], size: 1)` で根拠を 1 件取り、その `datatypeId` が `genetic_association` のものだけを残す。datasource の一覧をコードに書かないのは、版が上がって datasource が増えても取りこぼさないためである。API の `associationDatasources` は版 26.09 で空の配列を返したので、使わない
5. `snapshot.json` の `expression` に無い遺伝子について、`refresh.py` の `fetch_expression` を使って発現データを取得する
6. 終了時に版をもう一度取得し、開始時と違えば失敗にする。全件成功したときだけ、一時ファイルから `genetics.json` を置き換える

`orderByScore: "genetic_association"` の並びは、genetic association のスコアを持つ遺伝子を降順に置き、スコアを持たない遺伝子をその後ろに置く。
関節リウマチで確かめたところ、全件を取って数えた 0.1 以上の件数と、並べて 0.1 を下回るまで数えた件数が一致した。

取得の量は、版 26.09 に対して測った。
スナップショットの 370 疾患のうち、0.1 以上の遺伝子を 1 件以上持つ疾患はおよそ半数で、中央値は 1 件、最大は 742 件である。
0.1 以上の遺伝子の和集合はおよそ 4,200 件で、そのうち `snapshot.json` に発現データがあるのは約 250 件にとどまる。
発現データの取得件数は、薬剤の標的 764 件の約 5 倍になる。
`refresh` の発現の取得にかかる時間を測り、その 5 倍を目安にする。
版が違うので、この数は次の版では変わる。

`genetics.json` の形は次のとおり。

| キー | 中身 |
| --- | --- |
| `schema` | 1 |
| `data_version` | 取得時の Open Targets の版 |
| `retrieved_at` | 取得日時 |
| `source` | API の URL |
| `score_floor` | 保存した下限スコア。0.1 |
| `datasources` | genetic association に属すると判定した datasource の ID の一覧。表示名は画面側で ID から作る |
| `associations` | 疾患 ID をキーに、遺伝子の配列。各要素は `target_id`、`target`、`score`、`datasource_scores`（datasource ID をキーにしたスコア） |
| `expression` | `snapshot.json` に無い遺伝子の発現データ。形は `snapshot.json` の `expression` と同じ |

アプリは起動時に 2 つのファイルを読む。
`create_app()` には `genetics` を省略できる引数として足し、既存の呼び方とテストの fixture はそのまま動くようにする。
`genetics.json` が無いか、`data_version` が `snapshot.json` と一致しないときは、遺伝子ページに「Refresh the data with pixi run refresh-genetics」と表示し、比較図を出さない。
薬剤ページはこのとき影響を受けない。
発現データは 2 つの `expression` を結合して使い、同じ遺伝子が両方にあるときは `snapshot.json` を優先する。
同じ版から取得しているので中身は一致するはずだが、上書きせずに済む側を決めておく。

遺伝子ごとの代表的な根拠（GWAS の study、L2G スコア、変異）は初期版では保存しない。
根拠の行は疾患と遺伝子の組ごとに数十件あり、疾患全体では数万件になる。
代表 1 件を選ぶと、それが疾患関連スコアの唯一の根拠と読まれる。
詳細の表からは Open Targets の evidence ページ（`/evidence/<遺伝子 ID>/<疾患 ID>`）へリンクし、根拠はそこで読む。

## 集計は発現判定と細胞分類だけを共有する

`autoimmune_atlas/genetics.py` に、遺伝子ページ用の集計を置く。
`aggregation.summarize()` は入口で薬剤専用の絞り込みを呼び、結果にも薬剤数と標的不明薬剤数を持つので、一般化せず、遺伝子を薬剤記録の形に見せかけることもしない。
共有するのは `expression_metadata()`、`expression_state()`、`cell_catalog()` と、三値の集約である。

`summarize_genes(snapshot, genetics, score_threshold, threshold, *, method, specificity_threshold, level, cell_ids, disease_ids)` は、疾患と細胞の組ごとに次を返す。

| キー | 中身 |
| --- | --- |
| `disease_id`、`disease`、`cell_id`、`cell`、`ontology_id`、`cell_level`、`member_cell_ids` | 薬剤ページの `SummaryRow` と同じ |
| `count` | 発現判定が陽性の遺伝子数。疾患が未取得なら `None` |
| `percent` | 陽性の遺伝子数を分母で割った割合。分母が 0、または陽性が無く未判定が残るなら `None` |
| `denominator` | 閾値以上の遺伝子数 |
| `unknown` | 未判定の遺伝子数 |
| `status` | `unavailable`、`partial`、`complete`。薬剤ページと同じ意味 |
| `genes` | 陽性の遺伝子の ID と名前と genetic association のスコア |

三値判定、大分類の集約、下限値の扱い、割合の NA の規則は、薬剤ページの標的の規則をそのまま適用する。
薬剤の網羅性（`unmapped_drugs`、`total_drugs`）に当たるものは無い。

閾値の判定はスコアが閾値以上（等号を含む）とする。
Open Targets の表示と等号まで一致すると断定しない。

## ページは URL で切り替え、設定は両方を保持する

`dcc.Location` を 1 つ置き、パス名で 2 つのページを切り替える。
`/` が薬剤ページ、`/genetics` が遺伝子ページである。
Dash の pages 機能は使わない。
pages の登録先はプロセス全体で 1 つなので、fixture ごとに別のアプリを作る今のテストと相性が悪い。

2 つのページのレイアウトは起動時に両方組み立て、`dcc.Location` の callback が `hidden` を切り替える。
ページを切り替えても、もう一方の設定、適用済みの条件、展開状態、詳細の選択が失われないようにするためである。
遺伝子ページの部品の ID には `genetics-` を前に付け、薬剤ページの ID と重ねない。

Plotly の図は `hidden` の中では幅 0 で描かれ、表示に戻しても窓の大きさが変わるまで描き直されない。
ページを切り替えたときに `window` へ `resize` イベントを送る小さな JavaScript を `assets/` に置き、既存のツールチップの JavaScript と同じ形で扱う。
これで足りるかは実ブラウザーで確かめる。

上部バーのアプリ名の右に、`Drug targets` と `Genetic associations` の 2 つのリンクを置き、表示中のページに `aria-current="page"` を付ける。
Open Targets へのリンクはその右に残す。

## 遺伝子ページの操作と表示

上部パネルの区分は次のとおり。
薬剤ページと同じ順序で、Drug evidence だけが Genetic evidence に変わる。

| 区分 | 項目 |
| --- | --- |
| Comparison scope | 疾患。薬剤ページと同じ選択欄の部品を、別の ID で置く |
| Genetic evidence | genetic association のスコアの閾値。0.1 から 1 の数値入力、初期値 0.5。空欄と範囲外は初期値に戻し、適用した値を画面に記録する。既存の数値の検証は 0 未満しか弾かないので、下限 0.1 の検証を足す |
| Expression criteria | 発現基準、最低 CPM、CELLEX の閾値。薬剤ページと同じ |
| Display | Heatmap の実数／割合。標的／薬剤の切り替えは無い |

`Update` で上部パネルの設定をまとめて適用し、未反映の変更があることをボタンの隣に表示する規則も同じにする。
初期表示の疾患は薬剤ページと同じ 10 疾患とする。
閾値以上の遺伝子が 0 件の疾患が含まれるが、該当なしの表示を初期状態で確かめられる。

概要欄には、`genetics.json` 全体の値として、取得済みの疾患数、閾値 0.1 以上の関連の数、発現データを持つ遺伝子数、取得日時、版を示す。

見出しは `Genetically associated genes across diseases and cell types` とし、ⓘ の説明に次を書く。
差は疾患ごとの関連遺伝子の構成の違いであり、健常参照データを使うので疾患別の発現変化ではないこと。
関連スコアは Open Targets が複数の datasource から合成した値であり、因果や薬効を表さないこと。

比較図は Heatmap と Dot plot を切り替えられ、疾患 × 細胞の大分類に、閾値以上の遺伝子のうち発現判定が陽性の数と割合を示す。
大分類の展開、色範囲、円の面積、欠測の表示は薬剤ページの標的の図の規則を使う。
ホバーには疾患名、細胞名、実数、割合、分母、未判定の数を示し、薬剤の網羅性の行は出さない。
図を組み立てる関数は、`SummaryRow` の標的側の列だけを読む形で薬剤ページと共有できるかを実装時に確かめ、共有できなければ遺伝子用の関数を `ui/figures.py` に別に置く。

詳細欄は、比較図のクリックか選択欄で疾患を選ぶと更新する。
上から、遺伝子の根拠の表、遺伝子の発現図の順に置く。

根拠の表は AgGrid で、閾値以上の遺伝子 1 件を 1 行にする。

| 列 | 中身 |
| --- | --- |
| Gene | 遺伝子記号と Ensembl ID |
| Genetic association score | datatype のスコア。降順を初期の並びにする |
| 各 datasource | datasource ごとのスコア。無ければ空欄。列は `genetics.json` の `datasources` から作る |
| Open Targets | その遺伝子と疾患の evidence ページへのリンク |

列見出しのクリックで並べ替え、見出しの下の入力欄で絞り込めるようにする点は薬剤の表と同じにする。
発現判定はこの表に含めない。

発現図は、閾値以上の遺伝子について全元細胞型の発現値を示す。
z スコア、クラスタリング、大分類の平均、白丸の規則は薬剤ページの発現図と同じである。
遺伝子が 100 件を超える疾患では横幅が長くなるので、図内に横スクロールを設ける既存の枠をそのまま使う。
読めるかどうかは、実データで画面を確かめてから判断する。
関節リウマチで閾値 0.5 なら 140 件、0.1 なら 420 件である。

## 説明と用語

画面では、Open Targets の呼び方に合わせて `genetic association score` と書き、`genetic evidence` は区分の名前にだけ使う。
根拠の表の datasource 名は、`genetics.json` に保存した表示名を使う。
値が表示できない理由は、薬剤ページと同じく missing、unresolved、not loaded で呼び分け、閾値以上の遺伝子が無い場合を `no genes at or above the score threshold` とする。

## 初期版の範囲

datasource による絞り込み、上位 N 件の選択、下位語を含めた集計、代表的な根拠の保存、薬剤ページと遺伝子ページを重ねた図は作らない。
薬剤の標的と遺伝学的関連遺伝子の重なりを見る比較は、2 つのページが揃ってから、別に設計する。

## 検証

- `tests/test_refresh_genetics.py`：偽の query 関数を渡し、ページングが閾値で止まること、版の不一致で止まること、既存の発現データにある遺伝子を再取得しないこと、全件成功時だけ保存することを確かめる
- `tests/test_genetics.py`：閾値の等号、分母 0、未判定が残るときの NA、大分類の集約を、薬剤ページの集計と同じ fixture で確かめる
- `tests/test_ui.py`：`/genetics` の初期レイアウトに遺伝子ページの部品があること、`dcc.Location` の切り替えで `hidden` が入れ替わること、遺伝子ページの `Update` と詳細の callback が応答すること、薬剤ページの既存のテストが変わらず通ることを確かめる
- 実ブラウザーで、ページの往復で設定が保持されること、該当なしの疾患の表示、遺伝子数の多い疾患の発現図を確かめる

## 根拠

- [Open Targets: association scores](https://platform-docs.opentargets.org/associations)：datatype と datasource のスコアの合成
- [Open Targets: genetic association evidence](https://platform-docs.opentargets.org/evidence#genetic-association)：`genetic_association` に属する datasource
- `autoimmune_atlas/refresh.py`：既存の取得と保存の流れ。`fetch_expression` を再利用する
- `autoimmune_atlas/aggregation.py`：共有する発現判定と細胞分類
