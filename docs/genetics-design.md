# 遺伝学的関連遺伝子の発現から細胞への支持を比較する

薬剤の標的の代わりに、Open Targets の genetic association で疾患と結びついた遺伝子を使って、同じ細胞型別の比較を別ページに表示する。
この文書は、遺伝学的関連遺伝子のページの設計である。
[現行の設計](design.md)は、遺伝子集合の定め方と表示の規則についてこの文書を参照し、保存の形は[保存データの形と取得の手順](data.md)にある。

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

保存する関連遺伝子に下限は設けず、スコアを持つ遺伝子を全部保存する。
画面の閾値は 0 から 1 の範囲で変えられる。
初期版では下限 0.1 を置いていたが、弱い関連が見えないことが分かったので無くした。

閾値以上の遺伝子が 0 件の疾患は、取得済みだが該当なしとして扱う。
比較図では実数を 0、割合を NA とし、ホバーで「no genes at or above the score threshold」と示す。
未取得（not loaded）とは理由を分ける。
確かめた範囲では、pemphigus と dermatomyositis がこれにあたる。

## 関連遺伝子は snapshot.json から読み、代表の根拠は持たない

関連遺伝子は `pixi run refresh` が薬剤の記録と一緒に取得し、`data/snapshot.json` の `associations` と `datasources` に保存する。
取得の規則と保存の形は[保存データの形と取得の手順](data.md)にある。
`snapshot.json` が無いか読めないときは、遺伝子ページに「Refresh the data with pixi run refresh, then restart the app.」と表示し、比較図を出さない。

遺伝子ごとの代表的な根拠（GWAS の study、L2G スコア、変異）は初期版では保存しない。
根拠の行は疾患と遺伝子の組ごとに数十件あり、疾患全体では数万件になる。
代表 1 件を選ぶと、それが疾患関連スコアの唯一の根拠と読まれる。
詳細の表からは Open Targets の evidence ページ（`/evidence/<遺伝子 ID>/<疾患 ID>`）へリンクし、根拠はそこで読む。

## 集計は発現判定と細胞分類だけを共有する

`autoimmune_atlas/genetics.py` に、遺伝子ページ用の集計を置く。
`aggregation.summarize()` は入口で薬剤専用の絞り込みを呼び、結果にも薬剤数と標的不明薬剤数を持つので、一般化せず、遺伝子を薬剤記録の形に見せかけることもしない。
共有するのは `expression_metadata()`、`expression_state()`、`cell_catalog()` と、三値の集約である。

`summarize_genes(snapshot, score_threshold, threshold, *, method, specificity_threshold, level, cell_ids, disease_ids)` は、疾患と細胞の組ごとに次を返す。

| キー | 中身 |
| --- | --- |
| `disease_id`、`disease`、`cell_id`、`cell`、`ontology_id`、`cell_level`、`member_cell_ids` | 薬剤ページの `SummaryRow` と同じ |
| `count` | 発現判定が陽性の遺伝子数。疾患が未取得なら `None` |
| `percent` | 陽性の遺伝子数を分母で割った割合。分母が 0、または陽性が無く未判定が残るなら `None` |
| `denominator` | 閾値以上の遺伝子数 |
| `unknown` | 未判定の遺伝子数 |
| `status` | `unavailable`、`partial`、`complete`。薬剤ページと同じ意味 |

返す行は薬剤ページの `SummaryRow` と同じ型で、薬剤の列は空の値（0、`None`、空の配列）で埋める。
比較図を `kind="gene"` で共有するためである。

三値判定、大分類の集約、下限値の扱い、割合の NA の規則は、薬剤ページの標的の規則をそのまま適用する。
薬剤の網羅性（`unmapped_drugs`、`total_drugs`）に当たるものは無い。

閾値の判定はスコアが閾値以上（等号を含む）とする。
Open Targets の表示と等号まで一致すると断定しない。

発現の索引（`expression_metadata()` の結果）と細胞の分類は、遺伝子ページの callback を登録するときに 1 回だけ作り、`summarize_genes` と発現図に渡す。
結合した発現データは薬剤ページの数倍の遺伝子を持つので、callback のたびに作り直すと 1 回の操作に十数秒かかった。
`summarize_genes` はこの 2 つを省略できる引数として受け取り、省略したときは自分で作る。

## ページは URL で切り替え、設定は両方を保持する

`dcc.Location` を 1 つ置き、パス名で 2 つのページを切り替える。
`/` が薬剤ページ、`/genetics` が遺伝子ページである。
Dash の pages 機能は使わない。
pages の登録先はプロセス全体で 1 つなので、fixture ごとに別のアプリを作る今のテストと相性が悪い。

2 つのページのレイアウトは起動時に両方組み立て、`dcc.Location` の callback が `hidden` を切り替える。
ページを切り替えても、もう一方の設定、適用済みの条件、展開状態、詳細の選択が失われないようにするためである。
遺伝子ページの部品の ID には `genetics-` を前に付け、薬剤ページの ID と重ねない。

Plotly の図は `hidden` の中では幅 0 で描かれ、表示に戻しても窓の大きさが変わるまで描き直されない。
`assets/router.js` が 2 つのページの `hidden` 属性を `MutationObserver` で監視し、変わったときに `window` へ `resize` イベントを送る。
`dcc.Location` によるパスの変化は `popstate` では検出できないので、`hidden` の変化を手がかりにする。
これで足りるかは実ブラウザーで確かめる。

上部バーのアプリ名のすぐ右に、`Drug targets` と `Genetic associations` の 2 つのリンクを左寄せで置く。
2 つのリンクは `dcc.Link` なので、切り替えてもページを再読み込みせず、もう一方のページの状態が残る。
`dcc.Link` は `aria-current` を受け取らないので、表示中のページは `current` クラスで示す。
Open Targets へのリンクは、外部へのリンクだと分かるようにバーの右端に置く。
上部バーは両ページで使うので、`page_nav` として `ui/components.py` に置く。

## 遺伝子ページの操作と表示

上部パネルの区分は次のとおり。
薬剤ページと同じ順序で、Drug evidence だけが Genetic evidence に変わる。

| 区分 | 項目 |
| --- | --- |
| Comparison scope | 疾患。薬剤ページと同じ選択欄の部品を、別の ID で置く |
| Genetic evidence | genetic association のスコアの閾値。0 から 1 の数値入力、初期値 0.5。空欄と範囲外は初期値に戻し、適用した値を画面に記録する |
| Expression criteria | 発現基準、最低 CPM、CELLEX の閾値。薬剤ページと同じ |
| Display | Heatmap の実数／割合。標的／薬剤の切り替えは無い |

`Update` で上部パネルの設定をまとめて適用し、未反映の変更があることをボタンの隣に表示する規則も同じにする。
初期表示の疾患は薬剤ページと同じ 10 疾患とする。
閾値以上の遺伝子が 0 件の疾患が含まれるが、該当なしの表示を初期状態で確かめられる。

概要欄には、`snapshot.json` 全体の値として、取得済みの疾患数、関連の数、発現データを持つ遺伝子数、取得日時、版を示す。

見出しは `Genetically associated genes across diseases and cell types` とし、ⓘ の説明に次を書く。
差は疾患ごとの関連遺伝子の構成の違いであり、健常参照データを使うので疾患別の発現変化ではないこと。
関連スコアは Open Targets が複数の datasource から合成した値であり、因果や薬効を表さないこと。

比較図は Heatmap と Dot plot を切り替えられ、疾患 × 細胞の大分類に、閾値以上の遺伝子のうち発現判定が陽性の数と割合を示す。
大分類の展開、色範囲、円の面積、欠測の表示は薬剤ページの標的の図の規則を使う。
ホバーには疾患名、細胞名、実数、割合、分母、未判定の数を示し、薬剤の網羅性の行は出さない。
図を組み立てる関数は薬剤ページと共有し、`kind="gene"` のときは `SummaryRow` の標的側の列だけを読み、ホバーの件数の名前を Genes にする。

詳細欄は、選択欄で疾患を選ぶと更新する。比較図のマスのクリックでは切り替えない。
上から、遺伝子の根拠の表、遺伝子の発現図の順に置く。

根拠の表は AgGrid で、閾値以上の遺伝子 1 件を 1 行にする。

| 列 | 中身 |
| --- | --- |
| Gene | 遺伝子記号と、括弧に入れた Ensembl ID |
| Genetic association score | datatype のスコア。降順を初期の並びにする |
| datasource ごとの列 | datasource ごとのスコア。根拠が無ければダッシュ記号を表示する。列は `snapshot.json` の `datasources` から作り、見出しは ID のアンダースコアを空白に置き換えたものにする |
| Open Targets | その遺伝子と疾患の evidence ページへの `Evidence` リンク |

列見出しのクリックで並べ替え、見出しの下の入力欄で絞り込めるようにする点は薬剤の表と同じにする。
発現判定はこの表に含めない。

発現図は、閾値以上の遺伝子について全元細胞型の発現値を示す。
z スコア、クラスタリング、大分類の平均、白丸の規則は薬剤ページの発現図と同じである。
遺伝子が 100 件を超える疾患では横幅が長くなるので、図内に横スクロールを設ける既存の枠をそのまま使う。
読めるかどうかは、実データで画面を確かめてから判断する。
関節リウマチで閾値 0.5 なら 140 件、0.1 なら 420 件である。

## 説明と用語

画面では、Open Targets の呼び方に合わせて `genetic association score` と書き、`genetic evidence` は区分の名前にだけ使う。
根拠の表の datasource 名は、`snapshot.json` に保存した ID のアンダースコアを空白に置き換えて作る。
値が表示できない理由は、薬剤ページと同じく missing、unresolved、not loaded で呼び分け、閾値以上の遺伝子が無い場合を `no genes at or above the score threshold` とする。

## 初期版の範囲

datasource による絞り込み、上位 N 件の選択、下位語を含めた集計、代表的な根拠の保存、薬剤ページと遺伝子ページを重ねた図は作らない。
薬剤の標的と遺伝学的関連遺伝子の重なりを見る比較は、2 つのページが揃ってから、別に設計する。

## 検証

- `tests/test_refresh.py`：偽の query 関数を渡し、ページングがスコアを持たない遺伝子で止まること、途中の空ページで失敗すること、同じ版の発現を再利用すること、全件成功時だけ保存することを確かめる
- `tests/test_genetics.py`：閾値の等号、分母 0、未判定が残るときの NA、大分類の集約を、薬剤ページの集計と同じ fixture で確かめる
- `tests/test_genetics_ui.py`：`/genetics` の初期レイアウトに遺伝子ページの部品があること、`dcc.Location` の切り替えで `hidden` が入れ替わること、遺伝子ページの `Update` と詳細の callback が応答すること、薬剤ページの既存のテストが変わらず通ることを確かめる
- 実ブラウザーで、ページの往復で設定が保持されること、該当なしの疾患の表示、遺伝子数の多い疾患の発現図を確かめる

実ブラウザーでは、取得元の版 26.09 のデータで次を確かめた。
`/genetics` を直接開いても、上部バーのリンクで切り替えても、比較図は幅いっぱいに描かれる。
薬剤ページで Dot plot に切り替えてから遺伝子ページへ移って戻ると、Dot plot のままである。
関節リウマチは閾値 0.5 で 140 遺伝子、0.3 で 250 遺伝子が分母になり、閾値を変えて `Update` を押してから比較図が変わるまで 2 秒ほどだった。
関節リウマチの発現図は 140 列で、横スクロールで 1 列ずつ読める。
pemphigus のマスは 0 と表示され、ホバーに「No genes at or above the score threshold」が出る。
詳細で該当なしの疾患を選ぶと、表の代わりに案内文が出て、発現図は「No targets in the selected scope」になる。
比較図のマスをクリックして詳細の疾患を切り替える操作は、確かめた結果、実際のブラウザーでも動かなかったので、両ページから外した。

## 根拠

- [Open Targets: association scores](https://platform-docs.opentargets.org/associations)：datatype と datasource のスコアの合成
- [Open Targets: genetic association evidence](https://platform-docs.opentargets.org/evidence#genetic-association)：`genetic_association` に属する datasource
- `autoimmune_atlas/refresh.py`：関連遺伝子と発現の取得と保存
- `autoimmune_atlas/aggregation.py`：共有する発現判定と細胞分類
