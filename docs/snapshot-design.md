# 保存データを 1 ファイルにまとめ取得を一括にする

薬剤の記録、genetic association の関連遺伝子、細胞型別の発現を、1 つの `data/snapshot.json`（schema 3）に保存し、`pixi run refresh` 1 回で取得する。
この文書は、その保存の形と取得の手順の設計である。
遺伝子ページの集計と表示の規則は[遺伝学的関連遺伝子のページの設計](genetics-design.md)に、薬剤ページの規則は[現行の設計](design.md)にあり、どちらもデータの形についてはこの文書を参照する。

遺伝子ページを作ったときは、`snapshot.json`（schema 2）と別に `genetics.json` を置き、関連遺伝子はスコア 0.1 以上だけを保存した。
2 つのファイルの版を突き合わせる手間と、下限のせいで弱い関連が見えないことが分かったので、1 ファイルに戻す。
そのついでに、発現の行が細胞の定数を繰り返している無駄を無くす。

## 下限を設けず、スコアを持つ遺伝子を全部保存する

疾患ごとに、datatype `genetic_association` のスコアを持つ遺伝子を全部保存する。
下限は設けない。
画面の閾値は 0 から 1 で変えられ、初期値 0.5 は変えない。

取得元の版 26.09 で数えたところ、スコアを持つ遺伝子の和集合は約 6,800 件で、薬剤の標的と合わせると約 7,300 件だった。
下限 0.1 のときの約 4,400 件から 6 割ほど増える。
1 疾患あたりの最大は約 1,200 件で、スコアを持つ遺伝子が 1 件も無い疾患が半数近くある。
この数は版が変わると変わる。

## 細胞の定数は 1 回だけ持つ

schema 2 の発現の行は、遺伝子 × 細胞ごとに `cell` 名、`parent`、`parent_id`、`ancestor_ids` を繰り返していた。
170 細胞について、これらの値は遺伝子が違っても同じである。
版 26.09 の実データで、全行を照らし合わせて食い違いが無いことを確かめた。

schema 3 では、細胞の定数を `cells` の表に 1 回だけ置き、発現の行は細胞 ID と数値だけにする。

```text
schema: 3
roots, data_version, retrieved_at, source: schema 2 と同じ
diseases: schema 2 と同じ
records: 薬剤の記録。schema 2 と同じ
associations: {疾患 ID: [{target_id, target, score, datasource_scores}, ...]}
datasources: genetic association に属する datasource の ID の一覧
cells: {細胞 ID: {name, parent_id, parent, ancestor_ids}}
expression: {遺伝子 ID: [{cell_id, median, specificity_score}, ...]}
```

`associations` の 1 件は、`genetics.json` の schema 1 と同じ形である。
スコアの降順に並べて保存する。

`cells` には `parent` の表示名も入れる。
大分類の名前は `parent` から作っており、`parent_id` だけでは同じ表示にならない。
取得のときに、同じ細胞 ID の名前、親、祖先が遺伝子のあいだで一致することを確かめ、食い違えば失敗にする。
schema 2 の集計が行ごとに行っていた検査を、保存の時点に移す。

発現の行は名前付きの dict のままにする。
配列にすればさらに小さくなるが、読むときにキーの順序を覚える必要がある。
7,300 遺伝子 × 170 細胞で 1 行 60 バイト程度とみて、ファイルは 100 MB 前後になる見込みである。
この見込みは試作のファイルを読み込んで測り、「検証」の節に実測を書く。

## 集計は細胞の表を読む

`aggregation._cell_memberships()` は、発現の行ではなく `snapshot["cells"]` から細胞の名前、大分類、祖先を作る。
T 細胞の扱い（`CL_0000084` を祖先に持つ細胞を T 細胞にまとめる）は変えない。

`expression_metadata()` が返す行には、`cells` から引いた `cell` 名を付ける。
薬剤ページの根拠の行と詳細の表示が、metadata の `cell` 名を読んでいるためである。
保存の形と実行時の型は分ける。
`ExpressionRow`（保存）は 3 キーになり、`ExpressionMetadata`（実行時）は `cell` と `target_median` を持つ。

`load_snapshot()` は schema 3 だけを受け付ける。
schema 2 と `genetics.json` の読み込みは残さない。
`data/` は git の追跡外で、環境ごとに取り直す規約なので、アプリの側に移行の処理は書かない。
schema 2 のファイルがあれば、これまでどおり `pixi run refresh` を案内する。

fixture（`tests/core_fixture.py`、`tests/test_ui.py`、`tests/genetics_fixture.py`）は schema 3 の形に書き換える。
`genetics_fixture.py` は独立したファイルではなくなり、`associations` と `datasources` を持つ schema 3 の fixture に統合する。

## 取得は 1 本にし、同じ版の発現は取り直さない

`autoimmune_atlas/refresh.py` の `main()` を 1 本にする。

1. 版を取り、疾患の範囲を解決する
2. 疾患ごとに薬剤と作用機序を取る（schema 2 と同じ）
3. 疾患ごとに関連遺伝子を取る。`orderByScore: "genetic_association"` で並べ、genetic association のスコアを持たない遺伝子に当たった時点でその疾患のページングを止める
4. 現れた datasource ごとに根拠を 1 件取り、`genetic_association` に属するものだけを残す（`genetics.json` のときと同じ）
5. 薬剤の標的と関連遺伝子の和集合について、発現を取る
6. 終了時に版をもう一度取り、開始時と違えば失敗にする
7. 全件成功したときだけ、一時ファイルから `snapshot.json` を置き換える

**同じ版の発現は取り直さない。**
取得を始めるとき、既存の `data/snapshot.json` と `data/genetics.json` を読み、`data_version` が今回と同じなら、そこにある遺伝子の発現の行を再利用する。
schema 2 の行からは `cell_id`、`median`、`specificity_score` だけを取り出す。
版が違えば全部取り直す。
`genetics.json` を読むのは、schema 3 への切り替えのときに約 4,000 遺伝子ぶんの取得を省くためで、切り替えが済んだら消してよい。

再利用が効くのは、完成した保存データにある遺伝子だけである。
途中で失敗した取得は何も保存しないので、初回の長い取得を途中から再開することはできない。
再開が要るほど失敗が続くなら、取得済みの遺伝子を一時ファイルに書き出す形を別に考える。

取得の時間は、発現が 4,148 遺伝子で約 55 分（2 スレッド）だった実績から、和集合 7,300 遺伝子を全部取ると 1.5 時間から 2 時間とみる。
既存の 2 ファイルから再利用できれば、初回は約 2,400 遺伝子ぶんの 30 分程度で済む。
並列数を増やせば短くなるが、API のエラー率は測っていないので、2 のままにする。

`autoimmune_atlas/refresh_genetics.py` と `pixi run refresh-genetics` は消す。
`genetics.py` の `genes_for_disease()` と `summarize_genes()` は `snapshot["associations"]` を読む。
`create_app()` の `genetics` と `genetics_error` の引数、版の一致の確認、`genetics.json` が無いときの案内は消す。
遺伝子ページが出ないのは、`snapshot.json` そのものが無いか読めないときだけになる。

## 画面の変更

- 閾値の入力欄の `min` を 0 にし、`step` を 0.01 にする。0.05 のままだと 0.02 のような値を入力できない
- 入力値の検証も 0 以上 1 以下にする。`SCORE_FLOOR` の定数は消す
- 概要の「Gene–disease associations (score ≥ 0.1)」を「Gene–disease associations」にし、閾値の ⓘ から「0.1 未満は保存しない」の文を消す
- `genetics.json` の取得を促す画面と文言を消す

集計と表示の規則は変えない。

## 初期版の範囲

途中から再開できるチェックポイント、取得の並列数の調整、JSON 以外の保存形式は作らない。
ファイルの分割もしない。
節ごとに互いの ID で参照する形にしてあるので、一部だけ取り直す必要が出たら、そのときに節を別ファイルへ切り出す。

## 検証

- `tests/test_refresh.py`：関連遺伝子のページングがスコアの無い遺伝子で止まること、`cells` の一貫性の検査、同じ版の発現の再利用と版が違うときの全件取得、全件成功時だけ保存することを、偽の query 関数で確かめる
- `tests/test_aggregation.py`：`cells` の表から作った大分類と祖先が、schema 2 の fixture で得ていた結果と同じであること
- `tests/test_genetics.py`、`tests/test_genetics_ui.py`：schema 3 の fixture で既存のテストが通ること。閾値 0 と 0.02 が受け付けられること
- 実データで、schema 3 の `snapshot.json` の大きさ、起動時のメモリ、起動時間を測り、ここに書く
- 実ブラウザーで、閾値 0 のときに関節リウマチの分母がスコアを持つ遺伝子の数（版 26.09 で 697）になること

版 26.09 で測った結果は次のとおりである。
`snapshot.json` は 96 MB で、schema 2 の 2 ファイル（`snapshot.json` 101 MB と `genetics.json` 521 MB）の 6 分の 1 になった。
読み込みに 0.2 秒、画面の構成に 1.1 秒かかり、プロセスのメモリは 0.97 GB だった。
schema 2 のときは、2 つのファイルを読むとメモリが約 3 GB、画面の構成を作り始めるまでに約 9 秒かかっていた。
再利用が働き、7,299 遺伝子のうち 4,913 遺伝子を既存の 2 ファイルから取り、取得したのは 2,386 遺伝子だった。
関連遺伝子と発現の取得を合わせて 18 分で終わった。
genetic association に属する datasource は 5 つ（`eva`、`gene_burden`、`gwas_credible_sets`、`orphanet`、`uniprot_variants`）だった。
関節リウマチの閾値 0 の分母は 697 で、`associations` の配列の長さと一致した。
閾値 0.02 でも同じ 697 で、入力欄は 0.01 刻みの値を受け付けた。

## 根拠

- [Open Targets: association scores](https://platform-docs.opentargets.org/associations)：genetic association の datatype スコアと datasource スコア
- [Open Targets: baseline expression](https://github.com/opentargets/platform-docs/blob/main/target/baseline-expression.md)：細胞型別 pseudobulk と CELLEX
- `autoimmune_atlas/refresh.py`：取得と保存、`autoimmune_atlas/aggregation.py`：細胞の表の読み込みと集計
