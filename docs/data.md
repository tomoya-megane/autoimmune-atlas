# 保存データの形と取得の手順

アプリが読むデータは、`data/snapshot.json`（schema 3）の 1 ファイルである。
このファイルは git で追跡せず、環境ごとに `pixi run refresh` で Open Targets の公開 API から取得する。
アプリは起動時に `load_snapshot()` でこのファイルを読み、実行中は取り直さない。
この文書は、ファイルが持つキーの全部と、`pixi run refresh` がそれを作る手順を述べる。
集計の規則は[現行の設計](design.md)に、遺伝子ページの規則は[遺伝学的関連遺伝子のページの設計](genetics-design.md)に、疾患の起点は[対象疾患の起点](disease-roots.md)にある。

## ファイルは 12 のキーを持つ

最上位のキーは次の 12 個で、型は `autoimmune_atlas/models.py` の `Snapshot` が定める。

| キー | 型 | 中身 |
| --- | --- | --- |
| `schema` | `int` | 保存の形の版。`refresh` は `3` を書く |
| `root` | `str` か `RootIdentity` | 最初の起点の疾患 ID。`refresh` は `SCOPE_ROOTS` の先頭（`MONDO_0007179`）を文字列で書く。型は `{id, name}` の形も受け付ける |
| `roots` | `list[SnapshotRoot]` | 起点ごとの `id`、`name`、`include_descendants`、`count`。`count` は、その起点から数えた疾患の数（下位語を含めるなら下位語と起点の合計）である |
| `data_version` | `DataVersion` | Open Targets のデータの版。`year` と `month` は文字列、`iteration` は文字列か整数か `null` |
| `retrieved_at` | `str` | 保存した時刻。UTC の ISO 8601 形式 |
| `source` | `str` | 取得元の URL（`https://api.platform.opentargets.org/api/v4/graphql`） |
| `diseases` | `list[Disease]` | 対象疾患。名前の順に並ぶ |
| `records` | `list[DrugRecord]` | 疾患、薬剤、標的の組ごとの薬剤の記録 |
| `associations` | `dict[str, list[GeneAssociation]]` | 疾患 ID ごとの、genetic association の関連遺伝子 |
| `datasources` | `list[str]` | genetic association に属する datasource の ID |
| `cells` | `dict[str, CellDefinition]` | 細胞 ID ごとの名前、親分類、祖先 |
| `expression` | `dict[str, list[ExpressionRow]]` | 遺伝子 ID ごとの、細胞型別の発現の行 |

`roots` と `data_version` は、型の上では `NotRequired` である。
`refresh` は、この 2 つを含む 12 個を必ず書く。

## 疾患と薬剤の記録は Open Targets の値をそのまま持つ

`diseases` の 1 件は、`Disease` の 5 つのキーを持つ。

| キー | 中身 |
| --- | --- |
| `id` | 疾患 ID（`MONDO_0008383` など） |
| `name` | Open Targets の疾患名。書き換えない |
| `parent_ids` | 直接の親の疾患 ID。Open Targets の `parents` をそのまま写す。疾患選択の入れ子に使う |
| `status` | `refresh` は常に `ready` を書く。集計は `ready` 以外の疾患を未取得（not loaded）として扱う |
| `unclassified_stages` | その疾患の臨床候補のうち、`maxClinicalStage` が `UNKNOWN` か `WITHDRAWAL` の行の数。段階で絞る前に数える |

`parent_ids` と `unclassified_stages` は型の上では `NotRequired` だが、`refresh` は必ず書く。

`records` の 1 件は、疾患、元の薬剤、標的の組に対応する。
キーは `DrugRecord` の 14 個である。

| キー | 中身 |
| --- | --- |
| `disease_id`、`disease` | 疾患 ID と疾患名 |
| `drug_id`、`drug` | 元の薬剤の ChEMBL ID と名前。名前は小文字にする |
| `canonical_drug_id`、`canonical_drug` | 同じ有効成分の ID と名前。`parentMolecule` があればその値、無ければ元の薬剤の値。名前は小文字にする |
| `drug_type` | Open Targets の `drugType`（`Small molecule` など） |
| `modality` | `drug_type` から決めた区分（`small_molecule` など）。集計は読み込み時に `drug_type` から作り直す |
| `stage` | その疾患でその薬剤が到達した最高段階（`maxClinicalStage`） |
| `target_id`、`target` | 標的の Ensembl ID と遺伝子記号 |
| `mechanism` | 標的ごとの作用機序。複数あれば名前の順に並べ、`;` と空白 1 つでつなぐ |
| `action_types` | 作用の型（`INHIBITOR` など）。重複を除き、名前の順に並べる |
| `references` | 作用機序の出典。1 件は `source`、`ids`、`urls` の 3 キーで、同じ組は 1 回だけ持つ。`source` が無ければ空文字列にする |

`modality` と `action_types` は型の上では `NotRequired` だが、`refresh` は必ず書く。

作用機序の標的を持たない薬剤も、記録から外さない。
`target_id` を空文字列、`target` を `Unknown`、`mechanism` を空文字列、`action_types` と `references` を空の配列にした 1 行を置く。
型は `target_id`、`target`、`mechanism` に `null` も許すが、`refresh` は `null` を書かない。

保存するのは、段階が Phase I 以降（`PHASE_1`、`PHASE_1_2`、`PHASE_2`、`PHASE_2_3`、`PHASE_3`、`PREAPPROVAL`、`APPROVAL`、`PHASE_4`）の行だけである。
`UNKNOWN`、`WITHDRAWAL`、`EARLY_PHASE_1`、`IND`、`PRECLINICAL`、`PHASE_0` は既知の段階として読み、保存しない。
この 14 個のどれでもない段階が現れたら、取得を失敗にする。
Phase II 以降などの絞り込みは、集計が `stage` を読んで行う。

同じ疾患の中で同じ薬剤 ID が 2 回現れたとき、段階を満たす行に薬剤の情報が無いとき、`drugType` が既知の 11 種のどれでもないときも、取得を失敗にする。
同じ有効成分の記録をまとめるのは集計の仕事で、保存の時点では元の薬剤ごとに行を分けたままにする。

## 関連遺伝子は下限なしでスコアの降順に持つ

`associations` は、`diseases` にある全部の疾患 ID をキーに持つ。
genetic association のスコアを持つ遺伝子が 1 件も無い疾患は、空の配列になる。
キーが無い疾患は無いので、空の配列は「取得したが該当なし」を意味する。

1 件の `GeneAssociation` は 4 つのキーを持つ。

| キー | 中身 |
| --- | --- |
| `target_id` | 遺伝子の Ensembl ID |
| `target` | 遺伝子記号 |
| `score` | datatype `genetic_association` のスコア |
| `datasource_scores` | datasource ID ごとのスコア。`datasources` にある datasource だけを残す |

スコアを持つ遺伝子は全部保存し、下限を設けない。
遺伝子ページの閾値は 0 から 1 で変えられるので、どの閾値でも分母を保存データから数えられる。
並びはスコアの降順で、同じスコアなら `target_id` の順である。
API も降順で返すが、保存の形として並びを保証するために、保存の前に並べ直す。

`datasources` は、genetic association に属する datasource の ID を名前の順に並べたものである。
一覧はコードに書かず、取得のたびに決める。
版が上がって datasource が増えても、取りこぼさないためである。
決め方は[取得の手順](#取得は-1-本で全件成功したときだけ保存する)にある。
版 26.09 では、`eva`、`gene_burden`、`gwas_credible_sets`、`orphanet`、`uniprot_variants` の 5 つだった。

## 細胞の定数は 1 回だけ持ち、発現の行は 3 キーである

細胞の名前、親分類、祖先は、遺伝子が違っても同じ値になる。
そのため、発現の行には持たせず、`cells` に細胞ごとに 1 回だけ置く。

`cells` の値の `CellDefinition` は 4 つのキーを持つ。

| キー | 中身 |
| --- | --- |
| `name` | Tabula Sapiens の細胞型の名前（`biosampleName`） |
| `parent_id` | 取得元の親分類の ID。無ければ `null` |
| `parent` | 取得元の親分類の名前。無ければ `null` |
| `ancestor_ids` | Cell Ontology の祖先の ID。無ければ空の配列 |

`parent` は、`parent_id` から作れないので保存する。
集計は大分類の表示名を `parent` から作る。

`expression` の値は、1 遺伝子の細胞型別の発現の行の配列で、1 行の `ExpressionRow` は 3 つのキーを持つ。

| キー | 中身 |
| --- | --- |
| `cell_id` | 細胞 ID。`cells` のキーを指す |
| `median` | ドナー間中央値。単位は `CPM(pseudobulk sum[counts])`。欠測なら `null` |
| `specificity_score` | CELLEX の特異性スコア（0 から 1）。欠測なら `null` |

`expression` は、薬剤の標的と関連遺伝子の和集合の全部をキーに持つ。
Tabula Sapiens の細胞型別の値が無い遺伝子は、空の配列になる。

`cells` と `expression` のあいだには、2 つの不変条件がある。

- `expression` の全部の `cell_id` が、`cells` のキーにある
- `cells` は、保存した `expression` が使う細胞だけを持つ。細胞 ID の順に並ぶ

集計はこの 2 つを前提にしている。
`aggregation._cell_memberships()` は `cells` だけから細胞の名前、大分類、祖先を作る。
`CL_0000084` 自身と、それを `ancestor_ids` に持つ細胞は、取得元の親分類に関わらず T 細胞の大分類にまとめる。
`expression_metadata()` は、`cells` に無い `cell_id` を見つけると失敗する。
Target-relative median の基準になる標的内中央値は、その遺伝子の行が `cells` の全細胞を覆い、`median` に `null` が無いときだけ計算し、それ以外は `null` にする。
`cells` に使われない細胞が残っていると、全部の遺伝子で標的内中央値が `null` になる。
2 つ目の不変条件は、これを防ぐためにある。

発現の行は、名前付きの dict のままにする。
配列にすればファイルはさらに小さくなるが、読むときにキーの順序を覚えておく必要がある。
dict のままでも 96 MB に収まったので、読みやすさを優先する。

## 読み込みは msgspec で型ごと検証する

`autoimmune_atlas/snapshot.py` の `load_snapshot()` は、`msgspec.json.decode()` で `Snapshot` の型を指定してファイルを読む。
入れ子の `Disease`、`DrugRecord`、`GeneAssociation`、`CellDefinition`、`ExpressionRow` まで、キーと値の型をこの時点で確かめる。
結果は次の 4 通りである。

| ファイルの状態 | 結果 |
| --- | --- |
| 無い | `None` を返す |
| JSON として読めない | `ValueError("snapshot.json is not valid JSON: ...")` |
| 型に合わない | `ValueError("Unsupported snapshot format (...). Refresh it with: pixi run refresh")` |
| 型に合うが `schema` が 3 でない | `ValueError("Unsupported snapshot schema N. Refresh it with: pixi run refresh")` |

schema 1 と 2 のファイルは `cells` などのキーを持たないので、型に合わない場合にあたる。
古い形を読み替える処理は、アプリに置かない。
`data/` は環境ごとに取り直す規約なので、`pixi run refresh` を 1 回打てば済む。

`app.py` は、`None` のときも `ValueError` のときも、データ無しでアプリを組み立てる。
薬剤ページには `pixi run refresh` を打つ案内を表示する。
見出しは、ファイルが無ければ「No data loaded yet」、読めなければ「Data refresh required」で、後者ではエラーの文面も表示する。
遺伝子ページには「Genetic associations are not available」と、`pixi run refresh` のあとにアプリを再起動する案内を表示する。

## 取得は 1 本で、全件成功したときだけ保存する

`pixi run refresh` は `autoimmune_atlas/refresh.py` の `main()` を実行する。
`main()` は次の順に進み、どこかで失敗すれば何も保存せずに止まる。

1. `meta.dataVersion` でデータの版を取る
2. `SCOPE_ROOTS` の起点ごとに疾患と下位語を取り、和集合を対象疾患にする。起点が 1 つでも見つからなければ失敗にする
3. 対象疾患を 5 件ずつ `diseases(efoIds:)` で照会し、疾患名、親 ID、臨床候補を取る。返った疾患の集合が照会した 5 件と一致しないとき、臨床候補の行数が `count` と一致しないときは失敗にする。臨床候補は[前の節](#疾患と薬剤の記録は-open-targets-の値をそのまま持つ)の規則で段階を選ぶ
4. 段階を満たした薬剤を 20 件ずつ `drugs(chemblIds:)` で照会し、作用機序、作用の型、標的、出典を取る。返った薬剤の集合が照会した 20 件と一致しなければ失敗にする。ここで `records` を組み立てる
5. 疾患ごとに `associatedTargets(enableIndirect: false, orderByScore: "genetic_association")` を 500 件（`PAGE_SIZE`）ずつ取る。datatype `genetic_association` のスコアを持たない遺伝子に当たった時点で、その疾患のページングを止める。総件数に達する前に空のページが返れば、失敗にする
6. 関連遺伝子の `datasourceScores` に現れた datasource ごとに、それを持つ疾患と遺伝子の組を 1 つ選び、`evidences(size: 1)` で根拠を 1 件取る。その `datatypeId` が `genetic_association` の datasource だけを `datasources` に残し、`datasource_scores` もそれに絞る。根拠が 1 件も取れない datasource があれば失敗にする
7. 薬剤の標的（空文字列を除く）と関連遺伝子の和集合を、発現を取る遺伝子にする
8. 同じ版の既存ファイルから、発現の行と細胞の定数を再利用する。規則は[次の節](#同じ版の発現は既存のファイルから再利用する)にある
9. 残りの遺伝子の発現を、2 スレッド（`max_workers=2`）で取る。1 遺伝子につき `baselineExpression` を 3,000 行ずつ全ページ取り、件数が `count` と一致しなければ失敗にする
10. 発現に現れた細胞 ID が全部 `cells` にあることを確かめ、`cells` をその細胞だけに絞る
11. 版をもう一度取り、手順 1 と違えば「取得中にデータの版が変わりました。再取得してください」で失敗にする
12. `cell_catalog()` を 1 回呼び、同じ親分類 ID に別の名前が付いていないかを確かめる
13. 一時ファイルに書いてから `data/snapshot.json` を置き換える

手順 5 のページングは、`orderByScore: "genetic_association"` の並びに頼っている。
この並びでは、genetic association のスコアを持つ遺伝子が降順に先に並び、スコアを持たない遺伝子がその後ろに続く。
関節リウマチで、全件を取って数えた件数と、並べて途中で止めて数えた件数が一致することを確かめた。
HTTP と GraphQL が成功していても、ページが途中で空になることがあるので、空のページは失敗として扱う。

手順 6 で datasource の区分を API の一覧から取らないのは、`associationDatasources` が版 26.09 で空の配列を返したためである。

発現の行は、Tabula Sapiens（datasource `tabula_sapiens`）の細胞型別 pseudobulk だけを残す。
組織を持つ行と、細胞型を持たない行は捨てる。
次のどれかに当たれば、取得を失敗にする。

- 単位が `CPM(pseudobulk sum[counts])` でない
- `median` が `null` でも 0 以上の有限の数でもない
- `specificity_score` が `null` でも 0 以上 1 以下の有限の数でもない
- 1 遺伝子の中で同じ細胞 ID の行が 2 つある
- 同じ細胞 ID の名前、親分類、祖先が、遺伝子のあいだで食い違う（`merge_cells()`）

API への照会は、すべて `query_api()` を通る。
`httpx` でタイムアウト 45 秒の POST を送り、proxy は環境変数（`HTTPS_PROXY`、`https_proxy`）から `httpx` が読む。
HTTP の状態が 200 でないとき、応答に `errors` があるとき、`data` が無いときは失敗とみなし、1 秒、2 秒と待って最大 3 回まで試す。
3 回目も失敗すれば、例外をそのまま上げる。
呼び出し元は応答を `msgspec.convert()` で照会ごとの型へ変換するので、キーの欠けや型の違いもここで失敗になる。

`save_snapshot()` は、`json.dumps(allow_nan=False)` で書く。
`NaN` や `Infinity` が混ざっていれば、書く前に失敗する。
一時ファイルは `data/` の中に作り、書き終えてから `Path.replace()` で置き換える。
途中で失敗したときは一時ファイルを消すので、前回の `snapshot.json` はそのまま残る。

`main()` は、起点ごとの疾患数、チャンクごとの進み、関連遺伝子を持つ疾患の数、genetic association の datasource、発現の再利用と取得の件数を表示する。
最後に、保存先と、疾患、薬剤、標的、関連遺伝子の件数を表示する。

## 同じ版の発現は既存のファイルから再利用する

発現の取得は、全体の時間の大半を占める。
そこで `reusable_expression()` が、`data/snapshot.json` と `data/genetics.json` をこの順に読み、`data_version` が今回と同じなら発現の行を再利用する。
版が違えば、そのファイルは使わない。

- ファイルが無い、JSON として読めない、緩い型（`_StoredFile`）に合わない、版が違う、のどれかに当たれば、そのファイルを飛ばす
- schema 3 の `snapshot.json` は、`cells` と 3 キーの行をそのまま読む
- schema 2 の `snapshot.json` と schema 1 の `genetics.json` の行は、`cell`、`parent_id`、`parent`、`ancestor_ids` を持つ。そこから細胞の定数を `cells` へ移し、行は `cell_id`、`median`、`specificity_score` の 3 キーにする
- 同じ遺伝子が両方のファイルにあれば、先に読んだファイルの行を使う
- 細胞の定数は、後のファイルにしか無い遺伝子の行も含めて、全部の行から集めて比べる。食い違えば取得を失敗にする

`load_snapshot()` は schema 3 しか受け付けないので、再利用では別の緩い型で読む。
再利用した遺伝子のうち、今回の和集合に無いものは保存しない。

`genetics.json` を読むのは、schema 2 から schema 3 へ切り替える最初の取得で、約 4,000 遺伝子ぶんの取得を省くためである。
schema 3 の `snapshot.json` が一度できれば、`data/genetics.json` は消してよい。

再利用できるのは、完成した保存データにある遺伝子だけである。
途中で失敗した取得は何も保存しないので、長い取得を途中から再開する手段は無い。
再開が要るほど失敗が続くなら、取得済みの遺伝子を一時ファイルに書き出す形を別に考える。

## 実測

取得元の版 26.09 で、schema 3 の `snapshot.json` を作って測った。

| 項目 | 値 |
| --- | --- |
| ファイルの大きさ | 96 MB |
| 読み込み | 0.2 秒 |
| 画面の構成 | 1.1 秒 |
| プロセスのメモリ | 0.97 GB |
| 発現を持つ遺伝子 | 7,299。うち 4,913 を既存の 2 ファイルから再利用し、2,386 を取得した |
| 関連遺伝子と発現の取得時間 | 18 分 |
| 関節リウマチの関連遺伝子 | 697（閾値 0 と 0.02 のどちらでも分母が 697 で、`associations` の配列の長さと一致した） |
| genetic association の datasource | 5 つ |

schema 2 のときは、`snapshot.json` が 101 MB、`genetics.json` が 521 MB だった。
`genetics.json` は `snapshot.json` のおよそ 5 倍で、発現の行ごとに `cell`、`parent`、`ancestor_ids` を繰り返していたためである。
2 つのファイルを読むとプロセスのメモリは約 3 GB まで増え、画面の構成を作り始めるまでに約 9 秒かかっていた。
schema 3 の 1 ファイルは、2 ファイルの合計の 6 分の 1 になった。

発現の取得時間は、schema 2 のときに 4,148 遺伝子で約 55 分（2 スレッド）だった。
この実績から、再利用できるファイルが無い環境で約 7,300 遺伝子を全部取ると、1.5 時間から 2 時間かかると見積もっている。
並列数を増やせば短くなるが、API のエラー率を測っていないので、2 のままにしている。

これらの値は、版が変わると変わる。

## 経緯

遺伝子ページを作ったとき、`snapshot.json`（schema 2）とは別に `data/genetics.json` を置き、関連遺伝子はスコア 0.1 以上だけを保存した。
使ってみて、3 つの問題が分かった。
2 つのファイルの版が一致しているかを、読み込みのたびに確かめる必要があった。
下限 0.1 のせいで、それより弱い関連が画面から見えなかった。
発現の行が細胞の定数を繰り返し、2 つのファイルを読むとメモリが約 3 GB になった。

schema 3 で、2 つのファイルを 1 つにまとめ、下限を無くし、細胞の定数を `cells` に 1 回だけ置いた。
版 26.09 では、スコアを持つ遺伝子の和集合は約 6,800 件で、薬剤の標的と合わせると約 7,300 件になった。
下限 0.1 のときの約 4,400 件から、6 割ほど増えた。
`autoimmune_atlas/refresh_genetics.py` と `pixi run refresh-genetics` は、このときに消した。

初期版では、次のものを作らない。

- 途中から再開できるチェックポイント
- 取得の並列数を変える設定
- JSON 以外の保存形式
- ファイルの分割

節どうしは互いの ID で参照する形にしてあるので、一部だけ取り直す必要が出たら、そのときに節を別ファイルへ切り出す。

## 根拠

- [Open Targets: association scores](https://platform-docs.opentargets.org/associations)：genetic association の datatype スコアと datasource スコア
- [Open Targets: baseline expression](https://github.com/opentargets/platform-docs/blob/main/target/baseline-expression.md)：細胞型別 pseudobulk と CELLEX
- `autoimmune_atlas/models.py`：保存の形の型
- `autoimmune_atlas/refresh.py`：取得、検証、再利用、保存
- `autoimmune_atlas/snapshot.py`：読み込みと検証
