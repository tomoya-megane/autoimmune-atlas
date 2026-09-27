# Autoimmune Atlas

自己免疫疾患ごとに、臨床開発された薬剤の標的遺伝子がどの細胞型で発現するかを比較するローカルアプリ。
同じ比較を Open Targets の genetic association で疾患と結びついた遺伝子について行う、遺伝学的関連遺伝子のページ（`/genetics`）も持つ。
Python の Dash と Plotly を使う。
実行環境には Python 3.13 を使う。
集計規則と判断理由は[設計ノート](docs/design.md)に記録する。
文書の一覧は[ドキュメント索引](docs/README.md)にある。

画面構成は Open Targets の [Associations on the Fly](https://platform-docs.opentargets.org/web-interface/associations-on-the-fly) を参考にしている。
小さな概要欄、比較表の近くのフィルター、根拠を表示する詳細欄を採用した。
配色と余白はこのアプリ向けの調整であり、Open Targets のテーマをそのまま移植したものではない。

## コード構成

起動の入口 `app.py` だけをリポジトリ直下に置き、処理の実体は `autoimmune_atlas/` パッケージにまとめている。
データ更新は `pixi run refresh` が `autoimmune_atlas.refresh` をモジュールとして実行する。
集計、データ取得、画面表示のどこを読むべきかを、ファイル名から判断できる構成である。

```text
app.py                         Dash アプリの起動
autoimmune_atlas/
├── aggregation.py             薬剤、標的、細胞型の集計
├── disease_catalog.py         対象疾患と表示順の定義
├── genetics.py                遺伝学的関連遺伝子の集計
├── models.py                  スナップショットと集計結果の共有データ型
├── refresh.py                 公開 API の取得、正規化、保存（python -m で実行）
├── snapshot.py                保存済みデータの読み込みと検証
└── ui/
    ├── application.py         Dash アプリの組み立てとページの切り替え
    ├── components.py          再利用する画面部品と上部バー
    ├── config.py              画面で共有する固定値
    ├── controls.py            両ページの callback が共有する入力値の検証と行の切り替え
    ├── drugs_callbacks.py     薬剤ページの画面操作への応答
    ├── drugs_layout.py        薬剤ページの初期画面と詳細欄の配置
    ├── figures.py             比較図と発現図
    ├── genetics_callbacks.py  遺伝子ページの画面操作への応答
    └── genetics_layout.py     遺伝子ページの初期画面と詳細欄の配置
assets/                        CSS、JavaScript、アイコン
├── help.js                    ツールチップの位置調整と Escape キー
└── router.js                  ページ切り替え後の図の再描画
data/                          取得済みスナップショット
docs/                          設計ノートと対象疾患の根拠
tests/
├── test_aggregation.py        集計規則
├── test_refresh.py            公開 API の取得と保存
├── test_disease_catalog.py    対象疾患と表示順
├── test_genetics.py           遺伝子ページの集計規則
├── test_genetics_ui.py        遺伝子ページの図、画面、callback
├── ui_fixture.py              UI のテストが共有する型とスナップショット
├── test_disease_tree.py       疾患の選択ツリー
├── test_figures.py            ヒートマップとドットプロットの図
├── test_evidence.py           根拠の行と詳細パネル
├── test_inputs.py             閾値の入力と初期選択
├── test_snapshot.py           スナップショットの読み込み
├── test_application.py        アプリの入口と assets のパス
├── test_callbacks.py          画面と callback
└── test_help.cjs              ツールチップの JavaScript
```

## 開発環境

clone 後に、このリポジトリだけで環境を構築できる。

```bash
pixi install
```

作業の完了前に、整形、lint、型検査、Python と JavaScript のテストをまとめて実行する。

```bash
pixi run check
```

すべてのタスクが成功した状態を、変更完了の基準とする。
共通の検査基準は、それぞれのリポジトリ内に保持し、隣のリポジトリを実行時に参照しない。
基準を変えるときは、もう一方にも適用できる変更かを確認し、各リポジトリで個別に更新して検査する。
もう一方を参照できない場合は、こちらの作業を進め、未反映の項目を変更の説明に残す。

| 共通に保つもの | このリポジトリの保存先 |
| --- | --- |
| Ruff の基本ルール、Python の対象版、整形幅 | `ruff.toml` |
| basedpyright の診断設定 | `pyrightconfig.json` |
| Markdown のルール | `.markdownlint-cli2.jsonc` の `config` |

型検査の共通項目は `typeCheckingMode = recommended`、`reportUnusedCallResult = false`、`reportImplicitStringConcatenation = false` である。
仮想環境の指定、検索パス、除外対象は各リポジトリで定める。
Markdown を変更したときは `pixi run lint-markdown` で点検する。

実行環境、設定、基本点検はこのリポジトリだけで完結する。
スキルは複製せず、`../../mycompany/` を参照できるときだけ、作業に該当する `SKILL.md` を読んで使う。
参照できない場合も起動と基本点検は行えるが、スキルが必要で代わりの方法がない作業には制約が残る。
その場合は、利用できないスキルと制約を報告する。

## 起動

`data/` は大容量のため Git で追跡しない。
別環境で clone した後は、このディレクトリで公開 API からデータを取得する。

```bash
pixi run refresh
```

取得が完了してからアプリを起動する。

```bash
pixi run start
```

[アプリを開く](http://127.0.0.1:8050)。
起動時は保存済みデータを使い、外部 API を呼ばない。
停止は起動したターミナルで `Ctrl+C`。

コードを編集するときは開発モードで起動する。
通常モードが起動中なら、先に `Ctrl+C` で停止する。

```bash
pixi run dev
```

Python ファイルを保存するとアプリとブラウザーが再読み込みされ、`assets/` 内の CSS を保存するとスタイルが更新される。
Python の変更による再読み込みでは、選択したフィルターは初期値に戻る。
[Dash の開発モード](https://dash.plotly.com/devtools)を使用している。

## 集計の意味

対象は Open Targets の autoimmune disease (`MONDO_0007179`) と、[起点の一覧](docs/disease-roots.md)に挙げた疾患の下位語である。
MONDO が自己免疫疾患の枝の外に置く疾患（筋炎、血管炎、腎炎、炎症性腸疾患など）と、自己炎症性疾患を起点として足している。
疾患の階層に従うため、親概念や疾患サブタイプが混在する。
`Browse disease groups` から疾患群と疾患ファミリーを開き、比較する疾患を個別に選べる。
疾患本体のチェック欄の下にある `details` を開くと、病型や関連用語を個別に選べる。
名前で検索する選択欄とチェック欄は同期し、ヒートマップの列も群とファミリーの名称順に並ぶ。
親疾患のチェックはその疾患だけを選び、子疾患の選択や集計値の合算は行わない。
分類先が確認できない用語も `Other / unclassified` から選べる。
臨床段階は Phase I 以降、Phase II 以降、Phase III 以降、承認到達から選ぶ。
初期表示では、臨床段階に Phase III 以降、モダリティに All、発現基準に Fixed CPM + CELLEX specificity、表示値に Percent を使う。
各疾患の `drugAndClinicalCandidates.maxClinicalStage` を使い、薬剤全体の最高段階は使わない。
Phase I/II は Phase II 以降に、Phase II/III は Phase III 以降に含めない。
承認到達には歴史的な記録を含み、現在の販売継続や承認の有効性を保証しない。
症状や併存症に対する治療も含め、治療目的での分類は行わない。
[臨床段階の定義](https://platform-docs.opentargets.org/drug/clinical-report)を参照。

上部パネルで条件を変更し、`Update` を押すと結果へまとめて反映する。
入力中はヒートマップと詳細を保持し、未反映の変更があることをボタンの隣に表示する。

比較図は `Distinct targets`（標的数）と `Canonical drugs`（薬剤数）を選び、`Update` で切り替える。
初期表示は `Distinct targets` とし、切り替えても疾患・細胞の並びとフィルターを保つ。
図の形式は Heatmap と Dot plot から選び、`Update` を待たずに切り替える。
Heatmap は選択した Percent または Count を色で示す。
Dot plot は Percent を色、Count を圧縮した円の面積で同時に示すため、Measure の選択を使わない。
標的数は Ensembl 遺伝子 ID の重複を除いた数である。
薬剤数は `parentMolecule.id` で確認できる塩や水和物を同じ成分としてまとめ、一つでも標的が条件を満たせば一剤と数える。
元の薬剤 ID と成分 ID は詳細に残す。
標的割合の分母は、その疾患で段階とモダリティの条件を満たす薬剤の既知標的数である。
薬剤割合の分母は、少なくとも一つの既知標的を持つ成分数である。
どちらも細胞型ごとに分母を変えず、標的が判明した薬剤数と全薬剤数を併記する。
複合体の構成遺伝子はそれぞれ数えるため、標的複合体や作用機序の個数とは一致しない。

細胞への対応付けは、標的遺伝子の発現だけで判定する。
比較図の色は、発現条件を満たした標的・薬剤の割合または件数を表し、発現量そのものや薬効の強さを表さない。
疾患間の差は治療薬と標的の構成の差であり、同じ細胞型における疾患別の発現変化ではない。
発現値がない標的は未判定として扱い、測定値が閾値以下の標的と区別する。

発現データは Open Targets に統合された Tabula Sapiens のうち、細胞型のみで集約した値を使う。
組織別や組織と細胞型の組合せの行は混ぜない。
単位は `CPM(pseudobulk sum[counts])`、値はドナー間中央値である。
判定方法は次から選べる。

| 方法 | 条件 |
| --- | --- |
| Fixed CPM | 最低 CPM 以上。初期値は 0.5 |
| Fixed CPM + Target-relative median | 最低 CPM 以上、かつその標的の全参照細胞型にわたる中央値以上 |
| Fixed CPM + CELLEX specificity | 最低 CPM 以上、かつ CELLEX スコアが設定値以上。初期値は 0.5 |

閾値は画面から変更でき、入力を空にすると初期値を使う。
標的内中央値の参照細胞型は、表示する疾患や細胞を絞っても変えない。
CELLEX の閾値はアプリの判定規則であり、取得元の high ラベルと等号まで一致するとは限らない。
DICE の TPM と単位を混ぜず、単一細胞の生カウントも扱わない。
健常時の発現なので、疾患組織での発現や薬効を担う細胞を確定するものではない。
[発現データの作成方法](https://platform-docs.opentargets.org/target/baseline-expression)を参照。

ヒートマップの各マスには、確認できた支持の数または割合を記号なしで表示する。
ホバーと詳細の表示値に付く `≥` は、その値が下限であることを表す。
実数では未判定の標的や薬剤、標的未判明の薬剤を考慮する。
割合では分母内の未判定だけを考慮し、標的未判明の薬剤があるだけでは `≥` を付けない。
陽性が確認できず未判定が残るセルは、実数を下限値の 0、割合を NA とする。
疾患データ自体が利用できないセルは、実数と割合のどちらも NA とする。
条件を満たす薬剤がない疾患は標的数 0、分母がない割合は NA とする。
比較ヒートマップの欠測・未判定は、Percent では 0%、Count では 0 と表示し、ホバーで実測の 0 と区別する。
Dot plot は、Count と Percent がともに 0 の組み合わせと、どちらかが欠測・未判定の組み合わせを空白にする。
正の値の円はクリックとホバーに対応し、ホバーには Count、Percent、分母、未判定数を併記する。
発現詳細は Heatmap と Dot plot を即時に切り替え、選択は疾患を変えても保持する。
Heatmap の欠測は 0 と表示し、内部の欠測値や計算対象は変更しない。
発現 Dot plot は標的別 z スコアを色、CELLEX specificity を円の面積で示す。
円の尺度は 0 から 1 に固定する。
円の最小直径は 3 px とし、それを上回る範囲で面積を CELLEX に比例させる。
CPM か CELLEX が欠測なら空白にする。
大分類の CELLEX は、値がある元細胞を等重みで算術平均し、ホバーに有効セル数と総数を示す。
細胞型間の排他的な割合ではなく、複数の細胞型で同じ標的が数えられる。

細胞はヒートマップですべての大分類を表示し、細分類の展開と折りたたみで比較範囲を変える。
ヒートマップの大分類名を押すと、その直下に元細胞の行を展開する。
もう一度押すと元細胞の行を折りたたむ。
展開は `Update` を待たずに反映する。
展開前後で色が変わらないよう、色の範囲には折りたたみ中の元細胞も含める。
展開状態を変えても大分類の集計対象は変わらない。
大分類では所属細胞型で陽性になった標的・薬剤の集合を数え、サブタイプ数だけ重複して数えない。
T cell には CD4・CD8 系を含める。
表示は免疫系を先頭にした固定順とし、細分類では Cell Ontology の親子関係を使って近い細胞型をまとめる。
フィルターを変えても並び順を保ち、CD4 系や CD8 系の比較を追いやすくする。
分類ごとの細胞型数の違いは結果に影響し得る。
詳細の連続発現図は、選択した疾患と薬剤条件の標的を全元細胞型について示す。
色は標的ごとに全元細胞型の `log2(1 + CPM)` から計算した z スコアとし、ホバーでスコアと元の CPM を示す。
白い丸印は、現在適用中の発現条件を満たす標的と元細胞型の組み合わせを示す。
標的の列は全元細胞型の z スコアの類似性で並べ、発現が欠測の標的は末尾に置く。
条件を満たさなかった標的も確認できる。
初期表示は大分類とし、細胞名をクリックすると所属する元細胞を展開する。
大分類には、欠測を除いた所属細胞型の CPM の算術平均を表示する。
各細胞型を等しい重みで扱い、ホバーに平均 CPM と観測できた細胞型数を示す。
全所属細胞が欠測なら平均も欠測とする。
大分類の色は平均 CPM を `log2(1 + CPM)` に変換し、元細胞だけで求めた平均と標準偏差で標準化する。
大分類の値は表示専用で、標的内中央値、標準化の基準、標的の並び順、発現判定、比較図の集計には使わない。
白丸は元細胞だけに表示し、展開しても色の範囲と標的の並び順を保つ。

ヒートマップをクリックすると、詳細の疾患選択欄も更新する。
元記録表には、選択した疾患と薬剤条件に合う記録を表示する。
比較図の標的と薬剤は、いずれかの元細胞型で条件を満たせば1回だけ数える。
条件を変えても、表示対象に残っている選択は維持する。
元記録は10行ずつ表示し、元薬剤と標的の組み合わせごとに1行とする。
同じ薬剤でも標的が複数あれば複数行に分かれ、発現条件を満たさない標的の記録も残す。
表の下の `Page` で続きを表示できる。

モダリティは Open Targets の薬剤型で区分する。
`Small molecule`、`Antibody`、`Protein`、`Cell`、`Gene`、`Enzyme`、`Oligonucleotide`、`Antibody drug conjugate`、`Vaccine component`、`Oligosaccharide`、`Unknown` を個別に選べる。

## データ更新と検証

公開 API から取り直す。
初回は、発現データの取得に 1.5 時間から 2 時間かかる。
同じ版の既存ファイルから再利用できるときは 30 分程度で済む。

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
保存の形と取得の手順は[保存データの形と取得の手順](docs/data.md)にある。

```bash
pixi run test
```

テストは成分と標的の重複除去、疾患内の臨床段階、3種類の発現基準、大分類、欠測、画面の応答を確認する。
型は basedpyright で検査する。
ツールチップの位置調整と Escape キーの処理は、`pixi run check` が Node.js で検証する。
