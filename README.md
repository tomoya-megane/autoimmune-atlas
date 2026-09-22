# Autoimmune Atlas

自己免疫疾患ごとに、Phase III 以降の薬剤の分子標的を細胞型に対応付けるローカルアプリ。
Python の Dash と Plotly を使う。

画面構成は Open Targets の [Associations on the Fly](https://platform-docs.opentargets.org/web-interface/associations-on-the-fly) を参考にしている。
小さな概要欄、表示モードのタブ、比較表の近くのフィルター、根拠を表示する詳細欄を採用した。
配色と余白はこのアプリ向けの調整であり、Open Targets のテーマをそのまま移植したものではない。

## 起動

このディレクトリで実行する。

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

対象は Open Targets の autoimmune disease (`MONDO_0007179`) の下位語である。
疾患の階層に従うため、親概念や疾患サブタイプが混在する。
各疾患の `drugAndClinicalCandidates.maxClinicalStage` が `PHASE_3`、`PREAPPROVAL`、`PHASE_4`、`APPROVAL` の薬剤を含める。
薬剤全体の最高段階は使わない。
`PHASE_2_3`、段階不明、`WITHDRAWAL` は含めない。
これは開発段階の記録による選択であり、現在の開発継続や臨床的成功を保証しない。
[臨床段階の定義](https://platform-docs.opentargets.org/drug/clinical-report)を参照。

分子標的数は Ensembl 遺伝子 ID の重複を除いた数である。
同じ標的を持つ複数薬剤は一度だけ数える。
割合の分母は、その疾患でモダリティ条件を満たす薬剤の**既知の標的数**であり、細胞型ごとには変えない。
標的未判明の薬剤は分母に入れず、件数を別に表示する。
複合体の構成遺伝子はそれぞれ数えるため、標的複合体や作用機序の個数とは一致しない。

| 細胞への対応付け | 判定 | 解釈の限界 |
| --- | --- | --- |
| 標的発現 | Tabula Sapiens の細胞型別 pseudobulk のドナー中央値が閾値を超える | 健常時の発現であり、薬効や疾患時の発現を示さない |
| 薬効に関わる細胞 | 疾患、薬剤、標的、細胞を指定した出典付き注釈 | 注釈は一部のみ。未調査を陰性にしない |
| 未判定 | 発現値や薬効注釈がない | 0 と区別して NA または下限値として表示する |

発現データは Open Targets に統合された Tabula Sapiens のうち、細胞型のみで集約した値を使う。
組織別や組織と細胞型の組合せの行は混ぜない。
単位は `CPM(pseudobulk sum[counts])`、初期閾値は中央値 `> 0.5 CPM`。
閾値は画面から変更できる。
DICE の TPM と単位を混ぜず、単一細胞の生カウントも扱わない。
[発現データの作成方法](https://platform-docs.opentargets.org/target/baseline-expression)を参照。

`≥` は少なくとも確認できた標的数または割合を表す。
標的数では、標的の未判定や標的未判明の薬剤があるときに表示する。
割合では、分母に含む既知標的に未判定が残る場合だけ表示する。
標的未判明の薬剤があるだけでは、既知標的内の割合に `≥` を付けない。
割合はあくまで既知標的内の値で、全薬剤の標的を網羅した割合ではない。
陽性が確認できず未判定が残るセルは NA、全対象を判定できた陰性は 0 とする。
条件を満たす薬剤がない疾患は標的数 0、分母がない割合は NA とする。
細胞型間の排他的な割合ではなく、複数の細胞型で同じ標的が数えられる。

モダリティは Open Targets の薬剤型で区分する。
`Small molecule` は低分子、`Antibody` は抗体、既知のその他の型はその他とする。
型不明は「すべて」にのみ含める。

## 薬効注釈の追加

`data/cell_annotations.csv` に、疾患 ID、薬剤 ID、標的 ID、細胞 ID、細胞名、`status`、出典 URL、根拠の要約を記す。
`status` は `yes` または `no` とし、未調査は行を作らない。
ある薬剤や疾患の注釈を、他の薬剤や疾患へ自動で広げない。
編集後はアプリを再起動する。

初期注釈は EMA の [MabThera](https://www.ema.europa.eu/en/medicines/human/EPAR/mabthera)、[Orencia](https://www.ema.europa.eu/en/medicines/human/EPAR/orencia)、[Benlysta](https://www.ema.europa.eu/en/medicines/human/EPAR/benlysta) に基づく。
RA の rituximab と abatacept、SLE の belimumab に限って登録してあり、薬効モードは網羅的な比較には使えない。

## データ更新と検証

公開 API から取り直す。
通信量と API の応答により数分から十数分かかる。

```bash
pixi run refresh
```

全取得と件数の検証が成功したときだけ `data/snapshot.json` を置き換える。
途中で失敗した場合は前回のデータを保つ。
取得日時と Open Targets のデータ版はスナップショットに保存する。
更新後はアプリを再起動する。

```bash
pixi run test
```

テストは標的の重複除去、モダリティによる分母変更、欠測、薬効注釈の適用範囲、CSV、画面の応答を確認する。
