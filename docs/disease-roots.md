# 対象疾患の起点

対象疾患は、Open Targets の `autoimmune disease`（`MONDO_0007179`）と、表 A から表 C で採用した追加起点から作る。
下位語を含める起点では、その下位語も対象にする。
MONDO は idiopathic inflammatory myopathy、lupus nephritis、ANCA 関連血管炎、primary biliary cholangitis などを別の枝に置いているため、基底の 1 語だけでは拾えない。
このノートは、起点の一覧と、その根拠の記録である。

基準は 3 層にする。

1. 基底の起点は `MONDO_0007179` で、下位語をすべて含める。
2. 追加の起点は、MeSH の `Autoimmune Diseases`（D001327）と Disease Ontology の `autoimmune disease`（DOID:417）が自己免疫に置く語、臨床の分類で自己免疫疾患とされる語、そして自己炎症性疾患のうち、MONDO の基底の下に無いものとする。起点 1 件ごとに、どの分類が根拠かを表に残す。
3. 起点にしないものを名指しする。アレルギー性疾患、非免疫性の疾患を多く含む包括語、動物モデル、Open Targets に語が無いものがこれにあたる。

自己炎症性疾患を含めるのは利用者の判断である。
治療薬の標的が自己免疫疾患と重なり、比較の材料になるためで、対象は「自己免疫疾患」より「免疫介在性の炎症疾患」に近くなる。
アレルギー性疾患は含めない。

候補は 4 つの源から機械的に集めた。
Open Targets 26.06 の基底の下位語が 170 語、MeSH の枝が 74 語、DOID の枝が 96 語、これに加えて臨床の分類で自己免疫疾患とされる 80 語を手で挙げた。
MeSH と DOID の語は、MONDO が公開している ID の対応表（`mondo.sssom.tsv`、Monarch Initiative の GitHub リポジトリの `src/ontology/mappings/`、2026-09-01 版）で MONDO の ID に対応づけた。
170 語のうち 156 語が対応表で決まり、決まらなかった 14 語と、対応表の MONDO の ID が Open Targets に無い 24 語は、Open Targets の名前検索で補った。
手書きの語は ID を持たないので、すべて名前検索で対応づけた。
Open Targets の疾患 1 語ごとの `dbXRefs` にも同じ対応が入っており、一括の `disease.parquet`（約 7 MB）から逆引きの表を作れる。
parquet を読むパッケージが環境に無いので、今回は MONDO の対応表を使った。
参照する予定の Hayter & Cook 2012（81 疾患）は本文を入手できていないので、入手できたらこの一覧と差分を取る。

判定は 2026-09-24 に利用者が行った。
表 A と表 B の語は全部入れ、自己炎症は MONDO の `autoinflammatory syndrome` を起点にし、アレルギー性の 3 語だけを外した。
基底と追加の起点を合わせた対象は 370 語になる。

## 起点の書き方

下位語の列は、その語を起点にしたときに一緒に入る語の数である。
下位語に非免疫性の疾患が混ざる語は、「本体のみ」として下位語を含めずに入れる。
別の起点の下位語になっている語は、起点にしない。
版が上がって親子関係が変わると入らなくなるので、表には残して「含まれる先」を書く。

## 表 A：MeSH か DOID が自己免疫に置く語

MONDO の基底の下に無く、MeSH か DOID が自己免疫の枝に置いている語である。

| 起点 ID | 名前 | 自己免疫に置く分類 | MONDO の親 | 下位語 | 扱い |
| --- | --- | --- | --- | --- | --- |
| MONDO_0019100 | neuromyelitis optica | MeSH、DOID、手書き | nervous system disorder | 2 | 入れる |
| MONDO_0006702 | chronic inflammatory demyelinating polyradiculoneuropathy | MeSH、DOID、手書き | demyelinating polyneuropathy | 0 | 入れる |
| EFO_0020094 | Lambert-Eaton myasthenic syndrome | MeSH、DOID、手書き | paraneoplastic neurologic syndrome | 0 | 入れる |
| MONDO_0008491 | stiff-person syndrome | MeSH、手書き | nervous system disorder | 0 | 入れる |
| MONDO_0021081 | anti-NMDA receptor encephalitis | MeSH | encephalitis | 0 | 入れる |
| MONDO_0019383 | acute disseminated encephalomyelitis | MeSH | postinfectious encephalitis | 1 | 入れる |
| MONDO_0011716 | acute hemorrhagic leukoencephalitis | MeSH | acute disseminated encephalomyelitis | 0 | acute disseminated encephalomyelitis に含まれる |
| MONDO_0015342 | acute transverse myelitis | MeSH | myelitis | 1 | 入れる |
| MONDO_0005340 | alopecia areata | DOID、手書き | alopecia | 0 | 入れる |
| MONDO_0007191 | Behcet disease | DOID、手書き | skin vascular disease | 0 | 入れる |
| MONDO_0015492 | anti-neutrophil cytoplasmic antibody-associated vasculitis | MeSH | necrotizing vasculitis | 3 | 入れる |
| MONDO_0008538 | temporal arteritis（giant cell arteritis） | MeSH、手書き | central nervous system vasculitis | 0 | 入れる |
| MONDO_0003346 | central nervous system vasculitis | MeSH | vasculitis | 4 | 本体のみ入れる |
| MONDO_0005556 | lupus nephritis | MeSH、手書き | glomerulonephritis | 0 | 入れる |
| EFO_1001363 | Lupus Vasculitis, Central Nervous System | MeSH | meningoencephalitis | 0 | 入れる |
| MONDO_0005342 | IgA glomerulonephritis | MeSH、手書き | glomerulonephritis | 0 | 入れる |
| MONDO_0005376 | membranous glomerulonephritis | MeSH、手書き | glomerulonephritis | 1 | 入れる |
| MONDO_0011429 | juvenile idiopathic arthritis | MeSH、手書き | arthritic joint disease | 10 | 入れる |
| MONDO_0019355 | adult-onset Still disease | MeSH、手書き | autoinflammatory syndrome | 0 | autoinflammatory syndrome に含まれる |
| MONDO_0018092 | Vogt-Koyanagi-Harada disease | MeSH、DOID、手書き | panuveitis | 0 | 入れる |
| MONDO_0011599 | birdshot chorioretinopathy | MeSH、手書き | posterior uveitis | 0 | 入れる |
| MONDO_0019198 | sympathetic ophthalmia | MeSH、手書き | panuveitis | 0 | 入れる |
| MONDO_0015129 | chronic primary adrenal insufficiency（Addison disease） | MeSH、手書き | primary adrenal insufficiency | 16 | 本体のみ入れる |
| MONDO_0014629 | autoimmune interstitial lung disease-arthritis syndrome | DOID | type 1 interferonopathy | 0 | autoinflammatory syndrome に含まれる |

対応表に無く、名前検索で補った行が 4 つある。
Lambert-Eaton myasthenic syndrome と Lupus Vasculitis, Central Nervous System は、対応表の MONDO の ID（`MONDO:0018556`、`MONDO:0043985`）が Open Targets に無く、同じ名前の EFO の語を使った。
temporal arteritis と birdshot chorioretinopathy は、MeSH の ID が対応表に無かった。

本体のみにした理由は 2 つある。
central nervous system vasculitis の下位語 4 語には central nervous system AIDS arteritis が入る。
Addison disease の下位語 16 語は、先天性副腎過形成などの遺伝性疾患だけである。

antisynthetase syndrome（`MONDO_0019344`）も DOID と手書きの両方にあるが、idiopathic inflammatory myopathy の下位語なので、表 B の起点を入れれば一緒に入る。
同じく granulomatosis with polyangiitis、microscopic polyangiitis、eosinophilic granulomatosis with polyangiitis は、ANCA 関連血管炎の起点の下位語 3 語そのものである。

## 表 B：臨床の分類で自己免疫疾患とされる語

MeSH にも DOID にも無く、利用者の判定で入れる語である。

| 起点 ID | 名前 | MONDO の親 | 下位語 | 扱い |
| --- | --- | --- | --- | --- |
| MONDO_0600023 | idiopathic inflammatory myopathy | myositis disease | 12 | 入れる。分類基準（ACR/EULAR 2017） |
| MONDO_0007827 | inclusion body myositis | myositis disease | 4 | 入れる。分類基準では IIM に含まれる |
| MONDO_0005388 | primary biliary cholangitis | biliary tract disorder | 0 | 入れる |
| MONDO_0018646 | sclerosing cholangitis（primary sclerosing cholangitis） | cholangitis | 1 | 本体のみ入れる |
| MONDO_0008228 | pernicious anemia | megaloblastic anemia | 0 | 入れる |
| MONDO_0019740 | acquired thrombotic thrombocytopenic purpura | thrombotic thrombocytopenic purpura | 0 | 入れる |
| MONDO_0019125 | relapsing polychondritis | chondromalacia | 0 | 入れる |
| MONDO_0005083 | psoriasis | immune system disorder、dermatitis | 7 | 入れる |
| MONDO_0005011 | Crohn disease | inflammatory bowel disease | 10 | 入れる |
| MONDO_0005101 | ulcerative colitis | colitis | 2 | 入れる |
| MONDO_0000702 | microscopic colitis | colitis | 2 | 入れる |
| MONDO_0017991 | Takayasu arteritis | arteritis | 0 | 入れる |
| MONDO_0019170 | polyarteritis nodosa | arteritis | 3 | 入れる |
| EFO_1000965 | Henoch-Schoenlein purpura（IgA vasculitis） | hypersensitivity reaction disease | 0 | 入れる |
| MONDO_0012727 | Kawasaki disease | vasculitis | 0 | 入れる |
| MONDO_0019735 | polymyalgia rheumatica | rheumatic disorder | 0 | 入れる |
| MONDO_0016158 | narcolepsy-cataplexy syndrome | narcolepsy | 0 | 入れる |
| MONDO_0006572 | lichen planus | dermatitis | 5 | 入れる |
| MONDO_0007899 | lichen sclerosus et atrophicus | dermatitis | 1 | 入れる |
| MONDO_0044212 | chronic idiopathic urticaria | idiopathic urticaria | 0 | 入れる |
| MONDO_0006835 | lipoid nephrosis（minimal change disease） | glomerulonephritis | 0 | 入れる |
| MONDO_0006559 | hidradenitis suppurativa | hidradenitis | 3 | 入れる |
| MONDO_0019338 | sarcoidosis | autoinflammatory syndrome | 8 | autoinflammatory syndrome に含まれる |
| MONDO_0009813 | chronic recurrent multifocal osteomyelitis | autoinflammatory syndrome | 3 | autoinflammatory syndrome に含まれる |
| MONDO_0005361 | eosinophilic esophagitis | eosinophilic gastrointestinal disease | 1 | 入れない。アレルギー性 |
| MONDO_0004980 | atopic eczema | allergic disease | 8 | 入れない。アレルギー性 |
| MONDO_0004979 | asthma | bronchial disorder | 10 | 入れない。アレルギー性 |

sclerosing cholangitis の下位語 1 語は neonatal ichthyosis-sclerosing cholangitis syndrome という遺伝性疾患なので、本体のみにする。
lichen planus の下位語には lichenoid drug reaction が、psoriasis の下位語には遺伝子座ごとの語が入る。
どちらも本体と同じ枠で扱い、除かない。

## 表 C：自己炎症性疾患の起点

| 起点 ID | 名前 | MONDO の親 | 下位語 | 扱い |
| --- | --- | --- | --- | --- |
| MONDO_0019751 | autoinflammatory syndrome | rheumatic disorder | 97 | 入れる |

下位語には adult-onset Still disease、sarcoidosis、chronic recurrent multifocal osteomyelitis、CINCA syndrome、PFAPA syndrome、STING-associated vasculopathy、familial Mediterranean fever、cryopyrin-associated periodic syndrome、VEXAS syndrome、I 型インターフェロノパチー、systemic-onset juvenile idiopathic arthritis が入る。
周期熱症候群のような単一遺伝子疾患も入る。
1 つの起点にしたのは、版が上がって新しい疾患がこの枝に入ったときに、自動で対象に入るようにするためである。

## 表 D：起点にしない候補

機械で集めた候補のうち、基準の 3 層目で除くものと、対応づけが取れないものである。
次に版が上がったときに同じ語がまた候補に出るので、理由を残す。

| 語 | 出どころ | 除く理由 |
| --- | --- | --- |
| eosinophilic esophagitis、atopic eczema、asthma | 手書き | アレルギー性疾患 |
| uveitis（31）、interstitial lung disease（51）、aplastic anemia（42） | 手書き | 包括語。非免疫性の疾患を多く含む。個別の疾患を起点にする |
| experimental autoimmune encephalomyelitis、experimental autoimmune myasthenia gravis、experimental autoimmune neuritis | MeSH | 動物モデル |
| Caplan syndrome、cauda equina syndrome、central nervous system AIDS arteritis、Wolfram syndrome | MeSH | MeSH の枝に入っているが自己免疫疾患ではない |
| Rheumatoid Nodule、Rheumatoid Vasculitis、autoimmune cholangitis、Hirata disease、autoimmune polyendocrine syndrome type 2、Balo concentric sclerosis | MeSH、DOID | 対応表の MONDO の ID が Open Targets に無い。起点にできない |
| Polyradiculoneuropathy、Polyradiculopathy | MeSH | 対応表で MONDO の一般語に決まるが、免疫性でないものを含む包括語 |
| Diffuse Cerebral Sclerosis of Schilder、Latent Autoimmune Diabetes in Adults、Undifferentiated Connective Tissue Diseases、MOG antibody-associated disease | MeSH | 対応表に MeSH の ID が無く、名前検索でも Open Targets に対応する語が見つからない |
| cryoglobulinemic vasculitis、idiopathic pulmonary hemosiderosis、autoimmune neutropenia、Cogan syndrome | 手書き | Open Targets に対応する疾患語が無い。起点にできない |
| premature ovarian insufficiency、pure red cell aplasia | 手書き | 対応づけ先が HP の表現型で、疾患語が無い |
| autoimmune disease of urogenital tract、infantile onset multisystem autoimmune disease 1 から 5、ankylosing spondylitis 1 から 3 | DOID | DOID だけの群の名前と遺伝子座ごとの語。対応する MONDO の語は基底の下にある |

## 起点はコードとスナップショットで管理する

起点と下位語を含めるかどうかは、`backend/disease_catalog.py` の `SCOPE_ROOTS` に定義している。
`backend/refresh.py` は起点ごとの対象語を集めて和集合を作り、起点の一覧とともにスナップショットへ保存する。
閲覧用の群と各群の起点は、`backend/disease_catalog.py` の `DISEASE_GROUPS` に定義している。
対象を保守するときは、表の採否と `SCOPE_ROOTS`、全疾患の配置先、独立して表示する疾患を確認し、`tests/test_disease_catalog.py` の回帰テストも更新する。
