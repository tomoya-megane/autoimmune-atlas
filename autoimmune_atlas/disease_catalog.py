"""対象疾患の起点と、閲覧用の疾患群と、Open Targets の親子関係に基づく疾患ファミリー。

対象疾患は SCOPE_ROOTS の各起点とその下位語の和集合である。起点の根拠は docs/disease-roots.md にある。
群は表示上の整理であり、医学的な分類体系や集計単位を定義しない。
指定した疾患に最も近い親をファミリーに選び、複数の親がある場合は定義順で決める。
乾癬性関節炎と強直性脊椎炎は独立した起点とし、関節リウマチの下に表示しない。
根拠と保守方法は docs/design.md の疾患選択に記載する。
"""

from collections import deque

from autoimmune_atlas.models import (
    CatalogDisease,
    DiseaseCatalogGroup,
    DiseaseCatalogInput,
    DiseaseFamily,
    Snapshot,
)

# 対象疾患の起点。(ID, 下位語を含めるか)。先頭が基底の autoimmune disease。
# 下位語を含めない語は、下位語に遺伝性や感染性の疾患が混ざるので本体だけを入れる。
SCOPE_ROOTS = (
    ("MONDO_0007179", True),  # autoimmune disease
    ("MONDO_0019751", True),  # autoinflammatory syndrome
    ("MONDO_0019100", True),  # neuromyelitis optica
    ("MONDO_0006702", True),  # CIDP
    ("EFO_0020094", True),  # Lambert-Eaton myasthenic syndrome
    ("MONDO_0008491", True),  # stiff-person syndrome
    ("MONDO_0021081", True),  # anti-NMDA receptor encephalitis
    ("MONDO_0019383", True),  # acute disseminated encephalomyelitis
    ("MONDO_0015342", True),  # acute transverse myelitis
    ("MONDO_0016158", True),  # narcolepsy-cataplexy syndrome
    ("MONDO_0600023", True),  # idiopathic inflammatory myopathy
    ("MONDO_0007827", True),  # inclusion body myositis
    ("MONDO_0011429", True),  # juvenile idiopathic arthritis
    ("MONDO_0019735", True),  # polymyalgia rheumatica
    ("MONDO_0019125", True),  # relapsing polychondritis
    ("MONDO_0015492", True),  # ANCA-associated vasculitis
    ("MONDO_0008538", True),  # temporal arteritis (giant cell arteritis)
    ("MONDO_0003346", False),  # central nervous system vasculitis
    ("EFO_1001363", True),  # Lupus Vasculitis, Central Nervous System
    ("MONDO_0007191", True),  # Behcet disease
    ("MONDO_0017991", True),  # Takayasu arteritis
    ("MONDO_0019170", True),  # polyarteritis nodosa
    ("EFO_1000965", True),  # Henoch-Schoenlein purpura (IgA vasculitis)
    ("MONDO_0012727", True),  # Kawasaki disease
    ("MONDO_0005556", True),  # lupus nephritis
    ("MONDO_0005342", True),  # IgA glomerulonephritis
    ("MONDO_0005376", True),  # membranous glomerulonephritis
    ("MONDO_0006835", True),  # lipoid nephrosis (minimal change disease)
    ("MONDO_0018092", True),  # Vogt-Koyanagi-Harada disease
    ("MONDO_0011599", True),  # birdshot chorioretinopathy
    ("MONDO_0019198", True),  # sympathetic ophthalmia
    ("MONDO_0015129", False),  # chronic primary adrenal insufficiency (Addison disease)
    ("MONDO_0005388", True),  # primary biliary cholangitis
    ("MONDO_0018646", False),  # sclerosing cholangitis
    ("MONDO_0005011", True),  # Crohn disease
    ("MONDO_0005101", True),  # ulcerative colitis
    ("MONDO_0000702", True),  # microscopic colitis
    ("MONDO_0008228", True),  # pernicious anemia
    ("MONDO_0019740", True),  # acquired thrombotic thrombocytopenic purpura
    ("MONDO_0005083", True),  # psoriasis
    ("MONDO_0005340", True),  # alopecia areata
    ("MONDO_0006572", True),  # lichen planus
    ("MONDO_0007899", True),  # lichen sclerosus et atrophicus
    ("MONDO_0044212", True),  # chronic idiopathic urticaria
    ("MONDO_0006559", True),  # hidradenitis suppurativa
)

# 各 ID の名称は snapshot から取得する。ここで決めるのは閲覧用の群と起点だけ。
DISEASE_GROUPS = (
    (
        "systemic",
        "Systemic and connective tissue",
        (
            "MONDO_0007915",
            "MONDO_0004670",
            "MONDO_0005100",
            "MONDO_0019562",
            "MONDO_0019340",
            "MONDO_0010030",
            "MONDO_0016663",
            "MONDO_0005854",
            "MONDO_0017287",
            "MONDO_8000010",
            "MONDO_0019557",
            "MONDO_0019191",
            "MONDO_0030703",
            "MONDO_0005435",
            "MONDO_0005576",
            "MONDO_0017767",
        ),
    ),
    (
        "joints",
        "Joint and spine",
        (
            "MONDO_0008383",
            "MONDO_0011849",
            "MONDO_0005306",
            "MONDO_0000589",
            "MONDO_0011429",
            "MONDO_0019735",
            "MONDO_0019125",
        ),
    ),
    (
        "muscle",
        "Muscle",
        ("MONDO_0600023", "MONDO_0007827"),
    ),
    (
        "nervous",
        "Nervous system",
        (
            "EFO_0803536",
            "MONDO_0009688",
            "MONDO_0016218",
            "MONDO_0020640",
            "MONDO_0000774",
            "MONDO_0005851",
            "MONDO_0002977",
            "MONDO_0006704",
            "EFO_0803379",
            "MONDO_0019390",
            "MONDO_0019100",
            "MONDO_0006702",
            "EFO_0020094",
            "MONDO_0008491",
            "MONDO_0021081",
            "MONDO_0019383",
            "MONDO_0015342",
            "MONDO_0016158",
        ),
    ),
    (
        "endocrine",
        "Endocrine system",
        (
            "MONDO_0005147",
            "MONDO_0005623",
            "MONDO_0005364",
            "MONDO_0007699",
            "MONDO_0017278",
            "MONDO_0019835",
            "MONDO_0000569",
            "MONDO_0015129",
        ),
    ),
    (
        "skin",
        "Skin",
        (
            "MONDO_0006594",
            "MONDO_0019082",
            "MONDO_0018746",
            "MONDO_0015614",
            "MONDO_0018974",
            "MONDO_0006558",
            "MONDO_0019337",
            "MONDO_0008661",
            "MONDO_0025513",
            "MONDO_0005083",
            "MONDO_0005340",
            "MONDO_0006572",
            "MONDO_0007899",
            "MONDO_0044212",
            "MONDO_0006559",
        ),
    ),
    (
        "digestive",
        "Digestive system and liver",
        (
            "MONDO_0005130",
            "MONDO_0016264",
            "MONDO_0019787",
            "MONDO_0015175",
            "MONDO_0031014",
            "MONDO_0000588",
            "MONDO_0005011",
            "MONDO_0005101",
            "MONDO_0000702",
            "MONDO_0005388",
            "MONDO_0018646",
        ),
    ),
    (
        "blood",
        "Blood",
        (
            "MONDO_0020108",
            "MONDO_0019098",
            "MONDO_0016030",
            "MONDO_0000602",
            "MONDO_0008228",
            "MONDO_0019740",
        ),
    ),
    (
        "cardiovascular",
        "Cardiovascular system",
        ("MONDO_0030701", "MONDO_0022519", "MONDO_0000603"),
    ),
    (
        "vasculitis",
        "Vasculitis",
        (
            "MONDO_0015492",
            "MONDO_0008538",
            "MONDO_0003346",
            "EFO_1001363",
            "MONDO_0007191",
            "MONDO_0017991",
            "MONDO_0019170",
            "EFO_1000965",
            "MONDO_0012727",
        ),
    ),
    (
        "kidney",
        "Kidney",
        (
            "MONDO_0030700",
            "MONDO_0009303",
            "MONDO_0005556",
            "MONDO_0005342",
            "MONDO_0005376",
            "MONDO_0006835",
        ),
    ),
    ("lung", "Lung", ("MONDO_0012579",)),
    (
        "eye-ear-exocrine",
        "Eye, ear and exocrine glands",
        (
            "MONDO_0000587",
            "MONDO_0000586",
            "MONDO_0100014",
            "MONDO_0031012",
            "MONDO_0018092",
            "MONDO_0011599",
            "MONDO_0019198",
        ),
    ),
    ("autoinflammatory", "Autoinflammatory syndromes", ("MONDO_0019751",)),
    (
        "immune-dysregulation",
        "Immune dysregulation syndromes",
        (
            "MONDO_0000213",
            "MONDO_0017979",
            "MONDO_0010580",
            "EFO_0010647",
        ),
    ),
)

# 元の自己免疫疾患の範囲には残すが、自己免疫性の病型としては配置しない。
UNCLASSIFIED_IDS = {
    "MONDO_0008218",
    "MONDO_0011431",
}  # Hailey-Hailey disease, MASS syndrome


def disease_catalog(
    snapshot: Snapshot | DiseaseCatalogInput,
) -> list[DiseaseCatalogGroup]:
    """全疾患を一度ずつ収める。親情報や分類先がない用語も失わない。"""
    diseases: dict[str, CatalogDisease] = {
        disease["id"]: disease for disease in snapshot["diseases"]
    }
    roots = [root for _, _, ids in DISEASE_GROUPS for root in ids if root in diseases]
    rank = {root: i for i, root in enumerate(roots)}
    families: dict[str, list[CatalogDisease]] = {root: [] for root in roots}
    other: list[CatalogDisease] = []
    for disease in diseases.values():
        if disease["id"] in UNCLASSIFIED_IDS:
            other.append(disease)
            continue
        distances: dict[str, int] = {}
        queue: deque[tuple[str, int]] = deque([(disease["id"], 0)])
        while queue:
            node, distance = queue.popleft()
            if node in distances:
                continue
            distances[node] = distance
            parent_ids = (
                diseases[node].get("parent_ids", []) if node in diseases else []
            )
            queue.extend((parent, distance + 1) for parent in parent_ids)
        candidates = [root for root in roots if root in distances]
        if candidates:
            root = min(candidates, key=lambda item: (distances[item], rank[item]))
            families[root].append(disease)
        else:
            other.append(disease)
    groups: list[DiseaseCatalogGroup] = []
    for group_id, label, ids in DISEASE_GROUPS:
        entries: list[DiseaseFamily] = [
            {
                "id": root,
                "label": diseases[root]["name"],
                "diseases": sorted(
                    families[root],
                    key=lambda d: (d["id"] != root, d["name"].casefold(), d["id"]),
                ),
            }
            for root in ids
            if families.get(root)
        ]
        if entries:
            groups.append(
                {
                    "id": group_id,
                    "label": label,
                    "families": sorted(
                        entries, key=lambda f: (f["label"].casefold(), f["id"])
                    ),
                }
            )
    groups.sort(key=lambda group: (group["label"].casefold(), group["id"]))
    if other:
        groups.append(
            {
                "id": "other",
                "label": "Other / unclassified",
                "families": [
                    {
                        "id": "other",
                        "label": "Other terms",
                        "diseases": sorted(
                            other, key=lambda d: (d["name"].casefold(), d["id"])
                        ),
                    }
                ],
            }
        )
    return groups


def ordered_disease_ids(
    snapshot: Snapshot | DiseaseCatalogInput, selected: list[str] | None
) -> list[str]:
    """チェック順に依存せず、閲覧用の群とファミリーの順に並べる。"""
    selected_ids = set(selected or [])
    return [
        d["id"]
        for group in disease_catalog(snapshot)
        for family in group["families"]
        for d in family["diseases"]
        if d["id"] in selected_ids
    ]
