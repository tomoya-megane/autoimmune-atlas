"""閲覧用の疾患群と、Open Targets の親子関係に基づく疾患ファミリー。

群は表示上の整理であり、医学的な分類体系や集計単位を定義しない。
指定した疾患に最も近い親をファミリーに選び、複数の親がある場合は定義順で決める。
乾癬性関節炎と強直性脊椎炎は独立した起点とし、関節リウマチの下に表示しない。
根拠と保守方法は docs/design.md の疾患選択に記載する。
"""

from collections import deque


# 各 ID の名称は snapshot から取得する。ここで決めるのは閲覧用の群と起点だけ。
DISEASE_GROUPS = (
    ("systemic", "Systemic and connective tissue", (
        "MONDO_0007915", "MONDO_0004670", "MONDO_0005100", "MONDO_0019562", "MONDO_0019340",
        "MONDO_0010030", "MONDO_0016663", "MONDO_0005854", "MONDO_0017287", "MONDO_8000010",
        "MONDO_0019557", "MONDO_0019191", "MONDO_0030703", "MONDO_0005435", "MONDO_0005576", "MONDO_0017767",
    )),
    ("joints", "Joint and spine", (
        "MONDO_0008383", "MONDO_0011849", "MONDO_0005306", "MONDO_0000589",
    )),
    ("nervous", "Nervous system", (
        "MONDO_0005301", "MONDO_0009688", "MONDO_0016218", "MONDO_0020640", "MONDO_0000774",
        "MONDO_0005851", "MONDO_0002977", "MONDO_0006704", "EFO_0803379", "MONDO_0019390",
    )),
    ("endocrine", "Endocrine system", (
        "MONDO_0005147", "MONDO_0005623", "MONDO_0005364", "MONDO_0007699",
        "MONDO_0017278", "MONDO_0019835", "MONDO_0000569",
    )),
    ("skin", "Skin", (
        "MONDO_0006594", "MONDO_0019082", "MONDO_0018746", "MONDO_0015614",
        "MONDO_0018974", "MONDO_0006558", "MONDO_0019337", "MONDO_0008661", "MONDO_0025513",
    )),
    ("digestive", "Digestive system and liver", (
        "MONDO_0005130", "MONDO_0016264", "MONDO_0019787", "MONDO_0015175", "MONDO_0031014", "MONDO_0000588",
    )),
    ("blood", "Blood", (
        "MONDO_0020108", "MONDO_0019098", "MONDO_0016030", "MONDO_0000602",
    )),
    ("cardiovascular", "Cardiovascular system", ("MONDO_0030701", "MONDO_0022519", "MONDO_0000603")),
    ("kidney", "Kidney", ("MONDO_0030700", "MONDO_0009303")),
    ("lung", "Lung", ("MONDO_0012579",)),
    ("eye-ear-exocrine", "Eye, ear and exocrine glands", (
        "MONDO_0000587", "MONDO_0000586", "MONDO_0100014", "MONDO_0031012",
    )),
    ("immune-dysregulation", "Immune dysregulation syndromes", (
        "MONDO_0000213", "MONDO_0017979", "MONDO_0010580", "EFO_0010647",
    )),
)

# 元の自己免疫疾患の範囲には残すが、自己免疫性の病型としては配置しない。
UNCLASSIFIED_IDS = {"MONDO_0008218", "MONDO_0011431"}  # Hailey-Hailey disease, MASS syndrome


def disease_catalog(snapshot):
    """全疾患を一度ずつ収める。親情報や分類先がない用語も失わない。"""
    diseases = {d["id"]: d for d in snapshot["diseases"]}
    roots = [root for _, _, ids in DISEASE_GROUPS for root in ids if root in diseases]
    rank = {root: i for i, root in enumerate(roots)}
    families = {root: [] for root in roots}
    other = []
    for disease in diseases.values():
        if disease["id"] in UNCLASSIFIED_IDS:
            other.append(disease)
            continue
        distances, queue = {}, deque([(disease["id"], 0)])
        while queue:
            node, distance = queue.popleft()
            if node in distances:
                continue
            distances[node] = distance
            queue.extend((parent, distance + 1) for parent in diseases.get(node, {}).get("parent_ids", []))
        candidates = [root for root in roots if root in distances]
        if candidates:
            root = min(candidates, key=lambda item: (distances[item], rank[item]))
            families[root].append(disease)
        else:
            other.append(disease)
    groups = []
    for group_id, label, ids in DISEASE_GROUPS:
        entries = [{"id": root, "label": diseases[root]["name"], "diseases": sorted(families[root], key=lambda d: (d["id"] != root, d["name"].casefold(), d["id"]))} for root in ids if families.get(root)]
        if entries:
            groups.append({"id": group_id, "label": label, "families": sorted(entries, key=lambda f: (f["label"].casefold(), f["id"]))})
    groups.sort(key=lambda group: (group["label"].casefold(), group["id"]))
    if other:
        groups.append({"id": "other", "label": "Other / unclassified", "families": [{"id": "other", "label": "Other terms", "diseases": sorted(other, key=lambda d: (d["name"].casefold(), d["id"]))}]})
    return groups


def ordered_disease_ids(snapshot, selected):
    """チェック順に依存せず、閲覧用の群とファミリーの順に並べる。"""
    selected = set(selected or [])
    return [d["id"] for group in disease_catalog(snapshot) for family in group["families"] for d in family["diseases"] if d["id"] in selected]
