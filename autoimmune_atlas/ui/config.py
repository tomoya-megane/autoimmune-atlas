"""画面表示とコールバックで共有する固定値。"""

DEFAULT_EXPRESSION_THRESHOLD = 0.5
DEFAULT_SPECIFICITY_THRESHOLD = 0.5
SOURCE_PAGE_SIZE = 10
PARAMETER_IDS = (
    "measure",
    "modality",
    "target-class",
    "target-location",
    "stage",
    "method",
    "threshold",
    "specificity",
    "diseases",
    "heatmap-view",
)
STAGE_LABELS = {
    "phase1": "Phase I or later",
    "phase2": "Phase II or later",
    "phase3": "Phase III or later",
    "approved": "Approval reached",
}
METHOD_LABELS = {
    "fixed": "Fixed CPM",
    "relative": "Fixed CPM + Target-relative median",
    "specificity": "Fixed CPM + CELLEX specificity",
}
# 分類の無い標的の表示名。ChEMBL の実ラベル Unclassified protein と取り違えないよう別の語にする。
NOT_ANNOTATED_LABEL = "Not annotated"
DEFAULT_SCORE_THRESHOLD = 0.5
GENETICS_PARAMETER_IDS = (
    "genetics-measure",
    "genetics-score",
    "genetics-target-class",
    "genetics-target-location",
    "genetics-method",
    "genetics-threshold",
    "genetics-specificity",
    "genetics-diseases",
)
# 薬剤ページと遺伝子ページが初期選択に使う 10 疾患。
DEFAULT_DISEASE_TERMS = (
    "systemic lupus erythematosus",
    "systemic sclerosis",
    "Sjogren syndrome",
    "rheumatoid arthritis",
    "myasthenia gravis",
    "dermatomyositis",
    "type 1 diabetes mellitus",
    "anti-neutrophil cytoplasmic antibody-associated vasculitis",
    "pemphigus",
    "autoimmune hepatitis",
)
