"""画面表示とコールバックで共有する固定値。"""

# 下限は genetics.py が持ち、既存の import のためにここからも読めるようにする。
from autoimmune_atlas.genetics import SCORE_FLOOR as SCORE_FLOOR

DEFAULT_EXPRESSION_THRESHOLD = 0.5
DEFAULT_SPECIFICITY_THRESHOLD = 0.5
SOURCE_PAGE_SIZE = 10
PARAMETER_IDS = (
    "measure",
    "modality",
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
DEFAULT_SCORE_THRESHOLD = 0.5
GENETICS_PARAMETER_IDS = (
    "genetics-measure",
    "genetics-score",
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
