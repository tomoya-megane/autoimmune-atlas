"""画面表示とコールバックで共有する固定値。"""

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
