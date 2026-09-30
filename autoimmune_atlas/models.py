"""スナップショットと集計結果の共有データ型。"""

from typing import Literal, NotRequired, TypedDict


class Reference(TypedDict):
    source: str
    ids: list[str]
    urls: list[str]


class DrugRecord(TypedDict):
    disease_id: str
    disease: str
    drug_id: str
    drug: str
    canonical_drug_id: str
    canonical_drug: str
    modality: NotRequired[str]
    drug_type: str
    stage: str
    target_id: str | None
    target: str | None
    mechanism: str | None
    action_types: NotRequired[list[str]]
    references: list[Reference]


class FilteredRecord(DrugRecord):
    canonical_stage: str
    target_class: NotRequired[str]
    location_classes: NotRequired[list[str]]


class EvidenceRecord(FilteredRecord):
    cell_id: str
    cell: str
    cpm: float
    specificity_score: float | None
    target_median: float | None
    evidence: str
    note: str


class TargetRecord(TypedDict):
    target_id: NotRequired[str | None]
    target: NotRequired[str | None]


class EvidenceBase(FilteredRecord):
    group_cell_id: str
    group_cell: str
    references_json: str


class EvidenceTableRow(EvidenceBase):
    evidence_cell_id: str
    evidence_cell: str
    cpm: float | None
    specificity_score: float | None
    target_median: float | None
    support_state: Literal["positive", "negative", "unknown", "unmapped"]
    contributing: bool


class Disease(TypedDict):
    id: str
    name: str
    parent_ids: NotRequired[list[str]]
    status: str
    unclassified_stages: NotRequired[int]


class CatalogDisease(TypedDict):
    id: str
    name: str
    parent_ids: NotRequired[list[str]]


class ExpressionRow(TypedDict):
    """保存する発現の行。細胞の定数は Snapshot の cells にある。"""

    cell_id: str
    median: float | None
    specificity_score: float | None


class CellDefinition(TypedDict):
    """細胞 1 つの定数。遺伝子が違っても同じなので 1 回だけ保存する。"""

    name: str
    parent_id: str | None
    parent: str | None
    ancestor_ids: list[str]


class ExpressionMetadata(ExpressionRow):
    """実行時の行。cells から引いた cell 名と、標的内の中央値を持つ。"""

    cell: str
    target_median: float | None


class ExpressionStateInput(TypedDict):
    median: float | None
    target_median: NotRequired[float | None]
    specificity_score: NotRequired[float | None]


class RootIdentity(TypedDict):
    id: str
    name: str


class SnapshotRoot(RootIdentity):
    include_descendants: bool
    count: int


class DataVersion(TypedDict):
    year: str
    month: str
    iteration: str | int | None


class GeneAssociation(TypedDict):
    target_id: str
    target: str
    score: float
    datasource_scores: dict[str, float]


class TargetLocation(TypedDict):
    location: str
    source: str


class TargetAnnotation(TypedDict):
    """標的 1 つの分類。Open Targets の targetClass と subcellularLocations から作る。"""

    target_class: str | None
    locations: list[TargetLocation]


class Snapshot(TypedDict):
    schema: int
    root: str | RootIdentity
    roots: NotRequired[list[SnapshotRoot]]
    data_version: NotRequired[DataVersion]
    retrieved_at: str
    source: str
    diseases: list[Disease]
    records: list[DrugRecord]
    associations: dict[str, list[GeneAssociation]]
    datasources: list[str]
    cells: dict[str, CellDefinition]
    expression: dict[str, list[ExpressionRow]]
    targets: dict[str, TargetAnnotation]


class AggregationSnapshot(TypedDict):
    schema: int
    diseases: NotRequired[list[Disease]]
    records: NotRequired[list[DrugRecord]]
    associations: NotRequired[dict[str, list[GeneAssociation]]]
    datasources: NotRequired[list[str]]
    cells: NotRequired[dict[str, CellDefinition]]
    expression: NotRequired[dict[str, list[ExpressionRow]]]
    targets: NotRequired[dict[str, TargetAnnotation]]


class CoreSnapshot(TypedDict):
    schema: int
    diseases: list[Disease]
    records: list[DrugRecord]
    associations: dict[str, list[GeneAssociation]]
    datasources: list[str]
    cells: dict[str, CellDefinition]
    expression: dict[str, list[ExpressionRow]]
    targets: dict[str, TargetAnnotation]


class DiseaseCatalogInput(TypedDict):
    diseases: list[CatalogDisease]


class CellMembership(TypedDict):
    id: str
    name: str
    group_id: str
    group_name: str
    ancestor_ids: tuple[str, ...]


class CellCatalogEntry(TypedDict):
    id: str
    name: str
    members: list[str]
    ontology_id: NotRequired[str]
    cell_level: NotRequired[str]


class DiseaseFamily(TypedDict):
    id: str
    label: str
    diseases: list[CatalogDisease]


class DiseaseCatalogGroup(TypedDict):
    id: str
    label: str
    families: list[DiseaseFamily]


class SummaryRow(TypedDict):
    disease_id: str
    disease: str
    cell_id: str
    cell: str
    ontology_id: str
    cell_level: str
    count: int | None
    percent: float | None
    denominator: int
    unknown: int
    unmapped_drugs: int
    status: str
    records: list[EvidenceRecord]
    drug_count: int | None
    drug_percent: float | None
    drug_denominator: int
    unknown_drugs: int
    total_drugs: int
    mapped_drugs: int
    member_cell_ids: list[str]
