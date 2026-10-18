"""Parse and validate HPO disease-phenotype annotation files."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from rare_disease_ai.hpo_ontology import HPOOntology

EXPECTED_COLUMNS = (
    "database_id",
    "disease_name",
    "qualifier",
    "hpo_id",
    "reference",
    "evidence",
    "onset",
    "frequency",
    "sex",
    "modifier",
    "aspect",
    "biocuration",
)


@dataclass(frozen=True)
class HPOAnnotation:
    """One disease-phenotype annotation, preserving its source fields."""

    database_id: str
    disease_name: str
    qualifier: str
    hpo_id: str
    reference: str
    evidence: str
    onset: str
    frequency: str
    sex: str
    modifier: str
    aspect: str
    biocuration: str


@dataclass(frozen=True)
class HPOAnnotationValidationReport:
    """Summary of annotation-file integrity checks."""

    annotation_count: int
    unique_disease_count: int
    unknown_hpo_ids: tuple[str, ...]
    missing_disease_ids: int
    missing_hpo_ids: int

    @property
    def is_valid(self) -> bool:
        """Return whether every annotation has known required identifiers."""
        return (
            not self.unknown_hpo_ids and self.missing_disease_ids == 0 and self.missing_hpo_ids == 0
        )


def _annotation_from_row(row: dict[str, str | None]) -> HPOAnnotation:
    """Convert one tab-separated row into a typed annotation."""
    values = {column: (row.get(column) or "").strip() for column in EXPECTED_COLUMNS}

    return HPOAnnotation(**values)


def read_hpo_annotations(
    path: str | Path,
) -> tuple[HPOAnnotation, ...]:
    """Read HPOA annotations, skipping comment and metadata lines."""
    annotations: list[HPOAnnotation] = []

    with Path(path).open("r", encoding="utf-8-sig", newline="") as file:
        annotations.extend(_read_annotation_stream(file))

    return tuple(annotations)


def _read_annotation_stream(file: TextIO) -> list[HPOAnnotation]:
    """Read annotation rows from an open text stream."""
    data_lines = (line for line in file if line.strip() and not line.startswith("#"))
    reader = csv.DictReader(data_lines, delimiter="\t")

    if reader.fieldnames is None:
        raise ValueError("The HPO annotation file is empty.")

    reader.fieldnames = [name.strip() for name in reader.fieldnames]

    if tuple(reader.fieldnames) != EXPECTED_COLUMNS:
        raise ValueError(
            "Unexpected HPO annotation columns. "
            f"Expected {EXPECTED_COLUMNS}, got {tuple(reader.fieldnames)}."
        )

    annotations: list[HPOAnnotation] = []

    for line_number, row in enumerate(reader, start=2):
        if None in row:
            raise ValueError(
                f"Malformed annotation row near data row {line_number}: too many columns."
            )

        annotation = _annotation_from_row(row)

        if not annotation.database_id:
            raise ValueError(f"Missing database_id near data row {line_number}.")

        if not annotation.hpo_id:
            raise ValueError(f"Missing hpo_id near data row {line_number}.")

        annotations.append(annotation)

    return annotations


def validate_hpo_annotations(
    annotations: tuple[HPOAnnotation, ...],
    ontology: HPOOntology,
) -> HPOAnnotationValidationReport:
    """Check required identifiers and HPO IDs against the loaded ontology."""
    unknown_hpo_ids = sorted(
        {annotation.hpo_id for annotation in annotations if annotation.hpo_id not in ontology.terms}
    )

    return HPOAnnotationValidationReport(
        annotation_count=len(annotations),
        unique_disease_count=len({annotation.database_id for annotation in annotations}),
        unknown_hpo_ids=tuple(unknown_hpo_ids),
        missing_disease_ids=sum(not annotation.database_id for annotation in annotations),
        missing_hpo_ids=sum(not annotation.hpo_id for annotation in annotations),
    )
