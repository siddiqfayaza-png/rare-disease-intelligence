"""Import Gene Ontology Annotation Files (GAF 2.2)."""

import gzip
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path

from rare_disease_ai.mechanism_evidence import MechanismEvidence

ASPECT_TO_ANNOTATION_TYPE = {
    "F": "molecular_function",
    "P": "biological_process",
    "C": "cellular_component",
}


@dataclass(frozen=True)
class GOAnnotationImportResult:
    """Summary of imported GO annotations and data-quality issues."""

    annotations: tuple[MechanismEvidence, ...]
    rows_read: int
    row_errors: tuple[str, ...]
    rows_filtered_by_taxon: int


def iter_gaf_annotations(
    path: str | Path,
    *,
    source_release: str,
    source_url: str,
    term_labels: Mapping[str, str],
    taxon_id: str = "taxon:9606",
) -> Iterator[MechanismEvidence]:
    """Yield valid GAF annotations one at a time.

    This iterator skips malformed rows, unresolved GO labels, and rows whose
    primary taxon differs from the requested taxon. Unlike read_gaf_annotations,
    it does not accumulate annotations or row errors in memory.

    NOT qualifiers are preserved, not interpreted as positive evidence.
    For multi-taxon fields, only the first taxon is used for filtering.
    """
    if not isinstance(source_release, str) or not source_release.strip():
        raise ValueError("source_release must be non-empty")

    if not isinstance(source_url, str) or not source_url.strip():
        raise ValueError("source_url must be non-empty")

    if not isinstance(taxon_id, str) or not taxon_id.startswith("taxon:"):
        raise ValueError("taxon_id must use the taxon:NNNN format")

    file_path = Path(path)
    opener = gzip.open if file_path.suffix == ".gz" else open

    with opener(file_path, "rt", encoding="utf-8") as handle:
        for _, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("!"):
                continue

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 17:
                continue

            database, object_id, symbol = fields[0:3]
            qualifier = fields[3]
            go_id = fields[4]
            references = fields[5]
            evidence_code = fields[6]
            aspect = fields[8]
            taxon_field = fields[12]
            assigned_by = fields[14]

            required = {
                "database": database,
                "object_id": object_id,
                "gene_symbol": symbol,
                "GO ID": go_id,
                "reference": references,
                "evidence code": evidence_code,
                "aspect": aspect,
                "taxon": taxon_field,
                "assigned_by": assigned_by,
            }
            if any(not value.strip() for value in required.values()):
                continue

            primary_taxon = taxon_field.split("|")[0].strip()
            if primary_taxon != taxon_id:
                continue

            annotation_type = ASPECT_TO_ANNOTATION_TYPE.get(aspect)
            if annotation_type is None:
                continue

            term_label = term_labels.get(go_id)
            if not term_label or not term_label.strip():
                continue

            yield MechanismEvidence(
                gene_symbol=symbol,
                annotation_id=go_id,
                annotation_label=term_label,
                annotation_type=annotation_type,
                source_name=f"Gene Ontology Annotation ({database})",
                source_release=source_release,
                source_url=source_url,
                evidence_code=evidence_code,
                reference_ids=tuple(references.split("|")),
                qualifier=qualifier or None,
                taxon_id=primary_taxon,
                source_record_fields=tuple(fields),
            )


def read_gaf_annotations(
    path: str | Path,
    *,
    source_release: str,
    source_url: str,
    term_labels: Mapping[str, str],
    taxon_id: str = "taxon:9606",
) -> GOAnnotationImportResult:
    """Read GAF annotations and retain the existing diagnostic result API.

    This compatibility wrapper collects annotations and reports malformed
    rows, unsupported aspects, missing labels, and taxon-filter counts.
    Use iter_gaf_annotations for memory-efficient iteration.
    """
    if not isinstance(source_release, str) or not source_release.strip():
        raise ValueError("source_release must be non-empty")

    if not isinstance(source_url, str) or not source_url.strip():
        raise ValueError("source_url must be non-empty")

    if not isinstance(taxon_id, str) or not taxon_id.startswith("taxon:"):
        raise ValueError("taxon_id must use the taxon:NNNN format")

    file_path = Path(path)
    opener = gzip.open if file_path.suffix == ".gz" else open

    annotations: list[MechanismEvidence] = []
    errors: list[str] = []
    rows_read = 0
    rows_filtered_by_taxon = 0

    with opener(file_path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("!"):
                continue

            rows_read += 1
            fields = line.rstrip("\r\n").split("\t")

            if len(fields) != 17:
                errors.append(f"Line {line_number}: expected 17 GAF columns, found {len(fields)}")
                continue

            database, object_id, symbol = fields[0:3]
            qualifier = fields[3]
            go_id = fields[4]
            references = fields[5]
            evidence_code = fields[6]
            aspect = fields[8]
            taxon_field = fields[12]
            assigned_by = fields[14]

            required = {
                "database": database,
                "object_id": object_id,
                "gene_symbol": symbol,
                "GO ID": go_id,
                "reference": references,
                "evidence code": evidence_code,
                "aspect": aspect,
                "taxon": taxon_field,
                "assigned_by": assigned_by,
            }
            missing = [name for name, value in required.items() if not value.strip()]
            if missing:
                errors.append(f"Line {line_number}: missing required fields: {', '.join(missing)}")
                continue

            primary_taxon = taxon_field.split("|")[0].strip()
            if primary_taxon != taxon_id:
                rows_filtered_by_taxon += 1
                continue

            annotation_type = ASPECT_TO_ANNOTATION_TYPE.get(aspect)
            if annotation_type is None:
                errors.append(f"Line {line_number}: unsupported GO aspect {aspect!r}")
                continue

            term_label = term_labels.get(go_id)
            if not term_label or not term_label.strip():
                errors.append(f"Line {line_number}: no ontology label for {go_id}")
                continue

            annotations.append(
                MechanismEvidence(
                    gene_symbol=symbol,
                    annotation_id=go_id,
                    annotation_label=term_label,
                    annotation_type=annotation_type,
                    source_name=f"Gene Ontology Annotation ({database})",
                    source_release=source_release,
                    source_url=source_url,
                    evidence_code=evidence_code,
                    reference_ids=tuple(references.split("|")),
                    qualifier=qualifier or None,
                    taxon_id=primary_taxon,
                    source_record_fields=tuple(fields),
                )
            )

    return GOAnnotationImportResult(
        annotations=tuple(annotations),
        rows_read=rows_read,
        row_errors=tuple(errors),
        rows_filtered_by_taxon=rows_filtered_by_taxon,
    )
