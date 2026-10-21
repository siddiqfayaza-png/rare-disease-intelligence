"""Gene-indexed biological mechanism investigation with provenance preservation."""

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass

from rare_disease_ai.mechanism_evidence import MechanismEvidence


@dataclass(frozen=True)
class GeneMechanismResult:
    """GO mechanism annotations retrieved for one requested gene."""

    gene_symbol: str
    annotations: tuple[MechanismEvidence, ...]

    @property
    def annotation_count(self) -> int:
        return len(self.annotations)


@dataclass(frozen=True)
class MechanismInvestigationReport:
    """Structured, research-only mechanism lookup report."""

    results: tuple[GeneMechanismResult, ...]
    unmatched_gene_symbols: tuple[str, ...]
    requested_gene_count: int
    annotations_returned: int
    research_only: bool = True


def _gene_key(gene_symbol: str) -> str:
    """Create a case-insensitive lookup key without altering displayed symbols."""
    return gene_symbol.strip().casefold()


def _annotation_sort_key(
    annotation: MechanismEvidence,
) -> tuple[str, str, str, str, str, str]:
    """Return a stable ordering while retaining every annotation record."""
    return (
        annotation.annotation_type,
        annotation.annotation_id,
        annotation.source_name,
        annotation.source_release,
        annotation.evidence_code or "",
        "|".join(annotation.reference_ids),
    )


def build_mechanism_index(
    annotations: Iterable[MechanismEvidence],
) -> dict[str, tuple[MechanismEvidence, ...]]:
    """Index annotations by normalized gene symbol.

    The index key is case-insensitive. Annotation objects are retained
    unchanged, including qualifiers, evidence codes, references and sources.
    No evidence is combined into a score.
    """
    grouped: dict[str, list[MechanismEvidence]] = defaultdict(list)

    for position, annotation in enumerate(annotations):
        if not isinstance(annotation, MechanismEvidence):
            raise TypeError(f"Annotation at position {position} must be MechanismEvidence")

        key = _gene_key(annotation.gene_symbol)
        if not key:
            raise ValueError(f"Annotation at position {position} has an empty gene symbol")

        grouped[key].append(annotation)

    return {
        key: tuple(sorted(records, key=_annotation_sort_key))
        for key, records in sorted(grouped.items())
    }


def investigate_gene_mechanisms(
    gene_symbols: Iterable[str],
    mechanism_index: dict[str, tuple[MechanismEvidence, ...]],
) -> MechanismInvestigationReport:
    """Retrieve source-backed GO annotations for requested genes.

    Genes without annotations are reported explicitly. Duplicate requested
    gene symbols are collapsed case-insensitively, preserving first occurrence.
    This function retrieves annotations; it does not infer disease causality.
    """
    requested: list[str] = []
    seen: set[str] = set()

    for position, gene_symbol in enumerate(gene_symbols):
        if not isinstance(gene_symbol, str) or not gene_symbol.strip():
            raise ValueError(f"Gene symbol at position {position} must be a non-empty string")

        key = _gene_key(gene_symbol)
        if key not in seen:
            seen.add(key)
            requested.append(gene_symbol.strip())

    results: list[GeneMechanismResult] = []
    unmatched: list[str] = []
    total_annotations = 0

    for gene_symbol in requested:
        key = _gene_key(gene_symbol)
        annotations = mechanism_index.get(key, ())

        if not annotations:
            unmatched.append(gene_symbol)
            continue

        results.append(
            GeneMechanismResult(
                gene_symbol=gene_symbol,
                annotations=annotations,
            )
        )
        total_annotations += len(annotations)

    return MechanismInvestigationReport(
        results=tuple(results),
        unmatched_gene_symbols=tuple(unmatched),
        requested_gene_count=len(requested),
        annotations_returned=total_annotations,
    )


def mechanism_report_to_dict(
    report: MechanismInvestigationReport,
) -> dict[str, object]:
    """Convert a mechanism report to a JSON-serializable dictionary."""
    return {
        "report_type": "gene_mechanism_investigation",
        "research_only": report.research_only,
        "requested_gene_count": report.requested_gene_count,
        "matched_gene_count": len(report.results),
        "unmatched_gene_symbols": list(report.unmatched_gene_symbols),
        "annotations_returned": report.annotations_returned,
        "results": [
            {
                "gene_symbol": result.gene_symbol,
                "annotation_count": result.annotation_count,
                "annotations": [
                    {
                        "annotation_id": annotation.annotation_id,
                        "annotation_label": annotation.annotation_label,
                        "annotation_type": annotation.annotation_type,
                        "source_name": annotation.source_name,
                        "source_release": annotation.source_release,
                        "source_url": annotation.source_url,
                        "evidence_code": annotation.evidence_code,
                        "evidence_label": annotation.evidence_label,
                        "reference_ids": list(annotation.reference_ids),
                        "qualifier": annotation.qualifier,
                        "taxon_id": annotation.taxon_id,
                    }
                    for annotation in result.annotations
                ],
            }
            for result in report.results
        ],
        "interpretation_limits": [
            "GO annotations describe recorded gene-associated functions or processes.",
            "Annotations do not establish that a variant causes a disease.",
            "Evidence codes and qualifiers must be considered in context.",
            "Missing annotations do not establish absence of a biological mechanism.",
        ],
    }
