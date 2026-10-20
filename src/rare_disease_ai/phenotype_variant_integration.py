"""Integrate phenotype rankings with explicit gene-disease evidence.

Research prototype only. This module does not combine heuristic scores or
make clinical diagnoses.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from rare_disease_ai.phenotype_ranking import DiseaseRankingResult
from rare_disease_ai.variant_prioritization import VariantRankingResult

EVIDENCE_CATEGORIES = (
    "supporting",
    "disputed",
    "refuting",
    "no_known_relationship",
    "unspecified",
    "unrecognized",
)


@dataclass(frozen=True)
class GeneDiseaseEvidence:
    """A traceable gene-disease assertion with structured provenance."""

    disease_id: str
    gene_symbol: str
    source: str
    evidence_id: str
    classification: str | None = None
    submitter: str | None = None
    canonical_disease_id: str | None = None
    original_disease_id: str | None = None
    pmids: tuple[str, ...] = ()
    public_report_url: str | None = None
    assertion_criteria_url: str | None = None

    @property
    def evidence_category(self) -> str:
        """Categorize evidence without making a clinical conclusion."""
        categories = {
            "definitive": "supporting",
            "strong": "supporting",
            "moderate": "supporting",
            "limited": "supporting",
            "supportive": "supporting",
            "disputed evidence": "disputed",
            "refuted evidence": "refuting",
            "no known disease relationship": "no_known_relationship",
        }

        if self.classification is None or not self.classification.strip():
            return "unspecified"

        return categories.get(
            self.classification.strip().casefold(),
            "unrecognized",
        )


@dataclass(frozen=True)
class VariantEvidenceMatch:
    """A ranked variant linked through a supplied gene-disease relationship."""

    variant_rank: int
    gene_symbol: str
    variant_score: float
    consequence: str
    evidence_source: str
    evidence_id: str
    classification: str | None = None
    evidence_category: str = "unspecified"
    submitter: str | None = None
    canonical_disease_id: str | None = None
    original_disease_id: str | None = None
    pmids: tuple[str, ...] = ()
    public_report_url: str | None = None
    assertion_criteria_url: str | None = None


@dataclass(frozen=True)
class IntegratedDiseaseResult:
    """Disease phenotype ranking plus separately reported variant evidence."""

    disease_id: str
    disease_name: str
    phenotype_rank: int
    phenotype_score: float
    matched_variants: tuple[VariantEvidenceMatch, ...]
    evidence_summary: tuple[tuple[str, int], ...] = ()


def integrate_ranked_evidence(
    disease_rankings: tuple[DiseaseRankingResult, ...],
    variant_rankings: tuple[VariantRankingResult, ...],
    gene_disease_evidence: tuple[GeneDiseaseEvidence, ...],
) -> tuple[IntegratedDiseaseResult, ...]:
    """Attach variants through explicit gene-disease mappings.

    Phenotype and variant scores remain separate. Each evidence assertion is
    preserved independently, including conflicting classifications.
    Evidence summaries count assertions by category; they do not resolve
    conflicts, establish variant pathogenicity, or make clinical conclusions.
    """
    evidence_by_pair: dict[tuple[str, str], list[GeneDiseaseEvidence]] = {}

    for evidence in gene_disease_evidence:
        required_values = (
            evidence.disease_id,
            evidence.gene_symbol,
            evidence.source,
            evidence.evidence_id,
        )

        if not all(isinstance(value, str) and value.strip() for value in required_values):
            raise ValueError(
                "Gene-disease evidence requires disease_id, gene_symbol, source, and evidence_id."
            )

        pair = (
            evidence.disease_id.strip(),
            evidence.gene_symbol.strip().upper(),
        )
        evidence_by_pair.setdefault(pair, []).append(evidence)

    results: list[IntegratedDiseaseResult] = []

    for disease_rank, disease in enumerate(disease_rankings, start=1):
        matches: list[VariantEvidenceMatch] = []

        for variant_rank, ranked_variant in enumerate(variant_rankings, start=1):
            variant = ranked_variant.variant
            pair = (
                disease.database_id,
                variant.gene_symbol.strip().upper(),
            )

            for evidence in evidence_by_pair.get(pair, []):
                matches.append(
                    VariantEvidenceMatch(
                        variant_rank=variant_rank,
                        gene_symbol=variant.gene_symbol,
                        variant_score=ranked_variant.score,
                        consequence=variant.consequence,
                        evidence_source=evidence.source,
                        evidence_id=evidence.evidence_id,
                        classification=evidence.classification,
                        evidence_category=evidence.evidence_category,
                        submitter=evidence.submitter,
                        canonical_disease_id=evidence.canonical_disease_id,
                        original_disease_id=evidence.original_disease_id,
                        pmids=evidence.pmids,
                        public_report_url=evidence.public_report_url,
                        assertion_criteria_url=evidence.assertion_criteria_url,
                    )
                )

        category_counts = {
            category: sum(match.evidence_category == category for match in matches)
            for category in EVIDENCE_CATEGORIES
        }

        results.append(
            IntegratedDiseaseResult(
                disease_id=disease.database_id,
                disease_name=disease.disease_name,
                phenotype_rank=disease_rank,
                phenotype_score=disease.score,
                matched_variants=tuple(matches),
                evidence_summary=tuple(category_counts.items()),
            )
        )

    return tuple(results)


def integration_report_to_dict(
    results: tuple[IntegratedDiseaseResult, ...],
) -> dict[str, object]:
    """Convert integration results to a JSON-serializable research report."""
    diseases = []

    for result in results:
        disease = asdict(result)
        disease["matched_variants"] = [asdict(match) for match in result.matched_variants]
        disease["evidence_summary"] = dict(result.evidence_summary)
        diseases.append(disease)

    return {
        "research_only": True,
        "score_policy": (
            "Phenotype and variant heuristic scores are reported separately; "
            "they are not combined and are not clinical probabilities. "
            "GenCC classifications are preserved as evidence categories, "
            "not clinical conclusions. Conflicting assertions are not resolved. "
            "Evidence summaries count assertions by category and do not "
            "establish causality or variant pathogenicity."
        ),
        "evidence_category_definitions": {
            "supporting": (
                "Assertions classified as Definitive, Strong, Moderate, Limited, or Supportive."
            ),
            "disputed": "Assertions classified as Disputed Evidence.",
            "refuting": "Assertions classified as Refuted Evidence.",
            "no_known_relationship": ("Assertions classified as No Known Disease Relationship."),
            "unspecified": "Assertions without a classification.",
            "unrecognized": ("Assertions whose classification is not mapped by this module."),
        },
        "diseases": diseases,
    }
