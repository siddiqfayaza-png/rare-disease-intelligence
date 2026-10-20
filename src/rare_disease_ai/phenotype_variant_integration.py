"""Integrate phenotype rankings with explicit gene-disease evidence.

Research prototype only. This module does not combine heuristic scores or
make clinical diagnoses.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from rare_disease_ai.phenotype_ranking import DiseaseRankingResult
from rare_disease_ai.variant_prioritization import VariantRankingResult


@dataclass(frozen=True)
class GeneDiseaseEvidence:
    """An explicit, traceable gene-disease relationship."""

    disease_id: str
    gene_symbol: str
    source: str
    evidence_id: str


@dataclass(frozen=True)
class VariantEvidenceMatch:
    """A ranked variant linked through a supplied gene-disease relationship."""

    variant_rank: int
    gene_symbol: str
    variant_score: float
    consequence: str
    evidence_source: str
    evidence_id: str


@dataclass(frozen=True)
class IntegratedDiseaseResult:
    """Disease phenotype ranking plus separately reported variant evidence."""

    disease_id: str
    disease_name: str
    phenotype_rank: int
    phenotype_score: float
    matched_variants: tuple[VariantEvidenceMatch, ...]


def integrate_ranked_evidence(
    disease_rankings: tuple[DiseaseRankingResult, ...],
    variant_rankings: tuple[VariantRankingResult, ...],
    gene_disease_evidence: tuple[GeneDiseaseEvidence, ...],
) -> tuple[IntegratedDiseaseResult, ...]:
    """Attach variants to diseases only through explicit gene-disease mappings.

    Phenotype and variant scores remain separate. Disease ordering follows
    the phenotype ranker; variant ordering follows the variant ranker.
    """
    evidence_by_pair: dict[tuple[str, str], list[GeneDiseaseEvidence]] = {}

    for evidence in gene_disease_evidence:
        if not all(
            value.strip()
            for value in (
                evidence.disease_id,
                evidence.gene_symbol,
                evidence.source,
                evidence.evidence_id,
            )
        ):
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
            pair = (disease.database_id, variant.gene_symbol.strip().upper())

            for evidence in evidence_by_pair.get(pair, []):
                matches.append(
                    VariantEvidenceMatch(
                        variant_rank=variant_rank,
                        gene_symbol=variant.gene_symbol,
                        variant_score=ranked_variant.score,
                        consequence=variant.consequence,
                        evidence_source=evidence.source,
                        evidence_id=evidence.evidence_id,
                    )
                )

        results.append(
            IntegratedDiseaseResult(
                disease_id=disease.database_id,
                disease_name=disease.disease_name,
                phenotype_rank=disease_rank,
                phenotype_score=disease.score,
                matched_variants=tuple(matches),
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
        diseases.append(disease)

    return {
        "research_only": True,
        "score_policy": (
            "Phenotype and variant heuristic scores are reported separately; "
            "they are not combined and are not clinical probabilities."
        ),
        "diseases": diseases,
    }
