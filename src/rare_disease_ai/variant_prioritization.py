"""Validate and prioritize candidate genetic variants.

Research prototype only. Scores are heuristic and are not clinical
pathogenicity classifications or diagnostic conclusions.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class VariantRecord:
    """A normalized candidate variant record."""

    chromosome: str
    position: int
    reference: str
    alternate: str
    gene_symbol: str
    consequence: str
    zygosity: str
    allele_frequency: float | None = None


@dataclass(frozen=True)
class VariantValidationResult:
    """Validation outcome for a candidate variant."""

    is_valid: bool
    errors: tuple[str, ...]


@dataclass(frozen=True)
class VariantRankingResult:
    """Transparent heuristic ranking result."""

    variant: VariantRecord
    score: float
    reasons: tuple[str, ...]


# These are configurable prototype heuristics, not clinical standards.
CONSEQUENCE_SCORES = {
    "transcript_ablation": 5.0,
    "frameshift_variant": 5.0,
    "stop_gained": 5.0,
    "splice_acceptor_variant": 5.0,
    "splice_donor_variant": 5.0,
    "start_lost": 4.0,
    "stop_lost": 3.0,
    "missense_variant": 2.0,
    "inframe_insertion": 2.0,
    "inframe_deletion": 2.0,
    "synonymous_variant": 0.5,
    "intron_variant": 0.5,
    "intergenic_variant": 0.0,
}


def validate_variant(variant: VariantRecord) -> VariantValidationResult:
    """Check required fields and basic variant data consistency."""
    errors: list[str] = []

    if not variant.chromosome.strip():
        errors.append("Chromosome is required.")

    if variant.position < 1:
        errors.append("Position must be a positive integer.")

    if not variant.reference or not variant.alternate:
        errors.append("Reference and alternate alleles are required.")
    elif variant.reference.upper() == variant.alternate.upper():
        errors.append("Reference and alternate alleles must differ.")

    if not variant.gene_symbol.strip():
        errors.append("Gene symbol is required.")

    if variant.consequence not in CONSEQUENCE_SCORES:
        errors.append("Consequence is not supported by the prototype.")

    if variant.zygosity.lower() not in {
        "heterozygous",
        "homozygous",
        "hemizygous",
        "unknown",
    }:
        errors.append("Zygosity must be heterozygous, homozygous, hemizygous, or unknown.")

    frequency = variant.allele_frequency
    if frequency is not None and (not math.isfinite(frequency) or not 0.0 <= frequency <= 1.0):
        errors.append("Allele frequency must be between 0 and 1.")

    return VariantValidationResult(
        is_valid=not errors,
        errors=tuple(errors),
    )


def prioritize_variants(
    variants: tuple[VariantRecord, ...],
    *,
    max_allele_frequency: float | None = 0.01,
) -> tuple[VariantRankingResult, ...]:
    """Rank valid variants using transparent consequence/frequency heuristics.

    Variants with missing allele frequency are retained but receive no
    rarity bonus. No variant is classified as pathogenic or benign.
    """
    if max_allele_frequency is not None and not (
        math.isfinite(max_allele_frequency) and 0.0 <= max_allele_frequency <= 1.0
    ):
        raise ValueError("max_allele_frequency must be between 0 and 1.")

    results: list[VariantRankingResult] = []

    for variant in variants:
        validation = validate_variant(variant)
        if not validation.is_valid:
            continue

        score = CONSEQUENCE_SCORES[variant.consequence]
        reasons = [f"Consequence heuristic: {variant.consequence}"]

        frequency = variant.allele_frequency
        if frequency is None:
            reasons.append("Allele frequency unavailable; no rarity bonus.")
        elif max_allele_frequency is not None and frequency <= max_allele_frequency:
            score += 1.0
            reasons.append("Allele frequency meets the configured rarity threshold.")
        else:
            reasons.append("No rarity bonus applied.")

        results.append(
            VariantRankingResult(
                variant=variant,
                score=score,
                reasons=tuple(reasons),
            )
        )

    return tuple(
        sorted(
            results,
            key=lambda result: (
                -result.score,
                result.variant.gene_symbol,
                result.variant.chromosome,
                result.variant.position,
                result.variant.reference,
                result.variant.alternate,
            ),
        )
    )
