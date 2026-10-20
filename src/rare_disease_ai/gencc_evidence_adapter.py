"""Convert GenCC assertions to explicit, identifier-matched evidence mappings."""

from __future__ import annotations

from dataclasses import dataclass

from rare_disease_ai.gencc_import import GenCCAssertion
from rare_disease_ai.phenotype_variant_integration import GeneDiseaseEvidence


@dataclass(frozen=True)
class GenCCAdapterResult:
    """Summary of GenCC assertions matched to indexed disease identifiers."""

    evidence: tuple[GeneDiseaseEvidence, ...]
    matched_assertions: int
    unmatched_assertions: int
    matched_by_original_id: int
    matched_by_canonical_id: int


def adapt_gencc_assertions(
    assertions: tuple[GenCCAssertion, ...],
    indexed_disease_ids: set[str] | frozenset[str],
) -> GenCCAdapterResult:
    """Match exact disease identifiers, preferring the original identifier.

    Disease names are never used for matching. Each matched assertion is
    preserved separately, including its classification and provenance.
    No assertions are merged, no scores are combined, and no clinical
    conclusions are inferred.
    """
    evidence: list[GeneDiseaseEvidence] = []
    unmatched = 0
    original_matches = 0
    canonical_matches = 0

    for assertion in assertions:
        original_id = assertion.disease_original_curie.strip()
        canonical_id = assertion.disease_curie.strip()

        if original_id and original_id in indexed_disease_ids:
            matched_id = original_id
            original_matches += 1
        elif canonical_id and canonical_id in indexed_disease_ids:
            matched_id = canonical_id
            canonical_matches += 1
        else:
            unmatched += 1
            continue

        evidence.append(
            GeneDiseaseEvidence(
                disease_id=matched_id,
                gene_symbol=assertion.gene_symbol,
                source="GenCC",
                evidence_id=assertion.submission_id,
                classification=assertion.classification_title,
                submitter=assertion.submitter_title,
                canonical_disease_id=canonical_id or None,
                original_disease_id=original_id or None,
                pmids=assertion.pmids,
                public_report_url=assertion.public_report_url or None,
                assertion_criteria_url=assertion.assertion_criteria_url or None,
            )
        )

    return GenCCAdapterResult(
        evidence=tuple(evidence),
        matched_assertions=len(evidence),
        unmatched_assertions=unmatched,
        matched_by_original_id=original_matches,
        matched_by_canonical_id=canonical_matches,
    )
