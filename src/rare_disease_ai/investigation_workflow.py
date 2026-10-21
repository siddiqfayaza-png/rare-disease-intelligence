"""End-to-end research workflow for rare-disease evidence investigation.

This module orchestrates existing components. It does not diagnose patients,
infer variant pathogenicity, or combine phenotype and variant scores.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import asdict, dataclass

from rare_disease_ai.case_validation import PatientCase
from rare_disease_ai.disease_phenotypes import DiseasePhenotypeProfile
from rare_disease_ai.gencc_evidence_adapter import (
    GenCCAdapterResult,
    adapt_gencc_assertions,
)
from rare_disease_ai.gencc_import import GenCCAssertion
from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology
from rare_disease_ai.mechanism_evidence import MechanismEvidence
from rare_disease_ai.mechanism_investigation import (
    MechanismInvestigationReport,
    build_mechanism_index,
    investigate_gene_mechanisms,
    mechanism_report_to_dict,
)
from rare_disease_ai.phenotype_normalization import (
    NormalizedPatientCase,
    normalize_patient_case,
)
from rare_disease_ai.phenotype_ranking import (
    DiseaseRankingResult,
    rank_diseases,
)
from rare_disease_ai.phenotype_variant_integration import (
    IntegratedDiseaseResult,
    integrate_ranked_evidence,
    integration_report_to_dict,
)
from rare_disease_ai.variant_prioritization import (
    VariantRankingResult,
    VariantRecord,
    prioritize_variants,
)


@dataclass(frozen=True)
class InvestigationWorkflowResult:
    """Structured outputs and audit information from one investigation."""

    case_id: str
    normalized_case: NormalizedPatientCase
    disease_rankings: tuple[DiseaseRankingResult, ...]
    variant_rankings: tuple[VariantRankingResult, ...]
    gencc_adapter_result: GenCCAdapterResult
    integrated_diseases: tuple[IntegratedDiseaseResult, ...]
    mechanism_report: MechanismInvestigationReport
    research_only: bool = True


def run_investigation(
    case: PatientCase,
    ontology: HPOOntology,
    disease_index: dict[str, DiseasePhenotypeProfile],
    hpo_graph: HPOGraph,
    variants: tuple[VariantRecord, ...],
    gencc_assertions: tuple[GenCCAssertion, ...],
    mechanism_annotations: Iterable[MechanismEvidence],
    *,
    top_k: int | None = 20,
    max_allele_frequency: float | None = 0.01,
) -> InvestigationWorkflowResult:
    """Run phenotype ranking, variant prioritization, and evidence lookup.

    Inputs are already loaded into memory. Data ingestion and file handling
    remain the responsibility of the existing importer modules.
    """
    if top_k is not None and top_k < 1:
        raise ValueError("top_k must be at least 1 or None.")

    normalized_case = normalize_patient_case(case, ontology)

    # Only unambiguous, resolved HPO IDs enter phenotype ranking.
    resolved_hpo_ids = tuple(
        dict.fromkeys(
            record.hpo_id
            for record in normalized_case.records
            if record.status == "resolved" and record.hpo_id is not None
        )
    )

    if resolved_hpo_ids:
        disease_rankings = rank_diseases(
            resolved_hpo_ids,
            disease_index,
            hpo_graph,
            top_k=top_k,
        )
    else:
        disease_rankings = ()

    variant_rankings = prioritize_variants(
        variants,
        max_allele_frequency=max_allele_frequency,
    )

    gencc_result = adapt_gencc_assertions(
        gencc_assertions,
        set(disease_index),
    )

    integrated_diseases = integrate_ranked_evidence(
        disease_rankings,
        variant_rankings,
        gencc_result.evidence,
    )

    mechanism_index = build_mechanism_index(mechanism_annotations)
    genes_to_investigate = tuple(
        dict.fromkeys(result.variant.gene_symbol for result in variant_rankings)
    )
    mechanism_report = investigate_gene_mechanisms(
        genes_to_investigate,
        mechanism_index,
    )

    return InvestigationWorkflowResult(
        case_id=case.case_id,
        normalized_case=normalized_case,
        disease_rankings=disease_rankings,
        variant_rankings=variant_rankings,
        gencc_adapter_result=gencc_result,
        integrated_diseases=integrated_diseases,
        mechanism_report=mechanism_report,
    )


def investigation_report_to_dict(
    result: InvestigationWorkflowResult,
) -> dict[str, object]:
    """Create a JSON-serializable report with provenance and audit counts."""
    normalization_counts: dict[str, int] = {}
    for record in result.normalized_case.records:
        normalization_counts[record.status] = normalization_counts.get(record.status, 0) + 1

    return {
        "case_id": result.case_id,
        "research_only": True,
        "limitations": [
            "Research prototype; not for clinical diagnosis or treatment decisions.",
            "Phenotype scores are heuristic ranking signals, not probabilities.",
            "Variant scores are heuristic prioritization signals, not pathogenicity estimates.",
            "Phenotype and variant scores are kept separate and are not combined.",
            "Gene-disease assertions are preserved as evidence, not treated as proof of causality.",
            "Gene Ontology annotations describe gene-level annotations and do not establish "
            "that a variant causes a disease.",
            "Conflicting assertions are not automatically resolved.",
            "Unresolved, ambiguous, and obsolete phenotype inputs are retained in the audit "
            "but are not passed to phenotype ranking.",
        ],
        "normalization": {
            "records": [asdict(record) for record in result.normalized_case.records],
            "counts_by_status": normalization_counts,
            "resolved_unique_hpo_ids": list(
                dict.fromkeys(
                    record.hpo_id
                    for record in result.normalized_case.records
                    if record.status == "resolved" and record.hpo_id is not None
                )
            ),
        },
        "phenotype_rankings": [asdict(item) for item in result.disease_rankings],
        "variant_rankings": [
            {
                "variant": asdict(item.variant),
                "score": item.score,
                "reasons": list(item.reasons),
            }
            for item in result.variant_rankings
        ],
        "gencc_audit": {
            "matched_assertions": result.gencc_adapter_result.matched_assertions,
            "unmatched_assertions": result.gencc_adapter_result.unmatched_assertions,
            "matched_by_original_id": result.gencc_adapter_result.matched_by_original_id,
            "matched_by_canonical_id": result.gencc_adapter_result.matched_by_canonical_id,
            "evidence_assertions_retained": len(result.gencc_adapter_result.evidence),
        },
        "integrated_evidence": integration_report_to_dict(result.integrated_diseases),
        "mechanisms": mechanism_report_to_dict(result.mechanism_report),
    }
