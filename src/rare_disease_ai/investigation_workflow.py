"""End-to-end research workflow for rare-disease evidence investigation.

This module orchestrates existing components. It does not diagnose patients,
infer variant pathogenicity, or combine phenotype and variant scores.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import asdict, dataclass

from rare_disease_ai.biomedical_evidence_retrieval import (
    EvidenceRetrievalReport,
)
from rare_disease_ai.case_validation import PatientCase
from rare_disease_ai.disease_phenotypes import DiseasePhenotypeProfile
from rare_disease_ai.evidence_retrieval_service import EvidenceRetrievalService
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
    retrieved_evidence: EvidenceRetrievalReport | None = None


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
    evidence_retrieval_service: EvidenceRetrievalService | None = None,
    evidence_query: str | None = None,
    evidence_top_k: int = 5,
    evidence_entity_ids: tuple[str, ...] = (),
    evidence_types: tuple[str, ...] = (),
) -> InvestigationWorkflowResult:
    """Run phenotype ranking, variant prioritization, and evidence lookup.

    Inputs are already loaded into memory. Data ingestion and file handling
    remain the responsibility of the existing importer modules.

    Evidence retrieval is optional. When enabled, it runs an explicit query
    through the supplied retrieval service and preserves source provenance.
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

    # Keep general evidence retrieval separate from the typed GO mechanism
    # investigation. The query is explicit; no entity-ID format is assumed.
    if evidence_retrieval_service is None:
        if evidence_query is not None:
            raise ValueError(
                "evidence_retrieval_service is required when evidence_query is provided"
            )
        retrieved_evidence = None
    else:
        if evidence_query is None or not evidence_query.strip():
            raise ValueError(
                "A non-empty evidence_query is required when evidence retrieval is enabled"
            )

        retrieved_evidence = evidence_retrieval_service.search(
            evidence_query,
            top_k=evidence_top_k,
            entity_ids=evidence_entity_ids,
            evidence_types=evidence_types,
        )

    return InvestigationWorkflowResult(
        case_id=case.case_id,
        normalized_case=normalized_case,
        disease_rankings=disease_rankings,
        variant_rankings=variant_rankings,
        gencc_adapter_result=gencc_result,
        integrated_diseases=integrated_diseases,
        mechanism_report=mechanism_report,
        retrieved_evidence=retrieved_evidence,
    )


def investigation_report_to_dict(
    result: InvestigationWorkflowResult,
) -> dict[str, object]:
    """Create a JSON-serializable report with provenance and audit counts."""
    normalization_counts: dict[str, int] = {}
    for record in result.normalized_case.records:
        normalization_counts[record.status] = normalization_counts.get(record.status, 0) + 1

    retrieved_evidence_report: dict[str, object] | None = None

    if result.retrieved_evidence is not None:
        evidence_report = result.retrieved_evidence

        retrieved_evidence_report = {
            "research_only": evidence_report.research_only,
            "query": evidence_report.query,
            "total_documents": evidence_report.total_documents,
            "eligible_documents": evidence_report.eligible_documents,
            "returned_count": len(evidence_report.hits),
            "score_interpretation": (
                "Lexical query-token coverage; not diagnostic confidence or evidence quality."
            ),
            "hits": [
                {
                    "evidence_id": hit.document.evidence_id,
                    "title": hit.document.title,
                    "text": hit.document.text,
                    "evidence_type": hit.document.evidence_type,
                    "entity_ids": list(hit.document.entity_ids),
                    "reference_ids": list(hit.document.reference_ids),
                    "source_name": hit.document.source_name,
                    "source_release": hit.document.source_release,
                    "source_url": hit.document.source_url,
                    "score": hit.score,
                    "matched_terms": list(hit.matched_terms),
                }
                for hit in evidence_report.hits
            ],
        }

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
            "Retrieved evidence is research material and requires source-level interpretation.",
            "Evidence retrieval scores reflect lexical query-token coverage, not diagnostic "
            "confidence, evidence quality, or causal support.",
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
        "retrieved_evidence": retrieved_evidence_report,
    }
