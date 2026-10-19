"""Rank candidate diseases using patient HPO phenotypes."""

from __future__ import annotations

from dataclasses import dataclass

from rare_disease_ai.disease_phenotypes import DiseasePhenotypeProfile
from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.phenotype_information_content import (
    calculate_information_content,
)
from rare_disease_ai.phenotype_semantic_similarity import (
    best_match_average_similarity,
)

EXACT_MATCH_SCORE = 1.0
RELATED_MATCH_SCORE = 0.5


@dataclass(frozen=True)
class DiseaseRankingResult:
    """Auditable ranking result for one candidate disease."""

    database_id: str
    disease_name: str
    score: float
    semantic_score: float | None = None
    exact_matches: tuple[str, ...] = ()
    related_matches: tuple[str, ...] = ()
    unmatched_patient_hpo_ids: tuple[str, ...] = ()
    negative_matches: tuple[str, ...] = ()


def rank_diseases(
    patient_hpo_ids: tuple[str, ...],
    disease_index: dict[str, DiseasePhenotypeProfile],
    hpo_graph: HPOGraph,
    *,
    top_k: int | None = None,
    use_information_content: bool = False,
    use_semantic_similarity: bool = False,
) -> tuple[DiseaseRankingResult, ...]:
    """Rank diseases using exact and hierarchical phenotype matches.

    Exact matches receive a score of 1.0, and related ancestor/descendant
    matches receive 0.5. When enabled, information content weights patient
    phenotypes by their annotation-based specificity.

    When semantic similarity is enabled, an additional ontology-based score
    is calculated and reported. It does not change the default sort order.

    Negative annotations are reported separately and do not alter the score.
    All scores are research-ranking signals, not clinical probabilities.
    """
    unique_patient_hpo_ids = tuple(dict.fromkeys(patient_hpo_ids))

    if not unique_patient_hpo_ids:
        raise ValueError("At least one patient HPO identifier is required.")

    unknown_ids = sorted(set(unique_patient_hpo_ids) - set(hpo_graph.ontology.terms))
    if unknown_ids:
        raise ValueError(f"Unknown patient HPO identifiers: {unknown_ids}")

    if top_k is not None and top_k < 1:
        raise ValueError("top_k must be at least 1 or None.")

    needs_information_content = use_information_content or use_semantic_similarity

    information_content = (
        calculate_information_content(disease_index, hpo_graph) if needs_information_content else {}
    )

    total_patient_information = sum(
        information_content.get(hpo_id, 0.0) for hpo_id in unique_patient_hpo_ids
    )

    patient_ancestors = {
        hpo_id: set(hpo_graph.get_ancestors(hpo_id)) for hpo_id in unique_patient_hpo_ids
    }

    # Cache disease-term ancestors to avoid repeated graph traversal.
    disease_term_ancestors: dict[str, set[str]] = {}

    results: list[DiseaseRankingResult] = []

    for database_id, profile in disease_index.items():
        exact_matches: list[str] = []
        related_matches: list[str] = []
        unmatched: list[str] = []
        total_score = 0.0

        for patient_hpo_id in unique_patient_hpo_ids:
            patient_weight = information_content[patient_hpo_id] if use_information_content else 1.0

            if patient_hpo_id in profile.positive_hpo_ids:
                exact_matches.append(patient_hpo_id)
                total_score += patient_weight * EXACT_MATCH_SCORE
                continue

            patient_ancestors_for_term = patient_ancestors[patient_hpo_id]
            related = False

            for disease_hpo_id in profile.positive_hpo_ids:
                if disease_hpo_id not in hpo_graph.ontology.terms:
                    continue

                if disease_hpo_id not in disease_term_ancestors:
                    disease_term_ancestors[disease_hpo_id] = set(
                        hpo_graph.get_ancestors(disease_hpo_id)
                    )

                disease_ancestors = disease_term_ancestors[disease_hpo_id]

                if (
                    disease_hpo_id in patient_ancestors_for_term
                    or patient_hpo_id in disease_ancestors
                ):
                    related = True
                    break

            if related:
                related_matches.append(patient_hpo_id)
                total_score += patient_weight * RELATED_MATCH_SCORE
            else:
                unmatched.append(patient_hpo_id)

        denominator = (
            total_patient_information if use_information_content else len(unique_patient_hpo_ids)
        )

        score = total_score / denominator if denominator > 0 else 0.0

        negative_matches = tuple(sorted(set(unique_patient_hpo_ids) & profile.negative_hpo_ids))

        semantic_score = None
        if use_semantic_similarity:
            semantic_score = best_match_average_similarity(
                unique_patient_hpo_ids,
                profile.positive_hpo_ids,
                information_content,
                hpo_graph,
            )

        results.append(
            DiseaseRankingResult(
                database_id=database_id,
                disease_name=profile.disease_name,
                score=score,
                semantic_score=semantic_score,
                exact_matches=tuple(exact_matches),
                related_matches=tuple(related_matches),
                unmatched_patient_hpo_ids=tuple(unmatched),
                negative_matches=negative_matches,
            )
        )

    # Preserve the established ranking until semantic-score evaluation.
    if use_semantic_similarity:
        results.sort(
            key=lambda result: (
                -(result.semantic_score or 0.0),
                result.database_id,
            )
        )
    else:
        results.sort(key=lambda result: (-result.score, result.database_id))

    if top_k is not None:
        return tuple(results[:top_k])

    return tuple(results)
