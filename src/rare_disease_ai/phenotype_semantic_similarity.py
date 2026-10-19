"""Ontology-based semantic similarity between patient and disease phenotypes."""

from __future__ import annotations

from rare_disease_ai.hpo_graph import HPOGraph


def best_match_similarity(
    first_hpo_id: str,
    second_hpo_id: str,
    information_content: dict[str, float],
    hpo_graph: HPOGraph,
) -> float:
    """Return the information content of the most informative shared ancestor."""
    known_terms = hpo_graph.ontology.terms
    if first_hpo_id not in known_terms or second_hpo_id not in known_terms:
        raise ValueError("Both HPO identifiers must exist in the ontology.")

    first_ancestors = set(hpo_graph.get_ancestors(first_hpo_id))
    first_ancestors.add(first_hpo_id)

    second_ancestors = set(hpo_graph.get_ancestors(second_hpo_id))
    second_ancestors.add(second_hpo_id)

    shared_terms = first_ancestors & second_ancestors
    if not shared_terms:
        return 0.0

    return max(information_content.get(term, 0.0) for term in shared_terms)


def best_match_average_similarity(
    patient_hpo_ids: tuple[str, ...],
    disease_hpo_ids: frozenset[str],
    information_content: dict[str, float],
    hpo_graph: HPOGraph,
) -> float:
    """Average best-match similarities in both patient-to-disease directions.

    This is an experimental ranking metric, not a diagnostic probability.
    """
    patient_terms = tuple(dict.fromkeys(patient_hpo_ids))
    disease_terms = tuple(term for term in disease_hpo_ids if term in hpo_graph.ontology.terms)

    if not patient_terms or not disease_terms:
        return 0.0

    patient_to_disease = [
        max(
            best_match_similarity(patient_term, disease_term, information_content, hpo_graph)
            for disease_term in disease_terms
        )
        for patient_term in patient_terms
    ]

    disease_to_patient = [
        max(
            best_match_similarity(patient_term, disease_term, information_content, hpo_graph)
            for patient_term in patient_terms
        )
        for disease_term in disease_terms
    ]

    total_information = sum(information_content.get(term, 0.0) for term in patient_terms)
    total_information += sum(information_content.get(term, 0.0) for term in disease_terms)

    if total_information <= 0:
        return 0.0

    matched_information = sum(patient_to_disease) + sum(disease_to_patient)
    return matched_information / total_information
