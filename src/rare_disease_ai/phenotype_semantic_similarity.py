"""Ontology-based semantic similarity between patient and disease phenotypes."""

from __future__ import annotations

from rare_disease_ai.hpo_graph import HPOGraph


def _term_ancestors(
    hpo_id: str,
    hpo_graph: HPOGraph,
    ancestor_cache: dict[str, frozenset[str]] | None = None,
) -> frozenset[str]:
    """Return a term and its ancestors, optionally reusing cached results."""
    if ancestor_cache is not None and hpo_id in ancestor_cache:
        return ancestor_cache[hpo_id]

    ancestors = frozenset((*hpo_graph.get_ancestors(hpo_id), hpo_id))

    if ancestor_cache is not None:
        ancestor_cache[hpo_id] = ancestors

    return ancestors


def _cached_pair_similarity(
    first_hpo_id: str,
    second_hpo_id: str,
    information_content: dict[str, float],
    hpo_graph: HPOGraph,
    ancestor_cache: dict[str, frozenset[str]] | None,
    similarity_cache: dict[tuple[str, str], float] | None,
) -> float:
    """Compute symmetric term similarity while reusing prior calculations."""
    pair_key = tuple(sorted((first_hpo_id, second_hpo_id)))

    if similarity_cache is not None and pair_key in similarity_cache:
        return similarity_cache[pair_key]

    first_ancestors = _term_ancestors(first_hpo_id, hpo_graph, ancestor_cache)
    second_ancestors = _term_ancestors(second_hpo_id, hpo_graph, ancestor_cache)

    shared_terms = first_ancestors & second_ancestors
    score = max(
        (information_content.get(term, 0.0) for term in shared_terms),
        default=0.0,
    )

    if similarity_cache is not None:
        similarity_cache[pair_key] = score

    return score


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

    return _cached_pair_similarity(
        first_hpo_id,
        second_hpo_id,
        information_content,
        hpo_graph,
        ancestor_cache=None,
        similarity_cache=None,
    )


def best_match_average_similarity(
    patient_hpo_ids: tuple[str, ...],
    disease_hpo_ids: frozenset[str],
    information_content: dict[str, float],
    hpo_graph: HPOGraph,
    *,
    ancestor_cache: dict[str, frozenset[str]] | None = None,
    similarity_cache: dict[tuple[str, str], float] | None = None,
) -> float:
    """Average best-match similarities in both patient-to-disease directions.

    Optional caches can be shared across candidate diseases. The score formula
    is unchanged. This is an experimental ranking metric, not a probability.
    """
    patient_terms = tuple(dict.fromkeys(patient_hpo_ids))
    disease_terms = tuple(term for term in disease_hpo_ids if term in hpo_graph.ontology.terms)

    if not patient_terms or not disease_terms:
        return 0.0

    def pair_similarity(
        patient_term: str,
        disease_term: str,
    ) -> float:
        return _cached_pair_similarity(
            patient_term,
            disease_term,
            information_content,
            hpo_graph,
            ancestor_cache,
            similarity_cache,
        )

    patient_to_disease = [
        max(pair_similarity(patient_term, disease_term) for disease_term in disease_terms)
        for patient_term in patient_terms
    ]

    disease_to_patient = [
        max(pair_similarity(patient_term, disease_term) for patient_term in patient_terms)
        for disease_term in disease_terms
    ]

    total_information = sum(information_content.get(term, 0.0) for term in patient_terms)
    total_information += sum(information_content.get(term, 0.0) for term in disease_terms)

    if total_information <= 0:
        return 0.0

    matched_information = sum(patient_to_disease) + sum(disease_to_patient)
    return matched_information / total_information
