"""Calculate phenotype information content from disease annotations."""

from __future__ import annotations

from math import log

from rare_disease_ai.disease_phenotypes import DiseasePhenotypeProfile
from rare_disease_ai.hpo_graph import HPOGraph


def calculate_information_content(
    disease_index: dict[str, DiseasePhenotypeProfile],
    hpo_graph: HPOGraph,
) -> dict[str, float]:
    """Calculate smoothed information content for HPO terms.

    Disease annotations are propagated to their HPO ancestors.
    More specific terms generally receive higher information content.

    These values measure annotation-based information, not disease
    probabilities. The calculation is dependent on the dataset.
    """
    total_diseases = len(disease_index)
    if total_diseases == 0:
        raise ValueError("The disease index must not be empty.")

    disease_counts: dict[str, int] = {}

    for profile in disease_index.values():
        # Count each HPO term at most once per disease, including ancestors.
        terms_for_disease: set[str] = set()

        for hpo_id in profile.positive_hpo_ids:
            if hpo_id not in hpo_graph.ontology.terms:
                continue

            terms_for_disease.add(hpo_id)
            terms_for_disease.update(hpo_graph.get_ancestors(hpo_id))

        for hpo_id in terms_for_disease:
            disease_counts[hpo_id] = disease_counts.get(hpo_id, 0) + 1

    information_content: dict[str, float] = {}

    for hpo_id in hpo_graph.ontology.terms:
        count = disease_counts.get(hpo_id, 0)

        # Add-one smoothing avoids infinite values for unobserved terms.
        frequency = (count + 1) / (total_diseases + 1)
        information_content[hpo_id] = -log(frequency)

    return information_content
