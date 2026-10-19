import pytest

from rare_disease_ai.disease_phenotypes import DiseasePhenotypeProfile
from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology, HPOTerm
from rare_disease_ai.phenotype_information_content import (
    calculate_information_content,
)


def make_graph() -> HPOGraph:
    terms = {
        "HP:0000001": HPOTerm(
            hpo_id="HP:0000001",
            name="Root",
        ),
        "HP:0000002": HPOTerm(
            hpo_id="HP:0000002",
            name="Neurological phenotype",
            parents=("HP:0000001",),
        ),
        "HP:0000003": HPOTerm(
            hpo_id="HP:0000003",
            name="Rare phenotype",
            parents=("HP:0000002",),
        ),
    }
    ontology = HPOOntology(terms)
    return HPOGraph(ontology)


def make_profile(
    disease_id: str,
    positive_hpo_ids: frozenset[str],
) -> DiseasePhenotypeProfile:
    return DiseasePhenotypeProfile(
        database_id=disease_id,
        disease_name=disease_id,
        positive_hpo_ids=positive_hpo_ids,
        negative_hpo_ids=frozenset(),
        annotations=(),
    )


def test_specific_term_has_higher_information_content() -> None:
    graph = make_graph()
    disease_index = {
        "D1": make_profile("D1", frozenset({"HP:0000003"})),
        "D2": make_profile("D2", frozenset({"HP:0000002"})),
        "D3": make_profile("D3", frozenset({"HP:0000002"})),
    }

    ic = calculate_information_content(disease_index, graph)

    assert ic["HP:0000001"] == pytest.approx(0.0, abs=0.02)
    assert ic["HP:0000003"] > ic["HP:0000002"]


def test_empty_disease_index_is_rejected() -> None:
    with pytest.raises(ValueError, match="must not be empty"):
        calculate_information_content({}, make_graph())
