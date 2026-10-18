import pytest

from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology, HPOTerm


@pytest.fixture
def graph() -> HPOGraph:
    terms = {
        "HP:0000001": HPOTerm(
            hpo_id="HP:0000001",
            name="Root",
        ),
        "HP:0000002": HPOTerm(
            hpo_id="HP:0000002",
            name="Parent",
            parents=("HP:0000001",),
        ),
        "HP:0000003": HPOTerm(
            hpo_id="HP:0000003",
            name="Child",
            parents=("HP:0000002",),
        ),
        "HP:0000004": HPOTerm(
            hpo_id="HP:0000004",
            name="Other parent",
            parents=("HP:0000001",),
        ),
        "HP:0000005": HPOTerm(
            hpo_id="HP:0000005",
            name="Multiple inheritance child",
            parents=("HP:0000002", "HP:0000004"),
        ),
    }
    return HPOGraph(HPOOntology(terms))


def test_get_direct_parents(graph: HPOGraph) -> None:
    assert graph.get_parents("HP:0000003") == ("HP:0000002",)


def test_get_direct_children(graph: HPOGraph) -> None:
    assert graph.get_children("HP:0000002") == (
        "HP:0000003",
        "HP:0000005",
    )


def test_get_all_ancestors(graph: HPOGraph) -> None:
    assert graph.get_ancestors("HP:0000003") == (
        "HP:0000001",
        "HP:0000002",
    )


def test_ancestors_support_multiple_parents(graph: HPOGraph) -> None:
    assert graph.get_ancestors("HP:0000005") == (
        "HP:0000001",
        "HP:0000002",
        "HP:0000004",
    )


def test_unknown_term_raises_key_error(graph: HPOGraph) -> None:
    with pytest.raises(KeyError, match="Unknown HPO identifier"):
        graph.get_ancestors("HP:9999999")
