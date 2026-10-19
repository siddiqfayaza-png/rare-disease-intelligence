from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology, HPOTerm
from rare_disease_ai.phenotype_semantic_similarity import (
    best_match_average_similarity,
    best_match_similarity,
)


def make_graph() -> HPOGraph:
    terms = {
        "HP:0000001": HPOTerm("HP:0000001", "Root phenotype"),
        "HP:0000002": HPOTerm(
            "HP:0000002",
            "Neurological phenotype",
            parents=("HP:0000001",),
        ),
        "HP:0000003": HPOTerm(
            "HP:0000003",
            "Seizure",
            parents=("HP:0000002",),
        ),
        "HP:0000004": HPOTerm(
            "HP:0000004",
            "Hearing impairment",
            parents=("HP:0000001",),
        ),
    }
    return HPOGraph(HPOOntology(terms))


def test_exact_match_uses_term_information_content():
    graph = make_graph()
    ic = {
        "HP:0000001": 0.1,
        "HP:0000002": 0.5,
        "HP:0000003": 2.0,
        "HP:0000004": 1.0,
    }

    score = best_match_similarity("HP:0000003", "HP:0000003", ic, graph)

    assert score == 2.0


def test_related_terms_use_their_most_informative_shared_ancestor():
    graph = make_graph()
    ic = {
        "HP:0000001": 0.1,
        "HP:0000002": 0.5,
        "HP:0000003": 2.0,
        "HP:0000004": 1.0,
    }

    score = best_match_similarity("HP:0000003", "HP:0000002", ic, graph)

    assert score == 0.5


def test_unrelated_branches_share_only_root():
    graph = make_graph()
    ic = {
        "HP:0000001": 0.1,
        "HP:0000002": 0.5,
        "HP:0000003": 2.0,
        "HP:0000004": 1.0,
    }

    score = best_match_similarity("HP:0000003", "HP:0000004", ic, graph)

    assert score == 0.1


def test_unknown_hpo_term_is_rejected():
    graph = make_graph()

    try:
        best_match_similarity(
            "HP:9999999",
            "HP:0000003",
            {"HP:0000001": 0.1},
            graph,
        )
    except ValueError as error:
        assert "must exist in the ontology" in str(error)
    else:
        raise AssertionError("Expected unknown HPO identifier to be rejected.")


def test_empty_disease_profile_returns_zero():
    graph = make_graph()

    score = best_match_average_similarity(
        ("HP:0000003",),
        frozenset(),
        {"HP:0000001": 0.1, "HP:0000003": 2.0},
        graph,
    )

    assert score == 0.0


def test_exact_profile_match_has_positive_similarity():
    graph = make_graph()
    ic = {
        "HP:0000001": 0.1,
        "HP:0000002": 0.5,
        "HP:0000003": 2.0,
        "HP:0000004": 1.0,
    }

    score = best_match_average_similarity(
        ("HP:0000003",),
        frozenset({"HP:0000003"}),
        ic,
        graph,
    )

    assert score == 1.0
