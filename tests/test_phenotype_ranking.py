from rare_disease_ai.disease_phenotypes import DiseasePhenotypeProfile
from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology, HPOTerm
from rare_disease_ai.phenotype_ranking import rank_diseases


def make_profile(
    database_id: str,
    name: str,
    positive: set[str],
    negative: set[str] | None = None,
) -> DiseasePhenotypeProfile:
    return DiseasePhenotypeProfile(
        database_id=database_id,
        disease_name=name,
        positive_hpo_ids=frozenset(positive),
        negative_hpo_ids=frozenset(negative or set()),
        annotations=(),
    )


def make_graph() -> HPOGraph:
    terms = {
        "HP:0000001": HPOTerm(
            hpo_id="HP:0000001",
            name="Root phenotype",
        ),
        "HP:0000002": HPOTerm(
            hpo_id="HP:0000002",
            name="Neurological phenotype",
            parents=("HP:0000001",),
        ),
        "HP:0000003": HPOTerm(
            hpo_id="HP:0000003",
            name="Seizure",
            parents=("HP:0000002",),
        ),
        "HP:0000004": HPOTerm(
            hpo_id="HP:0000004",
            name="Hearing impairment",
            parents=("HP:0000001",),
        ),
    }
    return HPOGraph(HPOOntology(terms))


def test_exact_match_ranks_above_unmatched_disease():
    graph = make_graph()
    index = {
        "OMIM:1": make_profile("OMIM:1", "Seizure disease", {"HP:0000003"}),
        "OMIM:2": make_profile("OMIM:2", "Hearing disease", {"HP:0000004"}),
    }

    results = rank_diseases(("HP:0000003",), index, graph)

    assert results[0].database_id == "OMIM:1"
    assert results[0].score == 1.0
    assert results[0].exact_matches == ("HP:0000003",)


def test_related_match_receives_partial_score():
    graph = make_graph()
    index = {
        "OMIM:1": make_profile("OMIM:1", "Neurological disease", {"HP:0000002"}),
    }

    results = rank_diseases(("HP:0000003",), index, graph)

    assert results[0].score == 0.5
    assert results[0].related_matches == ("HP:0000003",)


def test_unmatched_phenotype_receives_zero():
    graph = make_graph()
    index = {
        "OMIM:1": make_profile("OMIM:1", "Hearing disease", {"HP:0000004"}),
    }

    results = rank_diseases(("HP:0000003",), index, graph)

    assert results[0].score == 0.0
    assert results[0].unmatched_patient_hpo_ids == ("HP:0000003",)


def test_negative_matches_are_reported_separately():
    graph = make_graph()
    index = {
        "OMIM:1": make_profile(
            "OMIM:1",
            "Seizure disease",
            {"HP:0000003"},
            {"HP:0000003"},
        ),
    }

    results = rank_diseases(("HP:0000003",), index, graph)

    assert results[0].score == 1.0
    assert results[0].negative_matches == ("HP:0000003",)


def test_duplicate_patient_terms_are_counted_once():
    graph = make_graph()
    index = {
        "OMIM:1": make_profile("OMIM:1", "Seizure disease", {"HP:0000003"}),
    }

    results = rank_diseases(("HP:0000003", "HP:0000003"), index, graph)

    assert results[0].score == 1.0
    assert results[0].exact_matches == ("HP:0000003",)


def test_unknown_patient_hpo_id_is_rejected():
    graph = make_graph()

    try:
        rank_diseases(("HP:9999999",), {}, graph)
    except ValueError as error:
        assert "Unknown patient HPO identifiers" in str(error)
    else:
        raise AssertionError("Expected unknown HPO identifier to be rejected.")


def test_empty_patient_phenotypes_are_rejected():
    graph = make_graph()

    try:
        rank_diseases((), {}, graph)
    except ValueError as error:
        assert "At least one patient HPO identifier" in str(error)
    else:
        raise AssertionError("Expected empty input to be rejected.")


def test_top_k_limits_number_of_results():
    graph = make_graph()
    index = {
        "OMIM:1": make_profile("OMIM:1", "Seizure", {"HP:0000003"}),
        "OMIM:2": make_profile("OMIM:2", "Hearing", {"HP:0000004"}),
    }

    results = rank_diseases(("HP:0000003",), index, graph, top_k=1)

    assert len(results) == 1
    assert results[0].database_id == "OMIM:1"


def test_semantic_score_is_none_by_default():
    graph = make_graph()
    index = {
        "OMIM:1": make_profile("OMIM:1", "Seizure disease", {"HP:0000003"}),
    }

    results = rank_diseases(("HP:0000003",), index, graph)

    assert results[0].semantic_score is None


def test_semantic_similarity_calculates_a_score():
    graph = make_graph()
    index = {
        "OMIM:1": make_profile(
            "OMIM:1",
            "Seizure disease",
            {"HP:0000003"},
        ),
        "OMIM:2": make_profile(
            "OMIM:2",
            "Hearing disease",
            {"HP:0000004"},
        ),
    }

    results = rank_diseases(
        ("HP:0000003",),
        index,
        graph,
        use_semantic_similarity=True,
    )

    seizure_result = next(result for result in results if result.database_id == "OMIM:1")

    assert seizure_result.semantic_score is not None
    assert seizure_result.semantic_score == 1.0


def test_semantic_similarity_sorts_candidates_by_score():
    """Rank candidates by semantic score when semantic ranking is enabled."""
    graph = make_graph()
    index = {
        "OMIM:1": make_profile(
            "OMIM:1",
            "Seizure disease",
            {"HP:0000003"},
        ),
        "OMIM:2": make_profile(
            "OMIM:2",
            "Neurological disease",
            {"HP:0000002"},
        ),
        "OMIM:3": make_profile(
            "OMIM:3",
            "Hearing disease",
            {"HP:0000004"},
        ),
    }

    results = rank_diseases(
        ("HP:0000003",),
        index,
        graph,
        use_semantic_similarity=True,
    )

    semantic_scores = [result.semantic_score for result in results]

    assert all(score is not None for score in semantic_scores)
    assert semantic_scores == sorted(semantic_scores, reverse=True)
    assert results[0].database_id == "OMIM:1"
