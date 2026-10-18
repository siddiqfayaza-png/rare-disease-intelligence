
from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology, HPOTerm
from rare_disease_ai.ontology_validation import validate_ontology


def make_ontology(
    terms: dict[str, HPOTerm],
) -> HPOOntology:
    return HPOOntology(terms)


def test_valid_ontology_passes_integrity_checks() -> None:
    ontology = make_ontology(
        {
            "HP:0000001": HPOTerm(
                hpo_id="HP:0000001",
                name="Root",
            ),
            "HP:0000002": HPOTerm(
                hpo_id="HP:0000002",
                name="Child",
                parents=("HP:0000001",),
            ),
        }
    )

    report = validate_ontology(ontology, HPOGraph(ontology))

    assert report.is_valid
    assert report.term_count == 2
    assert not report.cycle_detected


def test_missing_parent_reference_is_detected() -> None:
    ontology = make_ontology(
        {
            "HP:0000002": HPOTerm(
                hpo_id="HP:0000002",
                name="Child",
                parents=("HP:9999999",),
            ),
        }
    )

    report = validate_ontology(ontology, HPOGraph(ontology))

    assert not report.is_valid
    assert report.missing_parent_references == (("HP:0000002", "HP:9999999"),)


def test_missing_replacement_reference_is_detected() -> None:
    ontology = make_ontology(
        {
            "HP:0000002": HPOTerm(
                hpo_id="HP:0000002",
                name="Obsolete term",
                is_obsolete=True,
                replaced_by=("HP:9999999",),
            ),
        }
    )

    report = validate_ontology(ontology, HPOGraph(ontology))

    assert not report.is_valid
    assert report.missing_replacement_references == (("HP:0000002", "HP:9999999"),)


def test_cycle_is_detected() -> None:
    ontology = make_ontology(
        {
            "HP:0000001": HPOTerm(
                hpo_id="HP:0000001",
                name="First",
                parents=("HP:0000002",),
            ),
            "HP:0000002": HPOTerm(
                hpo_id="HP:0000002",
                name="Second",
                parents=("HP:0000001",),
            ),
        }
    )

    report = validate_ontology(ontology, HPOGraph(ontology))

    assert not report.is_valid
    assert report.cycle_detected


def test_unknown_parent_child_link_is_not_reported_as_inconsistent() -> None:
    ontology = make_ontology(
        {
            "HP:0000001": HPOTerm(
                hpo_id="HP:0000001",
                name="Root",
            ),
            "HP:0000002": HPOTerm(
                hpo_id="HP:0000002",
                name="Child",
                parents=("HP:0000001",),
            ),
        }
    )

    report = validate_ontology(ontology, HPOGraph(ontology))

    assert report.inconsistent_parent_child_links == ()
