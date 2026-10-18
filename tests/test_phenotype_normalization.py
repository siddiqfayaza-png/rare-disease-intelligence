from rare_disease_ai.case_validation import PatientCase
from rare_disease_ai.hpo_ontology import HPOOntology, HPOTerm
from rare_disease_ai.phenotype_normalization import normalize_patient_case


def make_ontology() -> HPOOntology:
    terms = {
        "HP:0001250": HPOTerm(
            hpo_id="HP:0001250",
            name="Seizure",
            synonyms=("Epileptic seizure", "Seizures"),
        ),
        "HP:0000002": HPOTerm(
            hpo_id="HP:0000002",
            name="Ambiguous feature",
        ),
        "HP:0000003": HPOTerm(
            hpo_id="HP:0000003",
            name="Another feature",
            synonyms=("Ambiguous feature",),
        ),
        "HP:0000004": HPOTerm(
            hpo_id="HP:0000004",
            name="Old phenotype",
            is_obsolete=True,
            replaced_by=("HP:0001250",),
        ),
    }
    return HPOOntology(terms)


def test_normalizes_preferred_label_and_synonym() -> None:
    case = PatientCase(
        case_id="case-001",
        phenotype_terms=("Seizure", "Seizures"),
    )

    result = normalize_patient_case(case, make_ontology())

    assert result.case_id == "case-001"
    assert [record.status for record in result.records] == [
        "resolved",
        "resolved",
    ]
    assert result.normalized_hpo_ids == ("HP:0001250",)
    assert result.records[1].matched_synonym == "Seizures"


def test_preserves_unresolved_phenotype() -> None:
    case = PatientCase(
        case_id="case-002",
        phenotype_terms=("Unknown feature",),
    )

    result = normalize_patient_case(case, make_ontology())

    assert result.records[0].status == "unresolved"
    assert result.records[0].original_text == "Unknown feature"
    assert result.normalized_hpo_ids == ()


def test_marks_ambiguous_phenotype() -> None:
    case = PatientCase(
        case_id="case-003",
        phenotype_terms=("Ambiguous feature",),
    )

    result = normalize_patient_case(case, make_ontology())

    record = result.records[0]
    assert record.status == "ambiguous"
    assert set(record.candidate_hpo_ids) == {
        "HP:0000002",
        "HP:0000003",
    }


def test_marks_obsolete_phenotype() -> None:
    case = PatientCase(
        case_id="case-004",
        phenotype_terms=("HP:0000004",),
    )

    result = normalize_patient_case(case, make_ontology())

    record = result.records[0]
    assert record.status == "obsolete"
    assert record.candidate_hpo_ids == ("HP:0001250",)
