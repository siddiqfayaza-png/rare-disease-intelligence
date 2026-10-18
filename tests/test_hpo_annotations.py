import pytest

from rare_disease_ai.hpo_annotations import (
    HPOAnnotation,
    read_hpo_annotations,
    validate_hpo_annotations,
)


def test_read_hpo_annotations_preserves_all_fields(tmp_path):
    path = tmp_path / "phenotype.hpoa"
    path.write_text(
        "# HPO annotation metadata\n"
        "database_id\tdisease_name\tqualifier\thpo_id\treference\t"
        "evidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration\n"
        "OMIM:123456\tExample disease\tNOT\tHP:0001250\tPMID:123\t"
        "IEA\t\t1/2\t\t\tP\tCurator[123]\n",
        encoding="utf-8",
    )

    annotations = read_hpo_annotations(path)

    assert len(annotations) == 1
    annotation = annotations[0]
    assert annotation.database_id == "OMIM:123456"
    assert annotation.disease_name == "Example disease"
    assert annotation.qualifier == "NOT"
    assert annotation.hpo_id == "HP:0001250"
    assert annotation.reference == "PMID:123"
    assert annotation.evidence == "IEA"
    assert annotation.frequency == "1/2"
    assert annotation.aspect == "P"
    assert annotation.onset == ""


def test_read_hpo_annotations_rejects_wrong_header(tmp_path):
    path = tmp_path / "invalid.hpoa"
    path.write_text("disease\thpo_id\nExample\tHP:0001250\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Unexpected HPO annotation columns"):
        read_hpo_annotations(path)


def test_validation_reports_unknown_hpo_ids():
    known_term = "HP:0001250"
    ontology = type("FakeOntology", (), {"terms": {known_term: object()}})()

    annotations = (
        HPOAnnotation(
            database_id="OMIM:123456",
            disease_name="Example disease",
            qualifier="",
            hpo_id=known_term,
            reference="PMID:123",
            evidence="IEA",
            onset="",
            frequency="",
            sex="",
            modifier="",
            aspect="P",
            biocuration="",
        ),
        HPOAnnotation(
            database_id="OMIM:654321",
            disease_name="Another disease",
            qualifier="",
            hpo_id="HP:9999999",
            reference="PMID:456",
            evidence="IEA",
            onset="",
            frequency="",
            sex="",
            modifier="",
            aspect="P",
            biocuration="",
        ),
    )

    report = validate_hpo_annotations(annotations, ontology)

    assert report.annotation_count == 2
    assert report.unique_disease_count == 2
    assert report.unknown_hpo_ids == ("HP:9999999",)
    assert report.is_valid is False


def test_validation_accepts_known_hpo_ids():
    ontology = type("FakeOntology", (), {"terms": {"HP:0001250": object()}})()
    annotation = HPOAnnotation(
        database_id="OMIM:123456",
        disease_name="Example disease",
        qualifier="",
        hpo_id="HP:0001250",
        reference="PMID:123",
        evidence="IEA",
        onset="",
        frequency="",
        sex="",
        modifier="",
        aspect="P",
        biocuration="",
    )

    report = validate_hpo_annotations((annotation,), ontology)

    assert report.is_valid is True
    assert report.annotation_count == 1
    assert report.unique_disease_count == 1
