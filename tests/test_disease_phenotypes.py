
from rare_disease_ai.disease_phenotypes import (
    build_disease_phenotype_index,
)
from rare_disease_ai.hpo_annotations import HPOAnnotation


def make_annotation(
    database_id="OMIM:123456",
    disease_name="Example disease",
    hpo_id="HP:0001250",
    qualifier="",
):
    return HPOAnnotation(
        database_id=database_id,
        disease_name=disease_name,
        qualifier=qualifier,
        hpo_id=hpo_id,
        reference="PMID:123",
        evidence="IEA",
        onset="",
        frequency="",
        sex="",
        modifier="",
        aspect="P",
        biocuration="",
    )


def test_index_groups_annotations_by_disease():
    annotations = (
        make_annotation(hpo_id="HP:0001250"),
        make_annotation(hpo_id="HP:0001263"),
        make_annotation(
            database_id="OMIM:654321",
            disease_name="Another disease",
            hpo_id="HP:0004322",
        ),
    )

    index = build_disease_phenotype_index(annotations)

    assert len(index) == 2
    assert index["OMIM:123456"].positive_hpo_ids == frozenset({"HP:0001250", "HP:0001263"})
    assert index["OMIM:654321"].positive_hpo_ids == frozenset({"HP:0004322"})


def test_index_separates_negative_annotations():
    annotations = (
        make_annotation(hpo_id="HP:0001250"),
        make_annotation(hpo_id="HP:0001263", qualifier="NOT"),
    )

    profile = build_disease_phenotype_index(annotations)["OMIM:123456"]

    assert profile.positive_hpo_ids == frozenset({"HP:0001250"})
    assert profile.negative_hpo_ids == frozenset({"HP:0001263"})
    assert len(profile.annotations) == 2


def test_index_preserves_annotation_evidence_and_references():
    annotation = make_annotation()
    index = build_disease_phenotype_index((annotation,))

    profile = index["OMIM:123456"]

    assert profile.annotations == (annotation,)
    assert profile.annotations[0].reference == "PMID:123"
    assert profile.annotations[0].evidence == "IEA"


def test_index_preserves_disease_name_variants():
    annotations = (
        make_annotation(disease_name="Parietal foramina"),
        make_annotation(disease_name="Parietal foramina 1"),
    )

    profile = build_disease_phenotype_index(annotations)["OMIM:123456"]

    assert profile.disease_name_variants == (
        "Parietal foramina",
        "Parietal foramina 1",
    )
    assert len(profile.annotations) == 2


def test_index_uses_identifier_when_name_is_empty():
    annotation = make_annotation(disease_name="")
    profile = build_disease_phenotype_index((annotation,))["OMIM:123456"]

    assert profile.disease_name == "OMIM:123456"
    assert profile.disease_name_variants == ()


def test_index_handles_empty_input():
    assert build_disease_phenotype_index(()) == {}
