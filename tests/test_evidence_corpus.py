from rare_disease_ai.evidence_corpus import (
    gencc_assertion_to_document,
    hpo_annotation_to_document,
    make_evidence_id,
    mechanism_evidence_to_document,
)
from rare_disease_ai.gencc_import import GenCCAssertion
from rare_disease_ai.hpo_annotations import HPOAnnotation
from rare_disease_ai.mechanism_evidence import MechanismEvidence


def make_annotation(**overrides):
    values = {
        "database_id": "OMIM:123456",
        "disease_name": "Example disorder",
        "qualifier": "NOT",
        "hpo_id": "HP:0000001",
        "reference": "PMID:123456",
        "evidence": "IEA",
        "onset": "",
        "frequency": "",
        "sex": "",
        "modifier": "",
        "aspect": "P",
        "biocuration": "HPO:curator[2026-01-01]",
    }
    values.update(overrides)
    return HPOAnnotation(**values)


def test_hpoa_conversion_preserves_provenance_and_qualifier():
    document = hpo_annotation_to_document(
        make_annotation(),
        release="2026-09-02",
        source_url="https://example.org/phenotype.hpoa",
        hpo_label="All",
    )

    assert document.title == "Example disorder: All"
    assert document.source_release == "2026-09-02"
    assert document.source_url == "https://example.org/phenotype.hpoa"
    assert document.evidence_type == "hpo_disease_phenotype_annotation"
    assert document.entity_ids == ("OMIM:123456", "HP:0000001")
    assert document.reference_ids == ("PMID:123456",)
    assert "qualifier=NOT" in document.text


def test_hpoa_conversion_uses_hpo_id_when_label_is_missing():
    document = hpo_annotation_to_document(
        make_annotation(),
        release="2026-09-02",
        source_url="https://example.org/phenotype.hpoa",
        hpo_label=None,
    )

    assert document.title == "Example disorder: HP:0000001"


def test_evidence_id_is_deterministic():
    first = make_evidence_id("hpoa", "OMIM:123456", "HP:0000001")
    second = make_evidence_id("hpoa", "OMIM:123456", "HP:0000001")

    assert first == second


def test_evidence_id_rejects_empty_identifiers():
    import pytest

    with pytest.raises(ValueError):
        make_evidence_id("hpoa", "OMIM:123456", "")




def make_gencc_assertion(**overrides):
    values = {
        "submission_id": "SGC-1",
        "version_number": 1,
        "gene_curie": "HGNC:1",
        "gene_symbol": "GENE1",
        "disease_curie": "MONDO:0000001",
        "disease_title": "Example disease",
        "disease_original_curie": "OMIM:123456",
        "classification_curie": "GENCC:CLASS",
        "classification_title": "Disputed Evidence",
        "moi_curie": "HP:0000006",
        "moi_title": "Autosomal dominant",
        "submitter_curie": "GENCC:000101",
        "submitter_title": "Example submitter",
        "submitted_as_date": "2025-01-01",
        "public_report_url": "https://example.org/report",
        "pmids": ("PMID:12345",),
        "assertion_criteria_url": "https://example.org/criteria",
    }
    values.update(overrides)
    return GenCCAssertion(**values)


def test_gencc_conversion_preserves_assertion_provenance():
    document = gencc_assertion_to_document(
        make_gencc_assertion(),
        release="2026-10-20",
        source_url="https://thegencc.org/download/action/submissions-export-csv?format=new",
    )

    assert document.source_release == "2026-10-20"
    assert document.evidence_type == "gene_disease_assertion"
    assert "classification=Disputed Evidence" in document.text
    assert "mode_of_inheritance=Autosomal dominant" in document.text
    assert "original_disease_id=OMIM:123456" in document.text
    assert "canonical_disease_id=MONDO:0000001" in document.text
    assert "submitter=Example submitter" in document.text
    assert document.reference_ids == ("PMID:12345",)


def test_gencc_conversion_is_deterministic():
    assertion = make_gencc_assertion()
    kwargs = {
        "release": "2026-10-20",
        "source_url": "https://example.org/gencc.csv",
    }

    first = gencc_assertion_to_document(assertion, **kwargs)
    second = gencc_assertion_to_document(assertion, **kwargs)

    assert first.evidence_id == second.evidence_id
    assert "not an independently verified clinical conclusion" in first.text


def test_gencc_conversion_keeps_submission_versions_distinct():
    kwargs = {
        "release": "2026-10-20",
        "source_url": "https://example.org/gencc.csv",
    }
    first = gencc_assertion_to_document(
        make_gencc_assertion(version_number=1), **kwargs
    )
    second = gencc_assertion_to_document(
        make_gencc_assertion(version_number=2), **kwargs
    )

    assert first.evidence_id != second.evidence_id
def make_mechanism_evidence(**overrides):
    values = {
        "gene_symbol": "BRCA1",
        "annotation_id": "GO:0006281",
        "annotation_label": "DNA repair",
        "annotation_type": "biological_process",
        "source_name": "Gene Ontology",
        "source_release": "2026-08-05",
        "source_url": "https://release.geneontology.org/2026-08-05/annotations/gaf/HUMAN-uniprot.gaf.gz",
        "evidence_code": "IDA",
        "evidence_label": "Inferred from Direct Assay",
        "reference_ids": ("PMID:12345",),
        "qualifier": "NOT|involved_in",
        "taxon_id": "taxon:9606",
    }
    values.update(overrides)
    return MechanismEvidence(**values)


def test_go_conversion_preserves_provenance_and_qualifier():
    document = mechanism_evidence_to_document(make_mechanism_evidence())

    assert document.title == "BRCA1: DNA repair (GO:0006281)"
    assert document.source_name == "Gene Ontology"
    assert document.source_release == "2026-08-05"
    assert document.source_url.startswith("https://release.geneontology.org/")
    assert document.evidence_type == "go_mechanism_annotation"
    assert "qualifier=NOT|involved_in" in document.text
    assert "evidence_code=IDA" in document.text
    assert document.entity_ids == ("BRCA1", "GO:0006281", "taxon:9606")
    assert document.reference_ids == ("PMID:12345",)


def test_go_conversion_is_deterministic():
    evidence = make_mechanism_evidence()

    first = mechanism_evidence_to_document(evidence)
    second = mechanism_evidence_to_document(evidence)

    assert first.evidence_id == second.evidence_id


def test_go_conversion_distinguishes_evidence_records():
    first = mechanism_evidence_to_document(make_mechanism_evidence())
    second = mechanism_evidence_to_document(
        make_mechanism_evidence(evidence_code="IEA")
    )

    assert first.evidence_id != second.evidence_id


def test_go_conversion_preserves_multiple_references():
    document = mechanism_evidence_to_document(
        make_mechanism_evidence(
            reference_ids=("PMID:12345", "PMID:67890")
        )
    )

    assert document.reference_ids == ("PMID:12345", "PMID:67890")
    assert "reference_ids=PMID:12345;PMID:67890" in document.text
