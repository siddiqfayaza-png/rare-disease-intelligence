from dataclasses import replace

import pytest

from rare_disease_ai.mechanism_evidence import (
    MechanismEvidence,
    validate_mechanism_evidence,
)


@pytest.fixture
def evidence() -> MechanismEvidence:
    return MechanismEvidence(
        gene_symbol="BRCA1",
        annotation_id="GO:0006281",
        annotation_label="DNA repair",
        annotation_type="biological_process",
        source_name="Example Ontology",
        source_release="2026-01",
        source_url="https://example.org/annotation/GO:0006281",
        evidence_code="EXP",
        evidence_label="experimental evidence",
        reference_ids=("PMID:12345678",),
        taxon_id="NCBITaxon:9606",
    )


def test_valid_evidence_preserves_provenance(
    evidence: MechanismEvidence,
) -> None:
    assert evidence.gene_symbol == "BRCA1"
    assert evidence.annotation_id == "GO:0006281"
    assert evidence.source_release == "2026-01"
    assert evidence.reference_ids == ("PMID:12345678",)
    assert evidence.evidence_code == "EXP"


@pytest.mark.parametrize(
    ("field_name", "value"),
    [
        ("gene_symbol", ""),
        ("annotation_id", " "),
        ("annotation_label", ""),
        ("source_name", ""),
        ("source_release", ""),
        ("source_url", ""),
    ],
)
def test_required_fields_cannot_be_blank(
    evidence: MechanismEvidence,
    field_name: str,
    value: str,
) -> None:
    with pytest.raises(ValueError):
        replace(evidence, **{field_name: value})


def test_unsupported_annotation_type_is_rejected(
    evidence: MechanismEvidence,
) -> None:
    with pytest.raises(ValueError, match="Unsupported annotation_type"):
        replace(evidence, annotation_type="clinical_diagnosis")


@pytest.mark.parametrize(
    "url",
    [
        "not-a-url",
        "ftp://example.org/annotation",
        "https:///missing-host",
    ],
)
def test_invalid_source_urls_are_rejected(
    evidence: MechanismEvidence,
    url: str,
) -> None:
    with pytest.raises(ValueError, match="source_url"):
        replace(evidence, source_url=url)


def test_reference_ids_must_be_a_tuple(
    evidence: MechanismEvidence,
) -> None:
    with pytest.raises(ValueError, match="tuple"):
        replace(evidence, reference_ids=["PMID:12345678"])


def test_blank_reference_ids_are_rejected(
    evidence: MechanismEvidence,
) -> None:
    with pytest.raises(ValueError, match="non-empty strings"):
        replace(evidence, reference_ids=("PMID:12345678", " "))


def test_optional_evidence_fields_reject_blank_strings(
    evidence: MechanismEvidence,
) -> None:
    with pytest.raises(ValueError, match="evidence_code"):
        replace(evidence, evidence_code=" ")


def test_duplicate_identity_is_reported(
    evidence: MechanismEvidence,
) -> None:
    duplicate = replace(evidence, annotation_label="Alternative label")

    issues = validate_mechanism_evidence((evidence, duplicate))

    assert len(issues) == 1
    assert "duplicate" in issues[0].lower()


def test_distinct_sources_are_not_duplicates(
    evidence: MechanismEvidence,
) -> None:
    second_source = replace(
        evidence,
        source_name="Another Ontology",
        source_url="https://other.example.org/annotation/GO:0006281",
    )

    assert validate_mechanism_evidence((evidence, second_source)) == ()


def test_invalid_record_type_is_reported() -> None:
    issues = validate_mechanism_evidence(("not evidence",))  # type: ignore[arg-type]

    assert len(issues) == 1
    assert "unsupported record type" in issues[0]
