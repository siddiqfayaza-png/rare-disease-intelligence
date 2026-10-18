"""Tests for patient case validation."""

import pytest

from rare_disease_ai.case_validation import PatientCase


def test_valid_case() -> None:
    case = PatientCase(
        case_id=" CASE-001 ",
        phenotype_terms=(" Seizures ", "Developmental delay"),
    )

    assert case.case_id == "CASE-001"
    assert case.phenotype_terms == ("Seizures", "Developmental delay")


def test_empty_case_id_is_rejected() -> None:
    with pytest.raises(ValueError):
        PatientCase(case_id=" ", phenotype_terms=("Seizures",))


def test_empty_phenotype_list_is_rejected() -> None:
    with pytest.raises(ValueError):
        PatientCase(case_id="CASE-001", phenotype_terms=())


def test_empty_phenotype_term_is_rejected() -> None:
    with pytest.raises(ValueError):
        PatientCase(case_id="CASE-001", phenotype_terms=("Seizures", " "))


def test_list_is_rejected_instead_of_tuple() -> None:
    with pytest.raises(TypeError):
        PatientCase(
            case_id="CASE-001",
            phenotype_terms=["Seizures"],  # type: ignore[arg-type]
        )
