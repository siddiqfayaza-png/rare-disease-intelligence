"""Tests for HPO OBO parsing and term normalization."""

from pathlib import Path

import pytest

from rare_disease_ai.hpo_ontology import HPOOntology


@pytest.fixture
def ontology_file(tmp_path: Path) -> Path:
    content = """format-version: 1.2

[Term]
id: HP:0000001
name: All
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0001250
name: Seizures
synonym: "Convulsions" EXACT []
is_a: HP:0000707 ! Abnormality of the nervous system

[Term]
id: HP:0009998
name: Old phenotype
is_obsolete: true
replaced_by: HP:0001250

[Typedef]
id: part_of
name: part of
"""
    path = tmp_path / "mini-hpo.obo"
    path.write_text(content, encoding="utf-8")
    return path


def test_loads_hpo_terms(ontology_file: Path) -> None:
    ontology = HPOOntology.from_obo(ontology_file)
    term = ontology.get_term("HP:0001250")

    assert term is not None
    assert term.name == "Seizures"
    assert term.parents == ("HP:0000707",)
    assert term.synonyms == ("Convulsions",)


def test_normalizes_label_case_insensitively(ontology_file: Path) -> None:
    result = HPOOntology.from_obo(ontology_file).normalize("  seizures ")

    assert result.is_resolved
    assert result.term is not None
    assert result.term.hpo_id == "HP:0001250"


def test_normalizes_synonym(ontology_file: Path) -> None:
    result = HPOOntology.from_obo(ontology_file).normalize("Convulsions")

    assert result.is_resolved
    assert result.term is not None
    assert result.term.hpo_id == "HP:0001250"


def test_unknown_term_is_not_guessed(ontology_file: Path) -> None:
    result = HPOOntology.from_obo(ontology_file).normalize("Unrecognized phenotype")

    assert not result.is_resolved
    assert not result.is_ambiguous
    assert result.term is None


def test_obsolete_term_is_not_resolved_as_active(ontology_file: Path) -> None:
    result = HPOOntology.from_obo(ontology_file).normalize("HP:0009998")

    assert not result.is_resolved
    assert result.term is None
    assert result.matched_terms[0].is_obsolete
    assert result.matched_terms[0].replaced_by == ("HP:0001250",)


def test_empty_query_is_rejected(ontology_file: Path) -> None:
    ontology = HPOOntology.from_obo(ontology_file)

    with pytest.raises(ValueError):
        ontology.normalize("  ")


def test_missing_file_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        HPOOntology.from_obo(tmp_path / "missing.obo")
