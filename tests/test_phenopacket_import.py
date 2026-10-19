"""Unit tests for the Phenopacket Store benchmark importer."""

from __future__ import annotations

import csv
import io
import json
import zipfile
from pathlib import Path
from typing import Any

import pytest

from rare_disease_ai.phenopacket_import import (
    PhenopacketImportResult,
    import_phenopackets,
)


def _create_synthetic_archive(
    tmp_path: Path,
    summary_rows: list[dict[str, str]],
    json_cases: dict[str, dict[str, Any]],
    archive_name: str = "phenopackets.zip",
    summary_path_in_zip: str = "0.1.0/phenopacket_store.summary.tsv",
) -> Path:
    """Helper to build a synthetic zip archive in a temporary directory."""
    zip_path = tmp_path / archive_name

    summary_columns = [
        "patient_id",
        "cohort",
        "disease_id",
        "disease",
        "gene",
        "allele_1",
        "allele_2",
        "PMID",
        "filename",
    ]

    with zipfile.ZipFile(zip_path, mode="w") as zf:
        # Write summary TSV
        summary_buffer = io.StringIO()
        writer = csv.DictWriter(summary_buffer, fieldnames=summary_columns, delimiter="\t")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow({col: row.get(col, "") for col in summary_columns})

        zf.writestr(summary_path_in_zip, summary_buffer.getvalue().encode("utf-8"))

        # Write each JSON file
        for rel_path, data in json_cases.items():
            zf.writestr(rel_path, json.dumps(data).encode("utf-8"))

    return zip_path


def test_import_valid_cases_unique_path_matching(tmp_path: Path) -> None:
    """Test standard import with unique path matching and field retention."""
    summary_rows = [
        {
            "patient_id": "Indiv_1",
            "cohort": "COHORT_A",
            "disease_id": "OMIM:100100",
            "disease": "Syndrome A",
            "PMID": "PMID:12345678",
            "filename": "COHORT_A/phenopackets/Case_1",
        },
        {
            "patient_id": "Indiv_2",
            "cohort": "COHORT_B",
            "disease_id": "OMIM:200200",
            "disease": "Syndrome B",
            "PMID": "PMID:87654321",
            "filename": "COHORT_B/phenopackets/Case_2",
        },
    ]

    json_cases = {
        "0.1.0/COHORT_A/Case_1.json": {
            "id": "Case_1",
            "subject": {"id": "Indiv_1"},
            "diseases": [{"term": {"id": "OMIM:100100", "label": "Syndrome A"}}],
            "phenotypicFeatures": [
                {"type": {"id": "HP:0001250", "label": "Seizures"}},
                {"type": {"id": "HP:0000256", "label": "Macrocephaly"}, "excluded": False},
                {"type": {"id": "HP:0000316", "label": "Hypertelorism"}, "excluded": True},
            ],
        },
        "0.1.0/COHORT_B/Case_2.json": {
            "id": "Case_2",
            "subject": {"id": "Indiv_2"},
            "diseases": [{"term": {"id": "OMIM:200200", "label": "Syndrome B"}}],
            "phenotypicFeatures": [
                {"type": {"id": "HP:0001250", "label": "Seizures"}, "excluded": False},
            ],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path)

    assert isinstance(result, PhenopacketImportResult)
    assert result.total_cases == 2
    assert result.eligible_count == 2
    assert result.ineligible_count == 0
    assert len(result.unmatched_summary_rows) == 0
    assert len(result.unmatched_json_files) == 0

    case1 = result.cases[0]
    assert case1.case_id == "0.1.0/COHORT_A/Case_1.json"
    assert case1.archive_path == "0.1.0/COHORT_A/Case_1.json"
    assert case1.cohort == "COHORT_A"
    assert case1.patient_id == "Indiv_1"
    assert case1.publication_id == "PMID:12345678"
    assert case1.reference_disease_id == "OMIM:100100"
    assert case1.disease_name == "Syndrome A"
    assert case1.positive_hpo_ids == ("HP:0001250", "HP:0000256")
    assert case1.excluded_hpo_ids == ("HP:0000316",)
    assert case1.target_disease_ids == frozenset({"OMIM:100100"})
    assert case1.is_eligible is True

    patient_case = case1.to_patient_case()
    assert patient_case.case_id == case1.case_id
    assert patient_case.phenotype_terms == ("HP:0001250", "HP:0000256")


def test_unmatched_and_ambiguous_summary_filenames(tmp_path: Path) -> None:
    """Test unmatched summary filenames and ambiguous multiple JSON matches."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/Case_Missing",
        },
        {
            "patient_id": "P2",
            "cohort": "C2",
            "disease_id": "OMIM:200200",
            "filename": "C2/phenopackets/Case_Ambiguous",
        },
    ]

    json_cases = {
        # Two files with the same stem 'Case_Ambiguous'
        "0.1.0/dir1/Case_Ambiguous.json": {
            "id": "Case_Ambiguous",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
            "diseases": [{"term": {"id": "OMIM:200200"}}],
        },
        "0.1.0/dir2/Case_Ambiguous.json": {
            "id": "Case_Ambiguous",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
            "diseases": [{"term": {"id": "OMIM:200200"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path)

    assert result.total_cases == 2
    assert result.eligible_count == 0
    assert result.ineligible_count == 2

    assert "C1/phenopackets/Case_Missing" in result.unmatched_summary_rows
    assert len(result.ambiguous_matches) == 1
    assert result.ambiguous_matches[0][0] == "C2/phenopackets/Case_Ambiguous"

    assert result.cases[0].is_eligible is False
    assert any("No matching JSON file" in r for r in result.cases[0].quality_reasons)
    assert result.cases[1].is_eligible is False
    assert any("Ambiguous match" in r for r in result.cases[1].quality_reasons)


def test_unmatched_json_files_reported(tmp_path: Path) -> None:
    """Test that JSON files not referenced in the summary TSV are tracked."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/Case_1",
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_1.json": {
            "id": "Case_1",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
            "diseases": [{"term": {"id": "OMIM:100100"}}],
        },
        "0.1.0/C1/Unreferenced_Case.json": {
            "id": "Unreferenced_Case",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
            "diseases": [{"term": {"id": "OMIM:100100"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path)

    assert result.eligible_count == 1
    assert result.unmatched_json_files == ("0.1.0/C1/Unreferenced_Case.json",)


def test_excluded_phenotypes_strictly_separated_from_positive_ranking_inputs(
    tmp_path: Path,
) -> None:
    """Test that excluded: true phenotypes are isolated from ranking inputs."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/Case_1",
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_1.json": {
            "id": "Case_1",
            "phenotypicFeatures": [
                {"type": {"id": "HP:0001250"}, "excluded": False},
                {"type": {"id": "HP:0000256"}, "excluded": True},
                {"type": {"id": "HP:0001263"}, "excluded": True},
            ],
            "diseases": [{"term": {"id": "OMIM:100100"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path)

    case = result.cases[0]
    assert case.positive_hpo_ids == ("HP:0001250",)
    assert case.excluded_hpo_ids == ("HP:0000256", "HP:0001263")
    assert case.ranking_input.positive_hpo_ids == ("HP:0001250",)
    assert case.ranking_input.excluded_hpo_ids == ("HP:0000256", "HP:0001263")
    assert case.to_patient_case().phenotype_terms == ("HP:0001250",)


def test_missing_excluded_field_defaults_to_false_and_is_recorded(tmp_path: Path) -> None:
    """Test that missing excluded field applies schema default and records audit note."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/Case_1",
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_1.json": {
            "id": "Case_1",
            "phenotypicFeatures": [
                {"type": {"id": "HP:0001250"}},  # excluded omitted
                {"type": {"id": "HP:0000256"}},  # excluded omitted
                {"type": {"id": "HP:0000316"}, "excluded": False},
            ],
            "diseases": [{"term": {"id": "OMIM:100100"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path)

    case = result.cases[0]
    assert case.positive_hpo_ids == ("HP:0001250", "HP:0000256", "HP:0000316")
    assert case.applied_default_excluded_count == 2
    assert any(
        "Applied schema default (excluded=False) to 2 phenotype(s)" in reason
        for reason in case.quality_reasons
    )
    assert any("not independently confirmed" in reason for reason in case.quality_reasons)


def test_rejects_malformed_and_unknown_hpo_ids_with_reasons(tmp_path: Path) -> None:
    """Test that malformed HPO IDs and IDs absent from ontology are rejected."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/Case_1",
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_1.json": {
            "id": "Case_1",
            "phenotypicFeatures": [
                {"type": {"id": "HP:0001250"}},  # Valid & in ontology
                {"type": {"id": "HP:123"}},  # Malformed (not 7 digits)
                {"type": {"id": "NOT_AN_HPO_ID"}},  # Malformed
                {"type": {"id": "HP:9999999"}},  # Valid format, but absent from ontology
            ],
            "diseases": [{"term": {"id": "OMIM:100100"}}],
        },
    }

    mock_ontology = type("FakeOntology", (), {"terms": {"HP:0001250": object()}})()

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path, ontology=mock_ontology)

    case = result.cases[0]
    assert case.positive_hpo_ids == ("HP:0001250",)
    assert case.is_eligible is True  # Still has one valid positive phenotype

    assert "HP:123" in result.malformed_hpo_ids
    assert "NOT_AN_HPO_ID" in result.malformed_hpo_ids
    assert "HP:9999999" in result.unknown_hpo_ids

    assert any("Malformed HPO ID rejected: 'HP:123'" in r for r in case.quality_reasons)
    assert any("Malformed HPO ID rejected: 'NOT_AN_HPO_ID'" in r for r in case.quality_reasons)
    assert any(
        "HPO ID absent from ontology rejected: 'HP:9999999'" in r for r in case.quality_reasons
    )


def test_cases_with_no_eligible_positive_phenotypes_are_ineligible(tmp_path: Path) -> None:
    """Test cases where no valid positive phenotypes exist are marked ineligible."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100",
            "filename": "C1/phenopackets/AllExcluded",
        },
        {
            "patient_id": "P2",
            "cohort": "C1",
            "disease_id": "OMIM:100",
            "filename": "C1/phenopackets/AllMalformed",
        },
        {
            "patient_id": "P3",
            "cohort": "C1",
            "disease_id": "OMIM:100",
            "filename": "C1/phenopackets/EmptyFeatures",
        },
    ]

    json_cases = {
        "0.1.0/C1/AllExcluded.json": {
            "id": "AllExcluded",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}, "excluded": True}],
            "diseases": [{"term": {"id": "OMIM:100"}}],
        },
        "0.1.0/C1/AllMalformed.json": {
            "id": "AllMalformed",
            "phenotypicFeatures": [{"type": {"id": "BAD_ID"}}],
            "diseases": [{"term": {"id": "OMIM:100"}}],
        },
        "0.1.0/C1/EmptyFeatures.json": {
            "id": "EmptyFeatures",
            "phenotypicFeatures": [],
            "diseases": [{"term": {"id": "OMIM:100"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path)

    assert result.eligible_count == 0
    assert result.ineligible_count == 3
    for case in result.cases:
        assert case.is_eligible is False
        assert any("no valid positive HPO phenotypes" in r for r in case.quality_reasons)
        with pytest.raises(ValueError, match="not benchmark-eligible"):
            case.to_patient_case()


def test_missing_summary_reference_diagnosis_is_ineligible(tmp_path: Path) -> None:
    """Test that a summary row without disease_id is marked ineligible."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "",  # Missing diagnosis
            "filename": "C1/phenopackets/Case_1",
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_1.json": {
            "id": "Case_1",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
            "diseases": [{"term": {"id": "OMIM:100100"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path)

    case = result.cases[0]
    assert case.is_eligible is False
    assert case.reference_disease_id == ""
    assert case.target_disease_ids == frozenset()
    assert any("missing reference disease ID" in r for r in case.quality_reasons)


def test_summary_json_diagnosis_mismatch_recorded(tmp_path: Path) -> None:
    """Test that summary vs JSON diagnosis mismatch is cross-checked and recorded."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",  # Summary disease
            "filename": "C1/phenopackets/Case_1",
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_1.json": {
            "id": "Case_1",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
            "diseases": [{"term": {"id": "OMIM:999999"}}],  # Different from summary
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)

    # By default, summary disease_id is authoritative reference; mismatch is recorded for audit
    result_default = import_phenopackets(archive_path, require_diagnosis_match=False)
    case_def = result_default.cases[0]
    assert case_def.reference_disease_id == "OMIM:100100"
    assert case_def.json_disease_ids == ("OMIM:999999",)
    assert case_def.diagnosis_mismatch is True
    assert case_def.is_eligible is True
    assert len(result_default.diagnosis_mismatches) == 1
    assert any("Diagnosis mismatch" in r for r in case_def.quality_reasons)

    # When require_diagnosis_match=True, mismatch marks case ineligible
    result_strict = import_phenopackets(archive_path, require_diagnosis_match=True)
    assert result_strict.cases[0].is_eligible is False
    assert any("mismatches JSON diseases" in r for r in result_strict.cases[0].quality_reasons)


def test_json_diseases_empty_recorded(tmp_path: Path) -> None:
    """Test that empty JSON diseases list is recorded in quality reasons."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/Case_1",
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_1.json": {
            "id": "Case_1",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
            "diseases": [],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path)

    case = result.cases[0]
    assert case.json_disease_ids == ()
    assert any("JSON diseases field is empty or missing" in r for r in case.quality_reasons)


def test_duplicate_summary_rows_flagged_and_ineligible(tmp_path: Path) -> None:
    """Test that duplicate summary rows targeting the same case are flagged as ineligible."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/Case_1",
        },
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/Case_1",  # Duplicate row
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_1.json": {
            "id": "Case_1",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
            "diseases": [{"term": {"id": "OMIM:100100"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path)

    assert result.total_cases == 2
    assert result.eligible_count == 1
    assert result.ineligible_count == 1

    assert result.cases[0].is_eligible is True
    assert result.cases[1].is_eligible is False

    assert "0.1.0/C1/Case_1.json" in result.duplicate_summary_rows
    assert any("Duplicate summary row" in r for r in result.cases[1].quality_reasons)


def test_repeated_patient_ids_reported_as_audit_observation_not_deduplicated(
    tmp_path: Path,
) -> None:
    """Test that repeated patient IDs across different case files are not deduplicated."""
    summary_rows = [
        {
            "patient_id": "Individual 1",
            "cohort": "Cohort_A",
            "disease_id": "OMIM:100100",
            "filename": "Cohort_A/phenopackets/Family1_Indiv1",
        },
        {
            "patient_id": "Individual 1",  # Same patient ID, different case file
            "cohort": "Cohort_B",
            "disease_id": "OMIM:200200",
            "filename": "Cohort_B/phenopackets/Family2_Indiv1",
        },
    ]

    json_cases = {
        "0.1.0/Cohort_A/Family1_Indiv1.json": {
            "id": "Family1_Indiv1",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
            "diseases": [{"term": {"id": "OMIM:100100"}}],
        },
        "0.1.0/Cohort_B/Family2_Indiv1.json": {
            "id": "Family2_Indiv1",
            "phenotypicFeatures": [{"type": {"id": "HP:0000256"}}],
            "diseases": [{"term": {"id": "OMIM:200200"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path)

    # Both cases must remain distinct and eligible!
    assert result.total_cases == 2
    assert result.eligible_count == 2
    assert result.cases[0].is_eligible is True
    assert result.cases[1].is_eligible is True
    assert result.cases[0].case_id != result.cases[1].case_id

    # Repeated patient ID is reported as audit observation
    assert ("Individual 1", 2) in result.repeated_patient_ids
    assert result.repeated_patient_id_counts["Individual 1"] == 2


def test_ranking_inputs_strictly_separated_from_ground_truth(tmp_path: Path) -> None:
    """Test that interpretations, genes, variants, and diagnoses are never in ranking inputs."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "disease": "Neurodevelopmental disorder",
            "gene": "ESAM",
            "allele_1": "c.35T>A",
            "filename": "C1/phenopackets/Case_1",
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_1.json": {
            "id": "Case_1",
            "subject": {"id": "P1"},
            "phenotypicFeatures": [
                {"type": {"id": "HP:0001250"}},
            ],
            "diseases": [{"term": {"id": "OMIM:100100", "label": "Neurodevelopmental disorder"}}],
            "interpretations": [{"id": "Interp_1", "diagnosis": {"genomicInterpretations": []}}],
            "genes": [{"id": "ESAM"}],
            "variants": [{"id": "var1"}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path)

    case = result.cases[0]
    # Only phenotypic features should be in positive_hpo_ids
    assert case.positive_hpo_ids == ("HP:0001250",)
    assert case.ranking_input.positive_hpo_ids == ("HP:0001250",)

    # Ground truth is isolated
    assert case.ground_truth.reference_disease_id == "OMIM:100100"
    assert case.ground_truth.disease_name == "Neurodevelopmental disorder"
    assert case.ground_truth.json_disease_ids == ("OMIM:100100",)


def test_missing_summary_file_raises_error(tmp_path: Path) -> None:
    """Test that archive missing the summary TSV raises FileNotFoundError."""
    zip_path = tmp_path / "empty_archive.zip"
    with zipfile.ZipFile(zip_path, mode="w") as zf:
        zf.writestr("test.txt", "hello")

    with pytest.raises(FileNotFoundError, match="Summary file .* not found in archive"):
        import_phenopackets(zip_path)


def test_corrupt_json_file_handled_gracefully(tmp_path: Path) -> None:
    """Test that a corrupt JSON file records reason and marks case ineligible without crashing."""
    zip_path = tmp_path / "corrupt.zip"
    with zipfile.ZipFile(zip_path, mode="w") as zf:
        zf.writestr(
            "0.1.0/phenopacket_store.summary.tsv",
            "patient_id\tcohort\tdisease_id\tdisease\tgene\tallele_1\tallele_2\tPMID\tfilename\n"
            "P1\tC1\tOMIM:100100\tDisease\t\t\t\t\tC1/phenopackets/CorruptCase\n",
        )
        zf.writestr("0.1.0/C1/CorruptCase.json", b"{broken json content")

    result = import_phenopackets(zip_path)
    assert result.total_cases == 1
    assert result.eligible_count == 0
    assert result.cases[0].is_eligible is False
    assert any("Failed to parse JSON file" in r for r in result.cases[0].quality_reasons)


def test_missing_archive_file_raises_file_not_found(tmp_path: Path) -> None:
    """Test that a non-existent archive path raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError, match="Phenopacket archive not found"):
        import_phenopackets(tmp_path / "non_existent.zip")


def test_strict_mode_empty_or_missing_json_diseases_is_ineligible(tmp_path: Path) -> None:
    """Test that strict mode marks cases with empty or missing JSON diseases as ineligible."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/Case_EmptyDiseases",
        },
        {
            "patient_id": "P2",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/Case_MissingDiseasesKey",
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_EmptyDiseases.json": {
            "id": "Case_EmptyDiseases",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
            "diseases": [],
        },
        "0.1.0/C1/Case_MissingDiseasesKey.json": {
            "id": "Case_MissingDiseasesKey",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path, require_diagnosis_match=True)

    assert result.total_cases == 2
    assert result.eligible_count == 0
    assert result.ineligible_count == 2

    for case in result.cases:
        assert case.is_eligible is False
        assert case.diagnosis_mismatch is False
        assert any(
            "JSON diseases list is empty or missing in strict mode" in r
            for r in case.quality_reasons
        )


def test_conflicting_positive_and_excluded_assertions_marked_ineligible(
    tmp_path: Path,
) -> None:
    """Test that conflicting positive and excluded assertions make a case ineligible."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/Case_Conflict",
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_Conflict.json": {
            "id": "Case_Conflict",
            "phenotypicFeatures": [
                {"type": {"id": "HP:0001250"}, "excluded": False},
                {"type": {"id": "HP:0001250"}, "excluded": True},
                {"type": {"id": "HP:0000256"}, "excluded": False},
            ],
            "diseases": [{"term": {"id": "OMIM:100100"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path)

    case = result.cases[0]
    # Original assertions are preserved for auditing
    assert "HP:0001250" in case.positive_hpo_ids
    assert "HP:0000256" in case.positive_hpo_ids
    assert "HP:0001250" in case.excluded_hpo_ids

    # Case is marked ineligible due to conflict
    assert case.is_eligible is False
    assert result.eligible_count == 0
    assert any("Conflicting positive and excluded assertions" in r for r in case.quality_reasons)
    assert any(
        "Ineligible: conflicting positive and excluded phenotype assertions" in r
        for r in case.quality_reasons
    )


def test_multiple_candidate_summary_files_without_exact_match_raises_value_error(
    tmp_path: Path,
) -> None:
    """Test that multiple candidate summary TSVs without an exact match raise ValueError."""
    zip_path = tmp_path / "ambiguous_summaries.zip"
    with zipfile.ZipFile(zip_path, mode="w") as zf:
        zf.writestr(
            "v1/phenopacket_store.summary.tsv",
            "patient_id\tcohort\tdisease_id\tdisease\tgene\tallele_1\tallele_2\tPMID\tfilename\n",
        )
        zf.writestr(
            "v2/phenopacket_store.summary.tsv",
            "patient_id\tcohort\tdisease_id\tdisease\tgene\tallele_1\tallele_2\tPMID\tfilename\n",
        )

    with pytest.raises(
        ValueError, match="Multiple candidate summary files found without an exact-path match"
    ):
        import_phenopackets(zip_path)

    # When exact root match is present alongside a subdirectory candidate, it resolves
    zip_exact_path = tmp_path / "exact_summary.zip"
    with zipfile.ZipFile(zip_exact_path, mode="w") as zf:
        zf.writestr(
            "phenopacket_store.summary.tsv",
            "patient_id\tcohort\tdisease_id\tdisease\tgene\tallele_1\tallele_2\tPMID\tfilename\n",
        )
        zf.writestr(
            "v2/phenopacket_store.summary.tsv",
            "patient_id\tcohort\tdisease_id\tdisease\tgene\tallele_1\tallele_2\tPMID\tfilename\n",
        )

    result = import_phenopackets(zip_exact_path)
    assert result.summary_path == "phenopacket_store.summary.tsv"


def test_permissive_mode_preserves_diagnosis_mismatch_and_empty_diseases_eligibility(
    tmp_path: Path,
) -> None:
    """Test that permissive mode preserves eligibility on mismatch or empty diseases."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/Case_Mismatch",
        },
        {
            "patient_id": "P2",
            "cohort": "C1",
            "disease_id": "OMIM:200200",
            "filename": "C1/phenopackets/Case_Empty",
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_Mismatch.json": {
            "id": "Case_Mismatch",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
            "diseases": [{"term": {"id": "OMIM:999999"}}],
        },
        "0.1.0/C1/Case_Empty.json": {
            "id": "Case_Empty",
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
            "diseases": [],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    result = import_phenopackets(archive_path, require_diagnosis_match=False)

    assert result.total_cases == 2
    assert result.eligible_count == 2
    assert result.cases[0].is_eligible is True
    assert result.cases[0].diagnosis_mismatch is True
    assert result.cases[1].is_eligible is True
    assert result.cases[1].diagnosis_mismatch is False
    assert any(
        "JSON diseases field is empty or missing" in r for r in result.cases[1].quality_reasons
    )
