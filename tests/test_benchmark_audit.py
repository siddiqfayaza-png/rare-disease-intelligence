"""Unit tests for the Phenopacket benchmark audit module using synthetic fixtures."""

from __future__ import annotations

import csv
import io
import json
import zipfile
from pathlib import Path
from typing import Any

import pytest

from rare_disease_ai.benchmark_audit import (
    audit_phenopacket_benchmark,
    compute_file_sha256,
    format_console_summary,
)
from rare_disease_ai.hpo_annotations import HPOAnnotation


def _create_synthetic_archive(
    tmp_path: Path,
    summary_rows: list[dict[str, str]],
    json_cases: dict[str, dict[str, Any]],
    archive_name: str = "phenopackets.zip",
) -> Path:
    """Helper to build a small synthetic zip archive."""
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
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=summary_columns, delimiter="\t")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow({col: row.get(col, "") for col in summary_columns})
        zf.writestr("0.1.0/phenopacket_store.summary.tsv", buffer.getvalue().encode("utf-8"))

        for path_str, data in json_cases.items():
            zf.writestr(path_str, json.dumps(data).encode("utf-8"))

    return zip_path


def _create_mock_ontology(known_terms: set[str]) -> Any:
    """Create a minimal mock ontology containing declared term identifiers."""
    terms = {term_id: type("Term", (), {"hpo_id": term_id})() for term_id in known_terms}
    return type("MockOntology", (), {"terms": terms, "get_term": lambda s, tid: terms.get(tid)})()


def _create_mock_annotations(disease_to_hpos: dict[str, list[str]]) -> tuple[HPOAnnotation, ...]:
    """Create synthetic disease phenotype annotations for index construction."""
    annotations = []
    for disease_id, hpo_ids in disease_to_hpos.items():
        for hpo_id in hpo_ids:
            annotations.append(
                HPOAnnotation(
                    database_id=disease_id,
                    disease_name=f"Disease {disease_id}",
                    qualifier="",
                    hpo_id=hpo_id,
                    reference="PMID:100",
                    evidence="IEA",
                    onset="",
                    frequency="",
                    sex="",
                    modifier="",
                    aspect="P",
                    biocuration="",
                )
            )
    return tuple(annotations)


def test_compute_file_sha256(tmp_path: Path) -> None:
    """Verify SHA-256 calculation for a local file."""
    test_file = tmp_path / "test.txt"
    test_file.write_text("benchmark audit payload\n", encoding="utf-8")
    digest = compute_file_sha256(test_file)
    assert isinstance(digest, str)
    assert len(digest) == 64


def test_audit_synthetic_benchmark_happy_path(tmp_path: Path) -> None:
    """Verify audit execution on a compliant synthetic benchmark archive."""
    summary_rows = [
        {
            "patient_id": "P_001",
            "cohort": "COHORT_ALPHA",
            "disease_id": "OMIM:100100",
            "disease": "Disorder Alpha",
            "PMID": "PMID:11111111",
            "filename": "COHORT_ALPHA/phenopackets/Case_1",
        },
        {
            "patient_id": "P_002",
            "cohort": "COHORT_BETA",
            "disease_id": "OMIM:200200",
            "disease": "Disorder Beta",
            "PMID": "PMID:22222222",
            "filename": "COHORT_BETA/phenopackets/Case_2",
        },
    ]

    json_cases = {
        "0.1.0/COHORT_ALPHA/Case_1.json": {
            "id": "Case_1",
            "subject": {"id": "P_001"},
            "diseases": [{"term": {"id": "OMIM:100100"}}],
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
        },
        "0.1.0/COHORT_BETA/Case_2.json": {
            "id": "Case_2",
            "subject": {"id": "P_002"},
            "diseases": [{"term": {"id": "OMIM:200200"}}],
            "phenotypicFeatures": [{"type": {"id": "HP:0000256"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    ontology = _create_mock_ontology({"HP:0001250", "HP:0000256"})
    annotations = _create_mock_annotations(
        {"OMIM:100100": ["HP:0001250"], "OMIM:200200": ["HP:0000256"]}
    )

    report = audit_phenopacket_benchmark(archive_path, ontology, annotations)

    counts = report["cohort_counts"]
    assert counts["total_cases_processed"] == 2
    assert counts["eligible_cases"] == 2
    assert counts["ineligible_cases"] == 0
    assert counts["summary_rows"] == 2
    assert counts["json_files"] == 2

    alignment = report["disease_index_alignment"]
    assert alignment["cases_with_reference_disease_in_index"] == 2
    assert alignment["cases_with_reference_disease_absent_from_index"] == 0
    assert alignment["absent_disease_ids"] == []

    grouping = report["grouping"]
    assert grouping["total_cohorts"] == 2
    assert grouping["total_publications"] == 2
    assert grouping["cases_per_cohort"] == {"COHORT_ALPHA": 1, "COHORT_BETA": 1}
    assert grouping["cases_per_publication"] == {"PMID:11111111": 1, "PMID:22222222": 1}

    # Verify absence of patient-level identifiers in JSON output
    serialized = json.dumps(report)
    assert "P_001" not in serialized
    assert "P_002" not in serialized


def test_audit_detects_unindexed_reference_disease(tmp_path: Path) -> None:
    """Verify detection and breakdown of reference disease IDs missing from index."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "COHORT_A",
            "disease_id": "OMIM:100100",
            "filename": "COHORT_A/phenopackets/Case_1",
        },
        {
            "patient_id": "P2",
            "cohort": "COHORT_A",
            "disease_id": "OMIM:999999",  # Unindexed disease
            "disease": "Novel Syndrome",
            "PMID": "PMID:33333333",
            "filename": "COHORT_A/phenopackets/Case_2",
        },
    ]

    json_cases = {
        "0.1.0/COHORT_A/Case_1.json": {
            "id": "Case_1",
            "diseases": [{"term": {"id": "OMIM:100100"}}],
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
        },
        "0.1.0/COHORT_A/Case_2.json": {
            "id": "Case_2",
            "diseases": [{"term": {"id": "OMIM:999999"}}],
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    ontology = _create_mock_ontology({"HP:0001250"})
    annotations = _create_mock_annotations({"OMIM:100100": ["HP:0001250"]})

    report = audit_phenopacket_benchmark(archive_path, ontology, annotations)

    alignment = report["disease_index_alignment"]
    assert alignment["cases_with_reference_disease_in_index"] == 1
    assert alignment["cases_with_reference_disease_absent_from_index"] == 1
    assert alignment["absent_disease_ids"] == ["OMIM:999999"]

    breakdown = alignment["absent_disease_distribution"]["OMIM:999999"]
    assert breakdown["disease_name"] == "Novel Syndrome"
    assert breakdown["case_count"] == 1
    assert breakdown["cohorts"] == ["COHORT_A"]
    assert breakdown["publications"] == ["PMID:33333333"]


def test_audit_detects_ineligible_cases_and_exclusion_reasons(tmp_path: Path) -> None:
    """Verify recording of ineligible cases due to missing phenotypes or conflicts."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/ZeroPheno",
        },
        {
            "patient_id": "P2",
            "cohort": "C1",
            "disease_id": "OMIM:100100",
            "filename": "C1/phenopackets/ConflictingPheno",
        },
    ]

    json_cases = {
        "0.1.0/C1/ZeroPheno.json": {
            "id": "ZeroPheno",
            "diseases": [{"term": {"id": "OMIM:100100"}}],
            "phenotypicFeatures": [
                {"type": {"id": "HP:0001250"}, "excluded": True},
            ],
        },
        "0.1.0/C1/ConflictingPheno.json": {
            "id": "ConflictingPheno",
            "diseases": [{"term": {"id": "OMIM:100100"}}],
            "phenotypicFeatures": [
                {"type": {"id": "HP:0001250"}, "excluded": False},
                {"type": {"id": "HP:0001250"}, "excluded": True},
            ],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    ontology = _create_mock_ontology({"HP:0001250"})
    annotations = _create_mock_annotations({"OMIM:100100": ["HP:0001250"]})

    report = audit_phenopacket_benchmark(archive_path, ontology, annotations)

    counts = report["cohort_counts"]
    assert counts["eligible_cases"] == 0
    assert counts["ineligible_cases"] == 2
    assert counts["cases_with_zero_positive_phenotypes"] == 1
    assert counts["cases_with_conflicting_phenotypes"] == 1

    reasons = report["exclusion_reasons"]
    assert len(reasons) >= 2


def test_audit_overlap_and_duplication_signals(tmp_path: Path) -> None:
    """Verify detection of shared phenotype profiles and repeated patient IDs."""
    summary_rows = [
        {
            "patient_id": "Patient 1",
            "cohort": "C1",
            "disease_id": "OMIM:100",
            "PMID": "PMID:100",
            "filename": "C1/phenopackets/Case_1",
        },
        {
            "patient_id": "Patient 1",  # Repeated patient ID
            "cohort": "C2",
            "disease_id": "OMIM:100",
            "PMID": "PMID:200",
            "filename": "C2/phenopackets/Case_2",
        },
        {
            "patient_id": "Patient 2",
            "cohort": "C1",
            "disease_id": "OMIM:100",
            "PMID": "PMID:100",
            "filename": "C1/phenopackets/Case_3",
        },
    ]

    json_cases = {
        # Case 1 & Case 3 share identical positive phenotype profiles within same pub
        "0.1.0/C1/Case_1.json": {
            "id": "Case_1",
            "diseases": [{"term": {"id": "OMIM:100"}}],
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
        },
        "0.1.0/C2/Case_2.json": {
            "id": "Case_2",
            "diseases": [{"term": {"id": "OMIM:100"}}],
            "phenotypicFeatures": [{"type": {"id": "HP:0000256"}}],
        },
        "0.1.0/C1/Case_3.json": {
            "id": "Case_3",
            "diseases": [{"term": {"id": "OMIM:100"}}],
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    ontology = _create_mock_ontology({"HP:0001250", "HP:0000256"})
    annotations = _create_mock_annotations({"OMIM:100": ["HP:0001250"]})

    report = audit_phenopacket_benchmark(archive_path, ontology, annotations)

    overlap = report["overlap_and_duplication_signals"]
    assert overlap["unique_repeated_patient_ids_count"] == 1
    assert overlap["total_cases_with_repeated_patient_ids"] == 2
    assert overlap["identical_phenotype_profile_clusters_count"] == 1
    assert overlap["cases_sharing_identical_phenotype_profile"] == 2
    assert overlap["same_publication_identical_phenotype_clusters_count"] == 1
    assert overlap["cases_sharing_phenotype_within_same_publication"] == 2

    # Verify no raw patient names leaked into report
    serialized = json.dumps(report)
    assert "Patient 1" not in serialized
    assert "Patient 2" not in serialized


def test_format_console_summary_output(tmp_path: Path) -> None:
    """Verify that console summary string contains all required section headers."""
    summary_rows = [
        {
            "patient_id": "P1",
            "cohort": "C1",
            "disease_id": "OMIM:100",
            "filename": "C1/phenopackets/Case_1",
        },
    ]

    json_cases = {
        "0.1.0/C1/Case_1.json": {
            "id": "Case_1",
            "diseases": [{"term": {"id": "OMIM:100"}}],
            "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
        },
    }

    archive_path = _create_synthetic_archive(tmp_path, summary_rows, json_cases)
    ontology = _create_mock_ontology({"HP:0001250"})
    annotations = _create_mock_annotations({"OMIM:100": ["HP:0001250"]})

    report = audit_phenopacket_benchmark(archive_path, ontology, annotations)
    summary_text = format_console_summary(report)

    assert "COHORT ELIGIBILITY & INTEGRITY" in summary_text
    assert "INELIGIBILITY REASONS" in summary_text
    assert "DISEASE INDEX ALIGNMENT" in summary_text
    assert "PHENOTYPE DATA QUALITY" in summary_text
    assert "OVERLAP & DUPLICATION SIGNALS" in summary_text
    assert "COHORTS & PUBLICATIONS" in summary_text


def test_audit_missing_archive_file_raises(tmp_path: Path) -> None:
    """Verify FileNotFoundError when archive file is missing."""
    ontology = _create_mock_ontology({"HP:0001250"})
    with pytest.raises(FileNotFoundError, match="Archive file not found"):
        audit_phenopacket_benchmark(tmp_path / "missing.zip", ontology, ())
