"""Reproducible audit of the Phenopacket Store benchmark cohort.

Audits case eligibility, data quality, ontology and disease index alignment,
and potential overlap signals for rare-disease candidate ranking benchmarks.
Ranking inputs are strictly separated from ground-truth diagnoses throughout.
No patient-level identifiers are exported in aggregate reports.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from rare_disease_ai.disease_phenotypes import build_disease_phenotype_index
from rare_disease_ai.hpo_ontology import HPOOntology
from rare_disease_ai.phenopacket_import import (
    ImportedPhenopacketCase,
    PhenopacketImportResult,
    import_phenopackets,
)


def compute_file_sha256(path: Path) -> str:
    """Compute the SHA-256 hexadecimal digest of a local file."""
    sha = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(1024 * 1024):
            sha.update(chunk)
    return sha.hexdigest()


def _load_optional_metadata(metadata_path: Path | None) -> dict[str, Any]:
    """Load metadata JSON file if present, returning an empty dict otherwise."""
    if metadata_path is not None and metadata_path.is_file():
        try:
            return json.loads(metadata_path.read_text(encoding="utf-8-sig"))
        except Exception:  # noqa: BLE001
            return {}
    return {}


def _compute_overlap_signals(cases: tuple[ImportedPhenopacketCase, ...]) -> dict[str, Any]:
    """Compute duplicate and overlap audit signals without exposing patient-level IDs."""
    patient_counts = Counter(c.patient_id for c in cases if c.patient_id)
    repeated_pids = {pid: count for pid, count in patient_counts.items() if count > 1}

    pub_pid_counts = Counter(
        (c.publication_id, c.patient_id)
        for c in cases
        if c.publication_id and c.patient_id
    )
    same_pub_same_pid = {pair: count for pair, count in pub_pid_counts.items() if count > 1}

    # Profile based strictly on sorted positive HPO ranking inputs
    pheno_counts = Counter(
        tuple(sorted(c.positive_hpo_ids))
        for c in cases
        if c.positive_hpo_ids
    )
    shared_pheno = {prof: count for prof, count in pheno_counts.items() if count > 1}

    pub_pheno_counts = Counter(
        (c.publication_id, tuple(sorted(c.positive_hpo_ids)))
        for c in cases
        if c.publication_id and c.positive_hpo_ids
    )
    same_pub_same_pheno = {pair: count for pair, count in pub_pheno_counts.items() if count > 1}

    disease_pheno_counts = Counter(
        (c.reference_disease_id, tuple(sorted(c.positive_hpo_ids)))
        for c in cases
        if c.reference_disease_id and c.positive_hpo_ids
    )
    same_disease_same_pheno = {
        pair: count for pair, count in disease_pheno_counts.items() if count > 1
    }

    return {
        "unique_repeated_patient_ids_count": len(repeated_pids),
        "total_cases_with_repeated_patient_ids": sum(repeated_pids.values()),
        "same_publication_same_patient_id_clusters": len(same_pub_same_pid),
        "identical_phenotype_profile_clusters_count": len(shared_pheno),
        "cases_sharing_identical_phenotype_profile": sum(shared_pheno.values()),
        "same_publication_identical_phenotype_clusters_count": len(same_pub_same_pheno),
        "cases_sharing_phenotype_within_same_publication": sum(same_pub_same_pheno.values()),
        "same_disease_identical_phenotype_clusters_count": len(same_disease_same_pheno),
        "cases_sharing_phenotype_with_same_disease": sum(same_disease_same_pheno.values()),
        "audit_interpretation": (
            "Identical phenotype profiles do not prove identical patients. In published "
            "cohorts, shared profiles often represent family members, minimal recorded "
            "feature sets, or standard syndromic presentations."
        ),
    }


def _compute_disease_index_alignment(
    cases: tuple[ImportedPhenopacketCase, ...],
    disease_index: dict[str, Any],
) -> dict[str, Any]:
    """Audit case reference diagnoses against the indexed disease annotations."""
    cases_in_index = sum(1 for c in cases if c.reference_disease_id in disease_index)
    cases_absent = sum(
        1 for c in cases if c.reference_disease_id and c.reference_disease_id not in disease_index
    )

    eligible_in_index = sum(
        1 for c in cases if c.is_eligible and c.reference_disease_id in disease_index
    )
    eligible_absent = sum(
        1 for c in cases
        if c.is_eligible and c.reference_disease_id and c.reference_disease_id not in disease_index
    )

    absent_disease_ids = sorted(
        {
            c.reference_disease_id
            for c in cases
            if c.reference_disease_id and c.reference_disease_id not in disease_index
        }
    )

    absent_breakdown: dict[str, dict[str, Any]] = {}
    for did in absent_disease_ids:
        matching_cases = [c for c in cases if c.reference_disease_id == did]
        cohorts = sorted({c.cohort for c in matching_cases if c.cohort})
        publications = sorted({c.publication_id for c in matching_cases if c.publication_id})
        disease_name = next(
            (c.disease_name for c in matching_cases if c.disease_name), did
        )
        absent_breakdown[did] = {
            "disease_name": disease_name,
            "case_count": len(matching_cases),
            "cohorts": cohorts,
            "publications": publications,
        }

    return {
        "total_diseases_in_index": len(disease_index),
        "cases_with_reference_disease_in_index": cases_in_index,
        "cases_with_reference_disease_absent_from_index": cases_absent,
        "eligible_cases_with_disease_in_index": eligible_in_index,
        "eligible_cases_with_disease_absent_from_index": eligible_absent,
        "absent_disease_ids": absent_disease_ids,
        "absent_disease_distribution": absent_breakdown,
    }


def _extract_exclusion_reasons(
    cases: tuple[ImportedPhenopacketCase, ...],
) -> dict[str, int]:
    """Summarize reasons why cases were marked ineligible for benchmark ranking."""
    reason_counts: Counter[str] = Counter()
    for case in cases:
        if not case.is_eligible:
            ineligibility_reasons = [r for r in case.quality_reasons if r.startswith("Ineligible:")]
            if ineligibility_reasons:
                for r in ineligibility_reasons:
                    reason_counts[r] += 1
            else:
                reason_counts["Ineligible: unspecified reason"] += 1
    return dict(sorted(reason_counts.items(), key=lambda x: (-x[1], x[0])))


def audit_phenopacket_benchmark(
    archive_path: str | Path,
    ontology: HPOOntology,
    disease_annotations: Any,
    *,
    archive_metadata_path: str | Path | None = None,
    ontology_metadata_path: str | Path | None = None,
    hpoa_metadata_path: str | Path | None = None,
) -> dict[str, Any]:
    """Audit the benchmark cohort against ontology and disease annotations.

    Runs import with strict diagnosis validation (require_diagnosis_match=True)
    and produces an auditable dictionary containing no patient-level identifiers.
    """
    archive_file = Path(archive_path)
    if not archive_file.is_file():
        raise FileNotFoundError(f"Archive file not found: {archive_file}")

    archive_sha256 = compute_file_sha256(archive_file)
    archive_size = archive_file.stat().st_size

    # Build disease index
    disease_index = build_disease_phenotype_index(disease_annotations)

    # Import archive with strict diagnosis validation
    import_result: PhenopacketImportResult = import_phenopackets(
        archive_file,
        ontology=ontology,
        require_diagnosis_match=True,
    )

    cases = import_result.cases

    # Load metadata descriptors if available
    ont_meta = _load_optional_metadata(
        Path(ontology_metadata_path) if ontology_metadata_path else None
    )
    hpoa_meta = _load_optional_metadata(
        Path(hpoa_metadata_path) if hpoa_metadata_path else None
    )

    disease_alignment = _compute_disease_index_alignment(cases, disease_index)
    overlap_signals = _compute_overlap_signals(cases)
    exclusion_reasons = _extract_exclusion_reasons(cases)

    # Phenotype data quality
    cases_with_defaulted = sum(1 for c in cases if c.applied_default_excluded_count > 0)
    total_defaulted_terms = sum(c.applied_default_excluded_count for c in cases)
    cases_with_zero_positives = sum(1 for c in cases if len(c.positive_hpo_ids) == 0)
    cases_with_conflicting = sum(
        1 for c in cases if set(c.positive_hpo_ids) & set(c.excluded_hpo_ids)
    )
    cases_missing_json_diagnoses = sum(1 for c in cases if not c.json_disease_ids)

    # Grouping counts
    cohort_counts = Counter(c.cohort for c in cases if c.cohort)
    pub_counts = Counter(c.publication_id for c in cases if c.publication_id)

    report: dict[str, Any] = {
        "report_metadata": {
            "generated_at": datetime.now(UTC).isoformat(),
            "benchmark_resource": "GA4GH Phenopacket Store",
            "benchmark_release": "0.1.27",
            "evaluation_readiness": (
                "Audit completed; research prototype validation only. Not for clinical use."
            ),
            "privacy_guarantee": (
                "No patient-level identifiers (names, MRNs, or raw subject IDs) "
                "are included in this aggregate audit report."
            ),
        },
        "source_metadata": {
            "archive": {
                "path": str(archive_file),
                "sha256": archive_sha256,
                "size_bytes": archive_size,
            },
            "ontology": {
                "release": ont_meta.get("release", "unknown"),
                "term_count": len(ontology.terms),
                "sha256": ont_meta.get("sha256", "unknown"),
            },
            "annotations": {
                "release": hpoa_meta.get("release", "unknown"),
                "annotation_count": len(disease_annotations),
                "sha256": hpoa_meta.get("sha256", "unknown"),
            },
        },
        "cohort_counts": {
            "summary_rows": import_result.total_summary_rows,
            "json_files": import_result.total_json_files,
            "total_cases_processed": import_result.total_cases,
            "eligible_cases": import_result.eligible_count,
            "ineligible_cases": import_result.ineligible_count,
            "unmatched_summary_rows": len(import_result.unmatched_summary_rows),
            "unmatched_json_files": len(import_result.unmatched_json_files),
            "ambiguous_matches": len(import_result.ambiguous_matches),
            "duplicate_summary_rows": len(import_result.duplicate_summary_rows),
            "diagnosis_mismatches": len(import_result.diagnosis_mismatches),
            "cases_with_missing_json_diagnoses": cases_missing_json_diagnoses,
            "cases_with_conflicting_phenotypes": cases_with_conflicting,
            "cases_with_zero_positive_phenotypes": cases_with_zero_positives,
        },
        "exclusion_reasons": exclusion_reasons,
        "disease_index_alignment": disease_alignment,
        "phenotype_quality": {
            "unique_malformed_hpo_ids": len(import_result.malformed_hpo_ids),
            "unique_unknown_hpo_ids": len(import_result.unknown_hpo_ids),
            "cases_with_defaulted_excluded_phenotypes": cases_with_defaulted,
            "total_phenotype_terms_defaulted_excluded": total_defaulted_terms,
        },
        "overlap_and_duplication_signals": overlap_signals,
        "grouping": {
            "total_cohorts": len(cohort_counts),
            "total_publications": len(pub_counts),
            "cases_per_cohort": dict(sorted(cohort_counts.items())),
            "cases_per_publication": dict(sorted(pub_counts.items())),
            "evaluation_split_guidance": (
                "All cases from the same publication must remain together in any future "
                "evaluation split to prevent cross-split leakage."
            ),
        },
    }

    return report


def format_console_summary(report: dict[str, Any]) -> str:
    """Format audit results as a concise, readable console summary."""
    cohort = report["cohort_counts"]
    alignment = report["disease_index_alignment"]
    quality = report["phenotype_quality"]
    overlap = report["overlap_and_duplication_signals"]
    grouping = report["grouping"]
    src = report["source_metadata"]

    total_cases = cohort["total_cases_processed"]
    el_pct = (cohort["eligible_cases"] / total_cases * 100) if total_cases else 0.0
    inel_pct = (cohort["ineligible_cases"] / total_cases * 100) if total_cases else 0.0

    lines = [
        "=" * 80,
        "PHENOPACKET STORE BENCHMARK AUDIT REPORT (Release 0.1.27)",
        "=" * 80,
        f"Generated:   {report['report_metadata']['generated_at']}",
        f"Archive:     {src['archive']['path']} (SHA256: {src['archive']['sha256'][:16]}...)",
        (
            f"Ontology:    Release {src['ontology']['release']} "
            f"({src['ontology']['term_count']:,} terms)"
        ),
        (
            f"Annotations: Release {src['annotations']['release']} "
            f"({src['annotations']['annotation_count']:,} records)"
        ),
        "",
        "1. COHORT ELIGIBILITY & INTEGRITY",
        f"   Total Summary Rows:             {cohort['summary_rows']:,}",
        f"   Total JSON Files:               {cohort['json_files']:,}",
        f"   Eligible Cases:                 {cohort['eligible_cases']:,} ({el_pct:.2f}%)",
        f"   Ineligible Cases:               {cohort['ineligible_cases']:,} ({inel_pct:.2f}%)",
        f"   Unmatched Summary Rows:         {cohort['unmatched_summary_rows']:,}",
        f"   Unmatched JSON Files:           {cohort['unmatched_json_files']:,}",
        f"   Duplicate Summary Rows:         {cohort['duplicate_summary_rows']:,}",
        f"   Ambiguous Path Matches:         {cohort['ambiguous_matches']:,}",
        f"   Diagnosis Mismatches:           {cohort['diagnosis_mismatches']:,}",
        f"   Missing JSON Diagnoses:         {cohort['cases_with_missing_json_diagnoses']:,}",
        f"   Conflicting Phenotypes:         {cohort['cases_with_conflicting_phenotypes']:,}",
        f"   Zero-Phenotype Cases:           {cohort['cases_with_zero_positive_phenotypes']:,}",
        "",
        "2. INELIGIBILITY REASONS",
    ]

    if report["exclusion_reasons"]:
        for reason, count in report["exclusion_reasons"].items():
            lines.append(f"   - {reason}: {count:,} case(s)")
    else:
        lines.append("   - None (100% of cases eligible)")

    in_idx = alignment["cases_with_reference_disease_in_index"]
    abs_idx = alignment["cases_with_reference_disease_absent_from_index"]

    lines.extend([
        "",
        "3. DISEASE INDEX ALIGNMENT",
        f"   Unique Diseases Indexed:        {alignment['total_diseases_in_index']:,}",
        f"   Cases with Disease in Index:    {in_idx:,}",
        f"   Cases Missing from Index:       {abs_idx:,}",
    ])

    if alignment["absent_disease_ids"]:
        lines.append("   Unindexed Reference Diagnoses:")
        for did, info in alignment["absent_disease_distribution"].items():
            cohort_str = ", ".join(info["cohorts"])
            pub_str = ", ".join(info["publications"])
            lines.append(
                f"     * {did} ({info['disease_name']}): {info['case_count']} case(s) "
                f"[cohort: {cohort_str}; pub: {pub_str}]"
            )

    defaulted_cases = quality["cases_with_defaulted_excluded_phenotypes"]
    defaulted_terms = quality["total_phenotype_terms_defaulted_excluded"]
    rep_pids = overlap["unique_repeated_patient_ids_count"]
    rep_cases = overlap["total_cases_with_repeated_patient_ids"]
    shared_clusters = overlap["identical_phenotype_profile_clusters_count"]
    shared_cases = overlap["cases_sharing_identical_phenotype_profile"]
    pub_clusters = overlap["same_publication_identical_phenotype_clusters_count"]
    pub_cases = overlap["cases_sharing_phenotype_within_same_publication"]
    same_pub_pids = overlap["same_publication_same_patient_id_clusters"]

    lines.extend([
        "",
        "4. PHENOTYPE DATA QUALITY",
        f"   Malformed HPO IDs:              {quality['unique_malformed_hpo_ids']:,}",
        f"   Unknown HPO IDs (not in OBO):   {quality['unique_unknown_hpo_ids']:,}",
        f"   Defaulted Excluded Cases:       {defaulted_cases:,} ({defaulted_terms:,} terms)",
        "",
        "5. OVERLAP & DUPLICATION SIGNALS",
        f"   Unique Repeated Patient IDs:    {rep_pids:,} ({rep_cases:,} cases)",
        f"   Same-Pub Same-Patient ID Dups:  {same_pub_pids:,}",
        f"   Identical Phenotype Profiles:   {shared_clusters:,} clusters ({shared_cases:,} cases)",
        f"   Same-Pub Identical Phenotypes:  {pub_clusters:,} clusters ({pub_cases:,} cases)",
        "   * Note: Identical phenotypes do not prove identical patients; do not deduplicate.",
        "",
        "6. COHORTS & PUBLICATIONS",
        f"   Total Unique Cohorts:           {grouping['total_cohorts']:,}",
        f"   Total Unique Publications:      {grouping['total_publications']:,}",
        "   * Guidance: All cases from same publication must remain together in evaluation splits.",
        "=" * 80,
    ])

    return "\n".join(lines)
