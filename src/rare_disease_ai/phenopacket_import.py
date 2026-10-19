"""Benchmark importer for Phenopacket Store case archives.

Extracts patient phenotypic profiles and reference disease diagnoses from
Phenopacket Store releases for disease-ranking benchmark evaluation.
Ranking inputs are strictly separated from ground-truth diagnoses.
"""

from __future__ import annotations

import csv
import io
import json
import re
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import TYPE_CHECKING, Any

from rare_disease_ai.case_validation import PatientCase

if TYPE_CHECKING:
    from rare_disease_ai.hpo_ontology import HPOOntology

HPO_ID_PATTERN = re.compile(r"^HP:\d{7}$")


@dataclass(frozen=True)
class PhenopacketRankingInput:
    """Phenotypic inputs strictly permitted for candidate ranking."""

    positive_hpo_ids: tuple[str, ...]
    excluded_hpo_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class PhenopacketGroundTruth:
    """Reference diagnostic targets for benchmark evaluation."""

    reference_disease_id: str
    disease_name: str = ""
    json_disease_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class ImportedPhenopacketCase:
    """An imported Phenopacket benchmark case with separated ranking and ground-truth data."""

    case_id: str
    archive_path: str
    cohort: str
    patient_id: str
    publication_id: str
    positive_hpo_ids: tuple[str, ...]
    excluded_hpo_ids: tuple[str, ...]
    reference_disease_id: str
    disease_name: str
    json_disease_ids: tuple[str, ...]
    is_eligible: bool
    quality_reasons: tuple[str, ...] = ()
    applied_default_excluded_count: int = 0
    diagnosis_mismatch: bool = False

    @property
    def ranking_input(self) -> PhenopacketRankingInput:
        """Return strictly isolated ranking inputs."""
        return PhenopacketRankingInput(
            positive_hpo_ids=self.positive_hpo_ids,
            excluded_hpo_ids=self.excluded_hpo_ids,
        )

    @property
    def ground_truth(self) -> PhenopacketGroundTruth:
        """Return evaluation ground truth."""
        return PhenopacketGroundTruth(
            reference_disease_id=self.reference_disease_id,
            disease_name=self.disease_name,
            json_disease_ids=self.json_disease_ids,
        )

    @property
    def target_disease_ids(self) -> frozenset[str]:
        """Return target disease ID set for benchmark evaluation."""
        if self.reference_disease_id:
            return frozenset({self.reference_disease_id})
        return frozenset()

    def to_patient_case(self) -> PatientCase:
        """Export as a validated PatientCase for phenotype ranking."""
        if not self.is_eligible or not self.positive_hpo_ids:
            raise ValueError(
                f"Case '{self.case_id}' is not benchmark-eligible or has no positive HPO IDs."
            )
        return PatientCase(case_id=self.case_id, phenotype_terms=self.positive_hpo_ids)


@dataclass(frozen=True)
class PhenopacketImportResult:
    """Auditable summary and case collection from a Phenopacket Store import."""

    cases: tuple[ImportedPhenopacketCase, ...]
    archive_path: str
    summary_path: str
    total_summary_rows: int
    total_json_files: int
    unmatched_summary_rows: tuple[str, ...] = ()
    unmatched_json_files: tuple[str, ...] = ()
    ambiguous_matches: tuple[tuple[str, tuple[str, ...]], ...] = ()
    duplicate_summary_rows: tuple[str, ...] = ()
    repeated_patient_ids: tuple[tuple[str, int], ...] = ()
    malformed_hpo_ids: tuple[str, ...] = ()
    unknown_hpo_ids: tuple[str, ...] = ()
    diagnosis_mismatches: tuple[tuple[str, str, tuple[str, ...]], ...] = ()

    @property
    def total_cases(self) -> int:
        """Total number of imported case records."""
        return len(self.cases)

    @property
    def eligible_cases(self) -> tuple[ImportedPhenopacketCase, ...]:
        """Cases meeting benchmark eligibility criteria."""
        return tuple(case for case in self.cases if case.is_eligible)

    @property
    def ineligible_cases(self) -> tuple[ImportedPhenopacketCase, ...]:
        """Cases failing benchmark eligibility criteria."""
        return tuple(case for case in self.cases if not case.is_eligible)

    @property
    def eligible_count(self) -> int:
        """Count of benchmark-eligible cases."""
        return len(self.eligible_cases)

    @property
    def ineligible_count(self) -> int:
        """Count of ineligible cases."""
        return len(self.ineligible_cases)

    @property
    def repeated_patient_id_counts(self) -> dict[str, int]:
        """Mapping of repeated patient IDs to occurrence counts."""
        return dict(self.repeated_patient_ids)


def _find_summary_path(
    namelist: list[str],
    summary_filename: str = "phenopacket_store.summary.tsv",
) -> str:
    """Find the summary TSV entry in an archive namelist."""
    candidates = [
        name
        for name in namelist
        if name == summary_filename or name.endswith("/" + summary_filename)
    ]
    if not candidates:
        raise FileNotFoundError(f"Summary file '{summary_filename}' not found in archive.")

    exact_matches = [name for name in candidates if name == summary_filename]
    if len(exact_matches) == 1:
        return exact_matches[0]

    if len(candidates) > 1:
        raise ValueError(
            "Multiple candidate summary files found without an exact-path match: "
            f"{tuple(sorted(candidates))}"
        )

    return candidates[0]


def _build_json_index(namelist: list[str]) -> tuple[dict[str, list[str]], tuple[str, ...]]:
    """Index JSON files by basename without .json extension."""
    json_paths = [
        name
        for name in namelist
        if name.endswith(".json") and not name.endswith("/")
    ]
    json_by_stem: dict[str, list[str]] = defaultdict(list)
    for path_str in json_paths:
        stem = PurePosixPath(path_str).stem
        json_by_stem[stem].append(path_str)

    return json_by_stem, tuple(sorted(json_paths))


def _extract_disease_ids_from_json(json_data: dict[str, Any]) -> tuple[str, ...]:
    """Extract disease identifiers from the diseases array in a phenopacket JSON."""
    raw_diseases = json_data.get("diseases", [])
    if not isinstance(raw_diseases, list):
        return ()

    extracted: list[str] = []
    for entry in raw_diseases:
        if not isinstance(entry, dict):
            continue
        term = entry.get("term")
        if isinstance(term, dict) and term.get("id"):
            disease_id = str(term["id"]).strip()
            if disease_id:
                extracted.append(disease_id)
        elif entry.get("id"):
            disease_id = str(entry["id"]).strip()
            if disease_id:
                extracted.append(disease_id)

    return tuple(dict.fromkeys(extracted))


def _is_term_in_ontology(hpo_id: str, ontology: Any) -> bool:
    """Check if an HPO identifier is known in the supplied ontology."""
    if hasattr(ontology, "get_term"):
        return ontology.get_term(hpo_id) is not None
    if hasattr(ontology, "terms"):
        return hpo_id in ontology.terms
    return False


def _process_phenotypic_features(
    raw_features: Any,
    ontology: HPOOntology | Any | None,
) -> tuple[
    tuple[str, ...],  # positive_hpo_ids
    tuple[str, ...],  # excluded_hpo_ids
    int,              # applied_default_count
    list[str],        # quality_reasons
    list[str],        # malformed_ids
    list[str],        # unknown_ids
]:
    """Extract and validate positive and excluded HPO terms from phenotypicFeatures."""
    positive_ids: list[str] = []
    excluded_ids: list[str] = []
    applied_default_count = 0
    reasons: list[str] = []
    malformed_ids: list[str] = []
    unknown_ids: list[str] = []

    if not isinstance(raw_features, list):
        reasons.append("phenotypicFeatures field is missing or not a list")
        return (), (), 0, reasons, malformed_ids, unknown_ids

    for index, feature in enumerate(raw_features):
        if not isinstance(feature, dict):
            reasons.append(f"phenotypicFeatures[{index}] is not an object")
            continue

        type_obj = feature.get("type")
        if not isinstance(type_obj, dict) or "id" not in type_obj:
            reasons.append(f"phenotypicFeatures[{index}] missing 'type.id'")
            continue

        raw_id = str(type_obj["id"]).strip()

        # Reject malformed HPO IDs
        if not HPO_ID_PATTERN.match(raw_id):
            reasons.append(f"Malformed HPO ID rejected: '{raw_id}'")
            malformed_ids.append(raw_id)
            continue

        # Reject HPO IDs absent from supplied ontology
        if ontology is not None and not _is_term_in_ontology(raw_id, ontology):
            reasons.append(f"HPO ID absent from ontology rejected: '{raw_id}'")
            unknown_ids.append(raw_id)
            continue

        excluded = feature.get("excluded")
        if excluded is True:
            excluded_ids.append(raw_id)
        elif excluded is False:
            positive_ids.append(raw_id)
        elif excluded is None:
            # Default of false applied per Phenopackets schema
            applied_default_count += 1
            positive_ids.append(raw_id)
        else:
            reasons.append(
                f"Non-boolean 'excluded' value ({excluded!r}) on '{raw_id}'; "
                "applied schema default false"
            )
            applied_default_count += 1
            positive_ids.append(raw_id)

    unique_positives = tuple(dict.fromkeys(positive_ids))
    unique_excluded = tuple(dict.fromkeys(excluded_ids))

    conflicts = sorted(set(unique_positives) & set(unique_excluded))
    if conflicts:
        reasons.append(f"Conflicting positive and excluded assertions for: {conflicts}")

    if applied_default_count > 0:
        reasons.append(
            f"Applied schema default (excluded=False) to {applied_default_count} phenotype(s); "
            "not independently confirmed"
        )

    return (
        unique_positives,
        unique_excluded,
        applied_default_count,
        reasons,
        malformed_ids,
        unknown_ids,
    )


def import_phenopackets(
    archive_path: str | Path,
    ontology: HPOOntology | Any | None = None,
    *,
    summary_filename: str = "phenopacket_store.summary.tsv",
    require_diagnosis_match: bool = False,
) -> PhenopacketImportResult:
    """Import Phenopacket Store cases from a ZIP archive without full extraction.

    Ranking inputs are strictly separated from ground-truth diagnoses.
    Cases are validated for benchmark eligibility based on:
    1. A usable case identity (unique matched JSON archive path, no duplicate summary rows).
    2. At least one valid positive HPO ID.
    3. A valid reference disease ID from the summary TSV.
    """
    path = Path(archive_path)
    if not path.is_file():
        raise FileNotFoundError(f"Phenopacket archive not found: {path}")

    with zipfile.ZipFile(path, mode="r") as archive:
        namelist = archive.namelist()
        summary_path = _find_summary_path(namelist, summary_filename)
        json_by_stem, all_json_paths = _build_json_index(namelist)

        with archive.open(summary_path) as summary_file:
            reader = csv.DictReader(
                io.TextIOWrapper(summary_file, encoding="utf-8-sig"),
                delimiter="\t",
            )
            summary_rows = list(reader)

        cases: list[ImportedPhenopacketCase] = []
        matched_json_paths: set[str] = set()
        seen_archive_paths: set[str] = set()

        unmatched_summary_rows: list[str] = []
        ambiguous_matches: list[tuple[str, tuple[str, ...]]] = []
        duplicate_summary_rows: list[str] = []
        diagnosis_mismatches: list[tuple[str, str, tuple[str, ...]]] = []
        all_malformed_hpo_ids: list[str] = []
        all_unknown_hpo_ids: list[str] = []

        for row_index, row in enumerate(summary_rows, start=1):
            raw_filename = (row.get("filename") or "").strip()
            cohort = (row.get("cohort") or "").strip()
            patient_id = (row.get("patient_id") or "").strip()
            disease_id = (row.get("disease_id") or "").strip()
            disease_name = (row.get("disease") or "").strip()
            publication_id = (
                row.get("PMID")
                or row.get("pmid")
                or row.get("publication_id")
                or ""
            ).strip()

            reasons: list[str] = []
            stem = PurePosixPath(raw_filename.replace("\\", "/")).stem if raw_filename else ""
            matching_paths = json_by_stem.get(stem, []) if stem else []

            is_duplicate_summary_row = False
            matched_archive_path = ""

            if not raw_filename or not matching_paths:
                reasons.append(
                    f"No matching JSON file found in archive for filename '{raw_filename}'"
                )
                unmatched_summary_rows.append(raw_filename or f"<row_{row_index}_no_filename>")
                case_id = raw_filename or f"unmatched_row_{row_index}"
            elif len(matching_paths) > 1:
                reasons.append(
                    f"Ambiguous match: multiple JSON files found for stem '{stem}': "
                    f"{tuple(matching_paths)}"
                )
                ambiguous_matches.append((raw_filename, tuple(matching_paths)))
                case_id = raw_filename
            else:
                matched_archive_path = matching_paths[0]
                case_id = matched_archive_path
                matched_json_paths.add(matched_archive_path)

                if matched_archive_path in seen_archive_paths:
                    is_duplicate_summary_row = True
                    duplicate_summary_rows.append(matched_archive_path)
                    reasons.append(
                        f"Duplicate summary row for matched archive path: '{matched_archive_path}'"
                    )
                else:
                    seen_archive_paths.add(matched_archive_path)

            json_data: dict[str, Any] = {}
            if matched_archive_path and not is_duplicate_summary_row:
                try:
                    with archive.open(matched_archive_path) as json_file:
                        json_data = json.load(json_file)
                except Exception as exc:  # noqa: BLE001
                    reasons.append(f"Failed to parse JSON file '{matched_archive_path}': {exc}")
                    json_data = {}

            if not isinstance(json_data, dict):
                reasons.append(f"Invalid JSON content in '{matched_archive_path}': expected object")
                json_data = {}

            # Fall back to subject ID in JSON if patient_id was not in summary
            subject = json_data.get("subject")
            if not patient_id and isinstance(subject, dict) and subject.get("id"):
                patient_id = str(subject["id"]).strip()

            # Cross-check disease diagnosis
            json_disease_ids = _extract_disease_ids_from_json(json_data)
            diagnosis_mismatch = False
            if matched_archive_path and not is_duplicate_summary_row:
                if not json_disease_ids:
                    reasons.append("JSON diseases field is empty or missing")
                elif disease_id and disease_id not in json_disease_ids:
                    diagnosis_mismatch = True
                    reasons.append(
                        f"Diagnosis mismatch: summary disease_id '{disease_id}' "
                        f"not found in JSON diseases {json_disease_ids}"
                    )
                    diagnosis_mismatches.append(
                        (case_id, disease_id, json_disease_ids)
                    )

            # Extract phenotypes
            raw_features = json_data.get("phenotypicFeatures", [])
            (
                positive_hpo_ids,
                excluded_hpo_ids,
                applied_default_count,
                feature_reasons,
                malformed_ids,
                unknown_ids,
            ) = _process_phenotypic_features(raw_features, ontology)

            reasons.extend(feature_reasons)
            all_malformed_hpo_ids.extend(malformed_ids)
            all_unknown_hpo_ids.extend(unknown_ids)

            # Determine benchmark eligibility
            has_usable_identity = bool(matched_archive_path) and not is_duplicate_summary_row
            has_positive_phenotypes = len(positive_hpo_ids) > 0
            has_reference_disease = bool(disease_id)
            has_conflicting_phenotypes = bool(set(positive_hpo_ids) & set(excluded_hpo_ids))

            if not has_usable_identity:
                reasons.append("Ineligible: unusable case identity")
            if not has_positive_phenotypes:
                reasons.append(
                    "Ineligible: no valid positive HPO phenotypes eligible for candidate ranking"
                )
            if not has_reference_disease:
                reasons.append("Ineligible: missing reference disease ID in summary")
            if has_conflicting_phenotypes:
                reasons.append(
                    "Ineligible: conflicting positive and excluded phenotype assertions"
                )

            diagnosis_match_valid = True
            if require_diagnosis_match:
                if not has_reference_disease:
                    diagnosis_match_valid = False
                if not json_disease_ids:
                    diagnosis_match_valid = False
                    reasons.append(
                        "Ineligible: JSON diseases list is empty or missing in strict mode"
                    )
                if diagnosis_mismatch:
                    diagnosis_match_valid = False
                    reasons.append(
                        "Ineligible: reference diagnosis mismatches JSON diseases"
                    )

            is_eligible = (
                has_usable_identity
                and has_positive_phenotypes
                and has_reference_disease
                and not has_conflicting_phenotypes
                and diagnosis_match_valid
            )

            cases.append(
                ImportedPhenopacketCase(
                    case_id=case_id,
                    archive_path=matched_archive_path,
                    cohort=cohort,
                    patient_id=patient_id,
                    publication_id=publication_id,
                    positive_hpo_ids=positive_hpo_ids,
                    excluded_hpo_ids=excluded_hpo_ids,
                    reference_disease_id=disease_id,
                    disease_name=disease_name,
                    json_disease_ids=json_disease_ids,
                    is_eligible=is_eligible,
                    quality_reasons=tuple(reasons),
                    applied_default_excluded_count=applied_default_count,
                    diagnosis_mismatch=diagnosis_mismatch,
                )
            )

        unmatched_json_files = tuple(sorted(set(all_json_paths) - matched_json_paths))

        patient_counts = Counter(c.patient_id for c in cases if c.patient_id)
        repeated_patient_ids = tuple(
            sorted(
                ((pid, count) for pid, count in patient_counts.items() if count > 1),
                key=lambda item: (-item[1], item[0]),
            )
        )

        return PhenopacketImportResult(
            cases=tuple(cases),
            archive_path=str(path),
            summary_path=summary_path,
            total_summary_rows=len(summary_rows),
            total_json_files=len(all_json_paths),
            unmatched_summary_rows=tuple(unmatched_summary_rows),
            unmatched_json_files=unmatched_json_files,
            ambiguous_matches=tuple(ambiguous_matches),
            duplicate_summary_rows=tuple(duplicate_summary_rows),
            repeated_patient_ids=repeated_patient_ids,
            malformed_hpo_ids=tuple(sorted(set(all_malformed_hpo_ids))),
            unknown_hpo_ids=tuple(sorted(set(all_unknown_hpo_ids))),
            diagnosis_mismatches=tuple(diagnosis_mismatches),
        )
