"""Normalize patient phenotype labels against the HPO ontology."""

from dataclasses import dataclass

from rare_disease_ai.case_validation import PatientCase
from rare_disease_ai.hpo_ontology import HPOOntology


@dataclass(frozen=True)
class PhenotypeNormalizationRecord:
    """Auditable result for one submitted phenotype."""

    original_text: str
    status: str
    hpo_id: str | None = None
    preferred_label: str | None = None
    matched_synonym: str | None = None
    candidate_hpo_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class NormalizedPatientCase:
    """Normalized phenotype results for a patient case."""

    case_id: str
    records: tuple[PhenotypeNormalizationRecord, ...]

    @property
    def normalized_hpo_ids(self) -> tuple[str, ...]:
        """Return unique resolved HPO IDs in input order."""
        return tuple(
            dict.fromkeys(record.hpo_id for record in self.records if record.hpo_id is not None)
        )


def normalize_patient_case(
    case: PatientCase,
    ontology: HPOOntology,
) -> NormalizedPatientCase:
    """Resolve each phenotype while preserving unresolved inputs."""
    records: list[PhenotypeNormalizationRecord] = []

    for phenotype in case.phenotype_terms:
        result = ontology.normalize(phenotype)

        if result.is_resolved:
            term = result.matched_terms[0]
            matched_synonym = phenotype if phenotype.casefold() != term.name.casefold() else None

            records.append(
                PhenotypeNormalizationRecord(
                    original_text=phenotype,
                    status="resolved",
                    hpo_id=term.hpo_id,
                    preferred_label=term.name,
                    matched_synonym=matched_synonym,
                )
            )
        elif result.is_ambiguous:
            records.append(
                PhenotypeNormalizationRecord(
                    original_text=phenotype,
                    status="ambiguous",
                    candidate_hpo_ids=tuple(term.hpo_id for term in result.matched_terms),
                )
            )
        elif result.matched_terms and result.matched_terms[0].is_obsolete:
            term = result.matched_terms[0]
            records.append(
                PhenotypeNormalizationRecord(
                    original_text=phenotype,
                    status="obsolete",
                    candidate_hpo_ids=term.replaced_by,
                )
            )
        else:
            records.append(
                PhenotypeNormalizationRecord(
                    original_text=phenotype,
                    status="unresolved",
                )
            )

    return NormalizedPatientCase(
        case_id=case.case_id,
        records=tuple(records),
    )
