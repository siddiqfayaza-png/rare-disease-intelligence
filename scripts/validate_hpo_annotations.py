"""Validate the local HPO disease-phenotype annotation dataset."""

from pathlib import Path

from rare_disease_ai.hpo_annotations import (
    read_hpo_annotations,
    validate_hpo_annotations,
)
from rare_disease_ai.hpo_ontology import HPOOntology

ROOT = Path(__file__).resolve().parents[1]
ONTOLOGY_PATH = ROOT / "data" / "raw" / "hp-base.obo"
ANNOTATION_PATH = ROOT / "data" / "external" / "phenotype.hpoa"


def main() -> None:
    """Load both datasets, validate annotations, and print a summary."""
    for path in (ONTOLOGY_PATH, ANNOTATION_PATH):
        if not path.is_file():
            raise FileNotFoundError(f"Required data file not found: {path}")

    print("Loading HPO ontology...")
    ontology = HPOOntology.from_obo(ONTOLOGY_PATH)

    print("Reading HPO disease-phenotype annotations...")
    annotations = read_hpo_annotations(ANNOTATION_PATH)

    print("Validating HPO identifiers...")
    report = validate_hpo_annotations(annotations, ontology)

    print("\n--- HPO Annotation Validation Report ---")
    print(f"Ontology terms: {len(ontology.terms):,}")
    print(f"Annotation rows: {report.annotation_count:,}")
    print(f"Unique disease identifiers: {report.unique_disease_count:,}")
    print(f"Annotations with missing disease IDs: {report.missing_disease_ids:,}")
    print(f"Annotations with missing HPO IDs: {report.missing_hpo_ids:,}")
    print(f"Unknown HPO identifiers: {len(report.unknown_hpo_ids):,}")
    print(f"Validation passed: {report.is_valid}")

    if report.unknown_hpo_ids:
        print("\nFirst 20 unknown HPO identifiers:")
        for hpo_id in report.unknown_hpo_ids[:20]:
            print(f"  {hpo_id}")

    if not report.is_valid:
        raise SystemExit("HPO annotation validation failed.")


if __name__ == "__main__":
    main()
