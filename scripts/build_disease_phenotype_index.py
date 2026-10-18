"""Build and summarize the disease-to-phenotype index."""

from pathlib import Path

from rare_disease_ai.disease_phenotypes import (
    build_disease_phenotype_index,
)
from rare_disease_ai.hpo_annotations import read_hpo_annotations

ROOT = Path(__file__).resolve().parents[1]
ANNOTATION_PATH = ROOT / "data" / "external" / "phenotype.hpoa"


def main() -> None:
    """Build the index and print basic dataset statistics."""
    if not ANNOTATION_PATH.is_file():
        raise FileNotFoundError(f"HPO annotation file not found: {ANNOTATION_PATH}")

    print("Reading HPO annotations...")
    annotations = read_hpo_annotations(ANNOTATION_PATH)

    print("Building disease-to-phenotype index...")
    index = build_disease_phenotype_index(annotations)

    positive_annotation_count = sum(
        annotation.qualifier.strip().upper() != "NOT" for annotation in annotations
    )
    negative_annotation_count = sum(
        annotation.qualifier.strip().upper() == "NOT" for annotation in annotations
    )

    diseases_with_positive_terms = sum(bool(profile.positive_hpo_ids) for profile in index.values())
    diseases_with_negative_terms = sum(bool(profile.negative_hpo_ids) for profile in index.values())
    diseases_with_both = sum(
        bool(profile.positive_hpo_ids) and bool(profile.negative_hpo_ids)
        for profile in index.values()
    )

    print("\n--- Disease-to-Phenotype Index Report ---")
    print(f"Annotation records: {len(annotations):,}")
    print(f"Unique diseases indexed: {len(index):,}")
    print(f"Positive annotation records: {positive_annotation_count:,}")
    print(f"NOT annotation records: {negative_annotation_count:,}")
    print(f"Diseases with positive terms: {diseases_with_positive_terms:,}")
    print(f"Diseases with NOT terms: {diseases_with_negative_terms:,}")
    print(f"Diseases with both: {diseases_with_both:,}")

    print("\nSample indexed diseases:")
    for disease_id in sorted(index)[:5]:
        profile = index[disease_id]
        print(f"\n{profile.database_id}: {profile.disease_name}")
        print(f"  Positive HPO terms: {len(profile.positive_hpo_ids)}")
        print(f"  Negative HPO terms: {len(profile.negative_hpo_ids)}")


if __name__ == "__main__":
    main()
