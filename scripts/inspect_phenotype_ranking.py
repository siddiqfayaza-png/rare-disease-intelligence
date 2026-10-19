"""Inspect phenotype overlap behind disease-ranking ties."""

from pathlib import Path

from rare_disease_ai.disease_phenotypes import (
    build_disease_phenotype_index,
)
from rare_disease_ai.hpo_annotations import read_hpo_annotations

ROOT = Path(__file__).resolve().parents[1]
ANNOTATION_PATH = ROOT / "data" / "external" / "phenotype.hpoa"

TEST_HPO_IDS = (
    "HP:0001250",  # Seizure
    "HP:0001263",  # Global developmental delay
    "HP:0001252",  # Hypotonia
)


def main() -> None:
    annotations = read_hpo_annotations(ANNOTATION_PATH)
    disease_index = build_disease_phenotype_index(annotations)

    print(f"Total disease profiles: {len(disease_index)}")

    for hpo_id in TEST_HPO_IDS:
        count = sum(hpo_id in profile.positive_hpo_ids for profile in disease_index.values())
        print(f"{hpo_id}: {count} diseases with an exact annotation")

    exact_match_counts = {0: 0, 1: 0, 2: 0, 3: 0}

    for profile in disease_index.values():
        matches = sum(hpo_id in profile.positive_hpo_ids for hpo_id in TEST_HPO_IDS)
        exact_match_counts[matches] += 1

    print("\nNumber of diseases by exact-match count:")
    for matches, count in exact_match_counts.items():
        print(f"{matches}/3 phenotypes: {count} diseases")


if __name__ == "__main__":
    main()
