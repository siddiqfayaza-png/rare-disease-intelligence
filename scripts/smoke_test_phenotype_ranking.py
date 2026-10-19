"""Compare phenotype-ranking methods using local HPO datasets."""

from pathlib import Path

from rare_disease_ai.case_validation import PatientCase
from rare_disease_ai.disease_phenotypes import (
    build_disease_phenotype_index,
)
from rare_disease_ai.hpo_annotations import read_hpo_annotations
from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology
from rare_disease_ai.phenotype_normalization import normalize_patient_case
from rare_disease_ai.phenotype_ranking import rank_diseases

ROOT = Path(__file__).resolve().parents[1]
ONTOLOGY_PATH = ROOT / "data" / "raw" / "hp-base.obo"
ANNOTATION_PATH = ROOT / "data" / "external" / "phenotype.hpoa"
TOP_K = 10


def print_results(title: str, results: tuple) -> None:
    """Print candidate rankings and summary statistics."""
    print(f"\n--- {title} ---")

    if not results:
        print("No candidate diseases returned.")
        return

    for position, result in enumerate(results, start=1):
        semantic_text = (
            f" | semantic={result.semantic_score:.3f}" if result.semantic_score is not None else ""
        )

        print(
            f"{position}. {result.database_id} | {result.disease_name} "
            f"| score={result.score:.3f}"
            f"{semantic_text}"
            f" | exact={len(result.exact_matches)}"
            f" | related={len(result.related_matches)}"
            f" | unmatched={len(result.unmatched_patient_hpo_ids)}"
            f" | negative={len(result.negative_matches)}"
        )

    highest_score = results[0].score
    tied_top_count = sum(result.score == highest_score for result in results)

    print(f"Highest score: {highest_score:.6f}")
    print(f"Candidates tied at the highest score in this result set: {tied_top_count}")


def main() -> None:
    """Run a comparative end-to-end ranking check."""
    for path in (ONTOLOGY_PATH, ANNOTATION_PATH):
        if not path.is_file():
            raise FileNotFoundError(f"Required dataset not found: {path}")

    print("Loading HPO ontology...")
    ontology = HPOOntology.from_obo(ONTOLOGY_PATH)
    graph = HPOGraph(ontology)

    print("Loading HPO disease annotations...")
    annotations = read_hpo_annotations(ANNOTATION_PATH)
    disease_index = build_disease_phenotype_index(annotations)

    print(f"Ontology terms: {len(ontology.terms):,}")
    print(f"Annotation records: {len(annotations):,}")
    print(f"Disease profiles: {len(disease_index):,}")

    case = PatientCase(
        case_id="synthetic-smoke-test-001",
        phenotype_terms=(
            "Seizures",
            "Developmental delay",
            "Hypotonia",
        ),
    )

    normalized = normalize_patient_case(case, ontology)

    print("\n--- Phenotype normalization ---")
    for record in normalized.records:
        details = f" -> {record.hpo_id} ({record.preferred_label})" if record.hpo_id else ""
        print(f"{record.original_text}: {record.status}{details}")

    patient_hpo_ids = normalized.normalized_hpo_ids

    if not patient_hpo_ids:
        raise RuntimeError("No test phenotypes were successfully normalized.")

    print(f"\nResolved unique HPO terms: {len(patient_hpo_ids)}")

    print("\nCalculating baseline rankings...")
    baseline_results = rank_diseases(
        patient_hpo_ids,
        disease_index,
        graph,
        top_k=TOP_K,
    )
    print_results("Baseline: exact and hierarchical matching", baseline_results)

    print("\nCalculating information-content rankings...")
    ic_results = rank_diseases(
        patient_hpo_ids,
        disease_index,
        graph,
        top_k=TOP_K,
        use_information_content=True,
    )
    print_results("Information-content weighted ranking", ic_results)

    print("\nCalculating semantic similarity scores...")
    semantic_results = rank_diseases(
        patient_hpo_ids,
        disease_index,
        graph,
        top_k=TOP_K,
        use_information_content=True,
        use_semantic_similarity=True,
    )
    print_results(
        "Information-content ranking with semantic scores",
        semantic_results,
    )

    print("\nComparison completed.")
    print("This synthetic example is a software smoke test, not a diagnostic-accuracy evaluation.")


if __name__ == "__main__":
    main()
