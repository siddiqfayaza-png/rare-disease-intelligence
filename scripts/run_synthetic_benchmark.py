"""Run a synthetic ranking benchmark for software validation only."""

from __future__ import annotations

import json
from pathlib import Path

from build_synthetic_benchmark import build_synthetic_disease_index

from rare_disease_ai.evaluation import evaluate_rankings
from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology
from rare_disease_ai.phenotype_ranking import rank_diseases

ROOT = Path(__file__).resolve().parents[1]
CASES_PATH = ROOT / "benchmarks" / "synthetic_cases.json"
ONTOLOGY_PATH = ROOT / "data" / "raw" / "hp-base.obo"


def main() -> None:
    """Run all ranking variants and print their evaluation metrics."""
    benchmark = json.loads(CASES_PATH.read_text(encoding="utf-8"))

    if benchmark.get("clinical_performance_claim") is not False:
        raise ValueError("Synthetic benchmark must not claim clinical performance.")

    ontology = HPOOntology.from_obo(ONTOLOGY_PATH)
    graph = HPOGraph(ontology)
    disease_index = build_synthetic_disease_index()

    cases = benchmark["cases"]
    ranked_methods = {
        "baseline": {},
        "information_content": {"use_information_content": True},
        "semantic_similarity": {"use_semantic_similarity": True},
    }

    target_ids = tuple(frozenset(case["target_disease_ids"]) for case in cases)

    print(f"Benchmark: {benchmark['benchmark_id']}")
    print(f"Purpose: {benchmark['purpose']}")
    print(f"Cases: {len(cases)}")
    print("WARNING: synthetic software validation only; not clinical evidence.")

    for method_name, options in ranked_methods.items():
        rankings = []

        for case in cases:
            results = rank_diseases(
                tuple(case["patient_hpo_ids"]),
                disease_index,
                graph,
                **options,
            )
            rankings.append(tuple(result.database_id for result in results))

        metrics = evaluate_rankings(tuple(rankings), target_ids)
        ranked_methods[method_name]["metrics"] = metrics.__dict__

        print(f"\n{method_name}:")
        print(json.dumps(metrics.__dict__, indent=2))

    output_path = ROOT / "benchmarks" / "synthetic_results.json"
    output_path.write_text(
        json.dumps(
            {
                "benchmark_id": benchmark["benchmark_id"],
                "benchmark_type": benchmark["benchmark_type"],
                "clinical_performance_claim": False,
                "results": ranked_methods,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"\nResults saved to: {output_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
