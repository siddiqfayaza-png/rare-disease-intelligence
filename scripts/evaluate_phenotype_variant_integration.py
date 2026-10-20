"""Reproducible real-data audit and synthetic workflow test for phenotype/variant integration.

This script measures HPOA-to-GenCC evidence coverage using local source datasets,
then exercises the existing phenotype-ranking, variant-prioritization, and
evidence-integration interfaces with explicitly synthetic variant records.

It does NOT measure clinical diagnostic accuracy or variant pathogenicity.
"""

from __future__ import annotations

import hashlib
import json
import platform
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from rare_disease_ai.disease_phenotypes import build_disease_phenotype_index
from rare_disease_ai.gencc_evidence_adapter import adapt_gencc_assertions
from rare_disease_ai.gencc_import import read_gencc_csv
from rare_disease_ai.hpo_annotations import read_hpo_annotations
from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology
from rare_disease_ai.phenotype_ranking import rank_diseases
from rare_disease_ai.phenotype_variant_integration import (
    EVIDENCE_CATEGORIES,
    GeneDiseaseEvidence,
    integrate_ranked_evidence,
    integration_report_to_dict,
)
from rare_disease_ai.variant_prioritization import VariantRecord, prioritize_variants

ROOT = Path(__file__).resolve().parents[1]
HPOA_PATH = ROOT / "data" / "external" / "phenotype.hpoa"
GENCC_PATH = ROOT / "data" / "external" / "gencc-submissions.csv"
ONTOLOGY_PATH = ROOT / "data" / "raw" / "hp-base.obo"
OUTPUT_PATH = ROOT / "reports" / "phase4_integration_evaluation.json"


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 checksum for a local input file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_inputs(paths: tuple[Path, ...]) -> None:
    """Fail early with an actionable error when a required local input is absent."""
    missing = [str(path.relative_to(ROOT)) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Required local dataset(s) not found: "
            + ", ".join(missing)
            + ". Verify that the files exist at these paths before rerunning."
        )


def count_classification_conflicts(
    evidence: tuple[GeneDiseaseEvidence, ...],
) -> tuple[int, int]:
    """Return distinct matched pairs and pairs with multiple raw classifications."""
    grouped: dict[tuple[str, str], set[str]] = defaultdict(set)
    for item in evidence:
        classification = (item.classification or "").strip() or "Unspecified"
        grouped[(item.disease_id, item.gene_symbol.strip().upper())].add(classification)
    conflicts = sum(len(classifications) > 1 for classifications in grouped.values())
    return len(grouped), conflicts


def run_synthetic_workflow(
    disease_index: dict[str, Any],
    hpo_graph: HPOGraph,
    evidence: tuple[GeneDiseaseEvidence, ...],
) -> dict[str, Any]:
    """Exercise ranking and integration with synthetic variants, not patient data."""
    evidence_by_disease: dict[str, list[GeneDiseaseEvidence]] = defaultdict(list)
    for item in evidence:
        evidence_by_disease[item.disease_id].append(item)

    # Select deterministically from real indexed disease evidence, then use that
    # disease's phenotype annotations only to construct a synthetic smoke-test case.
    selected_disease_id: str | None = None
    selected_gene: str | None = None

    for disease_id in sorted(evidence_by_disease):
        profile = disease_index.get(disease_id)
        if profile is None or not profile.positive_hpo_ids:
            continue

        candidate_genes = sorted(
            {
                item.gene_symbol.strip().upper()
                for item in evidence_by_disease[disease_id]
                if item.gene_symbol.strip()
            }
        )
        if candidate_genes:
            selected_disease_id = disease_id
            selected_gene = candidate_genes[0]
            break

    if selected_disease_id is None or selected_gene is None:
        raise RuntimeError(
            "No disease with both positive HPO annotations and matched GenCC "
            "evidence was available for the synthetic integration smoke test."
        )

    profile = disease_index[selected_disease_id]
    # Limit work while keeping the synthetic phenotype input deterministic.
    synthetic_hpo_ids = tuple(sorted(profile.positive_hpo_ids)[:20])
    rankings = rank_diseases(synthetic_hpo_ids, disease_index, hpo_graph)

    if not rankings:
        raise RuntimeError("Phenotype ranking returned no candidates for the synthetic case.")

    synthetic_variants = (
        VariantRecord(
            chromosome="1",
            position=100001,
            reference="A",
            alternate="T",
            gene_symbol=selected_gene,
            consequence="missense_variant",
            zygosity="heterozygous",
            allele_frequency=0.0001,
        ),
        VariantRecord(
            chromosome="1",
            position=100002,
            reference="C",
            alternate="G",
            gene_symbol="SYNTHETIC_UNMATCHED_GENE",
            consequence="synonymous_variant",
            zygosity="heterozygous",
            allele_frequency=None,
        ),
    )
    variant_rankings = prioritize_variants(synthetic_variants)

    # Use only real GenCC assertions whose disease ID and gene match the selected
    # synthetic case. Their provenance is preserved; this is not causal evidence.
    selected_evidence = tuple(
        item
        for item in evidence
        if item.disease_id == selected_disease_id
        and item.gene_symbol.strip().upper() == selected_gene
    )
    integrated = integrate_ranked_evidence(
        rankings,
        variant_rankings,
        selected_evidence,
    )
    disease_result = next(
        (item for item in integrated if item.disease_id == selected_disease_id),
        None,
    )
    report = integration_report_to_dict(
        tuple(item for item in integrated if item.disease_id == selected_disease_id)
    )

    matched_assertion_count = (
        len(disease_result.matched_variants) if disease_result is not None else 0
    )
    categories = dict(disease_result.evidence_summary) if disease_result is not None else {}

    return {
        "test_type": "synthetic_workflow_smoke_test",
        "synthetic_input": True,
        "clinical_accuracy_evaluated": False,
        "selected_disease_id_for_test_construction": selected_disease_id,
        "selected_gene_symbol_for_test_construction": selected_gene,
        "synthetic_phenotype_count": len(synthetic_hpo_ids),
        "phenotype_candidates_returned": len(rankings),
        "synthetic_variants_supplied": len(synthetic_variants),
        "valid_variants_prioritized": len(variant_rankings),
        "selected_disease_present_in_ranking": disease_result is not None,
        "matched_assertions_for_selected_disease_gene": matched_assertion_count,
        "evidence_categories_for_selected_disease": {
            category: categories.get(category, 0) for category in EVIDENCE_CATEGORIES
        },
        "integration_report": report,
        "limitations": [
            "Synthetic variant records are software-test fixtures, not patient data.",
            (
                "The disease and gene are selected from existing indexed "
                "annotations and GenCC evidence to exercise interfaces."
            ),
            ("The test is not independent and must not be interpreted as a diagnostic benchmark."),
            ("Phenotype and variant heuristic scores remain separate and are not probabilities."),
            ("GenCC assertions and categories are not independently adjudicated by this workflow."),
        ],
    }


def main() -> int:
    """Run the data audit, synthetic workflow test, and write a reproducible report."""
    require_inputs((HPOA_PATH, GENCC_PATH, ONTOLOGY_PATH))

    print("Loading HPO ontology...", flush=True)
    ontology = HPOOntology.from_obo(ONTOLOGY_PATH)
    graph = HPOGraph(ontology)

    print("Loading HPOA disease annotations...", flush=True)
    annotations = read_hpo_annotations(HPOA_PATH)
    disease_index = build_disease_phenotype_index(annotations)

    print("Loading GenCC assertions...", flush=True)
    gencc_result = read_gencc_csv(GENCC_PATH)
    if gencc_result.row_errors:
        raise RuntimeError(
            f"GenCC import produced {len(gencc_result.row_errors)} row errors; "
            "resolve them before generating the evaluation report."
        )

    adapter = adapt_gencc_assertions(gencc_result.assertions, set(disease_index))
    distinct_pairs, conflict_pairs = count_classification_conflicts(adapter.evidence)

    print("Running synthetic integration workflow...", flush=True)
    workflow = run_synthetic_workflow(disease_index, graph, adapter.evidence)

    input_metadata = {}
    for name, path in (
        ("hpo_ontology", ONTOLOGY_PATH),
        ("hpoa_annotations", HPOA_PATH),
        ("gencc_export", GENCC_PATH),
    ):
        input_metadata[name] = {
            "path": str(path.relative_to(ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }

    report = {
        "report_schema_version": "1.0",
        "report_date_utc": datetime.now(UTC).date().isoformat(),
        "research_only": True,
        "clinical_validation": False,
        "environment": {
            "python_version": platform.python_version(),
            "platform": platform.platform(),
        },
        "inputs": input_metadata,
        "dataset_counts": {
            "hpo_ontology_terms": len(ontology.terms),
            "hpoa_annotation_records": len(annotations),
            "hpoa_indexed_disease_ids": len(disease_index),
            "gencc_rows_read": gencc_result.rows_read,
            "gencc_assertions_imported": len(gencc_result.assertions),
            "gencc_import_errors": len(gencc_result.row_errors),
        },
        "gencc_matching_audit": {
            "matched_assertions": adapter.matched_assertions,
            "unmatched_assertions": adapter.unmatched_assertions,
            "matched_by_original_id": adapter.matched_by_original_id,
            "matched_by_canonical_id": adapter.matched_by_canonical_id,
            "distinct_matched_disease_gene_pairs": distinct_pairs,
            "matched_pairs_with_multiple_raw_classifications": conflict_pairs,
            "matching_policy": (
                "Exact identifier equality only; original disease identifier "
                "preferred; disease-name matching is not used."
            ),
            "evidence_provenance_retained": True,
            "clinical_interpretation_performed": False,
        },
        "synthetic_workflow": workflow,
        "limitations": [
            (
                "GenCC match counts measure identifier integration coverage, "
                "not diagnostic performance."
            ),
            (
                "The GenCC export may omit source assertions, including "
                "assertions from data sources excluded by the export."
            ),
            (
                "Unmatched identifiers are not necessarily incorrect; "
                "identifier namespaces and coverage can differ."
            ),
            (
                "The synthetic workflow is a software integration check, "
                "not an independent evaluation dataset."
            ),
            (
                "No patient-level diagnostic accuracy, sensitivity, "
                "specificity, or clinical utility is estimated."
            ),
        ],
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print("\n--- Phase 4 Integration Evaluation ---")
    print(f"HPO ontology terms: {len(ontology.terms):,}")
    print(f"HPOA annotations: {len(annotations):,}")
    print(f"GenCC assertions: {len(gencc_result.assertions):,}")
    print(f"Matched GenCC assertions: {adapter.matched_assertions:,}")
    print(f"Unmatched GenCC assertions: {adapter.unmatched_assertions:,}")
    print(f"Matched disease-gene pairs: {distinct_pairs:,}")
    print(f"Pairs with conflicting raw classifications: {conflict_pairs:,}")
    print(
        "Synthetic workflow variants prioritized: "
        f"{workflow['valid_variants_prioritized']:,}/"
        f"{workflow['synthetic_variants_supplied']:,}"
    )
    print(f"Report written to: {OUTPUT_PATH.relative_to(ROOT)}")
    print("Clinical accuracy evaluated: No")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
