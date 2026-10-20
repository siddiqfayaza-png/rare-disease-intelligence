"""Run reproducible publication-grouped evaluation on Phenopacket Store cases."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from rare_disease_ai.benchmark_evaluation import (  # noqa: E402
    evaluate_case_set,
    filter_rankable_cases,
    publication_grouped_split,
    summarize_split_groups,
)
from rare_disease_ai.disease_phenotypes import build_disease_phenotype_index  # noqa: E402
from rare_disease_ai.hpo_annotations import read_hpo_annotations  # noqa: E402
from rare_disease_ai.hpo_graph import HPOGraph  # noqa: E402
from rare_disease_ai.hpo_ontology import HPOOntology  # noqa: E402
from rare_disease_ai.phenopacket_import import import_phenopackets  # noqa: E402

DEFAULT_ARCHIVE = ROOT / "data" / "external" / "all_phenopackets-0.1.27.zip"
DEFAULT_ONTOLOGY = ROOT / "data" / "raw" / "hp-base.obo"
DEFAULT_ANNOTATIONS = ROOT / "data" / "external" / "phenotype.hpoa"
DEFAULT_OUTPUT = ROOT / "reports" / "phenopacket_ranking_evaluation.json"
DEFAULT_DEVELOPMENT_OUTPUT = ROOT / "reports" / "phenopacket_development_evaluation.json"


def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--ontology", type=Path, default=DEFAULT_ONTOLOGY)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--test-fraction", type=float, default=0.2)
    parser.add_argument("--seed", default="rdai-v1")
    parser.add_argument(
        "--development-only",
        action="store_true",
        help=(
            "Evaluate only development publication groups. "
            "When using the default output path, save to a separate development report."
        ),
    )
    parser.add_argument(
        "--skip-semantic",
        action="store_true",
        help="Run only the exact/hierarchical baseline.",
    )
    parser.add_argument(
        "--include-overall",
        action="store_true",
        help=(
            "Also evaluate every rankable case for descriptive comparison. "
            "This may take a long time."
        ),
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=25,
        help="Print progress after this many cases.",
    )
    return parser.parse_args(args)


def main(args: list[str] | None = None) -> None:
    opts = parse_args(args)

    if not 0.0 < opts.test_fraction < 1.0:
        raise SystemExit("--test-fraction must be strictly between 0 and 1.")

    if opts.progress_every < 1:
        raise SystemExit("--progress-every must be at least 1.")

    for required in (opts.archive, opts.ontology, opts.annotations):
        if not required.is_file():
            raise SystemExit(f"Required file not found: {required}")

    print("Loading ontology and disease annotations...", flush=True)
    ontology = HPOOntology.from_obo(opts.ontology)
    graph = HPOGraph(ontology)
    annotations = read_hpo_annotations(opts.annotations)
    disease_index = build_disease_phenotype_index(annotations)

    print("Importing Phenopacket Store cases with strict diagnosis matching...", flush=True)
    imported = import_phenopackets(
        opts.archive,
        ontology,
        require_diagnosis_match=True,
    )
    rankable, exclusion_counts = filter_rankable_cases(
        imported.cases,
        disease_index,
    )

    if len(rankable) < 2:
        raise SystemExit(f"Not enough rankable cases for evaluation: {len(rankable)}")

    train, test = publication_grouped_split(
        rankable,
        test_fraction=opts.test_fraction,
        seed=opts.seed,
    )

    print(
        f"Imported={imported.total_cases}; "
        f"rankable={len(rankable)}; "
        f"excluded={len(imported.cases) - len(rankable)}; "
        f"development-side={len(train)}; "
        f"held-out={len(test)}",
        flush=True,
    )

    if opts.development_only:
        evaluation_cases = train
        evaluation_label = "Development"
        evaluation_mode = "development"
    else:
        evaluation_cases = test
        evaluation_label = "Held-out"
        evaluation_mode = "heldout"

    print(
        f"Evaluating {evaluation_label.lower()} publication groups; "
        "diagnosis labels are not ranking inputs...",
        flush=True,
    )
    evaluation_metrics = evaluate_case_set(
        evaluation_cases,
        disease_index,
        graph,
        include_semantic=not opts.skip_semantic,
        progress_every=opts.progress_every,
        label=evaluation_label,
    )

    # Full-cohort results are optional and descriptive only.
    overall_metrics = None
    if opts.include_overall:
        print(
            "Evaluating all rankable cases for descriptive comparison...",
            flush=True,
        )
        overall_metrics = evaluate_case_set(
            rankable,
            disease_index,
            graph,
            include_semantic=not opts.skip_semantic,
            progress_every=opts.progress_every,
            label="Overall",
        )
    else:
        print(
            "Skipping full-dataset evaluation. Use --include-overall to enable it.",
            flush=True,
        )

    report = {
        "report_schema_version": "1.1",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "research_only": True,
        "clinical_validation": False,
        "evaluation_mode": evaluation_mode,
        "source": {
            "archive_filename": opts.archive.name,
            "archive_sha256": _sha256(opts.archive),
            "ontology_filename": opts.ontology.name,
            "ontology_sha256": _sha256(opts.ontology),
            "annotations_filename": opts.annotations.name,
            "annotations_sha256": _sha256(opts.annotations),
        },
        "methodology": {
            "split": "deterministic publication-grouped split",
            "test_fraction_requested": opts.test_fraction,
            "seed": opts.seed,
            "candidate_disease_count": len(disease_index),
            "ranker_inputs": ["positive_hpo_ids"],
            "ground_truth_used_only_for_evaluation": True,
            "excluded_phenotypes_used_for_scoring": False,
            "semantic_method_included": not opts.skip_semantic,
            "tie_break": "disease ID ascending after score ordering",
            "metrics": [
                "case_count",
                "recall_at_1",
                "recall_at_5",
                "recall_at_10",
                "recall_at_20",
                "mean_reciprocal_rank",
                "median_rank",
                "unranked_count",
            ],
            "full_cohort_evaluation_requested": opts.include_overall,
            "development_only": opts.development_only,
        },
        "cohort": {
            "imported_cases": imported.total_cases,
            "importer_eligible_cases": imported.eligible_count,
            "rankable_cases": len(rankable),
            "excluded_from_ranking_metrics": len(imported.cases) - len(rankable),
            "exclusion_counts_nonexclusive_stage_reasons": exclusion_counts,
            "publication_split": {
                "development_side": summarize_split_groups(train),
                "heldout_test": summarize_split_groups(test),
            },
            "evaluated_partition": summarize_split_groups(evaluation_cases),
            "import_integrity": {
                "unmatched_summary_rows": len(imported.unmatched_summary_rows),
                "unmatched_json_files": len(imported.unmatched_json_files),
                "ambiguous_matches": len(imported.ambiguous_matches),
                "duplicate_summary_rows": len(imported.duplicate_summary_rows),
                "diagnosis_mismatches": len(imported.diagnosis_mismatches),
                "malformed_hpo_ids": len(imported.malformed_hpo_ids),
                "unknown_hpo_ids": len(imported.unknown_hpo_ids),
            },
        },
        "primary_heldout_metrics": (None if opts.development_only else evaluation_metrics),
        "primary_development_metrics": (evaluation_metrics if opts.development_only else None),
        "descriptive_all_rankable_metrics": overall_metrics,
        "limitations": [
            "This is retrospective benchmark evaluation, not prospective or clinical validation.",
            "Publication grouping reduces within-publication leakage but does not prove "
            "patient-level independence or eliminate all phenotype overlap.",
            "The ranker is fixed; this evaluation does not train a model.",
            "The development partition is intended for experimentation. Repeatedly "
            "consulting held-out results to choose methods compromises their role as "
            "an untouched final evaluation.",
            "The disease index is derived from external HPO annotations; diseases absent "
            "from that index cannot be ranked and are excluded from ranking metrics.",
            "Phenotypes whose excluded field uses the schema default may not be "
            "independently confirmed positives.",
            "The semantic method is experimental and must be interpreted alongside the baseline.",
            "Full-cohort metrics are descriptive only.",
        ],
    }

    output_path = opts.output
    if opts.development_only and output_path == DEFAULT_OUTPUT:
        output_path = DEFAULT_DEVELOPMENT_OUTPUT

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print(f"\n{evaluation_label} metrics:", flush=True)
    print(json.dumps(evaluation_metrics, indent=2), flush=True)
    print(f"\nReport saved to: {output_path.resolve()}", flush=True)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


if __name__ == "__main__":
    main()
