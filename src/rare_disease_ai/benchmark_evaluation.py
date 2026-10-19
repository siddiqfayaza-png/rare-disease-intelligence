"""Publication-grouped evaluation for the Phenopacket Store ranking benchmark.

Reference diagnoses are used only as evaluation targets, never as ranking inputs.
Metrics are research benchmark results, not clinical validation.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from collections.abc import Iterable
from dataclasses import asdict
from typing import Any

from rare_disease_ai.evaluation import RankingMetrics, evaluate_rankings
from rare_disease_ai.phenopacket_import import ImportedPhenopacketCase
from rare_disease_ai.phenotype_ranking import rank_diseases


def publication_grouped_split(
    cases: Iterable[ImportedPhenopacketCase],
    *,
    test_fraction: float = 0.2,
    seed: str = "rdai-v1",
) -> tuple[tuple[ImportedPhenopacketCase, ...], tuple[ImportedPhenopacketCase, ...]]:
    """Split deterministically by publication, keeping each group together.

    Cases with no publication ID receive unique synthetic groups. If there are
    at least two groups, both partitions are guaranteed to be nonempty.
    """
    if not 0.0 < test_fraction < 1.0:
        raise ValueError("test_fraction must be strictly between 0 and 1.")

    case_list = tuple(cases)
    grouped: dict[str, list[ImportedPhenopacketCase]] = {}

    for case in case_list:
        group = case.publication_id.strip() or f"unknown-publication:{case.case_id}"
        grouped.setdefault(group, []).append(case)

    if len(grouped) < 2:
        raise ValueError("Publication-grouped split requires at least two distinct groups.")

    ordered_groups = sorted(
        grouped,
        key=lambda group: hashlib.sha256(f"{seed}:{group}".encode()).digest(),
    )
    test_group_count = round(len(ordered_groups) * test_fraction)
    test_group_count = max(1, min(test_group_count, len(ordered_groups) - 1))
    test_groups = set(ordered_groups[:test_group_count])

    train = tuple(
        case
        for case in case_list
        if (case.publication_id.strip() or f"unknown-publication:{case.case_id}") not in test_groups
    )
    test = tuple(
        case
        for case in case_list
        if (case.publication_id.strip() or f"unknown-publication:{case.case_id}") in test_groups
    )

    return train, test


def filter_rankable_cases(
    cases: Iterable[ImportedPhenopacketCase],
    disease_index: dict[str, Any],
) -> tuple[tuple[ImportedPhenopacketCase, ...], dict[str, int]]:
    """Keep importer-eligible cases with nonempty phenotypes and indexed targets."""
    selected: list[ImportedPhenopacketCase] = []
    reasons: Counter[str] = Counter()

    for case in cases:
        if not case.is_eligible:
            reasons["importer_ineligible"] += 1
            continue

        if not case.positive_hpo_ids:
            reasons["no_positive_phenotypes"] += 1
            continue

        if not case.reference_disease_id:
            reasons["missing_reference_disease"] += 1
            continue

        if case.reference_disease_id not in disease_index:
            reasons["reference_disease_not_in_index"] += 1
            continue

        selected.append(case)

    return tuple(selected), dict(sorted(reasons.items()))


def metrics_to_dict(metrics: RankingMetrics) -> dict[str, Any]:
    """Convert metrics to JSON-safe primitives."""
    return asdict(metrics)


def evaluate_case_set(
    cases: tuple[ImportedPhenopacketCase, ...],
    disease_index: dict[str, Any],
    hpo_graph: Any,
    *,
    include_semantic: bool = True,
    progress_every: int = 25,
    label: str = "Evaluation",
) -> dict[str, Any]:
    """Evaluate baseline and optional semantic ranking with progress reporting.

    Progress is printed every ``progress_every`` cases and on the final case.
    The ranking algorithm and evaluation metrics are otherwise unchanged.
    """
    if not cases:
        raise ValueError("No eligible cases available for evaluation.")

    if progress_every < 1:
        raise ValueError("progress_every must be at least 1.")

    baseline_rankings: list[tuple[str, ...]] = []
    semantic_rankings: list[tuple[str, ...]] = []
    targets: list[frozenset[str]] = []
    total = len(cases)

    for index, case in enumerate(cases, start=1):
        # Intentionally pass only phenotype IDs to the ranking engine.
        ranking_input = case.ranking_input

        baseline = rank_diseases(
            ranking_input.positive_hpo_ids,
            disease_index,
            hpo_graph,
        )
        baseline_rankings.append(tuple(item.database_id for item in baseline))

        if include_semantic:
            semantic = rank_diseases(
                ranking_input.positive_hpo_ids,
                disease_index,
                hpo_graph,
                use_semantic_similarity=True,
            )
            semantic_rankings.append(tuple(item.database_id for item in semantic))

        targets.append(case.target_disease_ids)

        if index % progress_every == 0 or index == total:
            print(
                f"{label}: {index}/{total} cases evaluated",
                flush=True,
            )

    output: dict[str, Any] = {
        "case_count": total,
        "baseline": metrics_to_dict(
            evaluate_rankings(
                tuple(baseline_rankings),
                tuple(targets),
            )
        ),
    }

    if include_semantic:
        output["semantic"] = metrics_to_dict(
            evaluate_rankings(
                tuple(semantic_rankings),
                tuple(targets),
            )
        )

    return output


def summarize_split_groups(
    cases: Iterable[ImportedPhenopacketCase],
) -> dict[str, int]:
    """Count cases and publication groups without emitting individual case IDs."""
    case_list = tuple(cases)
    groups = {
        case.publication_id.strip() or f"unknown-publication:{case.case_id}" for case in case_list
    }

    return {
        "cases": len(case_list),
        "publication_groups": len(groups),
    }
