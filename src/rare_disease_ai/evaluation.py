"""Evaluation metrics for rare-disease candidate ranking."""

from __future__ import annotations

from dataclasses import dataclass
from statistics import median


@dataclass(frozen=True)
class RankingMetrics:
    """Summary metrics for a disease-ranking benchmark."""

    case_count: int
    recall_at_1: float
    recall_at_5: float
    recall_at_10: float
    recall_at_20: float
    mean_reciprocal_rank: float
    median_rank: float | None
    unranked_count: int


def evaluate_rankings(
    ranked_disease_ids: tuple[tuple[str, ...], ...],
    target_disease_ids: tuple[frozenset[str], ...],
) -> RankingMetrics:
    """Evaluate ranked candidates against one or more target diseases per case.

    Each ranked_disease_ids entry contains disease IDs ordered from highest
    to lowest priority. Each target_disease_ids entry contains the relevant
    disease IDs for the corresponding case.

    A case is a hit at k if any target disease appears in the first k results.
    MRR uses the reciprocal rank of the first relevant disease. If no target
    disease appears in the ranking, that case contributes zero to MRR and is
    counted as unranked.

    Ties are handled by the order provided by the ranking engine. Disease-ID
    tie-breaking should be performed consistently before evaluation.
    """
    if len(ranked_disease_ids) != len(target_disease_ids):
        raise ValueError("Ranked results and target disease sets must have equal lengths.")

    if not ranked_disease_ids:
        raise ValueError("At least one evaluation case is required.")

    reciprocal_ranks: list[float] = []
    target_ranks: list[int | None] = []

    for case_index, (ranked_ids, targets) in enumerate(
        zip(ranked_disease_ids, target_disease_ids, strict=True)
    ):
        if not targets:
            raise ValueError(f"Evaluation case {case_index} has no target disease IDs.")

        first_relevant_rank = next(
            (rank for rank, disease_id in enumerate(ranked_ids, start=1) if disease_id in targets),
            None,
        )

        target_ranks.append(first_relevant_rank)
        reciprocal_ranks.append(
            1.0 / first_relevant_rank if first_relevant_rank is not None else 0.0
        )

    case_count = len(target_ranks)
    unranked_count = sum(rank is None for rank in target_ranks)
    ranked_positions = [rank for rank in target_ranks if rank is not None]

    def recall_at(k: int) -> float:
        """Return the fraction of cases with a target in the top k."""
        hits = sum(rank is not None and rank <= k for rank in target_ranks)
        return hits / case_count

    return RankingMetrics(
        case_count=case_count,
        recall_at_1=recall_at(1),
        recall_at_5=recall_at(5),
        recall_at_10=recall_at(10),
        recall_at_20=recall_at(20),
        mean_reciprocal_rank=sum(reciprocal_ranks) / case_count,
        median_rank=median(ranked_positions) if ranked_positions else None,
        unranked_count=unranked_count,
    )
