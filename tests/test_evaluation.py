import pytest

from rare_disease_ai.evaluation import evaluate_rankings


def test_evaluates_recall_mrr_and_median_rank():
    ranked = (
        ("A", "B", "C"),
        ("B", "A", "C"),
        ("A", "B", "C"),
    )
    targets = (
        frozenset({"A"}),
        frozenset({"A"}),
        frozenset({"D"}),
    )

    metrics = evaluate_rankings(ranked, targets)

    assert metrics.case_count == 3
    assert metrics.recall_at_1 == pytest.approx(1 / 3)
    assert metrics.recall_at_5 == pytest.approx(2 / 3)
    assert metrics.recall_at_10 == pytest.approx(2 / 3)
    assert metrics.recall_at_20 == pytest.approx(2 / 3)
    assert metrics.mean_reciprocal_rank == pytest.approx(1 / 2)
    assert metrics.median_rank == 1.5
    assert metrics.unranked_count == 1


def test_multiple_target_diseases_use_first_relevant_rank():
    ranked = (("A", "B", "C"),)
    targets = (frozenset({"B", "C"}),)

    metrics = evaluate_rankings(ranked, targets)

    assert metrics.recall_at_1 == 0.0
    assert metrics.mean_reciprocal_rank == pytest.approx(1 / 2)
    assert metrics.median_rank == 2
    assert metrics.unranked_count == 0


def test_all_unranked_cases_have_zero_mrr_and_no_median():
    metrics = evaluate_rankings(
        (("A", "B"), ("C",)),
        (frozenset({"X"}), frozenset({"Y"})),
    )

    assert metrics.mean_reciprocal_rank == 0.0
    assert metrics.median_rank is None
    assert metrics.unranked_count == 2


def test_rejects_mismatched_case_counts():
    with pytest.raises(ValueError, match="equal lengths"):
        evaluate_rankings(
            (("A",),),
            (),
        )


def test_rejects_empty_evaluation_cohort():
    with pytest.raises(ValueError, match="At least one"):
        evaluate_rankings((), ())


def test_rejects_case_without_target_disease():
    with pytest.raises(ValueError, match="no target disease IDs"):
        evaluate_rankings(
            (("A",),),
            (frozenset(),),
        )
