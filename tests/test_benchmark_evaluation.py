from dataclasses import dataclass

import pytest

from rare_disease_ai.benchmark_evaluation import (
    filter_rankable_cases,
    publication_grouped_split,
)


@dataclass(frozen=True)
class Case:
    case_id: str
    publication_id: str
    positive_hpo_ids: tuple[str, ...]
    reference_disease_id: str
    is_eligible: bool = True


@pytest.fixture
def cases():
    return (
        Case("a", "PMID1", ("HP:0000001",), "D1"),
        Case("b", "PMID1", ("HP:0000002",), "D1"),
        Case("c", "PMID2", ("HP:0000003",), "D2"),
        Case("d", "PMID3", ("HP:0000004",), "D3"),
        Case("e", "", ("HP:0000005",), "D1"),
        Case("f", "", ("HP:0000006",), "D1"),
    )


def test_publication_grouped_split_is_reproducible_and_keeps_groups_together(cases):
    train, test = publication_grouped_split(cases, test_fraction=0.4, seed="test-seed")
    train_again, test_again = publication_grouped_split(cases, test_fraction=0.4, seed="test-seed")
    assert [c.case_id for c in train] == [c.case_id for c in train_again]
    assert [c.case_id for c in test] == [c.case_id for c in test_again]
    train_publications = {c.publication_id for c in train if c.publication_id}
    test_publications = {c.publication_id for c in test if c.publication_id}
    assert train_publications.isdisjoint(test_publications)


def test_publication_split_rejects_invalid_fraction(cases):
    with pytest.raises(ValueError):
        publication_grouped_split(cases, test_fraction=1.0)


def test_filter_excludes_unrankable_targets_and_ineligible_cases(cases):
    candidates = {"D1": object(), "D2": object()}
    rankable, counts = filter_rankable_cases(cases, candidates)
    assert [c.case_id for c in rankable] == ["a", "b", "c", "e", "f"]
    assert counts == {"reference_disease_not_in_index": 1}


def test_filter_reports_missing_target_and_no_phenotypes():
    cases = (
        Case("missing-target", "PMID1", ("HP:0000001",), "UNKNOWN"),
        Case("no-phenotype", "PMID2", (), "D1"),
        Case("ineligible", "PMID3", ("HP:0000001",), "D1", is_eligible=False),
    )
    rankable, counts = filter_rankable_cases(cases, {"D1": object()})
    assert rankable == ()
    assert counts == {
        "importer_ineligible": 1,
        "no_positive_phenotypes": 1,
        "reference_disease_not_in_index": 1,
    }
