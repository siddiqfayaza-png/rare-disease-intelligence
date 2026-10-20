
from rare_disease_ai.gencc_evidence_adapter import adapt_gencc_assertions
from rare_disease_ai.gencc_import import GenCCAssertion


def make_assertion(
    *,
    submission_id: str = "SGC-1",
    original_id: str = "OMIM:123456",
    canonical_id: str = "MONDO:0000001",
    classification: str = "Definitive",
) -> GenCCAssertion:
    return GenCCAssertion(
        submission_id=submission_id,
        version_number=1,
        gene_curie="HGNC:1",
        gene_symbol="GENE1",
        disease_curie=canonical_id,
        disease_title="Example disease",
        disease_original_curie=original_id,
        classification_curie="GENCC:CLASS",
        classification_title=classification,
        moi_curie="HP:0000006",
        moi_title="Autosomal dominant",
        submitter_curie="GENCC:000101",
        submitter_title="Example submitter",
        submitted_as_date="2025-01-01",
        public_report_url="https://example.org/report",
        pmids=("PMID:12345",),
        assertion_criteria_url="https://example.org/criteria",
    )


def test_prefers_exact_original_disease_identifier() -> None:
    result = adapt_gencc_assertions(
        (make_assertion(),),
        {"OMIM:123456", "MONDO:0000001"},
    )

    assert result.matched_assertions == 1
    assert result.matched_by_original_id == 1
    assert result.matched_by_canonical_id == 0
    assert result.evidence[0].disease_id == "OMIM:123456"
    assert result.evidence[0].evidence_id == "SGC-1"
    assert result.evidence[0].classification == "Definitive"
    assert result.evidence[0].evidence_category == "supporting"
    assert result.evidence[0].submitter == "Example submitter"
    assert result.evidence[0].pmids == ("PMID:12345",)


def test_uses_canonical_id_when_original_id_is_not_indexed() -> None:
    result = adapt_gencc_assertions(
        (make_assertion(),),
        {"MONDO:0000001"},
    )

    assert result.matched_by_canonical_id == 1
    assert result.evidence[0].disease_id == "MONDO:0000001"


def test_does_not_match_by_disease_name() -> None:
    result = adapt_gencc_assertions(
        (make_assertion(),),
        {"Example disease"},
    )

    assert result.matched_assertions == 0
    assert result.unmatched_assertions == 1
    assert result.evidence == ()


def test_preserves_disputed_classification() -> None:
    result = adapt_gencc_assertions(
        (make_assertion(classification="Disputed Evidence"),),
        {"OMIM:123456"},
    )

    assert result.matched_assertions == 1
    assert result.evidence[0].classification == "Disputed Evidence"
    assert result.evidence[0].evidence_category == "disputed"


def test_reports_mixed_matched_and_unmatched_assertions() -> None:
    assertions = (
        make_assertion(submission_id="SGC-1"),
        make_assertion(
            submission_id="SGC-2",
            original_id="OMIM:999999",
            canonical_id="MONDO:9999999",
        ),
    )

    result = adapt_gencc_assertions(assertions, {"OMIM:123456"})

    assert result.matched_assertions == 1
    assert result.unmatched_assertions == 1
    assert result.evidence[0].evidence_id == "SGC-1"
