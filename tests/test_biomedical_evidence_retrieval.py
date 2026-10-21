import pytest

from rare_disease_ai.biomedical_evidence_retrieval import (
    EvidenceDocument,
    retrieve_evidence,
)


def make_document(
    evidence_id: str,
    title: str,
    text: str,
    *,
    entity_ids: tuple[str, ...] = ("GENE:1",),
    evidence_type: str = "gene_function",
) -> EvidenceDocument:
    return EvidenceDocument(
        evidence_id=evidence_id,
        title=title,
        text=text,
        source_name="Synthetic test source",
        source_release="test-1",
        source_url="https://example.org/evidence",
        evidence_type=evidence_type,
        entity_ids=entity_ids,
        reference_ids=("PMID:12345",),
    )


def test_retrieval_ranks_by_query_token_coverage_and_preserves_provenance():
    partial = make_document("E2", "Gene function", "Supports a cellular process")
    exact = make_document("E1", "Gene cellular process", "Gene function cellular process")

    report = retrieve_evidence("gene cellular process", (partial, exact))

    assert [hit.document.evidence_id for hit in report.hits] == ["E1", "E2"]
    assert report.hits[0].score == 1.0
    assert report.hits[0].matched_terms == ("cellular", "gene", "process")
    assert report.hits[0].document.source_release == "test-1"
    assert report.hits[0].document.source_url == "https://example.org/evidence"
    assert report.research_only is True


def test_retrieval_supports_exact_entity_and_type_filters():
    first = make_document("E1", "Gene process", "cellular process", entity_ids=("GENE:1",))
    second = make_document(
        "E2",
        "Gene process",
        "cellular process",
        entity_ids=("GENE:2",),
        evidence_type="pathway",
    )

    report = retrieve_evidence(
        "gene process",
        (first, second),
        entity_ids=("GENE:2",),
        evidence_types=("pathway",),
    )

    assert [hit.document.evidence_id for hit in report.hits] == ["E2"]
    assert report.eligible_documents == 1


def test_retrieval_returns_empty_hits_when_nothing_matches():
    document = make_document("E1", "Gene function", "cellular process")

    report = retrieve_evidence("mitochondrial transport", (document,))

    assert report.hits == ()
    assert report.total_documents == 1
    assert report.eligible_documents == 1


def test_retrieval_obeys_top_k():
    documents = tuple(
        make_document(f"E{i}", "Gene process", "cellular process")
        for i in range(3)
    )

    report = retrieve_evidence("gene process", documents, top_k=2)

    assert len(report.hits) == 2
    assert [hit.document.evidence_id for hit in report.hits] == ["E0", "E1"]


@pytest.mark.parametrize("query", ["", "   ", "the and of"])
def test_retrieval_rejects_empty_or_unsearchable_queries(query):
    with pytest.raises(ValueError):
        retrieve_evidence(query, ())


@pytest.mark.parametrize("top_k", [0, -1, True, 1.5])
def test_retrieval_rejects_invalid_top_k(top_k):
    with pytest.raises(ValueError):
        retrieve_evidence("gene", (), top_k=top_k)


def test_evidence_document_requires_valid_source_url():
    with pytest.raises(ValueError, match="source_url"):
        EvidenceDocument(
            evidence_id="E1",
            title="Title",
            text="Evidence text",
            source_name="Source",
            source_release="v1",
            source_url="not-a-url",
            evidence_type="gene_function",
        )


def test_retrieval_rejects_duplicate_evidence_ids():
    document = make_document("E1", "Gene process", "cellular process")

    with pytest.raises(ValueError, match="unique"):
        retrieve_evidence("gene", (document, document))
