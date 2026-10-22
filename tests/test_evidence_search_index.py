import json
import sqlite3

import pytest

from rare_disease_ai.biomedical_evidence_retrieval import (
    EvidenceDocument,
    retrieve_evidence,
)
from rare_disease_ai.evidence_search_index import (
    EvidenceSearchIndex,
    build_evidence_index,
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


def write_corpus(path, documents):
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for document in documents:
            payload = {
                "evidence_id": document.evidence_id,
                "title": document.title,
                "text": document.text,
                "source_name": document.source_name,
                "source_release": document.source_release,
                "source_url": document.source_url,
                "evidence_type": document.evidence_type,
                "entity_ids": list(document.entity_ids),
                "reference_ids": list(document.reference_ids),
            }
            handle.write(json.dumps(payload, ensure_ascii=False) + "\n")


@pytest.fixture
def indexed_corpus(tmp_path):
    documents = (
        make_document("E2", "Gene function", "Supports a cellular process"),
        make_document("E1", "Gene cellular process", "Gene function cellular process"),
        make_document(
            "E3",
            "Mitochondrial transport",
            "Mitochondrial transport pathway",
            entity_ids=("GENE:2",),
            evidence_type="pathway",
        ),
    )
    corpus_path = tmp_path / "corpus.jsonl"
    index_path = tmp_path / "index.sqlite3"
    write_corpus(corpus_path, documents)
    build_evidence_index(corpus_path, index_path, batch_size=2)
    return documents, index_path


def test_indexed_search_matches_baseline_ranking_and_scores(indexed_corpus):
    documents, index_path = indexed_corpus
    index = EvidenceSearchIndex(index_path)

    expected = retrieve_evidence("gene cellular process", documents)
    actual = index.search("gene cellular process")

    assert [hit.document.evidence_id for hit in actual.hits] == [
        hit.document.evidence_id for hit in expected.hits
    ]
    assert [hit.score for hit in actual.hits] == [hit.score for hit in expected.hits]
    assert [hit.matched_terms for hit in actual.hits] == [
        hit.matched_terms for hit in expected.hits
    ]


def test_indexed_search_preserves_provenance(indexed_corpus):
    _, index_path = indexed_corpus
    report = EvidenceSearchIndex(index_path).search("cellular process")

    assert report.research_only is True
    assert report.hits
    assert report.hits[0].document.source_name == "Synthetic test source"
    assert report.hits[0].document.source_release == "test-1"
    assert report.hits[0].document.source_url == "https://example.org/evidence"


def test_indexed_search_applies_exact_entity_and_type_filters(indexed_corpus):
    _, index_path = indexed_corpus

    report = EvidenceSearchIndex(index_path).search(
        "mitochondrial transport",
        entity_ids=("GENE:2",),
        evidence_types=("pathway",),
    )

    assert [hit.document.evidence_id for hit in report.hits] == ["E3"]
    assert report.eligible_documents == 1


def test_indexed_search_has_deterministic_tie_breaking(tmp_path):
    corpus_path = tmp_path / "corpus.jsonl"
    index_path = tmp_path / "index.sqlite3"
    documents = (
        make_document("E2", "Gene process", "cellular process"),
        make_document("E1", "Gene process", "cellular process"),
    )
    write_corpus(corpus_path, documents)
    build_evidence_index(corpus_path, index_path)

    report = EvidenceSearchIndex(index_path).search("gene process")

    assert [hit.document.evidence_id for hit in report.hits] == ["E1", "E2"]


def test_index_retains_duplicate_evidence_ids(tmp_path):
    corpus_path = tmp_path / "corpus.jsonl"
    index_path = tmp_path / "index.sqlite3"
    documents = (
        make_document("DUP", "Gene process", "cellular process"),
        make_document("DUP", "Gene process", "different source wording"),
    )
    write_corpus(corpus_path, documents)

    result = build_evidence_index(corpus_path, index_path)
    report = EvidenceSearchIndex(index_path).search("gene process")

    assert result["document_count"] == 2
    assert result["duplicate_evidence_id_count"] == 1
    assert len(report.hits) == 2


@pytest.mark.parametrize("query", ["", "   ", "the and of"])
def test_indexed_search_rejects_unsearchable_queries(indexed_corpus, query):
    _, index_path = indexed_corpus

    with pytest.raises(ValueError):
        EvidenceSearchIndex(index_path).search(query)


@pytest.mark.parametrize("top_k", [0, -1, True, 1.5])
def test_indexed_search_rejects_invalid_top_k(indexed_corpus, top_k):
    _, index_path = indexed_corpus

    with pytest.raises(ValueError):
        EvidenceSearchIndex(index_path).search("gene", top_k=top_k)


def test_index_records_corpus_checksum_and_count(indexed_corpus):
    documents, index_path = indexed_corpus
    index = EvidenceSearchIndex(index_path)

    assert index.metadata["document_count"] == str(len(documents))
    assert len(index.metadata["corpus_sha256"]) == 64


def test_failed_rebuild_does_not_replace_existing_index(tmp_path):
    corpus_path = tmp_path / "corpus.jsonl"
    index_path = tmp_path / "index.sqlite3"

    write_corpus(corpus_path, (make_document("E1", "Gene process", "cellular process"),))
    build_evidence_index(corpus_path, index_path)
    before = EvidenceSearchIndex(index_path).search("gene").hits[0].document.evidence_id

    corpus_path.write_text('{"invalid": "record"}\n', encoding="utf-8")

    with pytest.raises(ValueError, match="Invalid evidence record"):
        build_evidence_index(corpus_path, index_path)

    after = EvidenceSearchIndex(index_path).search("gene").hits[0].document.evidence_id
    assert after == before


def test_indexed_search_report_matches_baseline_metadata(indexed_corpus):
    documents, index_path = indexed_corpus
    query = "  GENE cellular process  "

    expected = retrieve_evidence(query, documents, top_k=2)
    actual = EvidenceSearchIndex(index_path).search(query, top_k=2)

    assert actual.query == expected.query
    assert actual.total_documents == expected.total_documents
    assert actual.eligible_documents == expected.eligible_documents
    assert actual.research_only is True
    assert len(actual.hits) == len(expected.hits)


def test_indexed_search_is_case_insensitive_and_deduplicates_query_tokens(
    indexed_corpus,
):
    _, index_path = indexed_corpus
    index = EvidenceSearchIndex(index_path)

    first = index.search("GENE gene CELLULAR")
    second = index.search("gene cellular")

    assert [hit.document.evidence_id for hit in first.hits] == [
        hit.document.evidence_id for hit in second.hits
    ]
    assert [hit.score for hit in first.hits] == [hit.score for hit in second.hits]
    assert all(hit.score <= 1.0 for hit in first.hits)


def test_indexed_search_returns_no_hits_when_filters_exclude_matches(
    indexed_corpus,
):
    _, index_path = indexed_corpus

    report = EvidenceSearchIndex(index_path).search(
        "gene cellular process",
        entity_ids=("GENE:DOES_NOT_EXIST",),
    )

    assert report.hits == ()
    assert report.total_documents == 3
    assert report.eligible_documents == 0


def test_indexed_search_combines_filters_and_respects_top_k(indexed_corpus):
    _, index_path = indexed_corpus
    index = EvidenceSearchIndex(index_path)

    report = index.search(
        "mitochondrial transport",
        top_k=1,
        entity_ids=("GENE:2",),
        evidence_types=("pathway",),
    )

    assert len(report.hits) == 1
    assert report.hits[0].document.evidence_id == "E3"
    assert report.total_documents == 3
    assert report.eligible_documents == 1


def test_search_rejects_missing_index(tmp_path):
    missing_path = tmp_path / "missing.sqlite3"

    with pytest.raises(FileNotFoundError, match="Evidence index not found"):
        EvidenceSearchIndex(missing_path)


def test_search_rejects_unsupported_schema_version(indexed_corpus):
    _, index_path = indexed_corpus

    with sqlite3.connect(index_path) as connection:
        connection.execute(
            "UPDATE metadata SET value = ? WHERE key = ?",
            ("unsupported-version", "schema_version"),
        )

    with pytest.raises(ValueError, match="schema version"):
        EvidenceSearchIndex(index_path)
