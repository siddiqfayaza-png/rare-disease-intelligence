import json

import pytest

from rare_disease_ai.biomedical_evidence_retrieval import EvidenceDocument
from rare_disease_ai.evidence_retrieval_service import EvidenceRetrievalService
from rare_disease_ai.evidence_search_index import build_evidence_index


def make_document(evidence_id: str) -> EvidenceDocument:
    return EvidenceDocument(
        evidence_id=evidence_id,
        title="Gene cellular process",
        text="Evidence supports a cellular process",
        source_name="Synthetic test source",
        source_release="test-1",
        source_url="https://example.org/evidence",
        evidence_type="gene_function",
        entity_ids=("GENE:1",),
        reference_ids=("PMID:12345",),
    )


def test_service_uses_in_memory_retrieval():
    documents = (make_document("E1"), make_document("E2"))
    service = EvidenceRetrievalService(documents=documents)

    report = service.search("gene cellular", top_k=1)

    assert len(report.hits) == 1
    assert report.hits[0].document.evidence_id == "E1"
    assert report.research_only is True


def test_service_uses_sqlite_index(tmp_path):
    document = make_document("E1")
    corpus_path = tmp_path / "corpus.jsonl"
    index_path = tmp_path / "index.sqlite3"

    record = {
        **document.__dict__,
        "entity_ids": list(document.entity_ids),
        "reference_ids": list(document.reference_ids),
    }
    corpus_path.write_text(
        json.dumps(record) + "\n",
        encoding="utf-8",
    )
    build_evidence_index(corpus_path, index_path)

    service = EvidenceRetrievalService(index_path=index_path)
    report = service.search("gene cellular")

    assert [hit.document.evidence_id for hit in report.hits] == ["E1"]
    assert report.research_only is True


def test_service_requires_exactly_one_backend():
    with pytest.raises(ValueError, match="exactly one"):
        EvidenceRetrievalService()

    with pytest.raises(ValueError, match="exactly one"):
        EvidenceRetrievalService(documents=(make_document("E1"),), index_path="index.sqlite3")
