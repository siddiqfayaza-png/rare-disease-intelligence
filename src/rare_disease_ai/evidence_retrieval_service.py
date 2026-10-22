"""Unified interface for in-memory and indexed biomedical evidence retrieval."""

from __future__ import annotations

from pathlib import Path

from rare_disease_ai.biomedical_evidence_retrieval import (
    EvidenceDocument,
    EvidenceRetrievalReport,
    retrieve_evidence,
)
from rare_disease_ai.evidence_search_index import (
    DEFAULT_INDEX,
    EvidenceSearchIndex,
)


class EvidenceRetrievalService:
    """Route evidence queries to an in-memory corpus or a SQLite index."""

    def __init__(
        self,
        *,
        documents: tuple[EvidenceDocument, ...] | None = None,
        index_path: str | Path | None = None,
    ) -> None:
        if (documents is None) == (index_path is None):
            raise ValueError("Provide exactly one of documents or index_path")

        self._documents = documents
        self._index = (
            EvidenceSearchIndex(index_path or DEFAULT_INDEX) if index_path is not None else None
        )

    def search(
        self,
        query: str,
        *,
        top_k: int = 5,
        entity_ids: tuple[str, ...] = (),
        evidence_types: tuple[str, ...] = (),
    ) -> EvidenceRetrievalReport:
        """Search the configured evidence source."""
        if self._index is not None:
            return self._index.search(
                query,
                top_k=top_k,
                entity_ids=entity_ids,
                evidence_types=evidence_types,
            )

        assert self._documents is not None
        return retrieve_evidence(
            query,
            self._documents,
            top_k=top_k,
            entity_ids=entity_ids,
            evidence_types=evidence_types,
        )
