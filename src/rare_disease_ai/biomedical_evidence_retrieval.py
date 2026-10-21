from __future__ import annotations

import re
from dataclasses import dataclass
from urllib.parse import urlparse

_TOKEN_PATTERN = re.compile(r"[a-z0-9][a-z0-9._:-]*")
_STOP_WORDS = frozenset(
    {"a", "an", "and", "are", "for", "in", "is", "of", "on", "or", "the", "to", "with"}
)


def _tokens(text: str) -> frozenset[str]:
    return frozenset(
        token
        for token in _TOKEN_PATTERN.findall(text.casefold())
        if len(token) > 1 and token not in _STOP_WORDS
    )


@dataclass(frozen=True)
class EvidenceDocument:
    """One source-backed biomedical evidence record."""

    evidence_id: str
    title: str
    text: str
    source_name: str
    source_release: str
    source_url: str
    evidence_type: str
    entity_ids: tuple[str, ...] = ()
    reference_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        required = (
            "evidence_id",
            "title",
            "text",
            "source_name",
            "source_release",
            "evidence_type",
        )
        for field_name in required:
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{field_name} must be a non-empty string")

        parsed_url = urlparse(self.source_url)
        if parsed_url.scheme not in {"http", "https"} or not parsed_url.netloc:
            raise ValueError("source_url must be an absolute HTTP(S) URL")

        for field_name in ("entity_ids", "reference_ids"):
            values = getattr(self, field_name)
            if not isinstance(values, tuple):
                raise ValueError(f"{field_name} must be a tuple of strings")
            if any(not isinstance(value, str) or not value.strip() for value in values):
                raise ValueError(f"{field_name} must contain non-empty strings")


@dataclass(frozen=True)
class EvidenceHit:
    """A retrieved document and its transparent lexical match details."""

    document: EvidenceDocument
    score: float
    matched_terms: tuple[str, ...]


@dataclass(frozen=True)
class EvidenceRetrievalReport:
    """Traceable results from a deterministic evidence retrieval query."""

    query: str
    hits: tuple[EvidenceHit, ...]
    total_documents: int
    eligible_documents: int
    research_only: bool = True


def retrieve_evidence(
    query: str,
    documents: tuple[EvidenceDocument, ...],
    *,
    top_k: int = 5,
    entity_ids: tuple[str, ...] = (),
    evidence_types: tuple[str, ...] = (),
) -> EvidenceRetrievalReport:
    """Retrieve records by lexical coverage, retaining source provenance.

    Score = matched unique query tokens / unique query tokens.
    Entity and evidence-type filters are exact, case-sensitive filters.
    This baseline does not infer relationships or assess evidence quality.
    """
    if not isinstance(query, str) or not query.strip():
        raise ValueError("query must be a non-empty string")
    if not isinstance(top_k, int) or isinstance(top_k, bool) or top_k < 1:
        raise ValueError("top_k must be a positive integer")

    query_tokens = _tokens(query)
    if not query_tokens:
        raise ValueError("query must contain at least one searchable token")

    if len({document.evidence_id for document in documents}) != len(documents):
        raise ValueError("evidence_id values must be unique")

    entity_filter = set(entity_ids)
    type_filter = set(evidence_types)
    eligible: list[EvidenceDocument] = []

    for document in documents:
        if entity_filter and not entity_filter.intersection(document.entity_ids):
            continue
        if type_filter and document.evidence_type not in type_filter:
            continue
        eligible.append(document)

    hits: list[EvidenceHit] = []
    for document in eligible:
        searchable_text = " ".join(
            (document.title, document.text, *document.entity_ids, *document.reference_ids)
        )
        matched = tuple(sorted(query_tokens & _tokens(searchable_text)))
        if not matched:
            continue

        hits.append(
            EvidenceHit(
                document=document,
                score=len(matched) / len(query_tokens),
                matched_terms=matched,
            )
        )

    hits.sort(key=lambda hit: (-hit.score, hit.document.evidence_id))

    return EvidenceRetrievalReport(
        query=query.strip(),
        hits=tuple(hits[:top_k]),
        total_documents=len(documents),
        eligible_documents=len(eligible),
    )
