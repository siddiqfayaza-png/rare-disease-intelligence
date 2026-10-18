"""Load and query terms from an HPO OBO ontology file."""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class HPOTerm:
    """A phenotype term parsed from an HPO ontology."""

    hpo_id: str
    name: str
    synonyms: tuple[str, ...] = ()
    parents: tuple[str, ...] = ()
    is_obsolete: bool = False
    replaced_by: tuple[str, ...] = ()


@dataclass(frozen=True)
class NormalizationResult:
    """Outcome of resolving a label or identifier to HPO terms."""

    query: str
    matched_terms: tuple[HPOTerm, ...]

    @property
    def is_resolved(self) -> bool:
        """Return whether exactly one active term was matched."""
        return len(self.matched_terms) == 1 and not self.matched_terms[0].is_obsolete

    @property
    def is_ambiguous(self) -> bool:
        """Return whether multiple terms were matched."""
        return len(self.matched_terms) > 1

    @property
    def term(self) -> HPOTerm | None:
        """Return the unique active match, or None if unresolved."""
        if self.is_resolved:
            return self.matched_terms[0]
        return None


def _unescape_obo(value: str) -> str:
    """Decode common quoted-string escapes in OBO values."""
    return value.replace(r"\"", '"').replace(r"\\", "\\")


def _parse_synonym(line: str) -> str | None:
    """Extract a quoted synonym label from an OBO synonym stanza."""
    match = re.match(r'^synonym:\s+"((?:\\.|[^"\\])*)"', line)
    if match is None:
        return None
    return _unescape_obo(match.group(1))


def _parse_term(stanza: list[str]) -> HPOTerm | None:
    """Parse one OBO [Term] stanza, ignoring non-HPO identifiers."""
    fields: dict[str, list[str]] = defaultdict(list)

    for line in stanza:
        if not line or line.startswith("!"):
            continue
        key, separator, value = line.partition(":")
        if separator:
            fields[key].append(value.strip())

    ids = fields.get("id", [])
    names = fields.get("name", [])
    if not ids or not names or not ids[0].startswith("HP:"):
        return None

    synonyms = tuple(
        synonym
        for line in fields.get("synonym", [])
        if (synonym := _parse_synonym(f"synonym: {line}")) is not None
    )
    parents = tuple(value.split()[0] for value in fields.get("is_a", []) if value.split())
    replaced_by = tuple(
        value.split()[0] for value in fields.get("replaced_by", []) if value.split()
    )

    return HPOTerm(
        hpo_id=ids[0],
        name=_unescape_obo(names[0]),
        synonyms=synonyms,
        parents=parents,
        is_obsolete=any(value.lower() == "true" for value in fields.get("is_obsolete", [])),
        replaced_by=replaced_by,
    )


def _normalize_label(value: str) -> str:
    """Normalize a label for case-insensitive exact matching."""
    return " ".join(value.casefold().split())


class HPOOntology:
    """In-memory HPO index with deterministic exact-match normalization."""

    def __init__(self, terms: dict[str, HPOTerm]) -> None:
        self.terms = dict(terms)
        self._label_index: dict[str, set[str]] = defaultdict(set)

        for term in self.terms.values():
            self._label_index[_normalize_label(term.name)].add(term.hpo_id)
            for synonym in term.synonyms:
                self._label_index[_normalize_label(synonym)].add(term.hpo_id)

    @classmethod
    def from_obo(cls, path: str | Path) -> HPOOntology:
        """Load HPO terms from an OBO file."""
        obo_path = Path(path)
        if not obo_path.is_file():
            raise FileNotFoundError(f"HPO ontology file not found: {obo_path}")

        terms: dict[str, HPOTerm] = {}
        stanza: list[str] = []
        stanza_type: str | None = None

        def flush() -> None:
            if stanza_type != "Term":
                return
            term = _parse_term(stanza)
            if term is not None:
                if term.hpo_id in terms:
                    raise ValueError(f"Duplicate HPO identifier: {term.hpo_id}")
                terms[term.hpo_id] = term

        with obo_path.open(encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if line.startswith("[") and line.endswith("]"):
                    flush()
                    stanza = []
                    stanza_type = line[1:-1]
                elif stanza_type == "Term":
                    stanza.append(line)
            flush()

        if not terms:
            raise ValueError(f"No HPO terms were parsed from: {obo_path}")

        return cls(terms)

    def get_term(self, hpo_id: str) -> HPOTerm | None:
        """Look up a term by HPO identifier."""
        return self.terms.get(hpo_id.strip())

    def normalize(self, query: str) -> NormalizationResult:
        """Resolve an exact HPO ID, preferred label, or synonym."""
        if not isinstance(query, str) or not query.strip():
            raise ValueError("Phenotype query must be a non-empty string.")

        cleaned = query.strip()
        if cleaned in self.terms:
            matches = (self.terms[cleaned],)
        else:
            ids = sorted(self._label_index.get(_normalize_label(cleaned), set()))
            matches = tuple(self.terms[hpo_id] for hpo_id in ids)

        return NormalizationResult(query=cleaned, matched_terms=matches)
