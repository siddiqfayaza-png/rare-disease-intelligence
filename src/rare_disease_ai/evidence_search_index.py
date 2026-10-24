"""Disk-backed, provenance-preserving biomedical evidence search index."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sqlite3
from pathlib import Path
from typing import Any

from rare_disease_ai.biomedical_evidence_retrieval import (
    EvidenceDocument,
    EvidenceHit,
    EvidenceRetrievalReport,
    _tokens,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CORPUS = PROJECT_ROOT / "data" / "processed" / "biomedical_evidence.jsonl"
DEFAULT_INDEX = PROJECT_ROOT / "data" / "processed" / "evidence_search.sqlite3"
SCHEMA_VERSION = "1"


def _connect(path: Path) -> sqlite3.Connection:
    """Open a SQLite database with a reasonable integrity/performance balance."""
    connection = sqlite3.connect(path)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    connection.execute("PRAGMA busy_timeout = 30000")
    return connection


def _create_schema(connection: sqlite3.Connection) -> None:
    """Create the document store, exact-token postings, and FTS5 index."""
    connection.executescript(
        """
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );

        CREATE TABLE documents (
            row_id INTEGER PRIMARY KEY,
            evidence_id TEXT NOT NULL,
            evidence_type TEXT NOT NULL,
            document_json TEXT NOT NULL
        );

        CREATE INDEX idx_documents_evidence_id
            ON documents(evidence_id);
        CREATE INDEX idx_documents_evidence_type
            ON documents(evidence_type);

        CREATE TABLE document_entities (
            row_id INTEGER NOT NULL REFERENCES documents(row_id),
            entity_id TEXT NOT NULL,
            PRIMARY KEY (row_id, entity_id)
        );
        CREATE INDEX idx_document_entities_entity
            ON document_entities(entity_id, row_id);

        CREATE TABLE token_postings (
            token TEXT NOT NULL,
            row_id INTEGER NOT NULL REFERENCES documents(row_id),
            PRIMARY KEY (token, row_id)
        );
        CREATE INDEX idx_token_postings_row
            ON token_postings(row_id, token);

        CREATE VIRTUAL TABLE evidence_fts USING fts5(
            title,
            text,
            entity_ids,
            reference_ids,
            tokenize = 'unicode61'
        );
        """
    )


def build_evidence_index(
    corpus_path: str | Path = DEFAULT_CORPUS,
    index_path: str | Path = DEFAULT_INDEX,
    *,
    batch_size: int = 2000,
) -> dict[str, Any]:
    """Stream a JSONL corpus into a new index and atomically publish it.

    Identical normalized payloads with the same evidence ID are indexed once.
    Conflicting payloads with the same evidence ID fail the build.
    The destination database is replaced only after a successful build.
    """
    source = Path(corpus_path)
    destination = Path(index_path)

    if not source.is_file():
        raise FileNotFoundError(f"Evidence corpus not found: {source}")

    if not isinstance(batch_size, int) or isinstance(batch_size, bool) or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".building")

    if temporary.exists():
        temporary.unlink()

    connection: sqlite3.Connection | None = None
    digest = hashlib.sha256()
    record_count = 0
    corpus_record_count = 0
    duplicate_id_count = 0

    try:
        connection = _connect(temporary)
        connection.execute("PRAGMA journal_mode = DELETE")
        connection.execute("PRAGMA synchronous = NORMAL")
        _create_schema(connection)

        # Track IDs and payload hashes in SQLite rather than retaining the
        # complete corpus in Python memory.
        connection.execute(
            """
            CREATE TEMP TABLE seen_evidence_ids (
                evidence_id TEXT PRIMARY KEY,
                payload_sha256 TEXT NOT NULL,
                first_line_number INTEGER NOT NULL
            )
            """
        )

        with source.open("rb") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                # Hash original bytes, including line endings, so the checksum
                # represents the exact input file, including duplicate rows.
                digest.update(raw_line)

                if not raw_line.strip():
                    continue

                try:
                    payload = json.loads(raw_line.decode("utf-8"))

                    if not isinstance(payload, dict):
                        raise ValueError("record must be a JSON object")

                    payload["entity_ids"] = tuple(payload.get("entity_ids", ()))
                    payload["reference_ids"] = tuple(payload.get("reference_ids", ()))

                    document = EvidenceDocument(**payload)

                except (
                    UnicodeDecodeError,
                    json.JSONDecodeError,
                    TypeError,
                    ValueError,
                ) as exc:
                    raise ValueError(
                        f"Invalid evidence record on line {line_number}: {exc}"
                    ) from exc

                corpus_record_count += 1

                # Serialize all normalized document fields consistently.
                serialized = json.dumps(
                    {
                        "evidence_id": document.evidence_id,
                        "title": document.title,
                        "text": document.text,
                        "source_name": document.source_name,
                        "source_release": document.source_release,
                        "source_url": document.source_url,
                        "evidence_type": document.evidence_type,
                        "entity_ids": list(document.entity_ids),
                        "reference_ids": list(document.reference_ids),
                    },
                    ensure_ascii=False,
                    separators=(",", ":"),
                )

                payload_sha256 = hashlib.sha256(serialized.encode("utf-8")).hexdigest()

                existing = connection.execute(
                    """
                    SELECT payload_sha256, first_line_number
                    FROM seen_evidence_ids
                    WHERE evidence_id = ?
                    """,
                    (document.evidence_id,),
                ).fetchone()

                if existing is not None:
                    if existing["payload_sha256"] == payload_sha256:
                        # Exact normalized duplicate: preserve the original
                        # corpus, but index this evidence record only once.
                        duplicate_id_count += 1
                        continue

                    # Never silently discard or choose between conflicting
                    # records. Raising here prevents publication of the build.
                    raise ValueError(
                        "Conflicting payloads for evidence ID "
                        f"{document.evidence_id!r} on line {line_number}; "
                        f"first seen on line {existing['first_line_number']}"
                    )

                connection.execute(
                    """
                    INSERT INTO seen_evidence_ids (
                        evidence_id, payload_sha256, first_line_number
                    )
                    VALUES (?, ?, ?)
                    """,
                    (
                        document.evidence_id,
                        payload_sha256,
                        line_number,
                    ),
                )

                row_id = record_count + 1

                connection.execute(
                    """
                    INSERT INTO documents
                        (row_id, evidence_id, evidence_type, document_json)
                    VALUES (?, ?, ?, ?)
                    """,
                    (
                        row_id,
                        document.evidence_id,
                        document.evidence_type,
                        serialized,
                    ),
                )

                connection.executemany(
                    """
                    INSERT OR IGNORE INTO document_entities(row_id, entity_id)
                    VALUES (?, ?)
                    """,
                    ((row_id, entity) for entity in document.entity_ids),
                )

                searchable_text = " ".join(
                    (
                        document.title,
                        document.text,
                        *document.entity_ids,
                        *document.reference_ids,
                    )
                )

                connection.executemany(
                    """
                    INSERT INTO token_postings(token, row_id)
                    VALUES (?, ?)
                    """,
                    ((token, row_id) for token in _tokens(searchable_text)),
                )

                connection.execute(
                    """
                    INSERT INTO evidence_fts(
                        rowid, title, text, entity_ids, reference_ids
                    )
                    VALUES (?, ?, ?, ?, ?)
                    """,
                    (
                        row_id,
                        document.title,
                        document.text,
                        " ".join(document.entity_ids),
                        " ".join(document.reference_ids),
                    ),
                )

                record_count += 1

                if record_count % batch_size == 0:
                    connection.commit()

        corpus_checksum = digest.hexdigest()

        metadata = {
            "schema_version": SCHEMA_VERSION,
            "corpus_path": str(source.resolve()),
            "corpus_sha256": corpus_checksum,
            "document_count": str(record_count),
            "corpus_record_count": str(corpus_record_count),
            # Backward-compatible key: now counts skipped duplicate
            # occurrences rather than duplicate rows retained in the index.
            "duplicate_evidence_id_count": str(duplicate_id_count),
            "duplicate_records_skipped": str(duplicate_id_count),
            # Successful builds cannot contain conflicting IDs.
            "conflicting_evidence_id_count": "0",
        }

        connection.executemany(
            "INSERT INTO metadata(key, value) VALUES (?, ?)",
            metadata.items(),
        )
        connection.commit()
        connection.execute("PRAGMA optimize")
        connection.close()
        connection = None

        # Publish only a completely built database. If any prior step fails,
        # the existing destination index remains untouched.
        os.replace(temporary, destination)

        return {
            "build_status": "success",
            "index_path": str(destination.resolve()),
            "document_count": record_count,
            "corpus_record_count": corpus_record_count,
            "duplicate_evidence_id_count": duplicate_id_count,
            "duplicate_records_skipped": duplicate_id_count,
            "conflicting_evidence_id_count": 0,
            "corpus_sha256": corpus_checksum,
            "schema_version": SCHEMA_VERSION,
        }

    except Exception:
        if connection is not None:
            connection.close()

        if temporary.exists():
            temporary.unlink()

        raise


class EvidenceSearchIndex:
    """Search an existing index using exact query-token coverage ranking."""

    def __init__(self, index_path: str | Path = DEFAULT_INDEX) -> None:
        self.index_path = Path(index_path)

        if not self.index_path.is_file():
            raise FileNotFoundError(f"Evidence index not found: {self.index_path}")

        with _connect(self.index_path) as connection:
            metadata = dict(connection.execute("SELECT key, value FROM metadata").fetchall())

            if metadata.get("schema_version") != SCHEMA_VERSION:
                raise ValueError("Unsupported or missing evidence index schema version")

            self.metadata = metadata

    def search(
        self,
        query: str,
        *,
        top_k: int = 5,
        entity_ids: tuple[str, ...] = (),
        evidence_types: tuple[str, ...] = (),
    ) -> EvidenceRetrievalReport:
        """Search while preserving baseline scores and exact filters."""
        if not isinstance(query, str) or not query.strip():
            raise ValueError("query must be a non-empty string")

        if not isinstance(top_k, int) or isinstance(top_k, bool) or top_k < 1:
            raise ValueError("top_k must be a positive integer")

        query_tokens = _tokens(query)

        if not query_tokens:
            raise ValueError("query must contain at least one searchable token")

        entity_filter = tuple(sorted(set(entity_ids)))
        type_filter = tuple(sorted(set(evidence_types)))

        with _connect(self.index_path) as connection:
            total_documents = connection.execute("SELECT COUNT(*) FROM documents").fetchone()[0]

            conditions: list[str] = []
            parameters: list[Any] = []

            if entity_filter:
                placeholders = ",".join("?" for _ in entity_filter)
                conditions.append(
                    f"""
                    EXISTS (
                        SELECT 1
                        FROM document_entities de
                        WHERE de.row_id = d.row_id
                          AND de.entity_id IN ({placeholders})
                    )
                    """
                )
                parameters.extend(entity_filter)

            if type_filter:
                placeholders = ",".join("?" for _ in type_filter)
                conditions.append(f"d.evidence_type IN ({placeholders})")
                parameters.extend(type_filter)

            filter_sql = "WHERE " + " AND ".join(conditions) if conditions else ""

            eligible_documents = connection.execute(
                f"""
                SELECT COUNT(*)
                FROM documents d
                {filter_sql}
                """,
                parameters,
            ).fetchone()[0]

            token_placeholders = ",".join("?" for _ in query_tokens)

            sql = f"""
                SELECT
                    d.row_id,
                    d.document_json,
                    COUNT(DISTINCT p.token) AS matched_count
                FROM token_postings p
                JOIN documents d ON d.row_id = p.row_id
                WHERE p.token IN ({token_placeholders})
                  {"AND " + " AND ".join(conditions) if conditions else ""}
                GROUP BY d.row_id
                ORDER BY
                    matched_count DESC,
                    d.evidence_id ASC,
                    d.row_id ASC
                LIMIT ?
            """

            rows = connection.execute(
                sql,
                [*sorted(query_tokens), *parameters, top_k],
            ).fetchall()

            hits: list[EvidenceHit] = []

            # Use a temporary query-token table instead of JSON1's
            # json_each(), making matched-term retrieval independent of
            # the SQLite JSON1 extension.
            connection.execute(
                """
                CREATE TEMP TABLE query_tokens (
                    token TEXT PRIMARY KEY
                )
                """
            )
            connection.executemany(
                "INSERT INTO query_tokens(token) VALUES (?)",
                ((token,) for token in sorted(query_tokens)),
            )

            for row in rows:
                matched_rows = connection.execute(
                    """
                    SELECT p.token
                    FROM token_postings p
                    JOIN query_tokens q ON q.token = p.token
                    WHERE p.row_id = ?
                    ORDER BY p.token
                    """,
                    (row["row_id"],),
                ).fetchall()

                document_data = json.loads(row["document_json"])
                document_data["entity_ids"] = tuple(document_data["entity_ids"])
                document_data["reference_ids"] = tuple(document_data["reference_ids"])

                document = EvidenceDocument(**document_data)
                matched_terms = tuple(item["token"] for item in matched_rows)

                hits.append(
                    EvidenceHit(
                        document=document,
                        score=len(matched_terms) / len(query_tokens),
                        matched_terms=matched_terms,
                    )
                )

        return EvidenceRetrievalReport(
            query=query.strip(),
            hits=tuple(hits),
            total_documents=total_documents,
            eligible_documents=eligible_documents,
        )


def main() -> None:
    """Build the index from the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    build_parser = subparsers.add_parser(
        "build",
        help="Build the evidence index",
    )
    build_parser.add_argument(
        "--corpus",
        type=Path,
        default=DEFAULT_CORPUS,
    )
    build_parser.add_argument(
        "--index",
        type=Path,
        default=DEFAULT_INDEX,
    )
    build_parser.add_argument(
        "--batch-size",
        type=int,
        default=2000,
    )

    args = parser.parse_args()

    if args.command == "build":
        result = build_evidence_index(
            args.corpus,
            args.index,
            batch_size=args.batch_size,
        )
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
