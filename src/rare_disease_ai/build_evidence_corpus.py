"""Build a checksum-verified, provenance-preserving biomedical evidence corpus."""

from __future__ import annotations

import gzip
import json
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

from rare_disease_ai.evidence_corpus import (
    PROJECT_ROOT,
    document_to_dict,
    gencc_assertion_to_document,
    hpo_annotation_to_document,
    load_release_metadata,
    mechanism_evidence_to_document,
    sha256_file,
)
from rare_disease_ai.gencc_import import read_gencc_csv
from rare_disease_ai.go_annotation_import import iter_gaf_annotations
from rare_disease_ai.go_ontology import read_go_term_labels
from rare_disease_ai.hpo_annotations import read_hpo_annotations
from rare_disease_ai.hpo_ontology import HPOOntology

HPO_ONTOLOGY_PATH = "data/raw/hp-base.obo"
HPOA_PATH = "data/external/phenotype.hpoa"
GENCC_PATH = "data/external/gencc-submissions.csv"
GO_ONTOLOGY_PATH = "data/external/go-basic-2026-08-05.obo"
GOA_PATH = "data/external/HUMAN-uniprot-2026-08-05.gaf.gz"

HPO_METADATA_PATH = "data/metadata/hpo_release.json"
HPOA_METADATA_PATH = "data/metadata/hpoa_release.json"
GENCC_METADATA_PATH = "data/metadata/gencc_release.json"
GO_METADATA_PATH = "data/metadata/go_release.json"

OUTPUT_PATH = "data/processed/biomedical_evidence.jsonl"
REPORT_PATH = "reports/evidence_corpus_build.json"


def _resolve(relative_path: str) -> Path:
    """Resolve a project-relative path."""
    return PROJECT_ROOT / relative_path


def _verify_checksum(
    relative_path: str,
    metadata: dict[str, Any],
    checksum_key: str = "sha256",
) -> str:
    """Verify a local input against its recorded SHA-256 checksum."""
    expected = metadata.get(checksum_key)
    if not isinstance(expected, str) or not expected.strip():
        raise ValueError(f"Missing {checksum_key!r} in metadata for {relative_path}")

    actual = sha256_file(_resolve(relative_path))
    if actual.casefold() != expected.casefold():
        raise ValueError(f"SHA-256 mismatch for {relative_path}: expected {expected}, got {actual}")

    return actual


def _load_sources() -> dict[str, Any]:
    """Load source metadata and verify all five source files."""
    metadata = {
        "hpo": load_release_metadata(HPO_METADATA_PATH),
        "hpoa": load_release_metadata(HPOA_METADATA_PATH),
        "gencc": load_release_metadata(GENCC_METADATA_PATH),
        "go": load_release_metadata(GO_METADATA_PATH),
    }

    verified = {
        "hpo_ontology": _verify_checksum(HPO_ONTOLOGY_PATH, metadata["hpo"]),
        "hpoa_annotations": _verify_checksum(HPOA_PATH, metadata["hpoa"]),
        "gencc_assertions": _verify_checksum(GENCC_PATH, metadata["gencc"]),
        "go_ontology": _verify_checksum(
            GO_ONTOLOGY_PATH,
            metadata["go"],
            "ontology_sha256",
        ),
        "go_annotations": _verify_checksum(
            GOA_PATH,
            metadata["go"],
            "annotation_sha256",
        ),
    }

    return {
        "metadata": metadata,
        "verified_sha256": verified,
    }


def _diagnose_gaf(
    path: Path,
    term_labels: dict[str, str],
) -> dict[str, int]:
    """Count GAF filtering and rejection reasons independently."""
    counts: Counter[str] = Counter()
    opener = gzip.open if path.suffix == ".gz" else open

    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("!"):
                continue

            counts["rows_read"] += 1
            fields = line.rstrip("\r\n").split("\t")

            if len(fields) != 17:
                counts["malformed_column_count"] += 1
                continue

            required = {
                "database": fields[0],
                "object_id": fields[1],
                "gene_symbol": fields[2],
                "GO ID": fields[4],
                "reference": fields[5],
                "evidence_code": fields[6],
                "aspect": fields[8],
                "taxon": fields[12],
                "assigned_by": fields[14],
            }

            if any(not value.strip() for value in required.values()):
                counts["missing_required_fields"] += 1
                continue

            primary_taxon = fields[12].split("|")[0].strip()
            if primary_taxon != "taxon:9606":
                counts["filtered_by_taxon"] += 1
                continue

            if fields[8] not in {"F", "P", "C"}:
                counts["unsupported_aspect"] += 1
                continue

            if not term_labels.get(fields[4], "").strip():
                counts["missing_go_label"] += 1
                continue

            counts["eligible_for_emission"] += 1

    return dict(sorted(counts.items()))


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    """Write JSON to a temporary sibling, then atomically replace the target."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            json.dump(
                payload,
                handle,
                ensure_ascii=False,
                indent=2,
                sort_keys=True,
            )
            handle.write("\n")

        os.replace(temporary_path, path)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def build_corpus() -> dict[str, Any]:
    """Build the corpus and return its reproducibility and quality report."""
    sources = _load_sources()
    metadata = sources["metadata"]
    verified = sources["verified_sha256"]

    hpo_ontology = HPOOntology.from_obo(_resolve(HPO_ONTOLOGY_PATH))
    hpo_annotations = read_hpo_annotations(_resolve(HPOA_PATH))
    gencc_result = read_gencc_csv(_resolve(GENCC_PATH))
    go_labels = read_go_term_labels(_resolve(GO_ONTOLOGY_PATH))

    go_diagnostics = _diagnose_gaf(_resolve(GOA_PATH), go_labels)

    hpoa_metadata = metadata["hpoa"]
    gencc_metadata = metadata["gencc"]
    go_metadata = metadata["go"]

    hpoa_release = str(hpoa_metadata.get("release", "unknown"))
    gencc_release = str(
        gencc_metadata.get("release") or gencc_metadata.get("retrieved_date") or "unknown"
    )

    # Keep the annotation archive date distinct from the ontology version.
    go_release = str(
        go_metadata.get("annotation_archive_date")
        or go_metadata.get("annotation_generated_at")
        or "unknown"
    )
    go_source_url = str(go_metadata["annotation_download_url"])

    output_path = _resolve(OUTPUT_PATH)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path = _resolve(REPORT_PATH)

    counts: Counter[str] = Counter()
    seen_ids: set[str] = set()
    duplicate_id_count = 0
    total_documents = 0
    temporary_path: Path | None = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)

            def emit(document: Any) -> None:
                nonlocal total_documents, duplicate_id_count

                payload = document_to_dict(document)
                evidence_id = payload["evidence_id"]

                if evidence_id in seen_ids:
                    duplicate_id_count += 1
                else:
                    seen_ids.add(evidence_id)

                handle.write(
                    json.dumps(
                        payload,
                        ensure_ascii=False,
                        sort_keys=True,
                        separators=(",", ":"),
                    )
                    + "\n"
                )
                counts[payload["evidence_type"]] += 1
                total_documents += 1

            for annotation in hpo_annotations:
                term = hpo_ontology.get_term(annotation.hpo_id)
                emit(
                    hpo_annotation_to_document(
                        annotation,
                        release=hpoa_release,
                        source_url=str(hpoa_metadata["source_url"]),
                        hpo_label=term.name if term else None,
                    )
                )

            for assertion in gencc_result.assertions:
                emit(
                    gencc_assertion_to_document(
                        assertion,
                        release=gencc_release,
                        source_url=str(gencc_metadata["source_url"]),
                    )
                )

            for evidence in iter_gaf_annotations(
                _resolve(GOA_PATH),
                source_release=go_release,
                source_url=go_source_url,
                term_labels=go_labels,
            ):
                emit(mechanism_evidence_to_document(evidence))

        if counts["go_mechanism_annotation"] != go_diagnostics.get("eligible_for_emission", 0):
            raise ValueError(
                "GO diagnostic/emission count mismatch: "
                f"eligible={go_diagnostics.get('eligible_for_emission', 0)}, "
                f"emitted={counts['go_mechanism_annotation']}"
            )

        os.replace(temporary_path, output_path)
        temporary_path = None
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()

    report: dict[str, Any] = {
        "schema_version": "1.0",
        "research_only": True,
        "clinical_validation_status": "not_clinically_validated",
        "build_status": "success",
        "output": {
            "path": OUTPUT_PATH,
            "sha256": sha256_file(output_path),
            "documents_written": total_documents,
            "counts_by_evidence_type": dict(sorted(counts.items())),
            "duplicate_evidence_ids": duplicate_id_count,
            "duplicate_records_retained": True,
        },
        "sources": {
            "metadata": metadata,
            "verified_sha256": verified,
            "hpo": {
                "path": HPO_ONTOLOGY_PATH,
                "terms_loaded": len(hpo_ontology.terms),
            },
            "hpoa": {
                "path": HPOA_PATH,
                "rows_parsed": len(hpo_annotations),
                "hpo_ids_not_in_ontology": sum(
                    hpo_ontology.get_term(item.hpo_id) is None for item in hpo_annotations
                ),
            },
            "gencc": {
                "path": GENCC_PATH,
                "rows_read": gencc_result.rows_read,
                "assertions_emitted": len(gencc_result.assertions),
                "row_errors": len(gencc_result.row_errors),
                "row_error_examples": list(gencc_result.row_errors[:10]),
            },
            "go_ontology": {
                "path": GO_ONTOLOGY_PATH,
                "active_labels_loaded": len(go_labels),
            },
            "go_annotations": {
                "path": GOA_PATH,
                "diagnostics": go_diagnostics,
                "documents_emitted": counts["go_mechanism_annotation"],
            },
        },
        "limitations": [
            "Research prototype; not for clinical decision-making.",
            "Source assertions and annotations are not independently clinically validated.",
            (
                "GO annotations preserve source qualifiers, including NOT; "
                "retrieval is not a clinical interpretation."
            ),
            "Duplicate evidence IDs are reported and duplicate records are retained.",
            "Source checksums verify file identity, not scientific correctness.",
        ],
    }

    _write_json_atomic(report_path, report)
    return report


def main() -> None:
    """Run the corpus builder from the command line."""
    report = build_corpus()
    print(
        json.dumps(
            {
                "build_status": report["build_status"],
                "documents_written": report["output"]["documents_written"],
                "output_path": report["output"]["path"],
                "output_sha256": report["output"]["sha256"],
                "report_path": REPORT_PATH,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
