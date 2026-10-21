"""Build a versioned, provenance-preserving biomedical evidence corpus."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from rare_disease_ai.biomedical_evidence_retrieval import EvidenceDocument
from rare_disease_ai.gencc_import import GenCCAssertion
from rare_disease_ai.hpo_annotations import HPOAnnotation
from rare_disease_ai.mechanism_evidence import MechanismEvidence

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def load_release_metadata(path: str | Path) -> dict[str, Any]:
    """Load a JSON release metadata file from the project root."""
    metadata_path = PROJECT_ROOT / path
    with metadata_path.open("r", encoding="utf-8-sig") as handle:
        metadata = json.load(handle)

    if not isinstance(metadata, dict):
        raise ValueError(f"Release metadata must be a JSON object: {metadata_path}")

    return metadata


def sha256_file(path: str | Path) -> str:
    """Calculate a file's SHA-256 checksum incrementally."""
    digest = hashlib.sha256()

    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def make_evidence_id(source: str, *identifiers: str) -> str:
    """Create a stable ID from source and source-record identifiers."""
    if not source.strip():
        raise ValueError("source must be non-empty")
    if not identifiers or any(not item.strip() for item in identifiers):
        raise ValueError("identifiers must be non-empty strings")

    identity = json.dumps(
        [source, *identifiers],
        ensure_ascii=False,
        separators=(",", ":"),
    )
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()
    return f"{source.casefold()}:{digest}"


def document_to_dict(document: EvidenceDocument) -> dict[str, Any]:
    """Serialize an evidence document deterministically."""
    return {
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


def hpo_annotation_to_document(
    annotation: HPOAnnotation,
    *,
    release: str,
    source_url: str,
    hpo_label: str | None,
) -> EvidenceDocument:
    """Convert one HPOA annotation while preserving all source fields."""
    title = (
        f"{annotation.disease_name}: {hpo_label}"
        if hpo_label
        else f"{annotation.disease_name}: {annotation.hpo_id}"
    )

    details = {
        "disease_id": annotation.database_id,
        "disease_name": annotation.disease_name,
        "hpo_id": annotation.hpo_id,
        "hpo_label": hpo_label,
        "qualifier": annotation.qualifier,
        "evidence": annotation.evidence,
        "frequency": annotation.frequency,
        "onset": annotation.onset,
        "sex": annotation.sex,
        "modifier": annotation.modifier,
        "aspect": annotation.aspect,
        "biocuration": annotation.biocuration,
        "reference": annotation.reference,
    }

    text = "; ".join(
        f"{key}={value}" for key, value in details.items() if value is not None and value != ""
    )

    references = (annotation.reference,) if annotation.reference.strip() else ()

    # Preserve every original HPOA field, including empty values.
    # Field order matches EXPECTED_COLUMNS in hpo_annotations.py.
    identity_fields = (
        annotation.database_id,
        annotation.disease_name,
        annotation.qualifier,
        annotation.hpo_id,
        annotation.reference,
        annotation.evidence,
        annotation.onset,
        annotation.frequency,
        annotation.sex,
        annotation.modifier,
        annotation.aspect,
        annotation.biocuration,
    )

    # Serialize the tuple as one identifier. This avoids ambiguity and
    # preserves empty fields when generating a deterministic evidence ID.
    evidence_id = make_evidence_id(
        "hpoa",
        json.dumps(
            identity_fields,
            ensure_ascii=False,
            separators=(",", ":"),
        ),
    )

    return EvidenceDocument(
        evidence_id=evidence_id,
        title=title,
        text=text,
        source_name="Human Phenotype Ontology disease annotations",
        source_release=release,
        source_url=source_url,
        evidence_type="hpo_disease_phenotype_annotation",
        entity_ids=(annotation.database_id, annotation.hpo_id),
        reference_ids=references,
    )


def gencc_assertion_to_document(
    assertion: GenCCAssertion,
    *,
    release: str,
    source_url: str,
) -> EvidenceDocument:
    """Convert one GenCC assertion into a provenance-preserving document."""
    details = {
        "submission_id": assertion.submission_id,
        "version_number": assertion.version_number,
        "gene_curie": assertion.gene_curie,
        "gene_symbol": assertion.gene_symbol,
        "canonical_disease_id": assertion.disease_curie,
        "canonical_disease_title": assertion.disease_title,
        "original_disease_id": assertion.disease_original_curie,
        "classification_curie": assertion.classification_curie,
        "classification": assertion.classification_title,
        "mode_of_inheritance_curie": assertion.moi_curie,
        "mode_of_inheritance": assertion.moi_title,
        "submitter_curie": assertion.submitter_curie,
        "submitter": assertion.submitter_title,
        "submitted_as_date": assertion.submitted_as_date,
        "public_report_url": assertion.public_report_url,
        "assertion_criteria_url": assertion.assertion_criteria_url,
        "pmids": ";".join(assertion.pmids),
        "interpretation": (
            "GenCC submission; preserve as an assertion, "
            "not an independently verified clinical conclusion."
        ),
    }

    text = "; ".join(
        f"{key}={value}" for key, value in details.items() if value is not None and value != ""
    )

    identifiers = (
        assertion.submission_id,
        str(assertion.version_number),
        assertion.gene_curie,
        assertion.disease_curie,
        assertion.disease_original_curie,
        assertion.classification_curie,
        assertion.moi_curie,
        assertion.submitter_curie,
    )

    evidence_id = make_evidence_id(
        "gencc",
        *(value for value in identifiers if value.strip()),
    )

    return EvidenceDocument(
        evidence_id=evidence_id,
        title=(f"{assertion.gene_symbol} — {assertion.disease_title or assertion.disease_curie}"),
        text=text,
        source_name="Gene Curation Coalition (GenCC)",
        source_release=release,
        source_url=source_url,
        evidence_type="gene_disease_assertion",
        entity_ids=tuple(
            value
            for value in (
                assertion.gene_curie,
                assertion.gene_symbol,
                assertion.disease_curie,
                assertion.disease_original_curie,
            )
            if value and value.strip()
        ),
        reference_ids=assertion.pmids,
    )


def mechanism_evidence_to_document(
    evidence: MechanismEvidence,
) -> EvidenceDocument:
    """Convert one GO mechanism annotation while preserving provenance."""
    gaf_column_names = (
        "database",
        "object_id",
        "gene_symbol",
        "qualifier",
        "go_id",
        "references",
        "evidence_code",
        "with_from",
        "aspect",
        "object_name",
        "synonyms",
        "object_type",
        "taxon",
        "date",
        "assigned_by",
        "annotation_extension",
        "gene_product_form_id",
    )

    details = {
        "gene_symbol": evidence.gene_symbol,
        "annotation_id": evidence.annotation_id,
        "annotation_label": evidence.annotation_label,
        "annotation_type": evidence.annotation_type,
        "evidence_code": evidence.evidence_code,
        "evidence_label": evidence.evidence_label,
        "qualifier": evidence.qualifier,
        "taxon_id": evidence.taxon_id,
        "source_name": evidence.source_name,
        "source_release": evidence.source_release,
        "source_url": evidence.source_url,
        "reference_ids": ";".join(evidence.reference_ids),
        "gaf_source_fields": (
            json.dumps(
                dict(
                    zip(
                        gaf_column_names,
                        evidence.source_record_fields,
                        strict=True,
                    )
                ),
                ensure_ascii=False,
                sort_keys=True,
            )
            if evidence.source_record_fields
            else None
        ),
    }

    text = "; ".join(
        f"{key}={value}" for key, value in details.items() if value is not None and value != ""
    )

    identity = json.dumps(
        [
            evidence.gene_symbol,
            evidence.annotation_id,
            evidence.annotation_type,
            evidence.source_release,
            evidence.evidence_code,
            evidence.reference_ids,
            evidence.qualifier,
            evidence.taxon_id,
            evidence.source_record_fields,
        ],
        ensure_ascii=False,
        separators=(",", ":"),
    )

    evidence_id = make_evidence_id("goa", identity)

    entity_ids = tuple(
        value
        for value in (
            evidence.gene_symbol,
            evidence.annotation_id,
            evidence.taxon_id,
        )
        if value and value.strip()
    )

    return EvidenceDocument(
        evidence_id=evidence_id,
        title=(f"{evidence.gene_symbol}: {evidence.annotation_label} ({evidence.annotation_id})"),
        text=text,
        source_name=evidence.source_name,
        source_release=evidence.source_release,
        source_url=evidence.source_url,
        evidence_type="go_mechanism_annotation",
        entity_ids=entity_ids,
        reference_ids=evidence.reference_ids,
    )
