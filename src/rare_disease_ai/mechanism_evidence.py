"""Typed, provenance-preserving biological mechanism evidence models."""

from dataclasses import dataclass
from urllib.parse import urlparse

SUPPORTED_ANNOTATION_TYPES = frozenset(
    {
        "pathway",
        "molecular_function",
        "biological_process",
        "cellular_component",
        "molecular_mechanism",
    }
)


@dataclass(frozen=True)
class MechanismEvidence:
    """One source-backed biological annotation for a gene."""

    gene_symbol: str
    annotation_id: str
    annotation_label: str
    annotation_type: str
    source_name: str
    source_release: str
    source_url: str
    evidence_code: str | None = None
    evidence_label: str | None = None
    reference_ids: tuple[str, ...] = ()
    qualifier: str | None = None
    taxon_id: str | None = None

    def __post_init__(self) -> None:
        required_fields = {
            "gene_symbol": self.gene_symbol,
            "annotation_id": self.annotation_id,
            "annotation_label": self.annotation_label,
            "source_name": self.source_name,
            "source_release": self.source_release,
            "source_url": self.source_url,
        }
        for field_name, value in required_fields.items():
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{field_name} must be a non-empty string")

        if self.annotation_type not in SUPPORTED_ANNOTATION_TYPES:
            raise ValueError(f"Unsupported annotation_type: {self.annotation_type}")

        parsed_url = urlparse(self.source_url)
        if parsed_url.scheme not in {"http", "https"} or not parsed_url.netloc:
            raise ValueError("source_url must be an absolute HTTP(S) URL")

        if not isinstance(self.reference_ids, tuple):
            raise ValueError("reference_ids must be a tuple of strings")
        if any(
            not isinstance(reference, str) or not reference.strip()
            for reference in self.reference_ids
        ):
            raise ValueError("reference_ids must contain non-empty strings")

        for field_name in (
            "evidence_code",
            "evidence_label",
            "qualifier",
            "taxon_id",
        ):
            value = getattr(self, field_name)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ValueError(f"{field_name} must be None or a non-empty string")


def validate_mechanism_evidence(
    records: tuple[MechanismEvidence, ...],
) -> tuple[str, ...]:
    """Report exact duplicate assertions without conflating distinct evidence."""

    issues: list[str] = []
    seen: set[tuple[object, ...]] = set()

    for index, record in enumerate(records):
        if not isinstance(record, MechanismEvidence):
            issues.append(f"Record {index}: unsupported record type")
            continue

        identity = (
            record.gene_symbol,
            record.annotation_id,
            record.annotation_type,
            record.source_name,
            record.source_release,
            record.evidence_code,
            record.reference_ids,
            record.qualifier,
            record.taxon_id,
        )

        if identity in seen:
            issues.append(f"Record {index}: duplicate evidence assertion")
        seen.add(identity)

    return tuple(issues)
