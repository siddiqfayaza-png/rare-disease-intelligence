"""Import and validate local GenCC gene-disease assertion exports.

Research prototype only. GenCC submissions are evidence records, not
independently verified clinical conclusions.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

REQUIRED_COLUMNS = (
    "sgc_id",
    "version_number",
    "gene_curie",
    "gene_symbol",
    "disease_curie",
    "disease_title",
    "disease_original_curie",
    "classification_curie",
    "classification_title",
    "moi_curie",
    "moi_title",
    "submitter_curie",
    "submitter_title",
    "submitted_as_date",
    "submitted_as_public_report_url",
    "submitted_as_pmids",
    "submitted_as_assertion_criteria_url",
)


@dataclass(frozen=True)
class GenCCAssertion:
    submission_id: str
    version_number: int
    gene_curie: str
    gene_symbol: str
    disease_curie: str
    disease_title: str
    disease_original_curie: str
    classification_curie: str
    classification_title: str
    moi_curie: str
    moi_title: str
    submitter_curie: str
    submitter_title: str
    submitted_as_date: str
    public_report_url: str
    pmids: tuple[str, ...]
    assertion_criteria_url: str


@dataclass(frozen=True)
class GenCCImportResult:
    assertions: tuple[GenCCAssertion, ...]
    rows_read: int
    row_errors: tuple[str, ...]


def _split_pmids(value: str) -> tuple[str, ...]:
    """Parse semicolon/comma-separated PMID values without inventing IDs."""
    normalized = value.replace(",", ";")
    return tuple(item.strip() for item in normalized.split(";") if item.strip())


def read_gencc_csv(path: str | Path) -> GenCCImportResult:
    """Read GenCC's new-format CSV, reporting malformed rows individually."""
    assertions: list[GenCCAssertion] = []
    errors: list[str] = []
    rows_read = 0

    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise ValueError("GenCC CSV is missing its header.")

        missing = sorted(set(REQUIRED_COLUMNS) - set(reader.fieldnames))
        if missing:
            raise ValueError("GenCC CSV is missing required columns: " + ", ".join(missing))

        for line_number, row in enumerate(reader, start=2):
            rows_read += 1
            try:
                values = {key: (row.get(key) or "").strip() for key in REQUIRED_COLUMNS}

                required_values = (
                    "sgc_id",
                    "version_number",
                    "gene_curie",
                    "gene_symbol",
                    "disease_curie",
                    "classification_title",
                    "moi_curie",
                    "submitter_curie",
                )
                missing_values = [key for key in required_values if not values[key]]
                if missing_values:
                    raise ValueError("missing values: " + ", ".join(missing_values))

                version = int(values["version_number"])
                if version < 1:
                    raise ValueError("version_number must be positive")

                assertions.append(
                    GenCCAssertion(
                        submission_id=values["sgc_id"],
                        version_number=version,
                        gene_curie=values["gene_curie"],
                        gene_symbol=values["gene_symbol"],
                        disease_curie=values["disease_curie"],
                        disease_title=values["disease_title"],
                        disease_original_curie=values["disease_original_curie"],
                        classification_curie=values["classification_curie"],
                        classification_title=values["classification_title"],
                        moi_curie=values["moi_curie"],
                        moi_title=values["moi_title"],
                        submitter_curie=values["submitter_curie"],
                        submitter_title=values["submitter_title"],
                        submitted_as_date=values["submitted_as_date"],
                        public_report_url=values["submitted_as_public_report_url"],
                        pmids=_split_pmids(values["submitted_as_pmids"]),
                        assertion_criteria_url=values["submitted_as_assertion_criteria_url"],
                    )
                )
            except (TypeError, ValueError) as exc:
                errors.append(f"line {line_number}: {exc}")

    return GenCCImportResult(
        assertions=tuple(assertions),
        rows_read=rows_read,
        row_errors=tuple(errors),
    )
