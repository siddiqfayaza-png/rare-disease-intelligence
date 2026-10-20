"""Generate an aggregate validation report for the local GenCC CSV."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from datetime import date
from pathlib import Path

from rare_disease_ai.gencc_import import read_gencc_csv

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "data" / "external" / "gencc-submissions.csv"
METADATA_PATH = ROOT / "data" / "metadata" / "gencc_release.json"
OUTPUT_PATH = ROOT / "reports" / "gencc_validation_report.json"


def main() -> int:
    if not CSV_PATH.is_file():
        raise FileNotFoundError(f"GenCC CSV not found: {CSV_PATH}")
    if not METADATA_PATH.is_file():
        raise FileNotFoundError(f"GenCC metadata not found: {METADATA_PATH}")

    metadata = json.loads(METADATA_PATH.read_text(encoding="utf-8-sig"))
    checksum = hashlib.sha256(CSV_PATH.read_bytes()).hexdigest()
    expected_checksum = str(metadata.get("sha256", "")).lower()

    result = read_gencc_csv(CSV_PATH)
    assertions = result.assertions

    classifications = Counter(
        item.classification_title for item in assertions
    )
    modes_of_inheritance = Counter(
        item.moi_title for item in assertions if item.moi_title
    )

    report = {
        "report_date": date.today().isoformat(),
        "source": metadata.get("source"),
        "source_url": metadata.get("source_url"),
        "retrieved_date": metadata.get("retrieved_date"),
        "input_file": "data/external/gencc-submissions.csv",
        "checksum": {
            "algorithm": "SHA-256",
            "actual": checksum,
            "expected": expected_checksum,
            "matches_metadata": checksum == expected_checksum,
        },
        "import_summary": {
            "rows_read": result.rows_read,
            "assertions_imported": len(assertions),
            "rows_with_errors": len(result.row_errors),
            "unique_submission_ids": len(
                {item.submission_id for item in assertions}
            ),
            "unique_gene_identifiers": len(
                {item.gene_curie for item in assertions if item.gene_curie}
            ),
            "unique_gene_symbols": len(
                {item.gene_symbol for item in assertions if item.gene_symbol}
            ),
            "unique_canonical_disease_identifiers": len(
                {item.disease_curie for item in assertions if item.disease_curie}
            ),
            "unique_original_disease_identifiers": len(
                {
                    item.disease_original_curie
                    for item in assertions
                    if item.disease_original_curie
                }
            ),
            "unique_submitters": len(
                {item.submitter_curie for item in assertions if item.submitter_curie}
            ),
            "assertions_with_publications": sum(
                bool(item.pmids) for item in assertions
            ),
        },
        "classification_counts": dict(sorted(classifications.items())),
        "mode_of_inheritance_counts": dict(
            sorted(modes_of_inheritance.items())
        ),
        "row_errors": list(result.row_errors),
        "interpretation_note": (
            "Counts describe imported GenCC submission records. They do not "
            "establish independent clinical validity. Canonical and original "
            "disease identifiers are counted separately and are not merged "
            "by name."
        ),
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"Rows read: {result.rows_read:,}")
    print(f"Assertions imported: {len(assertions):,}")
    print(f"Rows with errors: {len(result.row_errors):,}")
    print(f"Unique gene symbols: {report['import_summary']['unique_gene_symbols']:,}")
    print(
        "Unique canonical disease IDs: "
        f"{report['import_summary']['unique_canonical_disease_identifiers']:,}"
    )
    print(
        "Unique original disease IDs: "
        f"{report['import_summary']['unique_original_disease_identifiers']:,}"
    )
    print(f"Checksum matches metadata: {checksum == expected_checksum}")
    print(f"Report written to: {OUTPUT_PATH.relative_to(ROOT)}")

    return 0 if checksum == expected_checksum and not result.row_errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
