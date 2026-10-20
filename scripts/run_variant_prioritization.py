"""Run CSV-based variant prioritization and save research reports."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rare_disease_ai.variant_prioritization import (  # noqa: E402
    prioritize_variants,
    read_variant_csv,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate and rank candidate variants from a CSV.")
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input CSV containing normalized variant records.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "reports" / "variant_prioritization",
        help="Directory for output reports.",
    )
    parser.add_argument(
        "--max-allele-frequency",
        type=float,
        default=0.01,
        help="Heuristic rarity threshold, from 0 to 1. Use 'none' to disable.",
    )
    args = parser.parse_args()

    threshold = args.max_allele_frequency
    if not 0 <= threshold <= 1:
        parser.error("--max-allele-frequency must be between 0 and 1.")

    input_path = args.input
    if not input_path.is_absolute():
        input_path = Path.cwd() / input_path

    if not input_path.is_file():
        parser.error(f"Input CSV not found: {input_path}")

    try:
        imported = read_variant_csv(input_path)
    except (OSError, UnicodeError, ValueError, csv.Error) as exc:
        parser.error(f"Could not read input CSV: {exc}")

    ranked = prioritize_variants(
        imported.variants,
        max_allele_frequency=threshold,
    )

    output_dir = args.output_dir
    if not output_dir.is_absolute():
        output_dir = Path.cwd() / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    results_path = output_dir / "variant_ranking.csv"
    with results_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=[
                "rank",
                "chromosome",
                "position",
                "reference",
                "alternate",
                "gene_symbol",
                "consequence",
                "zygosity",
                "allele_frequency",
                "heuristic_score",
                "reasons",
            ],
        )
        writer.writeheader()

        for rank, result in enumerate(ranked, start=1):
            variant = result.variant
            writer.writerow(
                {
                    "rank": rank,
                    "chromosome": variant.chromosome,
                    "position": variant.position,
                    "reference": variant.reference,
                    "alternate": variant.alternate,
                    "gene_symbol": variant.gene_symbol,
                    "consequence": variant.consequence,
                    "zygosity": variant.zygosity,
                    "allele_frequency": (
                        "" if variant.allele_frequency is None else variant.allele_frequency
                    ),
                    "heuristic_score": result.score,
                    "reasons": "; ".join(result.reasons),
                }
            )

    report = {
        "research_only": True,
        "input_file": input_path.name,
        "rows_read": imported.rows_read,
        "valid_variants": len(imported.variants),
        "invalid_rows": len(imported.row_errors),
        "ranked_variants": len(ranked),
        "max_allele_frequency": threshold,
        "validation_errors": list(imported.row_errors),
        "score_interpretation": (
            "Heuristic prioritization only; not a pathogenicity classification "
            "or clinical diagnosis."
        ),
    }

    report_path = output_dir / "variant_validation_report.json"
    report_path.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print(f"Rows read: {imported.rows_read}")
    print(f"Valid variants: {len(imported.variants)}")
    print(f"Invalid rows: {len(imported.row_errors)}")
    print(f"Ranked variants: {len(ranked)}")
    print(f"Ranking CSV: {results_path}")
    print(f"Validation report: {report_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
