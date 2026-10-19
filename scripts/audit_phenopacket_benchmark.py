"""Reproducible audit script for the Phenopacket Store 0.1.27 benchmark cohort.

Loads the local Phenopacket Store archive, HPO ontology, and disease annotations,
audits cohort integrity and disease index alignment, and writes a machine-readable
audit report to reports/phenopacket_benchmark_audit.json without patient identifiers.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from rare_disease_ai.benchmark_audit import (
    audit_phenopacket_benchmark,
    format_console_summary,
)
from rare_disease_ai.hpo_annotations import read_hpo_annotations
from rare_disease_ai.hpo_ontology import HPOOntology

ROOT = Path(__file__).resolve().parents[1]

DEFAULT_ARCHIVE = ROOT / "data" / "external" / "all_phenopackets-0.1.27.zip"
DEFAULT_ONTOLOGY = ROOT / "data" / "raw" / "hp-base.obo"
DEFAULT_ANNOTATIONS = ROOT / "data" / "external" / "phenotype.hpoa"
DEFAULT_ONTOLOGY_METADATA = ROOT / "data" / "metadata" / "hpo_release.json"
DEFAULT_HPOA_METADATA = ROOT / "data" / "metadata" / "hpoa_release.json"
DEFAULT_OUTPUT = ROOT / "reports" / "phenopacket_benchmark_audit.json"


def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for the benchmark audit."""
    parser = argparse.ArgumentParser(
        description="Audit the Phenopacket Store 0.1.27 benchmark cohort."
    )
    parser.add_argument(
        "--archive",
        type=Path,
        default=DEFAULT_ARCHIVE,
        help="Path to Phenopackets zip archive.",
    )
    parser.add_argument(
        "--ontology",
        type=Path,
        default=DEFAULT_ONTOLOGY,
        help="Path to HPO OBO file.",
    )
    parser.add_argument(
        "--annotations",
        type=Path,
        default=DEFAULT_ANNOTATIONS,
        help="Path to HPOA annotations file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output path for JSON audit report.",
    )
    return parser.parse_args(args)


def main(args: list[str] | None = None) -> None:
    """Run the benchmark audit and output report and console summary."""
    opts = parse_args(args)

    for required_file in (opts.archive, opts.ontology, opts.annotations):
        if not required_file.is_file():
            print(f"Error: Required file not found: {required_file}", file=sys.stderr)
            sys.exit(1)

    print(f"Loading HPO ontology from: {opts.ontology}")
    ontology = HPOOntology.from_obo(opts.ontology)

    print(f"Loading HPO annotations from: {opts.annotations}")
    annotations = read_hpo_annotations(opts.annotations)

    print(f"Auditing Phenopacket archive: {opts.archive} (strict diagnosis matching)...")
    report = audit_phenopacket_benchmark(
        archive_path=opts.archive,
        ontology=ontology,
        disease_annotations=annotations,
        ontology_metadata_path=DEFAULT_ONTOLOGY_METADATA,
        hpoa_metadata_path=DEFAULT_HPOA_METADATA,
    )

    opts.output.parent.mkdir(parents=True, exist_ok=True)
    opts.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    console_summary = format_console_summary(report)
    print("\n" + console_summary)
    print(f"\nAudit report saved to: {opts.output.resolve()}\n")


if __name__ == "__main__":
    main()
