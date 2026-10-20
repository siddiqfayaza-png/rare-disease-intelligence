"""Audit exact GenCC disease-ID matches against the HPOA disease index."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from rare_disease_ai.disease_phenotypes import build_disease_phenotype_index
from rare_disease_ai.gencc_evidence_adapter import adapt_gencc_assertions
from rare_disease_ai.gencc_import import read_gencc_csv
from rare_disease_ai.hpo_annotations import read_hpo_annotations

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    gencc_path = ROOT / "data" / "external" / "gencc-submissions.csv"
    hpoa_path = ROOT / "data" / "external" / "phenotype.hpoa"
    output_path = ROOT / "reports" / "gencc_integration_audit.json"

    for path in (gencc_path, hpoa_path):
        if not path.is_file():
            raise FileNotFoundError(f"Required dataset not found: {path}")

    gencc_result = read_gencc_csv(gencc_path)
    hpo_annotations = read_hpo_annotations(hpoa_path)
    disease_index = build_disease_phenotype_index(hpo_annotations)
    indexed_ids = set(disease_index)

    adapter_result = adapt_gencc_assertions(
        gencc_result.assertions,
        indexed_ids,
    )

    # Count distinct disease-gene pairs while preserving all assertion records.
    assertions_by_pair: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for evidence in adapter_result.evidence:
        assertions_by_pair[(evidence.disease_id, evidence.gene_symbol.strip().upper())].append(
            {
                "submission_id": evidence.evidence_id,
                "classification": evidence.classification,
                "evidence_category": evidence.evidence_category,
                "submitter": evidence.submitter,
                "canonical_disease_id": evidence.canonical_disease_id,
                "original_disease_id": evidence.original_disease_id,
            }
        )

    classification_conflicts = []
    for (disease_id, gene_symbol), records in sorted(assertions_by_pair.items()):
        classifications = sorted({record["classification"] or "Unspecified" for record in records})
        if len(classifications) > 1:
            classification_conflicts.append(
                {
                    "disease_id": disease_id,
                    "gene_symbol": gene_symbol,
                    "assertion_count": len(records),
                    "classifications": classifications,
                }
            )

    report = {
        "research_only": True,
        "hpoa_disease_ids": len(indexed_ids),
        "gencc_rows_read": gencc_result.rows_read,
        "gencc_import_errors": len(gencc_result.row_errors),
        "gencc_assertions_imported": len(gencc_result.assertions),
        "matched_assertions": adapter_result.matched_assertions,
        "unmatched_assertions": adapter_result.unmatched_assertions,
        "matched_by_original_id": adapter_result.matched_by_original_id,
        "matched_by_canonical_id": adapter_result.matched_by_canonical_id,
        "distinct_matched_disease_gene_pairs": len(assertions_by_pair),
        "matched_pairs_with_multiple_classifications": len(classification_conflicts),
        "classification_conflicts": classification_conflicts,
        "limitations": [
            "Matches use exact disease identifier equality only; disease names are not used.",
            "Original disease identifiers are preferred when present in the HPOA index.",
            "A disease-gene assertion does not establish that a patient's variant is causal.",
            "Multiple classifications are reported, not reconciled into one conclusion.",
            "Match counts measure integration coverage, not diagnostic performance.",
        ],
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"HPOA disease IDs: {len(indexed_ids):,}")
    print(f"GenCC assertions imported: {len(gencc_result.assertions):,}")
    print(f"Matched assertions: {adapter_result.matched_assertions:,}")
    print(f"Unmatched assertions: {adapter_result.unmatched_assertions:,}")
    print(f"Matches by original ID: {adapter_result.matched_by_original_id:,}")
    print(f"Matches by canonical ID: {adapter_result.matched_by_canonical_id:,}")
    print(f"Distinct matched disease-gene pairs: {len(assertions_by_pair):,}")
    print(f"Pairs with multiple classifications: {len(classification_conflicts):,}")
    print(f"Audit written to: {output_path.relative_to(ROOT)}")

    return 1 if gencc_result.row_errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
