from __future__ import annotations

import hashlib
import json
from pathlib import Path

from rare_disease_ai.case_validation import PatientCase
from rare_disease_ai.disease_phenotypes import build_disease_phenotype_index
from rare_disease_ai.gencc_import import read_gencc_csv
from rare_disease_ai.go_annotation_import import read_gaf_annotations
from rare_disease_ai.go_ontology import read_go_term_labels
from rare_disease_ai.hpo_annotations import read_hpo_annotations
from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology
from rare_disease_ai.investigation_workflow import (
    investigation_report_to_dict,
    run_investigation,
)
from rare_disease_ai.variant_prioritization import VariantRecord

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data"
REPORT = ROOT / "reports" / "real_data_integration_smoke_test.json"


def load_metadata(relative_path: str) -> dict:
    return json.loads((ROOT / relative_path).read_text(encoding="utf-8-sig"))


def verify_sha256(relative_path: str, expected: str) -> str:
    path = ROOT / relative_path
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest.lower() != expected.lower():
        raise ValueError(
            f"SHA-256 mismatch for {relative_path}: "
            f"expected {expected}, got {digest}"
        )
    return digest


def main() -> None:
    hpo_meta = load_metadata("data/metadata/hpo_release.json")
    hpoa_meta = load_metadata("data/metadata/hpoa_release.json")
    gencc_meta = load_metadata("data/metadata/gencc_release.json")
    go_meta = load_metadata("data/metadata/go_release.json")

    checksums = {
        "hpo_ontology": verify_sha256(
            "data/raw/hp-base.obo", hpo_meta["sha256"]
        ),
        "hpo_annotations": verify_sha256(
            "data/external/phenotype.hpoa", hpoa_meta["sha256"]
        ),
        "gencc": verify_sha256(
            gencc_meta["local_file"], gencc_meta["sha256"]
        ),
        "go_ontology": verify_sha256(
            go_meta["ontology_file"], go_meta["ontology_sha256"]
        ),
        "go_annotations": verify_sha256(
            go_meta["annotation_file"], go_meta["annotation_sha256"]
        ),
    }
    print("All five reference-data checksums verified.")

    ontology = HPOOntology.from_obo(RAW / "raw" / "hp-base.obo")
    hpo_graph = HPOGraph(ontology)

    hpo_annotations = read_hpo_annotations(
        RAW / "external" / "phenotype.hpoa"
    )
    disease_index = build_disease_phenotype_index(hpo_annotations)
    print(f"HPO annotations loaded: {len(hpo_annotations):,}")
    print(f"Diseases indexed: {len(disease_index):,}")

    gencc_result = read_gencc_csv(
        ROOT / gencc_meta["local_file"]
    )
    print(
        f"GenCC rows: {gencc_result.rows_read:,}; "
        f"valid assertions: {len(gencc_result.assertions):,}; "
        f"row errors: {len(gencc_result.row_errors):,}"
    )

    go_labels = read_go_term_labels(
        ROOT / go_meta["ontology_file"]
    )
    go_result = read_gaf_annotations(
        ROOT / go_meta["annotation_file"],
        source_release=(
            f"GO archive {go_meta['ontology_archive_date']}; "
            f"ontology {go_meta['ontology_data_version']}"
        ),
        source_url=go_meta["annotation_download_url"],
        term_labels=go_labels,
    )
    print(
        f"GO annotation rows: {go_result.rows_read:,}; "
        f"imported: {len(go_result.annotations):,}; "
        f"filtered by taxon: {go_result.rows_filtered_by_taxon:,}; "
        f"errors: {len(go_result.row_errors):,}"
    )

    go_genes = {item.gene_symbol.casefold() for item in go_result.annotations}

    # Pick the first compatible, real GenCC assertion deterministically.
    # The case and variant created below remain synthetic test inputs.
    selected = None
    selected_disease_id = None
    selected_hpo_id = None

    for assertion in gencc_result.assertions:
        possible_ids = (
            assertion.disease_original_curie,
            assertion.disease_curie,
        )
        disease_id = next(
            (value for value in possible_ids if value and value in disease_index),
            None,
        )
        if not disease_id:
            continue
        if assertion.gene_symbol.casefold() not in go_genes:
            continue

        profile = disease_index[disease_id]

        # Prefer a clinically descriptive phenotype over inheritance-mode
        # annotations. Keep selection deterministic and require an HPO label.
        available_hpo_ids = [
            value
            for value in sorted(profile.positive_hpo_ids)
            if ontology.get_term(value) is not None
        ]
        clinical_hpo_ids = [
            value
            for value in available_hpo_ids
            if "inheritance" not in ontology.get_term(value).name.casefold()
            and value != "HP:0000006"
        ]
        preferred_hpo_ids = (
            "HP:0001363",  # Craniosynostosis
            "HP:0001263",  # Global developmental delay
            "HP:0001252",  # Hypotonia
        )
        hpo_id = next(
            (
                value
                for value in preferred_hpo_ids
                if value in clinical_hpo_ids
            ),
            next(iter(clinical_hpo_ids), None),
        )

        if hpo_id:
            selected = assertion
            selected_disease_id = disease_id
            selected_hpo_id = hpo_id
            break

    if selected is None or selected_hpo_id is None:
        raise RuntimeError(
            "No compatible GenCC assertion, HPOA disease phenotype, "
            "and GO-annotated gene were found. Check data releases and IDs."
        )

    term = ontology.get_term(selected_hpo_id)
    phenotype_label = term.name
    print(
        "Selected public-data fixture: "
        f"disease ID={selected_disease_id}, "
        f"HPO={selected_hpo_id} ({phenotype_label}), "
        f"gene={selected.gene_symbol}"
    )

    synthetic_case = PatientCase(
        case_id="SYNTHETIC-INTEGRATION-SMOKE-001",
        phenotype_terms=(phenotype_label,),
    )
    synthetic_variant = VariantRecord(
        chromosome="1",
        position=12345,
        reference="A",
        alternate="T",
        gene_symbol=selected.gene_symbol,
        consequence="missense_variant",
        zygosity="heterozygous",
        allele_frequency=0.001,
    )

    result = run_investigation(
        case=synthetic_case,
        ontology=ontology,
        disease_index=disease_index,
        hpo_graph=hpo_graph,
        variants=(synthetic_variant,),
        gencc_assertions=gencc_result.assertions,
        mechanism_annotations=go_result.annotations,
        top_k=20,
    )
    report = investigation_report_to_dict(result)

    output = {
        "report_schema_version": "1.0",
        "status": "completed",
        "test_type": "real_reference_data_synthetic_input_smoke_test",
        "research_only": True,
        "clinical_validation": False,
        "benchmark_evaluation": False,
        "synthetic_inputs": {
            "case_id": synthetic_case.case_id,
            "phenotype_term": phenotype_label,
            "variant": {
                "gene_symbol": selected.gene_symbol,
                "consequence": "missense_variant",
                "position_is_synthetic": True,
                "allele_frequency_is_synthetic": True,
            },
            "fixture_selection_note": (
                "The synthetic phenotype was selected from the chosen "
                "disease's HPOA profile to exercise integration. This is "
                "target-informed test construction, not independent "
                "diagnostic evaluation."
            ),
        },
        "reference_data": {
            "hpo_release": hpo_meta["release"],
            "hpoa_release": hpoa_meta["release"],
            "gencc_retrieved_date": gencc_meta["retrieved_date"],
            "go_archive_date": go_meta["ontology_archive_date"],
            "go_ontology_data_version": go_meta["ontology_data_version"],
            "go_annotation_generated_at": go_meta["annotation_generated_at"],
            "sha256_verified": checksums,
            "counts": {
                "hpo_annotations": len(hpo_annotations),
                "diseases_indexed": len(disease_index),
                "gencc_rows_read": gencc_result.rows_read,
                "gencc_assertions_imported": len(gencc_result.assertions),
                "gencc_row_errors": len(gencc_result.row_errors),
                "go_rows_read": go_result.rows_read,
                "go_annotations_imported": len(go_result.annotations),
                "go_rows_filtered_by_taxon": go_result.rows_filtered_by_taxon,
                "go_row_errors": len(go_result.row_errors),
            },
        },
        "workflow_result": report,
        "limitations": [
            "Synthetic inputs are not patient data.",
            "The phenotype was selected using the target disease profile.",
            "This is a smoke test, not an independent benchmark.",
            "Gene-disease assertions and GO annotations do not establish causality.",
            "Variant scores are heuristic and are not clinical pathogenicity classifications.",
            "The system is not validated for clinical decision-making.",
        ],
    }

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(
        json.dumps(output, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"\nReport written: {REPORT.relative_to(ROOT)}")
    print("Research-only smoke test completed.")


if __name__ == "__main__":
    main()
