from rare_disease_ai.case_validation import PatientCase
from rare_disease_ai.disease_phenotypes import DiseasePhenotypeProfile
from rare_disease_ai.evidence_retrieval_service import EvidenceRetrievalService
from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology, HPOTerm
from rare_disease_ai.investigation_workflow import (
    investigation_report_to_dict,
    run_investigation,
)
from rare_disease_ai.variant_prioritization import CONSEQUENCE_SCORES, VariantRecord


def test_workflow_retrieves_brca1_evidence_from_sqlite():
    ontology = HPOOntology(
        {
            "HP:0001250": HPOTerm(
                hpo_id="HP:0001250",
                name="Seizures",
            )
        }
    )
    graph = HPOGraph(ontology)

    disease_index = {
        "OMIM:123456": DiseasePhenotypeProfile(
            database_id="OMIM:123456",
            disease_name="Synthetic test disorder",
            positive_hpo_ids=frozenset({"HP:0001250"}),
            negative_hpo_ids=frozenset(),
            annotations=(),
        )
    }

    case = PatientCase(
        case_id="sqlite-integration-test",
        phenotype_terms=("Seizures",),
    )

    variant = VariantRecord(
        chromosome="1",
        position=12345,
        reference="A",
        alternate="T",
        gene_symbol="BRCA1",
        consequence=next(iter(CONSEQUENCE_SCORES)),
        zygosity="heterozygous",
        allele_frequency=0.001,
    )

    service = EvidenceRetrievalService(index_path="data/processed/evidence_search.sqlite3")

    result = run_investigation(
        case=case,
        ontology=ontology,
        disease_index=disease_index,
        hpo_graph=graph,
        variants=(variant,),
        gencc_assertions=(),
        mechanism_annotations=(),
        evidence_retrieval_service=service,
        evidence_query="DNA repair",
        evidence_top_k=5,
        evidence_entity_ids=("BRCA1",),
        evidence_types=("go_mechanism_annotation",),
    )

    report = investigation_report_to_dict(result)
    retrieved = report["retrieved_evidence"]

    assert retrieved["research_only"] is True
    assert retrieved["query"] == "DNA repair"
    assert retrieved["returned_count"] == 5

    for hit in retrieved["hits"]:
        assert "BRCA1" in hit["entity_ids"]
        assert hit["evidence_type"] == "go_mechanism_annotation"
        assert hit["source_name"]
        assert hit["source_url"].startswith("https://")
        assert hit["evidence_id"]
        assert hit["matched_terms"]

    print("SQLite-backed workflow integration passed.")
    print("Retrieved hits:", retrieved["returned_count"])
    print("First hit:", retrieved["hits"][0]["title"])
