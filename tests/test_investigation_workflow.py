import json

from rare_disease_ai.biomedical_evidence_retrieval import EvidenceDocument
from rare_disease_ai.case_validation import PatientCase
from rare_disease_ai.disease_phenotypes import DiseasePhenotypeProfile
from rare_disease_ai.evidence_retrieval_service import EvidenceRetrievalService
from rare_disease_ai.gencc_import import GenCCAssertion
from rare_disease_ai.hpo_graph import HPOGraph
from rare_disease_ai.hpo_ontology import HPOOntology, HPOTerm
from rare_disease_ai.investigation_workflow import (
    investigation_report_to_dict,
    run_investigation,
)
from rare_disease_ai.mechanism_evidence import MechanismEvidence
from rare_disease_ai.variant_prioritization import (
    CONSEQUENCE_SCORES,
    VariantRecord,
)


def test_end_to_end_investigation_report_is_json_serializable():
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
        case_id="synthetic-case-001",
        phenotype_terms=("Seizures", "Unknown synthetic symptom"),
    )

    variant = VariantRecord(
        chromosome="1",
        position=12345,
        reference="A",
        alternate="T",
        gene_symbol="GENE1",
        consequence=next(iter(CONSEQUENCE_SCORES)),
        zygosity="heterozygous",
        allele_frequency=0.001,
    )

    assertion = GenCCAssertion(
        submission_id="synthetic-submission",
        version_number=1,
        gene_curie="HGNC:1",
        gene_symbol="GENE1",
        disease_curie="OMIM:123456",
        disease_title="Synthetic test disorder",
        disease_original_curie="OMIM:123456",
        classification_curie="GENCC:100001",
        classification_title="Definitive",
        moi_curie="GENCC:100001",
        moi_title="Autosomal dominant",
        submitter_curie="GENCC:SUBMITTER1",
        submitter_title="Synthetic submitter",
        submitted_as_date="2026-01-01",
        public_report_url="https://example.org/report",
        pmids=("12345678",),
        assertion_criteria_url="https://example.org/criteria",
    )

    mechanism = MechanismEvidence(
        gene_symbol="GENE1",
        annotation_id="GO:0000001",
        annotation_label="Synthetic biological process",
        annotation_type="biological_process",
        source_name="Gene Ontology",
        source_release="synthetic-test",
        source_url="https://geneontology.org/",
        evidence_code="IDA",
        evidence_label="Inferred from Direct Assay",
        reference_ids=("PMID:12345678",),
        qualifier=None,
        taxon_id="taxon:9606",
    )

    # Create synthetic evidence for testing the optional retrieval stage.
    evidence_document = EvidenceDocument(
        evidence_id="SYNTHETIC:GO:1",
        title="Synthetic biological process",
        text="Synthetic evidence describes a biological process",
        source_name="Synthetic evidence source",
        source_release="test-1",
        source_url="https://example.org/evidence",
        evidence_type="go_mechanism",
        entity_ids=("GENE1",),
        reference_ids=("PMID:12345678",),
    )

    retrieval_service = EvidenceRetrievalService(
        documents=(evidence_document,),
    )

    result = run_investigation(
        case=case,
        ontology=ontology,
        disease_index=disease_index,
        hpo_graph=graph,
        variants=(variant,),
        gencc_assertions=(assertion,),
        mechanism_annotations=(mechanism,),
        top_k=10,
        evidence_retrieval_service=retrieval_service,
        evidence_query="synthetic biological process",
        evidence_entity_ids=("GENE1",),
        evidence_types=("go_mechanism",),
    )

    report = investigation_report_to_dict(result)
    encoded = json.dumps(report)
    decoded = json.loads(encoded)

    # Existing workflow assertions.
    assert decoded["research_only"] is True
    assert decoded["case_id"] == "synthetic-case-001"
    assert decoded["normalization"]["counts_by_status"]["resolved"] == 1
    assert decoded["normalization"]["counts_by_status"]["unresolved"] == 1
    assert decoded["phenotype_rankings"]
    assert decoded["variant_rankings"]
    assert decoded["gencc_audit"]["matched_assertions"] == 1
    assert decoded["integrated_evidence"]["diseases"]
    assert decoded["mechanisms"]["results"]

    # New evidence retrieval assertions.
    retrieved = decoded["retrieved_evidence"]

    assert retrieved["research_only"] is True
    assert retrieved["query"] == "synthetic biological process"
    assert retrieved["total_documents"] == 1
    assert retrieved["eligible_documents"] == 1
    assert retrieved["returned_count"] == 1
    assert "not diagnostic confidence" in retrieved["score_interpretation"]

    hit = retrieved["hits"][0]

    assert hit["evidence_id"] == "SYNTHETIC:GO:1"
    assert hit["title"] == "Synthetic biological process"
    assert hit["evidence_type"] == "go_mechanism"
    assert hit["entity_ids"] == ["GENE1"]
    assert hit["reference_ids"] == ["PMID:12345678"]
    assert hit["source_name"] == "Synthetic evidence source"
    assert hit["source_release"] == "test-1"
    assert hit["source_url"] == "https://example.org/evidence"
    assert hit["matched_terms"]
    assert hit["score"] > 0


def test_no_resolved_phenotypes_produces_no_disease_rankings():
    ontology = HPOOntology(
        {
            "HP:0001250": HPOTerm(
                hpo_id="HP:0001250",
                name="Seizures",
            )
        }
    )

    result = run_investigation(
        case=PatientCase(
            case_id="synthetic-unresolved-case",
            phenotype_terms=("Unknown synthetic symptom",),
        ),
        ontology=ontology,
        disease_index={},
        hpo_graph=HPOGraph(ontology),
        variants=(),
        gencc_assertions=(),
        mechanism_annotations=(),
    )

    assert result.disease_rankings == ()
    assert result.integrated_diseases == ()
    assert result.retrieved_evidence is None

    report = investigation_report_to_dict(result)

    assert report["retrieved_evidence"] is None
