from rare_disease_ai.phenotype_ranking import DiseaseRankingResult
from rare_disease_ai.phenotype_variant_integration import (
    GeneDiseaseEvidence,
    integrate_ranked_evidence,
    integration_report_to_dict,
)
from rare_disease_ai.variant_prioritization import (
    VariantRankingResult,
    VariantRecord,
)


def make_disease(
    disease_id: str,
    name: str,
    score: float,
) -> DiseaseRankingResult:
    return DiseaseRankingResult(
        database_id=disease_id,
        disease_name=name,
        score=score,
        semantic_score=None,
        exact_matches=(),
        related_matches=(),
        unmatched_patient_hpo_ids=(),
        negative_matches=(),
    )


def make_variant(
    gene: str,
    consequence: str = "frameshift_variant",
) -> VariantRankingResult:
    variant = VariantRecord(
        chromosome="1",
        position=12345,
        reference="A",
        alternate="T",
        gene_symbol=gene,
        consequence=consequence,
        zygosity="heterozygous",
        allele_frequency=0.001,
    )
    return VariantRankingResult(
        variant=variant,
        score=6.0,
        reasons=("Synthetic test evidence",),
    )


def test_links_variant_only_with_explicit_gene_disease_evidence():
    diseases = (
        make_disease("OMIM:100001", "Disease A", 1.0),
        make_disease("OMIM:100002", "Disease B", 0.5),
    )
    variants = (make_variant("GENE_A"), make_variant("GENE_B"))
    evidence = (
        GeneDiseaseEvidence(
            disease_id="OMIM:100001",
            gene_symbol="GENE_A",
            source="synthetic-test-source",
            evidence_id="TEST-001",
        ),
    )

    results = integrate_ranked_evidence(diseases, variants, evidence)

    assert results[0].disease_id == "OMIM:100001"
    assert results[0].phenotype_rank == 1
    assert len(results[0].matched_variants) == 1
    assert results[0].matched_variants[0].gene_symbol == "GENE_A"
    assert results[0].matched_variants[0].evidence_id == "TEST-001"
    assert results[1].matched_variants == ()


def test_keeps_phenotype_and_variant_scores_separate():
    diseases = (make_disease("OMIM:100001", "Disease A", 0.75),)
    variants = (make_variant("GENE_A"),)
    evidence = (
        GeneDiseaseEvidence(
            disease_id="OMIM:100001",
            gene_symbol="GENE_A",
            source="synthetic-test-source",
            evidence_id="TEST-001",
        ),
    )

    result = integrate_ranked_evidence(diseases, variants, evidence)[0]

    assert result.phenotype_score == 0.75
    assert result.matched_variants[0].variant_score == 6.0


def test_rejects_gene_disease_evidence_without_provenance():
    diseases = (make_disease("OMIM:100001", "Disease A", 1.0),)
    invalid_evidence = (
        GeneDiseaseEvidence(
            disease_id="OMIM:100001",
            gene_symbol="GENE_A",
            source="",
            evidence_id="TEST-001",
        ),
    )

    try:
        integrate_ranked_evidence(diseases, (make_variant("GENE_A"),), invalid_evidence)
    except ValueError as exc:
        assert "source" in str(exc)
    else:
        raise AssertionError("Expected missing evidence provenance to be rejected")


def test_report_is_marked_research_only():
    results = integrate_ranked_evidence(
        (make_disease("OMIM:100001", "Disease A", 1.0),),
        (),
        (),
    )

    report = integration_report_to_dict(results)

    assert report["research_only"] is True
    assert "not clinical probabilities" in report["score_policy"]
    assert report["diseases"][0]["matched_variants"] == []
